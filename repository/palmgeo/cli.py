"""Command line interface.

The pipeline is split into resumable stages. Every stage writes a ``.done``
marker per unit of work, so re-running after a network failure only redoes what
is missing -- which matters when a stage moves tens of gigabytes.

    CFG=static_driver/<run>/config/<run>.yml

    palmgeo sources   $CFG    # what data exists, under what licence
    palmgeo plan      $CFG    # what would be downloaded, no I/O
    palmgeo fetch     $CFG    # optional: pull every source first,
                              #   so the rest can run offline
    palmgeo vectors   $CFG    # cadastre, roads, clipped to the AOI
    palmgeo rasters   $CFG    # terrain, nDSM, NDVI, per tile
    palmgeo buildings $CFG    # heights joined onto footprints
    palmgeo trees     $CFG    # single trees from CHM + orthophoto
    palmgeo compose   $CFG    # land use x land cover -> PALM types
    palmgeo domains   $CFG    # aggregate onto the PALM grids
    palmgeo csd       $CFG    # write the palm_csd YAML
    palmgeo inspect   $CFG    # QGIS products
    palmgeo all       $CFG    # all of the above, in order
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from palmgeo import __version__
from palmgeo.cache import Cache
from palmgeo.config import Config, ConfigError, load
from palmgeo.grid import alignment_report

logger = logging.getLogger("palmgeo")

STAGES = [
    "sources", "plan", "vectors", "rasters", "buildings", "trees", "compose",
    "domains", "csd", "inspect",
]

#: Stages that move real volume, and therefore must not start before the data
#: sources and their licences have been looked at. `plan` is exempt on purpose:
#: it is a dry run, and telling somebody to approve their licences before they
#: are allowed to find out how big the download is would be the wrong order.
GATED_STAGES = ("vectors", "rasters")


def _setup_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)-18s %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stdout,
    )
    logging.getLogger("rasterio").setLevel(logging.WARNING)
    logging.getLogger("pyogrio").setLevel(logging.WARNING)
    logging.getLogger("fiona").setLevel(logging.WARNING)


def _cache(cfg: Config) -> Cache:
    return Cache(
        cfg.cachedir,
        budget_bytes=int(float(cfg.runtime.get("cache_budget_gb", 40)) * 1e9),
        keep=bool(cfg.runtime.get("keep_downloads", False)),
        retries=int(cfg.runtime.get("retries", 4)),
        timeout=float(cfg.runtime.get("timeout_s", 180)),
        root=cfg.root,
        offline=bool(cfg.runtime.get("offline", False)),
    )


def _state_path(cfg: Config) -> Path:
    return cfg.workdir / "state.json"


def _save_state(cfg: Config, stage: str, payload: Dict[str, Any]) -> None:
    path = _state_path(cfg)
    path.parent.mkdir(parents=True, exist_ok=True)
    state = json.loads(path.read_text()) if path.exists() else {}
    state[stage] = payload
    path.write_text(json.dumps(state, indent=2, default=str), encoding="utf-8")


def _load_state(cfg: Config) -> Dict[str, Any]:
    path = _state_path(cfg)
    return json.loads(path.read_text()) if path.exists() else {}


def _vector_paths(cfg: Config) -> Dict[str, Path]:
    vdir = cfg.workdir / "vector"
    names = [
        # raw cadastre extracts (written by the `vectors` stage)
        "cadastre_landuse", "cadastre_buildings", "cadastre_structures",
        # assembled products (written by the `buildings` stage)
        "buildings", "bridges",
        # topographic products
        "street_type", "street_crossings",
    ]
    return {n: vdir / f"{n}.gpkg" for n in names if (vdir / f"{n}.gpkg").exists()}


def _ndsm_sources(cfg: Config) -> List[Path]:
    """One nDSM mosaic per work tier, finest first.

    A run with several tiers only covers the fine area at fine resolution, so
    buildings outside it have to be sampled from a coarser tier rather than
    falling through to the 3D model. Callers sample the list in order and stop
    at the first tier that yields a value.
    """
    from palmgeo.rasters import build_vrt

    sources: List[Path] = []
    for tier in cfg.tiers:  # tiers are sorted finest first
        parts = sorted((cfg.tier_dir(tier) / "tiles").glob("*/ndsm.tif"))
        if not parts:
            continue
        vrt_dir = cfg.tier_dir(tier) / "vrt"
        vrt_dir.mkdir(parents=True, exist_ok=True)
        sources.append(build_vrt(parts, vrt_dir / "ndsm_all.vrt"))
    return sources


#: Bytes per pixel of a source's raster payload, by role (float32 / 4x uint8).
_ROLE_BYTES = {"dtm": 4.0, "dsm": 4.0, "ortho": 4.0}

#: COG internal block edge. Reading part of a block still transfers the block,
#: which dominates the cost of a heavily decimated read and is why a flat
#: per-tile figure is wrong by an order of magnitude at coarse work resolutions.
_BLOCK = 512

#: Slack over the pure block arithmetic for HTTP headers, IFD reads and the
#: warper occasionally choosing a finer overview than strictly needed.
_STREAM_SLACK = 1.6


def _streamed_bytes(native: float, work: float, per_px: float, tile: float) -> float:
    """Bytes moved to read one work tile from a COG at ``work`` resolution."""
    # GDAL serves the request from the finest overview no coarser than needed;
    # overview levels are powers of two.
    level = max(1, 2 ** int(math.floor(math.log2(max(work / native, 1.0)))))
    eff = native * level
    span = tile / eff                                   # pixels actually wanted
    blocks = math.ceil(span / _BLOCK) ** 2              # ... rounded to blocks
    return blocks * _BLOCK**2 * per_px * _STREAM_SLACK


def _transfer_estimate(cfg: Config, summaries: List[Dict[str, Any]]) -> float:
    """Approximate GB pulled over the network, honouring stream vs download."""
    from palmgeo.grid import TILE_SIZE
    from palmgeo.steps.tiles import _access_mode

    total = 0.0
    for tier, summary in zip(cfg.tiers, summaries):
        n = int(summary["tiles"])
        for role in ("dtm", "dsm", "ortho"):
            try:
                _, src = cfg.source_for_role(role)
            except Exception:
                continue
            native = float(src.get("resolution", 1.0))
            per_px = _ROLE_BYTES.get(role, 4.0)
            epochs = max(1, len(summary["ortho"])) if role == "ortho" else 1
            if _access_mode(cfg, native, tier.resolution) == "stream":
                per_tile = _streamed_bytes(native, tier.resolution, per_px, TILE_SIZE)
            else:
                # the whole source tile, compressed; overviews add ~1/3
                src_tile = float(src.get("tile_size", TILE_SIZE))
                px = (src_tile / native) ** 2
                share = (TILE_SIZE / src_tile) ** 2
                per_tile = px * per_px * 1.33 * share
            total += n * epochs * per_tile
    return total / 1e9


# ---------------------------------------------------------------------------
# stages
# ---------------------------------------------------------------------------
def stage_plan(cfg: Config) -> Dict[str, Any]:
    """Report what the run would do, without downloading raster data."""
    from palmgeo.steps.tiles import plan_tiles

    cache = _cache(cfg)
    print(f"\nProject '{cfg.get('project.name')}'  --  {cfg.crs}, profile "
          f"'{cfg.runtime.get('profile')}'\n")
    print("Domains")
    for d in cfg.domains:
        tier = cfg.tier_of(d)
        exact, message = alignment_report(d, tier.resolution)
        parent = f", parent {d.domain_parent}" if d.domain_parent else ""
        print(f"  {d.describe()}{parent}")
        print(f"      work tier {tier.name}, vegetation_mode "
              f"{d.options.get('vegetation_mode')}, {'exact' if exact else 'APPROXIMATE'}"
              f" aggregation")
        if not exact:
            print(f"      ! {message}")

    print("\nWork tiers")
    total_tiles = 0
    summaries: List[Dict[str, Any]] = []
    for tier in cfg.tiers:
        plans, summary = plan_tiles(cfg, cache, tier)
        summaries.append(summary)
        total_tiles += len(plans)
        print(f"  {tier.name}: {len(plans)} tile(s), {summary['area_km2']:.0f} km2, "
              f"domains {[d.name for d in tier.domains]}")
        if summary["tiles_without_terrain"]:
            print(f"      ! {summary['tiles_without_terrain']} tile(s) without terrain data")
        for report in summary["ortho"]:
            uniform = report.get("uniform")
            note = (
                f"single campaign around {uniform['anchor']} "
                f"(+/-{uniform['spread_days']} d)"
                if uniform else "per-tile epoch choice (expect NDVI steps on the ortho tile lattice)"
            )
            print(f"      ortho '{report['role']}': {note}")
            if report.get("relaxed"):
                print(f"      ! {report['role']}: {report['relaxed']}")
            if report.get("missing_tiles"):
                print(f"      ! {report['role']}: {report['missing_tiles']} tile(s) uncovered")

    # Volume estimate, so the disk and the time budget can be sized before
    # starting. It has to follow the access decision: a source read far below
    # its native resolution is streamed from a COG overview, and then the
    # transfer scales with the *work* grid, not with the file. Over a 2600 km2
    # root that is the difference between ~190 GB and ~10 GB, so a flat
    # per-tile figure would be useless here.
    est_gb = _transfer_estimate(cfg, summaries)
    print(f"\nEstimated source volume to transfer: ~{est_gb:.0f} GB "
          f"across {total_tiles} work tile(s).")
    print(f"Cache budget: {float(cfg.runtime.get('cache_budget_gb', 40)):.0f} GB, "
          f"keep_downloads={cfg.runtime.get('keep_downloads')}, "
          f"access={cfg.runtime.get('access')}\n")
    return {"tiers": summaries, "estimated_transfer_gb": round(est_gb, 1)}


def stage_vectors(cfg: Config) -> Dict[str, Any]:
    from palmgeo.steps import fetch, vectors

    fetch.ensure_available(cfg, _cache(cfg), "vectors")

    cache = _cache(cfg)
    cadastre = vectors.fetch_cadastre(cfg, cache)
    roads = vectors.fetch_roads(cfg, cache)
    logger.info("%s", cache.summary())
    return {
        "cadastre": {k: str(v) for k, v in cadastre.items()},
        "roads": str(roads) if roads else None,
    }


def stage_rasters(cfg: Config) -> Dict[str, Any]:
    from palmgeo.steps import fetch

    fetch.ensure_available(cfg, _cache(cfg), "rasters")
    from palmgeo.steps import tiles

    cache = _cache(cfg)
    return {"tiers": [tiles.run(cfg, cache, tier) for tier in cfg.tiers]}


def stage_buildings(cfg: Config) -> Dict[str, Any]:
    from palmgeo.steps import fetch

    fetch.ensure_available(cfg, _cache(cfg), "buildings")
    from palmgeo.steps import vectors

    cache = _cache(cfg)
    paths = _vector_paths(cfg)
    cadastre = {
        role: paths[f"cadastre_{role}"]
        for role in ("buildings", "structures")
        if f"cadastre_{role}" in paths
    }
    if "buildings" not in cadastre:
        raise RuntimeError("No cadastre buildings layer. Run the `vectors` stage first.")
    ndsm = _ndsm_sources(cfg)
    if not ndsm:
        logger.warning("No nDSM available -- run the `rasters` stage first")
    heights = vectors.fetch_building_heights(cfg, cache)
    buildings = vectors.build_buildings(cfg, cadastre, heights, ndsm)
    bridges = vectors.build_bridges(cfg, cadastre, ndsm)
    logger.info("%s", cache.summary())
    return {
        "buildings": str(buildings) if buildings else None,
        "bridges": str(bridges) if bridges else None,
    }


def stage_trees(cfg: Config) -> Dict[str, Any]:
    from palmgeo.steps import trees

    return trees.run(cfg)


def stage_compose(cfg: Config) -> Dict[str, Any]:
    from palmgeo.steps import compose

    vecs = _vector_paths(cfg)
    return {"tiers": [compose.run(cfg, tier, vecs) for tier in cfg.tiers]}


def stage_domains(cfg: Config) -> Dict[str, Any]:
    from palmgeo.steps import domains

    return domains.run(cfg)


def stage_csd(cfg: Config) -> Dict[str, Any]:
    from palmgeo.steps import csd

    return {"config": str(csd.build(cfg))}


def stage_inspect(cfg: Config) -> Dict[str, Any]:
    from palmgeo.steps import inspect

    return inspect.run(cfg)


def stage_sources(cfg: Config, approve: bool = False,
                  probe_alternatives: bool = True) -> Dict[str, Any]:
    """Report what data exists for this area and what may be done with it.

    Running this stage does not approve anything. Reading a licence table and
    agreeing to it are different acts, and an approval that happened as a side
    effect of looking would be worthless -- so the approval is a separate,
    explicit `--approve`, which is also what makes it meaningful when an agent
    passes it on someone's behalf.
    """
    from palmgeo.steps import inventory

    cache = _cache(cfg)
    reports = inventory.scan(cfg, cache, probe_alternatives=probe_alternatives)
    print(inventory.render(cfg, reports))

    path = inventory.write_markdown(cfg, reports)
    print(f"  Written to {path.relative_to(cfg.root)}")

    digest = inventory.fingerprint(reports)
    blocking = [r.role for r in reports if r.required and r.chosen is None]
    if blocking:
        raise RuntimeError(
            f"No source for the required role(s) {', '.join(blocking)}. "
            f"A driver cannot be built without them; see the report above for "
            f"what each one becomes."
        )

    # Re-running the report must not revoke an approval that still holds. It
    # would otherwise be impossible to run `all` after approving: `all` starts
    # with this stage, which would clear the approval, and the gate two stages
    # later would then block on it. An approval is only spent when the thing it
    # covers changes, which is exactly what the fingerprint comparison says.
    previous = _load_state(cfg).get("sources") or {}
    still_valid = bool(previous.get("approved")) and previous.get("fingerprint") == digest
    approved = bool(approve) or still_valid

    if approve:
        print(f"  Approved: {digest}\n")
    elif still_valid:
        print(f"  Already approved: {digest} (unchanged since "
              f"{previous.get('approved_at')})\n")
    else:
        print(f"  NOT approved. To accept these sources and their licences:\n"
              f"    python -m palmgeo.cli sources {cfg.path} --approve\n")

    return {
        "fingerprint": digest,
        "approved": approved,
        "approved_at": _utcnow() if approve else (previous.get("approved_at") if still_valid else None),
        "roles": {
            r.role: {
                "source": r.chosen.name if r.chosen else None,
                "licences": r.chosen.licence_ids if r.chosen else [],
                "coverage": r.chosen.coverage if r.chosen else None,
            }
            for r in reports
        },
        "licences": inventory.used_licences(cfg, reports),
        "obligations": {
            k: sorted(set(v))
            for k, v in inventory.combined_obligations(cfg, reports).items() if v
        },
    }


def cmd_fetch(cfg: Config, check_only: bool = False) -> Dict[str, Any]:
    """Materialise every source the run needs, so the rest can run offline."""
    from palmgeo.steps import fetch

    cache = _cache(cfg)
    if cache.offline and not check_only:
        raise RuntimeError(
            "runtime.offline is set, so `fetch` has nothing it may do. Unset it "
            "for the acquisition run, then set it again for the offline one."
        )
    return fetch.run(cfg, cache, check_only=check_only)


def _utcnow() -> str:
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def require_source_approval(cfg: Config) -> None:
    """Refuse a bulk-download stage until the sources have been approved.

    The check re-derives the fingerprint from the configuration as it stands
    now, so an approval does not survive a change of source or a change of
    licence -- which is the whole point of gating on a fingerprint rather than
    on a flag file.

    It deliberately re-runs the scan, which costs the index reads again. Those
    are cached; being certain that the approval matches what is about to be
    downloaded is worth more than the seconds.
    """
    from palmgeo.steps import inventory

    state = _load_state(cfg).get("sources") or {}
    hint = f"    python -m palmgeo.cli sources {cfg.path} --approve"
    if not state.get("approved"):
        raise RuntimeError(
            "The data sources for this run have not been approved. Before "
            "downloading, look at what is available and under which licence:\n"
            f"{hint}"
        )
    current = inventory.fingerprint(inventory.scan(cfg, _cache(cfg), probe_alternatives=False))
    if current != state.get("fingerprint"):
        raise RuntimeError(
            f"The approved sources ({state.get('fingerprint')}) are not the ones "
            f"this run would now use ({current}). A source or a licence has "
            f"changed since the approval. Look again:\n{hint}"
        )


STAGE_FUNCS = {
    "sources": stage_sources,
    "plan": stage_plan,
    "vectors": stage_vectors,
    "rasters": stage_rasters,
    "buildings": stage_buildings,
    "trees": stage_trees,
    "compose": stage_compose,
    "domains": stage_domains,
    "csd": stage_csd,
    "inspect": stage_inspect,
}


# ---------------------------------------------------------------------------
# auxiliary commands
# ---------------------------------------------------------------------------
def cmd_calibrate(cfg: Config, sample_tiles: int = 12) -> None:
    """Print an NDVI histogram so the vegetation threshold can be chosen.

    NDVI here comes from 8-bit orthophoto digital numbers, not calibrated
    reflectance, so the threshold is campaign-specific. A bimodal histogram with
    a clear trough is the sign of a usable separation.
    """
    import numpy as np
    import rasterio

    tier = cfg.tiers[0]
    parts = sorted((cfg.tier_dir(tier) / "tiles").glob("*/ndvi_leaf_on.tif"))
    if not parts:
        print("No NDVI tiles yet -- run the `rasters` stage first.")
        return
    step = max(1, len(parts) // sample_tiles)
    values: List[np.ndarray] = []
    for path in parts[::step][:sample_tiles]:
        with rasterio.open(path) as src:
            data = src.read(1, out_shape=(min(src.height, 250), min(src.width, 250)))
        values.append(data[np.isfinite(data)].ravel())
    if not values:
        print("No valid NDVI values found.")
        return
    sample = np.concatenate(values)
    counts, edges = np.histogram(sample, bins=40, range=(-0.4, 0.9))
    peak = counts.max()
    print(f"\nNDVI over {len(parts[::step][:sample_tiles])} sampled tile(s), "
          f"n = {sample.size}\n")
    for count, lo, hi in zip(counts, edges[:-1], edges[1:]):
        bar = "#" * int(60 * count / peak) if peak else ""
        print(f"  {lo:+.3f}..{hi:+.3f} |{bar}")
    for q in (5, 25, 50, 75, 95):
        print(f"  p{q:<2d} = {np.percentile(sample, q):+.3f}")
    print(f"\n  currently configured classification.ndvi_threshold: "
          f"{cfg.get('classification.ndvi_threshold')}\n")


def cmd_report(cfg: Config, what: str) -> None:
    """Print the unmapped source vocabulary so mappings can be completed."""
    if what in {"landuse", "all"}:
        from palmgeo.steps.compose import collect_landuse_keys

        keys, unmapped = collect_landuse_keys(cfg)
        print(f"\nland-use keys in the AOI: {len(keys)}, without an explicit rule: "
              f"{len(unmapped)}")
        for key in unmapped:
            print(f"  {key}")
    if what in {"street_type", "all"}:
        import geopandas as gpd

        path = cfg.workdir / "vector" / "street_type.gpkg"
        if path.exists():
            gdf = gpd.read_file(path, ignore_geometry=True)
            print(f"\nstreet_type distribution over {len(gdf)} segment(s):")
            for value, count in gdf["street_type"].value_counts().items():
                print(f"  {value}: {count}")
        else:
            print("\nNo roads layer; run the `vectors` stage with "
                  "sources.street_type.enabled: true")


def cmd_lonlat(lon: float, lat: float, epsg: int) -> None:
    from pyproj import Transformer

    t = Transformer.from_crs(4326, epsg, always_xy=True)
    x, y = t.transform(lon, lat)
    print(f"lon/lat ({lon}, {lat})  ->  EPSG:{epsg}  origin_x: {x:.1f}  origin_y: {y:.1f}")


def cmd_clean(cfg: Config, what: str) -> None:
    import shutil

    targets = {
        "cache": [cfg.cachedir],
        "tiers": [cfg.workdir / "tiers"],
        "vector": [cfg.workdir / "vector"],
        "work": [cfg.workdir],
        "out": [cfg.outdir],
        # Composition only, keeping the expensive downloaded/derived tiles.
        # Normally unnecessary -- a changed threshold invalidates the affected
        # tiles by itself -- but it forces a rebuild after a code change.
        "surface": sorted((cfg.workdir / "tiers").glob("*/surface")),
    }
    if what not in targets:
        print(f"clean: choose one of {sorted(targets)}")
        return
    paths = [p for p in targets[what] if p.exists()]
    if not paths:
        print(f"nothing to remove for '{what}'")
        return
    size = sum(f.stat().st_size for p in paths for f in p.rglob("*") if f.is_file())
    for path in paths:
        shutil.rmtree(path)
    print(f"Removed {', '.join(str(p) for p in paths)} ({size / 1e9:.2f} GB)")


# ---------------------------------------------------------------------------
def set_source_in_config(cfg: Config, role: str, name: str) -> bool:
    """Write ``sources.<role>: <name>`` into the run configuration file.

    Edited line by line rather than by re-serialising the YAML, because a run
    configuration is at least half comments -- the reasoning for every choice
    in it -- and a round trip through a plain YAML dumper would delete all of
    it and reorder the rest. The cost is that this only understands the two
    shapes the `sources:` block actually takes.

    Returns True if the file changed.
    """
    import re as _re

    if role not in Config._ROLE_SELECTOR:
        raise ConfigError(f"'{role}' is not a role palmgeo knows: "
                          f"{', '.join(Config._ROLE_SELECTOR)}")
    if name not in (cfg.sources.get("sources") or {}):
        raise ConfigError(f"'{name}' is not a source in {cfg.get('project.sources')}")
    declared = cfg.source(name).get("role")
    if declared != role:
        raise ConfigError(f"Source '{name}' has role '{declared}', not '{role}'")

    lines = cfg.path.read_text(encoding="utf-8").split("\n")
    # `topo` is selected through sources.street_type.source, one level deeper
    # than the others, because street handling carries its own enable switch.
    nested = role == "topo"
    key, indent = ("source", "    ") if nested else (role, "  ")

    start = next((i for i, l in enumerate(lines) if _re.match(r"^sources:\s*$", l)), None)
    if start is None:
        raise ConfigError(f"{cfg.path.name} has no top-level `sources:` block to edit")
    end = next((i for i in range(start + 1, len(lines))
                if lines[i] and not lines[i][0].isspace()), len(lines))

    if nested:
        # An inline mapping -- `street_type: {enabled: true, source: bdlm}` --
        # is valid YAML and common in these configs, but it holds the key on
        # the same line as its parent, where the block-style search below
        # cannot see it. Inserting a second `street_type:` would then produce
        # a duplicate key and silently drop half the setting, so refuse and
        # say exactly what to do instead.
        inline = next((i for i in range(start + 1, end)
                       if _re.match(r"^  street_type:\s*\{", lines[i])), None)
        if inline is not None:
            raise ConfigError(
                f"{cfg.path.name}:{inline + 1} writes street_type as an inline "
                f"mapping, which this editor will not rewrite without risking "
                f"the rest of the line. Expand it to block style:\n"
                f"  street_type:\n    enabled: true\n    source: {name}"
            )
        st = next((i for i in range(start + 1, end)
                   if _re.match(r"^  street_type:\s*$", lines[i])), None)
        if st is None:
            lines.insert(end, "  street_type:")
            lines.insert(end + 1, f"    source: {name}")
            cfg.path.write_text("\n".join(lines), encoding="utf-8")
            return True
        start, end = st, next((i for i in range(st + 1, end)
                               if lines[i][:3].strip() and not lines[i].startswith("    ")), end)

    for i in range(start + 1, end):
        if _re.match(rf"^{indent}{_re.escape(key)}:\s", lines[i]) or \
           _re.match(rf"^{indent}{_re.escape(key)}:\s*$", lines[i]):
            old = lines[i]
            comment = old.split("#", 1)[1] if "#" in old else ""
            lines[i] = f"{indent}{key}: {name}" + (f"  #{comment}" if comment else "")
            if lines[i] == old:
                return False
            cfg.path.write_text("\n".join(lines), encoding="utf-8")
            return True

    lines.insert(end, f"{indent}{key}: {name}")
    cfg.path.write_text("\n".join(lines), encoding="utf-8")
    return True


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="palmgeo",
        description="Build palm_csd inputs from open geodata.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--version", action="version", version=f"palmgeo {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="command", required=True)

    for stage in [*STAGES, "all"]:
        p = sub.add_parser(stage, help=f"run the '{stage}' stage")
        p.add_argument("config", help="path to the run configuration YAML")
        if stage == "sources":
            p.add_argument(
                "--approve", action="store_true",
                help="accept these sources and their licences, unlocking the "
                     "download stages")
            p.add_argument(
                "--set", action="append", default=[], metavar="ROLE=SOURCE",
                help="choose the source for a role and write it into the run "
                     "configuration, e.g. --set dsm=<source>. Repeatable.")
            p.add_argument(
                "--no-alternatives", action="store_true",
                help="do not probe the coverage of sources that are not chosen "
                     "(faster; each probe may fetch an index)")

    p = sub.add_parser(
        "fetch",
        help="download every source the run needs, so the rest can run offline")
    p.add_argument("config")
    p.add_argument("--check", action="store_true",
                   help="only report what is missing; download nothing")

    p = sub.add_parser("calibrate", help="print an NDVI histogram for threshold tuning")
    p.add_argument("config")
    p.add_argument("--tiles", type=int, default=12)

    p = sub.add_parser("report", help="list source vocabulary without a mapping rule")
    p.add_argument("config")
    p.add_argument("--what", default="all", choices=["all", "landuse", "street_type"])

    p = sub.add_parser("lonlat", help="convert WGS84 to the projected origin")
    p.add_argument("lon", type=float)
    p.add_argument("lat", type=float)
    p.add_argument("--epsg", type=int, required=True,
                   help="the run's grid.epsg -- the CRS the region's products are in")

    p = sub.add_parser("clean", help="delete intermediate data")
    p.add_argument("config")
    p.add_argument(
        "what", choices=["cache", "tiers", "surface", "vector", "work", "out"]
    )

    args = parser.parse_args(argv)
    _setup_logging(args.verbose)

    if args.command == "lonlat":
        cmd_lonlat(args.lon, args.lat, args.epsg)
        return 0

    try:
        cfg = load(args.config)
    except ConfigError as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        return 2

    try:
        if args.command == "fetch":
            _save_state(cfg, "fetch", cmd_fetch(cfg, check_only=args.check))
        elif args.command == "calibrate":
            cmd_calibrate(cfg, args.tiles)
        elif args.command == "report":
            cmd_report(cfg, args.what)
        elif args.command == "clean":
            cmd_clean(cfg, args.what)
        elif args.command == "sources":
            for assignment in args.set:
                if "=" not in assignment:
                    raise ConfigError(f"--set expects ROLE=SOURCE, got '{assignment}'")
                role, _, name = assignment.partition("=")
                if set_source_in_config(cfg, role.strip(), name.strip()):
                    logger.info("%s: sources.%s -> %s", cfg.path.name, role.strip(), name.strip())
            if args.set:
                cfg = load(args.config)          # re-read what we just wrote
            _save_state(cfg, "sources", stage_sources(
                cfg, approve=args.approve,
                probe_alternatives=not args.no_alternatives))
        elif args.command == "all":
            for stage in STAGES:
                if stage in GATED_STAGES:
                    require_source_approval(cfg)
                logger.info("=== stage: %s", stage)
                _save_state(cfg, stage, STAGE_FUNCS[stage](cfg) or {})
        else:
            if args.command in GATED_STAGES:
                require_source_approval(cfg)
            _save_state(cfg, args.command, STAGE_FUNCS[args.command](cfg) or {})
    except (RuntimeError, ConfigError) as exc:
        logger.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        logger.warning("Interrupted; re-run the same stage to resume.")
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
