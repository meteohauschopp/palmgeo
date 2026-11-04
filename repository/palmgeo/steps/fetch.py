"""Enumerate, and optionally download, every source file a run needs.

palmgeo streams by default: a source read far below its native resolution is
read through ``/vsicurl/`` at the matching overview level, which for a 20 cm
orthophoto read on a 1 m grid moves about 1/25 of the file. That is the single
reason a city-scale run is affordable — measured over Lower Saxony, streaming a
run costs 11 GB where downloading the same sources whole costs 345 GB.

Streaming has one cost, and this module exists for it: **the network is needed
while the rasters are being processed, not only beforehand.** Where that is not
acceptable — a compute node with no outbound route, a field deployment, an
archive that must be reproducible in ten years — the sources have to be on disk
first, and then the 30x transfer is the price of independence.

So there are two boundaries, and a run picks one:

``after rasters`` (the default)
    Acquisition is ``sources``, ``plan``, ``vectors``, ``rasters`` and
    ``buildings``. Everything after them reads only the run's own working
    directory and needs no network at all. For a nested city run that is 3.5 GB
    of derived data to carry to an offline machine, against 345 GB of sources —
    two orders of magnitude cheaper, and it is what most people mean when they
    say the heavy part should run somewhere else.

``after fetch`` (``palmgeo fetch``)
    Every source file the run will touch is downloaded whole, first. After that
    the entire pipeline runs with ``runtime.offline`` set and never opens a
    socket. This is the expensive answer, and it is the right one when the
    machine that processes genuinely cannot reach the internet.

Both rest on the same enumeration below, which is also what an offline run uses
to fail usefully: it reports everything that is missing at once, rather than
discovering it one file at a time across a few thousand retry cycles.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

from palmgeo.cache import Cache, OfflineError
from palmgeo.config import Config
from palmgeo.sources import load_index, regions_for_bounds

logger = logging.getLogger("palmgeo.fetch")


@dataclass(frozen=True)
class Requirement:
    """One source file the run will open, and what wants it."""

    url: str
    role: str
    #: Source name from the region file, used as the cache sub-directory so a
    #: download lands where the stage that needs it will look for it.
    source: str
    #: Which stage would fetch this, so a failure can name the fix.
    stage: str
    #: What it covers, for the report: a tile key, or the area of interest.
    where: str = ""

    @property
    def is_local(self) -> bool:
        """Whether this is a staged file rather than something to download."""
        return "://" not in self.url or self.url.startswith("file://")


def _aoi(cfg: Config):
    xs0, ys0, xs1, ys1 = zip(*(t.bounds for t in cfg.tiers))
    return (min(xs0), min(ys0), max(xs1), max(ys1))


def required(cfg: Config, cache: Cache) -> List[Requirement]:
    """Every source file this run will open, deduplicated, in stage order.

    Deliberately built from the same planning code the stages themselves use --
    ``plan_tiles`` for the rasters, the index searches for the vectors -- rather
    than from a parallel description that could drift. If this list is wrong,
    the run is wrong in the same way.
    """
    from palmgeo.steps.tiles import plan_tiles

    seen: Dict[str, Requirement] = {}

    def add(req: Requirement) -> None:
        seen.setdefault(req.url, req)

    # ---- rasters: one dtm, one dsm and one ortho item per work tile -------
    for tier in cfg.tiers:
        plans, _ = plan_tiles(cfg, cache, tier)
        names = {role: cfg.source_name_for_role(role) for role in ("dtm", "dsm", "ortho")}
        for plan in plans:
            for role, item in (("dtm", plan.dtm), ("dsm", plan.dsm)):
                if item is not None:
                    add(Requirement(item.url, role, names[role] or role,
                                    "rasters", plan.tile.key))
            for epoch, item in plan.ortho.items():
                if item is not None:
                    add(Requirement(item.url, f"ortho:{epoch}", names["ortho"] or "ortho",
                                    "rasters", plan.tile.key))

    aoi = _aoi(cfg)

    # ---- vectors: the cadastre units covering the AOI ---------------------
    name, src = cfg.source_for_role("cadastre")
    index = load_index(name, src, cache, transform_to_epsg=cfg.epsg)
    for unit in regions_for_bounds(index, aoi):
        add(Requirement(unit.url, "cadastre", name, "vectors", unit.key))

    # ---- vectors: road centrelines, if street handling is on -------------
    if cfg.get("sources.street_type.enabled", False):
        try:
            topo_name, topo_src = cfg.source_for_role("topo")
        except Exception:  # noqa: BLE001 - a region without one is not an error
            topo_name, topo_src = None, {}
        url = topo_src.get("url") if topo_src else None
        if url and topo_name:
            add(Requirement(str(url), "topo", topo_name, "vectors", "whole AOI"))

    # ---- buildings: the 3D model, where it is used at all -----------------
    if "lod2" in {str(cfg.get("buildings.height_source", "ndsm")),
                  str(cfg.get("buildings.height_fallback", "none"))}:
        try:
            b3_name, b3_src = cfg.source_for_role("buildings3d")
        except Exception:  # noqa: BLE001
            b3_name, b3_src = None, {}
        if b3_name:
            b3_index = load_index(b3_name, b3_src, cache, transform_to_epsg=cfg.epsg)
            for item in b3_index.search(aoi):
                add(Requirement(item.url, "buildings3d", b3_name, "buildings", item.key))

    order = {"rasters": 0, "vectors": 1, "buildings": 2}
    return sorted(seen.values(), key=lambda r: (order.get(r.stage, 9), r.role, r.where))


def missing(cache: Cache, reqs: List[Requirement]) -> List[Requirement]:
    """Which of ``reqs`` are not on disk right now.

    A staged file is checked where it stands; a remote one is checked against
    the cache under the sub-directory the fetching stage would use, because
    that is where it would have been put.
    """
    out: List[Requirement] = []
    for req in reqs:
        if req.is_local:
            raw = req.url[len("file://"):] if req.url.startswith("file://") else req.url
            path = Path(raw) if raw.startswith("/") else (cache.root or Path(".")) / raw
            if not path.is_file():
                out.append(req)
            continue
        name = req.url.rsplit("/", 1)[-1].split("?")[0]
        if not (cache.dir / req.source / name).is_file():
            out.append(req)
    return out


def describe_missing(reqs: List[Requirement], limit: int = 12) -> str:
    """A report an offline run can abort with, naming the whole gap at once."""
    by_stage: Dict[str, List[Requirement]] = {}
    for r in reqs:
        by_stage.setdefault(r.stage, []).append(r)
    lines = [f"{len(reqs)} source file(s) are needed and not on disk:"]
    for stage in ("vectors", "rasters", "buildings"):
        group = by_stage.get(stage)
        if not group:
            continue
        roles = sorted({r.role for r in group})
        lines.append(f"  {stage:10} {len(group):5} file(s)   roles: {', '.join(roles)}")
        for r in group[:limit // 3]:
            lines.append(f"             {r.role} {r.where}: {r.url}")
        if len(group) > limit // 3:
            lines.append(f"             ... and {len(group) - limit // 3} more")
    lines.append("")
    lines.append("Either run the acquisition stages with network access, or fetch")
    lines.append("everything up front so the rest can run offline:")
    lines.append("    python -m palmgeo.cli fetch <config>")
    return "\n".join(lines)


def ensure_available(cfg: Config, cache: Cache, stage: str) -> None:
    """Abort an offline run before it starts, listing everything it would need.

    Called by the stages that would otherwise fetch. The alternative -- letting
    each file fail as it is reached -- turns one missing tile into a few
    thousand retry cycles and buries the useful part of the message.
    """
    if not cache.offline:
        return
    reqs = [r for r in required(cfg, cache) if r.stage == stage]
    gap = missing(cache, reqs)
    if gap:
        raise OfflineError(
            f"Stage '{stage}' cannot run offline.\n" + describe_missing(gap)
        )
    logger.info("offline: all %d source file(s) for '%s' are on disk", len(reqs), stage)


def run(cfg: Config, cache: Cache, check_only: bool = False) -> Dict[str, Any]:
    """Download every source file the run needs, or report what is absent."""
    reqs = required(cfg, cache)
    gap = missing(cache, reqs)
    local = sum(1 for r in reqs if r.is_local)

    print(f"\nSource files this run will open: {len(reqs)} "
          f"({local} staged locally, {len(reqs) - local} remote)")
    print(f"Not on disk: {len(gap)}\n")
    if check_only or not gap:
        if gap:
            print(describe_missing(gap))
        else:
            print("Everything is on disk; the whole pipeline can run offline.\n")
        return {"required": len(reqs), "missing": len(gap), "downloaded": 0}

    if not cache.keep:
        # Without this the stages release each file after reading it, which is
        # correct for a streaming run and exactly wrong here: the point of
        # fetching is that the file is still there afterwards.
        raise RuntimeError(
            "`fetch` needs runtime.keep_downloads: true, otherwise every file it "
            "downloads is deleted again by the stage that reads it."
        )

    done = 0
    for req in gap:
        if req.is_local:
            continue
        cache.fetch(req.url, subdir=req.source)
        done += 1
        if done % 25 == 0 or done == len(gap):
            logger.info("  %d/%d fetched", done, len(gap))
    logger.info("%s", cache.summary())
    return {"required": len(reqs), "missing": len(gap), "downloaded": done}
