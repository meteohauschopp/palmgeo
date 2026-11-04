"""Aggregation of the work grid onto the PALM domain grids.

Each domain gets its own directory of GeoTIFFs written on **exactly** its grid
(same CRS, same pixel size, same origin). That is not cosmetic: it is the only
way palm_csd takes its windowed read path instead of loading whole rasters into
memory -- see the module docstring of :mod:`palmgeo.grid`.

Aggregating surface types needs care. Taking the mode of ``vegetation_type``
and of ``pavement_type`` independently would let a cell that is 40 % pavement
and 60 % vegetation end up carrying both, and palm_csd's precedence would then
silently award it to pavement. Instead the *group* raster decides the winner by
majority first, and only the winning group's subtype is kept.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import rasterio

from palmgeo import palm_types as pt
from palmgeo.config import Config
from palmgeo.grid import Domain, Tier, alignment_report
from palmgeo.rasters import RESAMPLING, build_vrt, gdal_env, resample_to_domain, write_raster
from palmgeo.steps.compose import I32_NODATA, PRODUCTS, U8_NODATA

logger = logging.getLogger("palmgeo.domains")

#: Written into the palm_csd `files:` block. Anything not listed here is
#: internal bookkeeping (surface_group, surface_class).
CSD_FILES = [
    "zt",
    "vegetation_type",
    "pavement_type",
    "water_type",
    "soil_type",
    "vegetation_height",
    "buildings_2d",
    "building_id",
    "building_type",
    "bridges_2d",
    "bridges_id",
    "street_type",
    "street_crossings",
]


def build_mosaics(cfg: Config, tier: Tier) -> Dict[str, Path]:
    """Create one VRT per product over all composed tiles of a tier."""
    surface_dir = cfg.tier_dir(tier) / "surface"
    vrt_dir = cfg.tier_dir(tier) / "vrt"
    vrt_dir.mkdir(parents=True, exist_ok=True)
    mosaics: Dict[str, Path] = {}
    for name, meta in PRODUCTS.items():
        parts = sorted(surface_dir.glob(f"*/{name}.tif"))
        if not parts:
            continue
        nodata = meta["nodata"]
        mosaics[name] = build_vrt(
            parts, vrt_dir / f"{name}.vrt",
            nodata=None if isinstance(nodata, float) and np.isnan(nodata) else nodata,
        )
    if not mosaics:
        raise RuntimeError(
            f"Tier {tier.name} has no composed tiles. Run the `compose` stage first."
        )
    logger.info("Tier %s: %d mosaic(s) over %d tile(s)",
                tier.name, len(mosaics), len(list(surface_dir.glob("*/zt.tif"))))
    return mosaics


def _aggregate(
    cfg: Config, domain: Domain, name: str, mosaic: Path, out_dir: Path
) -> Optional[Path]:
    meta = PRODUCTS[name]
    return resample_to_domain(
        mosaic,
        out_dir / f"{name}.tif",
        domain.transform,
        (domain.height, domain.width),
        cfg.crs,
        RESAMPLING[str(meta["resampling"])],
        str(meta["dtype"]),
        None if isinstance(meta["nodata"], float) and np.isnan(meta["nodata"]) else meta["nodata"],
        cfg.get("output.cog", {}) or {},
    )


def _load(path: Path) -> Tuple[np.ndarray, Any]:
    with rasterio.open(path) as src:
        return src.read(1), src.nodata


def build_domain(cfg: Config, domain: Domain, tier: Tier, mosaics: Dict[str, Path]) -> Path:
    """Write the full set of palm_csd input rasters for one domain."""
    out_dir = cfg.domain_dir(domain)
    out_dir.mkdir(parents=True, exist_ok=True)
    exact, message = alignment_report(domain, tier.resolution)
    (logger.debug if exact else logger.warning)("Aggregation %s", message)

    produced: Dict[str, Path] = {}
    for name in list(PRODUCTS):
        if name not in mosaics:
            continue
        path = _aggregate(cfg, domain, name, mosaics[name], out_dir)
        if path is not None:
            produced[name] = path

    # ---- resolve the surface group -----------------------------------------
    group_path = produced.get("surface_group")
    if group_path is None:
        raise RuntimeError("surface_group is missing; re-run the compose stage")
    group, _ = _load(group_path)

    # The group raster and the three type rasters are aggregated INDEPENDENTLY,
    # each by majority over its own block. That looks like it could disagree
    # with itself -- a block that is 60 % vegetation split across three
    # vegetation types, against 40 % nodata, would take the group from the
    # majority and the type from the nodata plurality, leaving a vegetated cell
    # with no vegetation type for palm_csd to read.
    #
    # It cannot, because GDAL's mode resampling excludes nodata from the vote,
    # so each type raster's majority is taken over the valid pixels only. The
    # loop below therefore only ever has to mask a type OUT where the group
    # disagrees; it never has to fill one in. verify_driver's "exactly one
    # surface type per cell" check is what confirms this holds, and it passes
    # at 100 % on every driver built so far.
    #
    # Written down because the alternative reading is plausible enough that
    # somebody will eventually add a fill step here that is not needed.
    opts = cfg.get("output.cog", {}) or {}
    for name, want in (
        ("vegetation_type", pt.GROUP_VEGETATION),
        ("pavement_type", pt.GROUP_PAVEMENT),
        ("water_type", pt.GROUP_WATER),
    ):
        path = produced.get(name)
        if path is None:
            continue
        data, _ = _load(path)
        masked = np.where(group == want, data, U8_NODATA).astype("uint8")
        write_raster(path, masked, domain.transform, cfg.crs, U8_NODATA, opts, overviews=False)

    # ---- buildings ---------------------------------------------------------
    is_building = group == pt.GROUP_BUILDING
    if "buildings_2d" in produced:
        data, _ = _load(produced["buildings_2d"])
        data = np.where(is_building & np.isfinite(data), data, np.nan).astype("float32")
        write_raster(produced["buildings_2d"], data, domain.transform, cfg.crs, np.nan, opts, overviews=False)
        building_present = np.isfinite(data)
    else:
        building_present = np.zeros((domain.height, domain.width), dtype=bool)
    for name, nodata, dtype in (
        ("building_id", I32_NODATA, "int32"),
        ("building_type", U8_NODATA, "uint8"),
    ):
        if name not in produced:
            continue
        data, _ = _load(produced[name])
        data = np.where(building_present, data, nodata).astype(dtype)
        write_raster(produced[name], data, domain.transform, cfg.crs, nodata, opts, overviews=False)

    # ---- bridges -----------------------------------------------------------
    # palm_csd raises "Buildings and bridges are overlapping" if a cell carries
    # both, which aggregation can create even though the 1 m fields are disjoint.
    if "bridges_2d" in produced:
        data, _ = _load(produced["bridges_2d"])
        data = np.where(~building_present & np.isfinite(data), data, np.nan).astype("float32")
        write_raster(produced["bridges_2d"], data, domain.transform, cfg.crs, np.nan, opts, overviews=False)
        bridge_present = np.isfinite(data)
        if "bridges_id" in produced:
            ids, _ = _load(produced["bridges_id"])
            ids = np.where(bridge_present, ids, I32_NODATA).astype("int32")
            write_raster(produced["bridges_id"], ids, domain.transform, cfg.crs, I32_NODATA, opts, overviews=False)

    # ---- soil / vegetation height consistency ------------------------------
    if "soil_type" in produced:
        data, _ = _load(produced["soil_type"])
        keep = (group == pt.GROUP_VEGETATION) | (group == pt.GROUP_PAVEMENT)
        write_raster(
            produced["soil_type"], np.where(keep, data, U8_NODATA).astype("uint8"),
            domain.transform, cfg.crs, U8_NODATA, opts, overviews=False,
        )
    if "vegetation_height" in produced:
        data, _ = _load(produced["vegetation_height"])
        data = np.where((group == pt.GROUP_VEGETATION) & np.isfinite(data), data, np.nan)
        write_raster(
            produced["vegetation_height"], data.astype("float32"),
            domain.transform, cfg.crs, np.nan, opts, overviews=False,
        )

    # ---- terrain completeness ---------------------------------------------
    # palm_csd declares `zt` all_or_none_missing: a single masked cell aborts
    # the run, so gaps are filled here where they are small enough to be safe.
    if "zt" in produced:
        data, _ = _load(produced["zt"])
        gaps = ~np.isfinite(data)
        if gaps.any():
            fraction = gaps.mean()
            if fraction > 0.02:
                raise RuntimeError(
                    f"Domain {domain.name}: {fraction:.1%} of the terrain is missing. "
                    f"The domain probably extends beyond the data coverage "
                    f"(the edge of the source's coverage, or a gap in the terrain model)."
                )
            filled = _fill_nearest(data)
            logger.warning(
                "Domain %s: filled %d terrain gap cell(s) (%.3f %%) by nearest neighbour",
                domain.name, int(gaps.sum()), fraction * 100,
            )
            write_raster(produced["zt"], filled, domain.transform, cfg.crs, None, opts, overviews=False)

    # ---- hydro-flattening ---------------------------------------------------
    # A water surface is flat. A coarse terrain model does not know that: over
    # the Aswan domain the 30 m FABDEM bleeds the riverbank into the Nile
    # polygon, so the modelled river surface spanned 0-22 m and 18.5 % of its
    # cells sat more than 2 m above the water. Eroding the mask by 80 m
    # removes essentially all of it (max 5.2 m, p95 0.5 m), which is what
    # identifies the signal as an edge artefact rather than islands.
    #
    # PALM would read that relief as real topography sitting in the middle of
    # the river and force the flow over it. Each connected water body is
    # therefore levelled to a single elevation, taken as the median of its own
    # eroded interior -- per body, not globally, so a canal on a terrace is
    # not dragged down to the river it never meets.
    flatten = cfg.get("water.flatten", False)
    if flatten and "zt" in produced and "water_type" in produced:
        from scipy import ndimage

        zt, _ = _load(produced["zt"])
        water_type, _ = _load(produced["water_type"])
        wet = water_type != U8_NODATA
        if wet.any():
            erode_m = float(cfg.get("water.flatten_erosion_m", 80.0))
            radius = max(int(round(erode_m / domain.pixel_size)), 1)
            labels, n_bodies = ndimage.label(wet)
            zt = zt.astype("float32", copy=True)
            n_flat = 0
            for body in range(1, n_bodies + 1):
                mask = labels == body
                interior = ndimage.binary_erosion(mask, iterations=radius)
                # A body narrower than twice the erosion radius erodes away
                # entirely; fall back to its own median, which is still a far
                # better estimate of its level than the bank-contaminated
                # values scattered across it.
                sample = zt[interior] if interior.any() else zt[mask]
                sample = sample[np.isfinite(sample)]
                if sample.size == 0:
                    continue
                if isinstance(flatten, (int, float)) and not isinstance(flatten, bool):
                    level = float(flatten)
                else:
                    level = float(np.median(sample))
                before = zt[mask]
                spread = float(np.nanmax(before) - np.nanmin(before)) if before.size else 0.0
                zt[mask] = level
                n_flat += int(mask.sum())
                if mask.sum() >= 100:
                    logger.info(
                        "Domain %s: water body %d (%d cells) levelled to %.2f m "
                        "(was spread over %.2f m)",
                        domain.name, body, int(mask.sum()), level, spread,
                    )
            logger.info("Domain %s: hydro-flattened %d water cell(s) in %d body/bodies",
                        domain.name, n_flat, n_bodies)
            write_raster(produced["zt"], zt, domain.transform, cfg.crs, None, opts,
                         overviews=False)

    # internal-only products do not belong in the palm_csd input directory
    for name in ("surface_group", "surface_class"):
        produced.pop(name, None)
        (out_dir / f"{name}.tif").unlink(missing_ok=True)

    _write_manifest(cfg, domain, tier, out_dir, produced, exact, message)
    return out_dir


def _fill_nearest(data: np.ndarray) -> np.ndarray:
    """Fill non-finite cells from their nearest finite neighbour."""
    from scipy.ndimage import distance_transform_edt

    invalid = ~np.isfinite(data)
    if not invalid.any():
        return data.astype("float32")
    _, (rows, cols) = distance_transform_edt(invalid, return_indices=True)
    return data[rows, cols].astype("float32")


def _write_manifest(
    cfg: Config,
    domain: Domain,
    tier: Tier,
    out_dir: Path,
    produced: Dict[str, Path],
    exact: bool,
    message: str,
) -> None:
    import json

    # What the tiles of this tier were built from, and when it was flown.
    # Read back off disk rather than passed down, so the manifest describes the
    # data that actually exists rather than what the plan intended -- the two
    # differ after a partial rebuild, and the manifest is what somebody
    # consults months later to explain a building that is not there.
    acquired: Dict[str, str] = {}
    tile_dir = cfg.tier_dir(tier) / "tiles"
    for product, stem in (("dtm", "zt"), ("dsm", "ndsm")):
        dates = sorted({
            f.read_text(encoding="utf-8").strip()
            for f in tile_dir.glob(f"*/{stem}.date")
        } - {"", "unknown"})
        if dates:
            acquired[product] = dates[0] if len(dates) == 1 else f"{dates[0]} … {dates[-1]}"
    for f in sorted(tile_dir.glob("*/ndvi_*.date")):
        role = f.stem.replace("ndvi_", "ortho:")
        value = f.read_text(encoding="utf-8").strip()
        if value and value != "unknown":
            acquired.setdefault(role, value)

    stats: Dict[str, Any] = {}
    for name, path in sorted(produced.items()):
        with rasterio.open(path) as src:
            data = src.read(1, masked=True)
            valid = int(data.count())
            stats[name] = {
                "valid_cells": valid,
                "coverage": round(valid / data.size, 4) if data.size else 0.0,
            }
            if valid and np.issubdtype(data.dtype, np.floating):
                stats[name]["min"] = round(float(data.min()), 3)
                stats[name]["max"] = round(float(data.max()), 3)
    (out_dir / "manifest.json").write_text(
        json.dumps(
            {
                "domain": domain.name,
                "grid": {
                    "pixel_size": domain.pixel_size,
                    "nx": domain.nx,
                    "ny": domain.ny,
                    "width_px": domain.width,
                    "height_px": domain.height,
                    "origin_x": domain.origin_x,
                    "origin_y": domain.origin_y,
                    "epsg": cfg.epsg,
                },
                "work_tier": {"resolution": tier.resolution, "aggregation": message},
                "exact_aggregation": exact,
                # Acquisition dates of the sources this domain rests on. The
                # driver itself records none of them, and they are what makes a
                # missing building explicable rather than mysterious.
                "acquired": acquired,
                "products": stats,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def run(cfg: Config) -> Dict[str, Any]:
    """Build the palm_csd input rasters for every configured domain."""
    result: Dict[str, Any] = {"domains": {}}
    with gdal_env(int(cfg.runtime.get("gdal_cachemax_mb", 512))):
        for tier in cfg.tiers:
            mosaics = build_mosaics(cfg, tier)
            for domain in tier.domains:
                logger.info("Domain %s", domain.describe())
                out_dir = build_domain(cfg, domain, tier, mosaics)
                result["domains"][domain.name] = str(out_dir)
    return result
