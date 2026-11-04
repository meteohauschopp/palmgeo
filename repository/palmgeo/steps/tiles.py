"""Per-tile raster derivation: terrain, normalised surface model, vegetation indices.

This is where the bulk of the data volume is touched, and where the
download -> process -> delete cycle happens. One work tile is 1 km x 1 km; at
1 m work resolution that is a 1000 x 1000 array per product, so a worker's
peak memory is a few tens of megabytes regardless of how large the area of
interest is.

Two access modes are available per source:

``stream``
    Read directly from the server through ``/vsicurl/``. Only the overview
    level matching the requested work resolution is transferred, which for a
    20 cm orthophoto read at 1 m is roughly 1/25 of the file. Nothing is
    written to disk.
``download``
    Fetch the whole tile into the bounded cache, read it, then release it.
    Used when the work resolution is close to the source resolution, where an
    overview would not help anyway.

Provenance of the nDSM
----------------------
The normalised digital surface model written here -- ``ndsm = dsm - dtm``,
clamped at zero, nodata propagated from either input -- is the canopy height
model of `canopy_height_berechnen.py` by **Joana-Pauline Rentz**, from her
Bachelor's thesis (full reference in :mod:`palmgeo.steps.trees`). The formula is trivial; what
is not trivial, and what the original established, is that it only means
anything when the two elevation models are co-registered and share a tiling, so
that a pixel of the surface model and the pixel of the terrain model beneath it
describe the same ground. Where they do not, the difference is not a canopy
height but a registration error, and it looks exactly like a tree. The source
requirements in PORTING.md exist for this reason.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from rasterio.enums import Resampling

from palmgeo.cache import Cache
from palmgeo.config import Config, ConfigError, marker_up_to_date
from palmgeo.grid import Tier, Tile
from palmgeo.rasters import gdal_env, read_on_grid, write_raster
from palmgeo.sources import (
    EpochSpec,
    Index,
    Item,
    load_index,
    newest_per_tile,
    pick_paired,
    select_epochs,
)

logger = logging.getLogger("palmgeo.tiles")

#: Streaming pays off once an overview level is at least this much coarser
#: than the native resolution.
STREAM_RATIO = 4.0

#: Epoch role whose RGB+NIR bands are kept on disk. The leaf-on image is the
#: one that carries foliage, so it feeds both the greenness criterion of the
#: vegetation mask and the single-tree segmentation.
KEEP_BANDS_ROLE = "leaf_on"


def _to_byte_cube(
    stack: "np.ndarray", ortho_src: Dict[str, Any], ortho_name: str, tile_key: str
) -> "np.ndarray":
    """Convert the kept visible/NIR bands to the 8-bit cube the later stages read.

    Why 8 bit at all: everything downstream that reads these bands was
    calibrated on 8-bit imagery. The broadleaf/conifer vote in ``steps/trees``
    compares a green-minus-red difference against the absolute thresholds 9.0
    and 15.0, which are numbers in 8-bit digital counts and mean nothing in
    another range. Normalising here keeps one calibration valid for every
    region instead of making each region re-tune a published method.

    Why the source has to declare its range: a plain cast silently saturates.
    Sentinel-2 L2A carries surface reflectance scaled by 10 000, so a scene
    whose bands average 1700-4200 becomes a uniform white image under
    ``clip(0, 255)`` -- which is what happened, undetected, to every driver
    built from the global fallback sources. Saturated bands do not look broken
    downstream: they look like a place where the greenness criterion never
    matches, because ``g/(r+g+b)`` is exactly 1/3 in all three channels and
    the test needs more than 0.35.

    So the range is declared per source as ``value_max``, and a source that
    declares none may not exceed 255. Refusing is the point: an undeclared
    16-bit source is a configuration gap, and the alternative to an error
    message is a silently degraded classification nobody notices for months.
    """
    declared = ortho_src.get("value_max")
    peak = float(np.nanmax(stack)) if np.isfinite(stack).any() else 0.0

    if declared is None:
        if peak > 255.0:
            raise RuntimeError(
                f"Ortho source '{ortho_name}' delivers values up to {peak:.0f} on "
                f"tile {tile_key}, but declares no `value_max`, so palmgeo cannot "
                f"scale it to the 8-bit range the classification is calibrated "
                f"for. Add `value_max:` to the source in the region file -- 10000 "
                f"for Sentinel-2 L2A surface reflectance, 255 for ordinary 8-bit "
                f"aerial imagery, 65535 for raw 16-bit."
            )
        scale = 1.0
    else:
        scale = 255.0 / float(declared)
        if peak > float(declared) * 1.001:
            # Not fatal: reflectance legitimately exceeds 1.0 over bright
            # cloud and snow. Worth saying, because the other explanation is
            # a wrong `value_max`, and the two look identical here.
            logger.warning(
                "Tile %s: '%s' reaches %.0f, above its declared value_max %s. "
                "Values above the declared maximum are clipped. Bright cloud or "
                "snow does this legitimately; so does a wrong value_max.",
                tile_key, ortho_name, peak, declared,
            )
    return np.clip(np.nan_to_num(stack * scale, nan=0.0), 0, 255).astype("uint8")


@dataclass
class TilePlan:
    """Everything one work tile needs, resolved to concrete source items."""

    tile: Tile
    dtm: Optional[Item]
    dsm: Optional[Item]
    ortho: Dict[str, Optional[Item]]


def _access_mode(cfg: Config, source_resolution: float, work_resolution: float) -> str:
    """Stream from the server, or pull the whole file into the cache first.

    Streaming wins when the work grid is much coarser than the source, because
    then only the matching overview level is transferred -- a 20 cm orthophoto
    read on a 1 m grid moves roughly 1/25 of the file. Below that ratio the
    overview is no cheaper than the file and the extra range requests are pure
    overhead.

    An offline run has no third option: there is no server to stream from, so
    the only readable copy is a local one, and `download` is what looks for it.
    """
    if bool(cfg.runtime.get("offline", False)):
        return "download"
    mode = str(cfg.runtime.get("access", "auto"))
    if mode in {"stream", "download"}:
        return mode
    return "stream" if work_resolution / max(source_resolution, 1e-9) >= STREAM_RATIO else "download"


def _open_path(
    cfg: Config, cache: Cache, item: Item, source_name: str, mode: str
) -> Tuple[str, Optional[Path]]:
    """Where to read this item from, and what to release afterwards.

    Streaming is skipped when the file is already in the cache. Fetching the
    overview level of a file that is sitting on disk in full is slower than
    reading it and moves bytes for nothing -- and after `palmgeo fetch` that is
    the normal situation, which is what makes an offline run work at all
    despite `auto` preferring to stream.
    """
    if mode == "stream":
        cached = cache.cached_path(item.url, subdir=source_name)
        if cached is None:
            return item.url, None
        # Returned as the release target like any other cache entry: it was
        # downloaded into the cache, so the ordinary keep/release policy is the
        # right one. `fetch` sets keep_downloads, which makes release a no-op.
        return str(cached), cached
    local = cache.fetch(item.url, subdir=source_name)
    return str(local), local


def plan_tiles(cfg: Config, cache: Cache, tier: Tier) -> Tuple[List[TilePlan], Dict[str, Any]]:
    """Resolve, for every tile of a tier, which source tile and epoch to use."""
    dtm_name, dtm_src = cfg.source_for_role("dtm")
    dsm_name, dsm_src = cfg.source_for_role("dsm")
    ortho_name, ortho_src = cfg.source_for_role("ortho")

    dtm_index = load_index(dtm_name, dtm_src, cache, transform_to_epsg=cfg.epsg)
    dsm_index = load_index(dsm_name, dsm_src, cache, transform_to_epsg=cfg.epsg)
    ortho_index = load_index(ortho_name, ortho_src, cache, transform_to_epsg=cfg.epsg)

    tiles = tier.tiles
    bounds = {t.key: t.bounds for t in tiles}
    dtm_choice = newest_per_tile(dtm_index, bounds)

    pair_on_date = dsm_src.get("pair_with") == dtm_name
    dsm_choice: Dict[str, Optional[Item]] = {}
    for t in tiles:
        ref = dtm_choice.get(t.key)
        if ref is not None and pair_on_date:
            dsm_choice[t.key] = pick_paired(dsm_index, ref, pair_on_date=True)
        else:
            hits = dsm_index.search(t.bounds)
            dsm_choice[t.key] = (
                max(hits, key=lambda it: it.date or date.min) if hits else None
            )

    reports: List[Dict[str, Any]] = []
    ortho_choice: Dict[str, Dict[str, Optional[Item]]] = {t.key: {} for t in tiles}
    for role, raw in (cfg.get("sources.ortho_epochs", {}) or {}).items():
        spec = EpochSpec.from_dict(role, raw or {})
        chosen, report = select_epochs(
            ortho_index,
            bounds,
            spec,
            prefer_uniform=bool(cfg.get("sources.prefer_uniform_epoch", True)),
            max_spread_days=int(cfg.get("sources.max_epoch_spread_days", 30)),
        )
        missing = [k for k, v in chosen.items() if v.item is None]
        if missing and spec.required:
            raise RuntimeError(
                f"Ortho epoch '{role}' is required but {len(missing)} tile(s) have no "
                f"coverage at all, e.g. {missing[:5]}. Relax sources.ortho_epochs."
                f"{role}.min_year or shrink the domains."
            )
        report["missing_tiles"] = len(missing)
        reports.append(report)
        for key, choice in chosen.items():
            ortho_choice[key][role] = choice.item

    plans = [
        TilePlan(t, dtm_choice.get(t.key), dsm_choice.get(t.key), ortho_choice[t.key])
        for t in tiles
    ]
    no_terrain = [p.tile.key for p in plans if p.dtm is None]
    if no_terrain:
        logger.warning(
            "%d tile(s) have no terrain coverage (beyond the source's coverage?): %s",
            len(no_terrain), no_terrain[:5],
        )
    # The orthophoto epochs were already reported; the elevation vintage was
    # not recorded anywhere, and it is the one that decides whether a building
    # exists at all -- a footprint the cadastre added after the last laser
    # scan has no height and is dropped.
    def _dates(items) -> List[str]:
        return sorted({it.date.isoformat() for it in items if it and it.date})

    summary = {
        "tier": tier.name,
        "tiles": len(plans),
        "area_km2": round(tier.area_km2_covered(), 1),
        "ortho": reports,
        "tiles_without_terrain": len(no_terrain),
        "dtm_dates": _dates(p.dtm for p in plans),
        "dsm_dates": _dates(p.dsm for p in plans),
    }
    return plans, summary


def process_tile(cfg: Config, cache: Cache, tier: Tier, plan: TilePlan) -> Optional[Path]:
    """Derive zt, nDSM and the vegetation indices for one work tile."""
    outdir = cfg.tier_dir(tier) / "tiles" / plan.tile.key
    done = outdir / ".done"
    fingerprint = cfg.stage_fingerprint("tiles")
    if marker_up_to_date(done, fingerprint):
        return outdir
    if plan.dtm is None:
        return None

    res = tier.resolution
    transform = plan.tile.transform(res)
    shape = plan.tile.shape(res)
    crs = cfg.crs
    opts = cfg.get("output.cog", {}) or {}

    dtm_name, dtm_src = cfg.source_for_role("dtm")
    dsm_name, dsm_src = cfg.source_for_role("dsm")
    ortho_name, ortho_src = cfg.source_for_role("ortho")

    released: List[Path] = []
    try:
        # ---- terrain -------------------------------------------------------
        mode = _access_mode(cfg, float(dtm_src.get("resolution", 1.0)), res)
        path, local = _open_path(cfg, cache, plan.dtm, dtm_name, mode)
        if local:
            released.append(local)
        zt = read_on_grid(
            path, transform, shape, band=1,
            resampling=Resampling.average, dtype="float32", fill=np.nan,
        )
        nodata = dtm_src.get("nodata")
        if nodata is not None and np.isfinite(float(nodata)):
            zt = np.where(np.isclose(zt, float(nodata)), np.nan, zt)
        # A DTM value of exactly 0 over land is a fill artefact, but 0 m above
        # the vertical datum is legitimate on the coast, so only the documented
        # nodata is removed.
        write_raster(outdir / "zt.tif", zt, transform, crs, np.nan, opts, overviews=False)
        # Beside the raster, as the ortho epochs already were. Cheap, and it is
        # what lets somebody answer "why is this building missing" months later
        # without re-deriving the plan.
        (outdir / "zt.date").write_text(
            plan.dtm.date.isoformat() if plan.dtm.date else "unknown", encoding="utf-8")

        # ---- normalised surface model --------------------------------------
        ndsm = np.full(shape, np.nan, dtype="float32")
        if plan.dsm is not None:
            mode = _access_mode(cfg, float(dsm_src.get("resolution", 1.0)), res)
            path, local = _open_path(cfg, cache, plan.dsm, dsm_name, mode)
            if local:
                released.append(local)
            dsm = read_on_grid(
                path, transform, shape, band=1,
                resampling=Resampling.average, dtype="float32", fill=np.nan,
            )
            dsm_nodata = dsm_src.get("nodata")
            if dsm_nodata is not None and isinstance(dsm_nodata, (int, float)) and np.isfinite(float(dsm_nodata)):
                dsm = np.where(np.isclose(dsm, float(dsm_nodata)), np.nan, dsm)
            ndsm = dsm - zt
            lo = float(cfg.get("classification.ndsm.min", 0.0))
            hi = float(cfg.get("classification.ndsm.max", 120.0))

            # A surface model below the terrain beneath it is physically
            # impossible, so every negative value is an error of some kind.
            # Clamping them to zero is right for the small ones -- they are
            # interpolation noise between two independently gridded products --
            # and dangerous for the large ones, because a systematic vertical
            # or horizontal offset between the two models produces exactly the
            # same signature and the clamp turns it into a plausible-looking
            # field of zeros. That is the failure this module's docstring
            # warns about, and it is invisible once the clamp has run.
            #
            # So measure it first. A few per cent of shallow negatives is
            # normal; a large share, or deep ones, means the two models do not
            # describe the same ground and no height derived from them means
            # anything.
            finite = np.isfinite(ndsm)
            n_finite = int(finite.sum())
            if n_finite:
                negative = finite & (ndsm < 0.0)
                share = float(negative.sum()) / n_finite
                deepest = (
                    float(np.nanmin(np.where(negative, ndsm, np.nan)))
                    if negative.any() else 0.0
                )
                warn_share = float(
                    cfg.get("classification.ndsm.negative_warn_fraction", 0.05)
                )
                # Depth separates the two causes, and they need different
                # answers, so the message must not conflate them.
                #
                # SHALLOW and widespread is a property of the models, not a
                # fault: where a bare-earth model is *derived* from a surface
                # model by removing buildings and canopy (FABDEM from
                # Copernicus DEM, for instance), the removal over-corrects in
                # places and the "terrain" ends up above the "surface" by a
                # fraction of a metre. Nothing is misplaced; the two simply
                # disagree within their own accuracy.
                #
                # DEEP means the two grids do not describe the same ground --
                # a vertical datum difference or a horizontal shift -- and
                # then every height derived from them is wrong by that amount,
                # not noisy by it.
                deep = deepest < -float(cfg.get("classification.ndsm.negative_deep_m", 3.0))
                if share > warn_share and deep:
                    logger.warning(
                        "Tile %s: %.1f %% of the nDSM is negative, down to %.1f m. "
                        "That is too deep to be model disagreement: '%s' and '%s' "
                        "are probably offset from each other vertically or "
                        "horizontally, and every height derived from them carries "
                        "that offset.",
                        plan.tile.key, 100.0 * share, deepest, dsm_name, dtm_name,
                    )
                elif share > warn_share:
                    logger.info(
                        "Tile %s: %.1f %% of the nDSM is negative but only down to "
                        "%.1f m, and is clamped to zero. Shallow negatives are "
                        "normal where the terrain model is derived from the "
                        "surface model ('%s' from '%s'); they are the removal "
                        "over-correcting, not a registration error.",
                        plan.tile.key, 100.0 * share, deepest, dtm_name, dsm_name,
                    )
            ndsm = np.where(finite, np.clip(ndsm, lo, hi), np.nan)
            smooth = int(cfg.get("classification.ndsm.smooth_radius", 0))
            if smooth > 0:
                from scipy.ndimage import median_filter

                valid = np.isfinite(ndsm)
                filled = np.where(valid, ndsm, 0.0)
                ndsm = np.where(
                    valid, median_filter(filled, size=2 * smooth + 1), np.nan
                )
        write_raster(outdir / "ndsm.tif", ndsm, transform, crs, np.nan, opts, overviews=False)
        (outdir / "ndsm.date").write_text(
            plan.dsm.date.isoformat() if plan.dsm is not None and plan.dsm.date
            else "unknown", encoding="utf-8")

        # ---- vegetation indices --------------------------------------------
        # Required: band order differs between products (RGBN, CIR, RGB plus a
        # separate NIR), and guessing it yields a plausible, wrong NDVI.
        bands = ortho_src.get("bands")
        if not bands or "red" not in bands or "nir" not in bands:
            raise ConfigError(
                f"sources.{ortho_name}: `bands` must give the 1-based band numbers "
                f"of at least `red` and `nir`, e.g. {{red: 1, green: 2, blue: 3, nir: 4}}"
            )
        for role, item in plan.ortho.items():
            if item is None:
                continue
            mode = _access_mode(cfg, float(ortho_src.get("resolution", 0.2)), res)
            path, local = _open_path(cfg, cache, item, ortho_name, mode)
            if local:
                released.append(local)
            # NDVI is formed from block-averaged digital numbers rather than by
            # averaging per-pixel NDVI: averaging the radiances first is the
            # radiometrically better-behaved of the two and much less noisy.
            red = read_on_grid(
                path, transform, shape, band=int(bands["red"]),
                resampling=Resampling.average, dtype="float32", fill=np.nan,
            )
            nir = read_on_grid(
                path, transform, shape, band=int(bands["nir"]),
                resampling=Resampling.average, dtype="float32", fill=np.nan,
            )
            total = nir + red
            ndvi = np.where(total > 1e-6, (nir - red) / np.where(total > 1e-6, total, 1.0), np.nan)
            write_raster(
                outdir / f"ndvi_{role}.tif", ndvi.astype("float32"),
                transform, crs, np.nan, opts, overviews=False,
            )
            (outdir / f"ndvi_{role}.date").write_text(
                item.date.isoformat() if item.date else "unknown", encoding="utf-8"
            )

            # NDVI alone misses shaded and duller vegetation, and the tree
            # segmentation needs brightness and hue per crown, so the visible
            # bands are kept as well for the epoch that carries the foliage.
            if role == KEEP_BANDS_ROLE and cfg.get("output.keep_ortho_bands", True):
                stack = [red]
                for band_name in ("green", "blue"):
                    index = bands.get(band_name)
                    stack.append(
                        read_on_grid(
                            path, transform, shape, band=int(index),
                            resampling=Resampling.average, dtype="float32", fill=np.nan,
                        )
                        if index
                        else np.full(shape, np.nan, dtype="float32")
                    )
                stack.append(nir)
                cube = _to_byte_cube(np.stack(stack), ortho_src, ortho_name, plan.tile.key)
                write_raster(
                    outdir / f"ortho_{role}.tif", cube,
                    transform, crs, None, opts, overviews=False,
                )
    except Exception:
        # never leave a half-written tile behind; the next run must redo it
        for f in outdir.glob("*"):
            f.unlink(missing_ok=True)
        raise
    finally:
        for p in released:
            cache.release(p)

    done.write_text(fingerprint)
    return outdir


def run(cfg: Config, cache: Cache, tier: Tier) -> Dict[str, Any]:
    """Process every tile of one resolution tier."""
    plans, summary = plan_tiles(cfg, cache, tier)
    workers = int(cfg.runtime.get("workers", 4))
    cache_mb = int(cfg.runtime.get("gdal_cachemax_mb", 512))
    logger.info(
        "Tier %s: %d tile(s) over %.0f km2 at %.0f m, %d worker(s)",
        tier.name, len(plans), summary["area_km2"], tier.resolution, workers,
    )

    failures: List[Tuple[str, str]] = []
    completed = 0
    with gdal_env(cache_mb):
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_tile, cfg, cache, tier, p): p for p in plans
            }
            for future in as_completed(futures):
                plan = futures[future]
                try:
                    future.result()
                except Exception as exc:  # noqa: BLE001 - collect and report
                    failures.append((plan.tile.key, str(exc)))
                    logger.error("Tile %s failed: %s", plan.tile.key, exc)
                completed += 1
                if completed % 25 == 0 or completed == len(plans):
                    logger.info("  %d/%d tiles", completed, len(plans))

    summary["failed"] = len(failures)
    if failures:
        summary["failures"] = failures[:20]
        logger.error(
            "%d tile(s) failed. Re-running the stage retries only those.", len(failures)
        )
    logger.info("%s", cache.summary())
    return summary
