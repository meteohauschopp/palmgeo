"""Single-tree detection from the canopy height model and the orthophoto.

Implements the segmentation and classification of Rentz as a pipeline stage:
the canopy height model (here the nDSM) is smoothed, local maxima are taken as
tree tops, and crowns are grown around them by a watershed whose elevation
combines the inverted canopy height with the edges of the visible image. Each
crown is then classified as broadleaf or coniferous from a weighted vote over
brightness, green-red difference, mean NDVI, crown shape and crown area, and
allometric relations turn the measured height into trunk diameter and leaf area
index.

The output columns map directly onto palm_csd's ``trees`` vector input, so the
result feeds ``generate_single_trees`` without further translation.

Changes against the original script, all made for the surrounding pipeline
rather than for the method:

* tiles are read with a halo and trees are kept only where their crown centroid
  falls inside the tile proper, so a crown straddling a tile edge is segmented
  whole and counted once; a tree is placed at that centroid, where the original
  placed it at the detected tree top;
* detections on building footprints are dropped, because a roof that reads as
  green is a green roof, not a tree;
* brightness, green-red difference and NDVI per crown are averaged from the
  bands on the work grid, where the original re-read the orthophoto at its
  native resolution inside each crown;
* lengths and areas are given in metres and converted to pixels, which equals
  the original's pixel values on its 1 m grid.

Provenance
----------
The method is a port of `klassifizierung_segmentierung.py` by **Joana-Pauline
Rentz**, written for her Bachelor's thesis:

    Rentz, J.-P. (2026): "Generierung von dreidimensionalen
    Einzelbaumrepräsentationen basierend auf Fernerkundungsdaten". Bachelor's
    thesis, Institute of Meteorology and Climatology, Leibniz University Hannover.

The original stood alone as a per-tile script. It is reproduced here rather than
imported, because the pipeline needs it tiled, parallel and resumable; the
original directory has been removed and this docstring, with the constants
below, is what remains of it.

Every numeric constant of the original is carried over, with the one exception
marked below, and they are the substance of the method -- each one was tuned against manually delineated
crowns on a 1 km orthophoto tile, so changing one without re-doing that work
changes the result without a basis for the change:

    smoothing            gaussian(sigma = 1.2 m) on the canopy height model
    tree tops            peak_local_max(min_distance = 3 m, threshold_abs = 2 m)
    crown mask           canopy height > 1 m; the original drops crowns under
                         20 px, here under half the minimum vegetation patch
                         (15 m2 by default) -- the exception
    watershed elevation  inverted smoothed height + 0.4 x normalised Sobel edges
    broadleaf/conifer    a weighted vote, conifer first:
                           brightness <= 72        conifer +3.0
                           brightness >= 82      broadleaf +3.0
                           green-red <= 9.0        conifer +2.0
                           green-red >= 15.0     broadleaf +2.0
                           mean NDVI <= 0.48       conifer +1.5
                           mean NDVI >= 0.58     broadleaf +1.5
                           crown form >= 1.25      conifer +0.5
                           roundness  >= 0.20    broadleaf +0.5
                           crown area <= 20 m2     conifer +1.0
                           crown area >= 26 m2   broadleaf +1.0
    trunk diameter       broadleaf  dbh = 1.31811 sqrt(h-1.3)
                                          / (1 - 0.151515 sqrt(h-1.3))
                         conifer    dbh = 19.5 / (3.85 - ln(h-1.3))
    leaf area index      broadleaf  lai = 0.307279 dbh^1.803
                                          / (pi (0.078307 dbh + 0.905689)^2)
                         conifer    lai = 12.25 (1 - exp(-0.1 dbh))^1.6

The vote is deliberately soft: no single band decides a tree. A conifer in
strong shadow loses the brightness and NDVI terms but still wins on crown form
and area, which is why the original used ten weighted terms rather than a
threshold on one.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from palmgeo.config import Config, ConfigError
from palmgeo.grid import Tier, Tile
from palmgeo.rasters import build_vrt, gdal_env, read_on_grid

logger = logging.getLogger("palmgeo.trees")

#: Single-tree segmentation is meaningless once a pixel is wider than a crown.
MAX_USEFUL_RESOLUTION = 2.0



def _drop_small(array: np.ndarray, min_pixels: int) -> np.ndarray:
    """Remove connected components below ``min_pixels``.

    Written against ``scipy.ndimage`` rather than ``skimage.morphology`` because
    the latter's ``min_size`` parameter is deprecated as of 0.26 with changed
    semantics, and this stage has to keep working across versions.
    """
    from scipy import ndimage

    if min_pixels <= 1:
        return array
    if array.dtype == bool:
        labels, count = ndimage.label(array)
    else:
        labels, count = array, int(array.max())
    if count == 0:
        return array
    sizes = np.bincount(labels.ravel())
    too_small = sizes < min_pixels
    too_small[0] = True
    drop = too_small[labels]
    return np.where(drop, 0, array)


def _params(cfg: Config) -> Dict[str, float]:
    section = cfg.get("vegetation.trees.segmentation", {}) or {}
    return {
        "ndvi_min": float(section.get("ndvi_min", 0.20)),
        "green_norm_min": float(section.get("green_norm_min", 0.35)),
        "green_over_red": float(section.get("green_over_red", 0.05)),
        "green_over_blue": float(section.get("green_over_blue", 0.005)),
        "height_min": float(section.get("height_min", 2.0)),
        "smooth_sigma_m": float(section.get("smooth_sigma_m", 1.2)),
        "min_distance_m": float(section.get("min_distance_m", 3.0)),
        "min_patch_area_m2": float(section.get("min_patch_area_m2", 30.0)),
        "min_crown_area_m2": float(section.get("min_crown_area_m2", 3.0)),
        "edge_weight": float(section.get("edge_weight", 0.4)),
        "halo_m": float(section.get("halo_m", 24.0)),
    }


def segment_tile(
    chm: np.ndarray,
    bands: np.ndarray,
    transform: Any,
    resolution: float,
    par: Dict[str, float],
) -> Optional["Any"]:
    """Segment one (haloed) window into crown polygons with attributes."""
    import geopandas as gpd
    from rasterio.features import shapes as rio_shapes
    from shapely.geometry import shape as to_shape
    from skimage.feature import peak_local_max
    from skimage.filters import gaussian, sobel
    from skimage.segmentation import watershed

    red, green, blue, nir = bands[0], bands[1], bands[2], bands[3]
    ndvi = (nir - red) / (nir + red + 1e-6)
    total = red + green + blue + 1e-6
    g_n, r_n, b_n = green / total, red / total, blue / total

    veg = (
        (ndvi > par["ndvi_min"])
        & (g_n > par["green_norm_min"])
        & (g_n > r_n + par["green_over_red"])
        & (g_n > b_n + par["green_over_blue"])
    )

    px_area = resolution * resolution
    chm = np.nan_to_num(chm, nan=0.0, posinf=0.0, neginf=0.0)
    smooth = gaussian(chm, sigma=max(par["smooth_sigma_m"] / resolution, 0.3), preserve_range=True)
    veg &= smooth >= par["height_min"]
    min_patch_px = max(1, int(round(par["min_patch_area_m2"] / px_area)))
    veg = _drop_small(veg, min_patch_px).astype(bool)
    if not veg.any():
        return None
    smooth = np.where(veg, smooth, 0.0)

    min_distance_px = max(1, int(round(par["min_distance_m"] / resolution)))
    peaks = peak_local_max(
        smooth,
        min_distance=min_distance_px,
        threshold_abs=par["height_min"],
        exclude_border=False,
    )
    if len(peaks) == 0:
        return None
    markers = np.zeros(chm.shape, dtype=np.int32)
    for index, (row, col) in enumerate(peaks, start=1):
        markers[row, col] = index

    grey = (red + green + blue) / 3.0
    edges = sobel(grey)
    edges = np.where(veg, edges, 0.0)
    elevation = -smooth + par["edge_weight"] * (edges / (edges.max() + 1e-6))
    labels = watershed(elevation, markers, mask=smooth > 1.0)
    labels = _drop_small(labels, max(1, min_patch_px // 2))
    if not labels.any():
        return None

    polygons, ids = [], []
    for geom, value in rio_shapes(labels.astype(np.int32), mask=labels > 0, transform=transform):
        polygons.append(to_shape(geom))
        ids.append(int(value))
    if not polygons:
        return None

    crowns = gpd.GeoDataFrame({"crown_id": ids}, geometry=polygons)
    crowns["crown_area"] = crowns.geometry.area
    crowns = crowns[crowns["crown_area"] > par["min_crown_area_m2"]].copy()
    if crowns.empty:
        return None
    crowns["roundness"] = [
        (4 * np.pi * g.area / (g.length**2)) if g.length > 0 else 0.0 for g in crowns.geometry
    ]
    crowns["CROWN_D"] = 2.0 * np.sqrt(crowns["crown_area"] / np.pi)

    # per-crown statistics, taken from the arrays directly instead of
    # re-opening the raster once per polygon as the original script does
    stats = _crown_statistics(crowns, labels, transform, chm, red, green, blue, ndvi)
    for key, values in stats.items():
        crowns[key] = values

    _classify(crowns)
    return crowns


def _crown_statistics(
    crowns: Any,
    labels: np.ndarray,
    transform: Any,
    chm: np.ndarray,
    red: np.ndarray,
    green: np.ndarray,
    blue: np.ndarray,
    ndvi: np.ndarray,
) -> Dict[str, List[float]]:
    """Aggregate the per-pixel arrays over each crown label in one pass."""
    order = {int(cid): i for i, cid in enumerate(crowns["crown_id"])}
    n = len(order)
    acc = {
        k: np.zeros(n) for k in ("h_max", "h_sum", "count", "bright", "green_red", "ndvi")
    }
    acc["h_max"][:] = -np.inf

    flat_labels = labels.ravel()
    valid = flat_labels > 0
    idx = flat_labels[valid]
    # A lookup table from watershed label to row in `crowns`, so the whole
    # tile can be accumulated with one vectorised pass instead of a Python
    # loop per crown. -1 marks a label that survived the watershed but was
    # dropped from `crowns` (too small, or its centroid fell in the halo),
    # and those pixels are filtered out below rather than accumulated into
    # some other crown's statistics.
    lookup = np.full(int(flat_labels.max()) + 1, -1, dtype=int)
    for cid, position in order.items():
        if cid <= int(flat_labels.max()):
            lookup[cid] = position
    positions = lookup[idx]
    good = positions >= 0
    positions = positions[good]

    def add(target: str, source: np.ndarray) -> None:
        np.add.at(acc[target], positions, source.ravel()[valid][good])

    add("h_sum", chm)
    add("bright", (red + green + blue) / 3.0)
    add("green_red", green - red)
    add("ndvi", ndvi)
    np.add.at(acc["count"], positions, 1.0)
    np.maximum.at(acc["h_max"], positions, chm.ravel()[valid][good])

    count = np.maximum(acc["count"], 1.0)
    h_max = np.where(np.isfinite(acc["h_max"]), acc["h_max"], 0.0)
    h_mean = acc["h_sum"] / count
    return {
        "HEIGHT_TOP": list(h_max),
        "h_mean": list(h_mean),
        "kronenform": list(h_max / (h_mean + 1e-6)),
        "helligkeit": list(acc["bright"] / count),
        "green_red": list(acc["green_red"] / count),
        "ndvi_mean": list(acc["ndvi"] / count),
    }


def _classify(crowns: Any) -> None:
    """Broadleaf/conifer vote and allometry, after Rentz."""
    conifer = np.zeros(len(crowns))
    broadleaf = np.zeros(len(crowns))

    conifer += (crowns["helligkeit"] <= 72.0) * 3.0
    broadleaf += (crowns["helligkeit"] >= 82.0) * 3.0
    conifer += (crowns["green_red"] <= 9.0) * 2.0
    broadleaf += (crowns["green_red"] >= 15.0) * 2.0
    conifer += (crowns["ndvi_mean"] <= 0.48) * 1.5
    broadleaf += (crowns["ndvi_mean"] >= 0.58) * 1.5
    conifer += (crowns["kronenform"] >= 1.25) * 0.5
    broadleaf += (crowns["roundness"] >= 0.20) * 0.5
    conifer += (crowns["crown_area"] <= 20.0) * 1.0
    broadleaf += (crowns["crown_area"] >= 26.0) * 1.0

    is_conifer = conifer > broadleaf
    # Which palm_csd tree type each class becomes is the region's choice
    # (mapping `single_trees`), applied in run().
    crowns["is_conifer"] = is_conifer

    height = np.clip(crowns["HEIGHT_TOP"].to_numpy(dtype=float), 1.31, None)
    root = np.sqrt(np.maximum(height - 1.3, 1e-6))
    with np.errstate(divide="ignore", invalid="ignore"):
        dbh_broadleaf = (1.31811 * root) / (1.0 - 0.151515 * root)
        dbh_conifer = 19.5 / (-np.log(np.maximum(height - 1.3, 1e-6)) + 3.85)
    dbh = np.where(is_conifer, dbh_conifer, dbh_broadleaf)
    dbh = np.clip(np.where(np.isfinite(dbh), dbh, 20.0), 1.0, 200.0)
    crowns["stammdurch"] = dbh / 100.0

    with np.errstate(divide="ignore", invalid="ignore"):
        lai_broadleaf = (0.307279 * dbh**1.803) / (np.pi * (0.078307 * dbh + 0.905689) ** 2)
        lai_conifer = 12.25 * (1.0 - np.exp(-0.1 * dbh)) ** 1.6
    lai = np.where(is_conifer, lai_conifer, lai_broadleaf)
    crowns["LAI"] = np.clip(np.where(np.isfinite(lai), lai, 3.0), 0.0, 12.0)


def _tile_trees(
    cfg: Config, tier: Tier, tile: Tile, ndsm_vrt: Path, ortho_vrt: Path, par: Dict[str, float]
) -> Optional[Any]:
    """Detect the trees whose top falls inside one work tile."""
    import geopandas as gpd
    from rasterio.enums import Resampling
    from rasterio.transform import Affine

    res = tier.resolution
    halo = int(round(par["halo_m"] / res))
    xmin, ymin, xmax, ymax = tile.bounds
    wide = (xmin - halo * res, ymin - halo * res, xmax + halo * res, ymax + halo * res)
    n = int(round((wide[2] - wide[0]) / res))
    transform = Affine.translation(wide[0], wide[3]) * Affine.scale(res, -res)

    chm = read_on_grid(
        ndsm_vrt, transform, (n, n), band=1,
        resampling=Resampling.average, dtype="float32", fill=np.nan,
    )
    if not np.isfinite(chm).any():
        return None
    bands = np.stack(
        [
            read_on_grid(
                ortho_vrt, transform, (n, n), band=b,
                resampling=Resampling.average, dtype="float32", fill=np.nan,
            )
            for b in (1, 2, 3, 4)
        ]
    )
    if not np.isfinite(bands).any():
        return None
    bands = np.nan_to_num(bands, nan=0.0)

    crowns = segment_tile(chm, bands, transform, res, par)
    if crowns is None or crowns.empty:
        return None

    # Represent each tree by its crown centroid and keep only those inside the
    # tile proper, so a crown reaching across the border is not counted twice.
    crowns = crowns.set_crs(cfg.crs, allow_override=True)
    points = crowns.copy()
    points["geometry"] = crowns.geometry.centroid
    inside = (
        (points.geometry.x >= xmin)
        & (points.geometry.x < xmax)
        & (points.geometry.y >= ymin)
        & (points.geometry.y < ymax)
    )
    points = points[inside]
    return points if not points.empty else None


def run(cfg: Config) -> Dict[str, Any]:
    """Detect single trees over the finest tier and write ``trees.gpkg``."""
    import geopandas as gpd
    import pandas as pd

    mode = str(cfg.get("vegetation.trees.mode", "rentz"))
    if mode != "rentz":
        logger.info("vegetation.trees.mode is '%s' -- skipping segmentation", mode)
        return {"mode": mode, "trees": 0}
    kinds = cfg.mapping.single_trees
    if kinds is None:
        raise ConfigError(
            "vegetation.trees.mode is 'rentz', but the region mapping declares no "
            "`single_trees`: which palm_csd tree type a detected conifer or "
            "broadleaf crown becomes is a regional choice."
        )

    tier = cfg.tiers[0]  # tiers are ordered finest first
    if tier.resolution > MAX_USEFUL_RESOLUTION:
        raise RuntimeError(
            f"The finest work tier is {tier.resolution:g} m. Single-tree segmentation "
            f"needs {MAX_USEFUL_RESOLUTION:g} m or finer -- a pixel has to be smaller "
            f"than a crown. Either add a finer domain or set "
            f"vegetation.trees.mode to 'unresolved'."
        )

    tile_dir = cfg.tier_dir(tier) / "tiles"
    ndsm_parts = sorted(tile_dir.glob("*/ndsm.tif"))
    ortho_parts = sorted(tile_dir.glob("*/ortho_leaf_on.tif"))
    if not ndsm_parts:
        raise RuntimeError("No nDSM tiles. Run the `rasters` stage first.")
    if not ortho_parts:
        raise RuntimeError(
            "No stored orthophoto bands. Re-run the `rasters` stage with "
            "output.keep_ortho_bands enabled (the default)."
        )
    vrt_dir = cfg.tier_dir(tier) / "vrt"
    vrt_dir.mkdir(parents=True, exist_ok=True)
    ndsm_vrt = build_vrt(ndsm_parts, vrt_dir / "ndsm_trees.vrt")
    ortho_vrt = build_vrt(ortho_parts, vrt_dir / "ortho_trees.vrt")

    par = _params(cfg)
    tiles = tier.tiles
    workers = int(cfg.runtime.get("workers", 4))
    logger.info(
        "Single-tree segmentation over %d tile(s) at %.0f m, %d worker(s)",
        len(tiles), tier.resolution, workers,
    )

    frames: List[Any] = []
    done = 0
    with gdal_env(int(cfg.runtime.get("gdal_cachemax_mb", 512))):
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_tile_trees, cfg, tier, t, ndsm_vrt, ortho_vrt, par): t
                for t in tiles
            }
            for future in as_completed(futures):
                tile = futures[future]
                try:
                    result = future.result()
                except Exception as exc:  # noqa: BLE001
                    logger.error("Tree segmentation failed for tile %s: %s", tile.key, exc)
                    result = None
                if result is not None and not result.empty:
                    frames.append(result)
                done += 1
                if done % 25 == 0 or done == len(tiles):
                    logger.info("  %d/%d tiles", done, len(tiles))

    if not frames:
        logger.warning("No trees detected over the area of interest")
        return {"mode": mode, "trees": 0}

    trees = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs=cfg.crs)

    # A green roof is not a tree.
    buildings = cfg.workdir / "vector" / "buildings.gpkg"
    if buildings.exists():
        footprints = gpd.read_file(buildings)
        before = len(trees)
        trees = trees.reset_index(drop=True)
        # `inner` rather than masking an `outer` result: a point falling into
        # two overlapping footprints appears twice, so the join is longer than
        # the input and cannot be used as a boolean mask.
        hit = gpd.sjoin(
            trees[["geometry"]], footprints[["geometry"]], how="inner", predicate="within"
        )
        trees = trees.drop(index=hit.index.unique())
        removed = before - len(trees)
        if removed:
            logger.info("Dropped %d detection(s) standing on building footprints", removed)

    is_conifer = trees["is_conifer"].astype(bool)
    trees["tree_type"] = np.where(
        is_conifer, kinds["conifer"]["tree_type"], kinds["broadleaf"]["tree_type"])
    trees["species"] = np.where(
        is_conifer, kinds["conifer"]["species"], kinds["broadleaf"]["species"])

    columns = {
        "HEIGHT_TOP": "tree_height",
        "CROWN_D": "tree_crown_diameter",
        "stammdurch": "tree_trunk_diameter",
        "LAI": "tree_lai",
        "tree_type": "tree_type",
    }
    out = trees[[*columns, "species", "geometry"]].rename(columns=columns)
    out = out[out["tree_height"].between(par["height_min"], 120.0)]
    # palm_csd validates against tree_crown_diameter 0.01..65 m
    out = out[out["tree_crown_diameter"].between(0.01, 65.0)]

    path = cfg.workdir / "vector" / "trees.gpkg"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    out.to_file(path, driver="GPKG")

    conifer = int(is_conifer.loc[out.index].sum())
    broadleaf = len(out) - conifer
    logger.info(
        "Trees: %d detected (%d broadleaf, %d coniferous), median height %.1f m, "
        "median crown %.1f m -> %s",
        len(out), broadleaf, conifer,
        float(out["tree_height"].median()), float(out["tree_crown_diameter"].median()),
        path.name,
    )
    return {
        "mode": mode,
        "trees": len(out),
        "broadleaf": broadleaf,
        "conifer": conifer,
        "path": str(path),
    }
