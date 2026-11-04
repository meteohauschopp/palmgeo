#!/usr/bin/env python3
"""Stage Egypt open-geodata sources for palmgeo (case-study port).

Written for, and measured over, Aswan; nothing in it is tied to that site any
more -- the area and the run's work directory come from the command line.

Copernicus DEM GLO-30 and FABDEM (EPSG:4326) must be reprojected to the target
CRS (EPSG:32636 by default, the CRS sources_eg.yml declares). Sentinel-2 tiles
in Egypt are natively in UTM. OSM data is extracted from the Overpass API in
WGS84 and reprojected.

Run before ``palmgeo vectors`` / ``palmgeo rasters`` / ``palmgeo all``.
``palmgeo plan`` does not need staged data (it only reads index metadata).

Usage::

    STAGE="python repository/ports/eg/stage_sources.py \
        --workdir static_driver/<run>/work --bbox <lon0> <lat0> <lon1> <lat1>"
    $STAGE index      # build GeoJSON indices
    $STAGE dem        # stage FABDEM (terrain)
    $STAGE dsm        # stage Copernicus DEM (surface)
    $STAGE ortho      # find + index Sentinel-2
    $STAGE osm        # extract OSM cadastre+roads
    $STAGE water      # fuse water polygons
    $STAGE buildings --open-buildings-tile <s2 token>   # fuse Open Buildings
    $STAGE all --open-buildings-tile <s2 token>         # everything above, in order

``--workdir`` must be the run's ``project.workdir``: sources_eg.yml finds the
staged files through ``{workdir}/staging``. ``--bbox`` is the area to stage,
in WGS84 lon/lat, and should cover the union of all domains plus palmgeo's
``grid.buffer`` and the 1 km work-tile rounding -- a margin of ~500 m on every
side is ample. FABDEM and Copernicus are staged from the single 1-degree tile
the box lies in; a box across a degree line is refused rather than silently
cut.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import math
import tempfile
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# A stale PROJ_DATA/PROJ_LIB in the environment is dropped by the palmgeo
# package itself (palmgeo/__init__.py) -- importing it first is the fix.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from palmgeo import USER_AGENT  # noqa: E402

logger = logging.getLogger("stage_eg")

# ---------------------------------------------------------------------------
# Run settings. Set from the command line by _configure(); nothing here names
# a site, because the port is for Egypt and the area is the caller's.
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[3]
WORKDIR = Path()
STAGING = Path()
INDEX_DIR = Path()
#: The CRS sources_eg.yml declares its products in (UTM 36N).
DEFAULT_EPSG = 32636
TARGET_EPSG = DEFAULT_EPSG

#: The area to stage, WGS84 (lon0, lat0, lon1, lat1).
BBOX_WGS84: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)

#: The same area in the target CRS: the envelope of BBOX_WGS84.
AOI_TARGET_BOUNDS: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)


def _configure(workdir: Path, bbox: List[float], epsg: int) -> None:
    """Set the run settings above from the command line."""
    global WORKDIR, STAGING, INDEX_DIR, TARGET_EPSG, BBOX_WGS84, AOI_TARGET_BOUNDS
    from pyproj import Transformer

    lon0, lat0, lon1, lat1 = (float(v) for v in bbox)
    if not (-180 <= lon0 < lon1 <= 180 and -90 <= lat0 < lat1 <= 90):
        raise SystemExit(f"--bbox {bbox}: expected lon0 lat0 lon1 lat1 in WGS84 degrees")
    WORKDIR = (workdir if workdir.is_absolute() else Path.cwd() / workdir).resolve()
    STAGING = WORKDIR / "staging"
    INDEX_DIR = STAGING / "indices"
    TARGET_EPSG = int(epsg)
    BBOX_WGS84 = (lon0, lat0, lon1, lat1)
    transformer = Transformer.from_crs(4326, TARGET_EPSG, always_xy=True)
    xs, ys = transformer.transform([lon0, lon1, lon0, lon1], [lat0, lat0, lat1, lat1])
    AOI_TARGET_BOUNDS = (min(xs), min(ys), max(xs), max(ys))


def _hemi(value: int, pos: str, neg: str, width: int) -> str:
    return f"{pos if value >= 0 else neg}{abs(value):0{width}d}"


def _degree_tile() -> Tuple[int, int]:
    """South-west corner (lat, lon) of the one 1-degree tile holding the box."""
    lon0, lat0, lon1, lat1 = BBOX_WGS84
    south, west = math.floor(lat0), math.floor(lon0)
    if math.ceil(lat1) - 1 != south or math.ceil(lon1) - 1 != west:
        raise RuntimeError(
            f"--bbox {BBOX_WGS84} crosses a whole-degree line. FABDEM and "
            f"Copernicus are staged from one 1-degree tile; split the area or "
            f"mosaic the neighbouring tiles by hand."
        )
    return south, west


def _degree_tile_name() -> str:
    """The tile name both FABDEM and this script use, e.g. ``N24E032``."""
    south, west = _degree_tile()
    return _hemi(south, "N", "S", 2) + _hemi(west, "E", "W", 3)


def _ensure_dirs() -> None:
    STAGING.mkdir(parents=True, exist_ok=True)
    INDEX_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# 1. FABDEM (DTM) -- download + reproject. Copernicus DEM (DSM) is staged the
#    same way by stage_dsm(): palmgeo reads without reprojecting, so a
#    geographic source has to be brought into the target CRS first.
# ---------------------------------------------------------------------------
def stage_dem() -> None:
    """Download the FABDEM tile under the box, reproject to the target CRS."""
    import rasterio
    from rasterio.warp import calculate_default_transform, reproject, Resampling

    _ensure_dirs()
    tile_name = _degree_tile_name()
    staged_path = STAGING / f"fabdem_{tile_name}_{TARGET_EPSG}.tif"

    if staged_path.exists():
        logger.info("FABDEM already staged: %s", staged_path)
        return

    # FABDEM groups its 1-degree tiles into 10-degree folders.
    south, west = _degree_tile()
    s10, w10 = 10 * math.floor(south / 10), 10 * math.floor(west / 10)
    group = (_hemi(s10, "N", "S", 2) + _hemi(w10, "E", "W", 3) + "-"
             + _hemi(s10 + 10, "N", "S", 2) + _hemi(w10 + 10, "E", "W", 3))
    tile_url = (
        "https://huggingface.co/buckets/links-ads/fabdem/resolve/tiles/"
        f"{group}_FABDEM_V1-2/{tile_name}_FABDEM_V1-2.tif"
    )
    logger.info("Downloading + reprojecting FABDEM tile %s ...", tile_name)

    # Open remote source
    src_path = f"/vsicurl/{tile_url}"
    with rasterio.open(src_path) as src:
        # Calculate transform for target CRS + bounds
        dst_transform, dst_width, dst_height = calculate_default_transform(
            src.crs, f"EPSG:{TARGET_EPSG}",
            src.width, src.height, *src.bounds,
            dst_width=None, dst_height=None,
            resolution=30.0,
        )
        dst_crs = rasterio.crs.CRS.from_epsg(TARGET_EPSG)
        profile = src.profile.copy()
        profile.update(
            crs=dst_crs,
            transform=dst_transform,
            width=dst_width,
            height=dst_height,
            compress="deflate",
            predictor=2,
        )

        with rasterio.open(staged_path, "w", **profile) as dst:
            reproject(
                source=rasterio.band(src, 1),
                destination=rasterio.band(dst, 1),
                src_transform=src.transform,
                src_crs=src.crs,
                src_nodata=src.nodata,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                dst_nodata=-9999.0,
                resampling=Resampling.bilinear,
            )

    logger.info("FABDEM staged: %s", staged_path)


# ---------------------------------------------------------------------------
# 2. Sentinel-2 orthophoto -- STAC query + index
# ---------------------------------------------------------------------------
def stage_ortho() -> None:
    """Query Earth Search STAC for Sentinel-2 L2A over the box and build index."""
    import requests

    _ensure_dirs()
    index_path = INDEX_DIR / "ortho_eg.json"
    if index_path.exists():
        logger.info("Ortho index already exists: %s", index_path)
        return

    bbox = list(BBOX_WGS84)
    query_param = json.dumps({"eo:cloud_cover": {"lt": 10}})
    url = (
        "https://earth-search.aws.element84.com/v1/search"
        f"?collections=sentinel-2-l2a"
        f"&bbox={bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}"
        f"&limit=10"
        f"&query={query_param}"
    )
    logger.info("Querying Earth Search STAC: %s", url)
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    data = resp.json()

    features = data.get("features", [])
    if not features:
        logger.error("No Sentinel-2 scenes found over the box with <10%% cloud cover")
        logger.info("Retrying with <30%% cloud cover ...")
        query_param = json.dumps({"eo:cloud_cover": {"lt": 30}})
        url = (
            "https://earth-search.aws.element84.com/v1/search"
            f"?collections=sentinel-2-l2a"
            f"&bbox={bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}"
            f"&limit=10"
            f"&query={query_param}"
        )
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        data = resp.json()
        features = data.get("features", [])

    if not features:
        logger.error("No Sentinel-2 scenes found. Cannot build ortho index.")
        return

    # Pick the scene with lowest cloud cover
    features.sort(key=lambda f: f.get("properties", {}).get("eo:cloud_cover", 100))
    best = features[0]
    props = best.get("properties", {})
    geom = best.get("geometry", {})
    datetime_str = props.get("datetime", "")
    cloud_cover = props.get("eo:cloud_cover", "?")
    logger.info(
        "Selected scene: %s (cloud %.1f%%, date %s)",
        props.get("id", "?"), cloud_cover, datetime_str[:10],
    )

    # Get the red, green, blue, nir asset URLs
    assets = best.get("assets", {})
    # Earth Search v1 uses lowercase asset keys: red, green, blue, nir
    band_map = {"red": "red", "green": "green", "blue": "blue", "nir": "nir"}
    asset_urls = {}
    for role, band_key in band_map.items():
        asset = assets.get(band_key, {})
        href = asset.get("href", "")
        if href:
            asset_urls[role] = href

    if not asset_urls:
        logger.error("No band assets found in STAC item")
        return

    # Download the four 10 m bands, reproject to the target CRS, and merge into a
    # single 4-band RGBN COG covering the AOI (band order red, green, blue,
    # nir -- matching sources_eg.yml's `bands: {red:1,green:2,blue:3,nir:4}`).
    # palmgeo's process_tile reads red/green/blue/nir from ONE file, so this
    # merging is required before palmgeo ever sees it (same single-mosaic
    # approach as the NYC port).
    import rasterio
    from rasterio.warp import calculate_default_transform, reproject, Resampling
    from rasterio.merge import merge as rasterio_merge

    out_cog = STAGING / f"sentinel2_rgbn_{TARGET_EPSG}.tif"
    if not out_cog.exists():
        # Earth Search v1 uses lowercase asset keys: red, green, blue, nir
        band_map_merge = {"red": "red", "green": "green", "blue": "blue", "nir": "nir"}
        # 1) Stream each band COG via /vsicurl, reproject to the target CRS
        band_datasets = []
        for band_role, band_key in band_map_merge.items():
            href = assets.get(band_key, {}).get("href", "")
            if not href:
                raise RuntimeError(f"missing band asset {band_key}")
            src_path = f"/vsicurl/{href}"
            with rasterio.open(src_path) as src:
                dst_crs = rasterio.crs.CRS.from_epsg(TARGET_EPSG)
                transform, width, height = calculate_default_transform(
                    src.crs, dst_crs, src.width, src.height, *src.bounds,
                    resolution=10.0,
                )
                profile = src.profile.copy()
                profile.update(
                    crs=dst_crs, transform=transform,
                    width=width, height=height,
                )
                # Write reprojected single-band to a temp file
                tmp = STAGING / f"_tmp_{band_role}.tif"
                with rasterio.open(tmp, "w", **profile) as dst:
                    reproject(
                        source=rasterio.band(src, 1),
                        destination=rasterio.band(dst, 1),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        src_nodata=src.nodata,
                        dst_transform=transform,
                        dst_crs=dst_crs,
                        dst_nodata=0,
                        resampling=Resampling.bilinear,
                    )
                band_datasets.append(rasterio.open(tmp))

        # 2) Stack into a 4-band RGBN COG
        count = len(band_datasets)
        profile = band_datasets[0].profile.copy()
        profile.update(
            count=count, compress="deflate", predictor=2,
            width=band_datasets[0].width, height=band_datasets[0].height,
        )
        with rasterio.open(out_cog, "w", **profile) as dst:
            for i, ds in enumerate(band_datasets, start=1):
                dst.write(ds.read(1), i)
                dst.update_tags(i, description=list(band_map_merge.keys())[i-1])

        # Cleanup temp files
        for ds in band_datasets:
            ds.close()
        for band_role in band_map_merge:
            tmp = STAGING / f"_tmp_{band_role}.tif"
            if tmp.exists():
                tmp.unlink()

        logger.info("Ortho RGBN COG written: %s", out_cog)

    # Build a single-feature GeoJSON index pointing at the merged COG.
    from shapely.geometry import mapping, box as shapely_box
    from pyproj import Transformer

    xmin, ymin, xmax, ymax = AOI_TARGET_BOUNDS

    feature = {
        "type": "Feature",
        "geometry": mapping(shapely_box(xmin, ymin, xmax, ymax)),
        "properties": {
            "tile_id": props.get("id", "sentinel2"),
            "datetime": datetime_str,
            "cloud_cover": cloud_cover,
            "staged_path": str(out_cog),
        },
    }

    index = {
        "type": "FeatureCollection",
        "features": [feature],
    }
    INDEX_DIR.mkdir(parents=True, exist_ok=True)
    index_path.write_text(json.dumps(index, indent=2))
    logger.info("Ortho index written: %s (%d scene)", index_path, 1)
    logger.info("  staged RGBN COG: %s", out_cog)


# ---------------------------------------------------------------------------
# helper: palmgeo's fetch_cadastre/fetch_roads call cache.unzip(archive) and
# then rglob for *.gpkg inside the zip (steps/tvectors.py). The staged GPKG
# must therefore be delivered as a zip containing exactly it.
# ---------------------------------------------------------------------------
def _zip_gpkg(gpkg: Path) -> Path:
    """Bundle a single GeoPackage into a zip beside it and return the zip path.

    The zip is named ``<stem>_files.zip`` (not ``<name>.zip``) so that
    ``cache.unzip(archive)`` creates a directory ``<stem>_files`` that does
    not collide with the ``.gpkg`` file itself.
    """
    zip_path = gpkg.with_name(gpkg.stem + "_files.zip")
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(gpkg, arcname=gpkg.name)
    logger.info("Zipped: %s", zip_path)
    return zip_path


# ---------------------------------------------------------------------------
# 3. OSM cadastre (landuse + buildings) + topology (highways)
# ---------------------------------------------------------------------------
def _overpass_query(bbox: Tuple[float, float, float, float], query: str) -> Dict:
    """Run an Overpass API query and return the JSON result."""
    import time
    import requests

    s, w, n, e = bbox[1], bbox[0], bbox[3], bbox[2]
    # `out geom` (NOT `out body;>;out skel qt;`) -- two reasons, both of which
    # cost real features in the 2026-09-07 build:
    #   1. `>` recurses down into every relation's member ways, and those
    #      members carry no tags of their own. The old converter turned each
    #      one into a standalone polygon, so the Nile relation arrived as 61
    #      untagged fragments plus 3 real features.
    #   2. `out geom` attaches the coordinate list to each way AND to each
    #      relation member, which is what the multipolygon stitcher below
    #      needs to reassemble an outer boundary that OSM has split across
    #      many ways -- exactly how the Nile riverbank is mapped.
    full_query = f"[out:json][timeout:300];({query});out geom;"
    url = "https://overpass-api.de/api/interpreter"
    logger.info("Overpass query (%d chars) ...", len(full_query))

    for attempt in range(4):
        resp = requests.post(
            url,
            data={"data": full_query},
            timeout=600,
            headers={"User-Agent": USER_AGENT},
        )
        if resp.status_code == 200:
            return resp.json()
        if resp.status_code in (429, 400, 504):
            wait = 30 * (attempt + 1)
            logger.warning("Overpass %d (attempt %d), waiting %ds ...", resp.status_code, attempt + 1, wait)
            time.sleep(wait)
        else:
            resp.raise_for_status()
    resp.raise_for_status()
    return {}


def _stitch_rings(members: List[Dict], role: str) -> List[Any]:
    """Reassemble a multipolygon relation's member ways into closed rings.

    An OSM multipolygon's outer boundary is routinely SPLIT across many open
    ways; each member on its own is a line, not a ring. The Nile riverbank
    through Aswan is one such relation. Treating the members as independent
    polygons -- or, as the previous version did, keeping only the FIRST one
    and closing it by force -- fabricates a shape that is not the river:
    it produced a 9.4 km2 blob spanning 18 km of northing, against the
    11.5 km2 the properly stitched relation actually covers.

    linemerge joins the members end to end; polygonize closes what merges
    into rings. Members that do not close are dropped rather than forced.
    """
    from shapely.geometry import LineString
    from shapely.ops import linemerge, polygonize

    lines = [
        LineString([(q["lon"], q["lat"]) for q in m["geometry"]])
        for m in members
        if m.get("role", "") == role and len(m.get("geometry") or []) >= 2
    ]
    if not lines:
        return []
    merged = linemerge(lines)
    return list(polygonize(merged if merged.geom_type == "MultiLineString" else [merged]))


def _overpass_to_features(data: Dict, require_tags: bool = True) -> List[Dict]:
    """Convert Overpass ``out geom`` JSON to GeoJSON *polygon* features (WGS84).

    ``way`` elements carry their coordinate ring inline; ``relation`` elements
    carry one ring per member, stitched by :func:`_stitch_rings`. Inner rings
    are subtracted, so an island in the Nile is not water.

    ``require_tags`` skips elements with no tags at all. Under ``out geom`` a
    relation's members are not emitted as separate elements, so this mostly
    guards against stray untagged ways; it is the second half of the fix for
    the 61 phantom water polygons.
    """
    from shapely.geometry import Polygon, mapping
    from shapely.ops import unary_union

    features: List[Dict] = []

    def _emit(geom, el, tags):
        if geom.is_empty:
            return
        if not geom.is_valid:
            geom = geom.buffer(0)
        if geom.is_empty:
            return
        parts = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
        for part in parts:
            if part.is_empty or part.area <= 0:
                continue
            features.append({
                "type": "Feature",
                "geometry": mapping(part),
                "properties": {
                    "osm_id": el["id"],
                    "osm_type": el["type"],
                    **tags,
                },
            })

    for el in data.get("elements", []):
        tags = el.get("tags", {}) or {}
        if require_tags and not tags:
            continue

        if el["type"] == "way":
            ring = [(q["lon"], q["lat"]) for q in (el.get("geometry") or [])]
            # An area needs a closed ring. Ways that are genuinely linear
            # (a river centreline, a road) are handled by
            # _overpass_to_lines instead, not silently closed into a
            # polygon the way the previous version did.
            if len(ring) < 4:
                continue
            if ring[0] != ring[-1]:
                ring.append(ring[0])
            if len(ring) < 4:
                continue
            try:
                _emit(Polygon(ring), el, tags)
            except Exception:
                continue

        elif el["type"] == "relation":
            members = el.get("members", []) or []
            outer = _stitch_rings(members, "outer")
            if not outer:
                continue
            try:
                geom = unary_union(outer)
                inner = _stitch_rings(members, "inner")
                if inner:
                    geom = geom.difference(unary_union(inner))
                _emit(geom, el, tags)
            except Exception:
                continue

    return features


def _overpass_to_lines(data: Dict) -> List[Dict]:
    """Convert Overpass ``out geom`` JSON to GeoJSON *line* features (WGS84).

    Used for the highway network and for water centrelines. Reads the inline
    ``geometry`` list, so it needs no separate node table -- the previous
    highway loop rebuilt a nodes dict from the response for every single way,
    which was both O(n^2) and incompatible with ``out geom``.
    """
    from shapely.geometry import LineString, mapping

    features: List[Dict] = []
    for el in data.get("elements", []):
        if el["type"] != "way":
            continue
        coords = [(q["lon"], q["lat"]) for q in (el.get("geometry") or [])]
        if len(coords) < 2:
            continue
        try:
            line = LineString(coords)
            if not line.is_valid or line.is_empty:
                continue
        except Exception:
            continue
        features.append({
            "type": "Feature",
            "geometry": mapping(line),
            "properties": {
                "osm_id": el["id"],
                "osm_type": el["type"],
                **(el.get("tags", {}) or {}),
            },
        })
    return features


def _drop_duplicate_osm(gdf, what: str):
    """Drop repeated OSM elements, keeping the first of each ``osm_id``."""
    if "osm_id" not in gdf.columns:
        return gdf
    before = len(gdf)
    gdf = gdf.drop_duplicates(subset=["osm_type", "osm_id"], keep="first").reset_index(drop=True)
    if len(gdf) != before:
        logger.info("  %s: %d duplicate rows dropped (%d unique)", what, before - len(gdf), len(gdf))
    return gdf


# Storey height for Aswan's building stock. Egyptian residential concrete-frame
# construction runs ~3.0 m floor to floor; the ground floor of a shop-house is
# taller. 3.2 m is the value the Egyptian residential building code assumes for
# a typical apartment storey, and it is what `building:levels` is multiplied by
# where OSM records levels but no explicit height.
STOREY_HEIGHT_M = 3.2


def _normalise_building_height(tags: Dict) -> None:
    """Fold OSM height tags into numeric ``osm_height`` / ``osm_levels``.

    OSM records building height three ways -- ``height``, ``building:height``
    and ``building:levels`` -- with units that may be bare metres, a "12 m"
    string, or (rarely) feet. palmgeo reads a single numeric column, so the
    variants are resolved here rather than in the pipeline.

    COVERAGE IS THE POINT: across the Aswan AOI only 32 of 4 025 buildings
    carry `building:levels` and exactly 1 carries `height` -- 0.8 %, not the
    "~40 %" that sources_eg.yml claimed before this was measured. These tags
    are therefore a garnish on the height model, not its backbone; the
    `buildings.height_default` typological prior does the actual work.
    """
    def _num(value) -> Optional[float]:
        if value is None:
            return None
        text = str(value).strip().lower().replace(",", ".")
        if not text:
            return None
        factor = 1.0
        for unit, mult in (("m", 1.0), ("meter", 1.0), ("metre", 1.0),
                           ("ft", 0.3048), ("'", 0.3048)):
            if text.endswith(unit):
                text, factor = text[: -len(unit)].strip(), mult
                break
        try:
            return float(text) * factor
        except ValueError:
            return None

    height = _num(tags.get("height")) or _num(tags.get("building:height"))
    levels = _num(tags.get("building:levels"))
    if height is not None and 2.0 <= height <= 200.0:
        tags["osm_height"] = height
    if levels is not None and 0 < levels <= 60:
        tags["osm_levels"] = levels
        if "osm_height" not in tags:
            # roof:levels are habitable storeys too where tagged
            roof = _num(tags.get("roof:levels")) or 0.0
            tags["osm_height"] = (levels + roof) * STOREY_HEIGHT_M


def stage_osm() -> None:
    """Extract landuse, buildings, and highways from OSM Overpass API."""
    import geopandas as gpd
    from shapely.geometry import shape, mapping
    from pyproj import Transformer

    _ensure_dirs()
    cadastre_index_path = INDEX_DIR / "cadastre_eg.json"
    topo_index_path = INDEX_DIR / "topo_eg.json"

    transformer = Transformer.from_crs(4326, TARGET_EPSG, always_xy=True)

    bbox = BBOX_WGS84
    # Overpass expects (south,west,north,east) = (lat_min,lon_min,lat_max,lon_max)
    s, w, n, e = bbox[1], bbox[0], bbox[3], bbox[2]
    BBOX = f"{s},{w},{n},{e}"

    # --- Landuse ---
    logger.info("Querying OSM landuse ...")
    lu_q = "\n".join([
        f'way["landuse"]({BBOX});',
        f'relation["landuse"]({BBOX});',
        f'way["natural"~"wood|scrub|grassland|wetland|sand|bare_rock|heath|desert"]({BBOX});',
        f'relation["natural"~"wood|scrub|grassland|wetland|sand|bare_rock|heath|desert"]({BBOX});',
        f'way["leisure"~"park|garden|pitch|golf_course|stadium"]({BBOX});',
        f'relation["leisure"~"park|garden|pitch|golf_course|stadium"]({BBOX});',
        f'way["amenity"~"parking|school|hospital|university|place_of_worship"]({BBOX});',
        f'relation["amenity"~"parking|school|hospital|university|place_of_worship"]({BBOX});',
    ])
    lu_data = _overpass_query(bbox, lu_q)
    lu_features = _overpass_to_features(lu_data)
    logger.info("  %d landuse/natural/leisure/amenity polygons", len(lu_features))

    # Compose composite keys: "landuse:residential", "natural:wood", etc.
    for f in lu_features:
        tags = f["properties"]
        key = None
        for prefix in ("landuse", "natural", "leisure", "amenity"):
            if prefix in tags:
                key = f"{prefix}:{tags[prefix]}"
                break
        tags["osm_key"] = key or "unknown"
        # Reproject geometry to target CRS
        geom = shape(f["geometry"])
        geom_proj = gpd.GeoSeries([geom], crs="EPSG:4326").to_crs(f"EPSG:{TARGET_EPSG}").iloc[0]
        f["geometry"] = mapping(geom_proj)

    # --- Buildings ---
    logger.info("Querying OSM buildings ...")
    import time; time.sleep(15)  # rate-limit courtesy
    bld_q = f'way["building"]({BBOX});\nrelation["building"]({BBOX});'
    bld_data = _overpass_query(bbox, bld_q)
    bld_features = _overpass_to_features(bld_data)
    logger.info("  %d building polygons", len(bld_features))
    for f in bld_features:
        _normalise_building_height(f["properties"])
        geom = shape(f["geometry"])
        geom_proj = gpd.GeoSeries([geom], crs="EPSG:4326").to_crs(f"EPSG:{TARGET_EPSG}").iloc[0]
        f["geometry"] = mapping(geom_proj)

    # Write fused cadastre GeoPackage.
    #
    # The GPKG is DELETED first and each layer written with append=False.
    # The previous version passed `append=bool(lu_features)` when writing the
    # building layer, so every re-run of `ports/eg/stage_sources.py osm` appended a
    # second copy of every building to the layer that was already there. The
    # 2026-09-07 build ran it four times: 16 100 building rows for 4 025
    # distinct footprints, every one of them stacked 4 deep. Duplicate
    # coincident footprints do not change the rasterised result, but they make
    # every building count in the logs and the README a 4x overstatement.
    cadastre_path = STAGING / "cadastre_eg.gpkg"
    if cadastre_path.exists():
        cadastre_path.unlink()
    if lu_features:
        gdf_lu = gpd.GeoDataFrame.from_features(lu_features, crs=f"EPSG:{TARGET_EPSG}")
        gdf_lu = _drop_duplicate_osm(gdf_lu, "landuse")
        gdf_lu.to_file(cadastre_path, layer="landuse", driver="GPKG")
    if bld_features:
        gdf_bld = gpd.GeoDataFrame.from_features(bld_features, crs=f"EPSG:{TARGET_EPSG}")
        gdf_bld = _drop_duplicate_osm(gdf_bld, "buildings")
        gdf_bld.to_file(cadastre_path, layer="buildings", driver="GPKG")
    logger.info("Cadastre written: %s", cadastre_path)
    cadastre_zip = _zip_gpkg(cadastre_path) if cadastre_path.exists() else None

    # Cadastre index
    from shapely.geometry import box as shapely_box
    xmin, ymin = transformer.transform(bbox[0], bbox[1])
    xmax, ymax = transformer.transform(bbox[2], bbox[3])
    cad_index = {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "geometry": mapping(shapely_box(xmin, ymin, xmax, ymax)),
            "properties": {
                "tile_id": "eg_cadastre",
                "staged_path": str(cadastre_zip or cadastre_path),
            },
        }],
    }
    INDEX_DIR.mkdir(parents=True, exist_ok=True)
    cadastre_index_path.write_text(json.dumps(cad_index, indent=2))
    logger.info("Cadastre index: %s", cadastre_index_path)

    # --- Highways ---
    logger.info("Querying OSM highways ...")
    time.sleep(15)  # rate-limit courtesy
    hw_q = f'way["highway"]({BBOX});'
    hw_data = _overpass_query(bbox, hw_q)
    hw_features = []
    for f in _overpass_to_lines(hw_data):
        tags = f["properties"]
        geom = shape(f["geometry"])
        geom_proj = gpd.GeoSeries([geom], crs="EPSG:4326").to_crs(f"EPSG:{TARGET_EPSG}").iloc[0]
        hw_features.append({
            "type": "Feature",
            "geometry": mapping(geom_proj),
            "properties": {
                "osm_id": tags.get("osm_id"),
                "highway": tags.get("highway", "unclassified"),
                "name": tags.get("name", ""),
            },
        })
    logger.info("  %d highway segments", len(hw_features))

    if hw_features:
        roads_dir = STAGING / "topo_eg_shp"
        roads_dir.mkdir(parents=True, exist_ok=True)
        roads_base = roads_dir / "highways.shp"
        gdf_hw = gpd.GeoDataFrame.from_features(hw_features, crs=f"EPSG:{TARGET_EPSG}")
        gdf_hw.to_file(str(roads_base), driver="ESRI Shapefile")
        # fetch_roads() rglobs for *.shp inside a zip
        import shutil
        roads_zip = STAGING / "topo_eg.shp.zip"
        if roads_zip.exists():
            roads_zip.unlink()
        with zipfile.ZipFile(roads_zip, "w", zipfile.ZIP_DEFLATED) as zf:
            for side in roads_dir.glob("highways.*"):
                zf.write(side, arcname=side.name)
        shutil.rmtree(roads_dir, ignore_errors=True)
        logger.info("Roads written: %s", roads_zip)

        topo_index = {
            "type": "FeatureCollection",
            "features": [{
                "type": "Feature",
                "geometry": mapping(shapely_box(xmin, ymin, xmax, ymax)),
                "properties": {
                    "tile_id": "eg_topo",
                    "staged_path": str(roads_zip),
                },
            }],
        }
        topo_index_path.write_text(json.dumps(topo_index, indent=2))
        logger.info("Topo index: %s", topo_index_path)


# ---------------------------------------------------------------------------
# 4. Water polygons
# ---------------------------------------------------------------------------
# OSM `water=*` / `waterway=*` subtags -> the synthetic land-use codes that
# mapping_eg.yml turns into PALM water_type values. The distinction is not
# cosmetic: PALM gives river, lake and pond different default temperatures and
# roughness, and the Nile is emphatically a river.
_WATER_CODES = {
    # moving water
    "river": "WTR_RIV", "canal": "WTR_RIV", "stream": "WTR_RIV",
    "ditch": "WTR_RIV", "fairway": "WTR_RIV", "riverbank": "WTR_RIV",
    # small standing or artificial bodies
    "pond": "WTR_POND", "basin": "WTR_POND", "reservoir_covered": "WTR_POND",
    "wastewater": "WTR_POND", "reflecting_pool": "WTR_POND",
    "fountain": "WTR_POND", "moat": "WTR_POND", "salt_pond": "WTR_POND",
    # genuine lakes and reservoirs
    "lake": "WTR_LAKE", "reservoir": "WTR_LAKE", "lagoon": "WTR_LAKE",
    "oxbow": "WTR_LAKE",
}

# Fallback half-widths for river/canal centrelines that OSM maps as a line
# with no accompanying riverbank polygon. Only used where a centreline does
# NOT already fall inside a mapped water polygon, so the Nile itself (which
# has a proper riverbank relation) is never buffered.
_CENTRELINE_HALF_WIDTH_M = {"river": 15.0, "canal": 5.0, "stream": 2.0}


def _water_code(tags: Dict) -> str:
    """Classify one OSM water feature into a WTR_* code."""
    sub = str(tags.get("water", "") or "").lower()
    if sub in _WATER_CODES:
        return _WATER_CODES[sub]
    way = str(tags.get("waterway", "") or "").lower()
    if way in _WATER_CODES:
        return _WATER_CODES[way]
    if str(tags.get("landuse", "") or "").lower() == "reservoir":
        return "WTR_LAKE"
    # `natural=water` with no subtag. Through Aswan the only water body of any
    # size is the Nile, so river is the right default here; a blanket "lake"
    # would give the river lake thermodynamics.
    return "WTR_RIV"


def stage_water() -> None:
    """Extract OSM water polygons and FUSE them into the cadastre landuse layer.

    This must run after :func:`stage_osm`, because it appends to the landuse
    layer that function writes and then re-zips and re-indexes the cadastre.

    WHY THE FUSION MATTERS: palmgeo reads water through the cadastre's
    ``landuse`` layer, keyed on ``osm_key`` and resolved against
    mapping_eg.yml's WTR_* rules -- there is no separate "water" role in the
    source registry. The previous version wrote a standalone
    ``water_eg.gpkg`` that no stage of the pipeline ever opened, which is why
    the 2026-09-07 driver had water_type nodata over 100 % of the domain and
    the Nile -- 15.3 % of the domain by area -- was simply absent.
    """
    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import shape, mapping, box as shapely_box
    from shapely.ops import unary_union

    _ensure_dirs()
    bbox = BBOX_WGS84
    _s, _w, _n, _e = bbox[1], bbox[0], bbox[3], bbox[2]
    _BBOX = f"{_s},{_w},{_n},{_e}"

    # --- water areas ---
    logger.info("Querying OSM water areas ...")
    water_q = "\n".join([
        f'way["natural"="water"]({_BBOX});',
        f'relation["natural"="water"]({_BBOX});',
        f'way["waterway"~"^(riverbank|dock)$"]({_BBOX});',
        f'relation["waterway"~"^(riverbank|dock)$"]({_BBOX});',
        f'way["landuse"="reservoir"]({_BBOX});',
        f'relation["landuse"="reservoir"]({_BBOX});',
    ])
    area_features = _overpass_to_features(_overpass_query(bbox, water_q))
    logger.info("  %d water area polygons", len(area_features))

    geoms, codes = [], []
    for f in area_features:
        geoms.append(shape(f["geometry"]))
        codes.append(_water_code(f["properties"]))

    # --- water centrelines (fallback only) ---
    import time
    time.sleep(15)  # rate-limit courtesy
    logger.info("Querying OSM water centrelines ...")
    line_q = f'way["waterway"~"^(river|canal|stream)$"]({_BBOX});'
    line_features = _overpass_to_lines(_overpass_query(bbox, line_q))
    logger.info("  %d water centrelines", len(line_features))

    mapped = unary_union(geoms) if geoms else None
    n_buffered = 0
    for f in line_features:
        line = shape(f["geometry"])
        kind = str(f["properties"].get("waterway", "")).lower()
        # Skip any centreline already covered by a mapped polygon -- otherwise
        # the Nile's own centreline would be buffered on top of its riverbank
        # relation, which changes nothing but wastes work and muddies the log.
        if mapped is not None and line.within(mapped.buffer(1e-6)):
            continue
        half = _CENTRELINE_HALF_WIDTH_M.get(kind)
        if not half:
            continue
        # Buffer in metres, so project first.
        line_m = gpd.GeoSeries([line], crs="EPSG:4326").to_crs(f"EPSG:{TARGET_EPSG}").iloc[0]
        poly_m = line_m.buffer(half)
        poly = gpd.GeoSeries([poly_m], crs=f"EPSG:{TARGET_EPSG}").to_crs("EPSG:4326").iloc[0]
        geoms.append(poly)
        codes.append(_WATER_CODES.get(kind, "WTR_RIV"))
        n_buffered += 1
    if n_buffered:
        logger.info("  %d centrelines buffered (no riverbank polygon mapped)", n_buffered)

    if not geoms:
        logger.warning("No water features found -- water_type will be empty")
        return

    water = gpd.GeoDataFrame({"osm_key": codes}, geometry=geoms, crs="EPSG:4326")
    water = water.to_crs(f"EPSG:{TARGET_EPSG}")
    water = water[water.geometry.notna() & ~water.geometry.is_empty]
    water["geometry"] = water.geometry.buffer(0)

    # Clip to the staged AOI so the Nile's 11.5 km2 relation -- which runs far
    # north and south of the domain -- does not bloat the cadastre.
    aoi = shapely_box(*AOI_TARGET_BOUNDS)
    water = gpd.clip(water, aoi)
    water = water[~water.geometry.is_empty & water.geometry.notna()].reset_index(drop=True)
    if water.empty:
        logger.warning("No water features inside the AOI -- water_type will be empty")
        return

    area_km2 = float(water.geometry.area.sum()) / 1e6
    logger.info("Water inside AOI: %d polygons, %.2f km2", len(water), area_km2)
    for code, sub in water.groupby("osm_key"):
        logger.info("  %-9s %2d polygons, %.2f km2", code, len(sub),
                    float(sub.geometry.area.sum()) / 1e6)

    # Keep the standalone file for inspection / QGIS, but it is no longer what
    # the pipeline reads.
    water_path = STAGING / "water_eg.gpkg"
    if water_path.exists():
        water_path.unlink()
    water.to_file(water_path, layer="water", driver="GPKG")
    logger.info("Water written (inspection copy): %s", water_path)

    # --- fuse into the cadastre landuse layer -------------------------------
    cadastre_path = STAGING / "cadastre_eg.gpkg"
    if not cadastre_path.exists():
        logger.error("Cadastre %s missing -- run `ports/eg/stage_sources.py osm` first; "
                     "water cannot be fused and will NOT reach the driver", cadastre_path)
        return

    landuse = gpd.read_file(cadastre_path, layer="landuse")
    landuse = landuse[landuse.get("osm_key", "").astype(str) != "WTR"]      # legacy rows
    landuse = landuse[~landuse["osm_key"].astype(str).str.startswith("WTR_")]
    buildings = gpd.read_file(cadastre_path, layer="buildings")

    fused = pd.concat(
        [landuse, water.to_crs(landuse.crs if landuse.crs else f"EPSG:{TARGET_EPSG}")],
        ignore_index=True,
    )
    fused = gpd.GeoDataFrame(fused, geometry="geometry",
                             crs=landuse.crs or f"EPSG:{TARGET_EPSG}")

    cadastre_path.unlink()
    fused.to_file(cadastre_path, layer="landuse", driver="GPKG")
    buildings.to_file(cadastre_path, layer="buildings", driver="GPKG")
    logger.info("Cadastre landuse re-written with %d water polygons fused in "
                "(%d landuse + %d water rows)", len(water), len(landuse), len(water))

    # The pipeline reads the zip, not the .gpkg, so it must be rebuilt too --
    # and the index rewritten to point at the new zip.
    cadastre_zip = _zip_gpkg(cadastre_path)
    from shapely.geometry import mapping as _mapping
    from pyproj import Transformer
    transformer = Transformer.from_crs(4326, TARGET_EPSG, always_xy=True)
    xmin, ymin = transformer.transform(bbox[0], bbox[1])
    xmax, ymax = transformer.transform(bbox[2], bbox[3])
    (INDEX_DIR / "cadastre_eg.json").write_text(json.dumps({
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "geometry": _mapping(shapely_box(xmin, ymin, xmax, ymax)),
            "properties": {"tile_id": "eg_cadastre", "staged_path": str(cadastre_zip)},
        }],
    }, indent=2))
    logger.info("Cadastre index refreshed: %s", INDEX_DIR / "cadastre_eg.json")


# ---------------------------------------------------------------------------
# 4b. Google Open Buildings -- ML building footprints
# ---------------------------------------------------------------------------
# WHY THIS SOURCE EXISTS IN AN OTHERWISE OSM-ONLY PORT.
#
# OSM's building coverage over central Aswan is drastically incomplete, and
# the gap is not visible from the totals: the 4 x 4 km root domain contains
# 322 km of mapped OSM road but only 591 mapped OSM buildings covering 1.4 %
# of it. The streets of the old city are all there; the buildings between
# them are not. (The southern districts outside the domain are mapped almost
# completely -- 5 571 buildings in the best-covered 4 km window -- so this is
# uneven mapping effort, not a domain that misses the city. OSM's own
# `place=city` node for Aswan sits inside the domain.)
#
# Google Open Buildings v3 is an ML footprint extraction over Africa and
# South Asia, published as S2-cell CSVs under CC BY 4.0. Over this domain it
# yields 35 123 footprints covering 25.4 % -- a plan-area fraction in the
# range expected of a dense Egyptian city core.
#
# VALIDATED 2026-09-08 against the OSM footprints that do exist:
#   - 68 % of OSM buildings in the domain are >50 % covered by a v3 footprint
#     at confidence >= 0.65, so the two sources broadly agree where both see.
#   - 0.5 % of v3 footprints fall on water (feluccas, river islands). They
#     are dropped against the staged water mask rather than left to become
#     buildings standing in the Nile.
#
# WHAT IT DOES NOT CARRY: height, building function, or an ID stable across
# versions. Fusing it therefore buys realistic PLAN FORM -- street canyons,
# block structure, frontage -- at the cost of making building height almost
# entirely a typological prior. That trade is the right one for an urban
# microclimate run, where the obstacle field dominates, but it must be stated
# in the driver's README rather than implied.
#: Open Buildings v3 is published per S2 level-4 cell; the cell is passed as
#: --open-buildings-tile (its token, the file-name prefix on the v3 download
#: page -- "143" for Aswan).
GOB_TILE_URL = ("https://storage.googleapis.com/open-buildings-data/v3/"
                "polygons_s2_level_4_gzip/{token}_buildings.csv.gz")
GOB_TILE_ID = ""

# Confidence floor. Google publishes v3 down to ~0.65 and recommends >= 0.75
# for precision-critical work. 0.70 is the balance this port takes: recall
# matters more than precision for an obstacle field (an absent block distorts
# the flow far more than a spurious 40 m2 shed), while 0.70 still sheds the
# noisiest tail. Measured over the root domain: 0.65 -> 25.4 % plan area,
# 0.70 -> 22.7 %, 0.75 -> 18.0 %.
GOB_MIN_CONFIDENCE = 0.70

# Below this a detection is more likely a shed, a wall or a noise blob than a
# building that shapes the wind field, and at 2 m grid spacing it cannot be
# resolved anyway.
GOB_MIN_AREA_M2 = 15.0


def stage_buildings() -> None:
    """Fuse Google Open Buildings v3 footprints into the cadastre.

    Must run after :func:`stage_osm` (which creates the cadastre) and after
    :func:`stage_water` (whose water mask removes footprints on the Nile).

    FUSION RULE: OSM wins. An Open Buildings footprint that overlaps any OSM
    footprint is dropped, so OSM's surveyed geometry, its ``building=*``
    function and its tagged heights all survive intact; Open Buildings fills
    only the (very large) gaps between them.
    """
    import csv
    import gzip

    import geopandas as gpd
    import pandas as pd
    from shapely import wkt
    from shapely.geometry import box as shapely_box
    from shapely.ops import unary_union

    _ensure_dirs()
    cadastre_path = STAGING / "cadastre_eg.gpkg"
    if not cadastre_path.exists():
        logger.error("Cadastre %s missing -- run `ports/eg/stage_sources.py osm` first",
                     cadastre_path)
        return

    if not GOB_TILE_ID:
        raise RuntimeError(
            "`buildings` needs --open-buildings-tile <s2 token>: the S2 level-4 "
            "cell of Google Open Buildings v3 that holds the box."
        )

    # --- download (cached) ---
    tile = STAGING / f"open_buildings_v3_{GOB_TILE_ID}.csv.gz"
    if tile.exists() and tile.stat().st_size > 1_000_000:
        logger.info("Open Buildings tile already staged: %s (%.0f MB)",
                    tile, tile.stat().st_size / 1e6)
    else:
        logger.info("Downloading Open Buildings tile %s (~70 MB) ...", GOB_TILE_ID)
        import requests
        with requests.get(GOB_TILE_URL.format(token=GOB_TILE_ID), stream=True, timeout=900,
                          headers={"User-Agent": USER_AGENT}) as resp:
            resp.raise_for_status()
            with open(tile, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
        logger.info("  staged %s (%.0f MB)", tile, tile.stat().st_size / 1e6)

    # --- filter to the AOI while streaming (the tile covers ~500 x 500 km) ---
    lon0, lat0, lon1, lat1 = BBOX_WGS84
    geoms, confs = [], []
    n_read = n_lowconf = n_small = 0
    with gzip.open(tile, "rt") as fh:
        for rec in csv.DictReader(fh):
            lon = float(rec["longitude"])
            lat = float(rec["latitude"])
            if not (lon0 <= lon <= lon1 and lat0 <= lat <= lat1):
                continue
            n_read += 1
            if float(rec["confidence"]) < GOB_MIN_CONFIDENCE:
                n_lowconf += 1
                continue
            if float(rec["area_in_meters"]) < GOB_MIN_AREA_M2:
                n_small += 1
                continue
            try:
                geoms.append(wkt.loads(rec["geometry"]))
            except Exception:
                continue
            confs.append(float(rec["confidence"]))
    logger.info("Open Buildings in AOI: %d read, %d below confidence %.2f, "
                "%d below %.0f m2 -> %d kept",
                n_read, n_lowconf, GOB_MIN_CONFIDENCE, n_small, GOB_MIN_AREA_M2,
                len(geoms))
    if not geoms:
        logger.warning("No Open Buildings footprints in the AOI -- cadastre unchanged")
        return

    gob = gpd.GeoDataFrame({"confidence": confs}, geometry=geoms, crs="EPSG:4326")
    gob = gob.to_crs(f"EPSG:{TARGET_EPSG}")
    gob["geometry"] = gob.geometry.buffer(0)
    gob = gob[gob.geometry.notna() & ~gob.geometry.is_empty]
    gob = gpd.clip(gob, shapely_box(*AOI_TARGET_BOUNDS))
    gob = gob[~gob.geometry.is_empty & gob.geometry.notna()].reset_index(drop=True)

    # --- drop footprints standing in the Nile ---
    landuse = gpd.read_file(cadastre_path, layer="landuse")
    water_rows = landuse[landuse["osm_key"].astype(str).str.startswith("WTR")]
    if len(water_rows):
        water_union = unary_union(water_rows.geometry.values)
        on_water = gob.geometry.centroid.within(water_union)
        if on_water.any():
            logger.info("  %d footprints dropped for standing on water", int(on_water.sum()))
            gob = gob[~on_water].reset_index(drop=True)
    else:
        logger.warning("  no water polygons in the cadastre -- run `water` first; "
                       "footprints on the Nile will NOT be removed")

    # --- OSM wins on overlap ---
    osm = gpd.read_file(cadastre_path, layer="buildings")
    if len(osm):
        hit = gpd.sjoin(gob[["geometry"]], osm[["geometry"]],
                        how="inner", predicate="intersects")
        overlapping = set(hit.index)
        if overlapping:
            logger.info("  %d footprints dropped as already mapped in OSM", len(overlapping))
            gob = gob.drop(index=list(overlapping)).reset_index(drop=True)

    # Open Buildings carries no function and no stable id. `building=yes` is
    # honest -- it says "a building, type unrecorded" -- and mapping_eg.yml
    # sends it to the residential default, which is what the overwhelming
    # majority of Aswan's unmapped fabric is. Ids are synthetic and negative
    # so they can never collide with a real OSM id.
    gob["osm_id"] = -(np.arange(len(gob), dtype="int64") + 1)
    gob["osm_type"] = "open_buildings_v3"
    gob["building"] = "yes"

    keep = [c for c in ("osm_id", "osm_type", "building", "name", "osm_height",
                        "osm_levels") if c in osm.columns]
    fused = pd.concat(
        [osm[keep + ["geometry"]], gob.reindex(columns=keep + ["geometry"])
         .assign(osm_id=gob["osm_id"], osm_type=gob["osm_type"], building=gob["building"])],
        ignore_index=True,
    )
    fused = gpd.GeoDataFrame(fused, geometry="geometry",
                             crs=osm.crs or f"EPSG:{TARGET_EPSG}")
    logger.info("Cadastre buildings: %d OSM + %d Open Buildings = %d "
                "(%.3f km2 of footprint in the AOI)",
                len(osm), len(gob), len(fused), fused.geometry.area.sum() / 1e6)

    cadastre_path.unlink()
    landuse.to_file(cadastre_path, layer="landuse", driver="GPKG")
    fused.to_file(cadastre_path, layer="buildings", driver="GPKG")

    cadastre_zip = _zip_gpkg(cadastre_path)
    from shapely.geometry import mapping as _mapping
    from pyproj import Transformer
    xmin, ymin, xmax, ymax = AOI_TARGET_BOUNDS
    (INDEX_DIR / "cadastre_eg.json").write_text(json.dumps({
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "geometry": _mapping(shapely_box(xmin, ymin, xmax, ymax)),
            "properties": {"tile_id": "eg_cadastre", "staged_path": str(cadastre_zip)},
        }],
    }, indent=2))
    logger.info("Cadastre index refreshed: %s", INDEX_DIR / "cadastre_eg.json")


# ---------------------------------------------------------------------------
# 5. Build all indices
# ---------------------------------------------------------------------------
def stage_index() -> None:
    """Build GeoJSON tile indices for all staged sources."""
    _ensure_dirs()
    stage_dem_index()
    stage_dsm_index()
    stage_ortho_index()
    stage_cadastre_index()
    stage_topo_index()


def stage_dem_index() -> None:
    """Build the DEM tile index."""
    from shapely.geometry import mapping, box as shapely_box
    from pyproj import Transformer

    index_path = INDEX_DIR / "dtm_eg.json"
    if index_path.exists():
        logger.info("DEM index already exists")
        return

    xmin, ymin, xmax, ymax = AOI_TARGET_BOUNDS
    tile_name = _degree_tile_name()
    staged = STAGING / f"fabdem_{tile_name}_{TARGET_EPSG}.tif"
    index = {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "geometry": mapping(shapely_box(xmin, ymin, xmax, ymax)),
            "properties": {
                "tile_id": tile_name,
                "staged_path": str(staged),
            },
        }],
    }
    INDEX_DIR.mkdir(parents=True, exist_ok=True)
    index_path.write_text(json.dumps(index, indent=2))
    logger.info("DEM index: %s", index_path)


def stage_dsm() -> None:
    """Download + reproject Copernicus DEM GLO-30 (DSM) to the target CRS.

    palmgeo's ``read_on_grid`` does windowed reads with NO reprojection, so a
    geographic (EPSG:4326) source can never intersect a projected tile and
    yields all-nodata. The DSM must be reprojected to the target CRS, exactly
    like the FABDEM DTM, or every per-tile nDSM comes out empty and all
    buildings are dropped as "too low".
    """
    import rasterio
    from rasterio.warp import calculate_default_transform, reproject, Resampling

    _ensure_dirs()
    staged_path = STAGING / f"copernicus_dsm_{TARGET_EPSG}.tif"
    if staged_path.exists():
        logger.info("Copernicus DSM already staged: %s", staged_path)
        return

    south, west = _degree_tile()
    name = (f"Copernicus_DSM_COG_10_{_hemi(south, 'N', 'S', 2)}_00_"
            f"{_hemi(west, 'E', 'W', 3)}_00_DEM")
    s3_url = f"https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif"
    logger.info("Streaming + reprojecting Copernicus DSM tile %s ...", _degree_tile_name())

    src_path = f"/vsicurl/{s3_url}"
    with rasterio.open(src_path) as src:
        dst_transform, dst_width, dst_height = calculate_default_transform(
            src.crs, f"EPSG:{TARGET_EPSG}",
            src.width, src.height, *src.bounds,
            dst_width=None, dst_height=None,
            resolution=30.0,
        )
        dst_crs = rasterio.crs.CRS.from_epsg(TARGET_EPSG)
        profile = src.profile.copy()
        profile.update(
            crs=dst_crs,
            transform=dst_transform,
            width=dst_width,
            height=dst_height,
            compress="deflate",
            predictor=2,
            nodata=-9999.0,
        )
        with rasterio.open(staged_path, "w", **profile) as dst:
            reproject(
                source=rasterio.band(src, 1),
                destination=rasterio.band(dst, 1),
                src_transform=src.transform,
                src_crs=src.crs,
                src_nodata=src.nodata,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                dst_nodata=-9999.0,
                resampling=Resampling.bilinear,
            )

    logger.info("Copernicus DSM staged: %s", staged_path)


def stage_dsm_index() -> None:
    """Build the Copernicus DEM (DSM) index.

    plan_tiles() calls load_index() for the dsm role, so even though the DSM
    is a single source, it is exposed through a one-feature GeoJSON index
    whose ``staged_path`` is the local reprojected tif produced by
    :func:`stage_dsm`.
    """
    from shapely.geometry import mapping, box as shapely_box
    from pyproj import Transformer

    index_path = INDEX_DIR / "dsm_eg.json"
    if index_path.exists():
        logger.info("DSM index already exists")
        return

    stage_dsm()

    xmin, ymin, xmax, ymax = AOI_TARGET_BOUNDS

    index = {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "geometry": mapping(shapely_box(xmin, ymin, xmax, ymax)),
            "properties": {
                "tile_id": _degree_tile_name(),
                "staged_path": str(STAGING / f"copernicus_dsm_{TARGET_EPSG}.tif"),
            },
        }],
    }
    INDEX_DIR.mkdir(parents=True, exist_ok=True)
    index_path.write_text(json.dumps(index, indent=2))
    logger.info("DSM index: %s", index_path)


def stage_ortho_index() -> None:
    """Build the Sentinel-2 orthophoto index (calls STAC if needed)."""
    stage_ortho()


def stage_cadastre_index() -> None:
    """Build the cadastre index (calls OSM extraction if needed)."""
    stage_osm()


def stage_topo_index() -> None:
    """Build the topology index (already built by stage_osm)."""
    pass


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="stage_eg",
        description="Stage Egypt open-geodata for palmgeo (case-study port).",
    )
    parser.add_argument(
        "command",
        choices=["index", "dem", "dsm", "ortho", "osm", "water", "buildings", "all"],
        help="Which staging step to run",
    )
    parser.add_argument(
        "--workdir", type=Path, required=True,
        help="the run's project.workdir; staged files go to <workdir>/staging",
    )
    parser.add_argument(
        "--bbox", type=float, nargs=4, required=True,
        metavar=("LON0", "LAT0", "LON1", "LAT1"),
        help="area to stage in WGS84 degrees: all domains plus ~500 m margin",
    )
    parser.add_argument(
        "--epsg", type=int, default=DEFAULT_EPSG,
        help=f"target CRS; must match sources_eg.yml (default {DEFAULT_EPSG})",
    )
    parser.add_argument(
        "--open-buildings-tile", default="", metavar="TOKEN",
        help="S2 level-4 cell token of Google Open Buildings v3 (`buildings` only)",
    )
    args = parser.parse_args(argv)
    _configure(args.workdir, args.bbox, args.epsg)
    global GOB_TILE_ID
    GOB_TILE_ID = str(args.open_buildings_tile)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)-18s %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stdout,
    )

    steps = {
        "index": stage_index,
        "dem": stage_dem,
        "dsm": stage_dsm,
        "ortho": stage_ortho,
        "osm": stage_osm,
        # `water` rewrites the cadastre landuse layer and `buildings` rewrites
        # its building layer against that water mask, so this order is load
        # bearing -- dicts preserve insertion order and `all` relies on it.
        "water": stage_water,
        "buildings": stage_buildings,
    }

    if args.command == "all":
        for name, fn in steps.items():
            logger.info("=== %s ===", name)
            try:
                fn()
            except Exception as exc:
                logger.error("%s failed: %s", name, exc)
                return 1
    else:
        try:
            steps[args.command]()
        except Exception as exc:
            logger.error("%s failed: %s", args.command, exc)
            return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
