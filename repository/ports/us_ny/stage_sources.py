#!/usr/bin/env python3
"""Stage NYC/NYS open-geodata sources for palmgeo's US-NY port (case study).

WHY THIS EXISTS: sources_us_ny.yml cannot point palmgeo directly at the real
NYC/NYS endpoints. Two hard constraints force a staging step in between --
see that file's own header for the full reasoning:

  1. CRS/units. Native NYC/NYS rasters are State Plane feet (EPSG:2263 for
     the 2017 LiDAR DEM/DSM per their metadata, EPSG:6539 VERIFIED by opening
     a real ortho tile -- see sources_us_ny.yml's ortho notes). palm_csd and
     palmgeo's own grid arithmetic are metric throughout, and palmgeo refuses
     to reproject raster tiles at load time.
  2. Layer fusion. palmgeo's cadastre role needs ONE zipped GeoPackage with
     `landuse`/`buildings` sublayers (steps/vectors.py::fetch_cadastre is
     hard-coded to that shape). NYC publishes land use (MapPLUTO) and
     building footprints as two separate datasets.

This script downloads from the REAL, verified endpoints below, reprojects
rasters to the run's target CRS (EPSG:6347, the CRS sources_us_ny.yml
declares), and fuses the two cadastre datasets -- then writes the
*_us_ny.geojson tile indices that sources_us_ny.yml's `index_url` entries
point at.

It was written for a run over Manhattan and both river banks. The area and the
run's work directory come from the command line; nothing here names a site:

    STAGE="python repository/ports/us_ny/stage_sources.py \
        --workdir static_driver/<run>/work --bbox <lon0> <lat0> <lon1> <lat1>"
    $STAGE index          # all four staged indices (fast, no raster I/O)
    $STAGE ortho --all    # NYSDOP tiles, ~10 MB each
    $STAGE dtm            # citywide DEM, a 15.6 GB download
    $STAGE dsm            # citywide DSM, a 30.1 GB download
    $STAGE cadastre       # MapPLUTO + Building Footprints + water, fused

``--workdir`` must be the run's ``project.workdir``; sources_us_ny.yml finds
the staged files through ``{workdir}/staging``. ``--bbox`` (WGS84) should
cover every domain plus ~500 m for palmgeo's ``grid.buffer`` and the 1 km
work-tile rounding. Every raster is warped onto its envelope in the target
CRS, so stage_nj_supplement.py has to be given the same box.

STATUS: every step has been run for real and built into finished drivers
(PORTING.md, section 5b, is the write-up). Two things stay manual:
  - the ortho mosaic sources_us_ny.yml points at is assembled by hand from
    the staged tiles: gdalbuildvrt, gdalwarp -nosrcalpha at 0.5 m, then
    tools/fix_colorinterp.py and tools/check_ortho.py (PORTING.md 5b gives
    the traps each of those steps guards against);
  - the coastline water polygons (904 MB, osmdata.openstreetmap.de) must
    already sit under repository/.cache/global/water-polygons-split-4326/;
    see sources_global.yml, section "water".
  - where the box crosses the Hudson, the New Jersey OSM layers
    (<staging>/osm/nj_landuse.gpkg, nj_buildings.gpkg) are cut from the
    Geofabrik New Jersey extract; the rasters for that side come from
    stage_nj_supplement.py.
The numbers quoted in the comments below (sizes, elevations, feature counts)
were measured on that run and describe its area, not yours.

Nothing here is needed for `palmgeo plan` to succeed (it only reads index
metadata, never raster pixels -- see cli.py::stage_plan / steps/tiles.py::
plan_tiles). It IS needed before `palmgeo vectors`/`rasters`/`all` will
produce a correct driver.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # ports/
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # palmgeo
from staging_index import aoi_bounds  # noqa: E402
from palmgeo import USER_AGENT  # noqa: E402

logger = logging.getLogger("stage_us_ny")

PROJECT_ROOT = Path(__file__).resolve().parents[3]
#: The CRS sources_us_ny.yml declares (NAD83(2011) / UTM 18N, metres).
TARGET_EPSG = 6347

# Run settings, set from the command line by _configure().
WORKDIR = Path()
STAGING = Path()
#: The staged area, WGS84 (lon0, lat0, lon1, lat1).
BBOX_WGS84: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)
#: The same area in the target CRS -- gdalwarp's -te and every index geometry.
AOI_TARGET_BOUNDS: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)

ORTHO_QUERY_URL = (
    "https://orthos-np.its.ny.gov/arcgis/rest/services/vector/ortho_indexes/"
    "MapServer/5/query"
)
# The dataset id f888-ni5f used in an earlier draft of this script
# (data.cityofnewyork.us/api/geospatial/f888-ni5f?method=export...) turned
# out to be an "href"-type Socrata view -- a redirect page, not real hosted
# data (confirmed 2026-09-04: both the NYC and the NY State mirror return
# "Unexportable view type: href" / only 2 PDF attachments, no shapefile).
# MapPLUTO's real, live, queryable source is NYC DCP's own ArcGIS
# FeatureServer, verified 2026-09-04: real schema (LandUse, BldgClass,
# Borough, ... confirmed present), real query tested for Borough='MN',
# maxRecordCount=2000 (pagination required, implemented in
# fetch_mappluto() below).
MAPPLUTO_QUERY_URL = (
    "https://services5.arcgis.com/GfwWNkhOj9bNBqoJ/arcgis/rest/services/"
    "MAPPLUTO/FeatureServer/0/query"
)
# Real, VERIFIED sizes/range-support as of 2026-09-03 (see module docstring).
DEM_2017_URL = (
    "https://data.cityofnewyork.us/api/views/7sc8-jtbz/files/"
    "a43310be-e215-4f65-8fec-707be53bfc1c"
    "?filename=Bare_Earth_Topographic_Hydroflattened.zip"
)
DSM_2017_URL = (
    "https://data.cityofnewyork.us/api/views/7sc8-jtbz/files/"
    "70ab3695-bf4e-4150-b8b5-c3590042ba57?filename=Highest_Hit.zip"
)
# VERIFIED WORKING 2026-09-04 (unlike f888-ni5f above): a real Socrata
# geospatial export, 288 MB citywide, real fields (base_bbl, construction_year,
# ground_elevation, mappluto_bbl, feature_code, ...). Downloaded once and
# cached as <staging>/building_footprints_citywide_native.geojson rather than
# re-fetched by stage_cadastre() every run.
BUILDING_FOOTPRINTS_URL = (
    "https://data.cityofnewyork.us/api/geospatial/5zhs-2jue?method=export&format=GeoJSON"
)

# WATER GAP FIX (2026-09-04) -- see mapping_us_ny.yml's landuse."WTR" entry.
# MapPLUTO has no polygon at all for the Hudson/East/Harlem rivers (they are
# not tax lots), so there is nothing to override via a mapping rule -- the
# fix has to ADD polygons, not remap existing ones. Two real, GLOBAL sources
# combine to cover both cases (see repository/config/sources_global.yml,
# section "water", for the general recipe this implements):
#   - inland lakes/ponds/reservoirs: OSM `natural=water`, queried live via
#     Overpass and cached as <staging>/osm_water_natural_water.json (354
#     features over the Manhattan run, verified 2026-09-04 -- Central Park's
#     ponds/reservoir among them).
#   - the big tidal rivers themselves: NOT tagged `natural=water` in OSM
#     (verified -- a targeted name search for "Hudson River"/"East River"/
#     "Harlem River" found no matching NYC-scale water polygon). They rely on
#     the `natural=coastline` convention instead, which needs the dedicated
#     osmdata.openstreetmap.de "water polygons" product (coastline ways
#     assembled into closed, split polygons) to become an area at all.
#     Downloaded once (904 MB, real, global, ODbL) to
#     repository/.cache/global/water-polygons-split-4326/ -- a SHARED asset,
#     not staged per-run, since it is region-agnostic by construction (see
#     sources_global.yml). Filtered to the Manhattan run's box: 5 features,
#     real NY Harbor/Hudson/East River polygons, verified 2026-09-04.
GLOBAL_WATER_POLYGONS = (
    PROJECT_ROOT / "repository" / ".cache" / "global" / "water-polygons-split-4326"
    / "water_polygons.shp"
)
WATER_LANDUSE_CODE = "WTR"  # synthetic MapPLUTO-shaped code, see mapping_us_ny.yml

# Verified against the NYC DCP PLUTO Data Dictionary (May 2026, 26v1).
MAPPLUTO_LANDUSE_LABELS = {
    "01": "One & Two Family Buildings",
    "02": "Multi-Family Walk-Up Buildings",
    "03": "Multi-Family Elevator Buildings",
    "04": "Mixed Residential & Commercial Buildings",
    "05": "Commercial & Office Buildings",
    "06": "Industrial & Manufacturing",
    "07": "Transportation & Utility",
    "08": "Public Facilities & Institutions",
    "09": "Open Space & Outdoor Recreation",
    "10": "Parking Facilities",
    "11": "Vacant Land",
}


def _bbox_polygon(bbox: tuple[float, float, float, float]) -> dict[str, Any]:
    x0, y0, x1, y1 = bbox
    return {
        "type": "Polygon",
        "coordinates": [[[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]],
    }


# Caches under the staging directory, set together with it.
BUILDING_FOOTPRINTS_CACHE = Path()
MAPPLUTO_CACHE = Path()
OSM_NATURAL_WATER_CACHE = Path()


def _configure(workdir: Path, bbox: list[float]) -> None:
    """Set the run settings from the command line."""
    global WORKDIR, STAGING, BBOX_WGS84, AOI_TARGET_BOUNDS
    global BUILDING_FOOTPRINTS_CACHE, MAPPLUTO_CACHE, OSM_NATURAL_WATER_CACHE
    lon0, lat0, lon1, lat1 = (float(v) for v in bbox)
    if not (-180 <= lon0 < lon1 <= 180 and -90 <= lat0 < lat1 <= 90):
        raise SystemExit(f"--bbox {bbox}: expected lon0 lat0 lon1 lat1 in WGS84 degrees")
    WORKDIR = (workdir if workdir.is_absolute() else Path.cwd() / workdir).resolve()
    STAGING = WORKDIR / "staging"
    BBOX_WGS84 = (lon0, lat0, lon1, lat1)
    AOI_TARGET_BOUNDS = aoi_bounds(BBOX_WGS84, TARGET_EPSG)
    BUILDING_FOOTPRINTS_CACHE = STAGING / "building_footprints_citywide_native.geojson"
    MAPPLUTO_CACHE = STAGING / "mappluto_aoi_native.gpkg"
    OSM_NATURAL_WATER_CACHE = STAGING / "osm_water_natural_water.json"


# ---------------------------------------------------------------------------
# index builders -- fast, no raster I/O, safe to re-run any time
# ---------------------------------------------------------------------------
def build_ortho_index() -> None:
    """Write a ONE-FEATURE index pointing at the staged ortho MOSAIC.

    Why a mosaic and not the 413 individual tiles: NYSDOP tiles are 2500 ft
    (762 m), palmgeo's work tiles are 1 km (grid.py TILE_SIZE), and
    steps/tiles.py picks exactly ONE ortho item per work tile. A 762 m tile
    cannot cover a 1 km work tile, so each work tile came out ~80 % black and
    its NDVI was meaningless (measured: Central Park work tile had NDVI p90
    +0.055 and 2.7 % of pixels above the vegetation threshold, against +0.384
    read straight from the same imagery). This is exactly pitfall 4 in
    PORTING.md -- "Do the tile grids of different products align?".
    PORTING.md suggests adjusting TILE_SIZE; mosaicking instead keeps
    palmgeo untouched and matches how dtm/dsm are already staged here.

    The mosaic is assembled by hand from the staged tiles (see the module
    docstring) at 0.5 m -- four times finer than a 2 m work tier, so nothing
    downstream is starved of detail, while a native-0.152 m mosaic of the
    Manhattan run would have been ~38 GB.
    """
    (STAGING / "indices").mkdir(parents=True, exist_ok=True)
    mosaic = STAGING / "ortho" / f"nyc_ortho_epsg{TARGET_EPSG}.tif"
    bbox_target = AOI_TARGET_BOUNDS   # what is really staged, see note
    feature = {
        "type": "Feature",
        "geometry": _bbox_polygon(bbox_target),
        "properties": {
            "tile_id": "nyc_ortho_2024_nysdop_naip2023",
            # The mosaic blends TWO campaigns: NYSDOP 2024 (New York, on top)
            # and NAIP 2023-08-20 (New Jersey, underneath, where NYSDOP stops
            # at the state line). One index feature can carry only one date,
            # so this is the dominant campaign's year -- not a measured date,
            # and not uniform across the mosaic. It feeds
            # `ortho_epochs`/`max_epoch_spread_days`, which therefore cannot
            # see the real ~10-month spread between the two banks. Acceptable
            # where phenology is disabled, as it was for the Manhattan run; wire a per-tile
            # index (build_ortho_tile_index) if that ever changes.
            "acquired": "2024-06-01",
            "staged_path": str(mosaic),
        },
    }
    out = {"type": "FeatureCollection", "features": [feature]}
    path = STAGING / "indices" / "ortho_us_ny.geojson"
    path.write_text(json.dumps(out))
    logger.info("ortho: wrote 1-feature mosaic index to %s", path)


def build_ortho_tile_index() -> None:
    """The older per-tile index -- kept because stage_ortho_tiles() needs the
    per-tile download URLs. Not used as sources_us_ny.yml's `index_url` any
    more; see build_ortho_index() for why."""
    (STAGING / "indices").mkdir(parents=True, exist_ok=True)
    params = {
        "geometry": ",".join(str(v) for v in BBOX_WGS84),
        "geometryType": "esriGeometryEnvelope",
        "inSR": 4326,
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "*",
        "f": "geojson",
    }
    resp = requests.get(ORTHO_QUERY_URL, params=params, timeout=60)
    resp.raise_for_status()
    native = resp.json()
    native_path = STAGING / "nysdop_2024_aoi_native.geojson"
    native_path.write_text(json.dumps(native))
    logger.info("ortho: %d native features queried live, cached at %s",
                len(native["features"]), native_path)

    from pyproj import Transformer

    tr = Transformer.from_crs(4326, TARGET_EPSG, always_xy=True)
    out_features = []
    for feat in native["features"]:
        props = feat["properties"]
        tile_id = Path(str(props["FILENAME"])).stem
        coords = feat["geometry"]["coordinates"][0]
        reproj_ring = [list(tr.transform(x, y)) for x, y in coords]
        out_features.append({
            "type": "Feature",
            "geometry": {"type": "Polygon", "coordinates": [reproj_ring]},
            "properties": {
                "tile_id": tile_id,
                "acquired": "2024-06-01",
                "staged_path": str(STAGING / "ortho" / f"{tile_id}_epsg{TARGET_EPSG}.tif"),
                "native_direct_dl": props["DIRECT_DL"],
            },
        })
    index_path = STAGING / "indices" / "ortho_tiles_us_ny.geojson"
    index_path.write_text(json.dumps({"type": "FeatureCollection", "features": out_features}))
    logger.info("ortho: wrote %d-feature tile index to %s", len(out_features), index_path)


def build_dtm_dsm_index() -> None:
    """Write single-feature indices for the two citywide DEM/DSM archives.

    Unlike NYSDOP, these are NOT pre-tiled -- one zip covers all 5 boroughs
    (real, verified sizes: 15.57 GB / 30.07 GB, see module docstring). A
    palmgeo tile_index with exactly one feature, whose geometry covers this
    run's whole AOI, is a legitimate degenerate case (sources.py::load_index
    has no minimum feature count) and lets `newest_per_tile`/`pick_paired`
    (steps/tiles.py) work unmodified: every work tile's bounds fall inside
    this one feature's bounds, so it is always the (only, and therefore
    newest) match.
    """
    bbox_target = AOI_TARGET_BOUNDS   # what is really staged, see note
    geometry = _bbox_polygon(bbox_target)
    # staged_name must match what stage_dtm_dsm() actually writes: the
    # metre-converted raster (`_m` suffix), NOT the intermediate still in US
    # survey feet. An earlier version emitted "..._dem_epsg6347.tif" here
    # while the pipeline wrote "..._dtm_epsg6347_m.tif" -- the index pointed
    # at a file that never existed.
    for role, tile_id, url, staged_name in (
        ("dtm", "nyc_dem_2017", DEM_2017_URL, f"nyc_dtm_epsg{TARGET_EPSG}_mw"),
        ("dsm", "nyc_dsm_2017", DSM_2017_URL, f"nyc_dsm_epsg{TARGET_EPSG}_mw"),
    ):
        feature = {
            "type": "Feature",
            "geometry": geometry,
            "properties": {
                "tile_id": tile_id,
                # 2017-06-01: midpoint of the verified 2017-05-03..2017-07-26
                # topographic capture window, not a per-tile measured date.
                "acquired": "2017-06-01",
                "staged_path": str(STAGING / role / f"{staged_name}.tif"),
                "native_url": url,
            },
        }
        out = {"type": "FeatureCollection", "features": [feature]}
        path = STAGING / "indices" / f"{role}_us_ny.geojson"
        path.write_text(json.dumps(out))
        logger.info("%s: wrote 1-feature staged index to %s", role, path)


def build_cadastre_index() -> None:
    """Single-feature index pointing at the fused, zipped GeoPackage that
    stage_cadastre() below produces (landuse from MapPLUTO, buildings from
    NYC Building Footprints -- see sources_us_ny.yml's nyc_cadastre notes for
    why these must be fused into one file before palmgeo can read them)."""
    bbox_target = AOI_TARGET_BOUNDS   # what is really staged, see note
    feature = {
        "type": "Feature",
        "geometry": _bbox_polygon(bbox_target),
        "properties": {
            "tile_id": "nyc_cadastre",
            "acquired": "2026-09-03",
            "staged_path": str(STAGING / "cadastre" / "nyc_cadastre.zip"),
        },
    }
    out = {"type": "FeatureCollection", "features": [feature]}
    path = STAGING / "indices" / "cadastre_us_ny.geojson"
    path.write_text(json.dumps(out))
    logger.info("cadastre: wrote 1-feature staged index to %s", path)


# ---------------------------------------------------------------------------
# actual raster/vector staging -- real downloads, only run on request
# ---------------------------------------------------------------------------
def stage_ortho_tiles(tile_ids: list[str] | None, run_all: bool = False) -> None:
    """Download + extract + reproject real NYSDOP tiles to local COGs.

    Proven end-to-end on tile 000260 by hand before this script existed
    (real download, JP2OpenJPEG read, EPSG:6539 -> EPSG:6347 reproject,
    5000x5000 -> 5057x5057 px, bounds sanity-checked against the run's
    domain origin). This function formalises that, nothing more.
    """
    import rasterio
    from rasterio.warp import Resampling, calculate_default_transform, reproject

    index_path = STAGING / "indices" / "ortho_tiles_us_ny.geojson"
    features = json.loads(index_path.read_text())["features"]
    if not run_all:
        wanted = set(tile_ids or [])
        features = [f for f in features if f["properties"]["tile_id"] in wanted]
    (STAGING / "ortho").mkdir(parents=True, exist_ok=True)

    for feat in features:
        props = feat["properties"]
        tile_id = props["tile_id"]
        staged_path = Path(props["staged_path"])
        if staged_path.exists():
            logger.info("ortho %s: already staged, skipping", tile_id)
            continue
        zip_path = STAGING / "ortho" / f"{tile_id}_native.zip"
        if not zip_path.exists():
            # gisdata.ny.gov resets connections under concurrent load -- with
            # 4 parallel workers two of them died on
            # ChunkedEncodingError/ConnectionResetError(104) after ~34 tiles
            # each. Retry with backoff instead of losing the worker.
            last_exc: Exception | None = None
            for attempt in range(1, 6):
                try:
                    logger.info("ortho %s: downloading %s (attempt %d)",
                                tile_id, props["native_direct_dl"], attempt)
                    r = requests.get(props["native_direct_dl"], timeout=180)
                    r.raise_for_status()
                    zip_path.write_bytes(r.content)
                    last_exc = None
                    break
                except Exception as exc:  # noqa: BLE001 - any transport error
                    last_exc = exc
                    zip_path.unlink(missing_ok=True)
                    wait = min(60, 5 * 2 ** (attempt - 1))
                    logger.warning("ortho %s: %s -- retrying in %d s",
                                   tile_id, type(exc).__name__, wait)
                    time.sleep(wait)
            if last_exc is not None:
                # skip this tile rather than kill the whole worker; a rerun
                # picks it up because staging is idempotent
                logger.error("ortho %s: giving up after 5 attempts (%s)",
                             tile_id, last_exc)
                continue
        with zipfile.ZipFile(zip_path) as zf:
            jp2_name = next(n for n in zf.namelist() if n.lower().endswith(".jp2")
                             and not n.lower().endswith(".aux.xml"))
            zf.extract(jp2_name, STAGING / "ortho")
        jp2_path = STAGING / "ortho" / jp2_name

        with rasterio.open(jp2_path) as src:
            transform, width, height = calculate_default_transform(
                src.crs, f"EPSG:{TARGET_EPSG}", src.width, src.height, *src.bounds)
            kwargs = src.meta.copy()
            kwargs.update(crs=f"EPSG:{TARGET_EPSG}", transform=transform, width=width,
                           height=height, driver="GTiff", compress="DEFLATE", tiled=True)
            with rasterio.open(staged_path, "w", **kwargs) as dst:
                for i in range(1, src.count + 1):
                    reproject(source=rasterio.band(src, i), destination=rasterio.band(dst, i),
                              src_transform=src.transform, src_crs=src.crs,
                              dst_transform=transform, dst_crs=f"EPSG:{TARGET_EPSG}",
                              resampling=Resampling.bilinear)
        jp2_path.unlink()
        zip_path.unlink()
        logger.info("ortho %s: staged to %s", tile_id, staged_path)


def stage_dtm_dsm(role: str) -> None:
    """Download the citywide DEM/DSM zip, clip to the AOI, reproject.

    GDAL's /vsicurl/ could not open these specific Socrata-hosted URLs
    directly in the environment this was written in (a HEAD-style probe
    timed out, even though plain `curl` with an explicit Range: header
    worked fine against the same URL -- likely a Socrata-side quirk, not a
    general vsicurl problem: a public S3-hosted COG opened via /vsicurl/
    without issue in the same session). The workaround below avoids
    GDAL-level REMOTE access entirely: download with `requests`/`curl`
    (already proved to work against this host), then let rasterio/GDAL open
    the LOCAL zip via /vsizip/ -- which DOES work locally, VERIFIED
    2026-09-04 (real windowed read of the DTM's Manhattan window, 76161 x
    37074 px, 101 s, real elevations 0..266.7 ft incl. the correct float32
    nodata sentinel -3.4e38). Real internal structure, VERIFIED 2026-09-04
    via `unzip -l` on the completed download: ONE citywide GeoTIFF per role
    (Bare_Earth_Topographic_Hydroflattened/bh_NYC.tif, 154859x153378 px;
    Highest_Hit/ is expected to follow the same
    "<Folder>/<prefix>_NYC.tif" pattern but was not confirmed independently
    -- check with `unzip -l` before assuming). Native CRS VERIFIED via
    gdalinfo: EPSG:2263 (NAD83 / NY Long Island ftUS -- the ORIGINAL
    realization, unlike nysdop's EPSG:6539 -- see sources_us_ny.yml's ortho
    notes for why the two differ and why it does not matter for the target
    CRS either way).
    """
    url = {"dtm": DEM_2017_URL, "dsm": DSM_2017_URL}[role]
    # VERIFIED for real via `unzip -l` on both completed downloads
    # 2026-09-04 (the dsm name is "hh_NYC.tif", NOT "bh_NYC.tif" -- an
    # earlier version of this comment guessed wrong before checking).
    internal_tif = {
        "dtm": "Bare_Earth_Topographic_Hydroflattened/bh_NYC.tif",
        "dsm": "Highest_Hit/hh_NYC.tif",
    }[role]
    zip_path = STAGING / role / f"nyc_2017_{role}_native.zip"
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    if not zip_path.exists():
        logger.warning("%s: about to download %s (~%s GB) -- this will take a while",
                        role, url, "15.6" if role == "dtm" else "30.1")
        with requests.get(url, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(zip_path, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)

    # EXTRACT to a plain local file rather than reading through /vsizip/.
    # Windowed+decimated reads via /vsizip/ on this exact ~16/30 GB citywide
    # single-file zip took down the whole interactive host process
    # three times in a row (2026-09-04) even running alone with no other
    # load -- a much simpler, well-trodden GDAL code path (a plain local
    # GeoTIFF, no zip layer at all) is used instead, at the cost of ~16/30 GB
    # of extra disk (194 GB free at the time this was written, so not a
    # constraint) and a slower, but far more predictable, unzip step first.
    extracted_path = STAGING / role / Path(internal_tif).name
    if not extracted_path.exists():
        logger.info("%s: extracting %s from %s (this takes a few minutes)",
                    role, internal_tif, zip_path)
        with zipfile.ZipFile(zip_path) as zf:
            with zf.open(internal_tif) as src_f, open(extracted_path, "wb") as dst_f:
                while chunk := src_f.read(1 << 24):
                    dst_f.write(chunk)
        logger.info("%s: extracted to %s", role, extracted_path)

    # CLIP + REPROJECT + UNIT-CONVERT, via gdalwarp rather than a numpy
    # windowed read.
    #
    # HISTORY, because this exact step is a trap: the first implementation
    # did a rasterio windowed read with `out_shape` decimation and computed
    # the decimation factor as `round(1.0 / src.transform.a)` -- intending
    # "target 1 metre". But `src.transform.a` is in the SOURCE CRS's units,
    # and EPSG:2263 is in US survey FEET, so transform.a == 1.0 and the
    # factor came out as 1: no decimation at all. The read then materialised
    # the full-resolution window (76161 x 37074 px = 2.82e9 px, ~11 GB as
    # float32, ~25 GB once the averaging and warp buffers were added) on a
    # 30 GB machine with 511 MB of swap. The kernel OOM-killer took the
    # process out -- and because it ran inside the editor's process scope
    # (app-code-*.scope), systemd tore the whole scope down with it. That
    # crashed the editor four times before the kernel log was consulted:
    #   Out of memory: Killed process ... (python3) anon-rss:25020096kB
    #   app-code-3783.scope: Failed with result 'oom-kill'.
    # The same metre/foot confusion this port documents so carefully for the
    # CRS choice, walked into in the resampling arithmetic.
    #
    # gdalwarp streams block-wise with an explicit memory budget, uses the
    # source's own overviews (this GeoTIFF has them: 77430x76689 down to
    # 1210x1199) and never materialises the raster. Run it under a guard
    # that sets `ulimit -v` and oom_score_adj=1000, so a future mistake
    # fails inside this process instead of killing the editor.
    staged_ft = STAGING / role / f"nyc_{role}_epsg{TARGET_EPSG}_ft.tif"
    staged_path = STAGING / role / f"nyc_{role}_epsg{TARGET_EPSG}_m.tif"
    if not staged_ft.exists():
        cmd = [
            "gdalwarp",
            "-t_srs", f"EPSG:{TARGET_EPSG}",
            "-te", *[str(v) for v in AOI_TARGET_BOUNDS],
            "-tr", "1", "1",
            "-r", "average",
            "-of", "GTiff",
            "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES", "-co", "BIGTIFF=IF_SAFER",
            "-wm", "512", "--config", "GDAL_CACHEMAX", "512",
            "-overwrite", str(extracted_path), str(staged_ft),
        ]
        logger.info("%s: %s", role, " ".join(cmd))
        subprocess.run(cmd, check=True)

    # gdalwarp reprojects COORDINATES but never touches pixel VALUES, so the
    # warped raster still carries elevation in US survey feet. PALM/palm_csd
    # are metric -- without this step terrain would be 3.28x too high, and
    # silently so. ft_to_m.py works block-wise (bounded memory) and also
    # normalises the -3.4e38 fill to the -9999.0 this port's sources_us_ny.yml
    # declares as `nodata`.
    if not staged_path.exists():
        logger.info("%s: converting ftUS -> m", role)
        subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "repository" / "tools" / "ft_to_m.py"),
             str(staged_ft), str(staged_path)],
            check=True,
        )
    logger.info("%s: staged to %s", role, staged_path)


def _fetch_osm_natural_water() -> None:
    """Query Overpass for `natural=water` over the box, cached as raw JSON."""
    lon0, lat0, lon1, lat1 = BBOX_WGS84
    bbox = f"{lat0},{lon0},{lat1},{lon1}"
    query = (f'[out:json][timeout:180];(way["natural"="water"]({bbox});'
             f'relation["natural"="water"]({bbox}););out geom;')
    logger.info("water: querying Overpass for natural=water over the box")
    # Overpass answers the default python-requests User-Agent with HTTP 406.
    r = requests.post("https://overpass-api.de/api/interpreter", data=query,
                      timeout=300, headers={"User-Agent": USER_AGENT})
    r.raise_for_status()
    OSM_NATURAL_WATER_CACHE.parent.mkdir(parents=True, exist_ok=True)
    OSM_NATURAL_WATER_CACHE.write_text(r.text)


def build_water_layer() -> "Any":
    """Assemble the water-fix polygons (see the module-level comment above
    GLOBAL_WATER_POLYGONS) into a GeoDataFrame shaped like a MapPLUTO
    `landuse` extract: a `LandUse` column set to WATER_LANDUSE_CODE, so
    mapping_us_ny.yml's ordinary landuse-rule mechanism (the ONLY place
    compose.py actually reads a `water:` key -- see steps/compose.py:124)
    picks these polygons up with no other code or schema change.

    RUN 2026-09-04: both inputs are real (Overpass live query cached at
    OSM_NATURAL_WATER_CACHE; GLOBAL_WATER_POLYGONS downloaded and filtered to
    the Manhattan run's box, 5 features, 1028 km2 before clipping -- see the
    module comment). This function's own clip-to-AOI + concat was not yet
    run end-to-end as part of a full stage_cadastre() pass.
    """
    import geopandas as gpd
    from shapely.geometry import Polygon, box

    aoi = box(*AOI_TARGET_BOUNDS)

    # OSM_NATURAL_WATER_CACHE is raw Overpass "out geom" JSON, not GeoJSON --
    # convert the two element shapes it actually contains. `way` elements
    # carry an ordered `geometry` ring directly; `relation` (multipolygon)
    # elements carry `members`, each with its own `geometry` ring -- take the
    # `outer` ones (holes from `inner` rings are ignored, a real but minor
    # simplification: none of the Manhattan `natural=water` relations found
    # 2026-09-04 have documented inner rings, e.g. an island in a reservoir).
    if not OSM_NATURAL_WATER_CACHE.exists():
        _fetch_osm_natural_water()
    raw = json.loads(OSM_NATURAL_WATER_CACHE.read_text())
    polys, codes = [], []

    # OSM's `water=*` subtag carries exactly the distinction PALM's water_type
    # draws. Getting this table wrong is not cosmetic: a first version had no
    # `river` entry at all, so the East River and Newtown Creek fell through
    # to the LAKE default, and `water=wastewater` -- the Wards Island
    # treatment plant tanks on Randall's Island -- became a lake too.
    _WATER_CODES = {
        # tidal channels, creeks, canals: moving water
        "river": "WTR_RIV", "canal": "WTR_RIV", "stream": "WTR_RIV",
        "ditch": "WTR_RIV", "tidal_channel": "WTR_RIV", "fairway": "WTR_RIV",
        # small artificial or standing bodies
        "pond": "WTR_POND", "wastewater": "WTR_POND", "basin": "WTR_POND",
        "reflecting_pool": "WTR_POND", "fountain": "WTR_POND",
        "salt_pond": "WTR_POND", "moat": "WTR_POND",
        # genuine lakes and reservoirs
        "lake": "WTR_LAKE", "reservoir": "WTR_LAKE", "lagoon": "WTR_LAKE",
        "oxbow": "WTR_LAKE",
    }

    def _code_for(tags: dict) -> str:
        sub = str(tags.get("water", "")).lower()
        if sub in _WATER_CODES:
            return _WATER_CODES[sub]
        if str(tags.get("waterway", "")).lower() in {"river", "riverbank", "canal", "stream"}:
            return "WTR_RIV"
        # Untagged `natural=water` is decided geometrically below: touching
        # the tidal system makes it river, otherwise lake. A blanket LAKE
        # default put the East River's untagged parts in the wrong class.
        return "WTR_UNTAGGED"

    from shapely.geometry import LineString
    from shapely.ops import linemerge, polygonize, unary_union

    def _rings(members, role):
        """Stitch a multipolygon relation's member ways into closed rings.

        An OSM multipolygon's outer boundary is routinely SPLIT across
        several open ways; each member on its own is a line, not a ring.
        Treating members as independent polygons (as a first version did)
        turned Central Park's Lake -- a 10-member relation, ~80 000 m2 --
        into two slivers of 31 and 4 284 m2, and it then classified as
        vegetation in the driver. linemerge + polygonize reassembles them.
        """
        lines = [LineString([(q["lon"], q["lat"]) for q in m["geometry"]])
                 for m in members
                 if m.get("role", "") == role and len(m.get("geometry") or []) >= 2]
        if not lines:
            return []
        merged = linemerge(lines)
        return list(polygonize(merged if merged.geom_type == "MultiLineString"
                               else [merged]))

    for el in raw.get("elements", []):
        tags = el.get("tags", {}) or {}
        code = _code_for(tags)
        if el["type"] == "way" and el.get("geometry"):
            ring = [(q["lon"], q["lat"]) for q in el["geometry"]]
            if len(ring) >= 4:
                polys.append(Polygon(ring)); codes.append(code)
        elif el["type"] == "relation":
            members = el.get("members", [])
            outer = _rings(members, "outer")
            if not outer:
                continue
            inner = _rings(members, "inner")
            shape = unary_union(outer)
            if inner:
                # islands inside a water body are not water
                shape = shape.difference(unary_union(inner))
            if not shape.is_empty:
                for g in (shape.geoms if shape.geom_type == "MultiPolygon" else [shape]):
                    polys.append(g); codes.append(code)
    keep = [(g, c) for g, c in zip(polys, codes) if g.is_valid]
    natural = gpd.GeoDataFrame(
        {"LandUse": [c for _, c in keep]},
        geometry=[g for g, _ in keep], crs=4326,
    ).to_crs(epsg=TARGET_EPSG)

    coastline = gpd.read_file(GLOBAL_WATER_POLYGONS, bbox=BBOX_WGS84)
    coastline = coastline.to_crs(epsg=TARGET_EPSG)

    # coastline-derived polygons are the tidal rivers/harbour -> river
    coastline = coastline.copy()
    coastline["LandUse"] = "WTR_RIV"
    water = gpd.GeoDataFrame(
        {"LandUse": list(natural["LandUse"]) + list(coastline["LandUse"])},
        geometry=list(natural.geometry) + list(coastline.geometry),
        crs=f"EPSG:{TARGET_EPSG}",
    )
    water = gpd.clip(water, aoi)

    # Resolve untagged bodies by contact with the tidal system: the
    # coastline-derived polygons ARE the ocean-connected water, so anything
    # touching them is part of it (river), and anything isolated is standing
    # water (lake).
    unt = water["LandUse"] == "WTR_UNTAGGED"
    if unt.any():
        tidal = water.loc[water["LandUse"] == "WTR_RIV", "geometry"]
        tidal_union = unary_union(list(tidal)) if len(tidal) else None
        def _resolve(geom):
            if tidal_union is not None and geom.intersects(tidal_union):
                return "WTR_RIV"
            return "WTR_LAKE"
        water.loc[unt, "LandUse"] = [
            _resolve(g) for g in water.loc[unt, "geometry"]
        ]
        logger.info("water: resolved %d untagged bodies by contact with the "
                    "tidal system", int(unt.sum()))
    labels = {"WTR_LAKE": "Water: lake/reservoir (OSM natural=water)",
              "WTR_POND": "Water: pond (OSM natural=water)",
              "WTR_RIV": "Water: tidal river/harbour (OSM coastline polygons)"}
    water["LandUse_label"] = water["LandUse"].map(labels)
    counts = water["LandUse"].value_counts().to_dict()
    logger.info("water: %d natural=water + %d coastline-derived, %d after clipping; "
                "by type %s", len(natural), len(coastline), len(water), counts)
    return water[["LandUse", "LandUse_label", "geometry"]]


def fetch_mappluto() -> "Any":
    """Page through the real, live MAPPLUTO ArcGIS FeatureServer over the AOI
    and return a GeoDataFrame with (at least) LandUse + geometry.

    The query is spatial, by the staged box -- see the comment on `params`
    below for why not by borough.
    """
    import geopandas as gpd

    # Cache the paged result: the AOI-wide query walks >222 000 tax lots at
    # 2000 per page and takes the better part of an hour. Iterating on the
    # water layer or the fusion should not pay that cost again -- delete
    # MAPPLUTO_CACHE to force a refresh.
    if MAPPLUTO_CACHE.exists():
        logger.info("mappluto: reusing cached fetch at %s", MAPPLUTO_CACHE)
        return gpd.read_file(MAPPLUTO_CACHE)

    all_features: list[dict[str, Any]] = []
    offset = 0
    page_size = 2000
    while True:
        params = {
            # Spatially by AOI, NOT `where=Borough='MN'`. The domain reaches
            # across the East River, so a borough filter silently leaves
            # Brooklyn, Queens and Roosevelt Island without land use OR
            # buildings -- visible in the first driver as an empty east bank.
            "geometry": ",".join(str(v) for v in BBOX_WGS84),
            "geometryType": "esriGeometryEnvelope",
            "inSR": 4326,
            "spatialRel": "esriSpatialRelIntersects",
            "outFields": "LandUse,BldgClass,Borough,Block,Lot",
            "resultOffset": offset,
            "resultRecordCount": page_size,
            "f": "geojson",
        }
        # Retry each page: the AOI-wide query pages through >222 000 tax lots
        # (Brooklyn and Queens dwarf Manhattan in lot count), and a single
        # transient 502 used to abort the whole fetch and discard everything
        # already downloaded.
        data = None
        for attempt in range(1, 7):
            try:
                resp = requests.get(MAPPLUTO_QUERY_URL, params=params, timeout=120)
                resp.raise_for_status()
                data = resp.json()
                break
            except Exception as exc:  # noqa: BLE001 - transport or HTTP error
                wait = min(60, 5 * 2 ** (attempt - 1))
                logger.warning("mappluto offset %d: %s -- retry %d in %d s",
                               offset, type(exc).__name__, attempt, wait)
                time.sleep(wait)
        if data is None:
            raise RuntimeError(
                f"mappluto: giving up at offset {offset} after 6 attempts; "
                f"{len(all_features)} features fetched so far"
            )
        feats = data.get("features", [])
        all_features.extend(feats)
        logger.info("mappluto: fetched %d features (offset %d)", len(feats), offset)
        if len(feats) < page_size:
            break
        offset += page_size
    out = gpd.GeoDataFrame.from_features(all_features, crs=4326)
    MAPPLUTO_CACHE.parent.mkdir(parents=True, exist_ok=True)
    out.to_file(MAPPLUTO_CACHE, driver="GPKG")
    logger.info("mappluto: cached %d features to %s", len(out), MAPPLUTO_CACHE)
    return out



def _nj_layers():
    """New Jersey land use and buildings, from the Geofabrik OSM extract.

    MapPLUTO stops at the state line and NYC Building Footprints stop at the
    city line, so once the domain reaches the Hudson's west bank everything
    there is unclassified and building-free -- which is what the driver
    showed. OSM is the only source covering both banks uniformly.

    Land-use keys are composed as `<osm-key>:<value>` (landuse > natural >
    leisure, first present wins) to match the rules copied into
    mapping_us_ny.yml from mapping_global.yml.

    Building heights are NOT read from OSM tags. MEASURED 2026-09-07 on the
    64 428 footprints of this AOI: 0.1 % carry `height` and 3.0 %
    `building:levels` -- an earlier version of this comment claimed "about a
    third", which was wrong. They are measured from the 3DEP DSM/DTM pair
    instead, exactly as on the NYC side. Where that DSM carries no structures
    (the 2014 NJ_SdL5 survey over the North Bergen industrial flats) the
    buildings are simply lost; an OSM height fallback would recover 3 % of
    them and is not worth wiring.
    """
    import geopandas as gpd

    lu_path = STAGING / "osm" / "nj_landuse.gpkg"
    bl_path = STAGING / "osm" / "nj_buildings.gpkg"
    if not lu_path.exists() or not bl_path.exists():
        logger.warning("no New Jersey OSM extract -- run "
                       "repository/ports/us_ny/stage_nj_supplement.py first; the west "
                       "bank will stay empty")
        return None, None

    lu = gpd.read_file(lu_path, layer="landuse")
    key = None
    for col in ("landuse", "natural", "leisure"):
        if col not in lu.columns:
            continue
        part = lu[col].map(lambda v, c=col: f"{c}:{v}" if v else None)
        key = part if key is None else key.fillna(part)
    lu["LandUse"] = key
    lu = lu[lu["LandUse"].notna()].copy()
    # Water is owned by build_water_layer(), which types every body as
    # WTR_RIV / WTR_LAKE / WTR_POND. Leaving OSM's own water tags in the
    # land-use layer duplicates those polygons under an UNTYPED key -- 94.8 %
    # of `natural:water` inside this domain was a second copy of a WTR_*
    # polygon. Two hazards, both silent: an untyped key falls through to
    # `fallback.low_ndvi` (asphalt -- a river rendered as pavement), and if it
    # rasterises after the typed copy it overwrites the type, which is how the
    # Hudson became a "lake" once already. Drop them here; nothing is lost.
    # NOT natural:wetland -- that is a real, separately mapped class
    # (vegetation:bogs_marshes), not open water, and dropping it would lose
    # the Hackensack marshes.
    _water_keys = ("natural:water", "natural:coastline", "natural:bay",
                   "natural:strait", "waterway:riverbank")
    _drop = lu["LandUse"].isin(_water_keys)
    if _drop.any():
        logger.info("nj: dropping %d OSM water polygons from the land-use layer "
                    "(build_water_layer owns water typing)", int(_drop.sum()))
        lu = lu[~_drop].copy()
    lu["LandUse_label"] = "New Jersey (OSM) " + lu["LandUse"].astype(str)

    bl = gpd.read_file(bl_path, layer="buildings").reset_index(drop=True)
    # palmgeo wants an int building id and a function code; OSM has neither in
    # the cadastre's vocabulary. 2100 is NYC Building Footprints' code for an
    # ordinary building, which is what mapping_us_ny.yml's by_gfk keys on.
    bl["bin"] = -(bl.index.to_numpy() + 1)   # negative: cannot collide with a real BIN
    bl["feature_code"] = "2100"
    return lu[["LandUse", "LandUse_label", "geometry"]], bl[["bin", "feature_code", "geometry"]]


def stage_cadastre() -> None:
    """Fetch MapPLUTO + Building Footprints over the AOI, fuse with the water
    layer above into one zipped GeoPackage with `landuse`/`buildings` layers.
    """
    import geopandas as gpd

    (STAGING / "cadastre").mkdir(parents=True, exist_ok=True)
    landuse = fetch_mappluto()
    landuse = landuse.to_crs(epsg=TARGET_EPSG)
    # Human-readable label per code, from the NYC DCP PLUTO Data Dictionary
    # (May 2026, 26v1) -- sources_us_ny.yml declares `label_field:
    # LandUse_label`, and `palmgeo report` prints it when listing vocabulary
    # that has no mapping rule yet, so it has to actually exist in the layer.
    landuse["LandUse_label"] = (
        landuse["LandUse"].map(MAPPLUTO_LANDUSE_LABELS).fillna("(unmapped MapPLUTO code)")
    )
    water = build_water_layer()
    cols = ["LandUse", "LandUse_label", "geometry"]
    landuse = gpd.GeoDataFrame(
        pd.concat([landuse[cols], water[cols]], ignore_index=True),
        crs=f"EPSG:{TARGET_EPSG}",
    )

    nj_lu, nj_bld = _nj_layers()
    if nj_lu is not None:
        landuse = gpd.GeoDataFrame(
            pd.concat([landuse[cols], nj_lu[cols]], ignore_index=True),
            crs=f"EPSG:{TARGET_EPSG}")
        logger.info("cadastre: + %d New Jersey land-use polygons (OSM)", len(nj_lu))

    # Clip SPATIALLY to the AOI, not by borough. The previous
    # `base_bbl.startswith("1")` filter kept Manhattan only and left every
    # building across the East River out of the driver.
    src = BUILDING_FOOTPRINTS_CACHE if BUILDING_FOOTPRINTS_CACHE.exists() else BUILDING_FOOTPRINTS_URL
    buildings = gpd.read_file(src, bbox=BBOX_WGS84)
    buildings = buildings.to_crs(epsg=TARGET_EPSG)

    if nj_bld is not None:
        keep = ["bin", "feature_code", "geometry"]
        buildings = gpd.GeoDataFrame(
            pd.concat([buildings[keep], nj_bld[keep]], ignore_index=True),
            crs=f"EPSG:{TARGET_EPSG}")
        logger.info("cadastre: + %d New Jersey building footprints (OSM)", len(nj_bld))

    gpkg_path = STAGING / "cadastre" / "nyc_cadastre.gpkg"
    landuse.to_file(gpkg_path, layer="landuse", driver="GPKG")
    buildings.to_file(gpkg_path, layer="buildings", driver="GPKG")

    # NOTE the stem: "nyc_cadastre.zip", not "nyc_cadastre.gpkg.zip".
    # Cache.unzip() extracts into <archive_dir>/<archive_stem>/, so the latter
    # would try to create a directory named exactly like the GeoPackage.
    zip_path = STAGING / "cadastre" / "nyc_cadastre.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(gpkg_path, gpkg_path.name)
    # Drop the loose .gpkg: palmgeo's Cache.unzip() extracts into a directory
    # named after the archive stem -- i.e. exactly "nyc_cadastre.gpkg"
    # -- and fails with FileExistsError if a file of that name sits beside
    # the archive.
    gpkg_path.unlink()
    logger.info("cadastre: fused %d landuse (incl. %d water) + %d building features -> %s",
                len(landuse), len(water), len(buildings), zip_path)


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", type=Path, required=True,
                    help="the run's project.workdir; staged files go to <workdir>/staging")
    ap.add_argument("--bbox", type=float, nargs=4, required=True,
                    metavar=("LON0", "LAT0", "LON1", "LAT1"),
                    help="area to stage in WGS84 degrees: all domains plus ~500 m margin")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("index", help="build all four staged tile indices (fast, no raster I/O)")
    p_ortho = sub.add_parser("ortho", help="stage real ortho tiles")
    p_ortho.add_argument("--sample", nargs="+", metavar="TILE_ID")
    p_ortho.add_argument("--all", action="store_true")
    sub.add_parser("dtm", help="stage the citywide DEM (15.6 GB download)")
    sub.add_parser("dsm", help="stage the citywide DSM (30.1 GB download)")
    sub.add_parser("cadastre", help="fuse MapPLUTO + Building Footprints + water")
    args = ap.parse_args()
    _configure(args.workdir, args.bbox)

    if args.cmd == "index":
        build_ortho_tile_index()
        build_ortho_index()
        build_dtm_dsm_index()
        build_cadastre_index()
    elif args.cmd == "ortho":
        stage_ortho_tiles(args.sample, run_all=args.all)
    elif args.cmd == "dtm":
        stage_dtm_dsm("dtm")
    elif args.cmd == "dsm":
        stage_dtm_dsm("dsm")
    elif args.cmd == "cadastre":
        stage_cadastre()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
