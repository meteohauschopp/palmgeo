#!/usr/bin/env python3
"""Stage New Jersey sources where a US-NY area crosses the Hudson.

WHY: New York City's own data stops at the city line, which runs down the
middle of the Hudson. Once a domain reaches the New Jersey bank (Jersey
City, Hoboken, the Palisades) every role has to be filled from somewhere
else, or that part of the domain is empty -- which is exactly what the
first two drivers of the Manhattan run looked like.

Usage, after stage_sources.py has staged the New York side into the same
work directory with the same --bbox (both warp onto its envelope, and the
merge needs identical grids):

    NJ="python repository/ports/us_ny/stage_nj_supplement.py \
        --workdir static_driver/<run>/work --bbox <lon0> <lat0> <lon1> <lat1> \
        --nj-bbox <lon0> <lat0> <lon1> <lat1>"
    $NJ dtm && $NJ dsm      # 3DEP terrain and surface over the box
    $NJ merge               # fill NYC nodata, with the measured datum offset
    $NJ naip                # list the NAIP tiles for the ortho mosaic

``--nj-bbox`` is the New Jersey part of the box, for the OSM and NAIP
queries. Mind PORTING.md section 5b before drawing it: the state line follows
the river, and a rectangle cut at one longitude loses New Jersey further
north.

What this adds, and why each source:

  terrain   USGS 3DEP lidar DTM, 2 m   \\  same NJ_SdL5_2014 collection, so the
  surface   USGS 3DEP lidar DSM, 2 m   /  pair is date-matched the way NYC's
                                         2017 DEM/DSM pair is -- nDSM from a
                                         mixed-epoch pair is not trustworthy.
  imagery   NAIP 2023, 0.3 m, 4-band      NYSDOP is a New York State
                                          programme and simply stops at the
                                          state line. NAIP is national and
                                          carries the NIR band NDVI needs.
  buildings OpenStreetMap `building=*`    ~16 000 footprints on the NJ strip.
  land use  OpenStreetMap `landuse=*`     MapPLUTO is NYC tax lots only.

UNITS -- the trap this project already fell into once: NYC's rasters are US
survey FEET and need converting, but 3DEP and NAIP are METRES (verified:
EPSG:26918 + NAVD88, both in metre). Do NOT run ft_to_m.py on these.

Planetary Computer assets need a SAS token (unsigned reads return HTTP 409);
`_signed()` fetches one per container and caches it for the run.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # ports/
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # palmgeo
from staging_index import aoi_bounds  # noqa: E402
from palmgeo import USER_AGENT  # noqa: E402

logger = logging.getLogger("stage_nj")

#: The CRS sources_us_ny.yml declares; must match stage_sources.py.
TARGET_EPSG = 6347

# Run settings, set from the command line by _configure().
STAGING = Path()
#: The staged area in the target CRS -- identical to stage_sources.py's.
AOI_TARGET_BOUNDS: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)
#: The New Jersey part of the area, WGS84: the OSM and NAIP queries.
NJ_BBOX_WGS84: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)
#: Raster search box: the area widened by this many degrees on every side.
#: Deliberately far wider than the AOI. The 3DEP collection's catalogue
#: bboxes are wrong AND its tile NAMES are offset differently between the
#: dtm and dsm sub-collections (the tiles covering the Manhattan run were
#: 6-8..6-11 in dtm but 6-3..6-5 in dsm, ~60 km apart under identical ids).
#: Only the real bounds read from each file can be trusted, so search broadly
#: and let the overlap filter below do the work.
STAC_SEARCH_MARGIN_DEG = 0.5
NJ_STAC_SEARCH_BBOX: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)


def _configure(workdir: Path, bbox: list[float], nj_bbox: list[float]) -> None:
    """Set the run settings from the command line."""
    global STAGING, AOI_TARGET_BOUNDS, NJ_BBOX_WGS84, NJ_STAC_SEARCH_BBOX
    for name, box in (("--bbox", bbox), ("--nj-bbox", nj_bbox)):
        lon0, lat0, lon1, lat1 = box
        if not (-180 <= lon0 < lon1 <= 180 and -90 <= lat0 < lat1 <= 90):
            raise SystemExit(f"{name} {box}: expected lon0 lat0 lon1 lat1 in WGS84 degrees")
    workdir = (workdir if workdir.is_absolute() else Path.cwd() / workdir).resolve()
    STAGING = workdir / "staging"
    AOI_TARGET_BOUNDS = aoi_bounds(tuple(bbox), TARGET_EPSG)
    NJ_BBOX_WGS84 = tuple(float(v) for v in nj_bbox)
    m = STAC_SEARCH_MARGIN_DEG
    NJ_STAC_SEARCH_BBOX = (bbox[0] - m, bbox[1] - m, bbox[2] + m, bbox[3] + m)

STAC_SEARCH = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
SAS_TOKEN = "https://planetarycomputer.microsoft.com/api/sas/v1/token/{account}/{container}"
OVERPASS = "https://overpass-api.de/api/interpreter"
_token_cache: dict[str, str] = {}


def _signed(href: str) -> str:
    """Planetary Computer blobs are private; sign per container."""
    rest = href.split("://", 1)[1]
    host, path = rest.split("/", 1)
    account = host.split(".")[0]
    container = path.split("/", 1)[0]
    key = f"{account}/{container}"
    if key not in _token_cache:
        r = requests.get(SAS_TOKEN.format(account=account, container=container), timeout=60)
        r.raise_for_status()
        _token_cache[key] = r.json()["token"]
    return f"{href}?{_token_cache[key]}"


def _stac_items(collection: str, bbox, datetime=None) -> list[dict]:
    body = {"collections": [collection], "bbox": list(bbox), "limit": 100}
    if datetime:
        body["datetime"] = datetime
    r = requests.post(STAC_SEARCH, json=body, timeout=120)
    r.raise_for_status()
    return r.json().get("features", [])


def stage_elevation(role: str) -> None:
    """Warp the NJ 3DEP DTM/DSM onto the AOI grid (metres, no conversion).

    TWO TRAPS, both hit for real on 2026-09-07:

    1. **The STAC bboxes in this collection are wrong.** A search for a 1 km
       box over Hoboken returned a tile whose actual georeferencing put it
       40 km south, and the project name ("NJ_SdL5") is no guide either --
       several of its tiles sit over Manhattan. Tiles are therefore selected
       by opening each one and reading its real bounds, never by the
       catalogue's geometry.
    2. **Warping several tiles in one gdalwarp call produced an empty
       raster**, although each tile warped correctly on its own. Building a
       VRT first and warping that once is both faster and reliable -- the
       same pattern the ortho mosaic already uses.
    """
    import rasterio

    collection = {"dtm": "3dep-lidar-dtm", "dsm": "3dep-lidar-dsm"}[role]
    items = _stac_items(collection, NJ_STAC_SEARCH_BBOX)
    if not items:
        raise RuntimeError(f"no {collection} items returned for the NJ search box")

    x0, y0, x1, y1 = AOI_TARGET_BOUNDS
    keep = []
    for it in items:
        href = _signed(it["assets"]["data"]["href"])
        try:
            with rasterio.open(f"/vsicurl/{href}") as s:
                b = s.bounds
        except Exception as exc:  # noqa: BLE001
            logger.warning("%s: cannot open %s (%s)", role, it["id"], type(exc).__name__)
            continue
        if not (b.right < x0 or b.left > x1 or b.top < y0 or b.bottom > y1):
            keep.append((it["id"], href))
    if not keep:
        raise RuntimeError(
            f"{collection}: none of the {len(items)} candidate tiles actually "
            f"overlaps AOI_TARGET_BOUNDS -- check the tiles' real bounds"
        )
    logger.info("%s: %d of %d tiles really overlap the AOI", role, len(keep), len(items))

    out = STAGING / role / f"nj_{role}_epsg{TARGET_EPSG}.tif"
    out.parent.mkdir(parents=True, exist_ok=True)
    vrt = STAGING / role / f"nj_{role}.vrt"
    listing = STAGING / role / f"nj_{role}_tiles.txt"
    listing.write_text("\n".join(f"/vsicurl/{h}" for _, h in keep))
    subprocess.run(["gdalbuildvrt", "-input_file_list", str(listing), str(vrt)], check=True)
    subprocess.run([
        "gdalwarp", "-t_srs", f"EPSG:{TARGET_EPSG}",
        "-te", *[str(v) for v in AOI_TARGET_BOUNDS],
        "-tr", "1", "1", "-r", "average",
        "-of", "GTiff", "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES",
        "-co", "BIGTIFF=IF_SAFER", "-srcnodata", "-9999", "-dstnodata", "-9999",
        "-wm", "512", "--config", "GDAL_CACHEMAX", "512", "-overwrite",
        str(vrt), str(out),
    ], check=True)

    with rasterio.open(out) as s:
        import numpy as np
        a = s.read(1, out_shape=(600, 430))
        valid = float(((a != -9999) & np.isfinite(a)).mean())
    logger.info("%s: staged NJ elevation -> %s (%.1f %% of the AOI has data)",
                role, out.name, 100 * valid)
    if valid <= 0:
        raise RuntimeError(f"{role}: the staged NJ raster is empty -- refusing to "
                           f"pretend this worked")


MIN_OVERLAP_PX = 20_000       # below this the median offset is not trustworthy
MAX_OFFSET_IQR_M = 1.5        # above this it is not a constant shift


def measure_datum_offset() -> float:
    """Measure the vertical offset between the NYC and NJ terrain, in metres.

    WHY THIS IS NOT OPTIONAL: the two surveys do not share a vertical datum.
    NYC's 2017 product is NAVD88 orthometric; the 3DEP tiles staged here come
    out ~31.7 m lower, which is the geoid separation over New York -- i.e.
    they are ellipsoidal heights. Merged untouched, Hoboken sits at -32 m and
    the Hudson's west bank becomes a 32 m cliff in the PALM terrain, while
    every structural check still passes because the raster is perfectly
    well-formed.

    Measured, not assumed. The wide STAC search (see NJ_STAC_SEARCH_BBOX)
    pulls in tiles that also cover New York, so the two models genuinely
    overlap on land and the shift can be read off directly. Measured
    2026-09-07: median +31.65 m over 103,267 px with an IQR of 0.31 m -- a
    tight spread around a large median is the signature of a datum shift
    rather than noise.

    The offset is always taken from the DTM pair and then applied to the DSM
    too. It is a property of the datum, not of the surface, and the DSM's own
    overlap is contaminated by real change between the epochs (roofs, canopy):
    its IQR is 1.09 m against the DTM's 0.31 m. Using one shift for both also
    leaves nDSM = DSM - DTM untouched, so building heights stay correct.
    """
    import numpy as np
    import rasterio

    nyc = STAGING / "dtm" / f"nyc_dtm_epsg{TARGET_EPSG}_m.tif"
    nj = STAGING / "dtm" / f"nj_dtm_epsg{TARGET_EPSG}.tif"
    with rasterio.open(nyc) as a, rasterio.open(nj) as b:
        sc = 8
        h, w = a.height // sc, a.width // sc
        va = a.read(1, out_shape=(h, w)).astype("f8")
        vb = b.read(1, out_shape=(h, w)).astype("f8")
    ok = (va != -9999) & (vb != -9999) & np.isfinite(va) & np.isfinite(vb)
    n = int(ok.sum())
    if n < MIN_OVERLAP_PX:
        raise RuntimeError(
            f"only {n:,} px where both terrain models have data -- too few to "
            f"measure the vertical datum offset. Widen NJ_STAC_SEARCH_BBOX so "
            f"the NJ tiles genuinely overlap the NYC survey, or supply the "
            f"offset explicitly; do NOT merge unshifted."
        )
    d = va[ok] - vb[ok]
    q25, med, q75 = (float(v) for v in np.percentile(d, [25, 50, 75]))
    iqr = q75 - q25
    logger.info("datum offset NYC-NJ: median %+.2f m over %d px (IQR %.2f m)",
                med, n, iqr)
    if iqr > MAX_OFFSET_IQR_M:
        raise RuntimeError(
            f"the NYC-NJ difference has an IQR of {iqr:.2f} m around a median "
            f"of {med:+.2f} m -- that is not a constant datum shift. Check the "
            f"horizontal alignment before applying any vertical correction."
        )
    return med


def merge_into_nyc(role: str, offset_m: float | None = None) -> None:
    """Fill the NYC raster's nodata with the NJ raster, block by block.

    `offset_m` is added to every NJ value first -- see measure_datum_offset()
    for why this is mandatory and why it is measured on the DTM for both
    roles. Pass None to measure it here.
    """
    import numpy as np
    import rasterio

    if offset_m is None:
        offset_m = measure_datum_offset()

    nyc = STAGING / role / f"nyc_{role}_epsg{TARGET_EPSG}_m.tif"
    nj = STAGING / role / f"nj_{role}_epsg{TARGET_EPSG}.tif"
    out = STAGING / role / f"merged_{role}_epsg{TARGET_EPSG}_m.tif"
    with rasterio.open(nyc) as a, rasterio.open(nj) as b:
        if (a.width, a.height) != (b.width, b.height):
            raise RuntimeError(
                f"grids differ: NYC {a.width}x{a.height} vs NJ {b.width}x{b.height} "
                f"-- both must be warped to AOI_TARGET_BOUNDS at the same -tr"
            )
        profile = a.profile.copy()
        profile.update(compress="DEFLATE", tiled=True, BIGTIFF="IF_SAFER", nodata=-9999.0)
        filled = 0
        with rasterio.open(out, "w", **profile) as dst:
            for _, win in a.block_windows(1):
                va = a.read(1, window=win)
                gap = (va == -9999.0) | ~np.isfinite(va)
                if gap.any():
                    vb = b.read(1, window=win)
                    ok = gap & (vb != -9999.0) & np.isfinite(vb)
                    if ok.any():
                        va = np.where(ok, vb + np.float32(offset_m), va)
                        filled += int(ok.sum())
                dst.write(va.astype("float32"), 1, window=win)
    logger.info("%s: filled %d px of NYC nodata from the NJ model "
                "(shifted %+.2f m) -> %s", role, filled, offset_m, out.name)


def fetch_osm(kind: str) -> Path:
    """Buildings or land use for the NJ strip, as GeoJSON-ish Overpass JSON."""
    q = {
        "buildings": '(way["building"]({bbox});relation["building"]({bbox}););',
        "landuse": '(way["landuse"]({bbox});relation["landuse"]({bbox});'
                    'way["natural"]({bbox});way["leisure"]({bbox}););',
    }[kind]
    bbox = f"{NJ_BBOX_WGS84[1]},{NJ_BBOX_WGS84[0]},{NJ_BBOX_WGS84[3]},{NJ_BBOX_WGS84[2]}"
    query = f"[out:json][timeout:180];{q.format(bbox=bbox)}out geom;"
    out = STAGING / f"osm_nj_{kind}.json"
    logger.info("osm %s: querying Overpass for the NJ strip", kind)
    # Overpass rejects the default python-requests User-Agent with HTTP 406
    # Not Acceptable (verified: same query, same body, only the UA differs ->
    # 406 vs 200). Identify the client properly, as its usage policy asks.
    r = requests.post(OVERPASS, data=query, timeout=300,
                      headers={"User-Agent": USER_AGENT})
    r.raise_for_status()
    out.write_text(r.text)
    n = len(r.json().get("elements", []))
    logger.info("osm %s: %d elements -> %s", kind, n, out)
    return out


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["dtm", "dsm", "merge-dtm", "merge-dsm",
                                     "merge", "datum-offset",
                                     "osm-buildings", "osm-landuse", "naip"])
    ap.add_argument("--workdir", type=Path, required=True,
                    help="the run's project.workdir, as given to stage_sources.py")
    ap.add_argument("--bbox", type=float, nargs=4, required=True,
                    metavar=("LON0", "LAT0", "LON1", "LAT1"),
                    help="the whole staged area, exactly as given to stage_sources.py")
    ap.add_argument("--nj-bbox", type=float, nargs=4, required=True,
                    metavar=("LON0", "LAT0", "LON1", "LAT1"),
                    help="the New Jersey part of it, for the OSM and NAIP queries")
    args = ap.parse_args()
    _configure(args.workdir, args.bbox, args.nj_bbox)
    if args.step in ("dtm", "dsm"):
        stage_elevation(args.step)
    elif args.step == "datum-offset":
        print(f"{measure_datum_offset():+.3f}")
    elif args.step == "merge":
        # Both roles in one go, sharing ONE measured offset -- see
        # measure_datum_offset(). This is the intended entry point.
        dz = measure_datum_offset()
        for r in ("dtm", "dsm"):
            merge_into_nyc(r, dz)
    elif args.step.startswith("merge-"):
        merge_into_nyc(args.step.split("-", 1)[1])
    elif args.step == "osm-buildings":
        fetch_osm("buildings")
    elif args.step == "osm-landuse":
        fetch_osm("landuse")
    elif args.step == "naip":
        items = _stac_items("naip", NJ_BBOX_WGS84, "2022-01-01/2025-12-31")
        # newest first, so the 2023 30 cm NJ coverage wins over older tiles
        items.sort(key=lambda i: i["properties"].get("datetime") or "", reverse=True)
        hrefs = [_signed(i["assets"]["image"]["href"]) for i in items]
        listing = STAGING / "ortho" / "naip_nj_sources.txt"
        listing.parent.mkdir(parents=True, exist_ok=True)
        listing.write_text("\n".join(f"/vsicurl/{h}" for h in hrefs) + "\n")
        logger.info("naip: %d tiles listed for the NJ strip -> %s", len(hrefs), listing)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
