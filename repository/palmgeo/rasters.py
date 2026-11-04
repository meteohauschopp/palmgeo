"""Raster I/O helpers.

Everything here rests on one assumption: **source and destination share a
single projected, north-up CRS.** palmgeo enforces it rather than coping with
violations -- ``config.load`` refuses a run whose ``grid.epsg`` differs from
the CRS the region file declares, and a region whose products genuinely differ
is expected to reproject into a staging area first (see PORTING.md). It is not
a property of any particular data provider; it is a precondition this module
is built on.

The reason it is worth enforcing is a single performance property: under that
assumption a *decimated windowed read* is possible. Reading a 2 km / 20 cm
orthophoto tile -- 10000 x 10000 x 4 bands, 400 MB uncompressed -- down to a
1 km / 1 m work grid then touches only the matching overview level and moves a
few megabytes. Without it, every read would go through a reprojection that
materialises the whole raster, and a city-scale run would move terabytes
instead of tens of gigabytes.

That is also why the destination grid is decided before the read rather than
after: the caller passes the exact transform and shape it wants, and this
module either satisfies it with a windowed read or says so.
"""

from __future__ import annotations

import logging
import math
import os
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Iterator, Optional, Sequence, Tuple

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import Affine, array_bounds
from rasterio.windows import from_bounds

from palmgeo.grid import Bounds

logger = logging.getLogger("palmgeo.rasters")

DEFAULT_COG: Dict[str, Any] = {
    "compress": "DEFLATE",
    "predictor": 2,
    "blocksize": 512,
    "overviews": True,
}

RESAMPLING = {
    "nearest": Resampling.nearest,
    "bilinear": Resampling.bilinear,
    "cubic": Resampling.cubic,
    "average": Resampling.average,
    "mode": Resampling.mode,
    "max": Resampling.max,
    "min": Resampling.min,
    "med": Resampling.med,
    "q1": Resampling.q1,
    "q3": Resampling.q3,
}


@contextmanager
def gdal_env(cache_mb: int = 512, **extra: Any) -> Iterator[None]:
    """Configure GDAL for cloud-optimised access within the block."""
    # rasterio.Env passes values straight to GDAL's config API, which is typed:
    # booleans must be bool and counters must be int, not their string spellings.
    env: Dict[str, Any] = {
        "GDAL_CACHEMAX": int(cache_mb),
        # Never list a bucket "directory" just to open one object.
        "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
        "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif,.TIF",
        "VSI_CACHE": True,
        "VSI_CACHE_SIZE": 64 * 1024 * 1024,
        "GDAL_HTTP_MAX_RETRY": 3,
        "GDAL_HTTP_RETRY_DELAY": 2,
        "GDAL_NUM_THREADS": "ALL_CPUS",
    }
    env.update(extra)
    with rasterio.Env(**env):
        yield


def vsi_path(url_or_path: str | Path) -> str:
    """Turn an http(s) URL into a GDAL /vsicurl/ path; pass through local paths."""
    text = str(url_or_path)
    if text.startswith(("http://", "https://")):
        return f"/vsicurl/{text}"
    return text


def _predictor(dtype: Any, opts: Dict[str, Any]) -> Dict[str, Any]:
    """Pick a valid TIFF predictor for the data type.

    GDAL accepts horizontal differencing (2) only for integers and needs
    floating-point predictor (3) for float bands; passing 2 for a float band is
    silently dropped, losing a sizeable chunk of the compression ratio on
    smooth fields such as terrain height.
    """
    if str(opts.get("compress", "")).upper() not in {"DEFLATE", "LZW", "ZSTD"}:
        return {}
    if np.issubdtype(dtype, np.floating):
        return {"predictor": 3}
    if np.issubdtype(dtype, np.integer):
        return {"predictor": int(opts.get("predictor", 2))}
    return {}


def _intersect(a: Bounds, b: Bounds) -> Optional[Bounds]:
    out = (max(a[0], b[0]), max(a[1], b[1]), min(a[2], b[2]), min(a[3], b[3]))
    if out[2] <= out[0] or out[3] <= out[1]:
        return None
    return out


def read_on_grid(
    source: str | Path,
    dst_transform: Affine,
    dst_shape: Tuple[int, int],
    band: int = 1,
    resampling: Resampling = Resampling.average,
    dtype: str = "float32",
    fill: float = np.nan,
) -> np.ndarray:
    """Read ``source`` resampled onto an explicit destination grid.

    The destination block is snapped *inward* to whole destination pixels so
    that the read window lies strictly inside the source; the uncovered margin
    keeps ``fill``. No `boundless` read and no full-array materialisation.
    """
    height, width = dst_shape
    out = np.full((height, width), fill, dtype=dtype)
    res_x = dst_transform.a
    res_y = -dst_transform.e
    left, bottom, right, top = array_bounds(height, width, dst_transform)

    with rasterio.open(vsi_path(source)) as src:
        inter = _intersect((left, bottom, right, top), tuple(src.bounds))
        if inter is None:
            return out
        # inward snap to destination pixel edges
        c0 = int(math.ceil((inter[0] - left) / res_x - 1e-6))
        c1 = int(math.floor((inter[2] - left) / res_x + 1e-6))
        r0 = int(math.ceil((top - inter[3]) / res_y - 1e-6))
        r1 = int(math.floor((top - inter[1]) / res_y + 1e-6))
        c0, r0 = max(0, c0), max(0, r0)
        c1, r1 = min(width, c1), min(height, r1)
        if c1 <= c0 or r1 <= r0:
            return out

        block = (
            left + c0 * res_x,
            top - r1 * res_y,
            left + c1 * res_x,
            top - r0 * res_y,
        )
        window = from_bounds(*block, transform=src.transform)
        data = src.read(
            band,
            window=window,
            out_shape=(r1 - r0, c1 - c0),
            resampling=resampling,
            masked=True,
        )
        # Cast before filling, not after: the source keeps its own dtype, so an
        # 8-bit orthophoto that declares a nodata value comes back as a masked
        # uint8 array, and filling that with NaN raises. Whether a tile is
        # masked at all depends on the individual product -- of the DOP20
        # epochs over the Hannover root, some campaigns set nodata and some do
        # not, so the same code path works on one tile and fails on the next.
        if np.ma.isMaskedArray(data):
            values = np.ma.getdata(data).astype(dtype, copy=True)
            mask = np.ma.getmaskarray(data)
            if mask.any():
                values[mask] = fill
        else:
            values = np.asarray(data).astype(dtype, copy=False)
        out[r0:r1, c0:c1] = values
    return out


def write_raster(
    path: Path,
    data: np.ndarray,
    transform: Affine,
    crs: str,
    nodata: Optional[float] = None,
    options: Optional[Dict[str, Any]] = None,
    overviews: Optional[bool] = None,
) -> Path:
    """Write a tiled, compressed (optionally cloud-optimised) GeoTIFF."""
    opts = {**DEFAULT_COG, **(options or {})}
    path.parent.mkdir(parents=True, exist_ok=True)
    if data.ndim == 2:
        data = data[np.newaxis, ...]
    count, height, width = data.shape
    blocksize = int(opts.get("blocksize", 512))
    # GDAL requires tile dimensions to be multiples of 16 and <= the raster
    tiled = width >= blocksize and height >= blocksize
    profile: Dict[str, Any] = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": count,
        "dtype": data.dtype.name,
        "crs": crs,
        "transform": transform,
        "compress": opts.get("compress", "DEFLATE"),
        "tiled": tiled,
        "BIGTIFF": "IF_SAFER",
    }
    if tiled:
        profile["blockxsize"] = blocksize
        profile["blockysize"] = blocksize
    profile.update(_predictor(data.dtype, opts))
    if nodata is not None:
        profile["nodata"] = nodata

    with rasterio.open(path, "w", **profile) as dst:
        dst.write(data)
        want_ovr = opts.get("overviews", True) if overviews is None else overviews
        if want_ovr and min(width, height) >= 2 * blocksize:
            factors = [f for f in (2, 4, 8, 16, 32) if min(width, height) // f >= 64]
            if factors:
                dst.build_overviews(factors, Resampling.average)
                dst.update_tags(ns="rio_overview", resampling="average")
    return path


def build_vrt(sources: Sequence[Path], out: Path, nodata: Optional[float] = None) -> Path:
    """Build a plain GDAL VRT mosaic over ``sources``.

    Written by hand instead of shelling out to ``gdalbuildvrt`` so the pipeline
    has no dependency on GDAL command line tools being on PATH and matching the
    Python bindings' PROJ version.
    """
    sources = [Path(s) for s in sources if Path(s).exists()]
    if not sources:
        raise ValueError("build_vrt: no existing source rasters")

    with rasterio.open(sources[0]) as ref:
        crs = ref.crs
        dtype = ref.dtypes[0]
        count = ref.count
        res_x = ref.transform.a
        res_y = -ref.transform.e
        nodata = ref.nodata if nodata is None else nodata

    # One pass over the sources; a tier can hold thousands of tiles and
    # re-opening each of them per band would dominate the runtime.
    infos = []
    for src_path in sources:
        with rasterio.open(src_path) as src:
            infos.append((src_path, src.bounds, src.width, src.height))
    left = min(i[1].left for i in infos)
    bottom = min(i[1].bottom for i in infos)
    right = max(i[1].right for i in infos)
    top = max(i[1].top for i in infos)
    width = int(round((right - left) / res_x))
    height = int(round((top - bottom) / res_y))

    gdal_dtype = {
        "uint8": "Byte", "int8": "Int8", "uint16": "UInt16", "int16": "Int16",
        "uint32": "UInt32", "int32": "Int32", "float32": "Float32",
        "float64": "Float64",
    }[dtype]

    lines = [
        f'<VRTDataset rasterXSize="{width}" rasterYSize="{height}">',
        f"  <SRS>{crs.to_wkt()}</SRS>",
        f"  <GeoTransform>{left}, {res_x}, 0.0, {top}, 0.0, {-res_y}</GeoTransform>",
    ]
    for band in range(1, count + 1):
        lines.append(f'  <VRTRasterBand dataType="{gdal_dtype}" band="{band}">')
        if nodata is not None:
            lines.append(f"    <NoDataValue>{nodata}</NoDataValue>")
        for src_path, sb, sw, sh in infos:
            x_off = int(round((sb.left - left) / res_x))
            y_off = int(round((top - sb.top) / res_y))
            rel = os.path.relpath(src_path, out.parent)
            lines += [
                '    <ComplexSource>',
                f'      <SourceFilename relativeToVRT="1">{rel}</SourceFilename>',
                f"      <SourceBand>{band}</SourceBand>",
                f'      <SrcRect xOff="0" yOff="0" xSize="{sw}" ySize="{sh}"/>',
                f'      <DstRect xOff="{x_off}" yOff="{y_off}" xSize="{sw}" ySize="{sh}"/>',
            ]
            if nodata is not None:
                lines.append(f"      <NODATA>{nodata}</NODATA>")
            lines.append("    </ComplexSource>")
        lines.append("  </VRTRasterBand>")
    lines.append("</VRTDataset>")

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def resample_to_domain(
    source: str | Path,
    out_path: Path,
    dst_transform: Affine,
    dst_shape: Tuple[int, int],
    crs: str,
    resampling: Resampling,
    dtype: str,
    nodata: Optional[float],
    options: Optional[Dict[str, Any]] = None,
    chunk_rows: int = 2048,
) -> Path:
    """Resample a (possibly very large) mosaic onto a domain grid, row-chunked.

    Peak memory is ``chunk_rows * width`` destination cells plus whatever the
    matching source window needs, so a 1500 x 1500 domain and a 40 000 x 40 000
    work mosaic cost the same.
    """
    height, width = dst_shape
    fill = nodata if nodata is not None else (np.nan if dtype.startswith("float") else 0)
    profile: Dict[str, Any] = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": 1,
        "dtype": dtype,
        "crs": crs,
        "transform": dst_transform,
        "compress": (options or DEFAULT_COG).get("compress", "DEFLATE"),
        "tiled": width >= 512 and height >= 512,
        "BIGTIFF": "IF_SAFER",
    }
    if profile["tiled"]:
        profile["blockxsize"] = profile["blockysize"] = 512
    if nodata is not None:
        profile["nodata"] = nodata
    profile.update(_predictor(np.dtype(dtype), options or DEFAULT_COG))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **profile) as dst:
        for r0 in range(0, height, chunk_rows):
            r1 = min(height, r0 + chunk_rows)
            chunk_transform = dst_transform * Affine.translation(0, r0)
            block = read_on_grid(
                source,
                chunk_transform,
                (r1 - r0, width),
                resampling=resampling,
                dtype=dtype,
                fill=fill,
            )
            dst.write(block.astype(dtype, copy=False), 1, window=((r0, r1), (0, width)))
    return out_path
