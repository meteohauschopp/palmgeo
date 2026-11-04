#!/usr/bin/env python3
"""Fill terrain/surface nodata inside the water mask with a flat water surface.

WHY: a survey that stops at an administrative line leaves the water beyond it
without elevation. NYC's 2017 topobathymetric survey, for which this was
written, stops at the city line down the middle of the Hudson. Over WATER that
gap has a physically correct answer: a hydro-flattened surface at the water
level -- 0 m NAVD88 there, because the rivers around Manhattan are tidal and
essentially at sea level (an East River sample of the NYC product reads
-0.35 m). Inland water sits higher; pass its level with --level.

This fills ONLY nodata that falls inside the water mask. Land gaps (the New
Jersey side) are deliberately left as nodata: filling those would invent
terrain, and the Palisades are 60-100 m cliffs, not sea level.

Block-wise, so memory stays bounded regardless of raster size.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import geometry_mask

#: Default water level, in the raster's vertical datum: sea level.
WATER_SURFACE_M = 0.0
NODATA = -9999.0


def fill(src_path: Path, dst_path: Path, water_gpkg: Path, layer: str = None,
         level: float = WATER_SURFACE_M) -> None:
    lu = gpd.read_file(water_gpkg, layer=layer) if layer else gpd.read_file(water_gpkg)
    # Match any WTR* code: the water layer is split by type (WTR_RIV /
    # WTR_LAKE / WTR_POND) so PALM gets river/lake/pond rather than one
    # blanket "river". An exact "WTR" test silently matched nothing after
    # that split and filled zero pixels.
    if "LandUse" in lu.columns:
        water = lu[lu["LandUse"].astype(str).str.startswith("WTR")]
    else:
        water = lu
    geoms = list(water.geometry)
    print(f"water polygons: {len(geoms)}")

    with rasterio.open(src_path) as src:
        profile = src.profile.copy()
        profile.update(compress="DEFLATE", tiled=True, BIGTIFF="IF_SAFER", nodata=NODATA)
        filled = 0
        with rasterio.open(dst_path, "w", **profile) as dst:
            for _, window in src.block_windows(1):
                data = src.read(1, window=window)
                gap = (data == src.nodata) | ~np.isfinite(data) | (data < -1e30)
                if gap.any():
                    wtr = ~geometry_mask(
                        geoms, out_shape=data.shape,
                        transform=src.window_transform(window), invert=False,
                    )
                    sel = gap & wtr
                    if sel.any():
                        data = np.where(sel, np.float32(level), data)
                        filled += int(sel.sum())
                data = np.where((data == src.nodata) | ~np.isfinite(data) | (data < -1e30),
                                np.float32(NODATA), data).astype("float32")
                dst.write(data, 1, window=window)
    print(f"filled {filled} px with a {level} m water surface -> {dst_path.name}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("src", type=Path, help="elevation raster with nodata gaps")
    ap.add_argument("dst", type=Path, help="output raster (must not be the input)")
    ap.add_argument("water", type=Path,
                    help="vector layer carrying the water polygons")
    ap.add_argument("--layer", default=None,
                    help="layer name inside the water file, if it has several")
    ap.add_argument("--level", type=float, default=WATER_SURFACE_M,
                    help="water surface in the raster's vertical datum, m "
                         f"(default {WATER_SURFACE_M:g}: sea level)")
    args = ap.parse_args()
    if args.src.resolve() == args.dst.resolve():
        ap.error("src and dst must differ -- this reads and writes block by "
                 "block and cannot do so in place")
    fill(args.src, args.dst, args.water, args.layer, args.level)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
