#!/usr/bin/env python3
"""Convert a US-survey-foot elevation raster to metres, block by block.

WHY: NYC's 2017 LiDAR DEM/DSM store elevation in US survey FEET (their CRS,
EPSG:2263, is "NAD83 / New York Long Island (ftUS)"). gdalwarp reprojects
coordinates but never touches pixel VALUES, so a warped-to-EPSG:6347 raster
still carries feet. palm_csd and PALM are metric throughout -- feeding the
unconverted raster in would inflate terrain by 3.28x, silently.

Memory is bounded by construction: one internal block row at a time, never
the whole raster. This is not caution for its own sake -- an earlier attempt
at the clip/resample step confused metres with feet in its decimation factor,
asked for ~25 GB in one allocation and was killed by the kernel, taking the
interactive session that launched it along. repository/tools/run_guarded.sh
exists for the same reason.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

US_SURVEY_FOOT = 1200.0 / 3937.0  # 0.3048006096..., NOT the international 0.3048
OUT_NODATA = -9999.0


def convert(src_path: Path, dst_path: Path) -> None:
    with rasterio.open(src_path) as src:
        profile = src.profile.copy()
        profile.update(dtype="float32", nodata=OUT_NODATA, compress="DEFLATE",
                        tiled=True, BIGTIFF="IF_SAFER")
        src_nodata = src.nodata
        with rasterio.open(dst_path, "w", **profile) as dst:
            # iterate the source's own blocks -> bounded memory
            for _, window in src.block_windows(1):
                data = src.read(1, window=window)
                if src_nodata is not None:
                    mask = (data == src_nodata) | ~np.isfinite(data)
                else:
                    mask = ~np.isfinite(data)
                # the topobathy DEM also uses very large negative values as
                # fill; treat anything below -1e30 as nodata too
                mask |= data < -1e30
                out = data.astype("float32") * np.float32(US_SURVEY_FOOT)
                out[mask] = OUT_NODATA
                dst.write(out, 1, window=window)
    print(f"converted {src_path.name} -> {dst_path.name} (ftUS -> m)")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("src", type=Path, help="raster whose VALUES are US survey feet")
    ap.add_argument("dst", type=Path, help="output raster in metres")
    args = ap.parse_args()
    if args.src.resolve() == args.dst.resolve():
        ap.error("src and dst must differ -- this reads and writes block by "
                 "block and cannot do so in place")
    convert(args.src, args.dst)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
