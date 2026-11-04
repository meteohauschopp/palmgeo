#!/usr/bin/env python3
"""Sanity-check a staged ortho mosaic BEFORE it is classified.

WHY THIS EXISTS
---------------
Every ortho bug this project has hit was silent: the command exited 0 and
produced a plausible-looking GeoTIFF. Each one cost a full pipeline run plus a
driver review to notice.

  * **NIR replaced by a synthetic alpha band.** The staged tiles label band 4
    ColorInterp=Alpha although it is NIR. Read without `-nosrcalpha`, GDAL
    substitutes an opaque 255, NDVI saturates near +1 and *everything* -- roofs
    included -- classifies as vegetation. Central Park came out 97.6 %
    vegetated and looked, at a glance, correct.
  * **A whole source silently dropped.** `gdalbuildvrt` adopts the first
    input's CRS and merely *warns* while skipping every input that disagrees.
    Mixing NAIP (EPSG:26918) with tiles in EPSG:6347 dropped all 605 tiles and
    yielded a mosaic with no imagery over New York City whatsoever.
  * **Tile-grid mismatch.** 762 m source tiles against 1 km work tiles left
    ~80 % of the work tiles mostly black.

The three checks below correspond one-to-one to those failures. Run this on
the staged mosaic; do not wait for the driver to look wrong.

    python repository/tools/check_ortho.py \\
        --ortho static_driver/<run>/work/staging/ortho/<mosaic>.tif \\
        --landuse static_driver/<run>/work/staging/cadastre/... \\
        --veg-query "LandUse == '09'" --built-query "LandUse in ['03','04','05']"

Exit code 1 if any check FAILs, so it can gate a staging script.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import rasterio
from rasterio.enums import Resampling

# NDVI a real 4-band leaf-on scene should never exceed over dense built-up.
# The alpha-band failure drives this to ~+1; a genuine city median sits near 0.
BUILT_NDVI_ALARM = 0.30


def _read(src, band, shape):
    return src.read(band, out_shape=shape, resampling=Resampling.average).astype("f4")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ortho", required=True)
    ap.add_argument("--nir-band", type=int, default=4,
                    help="1-based NIR band (default 4: R, G, B, NIR order)")
    ap.add_argument("--red-band", type=int, default=1,
                    help="1-based red band (default 1)")
    ap.add_argument("--landuse", help="vector file with a land-use column (optional)")
    ap.add_argument("--veg-query", help="pandas .query() selecting parks/green")
    ap.add_argument("--built-query", help="pandas .query() selecting dense built-up")
    ap.add_argument("--grid", type=int, default=8,
                    help="split the extent into GRID x GRID cells for the "
                         "coverage check (default 8)")
    ap.add_argument("--max-empty-cell", type=float, default=0.50,
                    help="FAIL if any cell is emptier than this (default 0.50)")
    ap.add_argument("--sample", type=int, default=2000, help="analysis size [px]")
    args = ap.parse_args()

    fails: list[str] = []
    with rasterio.open(args.ortho) as src:
        sc = max(src.width, src.height) / args.sample
        shape = (max(1, int(src.height / sc)), max(1, int(src.width / sc)))
        red = _read(src, args.red_band, shape)
        nir = _read(src, args.nir_band, shape)
        transform = src.transform * src.transform.scale(src.width / shape[1],
                                                        src.height / shape[0])
        bounds = src.bounds
        print(f"{args.ortho}\n  {src.count} bands, {src.width} x {src.height} px, "
              f"analysed at {shape[1]} x {shape[0]}")

        # -- check 1: is band `nir` actually data, or a synthetic alpha plane? --
        uniq = np.unique(nir[::7, ::7])
        print(f"\n[1] NIR band {args.nir_band}: min {nir.min():.0f} max {nir.max():.0f}, "
              f"{uniq.size} distinct values in the sample")
        if uniq.size <= 2 and nir.max() >= 254:
            fails.append(f"band {args.nir_band} is constant/near-constant at "
                         f"{nir.max():.0f} -- this is an ALPHA plane, not NIR. "
                         f"Re-warp with -nosrcalpha.")
            print("    FAIL  alpha plane, not NIR")
        else:
            print("    ok    varies like real imagery")

    data = (red > 0) | (nir > 0)
    print(f"\n[2] coverage: {100 * data.mean():.1f} % of the extent has imagery")

    # -- check 2: is any REGION empty? a dropped source leaves a clean block --
    H, W = data.shape
    worst = []
    for r in range(args.grid):
        for c in range(args.grid):
            cell = data[r * H // args.grid:(r + 1) * H // args.grid,
                        c * W // args.grid:(c + 1) * W // args.grid]
            if cell.size:
                worst.append((float(cell.mean()), r, c))
    worst.sort()
    for frac, r, c in worst[:3]:
        x = bounds.left + (c + 0.5) * (bounds.right - bounds.left) / args.grid
        y = bounds.top - (r + 0.5) * (bounds.top - bounds.bottom) / args.grid
        flag = "FAIL" if frac < args.max_empty_cell else "ok  "
        print(f"    {flag}  emptiest cell [{r},{c}] near ({x:.0f}, {y:.0f}): "
              f"{100 * frac:.1f} % covered")
    if worst and worst[0][0] < args.max_empty_cell:
        fails.append(f"a {args.grid}x{args.grid} cell is only "
                     f"{100 * worst[0][0]:.1f} % covered -- a source is probably "
                     f"missing over that block (check gdalbuildvrt for "
                     f"'Skipping' warnings)")

    # -- check 3: does NDVI separate green from built-up the way it must? --
    if args.landuse and (args.veg_query or args.built_query):
        import geopandas as gpd
        from rasterio.features import rasterize

        with np.errstate(invalid="ignore", divide="ignore"):
            ndvi = np.where(data, (nir - red) / np.maximum(nir + red, 1e-6), np.nan)
        gdf = gpd.read_file(args.landuse)
        print("\n[3] NDVI by land use")
        stats = {}
        for name, q in (("green", args.veg_query), ("built-up", args.built_query)):
            if not q:
                continue
            sel = gdf.query(q)
            if sel.empty:
                print(f"    --    {name}: query matched nothing, skipped")
                continue
            m = rasterize(((g, 1) for g in sel.geometry), out_shape=data.shape,
                          transform=transform, fill=0, dtype="uint8").astype(bool)
            v = ndvi[m & np.isfinite(ndvi)]
            if v.size < 100:
                print(f"    --    {name}: only {v.size} px, skipped")
                continue
            p = np.percentile(v, [25, 50, 75])
            stats[name] = p
            print(f"    {name:9s} n={v.size:8,d}  p25 {p[0]:+.3f}  "
                  f"median {p[1]:+.3f}  p75 {p[2]:+.3f}")
        if "built-up" in stats and stats["built-up"][1] > BUILT_NDVI_ALARM:
            fails.append(f"built-up median NDVI is {stats['built-up'][1]:+.3f}, "
                         f"above {BUILT_NDVI_ALARM} -- roofs cannot be that green; "
                         f"the NIR band is almost certainly wrong")
            print("    FAIL  built-up is too green to be real")
        if "green" in stats and "built-up" in stats:
            gap = stats["green"][1] - stats["built-up"][1]
            print(f"    separation of medians: {gap:+.3f}"
                  + ("" if gap > 0.05 else "   FAIL  too small to threshold on"))
            if gap <= 0.05:
                fails.append(f"green and built-up medians differ by only {gap:+.3f}; "
                             f"no ndvi_threshold can separate them")

    print()
    if fails:
        for f in fails:
            print(f"FAIL: {f}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
