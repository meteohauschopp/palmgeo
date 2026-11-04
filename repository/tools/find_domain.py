#!/usr/bin/env python3
"""Find the largest PALM-legal domain that stays inside the data coverage.

WHY THIS EXISTS
---------------
Picking a domain box by eye fails in two directions at once, and both cost a
full rebuild to discover:

* Too generous, and `domains` aborts on missing terrain -- a survey stops at
  some administrative border you did not think about. NYC's LiDAR ends at the
  city line, which runs down the middle of the Hudson; a box drawn around
  Manhattan had 38 % of its terrain missing.
* Too cautious, and you silently cut away the thing you meant to simulate.
  The first Manhattan box lost all of Upper Manhattan, and its west edge cut
  into the island because the search range had been capped without checking.

This searches the actual coverage mask instead. For every 2-smooth point count
(PALM's multigrid solver wants counts that factor into 2s and 3s only -- see
the palmgeo README section 3.2) and every candidate origin on a coarse grid, it
computes the share of missing terrain and the share of a target polygon set
retained, and reports the boxes that satisfy the coverage limit.

    python repository/tools/find_domain.py \\
        --terrain static_driver/<run>/work/staging/dtm/<dtm>.tif \\
        --target  static_driver/<run>/work/vector/cadastre_landuse.gpkg \\
        --target-drop-prefix LandUse=WTR \\
        --dx 16 --max-missing 0.02

`--target` is what you are trying to keep inside the box (here: the island's
land parcels). Everything else is measured against the terrain raster's own
nodata.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize


def two_smooth(n: int) -> bool:
    """PALM's multigrid solver coarsens by halving/thirding; a count carrying
    any other prime factor stops the coarsening dead (README section 3.2)."""
    m = n
    for p in (2, 3):
        while m % p == 0:
            m //= p
    return m == 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--terrain", required=True, help="staged terrain raster (nodata = gaps)")
    ap.add_argument("--target", help="vector file whose area the box should retain")
    ap.add_argument("--target-query", help="pandas .query() expression selecting target rows")
    ap.add_argument("--target-drop-prefix", metavar="COL=PREFIX",
                    help="drop target rows whose COL starts with PREFIX, e.g. "
                         "LandUse=WTR to exclude water from 'land to retain'")
    ap.add_argument("--dx", type=float, required=True, help="grid spacing [m]")
    ap.add_argument("--max-missing", type=float, default=0.02,
                    help="palmgeo's `domains` stage aborts above 2 %% (default 0.02)")
    ap.add_argument("--counts", type=int, nargs="+",
                    default=[384, 512, 576, 768, 864, 1024, 1152, 1296],
                    help="candidate point counts (2-smooth ones are kept)")
    ap.add_argument("--step", type=float, default=128.0, help="origin search step [m]")
    ap.add_argument("--downsample", type=int, default=8, help="analysis downsampling")
    ap.add_argument("--top", type=int, default=10)
    args = ap.parse_args()

    counts = [c for c in args.counts if two_smooth(c)]
    dropped = sorted(set(args.counts) - set(counts))
    if dropped:
        print(f"skipping non-2-smooth counts: {dropped}")

    with rasterio.open(args.terrain) as src:
        sc = args.downsample
        H, W = src.height // sc, src.width // sc
        a = src.read(1, out_shape=(H, W), resampling=Resampling.nearest)
        tr = src.transform * src.transform.scale(src.width / W, src.height / H)
        nod = (a == src.nodata) | ~np.isfinite(a) | (a < -1e30)
        L, B, R, T = src.bounds.left, src.bounds.bottom, src.bounds.right, src.bounds.top

        if args.target:
            gdf = gpd.read_file(args.target)
            if args.target_query:
                gdf = gdf.query(args.target_query)
            if args.target_drop_prefix:
                col, prefix = args.target_drop_prefix.split("=", 1)
                gdf = gdf[~gdf[col].astype(str).str.startswith(prefix)]
            if gdf.empty:
                print("target selection is empty -- check --target-query / "
                      "--target-drop-prefix")
                return 2
            target = rasterize(((g, 1) for g in gdf.geometry), out_shape=(H, W),
                               transform=tr, fill=0, dtype="uint8").astype(bool)
        else:
            target = np.zeros((H, W), bool)

    px, py = (R - L) / W, (T - B) / H
    total = max(int(target.sum()), 1)
    print(f"terrain {W}x{H} @ ~{px:.1f} m, missing overall {nod.mean():.1%}, "
          f"target cells {total:,}")

    ii_n = np.cumsum(np.cumsum(nod.astype(np.int64), 0), 1)
    ii_t = np.cumsum(np.cumsum(target.astype(np.int64), 0), 1)

    def box_sum(ii, r0, r1, c0, c1):
        s = ii[r1 - 1, c1 - 1]
        if r0 > 0:
            s -= ii[r0 - 1, c1 - 1]
        if c0 > 0:
            s -= ii[r1 - 1, c0 - 1]
        if r0 > 0 and c0 > 0:
            s += ii[r0 - 1, c0 - 1]
        return s

    results = []
    for nxp in counts:
        for nyp in counts:
            wm, hm = nxp * args.dx, nyp * args.dx
            if wm > (R - L) or hm > (T - B):
                continue
            for ox in np.arange(L, R - wm + 1, args.step):
                for oy in np.arange(B, T - hm + 1, args.step):
                    c0, c1 = int((ox - L) / px), int(np.ceil((ox + wm - L) / px))
                    r0, r1 = H - int(np.ceil((oy + hm - B) / py)), H - int((oy - B) / py)
                    if c1 <= c0 or r1 <= r0 or r0 < 0:
                        continue
                    miss = box_sum(ii_n, r0, r1, c0, c1) / ((r1 - r0) * (c1 - c0))
                    if miss <= args.max_missing:
                        kept = box_sum(ii_t, r0, r1, c0, c1) / total
                        results.append((kept, miss, ox, oy, nxp, nyp))

    if not results:
        print(f"\nNo box under {args.max_missing:.0%} missing terrain. Either the "
              f"coverage genuinely cannot hold one, or --max-missing is too tight.")
        return 1

    results.sort(reverse=True)
    print(f"\n{'target':>7} {'missing':>8} {'origin_x':>10} {'origin_y':>11} "
          f"{'nx':>5} {'ny':>5} {'extent km':>14}   (nx/ny below are POINT counts; "
          f"PALM's nx = points-1)")
    for kept, miss, ox, oy, nxp, nyp in results[:args.top]:
        # snap the origin to the grid spacing -- cleaner, and costs nothing
        sox, soy = round(ox / args.dx) * args.dx, round(oy / args.dx) * args.dx
        print(f"{kept:7.1%} {miss:8.2%} {sox:10.0f} {soy:11.0f} {nxp:5d} {nyp:5d} "
              f"{nxp*args.dx/1000:6.2f}x{nyp*args.dx/1000:6.2f}")
    best = results[0]
    print(f"\nbest: nx: {best[4]-1}, ny: {best[5]-1}, "
          f"origin_x: {round(best[2]/args.dx)*args.dx:.1f}, "
          f"origin_y: {round(best[3]/args.dx)*args.dx:.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
