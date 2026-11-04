#!/usr/bin/env python3
"""Is any PART of the driver systematically empty?

WHY THIS EXISTS
---------------
Aggregate numbers hide clipped sources. A run with "56 478 New Jersey
buildings" sounded complete; binned by row, building cover was 26 % up to
y = 4 517 000 and 0-5 % north of y = 4 519 000. That step is the signature of a
fetch clipped to a rectangle while the feature it was meant to follow -- a
river, a state line, a coast -- runs diagonally. No total, and no check that
judges a class as a whole, can show it. A user opening the driver sees it
immediately, which is the wrong way to find out.

This bins the domain and reports presence per block, for buildings by default.
Blocks are compared against the domain's own land median, so a genuinely
low-density region is not a failure -- what it looks for is a *cliff*.

    python repository/tools/check_coverage.py static_driver/<run>/config/<run>.yml
    python repository/tools/check_coverage.py static_driver/<run>/config/<run>.yml --var vegetation_type

Exit code 1 if a block looks systematically empty.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--var", default="buildings_2d",
                    help="driver variable to test for presence (default buildings_2d)")
    ap.add_argument("--blocks", type=int, default=16,
                    help="bin the domain into BLOCKS x BLOCKS (default 16)")
    ap.add_argument("--min-land", type=float, default=0.25,
                    help="ignore blocks that are less than this share dry land "
                         "(default 0.25) -- an empty block of open water is not "
                         "a finding")
    ap.add_argument("--green", type=float, default=0.60,
                    help="a block whose land is at least this share vegetation "
                         "is treated as open space, not as a missing source "
                         "(default 0.60)")
    ap.add_argument("--cliff", type=float, default=0.15,
                    help="flag a land block whose cover is below this fraction "
                         "of the domain's land median (default 0.15)")
    args = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import netCDF4
    from palmgeo import config as C

    cfg = C.load(args.config)
    base = str(cfg.get("project.name", "palm"))
    findings = 0

    for domain in cfg.domains:
        # Same resolution as verify_driver.py: the root domain keeps the bare
        # project name once finalize.py has run, children are suffixed.
        driver = cfg.outdir / base
        if domain.name != "root" or not driver.is_file():
            driver = cfg.outdir / f"{base}_{domain.name}"
        if not driver.exists():
            print(f"{domain.name}: no driver at {driver}, skipped")
            continue
        ds = netCDF4.Dataset(driver)
        var = ds.variables[args.var][:]
        water = ds.variables["water_type"][:]
        veg = ds.variables["vegetation_type"][:]
        ox, oy, dx = ds.origin_x, ds.origin_y, domain.pixel_size
        ny, nx = var.shape
        have = ~np.ma.getmaskarray(var)
        # water_type is masked exactly where there is no water, so its
        # mask array IS the dry-land mask.
        land = np.ma.getmaskarray(water)
        green = ~np.ma.getmaskarray(veg)
        ds.close()

        B = args.blocks
        cover = np.full((B, B), np.nan)
        landsh = np.zeros((B, B))
        greensh = np.zeros((B, B))
        for r in range(B):
            for c in range(B):
                sl = (slice(r * ny // B, (r + 1) * ny // B),
                      slice(c * nx // B, (c + 1) * nx // B))
                lm = land[sl]
                landsh[r, c] = lm.mean()
                if lm.sum() > 0:
                    cover[r, c] = have[sl][lm].mean()
                    greensh[r, c] = green[sl][lm].mean()

        usable = landsh >= args.min_land
        if not usable.any():
            print(f"{domain.name}: no block is {args.min_land:.0%} land; nothing to test")
            continue
        med = float(np.nanmedian(cover[usable]))
        print(f"\n=== {domain.name}: '{args.var}' cover over land, "
              f"{B}x{B} blocks ({nx*dx/B/1000:.2f} x {ny*dx/B/1000:.2f} km each) ===")
        print(f"median over land blocks: {med:.1%}   "
              f"(a block under {args.cliff:.0%} of that is flagged)\n")
        print("  north at top; '.' = under 5 %, digits = cover in tens of %, "
              "'~' = mostly water")
        for r in range(B - 1, -1, -1):
            line = ""
            for c in range(B):
                if not usable[r, c]:
                    line += "~"
                elif np.isnan(cover[r, c]) or cover[r, c] < 0.05:
                    line += "."
                else:
                    line += str(min(9, int(cover[r, c] * 10)))
            print("   " + line)

        # A park has no buildings and that is not a defect. Blocks whose land
        # is mostly vegetation are named as such instead of failed -- otherwise
        # Central Park is reported every run and the reader learns to skip the
        # findings, which is worse than not checking.
        bad = usable & (cover < args.cliff * med)
        parky = greensh >= args.green
        if bad.any():
            for r, c in zip(*np.where(bad)):
                x = ox + (c + 0.5) * nx * dx / B
                y = oy + (r + 0.5) * ny * dx / B
                if parky[r, c]:
                    print(f"    open space near ({x:.0f}, {y:.0f}): "
                          f"{cover[r, c]:.1%} buildings on land that is "
                          f"{greensh[r, c]:.0%} vegetation -- expected, not a finding")
                    continue
                print(f"    FAIL block near ({x:.0f}, {y:.0f}): "
                      f"{cover[r, c]:.1%} cover on {landsh[r, c]:.0%} land, "
                      f"only {greensh[r, c]:.0%} of it vegetation")
                findings += 1

    print()
    if findings:
        print(f"{findings} block(s) look systematically empty -- a source is "
              f"probably clipped short of the domain")
        return 1
    print("no systematically empty region found")
    return 0


if __name__ == "__main__":
    sys.exit(main())
