#!/usr/bin/env python3
"""Relabel a 4-band image's band 4 from Alpha to Undefined, in place.

WHY: a 4-band aerial image is R,G,B,NIR, but many products tag the 4th band
ColorInterp=Alpha. Nothing about the pixels is wrong -- the *label* is -- and
it breaks two things downstream, both quietly:

  * `gdalwarp` without `-nosrcalpha` replaces the NIR values with a synthetic
    opaque 255, so NDVI saturates at +1 and everything classifies as
    vegetation (this produced a 97.6 % vegetated Central Park).
  * `gdalbuildvrt` refuses to mix inputs whose band colour interpretations
    differ -- "expected Undefined, got Alpha" -- and merely *warns* while
    skipping them, so a mosaic silently loses whole sources.

`gdal_edit.py -colorinterp_4 undefined` is the documented fix but is broken in
some environments (a conda `osgeo` against a system libspatialite fails to
import). rasterio can do it in place, header-only, in milliseconds per file.

    python repository/tools/fix_colorinterp.py static_driver/<run>/work/staging/ortho/*.tif
"""
from __future__ import annotations

import sys

import rasterio
from rasterio.enums import ColorInterp

WANT = [ColorInterp.red, ColorInterp.green, ColorInterp.blue, ColorInterp.undefined]


def main(paths: list[str]) -> int:
    # A shell glob normally supplies the paths, so there are no options to
    # parse -- but a bare call, or a habitual --help, must explain itself
    # rather than try to open the flag as a raster.
    if not paths or paths[0] in ("-h", "--help"):
        print(__doc__)
        return 0 if paths else 2
    changed = skipped = failed = 0
    for p in paths:
        try:
            with rasterio.open(p, "r+") as d:
                if d.count != 4:
                    skipped += 1
                    continue
                if list(d.colorinterp) == WANT:
                    skipped += 1
                    continue
                d.colorinterp = WANT
                changed += 1
        except Exception as exc:  # noqa: BLE001
            print(f"  FAILED {p}: {type(exc).__name__}: {exc}", file=sys.stderr)
            failed += 1
    print(f"colorinterp: {changed} relabelled, {skipped} already correct or not "
          f"4-band, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
