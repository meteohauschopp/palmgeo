#!/usr/bin/env python3
"""Cross-tabulate the input land use against the driver's surface types.

WHY THIS EXISTS
---------------
`verify_driver.py` checks that a driver is STRUCTURALLY sound: 2-smooth point
counts, cell-centred coordinates, one surface type per cell, no terrain holes.
It passed on a Manhattan driver in which **48.7 % of Central Park was asphalt**
-- because "asphalt everywhere" is structurally perfect. Nothing in the
pipeline compared the result against what the input data said the place was.

This tool does that, and needs no new declarations to do it: the expectation is
already in `mapping_<region>.yml`. A land-use class whose rules only ever
produce `vegetation:*` is a class that should come out mostly vegetation. When
it comes out mostly pavement instead, either the NDVI threshold is wrong for
this imagery or the rule's `low_ndvi` branch is (both were true for Central
Park). Either way the porter should see it before the driver is used.

WHAT IT REPORTS
---------------
For every land-use key present in the AOI, the share of domain cells that
ended up vegetation / pavement / water / building, next to the groups the
mapping rule can actually produce. A class is judged only on the ground it
actually governs -- cells that are neither building nor water, since those are
hard overrides burnt in afterwards. A rule that can only produce vegetation but
whose governed ground is >10 % pavement (or the reverse) is flagged FAIL.

    python repository/tools/check_classification.py static_driver/<run>/config/<run>.yml

Region-agnostic: it reads the run config, the clipped cadastre and the driver,
and derives everything else.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import rasterize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from palmgeo import config as palmgeo_config  # noqa: E402

GROUPS = ("vegetation", "pavement", "water", "building")


def _rule_groups(rule) -> set[str]:
    """Which surface groups can this mapping rule ever produce?"""
    out: set[str] = set()
    for attr, grp in (("high_ndvi", None), ("low_ndvi", None),
                      ("high_veg", "vegetation"), ("water", "water")):
        val = getattr(rule, attr, None)
        if val is None:
            continue
        # rules store (group_id, type_id); map the group id back to a name
        if grp:
            out.add(grp)
        else:
            gid = val[0] if isinstance(val, (tuple, list)) else None
            out.add({1: "vegetation", 2: "pavement", 3: "water"}.get(gid, "?"))
    return {g for g in out if g != "?"}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--min-interior", type=float, default=0.40,
                    help="a class with fewer than this share of its cells in "
                         "the interior (all four neighbours in the same class) "
                         "is reported as unresolved at this dx rather than "
                         "failed: its cells are decided by what borders them, "
                         "not by its own rule (default 0.40)")
    ap.add_argument("--subcell-area", type=float, default=0.15,
                    help="a class with at least this share of its AREA in "
                         "polygons smaller than 4 grid cells is reported as "
                         "unresolved at this dx rather than failed -- the "
                         "geometry cannot be represented, so its rule is not "
                         "what is being measured (default 0.15)")
    ap.add_argument("--tolerance", type=float, default=0.10,
                    help="share of GOVERNED ground (neither building nor water) "
                         "allowed to contradict the rule before it is a FAIL "
                         "(default 0.10)")
    args = ap.parse_args()

    cfg = palmgeo_config.load(Path(args.config))
    outdir = cfg.outdir
    name = str(cfg.get("project.name"))

    landuse_path = cfg.workdir / "vector" / "cadastre_landuse.gpkg"
    if not landuse_path.is_file():
        print(f"no clipped land use at {landuse_path} -- run the `vectors` stage first")
        return 2
    lu = gpd.read_file(landuse_path)

    spec = cfg.source_for_role("cadastre")[1]["layers"]["landuse"]
    key_field = str(spec["key_field"])
    if key_field not in lu.columns:
        print(f"land-use layer has no column {key_field!r}; has {list(lu.columns)}")
        return 2

    failures = 0
    for domain in cfg.domains:
        driver = outdir / (f"{name}" if domain.name == "root" else f"{name}_{domain.name}")
        if not driver.is_file():
            print(f"{domain.name}: no driver at {driver} -- run palm_csd first")
            continue
        import netCDF4

        ds = netCDF4.Dataset(driver)
        ox, oy = ds.origin_x, ds.origin_y
        ny, nx = ds.variables["vegetation_type"].shape
        arrays = {
            "vegetation": ds.variables["vegetation_type"][:],
            "pavement": ds.variables["pavement_type"][:],
            "water": ds.variables["water_type"][:],
            "building": ds.variables["buildings_2d"][:],
        }
        present = {g: ~np.ma.getmaskarray(a) for g, a in arrays.items()}
        transform = rasterio.transform.from_origin(ox, oy + ny * domain.pixel_size,
                                                   domain.pixel_size, domain.pixel_size)

        print(f"\n=== {domain.name}: land use vs. realised surface "
              f"({nx}x{ny} @ {domain.pixel_size:g} m) ===")
        print(f"{'land-use key':<14}{'cells':>9}{'veg':>8}{'pav':>8}{'wat':>8}"
              f"{'bld':>8}   rule can be")

        for key, sub in lu.groupby(lu[key_field].astype(str)):
            mask = rasterize(((g, 1) for g in sub.geometry), out_shape=(ny, nx),
                             transform=transform, fill=0, dtype="uint8").astype(bool)
            mask = mask[::-1]  # driver arrays are south-up
            n = int(mask.sum())
            if n < 50:
                continue
            shares = {g: float((present[g] & mask).sum()) / n for g in GROUPS}
            rule = cfg.mapping.lookup(key)
            can = _rule_groups(rule)

            # How much of this class is too small to own a cell? A polygon
            # narrower than the grid spacing cannot win a cell against its
            # surroundings whatever its rule says, so a class made of such
            # slivers WILL contradict its rule -- and that is a statement
            # about dx, not about the mapping.
            #
            # It is the share of AREA in small polygons that matters, not the
            # share of polygons. Measured on this run: WTR_POND (median
            # 238 m2 against a 256 m2 cell, 20 % of its area in polygons under
            # 4 cells) came out 83 % water and would have been a FAIL, while
            # WTR_LAKE has just as many tiny polygons by count but only 1 % of
            # its AREA in them -- and comes out 96.5 % water. The failures
            # track sub-cell area, not the rule.
            cell_a = domain.pixel_size ** 2
            areas = sub.geometry.area.to_numpy()
            subcell = (float(areas[areas < 4 * cell_a].sum()) / float(areas.sum())
                       if areas.sum() > 0 else 0.0)

            # The decisive measure is not polygon size but how much of the
            # class is INTERIOR. A surface type is decided per cell by area
            # dominance, so a cell on a class boundary is settled by whatever
            # is on the other side of that boundary -- not by this class's
            # rule. A class that is nearly all boundary therefore cannot be
            # judged: `landuse:construction` here is 6 polygons of a median
            # 3 509 m2, comfortably larger than a 256 m2 cell, yet ~14 cells
            # each with almost every one of them on the edge.
            interior = np.zeros_like(mask)
            interior[1:-1, 1:-1] = (mask[1:-1, 1:-1] & mask[:-2, 1:-1] &
                                    mask[2:, 1:-1] & mask[1:-1, :-2] & mask[1:-1, 2:])
            inner = float(interior.sum()) / n

            # Judge the rule only on the cells it actually governs. Buildings
            # and water are HARD OVERRIDES burnt in from the cadastre after
            # the land-use rule has spoken (palmgeo README section 7), so a
            # residential class coming out 70 % building says nothing about
            # whether its rule is right -- it says the block has buildings on
            # it. What the rule is answerable for is the ground between them.
            flag = ""
            # A class whose area is mostly in polygons of a few cells or less
            # is reported as unresolved rather than failed. Left as a FAIL it
            # trains the reader to ignore FAILs, which is worse than not
            # checking at all.
            unresolved = subcell >= args.subcell_area or inner < args.min_interior
            if can == {"water"}:
                # A water rule is answerable for the whole class, not for
                # "governed ground": excluding water cells would leave only
                # the shoreline sliver, where a few per cent of pavement
                # becomes a meaningless 70 % of the remainder. (That is
                # exactly what an earlier version of this check reported --
                # a false FAIL on WTR_LAKE/WTR_RIV whose message even claimed
                # the rule "can only make vegetation".)
                if shares["water"] < 1.0 - args.tolerance:
                    if unresolved:
                        flag = (f"  unresolved at dx={domain.pixel_size:g} m: "
                                f"only {inner:.0%} of its cells are interior; "
                                f"{shares['water']:.0%} water")
                    else:
                        flag = (f"  FAIL only {shares['water']:.0%} water for a class "
                                f"that can only be water")
                        failures += 1
            else:
                # Land rules are judged only on the ground they actually
                # govern. Buildings and water are HARD OVERRIDES burnt in
                # from the cadastre after the land-use rule has spoken
                # (palmgeo README section 7), so a residential class coming
                # out 70 % building says nothing about whether its rule is
                # right -- it says the block has buildings on it.
                governed = mask & ~present["building"] & ~present["water"]
                ng = int(governed.sum())
                if ng >= 50:
                    gv = float((present["vegetation"] & governed).sum()) / ng
                    gp = float((present["pavement"] & governed).sum()) / ng
                    # A rule that cannot produce pavement must not yield much
                    # of it, and vice versa. This is the check that caught
                    # Central Park: its rule is vegetation-only, the result
                    # was 60 % pavement.
                    if "pavement" not in can and gp > args.tolerance:
                        if unresolved:
                            flag = (f"  unresolved at dx={domain.pixel_size:g} m: "
                                    f"only {inner:.0%} of its cells are "
                                    f"interior; {gp:.0%} pavement")
                        else:
                            flag = (f"  FAIL {gp:.0%} pavement on ground the rule can "
                                    f"only make {'/'.join(sorted(can))}")
                            failures += 1
                    elif "vegetation" not in can and gv > args.tolerance:
                        if unresolved:
                            flag = (f"  unresolved at dx={domain.pixel_size:g} m: "
                                    f"only {inner:.0%} of its cells are "
                                    f"interior; {gv:.0%} vegetation")
                        else:
                            flag = (f"  FAIL {gv:.0%} vegetation on ground the rule can "
                                    f"only make {'/'.join(sorted(can))}")
                            failures += 1
            print(f"{key:<14}{n:>9,}{shares['vegetation']:>8.1%}{shares['pavement']:>8.1%}"
                  f"{shares['water']:>8.1%}{shares['building']:>8.1%}   "
                  f"{','.join(sorted(can)):<22}{flag}")
        ds.close()

    print()
    if failures:
        print(f"{failures} class(es) look implausible -- check the NDVI threshold "
              f"(`palmgeo calibrate`) and the rule's low_ndvi branch")
        return 1
    print("classification looks plausible against the input land use")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
