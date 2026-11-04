#!/usr/bin/env python3
"""Acceptance check for a finished PALM static driver set.

    python tools/verify_driver.py config/<run>.yml

Checks the properties PALM cannot recover from once the driver is written --
point counts, cell centring, nest alignment, a shared vertical datum, exactly
one surface type per cell -- and then quantifies how much a nest disagrees with
its parent along the shared boundary, which is the one requirement that cannot
be met exactly by construction.

Everything is derived from the run configuration, so this verifies the driver
that was actually asked for rather than a hard-coded expectation.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import netCDF4 as nc  # noqa: E402

from palmgeo import config as pgconfig  # noqa: E402
from palmgeo.grid import multigrid_leftover  # noqa: E402

#: Written by palm_csd for every domain; a missing one is a broken driver.
REQUIRED = [
    "zt", "buildings_2d", "building_id", "building_type", "vegetation_type",
    "pavement_type", "water_type", "soil_type", "surface_fraction",
]
#: Coordinate and container variables that legitimately carry no _FillValue.
NO_FILL_OK = {
    "x", "y", "z", "zlad", "zsoil", "crs", "nsurface_fraction", "nvegetation_pars",
    "nalbedo_pars", "nbuilding_pars", "nsoil_pars", "nwater_pars", "npavement_pars",
    "nvegetation_type", "lat", "lon", "E_UTM", "N_UTM",
}

fails: list[str] = []
warns: list[str] = []


def check(ok: bool, label: str, detail: str = "", warn_only: bool = False) -> None:
    """Report one check, and remember the failures with their reason.

    The label names the property being tested, so it reads correctly after
    both PASS and FAIL. The closing summary carries the DETAIL as well, which
    it did not: a warning recorded as "2000 points on x are 2-smooth" was
    printed at the end with no tag and no reason, where it read as a statement
    that the property held -- the exact opposite of what it meant.
    """
    tag = "PASS" if ok else ("WARN" if warn_only else "FAIL")
    print(f"  {tag}  {label}" + (f"   {detail}" if detail else ""))
    if not ok:
        (warns if warn_only else fails).append(f"{label} ({detail})" if detail else label)


def driver_path(cfg, domain) -> Path:
    """palm_csd writes <file_out>_<domain>; the root is renamed for PALM."""
    base = str(cfg.get("project.name", "palm"))
    if domain.domain_parent is None:
        root = cfg.outdir / base
        return root if root.is_file() else cfg.outdir / f"{base}_{domain.name}"
    return cfg.outdir / f"{base}_{domain.name}"


def main(argv: list[str]) -> int:
    # Deliberately not argparse: this takes exactly one positional argument and
    # every other spelling is a mistake. What it must NOT do is pass "--help"
    # to the configuration loader, which then fails with a parse error about a
    # file the user never meant to name.
    if len(argv) != 2 or argv[1] in ("-h", "--help"):
        print(__doc__)
        return 0 if len(argv) == 2 else 2
    cfg = pgconfig.load(argv[1])
    by_name = {d.name: d for d in cfg.domains}
    print(f"Run:    {cfg.path}")
    print(f"Driver: {cfg.outdir}\n")

    print("=== Grid geometry ===")
    ds = {}
    for d in cfg.domains:
        path = driver_path(cfg, d)
        if not path.is_file():
            check(False, f"{d.name}: file present", str(path))
            continue
        f = nc.Dataset(path, "r")
        ds[d.name] = f
        x, y = np.asarray(f.variables["x"][:]), np.asarray(f.variables["y"][:])

        check(len(x) == d.width and len(y) == d.height, f"{d.name}: point count",
              f"x {len(x)} (want {d.width}), y {len(y)} (want {d.height})")
        for axis, n in (("x", len(x)), ("y", len(y))):
            rest = multigrid_leftover(n)
            check(rest == 1,
                  f"{d.name}: {axis} point count coarsens for the multigrid solver",
                  f"{n} leaves the factor {rest}" if rest != 1
                  else f"{n} coarsens all the way",
                  warn_only=True)
        for axis, v in (("x", x), ("y", y)):
            want = (np.arange(len(v)) + 0.5) * d.pixel_size
            check(np.allclose(v, want, atol=1e-6), f"{d.name}: {axis} cell-centred",
                  f"max|delta| {np.max(np.abs(v - want)):.3g} m")
        check(abs(f.getncattr("origin_x") - d.origin_x) < 1e-6
              and abs(f.getncattr("origin_y") - d.origin_y) < 1e-6, f"{d.name}: origin_x/y",
              f"({f.getncattr('origin_x'):.1f}, {f.getncattr('origin_y'):.1f})")
        want_rot = float(cfg.get("csd.settings.rotation_angle", 0.0))
        check(abs(float(f.getncattr("rotation_angle")) - want_rot) < 1e-9,
              f"{d.name}: rotation_angle = {want_rot:g}")

    if not ds:
        print("\nNo driver files found -- run palm_csd first.")
        return 1

    print("\n=== Vertical datum ===")
    zs = {n: float(f.getncattr("origin_z")) for n, f in ds.items()}
    check(len({round(v, 9) for v in zs.values()}) == 1,
          "origin_z identical in every file", f"{sorted(set(zs.values()))}")

    print("\n=== Nesting ===")
    for d in cfg.domains:
        if d.domain_parent is None:
            continue
        p = by_name[d.domain_parent]
        offx, offy = d.origin_x - p.origin_x, d.origin_y - p.origin_y
        check(offx % p.pixel_size == 0 and offy % p.pixel_size == 0,
              f"{d.name} in {p.name}: offset is a multiple of {p.pixel_size:g} m",
              f"{offx:g}/{p.pixel_size:g}={offx / p.pixel_size:g}, "
              f"{offy:g}/{p.pixel_size:g}={offy / p.pixel_size:g}")
        px0, py0, px1, py1 = p.bounds
        cx0, cy0, cx1, cy1 = d.bounds
        check(cx0 >= px0 and cy0 >= py0 and cx1 <= px1 and cy1 <= py1,
              f"{d.name} fully inside {p.name}")
        ratio = p.pixel_size / d.pixel_size
        check(abs(ratio - round(ratio)) < 1e-9 and round(ratio) > 1,
              f"{d.name}/{p.name}: dx ratio is an integer > 1", f"{ratio:g}")

    print("\n=== Variables and fill values ===")
    for name, f in ds.items():
        missing = [v for v in REQUIRED if v not in f.variables]
        check(not missing, f"{name}: required variables present",
              f"missing: {missing}" if missing else "")
        bad = [v for v in f.variables
               if v not in NO_FILL_OK and "_FillValue" not in f.variables[v].ncattrs()]
        check(not bad, f"{name}: _FillValue set", f"without: {bad}" if bad else "")

    print("\n=== Coverage ===")
    for name, f in ds.items():
        zt = np.ma.filled(f.variables["zt"][:].astype(float), np.nan)
        cover = 100 * np.isfinite(zt).mean()
        check(cover > 99.999, f"{name}: zt covers the whole domain", f"{cover:.4f} %")
        g = np.zeros(zt.shape, int)
        for v in ("vegetation_type", "pavement_type", "water_type", "buildings_2d"):
            if v in f.variables:
                g += np.isfinite(np.ma.filled(f.variables[v][:].astype(float),
                                              np.nan)).astype(int)
        check(bool((g == 1).all()), f"{name}: exactly one surface type per cell",
              f"{100 * (g == 1).mean():.4f} %")

    print("\n=== Child/parent boundary ===")
    print("  Measured against the parent terrain rasterised onto the PARENT's dz.")
    print("  With interpolate_terrain on, palm_csd builds the child boundary on the")
    print("  parent's vertical grid (create_driver.py: zt_ip = height_to_z_grid(zt_ip,")
    print("  parent.dz)), because PALM discretises the parent terrain that way at")
    print("  runtime; measuring against the parent file's raw zt is then misleading")
    print("  while that file carries its terrain unquantised (use_palm_z_axis: false).")
    print("  With it off, a difference of up to one parent grid spacing is the expected")
    print("  result, not a defect -- the check reports it and does not fail.")
    for d in cfg.domains:
        if d.domain_parent is None or d.name not in ds or d.domain_parent not in ds:
            continue
        p = by_name[d.domain_parent]
        r = int(round(p.pixel_size / d.pixel_size))
        ct = np.ma.filled(ds[d.name].variables["zt"][:].astype(float), np.nan)
        pt = np.ma.filled(ds[p.name].variables["zt"][:].astype(float), np.nan)
        j0 = int((d.origin_y - p.origin_y) / p.pixel_size)
        i0 = int((d.origin_x - p.origin_x) / p.pixel_size)
        nyp, nxp = d.height // r, d.width // r
        sub = pt[j0:j0 + nyp, i0:i0 + nxp]
        agg = ct[:nyp * r, :nxp * r].reshape(nyp, r, nxp, r).mean(axis=(1, 3))
        delta = np.abs(agg - np.round(sub / p.dz) * p.dz)
        ring = np.zeros(sub.shape, bool)
        ring[0, :] = ring[-1, :] = ring[:, 0] = ring[:, -1] = True
        blended = bool(d.options.get("interpolate_terrain", False))
        detail = (f"median {np.nanmedian(delta[ring]):.2f} m, "
                  f"p95 {np.nanpercentile(delta[ring], 95):.2f} m, "
                  f"max {np.nanmax(delta[ring]):.2f} m"
                  + (f"   (dz {d.dz:g} vs {p.dz:g} m)" if d.dz != p.dz else ""))
        if blended:
            # Tolerance is one whole parent level, not half of one. Child and
            # parent carry their terrain on different vertical ladders whenever
            # their dz differs, and a true elevation of 16 m lands on 16 in a
            # dz=16 child and on 0 or 32 in a dz=32 parent. That offset is
            # arithmetic, not a blending failure, and half a level would flag it.
            check(float(np.nanpercentile(delta[ring], 95)) <= p.dz,
                  f"{d.name} vs {p.name}: boundary on the parent z grid", detail)
        else:
            # Not a criterion: interpolate_terrain is off, so the nest is meant to
            # carry its own terrain to the edge. Report the size of the step so it
            # is a decision rather than a surprise.
            print(f"  INFO  {d.name} vs {p.name}: no blending "
                  f"(interpolate_terrain false), step at the boundary   {detail}")

    for f in ds.values():
        f.close()
    print()
    if warns:
        print(f"{len(warns)} warning(s): {', '.join(warns)}")
    print("ALL CHECKS PASSED" if not fails else f"{len(fails)} FAILED: {', '.join(fails)}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
