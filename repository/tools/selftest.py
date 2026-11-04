#!/usr/bin/env python3
"""Check that palmgeo still produces what it says it produces.

    python repository/tools/selftest.py            # offline, seconds
    python repository/tools/selftest.py --full     # rebuilds the example

The committed example under ``static_driver/example/`` is the reference. It was
built by this pipeline, verified, and its numbers are written down in its
README -- so the useful question is not "did a run finish" but "did it produce
the same thing as last time". Anything else is a test that passes while the
classification quietly changes underneath it.

Two levels, because they answer different questions and cost two orders of
magnitude apart:

**quick** (the default, no network, a few seconds)
    Everything that can be checked against what is already in the repository:
    every shipped configuration loads, the committed drivers open and carry the
    variables they should, each domain's manifest still matches the rasters
    beside it, and the two dependency lists agree. This is what to run after
    editing code, and what a continuous-integration job can run anywhere.

**full** (``--full``, needs the network, a few minutes)
    Rebuilds the example into a scratch directory and compares the result
    against the committed one, variable by variable: cell counts, minima,
    maxima and means. A drift in the classification shows up here as a changed
    mean, which is exactly the failure the quick pass cannot see. The example
    is in Lower Saxony, so this needs the LGLN download servers -- and the
    LGLN updates its data, so a drift can also mean the data moved on since the
    example was committed, not the code.

Exit status is 0 if everything passed and 1 otherwise, so it can be the thing a
CI job or a pre-release check runs.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_CFG = ROOT / "repository" / "config" / "example.yml"
EXAMPLE_DIR = ROOT / "static_driver" / "example"

#: Relative tolerance when comparing a rebuilt driver against the reference.
#: Not zero: GDAL resampling and floating-point summation are not bit-stable
#: across library versions, and a test that fails on the tenth decimal is a
#: test people learn to ignore. Large enough to survive a library upgrade,
#: small enough that a changed threshold or a changed source moves it.
TOLERANCE = 1e-3

failures: List[str] = []
checks = 0


def check(ok: bool, label: str, detail: str = "") -> bool:
    global checks
    checks += 1
    print(f"  {'PASS' if ok else 'FAIL'}  {label}" + (f"   {detail}" if detail else ""))
    if not ok:
        failures.append(f"{label} ({detail})" if detail else label)
    return ok


# ---------------------------------------------------------------------------
# quick: everything answerable from the repository itself
# ---------------------------------------------------------------------------
def check_configs() -> None:
    """Every shipped configuration must load and validate."""
    from palmgeo.config import load, ConfigError

    print("\n=== Shipped configurations ===")
    for path in sorted((ROOT / "repository" / "config").glob("*.yml")):
        if path.name.startswith(("sources_", "mapping_", "licences")):
            continue
        try:
            cfg = load(path)
            check(True, f"{path.name} loads", f"{len(cfg.domains)} domain(s)")
        except ConfigError as exc:
            check(False, f"{path.name} loads", str(exc)[:90])


def check_region_files() -> None:
    """Region files must parse, declare licences, and name no absolute path."""
    import yaml

    print("\n=== Region files ===")
    registry = yaml.safe_load(
        (ROOT / "repository" / "config" / "licences.yml").read_text(encoding="utf-8")
    )["licences"]
    for path in sorted((ROOT / "repository" / "config").glob("sources_*.yml")):
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        sources = data.get("sources") or {}
        bad_lic = []
        absolute = []
        for name, src in sources.items():
            value = src.get("licence")
            ids = [value] if isinstance(value, str) else list(value or [])
            bad_lic += [f"{name}:{i}" for i in ids if str(i) not in registry]
            for key in ("index_url", "url"):
                v = str(src.get(key, ""))
                bare = v[len("file://"):] if v.startswith("file://") else v
                if bare.startswith("/"):
                    absolute.append(f"{name}.{key}")
        check(not bad_lic, f"{path.name}: licence ids are in the registry",
              ", ".join(bad_lic))
        check(not absolute, f"{path.name}: no absolute local paths",
              ", ".join(absolute))


def check_dependency_lists() -> None:
    """pyproject.toml and requirements.txt must not drift apart."""
    import tomllib

    print("\n=== Dependency lists ===")
    proj = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    declared = {d.split(">")[0].split("=")[0].strip().lower()
                for d in proj["project"]["dependencies"]}
    listed = set()
    for line in (ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines():
        line = line.split("#")[0].strip()
        if line:
            listed.add(line.split(">")[0].split("=")[0].strip().lower())
    check(declared == listed, "pyproject.toml and requirements.txt agree",
          f"only in one: {sorted(declared ^ listed)}" if declared != listed else "")


def driver_stats(path: Path) -> Dict[str, Dict[str, float]]:
    """Per-variable summary of a driver, robust enough to compare across runs."""
    import netCDF4 as nc
    import numpy as np

    out: Dict[str, Dict[str, float]] = {}
    with nc.Dataset(path) as ds:
        for name, var in ds.variables.items():
            if var.dtype.kind not in "fiu" or not var.dimensions:
                continue
            data = var[:]
            values = data.compressed() if hasattr(data, "compressed") else np.ravel(data)
            values = values[np.isfinite(values)] if values.dtype.kind == "f" else values
            if values.size == 0:
                out[name] = {"count": 0}
                continue
            out[name] = {
                "count": float(values.size),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
                "mean": float(np.mean(values)),
            }
    return out


def check_committed_example() -> Optional[Dict[str, Dict[str, Dict[str, float]]]]:
    """The reference drivers must open, carry their variables, and match manifests."""
    import rasterio

    print("\n=== Committed example ===")
    drivers = sorted(EXAMPLE_DIR.glob("example_static*"))
    drivers = [d for d in drivers if d.is_file() and d.suffix == ""]
    if not check(len(drivers) == 3, "all three example drivers are present",
                 ", ".join(d.name for d in drivers) or "none found"):
        return None

    reference: Dict[str, Dict[str, Dict[str, float]]] = {}
    for drv in drivers:
        stats = driver_stats(drv)
        reference[drv.name] = stats
        required = {"zt", "buildings_2d", "vegetation_type", "pavement_type",
                    "water_type", "soil_type", "surface_fraction"}
        check(required <= set(stats), f"{drv.name}: required variables present",
              f"missing {sorted(required - set(stats))}" if not required <= set(stats) else "")

    nest = reference.get("example_static_N03", {})
    check("tree_id" in nest and nest["tree_id"]["count"] > 0,
          "the nest carries individual trees",
          f"{int(nest.get('tree_id', {}).get('count', 0))} tree cells")

    # Each manifest is a claim about the rasters beside it; check it still holds.
    for manifest_path in sorted(EXAMPLE_DIR.glob("products/*/manifest.json")):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        domain = manifest["domain"]
        mismatched = []
        for name, claim in manifest.get("products", {}).items():
            raster = manifest_path.parent / f"{name}.tif"
            if not raster.is_file():
                mismatched.append(f"{name} missing")
                continue
            with rasterio.open(raster) as src:
                actual = int(src.read(1, masked=True).count())
            if actual != claim.get("valid_cells"):
                mismatched.append(f"{name} {claim.get('valid_cells')}!={actual}")
        check(not mismatched, f"{domain}: manifest matches the rasters",
              "; ".join(mismatched[:3]))
        check(bool(manifest.get("acquired")), f"{domain}: manifest records acquisition dates",
              ", ".join(f"{k}={v}" for k, v in (manifest.get("acquired") or {}).items()))
    return reference


# ---------------------------------------------------------------------------
# full: rebuild and compare
# ---------------------------------------------------------------------------
def _run(args: List[str], cwd: Path) -> Tuple[int, str]:
    env = {**dict(__import__("os").environ), "PYTHONPATH": str(ROOT / "repository")}
    proc = subprocess.run(args, cwd=cwd, capture_output=True, text=True, env=env)
    return proc.returncode, proc.stdout + proc.stderr


def check_rebuild(reference: Dict[str, Dict[str, Dict[str, float]]]) -> None:
    """Rebuild the example into a scratch tree and compare against the reference."""
    import yaml

    print("\n=== Rebuild and compare ===")
    scratch = Path(tempfile.mkdtemp(prefix="palmgeo-selftest-"))
    cfg_dir = ROOT / "static_driver" / "_selftest" / "config"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = cfg_dir / "_selftest.yml"
    raw = yaml.safe_load(EXAMPLE_CFG.read_text(encoding="utf-8"))
    raw["project"]["name"] = "selftest_static"
    raw["project"]["outdir"] = "static_driver/_selftest"
    raw["project"]["workdir"] = str(scratch)
    cfg_path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    print(f"  scratch: {scratch}")
    print("  approving the example's sources (see static_driver/example/docs/DATA-SOURCES.md)")

    try:
        for args in (
            [sys.executable, "-m", "palmgeo.cli", "sources", str(cfg_path), "--approve"],
            [sys.executable, "-m", "palmgeo.cli", "all", str(cfg_path)],
        ):
            code, out = _run(args, ROOT)
            if not check(code == 0, f"`{args[3]}` succeeded", out.strip()[-160:] if code else ""):
                return

        csd = ROOT / "static_driver" / "_selftest" / "config" / "selftest_static_csd.yml"
        code, out = _run(
            ["env", "-u", "PROJ_DATA", "-u", "PROJ_LIB", "MPLBACKEND=Agg",
             "palm_csd", str(csd)], ROOT)
        if not check(code == 0, "palm_csd succeeded", out.strip()[-160:] if code else ""):
            return
        code, out = _run(
            [sys.executable, str(ROOT / "repository" / "tools" / "finalize.py"),
             str(cfg_path)], ROOT)
        code, out = _run(
            [sys.executable, str(ROOT / "repository" / "tools" / "verify_driver.py"),
             str(cfg_path)], ROOT)
        check("ALL CHECKS PASSED" in out, "verify_driver passes on the rebuild",
              out.strip()[-160:] if "ALL CHECKS PASSED" not in out else "")

        for ref_name, ref in reference.items():
            built = ROOT / "static_driver" / "_selftest" / ref_name.replace(
                "example_static", "selftest_static")
            if not built.is_file():
                check(False, f"{ref_name}: rebuilt driver exists", str(built))
                continue
            now = driver_stats(built)
            drifted = []
            for var, want in ref.items():
                got = now.get(var)
                if got is None:
                    drifted.append(f"{var} absent")
                    continue
                for key, value in want.items():
                    other = got.get(key)
                    if other is None:
                        drifted.append(f"{var}.{key} absent")
                    elif abs(other - value) > TOLERANCE * max(1.0, abs(value)):
                        drifted.append(f"{var}.{key} {value:.4g}->{other:.4g}")
            check(not drifted, f"{ref_name}: matches the committed reference",
                  "; ".join(drifted[:4]))
    finally:
        import shutil
        shutil.rmtree(scratch, ignore_errors=True)
        shutil.rmtree(ROOT / "static_driver" / "_selftest", ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--full", action="store_true",
                    help="also rebuild the example and compare (needs the network)")
    args = ap.parse_args()

    check_configs()
    check_region_files()
    check_dependency_lists()
    reference = check_committed_example()

    if args.full:
        if reference is None:
            check(False, "rebuild comparison", "no reference to compare against")
        else:
            check_rebuild(reference)

    print()
    if failures:
        print(f"{len(failures)} of {checks} FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"ALL {checks} CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
