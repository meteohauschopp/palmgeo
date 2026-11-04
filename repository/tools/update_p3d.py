#!/usr/bin/env python3
"""Bring an existing set of PALM p3d namelists in line with a palmgeo grid.

    python tools/update_p3d.py config/<run>.yml --layout config/<run>_layout.yml \
           --src <dir with the old namelists> [--dst <dir>] [--dry-run]

Only geometry is touched: ``nx``/``ny``/``dx``/``dy`` per domain, the header
comment that repeats them, and ``domain_layouts`` in the root namelist. Anything
physical -- ``nz``, ``dz``, boundary conditions, output, timing -- is copied
byte-identical, because that is the part a namelist author chose deliberately.

The PE counts usually have to change with the grid. PALM splits a domain so that
``npex`` divides ``nx+1`` and ``npey`` divides ``ny+1``; 2-smooth point counts
such as 1024 = 2^10 or 1536 = 2^9*3 carry no factor 5, so every npes containing
one (200, 1600, 2400, 3600, ...) becomes impossible. Each value in the layout
file is checked and refused if it cannot be split, with the valid counts on
either side printed -- how many cores to spend is an allocation decision, so the
tool will not silently pick a different one.

The layout file names the PALM domains and their PE budgets:

    job: myrun                     # PALM job name -> <job>_p3d, <job>_static
    domains:
      root: {npes: 96,   palm_name: myrun}
      N02:  {npes: 256,  palm_name: myrun_pH}
      ...
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import yaml  # noqa: E402

from palmgeo import config as pgconfig  # noqa: E402


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def decompose(npes: int, nx: int, ny: int) -> tuple[int, int] | None:
    """npex * npey = npes with npex | nx+1 and npey | ny+1, closest to the aspect.

    The aspect is compared in log space so that npex/npey and npey/npex are
    equally far from a square domain; ties go to the wider decomposition.
    """
    want = math.log((nx + 1) / (ny + 1))
    best = None
    for npex in divisors(nx + 1):
        if npes % npex:
            continue
        npey = npes // npex
        if (ny + 1) % npey:
            continue
        score = (abs(math.log(npex / npey) - want), npey)
        if best is None or score < best[0]:
            best = (score, npex, npey)
    return (best[1], best[2]) if best else None


def suggest(npes: int, nx: int, ny: int, span: int = 8) -> list[int]:
    """Valid PE counts on either side of an impossible one."""
    below = [n for n in range(npes - 1, 0, -1) if decompose(n, nx, ny)][:span // 2]
    above = [n for n in range(npes + 1, npes * 16 + 16) if decompose(n, nx, ny)][:span // 2]
    return sorted(below + above)


def domain_layouts(cfg, layout, chosen) -> str:
    """The root namelist's nesting table; offsets are relative to the ROOT origin."""
    root = next(d for d in cfg.domains if d.domain_parent is None)
    ids = {d.name: i + 1 for i, d in enumerate(cfg.domains)}
    lines = [
        "!-- Domain layout: name, id, parent_id, npes, x_origin, y_origin, z_origin",
        "!-- Offsets are relative to the ROOT origin, not to the immediate parent.",
        "!-- npes satisfies npex | nx+1 and npey | ny+1; the split is noted per line.",
    ]
    width = max(len(layout["domains"][d.name]["palm_name"]) for d in cfg.domains) + 3
    for i, d in enumerate(cfg.domains):
        head = "   domain_layouts = " if i == 0 else "                    "
        name = f"'{layout['domains'][d.name]['palm_name']}',"
        npes, npex, npey = chosen[d.name]
        parent = ids[d.domain_parent] if d.domain_parent else -1
        lines.append(
            f"{head}{name:<{width}} {ids[d.name]:d},  {parent:2d}, "
            f"{npes:5d}, {d.origin_x - root.origin_x:9.1f}, "
            f"{d.origin_y - root.origin_y:9.1f},  0.0,   ! {npex} x {npey}"
        )
    return "\n".join(lines) + "\n"


def rewrite(text: str, d, is_root: bool, layouts: str) -> str:
    text = re.sub(r"^(\s*nx\s*=\s*)-?\d+,", rf"\g<1>{d.nx},", text, flags=re.M)
    text = re.sub(r"^(\s*ny\s*=\s*)-?\d+,", rf"\g<1>{d.ny},", text, flags=re.M)
    text = re.sub(r"^(\s*dx\s*=\s*)[\d.]+,", rf"\g<1>{d.pixel_size:.1f},", text, flags=re.M)
    text = re.sub(r"^(\s*dy\s*=\s*)[\d.]+,", rf"\g<1>{d.pixel_size:.1f},", text, flags=re.M)
    text = re.sub(r"^! nx = \d+, ny = \d+,", f"! nx = {d.nx}, ny = {d.ny},", text, flags=re.M)
    if is_root:
        text = re.sub(
            r"^!-- Domain layout:.*?\n(?:\s*(?:domain_layouts\s*=)?\s*'[^\n]*\n)+",
            layouts, text, flags=re.M | re.S)
    return text


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--layout", required=True)
    ap.add_argument("--src", required=True, help="directory holding the current namelists")
    ap.add_argument("--dst", help="defaults to the run's outdir")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cfg = pgconfig.load(args.config)
    layout = yaml.safe_load(Path(args.layout).read_text())
    job = str(layout["job"])
    src, dst = Path(args.src), Path(args.dst) if args.dst else cfg.outdir

    missing = [d.name for d in cfg.domains if d.name not in layout.get("domains", {})]
    if missing:
        raise SystemExit(f"{args.layout}: no entry for domain(s) {', '.join(missing)}")

    chosen, bad = {}, []
    for d in cfg.domains:
        npes = int(layout["domains"][d.name]["npes"])
        got = decompose(npes, d.nx, d.ny)
        if got is None:
            bad.append((d, npes))
            print(f"  {d.name:6s} {d.nx + 1:5d} x {d.ny + 1:5d} pts   "
                  f"npes {npes:5d}   IMPOSSIBLE")
            continue
        chosen[d.name] = (npes, got[0], got[1])
        print(f"  {d.name:6s} {d.nx + 1:5d} x {d.ny + 1:5d} pts   "
              f"npes {npes:5d} = {got[0]} x {got[1]}")

    if bad:
        print("\nPALM needs npex | nx+1 and npey | ny+1. These do not factorise:")
        for d, npes in bad:
            opts = ", ".join(str(n) for n in suggest(npes, d.nx, d.ny))
            print(f"  {d.name}: npes {npes} on {d.nx + 1} x {d.ny + 1} -> try {opts}")
        raise SystemExit("\nEdit the npes in the layout file; the choice is an "
                         "allocation decision, not something to guess at.")

    layouts = domain_layouts(cfg, layout, chosen)
    dst.mkdir(parents=True, exist_ok=True)
    written = []
    for d in cfg.domains:
        is_root = d.domain_parent is None
        suffix = "" if is_root else f"_{d.name}"
        for stem in (f"{job}_p3d", f"{job}_p3dr"):
            s = src / f"{stem}{suffix}"
            if not s.is_file():
                continue
            out = rewrite(s.read_text(), d, is_root, layouts)
            if not args.dry_run:
                target = dst / s.name
                if not target.exists() or not s.samefile(target):
                    shutil.copy2(s, target)   # keep mode and timestamps
                target.write_text(out)
            written.append(s.name)

    if not written:
        raise SystemExit(f"no namelists named {job}_p3d* found in {src}")
    print(f"\n{'would write' if args.dry_run else 'wrote'} {len(written)} namelist(s) to {dst}")
    for n in sorted(written):
        print("   ", n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
