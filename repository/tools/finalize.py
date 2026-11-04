#!/usr/bin/env python3
"""Bring a project directory into its final shape after palm_csd has run.

    python repository/tools/finalize.py static_driver/<run>/config/<run>.yml [--dry-run]

palm_csd knows nothing about palmgeo's directory convention. It writes every
driver and every plot side by side into one output directory, and it names the
root domain like a nest. This script is the single step that turns that into
the layout the rest of the project assumes:

    static_driver/<run>/
    ├── <run>_static            root driver     PALM reads exactly this name
    ├── <run>_static_N02 …      nests           the suffix is what marks a nest
    ├── config/                 what described the run
    ├── products/               what palm_csd read
    ├── docs/                   what a human looks at
    └── work/                   intermediate data, not published

It does four things, each of them idempotent -- running it twice changes
nothing, which matters because it is the last step of a long pipeline and gets
re-run after every partial failure.

1.  **Rename the root driver.** palm_csd names outputs ``<file_out>_<domain>``
    with no special case for the root, so a root domain called ``root`` comes
    out as ``<run>_static_root``. PALM looks for the root under the bare name
    ``<run>_static``; the ``_N0x`` suffix is precisely what identifies a nest.
2.  **Move the plots into docs/.** ``palm_csd --pdf`` drops them next to the
    drivers. They are documentation, not model input.
3.  **Record the region YAMLs in config/.** A driver whose land-use mapping is
    unknown cannot be interpreted, and the mapping is not recoverable from the
    NetCDF. The copies are a record of what was used; they are not what a
    re-run reads, and each carries a header saying so.
4.  **Report what is left over.** Anything at the top of the project directory
    that is neither a driver nor one of the four sub-directories is named, not
    moved. Guessing is worse than telling.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path
from typing import List

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from palmgeo import config as pgconfig  # noqa: E402

#: Prepended to every copied region YAML. Without it the copy is
#: indistinguishable from the original and invites edits that do nothing.
RECORD_HEADER = (
    "# ---------------------------------------------------------------------------\n"
    "# RECORD -- a copy of {src} as it stood when this\n"
    "# driver was built. It documents what was used; it is NOT what a re-run reads.\n"
    "# Editing this file changes nothing. Edit the original above instead.\n"
    "# ---------------------------------------------------------------------------\n"
)


def _say(dry: bool, verb: str, what: str) -> None:
    print(f"  {'would ' + verb if dry else verb + 'd':<12} {what}")


def rename_root(cfg, dry: bool) -> int:
    """Strip the ``_root`` suffix palm_csd gives the root domain."""
    root = next((d for d in cfg.domains if d.domain_parent is None), None)
    if root is None:
        print("  no root domain in this configuration -- nothing to rename")
        return 0
    base = str(cfg.get("project.name", "palm"))
    n = 0
    for suffix in ("", ".pdf", ".png"):
        src = cfg.outdir / f"{base}_{root.name}{suffix}"
        dst = cfg.outdir / f"{base}{suffix}"
        if not src.is_file():
            continue
        if dst.exists():
            # palm_csd has just written `src`, so it is by definition this
            # run's output and `dst` is a previous one. Skipping -- which is
            # what this did -- leaves the OLD driver in place under the name
            # PALM reads, with the new one beside it as `..._root`. The run
            # looks like it succeeded and the driver does not change, which is
            # the worst of the three possible outcomes.
            import filecmp

            if filecmp.cmp(src, dst, shallow=False):
                print(f"  identical    {dst.name} is unchanged by this run")
                if not dry:
                    src.unlink()
                continue
            _say(dry, "replace", f"{dst.name} with the newer {src.name}")
            if not dry:
                src.replace(dst)
            n += 1
            continue
        _say(dry, "move", f"{src.name} -> {dst.name}")
        if not dry:
            src.rename(dst)
        n += 1
    return n


def collect_plots(cfg, dry: bool) -> int:
    """Move the palm_csd sheets out of the driver directory into docs/."""
    n = 0
    for src in sorted(cfg.outdir.glob("*.pdf")) + sorted(cfg.outdir.glob("*.png")):
        dst = cfg.docsdir / src.name
        if dst.exists() and dst.stat().st_mtime >= src.stat().st_mtime:
            if not dry:
                src.unlink()          # an older duplicate of what docs/ holds
            continue
        _say(dry, "move", f"{src.name} -> docs/")
        if not dry:
            cfg.docsdir.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dst))
        n += 1
    return n


def record_region_files(cfg, dry: bool) -> int:
    """Snapshot the sources/mapping pair this run actually resolved."""
    n = 0
    for key in ("sources", "mapping"):
        rel = cfg.get(f"project.{key}")
        if not rel:
            continue
        src = cfg.root / str(rel)
        if not src.is_file():
            print(f"  missing      project.{key} -> {rel}")
            continue
        dst = cfg.configdir / src.name
        body = RECORD_HEADER.format(src=rel) + src.read_text(encoding="utf-8")
        if dst.is_file() and dst.read_text(encoding="utf-8") == body:
            continue
        _say(dry, "record", f"{rel} -> config/{src.name}")
        if not dry:
            cfg.configdir.mkdir(parents=True, exist_ok=True)
            dst.write_text(body, encoding="utf-8")
        n += 1
    return n


def report_strays(cfg) -> List[str]:
    """Name anything at the top of the project directory that does not belong.

    Deliberately does not move it. A stray file at this point is either a
    palmpaint autosave, a leftover from a renamed run, or something the author
    put there on purpose -- and this script cannot tell which.
    """
    base = str(cfg.get("project.name", "palm"))
    keep = {"config", "products", "docs", "work"}
    strays = []
    for entry in sorted(cfg.outdir.iterdir()):
        if entry.name in keep or entry.name.startswith("."):
            continue
        if entry.name == base or entry.name.startswith(f"{base}_N"):
            continue                                    # root driver or a nest
        strays.append(entry.name)
    return strays


def main(argv: List[str]) -> int:
    args = [a for a in argv[1:] if not a.startswith("-")]
    dry = "--dry-run" in argv
    if len(args) != 1:
        print(__doc__)
        return 2

    cfg = pgconfig.load(args[0])
    if not cfg.outdir.is_dir():
        print(f"{cfg.outdir} does not exist -- has palm_csd run?")
        return 1

    print(f"Finalising {cfg.outdir}")
    changed = (rename_root(cfg, dry)
               + collect_plots(cfg, dry)
               + record_region_files(cfg, dry))
    if changed == 0:
        print("  already in final shape")

    strays = report_strays(cfg)
    if strays:
        print("\n  Left alone, because only you know what they are:")
        for s in strays:
            print(f"    {s}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
