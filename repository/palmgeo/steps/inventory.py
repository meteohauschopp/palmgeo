"""What data exists for this area, what it may be used for, and what it becomes.

This stage runs before anything is downloaded in bulk. It answers three
questions that a static driver otherwise buries:

**What is actually available here?** A region file lists what a region
*publishes*; it does not say whether any of it covers the particular few square
kilometres you asked for. The scan loads each source's index and measures
coverage over the area of interest, so a source that stops at a state border
shows up as a number rather than as a puzzling gap two hours into a download.

**Under what licence?** Every source names licences by id from
``repository/config/licences.yml``, and a driver inherits all of them at once.
The interesting quantity is therefore not any single licence but their
combination -- the strictest obligation from any input governs the whole
driver. One CC BY-NC-SA input makes the driver non-commercial and share-alike
no matter how permissive the other five are, and that is worth knowing before
the download rather than before the publication.

**What part of the driver rests on it?** ``palmgeo/roles.py`` carries, per
role, which driver variables it produces and what the driver looks like without
it. Printing that beside each source turns "no topo source" from a warning into
a decision: this driver will have no street types.

The stage writes its result to ``docs/DATA-SOURCES.md`` and records an approval
fingerprint in ``state.json``. The fingerprint covers the chosen source per role
and the licences those sources carry, so swapping a surface model or a licence
changing under a source invalidates the approval and asks again. Everything
else about a run -- domains, resolutions, thresholds -- deliberately does not,
because none of it changes what may be done with the result.
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from palmgeo import roles as roles_mod
from palmgeo.cache import Cache
from palmgeo.config import Config
from palmgeo.grid import Bounds
from palmgeo.sources import Index, load_index

logger = logging.getLogger("palmgeo.inventory")

#: Resolution of the coverage probe, in cells per side over the AOI. 160x160
#: is 25 600 point-in-box tests against a spatially hashed index -- milliseconds
#: -- and resolves a hole of about 1/160 of the domain, which is far finer than
#: any decision made from this number.
_PROBE = 160

#: The obligations a licence can impose, in the order they are reported. The
#: order is by how much each one constrains what you may do with the driver.
_OBLIGATIONS = (
    ("non_commercial", "non-commercial use only"),
    ("share_alike", "share-alike: a distributed derivative carries this licence"),
    ("redistribution_restricted", "the data itself may not be redistributed"),
    ("attribution", "attribution required"),
)


# ---------------------------------------------------------------------------
# data model
# ---------------------------------------------------------------------------
@dataclass
class SourceReport:
    """One source, as the report sees it."""

    name: str
    role: str
    chosen: bool
    #: Fraction of the AOI covered, 0..1, or None where it cannot be measured
    #: from an index (a single statewide file has no tile geometry to test).
    coverage: Optional[float]
    #: How the coverage figure was arrived at, so a reader can weigh it.
    coverage_basis: str
    licence_ids: List[str]
    #: Rights holder per licence id -- see Config.licence_holder.
    holders: Dict[str, Optional[str]]
    resolution: Optional[float]
    #: Oldest and newest acquisition date the index offers over the work tiles,
    #: as ISO strings. None where the source states no date at all.
    #:
    #: This is the question "how old is what I am about to build from", and it
    #: has no other home: the chosen orthophoto epoch is reported by `plan`,
    #: but the elevation vintage -- which decides whether a building exists at
    #: all -- was recorded nowhere before.
    date_min: Optional[str] = None
    date_max: Optional[str] = None
    #: A source the region file ships but switches off. Listed anyway, because
    #: "this region also publishes X, and here is why it is not used" is the
    #: kind of thing that is impossible to reconstruct later.
    enabled: bool = True
    note: Optional[str] = None
    error: Optional[str] = None


@dataclass
class RoleReport:
    """One role: what fills it, what else could, and what depends on it."""

    role: str
    required: bool
    summary: str
    feeds: str
    without: str
    chosen: Optional[SourceReport]
    alternatives: List[SourceReport] = field(default_factory=list)
    problem: Optional[str] = None


# ---------------------------------------------------------------------------
# coverage
# ---------------------------------------------------------------------------
def _work_tiles(cfg: Config) -> List[Bounds]:
    """The work tiles this run would actually fetch, deduplicated.

    Not the bounding box over the domains. Two nests 19 km apart span a box two
    and a half times their combined area, and probing that box would count
    ground the run never touches -- reporting a source as covering 60 % of "the
    area" when it covers all of every tile that will be read. Coverage has to
    be measured over the same tiles the download stage will ask for, which is
    what ``Tier.tiles`` already enumerates for exactly this reason.
    """
    seen: Dict[Tuple[int, int], Bounds] = {}
    for tier in cfg.tiers:
        for tile in tier.tiles:
            seen[(tile.east_km, tile.north_km)] = tile.bounds
    return [seen[k] for k in sorted(seen)]


def _aoi(cfg: Config) -> Bounds:
    """Bounding box of the work tiles -- for display and for sizing the probe."""
    tiles = _work_tiles(cfg)
    if not tiles:
        xs0, ys0, xs1, ys1 = zip(*(t.bounds for t in cfg.tiers))
        return (min(xs0), min(ys0), max(xs1), max(ys1))
    xs0, ys0, xs1, ys1 = zip(*tiles)
    return (min(xs0), min(ys0), max(xs1), max(ys1))


def coverage_fraction(
    index: Index, tiles: Sequence[Bounds], probe: int = _PROBE
) -> float:
    """Fraction of the work tiles that at least one index item covers.

    Measured by probing points rather than by unioning geometry. The items
    overlap each other, they are bounding boxes rather than true footprints,
    and a union would imply a precision this figure does not have; a probe grid
    is honest about being an estimate and costs nothing.

    The probe budget is spread over the tiles rather than over their bounding
    box, so that a run with widely separated domains spends its samples on
    ground that will actually be read.

    Bounding boxes overstate coverage where a tile is partly nodata -- a
    coastal or border tile is a full box with half a raster in it -- so the
    result is an upper bound, and the report says so rather than presenting it
    as measured coverage.
    """
    if not tiles:
        return 0.0
    # Samples per tile side: enough that the whole run gets about `probe`^2
    # samples, but never fewer than a 4x4 grid on any single tile.
    per_side = max(4, int(round(probe / max(1.0, len(tiles) ** 0.5))))
    hits = total = 0
    for x0, y0, x1, y1 in tiles:
        dx, dy = (x1 - x0) / per_side, (y1 - y0) / per_side
        for iy in range(per_side):
            cy = y0 + (iy + 0.5) * dy
            for ix in range(per_side):
                total += 1
                if index.contains(x0 + (ix + 0.5) * dx, cy):
                    hits += 1
    return hits / float(total) if total else 0.0


def _dates_over(index: Index, tiles: Sequence[Bounds]) -> Tuple[Optional[str], Optional[str]]:
    """Oldest and newest acquisition date the index offers over these tiles.

    Taken over the tiles the run will actually read, not over the whole index:
    a statewide orthophoto index spans a decade, while the campaign covering
    one city may be a single day, and the second number is the one that
    describes the driver.
    """
    dates = sorted(
        {it.date for b in tiles for it in index.search(b) if it.date is not None}
    )
    if not dates:
        return None, None
    return dates[0].isoformat(), dates[-1].isoformat()


def _probe_source(
    cfg: Config, cache: Cache, name: str, src: Dict[str, Any],
    tiles: Sequence[Bounds],
) -> Tuple[Optional[float], str, Optional[str], Optional[str], Optional[str]]:
    """Measure what can be measured about one source's coverage of the AOI.

    Returns (fraction, how it was established, error). A source that cannot be
    probed is not a failure -- a single statewide file genuinely has no index
    to interrogate -- so the basis string carries that distinction instead of
    an exception.
    """
    kind = str(src.get("kind", "")).lower()
    if kind in ("tile_index", "region_index"):
        try:
            index = load_index(name, src, cache, transform_to_epsg=cfg.epsg)
        except Exception as exc:  # noqa: BLE001 - reported, not raised
            return None, "index could not be read", f"{type(exc).__name__}: {exc}", None, None
        if len(index) == 0:
            return 0.0, "index is empty", None, None, None
        frac = coverage_fraction(index, tiles)
        lo, hi = _dates_over(index, tiles)
        unit = "tiles" if kind == "tile_index" else "administrative units"
        return frac, (f"probed over {len(tiles)} work tile(s) against {len(index)} "
                      f"{unit} (bounding boxes, upper bound)"), None, lo, hi
    if kind == "single_file":
        return None, "single file, no index to test coverage against", None, None, None
    return None, f"unknown source kind '{kind}'", None, None, None


# ---------------------------------------------------------------------------
# the scan
# ---------------------------------------------------------------------------
def scan(cfg: Config, cache: Cache, probe_alternatives: bool = True) -> List[RoleReport]:
    """Build one RoleReport per role the pipeline knows about."""
    tiles = _work_tiles(cfg)
    reports: List[RoleReport] = []

    for role in roles_mod.ROLES:
        names = cfg.sources_for_role(role.name)
        chosen_name = cfg.source_name_for_role(role.name)
        problem: Optional[str] = None

        if not names:
            problem = (
                f"no source in this region file declares role '{role.name}'"
                if not role.optional else None
            )
        elif chosen_name is None:
            problem = (
                f"{len(names)} sources declare this role and nothing chooses "
                f"between them: {', '.join(sorted(names))}"
            )

        entries: List[SourceReport] = []
        for name in names:
            is_chosen = name == chosen_name
            # Probing an alternative costs an index download it may not be
            # worth; the chosen source is always probed, the rest on request.
            src = cfg.source(name)
            enabled = src.get("enabled", True) is not False
            # A disabled source is never downloaded, so probing it would spend
            # an index fetch on a number nobody can act on.
            lo = hi = None
            if enabled and (is_chosen or probe_alternatives):
                frac, basis, err, lo, hi = _probe_source(cfg, cache, name, src, tiles)
            elif not enabled:
                frac, basis, err = None, "switched off in the region file", None
            else:
                frac, basis, err = None, "not probed", None
            entries.append(SourceReport(
                name=name,
                role=role.name,
                chosen=is_chosen,
                coverage=frac,
                coverage_basis=basis,
                licence_ids=cfg.licence_ids_for_source(name),
                holders={lid: cfg.licence_holder(name, lid)
                         for lid in cfg.licence_ids_for_source(name)},
                resolution=src.get("resolution"),
                date_min=lo,
                date_max=hi,
                enabled=enabled,
                note=_first_line(src.get("notes")),
                error=err,
            ))

        chosen = next((e for e in entries if e.chosen), None)
        reports.append(RoleReport(
            role=role.name,
            required=not role.optional,
            summary=role.summary,
            feeds=roles_mod.describe(role.name),
            without=role.without,
            chosen=chosen,
            alternatives=[e for e in entries if not e.chosen],
            problem=problem,
        ))
    return reports


def _first_line(value: Any) -> Optional[str]:
    """The first sentence of a source's notes, for a one-line display."""
    if not value:
        return None
    text = " ".join(str(value).split())
    return text if len(text) <= 110 else text[:107] + "..."


# ---------------------------------------------------------------------------
# licence arithmetic
# ---------------------------------------------------------------------------
def combined_obligations(cfg: Config, reports: Sequence[RoleReport]) -> Dict[str, List[str]]:
    """Which obligations the finished driver carries, and which input imposed each.

    The combination rule is the only one that is safe: an obligation applies to
    the driver if *any* input imposes it. There is no averaging and no majority
    -- one share-alike input among six makes the driver share-alike.

    Only *chosen* sources count. An alternative that would have been
    non-commercial imposes nothing as long as it is not used, which is exactly
    the information that makes the alternatives worth listing.
    """
    result: Dict[str, List[str]] = {key: [] for key, _ in _OBLIGATIONS}
    for report in reports:
        if report.chosen is None:
            continue
        for lid in report.chosen.licence_ids:
            entry = cfg.licences.get(lid) or {}
            for key, _ in _OBLIGATIONS:
                if (entry.get("obligations") or {}).get(key):
                    result[key].append(f"{report.chosen.name} ({lid})")
    return result


def used_licences(cfg: Config, reports: Sequence[RoleReport]) -> List[str]:
    """Every distinct licence id the chosen sources carry, in report order."""
    seen: List[str] = []
    for report in reports:
        if report.chosen is None:
            continue
        for lid in report.chosen.licence_ids:
            if lid not in seen:
                seen.append(lid)
    return seen


def attribution_lines(cfg: Config, reports: Sequence[RoleReport]) -> Tuple[List[str], List[str]]:
    """The credits, split into the ones that are required and the ones that are not.

    The distinction matters. A US federal work imposes no attribution at all --
    crediting USGS is good practice, not an obligation -- and listing it under
    "must travel with the driver" would misstate the licence in the one
    direction that costs the reader something: it makes a courtesy look like a
    condition, and the next person treats every line in the list as optional.

    Returns (required, customary).

    Where a template wants a year this stage cannot know -- acquisition years
    differ per tile and are chosen later, by the epoch selection -- the
    placeholder is left visible rather than filled with a guess, because a
    wrong year in a credit is worse than a blank one.
    """
    required: List[str] = []
    customary: List[str] = []
    for report in reports:
        if report.chosen is None:
            continue
        for lid in report.chosen.licence_ids:
            entry = cfg.licences.get(lid) or {}
            template = entry.get("attribution")
            if not template:
                continue
            holder = report.chosen.holders.get(lid) or "unknown holder"
            line = str(template).replace("{holder}", holder).replace("{year}", "<year>")
            target = required if (entry.get("obligations") or {}).get("attribution") else customary
            if line not in target:
                target.append(line)
    return required, customary



def currency(reports: Sequence[RoleReport]) -> List[Tuple[str, str, str, str]]:
    """Per chosen source: role, name, the acquisition window, and its age.

    Reported, not judged. Whether nine years between the elevation model and
    the imagery matters depends entirely on the study: for a terrain-driven
    flow question it is irrelevant, and for a newly built quarter it decides
    whether the buildings exist at all. palmgeo states the numbers and leaves
    the conclusion to somebody who knows which one they are doing.
    """
    from datetime import date as _date

    today = _date.today()
    out: List[Tuple[str, str, str, str]] = []
    for r in reports:
        c = r.chosen
        if c is None:
            continue
        if not c.date_min:
            out.append((r.role, c.name, "no date stated", ""))
            continue
        window = c.date_min if c.date_min == c.date_max else f"{c.date_min} … {c.date_max}"
        newest = _date.fromisoformat(c.date_max or c.date_min)
        years = (today - newest).days / 365.25
        age = "under a year" if years < 1 else f"{years:.0f} years"
        out.append((r.role, c.name, window, age))
    return out


def currency_span(reports: Sequence[RoleReport]) -> Optional[Tuple[str, str, float]]:
    """Oldest and newest acquisition across all chosen sources, and the gap."""
    from datetime import date as _date

    dates = [
        _date.fromisoformat(d)
        for r in reports if r.chosen
        for d in (r.chosen.date_min, r.chosen.date_max) if d
    ]
    if len(dates) < 2:
        return None
    lo, hi = min(dates), max(dates)
    return lo.isoformat(), hi.isoformat(), (hi - lo).days / 365.25


# ---------------------------------------------------------------------------
# approval fingerprint
# ---------------------------------------------------------------------------
def fingerprint(reports: Sequence[RoleReport]) -> str:
    """A stable digest of the decisions this approval covers.

    Deliberately narrow: the chosen source per role and the licences it
    carries, and nothing else. Changing a domain origin or an NDVI threshold
    does not change what the data may be used for, so it must not invalidate a
    licence approval and make the gate cry wolf. Changing a source, or a
    licence changing under a source, does.
    """
    payload = {
        r.role: {
            "source": r.chosen.name if r.chosen else None,
            "licences": sorted(r.chosen.licence_ids) if r.chosen else [],
        }
        for r in reports
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


# ---------------------------------------------------------------------------
# rendering
# ---------------------------------------------------------------------------
def _pct(value: Optional[float]) -> str:
    if value is None:
        return "     -"
    return f"{value * 100:5.1f} %"


def render(cfg: Config, reports: Sequence[RoleReport]) -> str:
    """The terminal report."""
    out: List[str] = []
    region = cfg.sources.get("region", "?")
    out.append(f"\nData sources for '{cfg.get('project.name')}'  --  region {region}, {cfg.crs}")
    x0, y0, x1, y1 = _aoi(cfg)
    n_tiles = len(_work_tiles(cfg))
    out.append(f"Area of interest: {n_tiles} work tile(s) within "
               f"{(x1 - x0) / 1000:.1f} x {(y1 - y0) / 1000:.1f} km\n")

    for r in reports:
        flag = "REQUIRED" if r.required else "optional"
        head = f"  {r.role:<12} {flag:<9}"
        if r.chosen is None:
            out.append(f"{head} -- NOT AVAILABLE")
            out.append(f"                         {r.summary}")
            if r.problem:
                out.append(f"                         ! {r.problem}")
            if r.without:
                out.append(f"                         -> {r.without}")
            out.append("")
            continue

        c = r.chosen
        lic = ", ".join(c.licence_ids)
        out.append(f"{head} {c.name:<22} {_pct(c.coverage)}  {lic}")
        out.append(f"                         becomes: {r.feeds}")
        out.append(f"                         {c.coverage_basis}")
        if c.error:
            out.append(f"                         ! {c.error}")
        for alt in r.alternatives:
            label = "alternative" if alt.enabled else "alt (off) "
            out.append(f"      {label}        {alt.name:<22} {_pct(alt.coverage)}  "
                       f"{', '.join(alt.licence_ids)}")
            if alt.note:
                out.append(f"                           {alt.note}")
        out.append("")

    # -- how old the data is, which nothing else records ---------------------
    rows = currency(reports)
    if rows:
        out.append("  " + "-" * 70)
        out.append("  Currency of the data")
        out.append("  " + "-" * 70)
        out.append(f"    {'role':<12} {'source':<20} {'acquired':<26} age")
        for role, name, window, age in rows:
            out.append(f"    {role:<12} {name:<20} {window:<26} {age}")
        span = currency_span(reports)
        if span:
            out.append("")
            out.append(f"    Oldest to newest: {span[0]} to {span[1]}, "
                       f"a span of {span[2]:.0f} years.")
            out.append("    Building geometry comes from the surface model, so anything")
            out.append("    built after its acquisition has no height there.")
        out.append("")

    # -- the part that decides what may be done with the result --------------
    out.append("  " + "-" * 70)
    out.append("  Licence position of the resulting driver")
    out.append("  " + "-" * 70)
    for lid in used_licences(cfg, reports):
        entry = cfg.licences.get(lid) or {}
        out.append(f"    {lid:<24} {entry.get('name', '?')}")
        if entry.get("url"):
            out.append(f"    {'':<24} {entry['url']}")
    obligations = combined_obligations(cfg, reports)
    out.append("")
    any_ob = False
    for key, text in _OBLIGATIONS:
        who = obligations.get(key) or []
        if not who:
            continue
        any_ob = True
        out.append(f"    [x] {text}")
        out.append(f"        from: {', '.join(sorted(set(who)))}")
    if not any_ob:
        out.append("    No obligations from any chosen source.")

    caveats = [
        (lid, (cfg.licences.get(lid) or {}).get("caveat"))
        for lid in used_licences(cfg, reports)
    ]
    caveats = [(lid, c) for lid, c in caveats if c]
    if caveats:
        out.append("\n    Known grey areas, read before publishing:")
        for lid, caveat in caveats:
            out.append(f"      {lid}: {' '.join(str(caveat).split())}")

    required, customary = attribution_lines(cfg, reports)
    if required:
        out.append("\n    Attribution that must travel with the driver:")
        out += [f"      {line}" for line in required]
    if customary:
        out.append("\n    Customary credit, not required by any licence:")
        out += [f"      {line}" for line in customary]
    if any("<year>" in line for line in required + customary):
        out.append("      (<year> is the acquisition year, which the epoch "
                   "selection in `plan` establishes per source)")

    out.append("\n  This is a summary of licence terms, not legal advice. The URLs above")
    out.append("  are authoritative and the terms can change.\n")
    return "\n".join(out)


def markdown(cfg: Config, reports: Sequence[RoleReport]) -> str:
    """``docs/DATA-SOURCES.md`` -- the same content, for the project record.

    Written as a document rather than a log because it outlives the run: it is
    what a reader of the finished driver consults to find out whether they may
    use it, and it is the material the project README draws on.
    """
    region = cfg.sources.get("region", "?")
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    out: List[str] = [
        f"# Data sources — {cfg.get('project.name')}",
        "",
        f"Region `{region}`, {cfg.crs}. Established {now} by "
        f"`palmgeo sources`, approval `{fingerprint(reports)}`.",
        "",
        "This file records what this driver is built from and what that means "
        "for its use. It is generated; edit the run configuration, not this file.",
        "",
        "## Sources",
        "",
        "| Role | | Source | Coverage of AOI | Licence | Becomes |",
        "|---|---|---|---|---|---|",
    ]
    for r in reports:
        need = "required" if r.required else "optional"
        if r.chosen is None:
            out.append(f"| `{r.role}` | {need} | — | — | — | *not available* |")
            continue
        c = r.chosen
        cover = f"{c.coverage * 100:.1f} %" if c.coverage is not None else "not measurable"
        out.append(
            f"| `{r.role}` | {need} | `{c.name}` | {cover} | "
            f"{', '.join(c.licence_ids)} | {r.feeds} |"
        )
        for alt in r.alternatives:
            out.append(
                f"| | *{'alternative' if alt.enabled else 'alternative, switched off'}* "
            f"| `{alt.name}` | "
                f"{f'{alt.coverage * 100:.1f} %' if alt.coverage is not None else 'not measurable'} | "
                f"{', '.join(alt.licence_ids)} | |"
            )

    rows = currency(reports)
    if rows:
        out += ["", "## Currency", "",
                "How old the data is. Nothing else in the driver records this, and it "
                "is what explains a building that is missing rather than wrong: "
                "building geometry comes from the surface model, so anything built "
                "after that model was flown has no height and drops out.",
                "",
                "The date is whatever the source declares as its currency field, and "
                "what that means differs per product: for an elevation model it is "
                "normally the survey, but for a derived product it can be the "
                "publication of a model whose geometry is older -- a 3D building "
                "model can carry this year's date and still lack every building its "
                "years-old laser scan did not see. Treat the column as an upper "
                "bound on how new the content is.",
                "",
                "| Role | Source | Acquired | Age |", "|---|---|---|---|"]
        for role, name, window, age in rows:
            out.append(f"| `{role}` | `{name}` | {window} | {age} |")
        span = currency_span(reports)
        if span:
            out += ["", f"Oldest to newest: **{span[0]} to {span[1]}**, a span of "
                        f"{span[2]:.0f} years.", ""]

    missing = [r for r in reports if r.chosen is None]
    if missing:
        out += ["", "## What is missing, and what that costs", ""]
        for r in missing:
            out.append(f"**`{r.role}`** — {r.summary}")
            out.append("")
            out.append(f"{r.without}")
            out.append("")

    out += ["", "## Licence position", ""]
    for lid in used_licences(cfg, reports):
        entry = cfg.licences.get(lid) or {}
        out.append(f"### {entry.get('name', lid)} (`{lid}`)")
        out.append("")
        if entry.get("url"):
            out.append(f"<{entry['url']}>")
            out.append("")
        if entry.get("summary"):
            out.append(" ".join(str(entry["summary"]).split()))
            out.append("")
        if entry.get("caveat"):
            out.append(f"> **Grey area.** {' '.join(str(entry['caveat']).split())}")
            out.append("")

    out += ["## What this driver requires of you", ""]
    obligations = combined_obligations(cfg, reports)
    any_ob = False
    for key, text in _OBLIGATIONS:
        who = obligations.get(key) or []
        if who:
            any_ob = True
            out.append(f"- **{text}** — from {', '.join(sorted(set(who)))}")
    if not any_ob:
        out.append("- Nothing. No chosen source imposes an obligation.")
    out.append("")

    required, customary = attribution_lines(cfg, reports)
    if required:
        out += ["### Attribution (required)", "",
                "This credit must travel with the driver and with anything "
                "published from it:", ""]
        out += [f"> {line}" for line in required]
        out.append("")
    if customary:
        out += ["### Credit (customary, not required)", ""]
        out += [f"> {line}" for line in customary]
        out.append("")
    if any("<year>" in line for line in required + customary):
        out += ["`<year>` is left for you to fill in, because what it means is the "
                "licence holder's call: for most it is the year the data was "
                "obtained, for Copernicus the year of acquisition. A wrong year in "
                "a credit is worse than an open one.", ""]

    out += ["---", "",
            "*A summary of licence terms, not legal advice. The linked texts are "
            "authoritative and can change.*", ""]
    return "\n".join(out)


def write_markdown(cfg: Config, reports: Sequence[RoleReport]) -> Path:
    cfg.docsdir.mkdir(parents=True, exist_ok=True)
    path = cfg.docsdir / "DATA-SOURCES.md"
    path.write_text(markdown(cfg, reports), encoding="utf-8")
    return path
