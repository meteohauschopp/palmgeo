"""Source indices, tile lookup and orthophoto epoch selection.

A *source index* is how palmgeo finds out which files cover a given patch of
ground. The shape it expects is a GeoJSON ``FeatureCollection`` with one
feature per available file, carrying that file's footprint as geometry and its
download URL and acquisition date as properties. Which property holds which is
declared per source in ``sources_<region>.yml``; nothing here knows a product
name.

Regions that publish such an index (Lower Saxony's LGLN does, one per product)
are read directly. Regions that do not are served by a staging script under
``repository/tools/`` that writes an index in this shape over whatever the
region actually offers -- so the difference between a well-served region and a
badly served one is confined to one script, and never reaches this module.

**Streaming.** These files reach 78 MB, and materialising one as Python objects
costs several hundred megabytes of geometry that is discarded immediately. They
are therefore parsed one feature at a time and condensed into a compact local
index that later runs reuse.

**Lookup by geometry, not by identifier.** Tiles are located by bounding-box
intersection through a grid hash rather than by parsing tile ids. That is
deliberate and was learned the hard way: the LoD2 index ships a malformed
``tile_id`` (``"25985897"`` for ``"325985897"``), and even where ids are sound
they are not comparable across products, because ortho tiles are 2 km where
terrain tiles are 1 km. Geometry is the only thing every index agrees on.
"""

from __future__ import annotations

import json
import logging
import math
import re
import time
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from palmgeo.cache import Cache
from palmgeo.grid import Bounds, intersects

logger = logging.getLogger("palmgeo.sources")

#: Grid hash cell size for the spatial lookup, in metres.
_HASH_CELL = 2000.0

#: Refresh a condensed index if it is older than this (days).
INDEX_TTL_DAYS = 30


# ---------------------------------------------------------------------------
# streaming GeoJSON
# ---------------------------------------------------------------------------
def stream_features(path: Path, chunk_size: int = 1 << 20) -> Iterator[Dict[str, Any]]:
    """Yield the features of a GeoJSON ``FeatureCollection`` one at a time.

    Written by hand rather than delegating to OGR/pyogrio because those read the
    whole collection into memory, which for the 78 MB DOP20 index means several
    hundred MB of resident geometry we immediately throw away.

    The scanner tracks brace depth while honouring JSON string literals and
    escapes, so it does not depend on the file being one-feature-per-line.
    """
    with path.open("r", encoding="utf-8") as fh:
        buf = ""
        # Advance to the opening bracket of the "features" array.
        while '"features"' not in buf:
            chunk = fh.read(chunk_size)
            if not chunk:
                raise ValueError(f"{path}: no 'features' member found")
            buf += chunk
        pos = buf.index("[", buf.index('"features"')) + 1

        depth = 0
        in_string = False
        escaped = False
        feature_start = -1
        eof = False
        while True:
            if pos >= len(buf):
                if eof:
                    return
                # Drop the consumed prefix exactly once per refill; slicing per
                # feature would make this quadratic in the file size.
                keep = feature_start if feature_start >= 0 else pos
                buf = buf[keep:]
                pos -= keep
                if feature_start >= 0:
                    feature_start = 0
                chunk = fh.read(chunk_size)
                if not chunk:
                    eof = True
                    if pos >= len(buf):
                        return
                else:
                    buf += chunk
                continue

            ch = buf[pos]
            pos += 1
            if in_string:
                if escaped:
                    escaped = False
                elif ch == "\\":
                    escaped = True
                elif ch == '"':
                    in_string = False
                continue
            if ch == '"':
                in_string = True
            elif ch == "{":
                if depth == 0:
                    feature_start = pos - 1
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0 and feature_start >= 0:
                    yield json.loads(buf[feature_start:pos])
                    feature_start = -1
            elif ch == "]" and depth == 0:
                return


def geometry_bounds(geometry: Dict[str, Any]) -> Optional[Bounds]:
    """Bounding box of any GeoJSON geometry, without building shapely objects."""
    if not geometry:
        return None
    coords = geometry.get("coordinates")
    if coords is None:
        return None
    xmin = ymin = math.inf
    xmax = ymax = -math.inf

    def walk(node: Any) -> None:
        nonlocal xmin, ymin, xmax, ymax
        if not node:
            return
        first = node[0]
        if isinstance(first, (int, float)):
            x, y = float(node[0]), float(node[1])
            xmin, ymin = min(xmin, x), min(ymin, y)
            xmax, ymax = max(xmax, x), max(ymax, y)
            return
        for child in node:
            walk(child)

    walk(coords)
    if xmin is math.inf:
        return None
    return (xmin, ymin, xmax, ymax)


def _parse_date(value: Any, fmt: Optional[str] = None) -> Optional[date]:
    """An acquisition date: ISO 8601 by default, else the source's ``date_format``."""
    if not value:
        return None
    try:
        if fmt:
            return datetime.strptime(str(value).strip(), str(fmt)).date()
        return datetime.strptime(str(value)[:10], "%Y-%m-%d").date()
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# condensed index
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Item:
    """One entry of a source index: one file, and where and when it applies.

    Deliberately four fields. An earlier version also carried every remaining
    GeoJSON property in an ``extra`` dict, which cost a dict per feature --
    around 100 000 of them for a statewide orthophoto index -- and was then
    dropped by the condensing step, so it existed on a first run and was empty
    on every run after. Nothing ever read it. If a region port needs a further
    property, the honest way is to give it a named field here and a declared
    key in ``sources_<region>.yml``, so that it survives condensing and means
    the same thing on every run.
    """

    bounds: Bounds
    url: str
    date: Optional[date]
    #: Stable identifier for the tile, used to pair products across sources.
    #: Derived from the file name where the source declares a pattern, and
    #: otherwise from the lower-left corner in whole kilometres.
    key: str

    @property
    def year(self) -> int:
        return self.date.year if self.date else 0

    @property
    def month(self) -> int:
        return self.date.month if self.date else 0


class Index:
    """Spatially searchable, condensed source index."""

    def __init__(self, name: str, items: Sequence[Item], crs: Optional[str]) -> None:
        self.name = name
        self.items = list(items)
        self.crs = crs
        self._hash: Dict[Tuple[int, int], List[int]] = {}
        for i, item in enumerate(self.items):
            for cell in self._cells(item.bounds):
                self._hash.setdefault(cell, []).append(i)

    @staticmethod
    def _cells(b: Bounds) -> Iterator[Tuple[int, int]]:
        """Every hash cell a bounding box touches.

        The epsilon on the upper edge stops a box that ends exactly on a cell
        boundary from claiming the next cell along, which would put every
        tile-aligned item into twice as many buckets as it needs.

        ``max(x0, ...)`` is what makes that safe for a box with no extent. A
        zero-width query -- a single point, as the coverage probe in
        ``steps/inventory.py`` uses -- would otherwise compute an upper cell
        one below its lower cell whenever the point sits exactly on a boundary,
        giving an empty range and therefore no results at all. The search would
        report "nothing covers this location" for a location that is covered.
        """
        x0 = int(math.floor(b[0] / _HASH_CELL))
        x1 = max(x0, int(math.floor((b[2] - 1e-6) / _HASH_CELL)))
        y0 = int(math.floor(b[1] / _HASH_CELL))
        y1 = max(y0, int(math.floor((b[3] - 1e-6) / _HASH_CELL)))
        for yy in range(y0, y1 + 1):
            for xx in range(x0, x1 + 1):
                yield (xx, yy)

    def search(self, b: Bounds) -> List[Item]:
        """All items whose bounding box overlaps ``b``.

        Overlap is strict: ``intersects`` requires positive area, so a query
        box with no extent matches nothing regardless of where it sits. Callers
        that want to ask "what covers this point" must give the point a small
        extent; ``contains`` below does that for them.
        """
        seen: set[int] = set()
        for cell in self._cells(b):
            seen.update(self._hash.get(cell, ()))
        return [self.items[i] for i in sorted(seen) if intersects(self.items[i].bounds, b)]

    def contains(self, x: float, y: float, eps: float = 1e-6) -> bool:
        """Whether any item covers the point ``(x, y)``.

        Exists because the obvious spelling -- ``search((x, y, x, y))`` -- is
        wrong twice over: a zero-area box fails the strict overlap test, and
        on a hash-cell boundary it used to address no bucket at all. Both are
        handled here once rather than at each call site.
        """
        return bool(self.search((x - eps, y - eps, x + eps, y + eps)))

    def __len__(self) -> int:
        return len(self.items)


def load_index(
    name: str,
    definition: Dict[str, Any],
    cache: Cache,
    transform_to_epsg: Optional[int] = None,
) -> Index:
    """Download (if needed), condense and return a source index."""
    condensed = cache.dir / "indices" / f"{name}.json"
    if condensed.is_file() and not _condensed_matches(condensed, transform_to_epsg):
        # The cache directory is shared between runs, and two runs may use the
        # same source in different CRSs (a global source in two UTM zones).
        # Condensed bounds in the wrong CRS would place every tile elsewhere.
        if getattr(cache, "offline", False):
            raise RuntimeError(
                f"Index {name} is cached in another CRS than EPSG:{transform_to_epsg} "
                f"and cannot be rebuilt offline."
            )
        logger.info("Index %s is cached in another CRS, rebuilding", name)
    elif condensed.is_file():
        age_days = (time.time() - condensed.stat().st_mtime) / 86400
        if age_days < INDEX_TTL_DAYS:
            return _read_condensed(name, condensed)
        if getattr(cache, "offline", False):
            # The TTL exists so that a long-lived cache notices new
            # acquisitions upstream. Offline there is nothing to notice, and a
            # stale index is strictly better than no run at all: it describes
            # data that was there a month ago, which is what is on disk anyway.
            logger.warning(
                "Index %s is %.0f days old and cannot be refreshed offline. "
                "Using it as it stands; it may not list acquisitions published "
                "since.", name, age_days,
            )
            return _read_condensed(name, condensed)
        logger.info("Index %s is %.0f days old, refreshing", name, age_days)

    url = definition.get("index_url")
    if not url:
        raise ValueError(f"Source '{name}' has no index_url")
    raw = cache.fetch(url, subdir="indices")

    url_field = definition.get("url_field")
    url_fallback = definition.get("url_field_fallback")
    date_field = definition.get("date_field")
    date_format = definition.get("date_format")
    name_regex = definition.get("tile_from_filename")
    pattern = re.compile(name_regex) if name_regex else None
    # Required, and checked when the configuration is loaded: there is no CRS
    # an index could safely be assumed to be in.
    index_crs = definition.get("index_crs")
    if not index_crs:
        raise ValueError(f"Source '{name}' has no index_crs")
    index_crs = str(index_crs)
    reproject = None
    if transform_to_epsg is not None:
        from pyproj import CRS, Transformer

        src_crs = CRS.from_user_input(index_crs)
        if src_crs != CRS.from_epsg(transform_to_epsg):
            reproject = Transformer.from_crs(src_crs, transform_to_epsg, always_xy=True)

    items: List[Item] = []
    skipped = 0
    for feature in stream_features(raw):
        props = feature.get("properties") or {}
        href = props.get(url_field) or (props.get(url_fallback) if url_fallback else None)
        if not href:
            skipped += 1
            continue
        b = geometry_bounds(feature.get("geometry") or {})
        if b is None:
            skipped += 1
            continue
        if reproject is not None:
            xs, ys = reproject.transform([b[0], b[2]], [b[1], b[3]])
            b = (min(xs), min(ys), max(xs), max(ys))
        key = ""
        if pattern is not None:
            match = pattern.search(href.rsplit("/", 1)[-1])
            if match:
                key = "_".join(match.groupdict().values())
        if not key:
            key = f"{int(b[0]) // 1000}_{int(b[1]) // 1000}"
        items.append(
            Item(
                bounds=b,
                url=href,
                date=_parse_date(props.get(date_field), date_format) if date_field else None,
                key=key,
            )
        )
    if skipped:
        logger.debug("Index %s: skipped %d features without url or geometry", name, skipped)
    if date_field and items and all(it.date is None for it in items):
        # Every item then drops out of the epoch selection, and the run ends in
        # "no coverage" -- which points at the area instead of at the format.
        logger.warning(
            "Index %s: no value of '%s' parsed as a date%s. Set `date_format` "
            "(strptime syntax) in the region file.", name, date_field,
            f" with date_format '{date_format}'" if date_format else " (expected YYYY-MM-DD)",
        )
    cache.release(raw)

    condensed.parent.mkdir(parents=True, exist_ok=True)
    with condensed.open("w", encoding="utf-8") as fh:
        json.dump(
            {
                "name": name,
                "crs": f"EPSG:{transform_to_epsg}" if transform_to_epsg else index_crs,
                "items": [
                    {
                        "b": list(it.bounds),
                        "u": it.url,
                        "d": it.date.isoformat() if it.date else None,
                        "k": it.key,
                    }
                    for it in items
                ],
            },
            fh,
        )
    logger.info("Index %s: %d tiles", name, len(items))
    return Index(name, items, f"EPSG:{transform_to_epsg}" if transform_to_epsg else index_crs)


def _condensed_matches(path: Path, epsg: Optional[int]) -> bool:
    """Whether a cached condensed index is in the CRS a run asks for."""
    if epsg is None:
        return True
    try:
        with path.open(encoding="utf-8") as fh:
            return json.load(fh).get("crs") == f"EPSG:{epsg}"
    except (OSError, ValueError):
        return False


def _read_condensed(name: str, path: Path) -> Index:
    with path.open(encoding="utf-8") as fh:
        data = json.load(fh)
    items = [
        Item(
            bounds=tuple(row["b"]),  # type: ignore[arg-type]
            url=row["u"],
            date=_parse_date(row["d"]),
            key=row["k"],
        )
        for row in data["items"]
    ]
    return Index(name, items, data.get("crs"))


# ---------------------------------------------------------------------------
# epoch selection
# ---------------------------------------------------------------------------
@dataclass
class EpochSpec:
    """Selection rule for one orthophoto role (e.g. leaf-on)."""

    name: str
    months: Sequence[int]
    min_year: int = 0
    w_season: float = 1.0
    w_recency: float = 0.35
    required: bool = True

    @classmethod
    def from_dict(cls, name: str, raw: Dict[str, Any]) -> "EpochSpec":
        return cls(
            name=name,
            months=[int(m) for m in (raw.get("months") or [])],
            min_year=int(raw.get("min_year", 0)),
            w_season=float(raw.get("w_season", 1.0)),
            w_recency=float(raw.get("w_recency", 0.35)),
            required=bool(raw.get("required", True)),
        )


def _score(item: Item, spec: EpochSpec, year_lo: int, year_hi: int) -> float:
    season = 1.0 if item.month in spec.months else 0.0
    span = max(1, year_hi - year_lo)
    recency = min(1.0, max(0.0, (item.year - year_lo) / span))
    return spec.w_season * season + spec.w_recency * recency


@dataclass
class EpochChoice:
    """Result of epoch selection for one work tile."""

    item: Optional[Item]
    relaxed: str = ""


def select_epochs(
    index: Index,
    tile_bounds: Dict[str, Bounds],
    spec: EpochSpec,
    prefer_uniform: bool = True,
    max_spread_days: int = 30,
) -> Tuple[Dict[str, EpochChoice], Dict[str, Any]]:
    """Choose one source item per work tile for the given epoch role.

    Ortho tiles are typically larger than a work tile (2 km in Lower Saxony),
    so choosing per tile independently imprints the acquisition lattice onto
    the NDVI field as visible steps. When
    ``prefer_uniform`` is set, the selector first looks for a single acquisition
    campaign (a cluster of dates no wider than ``max_spread_days``) that covers
    every tile, and only falls back to per-tile choice if none does.

    Returns:
        ``(choices, report)`` where ``report`` carries the diagnostics that
        ``palmgeo plan`` prints.
    """
    candidates: Dict[str, List[Item]] = {}
    for key, b in tile_bounds.items():
        hits = [it for it in index.search(b) if it.date is not None]
        candidates[key] = hits

    years = [it.year for hits in candidates.values() for it in hits if it.year]
    year_lo, year_hi = (min(years), max(years)) if years else (0, 1)

    def filtered(strictness: int) -> Dict[str, List[Item]]:
        """strictness 2 = season + min_year, 1 = min_year only, 0 = anything."""
        out: Dict[str, List[Item]] = {}
        for key, hits in candidates.items():
            sel = hits
            if strictness >= 1:
                sel = [it for it in sel if it.year >= spec.min_year]
            if strictness >= 2 and spec.months:
                sel = [it for it in sel if it.month in spec.months]
            out[key] = sel
        return out

    relaxed = ""
    pool = filtered(2)
    if any(not v for v in pool.values()):
        missing = sum(1 for v in pool.values() if not v)
        pool_relaxed = filtered(1)
        if all(pool_relaxed.values()):
            relaxed = f"season constraint dropped for {missing} tile(s)"
            pool = {k: (v or pool_relaxed[k]) for k, v in pool.items()}
        else:
            pool_any = filtered(0)
            relaxed = f"season and min_year constraints dropped for {missing} tile(s)"
            pool = {k: (v or pool_relaxed[k] or pool_any[k]) for k, v in pool.items()}

    report: Dict[str, Any] = {"role": spec.name, "relaxed": relaxed, "uniform": None}

    chosen: Dict[str, EpochChoice] = {}
    if prefer_uniform:
        all_dates = sorted({it.date for hits in pool.values() for it in hits if it.date})
        best: Tuple[float, Optional[date]] = (-1.0, None)
        for anchor in all_dates:
            window = [
                d for d in all_dates if abs((d - anchor).days) <= max_spread_days
            ]
            wset = set(window)
            covered = all(
                any(it.date in wset for it in hits) for hits in pool.values() if hits
            )
            if not covered:
                continue
            picks = [
                max(
                    (it for it in hits if it.date in wset),
                    key=lambda it: _score(it, spec, year_lo, year_hi),
                )
                for hits in pool.values()
                if hits
            ]
            mean = sum(_score(it, spec, year_lo, year_hi) for it in picks) / max(1, len(picks))
            if mean > best[0]:
                best = (mean, anchor)
        if best[1] is not None:
            anchor = best[1]
            wset = {d for d in all_dates if abs((d - anchor).days) <= max_spread_days}
            report["uniform"] = {
                "anchor": anchor.isoformat(),
                "spread_days": max_spread_days,
                "mean_score": round(best[0], 3),
            }
            for key, hits in pool.items():
                sel = [it for it in hits if it.date in wset]
                chosen[key] = EpochChoice(
                    max(sel, key=lambda it: _score(it, spec, year_lo, year_hi))
                    if sel
                    else None,
                    relaxed,
                )
            return chosen, report

    for key, hits in pool.items():
        chosen[key] = EpochChoice(
            max(hits, key=lambda it: _score(it, spec, year_lo, year_hi)) if hits else None,
            relaxed,
        )
    return chosen, report


def pick_paired(
    partner_index: Index, reference: Item, pair_on_date: bool = True
) -> Optional[Item]:
    """Find the partner tile matching ``reference`` in extent and, if asked, date.

    DGM1 and DOM1 share tile identifiers *and* acquisition dates for all 49 708
    statewide tiles, so pairing on the date yields a temporally consistent nDSM
    (``DOM1 - DGM1``). Image-based surface models such as bDOM20 have no such
    pairing, hence ``pair_on_date=False`` for those.
    """
    hits = [
        it
        for it in partner_index.search(reference.bounds)
        # require essentially the same footprint, not just any overlap
        if abs(it.bounds[0] - reference.bounds[0]) < 1.0
        and abs(it.bounds[1] - reference.bounds[1]) < 1.0
    ]
    if not hits:
        return None
    if pair_on_date and reference.date is not None:
        same = [it for it in hits if it.date == reference.date]
        if same:
            return same[0]
        return min(
            hits,
            key=lambda it: abs((it.date - reference.date).days) if it.date else 10**6,
        )
    return max(hits, key=lambda it: it.date or date.min)


def newest_per_tile(index: Index, tile_bounds: Dict[str, Bounds]) -> Dict[str, Optional[Item]]:
    """Pick the most recent item covering each work tile."""
    out: Dict[str, Optional[Item]] = {}
    for key, b in tile_bounds.items():
        hits = index.search(b)
        out[key] = max(hits, key=lambda it: it.date or date.min) if hits else None
    return out


def regions_for_bounds(index: Index, b: Bounds) -> List[Item]:
    """Administrative units (e.g. per-county cadastre extracts) intersecting the AOI."""
    return index.search(b)


def iter_urls(items: Iterable[Optional[Item]]) -> Iterator[str]:
    for item in items:
        if item is not None:
            yield item.url
