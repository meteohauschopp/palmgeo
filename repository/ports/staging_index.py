"""Write a source index in the shape palmgeo reads.

A *port* is what it takes to build drivers for a region whose data palmgeo
cannot fetch directly. Some regions publish a download index per product and
need no port at all; others publish an archive per state, or a STAC API, or a
tile server with its own tiling and its own coordinate reference system. A
staging script absorbs that difference: it fetches, reprojects, mosaics and
fuses whatever the region offers into files on disk, then writes an index that
points at them.

The index is the contract between a port and the pipeline, and it is small:

    {"type": "FeatureCollection",
     "features": [{"type": "Feature",
                   "geometry": <polygon, in the run's CRS>,
                   "properties": {"tile_id": ..., "staged_path": ..., ...}}]}

One feature per staged file, its geometry covering the ground that file covers.
``sources_<region>.yml`` then names which property holds the path
(``url_field``) and which the acquisition date (``date_field``).

A single feature covering the whole area of interest is a legitimate degenerate
case, and the usual one for a region that ships a single mosaic rather than
tiles: ``sources.py`` imposes no minimum, and every work tile falls inside that
one feature, so it is always the match.

WHY THIS IS SHARED
------------------
The three staging scripts written so far each wrote this structure by hand,
and each got the same two things wrong at least once:

*The index named a file that did not exist.* The New York port emitted
``..._dem_epsg6347.tif`` while the staging step wrote
``..._dtm_epsg6347_m.tif``. Nothing complained until the raster stage
tried to open it, hours later.

*The index kept working while the paths under it moved.* Absolute paths written
into an index survive a reorganisation of the project silently, and fail only
when something finally opens one.

Both are the same failure: an index is a promise about the filesystem that
nothing checks at the time it is made. ``write_index`` checks it, at the time
it is made, and refuses to write a promise it cannot keep.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

logger = logging.getLogger("palmgeo.ports")

Bounds = Tuple[float, float, float, float]


@dataclass
class StagedFile:
    """One file a port has produced, and the ground it covers."""

    #: Stable identifier. Used to pair products across sources (a terrain tile
    #: with the surface tile above it), so it must be equal for files that
    #: describe the same ground and different otherwise.
    tile_id: str
    #: The file itself. Absolute; see write_index for why it is checked.
    path: Path
    #: Extent in the RUN's coordinate reference system, not the source's.
    #: Reprojecting the geometry is the port's job, because only the port
    #: knows what the source was in.
    bounds: Bounds
    #: Acquisition date, ISO 8601. None where the source states none -- which
    #: is honest, and better than inventing one: palmgeo's epoch selection
    #: skips dateless items rather than ranking them wrongly.
    date: Optional[str] = None
    #: Anything else worth carrying; written into the feature's properties.
    extra: Dict[str, Any] = field(default_factory=dict)


def bbox_polygon(bounds: Bounds) -> Dict[str, Any]:
    """A closed GeoJSON polygon ring from a bounding box."""
    xmin, ymin, xmax, ymax = bounds
    return {
        "type": "Polygon",
        "coordinates": [[
            [xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax], [xmin, ymin],
        ]],
    }


def reproject_bounds(bounds: Bounds, src_epsg: int, dst_epsg: int) -> Bounds:
    """Transform a bounding box between coordinate reference systems.

    The corners are transformed and re-bounded rather than assumed to stay
    corners: away from the projection's centre a rectangle does not map to a
    rectangle, and taking only two opposite corners understates the extent.
    All four are used, which is exact enough for an index whose purpose is to
    answer "does this file overlap that tile".
    """
    if src_epsg == dst_epsg:
        return bounds
    from pyproj import Transformer

    tf = Transformer.from_crs(src_epsg, dst_epsg, always_xy=True)
    xmin, ymin, xmax, ymax = bounds
    xs, ys = tf.transform(
        [xmin, xmax, xmin, xmax], [ymin, ymin, ymax, ymax]
    )
    return (min(xs), min(ys), max(xs), max(ys))


def aoi_bounds(bbox_wgs84: Bounds, dst_epsg: int) -> Bounds:
    """The staged area in the run's CRS: a WGS84 box, snapped out to whole metres.

    Every raster a port warps is written onto exactly this extent (gdalwarp's
    ``-te``), and every single-feature index takes its geometry from it too --
    so the index describes what is staged, not what was ordered (PORTING.md,
    "A staged index must describe what is staged"). Two scripts staging into
    the same work directory get identical grids by passing the same box.
    """
    import math

    xmin, ymin, xmax, ymax = reproject_bounds(bbox_wgs84, 4326, dst_epsg)
    return (float(math.floor(xmin)), float(math.floor(ymin)),
            float(math.ceil(xmax)), float(math.ceil(ymax)))


def write_index(
    path: Path,
    files: Sequence[StagedFile],
    *,
    root: Optional[Path] = None,
    date_property: str = "acquired",
    require_files: bool = True,
    indent: Optional[int] = 2,
) -> Path:
    """Write ``files`` as a palmgeo source index, refusing broken promises.

    ``require_files`` is on by default and is the reason this function exists.
    An index entry is a claim that a particular file is on disk; checking the
    claim costs a stat call and turns a failure that would otherwise surface
    hours later, inside a worker thread, into one that surfaces here, naming
    every missing file at once.

    Turn it off only to write an index ahead of the data it describes, and
    expect to be told so in the log.
    """
    if not files:
        raise ValueError(f"{path}: refusing to write an index with no features")

    missing = [f for f in files if require_files and not f.path.exists()]
    if missing:
        raise FileNotFoundError(
            f"{path}: {len(missing)} staged file(s) named by this index do not "
            f"exist. An index that points at a missing file fails much later "
            f"and much less clearly:\n  "
            + "\n  ".join(str(f.path) for f in missing[:10])
        )
    if not require_files:
        logger.warning(
            "%s: written without checking that the staged files exist", path.name
        )

    features: List[Dict[str, Any]] = []
    for f in files:
        # Relative to the project root wherever possible. An index is written
        # once and read for as long as the project lives, so it outlasts the
        # directory layout it was written under: absolute paths inside indices
        # are what broke every staged region when the project was reorganised,
        # silently, until a worker thread tried to open one. The cache resolves
        # a relative path against the project root (see cache.Cache.fetch).
        staged = f.path.resolve()
        if root is not None:
            try:
                staged_value = str(staged.relative_to(Path(root).resolve()))
            except ValueError:
                logger.warning(
                    "%s lies outside the project root, so this index has to name "
                    "it absolutely and will not survive a move of the project.",
                    staged,
                )
                staged_value = str(staged)
        else:
            staged_value = str(staged)
        properties: Dict[str, Any] = {
            "tile_id": f.tile_id,
            "staged_path": staged_value,
        }
        if f.date:
            properties[date_property] = f.date
        properties.update(f.extra)
        features.append({
            "type": "Feature",
            "geometry": bbox_polygon(f.bounds),
            "properties": properties,
        })

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({"type": "FeatureCollection", "features": features}, indent=indent),
        encoding="utf-8",
    )
    logger.info("%s: %d feature(s)", path, len(features))
    return path


def repair_paths(index_path: Path, old_root: Path, new_root: Path) -> int:
    """Rewrite absolute ``staged_path`` values after the data has moved.

    Exists because reorganising a project moves the staged files while the
    indices keep pointing at where they used to be, and nothing notices until
    a stage tries to open one. Returns the number of features changed.
    """
    data = json.loads(index_path.read_text(encoding="utf-8"))
    old, new, changed = str(old_root), str(new_root), 0
    for feature in data.get("features", []):
        props = feature.get("properties") or {}
        value = props.get("staged_path")
        if isinstance(value, str) and value.startswith(old):
            props["staged_path"] = new + value[len(old):]
            changed += 1
    if changed:
        index_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
        logger.info("%s: repointed %d path(s)", index_path.name, changed)
    return changed
