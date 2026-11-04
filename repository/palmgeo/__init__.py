"""palmgeo -- open geodata to palm_csd static driver inputs.

A staged, resumable pipeline that turns publicly available open geodata --
national, regional or municipal -- into the georeferenced raster and vector
inputs of PALM's `palm_csd`, plus the matching palm_csd YAML configuration.

Region specifics live entirely in two YAML files (`sources_*.yml`,
`mapping_*.yml`); the Python code addresses sources by role only.
"""

from __future__ import annotations

import os
import sqlite3
import warnings
from pathlib import Path
from typing import Optional, Tuple

__version__ = "1.2.0"

#: Where palmgeo is published.
PROJECT_URL = "https://github.com/meteohauschopp/palmgeo"

#: Sent on every request, by the pipeline and by the staging scripts under
#: repository/ports/. Some open-data portals rate-limit or block clients that
#: do not identify themselves, and an operator who sees this in a log should be
#: able to tell which tool and which version is calling. Built from the package
#: version rather than written out, because a hard-coded version silently
#: drifts and then misinforms exactly the person it is meant to inform -- it
#: read "palmgeo/1.0" throughout 0.1 and 0.2. A project URL rather than a
#: personal address: this file is published, and a mail address in it is
#: harvested.
USER_AGENT = f"palmgeo/{__version__} (PALM static driver preprocessing; +{PROJECT_URL})"


def _proj_db_layout(db: Path) -> Optional[Tuple[int, int]]:
    """Read the ``DATABASE.LAYOUT.VERSION`` of a PROJ database."""
    try:
        with sqlite3.connect(f"file:{db}?mode=ro", uri=True) as conn:
            rows = dict(
                conn.execute(
                    "SELECT key, value FROM metadata "
                    "WHERE key LIKE 'DATABASE.LAYOUT.VERSION.%'"
                ).fetchall()
            )
        return (
            int(rows["DATABASE.LAYOUT.VERSION.MAJOR"]),
            int(rows["DATABASE.LAYOUT.VERSION.MINOR"]),
        )
    except Exception:  # noqa: BLE001 - a missing or unreadable db is "unknown"
        return None


def _fix_proj_env() -> None:
    """Drop a ``PROJ_DATA``/``PROJ_LIB`` that is older than rasterio's own copy.

    A stale ``PROJ_DATA`` left over from another PROJ installation is a common
    and badly diagnosable failure: GDAL then reads a ``proj.db`` with an older
    schema than the PROJ it was linked against, and every CRS operation dies
    with "The WKT could not be parsed. OGR Error code 6" -- which reads like a
    data problem but is purely environmental.

    This has to run before rasterio is imported, hence its place in the package
    ``__init__``.
    """
    candidates = [v for v in ("PROJ_DATA", "PROJ_LIB") if os.environ.get(v)]
    if not candidates:
        return
    bundled: Optional[Path] = None
    try:
        import rasterio

        path = Path(rasterio.__file__).parent / "proj_data" / "proj.db"
        bundled = path if path.is_file() else None
    except Exception:  # noqa: BLE001
        bundled = None
    if bundled is None:
        try:
            import pyproj

            bundled = Path(pyproj.datadir.get_data_dir()) / "proj.db"
        except Exception:  # noqa: BLE001
            return
    reference = _proj_db_layout(bundled)
    if reference is None:
        return
    for var in candidates:
        external = _proj_db_layout(Path(os.environ[var]) / "proj.db")
        if external is not None and external < reference:
            warnings.warn(
                f"{var}={os.environ[var]} points at a PROJ database with layout "
                f"{external[0]}.{external[1]}, older than the "
                f"{reference[0]}.{reference[1]} this installation expects. "
                f"Ignoring it for this process; unset {var} in your shell to "
                f"silence this.",
                RuntimeWarning,
                stacklevel=2,
            )
            del os.environ[var]


_fix_proj_env()
