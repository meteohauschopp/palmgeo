"""Bounded download cache.

The pipeline is designed to run on far less disk than the data it touches: a
48 x 48 km area of interest pulls roughly 70 GB of 20 cm orthophotos per epoch.
Source files are therefore treated as *transient* -- fetched, consumed, and
dropped again -- with a hard cap on how much may sit on disk at any moment.

Eviction is LRU over unreferenced entries. A file that a worker is currently
reading is pinned via :meth:`Cache.hold` and can never be evicted underneath it.

Ownership
---------
The cache directory holds more than downloads. It also holds the condensed
source indices (rebuildable, but only by re-fetching and re-parsing a GeoJSON
of up to 78 MB) and, for some regions, data a human put there by hand -- the
global coastline water polygons, for instance, which no URL in any region file
can restore.

**The cache therefore evicts only what it downloaded itself.** Ownership is
recorded in a manifest beside the data, so it survives a restart; anything not
in the manifest is invisible to the eviction logic no matter where it sits.
Deciding ownership by directory convention was the earlier design and was
wrong: it made the recovery cost of an eviction depend on a path prefix, and
the one file the cache could not restore lived under the same root as the ones
it could.

The manifest stays small in practice because entries are transient by design:
a stage fetches a tile, reads it, and calls :meth:`release`, which both deletes
the file and drops it from the manifest. Only a run with ``keep_downloads``
accumulates.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import threading
import time
import urllib.error
import urllib.request
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional

logger = logging.getLogger("palmgeo.cache")

#: Sent on every request; see palmgeo.USER_AGENT for why.
from palmgeo import USER_AGENT  # noqa: E402  (kept next to its use)

#: File name of the ownership manifest, inside the cache directory.
MANIFEST = ".palmgeo-cache.json"


@dataclass
class _Entry:
    path: Path
    size: int
    last_used: float = field(default_factory=time.time)
    holds: int = 0


class DownloadError(RuntimeError):
    """A source file could not be retrieved after all retries."""


class OfflineError(RuntimeError):
    """A file would have to be fetched, and the run was told not to.

    Separate from DownloadError because the two call for opposite responses:
    a download failure is usually worth retrying, and this never is.
    """


class Cache:
    """Size-capped, thread-safe download cache."""

    def __init__(
        self,
        directory: Path,
        budget_bytes: int,
        keep: bool = False,
        retries: int = 4,
        timeout: float = 180.0,
        root: Optional[Path] = None,
        offline: bool = False,
    ) -> None:
        self.dir = Path(directory)
        #: Project root, used to resolve the relative local paths that a
        #: staged source index carries. Optional so that the cache stays
        #: usable on its own in a test; without it a relative path is an
        #: error rather than a guess against the working directory, which
        #: would depend on where the caller happened to be standing.
        self.root = Path(root) if root is not None else None
        #: Refuse the network entirely. Not a retry policy -- the point is to
        #: fail at once, with a name, instead of spending four backoff rounds
        #: per file discovering the same thing a few thousand times.
        self.offline = bool(offline)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.budget = int(budget_bytes)
        self.keep = bool(keep)
        self.retries = int(retries)
        self.timeout = float(timeout)
        self._entries: Dict[Path, _Entry] = {}
        self._lock = threading.RLock()
        self._manifest = self.dir / MANIFEST
        # One lock per target path. Tile workers run concurrently and a 2 km
        # orthophoto tile serves four 1 km work tiles, so the same URL is
        # routinely requested by several threads at once; without this they
        # would write the same .part file simultaneously.
        self._downloads: Dict[Path, threading.Lock] = {}
        self.bytes_downloaded = 0
        self._adopt_existing()

    def _download_lock(self, target: Path) -> threading.Lock:
        """The per-target lock, created on first use and dropped when unused.

        Without pruning this dictionary grows one lock per distinct URL for the
        life of the process: a city-scale run touches tens of thousands of
        tiles, and every one of them would leave a lock object behind long
        after its download finished.
        """
        with self._lock:
            lock = self._downloads.get(target)
            if lock is None:
                lock = self._downloads[target] = threading.Lock()
            return lock

    def _release_download_lock(self, target: Path) -> None:
        """Forget a target's lock once nobody is waiting on it.

        Checked with ``locked()`` rather than by reference counting: a lock
        that is free at this moment has no waiter that could be surprised by
        its disappearance, because any later caller creates a fresh one under
        ``self._lock`` before touching the target at all.
        """
        with self._lock:
            lock = self._downloads.get(target)
            if lock is not None and not lock.locked():
                del self._downloads[target]

    # ------------------------------------------------------------- internals --
    def _adopt_existing(self) -> None:
        """Re-register what an earlier run downloaded, and nothing else.

        Reusing a tile that is already on disk is the difference between a
        resumed run and a repeated one, so the entries have to survive a
        restart. What must NOT survive is the assumption that everything under
        the cache directory belongs to the cache: the condensed indices and any
        hand-placed regional data live here too, and adopting them would make
        them eviction candidates from the second run onward -- silently
        deleting, under a tight budget, the one file that cannot be downloaded
        again.

        A manifest entry whose file has since disappeared is dropped rather
        than recreated; that is the normal state after somebody has cleared
        the directory by hand.
        """
        claimed = self._read_manifest()
        for rel in claimed:
            path = self.dir / rel
            if not path.exists() or path.name.endswith(".part"):
                continue
            size = (
                sum(f.stat().st_size for f in path.rglob("*") if f.is_file())
                if path.is_dir() else path.stat().st_size
            )
            self._entries[path] = _Entry(path, size, path.stat().st_mtime)
        if len(self._entries) != len(claimed):
            self._write_manifest()          # prune vanished entries

    def _read_manifest(self) -> List[str]:
        """Paths, relative to the cache directory, that the cache owns."""
        if not self._manifest.is_file():
            return []
        try:
            data = json.loads(self._manifest.read_text(encoding="utf-8"))
            return [str(p) for p in data.get("owned", [])]
        except (OSError, ValueError) as exc:
            # A corrupt manifest must not stop a run. The cost of ignoring it
            # is re-downloading files that are already on disk; the cost of
            # trusting it blindly would be evicting files that are not ours.
            logger.warning("Ignoring unreadable cache manifest %s: %s", self._manifest, exc)
            return []

    def _write_manifest(self) -> None:
        """Persist the ownership list. Caller must hold ``self._lock``."""
        try:
            payload = {
                "note": "Files this cache downloaded and may delete. "
                        "Anything under the cache directory but absent here is "
                        "left alone.",
                "owned": sorted(
                    str(p.relative_to(self.dir)) for p in self._entries
                    if self.dir in p.parents or p.parent == self.dir
                ),
            }
            tmp = self._manifest.with_suffix(".tmp")
            tmp.write_text(json.dumps(payload, indent=1), encoding="utf-8")
            os.replace(tmp, self._manifest)
        except OSError as exc:  # pragma: no cover - best effort
            logger.debug("Could not write cache manifest: %s", exc)

    def _total(self) -> int:
        return sum(e.size for e in self._entries.values())

    def _evict(self, need: int) -> None:
        """Drop least-recently-used unheld entries until ``need`` bytes fit."""
        if self.keep or self.budget <= 0:
            return
        with self._lock:
            total = self._total()
            if total + need <= self.budget:
                return
            candidates = sorted(
                (e for e in self._entries.values() if e.holds == 0),
                key=lambda e: e.last_used,
            )
            for entry in candidates:
                if total + need <= self.budget:
                    break
                total -= entry.size
                self._drop(entry.path)
            if total + need > self.budget:
                logger.warning(
                    "Cache budget of %.1f GB is too small for the working set "
                    "(%.1f GB held). Raise runtime.profiles.<profile>.cache_budget_gb.",
                    self.budget / 1e9,
                    total / 1e9,
                )

    def _drop(self, path: Path) -> None:
        """Delete an owned entry and forget it. Caller must hold ``self._lock``."""
        entry = self._entries.pop(path, None)
        if entry is not None:
            self._write_manifest()
        try:
            if path.is_dir():
                shutil.rmtree(path, ignore_errors=True)
            elif path.exists():
                path.unlink()
        except OSError as exc:  # pragma: no cover - best effort cleanup
            logger.debug("Could not remove %s: %s", path, exc)
        if entry is not None:
            logger.debug("Evicted %s (%.1f MB)", path.name, entry.size / 1e6)

    # ---------------------------------------------------------------- public --
    def fetch(self, url: str, subdir: str = "", filename: Optional[str] = None) -> Path:
        """Return a local path for ``url``, downloading it if necessary.

        A source that is already a local file is used IN PLACE -- no copy, no
        cache accounting, no eviction. Region ports that pre-process their
        sources (reprojecting into the run's CRS, fusing several datasets into
        one GeoPackage -- see repository/ports/us_ny/stage_sources.py) point
        their index `url_field` at staged files on disk. Without this branch
        those paths reach urllib as a schemeless URL and raise
        "unknown url type", and the `file://` workaround would copy the whole
        staged set (~24 GB for the Manhattan ortho tiles) into the cache for
        no benefit.
        """
        if url.startswith("file://"):
            url = url[len("file://"):]
        if url.startswith("/"):
            local = Path(url)
            if local.is_file():
                return local
            raise FileNotFoundError(
                f"staged source {local} does not exist -- run the region's "
                f"staging script before this stage"
            )
        if "://" not in url:
            # A relative local path, which is how a staged source index names
            # its files: relative to the project root, so that moving or
            # handing on a project does not break every index inside it.
            # Resolved here rather than by the caller because this is the one
            # place every local source passes through.
            if self.root is not None:
                local = (self.root / url).resolve()
                if local.is_file():
                    return local
                raise FileNotFoundError(
                    f"staged source {local} does not exist (named as '{url}', "
                    f"relative to the project root). Run the region's staging "
                    f"script, or check that the index still matches what it "
                    f"produced."
                )
            raise ValueError(
                f"'{url}' is neither a URL nor an absolute path, and this cache "
                f"was created without a project root to resolve it against. "
                f"A relative path resolved against the working directory would "
                f"depend on where the pipeline happened to be started."
            )
        name = filename or url.rsplit("/", 1)[-1].split("?")[0]
        target = self.dir / subdir / name
        with self._lock:
            entry = self._entries.get(target)
            if entry is not None and target.exists():
                entry.last_used = time.time()
                return target
        target.parent.mkdir(parents=True, exist_ok=True)

        try:
            with self._download_lock(target):
                # Another thread may have finished this download while we
                # waited on the lock, so the check has to be repeated inside
                # it -- the one outside only saves the common uncontended case.
                with self._lock:
                    entry = self._entries.get(target)
                    if entry is not None and target.exists():
                        entry.last_used = time.time()
                        return target
                return self._download(url, target, name)
        finally:
            self._release_download_lock(target)

    def cached_path(self, url: str, subdir: str = "",
                    filename: Optional[str] = None) -> Optional[Path]:
        """The local copy of ``url`` if there is one, without fetching.

        Lets a caller that would otherwise stream check whether the whole file
        is already here. Deliberately does not touch ``last_used``: asking is
        not using, and an entry should not be protected from eviction by
        somebody who then decided to stream instead.
        """
        if url.startswith("file://"):
            url = url[len("file://"):]
        if url.startswith("/"):
            local = Path(url)
            return local if local.is_file() else None
        if "://" not in url:
            if self.root is None:
                return None
            local = (self.root / url).resolve()
            return local if local.is_file() else None
        name = filename or url.rsplit("/", 1)[-1].split("?")[0]
        target = self.dir / subdir / name
        return target if target.is_file() else None

    def _download(self, url: str, target: Path, name: str) -> Path:
        if self.offline:
            raise OfflineError(
                f"{name} is not on disk and this run is offline "
                f"(runtime.offline). It would have come from {url}"
            )
        self._evict(need=0)
        tmp = target.with_suffix(target.suffix + ".part")
        last_exc: Exception | None = None
        for attempt in range(1, self.retries + 1):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    declared = response.headers.get("Content-Length")
                    if declared:
                        self._evict(need=int(declared))
                    with tmp.open("wb") as fh:
                        shutil.copyfileobj(response, fh, length=1 << 20)
                size = tmp.stat().st_size
                if declared and size != int(declared):
                    raise DownloadError(
                        f"truncated download: got {size} of {declared} bytes"
                    )
                os.replace(tmp, target)
                with self._lock:
                    self._entries[target] = _Entry(target, size)
                    self.bytes_downloaded += size
                    # Ownership is recorded at the moment the file becomes
                    # complete, not when the request starts: a manifest entry
                    # for a half-written .part would license a later run to
                    # delete something it never verified.
                    self._write_manifest()
                logger.debug("Fetched %s (%.1f MB)", name, size / 1e6)
                return target
            except (urllib.error.URLError, OSError, DownloadError) as exc:
                last_exc = exc
                tmp.unlink(missing_ok=True)
                if attempt < self.retries:
                    delay = min(2.0**attempt, 30.0)
                    logger.warning(
                        "Download failed (%s/%s) for %s: %s -- retrying in %.0f s",
                        attempt, self.retries, name, exc, delay,
                    )
                    time.sleep(delay)
        raise DownloadError(f"Could not download {url}: {last_exc}")

    @contextmanager
    def hold(self, path: Path) -> Iterator[Path]:
        """Pin ``path`` for the duration of the block so it survives eviction."""
        with self._lock:
            entry = self._entries.get(path)
            if entry is not None:
                entry.holds += 1
                entry.last_used = time.time()
        try:
            yield path
        finally:
            with self._lock:
                entry = self._entries.get(path)
                if entry is not None:
                    entry.holds = max(0, entry.holds - 1)

    def discard(self, path: Path) -> None:
        """Remove an entry unconditionally, keeping the accounting straight.

        Used for extraction directories, which are always transient regardless
        of ``keep_downloads``: the clipped result has already been written.
        """
        with self._lock:
            self._drop(path)

    def release(self, path: Path) -> None:
        """Delete ``path`` unless the profile asks for downloads to be kept.

        Only cache-owned entries (i.e. files that were actually downloaded and
        registered in ``self._entries``) are deleted. In-place staged source
        files that ``fetch`` handed straight back (the local-file fast path)
        are never registered and are therefore left untouched -- the fetch
        contract says they are used in place with "no cache accounting, no
        eviction", and a stage boundary must not destroy the region's staged
        sources.
        """
        if self.keep:
            return
        with self._lock:
            entry = self._entries.get(path)
            if entry is None:
                return
            if entry.holds > 0:
                return
            self._drop(path)

    def unzip(self, archive: Path, member_suffixes: tuple[str, ...] = ()) -> Path:
        """Extract ``archive`` next to itself and return the extraction directory.

        Args:
            archive: The ``.zip`` to extract.
            member_suffixes: If given, only members whose name ends with one of
                these suffixes are extracted (e.g. ``(".shp", ".dbf", ".shx",
                ".prj", ".cpg")`` to skip unrelated layers of a large archive).
        """
        out = archive.with_suffix("")
        marker = out / ".complete"
        if marker.exists():
            return out
        out.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(archive) as zf:
            members = [
                m
                for m in zf.namelist()
                if not m.endswith("/")
                and (not member_suffixes or m.lower().endswith(member_suffixes))
                # refuse absolute paths and traversal
                and not m.startswith(("/", "\\"))
                and ".." not in Path(m).parts
            ]
            zf.extractall(out, members=members)
        marker.touch()
        with self._lock:
            size = sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
            self._entries[out] = _Entry(out, size)
            self._write_manifest()
        return out

    def summary(self) -> str:
        with self._lock:
            return (
                f"cache: {len(self._entries)} entries, {self._total() / 1e9:.2f} GB on "
                f"disk, {self.bytes_downloaded / 1e9:.2f} GB downloaded this run "
                f"(budget {self.budget / 1e9:.0f} GB, keep={self.keep})"
            )
