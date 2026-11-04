"""Domain, resolution-tier and work-tile geometry.

Two invariants drive this module:

1.  **PALM's ``nx`` convention.** ``nx`` is the largest grid index, so a domain
    raster has ``nx + 1`` columns and ``ny + 1`` rows.

2.  **palm_csd's memory behaviour.** ``palm_csd/geo_converter.py`` only takes
    the windowed read path (``cut_aligned_raster_to_dst``) when the input CRS,
    the pixel size *and* the grid origin match the target domain to within
    0.2 px. Any mismatch falls through to ``reproject_raster_to_dst``, which
    calls ``src.read(masked=True)`` and pulls the **entire** raster into RAM.
    Every raster this pipeline hands to palm_csd is therefore written on
    exactly the domain grid.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Sequence, Tuple

from rasterio.transform import Affine

Bounds = Tuple[float, float, float, float]  # (xmin, ymin, xmax, ymax)

#: Side length of a work tile in metres.
#:
#: A work tile is palmgeo's own unit of work, not a source tile: it is the
#: square that one worker thread reads, derives and writes independently of
#: every other, which is what makes the pipeline parallel and resumable.
#:
#: 1 km was chosen because it coincides with the tiling of the reference
#: region's elevation products (Lower Saxony's DGM1/DOM1/LoD2 are 1 km, its
#: orthophotos 2 km), so one work tile maps onto exactly one source tile and
#: four work tiles onto one ortho tile -- no source file is opened twice for
#: the same work, and none is opened for a fraction of its content.
#:
#: It is NOT a claim about other regions. Manhattan and Aswan are served from
#: staged mosaics with no tiling of their own, and there 1 km simply sets the
#: memory ceiling: at 1 m work resolution a tile is a 1000x1000 array per
#: product, so a worker's peak footprint is tens of megabytes however large
#: the area of interest is. That property is the reason to keep the value
#: small, and the reason not to raise it for a region with coarser sources.
TILE_SIZE = 1000.0


def multigrid_leftover(points: int) -> int:
    """The part of a point count PALM's pressure solver cannot coarsen away.

    PALM's multigrid solver halves the grid repeatedly. Each halving needs an
    even count, so the coarsening stops at the first odd factor other than the
    3 that a final 3x3 level tolerates. This returns what is left after
    dividing out every 2 and 3 -- 1 for a count that coarsens all the way,
    and the obstructing factor otherwise.

        1536 = 2^9 * 3      -> 1     nine levels
        1537 = 29 * 53      -> 1537  no coarsening at all
        1500 = 2^2 * 3 * 5^3-> 125   two levels

    One point off is not a rounding error; it is the difference between a fast
    solver and an unusable one.
    """
    left = int(points)
    for factor in (2, 3):
        while left % factor == 0:
            left //= factor
    return left


def nearest_smooth(points: int, span: int = 64) -> List[int]:
    """Point counts near ``points`` that the multigrid solver can coarsen.

    Offered in error and warning messages so that the reader is not left to
    factorise numbers by hand. Searched outward from ``points`` so that the
    closest usable counts come first, which is what somebody who has already
    chosen an extent actually wants.
    """
    out: List[int] = []
    for delta in range(0, span + 1):
        for candidate in ({points - delta, points + delta} if delta else {points}):
            if candidate > 0 and multigrid_leftover(candidate) == 1:
                out.append(candidate)
        if len(out) >= 3:
            break
    return sorted(set(out), key=lambda c: (abs(c - points), c))[:3]


@dataclass(frozen=True)
class Domain:
    """One PALM model domain."""

    name: str
    pixel_size: float
    nx: int
    ny: int
    dz: float
    origin_x: float
    origin_y: float
    domain_parent: str | None = None
    options: Dict[str, object] = field(default_factory=dict)

    @property
    def width(self) -> int:
        """Number of raster columns (``nx + 1``)."""
        return self.nx + 1

    @property
    def height(self) -> int:
        """Number of raster rows (``ny + 1``)."""
        return self.ny + 1

    @property
    def bounds(self) -> Bounds:
        return (
            self.origin_x,
            self.origin_y,
            self.origin_x + self.width * self.pixel_size,
            self.origin_y + self.height * self.pixel_size,
        )

    @property
    def transform(self) -> Affine:
        """North-up affine transform; rasterio anchors on the upper-left corner."""
        _, _, _, ymax = self.bounds
        return Affine.translation(self.origin_x, ymax) * Affine.scale(
            self.pixel_size, -self.pixel_size
        )

    def describe(self) -> str:
        w = self.width * self.pixel_size
        h = self.height * self.pixel_size
        return (
            f"{self.name}: {self.width}x{self.height} px @ {self.pixel_size:g} m "
            f"= {w / 1000:.3f} x {h / 1000:.3f} km, origin ({self.origin_x:.1f}, "
            f"{self.origin_y:.1f})"
        )


@dataclass(frozen=True)
class Tile:
    """A square work tile, keyed by its lower-left corner in whole kilometres."""

    east_km: int
    north_km: int

    @property
    def key(self) -> str:
        return f"{self.east_km}_{self.north_km}"

    @property
    def bounds(self) -> Bounds:
        x = self.east_km * TILE_SIZE
        y = self.north_km * TILE_SIZE
        return (x, y, x + TILE_SIZE, y + TILE_SIZE)

    def transform(self, resolution: float) -> Affine:
        xmin, _, _, ymax = self.bounds
        return Affine.translation(xmin, ymax) * Affine.scale(resolution, -resolution)

    def shape(self, resolution: float) -> Tuple[int, int]:
        n = int(round(TILE_SIZE / resolution))
        return n, n


@dataclass
class Tier:
    """A group of domains processed at a common work resolution.

    Rasterising a 48 km root domain at 1 m only to average it down to 32 m
    costs ~1000x more ortho data than it can possibly resolve. Domains are
    therefore binned by the resolution they actually need.
    """

    resolution: float
    domains: List[Domain]
    buffer: float

    @property
    def name(self) -> str:
        return f"r{self.resolution:g}m"

    @property
    def bounds(self) -> Bounds:
        """Union of member domain bounds, buffered and snapped outward to tiles."""
        xs0, ys0, xs1, ys1 = zip(*(d.bounds for d in self.domains))
        b = (
            min(xs0) - self.buffer,
            min(ys0) - self.buffer,
            max(xs1) + self.buffer,
            max(ys1) + self.buffer,
        )
        return snap_bounds(b, TILE_SIZE)

    @property
    def tiles(self) -> List[Tile]:
        """Every work tile touching a member domain -- not the bounding box.

        Domains sharing a resolution need not be adjacent: two 8 m nests 19 km
        apart span a bounding box two and a half times their combined area, and
        filling it costs that ratio in downloads and processing. Tiles are
        therefore taken per domain and unioned.
        """
        seen: Dict[Tuple[int, int], Tile] = {}
        for domain in self.domains:
            xmin, ymin, xmax, ymax = domain.bounds
            padded = snap_bounds(
                (xmin - self.buffer, ymin - self.buffer,
                 xmax + self.buffer, ymax + self.buffer),
                TILE_SIZE,
            )
            for tile in iter_tiles(padded):
                seen.setdefault((tile.east_km, tile.north_km), tile)
        return [seen[k] for k in sorted(seen)]

    def area_km2_covered(self) -> float:
        """Area actually processed, i.e. the tiles, not the bounding box."""
        return len(self.tiles) * (TILE_SIZE / 1000.0) ** 2

    def area_km2(self) -> float:
        xmin, ymin, xmax, ymax = self.bounds
        return (xmax - xmin) * (ymax - ymin) / 1e6


def snap_bounds(b: Bounds, step: float) -> Bounds:
    """Expand ``b`` outward to the next multiple of ``step``."""
    xmin, ymin, xmax, ymax = b
    return (
        math.floor(xmin / step) * step,
        math.floor(ymin / step) * step,
        math.ceil(xmax / step) * step,
        math.ceil(ymax / step) * step,
    )


def iter_tiles(b: Bounds) -> Iterator[Tile]:
    """Yield every work tile intersecting ``b`` (which must be tile-snapped)."""
    xmin, ymin, xmax, ymax = b
    e0 = int(math.floor(xmin / TILE_SIZE))
    e1 = int(math.ceil(xmax / TILE_SIZE))
    n0 = int(math.floor(ymin / TILE_SIZE))
    n1 = int(math.ceil(ymax / TILE_SIZE))
    for north in range(n0, n1):
        for east in range(e0, e1):
            yield Tile(east, north)


def intersects(a: Bounds, b: Bounds) -> bool:
    """True if the two bounding boxes overlap with positive area."""
    return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])


def choose_work_resolution(
    pixel_size: float,
    detail_factor: float,
    minimum: float,
    maximum: float,
    snap: Sequence[float],
) -> float:
    """Pick the work resolution for a domain of the given grid spacing.

    The target is ``pixel_size / detail_factor`` -- fine enough that the
    majority filter over each PALM cell sees ``detail_factor**2`` samples --
    clamped to the configured range and snapped down to an allowed value.
    """
    target = max(minimum, min(maximum, pixel_size / float(detail_factor)))
    allowed = sorted(float(s) for s in snap)
    candidates = [s for s in allowed if s <= target + 1e-9]
    return candidates[-1] if candidates else allowed[0]


def build_tiers(
    domains: Sequence[Domain],
    buffer: float,
    detail_factor: float,
    minimum: float,
    maximum: float,
    snap: Sequence[float],
    fixed: float | None = None,
) -> List[Tier]:
    """Bin ``domains`` into resolution tiers, finest first."""
    binned: Dict[float, List[Domain]] = {}
    for d in domains:
        res = (
            float(fixed)
            if fixed
            else choose_work_resolution(d.pixel_size, detail_factor, minimum, maximum, snap)
        )
        binned.setdefault(res, []).append(d)
    return [Tier(res, binned[res], buffer) for res in sorted(binned)]


def alignment_report(domain: Domain, work_resolution: float) -> Tuple[bool, str]:
    """Check whether a domain cell is an exact block of work pixels.

    Two conditions, and they are different questions:

    *ratio*  -- the domain's pixel size is a whole multiple of the work
    resolution, so each domain cell covers a whole number of work pixels.

    *phase*  -- the domain origin sits on the work lattice, which is anchored
    at the CRS origin (a work tile starts at a whole kilometre, see
    ``Tile.bounds``), so the block boundaries line up with pixel boundaries
    rather than cutting through them.

    Both must hold for aggregation (``mode`` / ``average``) to be a clean block
    statistic. When one does not, GDAL still produces a sensible-looking
    result -- it resamples -- but each domain cell is then a weighted blend of
    partially covered work pixels, and a majority vote over a blend is not a
    majority vote. This is reported rather than enforced because the result is
    usable and the fix (moving an origin) changes what is simulated, which is
    the user's decision and not palmgeo's.
    """
    ratio = domain.pixel_size / work_resolution
    ratio_ok = abs(ratio - round(ratio)) < 1e-6
    off_x = domain.origin_x / work_resolution
    off_y = domain.origin_y / work_resolution
    phase_ok = (
        abs(off_x - round(off_x)) < 1e-6 and abs(off_y - round(off_y)) < 1e-6
    )
    if ratio_ok and phase_ok:
        return True, f"{domain.name}: exact ({int(round(ratio))}x{int(round(ratio))} work px per cell)"
    reasons = []
    if not ratio_ok:
        reasons.append(
            f"pixel_size {domain.pixel_size:g} is not a multiple of work resolution "
            f"{work_resolution:g}"
        )
    if not phase_ok:
        reasons.append(
            f"origin ({domain.origin_x:g}, {domain.origin_y:g}) is not on the "
            f"{work_resolution:g} m work lattice"
        )
    return False, f"{domain.name}: approximate -- " + "; ".join(reasons)
