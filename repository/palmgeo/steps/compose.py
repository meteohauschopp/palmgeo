"""Fusion of land use and land cover into PALM surface types.

The classification rule per pixel, in the order it is applied:

1.  **Land use** from the cadastre gives the *context* (a residential block, a
    road, an arable field).
2.  **Land cover** from NDVI decides whether that pixel is actually vegetated,
    selecting the class's ``high_ndvi`` or ``low_ndvi`` outcome. This is what
    turns a single "Wohnbauflaeche" polygon into lawn, driveway and hedge at
    1 m resolution.
3.  **Height** from the nDSM promotes vegetated pixels above
    ``height_high_vegetation_lower_threshold`` to a tall vegetation type.
4.  **Phenology** from a second, leaf-off orthophoto epoch splits tall
    vegetation into evergreen / deciduous / mixed.
5.  **Hard overrides** are burnt in last: user polygons, then water, then
    buildings and bridges.

palm_csd applies its own precedence on top of this
(``create_driver.py::process_surface_types``: building > water > pavement >
vegetation), so overlaps are resolved consistently rather than being an error.

Provenance of the greenness criterion
-------------------------------------
``greenness_mask`` below is the colour criterion of the vegetation mask in
`klassifizierung_segmentierung.py` by **Joana-Pauline Rentz**, from her
Bachelor's thesis (full reference in :mod:`palmgeo.steps.trees`). Its point is that NDVI and hue fail in different places: NDVI is a
*ratio* against the near infrared, and in deep shadow both bands collapse
towards the sensor noise floor, so the ratio becomes meaningless while the
*hue* of the shaded leaf is still unmistakably green. The original used it to
mask vegetation before crown segmentation; here it is the second chance a pixel
gets after the NDVI threshold has rejected it.

The three thresholds are the original's, on the normalised bands
g_n = G/(R+G+B) and likewise r_n, b_n:

    g_n > 0.35            green must carry more than a third of the brightness
    g_n > r_n + 0.05      and beat red by a clear margin -- this is what
                          separates foliage from bare soil and roof tile
    g_n > b_n + 0.005     and merely beat blue, a much weaker condition,
                          because blue is low in almost everything outdoors

Greenness alone is not allowed to classify a pixel as vegetation anywhere in
this module: it is always combined with a relaxed NDVI floor or a height
criterion, so that an algal pond or a moss-covered roof cannot pass.
"""

from __future__ import annotations

import logging
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import ndimage

# rasterize() warns that the pre-allocated `out` array carries no geotransform.
# The transform is passed explicitly on every call, so the warning is noise.
try:  # pragma: no cover - depends on the rasterio version
    from rasterio.errors import NotGeoreferencedWarning

    import warnings as _warnings

    _warnings.filterwarnings("ignore", category=NotGeoreferencedWarning)
except ImportError:  # pragma: no cover
    pass

from palmgeo import palm_types as pt
from palmgeo.config import Config, Mapping, SurfaceRule, marker_up_to_date
from palmgeo.grid import Tier, Tile
from palmgeo.rasters import gdal_env, write_raster

logger = logging.getLogger("palmgeo.compose")

U8_NODATA = 255
I32_NODATA = -9999

#: Vector layers compose_tile reads. Their content is part of the stage
#: fingerprint, so re-running `vectors` or `buildings` recomposes the tiles.
COMPOSE_VECTORS = ("cadastre_landuse", "buildings", "bridges", "street_type", "street_crossings")

#: Products written per tile, with the resampling rule used later when they are
#: aggregated onto a PALM domain grid. `mode` for anything categorical,
#: `average` for continuous fields; nodata is honoured by both, so a domain
#: cell only averages the pixels that actually carry the quantity.
PRODUCTS: Dict[str, Dict[str, Any]] = {
    "zt":                {"dtype": "float32", "nodata": np.nan,     "resampling": "average"},
    "surface_group":     {"dtype": "uint8",   "nodata": U8_NODATA,  "resampling": "mode"},
    "vegetation_type":   {"dtype": "uint8",   "nodata": U8_NODATA,  "resampling": "mode"},
    "pavement_type":     {"dtype": "uint8",   "nodata": U8_NODATA,  "resampling": "mode"},
    "water_type":        {"dtype": "uint8",   "nodata": U8_NODATA,  "resampling": "mode"},
    "soil_type":         {"dtype": "uint8",   "nodata": U8_NODATA,  "resampling": "mode"},
    "vegetation_height": {"dtype": "float32", "nodata": np.nan,     "resampling": "average"},
    "buildings_2d":      {"dtype": "float32", "nodata": np.nan,     "resampling": "average"},
    "building_id":       {"dtype": "int32",   "nodata": I32_NODATA, "resampling": "mode"},
    "building_type":     {"dtype": "uint8",   "nodata": U8_NODATA,  "resampling": "mode"},
    "bridges_2d":        {"dtype": "float32", "nodata": np.nan,     "resampling": "average"},
    "bridges_id":        {"dtype": "int32",   "nodata": I32_NODATA, "resampling": "mode"},
    "street_type":       {"dtype": "uint8",   "nodata": U8_NODATA,  "resampling": "mode"},
    "street_crossings":  {"dtype": "uint8",   "nodata": U8_NODATA,  "resampling": "mode"},
    "surface_class":     {"dtype": "int16",   "nodata": 0,          "resampling": "mode"},
}


# ---------------------------------------------------------------------------
# rule table
# ---------------------------------------------------------------------------
@dataclass
class RuleTable:
    """Land-use rules flattened into arrays for vectorised lookup.

    Index 0 is always the fallback rule, which is what pixels outside any
    land-use polygon get.
    """

    keys: List[str]
    index_of: Dict[str, int]
    water_type: np.ndarray
    high_group: np.ndarray
    high_type: np.ndarray
    low_group: np.ndarray
    low_type: np.ndarray
    hv_group: np.ndarray
    hv_type: np.ndarray
    hv_fixed: np.ndarray
    #: Where NDVI is unavailable, fall back to the `high` outcome for classes
    #: that are vegetation either way (fields, forest) and to `low` otherwise.
    prefer_high_without_ndvi: np.ndarray

    @classmethod
    def build(cls, mapping: Mapping, observed_keys: Sequence[str]) -> "RuleTable":
        rules: List[Tuple[str, SurfaceRule]] = [("<fallback>", mapping.fallback)]
        index_of: Dict[str, int] = {}
        for key in sorted(set(observed_keys)):
            index_of[key] = len(rules)
            rules.append((key, mapping.lookup(key)))

        n = len(rules)
        water_type = np.zeros(n, dtype="uint8")
        high_group = np.zeros(n, dtype="uint8")
        high_type = np.zeros(n, dtype="uint8")
        low_group = np.zeros(n, dtype="uint8")
        low_type = np.zeros(n, dtype="uint8")
        hv_group = np.zeros(n, dtype="uint8")
        hv_type = np.zeros(n, dtype="uint8")
        hv_fixed = np.zeros(n, dtype=bool)
        prefer_high = np.zeros(n, dtype=bool)

        default_hv = mapping.default_high_vegetation
        for i, (_, rule) in enumerate(rules):
            if rule.water is not None:
                water_type[i] = rule.water[1]
                continue
            high_group[i], high_type[i] = rule.high_ndvi  # type: ignore[misc]
            low_group[i], low_type[i] = rule.low_ndvi     # type: ignore[misc]
            hv = rule.high_veg or default_hv
            hv_group[i], hv_type[i] = hv
            hv_fixed[i] = rule.high_veg_fixed
            prefer_high[i] = low_group[i] == pt.GROUP_VEGETATION

        return cls(
            keys=[k for k, _ in rules],
            index_of=index_of,
            water_type=water_type,
            high_group=high_group,
            high_type=high_type,
            low_group=low_group,
            low_type=low_type,
            hv_group=hv_group,
            hv_type=hv_type,
            hv_fixed=hv_fixed,
            prefer_high_without_ndvi=prefer_high,
        )


def landuse_key_field(src: Dict[str, Any]) -> str:
    """The cadastre land-use layer's key column, as the region file names it."""
    from palmgeo.steps.vectors import declared_field

    spec = (src.get("layers") or {}).get("landuse") or {}
    return str(declared_field(spec, "key_field", "cadastre layers.landuse"))


def collect_landuse_keys(cfg: Config) -> Tuple[List[str], List[str]]:
    """Distinct land-use keys present in the AOI, and those with no rule."""
    import geopandas as gpd

    from palmgeo.steps.vectors import resolve_field

    path = cfg.workdir / "vector" / "cadastre_landuse.gpkg"
    if not path.exists():
        return [], []
    _, src = cfg.source_for_role("cadastre")
    configured = landuse_key_field(src)
    gdf = gpd.read_file(path, ignore_geometry=True)
    field = resolve_field(gdf.columns, configured)
    if field is None:
        logger.error(
            "Land-use layer has no '%s' column (available: %s); every pixel would "
            "use the fallback rule.", configured, ", ".join(sorted(gdf.columns)),
        )
        return [], []
    keys = sorted({str(v).strip() for v in gdf[field].dropna()})
    unmapped = [k for k in keys if cfg.mapping.lookup(k) is cfg.mapping.fallback]
    return keys, unmapped


# ---------------------------------------------------------------------------
# per-tile composition
# ---------------------------------------------------------------------------
def _read(path: Path, shape: Tuple[int, int]) -> Optional[np.ndarray]:
    import rasterio

    if not path.exists():
        return None
    with rasterio.open(path) as src:
        data = src.read(1, masked=True)
        arr = data.filled(np.nan).astype("float32")
    if arr.shape != shape:
        return None
    return arr


def _read_bands(path: Path, shape: Tuple[int, int]) -> Optional[np.ndarray]:
    """Read the stored RGB+NIR cube of a tile, or ``None`` if it was not kept."""
    import rasterio

    if not path.exists():
        return None
    with rasterio.open(path) as src:
        cube = src.read().astype("float32")
    if cube.shape[1:] != shape or cube.shape[0] < 3:
        return None
    return cube


def greenness_mask(
    bands: np.ndarray,
    green_norm_min: float,
    green_over_red: float,
    green_over_blue: float,
) -> np.ndarray:
    """Hue-based vegetation indicator from the visible bands.

    Normalising by the total brightness makes the criterion largely independent
    of illumination, which is exactly where NDVI computed from uncalibrated
    8-bit digital numbers becomes unreliable: in shadow both the red and the
    near-infrared band collapse towards the noise floor and their ratio stops
    being meaningful, while the *hue* of a shaded leaf is still green.

    After Rentz; see "Provenance of the greenness criterion" above.
    """
    red, green, blue = bands[0], bands[1], bands[2]
    total = red + green + blue + 1e-6
    g_n, r_n, b_n = green / total, red / total, blue / total
    return (
        (g_n > green_norm_min)
        & (g_n > r_n + green_over_red)
        & (g_n > b_n + green_over_blue)
    )


def vegetation_mask(
    cfg: Config,
    ndvi: Optional[np.ndarray],
    ndsm: np.ndarray,
    bands: Optional[np.ndarray],
    height_threshold: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Decide per pixel whether it is vegetated, and where that is unknowable.

    A single NDVI threshold systematically loses two things, both verified over
    the Hannover test area: the shaded halves of tree crowns, and street trees
    standing over paved surfaces. Two further criteria recover them without
    inventing vegetation on open pavement -- of the pixels the combined rule
    adds, 79 % stand more than 2 m tall:

    ``ndvi >= ndvi_threshold``
        the primary criterion.
    ``greenness and ndvi >= ndvi_threshold_relaxed``
        shaded or dull foliage, where the hue still resolves what the ratio
        cannot.
    ``greenness and ndsm >= height_high_vegetation_lower_threshold``
        a tall, green, non-building object is a tree whatever its NDVI.

    Buildings and water are burnt in afterwards and win by precedence, so a
    greenness false positive on an algal pond or a mossy roof cannot survive.
    """
    shape = ndsm.shape
    if ndvi is None:
        return np.zeros(shape, dtype=bool), np.ones(shape, dtype=bool)

    section = cfg.get("classification.vegetation", {}) or {}
    threshold = float(
        section.get("ndvi_threshold", cfg.get("classification.ndvi_threshold", 0.12))
    )
    unknown = ~np.isfinite(ndvi)
    vegetated = np.where(unknown, False, ndvi >= threshold)

    green_cfg = section.get("greenness", {}) or {}
    if bands is not None and green_cfg.get("enabled", True):
        green = greenness_mask(
            bands,
            float(green_cfg.get("green_norm_min", 0.35)),
            float(green_cfg.get("green_over_red", 0.05)),
            float(green_cfg.get("green_over_blue", 0.005)),
        )
        relaxed = float(section.get("ndvi_threshold_relaxed", 0.05))
        vegetated |= green & np.where(unknown, False, ndvi >= relaxed)
        if section.get("tall_green_is_vegetation", True):
            vegetated |= green & np.isfinite(ndsm) & (ndsm >= height_threshold)
    return vegetated, unknown


def _rasterize(
    path: Path,
    tile: Tile,
    transform,
    shape: Tuple[int, int],
    value_field: str,
    dtype: str,
    fill: float,
    extra_fields: Sequence[str] = (),
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Burn a vector layer onto the tile grid, reading only the tile's bbox."""
    import geopandas as gpd
    from rasterio.features import rasterize

    out = np.full(shape, fill, dtype=dtype)
    extras = {f: np.full(shape, fill, dtype=dtype) for f in extra_fields}
    if not path.exists():
        return out, extras
    # small halo so polygons only touching the edge still paint their sliver
    b = tile.bounds
    halo = 2.0
    gdf = gpd.read_file(path, bbox=(b[0] - halo, b[1] - halo, b[2] + halo, b[3] + halo))
    if gdf.empty:
        return out, extras
    for field, target in [(value_field, out), *((f, extras[f]) for f in extra_fields)]:
        if field not in gdf.columns:
            continue
        pairs = [
            (geom, value)
            for geom, value in zip(gdf.geometry, gdf[field])
            if geom is not None and not geom.is_empty and value is not None and np.isfinite(float(value))
        ]
        if not pairs:
            continue
        rasterize(
            pairs, out_shape=shape, transform=transform, fill=fill,
            dtype=dtype, out=target, all_touched=False,
        )
    return out, extras


def _fill_from_neighbours(
    values: np.ndarray, ids: np.ndarray, inside: np.ndarray
) -> np.ndarray:
    """Fill gaps in ``values`` from the nearest measured pixel of the same building.

    Used where vegetation hides the roof: the structure below a crown is assumed
    to continue from the nearest structure that *was* measured on the same
    building. On a roof the nearest donor is the roof itself, in a courtyard it
    is the courtyard floor -- which is the whole point, because the two need
    opposite answers and a single per-building height can only serve one.

    Gaps whose nearest donor belongs to a different building, and buildings with
    no measured pixel at all, are left untouched; the caller falls back to the
    per-building height there.
    """
    gaps = inside & ~np.isfinite(values)
    donors = inside & np.isfinite(values)
    if not gaps.any() or not donors.any():
        return values

    # EDT of the donor complement yields, per pixel, the index of the nearest
    # donor -- one O(n) pass for the whole tile, no per-building loop.
    _, (yy, xx) = ndimage.distance_transform_edt(
        ~donors, return_indices=True
    )
    src_y, src_x = yy[gaps], xx[gaps]
    same = ids[src_y, src_x] == ids[gaps]

    out = values.copy()
    target = np.flatnonzero(gaps.ravel())[same]
    out.ravel()[target] = values[src_y[same], src_x[same]]
    return out


def _deck_heights(
    ndsm: np.ndarray,
    ids: np.ndarray,
    deck: np.ndarray,
    fallback: np.ndarray,
    min_height: float,
    radius: int,
    hidden: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Per-pixel height of a bridge deck above the terrain model.

    PALM stacks a bridge on the *local* terrain of each grid cell
    (topography_mod.f90: building type 7 is not lifted to the highest terrain
    under it, unlike other buildings), so the deck top it builds is
    ``zt + bridges_2d`` cell by cell. A single height per deck gets that wrong
    twice. Arched decks and their ramps vary by metres along the span, and a
    median places the ends too high and the middle too low -- 19-31 % of deck
    pixels by more than 1 m over five test projects. And where the terrain
    model was not cut below the deck, ``zt`` already *is* the deck, so adding
    the deck's height on top builds it twice as high.

    The nDSM per pixel is exactly ``deck top - zt``, which answers both:
    pixels where the terrain model carries the deck measure close to zero and
    drop below ``min_height``; they stay terrain instead of becoming a second
    deck. Vegetation (``hidden``) overhanging a deck is refilled from the
    nearest measured pixel of the same bridge, and a median over a
    ``(2 * radius + 1)``-pixel window within the same bridge removes vehicles,
    railings and lamp posts. Pixels without any measurement keep the per-deck
    ``fallback``.
    """
    measured = np.where(deck & np.isfinite(ndsm), ndsm, np.nan)
    if hidden is not None:
        measured = np.where(deck & hidden, np.nan, measured)
        measured = _fill_from_neighbours(measured, ids, deck)

    if radius > 0 and deck.any():
        size = 2 * radius + 1
        win_v = np.lib.stride_tricks.sliding_window_view(
            np.pad(measured, radius, constant_values=np.nan), (size, size))
        win_i = np.lib.stride_tricks.sliding_window_view(
            np.pad(ids, radius, constant_values=I32_NODATA), (size, size))
        yy, xx = np.nonzero(deck)
        values = win_v[yy, xx].reshape(len(yy), -1)
        same = win_i[yy, xx].reshape(len(yy), -1) == ids[yy, xx][:, None]
        with warnings.catch_warnings():
            # all-NaN windows are expected (no DSM there); they fall back below
            warnings.simplefilter("ignore", RuntimeWarning)
            smoothed = np.nanmedian(np.where(same, values, np.nan), axis=1)
        measured = measured.copy()
        measured[yy, xx] = smoothed

    height = np.where(np.isfinite(measured), measured, fallback)
    keep = deck & np.isfinite(height) & (height >= min_height)
    return np.where(keep, height, np.nan).astype("float32")


def compose_tile(
    cfg: Config, tier: Tier, tile: Tile, rules: RuleTable, vectors: Dict[str, Path],
    fingerprint: str,
) -> Optional[Path]:
    """Turn one tile's derived rasters plus the vectors into PALM surface types."""
    tiledir = cfg.tier_dir(tier) / "tiles" / tile.key
    if not (tiledir / ".done").exists():
        return None
    outdir = cfg.tier_dir(tier) / "surface" / tile.key
    if marker_up_to_date(outdir / ".done", fingerprint):
        return outdir

    res = tier.resolution
    transform = tile.transform(res)
    shape = tile.shape(res)
    crs = cfg.crs
    opts = cfg.get("output.cog", {}) or {}

    zt = _read(tiledir / "zt.tif", shape)
    if zt is None:
        return None
    ndsm = _read(tiledir / "ndsm.tif", shape)
    if ndsm is None:
        ndsm = np.full(shape, np.nan, dtype="float32")
    ndvi_on = _read(tiledir / "ndvi_leaf_on.tif", shape)
    ndvi_off = _read(tiledir / "ndvi_leaf_off.tif", shape)

    # ---- 1. land use ------------------------------------------------------
    _, src = cfg.source_for_role("cadastre")
    configured_key = landuse_key_field(src)
    rule_id = np.zeros(shape, dtype="int32")
    lu_path = vectors.get("cadastre_landuse")
    if lu_path is not None and lu_path.exists():
        import geopandas as gpd
        from rasterio.features import rasterize

        from palmgeo.steps.vectors import resolve_field

        b = tile.bounds
        gdf = gpd.read_file(lu_path, bbox=(b[0] - 2, b[1] - 2, b[2] + 2, b[3] + 2))
        key_field = resolve_field(gdf.columns, configured_key) if not gdf.empty else None
        if not gdf.empty and key_field:
            pairs = [
                (geom, rules.index_of.get(str(value).strip(), 0))
                for geom, value in zip(gdf.geometry, gdf[key_field])
                if geom is not None and not geom.is_empty
            ]
            pairs = [(g, v) for g, v in pairs if v > 0]
            if pairs:
                rasterize(
                    pairs, out_shape=shape, transform=transform, fill=0,
                    dtype="int32", out=rule_id, all_touched=False,
                )

    # ---- 2./3. land cover and height --------------------------------------
    h_thr = float(cfg.get("csd.settings.height_high_vegetation_lower_threshold", 2.0))
    bands = _read_bands(tiledir / "ortho_leaf_on.tif", shape)
    vegetated, unknown = vegetation_mask(cfg, ndvi_on, ndsm, bands, h_thr)
    vegetated = np.where(unknown, rules.prefer_high_without_ndvi[rule_id], vegetated)

    group = np.where(vegetated, rules.high_group[rule_id], rules.low_group[rule_id])
    stype = np.where(vegetated, rules.high_type[rule_id], rules.low_type[rule_id])

    tall = vegetated & np.isfinite(ndsm) & (ndsm >= h_thr) & (group == pt.GROUP_VEGETATION)
    group = np.where(tall, rules.hv_group[rule_id], group)
    stype = np.where(tall, rules.hv_type[rule_id], stype)

    # ---- 4. phenology ------------------------------------------------------
    if (
        cfg.get("classification.phenology.enabled", False)
        and ndvi_off is not None
        and np.isfinite(ndvi_off).any()
    ):
        ph = cfg.mapping.phenology
        evergreen = ph.get("evergreen")
        deciduous = ph.get("deciduous")
        mixed = ph.get("mixed")
        ever_min = float(cfg.get("classification.phenology.evergreen_leaf_off_min", 0.30))
        dec_max = float(cfg.get("classification.phenology.deciduous_leaf_off_max", 0.15))
        if evergreen and deciduous and mixed:
            pheno = np.where(
                ndvi_off >= ever_min, evergreen[1],
                np.where(ndvi_off < dec_max, deciduous[1], mixed[1]),
            ).astype("uint8")
            apply = tall & ~rules.hv_fixed[rule_id] & np.isfinite(ndvi_off)
            stype = np.where(apply, pheno, stype)

    # ---- water -------------------------------------------------------------
    is_water = rules.water_type[rule_id] > 0
    group = np.where(is_water, pt.GROUP_WATER, group).astype("uint8")
    stype = np.where(is_water, rules.water_type[rule_id], stype).astype("uint8")

    # ---- 5a. user-supplied surface polygons --------------------------------
    for extra in cfg.get("extra_inputs.surfaces", []) or []:
        path = cfg.root / str(extra["path"])
        if not path.exists():
            logger.warning("extra_inputs.surfaces: %s does not exist, skipping", path)
            continue
        g_code, t_code = pt.parse_type(str(extra["type"]))
        import geopandas as gpd
        from rasterio.features import rasterize as _rio_rasterize

        b = tile.bounds
        gdf = gpd.read_file(path, layer=extra.get("layer"), bbox=(b[0] - 2, b[1] - 2, b[2] + 2, b[3] + 2))
        if gdf.empty:
            continue
        shapes = [(geom, 1) for geom in gdf.geometry if geom is not None and not geom.is_empty]
        if not shapes:
            continue
        burn = np.zeros(shape, dtype="uint8")
        _rio_rasterize(
            shapes, out_shape=shape, transform=transform, fill=0, dtype="uint8", out=burn,
        )
        group = np.where(burn == 1, g_code, group).astype("uint8")
        stype = np.where(burn == 1, t_code, stype).astype("uint8")

    # ---- 5b. buildings and bridges ----------------------------------------
    buildings_2d = np.full(shape, np.nan, dtype="float32")
    building_id = np.full(shape, I32_NODATA, dtype="int32")
    building_type = np.full(shape, U8_NODATA, dtype="uint8")
    bpath = vectors.get("buildings")
    if bpath is not None and bpath.exists():
        height, _ = _rasterize(
            bpath, tile, transform, shape, value_field="height",
            dtype="float32", fill=np.nan,
        )
        cap, _ = _rasterize(
            bpath, tile, transform, shape, value_field="height_cap",
            dtype="float32", fill=np.nan,
        )
        ids, _ = _rasterize(
            bpath, tile, transform, shape, value_field="building_id",
            dtype="int32", fill=I32_NODATA,
        )
        types, _ = _rasterize(
            bpath, tile, transform, shape, value_field="building_type",
            dtype="uint8", fill=U8_NODATA,
        )
        inside = ids != I32_NODATA
        min_height = float(cfg.get("buildings.min_height", 2.5))
        open_area = np.zeros(shape, dtype=bool)

        if cfg.get("buildings.per_pixel_height", True) and np.isfinite(ndsm).any():
            # A single height per footprint throws away a lot: over the Hannover
            # test area half the buildings span more than 5 m between the 5th
            # and 95th percentile of their own nDSM (median 5.0 m, p90 20.9 m).
            # Taking the measured surface per pixel keeps roof ridges, setbacks
            # and towers.
            measured = np.where(inside & np.isfinite(ndsm), ndsm, np.nan)

            # Vegetation inside a footprint is a tree crown or a green roof, not
            # structure, so the nDSM there measures leaves rather than the
            # building. Those pixels are dropped and refilled from the nearest
            # measured pixel of the *same* building (see _fill_from_neighbours).
            #
            # Substituting the building's own height instead -- which is what
            # this did before -- is only right when the vegetation sits on the
            # roof. In a courtyard, a light well or a setback the vegetation
            # stands well below the roof, and the substitution lifts the pixel
            # to the full building height: over one 1 m test tier it
            # inflated 0.4 % of all footprint area by more than 5 m, and turned
            # a 7 m courtyard floor inside a 15 m building into 23 m of solid
            # masonry. The nearest structural neighbour is right in both cases.
            #
            # The test uses NDVI alone rather than the combined vegetation mask.
            # The hue criterion cannot tell chlorophyll from green paint or
            # verdigris: 8 % of the Neues Rathaus dome passes it, because
            # patinated copper is green in the visible -- but its NDVI is -0.39,
            # because copper has no near-infrared signature. Using the combined
            # mask here would flatten the dome from 97.6 m to 92.9 m.
            if cfg.get("buildings.ignore_vegetation_on_roofs", True) and ndvi_on is not None:
                threshold = float(
                    (cfg.get("classification.vegetation", {}) or {}).get(
                        "ndvi_threshold", cfg.get("classification.ndvi_threshold", 0.12)
                    )
                )
                chlorophyll = np.isfinite(ndvi_on) & (ndvi_on >= threshold)
                measured = np.where(inside & chlorophyll, np.nan, measured)
                measured = _fill_from_neighbours(measured, ids, inside)

            # Optional second net against vegetation the mask missed: nothing
            # may exceed the building's own 99th percentile by more than the
            # tolerance. Off by default, because a genuine tower is exactly the
            # kind of feature that lives above its building's p99 -- with 1.15
            # the Neues Rathaus loses 12 m off its dome (97.6 -> 85.0 m).
            tolerance = cfg.get("buildings.per_pixel_tolerance")
            if tolerance:
                limit = np.where(np.isfinite(cap), cap * float(tolerance), np.inf)
                measured = np.minimum(measured, limit)

            if cfg.get("buildings.carve_open_areas", True):
                # Courtyards, light wells and passages that the cadastre still
                # maps as building but the surface model measures as open.
                # 2.5 % of footprint area over the test domain, and 5 of 82
                # buildings are more than a quarter open.
                open_area = inside & np.isfinite(ndsm) & (ndsm < min_height)

            # where the surface model has no value, keep the per-building height
            height = np.where(np.isfinite(measured), measured, height)

        valid = inside & ~open_area & np.isfinite(height) & (height >= min_height)
        buildings_2d = np.where(valid, height, np.nan).astype("float32")
        building_id = np.where(valid, ids, I32_NODATA).astype("int32")
        building_type = np.where(valid, types, U8_NODATA).astype("uint8")
        group = np.where(valid, pt.GROUP_BUILDING, group).astype("uint8")

    bridges_2d = np.full(shape, np.nan, dtype="float32")
    bridges_id = np.full(shape, I32_NODATA, dtype="int32")
    brpath = vectors.get("bridges")
    if brpath is not None and brpath.exists():
        bh, _ = _rasterize(brpath, tile, transform, shape, "height", "float32", np.nan)
        bi, _ = _rasterize(brpath, tile, transform, shape, "bridge_id", "int32", I32_NODATA)
        # a bridge deck is only a bridge where no building already stands
        free = (bi != I32_NODATA) & np.isfinite(bh) & (building_id == I32_NODATA)
        if cfg.get("buildings.bridges.per_pixel_height", True) and np.isfinite(ndsm).any():
            hidden = None
            if ndvi_on is not None:
                threshold = float(
                    (cfg.get("classification.vegetation", {}) or {}).get(
                        "ndvi_threshold", cfg.get("classification.ndvi_threshold", 0.12)
                    )
                )
                hidden = np.isfinite(ndvi_on) & (ndvi_on >= threshold)
            bridges_2d = _deck_heights(
                ndsm, bi, free, bh,
                min_height=float(cfg.get("buildings.bridges.min_height", 0.5)),
                radius=int(cfg.get("buildings.bridges.smooth_radius", 1)),
                hidden=hidden,
            )
        else:
            bridges_2d = np.where(free, bh, np.nan).astype("float32")
        bridges_id = np.where(np.isfinite(bridges_2d), bi, I32_NODATA).astype("int32")

    # ---- derived fields ----------------------------------------------------
    is_veg = group == pt.GROUP_VEGETATION
    is_pav = group == pt.GROUP_PAVEMENT
    is_wat = group == pt.GROUP_WATER

    vegetation_type = np.where(is_veg, stype, U8_NODATA).astype("uint8")
    pavement_type = np.where(is_pav, stype, U8_NODATA).astype("uint8")
    water_type = np.where(is_wat, stype, U8_NODATA).astype("uint8")

    vegetation_height = np.where(
        is_veg & np.isfinite(ndsm), np.maximum(ndsm, 0.0), np.nan
    ).astype("float32")

    soil = pt.SOIL_TYPE[str(cfg.mapping.soil.get("constant", "medium_fine"))]
    soil_type = np.where(is_veg | is_pav, soil, U8_NODATA).astype("uint8")
    for prefix, name in (cfg.mapping.soil.get("by_landuse_prefix") or {}).items():
        # `rule_ids`, not `ids`: `ids` is the rasterised building_id array a
        # hundred lines above, and reusing the name here reads like the two
        # are related. They are not.
        rule_ids = [i for k, i in rules.index_of.items() if k.startswith(str(prefix))]
        if rule_ids:
            hit = np.isin(rule_id, rule_ids) & (soil_type != U8_NODATA)
            soil_type = np.where(hit, pt.SOIL_TYPE[str(name)], soil_type).astype("uint8")

    street_type = np.full(shape, U8_NODATA, dtype="uint8")
    spath = vectors.get("street_type")
    if spath is not None and spath.exists():
        st, _ = _rasterize(spath, tile, transform, shape, "street_type", "uint8", U8_NODATA)
        street_type = st
    street_crossings = np.full(shape, U8_NODATA, dtype="uint8")
    cpath = vectors.get("street_crossings")
    if cpath is not None and cpath.exists():
        cr, _ = _rasterize(cpath, tile, transform, shape, "crossing", "uint8", U8_NODATA)
        street_crossings = cr

    # single-band, human-readable classification for QGIS
    surface_class = np.zeros(shape, dtype="int16")
    for g in (pt.GROUP_VEGETATION, pt.GROUP_PAVEMENT, pt.GROUP_WATER, pt.GROUP_BUILDING):
        sel = group == g
        if g == pt.GROUP_BUILDING:
            surface_class = np.where(sel, pt.INSPECTION_OFFSET[g], surface_class).astype("int16")
        else:
            surface_class = np.where(
                sel, pt.INSPECTION_OFFSET[g] + stype.astype("int16"), surface_class
            ).astype("int16")

    arrays = {
        "zt": zt,
        "surface_group": np.where(group == 0, U8_NODATA, group).astype("uint8"),
        "vegetation_type": vegetation_type,
        "pavement_type": pavement_type,
        "water_type": water_type,
        "soil_type": soil_type,
        "vegetation_height": vegetation_height,
        "buildings_2d": buildings_2d,
        "building_id": building_id,
        "building_type": building_type,
        "bridges_2d": bridges_2d,
        "bridges_id": bridges_id,
        "street_type": street_type,
        "street_crossings": street_crossings,
        "surface_class": surface_class,
    }
    try:
        for name, arr in arrays.items():
            meta = PRODUCTS[name]
            write_raster(
                outdir / f"{name}.tif", arr.astype(meta["dtype"], copy=False),
                transform, crs, meta["nodata"], opts, overviews=False,
            )
    except Exception:
        for f in outdir.glob("*"):
            f.unlink(missing_ok=True)
        raise
    (outdir / ".done").write_text(fingerprint)
    return outdir


def run(cfg: Config, tier: Tier, vectors: Dict[str, Path]) -> Dict[str, Any]:
    """Compose every tile of one tier."""
    keys, unmapped = collect_landuse_keys(cfg)
    if unmapped:
        logger.warning(
            "%d land-use key(s) in the AOI have no explicit rule and use the "
            "fallback: %s", len(unmapped), ", ".join(unmapped[:15]),
        )
    rules = RuleTable.build(cfg.mapping, keys)

    # Hashed once here rather than per tile: the layers can be large.
    fingerprint = cfg.stage_fingerprint(
        "compose", inputs={name: vectors.get(name) for name in COMPOSE_VECTORS}
    )

    tiles = tier.tiles
    workers = int(cfg.runtime.get("workers", 4))
    done = 0
    failures: List[Tuple[str, str]] = []
    with gdal_env(int(cfg.runtime.get("gdal_cachemax_mb", 512))):
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(compose_tile, cfg, tier, t, rules, vectors, fingerprint): t
                for t in tiles
            }
            for future in as_completed(futures):
                tile = futures[future]
                try:
                    future.result()
                except Exception as exc:  # noqa: BLE001
                    failures.append((tile.key, str(exc)))
                    logger.error("Compose failed for tile %s: %s", tile.key, exc)
                done += 1
                if done % 25 == 0 or done == len(tiles):
                    logger.info("  %d/%d tiles composed", done, len(tiles))
    return {
        "tier": tier.name,
        "tiles": len(tiles),
        "failed": len(failures),
        "unmapped_landuse_keys": unmapped,
    }
