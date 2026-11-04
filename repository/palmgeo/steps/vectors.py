"""Cadastre, 3D building models and topographic vectors.

Vector sources are handled differently from rasters for one reason: a building
footprint, a land-use polygon or a road centreline can straddle work-tile
boundaries. Anything computed per feature (a building's height, for instance)
must therefore be computed once over the whole area of interest, not per tile.

The stage output is a small set of local GeoPackages, each carrying an R-tree
spatial index, from which the tile workers later read only their own bounding
box. The multi-hundred-megabyte source archives are deleted as soon as their
content has been clipped.
"""

from __future__ import annotations

import logging
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from palmgeo import palm_types as pt
from palmgeo.cache import Cache
from palmgeo.config import Config, ConfigError
from palmgeo.grid import Bounds
from palmgeo.sources import Index, Item, load_index, newest_per_tile, regions_for_bounds

logger = logging.getLogger("palmgeo.vectors")

#: Read the source in stripes of this width so peak memory stays bounded even
#: for a 48 km area of interest against a whole-Landkreis GeoPackage.
CLIP_STRIPE_M = 4000.0


def resolve_field(columns: Iterable[str], name: str) -> Optional[str]:
    """Match a configured field name against actual columns, ignoring case.

    Attribute spellings differ between deliveries of the same product -- the
    Basis-DLM ships ``WDM`` where the specification writes ``wdm`` -- and a
    silent miss degrades into "everything is the default class" rather than an
    error. Resolving case-insensitively keeps the mapping files portable.
    """
    columns = list(columns)
    if name in columns:
        return name
    lowered = {c.lower(): c for c in columns}
    return lowered.get(name.lower())


def _aoi(cfg: Config) -> Bounds:
    """Union of all tier bounds -- everything the run can possibly need."""
    xs0, ys0, xs1, ys1 = zip(*(t.bounds for t in cfg.tiers))
    return (min(xs0), min(ys0), max(xs1), max(ys1))


def _write(gdf: gpd.GeoDataFrame, path: Path, append: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(path, driver="GPKG", mode="a" if append and path.exists() else "w")


def clip_layer(
    source: Path,
    layer: str,
    aoi: Bounds,
    out_path: Path,
    columns: Optional[Sequence[str]] = None,
    crs: Optional[str] = None,
) -> int:
    """Clip one layer to the area of interest, striped to bound memory."""
    if out_path.exists():
        out_path.unlink()
    total = 0
    first = True
    x = aoi[0]
    while x < aoi[2]:
        stripe = (x, aoi[1], min(x + CLIP_STRIPE_M, aoi[2]), aoi[3])
        try:
            gdf = gpd.read_file(source, layer=layer, bbox=stripe)
        except Exception as exc:  # noqa: BLE001 - surface as a clear message
            raise RuntimeError(f"Cannot read layer '{layer}' from {source}: {exc}") from exc
        x += CLIP_STRIPE_M
        if gdf.empty:
            continue
        if columns:
            keep = [c for c in (resolve_field(gdf.columns, name) for name in columns) if c]
            gdf = gdf[[*keep, gdf.geometry.name]]
        if crs and gdf.crs is None:
            gdf = gdf.set_crs(crs)
        # Drop geometries that are broken beyond repair. `buffer(0)` is the
        # standard trick for closing a self-intersecting ring, and asking
        # whether the repaired version is valid is how we tell "fixable" from
        # "not a polygon at all".
        #
        # Note what does NOT happen here: features are not cut to the stripe or
        # to the AOI. A bbox filter returns whole features, and cutting them
        # would split every polygon that crosses a stripe edge into pieces --
        # turning one building into two. The consequence is that a feature
        # straddling a stripe edge is returned by both stripes and written
        # twice; that is dealt with once, on merge, by dropping exact
        # duplicates rather than by mutilating geometry here.
        gdf = gdf[gdf.geometry.is_valid | gdf.geometry.buffer(0).is_valid]
        _write(gdf, out_path, append=not first)
        first = False
        total += len(gdf)
    logger.info("  %-12s -> %s (%d features)", layer, out_path.name, total)
    return total


def declared_field(spec: Dict[str, Any], key: str, where: str, required: bool = True) -> Optional[str]:
    """A column name the region file declares -- never a built-in default.

    The defaults that used to stand here were Lower Saxony's column names
    (ALKIS ``objid``/``gfk``/``tn``, LoD2 ``measHeight``, Basis-DLM ``WDM``).
    In any other region they match nothing, and every consequence of that is
    silent: no heights, every building the default type, every road
    unclassified. So a missing key is an error that names the key. A source
    that genuinely lacks the column says so with an explicit ``null``, which
    is accepted unless ``required``.
    """
    if key not in spec:
        raise ConfigError(
            f"{where}: `{key}` is not set. Name the column in the region file"
            + ("." if required else ", or set it to null if the source has none.")
        )
    value = spec[key]
    if value is None or str(value) == "":
        if required:
            raise ConfigError(f"{where}: `{key}` is required and cannot be null")
        return None
    return str(value)


# ---------------------------------------------------------------------------
# cadastre (e.g. ALKIS in Lower Saxony)
# ---------------------------------------------------------------------------
def fetch_cadastre(cfg: Config, cache: Cache) -> Dict[str, Path]:
    """Download the administrative units covering the AOI and clip their layers."""
    name, src = cfg.source_for_role("cadastre")
    aoi = _aoi(cfg)
    outdir = cfg.workdir / "vector"
    outdir.mkdir(parents=True, exist_ok=True)

    index = load_index(name, src, cache, transform_to_epsg=cfg.epsg)
    units = regions_for_bounds(index, aoi)
    if not units:
        raise RuntimeError(
            f"No {name} administrative unit intersects the area of interest {aoi}. "
            f"Check grid.epsg and the domain origins."
        )
    logger.info("Cadastre: %d administrative unit(s) intersect the AOI", len(units))

    layers = src.get("layers") or {}
    parts: Dict[str, List[Path]] = {role: [] for role in layers}
    for unit in units:
        archive = cache.fetch(unit.url, subdir=name)
        with cache.hold(archive):
            folder = cache.unzip(archive)
            gpkgs = sorted(folder.rglob("*.gpkg"))
            if not gpkgs:
                raise RuntimeError(f"No GeoPackage inside {archive}")
            for role, spec in layers.items():
                out = outdir / f"_{role}_{unit.key}_{archive.stem}.gpkg"
                n = clip_layer(
                    gpkgs[0],
                    str(spec["name"]),
                    aoi,
                    out,
                    columns=[
                        str(v)
                        for k, v in spec.items()
                        if k.endswith("_field") and isinstance(v, str)
                    ],
                    crs=cfg.crs,
                )
                if n:
                    parts[role].append(out)
            cache.discard(folder)
        cache.release(archive)

    merged: Dict[str, Path] = {}
    for role, files in parts.items():
        if not files:
            logger.warning("Cadastre layer '%s' is empty over the AOI", role)
            continue
        # prefixed so the raw cadastre extract is never overwritten by the
        # assembled `buildings.gpkg` that the buildings stage produces
        target = outdir / f"cadastre_{role}.gpkg"
        if target.exists():
            target.unlink()
        frames = []
        for f in files:
            frames.append(gpd.read_file(f))
            f.unlink()
        gdf = _dedupe(pd.concat(frames, ignore_index=True), role)
        _write(gpd.GeoDataFrame(gdf, geometry=gdf.geometry.name), target, append=False)
        merged[role] = target
    return merged


def _dedupe(gdf: "gpd.GeoDataFrame", label: str) -> "gpd.GeoDataFrame":
    """Drop rows that repeat a feature already present, geometry and all.

    The same feature reaches this point twice for two independent reasons, and
    both are properties of the data rather than mistakes:

    *Administrative extracts overlap.* A cadastre is published per district,
    and a building on a district boundary is delivered by both districts. A
    small area of interest can intersect dozens of them -- 45 for a 2 km
    domain over Hannover -- and 8 % of its footprints arrive twice.

    *Read stripes overlap.* ``clip_layer`` reads in 4 km stripes to bound
    memory, and a bbox filter returns whole features, so anything crossing a
    stripe edge comes back from both stripes.

    Measured across the projects that exist, this was 0.4 % to 3.5 % of every
    cadastre layer. The duplicates were exact -- identical geometry AND
    identical attributes in every group checked -- so the driver was correct;
    what was wrong was every count reported to the user, and roughly 2.5 % of
    the height sampling was done twice for nothing.

    Only exact duplicates are dropped. Two rows with the same geometry but
    different attributes are a genuine disagreement between two deliveries of
    the same feature, and choosing a winner is a data decision palmgeo has no
    basis for making. Those are kept and reported, so that somebody can look.
    """
    if gdf.empty:
        return gdf
    geom_name = gdf.geometry.name
    attrs = [c for c in gdf.columns if c != geom_name]
    key = gdf.geometry.to_wkb()

    exact = gdf.assign(_key=key).duplicated(subset=["_key", *attrs], keep="first")
    conflicting = 0
    repeated = key.duplicated(keep=False)
    if repeated.any():
        grouped = gdf.assign(_key=key)[repeated].groupby("_key")[attrs]
        conflicting = int(sum(len(g.drop_duplicates()) > 1 for _, g in grouped))

    if exact.any():
        logger.info(
            "  %-12s -> dropped %d exact duplicate(s) of %d feature(s)",
            label, int(exact.sum()), len(gdf),
        )
    if conflicting:
        logger.warning(
            "%s: %d geometry/geometries appear more than once with DIFFERENT "
            "attributes. These are kept as they are -- two deliveries of the "
            "same feature disagree, and palmgeo has no basis for picking a "
            "winner. Check the source extracts if the classification looks odd "
            "in those places.", label, conflicting,
        )
    return gdf[~exact.to_numpy()]


# ---------------------------------------------------------------------------
# 3D building models (role `buildings3d`; CityGML LoD2 in Lower Saxony)
# ---------------------------------------------------------------------------
def fetch_building_heights(cfg: Config, cache: Cache) -> Optional[pd.DataFrame]:
    """Collect per-building heights from the 3D building model.

    LoD2 geometry is a set of individual roof and wall faces (``MULTIPOLYGON Z``),
    so several rows share one ``gml_id`` and the geometry is unusable as a
    footprint. Only the attributes are harvested here; they are joined onto the
    clean 2D cadastre footprints in :func:`build_buildings`.
    """
    # Ask first whether the 3D model is wanted at all -- this needs no source
    # lookup, and a region that has no buildings3d source (the US-NY port, for
    # instance, where the DSM covers the whole AOI and NYC Building Footprints
    # carry their own height_roof) must not fail merely by resolving the role.
    if "lod2" not in {
        str(cfg.get("buildings.height_source", "ndsm")),
        str(cfg.get("buildings.height_fallback", "none")),
    }:
        return None
    try:
        name, src = cfg.source_for_role("buildings3d")
    except (KeyError, ConfigError):
        logger.warning(
            "No buildings3d source is defined for this region -- building "
            "heights come from the surface model alone"
        )
        return None

    aoi = _aoi(cfg)
    index = load_index(name, src, cache, transform_to_epsg=cfg.epsg)
    tiles = {"aoi": aoi}
    items = [it for it in index.search(aoi)]
    if not items:
        logger.warning("No %s tiles cover the AOI -- falling back to nDSM heights", name)
        return None

    # one item per tile key, most recent wins
    by_key: Dict[str, Item] = {}
    for it in items:
        prev = by_key.get(it.key)
        if prev is None or (it.date or prev.date) and (it.date and prev.date and it.date > prev.date):
            by_key[it.key] = it
    logger.info("Building models: %d %s tile(s)", len(by_key), name)

    f = src.get("fields") or {}
    where = f"sources.{name}.fields"
    key_field = declared_field(f, "building_key", where)
    height_field = declared_field(f, "height", where)
    # The reference from a model building to its cadastre footprint, and how
    # to cut the cadastre id out of it; and the building function, with the
    # separator in front of the code. All optional: a model without them
    # simply joins by nothing and carries no function.
    ext_field = f.get("extern_ref") or None
    ext_split = f.get("extern_ref_split") or None
    area_field = f.get("roof_area") or None
    func_field = f.get("function") or None
    func_split = f.get("function_split") or None

    frames: List[pd.DataFrame] = []
    for item in by_key.values():
        archive = cache.fetch(item.url, subdir=name)
        with cache.hold(archive):
            folder = cache.unzip(
                archive, member_suffixes=(".shp", ".dbf", ".shx", ".prj", ".cpg")
            )
            shapes = sorted(folder.rglob("*.shp"))
            if not shapes:
                logger.warning("No shapefile inside %s", archive.name)
            for shp in shapes:
                # ignore_geometry is essential here, not just an optimisation:
                # LoD2 ships PolyhedralSurface/TIN (WKB type 16), which shapely
                # cannot parse at all.
                df = gpd.read_file(shp, ignore_geometry=True)
                wanted = (key_field, ext_field, height_field, area_field, func_field)
                cols = [c for c in wanted if c and c in df.columns]
                if height_field not in cols:
                    continue
                frames.append(df[cols].copy())
            cache.discard(folder)
        cache.release(archive)

    if not frames:
        return None
    table = pd.concat(frames, ignore_index=True)

    # One row per LoD2 building part: measHeight is constant across its faces,
    # while the roof area has to be summed over them.
    agg: Dict[str, Any] = {height_field: "max"}
    if area_field and area_field in table.columns:
        table[area_field] = pd.to_numeric(table[area_field], errors="coerce").fillna(0.0)
        agg[area_field] = "sum"
    if ext_field in table.columns:
        agg[ext_field] = "first"
    if func_field in table.columns:
        agg[func_field] = "first"
    table = table.groupby(key_field, as_index=False).agg(agg)

    # Join key: the cadastre object id the model refers to -- in Lower Saxony
    # the ALKIS objid hidden behind "...$$$DENIAL870000ok0W".
    if ext_field and ext_field in table.columns:
        ref = table[ext_field].astype(str)
        table["cadastre_id"] = ref.str.rsplit(ext_split, n=1).str[-1] if ext_split else ref
    else:
        table["cadastre_id"] = ""
    # Where that reference is missing or not a cadastre id, a region may derive
    # the id from the model's own (`fields.id_rewrite`); Lower Saxony's two
    # differ only by namespace prefix (DENILD... vs DENIAL...).
    # `fields.cadastre_id_prefix` says what a valid reference starts with.
    prefix = f.get("cadastre_id_prefix")
    valid = (table["cadastre_id"].str.startswith(str(prefix)) if prefix
             else table["cadastre_id"].astype(bool))
    rewrite = f.get("id_rewrite")
    if rewrite:
        fallback = table[key_field].astype(str).str.replace(
            str(rewrite["pattern"]), str(rewrite["repl"]), regex=True
        )
        table["cadastre_id"] = table["cadastre_id"].where(valid, fallback)
    else:
        table["cadastre_id"] = table["cadastre_id"].where(valid, "")

    table["lod_height"] = pd.to_numeric(table[height_field], errors="coerce")
    if func_field and func_field in table.columns:
        values = table[func_field].astype(str)
        table["lod_function"] = values.str.rsplit(func_split, n=1).str[-1] if func_split else values
    else:
        table["lod_function"] = None
    table = table.dropna(subset=["lod_height"])
    table = table[table["cadastre_id"].astype(bool)]
    if table.empty:
        return None

    # Several LoD2 parts can share one cadastre footprint -- the Neues Rathaus
    # is 84 parts spanning 20-89 m against a handful of ALKIS polygons. Picking
    # any single part (or the maximum) would extrude the whole footprint to the
    # tower, so the parts are reduced with an explicit, configurable statistic.
    mode = str(cfg.get("buildings.height_aggregate", "area_weighted"))
    if area_field:
        table["lod_area"] = pd.to_numeric(table.get(area_field, 0.0), errors="coerce").fillna(0.0)
    else:
        table["lod_area"] = 0.0

    parts_per_building = table.groupby("cadastre_id").size()
    multi = int((parts_per_building > 1).sum())

    if mode == "area_weighted" and table["lod_area"].sum() > 0:
        def _weighted(group: pd.DataFrame) -> float:
            weights = group["lod_area"].to_numpy(dtype="float64")
            heights = group["lod_height"].to_numpy(dtype="float64")
            if weights.sum() <= 0:
                return float(np.nanmax(heights))
            return float(np.average(heights, weights=weights))

        heights = table.groupby("cadastre_id").apply(_weighted, include_groups=False)
        out = heights.rename("lod_height").reset_index()
    elif mode == "median":
        out = table.groupby("cadastre_id", as_index=False)["lod_height"].median()
    else:  # "max"
        if mode not in {"max", "area_weighted"}:
            logger.warning(
                "Unknown buildings.height_aggregate '%s', falling back to 'max'", mode
            )
        out = table.groupby("cadastre_id", as_index=False)["lod_height"].max()

    functions = table.drop_duplicates("cadastre_id")[["cadastre_id", "lod_function"]]
    out = out.merge(functions, on="cadastre_id", how="left")
    logger.info(
        "Building models: %d buildings with a measured height (%d of them assembled "
        "from several model parts, reduced with '%s')", len(out), multi, mode,
    )
    return out


# ---------------------------------------------------------------------------
# building assembly
# ---------------------------------------------------------------------------
def _building_type_code(gfk: Any, mapping: Dict[str, Any]) -> int:
    """Map a cadastre building function onto a PALM ``building_type``.

    PALM only knows residential/office x three construction periods. Few
    cadastres carry a construction year (neither ALKIS nor LoD2 in Lower
    Saxony does), so the period comes from
    ``mapping_*.yml: building_type.default_age`` and is uniform. See the
    assumptions section of the README.
    """
    age = str(mapping.get("default_age", "1951_2000"))
    base = str(mapping.get("default", "residential"))
    by_gfk = mapping.get("by_gfk") or {}
    text = "" if gfk is None else str(gfk)
    best = ""
    for prefix in by_gfk:
        p = str(prefix)
        if text.startswith(p) and len(p) > len(best):
            best, base = p, str(by_gfk[prefix])
    return pt.BUILDING_TYPE[f"{base}_{age}"]


def build_buildings(
    cfg: Config,
    cadastre: Dict[str, Path],
    heights: Optional[pd.DataFrame],
    ndsm_sources: Sequence[Path],
) -> Optional[Path]:
    """Assemble the definitive building layer: footprint + height + type + id.

    Height precedence follows ``buildings.height_source`` / ``height_fallback``:
    the 3D model's measured height where available and plausible, otherwise a
    percentile of the nDSM inside the (inward-buffered) footprint.
    """
    src_path = cadastre.get("buildings")
    if src_path is None:
        logger.warning("No cadastre building layer -- the driver will contain no buildings")
        return None

    name, src = cfg.source_for_role("cadastre")
    spec = (src.get("layers") or {}).get("buildings", {})

    gdf = gpd.read_file(src_path)
    if gdf.empty:
        return None
    if gdf.crs is None:
        gdf = gdf.set_crs(cfg.crs)
    gdf = gdf.to_crs(cfg.crs)
    gdf["geometry"] = gdf.geometry.buffer(0)              # repair self-intersections
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()]

    where = f"sources.{name}.layers.buildings"
    id_name = declared_field(spec, "id_field", where, required=False)
    fn_name = declared_field(spec, "function_field", where, required=False)
    id_field = resolve_field(gdf.columns, id_name) if id_name else None
    fn_field = resolve_field(gdf.columns, fn_name) if fn_name else None
    if id_name and id_field is None:
        logger.warning(
            "Cadastre buildings have no '%s' column; heights cannot be joined "
            "from the 3D model.", id_name,
        )
    gdf["cadastre_id"] = gdf[id_field].astype(str) if id_field else ""
    gdf["building_type"] = [
        _building_type_code(v, cfg.mapping.building_type)
        for v in (gdf[fn_field] if fn_field else [None] * len(gdf))
    ]

    min_h = float(cfg.get("buildings.min_height", 2.5))
    max_h = float(cfg.get("buildings.max_height", 300.0))
    source = str(cfg.get("buildings.height_source", "ndsm"))
    fallback = str(cfg.get("buildings.height_fallback", "none"))

    # The 3D model's function is a useful fallback where the cadastre lacks one,
    # and is taken whether or not the model supplies the height.
    lod = None
    if heights is not None and not heights.empty:
        lod = heights.set_index("cadastre_id")
        if "lod_function" in heights.columns and fn_field:
            missing_fn = gdf[fn_field].isna()
            mapped = gdf.loc[missing_fn, "cadastre_id"].map(lod["lod_function"])
            for idx, value in mapped.dropna().items():
                gdf.at[idx, "building_type"] = _building_type_code(
                    value, cfg.mapping.building_type
                )

    gdf["height"] = np.nan
    gdf["height_cap"] = np.nan
    n_ndsm = n_lod = 0

    def take_lod(mask: "pd.Series") -> int:
        if lod is None or not mask.any():
            return 0
        values = gdf.loc[mask, "cadastre_id"].map(lod["lod_height"])
        ok = values.between(min_h, max_h)
        gdf.loc[values.index[ok], "height"] = values[ok]
        return int(ok.sum())

    def take_ndsm(mask: "pd.Series") -> int:
        # Sample the finest tier first and let coarser ones fill what it does
        # not cover, instead of dropping straight through to the 3D model.
        taken = 0
        # Named `mosaic`, not `source`: the enclosing scope already has a
        # `source` holding the configured height_source ("ndsm"/"lod2"), and a
        # loop variable of the same name here reads like a reassignment of it.
        for mosaic in ndsm_sources:
            todo = mask & gdf["height"].isna()
            if not todo.any():
                break
            stats = sample_polygon_statistic(
                gdf.loc[todo, "geometry"],
                mosaic,
                percentiles=(float(cfg.get("buildings.ndsm_percentile", 75)), 99.0),
                erosion=float(cfg.get("buildings.footprint_erosion", 1.0)),
            )
            got = np.isfinite(stats[0])
            index = gdf.index[todo]
            gdf.loc[index[got], "height"] = stats[0][got]
            gdf.loc[index[got], "height_cap"] = stats[1][got]
            taken += int(got.sum())
        return taken

    def take_field(mask: "pd.Series") -> int:
        """Take a height carried by the cadastre itself, e.g. an OSM tag.

        ``buildings.height_field`` names a numeric column on the cadastre
        building layer (``osm_height`` for the OSM ports). A surveyed or
        tagged height beats any surface-model estimate, so this runs first
        where it is configured.
        """
        # The source spec's `height_field` is what keeps the column alive
        # through fetch_cadastre's clip (it selects every `*_field` value),
        # so honour it as the default and let the run config override.
        field = cfg.get("buildings.height_field") or spec.get("height_field")
        if not field or not mask.any():
            return 0
        column = resolve_field(gdf.columns, str(field))
        if column is None:
            logger.warning(
                "buildings.height_field '%s' is not a column on the cadastre "
                "building layer -- no tagged heights are available", field,
            )
            return 0
        values = pd.to_numeric(gdf.loc[mask, column], errors="coerce")
        ok = values.between(min_h, max_h)
        gdf.loc[values.index[ok], "height"] = values[ok]
        return int(ok.sum())

    n_field = take_field(pd.Series(True, index=gdf.index))

    # "none" disables a measured height source outright. It exists for regions
    # whose only surface model is far coarser than the buildings it would have
    # to measure: over Aswan a 30 m nDSM inside a 100 m2 ML footprint reports
    # the neighbourhood, not the building, and mixing those numbers into the
    # same field as a typological prior makes the field's accuracy
    # unstatable. Better to use the prior alone and say so.
    everything = gdf["height"].isna()
    if source == "lod2":
        n_lod = take_lod(everything)
        if fallback == "ndsm":
            n_ndsm = take_ndsm(gdf["height"].isna())
    elif source == "ndsm":
        n_ndsm = take_ndsm(everything)
        if fallback == "lod2":
            n_lod = take_lod(gdf["height"].isna())
    elif source not in ("none", ""):
        raise ConfigError(
            f"buildings.height_source '{source}' unknown (lod2, ndsm, none)"
        )

    # ---- typological default -------------------------------------------------
    # Last resort for footprints that no height source could measure. WHY THIS
    # EXISTS: over a region whose only surface model is a 30 m DEM, the nDSM
    # inside a 250 m2 footprint is routinely below min_height -- the 30 m cell
    # averages the building with the street beside it. Without a default those
    # footprints are dropped entirely, and 41 % of Aswan's buildings vanished
    # from the 2026-09-07 driver for want of a height rather than for want of a
    # building.
    #
    # A plausible typological height is strictly better for PALM than no
    # building at all: the flow field is dominated by the presence and plan
    # form of the obstacle, and a missing block is a far larger error than a
    # block that is a metre too short. It is, however, a PRIOR AND NOT A
    # MEASUREMENT -- `n_default` is logged separately for exactly that reason,
    # and any driver relying on it must say so.
    #
    # Accepts either a scalar (metres, applied to every building) or a mapping
    # keyed on PALM building_type name:
    #     height_default: 9.0
    #     height_default: {residential_1951_2000: 9.0, office_1951_2000: 12.0}
    n_default = 0
    default_spec = cfg.get("buildings.height_default")
    if default_spec is not None:
        todo = gdf["height"].isna()
        if todo.any():
            if isinstance(default_spec, dict):
                by_code = {}
                for key, value in default_spec.items():
                    code = pt.BUILDING_TYPE.get(str(key))
                    if code is None:
                        raise ValueError(
                            f"buildings.height_default: unknown building_type "
                            f"'{key}' (known: {sorted(pt.BUILDING_TYPE)})"
                        )
                    by_code[code] = float(value)
                values = gdf.loc[todo, "building_type"].map(by_code)
            else:
                values = pd.Series(float(default_spec), index=gdf.index[todo])
            ok = values.between(min_h, max_h)
            gdf.loc[values.index[ok], "height"] = values[ok]
            n_default = int(ok.sum())

    before = len(gdf)
    gdf = gdf[gdf["height"].between(min_h, max_h)].copy()
    # The cap bounds per-pixel heights later on; keep it above the robust height.
    gdf["height_cap"] = np.where(
        np.isfinite(gdf["height_cap"]),
        np.maximum(gdf["height_cap"], gdf["height"]),
        gdf["height"],
    )
    logger.info(
        "Buildings: %d kept (%d from cadastre height tags, %d from the surface "
        "model, %d from the 3D model, %d from the typological default, "
        "%d dropped as implausible or too low)",
        len(gdf), n_field, n_ndsm, n_lod, n_default, before - len(gdf),
    )
    if n_default:
        logger.warning(
            "%d of %d buildings (%.0f %%) carry buildings.height_default, not a "
            "measured height -- their heights are a typological prior",
            n_default, len(gdf), 100.0 * n_default / max(len(gdf), 1),
        )
    if gdf.empty:
        return None

    # palm_csd needs an int32 building_id; the cadastre id is a string.
    gdf = gdf.reset_index(drop=True)
    gdf["building_id"] = np.arange(1, len(gdf) + 1, dtype="int32")

    out = cfg.workdir / "vector" / "buildings.gpkg"
    if out.exists():
        out.unlink()
    gdf[
        ["building_id", "building_type", "height", "height_cap", "cadastre_id", "geometry"]
    ].to_file(out, driver="GPKG")
    return out


def build_bridges(cfg: Config, cadastre: Dict[str, Path], ndsm_sources: Sequence[Path]) -> Optional[Path]:
    """Extract bridge decks and give them a height above terrain."""
    if not cfg.get("buildings.bridges.enabled", True):
        return None
    src_path = cadastre.get("structures")
    if src_path is None:
        return None
    name, src = cfg.source_for_role("cadastre")
    spec = (src.get("layers") or {}).get("structures", {})
    type_field = declared_field(spec, "type_field", f"sources.{name}.layers.structures")
    values = {str(v) for v in (cfg.mapping.bridges.get("values") or spec.get("bridge_values") or [])}
    if not values:
        return None

    gdf = gpd.read_file(src_path)
    if gdf.empty:
        return None
    resolved = resolve_field(gdf.columns, type_field)
    if resolved is None:
        logger.warning("Bridges: no '%s' column in %s", type_field, src_path.name)
        return None
    gdf = gdf[gdf[resolved].astype(str).str.strip().isin(values)].copy()
    if gdf.empty:
        logger.info("Bridges: none in the AOI")
        return None
    if gdf.crs is None:
        gdf = gdf.set_crs(cfg.crs)
    gdf["geometry"] = gdf.geometry.buffer(0)
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()].reset_index(drop=True)

    if not ndsm_sources:
        logger.warning("Bridges found but no nDSM available -- skipping bridges")
        return None
    gdf["height"] = np.nan
    for source in ndsm_sources:
        todo = gdf["height"].isna()
        if not todo.any():
            break
        sampled = sample_polygon_statistic(
            gdf.loc[todo, "geometry"], source, percentiles=(50.0,), erosion=0.5
        )[0]
        got = np.isfinite(sampled)
        gdf.loc[gdf.index[todo][got], "height"] = sampled[got]
    gdf = gdf[gdf["height"].between(0.5, 100.0)].reset_index(drop=True)
    if gdf.empty:
        return None
    gdf["bridge_id"] = np.arange(1, len(gdf) + 1, dtype="int32")

    out = cfg.workdir / "vector" / "bridges.gpkg"
    if out.exists():
        out.unlink()
    gdf[["bridge_id", "height", "geometry"]].to_file(out, driver="GPKG")
    logger.info("Bridges: %d deck(s)", len(gdf))
    return out


def sample_polygon_statistic(
    geometries: gpd.GeoSeries,
    raster: Path,
    percentiles: Sequence[float] = (75.0,),
    erosion: float = 1.0,
) -> np.ndarray:
    """Percentiles of ``raster`` inside each polygon, after an inward buffer.

    The inward buffer keeps the mixed edge pixels of a roof -- half roof, half
    ground -- out of the statistic. Polygons that vanish under the buffer (small
    sheds, narrow wings) are sampled without it rather than dropped.

    Returns an array of shape ``(len(percentiles), len(geometries))``; several
    percentiles are read in one pass because the window read dominates the cost.
    """
    import rasterio
    from rasterio.features import geometry_mask
    from rasterio.windows import Window, from_bounds

    result = np.full((len(percentiles), len(geometries)), np.nan, dtype="float64")
    with rasterio.open(raster) as src:
        nodata = src.nodata
        for i, geom in enumerate(geometries.reset_index(drop=True)):
            if geom is None or geom.is_empty:
                continue
            shrunk = geom.buffer(-erosion) if erosion > 0 else geom
            if shrunk.is_empty:
                shrunk = geom
            b = shrunk.bounds
            if b[2] <= src.bounds.left or b[0] >= src.bounds.right:
                continue
            if b[3] <= src.bounds.bottom or b[1] >= src.bounds.top:
                continue
            window = from_bounds(*b, transform=src.transform).round_offsets().round_lengths()
            if window.width < 1 or window.height < 1:
                window = Window(window.col_off, window.row_off, 1, 1)
            # Clip to the dataset and paste into a NaN canvas by hand rather
            # than asking for a boundless read. rasterio implements boundless
            # by wrapping the dataset in a temporary VRT per call, and wrapping
            # a mosaic is priced by its source count, not by the window: over
            # a 2500-tile city-scale root that is 39.7 ms per polygon against
            # 0.16 ms here -- 9 hours versus 2 minutes for 663 000 footprints.
            col_off, row_off = int(window.col_off), int(window.row_off)
            w, h = int(window.width), int(window.height)
            c0, r0 = max(0, col_off), max(0, row_off)
            c1, r1 = min(src.width, col_off + w), min(src.height, row_off + h)
            if c1 <= c0 or r1 <= r0:
                continue
            data = np.full((h, w), np.nan, dtype="float32")
            data[r0 - row_off:r1 - row_off, c0 - col_off:c1 - col_off] = src.read(
                1, window=Window(c0, r0, c1 - c0, r1 - r0)
            )
            transform = src.window_transform(window)
            mask = geometry_mask(
                [shrunk], out_shape=data.shape, transform=transform, invert=True,
                all_touched=data.size <= 4,
            )
            values = data[mask]
            if nodata is not None:
                values = values[values != nodata]
            values = values[np.isfinite(values)]
            if values.size:
                result[:, i] = np.percentile(values, percentiles)
    return result


# ---------------------------------------------------------------------------
# topographic vectors (roads)
# ---------------------------------------------------------------------------
def fetch_roads(cfg: Config, cache: Cache) -> Optional[Path]:
    """Download, clip and classify road centrelines for ``street_type``."""
    if not cfg.get("sources.street_type.enabled", False):
        return None
    name, src = cfg.source_for_role("topo")
    aoi = _aoi(cfg)
    outdir = cfg.workdir / "vector"
    spec = (src.get("layers") or {}).get("roads") or {}
    layer_name = declared_field(spec, "name", f"sources.{name}.layers.roads")

    url = src.get("url")
    if not url:
        logger.warning("Topo source '%s' has no url -- skipping street_type", name)
        return None
    hint = src.get("size_hint")
    if hint:
        logger.info(
            "Topo source '%s': downloading %.0f MB once; it is clipped and deleted "
            "immediately afterwards.", name, float(hint) / 1e6,
        )
    archive = cache.fetch(url, subdir=name)
    out = outdir / "roads.gpkg"
    with cache.hold(archive):
        folder = cache.unzip(archive)
        matches = [p for p in folder.rglob("*.shp") if p.stem.lower() == layer_name.lower()]
        if not matches:
            available = sorted({p.stem for p in folder.rglob("*.shp")})
            logger.warning(
                "Layer '%s' not found in %s. Available: %s", layer_name, archive.name, available
            )
            cache.discard(folder)
            cache.release(archive)
            return None
        clip_layer(matches[0], matches[0].stem, aoi, out, crs=cfg.crs)
        cache.discard(folder)
    cache.release(archive)
    if not out.exists():
        return None

    configured = cfg.mapping.street_type.get("field") or spec.get("class_field")
    if not configured:
        raise ConfigError(
            f"street_type: name the road class column, as `street_type.field` in "
            f"the mapping or `class_field` under sources.{name}.layers.roads"
        )
    configured = str(configured)
    by_value = {str(k): str(v) for k, v in (cfg.mapping.street_type.get("by_value") or {}).items()}
    default = str(cfg.mapping.street_type.get("default", "unclassified"))

    gdf = gpd.read_file(out)
    if gdf.empty:
        return None
    field = resolve_field(gdf.columns, configured)
    if field is None:
        logger.warning(
            "street_type: no column '%s' in %s (available: %s). Everything becomes "
            "'%s' -- fix mapping.street_type.field.",
            configured, out.name, ", ".join(sorted(gdf.columns)[:20]), default,
        )
        raw = pd.Series([""] * len(gdf))
    else:
        raw = gdf[field].astype(str).str.strip()
    unmapped = sorted(set(raw) - set(by_value))
    if unmapped:
        logger.warning(
            "street_type: %d unmapped '%s' value(s), classified as '%s': %s",
            len(unmapped), configured, default, ", ".join(v or "<empty>" for v in unmapped[:15]),
        )
    gdf["street_type"] = [
        pt.STREET_TYPE[by_value.get(v, default)] for v in raw
    ]
    # PALM wants an area, not a centreline: buffer by half the carriageway width.
    # Optional: a source without a carriageway width gets the class default.
    width_name = spec.get("width_field")
    width_field = resolve_field(gdf.columns, str(width_name)) if width_name else None
    if width_field:
        # the Basis-DLM uses -9998 as "unknown"; any implausible value falls
        # through to the class default below
        width = pd.to_numeric(gdf[width_field], errors="coerce").fillna(0.0)
    else:
        width = pd.Series(np.zeros(len(gdf)))
    width = width.where(width > 1.0, _default_width(gdf["street_type"]))
    centrelines = gdf.geometry.copy()
    gdf["geometry"] = gdf.geometry.buffer(width.to_numpy() / 2.0, cap_style=2)

    out2 = outdir / "street_type.gpkg"
    if out2.exists():
        out2.unlink()
    gdf[["street_type", "geometry"]].to_file(out2, driver="GPKG")
    logger.info("Roads: %d segment(s) classified", len(gdf))

    _build_crossings(centrelines, outdir / "street_crossings.gpkg")
    out.unlink(missing_ok=True)
    return out2


def _build_crossings(centrelines: gpd.GeoSeries, out: Path, radius: float = 6.0) -> Optional[Path]:
    """Derive junction areas from the road network.

    PALM's ``street_crossings`` marks places a pedestrian agent may cross. A
    junction is approximated as a node shared by three or more centreline ends,
    buffered by ``radius``. This only feeds the multi-agent model, so a coarse
    approximation is adequate; supply a proper crossing data set through
    ``extra_inputs`` if the agent model matters for the study.
    """
    from collections import Counter

    from shapely.geometry import Point

    ends: Counter = Counter()
    for geom in centrelines:
        if geom is None or geom.is_empty:
            continue
        parts = geom.geoms if geom.geom_type.startswith("Multi") else [geom]
        for part in parts:
            coords = list(part.coords)
            if len(coords) < 2:
                continue
            for c in (coords[0], coords[-1]):
                ends[(round(c[0], 2), round(c[1], 2))] += 1
    nodes = [Point(x, y) for (x, y), n in ends.items() if n >= 3]
    if not nodes:
        return None
    gdf = gpd.GeoDataFrame(
        {"crossing": np.ones(len(nodes), dtype="uint8")},
        geometry=[p.buffer(radius) for p in nodes],
        crs=centrelines.crs,
    )
    if out.exists():
        out.unlink()
    gdf.to_file(out, driver="GPKG")
    logger.info("Street crossings: %d junction(s)", len(gdf))
    return out


def _default_width(street_type: pd.Series) -> pd.Series:
    """Fallback carriageway widths [m] where the source has none."""
    table = {
        pt.STREET_TYPE["motorway"]: 24.0,
        pt.STREET_TYPE["trunk"]: 14.0,
        pt.STREET_TYPE["primary"]: 12.0,
        pt.STREET_TYPE["secondary"]: 10.0,
        pt.STREET_TYPE["tertiary"]: 8.0,
        pt.STREET_TYPE["residential"]: 6.0,
        pt.STREET_TYPE["service"]: 4.0,
        pt.STREET_TYPE["unclassified"]: 5.0,
    }
    return street_type.map(table).fillna(5.0)
