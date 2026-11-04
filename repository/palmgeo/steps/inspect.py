"""Inspection products for QGIS.

The domain rasters are what palm_csd consumes, but they are already aggregated
onto the PALM grids and are therefore not the right artefact for judging the
classification itself. This stage writes the *unaggregated* work-resolution
surface classification as one single-band raster whose values read directly:

===========  ============================================
value        meaning
===========  ============================================
``0``        undefined
``1000+v``   vegetation, ``v`` = PALM ``vegetation_type``
``2000+p``   pavement, ``p`` = PALM ``pavement_type``
``3000+w``   water, ``w`` = PALM ``water_type``
``4000``     building
===========  ============================================

A QGIS layer style with the full legend is written next to it, so the file can
be dropped into a project and read without a lookup table. Polygonisation is
available but off by default -- a 20 x 20 km area at 1 m easily exceeds ten
million polygons.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import rasterio
from rasterio.enums import Resampling

from palmgeo import palm_types as pt
from palmgeo.config import Config
from palmgeo.grid import Tier
from palmgeo.rasters import build_vrt, gdal_env

logger = logging.getLogger("palmgeo.inspect")

#: A colour per surface group, varied by type so classes stay distinguishable.
GROUP_COLOURS: Dict[int, Tuple[int, int, int]] = {
    pt.GROUP_VEGETATION: (86, 148, 62),
    pt.GROUP_PAVEMENT: (140, 140, 145),
    pt.GROUP_WATER: (62, 122, 186),
    pt.GROUP_BUILDING: (176, 96, 84),
}


def _legend() -> List[Tuple[int, str, Tuple[int, int, int]]]:
    entries: List[Tuple[int, str, Tuple[int, int, int]]] = [(0, "undefined", (0, 0, 0))]
    for group, table in (
        (pt.GROUP_VEGETATION, pt.VEGETATION_TYPE),
        (pt.GROUP_PAVEMENT, pt.PAVEMENT_TYPE),
        (pt.GROUP_WATER, pt.WATER_TYPE),
    ):
        base = np.array(GROUP_COLOURS[group], dtype=float)
        values = sorted(set(table.values()))
        for i, value in enumerate(values):
            name = next(k for k, v in table.items() if v == value)
            # spread lightness across the group's members
            factor = 0.62 + 0.7 * (i / max(1, len(values) - 1))
            rgb = tuple(int(min(255, max(0, c * factor))) for c in base)
            entries.append((pt.INSPECTION_OFFSET[group] + value, f"{pt.NAME_OF_GROUP[group]}: {name}", rgb))
    entries.append((pt.INSPECTION_OFFSET[pt.GROUP_BUILDING], "building", GROUP_COLOURS[pt.GROUP_BUILDING]))
    return entries


def write_qml(path: Path) -> Path:
    """Write a QGIS paletted-raster style beside the inspection raster."""
    items = "\n".join(
        f'        <paletteEntry value="{value}" alpha="255" '
        f'color="#{r:02x}{g:02x}{b:02x}" label="{label}"/>'
        for value, label, (r, g, b) in _legend()
        if value != 0
    )
    xml = f"""<!DOCTYPE qgis>
<qgis version="3.28" styleCategories="AllStyleCategories">
  <pipe>
    <rasterrenderer type="paletted" band="1" opacity="1" alphaBand="-1" nodataColor="">
      <rasterTransparency/>
      <colorPalette>
{items}
      </colorPalette>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0"/>
    <huesaturation saturation="0" grayscaleMode="0"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
"""
    path.write_text(xml, encoding="utf-8")
    return path


def build_inspection(cfg: Config, tier: Tier) -> Optional[Path]:
    """Mosaic the composed ``surface_class`` tiles into one browsable COG."""
    parts = sorted((cfg.tier_dir(tier) / "surface").glob("*/surface_class.tif"))
    if not parts:
        logger.warning("Tier %s has no composed tiles -- nothing to inspect", tier.name)
        return None

    # These layers exist to be looked at while building a driver, not to be
    # shipped with it -- they are neither a palm_csd input nor documentation.
    # They therefore live in the working directory, which is not published.
    out_dir = cfg.workdir / "inspect"
    out_dir.mkdir(parents=True, exist_ok=True)
    vrt = build_vrt(parts, out_dir / f"surface_{tier.name}.vrt", nodata=0)
    out = out_dir / f"surface_{tier.name}.tif"

    opts = cfg.get("output.cog", {}) or {}
    with rasterio.open(vrt) as src:
        profile = src.profile.copy()
        profile.update(
            driver="GTiff",
            compress=opts.get("compress", "DEFLATE"),
            predictor=2,
            tiled=True,
            blockxsize=512,
            blockysize=512,
            BIGTIFF="IF_SAFER",
            nodata=0,
        )
        with rasterio.open(out, "w", **profile) as dst:
            # copy block-wise: a full-state tier can be tens of gigapixels
            for _, window in src.block_windows(1):
                dst.write(src.read(1, window=window), 1, window=window)
            if min(src.width, src.height) >= 1024:
                dst.build_overviews([2, 4, 8, 16, 32], Resampling.mode)
                dst.update_tags(ns="rio_overview", resampling="mode")

    write_qml(out.with_suffix(".qml"))
    logger.info("Inspection raster: %s (+ .qml style)", out)
    return out


def polygonize(cfg: Config, raster: Path) -> Optional[Path]:
    """Vectorise the inspection raster into a GeoPackage of typed polygons."""
    from rasterio.features import shapes as rio_shapes
    import geopandas as gpd
    from shapely.geometry import shape as to_shape

    out = raster.with_suffix(".gpkg")
    if out.exists():
        out.unlink()

    label_of = {value: label for value, label, _ in _legend()}
    written = 0
    first = True
    with rasterio.open(raster) as src:
        total = sum(1 for _ in src.block_windows(1))
        logger.info("Polygonising %d block(s) of %s", total, raster.name)
        for i, (_, window) in enumerate(src.block_windows(1), start=1):
            data = src.read(1, window=window)
            if not data.any():
                continue
            transform = src.window_transform(window)
            records = [
                {"code": int(value), "geometry": to_shape(geom)}
                for geom, value in rio_shapes(data, mask=data != 0, transform=transform)
            ]
            if not records:
                continue
            gdf = gpd.GeoDataFrame(records, crs=src.crs)
            gdf["group"] = [pt.NAME_OF_GROUP.get(c // 1000, "?") for c in gdf["code"]]
            gdf["type"] = [label_of.get(c, "?") for c in gdf["code"]]
            gdf["palm_value"] = [c % 1000 for c in gdf["code"]]
            gdf.to_file(out, driver="GPKG", layer="surface", mode="w" if first else "a")
            first = False
            written += len(gdf)
            if i % 50 == 0:
                logger.info("  block %d/%d, %d polygon(s)", i, total, written)
    if written == 0:
        return None
    logger.info(
        "Polygon layer: %s (%d polygons; note that these are block-wise and not "
        "dissolved across block borders)", out, written,
    )
    return out


def run(cfg: Config) -> Dict[str, Any]:
    """Build the inspection products for every tier."""
    result: Dict[str, Any] = {"rasters": [], "polygons": []}
    if not cfg.get("output.inspection_raster", True):
        return result
    with gdal_env(int(cfg.runtime.get("gdal_cachemax_mb", 512))):
        for tier in cfg.tiers:
            raster = build_inspection(cfg, tier)
            if raster is None:
                continue
            result["rasters"].append(str(raster))
            if cfg.get("output.inspection_polygons", False):
                poly = polygonize(cfg, raster)
                if poly:
                    result["polygons"].append(str(poly))
    return result
