"""Generation of the palm_csd configuration file.

palm_csd validates its YAML eagerly with pydantic and refuses to start if a
listed input file does not exist, so only products that were actually written
are referenced. Two details of palm_csd 25.10 are easy to get wrong because the
shipped documentation still describes the previous interface:

* the domain switch is ``generate_buildings_3d``; plain ``buildings_3d`` was
  removed and now raises;
* the ``file_*`` input keys were removed in favour of the ``files:`` mapping.

Both are verified against ``palm_csd/csd_config.py`` of the release rather than
against ``docs/yaml.md``.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import rasterio
import yaml

from palmgeo.config import Config
from palmgeo.grid import Domain
from palmgeo.steps.domains import CSD_FILES

logger = logging.getLogger("palmgeo.csd")

#: Domain keys palm_csd accepts; anything else in a domain block is dropped
#: with a warning rather than making palm_csd fail on an unknown field
#: (its models use ``extra="forbid"``).
CSD_DOMAIN_KEYS = {
    "pixel_size", "nx", "ny", "dz", "origin_x", "origin_y", "origin_lon", "origin_lat",
    "lower_left_x", "lower_left_y", "input_lower_left_x", "input_lower_left_y",
    "domain_parent", "input", "bridge_depth", "building_free_border_width",
    "building_free_border_pavement_type", "generate_buildings_3d",
    "estimate_lai_from_vegetation_height", "generate_single_trees",
    "generate_vegetation_patches", "overhanging_trees", "remove_low_lai_tree",
    "replace_high_vegetation_types", "interpolate_terrain", "use_palm_z_axis",
    "lcz_input", "dcep", "z_uhl", "udir", "lai_per_vegetation_height",
    "building_albedo_type", "building_emissivity", "building_fraction",
    "building_general_pars", "building_heat_capacity", "building_heat_conductivity",
    "building_indoor_pars", "building_lai", "building_roughness_length",
    "building_roughness_length_qh", "building_thickness", "building_transmissivity",
    "water_temperature",
}

CSD_SETTINGS_KEYS = {
    "epsg", "rotation_angle", "season", "downscaling_method", "upscaling_method",
    "ignore_input_georeferencing", "height_high_vegetation_lower_threshold",
    "height_rel_resolved_vegetation_lower_threshold", "lai_high_vegetation_default",
    "lai_low_vegetation_default", "lai_tree_lower_threshold", "lad_method",
    "lad_alpha", "lad_beta", "lad_z_max_rel", "patch_height_default",
    "replace_invalid_input_values", "soil_type_default", "use_lai_for_roofs",
    "use_lai_for_trees", "use_vegetation_height_for_trees", "vegetation_type_below_trees",
    "lai_per_vegetation_height",
}


def _has_data(path: Path) -> bool:
    """True if the raster holds at least one unmasked cell."""
    try:
        with rasterio.open(path) as src:
            return int(src.read(1, masked=True).count()) > 0
    except Exception:  # noqa: BLE001
        return False


def tree_source(cfg: Config) -> Tuple[Optional[Path], Dict[str, str]]:
    """Resolve where single trees come from, and how their columns are named.

    ``rentz`` writes palm_csd's own column names, so no mapping is needed;
    ``external`` data keeps whatever the supplier used and is translated by
    ``extra_inputs.trees.columns``.
    """
    mode = str(cfg.get("vegetation.trees.mode", "rentz"))
    if mode == "unresolved":
        return None, {}
    if mode == "rentz":
        path = cfg.workdir / "vector" / "trees.gpkg"
        if not path.exists():
            logger.warning(
                "vegetation.trees.mode is 'rentz' but %s does not exist -- "
                "run the `trees` stage first", path,
            )
            return None, {}
        return path, {}
    spec = cfg.get("extra_inputs.trees") or {}
    src = spec.get("path")
    if not src:
        return None, {}
    path = (cfg.root / str(src)).resolve()
    if not path.exists():
        logger.warning("extra_inputs.trees.path %s does not exist -- no single trees", path)
        return None, {}
    columns = {str(k): str(v) for k, v in (spec.get("columns") or {}).items()}
    return path, columns


def prepare_trees(cfg: Config, domain: Domain, src_path: Optional[Path]) -> Optional[str]:
    """Copy the tree data set into the domain's input directory, clipped.

    palm_csd prepends the section's ``path`` to every file name, so the tree
    file has to live underneath it. It is copied rather than referenced through
    ``../`` so that a domain directory stays self-contained and can be shipped
    to the cluster on its own.
    """
    if src_path is None:
        return None
    out_dir = cfg.domain_dir(domain)
    layer = (cfg.get("extra_inputs.trees") or {}).get("layer")
    target = out_dir / f"trees{src_path.suffix or '.gpkg'}"
    if src_path.suffix.lower() == ".shp":
        # a shapefile is a family of files
        for sidecar in src_path.parent.glob(src_path.stem + ".*"):
            shutil.copy2(sidecar, out_dir / f"trees{sidecar.suffix}")
        return target.name

    import geopandas as gpd

    gdf = gpd.read_file(src_path, layer=layer) if layer else gpd.read_file(src_path)
    if gdf.crs is None:
        gdf = gdf.set_crs(cfg.crs)
    gdf = gdf.to_crs(cfg.crs)
    # trim to the domain so every nest only carries its own trees
    xmin, ymin, xmax, ymax = domain.bounds
    gdf = gdf.cx[xmin:xmax, ymin:ymax]
    if gdf.empty:
        logger.warning("No trees inside domain %s", domain.name)
        return None
    if target.exists():
        target.unlink()
    gdf.to_file(target, driver="GPKG")
    logger.info("Domain %s: %d single tree(s)", domain.name, len(gdf))
    return target.name


def _config_path(cfg: Config, target: Path) -> str:
    """How a path is written into the generated palm_csd configuration.

    Relative to the project root, not absolute. palm_csd resolves what it finds
    here with a plain ``os.path.join`` against its own working directory
    (``csd_config.py``), so a relative path means "run palm_csd from the
    project root" -- which is what the README and static_driver/README have
    always documented, and what the printed command line does.

    Absolute would be robust against the working directory and wrong in every
    other respect. The generated configuration lives in the project's
    ``config/`` directory and is part of what a published driver carries, so an
    absolute path there is machine-specific -- a cluster has no
    ``/home/<somebody>`` -- and names the author's account in an artefact meant
    to be shared. A relative path that is resolved from the wrong directory
    fails loudly with a missing file; an absolute one that no longer exists
    fails the same way, having leaked a username first.

    A target outside the project tree cannot be expressed relatively and is
    written absolute, with a warning, because that is a deliberate choice
    somebody made by pointing ``outdir`` elsewhere.
    """
    try:
        return str(target.resolve().relative_to(cfg.root.resolve()))
    except ValueError:
        logger.warning(
            "%s lies outside the project root %s, so the generated palm_csd "
            "configuration has to name it absolutely. That file is no longer "
            "portable to another machine.", target, cfg.root,
        )
        return str(target.resolve())


def build(cfg: Config) -> Path:
    """Write ``<outdir>/config/<name>_csd.yml`` and return its path."""
    attributes = dict(cfg.get("csd.attributes", {}) or {})
    settings: Dict[str, Any] = {"epsg": cfg.epsg}
    for key, value in (cfg.get("csd.settings", {}) or {}).items():
        if key in CSD_SETTINGS_KEYS:
            settings[key] = value
        else:
            logger.warning("csd.settings: dropping unknown key '%s'", key)

    doc: Dict[str, Any] = {
        "attributes": attributes,
        "settings": settings,
        "output": {
            "path": _config_path(cfg, cfg.outdir),
            "file_out": str(cfg.get("project.name", "palm")),
            "version": 1,
        },
    }

    tree_path, tree_columns = tree_source(cfg)

    missing_report: Dict[str, List[str]] = {}
    for domain in cfg.domains:
        out_dir = cfg.domain_dir(domain)
        files: Dict[str, str] = {}
        missing: List[str] = []
        # palm_csd ignores bridges whose depth is below half the vertical grid
        # spacing, then warns about an unused input. Leave them out instead.
        bridge_depth = float(cfg.get("buildings.bridges.depth", 3.0))
        skip_bridges = bridge_depth < 0.5 * domain.dz
        for name in CSD_FILES:
            if skip_bridges and name.startswith("bridges"):
                continue
            path = out_dir / f"{name}.tif"
            if path.exists() and _has_data(path):
                files[name] = path.name
            else:
                missing.append(name)
        if skip_bridges and (out_dir / "bridges_2d.tif").exists():
            logger.info(
                "Domain %s: bridges omitted, bridge_depth %.1f m < dz/2 = %.1f m "
                "(palm_csd would discard them anyway)",
                domain.name, bridge_depth, 0.5 * domain.dz,
            )
        if "zt" not in files:
            raise RuntimeError(
                f"Domain {domain.name}: zt.tif is missing or empty. "
                f"Run the `domains` stage first."
            )
        if missing:
            missing_report[domain.name] = missing

        trees = (
            prepare_trees(cfg, domain, tree_path)
            if domain.options.get("generate_single_trees")
            else None
        )
        if trees:
            files["trees"] = trees
        elif domain.options.get("generate_single_trees"):
            logger.warning(
                "Domain %s asks for single trees but no tree data is available; "
                "disabling generate_single_trees for this domain.", domain.name,
            )
            domain.options["generate_single_trees"] = False

        section: Dict[str, Any] = {"path": _config_path(cfg, out_dir), "files": files}
        if trees and tree_columns:
            section["columns"] = tree_columns
        doc[f"input_{domain.name}"] = section

        block: Dict[str, Any] = {
            "pixel_size": float(domain.pixel_size),
            "nx": int(domain.nx),
            "ny": int(domain.ny),
            "dz": float(domain.dz),
            "origin_x": float(domain.origin_x),
            "origin_y": float(domain.origin_y),
            "input": domain.name,
        }
        if domain.domain_parent:
            block["domain_parent"] = domain.domain_parent
        if "bridges_2d" in files:
            block["bridge_depth"] = float(cfg.get("buildings.bridges.depth", 3.0))
        for key, value in domain.options.items():
            if key == "vegetation_mode":
                continue
            if key in CSD_DOMAIN_KEYS:
                block[key] = value
            else:
                logger.warning(
                    "Domain %s: dropping key '%s', which palm_csd does not accept",
                    domain.name, key,
                )
        doc[f"domain_{domain.name}"] = block

    for name, fields in missing_report.items():
        logger.info("Domain %s: no data for %s", name, ", ".join(fields))

    # The configuration is a description of the run and lives in config/; the
    # drivers palm_csd writes from it stay at the top of the project directory,
    # which is why output.path above is cfg.outdir and not this.
    cfg.outdir.mkdir(parents=True, exist_ok=True)
    cfg.configdir.mkdir(parents=True, exist_ok=True)
    out_path = cfg.configdir / f"{cfg.get('project.name', 'palm')}_csd.yml"
    base = str(cfg.get("project.name", "palm"))
    root = next((d for d in cfg.domains if d.domain_parent is None), None)
    header = (
        "# Generated by palmgeo -- do not edit by hand, edit the run config instead.\n"
        f"# Source run config: {_config_path(cfg, cfg.path)}\n"
        f"# Domains: {', '.join(d.name for d in cfg.domains)}\n"
        "#\n"
        "# Run with:  env -u PROJ_DATA -u PROJ_LIB palm_csd --pdf " + out_path.name + "\n"
    )
    if root is not None:
        # palm_csd names every output <file_out>_<domain>; PALM reads the root
        # driver without a suffix, because the suffix is what marks a nest.
        header += (
            f"# Then:      python repository/tools/finalize.py {_config_path(cfg, cfg.path)}"
            f"    # {base}_{root.name} -> {base}\n"
        )
    out_path.write_text(
        header + yaml.safe_dump(doc, sort_keys=False, allow_unicode=True, default_flow_style=False),
        encoding="utf-8",
    )
    logger.info("Wrote %s", out_path)
    return out_path
