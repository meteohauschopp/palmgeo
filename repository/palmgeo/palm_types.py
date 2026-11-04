"""PALM surface type vocabularies.

The tables mirror ``palm_csd/docs/types.md`` and the valid ranges in
``palm_csd/data/value_defaults.csv`` of the PALM model system. They are the
single source of truth for translating the human readable names used in
``mapping_*.yml`` into the integer codes written to the GeoTIFFs.
"""

from __future__ import annotations

from typing import Dict, Tuple

VEGETATION_TYPE: Dict[str, int] = {
    "user_defined": 0,
    "bare_soil": 1,
    "crops_mixed_farming": 2,
    "short_grass": 3,
    "evergreen_needleleaf_trees": 4,
    "deciduous_needleleaf_trees": 5,
    "evergreen_broadleaf_trees": 6,
    "deciduous_broadleaf_trees": 7,
    "tall_grass": 8,
    "desert": 9,
    "tundra": 10,
    "irrigated_crops": 11,
    "semidesert": 12,
    "ice_caps_glaciers": 13,
    "bogs_marshes": 14,
    "evergreen_shrubs": 15,
    "deciduous_shrubs": 16,
    "mixed_forest_woodland": 17,
    "interrupted_forest": 18,
}

PAVEMENT_TYPE: Dict[str, int] = {
    "asphalt_concrete_mix": 1,
    "asphalt": 2,
    "concrete": 3,
    "sett": 4,
    "paving_stones": 5,
    "cobblestone": 6,
    "metal": 7,
    "wood": 8,
    "gravel": 9,
    "fine_gravel": 10,
    "pebblestone": 11,
    "woodchips": 12,
    "tartan": 13,
    "artificial_turf": 14,
    "clay": 15,
}

WATER_TYPE: Dict[str, int] = {
    "lake": 1,
    "river": 2,
    "ocean": 3,
    "pond": 4,
    "fountain": 5,
}

SOIL_TYPE: Dict[str, int] = {
    "coarse": 1,
    "medium": 2,
    "medium_fine": 3,
    "fine": 4,
    "very_fine": 5,
    "organic": 6,
}

BUILDING_TYPE: Dict[str, int] = {
    "residential_1950": 1,
    "residential_1951_2000": 2,
    "residential_2001": 3,
    "office_1950": 4,
    "office_1951_2000": 5,
    "office_2001": 6,
    "bridges": 7,
}

STREET_TYPE: Dict[str, int] = {
    "unclassified": 1,
    "cycleway": 2,
    "footway_pedestrian": 3,
    "path": 4,
    "track": 5,
    "living_street": 6,
    "service": 7,
    "residential": 8,
    "tertiary": 9,
    "tertiary_link": 10,
    "secondary": 11,
    "secondary_link": 12,
    "primary": 13,
    "primary_link": 14,
    "trunk": 15,
    "trunk_link": 16,
    "motorway": 17,
    "motorway_link": 18,
    "raceway": 19,
}

#: palm_csd treats these vegetation types as high (grown) vegetation with large
#: roughness lengths (docs/types.md, section "High vegetation").
HIGH_VEGETATION = frozenset(
    VEGETATION_TYPE[n]
    for n in (
        "evergreen_needleleaf_trees",
        "deciduous_needleleaf_trees",
        "evergreen_broadleaf_trees",
        "deciduous_broadleaf_trees",
        "mixed_forest_woodland",
        "interrupted_forest",
    )
)

#: Largest roughness length PALM's land-surface model assigns to any of the
#: HIGH_VEGETATION types: 2.0 m for types 4-7 and 17, 1.1 m for type 18
#: (``vegetation_pars`` in land_surface_model_mod.f90, PALM 25.10). PALM stops
#: with LSM0048 when z0 >= 0.5 * z_mo, and z_mo, the first grid level, is dz/2;
#: so a domain that keeps these types needs dz > 4 * HIGH_VEGETATION_Z0.
HIGH_VEGETATION_Z0 = 2.0

#: Valid ranges from palm_csd/data/value_defaults.csv. Values outside are
#: rejected at config load time rather than at palm_csd runtime.
VALID_RANGE: Dict[str, Tuple[int, int]] = {
    "vegetation_type": (1, 18),
    "pavement_type": (0, 15),
    "water_type": (0, 5),
    "soil_type": (1, 6),
    "building_type": (1, 7),
    "street_type": (1, 19),
}

TABLES: Dict[str, Dict[str, int]] = {
    "vegetation": VEGETATION_TYPE,
    "pavement": PAVEMENT_TYPE,
    "water": WATER_TYPE,
    "soil": SOIL_TYPE,
    "building": BUILDING_TYPE,
    "street": STREET_TYPE,
}

#: Surface groups as encoded in the ``surface_group`` work raster. The order is
#: the precedence palm_csd applies in create_driver.py::process_surface_types:
#: building > water > pavement > vegetation.
GROUP_VEGETATION = 1
GROUP_PAVEMENT = 2
GROUP_WATER = 3
GROUP_BUILDING = 4
GROUP_NODATA = 0

GROUP_OF_NAME = {
    "vegetation": GROUP_VEGETATION,
    "pavement": GROUP_PAVEMENT,
    "water": GROUP_WATER,
    "building": GROUP_BUILDING,
}
NAME_OF_GROUP = {v: k for k, v in GROUP_OF_NAME.items()}

#: Offsets used by the single-band inspection raster so that a human (or QGIS)
#: can read group and type off one Int16 value: 1003 = vegetation/short_grass.
INSPECTION_OFFSET = {
    GROUP_VEGETATION: 1000,
    GROUP_PAVEMENT: 2000,
    GROUP_WATER: 3000,
    GROUP_BUILDING: 4000,
}


def parse_type(spec: str) -> Tuple[int, int]:
    """Resolve a ``"<group>:<name>"`` mapping entry.

    Args:
        spec: e.g. ``"vegetation:short_grass"``.

    Returns:
        ``(group_code, type_code)``.

    Raises:
        ValueError: Unknown group, unknown name, or value out of PALM's range.
    """
    if not isinstance(spec, str) or ":" not in spec:
        raise ValueError(f"Type spec must be '<group>:<name>', got {spec!r}")
    group, name = (part.strip() for part in spec.split(":", 1))
    if group not in GROUP_OF_NAME:
        raise ValueError(
            f"Unknown surface group {group!r} in {spec!r}; "
            f"expected one of {sorted(GROUP_OF_NAME)}"
        )
    table = TABLES[group]
    if name not in table:
        raise ValueError(
            f"Unknown {group} type {name!r}; expected one of {sorted(table)}"
        )
    value = table[name]
    lo, hi = VALID_RANGE[f"{group}_type"]
    if not lo <= value <= hi:
        raise ValueError(
            f"{group}_type {name!r} = {value} is outside PALM's valid range [{lo}, {hi}]"
        )
    return GROUP_OF_NAME[group], value


def inspection_code(group: int, type_value: int) -> int:
    """Encode a (group, type) pair for the single-band inspection raster."""
    return INSPECTION_OFFSET[group] + int(type_value)
