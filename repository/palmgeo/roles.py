"""What each data role is for, and what breaks without it.

palmgeo addresses every data source by a **role** -- ``dtm``, ``ortho``,
``cadastre`` and so on -- and never by product name. That indirection is what
lets one body of code build a driver for Lower Saxony, Manhattan and Aswan; it
is also what makes a source list unreadable, because a role name says nothing
about what the data ends up as.

This module is the missing half: for each role, which parts of the finished
static driver rest on it, whether the pipeline can run without it, and what the
driver looks like when it is absent. The `sources` stage prints this beside
each source so that the person approving a download can see what they are
approving, and so that "no `topo` source for this region" reads as "the driver
will have no street types" rather than as a warning to be clicked past.

It was a comment in the header of `sources_ni.yml` before it was code. A
comment cannot be printed next to the data it describes, and it was already
out of date.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass(frozen=True)
class Role:
    """One kind of input data, described by what depends on it."""

    name: str
    #: One line, in the vocabulary of somebody who knows cities, not GIS.
    summary: str
    #: Driver variables that cannot exist at all without this role.
    produces: Tuple[str, ...]
    #: Driver variables this role improves but does not solely determine.
    contributes_to: Tuple[str, ...] = ()
    #: False means the pipeline refuses to start without it.
    optional: bool = False
    #: What the driver looks like when this role has no source. Written for the
    #: report, so it must be a consequence, not a restatement of the absence.
    without: str = ""
    #: Roles whose data this one is measured against; a mismatch here is a
    #: silent error rather than a loud one, so it is worth showing.
    pairs_with: Tuple[str, ...] = ()


#: Ordered as the pipeline consumes them, which is also roughly the order of
#: how much a driver suffers when one is missing.
ROLES: Tuple[Role, ...] = (
    Role(
        name="dtm",
        summary="Bare-earth terrain: the height of the ground with everything on it removed.",
        produces=("zt",),
        contributes_to=("vegetation_height", "buildings_2d"),
        optional=False,
        without="Nothing can be built. Every height in a driver is measured "
                "from the terrain, so PALM has no ground to stand the city on.",
    ),
    Role(
        name="dsm",
        summary="Surface model: the height of whatever the laser or camera hit first "
                "-- roof, treetop or ground.",
        produces=("vegetation_height",),
        contributes_to=("buildings_2d",),
        optional=False,
        pairs_with=("dtm",),
        without="No object heights. Buildings fall back to the 3D model or to a "
                "typological guess, and vegetation height is unavailable, which "
                "means no resolved trees and no canopy.",
    ),
    Role(
        name="ortho",
        summary="Aerial imagery with a near-infrared band, from which vegetation "
                "is told apart from everything else.",
        produces=(),
        contributes_to=("vegetation_type", "pavement_type", "soil_type", "trees"),
        optional=True,
        without="The land-use / land-cover fusion collapses to plain land use: a "
                "residential block becomes one uniform surface instead of lawn, "
                "driveway and hedge, and tree species classification is lost.",
    ),
    Role(
        name="cadastre",
        summary="The official land register: what each parcel is *for*, plus building "
                "footprints and bridges.",
        produces=("building_id", "building_type", "bridges_2d", "bridges_id"),
        contributes_to=("vegetation_type", "pavement_type", "water_type", "soil_type"),
        optional=False,
        without="Nothing can be built. Land use is the context every pixel is "
                "classified in, and building footprints come from here.",
    ),
    Role(
        name="buildings3d",
        summary="Modelled buildings carrying a height attribute -- LoD2 and equivalents.",
        produces=(),
        contributes_to=("buildings_2d",),
        optional=True,
        without="Building heights rest entirely on the surface model. That is "
                "usually the better measurement anyway; the 3D model is the "
                "fallback where the surface model has no coverage.",
    ),
    Role(
        name="topo",
        summary="Road centrelines with a class attribute.",
        produces=("street_type", "street_crossings"),
        optional=True,
        without="The driver carries no street types and no crossings. PALM runs; "
                "anything keyed to street class -- traffic emissions, street "
                "canyon analysis -- does not.",
    ),
    Role(
        name="water",
        summary="Water bodies, where the cadastre does not carry them itself.",
        produces=(),
        contributes_to=("water_type",),
        optional=True,
        without="Water comes from the cadastre alone. Where that is incomplete -- "
                "typically outside Europe -- rivers and lakes classify as land, "
                "which changes the surface energy balance of the whole domain.",
    ),
)

ROLE_BY_NAME: Dict[str, Role] = {r.name: r for r in ROLES}

#: Roles the pipeline refuses to start without.
REQUIRED: Tuple[str, ...] = tuple(r.name for r in ROLES if not r.optional)


def get(name: str) -> Optional[Role]:
    return ROLE_BY_NAME.get(name)


def feeds(name: str) -> List[str]:
    """Every driver variable a role produces or contributes to, for display."""
    role = ROLE_BY_NAME.get(name)
    if role is None:
        return []
    return list(role.produces) + list(role.contributes_to)


def describe(name: str) -> str:
    """A one-line 'what it becomes' for the report's rightmost column."""
    role = ROLE_BY_NAME.get(name)
    if role is None:
        return "unknown role"
    parts = list(role.produces)
    if role.contributes_to:
        parts.append("+" + ", ".join(role.contributes_to))
    return ", ".join(parts) if parts else "(no direct driver variable)"
