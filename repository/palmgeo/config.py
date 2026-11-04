"""Configuration loading and validation.

Three YAML files make up a run:

``sources_<region>.yml``
    Where the data comes from. Addressed by *role*, never by name.
``mapping_<region>.yml``
    What the region's attribute vocabularies mean in PALM terms.
``<run>.yml``
    Domains, thresholds, output options.

Only the first two are region specific; porting the workflow means replacing
them, not editing Python. Validation is deliberately eager -- a typo in a PALM
type name should fail in the first second, not after a 70 GB download.
"""

from __future__ import annotations

import copy
import hashlib
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import yaml

from palmgeo import palm_types as pt
from palmgeo.grid import (
    Domain, Tier, build_tiers, multigrid_leftover, nearest_smooth,
)

#: Config subtrees each resumable stage actually reads. A stage's ``.done``
#: marker carries a hash of these, so editing a threshold invalidates exactly
#: the work that depended on it -- without this, a re-run silently keeps the
#: old tiles and the change appears to have no effect.
STAGE_KEYS: Dict[str, Tuple[str, ...]] = {
    "tiles": (
        "classification.ndsm",
        "sources.dsm",
        "sources.ortho_epochs",
        "sources.prefer_uniform_epoch",
        "sources.max_epoch_spread_days",
        "output.keep_ortho_bands",
    ),
    "compose": (
        "classification",
        "buildings",
        "vegetation",
        "soil",
        "extra_inputs",
        "csd.settings.height_high_vegetation_lower_threshold",
    ),
}

#: Generation of the fingerprint scheme per stage, the ``n`` in ``v<n>:``. A
#: marker of an older generation is trusted once and re-stamped (see
#: marker_up_to_date). ``compose`` is at 3 since its vector inputs are hashed.
FINGERPRINT_GENERATION: Dict[str, int] = {"tiles": 2, "compose": 3}

#: `vegetation_mode` shorthands expanded into palm_csd domain switches.
VEGETATION_MODES: Dict[str, Dict[str, bool]] = {
    # Keep tall vegetation as a rough LSM surface type. Correct for coarse
    # parent domains, where a canopy cannot be resolved vertically anyway --
    # and only there: forest keeps z0 = 2 m, PALM needs z0 < 0.5 * z_mo with
    # z_mo = dz/2, so dz must exceed 8 m, and the log law behind it wants the
    # first grid level well above the canopy, so in practice dz >= 16 m.
    # Below that, use `mixed`. _validate warns about the hard limit.
    "unresolved": {
        "generate_vegetation_patches": False,
        "replace_high_vegetation_types": False,
        "generate_single_trees": False,
    },
    # Resolve canopy patches into a 3D LAD field, but no individual trees.
    "mixed": {
        "generate_vegetation_patches": True,
        "replace_high_vegetation_types": True,
        "generate_single_trees": False,
    },
    # Everything resolved, including individual trees.
    "resolved": {
        "generate_vegetation_patches": True,
        "replace_high_vegetation_types": True,
        "generate_single_trees": True,
    },
}

#: Where individual trees come from.
#:
#: ``rentz``      detect them from the canopy height model and the orthophoto
#:               (see :mod:`palmgeo.steps.trees`; named after Joana-Pauline
#:               Rentz, whose method it is); the default, because the open
#:               geodata already contain everything the method needs.
#: ``external``   take them from ``extra_inputs.trees.path``, e.g. a municipal
#:               tree register.
#: ``unresolved`` no individual trees at all; canopy is only ever a patch.
TREE_MODES = ("rentz", "external", "unresolved")


logger = logging.getLogger("palmgeo.config")


class ConfigError(Exception):
    """Raised for any invalid or inconsistent configuration."""


def _read_yaml(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        raise ConfigError(f"Configuration file not found: {path}")
    with path.open(encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ConfigError(f"{path} must contain a YAML mapping at top level")
    return data


def _get(d: Dict[str, Any], dotted: str, default: Any = None) -> Any:
    """Look up ``a.b.c`` in nested dicts, returning ``default`` if absent."""
    node: Any = d
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return node


@dataclass
class SurfaceRule:
    """How one land-use class becomes a PALM surface type."""

    key: str
    water: Tuple[int, int] | None = None
    high_ndvi: Tuple[int, int] | None = None
    low_ndvi: Tuple[int, int] | None = None
    high_veg: Tuple[int, int] | None = None
    high_veg_fixed: bool = False


@dataclass
class Mapping:
    """Parsed and validated ``mapping_<region>.yml``."""

    region: str
    fallback: SurfaceRule
    landuse: Dict[str, SurfaceRule]
    default_high_vegetation: Tuple[int, int]
    phenology: Dict[str, Tuple[int, int]]
    building_type: Dict[str, Any]
    bridges: Dict[str, Any]
    street_type: Dict[str, Any]
    soil: Dict[str, Any]
    #: What a detected tree becomes in palm_csd, per crown class
    #: (``conifer``/``broadleaf`` -> ``tree_type`` index and ``species`` label),
    #: or None if the region declares none. See :mod:`palmgeo.steps.trees`.
    single_trees: Optional[Dict[str, Dict[str, Any]]] = None
    #: Land-use keys sorted longest-first, so prefix lookup is greedy.
    _ordered_keys: List[str] = field(default_factory=list)

    def lookup(self, key: str | None) -> SurfaceRule:
        """Longest-prefix lookup of a land-use key.

        Prefix matching suits hierarchical codes -- ALKIS's ``tn`` in Lower
        Saxony, where ``43001`` covers every ``43001 ...`` subclass. With flat
        vocabularies such as OSM tags it matches like an exact lookup, except
        that a short key also catches longer ones it is a prefix of.
        """
        if key:
            key = key.strip()
            for candidate in self._ordered_keys:
                if key.startswith(candidate):
                    return self.landuse[candidate]
        return self.fallback


def _rule(key: str, raw: Dict[str, Any], inherit: SurfaceRule | None) -> SurfaceRule:
    def conv(name: str) -> Tuple[int, int] | None:
        if name not in raw or raw[name] is None:
            return None
        try:
            return pt.parse_type(str(raw[name]))
        except ValueError as exc:
            raise ConfigError(f"landuse['{key}'].{name}: {exc}") from exc

    rule = SurfaceRule(
        key=key,
        water=conv("water"),
        high_ndvi=conv("high_ndvi"),
        low_ndvi=conv("low_ndvi"),
        high_veg=conv("high_veg"),
        high_veg_fixed=bool(raw.get("high_veg_fixed", False)),
    )
    if rule.water is None and inherit is not None:
        rule.high_ndvi = rule.high_ndvi or inherit.high_ndvi
        rule.low_ndvi = rule.low_ndvi or inherit.low_ndvi
        rule.high_veg = rule.high_veg or inherit.high_veg
    if rule.water is None and (rule.high_ndvi is None or rule.low_ndvi is None):
        raise ConfigError(
            f"landuse['{key}'] must define either `water` or both "
            f"`high_ndvi` and `low_ndvi`"
        )
    return rule


def _single_trees(raw: Any, path: Path) -> Optional[Dict[str, Dict[str, Any]]]:
    """Validate ``single_trees``: one palm_csd tree type per crown class.

    Which species a detected crown stands for is a property of the region, not
    of the detection: the same broadleaf crown is a beech in Lower Saxony and
    something else entirely on the Nile. So the region names it, and there is
    no default to fall back on.
    """
    if raw is None:
        return None
    if not isinstance(raw, dict) or set(raw) != {"conifer", "broadleaf"}:
        raise ConfigError(
            f"{path}: `single_trees` needs exactly the two crown classes "
            f"`conifer` and `broadleaf`, each with `tree_type` and `species`"
        )
    out: Dict[str, Dict[str, Any]] = {}
    for cls, entry in raw.items():
        if not isinstance(entry, dict) or "tree_type" not in entry or "species" not in entry:
            raise ConfigError(f"{path}: single_trees.{cls} needs `tree_type` and `species`")
        try:
            tree_type = int(entry["tree_type"])
        except (TypeError, ValueError):
            raise ConfigError(
                f"{path}: single_trees.{cls}.tree_type must be a palm_csd tree "
                f"type index (palm_csd/data/tree_defaults.csv)"
            ) from None
        if tree_type < 0:
            raise ConfigError(f"{path}: single_trees.{cls}.tree_type must be >= 0")
        out[str(cls)] = {"tree_type": tree_type, "species": str(entry["species"])}
    return out


def load_mapping(path: Path) -> Mapping:
    """Read and validate a region mapping file."""
    raw = _read_yaml(path)
    fb_raw = raw.get("fallback")
    if not isinstance(fb_raw, dict):
        raise ConfigError(f"{path}: `fallback` section is required")
    fallback = _rule("fallback", fb_raw, None)

    landuse: Dict[str, SurfaceRule] = {}
    for key, value in (raw.get("landuse") or {}).items():
        if not isinstance(value, dict):
            raise ConfigError(f"{path}: landuse['{key}'] must be a mapping")
        landuse[str(key)] = _rule(str(key), value, fallback)

    try:
        default_high = pt.parse_type(
            str(raw.get("default_high_vegetation", "vegetation:deciduous_broadleaf_trees"))
        )
        phen_raw = raw.get("phenology") or {}
        phenology = {k: pt.parse_type(str(v)) for k, v in phen_raw.items()}
    except ValueError as exc:
        raise ConfigError(f"{path}: {exc}") from exc

    for name, (group, _) in [("default_high_vegetation", default_high), *(
        (f"phenology.{k}", v) for k, v in phenology.items()
    )]:
        if group != pt.GROUP_VEGETATION:
            raise ConfigError(f"{path}: {name} must be a `vegetation:` type")

    mapping = Mapping(
        region=str(raw.get("region", "?")),
        fallback=fallback,
        landuse=landuse,
        default_high_vegetation=default_high,
        phenology=phenology,
        building_type=raw.get("building_type") or {},
        bridges=raw.get("bridges") or {},
        street_type=raw.get("street_type") or {},
        soil=raw.get("soil") or {},
        single_trees=_single_trees(raw.get("single_trees"), path),
    )
    mapping._ordered_keys = sorted(landuse, key=len, reverse=True)

    # Validate the building/street/soil vocabularies eagerly, too.
    age = str(mapping.building_type.get("default_age", "1951_2000"))
    for base in set((mapping.building_type.get("by_gfk") or {}).values()) | {
        str(mapping.building_type.get("default", "residential"))
    }:
        name = f"{base}_{age}"
        if name not in pt.BUILDING_TYPE:
            raise ConfigError(
                f"{path}: building_type '{base}' + default_age '{age}' -> '{name}' "
                f"is not a PALM building type ({sorted(pt.BUILDING_TYPE)})"
            )
    for value in set((mapping.street_type.get("by_value") or {}).values()) | {
        str(mapping.street_type.get("default", "unclassified"))
    }:
        if value not in pt.STREET_TYPE:
            raise ConfigError(f"{path}: unknown street_type '{value}'")
    for value in set((mapping.soil.get("by_landuse_prefix") or {}).values()) | {
        str(mapping.soil.get("constant", "medium_fine"))
    }:
        if value not in pt.SOIL_TYPE:
            raise ConfigError(f"{path}: unknown soil_type '{value}'")
    return mapping


@dataclass
class Config:
    """Everything one pipeline invocation needs."""

    path: Path
    root: Path
    raw: Dict[str, Any]
    sources: Dict[str, Any]
    #: The licence registry, keyed by id. Loaded once here so that no stage
    #: has to know where it lives or what to do when it is missing.
    licences: Dict[str, Any]
    mapping: Mapping
    domains: List[Domain]
    tiers: List[Tier]
    runtime: Dict[str, Any]

    # ---------------------------------------------------------------- paths --
    @property
    def workdir(self) -> Path:
        """Intermediates. Defaults to ``<outdir>/work`` so that everything
        belonging to one run stays under one directory; set
        ``project.workdir`` to an absolute path to put the tens of gigabytes
        on another disk."""
        custom = _get(self.raw, "project.workdir")
        return (self.root / str(custom)) if custom else (self.outdir / "work")

    @property
    def outdir(self) -> Path:
        return self.root / str(_get(self.raw, "project.outdir", "out"))

    @property
    def cachedir(self) -> Path:
        """Download cache, shared between runs when ``project.cachedir`` is set.

        The condensed source indices alone are ~150 MB of GeoJSON, so pointing
        several projects at one cache avoids re-fetching them per project. The
        cache is content-addressed by file name and safe to share.
        """
        custom = _get(self.raw, "project.cachedir")
        return (self.root / str(custom)) if custom else (self.workdir / "cache")

    def tier_dir(self, tier: Tier) -> Path:
        return self.workdir / "tiers" / tier.name

    def domain_dir(self, domain: Domain) -> Path:
        """Where this domain's palm_csd input rasters are written.

        ``products/`` rather than the output root, because the driver files
        themselves have to sit at the top of the project directory: PALM reads
        the root domain under the bare name ``<jobname>_static`` and will not
        look into a sub-directory for it.
        """
        return self.outdir / "products" / domain.name

    @property
    def configdir(self) -> Path:
        """Everything that *describes* the run: run config, palm_csd config,
        the region YAMLs that were used, PALM namelists."""
        return self.outdir / "config"

    @property
    def docsdir(self) -> Path:
        """Everything a human looks at to judge the driver: the palm_csd
        sheets and the project README."""
        return self.outdir / "docs"

    # -------------------------------------------------------------- helpers --
    def get(self, dotted: str, default: Any = None) -> Any:
        return _get(self.raw, dotted, default)

    @property
    def epsg(self) -> int:
        # Always present once load() has run: it fills grid.epsg from the
        # region file when the run config leaves it out.
        value = _get(self.raw, "grid.epsg")
        if value is None:
            raise ConfigError("grid.epsg is not set")
        return int(value)

    @property
    def crs(self) -> str:
        return f"EPSG:{self.epsg}"

    def source(self, name: str) -> Dict[str, Any]:
        src = (self.sources.get("sources") or {}).get(name)
        if src is None:
            raise ConfigError(
                f"Source '{name}' is not defined in {_get(self.raw, 'project.sources')}"
            )
        return src

    #: Where a run configuration names the source for each role. ``topo`` is
    #: the odd one out because street handling has its own enable switch.
    _ROLE_SELECTOR = {
        "dtm": "sources.dtm",
        "dsm": "sources.dsm",
        "ortho": "sources.ortho",
        "buildings3d": "sources.buildings3d",
        "cadastre": "sources.cadastre",
        "topo": "sources.street_type.source",
        "water": "sources.water",
    }

    def source_name_for_role(self, role: str) -> Optional[str]:
        """Which source fills ``role``, or None if the region offers none.

        Three levels, most specific first:

        1.  the run configuration, e.g. ``sources.dsm: <source name>``;
        2.  the region file's ``defaults:`` block;
        3.  the only source declaring that role, if there is exactly one.

        Level 3 exists so that a minimal region file needs no ``defaults:``
        block at all. Where several sources compete for a role and nothing
        chooses between them, this returns None rather than picking one:
        silently taking whichever came first in the file is how a driver ends
        up built from the wrong surface model.
        """
        key = self._ROLE_SELECTOR.get(role)
        chosen = _get(self.raw, key) if key else None
        if chosen:
            return str(chosen)
        chosen = (self.sources.get("defaults") or {}).get(role)
        if chosen:
            return str(chosen)
        candidates = [
            name for name, src in (self.sources.get("sources") or {}).items()
            if src.get("role") == role
        ]
        return candidates[0] if len(candidates) == 1 else None

    def sources_for_role(self, role: str) -> List[str]:
        """Every source in the region file that declares ``role``.

        The alternatives, in other words -- what the ``sources`` stage lists
        beneath the chosen one so that the choice is visible rather than
        implicit.
        """
        return [
            name for name, src in (self.sources.get("sources") or {}).items()
            if src.get("role") == role
        ]

    def source_for_role(self, role: str) -> Tuple[str, Dict[str, Any]]:
        """Return the (name, definition) of the source configured for ``role``."""
        name = self.source_name_for_role(role)
        if name is None:
            candidates = self.sources_for_role(role)
            where = str(_get(self.raw, "project.sources"))
            if not candidates:
                raise ConfigError(
                    f"No source in {where} declares role '{role}'. Either add "
                    f"one, or switch the feature off in the run configuration."
                )
            raise ConfigError(
                f"Role '{role}' is filled by {len(candidates)} sources in "
                f"{where} ({', '.join(sorted(candidates))}) and nothing chooses "
                f"between them. Set `sources.{role}:` in the run configuration, "
                f"or add a `defaults:` entry to {Path(where).name}."
            )
        src = self.source(name)
        if src.get("role") not in (role, None):
            raise ConfigError(
                f"Source '{name}' has role '{src.get('role')}' but was selected "
                f"for '{role}'"
            )
        return name, src

    def licence_ids_for_source(self, name: str) -> List[str]:
        """The licence ids a source declares, as a list even when there is one.

        A source may name several: the New York cadastre is city open data on
        one side of the Hudson and OpenStreetMap on the other, and both
        licences travel with the driver. A source that names none is reported
        as ``unknown``, which the registry deliberately models as the most
        restrictive position rather than as no position.
        """
        value = self.source(name).get("licence")
        if value is None:
            return ["unknown"]
        return [str(value)] if isinstance(value, str) else [str(v) for v in value]

    def licence_holder(self, name: str, licence_id: str) -> Optional[str]:
        """Who holds the rights in ``name`` under ``licence_id``.

        ``licence_holder`` is a plain string for the ordinary case of one
        source, one rights holder. Where a source is fused from several -- the
        New York cadastre is the city's on one side of the Hudson and
        OpenStreetMap's on the other -- it is a mapping keyed by licence id,
        because otherwise one holder gets credited under a licence that is not
        theirs, which is worse than crediting nobody.
        """
        value = self.source(name).get("licence_holder")
        if isinstance(value, dict):
            return value.get(licence_id)
        return str(value) if value else None

    def tier_of(self, domain: Domain) -> Tier:
        for tier in self.tiers:
            if domain in tier.domains:
                return tier
        raise ConfigError(f"Domain {domain.name} belongs to no tier")

    # --------------------------------------------------------- staleness --
    def stage_fingerprint(
        self, stage: str, inputs: Optional[Dict[str, Optional[Path]]] = None
    ) -> str:
        """Hash of everything ``stage`` reads, used to invalidate ``.done``.

        A stage's ``.done`` marker carries this value. When a setting the stage
        depends on changes, the hash changes, the marker no longer matches and
        exactly the affected work is redone -- without it, a re-run silently
        keeps the old tiles and the edit appears to have had no effect, which
        is the worst possible failure mode because it looks like success.

        What goes in is therefore not "the configuration" but "everything this
        stage reads", from all three files:

        *the run config*, through ``STAGE_KEYS``;

        *the region mapping*, hashed whole for ``compose``, because a changed
        PALM type assignment is as much a reason to recompose as a changed
        threshold;

        *the region sources*, hashed whole for ``tiles``, plus the source name
        actually resolved for each raster role. Both halves are needed and for
        different reasons. The file hash catches an edited ``index_url``,
        ``resolution`` or ``nodata`` -- changes that alter the pixels while
        leaving the run config untouched. The resolved names catch a change of
        source that the run config does not mention either, because since the
        ``defaults:`` block moved into the region file a run may name no source
        at all and still switch surface models when that file is edited.

        *the input files* the caller names in ``inputs``, hashed by content.
        ``compose`` reads the vector layers of earlier stages, and re-running
        ``vectors`` or ``buildings`` rewrites them -- with new building IDs at
        the least -- while every setting stays as it was. Without their hashes
        the old composition survived such a re-run: in the example, a root built
        from the new vectors sat on top of nests composed from the old ones. The
        caller hashes once per stage run, not per tile; the layers can be
        hundreds of megabytes.

        The ``v<n>:`` prefix marks the generation of this scheme, per stage
        (``FINGERPRINT_GENERATION``). See :func:`marker_up_to_date` for why that
        matters.
        """
        payload: Dict[str, Any] = {k: _get(self.raw, k) for k in STAGE_KEYS[stage]}
        if inputs:
            payload["_inputs"] = {
                name: _file_digest(Path(path))
                for name, path in sorted(inputs.items())
                if path is not None and Path(path).is_file()
            }
        if stage == "compose":
            mapping_path = self.root / str(_get(self.raw, "project.mapping", ""))
            if mapping_path.is_file():
                payload["_mapping"] = hashlib.sha1(mapping_path.read_bytes()).hexdigest()
        if stage == "tiles":
            sources_path = self.root / str(_get(self.raw, "project.sources", ""))
            if sources_path.is_file():
                payload["_sources"] = hashlib.sha1(sources_path.read_bytes()).hexdigest()
            payload["_resolved"] = {
                role: self.source_name_for_role(role)
                for role in ("dtm", "dsm", "ortho")
            }
        blob = json.dumps(payload, sort_keys=True, default=str).encode()
        generation = FINGERPRINT_GENERATION.get(stage, 2)
        return f"v{generation}:{hashlib.sha1(blob).hexdigest()[:16]}"


def marker_up_to_date(path: Path, fingerprint: str) -> bool:
    """Whether a stage's ``.done`` marker still matches the current settings.

    Markers from an older palmgeo are **trusted once and re-stamped** rather
    than invalidated. Two generations qualify: the earliest markers carried no
    fingerprint at all, and the next carried a bare hex digest computed before
    the region sources file was part of the hash.

    The alternative would be to invalidate them, which is the theoretically
    clean answer and the practically wrong one: it would discard tens of
    gigabytes of correctly derived tiles per project -- around 160 GB across
    the projects that exist -- to re-derive work that is almost certainly
    identical, because the region files did not change when the hashing scheme
    did. What the strictness is *for* is catching a change made from now on,
    and re-stamping preserves that: every marker written after this point
    carries the current generation and is compared in full.

    The generation is encoded as a ``v<n>:`` prefix rather than inferred from
    the hash, so that the next change to the scheme is a one-line decision
    instead of an archaeology exercise. The same rule covers every later
    generation: a marker from an older one is accepted once and re-stamped.
    When ``compose`` moved to generation 3, which hashes its vector inputs, that
    kept the compositions of every existing project -- consistent unless their
    vectors were rebuilt after composing, which is what generation 3 catches
    from now on. ``palmgeo clean <cfg> surface`` recomposes one that is not.
    """
    if not path.exists():
        return False
    stored = path.read_text().strip()
    if not stored or not stored.startswith("v") or _generation(stored) < _generation(fingerprint):
        # An older generation (0: empty, 1: bare digest, or a lower v<n>):
        # accept once, re-stamp, and check properly from here on.
        path.write_text(fingerprint)
        return True
    return stored == fingerprint


def _generation(fingerprint: str) -> int:
    """The ``n`` of a ``v<n>:...`` fingerprint, 0 if there is none."""
    head, sep, _ = fingerprint.partition(":")
    if not sep or not head.startswith("v") or not head[1:].isdigit():
        return 0
    return int(head[1:])


def _file_digest(path: Path) -> str:
    """SHA-1 of a file's content, read in 1 MB chunks."""
    digest = hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_project_root(start: Path) -> Path:
    """The directory every path in a run configuration is resolved against.

    A run configuration does not have to live in ``repository/config``. Configs
    that describe one particular site belong next to that site's data, in
    ``static_driver/<run>/config/``, so that publishing the repository does not
    publish them. Both locations have to resolve
    ``repository/config/sources_<region>.yml`` and ``static_driver/<run>`` to the
    same place, which they do as long as every path is read relative to the project
    root rather than to the config file.

    The root is the nearest ancestor holding ``repository/config``. Failing
    that -- a configuration kept somewhere else entirely -- the caller must say
    where the project is via ``project.root``.
    """
    for candidate in [start, *start.parents]:
        if (candidate / "repository" / "config").is_dir():
            return candidate
    raise ConfigError(
        f"Cannot locate the palmgeo project root above {start}: no ancestor "
        f"contains repository/config. Either move the configuration inside the "
        f"project tree or set `project.root` to the project directory."
    )


#: Placeholders a region file may use where it has to name a local path.
#:
#: A region file describes a REGION, and is meant to serve every run in it.
#: A staged source breaks that: the index it points at is produced per run and
#: lives in that run's working directory, so writing the path out ties
#: `sources_eg.yml` to one Egyptian project and makes the next one edit a file
#: that is supposed to be shared. Worse, an absolute path also ties it to one
#: machine and carries the author's account name into a published file.
#:
#: Expanded once at load time, against the run being loaded:
#:
#:     index_url: "{workdir}/staging/indices/dtm_eg.json"
#:
#: so Aswan and Luxor read the same region file and each finds its own index.
SOURCE_PLACEHOLDERS = ("workdir", "outdir", "project_root")


def _expand_placeholders(cfg: "Config") -> None:
    """Resolve ``{workdir}`` and friends in every string in the sources file.

    Applied to the whole ``sources`` subtree rather than to the two fields that
    need it today, because a future field naming a staged artefact would
    otherwise be silently left unexpanded -- and a literal ``{workdir}`` in a
    path fails in a way that reads like a corrupt configuration rather than a
    missing feature.
    """
    values = {
        "workdir": str(cfg.workdir),
        "outdir": str(cfg.outdir),
        "project_root": str(cfg.root),
    }

    def walk(node: Any) -> Any:
        if isinstance(node, str):
            if "{" not in node:
                return node
            out = node
            for key, value in values.items():
                out = out.replace("{" + key + "}", value)
            return out
        if isinstance(node, dict):
            return {k: walk(v) for k, v in node.items()}
        if isinstance(node, list):
            return [walk(v) for v in node]
        return node

    cfg.sources = walk(cfg.sources)


def _warn_absolute_local_paths(sources: Dict[str, Any], sources_path: Path) -> None:
    """Warn when a versioned region file hard-codes an absolute local path.

    This is the failure mode that survives longest without being noticed,
    because it works perfectly on the machine it was written on. It only
    surfaces when somebody else clones the repository, when the run moves to a
    cluster, or when the file is published with a home directory in it.

    A warning rather than an error: an absolute path is legitimate for data
    genuinely kept outside the project, and palmgeo cannot tell that case from
    a forgotten one. Naming the field is enough to make it a decision.
    """
    def scan(node: Any, trail: str) -> None:
        if isinstance(node, str):
            bare = node[len("file://"):] if node.startswith("file://") else node
            if bare.startswith("/") and "://" not in bare:
                logger.warning(
                    "%s: %s is an absolute local path (%s). It will not resolve "
                    "on another machine, and it carries this account's home "
                    "directory into a versioned file. Use a placeholder instead: "
                    "{workdir}, {outdir} or {project_root}.",
                    sources_path.name, trail, bare,
                )
        elif isinstance(node, dict):
            for key, value in node.items():
                scan(value, f"{trail}.{key}" if trail else str(key))
        elif isinstance(node, list):
            for i, value in enumerate(node):
                scan(value, f"{trail}[{i}]")

    scan(sources.get("sources") or {}, "sources")


def _validate_licence_ids(
    sources: Dict[str, Any],
    licences: Dict[str, Any],
    sources_path: Path,
    licences_path: Path,
) -> None:
    """Fail on a licence id the registry does not define.

    A typo here would otherwise resolve to nothing and be reported as "no
    licence stated", which reads like a gap in the data rather than a gap in
    this file -- and would send somebody looking in the wrong place.
    """
    bad: List[str] = []
    for name, src in (sources.get("sources") or {}).items():
        value = src.get("licence")
        if value is None:
            continue
        ids = [value] if isinstance(value, str) else list(value)
        for lid in ids:
            if str(lid) not in licences:
                bad.append(f"{name}: '{lid}'")
    if bad:
        raise ConfigError(
            f"{sources_path.name} names licence ids that {licences_path.name} "
            f"does not define: {', '.join(bad)}. Add them to the registry, or "
            f"correct the spelling -- do not leave the field empty, which the "
            f"registry models as 'unknown' and therefore as most restrictive."
        )


def _resolve_epsg(raw: Dict[str, Any], sources: Dict[str, Any], sources_path: Path) -> int:
    """The run's CRS: ``grid.epsg``, else the region file's ``crs:``.

    There is deliberately nothing after those two. A default CRS is a default
    *place* -- this used to fall back to UTM 32N, and a run anywhere else that
    forgot the key had every index reprojected into Germany and found no
    coverage, or was told its region file disagreed with a CRS it never wrote.
    The resolved value is written back into the run config, so everything
    downstream reads one number.
    """
    from pyproj import CRS

    declared = _get(raw, "grid.epsg")
    src_crs = str(sources.get("crs") or "").strip()
    if declared is None:
        if not src_crs:
            raise ConfigError(
                f"grid.epsg is not set and {sources_path.name} declares no `crs:`, "
                f"so the run's coordinate reference system is unknown. Set "
                f"grid.epsg to the projected CRS the region's products are in."
            )
        code = CRS.from_user_input(src_crs).to_epsg()
        if code is None:
            raise ConfigError(
                f"{sources_path.name}: `crs: {src_crs}` has no EPSG code; set grid.epsg"
            )
        epsg = int(code)
        raw.setdefault("grid", {})["epsg"] = epsg
    else:
        epsg = int(declared)
    if src_crs and src_crs.upper() != f"EPSG:{epsg}":
        raise ConfigError(
            f"grid.epsg is EPSG:{epsg} but {sources_path.name} declares its products "
            f"in {src_crs}. Reprojecting every tile would be lossy and slow; either "
            f"align grid.epsg with the source CRS or add an explicit reprojection step."
        )

    # Every distance in palmgeo -- pixel sizes, buffers, the 1 km work tiles,
    # tree spacing -- is in metres of the run's CRS. A geographic or a
    # foot-based CRS would run through without complaint and be wrong
    # everywhere at once.
    crs = CRS.from_epsg(epsg)
    units = {(axis.unit_name or "").lower() for axis in crs.axis_info}
    if not crs.is_projected or not units <= {"metre", "meter"}:
        raise ConfigError(
            f"EPSG:{epsg} ({crs.name}) is not a projected CRS in metres "
            f"(axis units: {', '.join(sorted(units)) or 'none'}). palmgeo, PALM "
            f"and palm_csd all work in metres; use a metric projected CRS, e.g. "
            f"the UTM zone of the site, and reproject feet-based products while "
            f"staging them (PORTING.md, section 5b)."
        )
    return epsg


def _validate_index_crs(sources: Dict[str, Any], sources_path: Path) -> None:
    """Every index must say what CRS its geometries are in.

    There is no safe default. GeoJSON is WGS84 by specification, a staged
    index is usually in the run's CRS, and the Lower Saxony indices are in
    UTM 32N -- guessing any one of them puts the others' tiles somewhere else
    and reports the area as uncovered.
    """
    from pyproj import CRS
    from pyproj.exceptions import CRSError

    missing, bad = [], []
    for name, src in (sources.get("sources") or {}).items():
        if not isinstance(src, dict) or not src.get("index_url"):
            continue
        value = src.get("index_crs")
        if not value:
            missing.append(str(name))
            continue
        try:
            CRS.from_user_input(str(value))
        except CRSError:
            bad.append(f"{name}: {value}")
    if missing:
        raise ConfigError(
            f"{sources_path.name}: source(s) {', '.join(missing)} have an "
            f"`index_url` but no `index_crs`. Name the CRS the index geometries "
            f"are in (EPSG:4326 for plain GeoJSON)."
        )
    if bad:
        raise ConfigError(f"{sources_path.name}: unknown index_crs: {'; '.join(bad)}")


def load(path: str | Path) -> Config:
    """Load a run configuration together with its sources and mapping files."""
    cfg_path = Path(path).resolve()
    raw = _read_yaml(cfg_path)
    declared = _get(raw, "project.root")
    root = (
        (cfg_path.parent / str(declared)).resolve() if declared
        else find_project_root(cfg_path.parent)
    )

    for key in ("project.sources", "project.mapping"):
        if not _get(raw, key):
            raise ConfigError(
                f"{cfg_path.name}: `{key}` is not set. It names the region file "
                f"to use, e.g. repository/config/{key.split('.')[1]}_<region>.yml; "
                f"there is no default region."
            )
    sources_path = root / str(_get(raw, "project.sources"))
    mapping_path = root / str(_get(raw, "project.mapping"))
    sources = _read_yaml(sources_path)
    mapping = load_mapping(mapping_path)
    src_region, map_region = sources.get("region"), mapping.region
    if src_region and map_region != "?" and str(src_region) != map_region:
        logger.warning(
            "%s is for region %s but %s is for region %s. Every land-use key the "
            "mapping does not know falls through to its fallback rule; check "
            "project.sources and project.mapping.",
            sources_path.name, src_region, mapping_path.name, map_region,
        )

    # The registry lives beside the region files and is not per-run, so it is
    # found relative to the project root rather than named in the run config.
    licences_path = root / "repository" / "config" / "licences.yml"
    if not licences_path.is_file():
        # Treated as an error rather than as "no licence information", because
        # the two are indistinguishable downstream and one of them is a lie:
        # the `sources` stage would report every source as `unknown` and every
        # driver as maximally restricted, which reads like a finding about the
        # data instead of a missing file in the installation.
        raise ConfigError(
            f"The licence registry {licences_path} is missing. It is part of "
            f"the repository and every source refers to it; without it palmgeo "
            f"cannot say what a driver built from these sources may be used for."
        )
    licences = _read_yaml(licences_path).get("licences") or {}
    _validate_licence_ids(sources, licences, sources_path, licences_path)

    epsg = _resolve_epsg(raw, sources, sources_path)
    _validate_index_crs(sources, sources_path)

    domains = _parse_domains(raw)
    profile_name = str(raw.get("profile", "local"))
    profiles = _get(raw, "runtime.profiles", {}) or {}
    if profile_name not in profiles:
        raise ConfigError(
            f"profile '{profile_name}' is not defined under runtime.profiles "
            f"({sorted(profiles)})"
        )
    runtime = copy.deepcopy(_get(raw, "runtime", {}) or {})
    runtime.update(copy.deepcopy(profiles[profile_name]))
    runtime["profile"] = profile_name

    wr = _get(raw, "grid.work_resolution", {}) or {}
    tiers = build_tiers(
        domains,
        buffer=float(_get(raw, "grid.buffer", 256.0)),
        detail_factor=float(wr.get("detail_factor", 8)),
        minimum=float(wr.get("min", 1.0)),
        maximum=float(wr.get("max", 8.0)),
        snap=wr.get("snap") or [1.0, 2.0, 4.0, 8.0],
        fixed=wr.get("fixed"),
    )

    cfg = Config(
        path=cfg_path,
        root=root,
        raw=raw,
        sources=sources,
        licences=licences,
        mapping=mapping,
        domains=domains,
        tiers=tiers,
        runtime=runtime,
    )
    # Placeholders first, then the audit: expanding turns "{workdir}/..." into
    # an absolute path, so scanning afterwards would flag every expansion as a
    # hard-coded one. The audit reads the file as written.
    _warn_absolute_local_paths(sources, sources_path)
    _expand_placeholders(cfg)
    _validate(cfg)
    return cfg


def _parse_domains(raw: Dict[str, Any]) -> List[Domain]:
    section = raw.get("domains")
    if not isinstance(section, dict) or not section:
        raise ConfigError("`domains` section is required and must not be empty")
    known = {
        "pixel_size", "nx", "ny", "dz", "origin_x", "origin_y",
        "origin_lon", "origin_lat", "domain_parent", "vegetation_mode",
    }
    interpolate_terrain_default = _get(raw, "csd.interpolate_terrain", False)
    if not isinstance(interpolate_terrain_default, bool):
        raise ConfigError(
            f"csd.interpolate_terrain must be true or false, "
            f"got {interpolate_terrain_default!r}"
        )
    domains: List[Domain] = []
    for name, spec in section.items():
        if not isinstance(spec, dict):
            raise ConfigError(f"domain '{name}' must be a mapping")
        for req in ("pixel_size", "nx", "ny", "dz"):
            if req not in spec:
                raise ConfigError(f"domain '{name}' is missing required key '{req}'")
        has_xy = spec.get("origin_x") is not None and spec.get("origin_y") is not None
        has_ll = spec.get("origin_lon") is not None and spec.get("origin_lat") is not None
        if not (has_xy or has_ll):
            raise ConfigError(
                f"domain '{name}' needs either origin_x/origin_y or origin_lon/origin_lat"
            )
        if not has_xy:
            raise ConfigError(
                f"domain '{name}' only has origin_lon/origin_lat. Convert them to "
                f"origin_x/origin_y in the target CRS first -- the tiling maths "
                f"needs projected coordinates. `palmgeo lonlat` does the conversion."
            )
        # everything not consumed by Domain is passed through to palm_csd
        options = {k: v for k, v in spec.items() if k not in known}
        mode = str(spec.get("vegetation_mode", "mixed"))
        if mode not in VEGETATION_MODES:
            raise ConfigError(
                f"domain '{name}': vegetation_mode '{mode}' unknown "
                f"({sorted(VEGETATION_MODES)})"
            )
        # explicit switches in the domain block win over the shorthand
        expanded = dict(VEGETATION_MODES[mode])
        expanded.update({k: v for k, v in options.items() if k in expanded})
        options.update(expanded)
        options["vegetation_mode"] = mode

        # Terrain blending at the nest boundary. Off by default: a nest is built
        # to resolve terrain its parent cannot, and blending replaces exactly
        # that over a 50-cell margin. Global, overridable per domain. The root
        # has no parent to blend with, so it is always false -- asking for it
        # there is a configuration error, not something to silence.
        is_root = spec.get("domain_parent") is None
        if is_root:
            if options.get("interpolate_terrain") is True:
                raise ConfigError(
                    f"domain '{name}' is the root and has no parent, so "
                    f"interpolate_terrain cannot be true: there is no coarser "
                    f"terrain to blend into. Remove the key or set it to false."
                )
            options["interpolate_terrain"] = False
        elif "interpolate_terrain" not in options:
            options["interpolate_terrain"] = interpolate_terrain_default
        domains.append(
            Domain(
                name=str(name),
                pixel_size=float(spec["pixel_size"]),
                nx=int(spec["nx"]),
                ny=int(spec["ny"]),
                dz=float(spec["dz"]),
                origin_x=float(spec["origin_x"]),
                origin_y=float(spec["origin_y"]),
                domain_parent=spec.get("domain_parent"),
                options=options,
            )
        )
    return domains


def _validate(cfg: Config) -> None:
    names = {d.name for d in cfg.domains}
    for d in cfg.domains:
        if d.domain_parent is not None and d.domain_parent not in names:
            raise ConfigError(
                f"domain '{d.name}': domain_parent '{d.domain_parent}' is not a "
                f"configured domain ({sorted(names)})"
            )
        if d.nx < 1 or d.ny < 1:
            raise ConfigError(f"domain '{d.name}': nx and ny must be >= 1")
        if d.pixel_size <= 0 or d.dz <= 0:
            raise ConfigError(f"domain '{d.name}': pixel_size and dz must be > 0")

    roots = [d for d in cfg.domains if d.domain_parent is None]
    if len(roots) != 1:
        raise ConfigError(
            f"exactly one domain must have no domain_parent, found {len(roots)}: "
            f"{[d.name for d in roots]}"
        )

    # ---------------------------------------------------------------- nesting --
    # PALM's nesting has three geometric preconditions, and it enforces them
    # itself -- by refusing to start, after the driver has been built and
    # submitted. Checking them here turns a failure that costs a queue slot and
    # a day into one that costs a second, which is the whole argument for
    # eager validation. All three are checked; none of them can be repaired
    # automatically, because every repair moves a boundary and therefore
    # changes what is simulated. The messages say what the nearest legal value
    # is and leave the choice where it belongs.
    by_name = {d.name: d for d in cfg.domains}
    for d in cfg.domains:
        if d.domain_parent is None:
            continue
        p = by_name[d.domain_parent]

        # 1. Containment. A nest that pokes out of its parent has no boundary
        #    condition on the protruding side; there is nothing to interpolate
        #    from.
        cb, pb = d.bounds, p.bounds
        if not (cb[0] >= pb[0] and cb[1] >= pb[1] and cb[2] <= pb[2] and cb[3] <= pb[3]):
            raise ConfigError(
                f"domain '{d.name}' {cb} is not contained in its parent "
                f"'{p.name}' {pb}"
            )

        # 2. Integer resolution ratio. The parent hands the child its boundary
        #    values, and a fractional ratio means a child cell straddles two
        #    parent cells with no consistent rule for which value it inherits.
        ratio = p.pixel_size / d.pixel_size
        if abs(ratio - round(ratio)) > 1e-9 or round(ratio) < 1:
            # Only ratios that yield a pixel size somebody would actually
            # write down. 32/3 = 10.6667 m is arithmetically valid and
            # nobody's grid; offering it as a suggestion is noise.
            options = [
                v for v in (p.pixel_size / n for n in (2, 3, 4, 5, 6, 8, 16))
                if v > 0 and abs(v * 100 - round(v * 100)) < 1e-9
            ]
            raise ConfigError(
                f"domain '{d.name}': pixel_size {d.pixel_size:g} m does not divide "
                f"its parent '{p.name}' pixel_size {p.pixel_size:g} m evenly "
                f"(ratio {ratio:.4f}). PALM needs a whole number of child cells "
                f"per parent cell. Nearby valid values: "
                f"{', '.join(f'{o:g}' for o in options)} m (4 is the usual choice)."
            )

        # 3. Origin phase. The child's corner must fall on a parent CELL
        #    boundary, not merely inside the parent. Off by half a parent cell
        #    and the nest border cuts through the cells it is supposed to read
        #    its boundary condition from.
        for axis, child_o, parent_o in (
            ("x", d.origin_x, p.origin_x), ("y", d.origin_y, p.origin_y)
        ):
            offset_cells = (child_o - parent_o) / p.pixel_size
            if abs(offset_cells - round(offset_cells)) > 1e-6:
                snapped = parent_o + round(offset_cells) * p.pixel_size
                raise ConfigError(
                    f"domain '{d.name}': origin_{axis} {child_o:g} is "
                    f"{offset_cells:.4f} parent cells from '{p.name}' "
                    f"origin_{axis} {parent_o:g}, so the nest border does not "
                    f"fall on a parent cell boundary. The nearest position that "
                    f"does is {snapped:g} ({snapped - child_o:+.3f} m). Moving it "
                    f"changes what is simulated, so palmgeo will not do it for you."
                )

    # ------------------------------------------------------- solver-friendly --
    # A warning, not an error. An unfavourable point count produces a correct
    # driver that PALM will happily read; what it costs is solver speed, which
    # is the user's trade to make. See grid.multigrid_leftover.
    for d in cfg.domains:
        bad = [
            (axis, count, multigrid_leftover(count))
            for axis, count in (("nx", d.width), ("ny", d.height))
            if multigrid_leftover(count) != 1
        ]
        if bad:
            detail = "; ".join(
                f"{axis}+1 = {count} leaves the factor {left}, so the multigrid "
                f"solver stops coarsening there (try "
                f"{', '.join(str(c) for c in nearest_smooth(count))})"
                for axis, count, left in bad
            )
            logger.warning("domain '%s': %s", d.name, detail)

    # ------------------------------------------------------ forest roughness --
    # A warning, not an error: whether the domain contains tall vegetation at
    # all is only known after `compose`, and PALM offers a way round the
    # abort. But that way round is easy to miss, and without it the driver is
    # rejected at model start -- after the job has waited in the queue.
    z0_limit_dz = 4.0 * pt.HIGH_VEGETATION_Z0
    for d in cfg.domains:
        if d.options.get("replace_high_vegetation_types"):
            continue
        if d.dz <= z0_limit_dz:
            logger.warning(
                "domain '%s': replace_high_vegetation_types is false (vegetation_mode "
                "'%s'), so tall vegetation stays an LSM surface type with z0 up to "
                "%g m. At dz = %g m the first grid level is %g m, and PALM stops with "
                "LSM0048 because z0 >= 0.5 * z_mo -- unless allow_roughness_limitation "
                "= .T., which clips z0 to %g m. Use vegetation_mode 'mixed' to resolve "
                "the canopy instead; 'unresolved' needs dz > %g m.",
                d.name, d.options.get("vegetation_mode"), pt.HIGH_VEGETATION_Z0,
                d.dz, d.dz / 2, d.dz / 4, z0_limit_dz,
            )

    # `resolved` promises single trees; they have to come from somewhere.
    mode = str(_get(cfg.raw, "vegetation.trees.mode", "rentz"))
    if mode not in TREE_MODES:
        raise ConfigError(
            f"vegetation.trees.mode '{mode}' unknown ({sorted(TREE_MODES)})"
        )
    wants_trees = [d.name for d in cfg.domains if d.options.get("generate_single_trees")]
    if wants_trees:
        if mode == "unresolved":
            raise ConfigError(
                f"domain(s) {wants_trees} use vegetation_mode 'resolved' "
                f"(generate_single_trees), but vegetation.trees.mode is 'unresolved'. "
                f"Set it to 'rentz' to detect trees from the open geodata, to "
                f"'external' plus extra_inputs.trees.path to supply your own, or use "
                f"vegetation_mode 'mixed', which resolves canopy patches without "
                f"individual trees."
            )
        if mode == "external" and not _get(cfg.raw, "extra_inputs.trees.path"):
            raise ConfigError(
                f"domain(s) {wants_trees} need single trees and "
                f"vegetation.trees.mode is 'external', but extra_inputs.trees.path "
                f"is not set."
            )
        if mode == "rentz" and cfg.mapping.single_trees is None:
            raise ConfigError(
                f"domain(s) {wants_trees} need single trees and "
                f"vegetation.trees.mode is 'rentz', but the region mapping "
                f"declares no `single_trees`: which palm_csd tree type a detected "
                f"conifer or broadleaf crown becomes is a regional choice."
            )
        if mode == "rentz":
            finest = min(t.resolution for t in cfg.tiers)
            if finest > 2.0:
                raise ConfigError(
                    f"domain(s) {wants_trees} need single trees and "
                    f"vegetation.trees.mode is 'rentz', but the finest work tier is "
                    f"{finest:g} m. Crown segmentation needs 2 m or finer. Lower "
                    f"grid.work_resolution.max (or set `fixed`), or use "
                    f"vegetation_mode 'mixed'."
                )

    if _get(cfg.raw, "classification.phenology.enabled", False):
        if not _get(cfg.raw, "sources.ortho_epochs.leaf_off"):
            raise ConfigError(
                "classification.phenology.enabled requires a "
                "sources.ortho_epochs.leaf_off definition"
            )
