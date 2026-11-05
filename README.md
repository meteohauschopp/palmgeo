# palmgeo — a PALM static driver from open geodata

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23973676.svg)](https://doi.org/10.5281/zenodo.23973676)

palmgeo turns publicly available open geodata — national, regional or
municipal — into a complete, verified
[PALM](https://palm.muk.uni-hannover.de/) static driver: the input rasters and
the YAML configuration for PALM's own `palm_csd`, for an arbitrary number of
nested domains with individually configurable extent, resolution, origin and
vegetation treatment.

**This repository is a workflow, not a finished product.** It is written so that
a competent agent — human or AI — can be handed a location and a domain layout
and produce a driver that PALM will actually run. The Python code is
region-agnostic; everything that knows about a particular country or state lives
in two declarative YAML files. The repository ships these regions:

| region | files | status |
|---|---|---|
| **Lower Saxony**, Germany | `sources_ni.yml`, `mapping_ni.yml` | the reference: works out of the box from LGLN OpenGeoData (`dl-de/by-2-0`); the committed example is built from it |
| **New York City** (+ New Jersey) | `sources_us_ny.yml`, `mapping_us_ny.yml` | case study: needs the staging scripts in `repository/ports/us_ny/` |
| **Egypt** | `sources_eg.yml`, `mapping_eg.yml` | case study on global data only: needs `repository/ports/eg/`; non-commercial by default (FABDEM) |
| **anywhere** | `sources_global.yml`, `mapping_global.yml` | not runnable by itself: a library of worldwide sources and OSM rules to build a new region from |

A new region means writing those two files — [`PORTING.md`](PORTING.md) says how.

> ### If you are an AI agent asked to build a driver with this repository
> Read **Appendix A** before you touch anything. It lists what you must have
> before starting, what to ask the user for when it was not in the prompt, what
> you may compute yourself, and what you must never decide alone.

---


## Contents

**[Part I — What palmgeo does, and why it is built this way](#part-i--what-palmgeo-does-and-why-it-is-built-this-way)**

| | |
|---|---|
| [1. The problem](#1-the-problem-open-geodata-to-a-static-driver) | open geodata to a static driver |
| [2. Repository layout](#2-repository-layout) | where everything lives |
| [3. Fifteen minutes to a driver](#3-fifteen-minutes-to-a-driver) | install, smoke test, first run |

**[Part II — The road through the pipeline](#part-ii--the-road-through-the-pipeline)**

| | |
|---|---|
| [4. Sources, licences, and what each becomes](#4-sources-licences-and-what-each-becomes) | what you may use, and for what |
| [5. Why rasters are written on the domain grid](#5-why-rasters-are-written-on-the-domain-grid) | tiers, streaming, transfer volume |
| [6. Classification: land use × land cover](#6-classification-land-use--land-cover) | how a pixel becomes a PALM type |
| [7. Vegetation](#7-vegetation-resolved-or-unresolved-per-nest) | patches, single trees, LAD |
| [8. Buildings](#8-buildings) | measured heights, courtyards, bridges |

**[Part III — What PALM will not forgive](#part-iii--what-palm-will-not-forgive)**

| | |
|---|---|
| [9. Grid, nesting, and terrain](#9-grid-nesting-and-terrain) | nx, 2-smoothness, cell centring, nesting |
| [10. Verification and worked examples](#10-verification-and-worked-examples) | proof it works |

**[Part IV — Operating it](#part-iv--operating-it)**

| | |
|---|---|
| [11. Running it](#11-running-it) | stages, resumability, cleaning |
| [12. Storage and memory](#12-storage-and-memory-behaviour) | running big on a small disk |
| [13. Assumptions](#13-assumptions) | what is decided, not measured |
| [14. Known limitations](#14-known-limitations) | what it gets wrong |
| [15. Troubleshooting](#15-troubleshooting) | the errors you will actually hit |
| [16. Versioning and attribution](#16-versioning-and-attribution) | branches, releases, credit |

**[Appendices](#appendices)**

| | |
|---|---|
| [A. Driving palmgeo with an agent](#appendix-a--driving-palmgeo-with-an-agent) | the intake checklist |
| [B. Case study: Lower Saxony](#appendix-b--case-study-lower-saxony) | one region's data, in detail |
| [C. Porting to another region](#appendix-c--porting-to-another-region) | new state, new country |


---

# Part I — What palmgeo does, and why it is built this way


## 1. The problem: open geodata to a static driver

To simulate the air over a city, PALM has to be told what the city *is*: a
number for every cell of a three-dimensional grid. How high the ground lies.
Where the buildings stand and how tall they are. Which square metre is asphalt
and which is lawn, where the trees are, how much leaf each one carries. That
description is one NetCDF file — the **static driver** — and PALM will not
start without it.

Almost all of it already exists as open government data: laser-scanned terrain
and surface models, aerial photography with a near-infrared band, the cadastre,
3D building models, road networks. What does not exist is the path from that
data to that file.

That path is longer than it looks. It means working out which of a dozen
products covers your patch of ground and which epoch of each to use; moving
tens of gigabytes without needing tens of gigabytes of disk; agreeing every
dataset onto one coordinate system and one grid; deciding what each cadastre
code means in PALM's vocabulary; measuring building heights rather than
trusting an attribute; finding individual trees in a canopy height model; and
getting a grid geometry right that PALM rejects outright if it is off by one.

**palmgeo is that path, written down and made repeatable.** It ends where
PALM's own `palm_csd` begins: palmgeo produces the GeoTIFFs and the
configuration that `palm_csd` turns into the driver. It never writes a NetCDF
driver itself and never runs PALM.

### What this document is for

The chapters below answer two questions, and the second one is the reason the
document is long. *What does each stage do* is a paragraph. *Why does it do it
that way* is the part that took the measurements — why building heights come
from the surface model rather than the 3D model, why NDVI alone loses half of
every tree crown, why a nest boundary shifts by metres when one switch is
turned on. Where a decision rests on a number, the number is here.

Three principles run through all of it, and they are worth stating once:

**No Python file names a data product.** Sources are addressed by *role* —
`dtm`, `ortho`, `cadastre` — and every product name, attribute name and code
value lives in two YAML files per region. Porting means writing those two
files. See [`PORTING.md`](PORTING.md).

**Nothing is discarded silently.** Where the pipeline drops a footprint, fills
a terrain gap, or falls back from a measurement to an assumption, it counts it
and says so. A run that quietly loses 3 % of its buildings is worse than one
that fails.

**Every stage is resumable.** A city-scale run moves tens of gigabytes over
hours. Each unit of work records a fingerprint of the settings it was built
with, so an interruption costs the tile it was on, and changing one threshold
rebuilds exactly the tiles that depended on it.

*If you are meeting palmgeo for the first time, [`HUMAN-README.md`](HUMAN-README.md)
is the shorter way in. If you are an AI agent asked to build a driver, read
**Appendix A** before touching anything.*


---


## 2. Repository layout

Two directories, and the difference between them is the whole idea: `repository/`
is the pipeline — generic, versioned, published — and `static_driver/` is what
you did with it — specific, private, enormous.

```
palmgeo/
├── README.md                 this file — the whole workflow, in detail
├── HUMAN-README.md           the short version, for someone meeting palmgeo
├── PORTING.md                how to adapt to another state or country
├── LICENSE
├── requirements.txt
│
├── repository/               ← everything that is versioned
│   ├── palmgeo/              the region-agnostic pipeline (importable package)
│   │   ├── cli.py            the stages, in order
│   │   ├── roles.py          what each data role becomes, and what breaks without it
│   │   ├── config.py grid.py sources.py rasters.py cache.py palm_types.py
│   │   └── steps/            one module per stage
│   ├── config/               only what is generic enough to publish
│   │   ├── sources_<rr>.yml  WHERE the data comes from     ─┐ replace these two
│   │   ├── mapping_<rr>.yml  WHAT its vocabulary means      ─┘ to port to a region
│   │   ├── licences.yml      the licence registry every source refers to
│   │   ├── template.yml      annotated starting point for a new run
│   │   └── example.yml       the worked example, and the end-to-end test
│   ├── tools/                things you run against any project
│   │   ├── selftest.py       has palmgeo stopped producing what it claims?
│   │   ├── finalize.py       bring a project directory into final shape
│   │   ├── verify_driver.py  structural acceptance check against the run config
│   │   ├── check_classification.py  plausibility: land use vs. realised surface
│   │   ├── check_coverage.py check for systematically empty parts of a driver
│   │   ├── check_ortho.py    sanity-check a staged mosaic before classifying it
│   │   ├── find_domain.py    largest PALM-legal domain inside the data coverage
│   │   ├── ft_to_m.py        elevation VALUES ft -> m (gdalwarp does not do this)
│   │   ├── fill_water.py     hydro-flatten nodata inside the water mask
│   │   ├── fix_colorinterp.py  relabel a NIR band mistagged as Alpha
│   │   ├── run_guarded.sh    memory-capped wrapper for a long job
│   │   └── update_p3d.py     bring existing p3d namelists in line with the grid
│   ├── ports/                one directory per region that needs staging
│   │   ├── staging_index.py  the contract: writing an index palmgeo can read
│   │   ├── eg/               Egypt, from the global fallback sources
│   │   └── us_ny/            New York City, plus a New Jersey supplement
│   └── .cache/               shared download cache (git-ignored)
│
└── static_driver/            ← one directory per project (git-ignored,
    │                             except `example/`, which is committed in full)
    └── <run>/
        ├── <run>_static          root domain — PALM reads exactly this name
        ├── <run>_static_N02 …    the nests — the suffix is what marks a nest
        ├── config/
        │   ├── <run>.yml             the run configuration
        │   ├── <run>_layout.yml      optional: PALM domain names and PE budget
        │   ├── <run>_static_csd.yml  the palm_csd configuration, generated
        │   ├── sources_<rr>.yml      a RECORD of the region files that were
        │   ├── mapping_<rr>.yml        used, copied here by finalize.py
        │   └── <run>_p3d[r][_N0x]    PALM namelists, if update_p3d.py was used
        ├── products/<domain>/*.tif   what palm_csd read, kept for inspection
        ├── docs/
        │   ├── <run>_static*.pdf     palm_csd plots, one sheet per domain
        │   ├── DATA-SOURCES.md       generated: what it is built from, and the
        │   │                           licence that governs using it
        │   └── README.md             what this run is, and on what data it rests
        └── work/                     intermediates (git-ignored, 50–100 GB)
            ├── state.json                the plan: tiers, tiles, ortho epochs
            ├── vector/*.gpkg             cadastre, roads, bridges, trees
            ├── tiers/<tier>/…            per-tile rasters, composed classification
            └── inspect/                  QGIS control layers with a class legend
```

### Why a project directory is shaped like that

The four sub-directories answer four different questions, and mixing them is
what makes a two-year-old driver uninterpretable:

| | |
|---|---|
| `config/` | **what was asked for.** Enough to rebuild the driver, and enough to know what its numbers mean — a driver whose land-use mapping is unknown cannot be read back. |
| `products/` | **what palm_csd was given.** The GeoTIFFs between palmgeo and palm_csd; where you look when a driver is wrong and you need to know which of the two produced the error. |
| `docs/` | **what a human looks at.** The palm_csd sheets, and a README saying what was built, what was assumed and what was dropped. The NetCDF records none of that. |
| `work/` | **what it cost to get there.** Gigabytes, reproducible, deletable. Nothing here is a deliverable. |

The drivers themselves stay at the top, uncontained, because PALM will not look
into a sub-directory for them: it reads the root domain under the bare job name
`<run>_static`.

### Where a run configuration lives, and why

**In its own project, at `static_driver/<run>/config/<run>.yml`** — not in
`repository/config/`. A run config names a site, an author, an institution and
a study; it is the private part of the work. `static_driver/` is git-ignored, so
publishing the repository publishes the pipeline and not the projects run
through it. Only the region files, the template and the smoke test are versioned.

The `sources_*.yml` and `mapping_*.yml` in a project's `config/` are **copies,
not inputs**. `finalize.py` puts them there so the project can be read on its
own; each carries a header saying so. A re-run reads the originals under
`repository/config/`.

### How paths resolve

**Every path in a run config is relative to the project root** — the directory
holding `repository/` — regardless of where the config file itself sits or what
the current working directory is. palmgeo finds that root by walking up from the
config file until it sees `repository/config`. So the same lines work from
anywhere inside the tree:

```yaml
project:
  sources:  repository/config/sources_ni.yml
  mapping:  repository/config/mapping_ni.yml
  outdir:   static_driver/<run>
  workdir:  static_driver/<run>/work    # the default; override to move it
  cachedir: repository/.cache
```

Absolute paths work too — point `workdir` at an external disk for a city-scale
run. For a config kept entirely outside the project tree, set `project.root` to
the project directory; palmgeo says so explicitly rather than guessing.


---


## 3. Fifteen minutes to a driver

```bash
git clone https://github.com/meteohauschopp/palmgeo.git
cd palmgeo                                   # the project root; stay here
pip install -r requirements.txt              # the same stack palm_csd needs
export PYTHONPATH=repository            # or prefix every command with it

CFG=repository/config/example.yml

python -m palmgeo.cli sources $CFG           # what data exists, under what licence
python -m palmgeo.cli sources $CFG --approve # accept it; unlocks the downloads
python -m palmgeo.cli plan    $CFG           # dry run: tiles, epochs, volume
python -m palmgeo.cli all     $CFG           # about two minutes

MPLBACKEND=Agg env -u PROJ_DATA -u PROJ_LIB \
    palm_csd --pdf static_driver/example/config/example_static_csd.yml
python repository/tools/finalize.py      $CFG
python repository/tools/verify_driver.py $CFG   # must end in ALL CHECKS PASSED
```

The `--approve` step is not a formality: `vectors` and `rasters` refuse to start
until the data sources and their licences have been accepted (§4).

`example.yml` covers 9 km² of Hannover city centre — a 2 km root at 32 m, a
1 km nest at 8 m and a 256 m nest at 2 m with individual trees — and takes
about two minutes against
the live services. It is the fastest way to prove that the network, the geo
stack and `palm_csd` all work before committing to a real domain.

**The result of that run is already in the repository**, at
[`static_driver/example/`](static_driver/example/): all three drivers, the palm_csd
inputs, the plots and a README describing what came out. So you can look at a
finished driver before running anything — and afterwards compare yours against
it.

If `palm_csd` is not on the `PATH`, run it out of the PALM installation:

```bash
CSD=/path/to/palm_model_system-v25.10/packages/static_driver/palm_csd
env -u PROJ_DATA -u PROJ_LIB MPLBACKEND=Agg PYTHONPATH=$CSD \
    python $CSD/palm_csd.py --pdf <run>_static_csd.yml
```

**palmgeo does not run `palm_csd`.** It writes its configuration and its input
rasters. That separation is deliberate: `palm_csd` is PALM's own tool, versioned
with PALM, and the driver it writes is the thing PALM supports.

### Requirements

Python ≥ 3.10, plus the stack in `requirements.txt` (geopandas, rasterio,
shapely, pyproj, numpy, pandas, scipy, PyYAML, and scikit-image for single-tree
segmentation). `palm_csd` needs the same libraries, so installing palmgeo into
the environment that already runs `palm_csd` is the least painful route.

Network access to the source services is required for the `rasters` and
`vectors` stages. Nothing else phones home.

---

# Part II — The road through the pipeline

Five chapters, in the order the pipeline runs them. Each says what the stage
produces and, more importantly, why it produces it that way — the alternatives
that were tried, and what they cost.


## 4. Sources, licences, and what each becomes

### What the pipeline produces

For every configured domain, a directory of GeoTIFFs on **exactly** that
domain's grid. The source chain is given by *role*, not by product name —
`dtm`, `ortho`, `cadastre` and the rest are what the code addresses, and which
product fills each one is a region's business (Appendix B works one region
through in full):

| File | PALM variable | Type | Built from |
|---|---|---|---|
| `zt.tif` | terrain height | float32 | `dtm` |
| `buildings_2d.tif` | building height | float32 | nDSM per pixel inside the `cadastre` footprint; `buildings3d` where the surface model has no cover; a typological prior only as a last resort, counted separately |
| `building_id.tif` | building id | int32 | `cadastre` footprint id, renumbered to int32 |
| `building_type.tif` | building type | uint8 | `cadastre` building function |
| `vegetation_type.tif` | vegetation type | uint8 | `cadastre` land use × `ortho` NDVI × nDSM × phenology |
| `pavement_type.tif` | pavement type | uint8 | `cadastre` land use × `ortho` NDVI |
| `water_type.tif` | water type | uint8 | `cadastre` land use, plus a `water` source where the cadastre carries none |
| `soil_type.tif` | soil type | uint8 | constant, refined by land use |
| `vegetation_height.tif` | vegetation height | float32 | nDSM = `dsm` − `dtm` |
| `bridges_2d.tif`, `bridges_id.tif` | bridges | float32 / int32 | `cadastre` structures + nDSM |
| `street_type.tif`, `street_crossings.tif` | street type / crossings | uint8 | `topo` |
| `trees.gpkg` | single trees | vector | `dsm`/`dtm` canopy height + `ortho`, per crown |

plus the `palm_csd` configuration ready to run, a provenance `manifest.json`
per domain, and a QGIS inspection layer at work resolution carrying the full
PALM classification of every pixel.

`palm_csd` then turns these into the NetCDF driver, adding `buildings_3d`,
`lad`, `bad`, `tree_id`, `tree_type`, `surface_fraction`, `albedo_pars` and the
rest of the PIDS static-driver vocabulary.

`palmgeo/roles.py` is the authoritative version of this mapping: it records,
per role, which of these variables depend on it and what the driver looks like
when the role has no source at all. The `sources` stage prints it.

### Which data, under which licence

A static driver is a derivative work of every dataset that went into it, and it
inherits every one of their licences at once. That is easy to lose sight of: by
the time the driver exists it is one NetCDF file with no visible seam between
the terrain that was public domain and the footprints that were share-alike.

`palmgeo sources` is the stage that keeps it visible. It runs before anything
is downloaded in bulk and answers three questions per data role.

**What is actually available here?** A region file lists what a region
*publishes*. Whether any of it covers your particular few square kilometres is
a different question, and the stage measures it: each source's index is probed
against the area of interest on a 160×160 grid. On the Manhattan run this is
what surfaces the fact that the city's elevation products cover 93.8 % of the
domain — the rest is New Jersey, where every New York source stops. Bounding
boxes overstate coverage where a tile is partly nodata, so the figure is
reported as an upper bound.

**Under what licence?** Every source names licences by id from
`repository/config/licences.yml`, a registry that carries, per licence, the
authoritative URL, an attribution template, any known grey area, and four
obligation flags:

| flag | what it decides |
|---|---|
| `attribution` | a credit must travel with the driver and with anything published from it |
| `share_alike` | a distributed derivative must carry the same licence |
| `non_commercial` | use primarily for commercial advantage is not permitted |
| `redistribution_restricted` | the data may not be passed on, only results derived from it |

The combination rule is the only safe one: **an obligation applies to the
driver if any single input imposes it.** There is no averaging. One CC BY-NC-SA
input among six makes the whole driver non-commercial and share-alike, which is
what the Aswan run reports — the FABDEM terrain alone determines it, and
nothing about the other five sources softens that.

A source that states no licence resolves to `unknown`, which the registry
deliberately models as the *most* restrictive position: attribution required,
redistribution restricted. An unexamined source is not a free source, and
modelling the gap as "no obligations" would be the one error that cannot be
caught later.

**What part of the driver rests on it?** `palmgeo/roles.py` records, per role,
which driver variables it produces and what the driver looks like without it.
That turns "no `topo` source for this region" from a warning into a statement:
this driver will have no street types and no crossings.

### The approval

Reading the report and agreeing to it are different acts, so running the stage
approves nothing. `--approve` records a fingerprint over the chosen source per
role and the licences those sources carry; `vectors` and `rasters` refuse to
start without it. `plan` is exempt on purpose — being made to approve licences
before being allowed to find out how large the download is would be the wrong
order.

The fingerprint is narrow by design. Changing a domain origin, a resolution or
an NDVI threshold does not change what may be done with the result, so it must
not invalidate the approval and make the gate cry wolf. Changing a source, or a
licence changing under a source, does.

The stage also writes `static_driver/<run>/docs/DATA-SOURCES.md`, which is what
a later reader of the driver consults to find out whether they may use it.


---


## 5. Why rasters are written on the domain grid

`palm_csd/geo_converter.py::read_raster_to_dst` chooses between two paths:

* `cut_aligned_raster_to_dst` — a **windowed** read, taken only when the source
  CRS, the pixel size *and* the grid origin (within 0.2 px) all match the target
  domain;
* `reproject_raster_to_dst` — everything else, which begins with
  `src_values = src.read(masked=True)`, i.e. **the entire raster into RAM**.

Handing palm_csd one shared 1 m master raster covering all domains would take the
second path for every domain and every variable. For a 48 × 48 km area that is
2.3 · 10⁹ cells — 9 GB per float variable. So each domain gets its own set of
rasters, at its own pixel size, anchored on its own origin. The remaining memory
cost is inherent to palm_csd:

```
RAM ≈ (nx+1) · (ny+1) · 4 bytes per 2D field
      20 km @ 2 m  = 10⁸ cells ≈ 400 MB per field
```

`palmgeo plan` prints, per domain, whether aggregation is *exact* (the domain
cell is a whole block of work pixels) or *approximate*, and why.

### Resolution tiers

Rasterising a 48 km root domain at 1 m only to average it down to 32 m would move
~1000× more orthophoto data than 32 m can resolve. Domains are therefore binned
into **work tiers**:

```
work_res(domain) = snap(clamp(pixel_size / detail_factor, min, max))
```

With the default `detail_factor: 8`, a 32 m root is processed at 4 m and a 2 m
nest at 1 m. Each PALM cell still sees 8 × 8 = 64 samples for its majority
statistic. The six-domain reference run then splits into

```
r1m   200 tile(s),  200 km²   domains N04, N05, N06
r4m  2500 tile(s), 2500 km²   domains root, N02, N03
```

Set `grid.work_resolution.fixed: 1.0` to force a single 1 m tier.

A tier's tiles are taken **per domain and unioned**, not from the bounding box of
the tier. Domains sharing a resolution need not be adjacent: the reference set
has two 8 m nests 19 km apart, whose bounding box is two and a half times their
combined area, and every tile in between would be downloaded and processed for
nothing. Per-domain planning brings that tier from ~280 tiles down to 200.

### What a run actually transfers

The decisive question for a large root is whether a source is *streamed* or
*downloaded*. `runtime.access: auto` streams whenever the work resolution is at
least `STREAM_RATIO` (4) times coarser than the source, because GDAL then serves
the read from a COG overview. Measured in Lower Saxony: a DOP20 orthophoto tile
is 565 MB on disk; read at 4 m it costs **3.3 MB**. Per 1 km work tile, one
epoch:

| source | native | @ 4 m | @ 1 m |
|---|---|---|---|
| DGM1 | 1 m | 1.04 MB (stream) | 5.65 MB (download) |
| DOM1 | 1 m | 1.10 MB (stream) | 6.46 MB (download) |
| DOP20 | 0.2 m | 1.16 MB (stream) | 39.5 MB (stream) |

`plan` estimates this from the COG block arithmetic rather than from a flat
per-tile figure — over a 2 500 km² root the two answers differ by a factor of
seven, which is the difference between "start it now" and "don't". For the
reference run it modelled 26 GB against 29.8 GB measured.


---


## 6. Classification: land use × land cover

The examples in this section are Lower Saxony's: the cadastre there is ALKIS,
and its land-use key `tn` is what `mapping_ni.yml` maps. Another region maps its
own vocabulary the same way — MapPLUTO codes in `mapping_us_ny.yml`, OSM tags in
`mapping_eg.yml` and `mapping_global.yml`.

ALKIS `nutzung` is a land **use** model. A single `41001 Wohnbaufläche` polygon
covers roofs, driveways, lawns and hedges alike, so mapping it 1:1 to a PALM
surface type would be wrong at 1 m. Each land-use class instead declares what it
becomes as a function of the observed land **cover**:

```yaml
"41001":                         # Wohnbaufläche
  high_ndvi: vegetation:short_grass
  low_ndvi:  pavement:paving_stones
"43002 VEG 1100":                # Wald, Laubholz
  high_ndvi: vegetation:short_grass
  high_veg:  vegetation:deciduous_broadleaf_trees
  high_veg_fixed: true           # never overridden by phenology
  low_ndvi:  vegetation:bare_soil
"44001": {water: water:river}
```

Applied per pixel, in order:

1. **Land use** → rule, by longest-prefix match on the land-use key (ALKIS `tn`
   in Lower Saxony), so `43002 VEG 1100` beats `43002`, and unknown keys fall
   through to `fallback`.
2. **Land cover** selects `high_ndvi` or `low_ndvi` — see below. This is what
   splits one residential polygon into lawn and driveway.
3. **nDSM** ≥ `height_high_vegetation_lower_threshold` promotes vegetated pixels
   to the class's `high_veg` type.
4. **Phenology** from the leaf-off epoch refines tall vegetation into
   evergreen / deciduous / mixed, except where `high_veg_fixed` says the cadastre
   already knows better.
5. **Hard overrides**, last and in this order: user polygons → water → buildings
   → bridges.

### Deciding what counts as vegetated

A single NDVI threshold is not enough, and the failure is systematic rather than
random. NDVI here comes from uncalibrated 8-bit digital numbers; in shadow both
the red and the near-infrared band collapse towards the noise floor and their
ratio stops meaning anything. Two things are lost as a result: the shaded halves
of tree crowns, and street trees standing over paved surfaces. Both matter for an
urban canopy study.

Three criteria are therefore combined, in an OR:

```
ndvi >= ndvi_threshold                                   # primary
greenness  and  ndvi >= ndvi_threshold_relaxed           # shaded or dull foliage
greenness  and  ndsm >= height_high_vegetation_threshold # a tall green object is a tree
```

where `greenness` is the hue criterion of Rentz's vegetation mask — green
normalised by total brightness, and green above red and blue by a margin.
Normalising by brightness is what makes it survive shadow: the *hue* of a shaded
leaf is still green even where its NDVI is not usable.

Measured over the Hannover test area, the combination adds 2.4 % of the domain,
and **79 % of what it adds stands above 2 m** — crown edges and street trees, not
invented lawn. Greenness alone also fires on 20 % of the water surface (algae,
reflections) and on 2 % of roofs; neither survives, because water and buildings
are burnt in from the cadastre afterwards and win by precedence.

All thresholds live under `classification.vegetation`, and `palmgeo calibrate`
prints the NDVI histogram of the actual imagery to set the primary one. The
criterion needs the visible bands, which the `rasters` stage keeps per tile when
`output.keep_ortho_bands` is on (the default).

> **Why not the BKG land cover model?** LBM-DE carries a sealing and a vegetation
> percentage per object, which sounds like the right prior. Its minimum mapping
> unit is 1 hectare, though, and everything this pipeline was missing is
> sub-hectare by nature — a street tree over a square, a lawn strip in a
> residential block, half a crown in shadow. For the coarse context it *can*
> resolve, ALKIS is finer and more current. It is wired into `sources_ni.yml` as
> `lbm_de` but disabled, and becomes useful when porting to a region without a
> detailed cadastral land-use layer.

Sealing is thereby a *per-pixel binary* decision rather than a fractional degree —
which is what PALM wants, since it assigns one discrete type per cell and derives
`surface_fraction` itself. Whether a non-vegetated pixel is *sealed* or *bare
soil* comes from the land-use context: unvegetated inside `43001 Ackerland` is
bare soil, unvegetated inside `41001 Wohnbaufläche` is pavement. That is the
whole point of fusing the two models.

### Aggregating types onto a coarse grid

Taking the mode of `vegetation_type` and of `pavement_type` independently would
let a cell that is 40 % pavement and 60 % vegetation carry both, and palm_csd's
precedence (`building > water > pavement > vegetation`, enforced in
`create_driver.py::process_surface_types`) would then silently award it to
pavement. Instead a `surface_group` raster decides the winner by majority first,
and only the winning group's subtype survives. Buildings are part of that same
majority vote, so a cell that is 45 % building becomes street or garden, not a
half-height block.


---


## 7. Vegetation: resolved or unresolved, per nest

`vegetation_mode` in a domain block expands into the palm_csd switches:

| mode | `generate_vegetation_patches` | `replace_high_vegetation_types` | `generate_single_trees` |
|---|---|---|---|
| `unresolved` | false | **false** | false |
| `mixed` | true | true | false |
| `resolved` | true | true | true |

`unresolved` deliberately keeps `replace_high_vegetation_types: false` so that
forest stays a rough LSM surface type instead of being flattened to short grass —
the physically right choice for a coarse parent domain, where a canopy cannot be
resolved vertically anyway. **Only there, though.** The tall-vegetation classes
keep PALM's forest roughness length, z0 = 2 m, and PALM's land-surface model
requires z0 < 0.5 · z_mo, where z_mo = dz/2 is the first grid level. At
`dz` ≤ 8 m it stops at model start with `LSM0048`. The namelist switch
`allow_roughness_limitation = .T.` lets it run by clipping z0 to dz/4, but at
8 m that clips nothing, and the log law behind the surface fluxes still wants
the first grid level well above the canopy. Rule of thumb: `unresolved` from
16 m upwards, `mixed` below. palmgeo warns at config load when a domain
breaks the hard limit.

Any individual switch written directly in the domain block overrides the
shorthand.

### Where individual trees come from

`vegetation.trees.mode` picks the source:

| mode | source |
|---|---|
| `rentz` (default) | detected from the canopy height model and the orthophoto by the `trees` stage |
| `external` | `extra_inputs.trees.path`, e.g. a municipal tree register |
| `unresolved` | none; canopy is only ever a patch |

The **`rentz`** method (`repository/palmgeo/steps/trees.py`) is named after
Joana-Pauline Rentz, whose Bachelor's thesis it is ported from (see
[Credits](#credits)). It needs nothing beyond what the pipeline already
downloads:

1. a vegetation mask from NDVI **and** the hue criterion, restricted to
   nDSM ≥ 2 m;
2. the canopy height model is smoothed and its local maxima are taken as tree
   tops, at least `min_distance_m` apart;
3. crowns are grown by a watershed whose elevation combines the inverted canopy
   height with the Sobel edges of the visible image, so touching crowns separate
   along the shading between them;
4. each crown is classified broadleaf or coniferous by a weighted vote over
   brightness, green−red, mean NDVI, crown shape and crown area;
5. allometry turns the measured height into trunk diameter and leaf area index,
   species-specific.

The output columns are palm_csd's own (`tree_height`, `tree_crown_diameter`,
`tree_trunk_diameter`, `tree_lai`, `tree_type`), and the types are the indices of
palm_csd's 87-species inventory, chosen per region (`single_trees` in the
mapping; in Lower Saxony 1 `Abies` for conifers, 33 `Fagus` for broadleaf). Two
departures from the original script, both forced by the
surrounding pipeline: tiles are read with a halo and a tree is kept only where
its crown centroid falls inside the tile proper, so a crown across a tile edge is
segmented whole and counted once; and detections on building footprints are
dropped, because a roof that reads as green is a green roof.

Segmentation is run on the **finest** work tier only and refuses to start above
2 m, where a pixel is wider than a crown. Whatever the source, the file is
reprojected, clipped per domain and copied into that domain's input directory, so
a domain directory stays self-contained.

Over the 4 km² around the Neues Rathaus the method finds 5231 trees (5117
broadleaf, 114 coniferous), median height 16.7 m, median crown 10.7 m; 419 of them
fall inside the 512 m domain. Their trunks are what gives the driver a non-zero
basal area density field.


---


## 8. Buildings

Footprint, identity and function come from the cadastre; **height comes from the
surface model**. In Lower Saxony's terms (other regions name their own columns
in `sources_<region>.yml`):

```
ALKIS gebaeude ──► footprint (2D, clean), objid → building_id, gfk → building_type
       │
       ├─ per pixel: nDSM inside the footprint          ──► buildings_2d
       ├─ per footprint: percentile(nDSM), eroded       ──► robust height
       └─ where the nDSM has no value: LoD2 measHeight  ──► fallback
```

### Why not the 3D building model

LoD2 is the obvious candidate and it is the wrong default, for two measured
reasons. Over the 512 m Hannover test domain:

* **One height per building throws away a lot.** Half the buildings span more
  than 5 m between the 5th and 95th percentile of their own footprint's nDSM —
  median 5.0 m, p90 20.9 m, max 59.2 m. That is roof ridges, setbacks, wings of
  different age, and towers, all flattened into one number. The driver now
  carries 28 256 distinct heights over its 42 866 building pixels, instead of one
  step per building.
* **The cadastre maps buildings that are no longer there.** 4.1 % of footprint
  area ends up carved out because the surface model measures it as open —
  courtyards, light wells, passages, demolitions. Five of 82 buildings are more
  than a quarter open.

For the Neues Rathaus specifically, LoD2 offers either 88.9 m (its tallest part,
which would extrude the whole 7513 m² footprint to the dome) or 35.7 m
area-weighted (which flattens the landmark). The measured surface reads p50 =
26.1 m with the dome reaching 97.6 m over 8.7 % of the footprint. Both LoD2
answers are wrong; the measurement is not.

Where a region has a 3D building model (role `buildings3d`; `lod2` in the
option names below, after the CityGML level of detail), it can still be fetched
as the fallback (`height_fallback: lod2`) for footprints the surface model does
not cover, as the example does, and `height_source: lod2` restores the old
behaviour. Without the key, there is no fallback. `height_aggregate` then decides how a building's several model parts
are reduced — `area_weighted` (default), `max` or `median`; over Hannover 6063 of
14 870 buildings consist of more than one part.

### Keeping the measurement honest

Three corrections are applied to the per-pixel height, each guarding against a
specific failure:

`ignore_vegetation_on_roofs`
: Vegetation inside a footprint is a crown or a green roof, so the nDSM there
  measures leaves, not the building. Those pixels are dropped and **refilled from
  the nearest measured pixel of the same building**.

  Substituting the building's robust height instead is only right when the
  vegetation sits on the roof. In a courtyard, a light well or a setback the
  vegetation stands well below the roof, and the substitution lifts the pixel to
  full building height — over the Welfengarten tier that inflated 0.4 % of all
  footprint area by more than 5 m, in the worst case turning a 7 m courtyard
  floor inside a 15 m building into 23 m of solid masonry. The nearest structural
  neighbour is the right answer in both cases: on a roof it *is* the roof, in a
  courtyard it is the courtyard floor. A single per-building number cannot serve
  both.

  **The test uses NDVI alone, not the combined vegetation mask** — the hue
  criterion cannot tell chlorophyll from green paint or verdigris. 8 % of the
  Neues Rathaus dome passes the hue test, because patinated copper is green in
  the visible; its NDVI is −0.39, because copper has no near-infrared signature.
  Using the combined mask here flattens the dome by 5 m.

`carve_open_areas`
: A footprint pixel whose nDSM is below `min_height` is not a building. This is
  what opens courtyards and passages the cadastre still maps as built.

  Note the reference: `min_height` is measured **above terrain**, so this opens
  courtyards whose floor is the ground, not those roofed at storey height. An
  inner court whose floor the surface model puts 7 m up is building to 7 m and
  open above — which is what it physically is. Both the ALS and the
  photogrammetric surface model have to agree on that floor for it to be trusted;
  at the Welfengarten they agree to 20 cm across nine years.

`per_pixel_tolerance`
: An optional cap at the building's own 99th percentile times a factor, against
  vegetation the NDVI test missed. **Off by default**, because a genuine tower is
  exactly what lives above its building's p99: at 1.15 the Rathaus dome loses
  12 m.

Set `per_pixel_height: false` for one robust height per building — appropriate at
coarse grid spacings, where per-pixel relief would be averaged away anyway.

Where several work tiers exist, footprints are sampled from the finest tier that
covers them and only fall through to the 3D model where no tier does. Over the
Welfenschloss area that is the difference between 18 465 and 4358 buildings
measured rather than modelled.

Bridges are extracted from ALKIS `bauwerke_f` (`bwf = 1800`) and given a deck
height from the nDSM. Because aggregation onto a coarse grid can put a building
and a bridge in the same cell — which palm_csd rejects — bridges are removed
wherever a building survives the majority vote.

**The deck height is taken per pixel**, for the same reason as roofs, and for one
more. PALM stacks a bridge on the *local* terrain of each cell
(`topography_mod.f90`: building type 7 is not lifted onto the highest terrain
under it, unlike every other building type), so the deck top it builds is
`zt + bridges_2d`, cell by cell. One height per deck, as before, got that wrong
twice. Arched decks and their ramps vary by metres along the span, and the
median put 19–31 % of deck pixels more than 1 m too high and 12–23 % more than
1 m too low, over five projects. And where the terrain model was not cut below
the deck — 1–9 % of deck pixels, and most of one Hannover bridge — `zt` already
*is* the deck, so the deck went in a second time on top of it, 4–5 m too high.
The nDSM per pixel is exactly deck top minus `zt`: pixels where the DTM carries
the deck measure near zero and stay terrain (`buildings.bridges.min_height`,
0.5 m), vegetation over a deck is refilled from the nearest pixel of the same
bridge, and a median over a 3 × 3 window within the bridge removes vehicles and
railings (`smooth_radius`). A deck the DTM carries over water still blocks the
flow below it — that would take cutting the DTM, which palmgeo does not do.

A domain containing bridges makes palm_csd print
`*_type arrays are not consistent! max: 2`. This is benign and internal to
palm_csd: it masks the surface types with `buildings_2d`, which excludes bridges,
but runs its final check against `building_type`, which includes them as type 7.
A bridge deck over a road legitimately carries both. Measured against the
criterion palm_csd actually uses, every cell carries exactly one type.

### Why maximum building heights fall as the grid coarsens

`buildings_2d` is averaged onto the domain grid over the building pixels of each
cell. The Neues Rathaus dome reads 97.6 m at 1 m but 39.5 m in a 32 m root,
because the dome covers ~50 m² of a 1024 m² cell. The alternative — taking the
maximum — would extrude the whole cell to dome height, which is aerodynamically
wrong. At 32 m the average is the honest representation; the fine nests keep the
detail. A coarse driver whose tallest building is *higher* than the fine one over
the same ground is a bug, not a feature.

---

# Part III — What PALM will not forgive

Everything so far has been about getting the data right. This part is about
getting the *geometry* right, which is a different problem: PALM enforces these
rules itself, by refusing to start, after the driver is built and the job is
queued. palmgeo checks them when the configuration loads.


## 9. Grid, nesting, and terrain

These are the rules that turn a plausible configuration into a driver PALM
rejects, or — worse — accepts and runs badly. palmgeo checks them when the
configuration loads, before any data is touched, and the distinction between
the two kinds matters:

* the **nesting** rules (§9.4) are hard errors, because a nest that violates
  them cannot be simulated at all, and the error names the nearest legal value;
* **2-smoothness** (§9.2) is a warning, because an unfavourable point count
  still produces a valid driver — it only makes the solver slow, and that is a
  trade you are allowed to make.

Nothing here is repaired automatically. Every repair moves a boundary, and
moving a boundary changes what is simulated.

### 9.1 `nx` is the last index, not the count

PALM counts grid points from 0, so a domain declared `nx = 1535` has **1536**
points and the driver's `x` dimension must be 1536. Off-by-one here is the most
common way to produce a driver that palm_csd writes happily and PALM refuses.

```
points = nx + 1          extent = (nx + 1) · dx
```

### 9.2 Point counts should be 2-smooth

PALM's multigrid pressure solver coarsens the grid by repeated halving. A point
count carrying a prime factor other than 2 or 3 stops the coarsening at the
level where that factor appears:

```
1536 = 2^9 · 3       coarsens to 3 x 3     — nine levels
1537 = 29 · 53       coarsens not at all   — one level
1500 = 2^2 · 3 · 5^3 coarsens to 375       — two levels
```

One point off is not a rounding error, it is the difference between a fast
solver and an unusable one. Useful counts: 512, 768, 1024, 1152, 1536, 2048,
3072, 4096.

This is a **warning, not an error** — unlike the nesting rules in §9.4. An
unfavourable count still produces a valid driver that PALM will read; what it
costs is solver speed, and that is a trade you are allowed to make. `load()`
warns and names the nearest usable counts; `verify_driver.py` reports the
leftover factor on the finished driver.

### 9.3 Coordinates are cell-centred

```
x(i) = (i + 0.5) · dx    for i = 0 … nx
```

`origin_x` / `origin_y` are the **lower-left corner** of the domain in the
projected CRS — the outer edge of the first cell, not its centre. The driver's
`x`/`y` variables are local offsets from that corner, so they always start at
`0.5·dx` regardless of where the domain sits.

### 9.4 Nesting alignment

For offline or online nesting, every one of these must hold. PALM enforces them
too — by refusing to start, after the driver has been built and the job
submitted — so palmgeo checks them when the configuration is loaded, which
turns a lost queue slot into a lost second:

* the child lies entirely inside its parent — *checked, hard error*;
* the parent-to-child resolution ratio is an integer (**4** is the usual choice,
  and what all the shipped examples use) — *checked, hard error*;
* `(origin_child − origin_parent) mod dx_parent = 0` in both x and y, so that
  nest borders fall on parent cell boundaries — *checked, hard error, and the
  message names the nearest origin that satisfies it*;
* every domain in one run shares one CRS, one `rotation_angle`, and one
  `origin_z` — *structural: there is one `grid.epsg` and one `csd.settings`
  block per run, so this cannot be violated from a configuration.*

None of the three is repaired automatically. Every repair moves a boundary and
therefore changes what is simulated, which is your decision and not palmgeo's;
the error says what the nearest legal value is and stops there.

`origin_z` is computed once over all domains in a single `palm_csd` invocation,
so building the whole set at once is what makes the vertical datum consistent.
Building a nest separately later will not reproduce it.

### 9.5 Terrain blending at nest boundaries

```yaml
csd:
  interpolate_terrain: false    # the default
```

A nest and its parent derive their terrain from the same source data but on
different grids, so they do not agree cell for cell at the shared boundary.
`interpolate_terrain` decides what palm_csd does about that.

**It is not only a boundary treatment, which is the thing to understand before
choosing.** `create_driver.py:479-488` shifts the *whole nest* vertically before
it blends anything:

```python
dz_mean = mean(zt_child) - mean(zt_ip)   # zt_ip = parent terrain, interpolated
domain.zt.values -= dz_mean              #        AND snapped to the parent's dz
```

Only then is the 50-cell margin merged. So turning it on moves every cell of the
nest, not just its edge.

Both effects were measured on a controlled A/B over a 32 / 16 / 4 m grid — same
work tiles, all 56 input GeoTIFFs bit-identical, one setting changed. `zt` was
the *only* variable that differed at all; buildings, surface types, LAD and
`surface_fraction` came out identical.

| | `false` (default) | `true` |
|---|---|---|
| terrain in the nest | as its own source data measures it | shifted bodily: **+1.14 m** (16 m nest), **+3.91 m**, **+5.12 m** (4 m nests) |
| outer 50 cells | its own orography | the parent's, flattened onto the parent's `dz` |
| step at the boundary | median 4.00 m, p95 8.00 m for the 4 m nests; up to 16 m at the 16 m nest | 0.00 m median and p95 everywhere |
| against raw, unquantised parent `zt` | median 3.62 m, p95 6.92 m | median 6.13 m, p95 15.00 m — *worse* |

**The default is `false` because the vertical shift is a global, silent change
to data you supplied.** A driver whose terrain sits five metres above what DGM1
measured — over half the total relief of that nest — should be something you
asked for, not something a default did to you. `zt` in a blended nest is not the
surveyed height, and any comparison against a levelled station height has to
account for it.

**That is not a recommendation to leave it off.** On the grid measured above,
blending was still judged the better choice: a step at the inflow sits on every
side of the domain for the whole simulation, in the place where the nest can
least absorb it, while the shift leaves the terrain *shape* exactly as measured
and costs a few metres of modelled atmosphere depth. Making the child's surface
sit at the parent's height is the whole point of the option, not a side effect.

Turn it on when the nest boundary crosses real orography. Leave it off when the
nests are flat relative to the parent's `dz`, or when absolute elevation matters.
If the shift comes out large, the actual culprit is usually the parent's `dz`:
in the A/B the 16 m parent quantised 8 × 8 km into three terrain levels, and both
nests inherited that step as their boundary condition. A finer parent `dz` shrinks
the shift and the step at once.

The switch is global and applies to every nest. **The root is always `false`**:
it has no parent to blend into, and asking for it there is refused rather than
silently ignored. A domain block may override the global value individually,
which is how you blend one difficult boundary and leave the others alone.

`verify_driver.py` measures the step either way. With blending on it is a
criterion — the reference run of §10 matches its parent to 0.00 m median and
0.00 m p95 at all five boundaries. With blending off the step is the expected
result, so the tool reports its size as `INFO` and does not fail.

### 9.6 The root driver's file name

`palm_csd` names every output `<file_out>_<domain>`, with no special case for
the root. PALM reads the root domain's driver as `<jobname>_static` — **without**
a suffix, because the suffix is what marks a nest. So a run named
`hannover_static` produces `hannover_static_root`, which has to become
`hannover_static`. `tools/finalize.py` does exactly that and nothing else.


---


## 10. Verification and worked examples

Three things are checked, at three different scopes, and they are not
interchangeable.

| | what it answers | cost |
|---|---|---|
| `tools/selftest.py` | has *palmgeo* stopped producing what it claims? | seconds, or minutes with `--full` |
| `tools/verify_driver.py` | is *this driver* structurally sound? | seconds |
| `tools/check_classification.py` | is this driver's classification *plausible*? | seconds |

### The self-test

```bash
python repository/tools/selftest.py           # offline, seconds
python repository/tools/selftest.py --full    # rebuilds the example, minutes
```

The committed example under `static_driver/example/` is the reference: it was
built by this pipeline, verified, and its numbers are written down. That makes
the useful question not *did a run finish* but *did it produce the same thing
as last time* — a test that only checks for a non-empty result passes happily
while the classification drifts underneath it.

The quick pass checks everything answerable from the repository itself: every
shipped configuration loads, every region file's licence ids exist in the
registry and name no absolute local path, `pyproject.toml` and
`requirements.txt` still agree, the committed drivers open and carry their
variables, the nest still contains individual trees, and each domain's
manifest still matches the rasters beside it. No network, so a CI job can run
it anywhere.

`--full` rebuilds the example into a scratch directory and compares the result
against the committed one variable by variable — cell counts, minima, maxima
and means, within a tolerance of 1e-3. The tolerance is deliberately not zero:
GDAL resampling and floating-point summation are not bit-stable across library
versions, and a test that fails on the tenth decimal is one people learn to
ignore. It is small enough that a changed threshold or a changed source moves
it.

Exit status is 0 or 1, so it can be the thing a release check runs.


### The example, which is also the smoke test

`repository/config/example.yml` covers 9 km² of Hannover city centre — 9 work tiles at
1 m, a 2 km root at 32 m from a second, 4 m work tier, a 1 km nest at 8 m and a
256 m nest at 2 m. It was run against the live
services and then through `palm_csd` 25.10 (the `vectors` and `buildings`
counts are from the build of 2026-11-04):

```
plan       9 tiles, 9 km², leaf-on campaign 2022-05-08, leaf-off 2025-03-06
rasters    9 tiles in 9 s, 0.07 GB transferred (orthophotos streamed, not stored)
vectors    ALKIS Region Hannover (330 MB) + Basis-DLM (536 MB) fetched,
           clipped to 4 606 land-use polygons / 9 787 footprints / 1 931 roads,
           archives deleted
buildings  8 485 kept: 9 126 heights measured from the nDSM, 16 from LoD2
           where the nDSM has none, 666 dropped as implausible; 32 bridge decks
compose    9 tiles, 18 distinct PALM classes, 100 % of pixels classified
domains    exactly one surface type on all 36 864 pixels of the three domains,
           every type value inside PALM's valid range, no terrain gaps
palm_csd   all three drivers written; N03 carries buildings_3d, lad, bad,
           tree_id, tree_type
```

Classification over the test area (measured on a build before v1.1; the shares
move by tenths of a percent as the LGLN updates its data), a plausible
signature for a German city centre with a lake and a park:

```
pavement: asphalt_concrete_mix          28.2 %
building                                25.9 %
vegetation: deciduous_broadleaf_trees   15.7 %
vegetation: short_grass                  7.3 %
pavement: paving_stones                  6.5 %
water: lake + river                      4.6 %
vegetation: mixed_forest_woodland        1.3 %
vegetation: evergreen_needleleaf_trees   0.3 %
…
```

The per-domain vegetation setting behaves as specified: the `unresolved` 32 m
root keeps tall vegetation as an LSM surface type (18.5 % of its cells), while
the `mixed` 8 m nest turns the same canopy into a 3D leaf area density field on
four levels and the `resolved` 2 m nest adds individual trees, both with short
grass beneath. Neither nest keeps forest as an LSM type; at 8 m PALM would
refuse that (§7).

### The reference run

Six nested domains at 32 / 8 / 2 m in the Hannover region, built to a grid
specified externally and verified against it. Its configuration belongs to a
study and is not published; the measurements are. To write a configuration of
your own, start from `repository/config/template.yml`.

| domain | role | parent | dx | extent [m] | points | nx, ny |
|---|---|---|---|---|---|---|
| `root` | root | — | 32 | 49152 × 49152 | 1536 × 1536 | 1535, 1535 |
| `N02` | parent | root | 8 | 8192 × 8192 | 1024 × 1024 | 1023, 1023 |
| `N03` | parent | root | 8 | 8192 × 8192 | 1024 × 1024 | 1023, 1023 |
| `N04` | child | N02 | 2 | 3072 × 3072 | 1536 × 1536 | 1535, 1535 |
| `N05` | child | N02 | 2 | 2048 × 2048 | 1024 × 1024 | 1023, 1023 |
| `N06` | child | N03 | 2 | 3072 × 2048 | 1536 × 1024 | 1535, 1023 |

Result, measured on the written NetCDF files:

```
domain   vegetation  pavement  buildings  water   max height  3d lvl  LAD lvl  trees   street_type
root         89.6 %     7.8 %      1.4 %   1.2 %     112.1 m       —        —      —        4.5 %
N02          45.5 %    35.4 %     15.8 %   3.4 %      92.0 m      13        6      —       10.9 %
N03          91.2 %     4.8 %      2.6 %   1.4 %      27.7 m       5        6      —        2.0 %
N04          40.5 %    39.2 %     17.0 %   3.3 %      78.1 m      41       26  16 622       6.6 %
N05          31.2 %    47.1 %     21.5 %   0.2 %      56.5 m      30       19   5 974       5.6 %
N06          80.7 %    12.0 %      6.7 %   0.7 %      27.2 m      15       24   6 448       3.1 %
```

Every cell of all six domains carries exactly one surface type; `zt` covers
100.0000 % of every domain; `origin_z` is identical (5.4454 m) in all six files.

**Nest-boundary terrain matches to 0.00 m median and 0.00 m p95 at all five
boundaries**, because this run sets `csd.interpolate_terrain: true` explicitly
(§9.5; the default is `false`) — but only when measured correctly. The setting
makes palm_csd interpolate the parent's `zt` onto the child grid and merge the
two over a 50-cell margin, and it deliberately snaps the interpolated terrain to
the **parent's** `dz` (`create_driver.py`: `zt_ip = height_to_z_grid(zt_ip,
parent.dz)`), because PALM will discretise the parent's terrain that way at run
time anyway. A root file carrying continuous terrain
(`use_palm_z_axis: false`) therefore appears to differ from its 8 m nest by a
whole 32 m step where there is none. `verify_driver.py` quantises the parent
before comparing, and says so in its output.

Buildings and land use are **not** blended across nest boundaries, and cannot be
without discarding the nest's resolution. Child and parent are aggregated
independently from the same work rasters using the same classification rules, so
they never disagree in *semantics* — but a building edge on a 2 m grid is
necessarily drawn differently from the same edge on an 8 m grid. That is what
nesting is for.

### Two small runs

Neither is in the repository: a run configuration stays with its project
(§2), and a driver is a few hundred megabytes to a gigabyte.

| domains | what it demonstrated |
|---|---|
| 512 m @ 1 m, single domain over the Neues Rathaus | per-pixel building height: the Rathaus dome resolved to 97.6 m, 28 256 distinct heights, 417 single trees |
| 2048 m @ 16 m / 512 m @ 4 m / 128 m @ 1 m | concentric three-level nest, one vegetation mode per level; the courtyard case of §8 |

---

# Part IV — Operating it


## 11. Running it

Run from the project root with `export PYTHONPATH=repository`. The working
directory does not actually matter — paths inside a config resolve against the
project root, not against the shell — but the commands read most clearly from
there. `CFG` below is `static_driver/<run>/config/<run>.yml`.

```bash
python -m palmgeo.cli sources    $CFG   # what data exists here, under what licence
python -m palmgeo.cli plan       $CFG   # dry run, no raster I/O
python -m palmgeo.cli vectors    $CFG   # cadastre, roads, bridges
python -m palmgeo.cli rasters    $CFG   # terrain, nDSM, NDVI per tile
python -m palmgeo.cli calibrate  $CFG   # NDVI histogram → threshold
python -m palmgeo.cli buildings  $CFG   # footprint heights
python -m palmgeo.cli trees      $CFG   # single-tree segmentation
python -m palmgeo.cli compose    $CFG   # per-pixel PALM classification
python -m palmgeo.cli report     $CFG   # unmapped source vocabulary
python -m palmgeo.cli domains    $CFG   # aggregate onto each domain grid
python -m palmgeo.cli csd        $CFG   # write the palm_csd YAML
python -m palmgeo.cli inspect    $CFG   # QGIS inspection layer

python -m palmgeo.cli all        $CFG   # all of the above, in order
```

`sources` and `plan` answer different questions and are both cheap. `sources`
asks *what may I use and what will it become*; `plan` asks *what will it cost*.
Neither downloads raster data.

**`vectors` and `rasters` refuse to start until the sources have been
approved** — `python -m palmgeo.cli sources $CFG --approve`. The approval is
recorded in `state.json` as a fingerprint over the chosen source per role and
the licences those sources carry, so changing a surface model or a licence
changing under a source invalidates it and asks again. Changing a domain
origin or a threshold does not, because neither changes what may be done with
the result.

then PALM's own tool and the two post-steps:

```bash
MPLBACKEND=Agg env -u PROJ_DATA -u PROJ_LIB palm_csd --pdf static_driver/<run>/config/<run>_static_csd.yml
python repository/tools/finalize.py         $CFG
python repository/tools/verify_driver.py       $CFG   # structure
python repository/tools/check_classification.py $CFG   # plausibility
```

`verify_driver.py` proves the driver is *structurally* sound. It cannot tell
you the classification is *right*: it passed on a Manhattan driver in which
half of Central Park had come out asphalt, because a uniformly wrong map is
structurally perfect. `check_classification.py` closes that gap by
cross-tabulating each land-use class against the surface types actually
produced, and failing a class whose rule can only produce vegetation but whose
ground came out pavement. Run both.

`palm_csd --pdf` (or `--png`) writes one sheet per domain next to the driver,
showing buildings, pavement, water and vegetation split into flat and resolved.
That is the artefact to look at when judging a driver, because it renders what
PALM will actually read rather than what went in. `-s` shows the plot
interactively and needs a display.

### Running without a network

palmgeo streams by default, and streaming needs the network *while* the rasters
are processed, not only beforehand. Where that is not acceptable — a compute
node with no outbound route, a field deployment, an archive that has to be
rebuildable years later — there are two boundaries to choose between, and they
differ by two orders of magnitude.

**The cheap boundary: after `rasters`.** Acquisition is `sources`, `plan`,
`vectors`, `rasters` and `buildings`; everything after them reads only the run's
own working directory. This needs no configuration at all — it is simply true —
and it is what most people mean when they say the heavy part should run
elsewhere. For a nested city run it means carrying **3.5 GB** of derived data to
the offline machine.

**The expensive boundary: after `fetch`.** Every source file the run will touch
is downloaded whole, first; afterwards the entire pipeline runs with no network
at all.

```bash
# once, with network. keep_downloads must be on, or the stages delete
# each file again after reading it.
python -m palmgeo.cli fetch $CFG --check    # what is missing, download nothing
python -m palmgeo.cli fetch $CFG

# then, anywhere
#   runtime.offline: true
python -m palmgeo.cli all $CFG
```

The price is the streaming saving, given up in full. Measured over Lower
Saxony:

| | streamed | fetched whole |
|---|---|---|
| Oldenburg, 1246 work tiles | 11 GB | **345 GB** |
| Hannover, 2700 work tiles | 25 GB | **747 GB** |

A region served by a staging port (`repository/ports/`) is a third case and the
happy one: its sources are already local files, so it is offline-capable from
the start. `palmgeo fetch --check` says which case a run is in.

**`runtime.offline: true` refuses the network rather than failing at it.** It
does not retry, and it does not discover the gap one file at a time: each
acquisition stage checks everything it will need before it starts and aborts
with the whole list. On a run whose sources are absent that is one message in
six seconds, where the retry-and-backoff path would have spent hours saying the
same thing a few thousand times.

Offline mode also stops preferring to stream a file that is already in the
cache — which is otherwise slower and moves bytes for nothing — and it accepts
a source index older than its 30-day refresh interval rather than aborting,
with a warning that it may not list recent acquisitions.

### Resumability and stage fingerprints

Every stage is resumable: work is marked complete per tile, so re-running after a
network failure redoes only what is missing.

That marker carries a **fingerprint of the settings the stage actually reads**
(`STAGE_KEYS` in `config.py`), so editing a threshold invalidates exactly the
tiles that depended on it and nothing else — change `buildings.min_height` and
the composition re-runs; change `csd.attributes.author` and nothing does. Without
this a re-run silently keeps the old tiles and a config change appears to have no
effect, which is a genuinely hard thing to notice: the pipeline reports success
and writes a driver, just not the one that was asked for. Markers written before
fingerprinting existed carry no hash; they are trusted once and re-stamped,
rather than forcing a multi-gigabyte re-download.

The composition also hashes the **content of the vector layers** it reads
(`COMPOSE_VECTORS` in `steps/compose.py`). Re-running `vectors` or `buildings`
rewrites them — with new building IDs at the least — while every setting stays
the same, and before this the old composition survived: the example once
carried a root built from new vectors on top of nests composed from the old
ones. A GeoPackage records when it was written, so any rewrite recomposes, even
of identical data; that is the safe direction, and composing is cheap next to
the tiles. Compose markers of the previous generation (`v2:`) are trusted once
and re-stamped like the oldest ones; `clean <cfg> surface` recomposes a project
whose vectors may have been rebuilt since.

### Cleaning

```bash
python -m palmgeo.cli clean $CFG cache|tiers|surface|vector|work|out
```

`surface` drops only the composition, keeping the expensive downloaded and
derived tiles — which is what a *code* change (as opposed to a config change)
calls for.

### The two empirical inputs

* `calibrate` prints an NDVI histogram over sampled tiles. NDVI is computed from
  8-bit digital numbers, not calibrated reflectance, so the threshold is
  campaign-specific; a bimodal histogram with a clear trough means the separation
  is usable.
* `report` lists every source vocabulary value that found no explicit mapping
  rule, so `mapping_*.yml` can be completed instead of guessed at.

### PALM namelists

`tools/update_p3d.py` brings an existing set of p3d namelists in line with a
palmgeo grid — `nx`, `ny`, `dx`, `dy`, the header comment, and `domain_layouts`
in the root namelist. Everything physical is copied byte-identical.

```bash
python repository/tools/update_p3d.py $CFG \
       --layout static_driver/<run>/config/<run>_layout.yml --src <dir with the old namelists>
```

The layout file names the PALM domains and their PE budget. PALM splits a domain
so that `npex` divides `nx+1` and `npey` divides `ny+1`; the tool refuses a count
it cannot split and prints the valid ones on either side, rather than choosing
for you. Note that 2-smooth point counts have no factor 5, so a PE count like 200
or 3600 becomes impossible the moment the grid is made multigrid-friendly.


---


## 12. Storage and memory behaviour

The pipeline is built to run on far less disk than it touches.

* **Work tile = 1 × 1 km.** At 1 m that is a 1000 × 1000 array per product, so a
  worker's peak memory is tens of megabytes no matter how large the AOI.
* **Bounded cache.** `cache_budget_gb` is a hard cap, LRU-evicted over
  unreferenced entries; files being read are pinned. Source tiles are released as
  soon as they have been consumed.
* **Streaming.** With `runtime.access: auto`, a source is read directly from the
  server through `/vsicurl/` whenever an overview level already satisfies the
  requested work resolution (≥ 4× coarser). Reading a 2 km / 20 cm orthophoto
  onto a 1 km / 1 m grid transfers roughly 1/25 of the file and writes nothing to
  disk. Measured: **1.7 s per km²** for two bands from the live server.
* **Region-wide monoliths** (in Lower Saxony: Basis-DLM 536 MB, whole-Landkreis
  ALKIS up to 330 MB) are downloaded once, clipped to the AOI in stripes, and
  deleted.
* **Aggregation is row-chunked**, so building a 1536 × 1536 domain from a
  49 000 × 49 000 work mosaic costs the same memory as a small one.

For scale: the six-domain reference run moved ≈ 30 GB, almost all of it streamed
from COG overviews and never written to disk, and produced 2.7 GB of driver in
1 h 30 min wall clock.


---


## 13. Assumptions

Every one of these changes the result and none is derivable from the open data.
**Report the relevant ones when handing over a driver.** Where a product or a
file is named, it is the Lower Saxony reference configuration's; the
assumption itself holds for any region.

1. **Building extent is the cadastre's, height is the measurement's.** A building
   the cadastre does not map does not exist for the driver, however clearly the
   surface model shows it; a footprint the cadastre still maps is carved out only
   where the surface model reads below `min_height`. Nothing creates buildings
   from the nDSM alone, because there is no reliable way to tell a building from
   a hedge or a lorry without the footprint.
2. **A courtyard is open only above its measured floor.** Openness is judged
   against terrain, so a roofed or raised inner court stays building up to the
   height the surface model gives it. This is the physically right answer when
   the surface model is right, and the failure mode is silent when it is not:
   neither DGM1 nor DOM1 can see into a court too narrow for the scan geometry,
   and both would then agree on a floor that is not there. Where a court matters,
   check it against a second, independent surface model — the photogrammetric
   `bdom20` is derived from different sensors and epochs than the ALS-based
   `dom1`.
3. **Building age.** PALM's `building_type` encodes residential/office × three
   construction periods. Few cadastres carry a construction year (neither ALKIS
   nor LoD2 does), so `mapping_<region>.yml: building_type.default_age` is
   applied uniformly (`1951_2000` in Lower Saxony). Supplying a building-age layer through `extra_inputs` is the
   only way to do better.
4. **NDVI and greenness thresholds.** Uncalibrated 8-bit DN, aggregated to work
   resolution before the ratio is formed (averaging radiances first is less noisy
   than averaging per-pixel NDVI). The defaults under `classification.vegetation`
   and the phenology cut-offs are starting points calibrated on one Hannover
   campaign — run `calibrate` for yours.
5. **Tree allometry.** Trunk diameter is derived from measured height and leaf
   area index from trunk diameter, with one relation per broad class. Species are
   not identified: every conifer and every broadleaf becomes the one palm_csd
   tree type the region declares for it (`single_trees` in
   `mapping_<region>.yml`; `Abies` and `Fagus` in Lower Saxony), which fixes
   crown shape and the LAD profile. Feed a real tree register through
   `vegetation.trees.mode: external` where species matter.
6. **Default tall vegetation** where phenology is unavailable is the region's
   `default_high_vegetation` — `deciduous_broadleaf_trees` in Lower Saxony, as
   the most common tall vegetation there.
7. **Soil type** is the constant `medium_fine` (PALM 3), refined to `organic`
   only for bog and marsh land use. No area-wide open soil texture layer is wired
   in; the LBEG BÜK50 WFS is the obvious upgrade.
8. **Pavement sub-types** per land-use class (asphalt for roads, paving stones
   for residential yards, gravel for railway ballast) are plausible defaults, not
   observations. They are one line each in `mapping_ni.yml`.
9. **`street_type` codes.** The Basis-DLM `WDM` values in `mapping_ni.yml` were
   verified against the delivered `ver01_l.dbf`; note the field is upper case in
   the delivery and lower case in the specification, which field resolution
   handles. PALM's `street_type` enum *is* the OSM `highway` taxonomy, so
   OpenStreetMap is an equally valid source; see [`PORTING.md`](PORTING.md).
10. **Bridge deck height** is the nDSM per pixel, smoothed within the deck
    (§8), and `bridge_depth` a uniform 3 m. Where the DTM was not cut below a
    deck, that part stays terrain: no flow passes under it.
11. **LAI** is left to palm_csd's height-based estimate. An NDVI→LAI relation is
    stubbed but disabled: a poorly calibrated one is worse than none.


---


## 14. Known limitations

* **NDVI seams on the orthophoto tile lattice** (2 km in Lower Saxony) where no
  single acquisition campaign covers the AOI. `plan` reports which case applies.
* **Aggregation is a per-class mode.** A cell that is 30 % grass, 30 % trees and
  40 % asphalt becomes asphalt although vegetation totals 60 %. The two-stage
  group-then-subtype vote fixes the group but not this within-group split.
* **Shadows** in the orthophoto depress NDVI, so north-facing façade strips and
  deep courtyards may classify as sealed. A shadow mask from the nDSM and the sun
  position at acquisition time would fix it; not implemented.
* **Vegetation under bridges and overhanging canopy** is not modelled; the nDSM
  is a 2.5D surface.
* **Water level.** The DGM over water is the water *surface* at flight time, so
  `zt` under rivers is not the bed.
* **Polygonised inspection output** is produced block-wise and not dissolved
  across block borders; it is for inspection, not for area statistics.
* **Coverage.** Domains reaching beyond the edge of the source coverage (a state
  or city border, usually) lose everything but terrain interpolation. `domains`
  aborts if more than 2 % of a domain's terrain is missing, and fills smaller
  gaps by nearest neighbour with a warning.
* **`albedo_type` is not written.**
* **Calibrated on Lower Saxony.** The NDVI, greenness and phenology thresholds
  and the broadleaf/conifer vote of the tree detection were tuned on 8-bit
  orthophotos of Lower Saxony, and the trunk-diameter and leaf-area equations
  for single trees are Central European. Elsewhere they are a starting point:
  `palmgeo calibrate` and `tools/check_classification.py` are how to check
  them (§11, §10).
* **Delivery formats.** A cadastre is read from a zipped GeoPackage, a 3D
  building model and a road network from zipped shapefiles, rasters from
  GeoTIFF. A region that publishes something else needs a staging step that
  converts it (`PORTING.md`).
* **Land-use rules match by prefix.** That suits hierarchical codes (ALKIS);
  with flat vocabularies such as OSM tags a short key also catches every
  longer key it is a prefix of. The mapping key `by_gfk` keeps its ALKIS name
  for any region's building function.
* **Seasons are configured, not inferred.** The orthophoto epochs
  (`sources.ortho_epochs.*.months`) in the template are for the northern
  hemisphere; south of the equator leaf-on and leaf-off swap.


---


## 15. Troubleshooting

**`The WKT could not be parsed. OGR Error code 6`**
A stale `PROJ_DATA` or `PROJ_LIB` in the shell shadows the PROJ database rasterio
ships with. It kills **both** palmgeo and palm_csd and reads like a data problem.
palmgeo detects and neutralises it at import and warns; palm_csd does not, so run
it as `MPLBACKEND=Agg env -u PROJ_DATA -u PROJ_LIB palm_csd …` or unset the variables.

**A stage finishes in a second and changes nothing.**
Its tiles are already marked complete. If a config change should have
invalidated them but did not, the key you edited is not in that stage's
`STAGE_KEYS` (`config.py`) — add it. `clean surface` forces recomposition after a
code change.

**`*_type arrays are not consistent! max: 2` from palm_csd.**
Benign, and only appears with bridges. See §8.

**palm_csd runs out of memory.**
A domain field is `(nx+1)·(ny+1)·4` bytes and palm_csd holds several at once. A
20 km domain at 2 m is 400 MB per field. Either coarsen, or split the area into
more, smaller nests.

**The driver looks right but PALM rejects the grid.**
Check §9.1 first: `nx` in the namelist must be one less than the driver's `x`
dimension.

**A nest boundary appears to jump by a whole parent `dz`.**
It is almost certainly the quantisation described in §10, not a defect. Measure against the
parent terrain rounded to the parent's `dz`; `verify_driver.py` does. If the jump
is real, check `csd.interpolate_terrain` (§9.5). It is **off by default**: the
nest then takes its terrain purely from the source data, and a step of up to one
parent grid spacing at the boundary is the expected result, not a defect. Turn it
on if that step is a problem for your case.

**`interpolate_terrain cannot be true` on the root.**
The root has no parent to blend into. Remove the key from the root's domain block;
the global switch already skips it.

**Downloads are slow or the transfer estimate looks enormous.**
Check that `runtime.access` is `auto` and that the work resolution is at least
4× coarser than the source, so COG overviews can serve the read. A `fixed: 1.0`
work resolution over a large root forces full-resolution downloads.


---


## 16. Versioning and attribution

Development happens on `main`. A state worth referring to gets a tag and a
GitHub release, and every release is archived on Zenodo with its own DOI
(see [Citing palmgeo](#citing-palmgeo)).
`repository/palmgeo/__init__.py` carries `__version__`, and it should match the
tag — it is also what identifies palmgeo to the servers it downloads from, so a
version that has drifted misinforms exactly the operator it is meant to inform.

Versions before 1.2.0 were developed in a private repository; their history
is not published, and the public repository starts at `v1.2.0`. They are
listed for what they introduced.

| version | what it was |
|---|---|
| `v0.1` | first working state: one region (Lower Saxony), the full pipeline |
| `v0.2` | two non-German regions, staging tools, plausibility checks |
| `v1.0` | reviewed end to end, documented for readers rather than for its author, and testable against itself |
| `v1.1` | bridge decks measured per pixel; `unresolved` vegetation restricted to coarse domains, with a warning before PALM's LSM0048 abort; `compose` re-runs when its vector inputs change; the example rebuilt as three domains under a 32 m root |
| `v1.2.0` | first public release, under the GPL and citable through Zenodo. No region default left in the code: the CRS, the cadastre and 3D-model column names and the single-tree species are declared by the region files; the New York and Egypt ports take their area on the command line |

### What 1.0 means here

Not "finished" — [section 14](#14-known-limitations) is still a list. It means
four things that were not true before:

**Every module has been read and argued with.** Four review passes over the
core, the acquisition, the derivation and the tools, carried out by AI coding
assistants (see [How palmgeo was developed](#how-palmgeo-was-developed)). They found a silent
corruption of any orthophoto that was not natively 8-bit, a cache that could
evict data it had never downloaded and could not restore, a stage fingerprint
blind to the region file, a cadastre merge that never de-duplicated, and a
generated configuration that named the author's home directory. Each fix was
recorded with what was wrong and how it was measured, in the private history
that preceded the first public release.

**Nothing in the repository points at a machine or a private study.** No
absolute path, no contact address, no research programme, no run configuration
of the author's. Who wrote palmgeo is stated where a citation needs it
([`CITATION.cff`](CITATION.cff)). What stays is measurements and public places,
because a justification without its number is not one.

**It says what it rests on.** Every driver carries its data basis, the licences
those sources impose, the conclusion drawn from them, and how old the data is —
the last of which explains more missing buildings than anything else.

**It can be checked.** `tools/selftest.py` compares a rebuild against the
committed example variable by variable, so a drift in the classification fails
a test instead of quietly shipping.


### Citing palmgeo

If you use palmgeo, or a driver built with it, in your work, please cite it.
[`CITATION.cff`](CITATION.cff) holds the reference; GitHub shows it under
"Cite this repository". Every release is archived on Zenodo. The concept DOI
[10.5281/zenodo.23973676](https://doi.org/10.5281/zenodo.23973676) always
resolves to the newest version; in a publication, cite the DOI of the version
you used, so a reader can retrieve exactly that code — for v1.2.0 that is
[10.5281/zenodo.23973677](https://doi.org/10.5281/zenodo.23973677).

### Licence

Copyright (C) 2026 Leibniz Universität Hannover, Institut für Meteorologie und
Klimatologie (IMUK).

palmgeo is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; see [`LICENSE`](LICENSE) for details.

### Data

**Data** is not covered by that licence. Every source carries its own, and a
driver inherits all of them. `palmgeo sources` reports them before anything is
downloaded and writes them, with the attribution each one requires, into
`docs/DATA-SOURCES.md` next to every driver. Keep that file with the driver, and
put the attribution into the driver's `source` attribute
(`csd.attributes.source` in the run config) when you pass a driver on.

The Lower Saxony sources are © LGLN under `dl-de/by-2-0`
(<https://www.govdata.de/dl-de/by-2-0>), obtained via
<https://opengeodata.lgln.niedersachsen.de/>. Any driver built from them, and any
publication using one, carries the attribution *"© LGLN \<year\> —
dl-de/by-2-0"*. The example under `static_driver/example/` is built from these
sources and carries *"© LGLN 2026 — dl-de/by-2-0"* in its `source` attribute.

`palm_csd` is part of the PALM model system (Leibniz Universität Hannover,
Technische Universität Berlin) and is GPLv3. palmgeo does not include or link it;
it writes its input and is run alongside it.

### Credits

The single-tree detection in `repository/palmgeo/steps/trees.py`, the
greenness criterion in `steps/compose.py` and the canopy height model in
`steps/tiles.py` are ports of scripts by **Joana-Pauline Rentz**, written for
her Bachelor's thesis at the Institute of Meteorology and Climatology, Leibniz
University Hannover:

> Rentz, J.-P. (2026): *Generierung von dreidimensionalen Einzelbaumrepräsentationen basierend auf Fernerkundungsdaten* [Generating three-dimensional single-tree
> representations from remote sensing data]. Bachelor's thesis, Institute of
> Meteorology and Climatology, Leibniz University Hannover.

The original scripts are not part of this repository; the method, every one of
its tuned constants and the reasoning behind them are documented in the
provenance section of each of those three modules. If your work relies on
palmgeo's single trees, please credit the thesis as well.

### How palmgeo was developed

palmgeo was developed with substantial help from AI coding assistants. They
wrote most of the code and carried out the module reviews described under
"What 1.0 means here". The goals, the design, the methodological decisions and
the checking of the results — the drivers, their plausibility, and their
behaviour in PALM — are the author's, and the author is responsible for the
result.

---

# Appendices


## Appendix A — Driving palmgeo with an agent

A static driver is not derivable from a place name. Roughly a dozen decisions
have to come from the person who will run the simulation, and getting one of
them wrong is expensive — a full run over a city takes hours and tens of
gigabytes of transfer.

### A.1 Ask for anything on this list that was not in the prompt

Ask **before** starting the pipeline, in one batch rather than one at a time.
Where the user has no opinion, offer the default in the right-hand column and
say what it implies.

| # | What you need | Why it cannot be guessed | Default if the user has no preference |
|---|---|---|---|
| 1 | **Region / country** | Decides which `sources_*.yml` applies. Lower Saxony works out of the box; New York City and Egypt need their staging scripts in `repository/ports/` first; anywhere else needs new source and mapping files — see Appendix C and `sources_global.yml`. | none — **blocking** |
| 2 | **Location** — a place name, a centre coordinate, or a lower-left corner, plus the CRS | `origin_x`/`origin_y` in a domain block are the **lower-left corner** in the projected CRS. A place name is not a domain. | ask; convert with `palmgeo lonlat <lon> <lat> --epsg <n>` |
| 3 | **Number of domains and the nesting tree** | Each nest costs its own rasterisation. | one root plus one nest |
| 4 | **Horizontal resolution per domain** (`pixel_size`) | Drives everything: work tier, transfer volume, runtime. | 32 m root, 8 m parent, 2 m child |
| 5 | **Extent per domain**, or the point count directly | See §9.2 — the point count is the thing PALM cares about. | derive from extent, then round **up** to the next 2-smooth count and say so |
| 6 | **Vertical spacing `dz` per domain** | Written into the driver for `use_palm_z_axis` and used for the terrain quantisation at nest boundaries. | `dz = pixel_size` |
| 7 | **Vegetation treatment per domain** | `unresolved` / `mixed` / `resolved` — see §7. `resolved` triples the runtime of the finest tier. | `unresolved` root (only at ≥ 16 m, §7), `mixed` parents, `resolved` at ≤ 2 m |
| 8 | **Terrain blending at nest boundaries** (`csd.interpolate_terrain`) | Not just a boundary treatment: turning it on shifts every nest bodily in the vertical, by metres. Off leaves a terrain step at the inflow instead. Neither is free — §9.5. | **`false`**; say you are using the default, and that `true` is the right choice where a nest boundary crosses real orography |
| 9 | **Season and `origin_time`** | Sets `season: summer\|winter` and the driver's time attribute; affects LAI and the leaf-on/leaf-off epoch choice. | summer, 12:00 UTC on 1 July of the most recent full year |
| 10 | **Where `palm_csd` is** and which Python runs it | palmgeo writes palm_csd's input; it does not run it. The two must share a working geo stack. | try `which palm_csd`; ask if absent |
| 11 | **Where `workdir` and `outdir` go** | A city-scale run wants 50–100 GB of scratch and can produce several GB of driver. | `static_driver/<run>` with `work/` inside it; override `workdir` to move the scratch to another disk |
| 12 | **Attribution metadata** — author, institution, contact, site name | Written into the driver's global attributes and into every file that leaves the machine. | ask; do not invent an affiliation |
| 13 | **PE budget per domain**, only if p3d namelists are wanted | An allocation decision. See `tools/update_p3d.py`. | ask; do not choose cores for the user |

### A.2 What you should work out yourself, without asking

* **Point counts from extent and resolution**, including the 2-smooth rounding
  (§9). Report the result; do not make the user do arithmetic.
* **Nest origin snapping** — a child origin must sit on a whole number of parent
  cells. If the user's requested origin is off by less than one parent cell,
  snap it and say by how much. If it is off by more, ask (§A.3).
* **Work tier resolution**, tile lists, epoch selection, cache budget: all
  computed by `plan`.
* **Which stages need re-running** after a config change — the stage
  fingerprints handle it (§11).
* **The whole classification chain**: thresholds have documented defaults, and
  `calibrate` closes the loop on the one that is campaign-specific.

### A.3 What you must never decide alone

* **Moving an origin or changing a point count the user specified.** Both change
  what is simulated. If the requested grid violates a hard rule from §9, say
  which rule, show the nearest grid that satisfies it, and ask.
* **Dropping data.** If footprints, tiles or trees are discarded, report the
  count and the reason. A run that silently drops 3 % of buildings is worse than
  one that fails.
* **The PE decomposition.** How many cores to spend is not a technical question.
* **Approving the data sources.** `palmgeo sources --approve` accepts, on the
  user's behalf, every licence obligation the finished driver will carry — a
  non-commercial clause, a share-alike clause, a restriction on redistributing
  the data. Show them the report and let them say yes. This is the one gate in
  the pipeline that exists for a reason outside the software.
* **Publishing or transmitting anything.** The driver may encode an author's
  name, a site, and a research programme.

### A.4 The loop

```
ask (§A.1)  ->  cp repository/config/template.yml static_driver/<run>/config/<run>.yml, then edit
            ->  palmgeo sources     # what data exists, under what licence, for what
            ->  show the user the licence position; get --approve from THEM
            ->  palmgeo plan        # read the output; it is a dry run for a reason
            ->  confirm tiles, transfer estimate and epochs with the user
            ->  palmgeo all
            ->  palm_csd --pdf
            ->  tools/finalize.py
            ->  tools/verify_driver.py     # must end in ALL CHECKS PASSED
            ->  look at the PDF plots
            ->  report: what was built, what was dropped, what is assumed
```

The run configuration belongs in `static_driver/<run>/config/`, not in
`repository/config/`. It names a site, an author and a study; `static_driver/`
is git-ignored, so a run stays private when the repository is published. Only
the region files (`sources_*.yml`, `mapping_*.yml`), the template and the smoke
test are versioned.

Every stage is resumable, so an interrupted run is not a lost run. `plan`
touches no raster and costs a minute; run it and read it before committing to
hours of downloads.


---


## Appendix B — Case study: Lower Saxony

*One region, worked through in full. It is here as an appendix rather than as a
chapter because nothing in it generalises: every product name, tile size and
attribute below is Lower Saxony's, and the pipeline knows none of them. What
does generalise is the shape of the questions — how many epochs does the
orthophoto have, are the terrain and surface models co-registered, does the
cadastre carry a building function — and those are worth reading before
answering them for a new region. `PORTING.md` asks them as a checklist.*

All products are reachable without registration and carry the same licence
(`dl-de/by-2-0`, attribution *"© LGLN \<year\> — dl-de/by-2-0"*). Each has a
machine-readable download index — a GeoJSON `FeatureCollection` with one feature
per tile — under
`https://arcgis-geojson.s3.eu-de.cloud-object-storage.appdomain.cloud/<product>/`.

| Product | Tiles | Resolution | Format | Epochs | Role here |
|---|---|---|---|---|---|
| **DGM1** | 1 × 1 km, 49 708 statewide | 1 m | COG Float32, nodata −9999 | 71 177 features ⇒ ⌀ 1.4/tile | terrain (`zt`) |
| **DOM1** | 1 × 1 km, 49 708 | 1 m | COG Float32 | 70 807 | nDSM |
| **bDOM20** | 1 × 1 km, 52 611 | 0.2 m | COG Float32, nodata NaN | 128 337, 2021–2026 | newer alternative DSM |
| **DOP20** | 2 × 2 km, 26 763 | 0.2 m | COG, **4 bands incl. NIR** | 93 952 ⇒ ⌀ 3.5/tile, 2012–2026 | NDVI, phenology |
| **LoD2** | 1 × 1 km | — | zipped 3D Shapefile / CityGML | 2024 | building heights |
| **ALKIS-Vektor** | per Landkreis (45 units, 19–330 MB) | — | zipped GeoPackage | 2026 | land use, footprints, bridges |
| **Basis-DLM** | statewide, 536 MB | — | zipped Shapefile | 2026 | roads → `street_type` |

Three properties of this catalogue shape the whole design.

### B.1 DGM1 and DOM1 are perfectly co-registered

All 49 708 tiles share both the tile identifier **and** the acquisition date
between the two products. Selecting the DOM1 tile with the same date as the
chosen DGM1 tile therefore yields a **temporally consistent** nDSM: the surface
and the ground come from one and the same airborne laser scan. bDOM20 is finer
and newer but has no such pairing, so `nDSM = bDOM20 − DGM1` mixes epochs — it is
offered (`sources.dsm: bdom20`) for areas where buildings changed after the last
ALS campaign, at the cost of that consistency.

### B.2 The orthophoto trade-off is real and quantifiable

Month histogram over all 93 952 DOP20 features:

```
Mar 29051  ████████████████████████████████████████████
Apr 16096  ████████████████████████
May 16911  █████████████████████████
Sep  9271  ██████████████
Jul  8535  █████████████
Jun  6989  ██████████
Oct  3521  █████
Aug  2834  ████
```

**The most recent orthophoto of a tile is more often than not a March or April
leaf-off image.** Naively taking the newest epoch would make NDVI collapse under
every deciduous tree. Epoch selection is therefore scored, not sorted:

```
score = w_season · 1[month ∈ months] + w_recency · normalised(year)
```

with `w_season` deliberately dominating `w_recency`. Because the average tile has
3.5 epochs, a leaf-on image usually exists; where none does, the selector relaxes
the season constraint, then the year floor, and **reports it** rather than
failing silently.

The same property is turned into an asset: requesting a second, *leaf-off* epoch
gives free phenology. Tall vegetation that is green in both epochs is coniferous,
green only in summer is deciduous, and in between is a mixed stand — exactly the
distinction PALM's `vegetation_type` draws between `evergreen_needleleaf_trees`,
`deciduous_broadleaf_trees` and `mixed_forest_woodland`.

Ortho tiles are 2 km while work tiles are 1 km, so a per-tile epoch choice
imprints the acquisition lattice on the NDVI field as visible 2 km steps.
`sources.prefer_uniform_epoch` first searches for a single acquisition campaign
(a date cluster no wider than `max_epoch_spread_days`) covering the whole area of
interest, and only falls back to per-tile choice if none exists.

### B.3 LoD2 geometry is not a footprint

A LoD2 tile holds one row per roof *and wall face*: 4 160 rows for 1 234 distinct
`gml_id` values in the sample tile, as `MULTIPOLYGON Z`. Projected to 2D, wall
faces degenerate to lines. LoD2 is therefore used for **attributes only**
(`measHeight`, `function`), joined onto the clean 2D ALKIS footprints.

The join key is hidden in `externRef`:
`http://…/art.htm#_9100$$$DENIAL870000ok0W` — the part after `$$$` is the ALKIS
`objid`. Where `externRef` is absent, the `gml_id` differs from the cadastre id
only by its namespace prefix (`DENILD…` vs `DENIAL…`), which the code uses as a
fallback.

The index's own `tile_id` field is **malformed** for LoD2 (`"25985897"` where
`"325985897"` is meant), so tiles are located by bounding box, never by
identifier.


---


## Appendix C — Porting to another region

See [`PORTING.md`](PORTING.md). In short: write a new `sources_<region>.yml` and
`mapping_<region>.yml`. No Python changes are needed as long as the new region
offers, in a projected CRS:

* a tiled terrain model and a tiled surface model (roles `dtm`, `dsm`),
* a tiled orthophoto **with a NIR band** (role `ortho`),
* a land-use polygon layer with a classification attribute (role `cadastre`),
* optionally building footprints, 3D building models and road centrelines.

The two hard requirements are the **NIR band** — without it there is no NDVI and
the whole classification chain collapses to land use alone — and a **projected**
CRS matching `grid.epsg`, since reprojecting every tile would be both lossy and
slow. palmgeo refuses to start if the sources declare a different CRS.

If you are an agent asked to port: this is a substantial piece of work
(cataloguing the services, mapping the land-use vocabulary, verifying attribute
names against a real delivery), not a config tweak. Say so before starting, and
plan to run `report` repeatedly to find the vocabulary you have not mapped yet.

