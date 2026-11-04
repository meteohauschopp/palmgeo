# Porting palmgeo to another region

The processing code never names a product, an attribute or a code value. It
addresses data by **role** and resolves vocabulary through the mapping file, so
a new region -- a state, a country, a city -- means writing two YAML files:

```
repository/config/sources_<region>.yml     endpoints, formats, field names, CRS
repository/config/mapping_<region>.yml     attribute values -> PALM types
```

and pointing a run config at them:

```yaml
project:
  sources: repository/config/sources_<region>.yml
  mapping: repository/config/mapping_<region>.yml
```

The code has no default region to fall back on: no default CRS, no default
column name, no default tree species. Whatever the region files leave out is
an error that names the missing key, not a silent Lower Saxony value. Start
from the nearest existing pair -- `sources_ni.yml` for a region with a
national or state geodata portal, `sources_global.yml` and `sources_eg.yml`
for one without.

---

## 1. What the new region must offer

| Role | Needed for | Hard requirement |
|---|---|---|
| `dtm` | `zt`, and the base of the nDSM | tiled raster, projected CRS |
| `dsm` | nDSM → vegetation height, building height fallback | tiled raster on the same tiling |
| `ortho` | NDVI, phenology, the whole land-cover half of the classification | **must have a NIR band** |
| `cadastre` | land use, building footprints, bridges | polygons with a classification attribute |
| `buildings3d` | building heights (optional but strongly preferred) | per-building height attribute |
| `topo` | `street_type`, `street_crossings` (optional) | road centrelines with a class attribute |
| `water` | water bodies where the cadastre carries none (optional) | polygons |

`palmgeo/roles.py` is the authoritative version of this table — it carries, per
role, which driver variables depend on it and what the driver looks like when
it is absent, and `palmgeo sources` prints it. Add a role there, not here.

### Every source must state a licence

This is not optional and not a formality. `sources_*.yml` gives each source a
`licence:` naming one or more ids from `repository/config/licences.yml`, plus a
`licence_holder:`. The `sources` stage combines them into the licence position
of the finished driver, and the combination rule is that any single input's
obligation governs the whole thing — one share-alike source makes the driver
share-alike.

```yaml
sources:
  my_terrain:
    role: dtm
    licence: dl-de-by-2.0             # an id from licences.yml
    licence_holder: National Mapping Agency
```

A source that is fused from several upstreams names several, and then
`licence_holder` becomes a mapping so that each holder is credited under their
own licence rather than under somebody else's:

```yaml
    licence: [us-municipal-open-data, ODbL-1.0]
    licence_holder:
      us-municipal-open-data: City Department of Planning
      ODbL-1.0: OpenStreetMap contributors
```

A licence not yet in the registry is added there once, with its URL, its
attribution template, its four obligation flags and — if it has one — the grey
area it is known to have. An id the registry does not define is a hard error at
load time, deliberately: a typo that resolved to nothing would be reported as
"no licence stated" and send the reader looking in the wrong place.

Omitting the field is allowed but resolves to `unknown`, which is modelled as
the most restrictive position there is. That is the correct default for a
source nobody has examined.

### Never write a local path into a region file

A region file describes a *region*, and is meant to serve every run in it. A
staged source works against that, because the index it points at is produced
per run and lives in that run's working directory. Writing the path out ties
`sources_eg.yml` to one Egyptian project, so the next Egyptian domain has to
edit a file that is supposed to be shared — and if the path is absolute it also
ties the file to one machine and carries an account name into a published
repository.

Use a placeholder instead. They are expanded once, against the run being
loaded:

| placeholder | expands to |
|---|---|
| `{workdir}` | the run's working directory, where its staged data lives |
| `{outdir}` | the run's project directory |
| `{project_root}` | the directory holding `repository/` |

```yaml
  fabdem:
    role: dtm
    kind: tile_index
    index_url: "{workdir}/staging/indices/dtm_eg.json"
```

The quotes matter: a YAML value starting with `{` is otherwise read as a flow
mapping. Aswan and Luxor then read the same region file and each finds its own
index.

palmgeo warns at load time when a versioned region file names an absolute local
path, because that is the mistake that survives longest — it works perfectly on
the machine it was written on, and only fails for the next person.

**Inside** the index, `staged_path` is written relative to the project root by
`ports/staging_index.py`, and the cache resolves it against that root. An
absolute path there is not published, but it does not survive the project being
moved or handed on, and it fails hours later inside a worker thread rather than
at the moment it becomes wrong.

### Declare which source fills each role

A region file ends with a `defaults:` block mapping role to source name. It is
only needed where several sources compete for one role, but writing it out is
cheaper than discovering later that the choice was made by dictionary order.

```yaml
defaults:
  dtm: my_terrain
  dsm: my_surface
  cadastre: my_cadastre
```

A run configuration overrides any of these with `sources.<role>:`, and
`palmgeo sources --set <role>=<name>` writes that override for you.

---

## 1a. Grid and coordinate reference

`grid.epsg` must equal the native CRS of the raster products. The loader
refuses a mismatch rather than reprojecting every tile, which would be both
lossy and slow. If your sources genuinely come in different CRSs, reproject
them into a staging area first and point `sources_*.yml` at that.

**Without a NIR band** the land-use × land-cover fusion collapses to plain land
use. The pipeline still runs — every pixel takes the `high_ndvi` or `low_ndvi`
branch according to each rule's `prefer_high_without_ndvi` heuristic — but the
1 m detail inside settlement polygons is lost. If only RGB orthophotos exist,
consider substituting a Sentinel-2 NDVI composite as the `ortho` source at
10 m; it is coarser but far better than nothing.

---

## 2. Source kinds

`sources_<region>.yml` describes each source with a `kind`:

### `tile_index`

A GeoJSON `FeatureCollection`, one feature per tile, whose properties carry the
download URL. Many mapping agencies publish exactly this; in Germany the LGLN
does, and NRW, Berlin, Brandenburg and others in near-identical form. The
Lower Saxony entry below is the worked example.

```yaml
dgm1:
  role: dtm
  kind: tile_index
  index_url: https://…/index.geojson
  index_crs: EPSG:25832        # CRS of the index geometries
  url_field: dgm1              # property holding the download URL
  url_field_fallback: null     # optional second choice
  date_field: Aktualitaet
  tile_field: tile_id          # or null if unusable
  tile_from_filename: 'dgm1_(?P<zone>\d{2})_(?P<east>\d{3})_(?P<north>\d{4})_'
  tile_size: 1000
  resolution: 1.0
  nodata: -9999.0
  format: cog
```

Tiles are located by **bounding-box intersection**, never by identifier, so a
broken or absent `tile_field` is harmless — set it to `null`. Index geometries
in a different CRS are reprojected once, at index build time.

`index_crs` is required: there is no CRS an index could safely be assumed to
be in (plain GeoJSON is EPSG:4326 by specification). `date_field` values are
read as ISO dates (`YYYY-MM-DD`); anything else needs `date_format` in
`strptime` syntax, e.g. `date_format: "%d.%m.%Y"`. An orthophoto source must
give its band order, `bands: {red: 1, green: 2, blue: 3, nir: 4}`.

The CRS of the whole region goes into `crs:` at the top of the file, and has
to be **projected and in metres**: every distance palmgeo computes is in the
units of that CRS. A run whose `grid.epsg` is unset takes it from there.

### `region_index`

The same structure, but one feature per administrative unit and the URL points
at a whole-unit archive. Only units intersecting the area of interest are
fetched, clipped and deleted.

### `single_file`

One statewide archive. Downloaded once, clipped, deleted. Give `size_hint` so
the user is told what is about to happen.

### What the files have to look like

The pipeline reads a fixed set of delivery formats, the ones Lower Saxony
publishes: a cadastre as a **zip holding a GeoPackage** with the layers named
under `layers:`, a 3D building model and a road network as **zipped
shapefiles**, rasters as **GeoTIFF** (ideally cloud-optimised). A region that
publishes something else -- a WFS, a STAC catalogue, separate land-use and
building datasets -- needs a staging step that converts into these, and an
index pointing at the result. That is what `repository/ports/` holds; its
README describes the index a port writes.

### Adding a new kind

If your region publishes through an OGC API, a WFS or a STAC catalogue, the
adapter belongs in `palmgeo/sources.py` next to `load_index`. It must yield
`Item(bounds, url, date, key, extra)`; everything downstream is unchanged.

---

## 3. Field names and case

Every column palmgeo reads is named in the region file; none has a built-in
default. The ones a pipeline stage needs, and fails without:

| where | key | what it is |
|---|---|---|
| `layers.landuse` (cadastre) | `key_field` | the land-use code the mapping keys on |
| `layers.buildings` (cadastre) | `id_field`, `function_field` | building id (joins the 3D model), building function (`building_type`); `null` where the source has none |
| `layers.structures` (cadastre, optional layer) | `type_field` | bridge type |
| `layers.roads` (topo) | `name`, and `class_field` unless the mapping's `street_type.field` names it | road layer and class; `width_field` is optional |
| `fields` (buildings3d) | `building_key`, `height` | model building id and height; `extern_ref`, `extern_ref_split`, `function`, `function_split`, `roof_area` are optional |

Where the 3D model's reference to its cadastre footprint is missing or not a
cadastre id, `fields.cadastre_id_prefix` says what a valid one starts with and
`fields.id_rewrite: {pattern, repl}` derives one from the model's own id --
Lower Saxony needs both, because its LoD2 and ALKIS ids differ only by a
namespace prefix.

The mapping declares `single_trees`: which palm_csd tree type (an index into
`palm_csd/data/tree_defaults.csv`) and species label a detected conifer and
broadleaf crown become. It is required wherever single trees are detected.

Attribute spellings vary between the specification and the actual delivery —
the Lower Saxony Basis-DLM ships `WDM` where the AAA specification writes `wdm`.
All field lookups go through `palmgeo.steps.vectors.resolve_field`, which
matches case-insensitively, so either spelling works. A declared field that
cannot be resolved produces a warning naming the columns that *do* exist,
rather than silently classifying everything as the default.

---

## 4. Writing the mapping

Each land-use class declares its outcome as a function of observed cover:

```yaml
landuse:
  "<key>":
    water:     water:<type>          # unconditional; omit for land classes
    high_ndvi: vegetation:<type>     # NDVI >= threshold
    low_ndvi:  pavement:<type>       # NDVI <  threshold
    high_veg:  vegetation:<type>     # vegetated AND nDSM >= 2 m
    high_veg_fixed: true             # never overridden by phenology
```

Lookup is **longest-prefix**, so a hierarchical code system needs no
enumeration of every leaf: define `43002` for forest and add `43002 VEG 1100`
only where the sub-class carries real information. A flat vocabulary works too
— the prefix rule degenerates to exact matching.

Type names and their valid ranges come from `palmgeo/palm_types.py` and are
checked when the config loads, so a typo fails in the first second.

### The empirical loop

Two commands close the gap between assumption and data:

```bash
CFG=static_driver/<run>/config/<run>.yml
python -m palmgeo.cli vectors   $CFG
python -m palmgeo.cli report    $CFG     # unmapped vocabulary
python -m palmgeo.cli rasters   $CFG
python -m palmgeo.cli calibrate $CFG     # NDVI histogram
```

`report` lists every land-use key and every road class present in the area of
interest that found no explicit rule. Work that list down to nothing before
trusting the result. `calibrate` prints the NDVI histogram of the actual
imagery, which is the only defensible way to choose `ndvi_threshold` for an
uncalibrated 8-bit orthophoto campaign.

---

## 5. Region-specific pitfalls worth checking

Each of these cost real debugging time in the Lower Saxony implementation and
are likely to have analogues elsewhere:

1. **Are the DTM and DSM from the same flight?** In Lower Saxony DGM1 and DOM1
   share tile identifiers *and* dates for all 49 708 tiles, so pairing on the
   date gives a temporally consistent nDSM (`pair_with` / `pair_on` in the
   source definition). Where no such pairing exists, buildings demolished or
   erected between the two epochs appear as artefacts.
2. **When were the orthophotos flown?** Check the month histogram of the index
   before trusting recency. In Lower Saxony 48 % of all epochs are March/April
   leaf-off. If your region flies only in leaf-off conditions, NDVI-based
   vegetation detection is not viable and you need a different cover source.
3. **Is the 3D building model's geometry a footprint?** LoD2 in the AdV
   standard is a set of roof and wall faces, not a footprint, and several rows
   share one building id. Use it for attributes and join onto cadastre
   footprints.
4. **Do the tile grids of different products align?** Ortho tiles are 2 km in
   Lower Saxony while terrain tiles are 1 km. The work tile is 1 km so that one
   work tile always falls inside exactly one ortho tile, avoiding seams *within*
   a tile. Adjust `TILE_SIZE` in `palmgeo/grid.py` if your region's grids differ.
5. **Are index identifiers trustworthy?** The Lower Saxony LoD2 index ships a
   truncated `tile_id`. Prefer geometry.

---

## 5b. What the first non-German port actually cost

Everything below was found the hard way while porting to New York City, for a
run over Manhattan and both river banks (`sources_us_ny.yml`,
`mapping_us_ny.yml`, `repository/ports/us_ny/`). None of it is NYC-specific;
all of it is the kind of thing that silently produces a driver that looks fine
and is wrong.

### Units are checked in two places, and only one of them is the CRS

`grid.epsg` matching the source CRS is necessary, not sufficient.

* **Horizontal.** NYC publishes in US survey feet (EPSG:2263). palm_csd is
  metric, so everything has to be reprojected into a metric CRS before palmgeo
  sees it. That much §1 already says.
* **Vertical.** `gdalwarp` reprojects *coordinates*; it never touches pixel
  *values*. A foot-based DEM warped into UTM still contains feet. Terrain came
  out 3.28x too high, and nothing in the pipeline noticed -- the numbers were
  plausible, just wrong. `repository/tools/ft_to_m.py` does the conversion;
  the check that caught it was comparing the max against a known landmark
  (Manhattan's highest natural point, 80.8 m).
* **In your own arithmetic.** A decimation factor computed as
  `round(1.0 / src.transform.a)` silently became 1 instead of 3, because
  `transform.a` is in the source CRS's units -- feet. The resulting read
  allocated 25 GB, the kernel OOM-killer took the process, and because it ran
  inside the editor's systemd scope the whole editor went down with it. Four
  times, before anyone looked at `journalctl -k`. **If a long job dies without
  a traceback, check the kernel log before theorising.**

### Two elevation sources will not share a vertical datum

Filling one survey's gap with another's data (see the next pitfall) fuses two
*vertical* reference frames, and nothing warns you. The 3DEP tiles staged for
the New Jersey side read ~31.7 m lower than New York City's 2017 product over
the same ground: NYC publishes NAVD88 orthometric heights, the 3DEP tiles are
ellipsoidal, and 31.7 m is simply the geoid separation over New York.

Merged untouched this puts Hoboken at -32 m and turns the Hudson's west bank
into a 32 m cliff -- while every structural check still passes, because the
raster is perfectly well-formed. It is the vertical twin of the feet/metres
bug above, and it is harder to see: the numbers are not implausible, they are
merely in the wrong frame.

**Measure the offset, do not assume it.** Search the second source's tiles
deliberately wider than the gap so the two models genuinely overlap on land,
then read the shift off the overlap:

    python repository/ports/us_ny/stage_nj_supplement.py datum-offset \
        --workdir static_driver/<run>/work --bbox ... --nj-bbox ...
    # median +31.65 m over 103,267 px, IQR 0.31 m

A *tight* spread around a *large* median is the signature of a datum shift; a
wide spread means the horizontal alignment is wrong and no vertical constant
will fix it. `ports/us_ny/stage_nj_supplement.py` refuses to merge if the overlap is too
small to measure or the IQR too wide.

Take the offset from the **DTM** pair and apply the same value to the DSM.
It is a property of the datum, not of the surface, the DSM's own overlap is
contaminated by real change between epochs (IQR 1.09 m against the DTM's
0.31 m), and one shared shift leaves nDSM = DSM - DTM untouched, so building
heights stay correct.

### A municipal survey stops at the municipal border

NYC's LiDAR ends at the city line -- which runs down the middle of the Hudson.
A domain drawn around Manhattan had 38 % of its terrain missing: the river's
west half plus New Jersey. Two different remedies, and they are not
interchangeable:

* **Water** has a physically correct fill: a hydro-flattened surface at local
  sea level (`repository/tools/fill_water.py`). This closed 14 of the 38 points.
* **Land** does not. The remaining 24 % was New Jersey, and nearest-neighbour
  filling it would have smeared Manhattan's terrain across the river. Either
  bring in a source that covers it (`sources_global.yml`) or move the domain.

Check the survey's actual footprint against your intended domain *before*
sizing the grid, not after `domains` aborts.

### A rectangle is the wrong shape for "the bit the other source is missing"

The New Jersey OSM extract was clipped to a strip ending at -73.985, on the
reasoning that "east of that the New York sources are authoritative". True of
the land -- but the state line runs down the middle of the Hudson, and the
Hudson bends east as it goes north: the New Jersey shore is near -74.02 at
40.77 and near -73.97 at 40.85. A rectangle cut at -73.985 therefore sliced
off New Jersey itself above ~40.82. Fort Lee and Englewood Cliffs got no
buildings from OSM, and none from NYC Building Footprints either, because
those stop at the city line. The result was a 5 km strip of land in the driver
with no buildings at all -- reported by the user, not by any check.

Two rules:

* **Clip the gap-filling source to the whole staged AOI, not to a hand-drawn
  strip.** Where the two sources genuinely overlap, decide by precedence at
  mosaic/fusion time, not by geometry at fetch time. Regional OSM extracts are
  already cut to the administrative boundary, so a wide clip cannot duplicate
  anything on the other side of it.
* **Check coverage as a function of position, not as a total.** 56 478 New
  Jersey buildings sounded like plenty. Binned by row, building cover was 26 %
  up to y = 4 517 000 and 0-5 % north of y = 4 519 000 -- a step function, which
  is what a clipped fetch looks like and what no aggregate number reveals.
  `check_ortho.py` does exactly this for imagery (its emptiest-cell test); the
  same idea applied to vector coverage is what caught this.

### Do not scope the cadastre by administrative unit

Fetching MapPLUTO with `where=Borough='MN'` and filtering buildings on the
Manhattan BBL prefix produced a driver with an empty east bank -- every
building in Brooklyn and Queens inside the domain was simply absent. Query the
cadastre **spatially, by AOI**. The administrative unit the site sits in is
almost never the same shape as the domain.

### Tile grids: the work tile must fit inside one source tile

§5.4 states the rule; here is the failure mode. NYSDOP orthophotos are 2500 ft
(762 m) while palmgeo's work tile is 1 km, and `steps/tiles.py` binds exactly
one ortho item per work tile. Every work tile therefore came out ~80 % black,
and its NDVI was noise. The fix used here was to mosaic the tiles into one
raster during staging (as the DEM/DSM already are) rather than to change
`TILE_SIZE` -- mosaicking at 0.5 m for a 2 m work tier keeps ample
oversampling and leaves palmgeo untouched.

### A 4-band image will lose its NIR band to alpha semantics

GDAL labels the 4th band of a 4-band byte image `ColorInterp=Alpha`.
`gdalwarp` then treats it as transparency and writes a synthetic opaque band
(255 everywhere) instead of the NIR data. NDVI silently becomes an inverted
red-brightness map. It *looks* right -- 97 % "vegetation" over a park -- and
classifies asphalt as vegetation too. Use `gdalwarp -nosrcalpha`, and verify
by comparing NDVI over a known built-up area, not only over a known park.

### Calibrate the NDVI threshold against references, not a random sample

`palmgeo calibrate` samples tiles at random. Over an AOI that is a fifth open
water, that sample is dominated by water and reports a median NDVI around
-0.5, which tells you nothing about where to put the threshold. Calibrate
against two known populations instead -- a park polygon and a dense built-up
polygon from the cadastre you already have:

    Central Park   p25 +0.004  median +0.101  p75 +0.204
    dense built-up p25 -0.111  median -0.077  p75 -0.042

The inherited Lower Saxony threshold of 0.12 put only 43 % of the park in
vegetation. Note also the absolute level: a leaf-on summer canopy sits at
+0.5..+0.7, so a park median of +0.10 is itself evidence that the imagery is
off-season -- pitfall 2 above, worth confirming before trusting any threshold.

### Structural verification is not quality verification

`verify_driver.py` passed on a driver in which **half of Central Park was
asphalt**, because "asphalt everywhere" is structurally immaculate: 2-smooth
grid, cell-centred coordinates, exactly one surface type per cell, no terrain
holes. Run `repository/tools/check_classification.py` as well. It
cross-tabulates each land-use class against the surface types that actually
came out, judging each class only on the ground its rule governs (buildings
and water are hard overrides burnt in afterwards), and fails a class whose
rule can only produce vegetation but whose ground came out pavement. It needs
no new declarations -- the expectation is already in `mapping_<region>.yml`.

    python repository/tools/check_classification.py static_driver/<run>/config/<run>.yml

### `gdalbuildvrt` only ever warns when it drops an input

Filling a gap with a second source (here: national NAIP where the state
programme NYSDOP stops at the state line) fuses two products that were never
meant to meet -- and the tool that fuses them drops inputs for several
unrelated reasons, each time with a warning and a zero exit code. Three were
hit for real on one afternoon:

* **Different CRS.** It adopts the first input's CRS and merely *warns* while
  skipping every input that disagrees. NAIP is `EPSG:26918`
  ("NAD83 / UTM 18N"), the staged tiles `EPSG:6347` ("NAD83(2011) / UTM 18N").
  Listing NAIP first dropped **all 605** state tiles -- a mosaic with no
  imagery over the actual city, built without an error. Reproject the odd one
  out into the majority CRS *before* the VRT, then grep the build's stderr for
  `Skipping` and abort on a hit.
* **The comparison is on the CRS *description*, not the EPSG code.** Tiles
  carrying a GeoTIFF-key CRS definition do not byte-match the EPSG registry's
  for the same code, so `-t_srs EPSG:6347` can still be "heterogeneous". Take
  the target with `gdalsrsinfo -o wkt1 <a real tile>` and warp to that file.

* **Different band colour interpretation.** "expected Undefined, got Alpha".
  The staged tiles label band 4 (NIR) as Alpha; a NAIP composite warped with
  `-nosrcalpha` has it Undefined. Same outcome: 605 tiles dropped, exit 0.
  `repository/tools/fix_colorinterp.py` relabels band 4 in place with rasterio
  (`gdal_edit.py -colorinterp_4 undefined` is the documented fix but fails to
  import in some conda environments).
* **Unreadable inputs.** Two tiles were 0 bytes from an earlier interrupted
  run and simply vanished from the mosaic. Check for empty files before
  building, not after.

Two counting traps when you try to verify the result:

* A VRT holds **one `<SourceFilename>` per band per input**, so grepping the
  tag reports 4x the input count for 4-band imagery. Compare *distinct* paths.
* `cat a b` where `a` lacks a trailing newline glues `a`'s last line to `b`'s
  first, quietly losing one input from each. Write generated file lists with a
  trailing newline.

`repository/tools/check_ortho.py` checks a staged mosaic for all of this before
anything is classified: NIR is real data and not a constant alpha plane, no
block of the extent is empty (a dropped source leaves a clean rectangle), and
NDVI actually separates green from built-up.

    python repository/tools/check_ortho.py --ortho static_driver/<run>/work/staging/ortho/<m>.tif \
        --landuse static_driver/<run>/work/vector/cadastre_landuse.gpkg \
        --veg-query "LandUse == '09'" --built-query "LandUse in ['03','04','05']"

### A staged index must describe what is staged

The single-feature indices for staged rasters had their geometry derived from
the WGS84 bbox the data was *ordered* for, which reprojected to a box 763 m
wider than the rasters actually written. palmgeo then claimed coverage for
work tiles that were pure nodata, and reported them downstream as missing
tiles -- with nothing pointing at the index as the thing that lied. Derive the
index geometry from the same constant the warp's `-te` uses.

### Do not trust a filename or a catalogue bbox to place a tile

In `3dep-lidar-dtm` / `3dep-lidar-dsm` the STAC bboxes are wrong (a search over
Hoboken returned a tile georeferenced 40 km south) **and** the tile *ids* are
offset differently between the two collections: the tiles covering this AOI are
`6-8..6-11` in dtm but `6-3..6-5` in dsm, ~60 km apart under otherwise
identical naming. Search deliberately wide, then open every candidate and
select on its real bounds. Warping several such tiles in one `gdalwarp` call
also produced an empty raster although each warped correctly alone -- build a
VRT and warp that once instead, and assert the result is non-empty.

## 6. Data that is not open geodata

Anything you produce yourself enters through `extra_inputs` and needs no code
change:

```yaml
extra_inputs:
  trees:                              # -> palm_csd `trees` + `columns`
    path: data/tree_register.gpkg
    columns: {height: tree_height, crown: tree_crown_diameter}
  surfaces:                           # burnt in after the cadastre classification
    - {path: data/green_roofs.gpkg, type: vegetation:short_grass}
  rasters: {}                         # role -> raster, replaces a derived product
```

`surfaces` entries are applied in order, after land use and cover but before
water and buildings, so they override the automatic classification without
overriding the cadastre's hard facts.
