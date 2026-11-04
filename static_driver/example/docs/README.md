# palmgeo example — Hannover city centre

| | |
|---|---|
| **Domains** | `root` 64×64 @ 32 m<br>`N02` 128×128 @ 8 m, parent `root`<br>`N03` 128×128 @ 2 m, parent `N02` |
| **Region** | `NI`, EPSG:25832 |
| **Built** | 2026-09-09; rebuilt 2026-11-04 from the then-current LGLN data, 7.6 MB of drivers, palmgeo 1.2.0 |
| **Data basis** | `dgm1`, `dom1`, `dop20`, `alkis`, `lod2`, `bdlm` |
| **Licences** | dl-de-by-2.0 — see [`DATA-SOURCES.md`](DATA-SOURCES.md) |

**What you may do with it.** Use it freely, including commercially, and
redistribute it, as long as the sources are credited: *© LGLN 2026 —
dl-de/by-2-0*. The run configuration puts that credit into every driver's
`source` attribute (`csd.attributes.source` in `repository/config/example.yml`);
leave it there.

**This is the only project committed to the repository.** Every other run is
git-ignored, because a city-scale driver is hundreds of megabytes. This one is
7.6 MB, so you can open a real PALM static driver without running anything —
and, having run it yourself, compare your result against a known one.

---

## Why this exists twice over

It is the worked example *and* the end-to-end test, deliberately the same
artefact. A test nobody looks at rots quietly; an example nobody runs rots the
same way. Keeping them one thing means every pipeline change is exercised on
what newcomers read, and every claim here is one the test would catch going
wrong.

## What it covers

9 km² of Hannover city centre, in about two minutes on a warm cache. Three
nested domains, one per vegetation mode, each at a resolution where that mode is
the right one:

* **`root`**, 2048 × 2048 m at 32 m, `vegetation_mode: unresolved` — tall
  vegetation stays a rough LSM surface type (z0 = 2 m), which is correct for a
  domain whose first grid level, at 16 m, sits well above the canopy. It is
  also the only domain built from a second, 4 m work tier.
* **`N02`**, 1024 × 1024 m at 8 m, `vegetation_mode: mixed` — canopy patches
  become a 3D leaf-area field on the 8 m grid, and the ground beneath them
  short grass. Not `unresolved`: at 8 m the first grid level is 4 m, and PALM
  refuses to start (`LSM0048`) when z0 is not below half of it.
* **`N03`**, 256 × 256 m at 2 m, `vegetation_mode: resolved` — canopy patches
  become a 3D leaf-area field *and* individual trees are segmented from the
  canopy height model and written into the driver.

That third mode is the point. The crown segmentation runs at any setting, but
only `resolved` writes its result, and an example that skipped it would skip
the part of palmgeo that is not available anywhere else.

## What came out

Measured on the build of 2026-11-04 (LGLN data as of that day; the cadastre's
currency date was 2026-10-30):

| | root (32 m) | N02 (8 m) | N03 (2 m) |
|---|---|---|---|
| vegetation | 22.4 % | 23.6 % | 15.4 % |
| … of it tall vegetation as LSM type | 18.5 % of the domain | — | — |
| pavement | 48.9 % | 45.6 % | 47.3 % |
| buildings | 23.9 % | 25.9 % | 37.3 % |
| water | 4.8 % | 4.9 % | 0.0 % |
| leaf-area-density field | — | 16.5 % of columns, 4 levels up to 28 m | yes, 13 levels up to 25 m |
| individual trees | — | — | **68** |

8485 buildings were kept over the area of interest: 9126 measured from the
surface model, 16 from the 3D building model where the surface model had no
value, and 666 dropped as implausible or below 2.5 m. 160 tree detections
standing on building footprints were discarded — a roof that reads as green is
a green roof, not a tree.

The columns cover different ground — the outer domains hold the lake and the
park, `N03` only the dense centre — so they are not comparable as they stand.
Over `N03`'s 256 m, the 8 m and 2 m grids agree closely: 37.3 % buildings at
2 m against 37.2 % at 8 m, 15.4 % vegetation against 14.6 %, which assigns the
smallest green patches to the pavement around them.

## A worked example of a data problem: the missing houses

*Counted on the build of 2026-10-02, on the cadastre as of 2026-09-02. The
drivers committed now were rebuilt on 2026-11-04 from newer LGLN data, so
single counts may differ; the mechanism does not.*

Open `N03`'s plot and look at its bottom-right corner. A row of houses that
exists on the ground is not in the driver. That is not a palmgeo bug, and
following it through is the most useful thing this example teaches.

**What is there.** The cadastre knows the buildings: 18 footprints intersect
that corner, and their currency date is 2026-09-02. Nine of them reach the
driver; nine do not — seven residential buildings between 323 and 533 m², plus
two small outbuildings.

**Why they drop out.** Building height is *measured*, from the surface model
minus the terrain model. Sample the nDSM inside those footprints:

| | 75th percentile of the nDSM inside the footprint |
|---|---|
| the nine that survived | 9.9, 10.4, 15.6, 15.7, 15.7, 15.8, 16.1, 16.2 m |
| the nine that dropped | 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.2, 2.3 m |

Six of them measure **exactly zero**, and the rest fall below the 2.5 m
`buildings.min_height`. Where the cadastre maps a house, the surface model sees
open ground.

**Why the surface model sees open ground.** It was flown on **2016-04-18** —
see the Currency table in [`DATA-SOURCES.md`](DATA-SOURCES.md). The cadastre is
ten years newer. Anything built in between exists in the register and not in
the elevation data.

**Why the fallback does not save them.** `buildings.height_fallback: lod2`
should supply a height where the surface model has none, and the LoD2 model is
dated 2024. It has heights for two of the nine, both 2.3 m — because LoD2's
*geometry* derives from the same 2016 laser scan, and its 2024 date is when the
model was published, not when the ground was measured. This is why the Currency
table says the date is an upper bound on how new the content is.

**What palmgeo did right.** It dropped them and counted them: the run reports
`666 dropped as implausible or too low` across the area of interest. A driver
with a nine-metre building where the data says bare ground would be worse than
a driver with a hole, and a silent hole would be worse than both.

**What you would do about it.** Nothing, if the study is about the rest of the
area. If those houses matter: supply them through `extra_inputs.surfaces`,
raise `buildings.height_default` to give unmeasured footprints a typological
height (and accept that it is a prior, not a measurement — palmgeo counts those
separately for exactly that reason), or find a newer surface model. What you
cannot do is not know, which is what the Currency table is for.

## Running it yourself

```bash
export PYTHONPATH=repository
CFG=repository/config/example.yml

python -m palmgeo.cli sources $CFG              # what it rests on, and its licence
python -m palmgeo.cli sources $CFG --approve    # accept, unlocking the downloads
python -m palmgeo.cli all     $CFG              # about two minutes

MPLBACKEND=Agg env -u PROJ_DATA -u PROJ_LIB \
    palm_csd --pdf static_driver/example/config/example_static_csd.yml
python repository/tools/finalize.py      $CFG
python repository/tools/verify_driver.py $CFG   # must end in ALL CHECKS PASSED
```

`MPLBACKEND=Agg` is not optional: without it matplotlib segfaults while
resolving its default backend, *after* palm_csd has written a perfectly good
driver — so the crash looks like a driver problem and is not one.

The run configuration lives in `repository/config/example.yml` rather than in
this directory's `config/`, because unlike every other project it is versioned
and shared. What `config/` holds here is what the run produced: the generated
palm_csd configuration, and a record of the region files as they stood.

## What it does not show

* **One region.** Lower Saxony publishes a download index per product, which is
  the easy case. A region without one needs a staging port — see
  [`../../../repository/ports/`](../../../repository/ports/).
* **Scale.** 9 km² at 1 m work resolution is 187 MB of intermediates. A city at
  the same detail is 50–100 GB, and the trade-offs in README part IV are about
  that regime, not this one.
* **A real author.** `csd.attributes` carries a placeholder institution on
  purpose: palmgeo does not invent an affiliation, and a real run replaces it.
