# palmgeo, in plain language

*The short version. [`README.md`](README.md) is the full workflow, and it is
long on purpose; this file is for someone meeting palmgeo for the first time.*

palmgeo was built to be used with an agentic AI: you describe what you want to
simulate, and the agent runs the pipeline and asks you for every decision that
is yours to make — see [Using palmgeo with an AI agent](#using-palmgeo-with-an-ai-agent).

---

## What problem this solves

To simulate the air in a city, [PALM](https://palm.muk.uni-hannover.de/) has to
be told what the city *is*. Not a picture of it — a number for every cell of a
three-dimensional grid: how high the ground is, where the buildings stand and
how tall they are, which patch is asphalt and which is lawn, where the trees
are and how much leaf each one carries. That description is one file, the
**static driver**, and PALM will not start without it.

Everything in it exists as open government data already: laser-scanned terrain
and surface models, aerial photographs with an infrared band, the cadastre, 3D
building models, road networks. What does not exist is the path from that data
to that file. Getting there means finding which of a dozen products covers your
patch of ground, downloading tens of gigabytes, agreeing every dataset to one
coordinate system and one grid, deciding what each cadastre code means in PALM's
vocabulary, measuring building heights, finding individual trees — and getting a
grid geometry right that PALM will reject outright if it is off by one.

**palmgeo is that path, written down and made repeatable.**

---

## Where it sits

```
   open government geodata
   terrain · surface model · aerial imagery · cadastre · 3D buildings · roads
              │
              │   ┌──────────────────────────────────────────────┐
              └──▶│  palmgeo            ← this repository        │
                  │  finds it, fetches it, agrees it onto one    │
                  │  grid, classifies every pixel, measures      │
                  │  every building, detects every tree          │
                  └──────────────────────────────────────────────┘
                       │  GeoTIFFs + a palm_csd configuration
                       ▼
                  ┌──────────────────────────────────────────────┐
                  │  palm_csd           ← part of PALM           │
                  │  turns them into the NetCDF static driver    │
                  └──────────────────────────────────────────────┘
                       │  <run>_static, <run>_static_N02 …
                       ├─────────────────────┐
                       ▼                     ▼
                  ┌──────────────┐    ┌──────────────────────────┐
                  │  PALM        │    │  PALMPaint  (optional)   │
                  │  simulates   │    │  edit a driver by hand:  │
                  └──────────────┘    │  add a planned building, │
                                      │  remove one, repaint a   │
                                      │  surface for a scenario  │
                                      └──────────────────────────┘
```

palmgeo comes **first**, before both of them. `palm_csd` is PALM's own tool and
already knows how to write a valid driver — but only if you hand it a complete,
consistent, correctly gridded set of GeoTIFFs and a configuration. Producing
*that* is the hard part, and it is the part palmgeo does. palmgeo never writes a
NetCDF driver itself and never runs PALM.

**[PALMPaint](https://doi.org/10.5281/zenodo.23192509)** comes last and is for a
different job: changing a driver by hand after it exists — a building that is planned but not yet built, a scenario where
a car park becomes a park. palmgeo gives you reality; PALMPaint lets you edit it.

---

## What it actually does

| Stage | In one sentence |
|---|---|
| `sources` | Shows what geodata actually exists for your area, under which licence you may use it, and which part of the driver each dataset becomes — **before** anything is downloaded. Nothing bulk downloads until you approve it. |
| `plan` | Works out which data tiles your domains need and how many gigabytes that is — **without downloading anything**. |
| `vectors` | Fetches the cadastre, building footprints, bridges and roads, and clips them to your area. |
| `rasters` | Streams terrain, surface model and aerial imagery tile by tile; computes object height and vegetation index; throws each tile away once its numbers are extracted. |
| `buildings` | Measures each building's height from the laser surface model rather than trusting an attribute — with the 3D model as a fallback. |
| `trees` | Finds individual trees: tree tops as peaks in the canopy height, crowns grown outward from them, each one sorted into broadleaf or conifer and given a trunk diameter and a leaf area. |
| `compose` | Decides what every single pixel *is*, by combining what the cadastre says it is *for* with what the imagery says is *growing* on it. A residential block becomes lawn, driveway and hedge. |
| `domains` | Aggregates that onto each PALM domain's own grid — every nest gets its own consistent set. |
| `csd` | Writes the palm_csd configuration that turns all of it into a driver. |

Then `palm_csd` runs, and two checks follow: one that the driver is
structurally valid, and one that the classification is *plausible* — because a
uniformly wrong map is structurally perfect.

## Two things worth knowing

**It is region-agnostic by construction.** The code relies on no data product,
no cadastre code, no column name and no coordinate system of any one place. All
of that lives in two YAML files per region — *where the data is* and *what its
vocabulary means*. The repository ships Lower Saxony as the reference, New York
City and Egypt as worked case studies, and a library of worldwide sources to
start any other region from. Porting to a new state or country means writing
those two files, not touching the code. See [`PORTING.md`](PORTING.md).

Some regions need one step more. A region that publishes a download index per
product can be read directly; one that publishes a single archive per state, or
a satellite API, or nothing but a map server, needs a *port* — a script under
`repository/ports/` that fetches and reprojects and fuses whatever the region
offers into files palmgeo can index. That is a few hundred lines of work per
region, and knowing which case you are in is the first question to answer.

**It tells you what you are allowed to do with the result.** A driver inherits
the licence of every dataset in it, and by the time it exists there is no
visible seam between the terrain that was public domain and the footprints that
were share-alike. The `sources` stage works out the combined position and says
it plainly, before anything is downloaded:

```
[x] non-commercial use only
    from: fabdem (CC-BY-NC-SA-4.0)
[x] share-alike: a distributed derivative carries this licence
    from: fabdem (CC-BY-NC-SA-4.0), osm_cadastre (ODbL-1.0), osm_highways (ODbL-1.0)
```

That is a real driver — an Egyptian one — and a single input decides it. The
rule is that an obligation applies if *any* source imposes it, so one
non-commercial dataset among six makes the whole driver non-commercial. The
stage also refuses to let the download stages start until somebody has accepted
that, which is the one thing in the pipeline that is gated for a reason outside
the software.

**It is resumable.** A city-scale run moves tens of gigabytes over hours. Every
tile records a fingerprint of the settings it was built with, so an interrupted
run continues where it stopped, and changing one threshold rebuilds exactly the
tiles that depended on it — not the other 900.

---

## Using palmgeo with an AI agent

palmgeo is written to be **driven by an agent**. The README's **Appendix A** is
addressed to one directly: what it must ask you before starting, what it should
work out by itself, and what it must never decide alone (moving your domain
origin, discarding data quietly, choosing how many cores you spend).

You do not need to know the stages, the flags or the YAML. You need to know
*what you want simulated*. Point an agent at this repository and describe that.

### An example prompt

> Build me a PALM static driver with palmgeo, in `/path/to/palmgeo`.
>
> **Where:** Hannover, Germany — centred on the Neues Rathaus.
>
> **What I want:** two domains. A root of about 2 km × 2 km at 8 m, and one
> nest of about 500 m × 500 m at 2 m centred on the building. Summer, midday.
>
> **What I care about:** the courtyard and the park to the south. I want
> individual trees resolved in the nest, not a vegetation patch.
>
> **Attribution:** author *<your name>*, institution *<your institution>*,
> contact *<your address>*.
>
> Read `README.md` Appendix A first. Ask me about anything on that list I have
> not told you, in one batch — do not guess. Then run `sources` and show me the
> licence position before approving anything, and `plan` with the tile count
> and transfer estimate. Wait for my go-ahead before downloading. Afterwards,
> run the verification and tell me what was built, what was dropped and what
> you assumed.

**Why that prompt is shaped like that:**

* **It names the region first.** Only regions with a `sources_*.yml` can be
  built at all. Anywhere else needs those two YAML files written first, which is
  a real piece of work — better to learn that in the first minute.
* **It gives extents and resolutions, not a grid.** Point counts, 2-smooth
  rounding and nest origin snapping are arithmetic, and the agent should do
  them and report the result.
* **It says what you care about.** "Individual trees in the nest" is a modelling
  decision with a real cost — it roughly triples the finest tier's runtime. The
  agent cannot infer it from a place name.
* **It withholds the licence approval.** `sources --approve` accepts, on your
  behalf, every obligation the driver will carry. An agent should show you the
  report and let you say yes; the pipeline enforces this by refusing to
  download until the approval exists.
* **It asks for a stop before the download.** `plan` does no I/O and takes a
  minute. A wrong grid found there costs a minute; found after the download, it
  costs a day.
* **It asks what was dropped.** A run that silently discards 3 % of the
  buildings is worse than one that fails, and it will not tell you unless asked.

---

## Licence, and what a driver may be used for

The **code** is free software under the GNU GPL, version 3 or later — see
[`LICENSE`](LICENSE). To cite palmgeo, use [`CITATION.cff`](CITATION.cff) or
the Zenodo DOI of the release you used.

The **data** is not, and this matters more than it sounds. Every source carries
its own licence, and the driver you build inherits them all. Lower Saxony's
products are © LGLN under `dl-de/by-2-0`, which requires attribution and permits
commercial use; other regions and other products differ, and some open datasets
forbid exactly the use you had in mind. `palmgeo sources` works out the combined position for
your particular run and writes it to `static_driver/<run>/docs/DATA-SOURCES.md`
— read that before publishing anything built from the driver.

`palm_csd` and PALM are part of the PALM model system (Leibniz Universität
Hannover, Technische Universität Berlin) and are GPLv3. palmgeo neither includes
nor links them; it writes their input and is run alongside them.

The single-tree detection is ported from the Bachelor's thesis of Joana-Pauline
Rentz (2026, Institute of Meteorology and Climatology, Leibniz University
Hannover); the README's Credits section has the full reference.
