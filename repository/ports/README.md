# repository/ports/

**One directory per region whose data palmgeo cannot fetch by itself.**

Most of palmgeo knows nothing about any region: the pipeline addresses data by
*role*, and a region declares its products in two YAML files. That works as
long as the region publishes something palmgeo can read directly — a download
index per product, listing one tile per feature with its geometry, its date and
its URL.

Plenty of regions do not. They publish one archive for a whole state, or a STAC
API, or a tile server in its own projection, or nothing but an overpass query
and a hope. A **port** absorbs that difference. It fetches what the region
offers, reprojects and mosaics and fuses it into files on disk, and writes an
index in the shape `sources.py` reads. Everything downstream then works
unchanged — which is the point, and the reason the ugliness is confined here.

```
repository/ports/
├── staging_index.py       the contract: writing an index palmgeo can read
├── eg/stage_sources.py    Egypt — global fallback data, no national products
└── us_ny/
    ├── stage_sources.py        New York City and State
    └── stage_nj_supplement.py  New Jersey, where the New York sources stop
```

Both are case studies, written for one run each (over Manhattan, and over
Aswan) and carried through to finished drivers. Neither names its site any
more: every script takes the run's work directory and the area to stage on its
command line (`--workdir`, `--bbox`; see each script's docstring).

`us_ny` has two scripts because the Manhattan domain crossed the Hudson and every
New York source ends at the state line. The supplement fills each role from a
national or global source for the western bank, and merges the two — including
the vertical datum offset between them, which it measures rather than assumes.

The scripts predate `staging_index.py` and still write their indices by hand;
new ports should use it.

## Why these are not in `tools/`

`tools/` holds things you run against any project. A port is the opposite: it
is specific to one region, it runs once, and its output is data rather than an
answer. Keeping them apart means `tools/` stays readable, and it makes the size
of a port visible — roughly 900 to 1300 lines each, which is the honest cost of
a region that publishes badly, and worth knowing before promising one.

## Writing a port

Read [`../../PORTING.md`](../../PORTING.md) first: it covers the two YAML files,
which roles are required, what each becomes, and the licence declarations. A
port is only needed when those two files cannot address the region's data
directly.

Then write a script that produces staged files and calls `write_index`:

```python
from staging_index import StagedFile, reproject_bounds, write_index

write_index(
    indices_dir / "dtm_<region>.json",
    [StagedFile(tile_id="N24E032", path=staged_tif,
                bounds=reproject_bounds(bbox_wgs84, 4326, target_epsg),
                date="2023-04-11")],
)
```

`write_index` refuses to write an entry whose file is not on disk. That check
is the whole reason the module exists: every port written so far has, at least
once, emitted an index naming a file that a later rename made wrong, and the
mistake surfaced hours later inside a worker thread rather than at the moment
it was made.

Point the region file at what you wrote:

```yaml
  my_terrain:
    role: dtm
    kind: tile_index
    index_url: <path to the index>
    url_field: staged_path
    date_field: acquired
```

## A port is a last resort

It is a few hundred lines of code that only one region will ever run, and it
has to be re-run whenever the upstream changes. Before writing one, check
whether the region publishes an index already — `sources_ni.yml` needs no port
at all — and whether the global fallback pair in `sources_global.yml` is good
enough for the resolution you need. It usually is not below about 16 m, and
`PORTING.md` says why.
