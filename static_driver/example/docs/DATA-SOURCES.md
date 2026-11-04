# Data sources — example_static

Region `NI`, EPSG:25832. Established 2026-11-04 by `palmgeo sources`, approval `67b187b787952b76`.

This file records what this driver is built from and what that means for its use. It is generated; edit the run configuration, not this file.

## Sources

| Role | | Source | Coverage of AOI | Licence | Becomes |
|---|---|---|---|---|---|
| `dtm` | required | `dgm1` | 100.0 % | dl-de-by-2.0 | zt, +vegetation_height, buildings_2d |
| `dsm` | required | `dom1` | 100.0 % | dl-de-by-2.0 | vegetation_height, +buildings_2d |
| | *alternative* | `bdom20` | 100.0 % | dl-de-by-2.0 | |
| `ortho` | optional | `dop20` | 100.0 % | dl-de-by-2.0 | +vegetation_type, pavement_type, soil_type, trees |
| `cadastre` | required | `alkis` | 100.0 % | dl-de-by-2.0 | building_id, building_type, bridges_2d, bridges_id, +vegetation_type, pavement_type, water_type, soil_type |
| | *alternative, switched off* | `lbm_de` | not measurable | CC-BY-4.0 | |
| | *alternative, switched off* | `hausumringe` | not measurable | dl-de-by-2.0 | |
| `buildings3d` | optional | `lod2` | 100.0 % | dl-de-by-2.0 | +buildings_2d |
| `topo` | optional | `bdlm` | not measurable | dl-de-by-2.0 | street_type, street_crossings |
| `water` | optional | — | — | — | *not available* |

## Currency

How old the data is. Nothing else in the driver records this, and it is what explains a building that is missing rather than wrong: building geometry comes from the surface model, so anything built after that model was flown has no height and drops out.

The date is whatever the source declares as its currency field, and what that means differs per product: for an elevation model it is normally the survey, but for a derived product it can be the publication of a model whose geometry is older -- a 3D building model can carry this year's date and still lack every building its years-old laser scan did not see. Treat the column as an upper bound on how new the content is.

| Role | Source | Acquired | Age |
|---|---|---|---|
| `dtm` | `dgm1` | 2016-04-18 | 11 years |
| `dsm` | `dom1` | 2016-04-18 | 11 years |
| `ortho` | `dop20` | 2013-06-18 … 2025-03-06 | 2 years |
| `cadastre` | `alkis` | 2026-10-30 | under a year |
| `buildings3d` | `lod2` | 2024-06-11 … 2024-06-12 | 2 years |
| `topo` | `bdlm` | no date stated |  |

Oldest to newest: **2013-06-18 to 2026-10-30**, a span of 13 years.


## What is missing, and what that costs

**`water`** — Water bodies, where the cadastre does not carry them itself.

Water comes from the cadastre alone. Where that is incomplete -- typically outside Europe -- rivers and lakes classify as land, which changes the surface energy balance of the whole domain.


## Licence position

### Datenlizenz Deutschland -- Namensnennung 2.0 (`dl-de-by-2.0`)

<https://www.govdata.de/dl-de/by-2-0>

The standard licence of German open government data. Free use including commercially, free redistribution and modification, on the single condition that the source, the holder and the licence are named, and that the notice is kept on any modification.

## What this driver requires of you

- **attribution required** — from alkis (dl-de-by-2.0), bdlm (dl-de-by-2.0), dgm1 (dl-de-by-2.0), dom1 (dl-de-by-2.0), dop20 (dl-de-by-2.0), lod2 (dl-de-by-2.0)

### Attribution (required)

This credit must travel with the driver and with anything published from it:

> (C) LGLN <year> -- dl-de/by-2-0

`<year>` is left for you to fill in, because what it means is the licence holder's call: for most it is the year the data was obtained, for Copernicus the year of acquisition. A wrong year in a credit is worse than an open one.

---

*A summary of licence terms, not legal advice. The linked texts are authoritative and can change.*
