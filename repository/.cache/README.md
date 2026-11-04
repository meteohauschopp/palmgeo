# .cache/

Shared download cache. **Not versioned** — it is a cache.

```
.cache/
├── indices/<product>.json   condensed source download indices
└── <product>/               source tiles being processed, LRU-evicted
```

The indices can be large -- for Lower Saxony ~150 MB of GeoJSON, one feature
per tile statewide: 71 177 DGM1 features, 93 952 DOP20, and so on -- which is
why several runs should point
`project.cachedir` at one directory rather than each fetching their own copy. The
cache is content-addressed by file name and safe to share between concurrent
runs; files being read are pinned against eviction.

`runtime.profiles.<profile>.cache_budget_gb` is a hard cap. Source tiles are
released as soon as they have been consumed, so with `keep_downloads: false` the
cache stays near the index size even over a run that moves 30 GB. Set
`keep_downloads: true` when the machine that computes has no internet and needs
everything prefetched.

```bash
# from the project root, the directory holding repository/
PYTHONPATH=repository python -m palmgeo.cli clean static_driver/<run>/config/<run>.yml cache
```

Deleting it costs a re-download of the indices (a few minutes), nothing else.
