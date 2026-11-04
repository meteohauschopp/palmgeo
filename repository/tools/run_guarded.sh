#!/bin/bash
# Run a long palmgeo/GDAL job so that a runaway allocation kills only itself.
#
# WHY: a single arithmetic slip in a staging step -- a factor confusing metres
# with feet in a decimation ratio -- turned a modest read into a ~25 GB
# allocation. The kernel's OOM killer resolved it, but the OOM killer chooses
# its victim by score, not by blame, and on a desktop the process it picks is
# often the interactive session that launched the job rather than the job.
# Losing an editor to a background task is a bad way to learn about a typo.
#
#   ulimit -v         caps address space, so an over-large allocation fails
#                     INSIDE this process, as an ordinary MemoryError with a
#                     traceback pointing at the line that asked for it
#   oom_score_adj     if the system runs out of memory anyway, marks this
#                     process as the one to kill first. Worth setting
#                     explicitly: a job inherits its launcher's score, which
#                     on a desktop session is deliberately low.
#
# Verify the limit took with: grep 'Max address space' /proc/<pid>/limits
#
# usage: run_guarded.sh <limit_gb> <command...>
#   e.g. run_guarded.sh 8 python -m palmgeo.cli rasters <config>
set -euo pipefail
LIMIT_GB="$1"; shift
ulimit -v $((LIMIT_GB * 1024 * 1024))
echo 1000 > /proc/self/oom_score_adj 2>/dev/null || true
unset PROJ_DATA PROJ_LIB   # a stale PROJ_DATA shadows rasterio's own copy
                           # and breaks every CRS operation; see palmgeo/__init__.py
exec nice -n 10 "$@"
