#!/bin/zsh
# Phase A of the census: 8 radars, one winter at a time, under the thermal guard.
#   tools/run_census_phase_a.sh WINTER_START WINTER_END [WORKERS]
# e.g. tools/run_census_phase_a.sh 2024-11-01 2025-03-31 2
# Days whose ADS-B is not downloaded yet are skipped and can be picked up by running the
# same command again later (finished radar-days are not redone).
cd "$(dirname "$0")/.." || exit 1
PY=.venv/bin/python
[[ -x $PY ]] || PY=../../../.venv/bin/python
SITES=KFTG,KLOT,KMPX,KDTX,KCLE,KBUF,KPBZ,KIND
mkdir -p logs
exec $PY -W ignore tools/thermal/guarded_run.py -- $PY -W ignore scripts/census.py "$1" "$2" \
    --workers "${3:-2}" --sites $SITES
