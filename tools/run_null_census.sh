#!/bin/zsh
# Warm-season null calibration: the census pipeline on summer days at the Phase A radars, with
# only tracks below 13 000 ft at points the HRRR puts above freezing (no seeding possible).
#   tools/run_null_census.sh START END [WORKERS] [STEP]
# e.g. tools/run_null_census.sh 2025-07-01 2025-07-31 2 3
cd "$(dirname "$0")/.." || exit 1
PY=.venv/bin/python
[[ -x $PY ]] || PY=../../../.venv/bin/python
SITES=KFTG,KLOT,KMPX,KDTX,KCLE,KBUF,KPBZ,KIND
exec $PY -W ignore tools/thermal/guarded_run.py -- $PY -W ignore scripts/census.py "$1" "$2" \
    --workers "${3:-2}" --step "${4:-3}" --null-alt 13000 --sites $SITES
