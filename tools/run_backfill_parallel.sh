#!/bin/zsh
# Add orientation / parallel-neighbour counts to match tables saved before the census
# recorded them: census days, then warm-season null days. Under the thermal guard.
#   tools/run_backfill_parallel.sh [WORKERS]
cd "$(dirname "$0")/.." || exit 1
PY=.venv/bin/python
[[ -x $PY ]] || PY=../../../.venv/bin/python
for ev in data/events data/events_null; do
  echo "=== $(date '+%F %T') backfill $ev"
  $PY -W ignore tools/thermal/guarded_run.py -- $PY -W ignore scripts/backfill_parallel.py \
      --events $ev --workers "${1:-2}"
done
echo "=== $(date '+%F %T') backfill done"
