#!/bin/zsh
# Per-streak descriptions and seeding-layer exposure for a list of radar-days.
#   tools/describe_days.sh "KMPX 2025-02-15" "KIND 2024-11-28" ...
# Days that already have both outputs are skipped.
cd "$(dirname "$0")/.." || exit 1
PY=.venv/bin/python
[[ -x $PY ]] || PY=../../../.venv/bin/python
for sd in "$@"; do
  set -- ${=sd}
  site=$1; day=$2
  if [[ -f data/events/$site/${day}_described.csv && -f data/events/$site/${day}_exposure.csv ]]; then
    echo "$site $day: already described"; continue
  fi
  $PY -W ignore scripts/describe_day.py $site $day 2>&1 | grep -v "Unknown message\|unexpected" | tail -1
  $PY -W ignore scripts/seeding_efficiency.py $site $day 2>&1 | grep -v "Unknown message\|unexpected" | tail -1
done
echo DESCRIBE_DONE
