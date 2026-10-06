#!/bin/zsh
# Run the Phase A census winter by winter, one at a time (the laptop's CPU is the limit).
# Waits for any census run already going, then runs each winter; days whose ADS-B is not
# downloaded yet are skipped and picked up on a later pass (finished days are not redone).
#   tools/run_census_queue.sh [WORKERS]
cd "$(dirname "$0")/.." || exit 1
W=${1:-3}
# wait for the winter 2024-25 run started separately
while pgrep -f "scripts/census.py 2024-11-01" > /dev/null; do sleep 60; done
for range in "2024-11-01 2025-03-31" "2023-11-01 2024-03-31" "2025-11-01 2026-03-31" "2023-03-10 2023-03-31"; do
  echo "=== $(date '+%F %T') census $range"
  tools/run_census_phase_a.sh ${=range} $W
done
echo "=== $(date '+%F %T') queue done"
