#!/bin/zsh
# A larger warm-season null: every day of July and August 2025 and of July 2024 at the eight
# radars (the first null used every third day of July 2025, 80 tested radar-days, too few to
# check calibration at the 1% level). Waits for the census queue to finish first.
#   tools/run_null_extension.sh [WORKERS]
cd "$(dirname "$0")/.." || exit 1
W=${1:-4}
while pgrep -f "run_census_queue.sh" > /dev/null; do sleep 120; done
for range in "2025-07-01 2025-08-31" "2024-07-01 2024-07-31"; do
  echo "=== $(date '+%F %T') null census $range"
  tools/run_null_census.sh ${=range} $W 1
done
echo "=== $(date '+%F %T') null extension done"
