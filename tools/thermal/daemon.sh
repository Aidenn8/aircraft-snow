#!/bin/bash
# Thermal guard for a fanless MacBook Air.
#
# Every 20 s it reads macOS's own thermal pressure (0 nominal, 1 fair, 2 serious,
# 3 critical) and the battery temperature, writes the verdict to ./state and
# appends a line to ./thermal.log. Heavy jobs are started through guarded_run.py,
# which registers their process group in ./pgids/; this daemon then enforces the
# verdict on them directly:
#   PAUSE (serious, or battery >= 40 C): SIGSTOP the job (frozen, no CPU use)
#   resume (fair or better and battery < 38 C): SIGCONT
#   STOP  (critical, or battery >= 45 C): SIGTERM the job and leave state=STOP
# STOP is sticky: nothing restarts until someone deletes ./state.
#
# Run it wrapped in `caffeinate -i -s` so the Mac does not sleep mid-job.

DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$DIR/pgids"
STATE_FILE="$DIR/state"
LOG="$DIR/thermal.log"
PAUSE_C=${PAUSE_C:-40.0}; RESUME_C=${RESUME_C:-38.0}; STOP_C=${STOP_C:-45.0}

verdict="OK"
[ -f "$STATE_FILE" ] && grep -q '^STOP' "$STATE_FILE" && verdict="STOP"

while true; do
  ts=$(date '+%Y-%m-%d %H:%M:%S')
  therm=$("$DIR/thermalstate" 2>/dev/null || echo -1)
  raw=$(ioreg -rn AppleSmartBattery | awk '/"Temperature" =/ {print $3; exit}')
  batt=$(awk -v r="${raw:-0}" 'BEGIN {printf "%.1f", r/100}')
  load=$(sysctl -n vm.loadavg | awk '{print $2}')
  ac=$(ioreg -rn AppleSmartBattery | awk '/"ExternalConnected" =/ {print $3; exit}')

  prev="$verdict"
  if [ "$verdict" != "STOP" ]; then
    if [ "$therm" -ge 3 ] || awk -v b="$batt" -v s="$STOP_C" 'BEGIN{exit !(b>=s)}'; then
      verdict="STOP"
    elif [ "$therm" -ge 2 ] || awk -v b="$batt" -v p="$PAUSE_C" 'BEGIN{exit !(b>=p)}'; then
      verdict="PAUSE"
    elif [ "$verdict" = "PAUSE" ]; then
      if [ "$therm" -le 1 ] && awk -v b="$batt" -v r="$RESUME_C" 'BEGIN{exit !(b<r)}'; then
        verdict="OK"
      fi
    else
      verdict="OK"
    fi
  fi

  echo "$verdict therm=$therm batt=${batt}C load=$load ac=$ac at $ts" > "$STATE_FILE"
  echo "$ts verdict=$verdict therm=$therm batt=${batt}C load1=$load ac=$ac" >> "$LOG"
  [ "$prev" != "$verdict" ] && echo "$ts TRANSITION $prev -> $verdict (therm=$therm batt=${batt}C)" >> "$LOG"

  for f in "$DIR"/pgids/*; do
    [ -e "$f" ] || continue
    pg=$(basename "$f")
    if ! kill -0 -"$pg" 2>/dev/null; then rm -f "$f"; continue; fi
    case "$verdict" in
      PAUSE) kill -STOP -"$pg" 2>/dev/null ;;
      OK)    kill -CONT -"$pg" 2>/dev/null ;;
      STOP)  kill -CONT -"$pg" 2>/dev/null; kill -TERM -"$pg" 2>/dev/null; rm -f "$f" ;;
    esac
  done
  sleep 20
done
