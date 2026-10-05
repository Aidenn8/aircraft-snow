"""Run a command under the thermal guard.

    python tools/thermal/guarded_run.py -- <command> [args...]

The command gets its own process group, registered in tools/thermal/pgids/, so
daemon.sh can freeze (SIGSTOP), resume (SIGCONT) or end (SIGTERM) it and all of
its worker processes together. It also runs at reduced priority (nice 10).
Refuses to start while the guard says PAUSE or STOP, or if the daemon is not running.
"""
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
STATE = HERE / "state"


def guard_verdict():
    try:
        text = STATE.read_text()
    except FileNotFoundError:
        return None, 0.0
    return text.split()[0], time.time() - STATE.stat().st_mtime


def main():
    args = sys.argv[1:]
    if args and args[0] == "--":
        args = args[1:]
    if not args:
        sys.exit(__doc__)
    verdict, age = guard_verdict()
    if verdict is None or age > 90:
        sys.exit("thermal guard daemon is not running (state file missing or stale); refusing to start")
    while verdict == "PAUSE":
        print("thermal guard says PAUSE; waiting to start", file=sys.stderr, flush=True)
        time.sleep(30)
        verdict, _ = guard_verdict()
    if verdict == "STOP":
        sys.exit("thermal guard says STOP; refusing to start")
    os.setpgrp()
    os.nice(10)
    (HERE / "pgids").mkdir(exist_ok=True)
    (HERE / "pgids" / str(os.getpgrp())).write_text(" ".join(args))
    os.execvp(args[0], args)


if __name__ == "__main__":
    main()
