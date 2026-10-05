"""Run the event tests for every census radar over a range of analysis days.

usage: census.py START END [--workers N] [--sites config/census_sites.csv] [--force]

For each analysis day D (12 UTC D to 12 UTC D+1) whose ADS-B days D and D+1 are present,
the continental files are split per radar (snow.adsb_split) and every radar-day is run
(snow.pipeline.run_day, streak test only) in a pool of worker processes. Results land in data/events/SITE/.
Days already done are skipped unless --force. A one-line summary per radar-day goes to
stdout; tracebacks to stderr.
"""
import argparse
import concurrent.futures as cf
import datetime as dt
import sys
import traceback

from snow import adsb, adsb_split


def one(site: str, day: dt.date, force: bool) -> str:
    from snow import pipeline
    try:
        s = pipeline.run_day(site, day, force=force, event_test=False)
    except Exception as e:
        traceback.print_exc(file=sys.stderr)
        return f"{site} {day}: FAILED {e!r}"
    if "screened_out" in s:
        return f"{site} {day}: screened out ({s['screened_out']})"
    sk = s.get("streaks", {})
    return (f"{site} {day}: flights {s['n_flights']:5d} streaks {sk.get('n_streaks', 0):4d} "
            f"matched {sk.get('n_matched', 0):3d} E={sk.get('E', 0):6.1f} p_day={sk.get('p_day', 1):.3f} "
            f"p_gamma={sk.get('p_day_gamma', float('nan')):.2g} "
            f"p_rev={s.get('streaks_reverse', {}).get('p_day', float('nan')):.3f} "
            f"(R={sk.get('n_rotations', 0)}) ({s['seconds']:.0f}s)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("start")
    ap.add_argument("end")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    sites = adsb_split.census_sites()
    d0, d1 = dt.date.fromisoformat(a.start), dt.date.fromisoformat(a.end)
    jobs = []
    d = d0
    while d <= d1:
        if adsb.day_path(d).exists() and adsb.day_path(d + dt.timedelta(days=1)).exists():
            adsb_split.split_day(d, sites)
            adsb_split.split_day(d + dt.timedelta(days=1), sites)
            jobs += [(s, d) for s in sites]
        else:
            print(f"{d}: no ADS-B yet", flush=True)
        d += dt.timedelta(days=1)
    print(f"{len(jobs)} radar-days queued", flush=True)
    with cf.ProcessPoolExecutor(a.workers) as ex:
        futs = [ex.submit(one, s, d, a.force) for s, d in jobs]
        for f in cf.as_completed(futs):
            print(f.result(), flush=True)


if __name__ == "__main__":
    main()
