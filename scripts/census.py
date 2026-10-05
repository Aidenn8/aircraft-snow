"""Run the event tests for every census radar over a range of analysis days.

usage: census.py START END [--workers N] [--sites KFTG,KLOT,...] [--force]

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


def one(site: str, day: dt.date, force: bool, null_alt: float | None = None) -> str:
    from snow import pipeline
    try:
        s = pipeline.run_day(site, day, force=force, event_test=False, null_alt_max_ft=null_alt)
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
    ap.add_argument("--sites", default="", help="comma-separated subset of the census radars")
    ap.add_argument("--null-alt", type=float, default=None,
                    help="null mode (warm-season calibration): tracks below this altitude, results "
                         "to data/events_null")
    ap.add_argument("--step", type=int, default=1, help="every STEP-th day")
    a = ap.parse_args()
    all_sites = adsb_split.census_sites()
    sites = [s for s in a.sites.split(",") if s] or all_sites
    unknown = set(sites) - set(all_sites)
    if unknown:
        raise SystemExit(f"not census radars: {sorted(unknown)}")
    d0, d1 = dt.date.fromisoformat(a.start), dt.date.fromisoformat(a.end)
    # days are split (in this process) and submitted one at a time, so workers start on the
    # first days while later ones are still being split
    n = 0
    futs = []
    with cf.ProcessPoolExecutor(a.workers) as ex:
        d = d0
        while d <= d1:
            if (d - d0).days % a.step:
                d += dt.timedelta(days=1)
                continue
            if adsb.day_path(d).exists() and adsb.day_path(d + dt.timedelta(days=1)).exists():
                # split for every census radar while the continental file is open (cheap,
                # and later phases can reuse it)
                adsb_split.split_day(d, all_sites)
                adsb_split.split_day(d + dt.timedelta(days=1), all_sites)
                for s in sites:
                    f = ex.submit(one, s, d, a.force, a.null_alt)
                    f.add_done_callback(lambda f: print(f.result(), flush=True))
                    futs.append(f)
                n += len(sites)
            else:
                print(f"{d}: no ADS-B yet", flush=True)
            d += dt.timedelta(days=1)
        print(f"{n} radar-days queued", flush=True)
        cf.wait(futs)

if __name__ == "__main__":
    main()
