"""Day-swap test: do this day's flights match its streaks better than other days' flights?

usage: dayswap_null.py SITE YYYY-MM-DD [--null-alt FT] [--span 14]
       dayswap_null.py --confirmed [--events data/events] [--null-alt FT] [--span 14]

Rotations about the radar break any link between where streak-like echo forms and where
aircraft fly: terrain-locked convergence lines and bands along approach corridors, runways
and bands that both follow the wind, routing around storms. Flights from other days at the
same time of day (within +-SPAN days, snow.pipeline.dayswap_sets) keep all of that: the
same corridors, runway patterns and timing, through this day's weather and wind, but they
cannot have seeded this day's cloud. The test is conservative: an other-day flight that
happens to fly a seeder's path just ahead of it can share the seeder's streak.
The streak table, day statistic and false discovery rates are computed exactly as for
rotations, with each other day as one placebo set. With --null-alt, null mode (tracks below
that altitude, warm points only), as for warm-season null days.

--confirmed runs it for every radar-day in --events that was confirmed with 200 rotations
(the census's candidate days) and has no day-swap result yet.
Writes data/events[_null]/SITE/YYYY-MM-DD_dayswap.json and _dayswap_streaks.parquet.
"""
import argparse
import datetime as dt
import json
import pathlib
import sys
import time
import traceback

import numpy as np

from snow import match, pipeline

ROOT = pathlib.Path(__file__).resolve().parents[1]


def run(site: str, day: dt.date, null_alt=None, span: int = 14) -> dict:
    tic = time.time()
    null = null_alt is not None
    prep = pipeline.prepare(site, day, alt_max_ft=null_alt)
    offsets = tuple(k for k in range(-span, span + 1) if k != 0)
    sets, used = pipeline.dayswap_sets(prep, site, day, offsets, alt_max_ft=null_alt)
    tmax = None if null else match.T_SEED_MAX_C
    table = pipeline.match_day(prep, max_temp_c=tmax, placebo_sets=sets)
    st = match.day_statistic(table)
    st.pop("E_null", None)
    out = {"site": site, "day": str(day), "null_mode": null, "offsets": used, "stats": st,
           "fdr": match.empirical_fdr(table).to_dict(orient="records") if not table.empty else [],
           "seconds": round(time.time() - tic, 1)}
    d = ROOT / "data" / ("events_null" if null else "events") / site
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{day}_dayswap.json").write_text(json.dumps(out, indent=1, default=float))
    if not table.empty:
        t = table.copy()
        t["placebo_cov"] = t.placebo_cov.map(lambda v: np.asarray(v, np.float32).tolist())
        t.to_parquet(d / f"{day}_dayswap_streaks.parquet")
    return out


def report(out: dict):
    st = out["stats"]
    print(f"{out['site']} {out['day']} day-swap ({len(out['offsets'])} other days): E = {st['E']:.1f}, "
          f"null median {st.get('E_null_median', 0):.1f}, q95 {st.get('E_null_q95', 0):.1f}, "
          f"rank p = {st['p_day']:.3f}, gamma p = {st.get('p_day_gamma', float('nan')):.3g}, "
          f"flights matched {st.get('n_flights_matched', 0)} ({out['seconds']:.0f}s)", flush=True)
    for f in out["fdr"]:
        print("   ", {k: round(v, 3) for k, v in f.items()}, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site", nargs="?")
    ap.add_argument("day", nargs="?")
    ap.add_argument("--null-alt", type=float, default=None)
    ap.add_argument("--span", type=int, default=14)
    ap.add_argument("--confirmed", action="store_true")
    ap.add_argument("--events", default="")
    a = ap.parse_args()
    if not a.confirmed:
        report(run(a.site, dt.date.fromisoformat(a.day), a.null_alt, a.span))
        return
    root = pathlib.Path(a.events) if a.events else ROOT / "data" / ("events_null" if a.null_alt else "events")
    for path in sorted(root.glob("*/*.json")):
        if "_" in path.stem:
            continue
        s = json.loads(path.read_text())
        if s.get("streaks", {}).get("n_rotations", 0) < 100:
            continue
        if path.with_name(f"{path.stem}_dayswap.json").exists():
            continue
        try:
            report(run(s["site"], dt.date.fromisoformat(s["day"]), a.null_alt, a.span))
        except Exception as e:
            print(f"{s['site']} {s['day']}: FAILED {e!r}", flush=True)
            traceback.print_exc(file=sys.stderr)


if __name__ == "__main__":
    main()
