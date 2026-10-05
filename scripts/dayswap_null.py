"""Day-swap test: do this day's flights match its streaks better than other days' flights?

usage: dayswap_null.py SITE YYYY-MM-DD [--null-alt FT] [--span 14]

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
Writes data/events[_null]/SITE/YYYY-MM-DD_dayswap.json and _dayswap_streaks.parquet.
"""
import argparse
import datetime as dt
import json
import pathlib
import time

import numpy as np

from snow import match, pipeline

ROOT = pathlib.Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("--null-alt", type=float, default=None)
    ap.add_argument("--span", type=int, default=14)
    a = ap.parse_args()
    tic = time.time()
    day = dt.date.fromisoformat(a.day)
    null = a.null_alt is not None
    prep = pipeline.prepare(a.site, day, alt_max_ft=a.null_alt)
    offsets = tuple(k for k in range(-a.span, a.span + 1) if k != 0)
    sets, used = pipeline.dayswap_sets(prep, a.site, day, offsets, alt_max_ft=a.null_alt)
    print(f"{len(sets)} other days as placebo sets: {used}", flush=True)
    tmax = None if null else match.T_SEED_MAX_C
    table = pipeline.match_day(prep, max_temp_c=tmax, placebo_sets=sets)
    st = match.day_statistic(table)
    st.pop("E_null", None)
    out = {"site": a.site, "day": str(day), "null_mode": null, "offsets": used, "stats": st,
           "fdr": match.empirical_fdr(table).to_dict(orient="records") if not table.empty else [],
           "seconds": round(time.time() - tic, 1)}
    d = ROOT / "data" / ("events_null" if null else "events") / a.site
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{day}_dayswap.json").write_text(json.dumps(out, indent=1, default=float))
    if not table.empty:
        t = table.copy()
        t["placebo_cov"] = t.placebo_cov.map(lambda v: np.asarray(v, np.float32).tolist())
        t.to_parquet(d / f"{day}_dayswap_streaks.parquet")
    print(f"{a.site} {day} day-swap: E = {st['E']:.1f}, null median {st.get('E_null_median', 0):.1f}, "
          f"q95 {st.get('E_null_q95', 0):.1f}, rank p = {st['p_day']:.3f}, gamma p = "
          f"{st.get('p_day_gamma', float('nan')):.3g}, flights matched {st.get('n_flights_matched', 0)}", flush=True)
    for f in out["fdr"]:
        print({k: round(v, 3) for k, v in f.items()}, flush=True)


if __name__ == "__main__":
    main()
