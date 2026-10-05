"""Robustness check of the rotation null with orientation-preserving placebos.

usage: shift_null.py SITE YYYY-MM-DD [--null-alt FT] [--inner-km 49]

Rotating flights about the radar changes their orientation. If streak-like echo has a
preferred orientation that real flight paths share (bands and runways both follow the
wind), rotations would favour the real flights and the rotation p-value would be too small.
Translating every flight by 8-16 km (match.shift_placebos) keeps each track's orientation,
and its orientation relative to the wind, exactly, and only moves it. Both placebo families
are run here on the same points: tracks within --inner-km of the radar, so that shifted
copies stay inside the analysed domain. With --null-alt, null mode (tracks below that
altitude, no temperature requirement), as for warm-season null days.
Writes data/events[_null]/SITE/YYYY-MM-DD_shiftnull.json.
"""
import argparse
import datetime as dt
import json
import pathlib

import numpy as np
import pandas as pd

from snow import match, pipeline

ROOT = pathlib.Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("--null-alt", type=float, default=None)
    ap.add_argument("--inner-km", type=float, default=49.0)
    a = ap.parse_args()
    day = dt.date.fromisoformat(a.day)
    null = a.null_alt is not None
    prep = pipeline.prepare(a.site, day, alt_max_ft=a.null_alt)
    G = prep["grids"]
    pts, dr = pipeline.matching_points(prep)
    inner = np.hypot(pts.x, pts.y).to_numpy() <= a.inner_km * 1e3
    pts, dr = pts[inner].reset_index(drop=True), dr[inner]
    kept, _, xs, ys = pipeline.streak_objects_day(prep)
    tmax = None if null else match.T_SEED_MAX_C
    out = {"site": a.site, "day": str(day), "null_mode": null, "inner_km": a.inner_km}
    for name, angles in (("rotation", pipeline.SCREEN_ANGLES), ("shift", match.shift_placebos())):
        rows = [match.match_streaks(G, G.t[k], lab, obj, xs, ys, pts, dr, angles=angles, max_temp_c=tmax)
                for k, obj, lab in kept]
        rows = [r for r in rows if not r.empty]
        table = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
        st = match.day_statistic(table)
        st.pop("E_null", None)
        out[name] = st
        print(f"{name:8s}: E = {st['E']:6.1f}, null median {st.get('E_null_median', 0):6.1f}, "
              f"q95 {st.get('E_null_q95', 0):6.1f}, rank p = {st['p_day']:.3f}, "
              f"gamma p = {st.get('p_day_gamma', float('nan')):.3g}, flights matched {st.get('n_flights_matched', 0)}",
              flush=True)
    d = ROOT / "data" / ("events_null" if null else "events") / a.site
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{day}_shiftnull.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
