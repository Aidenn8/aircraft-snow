"""Dual-polarisation signature of matched streaks versus other echo (Level II, lowest tilt).

usage: dualpol_streaks.py SITE YYYY-MM-DD START END [--alpha 0.005]

For each Level II volume in the window, the lowest sweep carrying ZDR and RHO is sampled
at the pixels of streak objects matched at p <= alpha (from the nearest Level III scan),
and at other echo pixels (>= 0 dBZ, outside any streak object, 10-60 km range).
Prints medians and quartiles; writes data/events/SITE/YYYY-MM-DD_dualpol.parquet.
"""
import argparse
import datetime as dt
import pathlib

import numpy as np
import pandas as pd

from snow import pipeline, radar, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("start")
    ap.add_argument("end")
    ap.add_argument("--alpha", type=float, default=0.005)
    a = ap.parse_args()
    day = dt.date.fromisoformat(a.day)
    UTC = dt.timezone.utc
    ta = dt.datetime.fromisoformat(a.start).replace(tzinfo=UTC)
    tb = dt.datetime.fromisoformat(a.end).replace(tzinfo=UTC)
    S = pd.read_parquet(ROOT / "data" / "events" / a.site / f"{day}_streaks.parquet")
    S = S[(S.coverage > 0) & (S.p_value <= a.alpha)]
    prep = pipeline.prepare(a.site, day)
    G = prep["grids"]
    xs = G.x0 + G.res * np.arange(G.z.shape[2])
    ys = G.y0 + G.res * np.arange(G.z.shape[1])
    X, Y = np.meshgrid(xs, ys)
    R = np.hypot(X, Y)
    vols = []
    d = ta.date()
    while d <= tb.date():
        vols += [v for v in radar.list_level2(a.site, d) if ta <= v[0] <= tb]
        d += dt.timedelta(days=1)
    rows = []
    for t, key, _ in vols:
        k = int(np.argmin(np.abs(G.t - t.timestamp())))
        if abs(G.t[k] - t.timestamp()) > 120:
            continue
        z = np.where(np.isnan(G.z[k]), streaks.FLOOR_DBZ, G.z[k])
        obj, lab = streaks.find_streaks(streaks.GridScan(G.t[k], xs, ys, z, np.isfinite(G.z[k])))
        matched = S[np.abs(S.t - G.t[k]) < 1.0].label.astype(int).tolist()   # same scan (epoch s)
        if not matched:
            continue
        in_streak = np.isin(lab, matched)
        any_streak = lab > 0
        sweeps = radar.read_level2(radar.fetch(radar.L2_BUCKET, key), moments=("REF", "ZDR", "RHO"))
        sw = next(s for s in sweeps if "ZDR" in s.data and "RHO" in s.data)
        vals = {m: radar.sample(sw, X, Y, m) for m in ("REF", "ZDR", "RHO")}
        bg = (~any_streak) & (vals["REF"] >= 0) & (R >= 10e3) & (R <= 60e3)
        st = in_streak & (vals["REF"] >= 0)
        for name, mask in (("streak", st), ("other echo", bg)):
            rows.append(pd.DataFrame({"time": t, "region": name, "REF": vals["REF"][mask],
                                      "ZDR": vals["ZDR"][mask], "RHO": vals["RHO"][mask]}))
    df = pd.concat(rows, ignore_index=True).dropna()
    df.to_parquet(ROOT / "data" / "events" / a.site / f"{day}_dualpol.parquet")
    for name, g in df.groupby("region"):
        q = g[["REF", "ZDR", "RHO"]].quantile([0.25, 0.5, 0.75]).round(2)
        print(f"{name}: {len(g)} gates\n{q.to_string()}")
    # same reflectivity range only, to separate the effect of intensity
    m = (df.REF >= 10) & (df.REF <= 20)
    print("restricted to 10-20 dBZ:")
    print(df[m].groupby("region")[["ZDR", "RHO"]].median().round(2))


if __name__ == "__main__":
    main()
