"""What share of aircraft crossing the seeding layer left a radar-visible streak?

usage: seeding_efficiency.py SITE YYYY-MM-DD [--alpha 0.01] [--layer-ft 1000] [--pad-min 60]

On an event day, the seeding layer and period are taken from the confident matches
(p <= alpha): altitudes within --layer-ft of the matched seeding altitudes, passage times
from --pad-min before the first to --pad-min after the last match. Every flight with a
track segment in that layer and period is "exposed". For each exposed flight, and only
where its air was echo-free before it arrived (first passage), the question is whether a
streak object later lay along its drifted track: the best single-flight F1 against all
streak objects in the scans 10-50 min after its passage (no rotations; this is a
description of the outcome, not a test).

Writes data/events/SITE/YYYY-MM-DD_exposure.csv with one row per exposed flight.
"""
import argparse
import datetime as dt
import pathlib

import numpy as np
import pandas as pd

from snow import flights, match, pipeline, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]
PLACEBO_ROT = (60.0, 120.0, 180.0, 240.0, 300.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--layer-ft", type=float, default=1000.0)
    ap.add_argument("--pad-min", type=float, default=60.0)
    ap.add_argument("--f1", type=float, default=0.5)
    a = ap.parse_args()
    day = dt.date.fromisoformat(a.day)
    ev = ROOT / "data" / "events" / a.site
    S = pd.read_parquet(ev / f"{day}_streaks.parquet")
    conf = S[(S.coverage > 0) & (S.p_value <= a.alpha)]
    desc = pd.read_csv(ev / f"{day}_described.csv")
    if conf.empty or desc.empty:
        raise SystemExit("no confident matches")
    alt_lo = desc.seed_alt_ft.min() - a.layer_ft
    alt_hi = desc.seed_alt_ft.max() + a.layer_ft
    t_lo = conf.t_pass.min() - a.pad_min * 60
    t_hi = conf.t_pass.max() + a.pad_min * 60
    prep = pipeline.prepare(a.site, day)
    G = prep["grids"]
    pts, dr = pipeline.matching_points(prep)
    calls = pd.concat([flights.load_calls(day, a.site), flights.load_calls(day + dt.timedelta(days=1), a.site)])
    ft = flights.flight_table(prep["positions"], calls)
    xs = G.x0 + G.res * np.arange(G.z.shape[2])
    ys = G.y0 + G.res * np.arange(G.z.shape[1])
    in_layer = (pts.alt_ft >= alt_lo) & (pts.alt_ft <= alt_hi) & (pts.t >= t_lo) & (pts.t <= t_hi)
    exposed = np.unique(pts.flight[in_layer])
    # streak objects of every scan in the period (fixed echoes removed), with skeletons
    kept, _, _, _ = pipeline.streak_objects_day(prep)
    skel = {}
    for k, obj, lab in kept:
        if t_lo + match.AGE_MIN_S <= G.t[k] <= t_hi + match.AGE_MAX_S:
            skel[k] = [match.skeleton_xy(lab, int(L), xs, ys) for L in obj.label]
    seeders = set(int(f) for f in conf.flight)
    rows = []
    for f in exposed:
        sel = (pts.flight.to_numpy() == f) & in_layer.to_numpy()
        p = pts[sel]
        d = dr[sel]
        t_pass = float(np.median(p.t))
        # first passage: the air along this flight's layer segment was echo-free beforehand
        before = np.nonzero((G.t >= t_pass - 30 * 60) & (G.t <= t_pass - 3 * 60))[0]
        quiet = []
        for k in before:
            age = G.t[k] - p.t.to_numpy(float)
            z = G.lookup(k, p.x.to_numpy() + d[:, 0] * age, p.y.to_numpy() + d[:, 1] * age)
            quiet.append(z[np.isfinite(z)])
        qv = np.concatenate(quiet) if quiet else np.array([])
        quiet_frac = float(np.mean(qv < match.QUIET_DBZ)) if qv.size else np.nan
        def best_f1(rot):
            b = 0.0
            for k, sks in skel.items():
                T = G.t[k]
                age = T - t_pass
                if not (match.AGE_MIN_S <= age <= match.AGE_MAX_S):
                    continue
                # the same physical requirements as the census matcher (fall speed, temperature)
                mk = match._physical_mask(p, T - p.t.to_numpy(float))
                if mk.sum() < 4:
                    continue
                for sx, sy in sks:
                    if len(sx) < 4:
                        continue
                    c = match._candidates(sx, sy, p.x.to_numpy()[mk], p.y.to_numpy()[mk],
                                          p.t.to_numpy(float)[mk], T, d[mk], p.flight.to_numpy()[mk], rot=rot)
                    if c:
                        b = max(b, c[0][0])
            return b

        best = best_f1(0.0)
        # chance level: the same flight rotated about the radar to other bearings
        placebo = [best_f1(r) for r in PLACEBO_ROT]
        prof = flights.interp_profile_to_points(p, prep["profiles"], fields=("temp",))
        info = ft.loc[f] if f in ft.index else None
        rows.append({"flight": int(f), "callsign": getattr(info, "callsign", ""),
                     "icao_type": getattr(info, "icao_type", ""), "engine_class": getattr(info, "engine_class", "unknown"),
                     "phase": getattr(info, "phase", ""), "t_pass": t_pass, "alt_ft": float(p.alt_ft.median()),
                     "temp_c": float(prof.temp.median()), "n_points": int(sel.sum()),
                     "quiet_before": quiet_frac, "first_passage": bool(quiet_frac >= match.NEW_MIN_QUIET),
                     "best_f1": best, "made_streak": best >= a.f1, "confident_match": int(f) in seeders,
                     "placebo_rate": float(np.mean(np.array(placebo) >= a.f1))})
    out = pd.DataFrame(rows)
    out.to_csv(ev / f"{day}_exposure.csv", index=False)
    fp = out[out.first_passage]
    print(f"layer {alt_lo:.0f}-{alt_hi:.0f} ft, period {pd.Timestamp(t_lo, unit='s'):%H:%M}-"
          f"{pd.Timestamp(t_hi, unit='s'):%H:%M} UTC: {len(out)} exposed flights, {len(fp)} first passages")
    g = fp.groupby("engine_class").agg(n=("flight", "size"), made_streak=("made_streak", "mean"),
                                       chance=("placebo_rate", "mean"), confident=("confident_match", "mean"),
                                       temp=("temp_c", "median"))
    # share above chance, scaled to what chance leaves room for
    g["excess"] = (g.made_streak - g.chance) / (1 - g.chance)
    print(g.round(3))
    bins = [-30, -20, -17, -14, -11, -8, 0]
    fp = fp.assign(tbin=pd.cut(fp.temp_c, bins))
    print(fp[fp.engine_class == "jet"].groupby("tbin", observed=True).agg(
        n=("flight", "size"), made_streak=("made_streak", "mean"), chance=("placebo_rate", "mean")).round(3))


if __name__ == "__main__":
    main()
