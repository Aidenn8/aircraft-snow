"""Do streak objects move with the wind, or stay where they are?

usage: streak_motion.py SITE YYYY-MM-DD [--null] [--alpha 0.01] [--lag-min 20]

Snow made by an aircraft drifts with the wind once the aircraft has gone; fixed echoes
(clutter, wind farms, terrain-locked features) stay put. For each streak object matched at
p <= alpha (and, for comparison, a sample of all objects), the skeleton is sampled in the
scan nearest LAG minutes later, (a) moved by the wind and (b) where it was. The wind is the
HRRR mean over the layer from the beam to 1.5 km above it at the object's position and time,
so it does not depend on which flight matched. Score = median moved - median stayed (dB);
displacement = |wind| x lag. Prints both and writes data/events[_null]/SITE/DAY_motion.csv.
"""
import argparse
import datetime as dt
import pathlib

import numpy as np
import pandas as pd

from snow import epoch, match, pipeline, radar, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]
LAYER_M = 1500.0


def object_wind(prep, x, y, T):
    """HRRR mean wind between the beam and LAYER_M above it, at (x, y) and time T."""
    sw0, profs = prep["sweep0"], prep["profiles"]
    runs = np.array(sorted(profs.run.unique()))
    run_s = np.array([pd.Timestamp(r).timestamp() for r in runs])
    prof = profs[profs.run == runs[np.abs(run_s - T).argmin()]].sort_values("p")
    hb, _ = radar.beam_height_ground_range(np.hypot(x, y), sw0.elevation, sw0.radar_height)
    hgt, p = prof.HGT.to_numpy(float), prof.p.to_numpy(float)
    ok = np.isfinite(hgt)
    order = np.argsort(hgt[ok])
    p_lo = float(np.exp(np.interp(hb, hgt[ok][order], np.log(p[ok][order]))))
    p_hi = float(np.exp(np.interp(hb + LAYER_M, hgt[ok][order], np.log(p[ok][order]))))
    u, v = epoch.layer_mean_wind(prof, p_hi, p_lo)
    return float(u[0]), float(v[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("--null", action="store_true")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--lag-min", type=float, default=20.0)
    ap.add_argument("--sample", type=int, default=60)
    a = ap.parse_args()
    day = dt.date.fromisoformat(a.day)
    ev = ROOT / "data" / ("events_null" if a.null else "events") / a.site
    S = pd.read_parquet(ev / f"{day}_streaks.parquet")
    prep = pipeline.prepare(a.site, day, alt_max_ft=13000 if a.null else None)
    G = prep["grids"]
    kept, _, xs, ys = pipeline.streak_objects_day(prep)
    labs = {float(G.t[k]): lab for k, _, lab in kept}
    matched = S[(S.coverage > 0) & (S.p_value <= a.alpha)].assign(group="matched")
    rest = S[~S.index.isin(matched.index)]
    others = rest.sample(min(a.sample, len(rest)), random_state=0).assign(group="other")
    rows = []
    for _, o in pd.concat([matched, others]).iterrows():
        lab = labs.get(min(labs, key=lambda T: abs(T - o.t)))
        sx, sy = match.skeleton_xy(lab, int(o.label), xs, ys)
        if len(sx) < 4:
            continue
        T2 = o.t + a.lag_min * 60
        k2 = int(np.argmin(np.abs(G.t - T2)))
        if abs(G.t[k2] - T2) > 300:
            continue
        lag = G.t[k2] - o.t
        u, v = object_wind(prep, float(o.x), float(o.y), float(o.t))
        moved = G.lookup(k2, sx + u * lag, sy + v * lag)
        stayed = G.lookup(k2, sx, sy)
        moved = np.where(np.isfinite(moved), moved, streaks.FLOOR_DBZ)
        stayed = np.where(np.isfinite(stayed), stayed, streaks.FLOOR_DBZ)
        rows.append({"group": o.group, "t": o.t, "label": int(o.label), "x_km": o.x / 1e3, "y_km": o.y / 1e3,
                     "p_value": o.p_value, "u": u, "v": v, "disp_km": np.hypot(u, v) * lag / 1e3,
                     "z_moved": float(np.median(moved)), "z_stayed": float(np.median(stayed)),
                     "score_db": float(np.median(moved) - np.median(stayed))})
    out = pd.DataFrame(rows)
    out.to_csv(ev / f"{day}_motion.csv", index=False)
    pd.set_option("display.width", 200)
    for g, d in out.groupby("group"):
        far = d[d.disp_km >= 2.5]
        print(f"{a.site} {day} {g}: n={len(d)}, displacement median {d.disp_km.median():.1f} km; "
              f"with >= 2.5 km: n={len(far)}, score median {far.score_db.median():+.1f} dB, "
              f"share moving (score > 0): {np.mean(far.score_db > 0):.2f}")
    print(out[out.group == "matched"].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
