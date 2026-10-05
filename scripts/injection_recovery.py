"""Recall of the streak-matching detector from synthetic streaks injected into real radar.

usage: injection_recovery.py SITE DAY [DAY ...] [--n N] [--peaks 5,10,15,20] [--seed S]

For each trial a fall speed is drawn from V_FALL (narrower than the detector's allowed
range, so the test does not assume what it checks) and a supercooled layer thickness from
LAYER_M. A real flight is drawn at random from those with track points that could
physically seed: between 6000 and 14000 ft, 10-60 km from the radar, HRRR temperature at or
below the detector's seeding limit, and high enough above the beam that their ice reaches it
10-50 min later at that fall speed. The layer is centred on the height above the beam of one
of the flight's points, chosen at random, and only the flight's points inside it seed; the
trial needs at least MIN_DENSE of them (10 km of track), otherwise another flight is drawn.
A flight crossing the layer level leaves a long streak, one descending through it a short
one, as in nature. A streak (snow.inject, physical timing) is added along the seeding
points' drifted track, with peak reflectivity from --peaks and a width when it appears
drawn from WIDTH_M (the widths of the detected Denver streaks; so this is recall for streaks
shaped like those). The streak detector and single-flight matcher (with rotation placebos)
are then run on the scans 10-60 min after the flight's passage. A trial is recovered_any
when a streak object overlapping the injected streak's footprint is matched to some flight
with p <= ALPHA (event recall), and recovered when it is matched to the injected flight
itself (attribution; in a stream of aircraft on one path, one just ahead can take the
match). Results go to data/injection/SITE_DAYS_sSEED.parquet.
"""
import argparse
import datetime as dt
import pathlib
import time

import numpy as np
import pandas as pd

from snow import inject, match, pipeline, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]
ALPHA = 0.05
MIN_DENSE = 40                      # 10 km of interpolated track
BAND = (6000.0, 14000.0)
V_FALL = (0.4, 1.2)                 # m/s, drawn uniformly per trial
LAYER_M = (300.0, 800.0)            # supercooled layer depth, drawn uniformly per trial
WIDTH_M = (1500.0, 3000.0)          # streak width above 3 dBZ when it appears (Denver: 1.8-3.1 km)
MAX_DRAWS = 50


def run_trials(site, day, n, peaks, rng):
    prep = pipeline.prepare(site, day)
    G = prep["grids"]
    pts, dr = pipeline.matching_points(prep)
    xs = G.x0 + G.res * np.arange(G.z.shape[2])
    ys = G.y0 + G.res * np.arange(G.z.shape[1])
    rng_pts = np.hypot(pts.x, pts.y)
    t0 = prep["t0"].timestamp()
    sel = pts[(pts.alt_ft >= BAND[0]) & (pts.alt_ft < BAND[1]) & (rng_pts > 10e3) & (rng_pts < 60e3)
              & (pts.t > t0 + 1800) & (pts.t < t0 + 86400 - 5400)
              & (pts.temp_c <= match.T_SEED_MAX_C) & (pts.dh > 0)]
    rows = []
    dh_all = sel.dh.to_numpy()
    fl_all = sel.flight.to_numpy()
    for trial in range(n):
        v = float(rng.uniform(*V_FALL))
        depth = float(rng.uniform(*LAYER_M))
        arrival = dh_all / v
        ok = (arrival >= match.AGE_MIN_S) & (arrival <= match.AGE_MAX_S)
        counts = pd.Series(fl_all[ok]).value_counts()
        candidates = counts[counts >= MIN_DENSE].index.to_numpy()
        fp = None
        for _ in range(MAX_DRAWS if len(candidates) else 0):
            f = rng.choice(candidates)
            mine = np.nonzero(ok & (fl_all == f))[0]
            centre = dh_all[rng.choice(mine)]
            layer = mine[np.abs(dh_all[mine] - centre) <= depth / 2]
            if len(layer) >= MIN_DENSE:
                fp = sel.iloc[layer]
                break
        if fp is None:
            continue
        idx = fp.index.to_numpy()
        peak = float(rng.choice(peaks))
        width = float(rng.uniform(*WIDTH_M))
        s0 = inject.sigma_for_width(width, peak)
        Gi, err = inject.inject(G, fp.x.to_numpy(), fp.y.to_numpy(), fp.t.to_numpy(float), dr[idx],
                                peak, rng=rng, dh=fp.dh.to_numpy(), v_fall=v, s0_m=s0)
        arr = fp.dh.to_numpy() / v
        t_pass = float(fp.t.median())
        scans = np.nonzero((Gi.t >= t_pass + 10 * 60) & (Gi.t <= t_pass + 60 * 60))[0]
        best_p, best_cov, found_obj, best_any = 1.0, 0.0, False, 1.0
        for k in scans:
            T = Gi.t[k]
            z = np.where(np.isnan(Gi.z[k]), streaks.FLOOR_DBZ, Gi.z[k])
            obj, lab = streaks.find_streaks(streaks.GridScan(T, xs, ys, z, np.isfinite(Gi.z[k])))
            if obj.empty:
                continue
            # objects overlapping the injected streak: the injection's own footprint at T
            age = T - fp.t.to_numpy(float)
            w = inject.envelope_fall(age, arr) > 0.3
            if not w.any():
                continue
            px = fp.x.to_numpy()[w] + (dr[idx][w, 0] + err[0]) * age[w]
            py = fp.y.to_numpy()[w] + (dr[idx][w, 1] + err[1]) * age[w]
            ix = np.clip(np.rint((px - xs[0]) / G.res).astype(int), 0, len(xs) - 1)
            iy = np.clip(np.rint((py - ys[0]) / G.res).astype(int), 0, len(ys) - 1)
            labs = set(np.unique(lab[iy, ix])) - {0}
            if not labs:
                continue
            found_obj = True
            o = obj[obj.label.isin(labs)]
            m_all = match.match_streaks(Gi, T, lab, o, xs, ys, pts, dr)
            # the streak found and matched to some flight (event recall) ...
            if not m_all.empty:
                best_any = min(best_any, float(m_all.p_value.min()))
            # ... and to the injected flight itself (attribution)
            m = m_all[m_all.flight == f] if not m_all.empty else m_all
            if m.empty:
                continue
            j = m.p_value.idxmin()
            if m.p_value[j] < best_p or (m.p_value[j] == best_p and m.coverage[j] > best_cov):
                best_p, best_cov = float(m.p_value[j]), float(m.coverage[j])
            if best_p <= ALPHA:
                break
        rows.append({"site": site, "day": str(day), "trial": trial, "flight": int(f), "peak_dbz": peak,
                     "n_points": len(fp), "alt_ft": float(fp.alt_ft.median()), "v_fall": v, "layer_m": depth, "width_m": width,
                     "alt_range_ft": float(fp.alt_ft.max() - fp.alt_ft.min()),
                     "dh_m": float(fp.dh.median()), "temp_c": float(fp.temp_c.median()),
                     "wind_err_u": err[0], "wind_err_v": err[1], "object_found": found_obj,
                     "best_p": best_p, "best_cov": best_cov, "recovered": best_p <= ALPHA,
                     "best_p_any": best_any, "recovered_any": best_any <= ALPHA})
        print(f"{site} {day} trial {trial}: peak {peak:4.0f} dBZ object {found_obj!s:5} p {best_p:.3f}",
              flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("days", nargs="+")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--peaks", default="5,10,15,20")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    peaks = [float(p) for p in a.peaks.split(",")]
    rows = []
    for d in a.days:
        tic = time.time()
        rows += run_trials(a.site, dt.date.fromisoformat(d), a.n, peaks, rng)
        print(f"{d}: {time.time() - tic:.0f} s", flush=True)
    out = ROOT / "data" / "injection" / f"{a.site}_{'_'.join(a.days)}_s{a.seed}.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_parquet(out)
    print(df.groupby("peak_dbz")[["object_found", "recovered_any", "recovered"]].mean())


if __name__ == "__main__":
    main()
