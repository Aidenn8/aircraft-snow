"""Composite life cycle of matched streaks: reflectivity across the drifted track versus age.

usage: streak_lifecycle.py SITE YYYY-MM-DD [--alpha 0.005]

For every streak matched at p <= alpha, the matched flight's seeding points (its points on
the streak at the matched scan) are followed with the matched drift through every scan from
-30 to +120 min of age. At each scan, reflectivity is sampled on lines parallel to the
track at cross-track offsets -8..+8 km. Each flight is used once (its earliest matched
streak). Writes data/events/SITE/YYYY-MM-DD_lifecycle.npz and a figure in review/.
"""
import argparse
import datetime as dt
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from snow import epoch, match, pipeline, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]
OFFSETS = np.arange(-8000, 8001, 500)
AGE_BINS = np.arange(-30 * 60, 120 * 60 + 1, 5 * 60)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("--alpha", type=float, default=0.005)
    a = ap.parse_args()
    day = dt.date.fromisoformat(a.day)
    ev = ROOT / "data" / "events" / a.site
    S = pd.read_parquet(ev / f"{day}_streaks.parquet")
    S = S[(S.coverage > 0) & (S.p_value <= a.alpha)].sort_values("t").drop_duplicates("flight")
    prep = pipeline.prepare(a.site, day)
    G = prep["grids"]
    pts, dr = pipeline.matching_points(prep)
    xs = G.x0 + G.res * np.arange(G.z.shape[2])
    ys = G.y0 + G.res * np.arange(G.z.shape[1])
    nA, nO = len(AGE_BINS) - 1, len(OFFSETS)
    acc = np.zeros((nA, nO))
    cnt = np.zeros((nA, nO))
    per_flight = []
    for _, o in S.iterrows():
        k = int(np.argmin(np.abs(G.t - o.t)))
        z = np.where(np.isnan(G.z[k]), streaks.FLOOR_DBZ, G.z[k])
        obj, lab = streaks.find_streaks(streaks.GridScan(o.t, xs, ys, z, np.isfinite(G.z[k])))
        sx, sy = match.skeleton_xy(lab, int(o.label), xs, ys)
        mp = match.matched_points(sx, sy, pts, dr, o.t, int(o.flight), o.du, o.dv)
        if len(mp) < 5:
            continue
        idx = mp.index.to_numpy()
        hd = epoch.heading(mp.assign(flight=0).reset_index(drop=True))
        perp = np.c_[-hd[:, 1], hd[:, 0]]
        u = dr[idx, 0] + o.du
        v = dr[idx, 1] + o.dv
        f_acc = np.zeros((nA, nO))
        f_cnt = np.zeros((nA, nO))
        for kk, T in enumerate(G.t):
            age = T - mp.t.to_numpy(float)
            m = (age >= AGE_BINS[0]) & (age < AGE_BINS[-1])
            if not m.any():
                continue
            b = np.digitize(age[m], AGE_BINS) - 1
            px = mp.x.to_numpy()[m] + u[m] * age[m]
            py = mp.y.to_numpy()[m] + v[m] * age[m]
            for j, off in enumerate(OFFSETS):
                zz = G.lookup(kk, px + off * perp[m, 0], py + off * perp[m, 1])
                ok = np.isfinite(zz)
                np.add.at(f_acc[:, j], b[ok], zz[ok])
                np.add.at(f_cnt[:, j], b[ok], 1)
        acc += f_acc
        cnt += f_cnt
        with np.errstate(invalid="ignore"):
            per_flight.append(f_acc / f_cnt)
    with np.errstate(invalid="ignore"):
        Z = acc / cnt
    np.savez(ev / f"{day}_lifecycle.npz", Z=Z, cnt=cnt, offsets=OFFSETS, age_bins=AGE_BINS,
             per_flight=np.array(per_flight), n_flights=len(per_flight))
    ages = 0.5 * (AGE_BINS[:-1] + AGE_BINS[1:]) / 60
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    im = axes[0].pcolormesh(OFFSETS / 1e3, ages, Z, cmap="turbo", vmin=-10, vmax=15, shading="nearest")
    axes[0].set_xlabel("cross-track offset (km)")
    axes[0].set_ylabel("minutes since aircraft passage")
    axes[0].set_title(f"{a.site} {day}: mean reflectivity around {len(per_flight)} seeding flights")
    fig.colorbar(im, ax=axes[0], label="dBZ (floor -10)")
    centre = Z[:, nO // 2]
    side = np.nanmean(Z[:, (np.abs(OFFSETS) >= 5000)], axis=1)
    axes[1].plot(ages, centre, label="on the drifted track")
    axes[1].plot(ages, side, label="5-8 km to either side")
    axes[1].axvline(0, color="k", lw=0.5)
    axes[1].set_xlabel("minutes since aircraft passage")
    axes[1].set_ylabel("mean dBZ")
    axes[1].legend()
    out = ROOT / "review" / f"{a.site}_{day}_lifecycle.png"
    out.parent.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=80)
    # summary numbers: half-width at half maximum above the side level, per age bin
    for i, ag in enumerate(ages):
        prof = Z[i] - side[i]
        if np.nanmax(prof) <= 2:
            continue
        half = np.nanmax(prof) / 2
        above = OFFSETS[prof >= half]
        print(f"age {ag:5.1f} min: peak {np.nanmax(Z[i]):5.1f} dBZ (+{np.nanmax(prof):4.1f} over sides), "
              f"FWHM {(above.max() - above.min() + 500) / 1e3:4.1f} km")
    print(out)


if __name__ == "__main__":
    main()
