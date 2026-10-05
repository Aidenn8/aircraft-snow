"""How far does the best-fitting drift of each confident match lie from the HRRR drift?

usage: drift_check.py SITE YYYY-MM-DD [--alpha 0.005] [--max 3.0] [--step 0.25]

Matching allows the HRRR drift to be wrong by up to +-1 m/s per component. If the best
perturbation of many matches sits at that edge, the true drift error may be larger and
the window is clipping matches. For each streak matched at p <= alpha, the matched
flight's F1 against its streak is mapped over a wider grid of perturbations (same
coverage / precision definitions as match._candidates, single flight).
Writes data/events/SITE/YYYY-MM-DD_drift.csv.
"""
import argparse
import datetime as dt
import pathlib

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from snow import match, pipeline

ROOT = pathlib.Path(__file__).resolve().parents[1]


def f1_single(tree, n_sk, sx, sy, px, py):
    """F1 of one flight's drifted points against one skeleton (as in match._candidates)."""
    lo_x, hi_x = sx.min() - match.NEAR_M, sx.max() + match.NEAR_M
    lo_y, hi_y = sy.min() - match.NEAR_M, sy.max() + match.NEAR_M
    box = (px > lo_x) & (px < hi_x) & (py > lo_y) & (py < hi_y)
    if not box.any():
        return 0.0, 0.0, 0.0
    d, nearest = tree.query(np.c_[px[box], py[box]], distance_upper_bound=match.NEAR_M)
    near = np.isfinite(d)
    if not near.any():
        return 0.0, 0.0, 0.0
    hit = d[near] <= match.TOL_M
    prec = hit.sum() / near.sum()
    cov = np.unique(nearest[near][hit]).size / n_sk
    f1 = 2 * cov * prec / (cov + prec) if cov + prec > 0 else 0.0
    return f1, cov, prec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("--alpha", type=float, default=0.005)
    ap.add_argument("--max", type=float, default=3.0)
    ap.add_argument("--step", type=float, default=0.25)
    a = ap.parse_args()
    day = dt.date.fromisoformat(a.day)
    ev = ROOT / "data" / "events" / a.site
    S = pd.read_parquet(ev / f"{day}_streaks.parquet")
    conf = S[(S.coverage > 0) & (S.p_value <= a.alpha)]
    prep = pipeline.prepare(a.site, day)
    G = prep["grids"]
    pts, dr = pipeline.matching_points(prep)
    kept, _, xs, ys = pipeline.streak_objects_day(prep)
    labs = {float(G.t[k]): lab for k, _, lab in kept}
    grid = np.round(np.arange(-a.max, a.max + 1e-9, a.step), 3)
    rows = []
    for _, o in conf.iterrows():
        lab = labs[min(labs, key=lambda T: abs(T - o.t))]
        sx, sy = match.skeleton_xy(lab, int(o.label), xs, ys)
        tree = cKDTree(np.c_[sx, sy])
        T = float(o.t)
        age = T - pts.t.to_numpy(float)
        m = ((pts.flight.to_numpy() == o.flight) & (age >= match.AGE_MIN_S) & (age <= match.AGE_MAX_S)
             & match._physical_mask(pts, age))
        x, y, ag, d = pts.x.to_numpy()[m], pts.y.to_numpy()[m], age[m], dr[m]
        F = np.zeros((len(grid), len(grid)))
        for i, du in enumerate(grid):
            for j, dv in enumerate(grid):
                F[i, j] = f1_single(tree, len(sx), sx, sy, x + (d[:, 0] + du) * ag, y + (d[:, 1] + dv) * ag)[0]
        i, j = np.unravel_index(np.argmax(F), F.shape)
        # how sharp is the fit: the perturbations within 0.05 of the best F1. Along a straight
        # streak the drift component parallel to it is barely constrained, so this region is
        # an elongated ridge; the perpendicular width is what constrains the wind.
        good = F >= F.max() - 0.05
        gu, gv = np.nonzero(good)
        i0 = j0 = int(np.argmin(np.abs(grid)))
        inner = np.abs(grid) <= max(match.DRIFT_PERTURB_MS) + 1e-9
        f1_in = float(F[np.ix_(inner, inner)].max())
        # nearest perturbation (to HRRR) that reaches within 0.05 of the best F1
        r = np.hypot(grid[gu], grid[gv])
        rows.append({"t": T, "flight": int(o.flight), "f1_match": float(o.coverage), "du_match": o.du,
                     "dv_match": o.dv, "f1_hrrr": float(F[i0, j0]), "f1_window": f1_in,
                     "f1_best": float(F.max()), "gain_outside": float(F.max()) - f1_in,
                     "du_best": grid[i], "dv_best": grid[j], "r_needed": float(r.min()),
                     "du_lo": grid[gu.min()], "du_hi": grid[gu.max()],
                     "dv_lo": grid[gv.min()], "dv_hi": grid[gv.max()],
                     "u_hrrr": float(np.median(d[:, 0])), "v_hrrr": float(np.median(d[:, 1]))})
    out = pd.DataFrame(rows)
    out.to_csv(ev / f"{day}_drift.csv", index=False)
    pd.set_option("display.width", 200)
    print(out.round(2).to_string(index=False))
    print(f"median best perturbation: du {out.du_best.median():+.2f}, dv {out.dv_best.median():+.2f} m/s; "
          f"|best| > 1 m/s in {int(((out.du_best.abs() > 1) | (out.dv_best.abs() > 1)).sum())} of {len(out)}")
    print(f"F1 gained by allowing up to {a.max:.1f} m/s instead of {max(match.DRIFT_PERTURB_MS):.1f}: "
          f"median {out.gain_outside.median():.3f}, max {out.gain_outside.max():.3f}; "
          f"smallest perturbation within 0.05 of the best: median {out.r_needed.median():.2f} m/s, "
          f"> 1 m/s in {int((out.r_needed > 1).sum())} of {len(out)}")


if __name__ == "__main__":
    main()
