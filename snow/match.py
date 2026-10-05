"""Match radar streak objects to the single flight whose drifted track they follow.

For a streak at scan time T with skeleton (centre-line) pixels s_k, and a flight f whose
points (x_i, y_i, t_i; tracks interpolated to ~250 m spacing) drift with (u_i + du, v_i + dv)
to
    p_i(T) = (x_i + (u_i + du)(T - t_i), y_i + (v_i + dv)(T - t_i)),   T - t_i in [AGE_MIN_S, AGE_MAX_S]
two quantities measure how well the flight explains the streak:
    coverage   fraction of skeleton pixels within TOL_M of some p_i(T)
    precision  fraction of the flight's drifted points near the streak (within NEAR_M of
               the skeleton) that are within TOL_M of it; a flight that merely crosses the
               streak has low precision
and the match score is their harmonic mean, F1 = 2 C P / (C + P). (du, dv) ranges over
DRIFT_PERTURB_MS to allow for HRRR wind error; the best perturbation is kept.

Newness: a streak made by flight f did not exist before f passed. The skeleton is moved
back with the matched drift to every scan between NEW_WINDOW_S[0] and NEW_WINDOW_S[1]
seconds relative to f's passage (the median time of its matching points); the fraction of
those samples below QUIET_DBZ must be at least NEW_MIN_QUIET, over at least NEW_MIN_SCANS
scans. Among flights with coverage and precision >= MIN_COVERAGE the best new one (highest
F1) is the match. (In a stream of aircraft on one path the streak predates the later
aircraft, so the match is normally the first one through new cloud.)

Chance: the whole procedure is repeated for placebo flight sets, every flight rotated about
the radar by an angle in PLACEBO_ANGLES (same ranges, shapes, timing and wind, wrong place).
The p-value is (1 + number of rotations whose best new F1 >= the real one) /
(1 + number of rotations).
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from skimage import morphology

# Streaks are first visible 15-40 min after passage (Denver case; Iwasaki & Fujii 2025 find
# hole-punch lines within 10-20 min), so only aircraft that passed 10-50 min before the scan
# are candidates; older streaks were already candidates when younger.
AGE_MIN_S = 10 * 60
AGE_MAX_S = 50 * 60
TOL_M = 1000.0
NEAR_M = 3000.0
DRIFT_PERTURB_MS = (-1.0, -0.5, 0.0, 0.5, 1.0)
PLACEBO_ANGLES = tuple(np.arange(20.0, 340.1, 3.2))     # 101 rotations
MIN_COVERAGE = 0.5
MAX_CANDIDATES = 6
NEW_WINDOW_S = (-30 * 60, -3 * 60)
QUIET_DBZ = 0.0
NEW_MIN_QUIET = 0.7
NEW_MIN_SCANS = 2


def skeleton_xy(label_img: np.ndarray, label: int, xs: np.ndarray, ys: np.ndarray):
    sk = morphology.skeletonize(label_img == label)
    iy, ix = np.nonzero(sk)
    return xs[ix], ys[iy]


def _candidates(sx, sy, x, y, t, T, drift, flight, rot=0.0, tree=None):
    """Flights explaining the skeleton, best first: list of (F1, flight, du, dv, t_pass,
    u_med, v_med, coverage, precision) for the best drift perturbation of each flight."""
    from scipy.spatial import cKDTree
    if rot:
        a = np.radians(rot)
        x, y = x * np.cos(a) - y * np.sin(a), x * np.sin(a) + y * np.cos(a)
    tree = tree or cKDTree(np.c_[sx, sy])
    n_sk = len(sx)
    age = T - t
    best = {}
    lo_x, hi_x = sx.min() - NEAR_M, sx.max() + NEAR_M
    lo_y, hi_y = sy.min() - NEAR_M, sy.max() + NEAR_M
    for du, dv in itertools.product(DRIFT_PERTURB_MS, DRIFT_PERTURB_MS):
        px = x + (drift[:, 0] + du) * age
        py = y + (drift[:, 1] + dv) * age
        box = (px > lo_x) & (px < hi_x) & (py > lo_y) & (py < hi_y)
        if not box.any():
            continue
        idx = np.nonzero(box)[0]
        d, nearest = tree.query(np.c_[px[idx], py[idx]], distance_upper_bound=NEAR_M)
        near = np.isfinite(d)
        if not near.any():
            continue
        idx, d, nearest = idx[near], d[near], nearest[near]
        fl = flight[idx]
        hit = d <= TOL_M
        for f in np.unique(fl[hit]):
            mf = fl == f
            prec = float(hit[mf].mean())
            if prec < MIN_COVERAGE:
                continue
            # skeleton pixels within TOL of any of the flight's points
            sel = idx[mf & hit]
            covered = tree.query_ball_point(np.c_[px[sel], py[sel]], r=TOL_M)
            cov = len(set().union(*covered)) / n_sk
            if cov < MIN_COVERAGE:
                continue
            f1 = 2 * cov * prec / (cov + prec)
            if f1 > best.get(f, (0,))[0]:
                close = idx[mf & hit]
                best[f] = (f1, int(f), du, dv, float(np.median(t[close])),
                           float(np.median(drift[close, 0]) + du), float(np.median(drift[close, 1]) + dv),
                           cov, prec)
    return sorted(best.values(), key=lambda c: -c[0])[:MAX_CANDIDATES]


def newness(grids, sx, sy, T, t_pass, u, v):
    """Fraction of back-drifted skeleton samples below QUIET_DBZ in scans between
    t_pass + NEW_WINDOW_S[0] and t_pass + NEW_WINDOW_S[1], and the number of scans used."""
    sel = np.nonzero((grids.t >= t_pass + NEW_WINDOW_S[0]) & (grids.t <= t_pass + NEW_WINDOW_S[1]))[0]
    vals = []
    for k in sel:
        dtk = T - grids.t[k]
        z = grids.lookup(k, sx - u * dtk, sy - v * dtk)
        vals.append(z[np.isfinite(z)])
    if not vals:
        return np.nan, 0
    allv = np.concatenate(vals)
    if allv.size == 0:
        return np.nan, 0
    return float(np.mean(allv < QUIET_DBZ)), len(sel)


def best_new_match(grids, sx, sy, x, y, t, T, drift, flight, rot=0.0, tree=None):
    """Best candidate flight that also passes the newness test: (F1, flight, du, dv, t_pass,
    quiet fraction, coverage, precision), or zeros if none."""
    for f1, f, du, dv, t_pass, u, v, cov, prec in _candidates(sx, sy, x, y, t, T, drift, flight, rot, tree):
        q, n = newness(grids, sx, sy, T, t_pass, u, v)
        if n >= NEW_MIN_SCANS and q >= NEW_MIN_QUIET:
            return f1, f, du, dv, t_pass, q, cov, prec
    return 0.0, -1, np.nan, np.nan, np.nan, np.nan, 0.0, 0.0


def match_streaks(grids, T: float, label_img, objects: pd.DataFrame, xs, ys,
                  points: pd.DataFrame, drift: np.ndarray, angles=PLACEBO_ANGLES) -> pd.DataFrame:
    """For every streak object: the best real flight that passes the newness test, its F1
    score (column `coverage` for compatibility, plus `cov` and `precision`), and the
    rotation p-value of that score."""
    from scipy.spatial import cKDTree
    age_all = T - points.t.to_numpy(float)
    m = (age_all >= AGE_MIN_S) & (age_all <= AGE_MAX_S) & np.isfinite(drift[:, 0])
    x, y = points.x.to_numpy()[m], points.y.to_numpy()[m]
    t, dr, fl = points.t.to_numpy(float)[m], drift[m], points.flight.to_numpy()[m]
    rows = []
    for _, o in objects.iterrows():
        sx, sy = skeleton_xy(label_img, o.label, xs, ys)
        if len(sx) < 4:
            continue
        tree = cKDTree(np.c_[sx, sy])
        f1, f, du, dv, t_pass, q, cov, prec = best_new_match(grids, sx, sy, x, y, t, T, dr, fl, tree=tree)
        plac = np.array([best_new_match(grids, sx, sy, x, y, t, T, dr, fl, rot=r, tree=tree)[0]
                         for r in angles])
        rows.append({**o.to_dict(), "flight": f, "coverage": f1, "cov": cov, "precision": prec,
                     "du": du, "dv": dv, "t_pass": t_pass, "age_s": T - t_pass if f >= 0 else np.nan,
                     "quiet_before": q, "placebo_mean": float(plac.mean()),
                     "placebo_q95": float(np.quantile(plac, 0.95)),
                     "placebo_frac_matched": float(np.mean(plac > 0)),
                     "p_value": float((1 + np.sum(plac >= f1)) / (1 + len(plac))) if f1 > 0 else 1.0,
                     "placebo_cov": plac.astype(np.float32)})
    return pd.DataFrame(rows)


def day_statistic(table: pd.DataFrame) -> dict:
    """Combine all streaks of a radar-day: E = sum of -log10 p over streaks.

    Null: under 'no aircraft effect' the real flight set is exchangeable with the rotated
    ones, so each rotation in turn is treated as the real set and scored against the other
    rotations plus the real set; its E is computed the same way. The day's p-value is
    (1 + number of rotations with E >= real E) / (1 + number of rotations)."""
    if table.empty:
        return {"E": 0.0, "p_day": 1.0, "n_streaks": 0, "n_matched": 0}
    C0 = table.coverage.to_numpy(float)                      # (n_streaks,)
    CP = np.stack(table.placebo_cov.to_numpy()).astype(float)  # (n_streaks, R)
    R = CP.shape[1]
    allc = np.concatenate([C0[:, None], CP], axis=1)          # column 0 = real

    def pvals(col):
        c = allc[:, col]
        ge = (allc >= c[:, None]).sum(axis=1) - 1              # others with coverage >= this
        p = (1 + ge) / (1 + R)
        return np.where(c > 0, p, 1.0)

    E = lambda col: float(-np.log10(pvals(col)).sum())
    E_real = E(0)
    E_null = np.array([E(j) for j in range(1, R + 1)])
    return {"E": E_real, "p_day": float((1 + np.sum(E_null >= E_real)) / (1 + R)),
            "E_null_median": float(np.median(E_null)), "E_null_q95": float(np.quantile(E_null, 0.95)),
            "n_streaks": int(len(table)), "n_matched": int((C0 > 0).sum())}
