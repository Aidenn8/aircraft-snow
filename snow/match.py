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

Fall speed: ice made at height h above the radar beam (beam centre h_b, half-width w) is
seen by the beam after falling h - h_b, so a point of age a can only explain echo at the
beam if its implied fall speed lies in [V_FALL_MIN_MS, V_FALL_MAX_MS]:
    V_FALL_MIN * a <= (h - h_b) + w   and   (h - h_b) - w <= V_FALL_MAX * a.
The range spans snow crystals, aggregates and lightly rimed particles (Locatelli & Hobbs
1974). Points that fail it are not candidates. Rotation about the radar keeps range, height
and time, so real and rotated flights are filtered alike. (Tables without the height
columns dh, bw are matched without this filter.)

Temperature: only points where the HRRR analysis is at or below T_SEED_MAX_C are candidates
(supercooled cloud is needed; the bound is far warmer than any reported seeding threshold,
so it removes only points that cannot seed). Null-mode runs on warm days switch it off.

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
V_FALL_MIN_MS = 0.2
V_FALL_MAX_MS = 2.0
BEAM_HALF_WIDTH_DEG = 0.475       # half the 0.95 deg half-power beamwidth
T_SEED_MAX_C = -4.0


def fall_ok(dh, bw, age):
    """Whether ice made dh metres above the beam centre (beam half-width bw) can be at the
    beam after |age| seconds, for a fall speed in [V_FALL_MIN_MS, V_FALL_MAX_MS]. The
    absolute age keeps the time-reversed diagnostic a mirror image of the forward test."""
    a = np.abs(np.asarray(age, float))
    dh = np.asarray(dh, float)
    bw = np.asarray(bw, float)
    return (dh + bw >= V_FALL_MIN_MS * a) & (dh - bw <= V_FALL_MAX_MS * a)


def _physical_mask(points: pd.DataFrame, age: np.ndarray, max_temp_c=T_SEED_MAX_C) -> np.ndarray:
    """Points that could have made echo seen at the beam at this age: fall speed in range,
    and (unless max_temp_c is None) cold enough for supercooled cloud."""
    m = np.ones(len(points), bool)
    if "dh" in points and "bw" in points:
        m &= fall_ok(points.dh.to_numpy(float), points.bw.to_numpy(float), age)
    if max_temp_c is not None and "temp_c" in points:
        m &= points.temp_c.to_numpy(float) <= max_temp_c
    return m


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
    # conservative prefilter: no perturbation moves a point by more than max|d| * sqrt(2) * age
    bx = x + drift[:, 0] * age
    by = y + drift[:, 1] * age
    reach = max(abs(p) for p in DRIFT_PERTURB_MS) * np.sqrt(2) * np.abs(age)
    pre = (bx > lo_x - reach) & (bx < hi_x + reach) & (by > lo_y - reach) & (by < hi_y + reach)
    if not pre.any():
        return []
    x, y, t, drift, flight, age = x[pre], y[pre], t[pre], drift[pre], flight[pre], age[pre]
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
        if not hit.any():
            continue
        # vectorised over flights: precision = hits / near points; coverage = distinct
        # skeleton pixels that are some hit point's nearest pixel (tracks are interpolated
        # to 250 m and the skeleton has 500 m pixels, so a track running along the streak is
        # nearest to every pixel it covers)
        uf, fi = np.unique(fl, return_inverse=True)
        n_near = np.bincount(fi, minlength=len(uf))
        n_hit = np.bincount(fi[hit], minlength=len(uf))
        covered = np.zeros(len(uf) * n_sk, bool)          # (flight, skeleton pixel) pairs
        covered[fi[hit] * n_sk + nearest[hit]] = True
        n_cov = covered.reshape(len(uf), n_sk).sum(axis=1)
        prec_all = n_hit / np.maximum(n_near, 1)
        cov_all = n_cov / n_sk
        for j in np.nonzero((prec_all >= MIN_COVERAGE) & (cov_all >= MIN_COVERAGE))[0]:
            f = uf[j]
            prec, cov = float(prec_all[j]), float(cov_all[j])
            mf = fl == f
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


def best_new_match(grids, sx, sy, x, y, t, T, drift, flight, rot=0.0, tree=None,
                   require_new=True):
    """Best candidate flight that also passes the newness test (unless require_new is
    False): (F1, flight, du, dv, t_pass, quiet fraction, coverage, precision), or zeros."""
    for f1, f, du, dv, t_pass, u, v, cov, prec in _candidates(sx, sy, x, y, t, T, drift, flight, rot, tree):
        if not require_new:
            return f1, f, du, dv, t_pass, np.nan, cov, prec
        q, n = newness(grids, sx, sy, T, t_pass, u, v)
        if n >= NEW_MIN_SCANS and q >= NEW_MIN_QUIET:
            return f1, f, du, dv, t_pass, q, cov, prec
    return 0.0, -1, np.nan, np.nan, np.nan, np.nan, 0.0, 0.0


def age_window(reverse: bool = False):
    """Ages (scan time minus passage time) of candidate points. Reverse: aircraft that will
    pass AGE_MIN-AGE_MAX after the scan (a diagnostic; they cannot have made the echo)."""
    return (-AGE_MAX_S, -AGE_MIN_S) if reverse else (AGE_MIN_S, AGE_MAX_S)


def matched_points(sx, sy, points: pd.DataFrame, drift: np.ndarray, T: float, flight: int,
                   du: float, dv: float, reverse: bool = False, max_temp_c=T_SEED_MAX_C) -> pd.DataFrame:
    """The matched flight's (interpolated) points whose drifted position at T lies within
    TOL_M of the streak's skeleton: these are where the aircraft seeded."""
    from scipy.spatial import cKDTree
    f = (points.flight.to_numpy() == flight)
    age = T - points.t.to_numpy(float)
    a0, a1 = age_window(reverse)
    m = f & (age >= a0) & (age <= a1) & _physical_mask(points, age, max_temp_c)
    p = points[m]
    px = p.x.to_numpy() + (drift[m, 0] + du) * age[m]
    py = p.y.to_numpy() + (drift[m, 1] + dv) * age[m]
    d, _ = cKDTree(np.c_[sx, sy]).query(np.c_[px, py], distance_upper_bound=TOL_M)
    return p[np.isfinite(d)]


def match_streaks(grids, T: float, label_img, objects: pd.DataFrame, xs, ys,
                  points: pd.DataFrame, drift: np.ndarray, angles=PLACEBO_ANGLES,
                  reverse: bool = False, max_temp_c=T_SEED_MAX_C) -> pd.DataFrame:
    """For every streak object: the best real flight that passes the newness test, its F1
    score (column `coverage` for compatibility, plus `cov` and `precision`), and the
    rotation p-value of that score. reverse=True runs the time-reversed diagnostic: aircraft
    passing 10-50 min after the scan, no newness test."""
    from scipy.spatial import cKDTree
    age_all = T - points.t.to_numpy(float)
    a0, a1 = age_window(reverse)
    m = ((age_all >= a0) & (age_all <= a1) & np.isfinite(drift[:, 0])
         & _physical_mask(points, age_all, max_temp_c))
    x, y = points.x.to_numpy()[m], points.y.to_numpy()[m]
    t, dr, fl = points.t.to_numpy(float)[m], drift[m], points.flight.to_numpy()[m]
    rows = []
    for _, o in objects.iterrows():
        sx, sy = skeleton_xy(label_img, o.label, xs, ys)
        if len(sx) < 4:
            continue
        tree = cKDTree(np.c_[sx, sy])
        f1, f, du, dv, t_pass, q, cov, prec = best_new_match(grids, sx, sy, x, y, t, T, dr, fl, tree=tree,
                                                             require_new=not reverse)
        pm = [best_new_match(grids, sx, sy, x, y, t, T, dr, fl, rot=r, tree=tree, require_new=not reverse)
              for r in angles]
        plac = np.array([m_[0] for m_ in pm])
        plac_flight = np.array([m_[1] for m_ in pm], np.int64)
        alt_med, p_med, n_seed, v_fall = np.nan, np.nan, 0, np.nan
        if f >= 0:
            mp = matched_points(sx, sy, points, drift, T, f, du, dv, reverse, max_temp_c)
            if len(mp):
                alt_med, p_med, n_seed = float(mp.alt_ft.median()), float(mp.p_hpa.median()), len(mp)
                if "dh" in mp:
                    v_fall = float(np.median(mp.dh.to_numpy(float) / np.abs(T - mp.t.to_numpy(float))))
        rows.append({**o.to_dict(), "flight": f, "coverage": f1, "cov": cov, "precision": prec,
                     "seed_alt_ft": alt_med, "seed_p_hpa": p_med, "n_seed_points": n_seed,
                     "v_fall": v_fall,
                     "du": du, "dv": dv, "t_pass": t_pass, "age_s": T - t_pass if f >= 0 else np.nan,
                     "quiet_before": q, "placebo_mean": float(plac.mean()),
                     "placebo_q95": float(np.quantile(plac, 0.95)),
                     "placebo_frac_matched": float(np.mean(plac > 0)),
                     "p_value": float((1 + np.sum(plac >= f1)) / (1 + len(plac))) if f1 > 0 else 1.0,
                     "placebo_cov": plac.astype(np.float32), "placebo_flight": plac_flight})
    return pd.DataFrame(rows)


def pseudo_pvalues(table: pd.DataFrame) -> np.ndarray:
    """p-values of every streak for the real set (column 0) and for each rotation treated
    as if it were real (columns 1..R), each against all the other sets. Shape (n, R + 1)."""
    C0 = table.coverage.to_numpy(float)
    CP = np.stack(table.placebo_cov.to_numpy()).astype(float)
    allc = np.concatenate([C0[:, None], CP], axis=1)
    R = CP.shape[1]
    out = np.ones_like(allc)
    for col in range(R + 1):
        c = allc[:, col]
        ge = (allc >= c[:, None]).sum(axis=1) - 1
        out[:, col] = np.where(c > 0, (1 + ge) / (1 + R), 1.0)
    return out


def flight_matrix(table: pd.DataFrame) -> np.ndarray:
    """Matched flight per streak for the real set (column 0) and each rotation (1..R);
    -1 where nothing matched."""
    F0 = table.flight.to_numpy(np.int64)
    if "placebo_flight" in table:
        FP = np.stack(table.placebo_flight.to_numpy()).astype(np.int64)
    else:     # older tables: treat every placebo match as a distinct flight
        CP = np.stack(table.placebo_cov.to_numpy())
        FP = np.where(CP > 0, -2 - np.arange(CP.size).reshape(CP.shape), -1)
    return np.concatenate([F0[:, None], FP], axis=1)


def per_flight_min_p(P: np.ndarray, F: np.ndarray, col: int) -> np.ndarray:
    """For one set (column), the smallest p over the streaks matched to each flight."""
    f = F[:, col]
    m = f != -1
    if not m.any():
        return np.array([])
    s = pd.Series(P[m, col]).groupby(f[m]).min()
    return s.to_numpy()


def empirical_fdr(table: pd.DataFrame, thresholds=(0.005, 0.01, 0.02, 0.05)) -> pd.DataFrame:
    """For each p-value threshold: the number of real flights with a streak match at or
    below it (each flight counted once, by its best streak), the mean number the rotations
    produce when treated as real (false matches expected under 'no aircraft effect'), the
    implied false discovery rate and the estimated number of genuinely seeding flights."""
    P = pseudo_pvalues(table)
    F = flight_matrix(table)
    R = P.shape[1] - 1
    mins = [per_flight_min_p(P, F, j) for j in range(R + 1)]
    rows = []
    for tau in thresholds:
        n_real = int((mins[0] <= tau).sum())
        null = np.array([(m <= tau).sum() for m in mins[1:]])
        exp_false = float(null.mean())
        rows.append({"threshold": tau, "n_real": n_real, "expected_false": exp_false,
                     "null_q95": float(np.quantile(null, 0.95)),
                     "fdr": exp_false / n_real if n_real else np.nan,
                     "excess": n_real - exp_false})
    return pd.DataFrame(rows)


def tail_pvalues(E_real: float, E_null: np.ndarray) -> dict:
    """Smooth tail estimates of P(E_null >= E_real) beyond the 1/(R+1) resolution of the rank
    p-value: a z-score under a normal approximation, and a gamma distribution fitted to the
    null scores (E is a non-negative sum of evidence, typically right-skewed). Both are
    approximations whose calibration is checked on null days."""
    from scipy import stats
    E_null = np.asarray(E_null, float)
    mu, sd = float(E_null.mean()), float(E_null.std(ddof=1)) if len(E_null) > 1 else 0.0
    out = {"z_day": (E_real - mu) / sd if sd > 0 else np.nan}
    out["p_day_z"] = float(stats.norm.sf(out["z_day"])) if sd > 0 else np.nan
    pos = E_null[E_null > 0]
    if len(pos) >= 10 and np.std(pos) > 0:
        # moment-matched gamma for the positive part, mixed with the point mass at zero
        frac0 = 1 - len(pos) / len(E_null)
        k = pos.mean() ** 2 / pos.var(ddof=1)
        theta = pos.var(ddof=1) / pos.mean()
        tail = stats.gamma.sf(E_real, k, scale=theta) if E_real > 0 else 1.0
        out["p_day_gamma"] = float((1 - frac0) * tail + (frac0 if E_real <= 0 else 0.0))
    else:
        out["p_day_gamma"] = np.nan
    return out


def day_statistic(table: pd.DataFrame) -> dict:
    """Combine a radar-day's evidence: E = sum over matched flights of -log10 of the
    flight's best streak p-value (a flight matched to the same echo in many scans counts
    once).

    Null: under 'no aircraft effect' the real flight set is exchangeable with the rotated
    ones, so each rotation in turn is treated as the real set and scored against the other
    rotations plus the real set; its E is computed the same way. The day's p-value is
    (1 + number of rotations with E >= real E) / (1 + number of rotations)."""
    if table.empty:
        return {"E": 0.0, "p_day": 1.0, "n_streaks": 0, "n_matched": 0}
    C0 = table.coverage.to_numpy(float)
    P = pseudo_pvalues(table)                     # column 0 real, 1..R rotations as if real
    F = flight_matrix(table)
    R = P.shape[1] - 1
    E_all = np.array([-np.log10(per_flight_min_p(P, F, j)).sum() for j in range(R + 1)])
    E_real, E_null = float(E_all[0]), E_all[1:]
    return {"E": E_real, "p_day": float((1 + np.sum(E_null >= E_real)) / (1 + R)),
            **tail_pvalues(E_real, E_null),
            "n_flights_matched": int(len(per_flight_min_p(P, F, 0))),
            "E_null_median": float(np.median(E_null)), "E_null_q95": float(np.quantile(E_null, 0.95)),
            # kept so that the tail estimates can be calibrated across days: each rotation,
            # scored against the others, gives a null p-value (see census_summary.py)
            "E_null": [round(float(e), 4) for e in E_null],
            "n_streaks": int(len(table)), "n_matched": int((C0 > 0).sum())}
