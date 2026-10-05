"""Associate radar streak objects with aircraft tracks drifted by the wind.

For a scan at time T, every track point passed at t_i and drifts to
    x_i + (u_i + du)(T - t_i),  y_i + (v_i + dv)(T - t_i)
for ages T - t_i in [AGE_MIN_S, AGE_MAX_S], with (u_i, v_i) the HRRR layer wind between the
radar beam and the point's level (epoch.drift_velocity) and (du, dv) a small perturbation
from DRIFT_PERTURB_MS that allows for wind error. The drifted points of one altitude band
are rasterised on the streak grid and dilated by TOLERANCE_M; a streak's coverage is the
fraction of its pixels inside that footprint. The best coverage over bands and perturbations
is kept.

The same is done for placebo track sets: every point rotated about the radar by an angle in
PLACEBO_ROTATIONS_DEG (same ranges, same shapes, same traffic, same wind; wrong place). The
p-value of a streak's association is (1 + number of placebos with coverage >= real) /
(1 + number of placebos).
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy import ndimage
from skimage import morphology

AGE_MIN_S = 5 * 60
AGE_MAX_S = 120 * 60
TOLERANCE_M = 1500.0
DRIFT_PERTURB_MS = (-1.5, 0.0, 1.5)
PLACEBO_ROTATIONS_DEG = tuple(range(30, 331, 27))     # 12 rotations, none near 0
ALT_BANDS_FT = tuple((a, a + 3000) for a in range(1000, 22001, 1500))


def _footprint(px, py, xs, ys, tol_m):
    res = xs[1] - xs[0]
    ix = np.rint((px - xs[0]) / res).astype(int)
    iy = np.rint((py - ys[0]) / res).astype(int)
    ok = (ix >= 0) & (ix < len(xs)) & (iy >= 0) & (iy < len(ys))
    img = np.zeros((len(ys), len(xs)), bool)
    img[iy[ok], ix[ok]] = True
    r = int(round(tol_m / res))
    return ndimage.binary_dilation(img, structure=morphology.disk(r))


def coverage_table(T: float, label_img: np.ndarray, labels, xs, ys, points: pd.DataFrame,
                   drift: np.ndarray, rotate_deg: float = 0.0) -> pd.DataFrame:
    """Best coverage per streak label over altitude bands and drift perturbations.
    Returns columns label, coverage, band_lo, du, dv."""
    x0 = points.x.to_numpy()
    y0 = points.y.to_numpy()
    if rotate_deg:
        a = np.radians(rotate_deg)
        x0, y0 = x0 * np.cos(a) - y0 * np.sin(a), x0 * np.sin(a) + y0 * np.cos(a)
    age = T - points.t.to_numpy(float)
    alt = points.alt_ft.to_numpy()
    in_age = (age >= AGE_MIN_S) & (age <= AGE_MAX_S) & np.isfinite(drift[:, 0])
    best = {lab: (0.0, np.nan, np.nan, np.nan) for lab in labels}
    if not in_age.any() or not len(labels):
        return pd.DataFrame([(k, *v) for k, v in best.items()],
                            columns=["label", "coverage", "band_lo", "du", "dv"])
    pix = {lab: np.nonzero(label_img == lab) for lab in labels}
    for lo, hi in ALT_BANDS_FT:
        m = in_age & (alt >= lo) & (alt < hi)
        if m.sum() < 5:
            continue
        for du, dv in itertools.product(DRIFT_PERTURB_MS, DRIFT_PERTURB_MS):
            px = x0[m] + (drift[m, 0] + du) * age[m]
            py = y0[m] + (drift[m, 1] + dv) * age[m]
            fp = _footprint(px, py, xs, ys, TOLERANCE_M)
            for lab in labels:
                yy, xx = pix[lab]
                c = float(fp[yy, xx].mean())
                if c > best[lab][0]:
                    best[lab] = (c, lo, du, dv)
    return pd.DataFrame([(k, *v) for k, v in best.items()],
                        columns=["label", "coverage", "band_lo", "du", "dv"])


def associate(T: float, label_img, objects: pd.DataFrame, xs, ys, points, drift,
              rotations=PLACEBO_ROTATIONS_DEG) -> pd.DataFrame:
    """Real coverage, placebo coverages and the association p-value for each streak."""
    labels = list(objects.label)
    real = coverage_table(T, label_img, labels, xs, ys, points, drift)
    plac = [coverage_table(T, label_img, labels, xs, ys, points, drift, rotate_deg=r)
            .set_index("label").coverage.rename(f"rot{r}") for r in rotations]
    P = pd.concat(plac, axis=1)
    out = objects.merge(real, on="label")
    Pm = P.reindex(out.label).to_numpy()
    out["placebo_mean"] = Pm.mean(axis=1)
    out["placebo_max"] = Pm.max(axis=1)
    out["p_value"] = (1 + (Pm >= out.coverage.to_numpy()[:, None]).sum(axis=1)) / (1 + Pm.shape[1])
    return out


def contributing_flights(T: float, label_img, label: int, xs, ys, points, drift, band_lo,
                         du, dv) -> pd.DataFrame:
    """Flights whose drifted points (in the streak's best band and drift) fall on the streak:
    number of points inside, their median age and altitude."""
    age = T - points.t.to_numpy(float)
    alt = points.alt_ft.to_numpy()
    m = ((age >= AGE_MIN_S) & (age <= AGE_MAX_S) & (alt >= band_lo) & (alt < band_lo + 3000)
         & np.isfinite(drift[:, 0]))
    px = points.x.to_numpy()[m] + (drift[m, 0] + du) * age[m]
    py = points.y.to_numpy()[m] + (drift[m, 1] + dv) * age[m]
    res = xs[1] - xs[0]
    ix = np.rint((px - xs[0]) / res).astype(int)
    iy = np.rint((py - ys[0]) / res).astype(int)
    ok = (ix >= 0) & (ix < len(xs)) & (iy >= 0) & (iy < len(ys))
    inside = np.zeros(m.sum(), bool)
    grown = ndimage.binary_dilation(label_img == label, structure=morphology.disk(int(TOLERANCE_M / res)))
    inside[ok] = grown[iy[ok], ix[ok]]
    sub = points[m][inside].assign(age=age[m][inside])
    if sub.empty:
        return pd.DataFrame(columns=["flight", "n_points", "age_med_s", "alt_med_ft"])
    return (sub.groupby("flight").agg(n_points=("t", "size"), age_med_s=("age", "median"),
                                      alt_med_ft=("alt_ft", "median")).reset_index())
