"""Calibration figure: how often a real streak is recovered, and how often chance passes.

usage: fig_calibration.py --inj A.parquet [B.parquet ...] --census data/census
       [--null-census data/census_null] --out FIG.png

(a) Injection-recovery: share of synthetic streaks (injected along real flights' drifted
    tracks into real radar on quiet days) found as a streak object, matched to some flight
    at p <= 0.05 (event recall), and matched to the injected flight (attribution), by peak
    reflectivity, with 95% Clopper-Pearson intervals.
(b) Day-test calibration: for nominal levels, the share of p-values at or below the level,
    for the pooled leave-one-out rotation nulls of the census (does the gamma tail mean what
    it says?) and for real flights on summer null days (are real and rotated flights
    exchangeable when no seeding is possible?). Calibrated tests lie on the diagonal.
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import beta

from snow import figstyle


def clopper_pearson(k, n, conf=0.95):
    a = (1 - conf) / 2
    lo = beta.ppf(a, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - a, k + 1, n - k) if k < n else 1.0
    return lo, hi


def share_below(p, levels):
    p = np.asarray(p, float)
    p = p[np.isfinite(p)]
    return np.array([np.mean(p <= L) if len(p) else np.nan for L in levels]), len(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inj", nargs="+", required=True)
    ap.add_argument("--census", default="")
    ap.add_argument("--null-census", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    figstyle.apply()
    df = pd.concat([pd.read_parquet(p) for p in a.inj], ignore_index=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 3.9), layout="constrained")

    peaks = np.sort(df.peak_dbz.unique())
    series = [("object_found", "found as a streak object", figstyle.MUTED, -0.4)]
    if "recovered_any" in df:
        series.append(("recovered_any", "matched to a flight, p ≤ 0.05", figstyle.ORANGE, 0.0))
    series.append(("recovered", "matched to the injected flight", figstyle.BLUE, 0.4))
    for col, name, colour, dx in series:
        rows = []
        for pk in peaks:
            d = df[df.peak_dbz == pk]
            k, n = int(d[col].sum()), len(d)
            lo, hi = clopper_pearson(k, n)
            rows.append((pk, k / n, lo, hi, n))
        r = np.array(rows)
        ax1.errorbar(r[:, 0] + dx, r[:, 1], yerr=[r[:, 1] - r[:, 2], r[:, 3] - r[:, 1]], fmt="o-",
                     color=colour, ms=6, lw=1.6, capsize=0, elinewidth=1.2, label=name)
    for pk in peaks:
        n = int((df.peak_dbz == pk).sum())
        ax1.text(pk, 0.02, f"n={n}", ha="center", fontsize=7.5, color=figstyle.MUTED)
    ax1.set_ylim(0, 1.05)
    ax1.set_xticks(peaks)
    ax1.set_xlabel("peak reflectivity of the injected streak (dBZ)")
    ax1.set_ylabel("share of injected streaks")
    ax1.set_title(f"(a) Recovery of {len(df)} synthetic streaks")
    ax1.legend(loc="center right")

    levels = np.array([0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001])
    ax2.plot(levels, levels, color=figstyle.AXIS, lw=1, ls="--")
    ax2.text(0.0012, 0.0012 * 1.6, "calibrated", rotation=45, fontsize=7.5, color=figstyle.MUTED,
             rotation_mode="anchor")
    drawn = False
    if a.census and (pathlib.Path(a.census) / "null_pvalues.npy").exists():
        null = np.load(pathlib.Path(a.census) / "null_pvalues.npy")
        s, n = share_below(null, levels)
        ok = s > 0
        ax2.plot(levels[ok], s[ok], "o-", color=figstyle.BLUE, ms=5, label=f"rotations scored as real (n={n:,})")
        drawn = True
    if a.null_census and (pathlib.Path(a.null_census) / "radar_days.csv").exists():
        rd = pd.read_csv(pathlib.Path(a.null_census) / "radar_days.csv")
        s, n = share_below(rd.p_gamma, levels)
        lo_hi = np.array([clopper_pearson(int(round(v * n)), n) for v in s])
        ok = s > 0
        ax2.errorbar(levels[ok], s[ok], yerr=[s[ok] - lo_hi[ok, 0], lo_hi[ok, 1] - s[ok]], fmt="s",
                     color=figstyle.ORANGE, ms=5, capsize=0, elinewidth=1.2,
                     label=f"real flights, summer null days (n={n})")
        drawn = True
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("nominal p-value level")
    ax2.set_ylabel("share at or below the level")
    ax2.set_title("(b) Day test under no seeding")
    if drawn:
        ax2.legend(loc="upper left")
    fig.savefig(a.out)
    print(a.out)


if __name__ == "__main__":
    main()
