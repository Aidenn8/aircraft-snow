"""Publication figure of the streak life cycle (from streak_lifecycle.py output).

usage: fig_lifecycle.py SITE YYYY-MM-DD --out FIG.png [--boot 2000]

(a) Reflectivity on the drifted track minus the mean 5-8 km to either side, against minutes
    since the aircraft passed, with a 90% bootstrap interval over flights. The 10-50 min
    matching window is marked: inside it the composite is shaped by selection (each flight
    was chosen because a streak lay on its track then), outside it is not.
(b) Mean cross-track profiles in four age ranges.
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from snow import figstyle, match

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROFILE_AGES = ((-30, -5), (15, 30), (30, 50), (55, 80))     # minutes


def contrast(per_flight, offsets):
    """Per-flight centre minus side (dBZ), shape (flights, ages)."""
    centre = per_flight[:, :, np.argmin(np.abs(offsets))]
    side = np.nanmean(per_flight[:, :, np.abs(offsets) >= 5000], axis=2)
    return centre - side


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("--out", required=True)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    figstyle.apply()
    L = np.load(ROOT / "data" / "events" / a.site / f"{a.day}_lifecycle.npz")
    off, ab, PF = L["offsets"], L["age_bins"], L["per_flight"]
    ages = 0.5 * (ab[:-1] + ab[1:]) / 60
    C = contrast(PF, off)
    n = C.shape[0]
    rng = np.random.default_rng(a.seed)
    boot = np.array([np.nanmean(C[rng.integers(0, n, n)], axis=0) for _ in range(a.boot)])
    mean = np.nanmean(C, axis=0)
    lo, hi = np.nanquantile(boot, [0.05, 0.95], axis=0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 3.9), layout="constrained",
                                   gridspec_kw={"width_ratios": [1.35, 1]})
    w0, w1 = match.AGE_MIN_S / 60, match.AGE_MAX_S / 60
    ax1.axvspan(w0, w1, color=figstyle.GRID, alpha=0.55, lw=0)
    ax1.text((w0 + w1) / 2, 0.97, "matching window\n(selected for echo)", transform=ax1.get_xaxis_transform(),
             ha="center", va="top", fontsize=8, color=figstyle.MUTED)
    ax1.axhline(0, color=figstyle.AXIS, lw=0.8)
    ax1.axvline(0, color=figstyle.AXIS, lw=0.8)
    ax1.fill_between(ages, lo, hi, color=figstyle.BLUE, alpha=0.18, lw=0)
    ax1.plot(ages, mean, color=figstyle.BLUE, lw=2)
    ax1.set_xlim(ages[0] - 2.5, ages[-1] + 2.5)
    ax1.set_xlabel("minutes since the aircraft passed")
    ax1.set_ylabel("track minus 5–8 km either side (dB)")
    ax1.set_title(f"(a) Echo on the drifted track, {n} seeding flights")
    ax1.text(ages[-1], hi[-1] + 0.4, "90% interval\nover flights", ha="right", va="bottom", fontsize=8,
             color=figstyle.INK_2)

    # before passage in grey (the baseline), then one hue darkening with age
    shades = [figstyle.MUTED, figstyle.BLUE_RAMP[3], figstyle.BLUE_RAMP[5], figstyle.BLUE_RAMP[7]]
    assert len(shades) == len(PROFILE_AGES)
    for (g0, g1), col in zip(PROFILE_AGES, shades):
        m = (ages >= g0) & (ages <= g1)
        prof = np.nanmean(np.nanmean(PF[:, m, :], axis=1), axis=0)
        ax2.plot(off / 1e3, prof, color=col, lw=2, label=f"{g0} to {g1} min")
    ax2.set_xlabel("distance across the drifted track (km)")
    ax2.set_ylabel("mean reflectivity (dBZ)")
    ax2.set_title("(b) Cross-track profile by age")
    ax2.legend(loc="upper right", title="since passage", title_fontsize=8)
    fig.savefig(a.out)
    print(a.out)
    for i, ag in enumerate(ages):
        print(f"{ag:6.1f} min: {mean[i]:5.1f} dB [{lo[i]:5.1f}, {hi[i]:5.1f}]")


if __name__ == "__main__":
    main()
