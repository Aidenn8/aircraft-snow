"""Census figure: where and how often, and at what temperatures, aircraft made snow.

usage: fig_census.py [--census data/census] --out FIG.png

(a) Share of radar-days with an event day, per radar and winter (from rates.csv).
(b) HRRR temperature at the seeding points of every catalogued seeding flight, by engine
    class (from seeders.csv): each flight a dot, with the median marked.
(c) Implied fall speed of the seeding flights' ice (height above the beam / delay).
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from snow import figstyle

ROOT = pathlib.Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", default=str(ROOT / "data" / "census"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    figstyle.apply()
    cen = pathlib.Path(a.census)
    rates = pd.read_csv(cen / "rates.csv")
    S = pd.read_csv(cen / "seeders.csv") if (cen / "seeders.csv").exists() else pd.DataFrame()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), layout="constrained",
                             gridspec_kw={"width_ratios": [1.5, 1, 1]})

    ax = axes[0]
    winters = sorted(rates.winter.unique())
    sites = (rates.groupby("site").event_days.sum().sort_values(ascending=False).index.tolist())
    width = 0.8 / max(len(winters), 1)
    for i, w in enumerate(winters):
        r = rates[rates.winter == w].set_index("site").reindex(sites)
        x = np.arange(len(sites)) + (i - (len(winters) - 1) / 2) * width
        ax.bar(x, 100 * r.event_days / r.radar_days, width * 0.92, color=figstyle.SERIES[i], label=w)
    ax.set_xticks(np.arange(len(sites)), sites)
    ax.set_ylabel("radar-days with aircraft snow (%)")
    ax.set_title("(a) Event days by radar and winter")
    ax.legend(title="winter (Nov-Mar)", title_fontsize=8)
    ax.grid(axis="x", visible=False)

    ax = axes[1]
    if not S.empty:
        classes = [c for c in ("jet", "turboprop", "piston") if (S.engine_class == c).any()]
        rng = np.random.default_rng(0)
        for i, c in enumerate(classes):
            v = S[S.engine_class == c].seed_temp_c.dropna().to_numpy()
            ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=14, color=figstyle.SERIES[i], alpha=0.7, lw=0)
            ax.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color=figstyle.INK, lw=2)
            ax.annotate(f"n={len(v)}", (i, np.nanmax(v) if len(v) else 0), xytext=(0, 6),
                        textcoords="offset points", ha="center", fontsize=8, color=figstyle.MUTED)
        ax.set_xticks(range(len(classes)), [f"{c}s" for c in classes])
        ax.invert_yaxis()
    ax.set_ylabel("HRRR temperature at the seeding points (°C)")
    ax.set_title("(b) Seeding temperature")
    ax.grid(axis="x", visible=False)

    ax = axes[2]
    if not S.empty and "v_fall" in S:
        v = S.v_fall.dropna().to_numpy()
        ax.hist(v, bins=np.arange(0, 2.05, 0.1), color=figstyle.BLUE, edgecolor=figstyle.SURFACE, lw=1)
        ax.axvline(np.median(v), color=figstyle.INK, lw=1.2)
        ax.annotate(f"median {np.median(v):.2f} m/s", (np.median(v), ax.get_ylim()[1] * 0.92), xytext=(5, 0),
                    textcoords="offset points", fontsize=8.5, color=figstyle.INK_2)
    ax.set_xlabel("implied fall speed (m/s)")
    ax.set_ylabel("seeding flights")
    ax.set_title("(c) How fast their ice fell")
    fig.savefig(a.out)
    print(a.out)


if __name__ == "__main__":
    main()
