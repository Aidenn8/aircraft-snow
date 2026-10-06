"""Census figure: how often, in what clouds, and at what temperatures aircraft made snow.

usage: fig_census.py [--census data/census] --out FIG.png

(a) Per radar, the radar-days whose calibrated p-value is at or below 0.01, against the
    number expected with no seeding at all (1% of the tested days)       (rates.csv)
(b) Share of tested radar-days at p <= 0.05 by hours of supercooled cloud water over the
    radar that day, with 90% intervals; 5% is what no seeding gives         (slw.csv)
(c) Jets crossing the seeding layer on catalogue days: share that left a streak, by HRRR
    temperature, arrivals and departures, with 90% intervals; the open markers are the
    chance level from rotated tracks                                        (efficiency.csv)
(d) Seeding temperature of every catalogue flight, by engine class, against the thresholds
    usually quoted (jets about -20 to -25 C, turboprops about -10 to -15 C) (seeders.csv)
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
SITE_ORDER = ["KFTG", "KMPX", "KLOT", "KIND", "KDTX", "KCLE", "KBUF", "KPBZ"]


def tbin_mid(s: str) -> float:
    lo, hi = (float(v) for v in s.strip("(]").split(","))
    return (lo + hi) / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", default=str(ROOT / "data" / "census"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    figstyle.apply()
    cen = pathlib.Path(a.census)
    names = pd.read_csv(ROOT / "data" / "sites.csv").drop_duplicates("radar").set_index("radar").municipality
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6), layout="constrained")

    # (a) how often, per radar
    ax = axes[0, 0]
    r = pd.read_csv(cen / "rates.csv")
    r = r[(r.score == "census") & (r.winter == "all") & (r.site != "all")].set_index("site")
    order = [s for s in SITE_ORDER if s in r.index]
    y = np.arange(len(order))[::-1]
    ax.barh(y, r.loc[order, "expected_p01"], height=0.55, color=figstyle.GRID, label="expected with no seeding")
    ax.plot(r.loc[order, "p01"], y, "o", color=figstyle.BLUE, ms=8, label="observed")
    ax.set_yticks(y, [f"{names.get(s, s)} ({s})" for s in order])
    ax.set_xlabel("radar-days with p ≤ 0.01")
    tot = pd.read_csv(cen / "rates.csv")
    tot = tot[(tot.score == "census") & (tot.site == "all") & (tot.winter == "all")].iloc[0]
    ax.set_title(f"(a) Days with streaks beyond chance: {tot.p01:.0f} observed, {tot.expected_p01:.1f} expected",
                 fontsize=10.5)
    ax.set_xlim(left=-0.3)
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right")

    # (b) supercooled cloud
    ax = axes[0, 1]
    s = pd.read_csv(cen / "slw.csv")
    s = s[s.stratum.astype(str) == "all"]
    k, n = s.p05.to_numpy(float), s.tested.to_numpy(float)
    z = 1.645
    ph = k / n
    c = (ph + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    x = np.arange(len(s))
    ax.bar(x, 100 * ph, width=0.6, color=figstyle.BLUE)
    ax.errorbar(x, 100 * ph, yerr=[100 * (ph - (c - h)), 100 * ((c + h) - ph)], fmt="none",
                ecolor=figstyle.INK_2, elinewidth=1.2, capsize=0)
    ax.axhline(5, color=figstyle.INK_2, lw=1, ls=(0, (4, 3)))
    ax.annotate("no seeding: 5%", (len(s) - 0.5, 5), xytext=(0, 4), textcoords="offset points",
                ha="right", fontsize=8.5, color=figstyle.INK_2)
    ax.set_xticks(x, [f"{lab}\n({int(m)} days)" for lab, m in zip(s.slw_class, s.tested)])
    ax.set_xlabel("hours of supercooled cloud water over the radar (HRRR)")
    ax.set_ylabel("radar-days with p ≤ 0.05 (%)")
    ax.set_title("(b) The excess comes on days with supercooled cloud", fontsize=10.5)
    ax.grid(axis="x", visible=False)

    # (c) efficiency by temperature, jets, arrivals and departures
    ax = axes[1, 0]
    e = pd.read_csv(cen / "efficiency.csv")
    e = e[(e.engine_class == "jet") & (e.tbin != "all") & e.phase.isin(["arrival", "departure"]) & (e.n >= 20)]
    for ph_, col, dx in (("arrival", figstyle.BLUE, -0.25), ("departure", figstyle.ORANGE, 0.25)):
        g = e[e.phase == ph_].copy()
        g["mid"] = g.tbin.map(tbin_mid)
        g = g.sort_values("mid")
        ax.errorbar(g.mid + dx, 100 * g.made_streak,
                    yerr=[100 * (g.made_streak - g.lo90), 100 * (g.hi90 - g.made_streak)],
                    fmt="o-", color=col, ms=6, lw=1.6, elinewidth=1.1, capsize=0, label=f"jet {ph_}s")
        ax.plot(g.mid + dx, 100 * g.chance, "o", mfc=figstyle.SURFACE, mec=col, ms=5, ls="none")
    ax.plot([], [], "o", mfc=figstyle.SURFACE, mec=figstyle.INK_2, ms=5, ls="none", label="chance (rotated tracks)")
    ax.set_xlabel("temperature where the aircraft crossed the layer (°C, HRRR)")
    ax.set_ylabel("first passages that left a streak (%)")
    ax.set_title("(c) Arrivals seed more often, and down to about −8 °C", fontsize=10.5)
    ax.invert_xaxis()
    ax.legend(loc="upper left")

    # (d) seeding temperatures against the quoted thresholds
    ax = axes[1, 1]
    S = pd.read_csv(cen / "seeders.csv")
    classes = [c for c in ("jet", "turboprop", "piston") if c in set(S.engine_class)]
    ax.axvspan(-25, -20, color=figstyle.BLUE, alpha=0.10, lw=0)
    ax.axvspan(-15, -10, color=figstyle.AQUA, alpha=0.12, lw=0)
    rng = np.random.default_rng(1)
    for i, cl in enumerate(classes):
        t = S.loc[S.engine_class == cl, "seed_temp_c"].to_numpy(float)
        yy = len(classes) - 1 - i + rng.uniform(-0.18, 0.18, len(t))
        ax.plot(t, yy, "o", ms=4.5, color=figstyle.SERIES[i], alpha=0.75, mec="none")
        ax.plot([np.median(t)] * 2, [len(classes) - 1 - i - 0.3, len(classes) - 1 - i + 0.3],
                color=figstyle.INK, lw=2)
    ax.set_yticks(range(len(classes))[::-1], [f"{c}s ({int((S.engine_class == c).sum())})" for c in classes])
    ax.text(-22.5, len(classes) - 0.45, "quoted for jets", ha="center", fontsize=8.5, color=figstyle.INK_2)
    ax.text(-12.5, len(classes) - 0.45, "quoted for\nturboprops", ha="center", fontsize=8.5, color=figstyle.INK_2)
    ax.set_ylim(-0.6, len(classes) - 0.1)
    ax.set_xlim(-27, -3)
    ax.invert_xaxis()
    ax.set_xlabel("seeding temperature of catalogue flights (°C, HRRR); bar = median")
    ax.set_title("(d) Catalogue jets seeded warmer than the usual jet threshold", fontsize=10.5)
    ax.grid(axis="y", visible=False)

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()
