"""How often an aircraft passing through the seeding layer left a streak, by temperature.

usage: fig_efficiency.py --exposure A_exposure.csv [B_exposure.csv ...] --out FIG.png
       [--bins -26,-20,-17,-14,-11,-8,-4] [--min-n 8]

Input: the per-flight tables written by seeding_efficiency.py for event days. Only first
passages count (the flight's air was echo-free before it arrived), since a later aircraft
on the same path flies through ice that is already there. For each temperature bin and
engine class: the share of flights whose drifted track carried a streak (F1 >= 0.5), and
the chance level, the same share for the flights rotated about the radar. Bars are 90%
Wilson intervals; bins with fewer than --min-n flights are not drawn. Also writes the table
next to the figure (.csv).
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from snow import figstyle


def wilson(k, n, z=1.645):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / d
    return c - h, c + h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exposure", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--bins", default="-26,-20,-17,-14,-11,-8,-4")
    ap.add_argument("--min-n", type=int, default=8)
    a = ap.parse_args()
    figstyle.apply()
    df = pd.concat([pd.read_csv(p).assign(source=pathlib.Path(p).name) for p in a.exposure], ignore_index=True)
    df = df[df.first_passage.astype(bool)]
    edges = [float(v) for v in a.bins.split(",")]
    df["tbin"] = pd.cut(df.temp_c, edges)
    rows = []
    for (cls, tb), g in df.groupby(["engine_class", "tbin"], observed=True):
        n, k = len(g), int(g.made_streak.sum())
        lo, hi = wilson(k, n)
        rows.append({"engine_class": cls, "t_lo": tb.left, "t_hi": tb.right, "t_mid": tb.mid, "n": n,
                     "made_streak": k, "share": k / n, "ci_lo": lo, "ci_hi": hi,
                     "chance": float(g.placebo_rate.mean())})
    tab = pd.DataFrame(rows)
    tab.to_csv(pathlib.Path(a.out).with_suffix(".csv"), index=False)
    print(tab.round(3).to_string(index=False))

    fig, ax = plt.subplots(figsize=(6.8, 4.0), layout="constrained")
    styles = {"jet": (figstyle.BLUE, "o", -0.25), "turboprop": (figstyle.ORANGE, "s", 0.25)}
    for cls, (col, mk, dx) in styles.items():
        t = tab[(tab.engine_class == cls) & (tab.n >= a.min_n)]
        if t.empty:
            continue
        ax.errorbar(t.t_mid + dx, t.share, yerr=[t.share - t.ci_lo, t.ci_hi - t.share], fmt=mk, color=col,
                    ms=7, capsize=0, elinewidth=1.4, label=f"{cls}s, left a streak")
        ax.plot(t.t_mid + dx, t.chance, mk, mfc="none", mec=col, ms=7, mew=1.3, alpha=0.8,
                label=f"{cls}s, chance (rotated)")
        for _, r in t.iterrows():
            ax.annotate(f"n={r.n}", (r.t_mid + dx, r.ci_hi), textcoords="offset points", xytext=(0, 4),
                        ha="center", fontsize=7.5, color=figstyle.MUTED)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, 1.0)
    ax.invert_xaxis()
    ax.set_xlabel("HRRR temperature where the aircraft crossed the layer (°C)")
    ax.set_ylabel("share of first passages")
    ax.set_title("Aircraft that left a streak, by temperature")
    ax.legend(loc="upper left")
    fig.savefig(a.out)
    print(a.out)


if __name__ == "__main__":
    main()
