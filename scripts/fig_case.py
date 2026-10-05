"""Case-study figure: radar snapshots with the drifted tracks of the aircraft matched to
streaks (p <= alpha), and recent traffic for context.

usage: fig_case.py SITE YYYY-MM-DD HH:MM [HH:MM ...] --out FIG.png
       [--box x0,x1,y0,y1] [--alt 8000,13000] [--alpha 0.01] [--airport LAT,LON,NAME]
"""
import argparse
import datetime as dt
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from snow import figstyle, geo, match, pipeline, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("times", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--box", default="-45,30,-48,27")
    ap.add_argument("--alt", default="8000,13000")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--airport", default="")
    ap.add_argument("--title", default="")
    ap.add_argument("--cols", type=int, default=2)
    a = ap.parse_args()
    figstyle.apply()
    day = dt.date.fromisoformat(a.day)
    x0, x1, y0, y1 = (float(v) for v in a.box.split(","))
    lo, hi = (float(v) for v in a.alt.split(","))
    prep = pipeline.prepare(a.site, day)
    G = prep["grids"]
    pts, dr = pipeline.matching_points(prep)
    pos = prep["positions"]
    S = pd.read_parquet(ROOT / "data" / "events" / a.site / f"{day}_streaks.parquet")
    S = S[(S.coverage > 0) & (S.p_value <= a.alpha)]
    xs = G.x0 + G.res * np.arange(G.z.shape[2])
    ys = G.y0 + G.res * np.arange(G.z.shape[1])
    t0 = prep["t0"]
    ap_xy = None
    if a.airport:
        la, lo_, name = a.airport.split(",")
        ax_, ay_ = geo.to_xy(float(la), float(lo_), prep["sweep0"].radar_lat, prep["sweep0"].radar_lon)
        ap_xy = (float(ax_) / 1e3, float(ay_) / 1e3, name)
    desc_path = ROOT / "data" / "events" / a.site / f"{day}_described.csv"
    desc = pd.read_csv(desc_path) if desc_path.exists() else pd.DataFrame(columns=["flight"])
    names = {int(r.flight): f"{r.callsign} ({r.icao_type})" for _, r in desc.iterrows()}
    n = len(a.times)
    cols = min(a.cols, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(5.0 * cols + 1.2, 4.7 * rows + 0.6), squeeze=False,
                             layout="constrained")
    for ax, hhmm in zip(axes.flat, a.times):
        h, m = (int(v) for v in hhmm.split(":"))
        T = t0.replace(hour=h, minute=m)
        if T < t0:
            T += dt.timedelta(days=1)
        k = int(np.argmin(np.abs(G.t - T.timestamp())))
        Tk = G.t[k]
        z = G.z[k]
        im = ax.pcolormesh(xs / 1e3, ys / 1e3, np.where(np.isnan(z), -99, z), cmap=figstyle.REFLECTIVITY,
                           vmin=0, vmax=30, shading="nearest", rasterized=True)
        # context: actual tracks of traffic in the altitude band over the previous 45 min
        w = pos[(pos.t > Tk - 2700) & (pos.t <= Tk) & (pos.alt_ft >= lo) & (pos.alt_ft <= hi)]
        for _, g in w.groupby("flight"):
            ax.plot(g.x / 1e3, g.y / 1e3, color=figstyle.MUTED, lw=0.5, alpha=0.55)
        # matched flights at this scan: their seeding points drifted to the scan time
        here = S[(np.abs(S.t - Tk) < 1.0)          # same scan (epoch seconds), inside the panel
                 & (S.x / 1e3 > x0) & (S.x / 1e3 < x1) & (S.y / 1e3 > y0) & (S.y / 1e3 < y1)]
        zz = np.where(np.isnan(z), streaks.FLOOR_DBZ, z)
        obj, lab = streaks.find_streaks(streaks.GridScan(Tk, xs, ys, zz, np.isfinite(z)))
        for _, o in here.iterrows():
            sx, sy = match.skeleton_xy(lab, int(o.label), xs, ys)
            mp = match.matched_points(sx, sy, pts, dr, Tk, int(o.flight), o.du, o.dv)
            if mp.empty:
                continue
            age = Tk - mp.t.to_numpy(float)
            idx = mp.index.to_numpy()
            px = mp.x.to_numpy() + (dr[idx, 0] + o.du) * age
            py = mp.y.to_numpy() + (dr[idx, 1] + o.dv) * age
            order = np.argsort(mp.t.to_numpy())
            ax.plot(px[order] / 1e3, py[order] / 1e3, color=figstyle.ORANGE, lw=1.6, alpha=0.95)
        ax.plot(0, 0, marker="^", color=figstyle.INK, ms=6, ls="none")
        if ap_xy:
            ax.plot(ap_xy[0], ap_xy[1], marker="s", color=figstyle.INK, ms=5, ls="none")
            ax.annotate(ap_xy[2], (ap_xy[0], ap_xy[1]), xytext=(5, 4), textcoords="offset points",
                        fontsize=8.5, color=figstyle.INK)
        ax.annotate("radar", (0, 0), xytext=(5, -10), textcoords="offset points", fontsize=8.5,
                    color=figstyle.INK_2)
        figstyle.map_axes(ax)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        local = pd.Timestamp(Tk, unit="s", tz="UTC")
        who = "; ".join(f"{names.get(int(o.flight), 'flight')}, {o.age_s / 60:.0f} min after passing"
                        for _, o in here.iterrows())
        ax.set_title(f"{local:%H:%M} UTC" + (f"  {who}" if who else ""), fontsize=10)
    for ax in axes.flat[n:]:
        ax.axis("off")
    cb = fig.colorbar(im, ax=axes, shrink=0.5, pad=0.01, extend="max")
    cb.set_label("reflectivity at the lowest radar tilt (dBZ)")
    cb.outline.set_visible(False)
    handles = [plt.Line2D([], [], color=figstyle.ORANGE, lw=1.6, label="matched aircraft: seeded air carried by the wind to this time"),
               plt.Line2D([], [], color=figstyle.MUTED, lw=0.8, label=f"other traffic {lo / 1000:.1f}-{hi / 1000:.0f} kft, previous 45 min")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2)
    if a.title:
        fig.suptitle(a.title, x=0.02, ha="left", fontsize=13, fontweight="semibold", color=figstyle.INK)
    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()
