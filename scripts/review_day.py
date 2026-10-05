"""Review sheet for one radar-day: the event test's z-map, the streak-match timeline, and
zoomed radar panels of the best-matched streaks with the matched flight's track.

usage: review_day.py SITE YYYY-MM-DD [--top N]
Writes review/SITE_YYYY-MM-DD.png
"""
import argparse
import datetime as dt
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from snow import adsb, flights, pipeline, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("--top", type=int, default=8)
    a = ap.parse_args()
    day = dt.date.fromisoformat(a.day)
    ev_dir = ROOT / "data" / "events" / a.site
    summ = json.loads((ev_dir / f"{day}.json").read_text())
    acc = np.load(ev_dir / f"{day}_acc.npz")
    st_path = ev_dir / f"{day}_streaks.parquet"
    S = pd.read_parquet(st_path) if st_path.exists() else pd.DataFrame()
    prep = pipeline.prepare(a.site, day)
    G, pts = prep["grids"], prep["points"]
    calls = pd.concat([adsb.load_calls(day), adsb.load_calls(day + dt.timedelta(days=1))])
    ft = flights.flight_table(prep["positions"], calls)
    xs = G.x0 + G.res * np.arange(G.z.shape[2])
    ys = G.y0 + G.res * np.arange(G.z.shape[1])

    top = S[S.coverage > 0].sort_values(["p_value", "coverage"], ascending=[True, False]).head(a.top) if len(S) else S
    n_top = len(top)
    cols = 4
    rows = 1 + int(np.ceil(max(n_top, 1) / cols))
    fig = plt.figure(figsize=(4.5 * cols, 4.4 * rows))
    # z-map of the primary statistic and of the new-echo statistic
    edges = acc["band_edges"]
    t0 = float(acc["t_start"])
    for i, stat in enumerate(("z", "new")):
        ax = fig.add_subplot(rows, cols, 1 + i)
        zm = acc[f"zmap_{stat}"]
        im = ax.imshow(zm, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-4, vmax=4,
                       extent=[0, zm.shape[1], edges[0] / 1000, edges[0] / 1000 + zm.shape[0]])
        ax.set_xlabel("passage hour (UTC)")
        ticks = np.arange(0, 24, 3)
        ax.set_xticks(ticks + 0.5)
        ax.set_xticklabels([f"{(dt.datetime.fromtimestamp(t0, dt.timezone.utc).hour + h) % 24:02d}" for h in ticks])
        ax.set_ylabel("altitude band start (kft)")
        p = summ["stats"][stat]["p_cluster"]
        ax.set_title(f"event test '{stat}': cluster p = {p:.3f}", fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.8, label="z")
    ax = fig.add_subplot(rows, cols, 3)
    if len(S):
        tt = pd.to_datetime(S.t, unit="s", utc=True)
        sc = ax.scatter(tt, S.coverage, c=-np.log10(S.p_value.clip(1e-3)), cmap="viridis", vmin=0, vmax=2, s=12)
        fig.colorbar(sc, ax=ax, shrink=0.8, label="-log10 p")
        ax.set_ylabel("best new single-flight coverage")
        ax.tick_params(axis="x", labelrotation=45, labelsize=7)
    sk = summ.get("streaks", {})
    ax.set_title(f"streak matches: {sk.get('n_matched', 0)}/{sk.get('n_streaks', 0)}  "
                 f"E={sk.get('E', 0):.1f} p_day={sk.get('p_day', 1):.3f}", fontsize=9)
    ax = fig.add_subplot(rows, cols, 4)
    ax.axis("off")
    ax.text(0, 1, f"{a.site} {day} (12 UTC to 12 UTC)\nflights {summ['n_flights']}, scans {summ['n_scans']}",
            va="top", fontsize=10)
    for j, (_, o) in enumerate(top.iterrows()):
        ax = fig.add_subplot(rows, cols, cols + 1 + j)
        k = int(np.argmin(np.abs(G.t - o.t)))
        ax.pcolormesh(xs / 1e3, ys / 1e3, np.where(np.isnan(G.z[k]), -10, G.z[k]), cmap="Greys",
                      vmin=-5, vmax=25, shading="nearest")
        z = np.where(np.isnan(G.z[k]), streaks.FLOOR_DBZ, G.z[k])
        obj, lab = streaks.find_streaks(streaks.GridScan(o.t, xs, ys, z, np.isfinite(G.z[k])))
        ax.contour(xs / 1e3, ys / 1e3, lab == o.label, levels=[0.5], colors="C2", linewidths=1)
        f = pts[pts.flight == o.flight]
        age = o.t - f.t.to_numpy()
        ok = (age >= 0) & (age <= 60 * 60)
        dx = (prep["drift"][f.index.to_numpy(), 0] + o.du) * age
        dy = (prep["drift"][f.index.to_numpy(), 1] + o.dv) * age
        ax.plot(f.x[ok] / 1e3, f.y[ok] / 1e3, ".", ms=1.5, color="C0", label="track")
        ax.plot((f.x[ok] + dx[ok]) / 1e3, (f.y[ok] + dy[ok]) / 1e3, ".", ms=1.5, color="C3", label="drifted")
        ax.set_xlim(o.x / 1e3 - 25, o.x / 1e3 + 25)
        ax.set_ylim(o.y / 1e3 - 25, o.y / 1e3 + 25)
        ax.set_aspect(1)
        info = ft.loc[o.flight] if o.flight in ft.index else None
        name = f"{info.callsign} {info.icao_type}" if info is not None else "?"
        ax.set_title(f"{pd.Timestamp(o.t, unit='s', tz='UTC'):%H:%M}Z {name}\ncov {o.coverage:.2f} "
                     f"p {o.p_value:.3f} age {o.age_s / 60:.0f} min  alt {f.alt_ft.median():.0f} ft", fontsize=8)
        if j == 0:
            ax.legend(fontsize=7, markerscale=5, loc="lower left")
    out = ROOT / "review" / f"{a.site}_{day}.png"
    out.parent.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=60)
    print(out)


if __name__ == "__main__":
    main()
