"""Exploratory panels: lowest-tilt radar reflectivity with recent low-level flight tracks.

usage: explore_case.py SITE START END OUT.png [--l2] [--alt0 FT --alt1 FT] [--radius KM]
Times are UTC ISO strings. Tracks shown are those within the pressure-altitude band
during the 30 minutes before each radar scan.
"""
import argparse
import datetime as dt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from snow import adsb, geo, radar


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("start")
    ap.add_argument("end")
    ap.add_argument("out")
    ap.add_argument("--l2", action="store_true")
    ap.add_argument("--alt0", type=float, default=0)
    ap.add_argument("--alt1", type=float, default=30000)
    ap.add_argument("--radius", type=float, default=60)
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--lookback", type=float, default=30)
    a = ap.parse_args()
    t0 = dt.datetime.fromisoformat(a.start).replace(tzinfo=dt.timezone.utc)
    t1 = dt.datetime.fromisoformat(a.end).replace(tzinfo=dt.timezone.utc)
    days = sorted({t0.date(), t1.date(), (t0 - dt.timedelta(hours=1)).date()})
    if a.l2:
        files = [f for d in days for f in radar.list_level2(a.site, d)]
    else:
        files = [f for d in days for f in radar.list_level3(a.site[1:], "N0B", d)]
    files = [f for f in files if t0 <= f[0] <= t1][:: a.every]
    pos = pd.concat([adsb.load_day(d) for d in days if adsb.day_path(d).exists()])
    n = len(files)
    cols = min(4, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 4.2 * rows), squeeze=False)
    R = a.radius * 1000
    for ax, (t, key, _) in zip(axes.flat, files):
        if a.l2:
            sw = radar.read_level2(radar.fetch(radar.L2_BUCKET, key), moments=("REF",))[0]
        else:
            sw = radar.read_level3(radar.fetch(radar.L3_BUCKET, key))
        x, y, h = geo.gate_geometry(sw)
        z = sw.data["REF"]
        m = (np.abs(x) < R) & (np.abs(y) < R)
        ax.pcolormesh(x / 1e3, y / 1e3, np.where(m, z, np.nan), cmap="turbo", vmin=-10, vmax=30,
                      shading="auto")
        ts = t.timestamp()
        w = pos[(pos.t > ts - a.lookback * 60) & (pos.t <= ts) & (pos.alt * 25.0 >= a.alt0)
                & (pos.alt * 25.0 <= a.alt1)]
        px, py = geo.to_xy(w.lat.to_numpy(), w.lon.to_numpy(), sw.radar_lat, sw.radar_lon)
        keep = (np.abs(px) < R) & (np.abs(py) < R)
        ax.scatter(px[keep] / 1e3, py[keep] / 1e3, s=0.4, c="k", alpha=0.6, lw=0)
        ax.set_xlim(-a.radius, a.radius)
        ax.set_ylim(-a.radius, a.radius)
        ax.set_aspect(1)
        ax.set_title(f"{t:%H:%M:%S}Z el {sw.elevation:.1f}", fontsize=9)
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.suptitle(f"{a.site} {t0:%Y-%m-%d} reflectivity, tracks {a.alt0:.0f}-{a.alt1:.0f} ft (prev {a.lookback:.0f} min)")
    fig.tight_layout()
    fig.savefig(a.out, dpi=80)


if __name__ == "__main__":
    main()
