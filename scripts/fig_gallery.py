"""Gallery of census cases: one panel per aircraft streak, from different radars and days.

usage: fig_gallery.py CASES.csv --out FIG.png [--cols 3] [--alpha 0.05] [--title TEXT]

CASES.csv has columns site, day (YYYY-MM-DD, the 12 UTC radar-day), time (HH:MM UTC, the
scan to show) and callsign, and optionally half_km (panel half-width, default 22). Each
panel shows the lowest-tilt reflectivity around the streak the flight was matched to in
that scan (p <= --alpha), the flight's own track over the previous hour, and the part of
it that matched: the seeding points carried by the wind to the scan time. The panel title
gives the radar and time, the flight, its type, the seeding altitude and temperature (HRRR)
and the time since the aircraft passed.
"""
import argparse
import datetime as dt
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from snow import aircraft, figstyle, flights, match, pipeline, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]


def site_names() -> dict[str, str]:
    s = pd.read_csv(ROOT / "data" / "sites.csv").drop_duplicates("radar")
    return dict(zip(s.radar, s.municipality))


# names the type database gets wrong or ambiguous for these purposes
TYPE_NAMES = {"B77L": "Boeing 777-200LR / 777F"}
MAKERS = {"MCDONNELL-DOUGLAS": "McDonnell Douglas", "DE HAVILLAND CANADA": "De Havilland Canada"}


def type_name(icao: str) -> str:
    if icao in TYPE_NAMES:
        return TYPE_NAMES[icao]
    t = aircraft.types()
    if icao in t.index and isinstance(t.loc[icao, "type_name"], str):
        name = t.loc[icao, "type_name"]
        # the database is upper case; keep model codes as they are, title-case the maker
        maker, _, model = name.partition(" ")
        return f"{MAKERS.get(maker, maker.title())} {model}".strip()
    return icao


def minus(v: float) -> str:
    return f"{v:.0f}".replace("-", "\u2212")


def scan_time(day: dt.date, hhmm: str) -> dt.datetime:
    t0, _ = pipeline.day_window(day)
    h, m = (int(v) for v in hhmm.split(":"))
    T = t0.replace(hour=h, minute=m)
    return T + dt.timedelta(days=1) if T < t0 else T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cols", type=int, default=3)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--title", default="")
    a = ap.parse_args()
    figstyle.apply()
    cases = pd.read_csv(a.cases, dtype={"time": str})
    names = site_names()
    n = len(cases)
    cols = min(a.cols, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4.6 * cols + 1.0, 4.9 * rows + 0.6), squeeze=False,
                             layout="constrained")
    preps: dict = {}
    im = None
    for ax, c in zip(axes.flat, cases.itertuples()):
        day = dt.date.fromisoformat(c.day)
        key = (c.site, c.day)
        if key not in preps:
            prep = pipeline.prepare(c.site, day)
            calls = pd.concat([flights.load_calls(day, c.site), flights.load_calls(day + dt.timedelta(days=1), c.site)])
            S = pd.read_parquet(ROOT / "data" / "events" / c.site / f"{day}_streaks.parquet")
            preps[key] = (prep, flights.flight_table(prep["positions"], calls), S, pipeline.matching_points(prep))
        prep, ft, S, (pts, dr) = preps[key]
        G = prep["grids"]
        xs = G.x0 + G.res * np.arange(G.z.shape[2])
        ys = G.y0 + G.res * np.arange(G.z.shape[1])
        k = int(np.argmin(np.abs(G.t - scan_time(day, c.time).timestamp())))
        Tk = G.t[k]
        ids = ft.index[ft.callsign.str.strip() == c.callsign.strip()]
        m = S[(np.abs(S.t - Tk) < 1.0) & S.flight.isin(ids) & (S.coverage > 0) & (S.p_value <= a.alpha)]
        if m.empty:
            raise SystemExit(f"{c.site} {c.day} {c.time} {c.callsign}: no match at p <= {a.alpha} in that scan")
        o = m.sort_values("p_value").iloc[0]
        info = ft.loc[o.flight]
        z = G.z[k]
        zz = np.where(np.isnan(z), streaks.FLOOR_DBZ, z)
        _, lab = streaks.find_streaks(streaks.GridScan(Tk, xs, ys, zz, np.isfinite(z)))
        sx, sy = match.skeleton_xy(lab, int(o.label), xs, ys)
        mp = match.matched_points(sx, sy, pts, dr, Tk, int(o.flight), o.du, o.dv)
        age = Tk - mp.t.to_numpy(float)
        idx = mp.index.to_numpy()
        px = mp.x.to_numpy() + (dr[idx, 0] + o.du) * age
        py = mp.y.to_numpy() + (dr[idx, 1] + o.dv) * age
        order = np.argsort(mp.t.to_numpy())
        half = float(getattr(c, "half_km", np.nan)) if hasattr(c, "half_km") else np.nan
        half = 22.0 if not np.isfinite(half) else half
        cx, cy = o.x / 1e3, o.y / 1e3
        im = ax.pcolormesh(xs / 1e3, ys / 1e3, np.where(np.isnan(z), -99, z), cmap=figstyle.REFLECTIVITY,
                           vmin=0, vmax=30, shading="nearest", rasterized=True)
        pos = prep["positions"]
        tr = pos[(pos.flight == o.flight) & (pos.t > Tk - 3600) & (pos.t <= Tk)]
        ax.plot(tr.x / 1e3, tr.y / 1e3, color=figstyle.INK_2, lw=0.9, ls=(0, (3, 2)))
        ax.plot(px[order] / 1e3, py[order] / 1e3, color=figstyle.ORANGE, lw=2.0)
        if abs(cx) < half and abs(cy) < half:
            ax.plot(0, 0, marker="^", color=figstyle.INK, ms=6, ls="none")
        figstyle.map_axes(ax)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xlim(cx - half, cx + half)
        ax.set_ylim(cy - half, cy + half)
        when = pd.Timestamp(Tk, unit="s", tz="UTC")
        alt = float(np.median(mp.alt_ft)) if "alt_ft" in mp else float(o.seed_alt_ft)
        temp = float(np.median(mp.temp_c)) if "temp_c" in mp else np.nan
        ax.set_title(f"{names.get(c.site, c.site)} ({c.site}), {when:%-d %b %Y %H:%M} UTC\n"
                     f"{info.callsign.strip()} · {type_name(info.icao_type)}\n"
                     f"{round(alt, -2):,.0f} ft · {minus(temp)} °C · {o.age_s / 60:.0f} min after passing",
                     fontsize=9.5, fontweight="normal", linespacing=1.35)
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.supxlabel("km east of the radar", fontsize=9.5, color=figstyle.INK_2)
    fig.supylabel("km north of the radar", fontsize=9.5, color=figstyle.INK_2)
    cb = fig.colorbar(im, ax=axes, shrink=0.45, pad=0.01, extend="max")
    cb.set_label("reflectivity at the lowest radar tilt (dBZ)")
    cb.outline.set_visible(False)
    handles = [plt.Line2D([], [], color=figstyle.ORANGE, lw=2.0,
                          label="where the aircraft flew, carried by the wind to the scan time"),
               plt.Line2D([], [], color=figstyle.INK_2, lw=0.9, ls=(0, (3, 2)),
                          label="the aircraft's track, previous hour")]
    fig.legend(handles=handles, loc="outside upper center", ncol=2)
    fig.get_layout_engine().set(h_pad=0.12, hspace=0.04)
    if a.title:
        fig.suptitle(a.title, x=0.02, ha="left", fontsize=13, fontweight="semibold", color=figstyle.INK)
    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()
