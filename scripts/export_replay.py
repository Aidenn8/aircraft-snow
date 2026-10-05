"""Export a radar-day window for the interactive replay: reflectivity frames as PNG (base64),
flight tracks and matched streaks, in one JSON file.

usage: export_replay.py SITE YYYY-MM-DD START END OUT.json [--alt 6000,14000] [--alpha 0.01]
START/END are UTC times (HH:MM, END may be on the next day if earlier than START).
"""
import argparse
import base64
import datetime as dt
import io
import json
import pathlib

import numpy as np
import pandas as pd
from PIL import Image

from snow import flights, geo, match, pipeline, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _finite(v):
    """A JSON-safe number: None for missing or non-finite values (JSON has no NaN)."""
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def encode(z):
    """8-bit greyscale frame: 0 = no echo (< 0 dBZ or no data), 1-255 = 0-35 dBZ linearly.
    The page colours it with a theme-specific single-hue ramp."""
    v = np.where(np.isfinite(z) & (z >= 0), 1 + np.clip(z, 0, 35) / 35 * 254, 0)
    return v.astype(np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("start")
    ap.add_argument("end")
    ap.add_argument("out")
    ap.add_argument("--alt", default="6000,14000")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--airport", default="", help="LAT,LON,NAME")
    ap.add_argument("--tz-offset", type=float, default=0.0, help="local time minus UTC, hours")
    ap.add_argument("--tz", default="UTC")
    a = ap.parse_args()
    day = dt.date.fromisoformat(a.day)
    lo_alt, hi_alt = (float(v) for v in a.alt.split(","))
    prep = pipeline.prepare(a.site, day)
    G = prep["grids"]
    t0 = prep["t0"]
    def at(hhmm):
        h, m = (int(v) for v in hhmm.split(":"))
        t = t0.replace(hour=h, minute=m)
        return t if t >= t0 else t + dt.timedelta(days=1)
    ta, tb = at(a.start), at(a.end)
    if tb <= ta:
        tb += dt.timedelta(days=1)
    ks = np.nonzero((G.t >= ta.timestamp()) & (G.t <= tb.timestamp()))[0]
    half = 50_000.0
    n = G.z.shape[1]
    xs = G.x0 + G.res * np.arange(n)
    sel = np.nonzero(np.abs(xs) <= half)[0]
    frames = []
    for k in ks:
        z = G.z[k][np.ix_(sel, sel)]
        img = Image.fromarray(encode(z[::-1]), "L")               # north up, greyscale code
        buf = io.BytesIO()
        img.save(buf, format="PNG", optimize=True)
        frames.append({"t": float(G.t[k]), "png": base64.b64encode(buf.getvalue()).decode()})
    pos = prep["positions"]
    w = pos[(pos.alt_ft >= lo_alt) & (pos.alt_ft <= hi_alt) & (pos.t >= ta.timestamp() - 3600)
            & (pos.t <= tb.timestamp()) & (np.abs(pos.x) <= half) & (np.abs(pos.y) <= half)]
    calls = pd.concat([flights.load_calls(day, a.site), flights.load_calls(day + dt.timedelta(days=1), a.site)])
    ft = flights.flight_table(pos, calls)
    S = pd.read_parquet(ROOT / "data" / "events" / a.site / f"{day}_streaks.parquet")
    S = S[(S.coverage > 0) & (S.p_value <= a.alpha)]
    seeders = set(int(f) for f in S.flight)
    # drift (m/s) of seeding flights' track points, from the matched perturbation of their
    # best streak, so the page can carry their air downwind
    pts_all, drift_all = prep["points"], prep["drift"]
    best = S.sort_values("p_value").drop_duplicates("flight").set_index("flight")
    tracks = []
    for fid, g in w.groupby("flight"):
        if len(g) < 3:
            continue
        info = ft.loc[fid]
        tr = {"id": int(fid), "callsign": info.callsign or info.reg or "",
              "type": info.icao_type if isinstance(info.icao_type, str) else "",
              "engine": info.engine_class, "seeder": int(fid) in seeders,
              "t": [float(v) for v in g.t], "x": [round(float(v) / 1000, 2) for v in g.x],
              "y": [round(float(v) / 1000, 2) for v in g.y], "alt": [int(v) for v in g.alt_ft]}
        if tr["seeder"]:
            sel = pts_all.flight.to_numpy() == fid
            p = pts_all[sel]
            d = drift_all[sel]
            o = best.loc[fid]
            tr["seed"] = {"t": [float(v) for v in p.t], "x": [round(float(v) / 1000, 3) for v in p.x],
                          "y": [round(float(v) / 1000, 3) for v in p.y],
                          "u": [round(float(v) + float(o.du), 2) for v in d[:, 0]],
                          "v": [round(float(v) + float(o.dv), 2) for v in d[:, 1]],
                          "alt": [int(v) for v in p.alt_ft], "t_pass": float(o.t_pass),
                          "seed_alt_ft": _finite(o.get("seed_alt_ft", None))}
        tracks.append(tr)
    matches = [{"t": float(o.t), "flight": int(o.flight), "x": round(o.x / 1000, 2), "y": round(o.y / 1000, 2),
                "age_min": round(o.age_s / 60, 1), "p": float(o.p_value)}
               for _, o in S.iterrows() if ta.timestamp() <= o.t <= tb.timestamp()]
    lat0, lon0 = prep["sweep0"].radar_lat, prep["sweep0"].radar_lon
    # per-seeder details from the described matches (seeding altitude and HRRR temperature)
    desc_path = ROOT / "data" / "events" / a.site / f"{day}_described.csv"
    desc = pd.read_csv(desc_path) if desc_path.exists() else pd.DataFrame()
    seeder_info = {}
    for fid, o in best.iterrows():
        d = desc[desc.flight == fid] if len(desc) else desc
        seeder_info[int(fid)] = {
            "t_pass": float(o.t_pass), "t_match": float(o.t), "age_min": round(float(o.age_s) / 60, 1),
            "p": float(o.p_value),
            "seed_alt_ft": _finite(d.seed_alt_ft.iloc[0]) if len(d) else None,
            "seed_temp_c": _finite(d.seed_temp_c.iloc[0]) if len(d) else None}
    airport = None
    if a.airport:
        la, lo_, name = a.airport.split(",")
        ax_, ay_ = geo.to_xy(float(la), float(lo_), lat0, lon0)
        airport = {"x": round(float(ax_) / 1000, 2), "y": round(float(ay_) / 1000, 2), "name": name}
    out = {"site": a.site, "day": a.day, "radar": {"lat": lat0, "lon": lon0}, "half_km": half / 1000,
           "res_km": G.res / 1000, "frames": frames, "tracks": tracks, "matches": matches,
           "alt_band_ft": [lo_alt, hi_alt], "seeder_info": seeder_info, "airport": airport,
           "local_offset_h": a.tz_offset, "local_tz": a.tz}
    pathlib.Path(a.out).write_text(json.dumps(out, separators=(",", ":"), allow_nan=False))
    print(f"{len(frames)} frames, {len(tracks)} tracks ({len(seeders & set(t['id'] for t in tracks))} seeders), "
          f"{len(matches)} matches, {pathlib.Path(a.out).stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
