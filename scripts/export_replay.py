"""Export a radar-day window for the interactive replay: reflectivity frames as PNG (base64),
flight tracks and matched streaks, in one JSON file.

usage: export_replay.py SITE YYYY-MM-DD START END OUT.json [--alt 6000,14000] [--alpha 0.005]
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

from snow import flights, match, pipeline, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]
# reflectivity colour table: transparent below 0 dBZ, then a perceptual ramp to 35 dBZ
STOPS = [(0, (210, 225, 240)), (5, (140, 190, 230)), (10, (70, 150, 210)), (15, (40, 170, 120)),
         (20, (240, 200, 60)), (25, (240, 130, 40)), (30, (215, 50, 40)), (35, (150, 20, 90))]


def colour(z):
    """RGBA image from floored reflectivity; below 0 dBZ fully transparent."""
    lv = np.array([s[0] for s in STOPS], float)
    rgb = np.array([s[1] for s in STOPS], float)
    zc = np.clip(z, lv[0], lv[-1])
    out = np.zeros(z.shape + (4,), np.uint8)
    for c in range(3):
        out[..., c] = np.interp(zc, lv, rgb[:, c]).astype(np.uint8)
    alpha = np.clip((z - 0.0) / 3.0, 0, 1) * 230
    out[..., 3] = np.where(np.isfinite(z), alpha, 0).astype(np.uint8)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("start")
    ap.add_argument("end")
    ap.add_argument("out")
    ap.add_argument("--alt", default="6000,14000")
    ap.add_argument("--alpha", type=float, default=0.005)
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
        img = Image.fromarray(colour(z[::-1]), "RGBA")          # north up
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
    tracks = []
    for fid, g in w.groupby("flight"):
        if len(g) < 3:
            continue
        info = ft.loc[fid]
        tracks.append({"id": int(fid), "callsign": info.callsign or info.reg or "",
                       "type": info.icao_type if isinstance(info.icao_type, str) else "",
                       "engine": info.engine_class, "seeder": int(fid) in seeders,
                       "t": [float(v) for v in g.t], "x": [round(float(v) / 1000, 2) for v in g.x],
                       "y": [round(float(v) / 1000, 2) for v in g.y], "alt": [int(v) for v in g.alt_ft]})
    matches = [{"t": float(o.t), "flight": int(o.flight), "x": round(o.x / 1000, 2), "y": round(o.y / 1000, 2),
                "age_min": round(o.age_s / 60, 1), "p": float(o.p_value)}
               for _, o in S.iterrows() if ta.timestamp() <= o.t <= tb.timestamp()]
    lat0, lon0 = prep["sweep0"].radar_lat, prep["sweep0"].radar_lon
    out = {"site": a.site, "day": a.day, "radar": {"lat": lat0, "lon": lon0}, "half_km": half / 1000,
           "res_km": G.res / 1000, "frames": frames, "tracks": tracks, "matches": matches,
           "alt_band_ft": [lo_alt, hi_alt]}
    pathlib.Path(a.out).write_text(json.dumps(out, separators=(",", ":")))
    print(f"{len(frames)} frames, {len(tracks)} tracks ({len(seeders & set(t['id'] for t in tracks))} seeders), "
          f"{len(matches)} matches, {pathlib.Path(a.out).stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
