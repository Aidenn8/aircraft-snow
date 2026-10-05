"""Run the event test for one radar and one analysis day.

An analysis day starts at 12:00 UTC (early morning in the US) and lasts 24 hours, so a US
evening is not split across two UTC dates. Passages in the day are tested; radar scans from
30 min before the day to 60 min after it are used.
"""
from __future__ import annotations

import datetime as dt
import json
import pathlib
import time

import numpy as np
import pandas as pd

from . import adsb, epoch, event, flights, hrrr, radar
from .radar import beam_height_ground_range

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "events"
UTC = dt.timezone.utc
DAY_START_HOUR = 12
MIN_RANGE_M = 8_000.0
MAX_RANGE_M = 65_000.0
ALT_RANGE_FT = (2_000.0, 24_000.0)
BEAM_MARGIN_M = 500.0
HRRR_BOX = 2


def day_window(day: dt.date):
    t0 = dt.datetime(day.year, day.month, day.day, DAY_START_HOUR, tzinfo=UTC)
    return t0, t0 + dt.timedelta(hours=24)


def load_scans(site: str, t_a: dt.datetime, t_b: dt.datetime):
    site3 = site[1:] if len(site) == 4 else site
    files = []
    d = t_a.date()
    while d <= t_b.date():
        files += [f for f in radar.list_level3(site3, "N0B", d) if t_a <= f[0] <= t_b]
        d += dt.timedelta(days=1)
    scans = []
    for t, key, _ in files:
        try:
            scans.append((t.timestamp(), radar.read_level3(radar.fetch(radar.L3_BUCKET, key))))
        except Exception as e:      # an occasional corrupt file should not stop the day
            print(f"  skip {key}: {e!r}")
    return scans


def prepare(site: str, day: dt.date, scans=None):
    """Everything the event test needs: grids, points with drift, the radar sweep header."""
    t0, t1 = day_window(day)
    if scans is None:
        scans = load_scans(site, t0 - dt.timedelta(minutes=30), t1 + dt.timedelta(minutes=60))
    if len(scans) < 20:
        raise RuntimeError(f"{site} {day}: only {len(scans)} radar scans")
    sw0 = scans[0][1]
    lat0, lon0 = sw0.radar_lat, sw0.radar_lon
    pos = flights.segment(flights.load_window(lat0, lon0, MAX_RANGE_M / 1000 + 5,
                                              t0 - dt.timedelta(minutes=5), t1))
    pts = pos[(pos.alt_ft >= ALT_RANGE_FT[0]) & (pos.alt_ft < ALT_RANGE_FT[1])].copy()
    rng = np.hypot(pts.x, pts.y).to_numpy()
    h_beam, _ = beam_height_ground_range(rng, sw0.elevation, sw0.radar_height)
    # pressure altitude ~ altitude above sea level to within a few hundred feet
    alt_m = pts.alt_ft.to_numpy() * 0.3048
    keep = (h_beam <= alt_m + BEAM_MARGIN_M) & (rng >= MIN_RANGE_M) & (rng <= MAX_RANGE_M)
    pts = pts[keep].reset_index(drop=True)
    runs = [t0.replace(tzinfo=None) + dt.timedelta(hours=k) for k in range(-1, 26)]
    profs = pd.concat([hrrr.profile(r, lat0, lon0, box=HRRR_BOX) for r in runs], ignore_index=True)
    drift = epoch.drift_velocity(pts, profs, sw0.radar_height, sw0.elevation)
    # every aircraft at any altitude (including ones on the ground) is masked out of the radar
    allpos = flights.load_window(lat0, lon0, MAX_RANGE_M / 1000 + 15, t0 - dt.timedelta(hours=1),
                                 t1 + dt.timedelta(hours=1.5), max_alt_ft=60_000)
    grids = event.Grids.from_sweeps(scans, aircraft=allpos[["hex", "t", "x", "y", "gs"]])
    return {"grids": grids, "points": pts, "drift": drift, "profiles": profs, "sweep0": sw0,
            "positions": pos, "t0": t0}


def run_day(site: str, day: dt.date, save: bool = True, force: bool = False) -> dict:
    out_dir = OUT / site
    summ_path = out_dir / f"{day:%Y-%m-%d}.json"
    if summ_path.exists() and not force:
        return json.loads(summ_path.read_text())
    tic = time.time()
    prep = prepare(site, day)
    res = event.run(prep["grids"], prep["points"], prep["drift"], prep["t0"].timestamp(), 24)
    summary = {"site": site, "day": f"{day}", "n_scans": int(len(prep["grids"].t)),
               "n_points": int(len(prep["points"])),
               "n_flights": int(prep["points"].flight.nunique()), "stats": {}}
    zmaps = {}
    for stat in event.STATS + ("new",):
        ev = event.evaluate_clusters(res, stat)
        zmaps[stat] = ev.pop("z_map")
        ev.pop("placebo_masses")
        summary["stats"][stat] = ev
    summary["seconds"] = round(time.time() - tic, 1)
    if save:
        out_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out_dir / f"{day:%Y-%m-%d}_acc.npz", acc_real=res["acc_real"],
                            acc_placebo=res["acc_placebo"], new_real=res["new_real"],
                            new_placebo=res["new_placebo"], band_edges=res["band_edges"],
                            angles=res["angles"], t_start=res["t_start"],
                            **{f"zmap_{k}": v for k, v in zmaps.items()})
        summ_path.write_text(json.dumps(summary, indent=1, default=float))
    return summary
