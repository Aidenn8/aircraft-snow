"""Split continental ADS-B day files into small per-radar files.

Each census radar gets data/adsb_radar/SITE/YYYY-MM-DD.parquet with every position (any
altitude, airborne or not) within RADIUS_KM of the radar, plus the matching callsign
records. Reading these instead of the continental file makes a radar-day ~20x cheaper.

usage: python -m snow.adsb_split DATES... [--sites census_sites.csv]
"""
from __future__ import annotations

import argparse
import datetime as dt
import pathlib

import numpy as np
import pandas as pd

from . import adsb, geo, sites

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "adsb_radar"
RADIUS_KM = 100.0


def split_day(day: dt.date, site_list, force=False):
    todo = [s for s in site_list if force or not (OUT / s / f"{day}.parquet").exists()]
    if not todo or not adsb.day_path(day).exists():
        return 0
    pos = adsb.load_day(day)
    calls = adsb.load_calls(day)
    lat = pos.lat.to_numpy()
    lon = pos.lon.to_numpy()
    for s in todo:
        lat0, lon0 = sites.location(s)
        dlat = RADIUS_KM / 111.0
        dlon = RADIUS_KM / (111.0 * np.cos(np.radians(lat0)))
        box = (lat > lat0 - dlat) & (lat < lat0 + dlat) & (lon > lon0 - dlon) & (lon < lon0 + dlon)
        sub = pos[box]
        d = geo.haversine(sub.lat.to_numpy(), sub.lon.to_numpy(), lat0, lon0)
        sub = sub[d <= RADIUS_KM * 1000]
        c = calls[calls.hex.isin(np.unique(sub.hex.to_numpy()))]
        (OUT / s).mkdir(parents=True, exist_ok=True)
        sub.to_parquet(OUT / s / f"{day}.parquet", index=False, compression="zstd")
        c.to_parquet(OUT / s / f"{day}.calls.parquet", index=False, compression="zstd")
    return len(todo)


def census_sites(path=ROOT / "config" / "census_sites.csv"):
    return list(pd.read_csv(path).radar)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dates", nargs="+")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    s_list = census_sites() + ["KFTG"]
    s_list = list(dict.fromkeys(s_list))
    for day in adsb.parse_dates(a.dates):
        n = split_day(day, s_list, a.force)
        if n:
            print(f"{day}: split for {n} radars", flush=True)


if __name__ == "__main__":
    main()
