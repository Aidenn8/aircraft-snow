"""Supercooled liquid water above each radar on each census day, from the HRRR analyses.

usage: slw_covariate.py [--census data/census] [--threshold 0.01]

Aircraft make ice only where they fly through cloud of supercooled liquid droplets, so days
on which seeding is genuinely detected should be days with HRRR cloud water at temperatures
of -4 C or colder at the levels flown (matching uses 2 000-24 000 ft). Artefacts of the
method (snow bands, clutter) have no reason to follow that. For each radar-day in
CENSUS/radar_days.csv this reads the hourly profiles the census already fetched (the 5 x 5
point average over the radar; analyses for 12 UTC to 11 UTC) and records

  slw_hours      hours (of 24) in which some level in 2 000-24 000 ft with T <= -4 C has cloud
                 water mixing ratio >= --threshold g/kg
  slw_hours_low  the same with a ceiling of 16 000 ft (where streaks are seen best)
  slw_max_gkg    largest such cloud water over the day
  cold_hours     hours in which some level in the layer is at -4 C or colder
  t700_c         mean 700 hPa temperature

and writes CENSUS/slw_days.csv. Missing profiles are fetched (and cached) as in the census.
"""
import argparse
import datetime as dt
import pathlib

import numpy as np
import pandas as pd

from snow import hrrr, pipeline

ROOT = pathlib.Path(__file__).resolve().parents[1]
FT = 0.3048
LAYER_FT = (2_000.0, 24_000.0)
LOW_CEILING_FT = 16_000.0
T_MAX_C = -4.0


def cached_locations() -> dict[str, tuple[float, float]]:
    """Radar site -> the (lat, lon) spelling most used in the profile cache (the Level III
    header's coordinates, written with differing precision in different years; all fall in
    the same HRRR box)."""
    sites = pd.read_csv(ROOT / "data" / "sites.csv").drop_duplicates("radar").set_index("radar")
    counts: dict[tuple[float, float], int] = {}
    for f in (hrrr.CACHE / "profiles").glob("*_2.parquet"):
        parts = f.stem.split("_")
        key = (float(parts[2]), float(parts[3]))
        counts[key] = counts.get(key, 0) + 1
    out = {}
    for radar, r in sites.iterrows():
        near = [(n, k) for k, n in counts.items()
                if abs(k[0] - r.radar_lat) < 0.01 and abs(k[1] - r.radar_lon) < 0.01]
        out[radar] = max(near)[1] if near else (round(r.radar_lat, 4), round(r.radar_lon, 4))
    return out


def hour_metrics(p: pd.DataFrame, threshold: float) -> dict:
    z = p.HGT.to_numpy() / FT
    t = p.temp.to_numpy()
    w = p.CLMR.to_numpy()
    lay = (z >= LAYER_FT[0]) & (z <= LAYER_FT[1]) & (t <= T_MAX_C)
    low = lay & (z <= LOW_CEILING_FT)
    wmax = float(np.nanmax(w[lay])) if lay.any() else 0.0
    wlow = float(np.nanmax(w[low])) if low.any() else 0.0
    t700 = p.loc[p.p == 700, "temp"]
    return {"slw": wmax >= threshold, "slw_low": wlow >= threshold, "w": wmax, "cold": bool(lay.any()),
            "t700": float(t700.iloc[0]) if len(t700) else np.nan}


def day_metrics(site: str, day: dt.date, loc: tuple[float, float], threshold: float) -> dict:
    t0, _ = pipeline.day_window(day)
    hours = [hour_metrics(hrrr.profile(t0.replace(tzinfo=None) + dt.timedelta(hours=k), loc[0], loc[1],
                                       box=pipeline.HRRR_BOX), threshold) for k in range(24)]
    return {"site": site, "day": str(day),
            "slw_hours": sum(h["slw"] for h in hours),
            "slw_hours_low": sum(h["slw_low"] for h in hours),
            "slw_max_gkg": max(h["w"] for h in hours),
            "cold_hours": sum(h["cold"] for h in hours),
            "t700_c": float(np.nanmean([h["t700"] for h in hours]))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", default=str(ROOT / "data" / "census"))
    ap.add_argument("--threshold", type=float, default=0.01, help="cloud water, g/kg")
    a = ap.parse_args()
    cen = pathlib.Path(a.census)
    rd = pd.read_csv(cen / "radar_days.csv")
    locs = cached_locations()
    rows = []
    for i, r in enumerate(rd[["site", "day"]].drop_duplicates().itertuples()):
        rows.append(day_metrics(r.site, dt.date.fromisoformat(r.day), locs[r.site], a.threshold))
        if i % 200 == 0:
            print(f"{i} radar-days", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(cen / "slw_days.csv", index=False)
    print(f"{len(out)} radar-days written; median slw_hours {out.slw_hours.median():.0f}, "
          f"share with none {np.mean(out.slw_hours == 0):.2f}")


if __name__ == "__main__":
    main()
