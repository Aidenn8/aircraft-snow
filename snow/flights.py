"""Turn raw ADS-B positions into flights near an airport, with identity and environment.

A flight is a run of positions from one aircraft with no gap longer than GAP_S. Each
flight gets the callsign broadcast during it, the registry entry (type, engine class) and
a phase label (arrival, departure, overflight) from its altitude trend near the airport.
"""
from __future__ import annotations

import datetime as dt

import numpy as np
import pandas as pd

from . import adsb, aircraft, geo
from .weather import pressure_from_pressure_altitude

GAP_S = 20 * 60


def _day_positions(d: dt.date, site: str | None) -> pd.DataFrame:
    """One day's positions: the per-radar split file when available (and radius allows),
    else the continental file."""
    if site is not None:
        path = adsb.ROOT / "data" / "adsb_radar" / site / f"{d}.parquet"
        if path.exists():
            return pd.read_parquet(path)
    if not adsb.day_path(d).exists():
        raise FileNotFoundError(f"ADS-B for {d} not downloaded")
    return adsb.load_day(d)


def load_calls(d: dt.date, site: str | None = None) -> pd.DataFrame:
    if site is not None:
        path = adsb.ROOT / "data" / "adsb_radar" / site / f"{d}.calls.parquet"
        if path.exists():
            return pd.read_parquet(path)
    return adsb.load_calls(d)


def load_window(lat0: float, lon0: float, radius_km: float, start: dt.datetime,
                end: dt.datetime, max_alt_ft: float = 30_000, site: str | None = None) -> pd.DataFrame:
    """Airborne positions within radius_km of (lat0, lon0) between start and end. With
    `site`, per-radar split files (snow.adsb_split, radius 100 km) are used when present."""
    if site is not None and radius_km > 100.0:
        site = None
    days = sorted({(start + dt.timedelta(days=k)).date()
                   for k in range((end.date() - start.date()).days + 1)})
    parts = []
    for d in days:
        p = _day_positions(d, site)
        p = p[(p.t >= start.timestamp()) & (p.t <= end.timestamp())]
        dlat = radius_km / 111.0
        dlon = radius_km / (111.0 * np.cos(np.radians(lat0)))
        p = p[(p.lat > lat0 - dlat) & (p.lat < lat0 + dlat) & (p.lon > lon0 - dlon) & (p.lon < lon0 + dlon)]
        parts.append(p)
    pos = pd.concat(parts, ignore_index=True)
    pos = pos[(pos.alt != adsb.ALT_GROUND) & (pos.alt != adsb.ALT_UNKNOWN)]
    pos = pos[pos.alt.astype(np.int32) * 25 <= max_alt_ft]
    d = geo.haversine(pos.lat.to_numpy(), pos.lon.to_numpy(), lat0, lon0)
    pos = pos[d <= radius_km * 1000].copy()
    pos["alt_ft"] = pos.alt.astype(np.float64) * 25.0
    pos["p_hpa"] = pressure_from_pressure_altitude(pos["alt_ft"])
    x, y = geo.to_xy(pos.lat.to_numpy(), pos.lon.to_numpy(), lat0, lon0)
    pos["x"], pos["y"] = x, y
    return pos.sort_values(["hex", "t"]).reset_index(drop=True)


def segment(pos: pd.DataFrame, gap_s: float = GAP_S) -> pd.DataFrame:
    """Add a flight id: a new flight starts after a gap longer than gap_s."""
    pos = pos.sort_values(["hex", "t"]).copy()
    new = (pos.hex.diff() != 0) | (pos.t.diff() > gap_s)
    pos["flight"] = np.cumsum(new.to_numpy()) - 1
    return pos


def densify(points: pd.DataFrame, step_m: float = 250.0, max_gap_s: float = 60.0,
            columns=("x", "y", "t", "alt_ft", "p_hpa", "gs")) -> pd.DataFrame:
    """Linearly interpolate each flight's track so consecutive points are at most step_m
    apart (heatmap positions are 10-30 s apart, i.e. up to several km at airliner speeds).
    Pairs more than max_gap_s apart are not bridged. Keeps the flight column and adds
    `src`, the index of the original point each new point starts from."""
    p = points.sort_values(["flight", "t"])
    fl = p.flight.to_numpy()
    cols = {c: p[c].to_numpy(float) for c in columns if c in p}
    src = p.index.to_numpy()
    same = (fl[1:] == fl[:-1]) & (np.diff(cols["t"]) <= max_gap_s)
    dist = np.hypot(np.diff(cols["x"]), np.diff(cols["y"]))
    nseg = np.where(same, np.maximum(np.ceil(dist / step_m), 1), 1).astype(int)
    # each original point i (except a flight's last) spawns nseg[i] points at fractions k/nseg
    reps = np.r_[nseg, 1]
    i0 = np.repeat(np.arange(len(p)), reps)
    k = np.concatenate([np.arange(n) for n in reps])
    frac = np.where(np.r_[same, False][i0], k / np.repeat(reps, reps), 0.0)
    i1 = np.minimum(i0 + 1, len(p) - 1)
    out = {c: v[i0] + frac * (v[i1] - v[i0]) for c, v in cols.items()}
    out["flight"] = fl[i0]
    out["src"] = src[i0]
    return pd.DataFrame(out)


def flight_table(pos: pd.DataFrame, calls: pd.DataFrame | None = None) -> pd.DataFrame:
    """One row per flight: identity, time span, altitude trend and engine class."""
    g = pos.groupby("flight")
    f = g.agg(hex=("hex", "first"), t0=("t", "min"), t1=("t", "max"), n=("t", "size"),
              alt_first=("alt_ft", "first"), alt_last=("alt_ft", "last"),
              alt_min=("alt_ft", "min"), alt_max=("alt_ft", "max"))
    f["phase"] = np.select(
        [f.alt_last < f.alt_first - 2000, f.alt_last > f.alt_first + 2000],
        ["arrival", "departure"], "level")
    if calls is not None and len(calls):
        cs = []
        for fid, row in f.iterrows():
            c = calls[(calls.hex == row.hex) & (calls.t1 >= row.t0 - 600) & (calls.t0 <= row.t1 + 600)]
            c = c[c.callsign.str.len() > 0]
            cs.append(c.sort_values("n").callsign.iloc[-1] if len(c) else "")
        f["callsign"] = cs
    else:
        f["callsign"] = ""
    reg = aircraft.describe(f.hex.to_numpy())
    f = f.join(reg[["reg", "icao_type", "desc", "engine_class", "wtc", "year", "owner"]], on="hex")
    f["hexs"] = f.hex.map(lambda h: f"{h & 0xFFFFFF:06x}")
    return f


def interp_profile_to_points(points: pd.DataFrame, profiles: pd.DataFrame,
                             fields=("temp", "RH", "CLMR", "UGRD", "VGRD")) -> pd.DataFrame:
    """Interpolate hourly HRRR profiles (columns run, p, fields) to each point's time and
    pressure: linear in time between analyses and in log-pressure between levels.
    The points' own time column is `t` (epoch seconds), so no field may be called `t`."""
    if "t" in fields:
        raise ValueError("field name 't' would overwrite the point times")
    out = points.copy()
    runs = np.array(sorted(profiles.run.unique()))
    run_s = np.array([pd.Timestamp(r).timestamp() for r in runs])
    tables = {r: profiles[profiles.run == r].sort_values("p") for r in runs}
    lp_pts = np.log(points.p_hpa.to_numpy())
    tt = points.t.to_numpy(float)
    i1 = np.clip(np.searchsorted(run_s, tt), 1, len(runs) - 1)
    i0 = i1 - 1
    w = np.clip((tt - run_s[i0]) / (run_s[i1] - run_s[i0]), 0, 1)
    for fld in fields:
        v0 = np.full(len(points), np.nan)
        v1 = np.full(len(points), np.nan)
        for k, r in enumerate(runs):
            tab = tables[r]
            lp = np.log(tab.p.to_numpy(float))
            f = tab[fld].to_numpy(float)
            interp = np.interp(lp_pts, lp, f, left=np.nan, right=np.nan)
            v0 = np.where(i0 == k, interp, v0)
            v1 = np.where(i1 == k, interp, v1)
        out[fld] = (1 - w) * v0 + w * v1
    return out
