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

from . import adsb, epoch, event, flights, hrrr, match, radar, screen, sites, streaks
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


def site_location(site: str) -> tuple[float, float]:
    """Radar latitude/longitude from the NCEI station list."""
    return sites.location(site)


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


def prepare(site: str, day: dt.date, scans=None, alt_max_ft: float | None = None):
    """Everything the event test needs: grids, points with drift, the radar sweep header."""
    t0, t1 = day_window(day)
    if scans is None:
        scans = load_scans(site, t0 - dt.timedelta(minutes=30), t1 + dt.timedelta(minutes=60))
    if len(scans) < 20:
        raise RuntimeError(f"{site} {day}: only {len(scans)} radar scans")
    sw0 = scans[0][1]
    lat0, lon0 = sw0.radar_lat, sw0.radar_lon
    pos = flights.segment(flights.load_window(lat0, lon0, MAX_RANGE_M / 1000 + 5,
                                              t0 - dt.timedelta(minutes=5), t1, site=site))
    hi = ALT_RANGE_FT[1] if alt_max_ft is None else min(alt_max_ft, ALT_RANGE_FT[1])
    pts = pos[(pos.alt_ft >= ALT_RANGE_FT[0]) & (pos.alt_ft < hi)].copy()
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
                                 t1 + dt.timedelta(hours=1.5), max_alt_ft=60_000, site=site)
    grids = event.Grids.from_sweeps(scans, aircraft=allpos[["hex", "t", "x", "y", "gs"]])
    # tracks interpolated to 250 m spacing for streak matching; each interpolated point
    # inherits the drift of the reported position it starts from
    dense = flights.densify(pts)
    drift_dense = drift[dense.src.to_numpy()]
    return {"grids": grids, "points": pts, "drift": drift, "profiles": profs, "sweep0": sw0,
            "positions": pos, "t0": t0, "points_dense": dense, "drift_dense": drift_dense}


PERSISTENT_FRACTION = 0.4      # pixels inside streak objects in >= this share of scans
MIN_MATCH_GS_KT = 80.0
PERSISTENT_OVERLAP = 0.5       # objects mostly on such pixels are fixed echo, not snow
SCREEN_ANGLES = tuple(np.linspace(20.0, 340.0, 24))
CONFIRM_ANGLES = tuple(np.linspace(20.0, 340.0, 200))


def streak_objects_day(prep: dict):
    """Streak objects of every scan in the day, after removing fixed echoes: returns
    a list of (scan index, objects, label image) and the persistence map."""
    G = prep["grids"]
    xs = G.x0 + G.res * np.arange(G.z.shape[2])
    ys = G.y0 + G.res * np.arange(G.z.shape[1])
    t0 = prep["t0"].timestamp()
    found, count, n = [], np.zeros(G.z.shape[1:], np.int32), 0
    for k, T in enumerate(G.t):
        if not (t0 <= T <= t0 + 86400 + 3600):
            continue
        z = np.where(np.isnan(G.z[k]), streaks.FLOOR_DBZ, G.z[k])
        obj, lab = streaks.find_streaks(streaks.GridScan(T, xs, ys, z, np.isfinite(G.z[k])))
        n += 1
        count += lab > 0
        if not obj.empty:
            found.append((k, obj, lab))
    persist = count >= max(PERSISTENT_FRACTION * n, 3)
    kept = []
    for k, obj, lab in found:
        frac = obj.label.map(lambda L: float(persist[lab == L].mean()))
        o = obj[frac.to_numpy() < PERSISTENT_OVERLAP]
        if not o.empty:
            kept.append((k, o, lab))
    return kept, persist, xs, ys


def match_day(prep: dict, angles=SCREEN_ANGLES) -> pd.DataFrame:
    """Streak objects in every scan of the day (fixed echoes removed), each matched to its
    best new flight with rotation placebos."""
    G, pts, dr = prep["grids"], prep["points_dense"], prep["drift_dense"]
    # slow points (helicopters, hovering, slow general aviation: < 80 kt) are excluded: their
    # tracks loop over small areas, so their drifted "lines" are blobs that fit any echo
    fast = ~((pts.gs.to_numpy() >= 0) & (pts.gs.to_numpy() < MIN_MATCH_GS_KT * 10))
    pts, dr = pts[fast].reset_index(drop=True), dr[fast]
    kept, _, xs, ys = streak_objects_day(prep)
    rows = [match.match_streaks(G, G.t[k], lab, obj, xs, ys, pts, dr, angles=angles)
            for k, obj, lab in kept]
    rows = [r for r in rows if not r.empty]
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def drop_radar_cache(site: str, day: dt.date):
    """Delete this radar-day's cached Level III files (they can be fetched again)."""
    t0, t1 = day_window(day)
    site3 = site[1:]
    lo, hi = t0 - dt.timedelta(minutes=30), t1 + dt.timedelta(minutes=60)
    n = 0
    for d in {lo.date(), hi.date()}:
        for p in (radar.CACHE / "l3").glob(f"{site3}_N0B_{d:%Y_%m_%d}_*"):
            t = dt.datetime.strptime(p.name[-19:], "%Y_%m_%d_%H_%M_%S").replace(tzinfo=UTC)
            if lo <= t <= hi:
                p.unlink()
                n += 1
    return n


def run_day(site: str, day: dt.date, save: bool = True, force: bool = False,
            null_alt_max_ft: float | None = None, keep_radar: bool | None = None) -> dict:
    """Run both tests for one radar-day (see below). Unless keep_radar is True, the day's
    cached radar files are deleted afterwards when neither test reaches p <= 0.05."""
    summary = _run_day(site, day, save, force, null_alt_max_ft)
    interesting = (summary.get("streaks", {}).get("p_day", 1.0) <= 0.05
                   or summary.get("stats", {}).get(event.PRIMARY, {}).get("p_cluster", 1.0) <= 0.05)
    if keep_radar is False or (keep_radar is None and not interesting):
        drop_radar_cache(site, day)
    return summary


def _run_day(site: str, day: dt.date, save: bool = True, force: bool = False,
             null_alt_max_ft: float | None = None) -> dict:
    """Run both tests for one radar-day. With null_alt_max_ft set (null mode, for
    false-alarm calibration on warm days), only tracks below that altitude are used, the
    screens are skipped and results go to data/events_null/."""
    null = null_alt_max_ft is not None
    out_dir = (OUT.parent / "events_null" if null else OUT) / site
    summ_path = out_dir / f"{day:%Y-%m-%d}.json"
    if summ_path.exists() and not force:
        return json.loads(summ_path.read_text())
    tic = time.time()
    t0, _ = day_window(day)
    lat0, lon0 = site_location(site)
    cm = screen.cold_moist(lat0, lon0, t0)
    if not cm["pass"] and not null:
        summary = {"site": site, "day": f"{day}", "screened_out": "no cold moist layer",
                   "cold_moist": cm, "seconds": round(time.time() - tic, 1)}
        if save:
            out_dir.mkdir(parents=True, exist_ok=True)
            summ_path.write_text(json.dumps(summary, indent=1, default=float))
        return summary
    prep = prepare(site, day, alt_max_ft=null_alt_max_ft)
    ep = screen.echo_present(prep["grids"])
    if not ep["pass"] and not null:
        summary = {"site": site, "day": f"{day}", "screened_out": "no echo",
                   "cold_moist": cm, "echo": ep, "seconds": round(time.time() - tic, 1)}
        if save:
            out_dir.mkdir(parents=True, exist_ok=True)
            summ_path.write_text(json.dumps(summary, indent=1, default=float))
        return summary
    res = event.run(prep["grids"], prep["points"], prep["drift"], prep["t0"].timestamp(), 24)
    summary = {"site": site, "day": f"{day}", "n_scans": int(len(prep["grids"].t)),
               "n_points": int(len(prep["points"])),
               "n_flights": int(prep["points"].flight.nunique()), "stats": {},
               "cold_moist": cm, "echo": ep}
    zmaps = {}
    for stat in event.STATS + ("new",):
        ev = event.evaluate_clusters(res, stat)
        zmaps[stat] = ev.pop("z_map")
        ev.pop("placebo_masses")
        summary["stats"][stat] = ev
    # streak objects matched to single flights (match.py): screen with few rotations,
    # confirm promising days with many
    streak_table = match_day(prep, SCREEN_ANGLES)
    summary["streaks"] = match.day_statistic(streak_table)
    summary["streaks"]["n_rotations"] = len(SCREEN_ANGLES)
    if summary["streaks"]["p_day"] <= 0.05 or summary["stats"][event.PRIMARY]["p_cluster"] <= 0.05:
        streak_table = match_day(prep, CONFIRM_ANGLES)
        summary["streaks"] = match.day_statistic(streak_table)
        summary["streaks"]["n_rotations"] = len(CONFIRM_ANGLES)
        if not streak_table.empty:
            summary["streaks"]["fdr"] = match.empirical_fdr(streak_table).to_dict(orient="records")
    if save and not streak_table.empty:
        out_dir.mkdir(parents=True, exist_ok=True)
        st = streak_table.copy()
        st["placebo_cov"] = st.placebo_cov.map(lambda a: np.asarray(a, np.float32).tolist())
        st.to_parquet(out_dir / f"{day:%Y-%m-%d}_streaks.parquet")
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
