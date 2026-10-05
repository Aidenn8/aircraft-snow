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


def load_scans(site: str, t_a: dt.datetime, t_b: dt.datetime, rho: bool = True):
    """Level III lowest-tilt reflectivity scans (N0B) between t_a and t_b, each with the same
    scan's correlation coefficient (N0C) attached as Sweep.rho when available."""
    site3 = site[1:] if len(site) == 4 else site
    files, rho_files = [], {}
    d = t_a.date()
    while d <= t_b.date():
        files += [f for f in radar.list_level3(site3, "N0B", d) if t_a <= f[0] <= t_b]
        if rho:
            rho_files.update({f[0]: f[1] for f in radar.list_level3(site3, "N0C", d) if t_a <= f[0] <= t_b})
        d += dt.timedelta(days=1)
    scans = []
    for t, key, _ in files:
        try:
            sw = radar.read_level3(radar.fetch(radar.L3_BUCKET, key))
        except Exception as e:      # an occasional corrupt file should not stop the day
            print(f"  skip {key}: {e!r}")
            continue
        if t in rho_files:
            try:
                sw.rho = radar.read_level3(radar.fetch(radar.L3_BUCKET, rho_files[t]), moment="RHO")
            except Exception as e:
                print(f"  no RHO for {key}: {e!r}")
        scans.append((t.timestamp(), sw))
    return scans


def track_points(pos: pd.DataFrame, profs: pd.DataFrame, sw0, alt_max_ft: float | None = None):
    """Track points that can be candidates (altitude band, range, not far below the beam),
    with their drift, height above the beam, beam half-width and HRRR temperature, and the
    same interpolated to 250 m spacing (each interpolated point inherits the drift of the
    reported position it starts from)."""
    hi = ALT_RANGE_FT[1] if alt_max_ft is None else min(alt_max_ft, ALT_RANGE_FT[1])
    pts = pos[(pos.alt_ft >= ALT_RANGE_FT[0]) & (pos.alt_ft < hi)].copy()
    rng = np.hypot(pts.x, pts.y).to_numpy()
    h_beam, _ = beam_height_ground_range(rng, sw0.elevation, sw0.radar_height)
    # pressure altitude ~ altitude above sea level to within a few hundred feet (a loose
    # preselection; the fall-speed test in match.py uses geometric height)
    alt_m = pts.alt_ft.to_numpy() * 0.3048
    keep = (h_beam <= alt_m + BEAM_MARGIN_M) & (rng >= MIN_RANGE_M) & (rng <= MAX_RANGE_M)
    pts = pts[keep].reset_index(drop=True)
    drift = epoch.drift_velocity(pts, profs, sw0.radar_height, sw0.elevation)
    # height above the beam and beam half-width, for the fall-speed test in match.py
    pts["dh"], pts["bw"] = epoch.height_above_beam(pts, profs, sw0.radar_height, sw0.elevation,
                                                   match.BEAM_HALF_WIDTH_DEG)
    # HRRR temperature at each point (time- and log-pressure-interpolated)
    pts["temp_c"] = flights.interp_profile_to_points(pts, profs, fields=("temp",))["temp"].to_numpy()
    dense = flights.densify(pts, columns=("x", "y", "t", "alt_ft", "p_hpa", "gs", "dh", "bw", "temp_c"))
    drift_dense = drift[dense.src.to_numpy()]
    return pts, drift, dense, drift_dense


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
    runs = [t0.replace(tzinfo=None) + dt.timedelta(hours=k) for k in range(-1, 26)]
    profs = pd.concat([hrrr.profile(r, lat0, lon0, box=HRRR_BOX) for r in runs], ignore_index=True)
    pts, drift, dense, drift_dense = track_points(pos, profs, sw0, alt_max_ft)
    # every aircraft at any altitude (including ones on the ground) is masked out of the radar
    allpos = flights.load_window(lat0, lon0, MAX_RANGE_M / 1000 + 15, t0 - dt.timedelta(hours=1),
                                 t1 + dt.timedelta(hours=1.5), max_alt_ft=60_000, site=site)
    grids = event.Grids.from_sweeps(scans, aircraft=allpos[["hex", "t", "x", "y", "gs"]])
    return {"grids": grids, "points": pts, "drift": drift, "profiles": profs, "sweep0": sw0,
            "positions": pos, "t0": t0, "points_dense": dense, "drift_dense": drift_dense}


PERSISTENT_FRACTION = 0.4      # pixels inside streak objects in >= this share of scans
MIN_MATCH_GS_KT = 80.0
PERSISTENT_OVERLAP = 0.5       # objects mostly on such pixels are fixed echo, not snow
SCREEN_ANGLES = tuple(np.linspace(20.0, 340.0, 24))
SCREEN_P = 0.05                # confirm a day whose screening rank p-value is <= this
SCREEN_P_GAMMA = 0.1           # ... or whose gamma-tail p-value is <= this
CONFIRM_ANGLES = tuple(np.linspace(20.0, 340.0, 200))


def streak_objects_day(prep: dict):
    """Streak objects of every scan in the day, after removing fixed echoes: returns
    a list of (scan index, objects, label image) and the persistence map. The result is
    cached in prep (it depends only on the radar)."""
    if "_streak_objects" in prep:
        return prep["_streak_objects"]
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
    prep["_streak_objects"] = (kept, persist, xs, ys)
    return kept, persist, xs, ys


def _fast_only(pts: pd.DataFrame, dr: np.ndarray):
    gs = pts.gs.to_numpy()
    fast = ~((gs >= 0) & (gs < MIN_MATCH_GS_KT * 10))
    return pts[fast].reset_index(drop=True), dr[fast]


def matching_points(prep: dict):
    """Interpolated track points used for streak matching, and their drift. Slow points
    (helicopters, hovering, slow general aviation: < MIN_MATCH_GS_KT) are excluded: their
    tracks loop over small areas, so their drifted "lines" are blobs that fit any echo."""
    return _fast_only(prep["points_dense"], prep["drift_dense"])


DAYSWAP_OFFSETS = tuple(k for k in range(-14, 15) if k != 0)


def dayswap_sets(prep: dict, site: str, day: dt.date, offsets=DAYSWAP_OFFSETS,
                 alt_max_ft: float | None = None):
    """Placebo flight sets from other days: the flights of day + k (same hours), moved in time
    by -k days onto this day, and prepared exactly like the real flights with this day's HRRR
    (drift, height above the beam, temperature). They fly the same corridors and runway
    patterns through the same weather but cannot have seeded this day's cloud. Days without
    ADS-B are skipped. Returns a list of (points, drift) and the offsets used."""
    t0, t1 = day_window(day)
    sw0, profs = prep["sweep0"], prep["profiles"]
    sets, used = [], []
    for k in offsets:
        d = day + dt.timedelta(days=k)
        if not (adsb.day_path(d).exists() and adsb.day_path(d + dt.timedelta(days=1)).exists()):
            continue
        shift = dt.timedelta(days=k)
        pos = flights.load_window(sw0.radar_lat, sw0.radar_lon, MAX_RANGE_M / 1000 + 5,
                                  t0 + shift - dt.timedelta(minutes=5), t1 + shift, site=site)
        pos = flights.segment(pos)
        pos["t"] = pos.t - k * 86400
        _, _, dense, drift_dense = track_points(pos, profs, sw0, alt_max_ft)
        p, d_ = _fast_only(dense, drift_dense)
        p["flight"] = p.flight + (1000 + k) * 1_000_000        # distinct from real flight ids
        sets.append((p, d_))
        used.append(k)
    return sets, used


def screen_passes(day_stat: dict) -> bool:
    """Whether a screened day goes on to confirmation with CONFIRM_ANGLES rotations."""
    pg = day_stat.get("p_day_gamma", np.nan)
    return bool(day_stat.get("p_day", 1.0) <= SCREEN_P or (np.isfinite(pg) and pg <= SCREEN_P_GAMMA))


def match_day(prep: dict, angles=SCREEN_ANGLES, reverse: bool = False,
              max_temp_c=match.T_SEED_MAX_C, placebo_sets=None) -> pd.DataFrame:
    """Streak objects in every scan of the day (fixed echoes removed), each matched to its
    best new flight with rotation placebos (reverse=True: the time-reversed diagnostic;
    max_temp_c=None: null mode, warm points only; placebo_sets: other flight sets instead of
    rotations, see dayswap_sets)."""
    G = prep["grids"]
    pts, dr = matching_points(prep)
    kept, _, xs, ys = streak_objects_day(prep)
    rows = [match.match_streaks(G, G.t[k], lab, obj, xs, ys, pts, dr, angles=angles, reverse=reverse,
                                max_temp_c=max_temp_c, placebo_sets=placebo_sets)
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
        for p in [*(radar.CACHE / "l3").glob(f"{site3}_N0B_{d:%Y_%m_%d}_*"),
                  *(radar.CACHE / "l3").glob(f"{site3}_N0C_{d:%Y_%m_%d}_*")]:
            t = dt.datetime.strptime(p.name[-19:], "%Y_%m_%d_%H_%M_%S").replace(tzinfo=UTC)
            if lo <= t <= hi:
                p.unlink()
                n += 1
    return n


def run_day(site: str, day: dt.date, save: bool = True, force: bool = False,
            null_alt_max_ft: float | None = None, keep_radar: bool | None = None,
            event_test: bool = True, reverse_all: bool = False) -> dict:
    """Run the tests for one radar-day (see below). Unless keep_radar is True, the day's
    cached radar files are deleted afterwards when no test reaches p <= 0.05."""
    summary = _run_day(site, day, save, force, null_alt_max_ft, event_test, reverse_all)
    interesting = (summary.get("streaks", {}).get("p_day", 1.0) <= 0.05
                   or summary.get("stats", {}).get(event.PRIMARY, {}).get("p_cluster", 1.0) <= 0.05)
    if keep_radar is False or (keep_radar is None and not interesting):
        drop_radar_cache(site, day)
    return summary


def _run_day(site: str, day: dt.date, save: bool = True, force: bool = False,
             null_alt_max_ft: float | None = None, event_test: bool = True,
             reverse_all: bool = False) -> dict:
    """Run the tests for one radar-day. With null_alt_max_ft set (null mode, for
    false-alarm calibration on warm days), only tracks below that altitude are used, the
    screens are skipped and results go to data/events_null/. event_test=False skips the
    parcel-following event test (event.py), which is not calibrated on winter days (see
    docs/methods.md); the census runs without it. The time-reversed diagnostic runs on
    candidate days (those sent to confirmation), or on every day with reverse_all."""
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
    summary = {"site": site, "day": f"{day}", "n_scans": int(len(prep["grids"].t)),
               "n_points": int(len(prep["points"])),
               "n_flights": int(prep["points"].flight.nunique()), "stats": {},
               "cold_moist": cm, "echo": ep}
    res, zmaps = None, {}
    if event_test:
        res = event.run(prep["grids"], prep["points"], prep["drift"], prep["t0"].timestamp(), 24)
        for stat in event.STATS + ("new",):
            ev = event.evaluate_clusters(res, stat)
            zmaps[stat] = ev.pop("z_map")
            ev.pop("placebo_masses")
            summary["stats"][stat] = ev
    # streak objects matched to single flights (match.py): screen with few rotations,
    # confirm promising days with many. With 24 rotations the rank p-value cannot go below
    # 0.04, so the smooth gamma-tail estimate also lets a borderline day through.
    tmax = None if null else match.T_SEED_MAX_C
    summary["max_temp_c"] = tmax
    streak_table = match_day(prep, SCREEN_ANGLES, max_temp_c=tmax)
    summary["streaks"] = match.day_statistic(streak_table, SCREEN_ANGLES)
    summary["streaks"]["n_rotations"] = len(SCREEN_ANGLES)
    # the screening result is kept for every day: it is the one statistic computed the same
    # way on all days, so census-level selection uses it (census_summary.py)
    summary["streaks_screen"] = dict(summary["streaks"])
    candidate = screen_passes(summary["streaks"])
    if candidate:
        streak_table = match_day(prep, CONFIRM_ANGLES, max_temp_c=tmax)
        summary["streaks"] = match.day_statistic(streak_table, CONFIRM_ANGLES)
        summary["streaks"]["n_rotations"] = len(CONFIRM_ANGLES)
        if not streak_table.empty:
            summary["streaks"]["fdr"] = match.empirical_fdr(streak_table, angles=CONFIRM_ANGLES).to_dict(orient="records")
    # time-reversed diagnostic, on candidate days only (it costs as much as the screen):
    # aircraft that pass 10-50 min AFTER each scan cannot have made its echoes; their excess
    # over rotations measures alignment with traffic corridors that is not caused by seeding
    if candidate or reverse_all:
        summary["streaks_reverse"] = match.day_statistic(match_day(prep, SCREEN_ANGLES, reverse=True,
                                                                   max_temp_c=tmax), SCREEN_ANGLES)
    if save and not streak_table.empty:
        out_dir.mkdir(parents=True, exist_ok=True)
        st = streak_table.copy()
        st["placebo_cov"] = st.placebo_cov.map(lambda a: np.asarray(a, np.float32).tolist())
        st.to_parquet(out_dir / f"{day:%Y-%m-%d}_streaks.parquet")
    summary["seconds"] = round(time.time() - tic, 1)
    if save and res is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out_dir / f"{day:%Y-%m-%d}_acc.npz", acc_real=res["acc_real"],
                            acc_placebo=res["acc_placebo"], new_real=res["new_real"],
                            new_placebo=res["new_placebo"], band_edges=res["band_edges"],
                            angles=res["angles"], t_start=res["t_start"],
                            **{f"zmap_{k}": v for k, v in zmaps.items()})
    if save:
        out_dir.mkdir(parents=True, exist_ok=True)
        summ_path.write_text(json.dumps(summary, indent=1, default=float))
    return summary
