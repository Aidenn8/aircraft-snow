"""Describe the matched streaks of a radar-day: aircraft, seeding altitude and temperature,
streak size and intensity, and a snowfall-rate estimate.

usage: describe_day.py SITE YYYY-MM-DD [--alpha 0.05]
Writes data/events/SITE/YYYY-MM-DD_described.csv

Snowfall rate (liquid equivalent, mm/h) from reflectivity with Ze = 110 S^2 (Wolfe &
Snider 2012, Cheyenne WY, a high-altitude Front Range WSR-88D); Ze-S relations vary
several-fold between storms (Rasmussen et al. 2003), so these are order-of-magnitude values.
"""
import argparse
import datetime as dt
import pathlib

import numpy as np
import pandas as pd

from snow import flights, match, pipeline, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]
ZS_ALPHA = 110.0


def snow_rate(dbz):
    return np.sqrt(10 ** (np.asarray(dbz, float) / 10.0) / ZS_ALPHA)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("day")
    ap.add_argument("--alpha", type=float, default=0.05)
    a = ap.parse_args()
    day = dt.date.fromisoformat(a.day)
    ev = ROOT / "data" / "events" / a.site
    S = pd.read_parquet(ev / f"{day}_streaks.parquet")
    S = S[(S.coverage > 0) & (S.p_value <= a.alpha)].copy()
    prep = pipeline.prepare(a.site, day)
    G = prep["grids"]
    pts, dr = pipeline.matching_points(prep)
    calls = pd.concat([flights.load_calls(day, a.site), flights.load_calls(day + dt.timedelta(days=1), a.site)])
    ft = flights.flight_table(prep["positions"], calls)
    xs = G.x0 + G.res * np.arange(G.z.shape[2])
    ys = G.y0 + G.res * np.arange(G.z.shape[1])
    rows = []
    for _, o in S.iterrows():
        k = int(np.argmin(np.abs(G.t - o.t)))
        z = np.where(np.isnan(G.z[k]), streaks.FLOOR_DBZ, G.z[k])
        obj, lab = streaks.find_streaks(streaks.GridScan(o.t, xs, ys, z, np.isfinite(G.z[k])))
        sx, sy = match.skeleton_xy(lab, int(o.label), xs, ys)
        # snowfall under the streak at this scan: sum of pixel rates x pixel area, as water
        # (1 mm over 1 km^2 = 1000 t); below 0 dBZ counts as nothing
        zpix = G.z[k][lab == int(o.label)]
        zpix = zpix[np.isfinite(zpix) & (zpix >= 0)]
        flux_t_h = float(snow_rate(zpix).sum() * (G.res / 1e3) ** 2 * 1000.0)
        mp = match.matched_points(sx, sy, pts, dr, o.t, int(o.flight), o.du, o.dv)
        mp = flights.interp_profile_to_points(mp, prep["profiles"], fields=("temp", "RH", "CLMR"))
        info = ft.loc[int(o.flight)]
        rows.append({
            "time_utc": pd.Timestamp(o.t, unit="s", tz="UTC"), "x_km": o.x / 1e3, "y_km": o.y / 1e3,
            "length_km": o.length_m / 1e3, "width_km": o.width_m / 1e3, "z_mean": o.z_mean,
            "z_max": o.z_max, "snow_mm_h_mean": float(snow_rate(o.z_mean)),
            "snow_mm_h_max": float(snow_rate(o.z_max)), "snow_flux_t_h": flux_t_h,
            "area_km2": o.area_km2,
            "f1": o.coverage, "coverage": o["cov"], "precision": o.precision, "p_value": o.p_value,
            "age_min": o.age_s / 60, "quiet_before": o.quiet_before,
            "callsign": info.callsign, "reg": info.reg, "icao_type": info.icao_type,
            "engine_class": info.engine_class, "phase": info.phase, "flight": int(o.flight),
            "seed_alt_ft": float(mp.alt_ft.median()) if len(mp) else np.nan,
            "seed_temp_c": float(mp.temp.median()) if len(mp) else np.nan,
            "seed_temp_min": float(mp.temp.min()) if len(mp) else np.nan,
            "seed_temp_max": float(mp.temp.max()) if len(mp) else np.nan,
            "seed_rh": float(mp.RH.median()) if len(mp) else np.nan,
            "seed_dh_m": float(mp.dh.median()) if len(mp) and "dh" in mp else np.nan,
            "v_fall": float(np.median(mp.dh / (o.t - mp.t))) if len(mp) and "dh" in mp else np.nan,
            "n_seed_points": len(mp),
        })
    out = pd.DataFrame(rows).sort_values("time_utc")
    path = ev / f"{day}_described.csv"
    out.to_csv(path, index=False)
    pd.set_option("display.width", 250)
    print(out[["time_utc", "callsign", "icao_type", "engine_class", "phase", "seed_alt_ft", "seed_temp_c",
               "seed_dh_m", "v_fall", "age_min", "length_km", "z_max", "snow_mm_h_max", "f1",
               "p_value"]].round(3).to_string())
    print(path)


if __name__ == "__main__":
    main()
