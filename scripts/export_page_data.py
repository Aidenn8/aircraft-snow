"""Collect the numbers the web page shows into site/page_data.json.

usage: export_page_data.py [--site KFTG] [--day 2025-01-25] [--census data/census]

From the Denver case: the matched aircraft (one row per flight, its best streak), the
day's statistics, the life-cycle composite with its 90% bootstrap band, the efficiency
table, and the dual-polarisation comparison. From the census, when present
(census_summary.py and census_analysis.py outputs): totals and per-radar counts against the
no-seeding expectation, the counts by hours of supercooled cloud, the efficiency table by
flight phase, the catalogue flights, and short notes on the band-filtered census and the
arrival/departure comparison. With --band, the band-filtered census_summary output.
"""
import argparse
import json
import pathlib
import warnings

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _clean(v):
    """JSON-safe: NaN and infinities become None, numpy scalars become Python numbers."""
    if isinstance(v, dict):
        return {k: _clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_clean(x) for x in v]
    if isinstance(v, (np.floating, float)):
        return float(v) if np.isfinite(v) else None
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.bool_):
        return bool(v)
    return v


def records(df: pd.DataFrame, cols, digits=3):
    out = df[cols].copy()
    for c in cols:
        if pd.api.types.is_float_dtype(out[c]):
            out[c] = out[c].round(6 if c.startswith("p_") or c.endswith("_p") else digits)
    return out.to_dict(orient="records")


def lifecycle(path: pathlib.Path, boot=2000, seed=1):
    L = np.load(path)
    off, ab, PF = L["offsets"], L["age_bins"], L["per_flight"]
    ages = 0.5 * (ab[:-1] + ab[1:]) / 60
    centre = PF[:, :, np.argmin(np.abs(off))]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)          # age bins with no samples
        side = np.nanmean(PF[:, :, np.abs(off) >= 5000], axis=2)
    C = centre - side
    rng = np.random.default_rng(seed)
    n = C.shape[0]
    bs = np.array([np.nanmean(C[rng.integers(0, n, n)], axis=0) for _ in range(boot)])
    lo, hi = np.nanquantile(bs, [0.05, 0.95], axis=0)
    return {"age_min": ages.tolist(), "mean_db": np.nanmean(C, axis=0).round(2).tolist(),
            "lo_db": lo.round(2).tolist(), "hi_db": hi.round(2).tolist(), "n_flights": int(n)}


def winters_text(days: pd.Series) -> str:
    """'three winters (November 2023 to March 2026) and March 2023': whole winters are those
    with at least four of their months in the census, the rest are named by month."""
    d = pd.to_datetime(days)
    y = pd.Series(np.where(d.dt.month >= 7, d.dt.year, d.dt.year - 1), index=d.index)
    months = d.dt.month.groupby(y).nunique()
    whole = sorted(int(w) for w, n in months.items() if n >= 4)
    words = {1: "one winter", 2: "two winters", 3: "three winters", 4: "four winters", 5: "five winters"}
    parts = [f"{words.get(len(whole), f'{len(whole)} winters')} (November {whole[0]} to March {whole[-1] + 1})"] if whole else []
    for w in sorted(int(w) for w in months.index if int(w) not in whole):
        sub = d[y == w]
        parts.append(f"{sub.min():%B %Y}" if sub.min().month == sub.max().month else f"{sub.min():%B} to {sub.max():%B %Y}")
    return " and ".join(parts)


def census_block(cen: pathlib.Path, band: pathlib.Path) -> dict:
    from snow import aircraft
    rd = pd.read_csv(cen / "radar_days.csv")
    rates = pd.read_csv(cen / "rates.csv")
    r = rates[rates.score == "census"]
    tot = r[(r.site == "all") & (r.winter == "all")].iloc[0]
    names = pd.read_csv(ROOT / "data" / "sites.csv").drop_duplicates("radar").set_index("radar").municipality
    sites = r[(r.site != "all") & (r.winter == "all")].copy()
    sites["name"] = sites.site.map(names)
    C = {"total": {"radar_days": int(len(rd)), "tested": int(tot.tested), "p01": int(tot.p01),
                   "expected_p01": float(tot.expected_p01), "p05": int(tot.p05), "expected_p05": float(tot.expected_p05)},
         "winters_text": winters_text(rd.day),
         "sites": records(sites, ["site", "name", "tested", "p01", "expected_p01", "p05", "expected_p05"], 2)}
    if (cen / "slw.csv").exists():
        sl = pd.read_csv(cen / "slw.csv")
        C["slw"] = records(sl[sl.stratum.astype(str) == "all"], ["slw_class", "tested", "p01", "p05"], 0)
    if (cen / "efficiency.csv").exists():
        C["efficiency"] = records(pd.read_csv(cen / "efficiency.csv"),
                                  ["engine_class", "phase", "tbin", "n", "k", "made_streak", "lo90", "hi90", "chance"], 4)
    if (cen / "seeders.csv").exists():
        S = pd.read_csv(cen / "seeders.csv")
        S["site_name"] = S.site.map(names)
        S["type_name"] = S.icao_type.map(aircraft.display_name)
        S["date"] = pd.to_datetime(S.time_utc, utc=True).dt.strftime("%-d %b %Y")
        S = S.sort_values("time_utc")
        C["seeders"] = records(S, ["date", "site", "site_name", "callsign", "icao_type", "type_name", "engine_class", "phase",
                                   "seed_alt_ft", "seed_temp_c", "age_min", "length_km", "snow_mm_h_max", "v_fall"], 2)
    cs = cen / "catalogue_summary.json"
    if cs.exists():
        k = json.loads(cs.read_text())
        C["catalogue"] = {"flights": k["catalogue_flights"], "days": k["catalogue_days"], "radars": k["catalogue_radars"],
                          "flights_p01_all": k["flights_p01_all"], "expected_false_p01": k["expected_false_p01"]}
    b = rates[(rates.score == "band-filtered") & (rates.site == "all") & (rates.winter == "all")]
    if len(b):
        b = b.iloc[0]
        C["band_note"] = (f"Snow bands lie along the wind, and so do runways, so a band can line up with a "
                          f"traffic corridor by chance more often than with rotated traffic. Scoring every day "
                          f"again without streaks that sit in families of parallel bands leaves "
                          f"<strong>{int(b.p01)}</strong> days at p&nbsp;&le;&nbsp;0.01 against "
                          f"{b.expected_p01:.1f} expected (and {int(b.p05)} at p&nbsp;&le;&nbsp;0.05 against "
                          f"{b.expected_p05:.1f}).")
    ph = cen / "phase.csv"
    if ph.exists():
        P = pd.read_csv(ph).set_index("comparison")
        if "by path length, temperature and day" in P.index:
            q = P.loc["by path length, temperature and day"]
            crude = P.loc["crude"]
            C["phase_note"] = (f"Arriving jets had about {crude['or']:.0f} times the odds of leaving a streak that "
                               f"departing ones had. Much of that is geometry: arrivals descend gently and fly level legs, so "
                               f"they spend about twice as long in a cloud layer as departures climbing through it. "
                               f"Compared at the same time in the layer, the same temperature and on the same day, "
                               f"arrivals still had <strong>{q['or']:.1f}</strong> times the odds of leaving a streak "
                               f"(90% interval {q['lo']:.1f}&ndash;{q['hi']:.1f}): consistent with the stronger cooling "
                               f"over a wing with gear and flaps out at low speed, which this census cannot separate "
                               f"from the slower speed itself.")
    return C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", default="KFTG")
    ap.add_argument("--day", default="2025-01-25")
    ap.add_argument("--census", default=str(ROOT / "data" / "census"))
    ap.add_argument("--band", default=str(ROOT / "data" / "census_band"))
    ap.add_argument("--out", default=str(ROOT / "site" / "page_data.json"))
    a = ap.parse_args()
    ev = ROOT / "data" / "events" / a.site
    s = json.loads((ev / f"{a.day}.json").read_text())
    desc = pd.read_csv(ev / f"{a.day}_described.csv")
    best = desc.sort_values("p_value").drop_duplicates("flight").copy()
    passed = pd.to_datetime(best.time_utc) - pd.to_timedelta(best.age_min, unit="min")
    best = best.assign(_passed=passed).sort_values("_passed")          # in the order they flew
    best["pass_utc"] = best._passed.dt.strftime("%H:%M")
    best["scan_utc"] = pd.to_datetime(best.time_utc).dt.strftime("%H:%M")
    fdr = {f["threshold"]: f for f in s["streaks"].get("fdr", [])}
    out = {
        "case": {
            "site": a.site, "day": a.day,
            "n_streak_objects": s["streaks"]["n_streaks"], "n_matched_objects": s["streaks"]["n_matched"],
            "E": s["streaks"]["E"], "E_null_median": s["streaks"]["E_null_median"],
            "E_null_q95": s["streaks"]["E_null_q95"], "p_day_gamma": s["streaks"]["p_day_gamma"],
            "n_rotations": s["streaks"]["n_rotations"],
            "flights_p01": fdr.get(0.01, {}).get("n_real"), "expected_false_p01": fdr.get(0.01, {}).get("expected_false"),
            "p_rev": s.get("streaks_reverse", {}).get("p_day"),
            "aircraft": records(best, ["callsign", "icao_type", "engine_class", "phase", "pass_utc", "scan_utc",
                                       "age_min", "seed_alt_ft", "seed_temp_c", "seed_dh_m", "v_fall", "length_km",
                                       "width_km", "z_max", "snow_mm_h_max", "snow_flux_t_h", "f1", "p_value"], 2),
        },
        "lifecycle": lifecycle(ev / f"{a.day}_lifecycle.npz"),
    }
    eff_path = ROOT / "figures" / "denver_efficiency.csv"
    if eff_path.exists():
        eff = pd.read_csv(eff_path)
        out["efficiency"] = records(eff, ["engine_class", "t_lo", "t_hi", "n", "made_streak", "share", "ci_lo", "ci_hi",
                                          "chance"], 3)
    dp = ev / f"{a.day}_dualpol_paired.csv"
    if dp.exists():
        out["dualpol"] = records(pd.read_csv(dp), ["dbz", "moment", "n_scans", "streak_median", "other_median",
                                                   "streak_higher", "sign_p"], 3)
    cen = pathlib.Path(a.census)
    if (cen / "rates.csv").exists():
        out["census"] = census_block(cen, pathlib.Path(a.band))
    pathlib.Path(a.out).write_text(json.dumps(_clean(out), separators=(",", ":"), allow_nan=False))
    print(a.out, f"{pathlib.Path(a.out).stat().st_size / 1e3:.0f} kB")


if __name__ == "__main__":
    main()
