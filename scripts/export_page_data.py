"""Collect the numbers the web page shows into site/page_data.json.

usage: export_page_data.py [--site KFTG] [--day 2025-01-25] [--census data/census]

From the Denver case: the matched aircraft (one row per flight, its best streak), the
day's statistics, the life-cycle composite with its 90% bootstrap band, the efficiency
table, and the dual-polarisation comparison. From the census, when present: rates per
radar and winter, the calibration tables and the seeding flights.
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", default="KFTG")
    ap.add_argument("--day", default="2025-01-25")
    ap.add_argument("--census", default=str(ROOT / "data" / "census"))
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
        out["census"] = {"rates": pd.read_csv(cen / "rates.csv").to_dict(orient="records")}
        for name in ("calibration", "seeders", "by_type", "efficiency"):
            if (cen / f"{name}.csv").exists():
                out["census"][name] = pd.read_csv(cen / f"{name}.csv").to_dict(orient="records")
        # false alarms: warm null days scored against the census's fitted null (census_summary.py)
        fit = cen / "null_fit.json"
        if fit.exists():
            nf = json.loads(fit.read_text())
            if "null_days" in nf:
                out["census"]["null"] = {"tested": nf["null_days"], "selected": nf["null_selected"],
                                         "cutoff": nf["p_cutoff"], "fit": nf["fit"]}
    pathlib.Path(a.out).write_text(json.dumps(_clean(out), separators=(",", ":"), allow_nan=False))
    print(a.out, f"{pathlib.Path(a.out).stat().st_size / 1e3:.0f} kB")


if __name__ == "__main__":
    main()
