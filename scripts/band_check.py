"""Does a day's evidence survive without streaks that belong to families of parallel bands?

usage: band_check.py SITE YYYY-MM-DD [SITE YYYY-MM-DD ...]
       band_check.py --census data/census --p 0.01

Rotating flights about the radar changes their orientation. Runways are chosen by the wind
and snow bands line up with it, so on a banded day the real arrival paths run along the
bands and rotated copies do not: band pieces match the real flights more often than the
rotated ones, and the rotation null is too easy to beat. Whether a streak object sits in a
family of parallel ones (streaks.parallel_counts: other objects in the same scan within
15 deg of its orientation and 60 km) depends on the radar echo alone, not on any flight, so
dropping such objects from the real and the rotated sets alike leaves the test valid and
removes the scenes where orientation coincidences happen. A busy arrival pattern can leave
a family of parallel aircraft streaks too (downwind legs flown again and again), so the
cost is some power; several thresholds are reported.

For each day the saved streak table (all 200 rotations, data/events/SITE/DAY_streaks.parquet)
is re-scored with match.day_statistic: first with every streak (it must reproduce the
saved E), then without objects with at least 8, 5 and 3 parallel neighbours. The p-value
is the rank of the real E among the rotations, each scored as if real against the sets at
least 20 deg away; z as in census_summary (s.d. floored at 0.7).
Writes data/events/SITE/DAY_bandcheck.json and, with --census, CENSUS/band_check.csv.
"""
import argparse
import datetime as dt
import json
import pathlib

import numpy as np
import pandas as pd

from snow import match, pipeline, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]
THRESHOLDS = (8, 5, 3)
SD_FLOOR = 0.7


def parallel_table(prep: dict) -> pd.DataFrame:
    """(t, label, parallel) for every streak object kept for the day."""
    kept, _, _, _ = pipeline.streak_objects_day(prep)
    rows = []
    for k, obj, lab in kept:
        par = streaks.parallel_counts(obj, streaks.orientations(lab, obj.label))
        rows.append(pd.DataFrame({"t": obj.t.to_numpy(float), "label": obj.label.to_numpy(int), "parallel": par}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["t", "label", "parallel"])


def score(table: pd.DataFrame, angles) -> dict:
    if table.empty:
        return {"E": 0.0, "p_rank": 1.0, "z": 0.0, "n_streaks": 0, "n_flights_matched": 0}
    st = match.day_statistic(table, angles)
    E_all = np.array([st["E"]] + list(st.get("E_null", [])), float)
    comp = match.comparable_sets(angles)
    ref = E_all[comp[0]]
    z = (E_all[0] - ref.mean()) / max(ref.std(ddof=1), SD_FLOOR)
    return {"E": float(st["E"]), "p_rank": float(st["p_day"]), "z": float(z),
            "E_null_median": float(st.get("E_null_median", 0.0)), "E_null_q95": float(st.get("E_null_q95", 0.0)),
            "n_streaks": int(len(table)), "n_flights_matched": int(st.get("n_flights_matched", 0))}


def check_day(site: str, day: dt.date) -> dict:
    ev = ROOT / "data" / "events" / site
    S = pd.read_parquet(ev / f"{day}_streaks.parquet")
    saved = json.loads((ev / f"{day}.json").read_text())["streaks"]
    R = len(S.placebo_cov.iloc[0])
    angles = {len(pipeline.CONFIRM_ANGLES): pipeline.CONFIRM_ANGLES,
              len(pipeline.SCREEN_ANGLES): pipeline.SCREEN_ANGLES}[R]
    if "parallel" not in S:                  # tables saved before the census recorded it
        par = parallel_table(pipeline.prepare(site, day))
        S = S.assign(_t=S.t.round(3)).merge(par.assign(_t=par.t.round(3)).drop(columns="t"),
                                            on=["_t", "label"], how="left").drop(columns="_t")
    missing = int(S.parallel.isna().sum())
    if missing:
        raise RuntimeError(f"{site} {day}: {missing} saved streaks not found among the day's objects")
    out = {"site": site, "day": str(day), "rotations": R, "all": score(S, angles)}
    if not np.isclose(out["all"]["E"], saved["E"], rtol=1e-6, atol=1e-6):
        raise RuntimeError(f"{site} {day}: re-scored E {out['all']['E']:.4f} != saved {saved['E']:.4f}")
    m = (S.coverage > 0) & (S.p_value <= 0.01)
    out["share_banded_objects"] = {str(n): float(np.mean(S.parallel >= n)) for n in THRESHOLDS}
    out["share_banded_matches"] = {str(n): float(np.mean(S.parallel[m] >= n)) if m.any() else np.nan
                                   for n in THRESHOLDS}
    for n in THRESHOLDS:
        out[f"max_parallel_{n - 1}"] = score(S[S.parallel < n].reset_index(drop=True), angles)
    (ev / f"{day}_bandcheck.json").write_text(json.dumps(out, indent=1, default=float))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("days", nargs="*", help="SITE YYYY-MM-DD pairs")
    ap.add_argument("--census", default=None)
    ap.add_argument("--p", type=float, default=0.01, help="with --census: days at or below this day p")
    a = ap.parse_args()
    pairs = list(zip(a.days[::2], a.days[1::2]))
    if a.census:
        rd = pd.read_csv(pathlib.Path(a.census) / "radar_days.csv")
        sel = rd[rd.screened_out.isna() & (rd.p_day <= a.p)].sort_values("p_day")
        pairs += list(zip(sel.site, sel.day))
    rows = []
    for site, d in pairs:
        try:
            o = check_day(site, dt.date.fromisoformat(d))
        except Exception as e:                  # noqa: BLE001 - report and continue
            print(f"{site} {d}: FAILED {e!r}", flush=True)
            continue
        line = [f"{site} {d}"] + [f"{k}: E {o[k]['E']:5.1f} p {o[k]['p_rank']:.3f} z {o[k]['z']:5.2f} n {o[k]['n_flights_matched']:2d}"
                                  for k in ["all"] + [f"max_parallel_{n - 1}" for n in THRESHOLDS]]
        print(" | ".join(line), flush=True)
        row = {"site": site, "day": d, "share_banded_matches_3": o["share_banded_matches"]["3"]}
        for k in ["all"] + [f"max_parallel_{n - 1}" for n in THRESHOLDS]:
            row.update({f"{k}_{c}": o[k][c] for c in ("E", "p_rank", "z", "n_flights_matched")})
        rows.append(row)
    if a.census and rows:
        pd.DataFrame(rows).to_csv(pathlib.Path(a.census) / "band_check.csv", index=False)


if __name__ == "__main__":
    main()
