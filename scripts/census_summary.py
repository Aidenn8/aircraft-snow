"""Collect radar-day results into tables, select event days, and build the catalogue.

usage: census_summary.py [--events data/events] [--out data/census] [--q 0.05] [--alpha 0.005]

Day selection. Every tested radar-day has a screening score E (24 rotations) and a
gamma-tail p-value for it. The gamma tail is an approximation, so it is calibrated
empirically: on every day, each rotation in turn is treated as the real flight set and
scored against the other rotations in the same way, giving a null p-value; pooled over all
tested days these show how often the gamma tail produces a p-value this small by chance.
A day's calibrated p-value is the share of pooled null p-values at or below its own,
(1 + count) / (1 + N). Event days are those selected by Benjamini-Hochberg at FDR --q on
the calibrated p-values.

This calibration checks the tail approximation, not whether the real flights are
exchangeable with rotated ones (terrain, beam blockage and traffic corridors could break
that); null days in summer, where aircraft cannot seed, check that (run this script on
data/events_null).

Outputs (in --out):
  radar_days.csv    one row per radar-day: screens, tests, reversed diagnostic, p-values
  null_pvalues.npy  the pooled null p-values
  calibration.csv   share of pooled null p-values at or below nominal levels
  event_days.csv    radar-days selected by BH
  catalogue.csv     for event days, each flight's best streak match with p <= --alpha
"""
import argparse
import json
import pathlib

import numpy as np
import pandas as pd

from snow import match

ROOT = pathlib.Path(__file__).resolve().parents[1]


def screen_stats(s: dict, path: pathlib.Path) -> dict:
    """The day's screening statistics, with the null scores. Days run before the screening
    result was stored separately fall back to the final result, and days run before the null
    scores were stored recompute them from the saved streak table."""
    sk = s.get("streaks_screen") or s.get("streaks", {})
    if "E_null" not in sk and sk.get("n_streaks", 0) > 0:
        tp = path.with_name(path.stem + "_streaks.parquet")
        if tp.exists():
            sk = {**sk, **match.day_statistic(pd.read_parquet(tp)), "recomputed": True}
    return sk


def radar_days(ev_root: pathlib.Path) -> tuple[pd.DataFrame, dict]:
    rows, nulls = [], {}
    for path in sorted(ev_root.glob("*/*.json")):
        s = json.loads(path.read_text())
        r = {"site": s["site"], "day": s["day"], "screened_out": s.get("screened_out", "")}
        if not r["screened_out"]:
            sc = screen_stats(s, path)
            sk = s.get("streaks", {})
            rv = s.get("streaks_reverse", {})
            st = s.get("stats", {})
            r.update({"n_flights": s.get("n_flights"), "n_scans": s.get("n_scans"),
                      "n_streaks": sk.get("n_streaks"), "n_matched": sk.get("n_matched"),
                      "n_flights_matched": sk.get("n_flights_matched"),
                      "E_screen": sc.get("E"), "p_screen": sc.get("p_day"),
                      "p_gamma": sc.get("p_day_gamma", np.nan), "R_screen": len(sc.get("E_null", [])),
                      "recomputed": bool(sc.get("recomputed", False)),
                      "E": sk.get("E"), "p_day": sk.get("p_day"), "n_rotations": sk.get("n_rotations"),
                      "p_day_gamma": sk.get("p_day_gamma", np.nan),
                      "p_event": st.get("z", {}).get("p_cluster", np.nan),
                      "p_rev": rv.get("p_day"), "E_rev": rv.get("E"),
                      "cold_moist_hours": len(s.get("cold_moist", {}).get("hours", [])),
                      "seconds": s.get("seconds")})
            fdr = {f["threshold"]: f for f in sk.get("fdr", [])}
            if 0.005 in fdr:
                r.update({"n_flights_p005": fdr[0.005]["n_real"], "exp_false_p005": fdr[0.005]["expected_false"]})
            if len(sc.get("E_null", [])):
                nulls[(s["site"], s["day"])] = np.asarray(sc["E_null"], float)
        rows.append(r)
    return pd.DataFrame(rows), nulls


def loo_null_pvalues(E_null: np.ndarray) -> np.ndarray:
    """Gamma-tail p-value of each rotation scored against the other rotations of its day."""
    out = np.empty(len(E_null))
    for j in range(len(E_null)):
        out[j] = match.tail_pvalues(E_null[j], np.delete(E_null, j))["p_day_gamma"]
    return out


def calibrated(p: np.ndarray, null: np.ndarray) -> np.ndarray:
    """(1 + number of pooled null p-values <= p) / (1 + number of null p-values)."""
    null = np.sort(null[np.isfinite(null)])
    p = np.asarray(p, float)
    k = np.searchsorted(null, np.where(np.isfinite(p), p, np.inf), side="right")
    return np.where(np.isfinite(p), (1 + k) / (1 + len(null)), 1.0)


def benjamini_hochberg(p: np.ndarray, q: float) -> np.ndarray:
    """Boolean mask of rejected hypotheses (BH step-up at FDR q)."""
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return np.zeros(0, bool)
    order = np.argsort(p)
    thresh = q * (np.arange(1, n + 1) / n)
    ok = p[order] <= thresh
    k = np.nonzero(ok)[0].max() + 1 if ok.any() else 0
    mask = np.zeros(n, bool)
    mask[order[:k]] = True
    return mask


LEVELS = (0.1, 0.05, 0.01, 1e-3, 1e-4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default=str(ROOT / "data" / "events"))
    ap.add_argument("--out", default=str(ROOT / "data" / "census"))
    ap.add_argument("--q", type=float, default=0.05)
    ap.add_argument("--alpha", type=float, default=0.005)
    a = ap.parse_args()
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rd, nulls = radar_days(pathlib.Path(a.events))
    tested = rd[rd.screened_out == ""].copy()
    if tested.empty:
        rd.to_csv(out / "radar_days.csv", index=False)
        print(f"radar-days: {len(rd)}, none tested")
        return
    # days with no streak objects at all: E = 0 for the real set and every rotation
    tested["p_gamma"] = tested.p_gamma.where(tested.n_streaks.fillna(0) > 0, 1.0)
    null = np.concatenate([loo_null_pvalues(v) for v in nulls.values()]) if nulls else np.array([])
    np.save(out / "null_pvalues.npy", null)
    cal = pd.DataFrame({"level": LEVELS, "null_share": [float(np.mean(null <= L)) if len(null) else np.nan
                                                        for L in LEVELS], "n_null": len(null)})
    cal.to_csv(out / "calibration.csv", index=False)
    tested["p_cal"] = calibrated(tested.p_gamma.to_numpy(float), null) if len(null) else tested.p_gamma
    tested["event"] = benjamini_hochberg(tested.p_cal.to_numpy(float), a.q)
    rd = rd.merge(tested[["site", "day", "p_cal", "event"]], on=["site", "day"], how="left")
    rd["event"] = rd.event.fillna(False).astype(bool)
    rd.to_csv(out / "radar_days.csv", index=False)
    ev = tested[tested.event]
    ev.to_csv(out / "event_days.csv", index=False)
    n_rec = int(tested.recomputed.sum())
    print(f"radar-days: {len(rd)}, screened out {int((rd.screened_out != '').sum())}, tested {len(tested)}"
          + (f" ({n_rec} with null scores recomputed from saved tables)" if n_rec else ""))
    print(f"pooled null p-values: {len(null)}; share at or below nominal level:")
    print(cal.to_string(index=False))
    print(f"real days: share with screening rank p <= 0.05: {np.mean(tested.p_screen <= 0.05):.3f}; "
          f"gamma p <= 0.05: {np.mean(tested.p_gamma <= 0.05):.3f}; calibrated p <= 0.05: "
          f"{np.mean(tested.p_cal <= 0.05):.3f}; reversed rank p <= 0.05: "
          f"{np.mean(tested.p_rev.dropna() <= 0.05):.3f}")
    print(f"event days (BH q={a.q} on calibrated p): {len(ev)}")
    rows = []
    for _, d in ev.iterrows():
        path = pathlib.Path(a.events) / d.site / f"{d.day}_streaks.parquet"
        if not path.exists():
            continue
        S = pd.read_parquet(path)
        S = S[(S.coverage > 0) & (S.p_value <= a.alpha)]
        S = S.sort_values("p_value").drop_duplicates("flight")
        S["site"], S["day"] = d.site, d.day
        rows.append(S.drop(columns=["placebo_cov", "placebo_flight"], errors="ignore"))
    if rows:
        cat = pd.concat(rows, ignore_index=True)
        cat.to_csv(out / "catalogue.csv", index=False)
        print(f"catalogue: {len(cat)} flights with a streak at p <= {a.alpha}")


if __name__ == "__main__":
    main()
