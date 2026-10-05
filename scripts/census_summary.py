"""Collect radar-day results into tables, select event days, and build the catalogue.

usage: census_summary.py [--events data/events] [--out data/census] [--q 0.05] [--alpha 0.005]

Day selection. Every tested radar-day has a screening score E (24 rotations) and a
gamma-tail p-value for it; days that screen at rank p <= 0.05 or gamma p <= 0.1 are re-run
with 200 rotations, which refines the same p-value (a sequential Monte Carlo test, Besag &
Clifford 1991). A day's p-value is the refined one if it was confirmed, else its screening
one. Event days are selected by Benjamini-Hochberg at FDR --q on these gamma-tail p-values.

The gamma tail is an approximation, so it is checked empirically: on every day, each
rotation in turn is treated as the real flight set and scored against the other rotations
in the same way, giving a null p-value; pooled over days (separately for 24- and
200-rotation runs) these show how often the gamma tail produces a p-value this small by
chance (calibration.csv). As a robustness check, each day also gets an empirically
calibrated p-value, the share of pooled null p-values at or below its own, (1 + count) /
(1 + N), and a second BH selection on those (column event_cal); it cannot resolve p-values
below 1/(1 + N).

This checks the tail approximation, not whether the real flights are exchangeable with
rotated ones. They are not always: on a warm null day the rotation test reached p = 0.01
(docs/methods.md), because rotation also breaks links between corridors and echo that are
not seeding. So a day selected on rotations is an event day only if it is confirmed by the
day-swap test (dayswap_null.py: other days' flights on this day's radar and wind), with
gamma-tail p <= --swap-alpha. This is an intersection-union rule: a confound that makes the
rotation p-value too small leaves the day-swap p-value valid (conservative). Selected days
without a day-swap result yet are reported as pending.

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


def radar_days(ev_root: pathlib.Path) -> tuple[pd.DataFrame, dict, dict]:
    rows, nulls24, nulls200 = [], {}, {}
    for path in sorted(ev_root.glob("*/*.json")):
        if "_" in path.stem:          # auxiliary results (_dayswap, _shiftnull), not radar-days
            continue
        s = json.loads(path.read_text())
        r = {"site": s["site"], "day": s["day"], "screened_out": s.get("screened_out", "")}
        if not r["screened_out"]:
            sc = screen_stats(s, path)
            sk = s.get("streaks", {})
            rv = s.get("streaks_reverse", {})
            st = s.get("stats", {})
            confirmed = sk.get("n_rotations", 0) > len(sc.get("E_null", [])) > 0
            r.update({"n_flights": s.get("n_flights"), "n_scans": s.get("n_scans"),
                      "n_streaks": sk.get("n_streaks"), "n_matched": sk.get("n_matched"),
                      "n_flights_matched": sk.get("n_flights_matched"),
                      "E_screen": sc.get("E"), "p_screen": sc.get("p_day"),
                      "p_gamma_screen": sc.get("p_day_gamma", np.nan), "R_screen": len(sc.get("E_null", [])),
                      "recomputed": bool(sc.get("recomputed", False)), "confirmed": confirmed,
                      "E": sk.get("E"), "p_day": sk.get("p_day"), "n_rotations": sk.get("n_rotations"),
                      "p_gamma_confirm": sk.get("p_day_gamma", np.nan) if confirmed else np.nan,
                      "p_event": st.get("z", {}).get("p_cluster", np.nan),
                      "p_rev": rv.get("p_day"), "E_rev": rv.get("E"),
                      "cold_moist_hours": len(s.get("cold_moist", {}).get("hours", [])),
                      "seconds": s.get("seconds")})
            sw = path.with_name(f"{path.stem}_dayswap.json")
            if sw.exists():
                ds = json.loads(sw.read_text())
                r.update({"p_swap": ds["stats"].get("p_day"), "p_swap_gamma": ds["stats"].get("p_day_gamma"),
                          "E_swap": ds["stats"].get("E"), "swap_sets": len(ds.get("offsets", []))})
            fdr = {f["threshold"]: f for f in sk.get("fdr", [])}
            if 0.005 in fdr:
                r.update({"n_flights_p005": fdr[0.005]["n_real"], "exp_false_p005": fdr[0.005]["expected_false"]})
            if len(sc.get("E_null", [])):
                nulls24[(s["site"], s["day"])] = np.asarray(sc["E_null"], float)
            if confirmed and len(sk.get("E_null", [])):
                nulls200[(s["site"], s["day"])] = np.asarray(sk["E_null"], float)
        rows.append(r)
    return pd.DataFrame(rows), nulls24, nulls200


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
    ap.add_argument("--swap-alpha", type=float, default=0.05)
    a = ap.parse_args()
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rd, nulls24, nulls200 = radar_days(pathlib.Path(a.events))
    tested = rd[rd.screened_out == ""].copy()
    if tested.empty:
        rd.to_csv(out / "radar_days.csv", index=False)
        print(f"radar-days: {len(rd)}, none tested")
        return
    # days with no streak objects at all: E = 0 for the real set and every rotation
    none = tested.n_streaks.fillna(0) == 0
    tested.loc[none, "p_gamma_screen"] = 1.0
    null24 = np.concatenate([loo_null_pvalues(v) for v in nulls24.values()]) if nulls24 else np.array([])
    null200 = np.concatenate([loo_null_pvalues(v) for v in nulls200.values()]) if nulls200 else np.array([])
    np.save(out / "null_pvalues.npy", null24)
    np.save(out / "null_pvalues_200.npy", null200)
    cal = pd.DataFrame([{"rotations": R, "level": L, "null_share": float(np.mean(nl <= L)) if len(nl) else np.nan,
                         "n_null": len(nl)} for R, nl in ((24, null24), (200, null200)) for L in LEVELS])
    cal.to_csv(out / "calibration.csv", index=False)
    conf = tested.confirmed.fillna(False).astype(bool)
    tested["p_final"] = np.where(conf, tested.p_gamma_confirm, tested.p_gamma_screen)
    tested["p_final"] = tested.p_final.fillna(1.0)
    ones = np.ones(len(tested))
    cal24 = calibrated(tested.p_gamma_screen.fillna(1.0).to_numpy(float), null24) if len(null24) else ones
    cal200 = calibrated(tested.p_gamma_confirm.fillna(1.0).to_numpy(float), null200) if len(null200) else cal24
    tested["p_cal"] = np.where(conf, cal200, cal24)
    tested["selected_rot"] = benjamini_hochberg(tested.p_final.to_numpy(float), a.q)
    tested["event_cal"] = benjamini_hochberg(tested.p_cal.to_numpy(float), a.q)
    if "p_swap_gamma" not in tested:
        tested["p_swap_gamma"] = np.nan
    tested["swap_pass"] = tested.p_swap_gamma <= a.swap_alpha
    tested["swap_pending"] = tested.selected_rot & tested.p_swap_gamma.isna()
    tested["event"] = tested.selected_rot & tested.swap_pass
    keep = ["site", "day", "p_final", "p_cal", "selected_rot", "swap_pass", "swap_pending", "event", "event_cal"]
    rd = rd.merge(tested[keep], on=["site", "day"], how="left")
    for c in ("selected_rot", "swap_pass", "swap_pending", "event", "event_cal"):
        rd[c] = rd[c].fillna(False).astype(bool)
    rd.to_csv(out / "radar_days.csv", index=False)
    ev = tested[tested.event]
    ev.to_csv(out / "event_days.csv", index=False)
    n_rec = int(tested.recomputed.sum())
    print(f"radar-days: {len(rd)}, screened out {int((rd.screened_out != '').sum())}, tested {len(tested)}, "
          f"confirmed with 200 rotations {int(conf.sum())}"
          + (f" ({n_rec} with null scores recomputed from saved tables)" if n_rec else ""))
    print("pooled null p-values (rotations scored as real): share at or below nominal level")
    print(cal.pivot(index="level", columns="rotations", values="null_share").to_string())
    print(f"real days: screening rank p <= 0.05: {np.mean(tested.p_screen <= 0.05):.3f}; "
          f"final gamma p <= 0.05: {np.mean(tested.p_final <= 0.05):.3f}; reversed rank p <= 0.05: "
          f"{np.mean(tested.p_rev.dropna() <= 0.05):.3f}")
    print(f"selected on rotations (BH q={a.q}): {int(tested.selected_rot.sum())} "
          f"({int(tested.event_cal.sum())} on empirically calibrated p); confirmed by day-swap "
          f"(gamma p <= {a.swap_alpha}): {len(ev)} event days; day-swap pending: {int(tested.swap_pending.sum())}")
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
