"""Collect radar-day results into tables, select event days, and build the catalogue.

usage: census_summary.py [--events data/events] [--q 0.10] [--alpha 0.005]

Outputs (data/census/):
  radar_days.csv   one row per radar-day: screens, both tests, reversed diagnostic
  event_days.csv   radar-days selected by Benjamini-Hochberg at --q on the streak test
  catalogue.csv    for event days, each flight's best streak match with p <= --alpha
"""
import argparse
import json
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]


def radar_days(ev_root: pathlib.Path) -> pd.DataFrame:
    rows = []
    for path in sorted(ev_root.glob("*/*.json")):
        s = json.loads(path.read_text())
        r = {"site": s["site"], "day": s["day"], "screened_out": s.get("screened_out", "")}
        if not r["screened_out"]:
            sk = s.get("streaks", {})
            rv = s.get("streaks_reverse", {})
            r.update({"n_flights": s.get("n_flights"), "n_scans": s.get("n_scans"),
                      "n_streaks": sk.get("n_streaks"), "n_matched": sk.get("n_matched"),
                      "n_flights_matched": sk.get("n_flights_matched"), "E": sk.get("E"),
                      "p_day": sk.get("p_day"), "n_rotations": sk.get("n_rotations"),
                      "p_event": s["stats"]["z"]["p_cluster"], "p_event_new": s["stats"]["new"]["p_cluster"],
                      "p_rev": rv.get("p_day"), "E_rev": rv.get("E"),
                      "cold_moist_hours": len(s.get("cold_moist", {}).get("hours", [])),
                      "seconds": s.get("seconds")})
        rows.append(r)
    return pd.DataFrame(rows)


def benjamini_hochberg(p: np.ndarray, q: float) -> np.ndarray:
    """Boolean mask of rejected hypotheses (BH step-up at FDR q)."""
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    thresh = q * (np.arange(1, n + 1) / n)
    ok = p[order] <= thresh
    k = np.nonzero(ok)[0].max() + 1 if ok.any() else 0
    mask = np.zeros(n, bool)
    mask[order[:k]] = True
    return mask


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default=str(ROOT / "data" / "events"))
    ap.add_argument("--q", type=float, default=0.10)
    ap.add_argument("--alpha", type=float, default=0.005)
    a = ap.parse_args()
    out = ROOT / "data" / "census"
    out.mkdir(parents=True, exist_ok=True)
    rd = radar_days(pathlib.Path(a.events))
    rd.to_csv(out / "radar_days.csv", index=False)
    tested = rd[rd.screened_out == ""].dropna(subset=["p_day"]).copy()
    tested["event"] = benjamini_hochberg(tested.p_day.to_numpy(), a.q)
    ev = tested[tested.event]
    ev.to_csv(out / "event_days.csv", index=False)
    print(f"radar-days: {len(rd)}, screened out {int((rd.screened_out != '').sum())}, tested {len(tested)}")
    print(f"event days (BH q={a.q}): {len(ev)}")
    print(f"share of tested days with p_day <= 0.05: {np.mean(tested.p_day <= 0.05):.3f}; "
          f"reversed test: {np.mean(tested.p_rev.dropna() <= 0.05):.3f}; "
          f"event test (z): {np.mean(tested.p_event <= 0.05):.3f}")
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
