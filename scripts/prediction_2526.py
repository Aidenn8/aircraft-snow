"""Score the prediction recorded in docs/methods.md before winter 2025-26 was analysed.

usage: prediction_2526.py [--census data/census]

The covariate slw_hours_cold (hours of HRRR supercooled cloud water at -10 C or colder over
the radar; slw_covariate.py) was chosen after seeing winters 2023-24 and 2024-25, so it is
tested on 2025-26 alone (November 2025 to March 2026, tested radar-days, calibrated p from
census_summary.py). The predictions:
  (1) days with no such cloud show no excess at p <= 0.05: the observed count is within the
      one-sided Poisson 95% range of 5% of them (P(X >= k) >= 0.05);
  (2) days with 4 or more such hours show an excess at p <= 0.05: P(X >= k) < 0.05 against
      5% of them;
  (3) the share at p <= 0.05 is larger with 4 or more hours than with none.
Writes CENSUS/prediction_2526.json and prints the outcome.
"""
import argparse
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", default=str(ROOT / "data" / "census"))
    a = ap.parse_args()
    cen = pathlib.Path(a.census)
    rd = pd.read_csv(cen / "radar_days.csv").merge(pd.read_csv(cen / "slw_days.csv"), on=["site", "day"])
    t = rd[rd.screened_out.fillna("") == ""]
    t = t[(t.day >= "2025-11-01") & (t.day <= "2026-03-31")]
    out = {"tested": int(len(t))}
    for name, sel in (("none", t.slw_hours_cold == 0), ("1-3 h", t.slw_hours_cold.between(1, 3)),
                      ("4+ h", t.slw_hours_cold >= 4)):
        g = t[sel]
        n, k = len(g), int((g.p_day <= 0.05).sum())
        out[name] = {"days": n, "p05": k, "expected": 0.05 * n, "share": k / n if n else np.nan,
                     "poisson_p": float(stats.poisson.sf(k - 1, 0.05 * n)) if n else np.nan,
                     "p01": int((g.p_day <= 0.01).sum())}
    out["prediction_1"] = bool(out["none"]["poisson_p"] >= 0.05)
    out["prediction_2"] = bool(out["4+ h"]["poisson_p"] < 0.05)
    out["prediction_3"] = bool(out["4+ h"]["share"] > out["none"]["share"])
    (cen / "prediction_2526.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
