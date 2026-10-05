"""Copy the small result tables into results/ (tracked in git; data/ is not).

usage: export_results.py

results/census/        census_summary.py and census_analysis.py tables for the winters
results/census_null/   the same for the warm-season null days
results/injection/     injection-recovery trials (all runs of the final model)
results/denver/        the Denver case: streak descriptions, exposure table, drift check,
                       dual-polarisation comparison, day summary
Large per-streak tables with placebo arrays stay in data/.
"""
import json
import pathlib
import shutil

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FINAL_INJECTION = ("KFTG_2025-01-01_2025-01-02_2025-01-03_s7.parquet",)


def copy_tables(src: pathlib.Path, dst: pathlib.Path, names):
    dst.mkdir(parents=True, exist_ok=True)
    n = 0
    for name in names:
        p = src / name
        if p.exists():
            shutil.copy2(p, dst / name)
            n += 1
    return n


def main():
    tables = ("radar_days.csv", "event_days.csv", "catalogue.csv", "calibration.csv", "rates.csv",
              "seeders.csv", "by_type.csv", "efficiency.csv")
    print("census:", copy_tables(ROOT / "data" / "census", OUT / "census", tables))
    print("census_null:", copy_tables(ROOT / "data" / "census_null", OUT / "census_null", tables))
    inj = OUT / "injection"
    inj.mkdir(parents=True, exist_ok=True)
    for name in FINAL_INJECTION:
        p = ROOT / "data" / "injection" / name
        if p.exists():
            pd.read_parquet(p).to_csv(inj / name.replace(".parquet", ".csv"), index=False)
            print("injection:", name)
    ev = ROOT / "data" / "events" / "KFTG"
    den = OUT / "denver"
    n = copy_tables(ev, den, ("2025-01-25_described.csv", "2025-01-25_exposure.csv", "2025-01-25_drift.csv",
                              "2025-01-25_dualpol_paired.csv"))
    s = json.loads((ev / "2025-01-25.json").read_text())
    for k in ("streaks", "streaks_screen", "streaks_reverse"):
        s.get(k, {}).pop("E_null", None)
    (den / "2025-01-25_summary.json").write_text(json.dumps(s, indent=1, default=float))
    for name in ("2025-01-25_dayswap.json",):
        if (ev / name).exists():
            shutil.copy2(ev / name, den / name)
    print("denver:", n + 1)


if __name__ == "__main__":
    main()
