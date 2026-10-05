"""Collect radar-day results into tables, select event days, and build the catalogue.

usage: census_summary.py [--events data/events] [--out data/census]
                         [--null-events data/events_null] [--q 0.05] [--alpha 0.01]

Why an empirical null. Each radar-day is scored by E, the evidence of its flights' streak
matches, against the same score for the flights rotated about the radar (24 rotations,
screening run; see snow/match.py). Two things make the rotation p-value unreliable in the
far tail that a census-wide false discovery rate needs: neighbouring rotations match the
same streaks, so there are only some 15-25 independent comparisons per day and no p-value
below about 1/20 means what it says; and on some days echo lines up with traffic corridors
for reasons other than seeding (fixed echoes, routing around storms), which rotations break.
Both act on most days alike. So, as is standard in large-scale testing when the theoretical
null is in doubt (Efron 2004, 2010), the null is estimated from the days themselves: each
day's standardised score z = (E - mean of the rotation scores) / their s.d. is computed, the
bulk of the z values across all tested days (most days have no seeding) is fitted by a
normal distribution using the median and the interquartile range, and each day's p-value is
its upper tail under that fitted null. Event days are those selected by Benjamini-Hochberg
at FDR --q on these p-values.

Checks: (1) warm-season null days (--null-events: tracks below 13 000 ft in July, where no
aircraft can seed), scored against the same fitted null; the share the selection would
pick is a direct false-alarm rate. (2) The theoretical selection (BH on the gamma-tail
p-value of each day's own rotations, refined with 200 rotations where the day was
confirmed) is reported alongside (column event_theory). (3) The null days' own z
distribution, centre versus tail, shows whether a normal null is too thin-tailed.

Outputs (in --out):
  radar_days.csv    one row per radar-day: screens, scores, p-values, selections
  event_days.csv    radar-days selected (empirical null, BH)
  null_fit.json     the fitted null (centre, width), and the null-day check
  catalogue.csv     for event days, each flight's best streak match with p <= --alpha
"""
import argparse
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

from snow import match

ROOT = pathlib.Path(__file__).resolve().parents[1]


def radar_days(ev_root: pathlib.Path) -> pd.DataFrame:
    rows = []
    for path in sorted(ev_root.glob("*/*.json")):
        if "_" in path.stem:          # auxiliary results (_dayswap, _shiftnull), not radar-days
            continue
        s = json.loads(path.read_text())
        r = {"site": s["site"], "day": s["day"], "screened_out": s.get("screened_out", "")}
        if not r["screened_out"]:
            sc = s.get("streaks_screen") or s.get("streaks", {})
            sk = s.get("streaks", {})
            rv = s.get("streaks_reverse", {})
            confirmed = sk.get("n_rotations", 0) > len(sc.get("E_null", [])) > 0
            p_rank = sk.get("p_day", 1.0)
            p_gamma = sk.get("p_day_gamma", np.nan)
            r.update({"n_flights": s.get("n_flights"), "n_scans": s.get("n_scans"),
                      "n_streaks": sk.get("n_streaks", 0), "n_matched": sk.get("n_matched", 0),
                      "n_flights_matched": sk.get("n_flights_matched", 0),
                      "E_screen": sc.get("E", 0.0), "z_screen": sc.get("z_day", np.nan),
                      "p_screen": sc.get("p_day", 1.0), "confirmed": confirmed,
                      "E": sk.get("E", 0.0), "n_rotations": sk.get("n_rotations", 0),
                      # the gamma tail is undefined with too few positive null scores: then the
                      # rank p-value is used (it is conservative there)
                      "p_theory": p_gamma if np.isfinite(p_gamma) else p_rank,
                      "p_rev": rv.get("p_day", np.nan),
                      "max_temp_c": s.get("max_temp_c"), "seconds": s.get("seconds")})
            fdr = {f["threshold"]: f for f in sk.get("fdr", [])}
            for tau in (0.005, 0.01, 0.05):
                if tau in fdr:
                    r[f"n_flights_p{tau}"] = fdr[tau]["n_real"]
                    r[f"exp_false_p{tau}"] = fdr[tau]["expected_false"]
        rows.append(r)
    return pd.DataFrame(rows)


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


def fit_null(z: np.ndarray) -> dict:
    """Normal null fitted to the bulk of the z values: centre = median, width = IQR / 1.349."""
    z = z[np.isfinite(z)]
    q25, q50, q75 = np.quantile(z, [0.25, 0.5, 0.75])
    return {"centre": float(q50), "width": float((q75 - q25) / 1.349), "n": int(len(z))}


def tail_p(z: np.ndarray, null: dict) -> np.ndarray:
    """Upper-tail p-value under the fitted null; days without a score get 1."""
    z = np.asarray(z, float)
    return np.where(np.isfinite(z), stats.norm.sf((z - null["centre"]) / null["width"]), 1.0)


def scored(rd: pd.DataFrame) -> pd.DataFrame:
    t = rd[rd.screened_out.fillna("") == ""].copy()
    # a day with no streak objects (or none matched by any set) has no score: no evidence
    t.loc[t.n_streaks.fillna(0) == 0, "z_screen"] = np.nan
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default=str(ROOT / "data" / "events"))
    ap.add_argument("--null-events", default=str(ROOT / "data" / "events_null"))
    ap.add_argument("--out", default=str(ROOT / "data" / "census"))
    ap.add_argument("--q", type=float, default=0.05)
    ap.add_argument("--alpha", type=float, default=0.01)
    a = ap.parse_args()
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rd = radar_days(pathlib.Path(a.events))
    tested = scored(rd)
    if tested.empty:
        rd.to_csv(out / "radar_days.csv", index=False)
        print(f"radar-days: {len(rd)}, none tested")
        return
    null = fit_null(tested.z_screen.to_numpy(float))
    tested["p_emp"] = tail_p(tested.z_screen.to_numpy(float), null)
    tested["event"] = benjamini_hochberg(tested.p_emp.to_numpy(float), a.q)
    tested["event_theory"] = benjamini_hochberg(tested.p_theory.fillna(1.0).to_numpy(float), a.q)
    cutoff = float(tested.loc[tested.event, "p_emp"].max()) if tested.event.any() else 0.0
    check = {"fit": null, "q": a.q, "p_cutoff": cutoff}
    nroot = pathlib.Path(a.null_events)
    if nroot.exists() and any(nroot.glob("*/*.json")):
        nd = scored(radar_days(nroot))
        nz = nd.z_screen.to_numpy(float)
        nd["p_emp"] = tail_p(nz, null)
        own = fit_null(nz) if np.isfinite(nz).sum() >= 10 else None
        check.update({"null_days": int(len(nd)), "null_selected": int((nd.p_emp <= cutoff).sum()) if cutoff else 0,
                      "null_share_p05": float(np.mean(nd.p_emp <= 0.05)),
                      "null_share_p01": float(np.mean(nd.p_emp <= 0.01)),
                      "null_own_fit": own})
        if own:
            # is the null's tail heavier than a normal fitted to its centre?
            zz = nz[np.isfinite(nz)]
            for k in (2.0, 3.0):
                check[f"null_beyond_{k:.0f}sd"] = {
                    "observed": int(np.sum(zz > own["centre"] + k * own["width"])),
                    "normal_expected": float(len(zz) * stats.norm.sf(k))}
        nd.to_csv(out / "null_days.csv", index=False)
    (out / "null_fit.json").write_text(json.dumps(check, indent=1, default=float))
    keep = ["site", "day", "z_screen", "p_emp", "event", "event_theory"]
    rd = rd.merge(tested[keep], on=["site", "day"], how="left")
    for c in ("event", "event_theory"):
        rd[c] = rd[c].fillna(False).astype(bool)
    rd.to_csv(out / "radar_days.csv", index=False)
    ev = tested[tested.event]
    ev.to_csv(out / "event_days.csv", index=False)
    print(f"radar-days: {len(rd)}, screened out {int((rd.screened_out.fillna('') != '').sum())}, "
          f"tested {len(tested)} ({int(np.isfinite(tested.z_screen).sum())} with a score), "
          f"confirmed with 200 rotations {int(tested.confirmed.sum())}")
    print(f"empirical null: centre {null['centre']:.2f}, width {null['width']:.2f} (theoretical: 0, 1)")
    print(f"event days (BH q={a.q}): {len(ev)} on the empirical null (p cutoff {cutoff:.2g}); "
          f"{int(tested.event_theory.sum())} on the theoretical gamma-tail p")
    if "null_days" in check:
        print(f"warm null days: {check['null_days']}; selected by the same cutoff: {check['null_selected']}; "
              f"share with p <= 0.05: {check['null_share_p05']:.3f}, <= 0.01: {check['null_share_p01']:.3f}")
        if check.get("null_own_fit"):
            print(f"  null days' own fit: centre {check['null_own_fit']['centre']:.2f}, "
                  f"width {check['null_own_fit']['width']:.2f}; beyond 2 / 3 widths: "
                  f"{check['null_beyond_2sd']['observed']} / {check['null_beyond_3sd']['observed']} "
                  f"(normal: {check['null_beyond_2sd']['normal_expected']:.1f} / "
                  f"{check['null_beyond_3sd']['normal_expected']:.1f})")
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
        print(f"catalogue: {len(cat)} flights with a streak at p <= {a.alpha} on event days")


if __name__ == "__main__":
    main()
