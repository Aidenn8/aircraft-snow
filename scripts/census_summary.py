"""Collect radar-day results into tables, estimate how many days had aircraft snow, and
select and rank candidate days.

usage: census_summary.py [--events data/events] [--out data/census]
                         [--null-events data/events_null] [--q 0.05] [--alpha 0.01] [--top 30]

The day's p-value. Each radar-day's evidence E (its flights' streak matches) is compared
with E for the flights rotated about the radar (snow/match.py: every set scored against
sets at least 20 deg away, p-values floored at the coarsest resolution). The p-value is the
gamma-tail estimate from the 200-rotation run where the day was confirmed, else from the
24-rotation screen; with too few positive null scores for a gamma fit, the rank p-value.
These are valid under the rotation null (checked by simulation and on warm null days),
but a day has only some 15-25 independent rotations, so p-values far below ~0.001 rest on
the gamma extrapolation, and naming individual days at a census-wide false discovery rate
needs more than that across hundreds of days.

So the census answers two questions. (1) How many days had aircraft snow? Under no
seeding the day p-values are uniform, so the days at or below p = 0.01 (or 0.05) beyond the
1% (5%) the null predicts are event days: for rare, strong events this excess is the sharp
estimate. Storey's estimate of all event days, 1 - pi0 with pi0 = #(p > 0.5) / (N / 2), is
also given, with a bootstrap interval, but it is too noisy to resolve a few per cent. (2) Which days? Benjamini-Hochberg at --q, Storey
q-values, and a ranked list of the --top candidates for inspection. Warm-season null days
(--null-events) are summarised the same way: their estimated event share should be ~0.
A robustness check fits an empirical null (median, IQR; Efron 2004) to the probit of the
day p-values of days with any evidence (p < 0.99), for both sets of days.

Outputs (in --out):
  radar_days.csv    one row per radar-day: screens, scores, p-values, selections, q-values
  event_days.csv    radar-days selected by BH at --q
  candidates.csv    the --top days ranked by p-value
  summary.json      event-share estimates, selections, null-day check, empirical-null fits
  catalogue.csv     for selected days, each flight's best streak match with p <= --alpha
"""
import argparse
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[1]
LAMBDA = 0.5


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
                      "p_rank": p_rank,
                      "p_day": float(p_gamma) if np.isfinite(p_gamma) else float(p_rank),
                      "p_rev": rv.get("p_day", np.nan),
                      "max_temp_c": s.get("max_temp_c"), "seconds": s.get("seconds")})
            fdr = {f["threshold"]: f for f in sk.get("fdr", [])}
            for tau in (0.01, 0.05):
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


def storey_pi0(p: np.ndarray, lam: float = LAMBDA) -> float:
    p = np.asarray(p, float)
    return float(min(1.0, np.mean(p > lam) / (1 - lam))) if len(p) else 1.0


def q_values(p: np.ndarray, pi0: float) -> np.ndarray:
    """Storey q-values: pi0 x the BH-adjusted p-values."""
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n)
    adj[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(pi0 * adj, 1.0)


def event_share(p: np.ndarray, n_boot: int = 2000, seed: int = 0) -> dict:
    """Estimated share and number of event days (1 - pi0), with a bootstrap 90% interval."""
    p = np.asarray(p, float)
    rng = np.random.default_rng(seed)
    pi0 = storey_pi0(p)
    boot = np.array([1 - storey_pi0(rng.choice(p, len(p))) for _ in range(n_boot)]) if len(p) else np.zeros(1)
    lo, hi = np.quantile(boot, [0.05, 0.95])
    out = {"n_days": int(len(p)), "pi0": pi0, "share": 1 - pi0, "share_lo": float(lo), "share_hi": float(hi),
           "n_events_storey": (1 - pi0) * len(p), "n_events_storey_lo": float(lo) * len(p),
           "n_events_storey_hi": float(hi) * len(p)}
    # for rare, strong events the excess at small p is far sharper than pi0: the days at or
    # below t, minus the t * N the null predicts (Poisson 90% interval on the null count)
    for t, key in ((0.01, "p01"), (0.05, "p05")):
        obs, exp = int(np.sum(p <= t)), t * len(p)
        out[f"observed_{key}"] = obs
        out[f"expected_{key}_null"] = exp
        out[f"excess_{key}"] = obs - exp
        out[f"excess_{key}_lo"] = obs - stats.poisson.ppf(0.95, exp) if exp > 0 else float(obs)
        out[f"excess_{key}_hi"] = obs - stats.poisson.ppf(0.05, exp) if exp > 0 else float(obs)
    return out


def fit_null(z: np.ndarray) -> dict:
    """Normal null fitted to the bulk of the z values: centre = median, width = IQR / 1.349."""
    z = z[np.isfinite(z)]
    if len(z) < 5:
        return {"centre": np.nan, "width": np.nan, "n": int(len(z))}
    q25, q50, q75 = np.quantile(z, [0.25, 0.5, 0.75])
    return {"centre": float(q50), "width": float((q75 - q25) / 1.349), "n": int(len(z))}


def tail_p(z: np.ndarray, null: dict) -> np.ndarray:
    """Upper-tail p-value under a fitted null; days without a score get 1."""
    z = np.asarray(z, float)
    return np.where(np.isfinite(z), stats.norm.sf((z - null["centre"]) / null["width"]), 1.0)


def probit(p: np.ndarray) -> np.ndarray:
    """z = Phi^-1(1 - p) for days with any evidence (p < 0.99), NaN otherwise."""
    p = np.asarray(p, float)
    return np.where(p < 0.99, stats.norm.isf(np.clip(p, 1e-12, 1.0)), np.nan)


def tested_days(rd: pd.DataFrame) -> pd.DataFrame:
    t = rd[rd.screened_out.fillna("") == ""].copy()
    t["p_day"] = t.p_day.fillna(1.0)
    return t


def summarise(t: pd.DataFrame, q: float) -> dict:
    p = t.p_day.to_numpy(float)
    return {"event_share": event_share(p), "bh_selected": int(benjamini_hochberg(p, q).sum()),
            "empirical_null_probit": fit_null(probit(p))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default=str(ROOT / "data" / "events"))
    ap.add_argument("--null-events", default=str(ROOT / "data" / "events_null"))
    ap.add_argument("--out", default=str(ROOT / "data" / "census"))
    ap.add_argument("--q", type=float, default=0.05)
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--top", type=int, default=30)
    a = ap.parse_args()
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rd = radar_days(pathlib.Path(a.events))
    tested = tested_days(rd)
    if tested.empty:
        rd.to_csv(out / "radar_days.csv", index=False)
        print(f"radar-days: {len(rd)}, none tested")
        return
    summ = {"census": summarise(tested, a.q)}
    pi0 = summ["census"]["event_share"]["pi0"]
    tested["q_value"] = q_values(tested.p_day.to_numpy(float), pi0)
    tested["event"] = benjamini_hochberg(tested.p_day.to_numpy(float), a.q)
    cutoff = float(tested.loc[tested.event, "p_day"].max()) if tested.event.any() else 0.0
    summ["census"]["p_cutoff"] = cutoff
    nroot = pathlib.Path(a.null_events)
    if nroot.exists() and any(nroot.glob("*/*.json")):
        nd = tested_days(radar_days(nroot))
        if len(nd):
            summ["null"] = summarise(nd, a.q)
            summ["null"]["selected_at_census_cutoff"] = int((nd.p_day <= cutoff).sum()) if cutoff else 0
            nd.to_csv(out / "null_days.csv", index=False)
    (out / "summary.json").write_text(json.dumps(summ, indent=1, default=float))
    keep = ["site", "day", "p_day", "q_value", "event"]
    rd = rd.merge(tested[keep], on=["site", "day"], how="left")
    rd["event"] = rd.event.fillna(False).astype(bool)
    rd.to_csv(out / "radar_days.csv", index=False)
    ev = tested[tested.event]
    ev.to_csv(out / "event_days.csv", index=False)
    cand = tested.sort_values("p_day").head(a.top)
    cand.to_csv(out / "candidates.csv", index=False)

    es = summ["census"]["event_share"]
    print(f"radar-days: {len(rd)}, screened out {int((rd.screened_out.fillna('') != '').sum())}, "
          f"tested {len(tested)}, confirmed with 200 rotations {int(tested.confirmed.sum())}")
    print(f"day p-values <= 0.01: {es['observed_p01']} (null: {es['expected_p01_null']:.1f}); "
          f"<= 0.05: {es['observed_p05']} (null: {es['expected_p05_null']:.1f})")
    print(f"excess of days at p <= 0.01: {es['excess_p01']:.1f} (90%: {es['excess_p01_lo']:.1f} to "
          f"{es['excess_p01_hi']:.1f}); at p <= 0.05: {es['excess_p05']:.1f}; Storey estimate of all event days "
          f"{es['n_events_storey']:.1f} ({es['n_events_storey_lo']:.1f}-{es['n_events_storey_hi']:.1f})")
    print(f"BH q={a.q}: {len(ev)} days selected (p cutoff {cutoff:.2g})")
    en = summ["census"]["empirical_null_probit"]
    print(f"empirical null of probit p (days with evidence, n={en['n']}): centre {en['centre']:.2f}, "
          f"width {en['width']:.2f}")
    if "null" in summ:
        ns = summ["null"]["event_share"]
        print(f"warm null days: {ns['n_days']}; p <= 0.01: {ns['observed_p01']} (null {ns['expected_p01_null']:.1f}), "
              f"<= 0.05: {ns['observed_p05']} (null {ns['expected_p05_null']:.1f}); excess at 0.01 "
              f"{ns['excess_p01']:.1f}; BH selected "
              f"{summ['null']['bh_selected']}; at the census cutoff {summ['null']['selected_at_census_cutoff']}")
    print("top candidates:")
    print(cand[["site", "day", "n_streaks", "n_flights_matched", "E", "p_day", "q_value", "p_rev"]]
          .head(10).to_string(index=False))
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
        print(f"catalogue: {len(cat)} flights with a streak at p <= {a.alpha} on selected days")


if __name__ == "__main__":
    main()
