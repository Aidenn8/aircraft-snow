"""Collect radar-day results into tables, estimate how many days had aircraft snow, and
select and rank candidate days.

usage: census_summary.py [--events data/events] [--out data/census]
                         [--null-events data/events_null] [--q 0.05] [--alpha 0.01] [--top 30]

The day's evidence. Each radar-day's flights are scored by E, the sum over matched flights
of -log10 of their streak p-values, and so is every copy of the day's flights rotated about
the radar, each scored against sets at least 20 deg away (snow/match.py). The day's score is
z = (E - mean) / s.d. of the sets it is compared with (s.d. floored), which keeps the margin
by which the real flights beat the rotations.

The day's p-value. A day has only some 15-25 independent rotations, so neither the rank of
E among them (which cannot say by how much a day beats them all) nor a normal tail for z
(far too optimistic on sparse days, where the rotations match almost nothing) is exact.
Instead every rotation's own z, computed the same way, is pooled over all tested days,
separately for 24- and 200-rotation runs and within strata of the day's chance-match level,
and a real day's p-value is the share of its pool at or above its z (self_calibrate). Under
no seeding the real set is one more rotation, so this is valid, and the pool is large.

Band-filtered census (--band-max N). Runways are chosen by the wind and snow bands line up
with it, so on a banded day pieces of bands match the real arrival paths more often than
rotated copies of them. Whether a streak object sits in a family of parallel ones (its
'parallel' count: other objects in the same scan within 15 deg and 60 km) depends on the
radar alone, so every day can be re-scored, real and rotated sets alike, without objects
having more than N parallel neighbours (from the saved match tables), and the whole
analysis below repeated on that score; write it to its own --out.

What is reported. (1) How many days had aircraft snow: the days at or below p = 0.01 (0.05)
beyond the 1% (5%) expected under no seeding, with an interval; Storey's estimate of all
event days as a noisier companion. (2) Which days: Benjamini-Hochberg at --q, Storey
q-values, and the --top candidates for inspection. (3) Checks: warm-season null days
(--null-events, tracks below 13 000 ft in July, where no aircraft can seed), calibrated the
same way against their own rotations, should show no excess; the false-alarm rate they show
at each level is also used to recompute the census excess ("calibrated"); the rank p-value
(p_rank) and a score counting places instead of flights (p_place) are kept as sensitivities.

Outputs (in --out):
  radar_days.csv    one row per radar-day: screens, scores, p-values, selections, q-values
  event_days.csv    radar-days selected by BH at --q
  candidates.csv    the --top days ranked by p-value
  null_days.csv     the warm null days, scored the same way
  summary.json      event-day estimates, selections, null-day check
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


def rescore(path: pathlib.Path, band_max: int) -> dict:
    """The day's streak statistic re-computed from its saved match table without objects
    having more than band_max parallel neighbours (same keys as the saved 'streaks')."""
    from snow import match, pipeline
    if not path.exists():
        return {"E": 0.0, "p_day": 1.0, "n_streaks": 0, "n_matched": 0}
    S = pd.read_parquet(path)
    if S.empty:
        return {"E": 0.0, "p_day": 1.0, "n_streaks": 0, "n_matched": 0}
    if "parallel" not in S:
        raise RuntimeError(f"{path}: no parallel counts (run scripts/backfill_parallel.py)")
    R = len(S.placebo_cov.iloc[0])
    angles = {len(pipeline.CONFIRM_ANGLES): pipeline.CONFIRM_ANGLES,
              len(pipeline.SCREEN_ANGLES): pipeline.SCREEN_ANGLES}[R]
    S = S[S.parallel <= band_max].reset_index(drop=True)
    if S.empty:
        return {"E": 0.0, "p_day": 1.0, "n_streaks": 0, "n_matched": 0, "n_rotations": R,
                "E_null": [0.0] * R}
    st = match.day_statistic(S, angles)
    st["n_rotations"] = R
    st["fdr"] = match.empirical_fdr(S, thresholds=(0.01, 0.05), angles=angles).to_dict(orient="records")
    return st


def radar_days(ev_root: pathlib.Path, band_max: int | None = None) -> pd.DataFrame:
    rows = []
    for path in sorted(ev_root.glob("*/*.json")):
        if "_" in path.stem:          # auxiliary results (_dayswap, _shiftnull), not radar-days
            continue
        s = json.loads(path.read_text())
        r = {"site": s["site"], "day": s["day"], "screened_out": s.get("screened_out", "")}
        if not r["screened_out"]:
            sc = s.get("streaks_screen") or s.get("streaks", {})
            sk = s.get("streaks", {})
            if band_max is not None:
                sk = rescore(path.parent / f"{path.stem}_streaks.parquet", band_max)
            rv = s.get("streaks_reverse", {})
            confirmed = sk.get("n_rotations", 0) > len(sc.get("E_null", [])) > 0
            p_rank = sk.get("p_day", 1.0)
            p_gamma = sk.get("p_day_gamma", np.nan)
            fin = sk if len(sk.get("E_null", [])) else sc
            r["_E_all"] = [float(fin.get("E", 0.0))] + [float(e) for e in fin.get("E_null", [])]
            r.update({"n_flights": s.get("n_flights"), "n_scans": s.get("n_scans"),
                      "n_streaks": sk.get("n_streaks", 0), "n_matched": sk.get("n_matched", 0),
                      "n_flights_matched": sk.get("n_flights_matched", 0),
                      "E_screen": sc.get("E", 0.0), "z_screen": sc.get("z_day", np.nan),
                      "p_screen": sc.get("p_day", 1.0), "confirmed": confirmed,
                      "E": sk.get("E", 0.0), "n_rotations": sk.get("n_rotations", 0),
                      "p_rank": p_rank,
                      "p_day": float(p_rank), "p_gamma": float(p_gamma),
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


SD_FLOOR = 0.7            # z uses max(s.d., this): about half one flight's evidence at p = 0.04
N_STRATA = 3              # days grouped by how many chance matches the rotations make


def set_z(E_all: np.ndarray) -> np.ndarray | None:
    """z of every set (column 0 real, 1..R rotations) against the sets it is comparable with
    (match.comparable_sets: at least 20 deg away), with the s.d. floored at SD_FLOOR."""
    from snow import match, pipeline
    R = len(E_all) - 1
    if R == len(pipeline.CONFIRM_ANGLES):
        angles = pipeline.CONFIRM_ANGLES
    elif R == len(pipeline.SCREEN_ANGLES):
        angles = pipeline.SCREEN_ANGLES
    else:
        return None
    comp = match.comparable_sets(angles)
    z = np.empty(R + 1)
    for j in range(R + 1):
        ref = E_all[comp[j]]
        z[j] = (E_all[j] - ref.mean()) / max(ref.std(ddof=1), SD_FLOOR)
    return z


def self_calibrate(t: pd.DataFrame) -> pd.DataFrame:
    """Calibrated day p-values from the rotations themselves. Every rotation, scored as if it
    were the real set, has a z like the real one's; under no seeding the real z is just one
    more of them. The rotations' z values are pooled over all tested days, separately for
    24- and 200-rotation runs and within strata of the day's chance-match level (median
    rotation E), and each real day's p-value is the share of its pool at or above its z,
    (1 + count) / (1 + N). Unlike the rank p-value this keeps the margin by which a day
    beats its rotations, and unlike the nominal z it is calibrated for sparse days."""
    t = t.copy()
    zs, pools = [], {}
    dens = t["_E_all"].map(lambda e: float(np.median(e[1:])) if len(e) > 1 else 0.0)
    t["chance_level"] = dens
    edges = np.quantile(dens[dens > 0], np.linspace(0, 1, N_STRATA + 1)[1:-1]) if (dens > 0).sum() > N_STRATA else []
    t["stratum"] = np.where(dens <= 0, 0, 1 + np.searchsorted(edges, dens, side="right"))
    for i, r in t.iterrows():
        e = np.asarray(r["_E_all"], float)
        z = set_z(e) if len(e) > 1 else None
        if z is None:
            zs.append(np.nan)
            continue
        zs.append(z[0])
        pools.setdefault((len(e) - 1, r["stratum"]), []).append(z[1:])
    t["z_day"] = zs
    pools = {k: np.sort(np.concatenate(v)) for k, v in pools.items()}
    pc, npool = [], []
    for _, r in t.iterrows():
        pool = pools.get((len(r["_E_all"]) - 1, r["stratum"]))
        if pool is None or not np.isfinite(r["z_day"]):
            pc.append(1.0)
            npool.append(0)
            continue
        k = len(pool) - np.searchsorted(pool, r["z_day"], side="left")      # pool values >= z
        pc.append((1 + k) / (1 + len(pool)))
        npool.append(len(pool))
    t["p_cal"] = pc
    t["pool_n"] = npool
    return t


def place_p(path: pathlib.Path) -> float:
    """Sensitivity: the day's rank p-value with evidence counted once per place
    (match.place_evidence) instead of once per flight, from the saved streak table."""
    from snow import match, pipeline
    if not path.exists():
        return 1.0
    S = pd.read_parquet(path)
    if S.empty:
        return 1.0
    R = len(S.placebo_cov.iloc[0])
    angles = pipeline.CONFIRM_ANGLES if R == len(pipeline.CONFIRM_ANGLES) else pipeline.SCREEN_ANGLES
    return float(match.day_statistic(S, angles, by_place=True)["p_day"])


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
    ap.add_argument("--band-max", type=int, default=None,
                    help="re-score days without streak objects with more parallel neighbours than this")
    a = ap.parse_args()
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rd = radar_days(pathlib.Path(a.events), a.band_max)
    tested = tested_days(rd)
    if tested.empty:
        rd.drop(columns=["_E_all"], errors="ignore").to_csv(out / "radar_days.csv", index=False)
        print(f"radar-days: {len(rd)}, none tested")
        return
    tested = self_calibrate(tested)
    tested["p_rank"] = tested.p_day
    tested["p_day"] = tested.p_cal                 # the calibrated p-value is the day's p-value
    summ = {"census": summarise(tested, a.q)}
    pi0 = summ["census"]["event_share"]["pi0"]
    tested["q_value"] = q_values(tested.p_day.to_numpy(float), pi0)
    tested["event"] = benjamini_hochberg(tested.p_day.to_numpy(float), a.q)
    cutoff = float(tested.loc[tested.event, "p_day"].max()) if tested.event.any() else 0.0
    summ["census"]["p_cutoff"] = cutoff
    nroot = pathlib.Path(a.null_events)
    if nroot.exists() and any(nroot.glob("*/*.json")):
        nd = tested_days(radar_days(nroot, a.band_max))
        if len(nd):
            nd = self_calibrate(nd)
            nd["p_rank"] = nd.p_day
            nd["p_day"] = nd.p_cal
            summ["null"] = summarise(nd, a.q)
            summ["null"]["selected_at_census_cutoff"] = int((nd.p_day <= cutoff).sum()) if cutoff else 0
            nd.drop(columns=["_E_all"]).to_csv(out / "null_days.csv", index=False)
            # the rotation p-value is not exact on sparse days (few independent rotations), so
            # the false-alarm rate at each level is taken from the null days themselves, and
            # the census excess is measured against it (binomial uncertainty on both counts)
            cal = {}
            for t, key in ((0.01, "p01"), (0.05, "p05")):
                k0, n0 = int((nd.p_day <= t).sum()), len(nd)
                k1, n1 = int((tested.p_day <= t).sum()), len(tested)
                f0 = (k0 + 0.5) / (n0 + 1)                       # shrunk a little towards 0.5/n
                lo0, hi0 = stats.beta.ppf([0.05, 0.95], k0 + 0.5, n0 - k0 + 0.5)
                cal[key] = {"null_rate": f0, "null_rate_90": [float(lo0), float(hi0)], "null_k": k0, "null_n": n0,
                            "census_k": k1, "census_n": n1, "expected_false": f0 * n1,
                            "excess": k1 - f0 * n1, "excess_lo": k1 - hi0 * n1, "excess_hi": k1 - lo0 * n1}
            summ["calibrated"] = cal
    (out / "summary.json").write_text(json.dumps(summ, indent=1, default=float))
    tested["p_place"] = ([place_p(pathlib.Path(a.events) / r.site / f"{r.day}_streaks.parquet")
                          for r in tested.itertuples()] if a.band_max is None else np.nan)
    keep = ["site", "day", "p_day", "z_day", "p_rank", "p_place", "chance_level", "stratum", "q_value", "event"]
    rd = rd.drop(columns=["p_day", "_E_all"], errors="ignore").merge(tested[keep], on=["site", "day"], how="left")
    tested = tested.drop(columns=["_E_all"])
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
    if "calibrated" in summ:
        for key, c in summ["calibrated"].items():
            print(f"calibrated by null days at {key}: null rate {c['null_rate']:.3f} ({c['null_k']}/{c['null_n']}; 90% "
                  f"{c['null_rate_90'][0]:.3f}-{c['null_rate_90'][1]:.3f}); census {c['census_k']}/{c['census_n']}, "
                  f"expected false {c['expected_false']:.1f}, excess {c['excess']:.1f} ({c['excess_lo']:.1f} to {c['excess_hi']:.1f})")
    print("top candidates:")
    print(cand[["site", "day", "n_streaks", "n_flights_matched", "E", "z_day", "p_rank", "p_day", "q_value", "p_rev"]]
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
