"""Census results: how often, where, which aircraft, at what temperatures.

usage: census_analysis.py [--census data/census] [--band data/census_band] [--events data/events]
                          [--review results/census/visual_review.csv] [--alpha 0.01]

Three layers, kept apart because they rest on different things.

1. How often (no human judgement). For each radar and winter and overall: tested radar-days
   and the days whose calibrated p-value (census_summary.py) is at or below 0.01 and 0.05,
   against the share expected under no seeding. With --band, the same for the census
   re-scored without streaks in families of parallel bands (census_summary.py --band-max).
   The counts are also split by hours of supercooled cloud water over the radar
   (slw_covariate.py), which real seeding needs and artefacts do not.          -> rates.csv, slw.csv

2. Which flights (the catalogue). Days whose matches were judged by eye to be aircraft
   streaks (--review, visual == convincing). On them, every flight whose best streak match
   has p <= --alpha and is not 'doubtful' by the objective rules of inspect_candidates.py
   (isolation, parallel bands, fall speed, altitude). The rotations' expected number of
   false matches on those days is summed for comparison.                       -> seeders.csv

3. What they show. Seeding temperature and altitude, aircraft type and engine class, fall
   speed, streak size and snowfall (describe_day.py, *_described.csv), and on catalogue days
   the share of aircraft crossing the seeding layer that left a streak, by engine class and
   temperature, against the rotation chance level (seeding_efficiency.py, *_exposure.csv),
   and the traffic mix of exposed flights against that of seeders.  -> by_type.csv,
                                                                       efficiency.csv, mix.csv
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[1]
TEMP_BINS = [-30, -20, -17, -14, -11, -8, -4]
SLW_BINS = [-1, 0, 3, 8, 24]
SLW_LABELS = ["0 h", "1-3 h", "4-8 h", "9-24 h"]


def winter_of(day: pd.Series) -> pd.Series:
    """Winter label: November-March of Y/Y+1 is 'Y-YY'; other months are labelled by year."""
    d = pd.to_datetime(day)
    y = np.where(d.dt.month >= 7, d.dt.year, d.dt.year - 1)
    return pd.Series([f"{a}-{str(a + 1)[2:]}" for a in y], index=day.index)


def counts(g: pd.DataFrame) -> pd.Series:
    """Tested days, days at p <= 0.01 / 0.05, the numbers expected under no seeding, and the
    one-sided Poisson probability of at least the observed count at 0.01."""
    n = len(g)
    k1, k5 = int((g.p_day <= 0.01).sum()), int((g.p_day <= 0.05).sum())
    return pd.Series({"tested": n, "p01": k1, "expected_p01": 0.01 * n, "p05": k5, "expected_p05": 0.05 * n,
                      "excess_p01": k1 - 0.01 * n, "excess_p05": k5 - 0.05 * n,
                      "poisson_p01": float(stats.poisson.sf(k1 - 1, 0.01 * n)) if n else np.nan})


def tested(cen: pathlib.Path) -> pd.DataFrame:
    rd = pd.read_csv(cen / "radar_days.csv")
    t = rd[rd.screened_out.fillna("") == ""].copy()
    t["winter"] = winter_of(t.day)
    return t


def how_often(cen: pathlib.Path, band: pathlib.Path | None) -> pd.DataFrame:
    out = []
    for name, d in (("census", cen), ("band-filtered", band)):
        if d is None or not (d / "radar_days.csv").exists():
            continue
        t = tested(d)
        by = t.groupby(["site", "winter"]).apply(counts, include_groups=False).reset_index()
        site = t.groupby("site").apply(counts, include_groups=False).reset_index().assign(winter="all")
        win = t.groupby("winter").apply(counts, include_groups=False).reset_index().assign(site="all")
        tot = counts(t).to_frame().T.assign(site="all", winter="all")
        out.append(pd.concat([by, site, win, tot], ignore_index=True).assign(score=name))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def by_slw(cen: pathlib.Path) -> pd.DataFrame:
    path = cen / "slw_days.csv"
    if not path.exists():
        return pd.DataFrame()
    t = tested(cen).merge(pd.read_csv(path), on=["site", "day"], how="left")
    t["slw_class"] = pd.cut(t.slw_hours, SLW_BINS, labels=SLW_LABELS)
    s = t.groupby("slw_class", observed=True).apply(counts, include_groups=False).reset_index()
    # within each chance-level stratum too: does supercooled cloud matter at equal echo?
    st = t.groupby(["stratum", "slw_class"], observed=True).apply(counts, include_groups=False).reset_index()
    return pd.concat([s.assign(stratum="all"), st], ignore_index=True)


def catalogue(cen: pathlib.Path, evr: pathlib.Path, review: pathlib.Path, alpha: float):
    rv = pd.read_csv(review)
    days = rv[rv.visual == "convincing"][["site", "day"]].drop_duplicates()
    insp = pd.read_csv(cen / "inspection_flights.csv")
    rd = pd.read_csv(cen / "radar_days.csv")
    rows, missing = [], []
    for d in days.itertuples():
        p = evr / d.site / f"{d.day}_described.csv"
        if not p.exists():
            missing.append(f"{d.site} {d.day} (no description)")
            continue
        s = pd.read_csv(p)
        s = s[s.p_value <= alpha].sort_values(["p_value", "f1"], ascending=[True, False]).drop_duplicates("flight")
        ins = insp[(insp.site == d.site) & (insp.day == d.day)][["flight", "isolation", "n_scans", "parallel", "class"]]
        if ins.empty:
            missing.append(f"{d.site} {d.day} (not inspected)")
            continue
        s = s.merge(ins, on="flight", how="left")
        rows.append(s.assign(site=d.site, day=d.day))
    S = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if len(S):
        S = S[S["class"].fillna("doubtful") != "doubtful"].copy()
        S["winter"] = winter_of(S.day)
    dd = rd.merge(days, on=["site", "day"])
    false = {"days": len(days), "flights_p01_all": float(dd.get("n_flights_p0.01", pd.Series(dtype=float)).sum()),
             "expected_false_p01": float(dd.get("exp_false_p0.01", pd.Series(dtype=float)).sum())}
    return S, false, missing


def wilson(k: np.ndarray, n: np.ndarray, z: float = 1.645) -> tuple[np.ndarray, np.ndarray]:
    """Wilson score interval (90% by default) for k successes in n trials."""
    k, n = np.asarray(k, float), np.asarray(n, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        ph = k / n
        c = (ph + z * z / (2 * n)) / (1 + z * z / n)
        h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def rate_table(E: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    g = (E.groupby(keys, observed=True)
         .agg(n=("flight", "size"), k=("made_streak", "sum"), chance=("placebo_rate", "mean")).reset_index())
    g["made_streak"] = g.k / g.n
    g["lo90"], g["hi90"] = wilson(g.k, g.n)
    g["excess"] = (g.made_streak - g.chance) / (1 - g.chance)
    return g


def mantel_haenszel(df: pd.DataFrame, exposed: str, outcome: str, strata: list[str], z: float = 1.645) -> dict:
    """Mantel-Haenszel common odds ratio of outcome for exposed vs not, over strata, with the
    Robins-Breslow-Greenland interval (90% by default)."""
    R = S = PR = PSQR = QS = 0.0
    for _, g in df.groupby(strata, observed=True):
        e, y = g[exposed].astype(bool), g[outcome].astype(bool)
        a, b = float((e & y).sum()), float((e & ~y).sum())
        c, d = float((~e & y).sum()), float((~e & ~y).sum())
        n = a + b + c + d
        if n == 0 or (a + b) == 0 or (c + d) == 0:
            continue
        r, s_, P, Q = a * d / n, b * c / n, (a + d) / n, (b + c) / n
        R, S = R + r, S + s_
        PR, PSQR, QS = PR + P * r, PSQR + P * s_ + Q * r, QS + Q * s_
    if R == 0 or S == 0:
        return {"or": np.nan, "lo": np.nan, "hi": np.nan}
    var = PR / (2 * R * R) + PSQR / (2 * R * S) + QS / (2 * S * S)
    lor = np.log(R / S)
    return {"or": float(np.exp(lor)), "lo": float(np.exp(lor - z * np.sqrt(var))),
            "hi": float(np.exp(lor + z * np.sqrt(var)))}


def phase_effect(E: pd.DataFrame) -> pd.DataFrame:
    """Jet arrivals against departures. Arrivals spend much longer in a layer (shallow descent,
    level segments) and a longer seeded line is likelier to make a detectable streak, so the
    comparison is stratified by the number of track points in the layer (quintiles) and by
    temperature. What remains mixes configuration (gear and flaps out) with airspeed (slower
    flight seeds more air per km); it does not separate them."""
    J = E[(E.engine_class == "jet") & E.phase.isin(["arrival", "departure"])].copy()
    if J.empty:
        return pd.DataFrame()
    J["path_q"] = pd.qcut(J.n_points, 5, duplicates="drop").astype(str)
    J["arrival"] = J.phase == "arrival"
    rows = [{"comparison": "crude", **mantel_haenszel(J.assign(_s=0), "arrival", "made_streak", ["_s"])},
            {"comparison": "by path length", **mantel_haenszel(J, "arrival", "made_streak", ["path_q"])},
            {"comparison": "by path length and temperature",
             **mantel_haenszel(J, "arrival", "made_streak", ["path_q", "tbin"])},
            # flights on one day share its clouds: compare arrivals with departures of the same day
            {"comparison": "by path length, temperature and day",
             **mantel_haenszel(J, "arrival", "made_streak", ["path_q", "tbin", "site", "day"])}]
    return pd.DataFrame(rows)


def efficiency(evr: pathlib.Path, days: pd.DataFrame):
    expo = []
    for d in days.itertuples():
        q = evr / d.site / f"{d.day}_exposure.csv"
        if q.exists():
            expo.append(pd.read_csv(q).assign(site=d.site, day=d.day))
    if not expo:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    E = pd.concat(expo, ignore_index=True)
    E = E[E.first_passage.astype(bool)].copy()
    E["tbin"] = pd.cut(E.temp_c, TEMP_BINS).astype(str)
    E["made_streak"] = E.made_streak.astype(bool)
    # by engine class and temperature, and (jets, which have the numbers) by flight phase:
    # gear and flaps out at low speed on approach cool the air over the wing most
    eff = pd.concat([rate_table(E, ["engine_class", "tbin"]).assign(phase="all"),
                     rate_table(E, ["engine_class"]).assign(tbin="all", phase="all"),
                     rate_table(E[E.engine_class == "jet"], ["engine_class", "phase", "tbin"]),
                     rate_table(E[E.engine_class == "jet"], ["engine_class", "phase"]).assign(tbin="all")],
                    ignore_index=True)
    # the traffic that crossed the layer against the traffic that seeded, by engine class
    mix = (E.groupby("engine_class").agg(exposed=("flight", "size"), seeded=("made_streak", "sum")).reset_index())
    mix["share_exposed"] = mix.exposed / mix.exposed.sum()
    mix["share_seeded"] = mix.seeded / max(mix.seeded.sum(), 1)
    return eff, mix, phase_effect(E)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", default=str(ROOT / "data" / "census"))
    ap.add_argument("--band", default=None, help="the band-filtered census_summary output")
    ap.add_argument("--events", default=str(ROOT / "data" / "events"))
    ap.add_argument("--review", default=str(ROOT / "results" / "census" / "visual_review.csv"))
    ap.add_argument("--alpha", type=float, default=0.01)
    a = ap.parse_args()
    cen, evr = pathlib.Path(a.census), pathlib.Path(a.events)
    band = pathlib.Path(a.band) if a.band else None

    rates = how_often(cen, band)
    rates.to_csv(cen / "rates.csv", index=False)
    tot = rates[(rates.site == "all") & (rates.winter == "all")]
    print("1. how often (all radars and winters):")
    print(tot.round(2).to_string(index=False))
    print(rates[(rates.site != "all") & (rates.winter == "all")].round(2).to_string(index=False))

    slw = by_slw(cen)
    if len(slw):
        slw.to_csv(cen / "slw.csv", index=False)
        print("\nby hours of supercooled cloud over the radar:")
        print(slw[slw.stratum == "all"].round(3).to_string(index=False))

    S, false, missing = catalogue(cen, evr, pathlib.Path(a.review), a.alpha)
    for m in missing:
        print(f"  catalogue day skipped: {m}")
    print(f"\n2. catalogue: {false['days']} days judged convincing; on them {false['flights_p01_all']:.0f} flights "
          f"matched at p <= 0.01 against {false['expected_false_p01']:.1f} expected by chance")
    if len(S):
        S.to_csv(cen / "seeders.csv", index=False)
        print(f"   {len(S)} catalogue flights (not doubtful) on {S[['site', 'day']].drop_duplicates().shape[0]} days, "
              f"{S.site.nunique()} radars")
        bt = (S.groupby(["engine_class", "icao_type"]).agg(flights=("flight", "size"),
                                                            temp_median=("seed_temp_c", "median"),
                                                            alt_median_ft=("seed_alt_ft", "median"))
              .reset_index().sort_values("flights", ascending=False))
        bt.to_csv(cen / "by_type.csv", index=False)
        print("\n3. seeding temperature (C) by engine class:")
        print(S.groupby("engine_class").seed_temp_c.describe().round(1).to_string())
        print("seeding altitude (ft):", S.seed_alt_ft.describe().round(0).to_dict())
        print("phase:", S.groupby("phase").size().to_dict())
        print("fall speed (m/s):", S.v_fall.describe().round(2).to_dict())
        print("age at match (min):", S.age_min.describe().round(0).to_dict())
        print("length (km):", S.length_km.describe().round(1).to_dict(), " width (km):",
              S.width_km.describe().round(1).to_dict())
        print("peak snowfall (mm/h):", S.snow_mm_h_max.describe().round(2).to_dict())
        print("snowfall under the streak (t/h water):", S.snow_flux_t_h.describe().round(0).to_dict())
    rv = pd.read_csv(a.review)
    eff, mix, ph = efficiency(evr, rv[rv.visual == "convincing"][["site", "day"]].drop_duplicates())
    if len(eff):
        eff.to_csv(cen / "efficiency.csv", index=False)
        mix.to_csv(cen / "mix.csv", index=False)
        if len(ph):
            ph.to_csv(cen / "phase.csv", index=False)
            print("\njet arrivals vs departures, odds ratio of leaving a streak (90%):")
            print(ph.round(2).to_string(index=False))
        print("\nfirst passages through the seeding layer on catalogue days:")
        print(eff.round(3).to_string(index=False))
        print(mix.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
