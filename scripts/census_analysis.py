"""Census results: how often, where, which aircraft, at what temperatures.

usage: census_analysis.py [--census data/census] [--events data/events] [--alpha 0.01]

Reads census_summary.py's tables and, for each event day, the per-streak descriptions
(describe_day.py, *_described.csv) and per-flight exposure tables (seeding_efficiency.py,
*_exposure.csv). Writes to --census:
  rates.csv          per radar and winter: radar-days, tested, event days, share
  seeders.csv        one row per seeding flight (best streak), with type, engine class,
                     phase, seeding temperature and altitude, implied fall speed, snowfall
  by_type.csv        seeding flights by ICAO type and engine class
  efficiency.csv     first-passage flights on event days, pooled by engine class and
                     temperature bin, with the rotation chance level
and prints a summary.
"""
import argparse
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
TEMP_BINS = [-30, -20, -17, -14, -11, -8, -4]


def winter_of(day: pd.Series) -> pd.Series:
    """Winter label: November-March of Y/Y+1 is 'Y-YY'; other months are labelled by year."""
    d = pd.to_datetime(day)
    y = np.where(d.dt.month >= 7, d.dt.year, d.dt.year - 1)
    return pd.Series([f"{a}-{str(a + 1)[2:]}" for a in y], index=day.index)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", default=str(ROOT / "data" / "census"))
    ap.add_argument("--events", default=str(ROOT / "data" / "events"))
    ap.add_argument("--alpha", type=float, default=0.01)
    a = ap.parse_args()
    cen, evr = pathlib.Path(a.census), pathlib.Path(a.events)
    rd = pd.read_csv(cen / "radar_days.csv")
    rd["winter"] = winter_of(rd.day)
    rd["tested"] = rd.screened_out.fillna("") == ""
    rates = rd.groupby(["site", "winter"]).agg(radar_days=("day", "size"), tested=("tested", "sum"),
                                               event_days=("event", "sum")).reset_index()
    rates["event_share_of_days"] = rates.event_days / rates.radar_days
    rates.to_csv(cen / "rates.csv", index=False)
    print(rates.to_string(index=False))
    print(f"\nall: {len(rd)} radar-days, {int(rd.tested.sum())} tested, {int(rd.event.sum())} event days")

    ev = rd[rd.event]
    seeders, expo = [], []
    for _, d in ev.iterrows():
        p = evr / d.site / f"{d.day}_described.csv"
        if p.exists():
            s = pd.read_csv(p)
            s = s[s.p_value <= a.alpha].sort_values("p_value").drop_duplicates("flight")
            seeders.append(s.assign(site=d.site, day=d.day, winter=d.winter))
        else:
            print(f"  {d.site} {d.day}: no description yet")
        q = evr / d.site / f"{d.day}_exposure.csv"
        if q.exists():
            expo.append(pd.read_csv(q).assign(site=d.site, day=d.day))
    if seeders:
        S = pd.concat(seeders, ignore_index=True)
        S.to_csv(cen / "seeders.csv", index=False)
        bt = (S.groupby(["engine_class", "icao_type"]).agg(flights=("flight", "size"),
                                                            temp_median=("seed_temp_c", "median"),
                                                            temp_max=("seed_temp_c", "max"))
              .reset_index().sort_values("flights", ascending=False))
        bt.to_csv(cen / "by_type.csv", index=False)
        print(f"\n{len(S)} seeding flights on {S[['site', 'day']].drop_duplicates().shape[0]} event days")
        print(S.groupby("engine_class").seed_temp_c.describe().round(1))
        print(S.groupby("phase").size())
        if "v_fall" in S:
            print("implied fall speed (m/s):", S.v_fall.describe().round(2).to_dict())
        if "snow_flux_t_h" in S:
            print("snowfall under the streak (t/h water):", S.snow_flux_t_h.describe().round(0).to_dict())
    if expo:
        E = pd.concat(expo, ignore_index=True)
        E = E[E.first_passage.astype(bool)].copy()
        E["tbin"] = pd.cut(E.temp_c, TEMP_BINS)
        eff = (E.groupby(["engine_class", "tbin"], observed=True)
               .agg(n=("flight", "size"), made_streak=("made_streak", "mean"), chance=("placebo_rate", "mean"))
               .reset_index())
        eff["excess"] = (eff.made_streak - eff.chance) / (1 - eff.chance)
        eff.to_csv(cen / "efficiency.csv", index=False)
        print("\nfirst passages on event days:")
        print(eff.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
