"""Objective inspection of candidate days: does each matched flight look like an aircraft
streak, or like a piece of natural precipitation that a flight path happened to follow?

usage: inspect_candidates.py [--census data/census] [--events data/events] [--top 30]
                             [--alpha 0.05] [--days SITE:YYYY-MM-DD,...]

For every flight matched at p <= alpha on a candidate day (census_summary.py's candidates,
or --days), at its best-matched scan:
  isolation    share of the ring 2-8 km from the streak's centre line with echo >= 5 dBZ
               (natural bands and widespread snow fill it; an aircraft streak in clear air
               does not)
  n_scans      scans in which the same flight matched a streak at p <= alpha (an aircraft
               streak persists and grows for half an hour or more)
  parallel     other streak objects in the same scan within 15 deg of its orientation and
               60 km of it: natural snow bands come in families of parallel lines (lake-effect
               and banded snow), an aircraft streak usually does not
  v_fall, seed altitude   from the match
and a rule-based class, fixed in advance:
  convincing   isolation <= 0.30, n_scans >= 2, 0.3 <= v_fall <= 1.5 m/s, seed below 16 000 ft
  doubtful     isolation > 0.60, or parallel >= 3, or v_fall outside 0.25-1.6 m/s, or seed at
               or above 16 000 ft
  plausible    otherwise
The day's class is that of its best flight. Writes CENSUS/inspection_flights.csv and
CENSUS/inspection_days.csv.
"""
import argparse
import datetime as dt
import pathlib

import numpy as np
import pandas as pd
from scipy import ndimage

from snow import match, pipeline, streaks

ROOT = pathlib.Path(__file__).resolve().parents[1]
RING_M = (2000.0, 8000.0)
ECHO_DBZ = 5.0


PARALLEL_DEG = 15.0
PARALLEL_KM = 60.0


def orientation(lab: np.ndarray, label: int) -> float:
    """Orientation (deg, 0-180) of an object's pixels from their principal axis."""
    yy, xx = np.nonzero(lab == label)
    if len(xx) < 3:
        return np.nan
    c = np.cov(np.vstack([xx - xx.mean(), yy - yy.mean()]))
    w, v = np.linalg.eigh(c)
    vx, vy = v[:, np.argmax(w)]
    return float(np.degrees(np.arctan2(vy, vx)) % 180.0)


def parallel_neighbours(lab: np.ndarray, objs: pd.DataFrame, label: int, res_m: float) -> int:
    """Other objects in the scan within PARALLEL_DEG of this one's orientation, PARALLEL_KM of it."""
    th0 = orientation(lab, label)
    me = objs[objs.label == label].iloc[0]
    n = 0
    for _, o in objs.iterrows():
        if o.label == label or np.hypot(o.x - me.x, o.y - me.y) > PARALLEL_KM * 1e3:
            continue
        th = orientation(lab, int(o.label))
        d = abs((th - th0 + 90.0) % 180.0 - 90.0)
        if np.isfinite(d) and d <= PARALLEL_DEG:
            n += 1
    return n


def classify(r) -> str:
    if (r.isolation > 0.60 or r.parallel >= 3 or not (0.25 <= r.v_fall <= 1.6)
            or r.seed_alt_ft >= 16000):
        return "doubtful"
    if r.isolation <= 0.30 and r.n_scans >= 2 and 0.3 <= r.v_fall <= 1.5 and r.seed_alt_ft < 16000:
        return "convincing"
    return "plausible"


def inspect_day(site: str, day: dt.date, ev: pathlib.Path, alpha: float) -> list[dict]:
    path = ev / site / f"{day}_streaks.parquet"
    if not path.exists():
        return []
    S = pd.read_parquet(path)
    M = S[(S.coverage > 0) & (S.p_value <= alpha)]
    if M.empty:
        return []
    prep = pipeline.prepare(site, day)
    G = prep["grids"]
    kept, _, xs, ys = pipeline.streak_objects_day(prep)
    labs = {float(G.t[k]): (k, lab, obj) for k, obj, lab in kept}
    res = G.res
    rows = []
    for f, g in M.groupby("flight"):
        o = g.sort_values(["p_value", "coverage"], ascending=[True, False]).iloc[0]
        k, lab, objs = labs[min(labs, key=lambda T: abs(T - o.t))]
        obj = lab == int(o.label)
        dist = ndimage.distance_transform_edt(~obj) * res            # m from the object
        ring = (dist >= RING_M[0]) & (dist <= RING_M[1])
        z = G.z[k]
        valid = ring & np.isfinite(z)
        iso = float(np.mean(z[valid] >= ECHO_DBZ)) if valid.any() else np.nan
        rows.append({"site": site, "day": str(day), "flight": int(f), "t": o.t, "p_value": o.p_value,
                     "f1": o.coverage, "isolation": iso, "n_scans": int(g.t.nunique()),
                     "parallel": parallel_neighbours(lab, objs, int(o.label), res),
                     "v_fall": float(o.v_fall) if "v_fall" in o and np.isfinite(o.v_fall) else np.nan,
                     "seed_alt_ft": float(o.seed_alt_ft), "length_km": o.length_m / 1e3, "z_max": o.z_max})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", default=str(ROOT / "data" / "census"))
    ap.add_argument("--events", default=str(ROOT / "data" / "events"))
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--days", default="")
    a = ap.parse_args()
    cen, ev = pathlib.Path(a.census), pathlib.Path(a.events)
    if a.days:
        days = [tuple(x.split(":")) for x in a.days.split(",")]
    else:
        cand = pd.read_csv(cen / "candidates.csv").head(a.top)
        days = list(zip(cand.site, cand.day))
    rows = []
    for site, d in days:
        try:
            rows += inspect_day(site, dt.date.fromisoformat(d), ev, a.alpha)
        except Exception as e:                         # a day whose radar is gone, for one
            print(f"{site} {d}: {e!r}")
        print(f"{site} {d}: done", flush=True)
    F = pd.DataFrame(rows)
    if F.empty:
        print("nothing to inspect")
        return
    F["class"] = F.apply(classify, axis=1)
    rank = {"convincing": 0, "plausible": 1, "doubtful": 2}
    D = (F.assign(r=F["class"].map(rank)).sort_values("r").groupby(["site", "day"])
         .agg(day_class=("class", "first"), n_flights=("flight", "size"),
              n_convincing=("class", lambda c: int((c == "convincing").sum())),
              n_plausible=("class", lambda c: int((c == "plausible").sum())))
         .reset_index())
    cen.mkdir(parents=True, exist_ok=True)
    F.to_csv(cen / "inspection_flights.csv", index=False)
    D.to_csv(cen / "inspection_days.csv", index=False)
    pd.set_option("display.width", 200)
    print(F.round(3).to_string(index=False))
    print(D.to_string(index=False))


if __name__ == "__main__":
    main()
