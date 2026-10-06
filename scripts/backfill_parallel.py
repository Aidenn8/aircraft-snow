"""Add streak orientation and parallel-neighbour counts to match tables saved before the
census recorded them (pipeline.streak_objects_day now does).

usage: backfill_parallel.py [--events data/events] [--workers 2]

For every SITE/DAY whose _streaks.parquet (and _screen_streaks.parquet) has no 'parallel'
column, the day's radar is loaded again (fetched if its cache was dropped), the streak
objects are found as in the census, and their orientation and parallel count are joined to
the saved rows by scan time and label; a saved row with no matching object is an error
(the objects depend on the radar alone, so they must be the same). Radar files fetched
here are deleted again afterwards.
"""
import argparse
import datetime as dt
import multiprocessing as mp
import pathlib
import signal

import pandas as pd
import pyarrow.parquet as pq

ROOT = pathlib.Path(__file__).resolve().parents[1]


def needs(path: pathlib.Path) -> bool:
    return path.exists() and "parallel" not in pq.read_schema(path).names


def cached_scans(site: str, day: dt.date) -> int:
    from snow import pipeline, radar
    t0, t1 = pipeline.day_window(day)
    return sum(1 for d in {t0.date(), t1.date()}
               for _ in (radar.CACHE / "l3").glob(f"{site[1:]}_N0B_{d:%Y_%m_%d}_*"))


def backfill(job: tuple[str, str, str]) -> str:
    from snow import pipeline, streaks
    events, site, day_s = job
    day = dt.date.fromisoformat(day_s)
    ev = pathlib.Path(events) / site
    paths = [p for p in (ev / f"{day}_streaks.parquet", ev / f"{day}_screen_streaks.parquet") if needs(p)]
    if not paths:
        return f"{site} {day}: nothing to do"
    had_radar = cached_scans(site, day) >= 20
    try:
        prep = pipeline.prepare(site, day)
        kept, _, _, _ = pipeline.streak_objects_day(prep)   # now carries orientation, parallel
        par = pd.concat([o[["t", "label", "orientation", "parallel"]] for _, o, _ in kept], ignore_index=True)
        par = par.assign(_t=par.t.round(3)).drop(columns="t")
        for p in paths:
            S = pd.read_parquet(p)
            S = S.assign(_t=S.t.round(3)).merge(par, on=["_t", "label"], how="left").drop(columns="_t")
            miss = int(S.parallel.isna().sum())
            if miss:
                return f"{site} {day}: FAILED {miss} of {len(S)} rows of {p.name} have no object"
            S["parallel"] = S.parallel.astype(int)
            tmp = p.with_suffix(".tmp")
            S.to_parquet(tmp)
            tmp.replace(p)
        return f"{site} {day}: {len(paths)} table(s), {len(par)} objects"
    except Exception as e:                     # noqa: BLE001 - report and continue
        return f"{site} {day}: FAILED {e!r}"
    finally:
        if not had_radar:
            pipeline.drop_radar_cache(site, day)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default=str(ROOT / "data" / "events"))
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    ev = pathlib.Path(a.events)
    jobs = []
    for js in sorted(ev.glob("K???/????-??-??.json")):
        site, day = js.parent.name, js.stem
        if any(needs(js.parent / f"{day}{s}") for s in ("_streaks.parquet", "_screen_streaks.parquet")):
            jobs.append((str(ev), site, day))
    # day by day across radars, never adjacent days of one radar at once: a day's cache
    # clean-up deletes the radar files of the hour it shares with its neighbours
    jobs.sort(key=lambda j: (j[2], j[1]))
    print(f"{len(jobs)} radar-days to backfill in {ev}", flush=True)
    with mp.get_context("spawn").Pool(a.workers) as pool:
        def stop(*_):
            pool.terminate()
            raise SystemExit(1)
        signal.signal(signal.SIGTERM, stop)
        for msg in pool.imap_unordered(backfill, jobs):
            print(msg, flush=True)


if __name__ == "__main__":
    main()
