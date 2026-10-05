"""Run the event test for one radar over a range of analysis days.

usage: run_days.py SITE START END [--force]
An analysis day D covers 12 UTC on D to 12 UTC on D+1, so ADS-B for D and D+1 is needed;
days without it are skipped and listed.
"""
import argparse
import datetime as dt
import json
import sys
import traceback

from snow import adsb, pipeline


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("start")
    ap.add_argument("end")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--null-alt", type=float, default=None,
                    help="null mode: only tracks below this altitude (ft), no screens")
    a = ap.parse_args()
    d = dt.date.fromisoformat(a.start)
    end = dt.date.fromisoformat(a.end)
    while d <= end:
        need = [d, d + dt.timedelta(days=1)]
        if not all(adsb.day_path(x).exists() for x in need):
            print(f"{a.site} {d}: skipped (no ADS-B)", flush=True)
            d += dt.timedelta(days=1)
            continue
        try:
            s = pipeline.run_day(a.site, d, force=a.force, null_alt_max_ft=a.null_alt)
            if "screened_out" in s:
                print(f"{a.site} {d}: screened out ({s['screened_out']})", flush=True)
                d += dt.timedelta(days=1)
                continue
            st = s["stats"]
            sk = s.get("streaks", {})
            print(f"{a.site} {d}: flights {s['n_flights']:5d} scans {s['n_scans']:3d}  "
                  f"streaks {sk.get('n_streaks', 0):4d} matched {sk.get('n_matched', 0):3d} "
                  f"E={sk.get('E', 0):6.1f} p_day={sk.get('p_day', 1):.3f}  "
                  + "  ".join(f"{k}: p={v['p_cluster']:.3f}" for k, v in st.items())
                  + f"  ({s.get('seconds', 0):.0f}s)", flush=True)
        except Exception as e:
            print(f"{a.site} {d}: FAILED {e!r}", flush=True)
            traceback.print_exc(file=sys.stderr)
        d += dt.timedelta(days=1)


if __name__ == "__main__":
    main()
