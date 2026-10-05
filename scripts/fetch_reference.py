"""Download the reference tables used to describe aircraft (tar1090-db, by wiedehopf).

aircraft.csv.gz      hex;registration;ICAO type;flags;description;year;owner
icao_aircraft_types  ICAO type designator -> [description, wake category] e.g. "L2J", "M"
"""
import pathlib

import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "raw" / "aircraft_db"
BASE = "https://raw.githubusercontent.com/wiedehopf/tar1090-db"
FILES = {
    "aircraft.csv.gz": f"{BASE}/csv/aircraft.csv.gz",
    "icao_aircraft_types.js": f"{BASE}/master/db/icao_aircraft_types.js",
    "icao_aircraft_types2.js": f"{BASE}/master/db/icao_aircraft_types2.js",
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        (OUT / name).write_bytes(r.content)
        print(f"{name}: {len(r.content):,} bytes")


if __name__ == "__main__":
    main()
