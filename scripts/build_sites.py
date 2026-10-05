"""Pair large US airports with their nearest NEXRAD radar -> data/sites.csv.

Sources: NOAA NCEI nexrad-stations.txt (fixed width) and OurAirports airports.csv.
"""
import pathlib

import numpy as np
import pandas as pd

from snow import sites
from snow.geo import haversine

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "sites"


def nexrad() -> pd.DataFrame:
    return sites.nexrad().reset_index()


def airports() -> pd.DataFrame:
    a = pd.read_csv(RAW / "airports.csv", keep_default_na=False)
    a = a[(a.iso_country == "US") & (a.type == "large_airport") & (a.scheduled_service == "yes")]
    return a[["ident", "iata_code", "name", "municipality", "iso_region", "latitude_deg",
              "longitude_deg", "elevation_ft"]]


def main():
    nx = nexrad()
    ap = airports()
    d = haversine(ap.latitude_deg.to_numpy()[:, None], ap.longitude_deg.to_numpy()[:, None],
                  nx.lat.to_numpy()[None, :], nx.lon.to_numpy()[None, :]) / 1000
    j = d.argmin(axis=1)
    ap = ap.assign(radar=nx.icao.to_numpy()[j], radar_km=d[np.arange(len(ap)), j],
                   radar_lat=nx.lat.to_numpy()[j], radar_lon=nx.lon.to_numpy()[j],
                   radar_elev_ft=nx.elev.to_numpy()[j])
    ap = ap.sort_values("radar_km")
    out = ROOT / "data" / "sites.csv"
    ap.to_csv(out, index=False)
    pd.set_option("display.width", 200)
    print(ap[["ident", "iata_code", "municipality", "iso_region", "latitude_deg", "radar", "radar_km",
              "radar_elev_ft", "elevation_ft"]].round(1).to_string())


if __name__ == "__main__":
    main()
