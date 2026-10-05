"""Weather observations: radiosondes (NOAA IGRA2), airport METARs and pilot reports (Iowa
Environmental Mesonet), and the thermodynamic helpers used to find supercooled cloud.

IGRA2 data record columns (igra2-data-format.txt, v2.2):
  LVLTYP1 1, LVLTYP2 2, ETIME 4-8, PRESS 10-15 (Pa), PFLAG 16, GPH 17-21 (m), ZFLAG 22,
  TEMP 23-27 (0.1 C), TFLAG 28, RH 29-33 (0.1 %), DPDP 35-39 (0.1 C), WDIR 41-45 (deg),
  WSPD 47-51 (0.1 m/s); -8888 and -9999 mark removed or missing values.
"""
from __future__ import annotations

import datetime as dt
import io
import pathlib
import zipfile

import numpy as np
import pandas as pd
import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "raw" / "weather"
IGRA = "https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive"
IEM = "https://mesonet.agron.iastate.edu"

_session = requests.Session()


# ----------------------------------------------------------------------------- thermodynamics

def esat_water(t_c):
    """Saturation vapour pressure over liquid water, Pa (Murphy & Koop 2005, eq. 10)."""
    t = np.asarray(t_c, float) + 273.15
    return np.exp(54.842763 - 6763.22 / t - 4.210 * np.log(t) + 0.000367 * t
                  + np.tanh(0.0415 * (t - 218.8))
                  * (53.878 - 1331.22 / t - 9.44523 * np.log(t) + 0.014025 * t))


def esat_ice(t_c):
    """Saturation vapour pressure over ice, Pa (Murphy & Koop 2005, eq. 7), T < 273.16 K."""
    t = np.asarray(t_c, float) + 273.15
    return np.exp(9.550426 - 5723.265 / t + 3.53068 * np.log(t) - 0.00728332 * t)


def rh_ice_from_rh_water(rh_w, t_c):
    return np.asarray(rh_w, float) * esat_water(t_c) / esat_ice(t_c)


def pressure_from_pressure_altitude(h_ft):
    """ISA pressure (hPa) for a pressure altitude in feet (valid below 36 089 ft)."""
    return 1013.25 * (1.0 - 6.8755856e-6 * np.asarray(h_ft, float)) ** 5.2558797


def pressure_altitude_from_pressure(p_hpa):
    return (1.0 - (np.asarray(p_hpa, float) / 1013.25) ** (1 / 5.2558797)) / 6.8755856e-6


# ----------------------------------------------------------------------------- radiosondes

def igra_stations() -> pd.DataFrame:
    path = CACHE / "igra2-station-list.txt"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        r = _session.get(f"{IGRA}/doc/igra2-station-list.txt", timeout=120)
        r.raise_for_status()
        path.write_bytes(r.content)
    # columns: ID 1-11, LAT 13-20, LON 22-30, ELEV 32-37, STATE 39-40, NAME 42-71,
    # FSTYEAR 73-76, LSTYEAR 78-81, NOBS 83-88 (igra2-list-format.txt)
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        try:
            rows.append({"id": line[0:11].strip(), "lat": float(line[12:20]),
                         "lon": float(line[21:30]), "elev": float(line[31:37]),
                         "state": line[38:40].strip(), "name": line[41:71].strip(),
                         "first_year": int(line[72:76]), "last_year": int(line[77:81])})
        except ValueError:
            continue  # a few entries lack coordinates
    return pd.DataFrame(rows)


def _igra_zip(station: str) -> pathlib.Path:
    path = CACHE / "igra" / f"{station}-data.txt.zip"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        r = _session.get(f"{IGRA}/access/data-por/{station}-data.txt.zip", timeout=600)
        r.raise_for_status()
        path.write_bytes(r.content)
    return path


def _num(s: str, scale=1.0):
    s = s.strip()
    if not s:
        return np.nan
    v = int(s)
    return np.nan if v in (-8888, -9999) else v * scale


def soundings(station: str, start: dt.datetime, end: dt.datetime) -> pd.DataFrame:
    """All levels of the soundings launched between start and end (UTC, nominal time)."""
    rows = []
    with zipfile.ZipFile(_igra_zip(station)) as z:
        name = z.namelist()[0]
        with z.open(name) as f:
            keep = False
            for raw in io.TextIOWrapper(f, encoding="ascii", errors="replace"):
                if raw.startswith("#"):
                    y, mo, d, hh = int(raw[13:17]), int(raw[18:20]), int(raw[21:23]), int(raw[24:26])
                    if y < start.year - 1:
                        keep = False
                        continue
                    hh = 0 if hh == 99 else hh
                    t = dt.datetime(y, mo, d, hh, tzinfo=dt.timezone.utc)
                    keep = start <= t <= end
                    if t > end + dt.timedelta(days=1):
                        break
                    rel = raw[27:31]
                    continue
                if not keep:
                    continue
                rows.append({
                    "time": t, "lvltyp": raw[0:2], "etime": raw[3:8].strip(),
                    "p": _num(raw[9:15], 0.01), "z": _num(raw[16:21]), "t": _num(raw[22:27], 0.1),
                    "rh": _num(raw[28:33], 0.1), "dpd": _num(raw[34:39], 0.1),
                    "wdir": _num(raw[40:45]), "wspd": _num(raw[46:51], 0.1),
                })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["td"] = df["t"] - df["dpd"]
    # relative humidity over water from the dewpoint where RH itself is not reported
    rh_td = 100.0 * esat_water(df["td"]) / esat_water(df["t"])
    df["rh"] = df["rh"].fillna(rh_td)
    df["rhi"] = rh_ice_from_rh_water(df["rh"], df["t"])
    wd = np.radians(df["wdir"])
    df["u"] = -df["wspd"] * np.sin(wd)
    df["v"] = -df["wspd"] * np.cos(wd)
    return df


# ----------------------------------------------------------------------------- METAR, PIREP

def metars(station: str, start: dt.datetime, end: dt.datetime) -> pd.DataFrame:
    """Routine and special METARs from the IEM ASOS archive (station without the K)."""
    key = f"{station}_{start:%Y%m%d%H}_{end:%Y%m%d%H}.csv"
    path = CACHE / "metar" / key
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        params = {"station": station, "data": "all", "sts": f"{start:%Y-%m-%dT%H:%MZ}",
                  "ets": f"{end:%Y-%m-%dT%H:%MZ}", "tz": "Etc/UTC", "format": "onlycomma",
                  "latlon": "no", "elev": "no", "missing": "M", "trace": "T", "direct": "no",
                  "report_type": ["3", "4"]}
        r = _session.get(f"{IEM}/cgi-bin/request/asos.py", params=params, timeout=300)
        r.raise_for_status()
        path.write_bytes(r.content)
    df = pd.read_csv(path, na_values=["M"], keep_default_na=True)
    if not df.empty:
        df["valid"] = pd.to_datetime(df["valid"], utc=True)
    return df


def pireps(start: dt.datetime, end: dt.datetime) -> pd.DataFrame:
    """All US pilot reports between start and end (IEM archive)."""
    key = f"pireps_{start:%Y%m%d%H}_{end:%Y%m%d%H}.csv"
    path = CACHE / "pirep" / key
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        params = {"sts": f"{start:%Y-%m-%dT%H:%MZ}", "ets": f"{end:%Y-%m-%dT%H:%MZ}",
                  "fmt": "csv", "artcc": "_ALL"}
        r = _session.get(f"{IEM}/cgi-bin/request/gis/pireps.py", params=params, timeout=300)
        r.raise_for_status()
        path.write_bytes(r.content)
    df = pd.read_csv(path)
    return df
