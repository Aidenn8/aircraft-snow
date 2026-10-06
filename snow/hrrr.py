"""Vertical profiles from the HRRR analysis (NOAA High-Resolution Rapid Refresh, 3 km),
read from the University of Utah's cloud-optimised Zarr copy on AWS (s3://hrrrzarr).

Layout: prs/YYYYMMDD/YYYYMMDD_HHz_anl.zarr/<level>/<VAR>/<level>/<VAR>/<j>.<i>, one
150 x 150 point tile per file (zarr v2, blosc-lz4, float32, fill -9999). The grid is the
HRRR Lambert conformal grid, 1799 x 1059 points at 3 km, on a sphere of radius 6 371 229 m
with standard parallels 38.5 N and central meridian 97.5 W. Grid coordinates (metres) are
read from the store itself and the projection is checked against the published corner.

Only the tile containing a point is downloaded, so one profile costs n_levels x n_vars
requests of ~10 kB each.
"""
from __future__ import annotations

import concurrent.futures as cf
import datetime as dt
import functools
import json
import pathlib

import numpy as np
import pandas as pd
import pyproj
import requests
from numcodecs import Blosc

ROOT = pathlib.Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "raw" / "hrrr"
BUCKET = "https://hrrrzarr.s3.amazonaws.com"
CHUNK = 150
NY, NX = 1059, 1799
LEVELS = tuple(range(1000, 275, -25))     # 1000 ... 300 hPa every 25 hPa
VARS = ("TMP", "RH", "CLMR", "CIMIXR", "SNMR", "UGRD", "VGRD", "VVEL", "HGT")

HRRR_CRS = pyproj.CRS.from_proj4(
    "+proj=lcc +lat_0=38.5 +lon_0=-97.5 +lat_1=38.5 +lat_2=38.5 +a=6371229 +b=6371229 +units=m")

_session = requests.Session()
_blosc = Blosc()


def _get(url: str) -> bytes | None:
    for attempt in range(4):
        try:
            r = _session.get(url, timeout=60)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.content
        except requests.RequestException:
            if attempt == 3:
                raise
    return None


@functools.lru_cache(maxsize=1)
def grid_xy() -> tuple[np.ndarray, np.ndarray]:
    """Projection x (1799) and y (1059) coordinates in metres, from the store."""
    path = CACHE / "grid_xy.npz"
    if path.exists():
        z = np.load(path)
        return z["x"], z["y"]
    base = f"{BUCKET}/prs/20250126/20250126_01z_anl.zarr/700mb/TMP"
    x = np.frombuffer(_blosc.decode(_get(f"{base}/projection_x_coordinate/0")), "<f8")
    y = np.frombuffer(_blosc.decode(_get(f"{base}/projection_y_coordinate/0")), "<f8")
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, x=x, y=y)
    return x, y


def grid_index(lat: float, lon: float) -> tuple[int, int]:
    """Nearest HRRR grid point (j row, i column) to a latitude/longitude."""
    t = pyproj.Transformer.from_crs("EPSG:4326", HRRR_CRS, always_xy=True)
    px, py = t.transform(lon, lat)
    x, y = grid_xy()
    i = int(np.argmin(np.abs(x - px)))
    j = int(np.argmin(np.abs(y - py)))
    return j, i


def grid_latlon(j: int, i: int) -> tuple[float, float]:
    x, y = grid_xy()
    t = pyproj.Transformer.from_crs(HRRR_CRS, "EPSG:4326", always_xy=True)
    lon, lat = t.transform(x[i], y[j])
    return float(lat), float(lon)


def _tile(run: dt.datetime, level: str, var: str, cj: int, ci: int) -> np.ndarray | None:
    url = (f"{BUCKET}/prs/{run:%Y%m%d}/{run:%Y%m%d_%H}z_anl.zarr/"
           f"{level}/{var}/{level}/{var}/{cj}.{ci}")
    raw = _get(url)
    if raw is None:
        return None
    return _decode_tile(_blosc.decode(raw))


def _decode_tile(buf: bytes) -> np.ndarray:
    """A decompressed chunk as float32 with fill values as NaN. The archive's data type
    changed: 2023 runs are float16 with fill -10000, 2025 runs float32 with fill -9999
    (each array's .zarray says which); the byte count tells them apart."""
    width = len(buf) // (CHUNK * CHUNK)
    if width not in (2, 4) or len(buf) != width * CHUNK * CHUNK:
        raise ValueError(f"HRRR chunk of {len(buf)} bytes is neither float16 nor float32 {CHUNK}x{CHUNK}")
    a = np.frombuffer(buf, "<f2" if width == 2 else "<f4").astype(np.float32).reshape(CHUNK, CHUNK)
    a[(a == -9999.0) | (a == -10000.0)] = np.nan
    return a


def profile(run: dt.datetime, lat: float, lon: float, box: int = 0,
            levels=LEVELS, variables=VARS, workers: int = 32) -> pd.DataFrame:
    """Analysis profile at the grid point nearest (lat, lon) for one hourly run.

    With box > 0 the values are averaged over the (2 box + 1)^2 neighbouring points that
    fall inside the same tile (a 3 km point is noisy for cloud water).
    Returns one row per pressure level; temperature (column temp) in C, mixing ratios in g/kg,
    winds in m/s, VVEL in Pa/s, HGT in m.
    """
    run = run.replace(minute=0, second=0, microsecond=0)
    cache = CACHE / "profiles" / f"{run:%Y%m%d_%H}_{lat:.4f}_{lon:.4f}_{box}.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    j, i = grid_index(lat, lon)
    cj, ci = j // CHUNK, i // CHUNK
    jj, ii = j - cj * CHUNK, i - ci * CHUNK
    j0, j1 = max(jj - box, 0), min(jj + box + 1, CHUNK)
    i0, i1 = max(ii - box, 0), min(ii + box + 1, CHUNK)
    jobs = [(f"{p}mb", v) for p in levels for v in variables]
    with cf.ThreadPoolExecutor(workers) as ex:
        tiles = list(ex.map(lambda lv: _tile(run, lv[0], lv[1], cj, ci), jobs))
    rows = {}
    for (lev, var), a in zip(jobs, tiles):
        p = int(lev[:-2])
        val = np.nan if a is None else float(np.nanmean(a[j0:j1, i0:i1]))
        rows.setdefault(p, {"p": p})[var] = val
    df = pd.DataFrame(sorted(rows.values(), key=lambda r: -r["p"]))
    if "TMP" in df:
        df["temp"] = df["TMP"] - 273.15
    for m in ("CLMR", "CIMIXR", "SNMR"):
        if m in df:
            df[m] = df[m] * 1000.0  # kg/kg -> g/kg
    df["run"] = pd.Timestamp(run)
    df["lat"], df["lon"] = lat, lon
    cache.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache, index=False)
    return df


def grid_projection_check() -> dict:
    """Compare the store's grid with the published HRRR first grid point (21.138123 N,
    122.719528 W) and grid spacing (3 km)."""
    x, y = grid_xy()
    lat, lon = grid_latlon(0, 0)
    return {"corner_lat": lat, "corner_lon": lon,
            "corner_error_m": float(np.hypot(*_err(lat, lon))),
            "dx": float(np.median(np.diff(x))), "dy": float(np.median(np.diff(y)))}


def _err(lat, lon):
    from .geo import haversine
    return (haversine(lat, lon, 21.138123, -122.719528), 0.0)
