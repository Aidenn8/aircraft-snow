"""Coordinate helpers: a local azimuthal-equidistant plane around a point (metres east and
north), great-circle distances, and the radar gate geometry."""
from __future__ import annotations

import functools

import numpy as np
import pyproj

from .radar import Sweep, beam_height_ground_range

EARTH_RADIUS = 6_371_000.0


@functools.lru_cache(maxsize=64)
def _transformer(lat0: float, lon0: float) -> pyproj.Transformer:
    aeqd = pyproj.CRS.from_proj4(f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +datum=WGS84 +units=m")
    return pyproj.Transformer.from_crs("EPSG:4326", aeqd, always_xy=True)


def to_xy(lat, lon, lat0, lon0):
    """(lat, lon) -> (x east, y north) in metres in the plane centred on (lat0, lon0)."""
    x, y = _transformer(float(lat0), float(lon0)).transform(np.asarray(lon), np.asarray(lat))
    return np.asarray(x), np.asarray(y)


def to_latlon(x, y, lat0, lon0):
    lon, lat = _transformer(float(lat0), float(lon0)).transform(
        np.asarray(x), np.asarray(y), direction="INVERSE")
    return np.asarray(lat), np.asarray(lon)


def haversine(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * EARTH_RADIUS * np.arcsin(np.sqrt(a))


def gate_geometry(sweep: Sweep, moment: str = "REF"):
    """Per-gate (x, y, height ASL) in metres, x/y in the plane centred on the radar.
    Arrays are (n_rays, n_gates). The plane's distance from the centre is the ground
    range s along the earth, which aeqd preserves exactly."""
    r = sweep.ranges[moment]
    h, s = beam_height_ground_range(r[None, :], sweep.elevation, sweep.radar_height)
    az = np.radians(sweep.azimuth)[:, None]
    x = s * np.sin(az)
    y = s * np.cos(az)
    return x, y, np.broadcast_to(h, x.shape)
