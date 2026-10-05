"""NEXRAD (WSR-88D) radar data from the public AWS archives, read with MetPy.

Level II  s3://unidata-nexrad-level2/YYYY/MM/DD/KSSS/KSSSYYYYMMDD_HHMMSS_V06
          full volume scans: every elevation, reflectivity plus Doppler and dual-pol moments.
Level III s3://unidata-nexrad-level3/SSS_PPP_YYYY_MM_DD_HH_MM_SS
          single-sweep products, e.g. N0B (0.5 deg super-resolution base reflectivity).

Beam geometry uses the standard 4/3 effective earth radius model (Doviak & Zrnic 1993):
  h = sqrt(r^2 + (k a)^2 + 2 r k a sin(theta)) - k a + h_radar
  s = k a asin(r cos(theta) / (k a + h - h_radar))
with r the slant range, theta the elevation angle, a the earth radius and k = 4/3.
"""
from __future__ import annotations

import dataclasses
import datetime as dt
import pathlib
import re
import xml.etree.ElementTree as ET

import numpy as np
import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "radar"
L2_BUCKET = "https://unidata-nexrad-level2.s3.amazonaws.com"
L3_BUCKET = "https://unidata-nexrad-level3.s3.amazonaws.com"
EARTH_RADIUS = 6_371_000.0
K_EFF = 4.0 / 3.0

_session = requests.Session()
_NS = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}


def _list(bucket: str, prefix: str) -> list[tuple[str, int]]:
    """All (key, size) under a prefix, following S3 continuation tokens."""
    out, token = [], None
    while True:
        params = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if token:
            params["continuation-token"] = token
        r = _session.get(bucket + "/", params=params, timeout=60)
        r.raise_for_status()
        root = ET.fromstring(r.content)
        for c in root.findall("s3:Contents", _NS):
            out.append((c.find("s3:Key", _NS).text, int(c.find("s3:Size", _NS).text)))
        if root.findtext("s3:IsTruncated", namespaces=_NS) == "true":
            token = root.findtext("s3:NextContinuationToken", namespaces=_NS)
        else:
            return out


_L2_RE = re.compile(r"([A-Z]{4})(\d{8})_(\d{6})(_V0\d)?(\.gz)?$")


def list_level2(site: str, day: dt.date) -> list[tuple[dt.datetime, str, int]]:
    """Volume scans for one radar and UTC day: (start time, key, bytes). Skips _MDM files."""
    keys = _list(L2_BUCKET, f"{day:%Y/%m/%d}/{site}/")
    out = []
    for key, size in keys:
        m = _L2_RE.search(key)
        if not m:
            continue  # *_MDM (model data message) and other auxiliary files
        t = dt.datetime.strptime(m.group(2) + m.group(3), "%Y%m%d%H%M%S").replace(tzinfo=dt.timezone.utc)
        out.append((t, key, size))
    return sorted(out)


def list_level3(site3: str, product: str, day: dt.date) -> list[tuple[dt.datetime, str, int]]:
    keys = _list(L3_BUCKET, f"{site3}_{product}_{day:%Y_%m_%d}")
    out = []
    for key, size in keys:
        t = dt.datetime.strptime(key[-19:], "%Y_%m_%d_%H_%M_%S").replace(tzinfo=dt.timezone.utc)
        out.append((t, key, size))
    return sorted(out)


def fetch(bucket: str, key: str) -> pathlib.Path:
    sub = "l2" if bucket == L2_BUCKET else "l3"
    path = CACHE / sub / key
    if path.exists() and path.stat().st_size > 0:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    r = _session.get(f"{bucket}/{key}", timeout=180)
    r.raise_for_status()
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_bytes(r.content)
    tmp.rename(path)
    return path


@dataclasses.dataclass
class Sweep:
    """One elevation sweep: moment arrays are (n_rays, n_gates), NaN where no data."""
    site: str
    time: dt.datetime            # time of the first ray
    elevation: float             # nominal elevation angle, degrees
    azimuth: np.ndarray          # degrees clockwise from north, per ray
    ray_time: np.ndarray         # seconds since `time`, per ray
    ranges: dict                 # moment -> gate centre ranges in metres
    data: dict                   # moment -> 2-D float array
    radar_lat: float
    radar_lon: float
    radar_height: float          # metres above sea level (site elevation + feedhorn)
    vcp: int | None = None
    duration_s: float | None = None   # sweep duration (Level III: product delta_time)
    az_start: float | None = None     # azimuth of the first radial; WSR-88D rotates clockwise

    def beam_time(self, az_deg):
        """Seconds after self.time at which the beam pointed at azimuth az_deg.
        Level II: from the recorded ray times; Level III: assumes uniform clockwise rotation
        from az_start over duration_s (checked against Level II: 71 s vs 71.8 s, KFTG)."""
        az = np.asarray(az_deg, float) % 360.0
        if self.duration_s is not None and self.az_start is not None:
            return ((az - self.az_start) % 360.0) / 360.0 * self.duration_s
        order = np.argsort(self.azimuth)
        a = self.azimuth[order]
        rt = self.ray_time[order]
        k = np.clip(np.searchsorted(a, az), 0, len(a) - 1)
        return rt[k]


def beam_height_ground_range(r, elev_deg, h0):
    """Height above sea level and ground range for slant range r (m) at elevation (deg)."""
    ka = K_EFF * EARTH_RADIUS
    th = np.radians(elev_deg)
    h = np.sqrt(r ** 2 + ka ** 2 + 2 * r * ka * np.sin(th)) - ka + h0
    s = ka * np.arcsin(r * np.cos(th) / (ka + h - h0))
    return h, s


def sample(sweep: Sweep, x, y, moment: str = "REF", max_az_gap: float = 1.0):
    """Nearest-gate value of a moment at plane coordinates (x east, y north, metres from
    the radar). Returns NaN outside the sweep. The ground range is converted to slant
    range by inverting the 4/3-earth geometry of the sweep's gates."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    s = np.hypot(x, y)
    az = np.degrees(np.arctan2(x, y)) % 360.0
    r = sweep.ranges[moment]
    _, s_gates = beam_height_ground_range(r, sweep.elevation, sweep.radar_height)
    dr = r[1] - r[0]
    gate = np.rint(np.interp(s, s_gates, np.arange(len(r)), left=-1, right=len(r))).astype(int)
    order = np.argsort(sweep.azimuth)
    a_sorted = sweep.azimuth[order]
    k = np.searchsorted(a_sorted, az)
    k0 = (k - 1) % len(a_sorted)
    k1 = k % len(a_sorted)
    d0 = np.abs(((az - a_sorted[k0]) + 180) % 360 - 180)
    d1 = np.abs(((az - a_sorted[k1]) + 180) % 360 - 180)
    ray = np.where(d0 <= d1, order[k0], order[k1])
    gap = np.minimum(d0, d1)
    data = sweep.data[moment]
    ok = (gate >= 0) & (gate < data.shape[1]) & (gap <= max_az_gap) & (s - s_gates[0] > -dr)
    out = np.full(x.shape, np.nan, np.float32)
    out[ok] = data[ray[ok], gate[ok]]
    return out


def to_grid(sweep: Sweep, half_width_m: float, res_m: float, moment: str = "REF"):
    """Resample a sweep onto a square grid centred on the radar: returns (xs, ys, Z) with
    Z[iy, ix] at (xs[ix], ys[iy])."""
    xs = np.arange(-half_width_m, half_width_m + res_m / 2, res_m)
    ys = xs.copy()
    X, Y = np.meshgrid(xs, ys)
    return xs, ys, sample(sweep, X, Y, moment)


def read_level2(path, moments=("REF", "ZDR", "RHO", "VEL", "SW", "PHI")) -> list[Sweep]:
    """Read a Level II volume into sweeps. Split-cut duplicates are kept as separate sweeps
    (at low elevations the volume coverage patterns scan once for reflectivity and once
    for Doppler; callers pick the sweep that carries the moment they need)."""
    from metpy.io import Level2File

    f = Level2File(str(path))
    site = f.stid.decode() if isinstance(f.stid, bytes) else str(f.stid)
    # each ray is (header, volume constants, elevation constants, radial constants, moments)
    vol = f.sweeps[0][0][1]
    lat, lon = float(vol.lat), float(vol.lon)
    h0 = float(vol.site_amsl) + float(vol.feedhorn_agl)
    vcp = int(vol.vcp)
    out = []
    for rays in f.sweeps:
        if not rays:
            continue
        az = np.array([r[0].az_angle for r in rays], float)
        # the median is the pointing angle; the first rays can lag while the antenna settles
        el = float(np.median([r[0].el_angle for r in rays]))
        # ray times: header has days since 1970-01-01 (1-based) and ms of day
        rt = np.array([(r[0].date - 1) * 86400.0 + r[0].time_ms / 1000.0 for r in rays])
        tstart = dt.datetime.fromtimestamp(rt[0], tz=dt.timezone.utc)
        data, ranges = {}, {}
        for mom in moments:
            key = mom.encode()
            if key not in rays[0][4]:
                continue
            hdr = rays[0][4][key][0]
            ngates = max(len(r[4][key][1]) for r in rays if key in r[4])
            arr = np.full((len(rays), ngates), np.nan, np.float32)
            for i, r in enumerate(rays):
                if key in r[4]:
                    v = r[4][key][1]
                    arr[i, : len(v)] = v
            data[mom] = arr
            ranges[mom] = (hdr.first_gate + hdr.gate_width * np.arange(ngates)) * 1000.0
        if not data:
            continue
        out.append(Sweep(site=site, time=tstart, elevation=el, azimuth=az, ray_time=rt - rt[0],
                         ranges=ranges, data=data, radar_lat=lat, radar_lon=lon,
                         radar_height=h0, vcp=vcp))
    return out


def read_level3(path) -> Sweep:
    """Read a Level III radial product (e.g. N0B) into a Sweep with moment 'REF'."""
    from metpy.io import Level3File

    f = Level3File(str(path))
    d = f.sym_block[0][0]
    start = np.asarray(d["start_az"], float)
    end = np.asarray(d["end_az"], float)
    width = (end - start) % 360.0
    az = (start + 0.5 * width) % 360.0
    values = np.asarray(f.map_data(np.asarray(d["data"])), np.float32)
    # bins of equal width from 0 km to max_range (N0B: 1840 x 0.25 km to 460 km), so bin 8
    # (2.0-2.25 km) lines up with the first Level II gate centred at 2.125 km
    ngates = values.shape[1]
    edges = np.linspace(0.0, float(f.max_range), ngates + 1)
    rng = 0.5 * (edges[:-1] + edges[1:]) * 1000.0
    t = f.metadata["vol_time"].replace(tzinfo=dt.timezone.utc)
    dur = f.metadata.get("delta_time")
    return Sweep(site=f.siteID, time=t, elevation=float(f.metadata.get("el_angle", np.nan)),
                 azimuth=az, ray_time=np.zeros(len(az)), ranges={"REF": rng},
                 data={"REF": values}, radar_lat=float(f.lat), radar_lon=float(f.lon),
                 radar_height=float(f.height) * 0.3048, vcp=f.prod_desc.vcp,
                 duration_s=float(dur) if dur else None, az_start=float(start[0]))
