"""Parcel-following superposed-epoch analysis of radar reflectivity around aircraft passages.

Idea: if an aircraft seeds a supercooled cloud, ice appears in the air it flew through and
moves with the wind. So follow each track point's air parcel forward and backward in time
and record the radar reflectivity there. A seeded parcel is echo-free before the aircraft
arrives and develops echo some 10-20 minutes after. Pre-existing precipitation shows echo
before and after; ground clutter and other fixed echoes do not move with the wind.

For a track point i at (x_i, y_i, t_i) and a radar scan at time T, the parcel is at
    x_i + u_i (T - t_i),  y_i + v_i (T - t_i)
with (u_i, v_i) the drift velocity: the mean HRRR wind between the radar beam and the
seeding level, since ice that reaches the beam has fallen through that layer.

Control parcels start on lines offset sideways from the track (perpendicular to the local
heading, by +-CONTROL_OFFSET_M) and drift the same way. Their before/after change measures
any background evolution, e.g. snowfall already intensifying over the whole area.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import radar

ECHO_FLOOR_DBZ = -10.0          # below-threshold gates count as this value
CONTROL_OFFSET_M = 6000.0


def heading(points: pd.DataFrame) -> np.ndarray:
    """Unit vectors along each flight's track (central differences within a flight)."""
    out = np.zeros((len(points), 2))
    for _, idx in points.groupby("flight").indices.items():
        x = points.x.to_numpy()[idx]
        y = points.y.to_numpy()[idx]
        if len(idx) < 2:
            continue
        dx = np.gradient(x)
        dy = np.gradient(y)
        n = np.hypot(dx, dy)
        n[n == 0] = 1.0
        out[idx, 0] = dx / n
        out[idx, 1] = dy / n
    return out


def layer_mean_wind(profile: pd.DataFrame, p_top, p_bottom):
    """Mean HRRR u, v over pressure levels between p_top and p_bottom (hPa) for one run,
    weighted equally by level (levels are equally spaced in pressure)."""
    p = profile.p.to_numpy(float)
    u = profile.UGRD.to_numpy(float)
    v = profile.VGRD.to_numpy(float)
    p_top = np.atleast_1d(np.asarray(p_top, float))
    p_bottom = np.atleast_1d(np.asarray(p_bottom, float))
    uu = np.empty(len(p_top))
    vv = np.empty(len(p_top))
    for k, (a, b) in enumerate(zip(np.minimum(p_top, p_bottom), np.maximum(p_top, p_bottom))):
        m = (p >= a - 1e-6) & (p <= b + 1e-6)
        if not m.any():
            j = np.argmin(np.abs(p - 0.5 * (a + b)))
            m = np.zeros_like(p, bool)
            m[j] = True
        uu[k] = np.nanmean(u[m])
        vv[k] = np.nanmean(v[m])
    return uu, vv


def drift_velocity(points: pd.DataFrame, profiles: pd.DataFrame, radar_height_m: float,
                   beam_elev_deg: float) -> np.ndarray:
    """Drift (u, v) per point: HRRR layer-mean wind between the radar beam height at the
    point's ground range and the point's own pressure level, from the nearest hourly run."""
    from .radar import beam_height_ground_range
    runs = np.array(sorted(profiles.run.unique()))
    run_s = np.array([pd.Timestamp(r).timestamp() for r in runs])
    rng = np.hypot(points.x.to_numpy(), points.y.to_numpy())
    h_beam, _ = beam_height_ground_range(rng, beam_elev_deg, radar_height_m)
    out = np.full((len(points), 2), np.nan)
    which = np.abs(points.t.to_numpy(float)[:, None] - run_s[None, :]).argmin(axis=1)
    for k, r in enumerate(runs):
        m = which == k
        if not m.any():
            continue
        prof = profiles[profiles.run == r].sort_values("p")
        hgt = prof.HGT.to_numpy(float)
        p = prof.p.to_numpy(float)
        ok = np.isfinite(hgt)
        # pressure at beam height by log-linear interpolation in height (hgt rises as p falls)
        order = np.argsort(hgt[ok])
        p_beam = np.exp(np.interp(h_beam[m], hgt[ok][order], np.log(p[ok][order])))
        u, v = layer_mean_wind(prof, points.p_hpa.to_numpy()[m], p_beam)
        out[m, 0], out[m, 1] = u, v
    return out


def sample_parcels(points: pd.DataFrame, scans: list[tuple[float, "radar.Sweep"]],
                   drift: np.ndarray, age_min_s: float = -1800, age_max_s: float = 5400,
                   offsets_m=(0.0, CONTROL_OFFSET_M, -CONTROL_OFFSET_M)) -> pd.DataFrame:
    """Radar reflectivity at each point's drifting parcel for every scan with an age in
    [age_min_s, age_max_s]. One row per (point, scan, offset)."""
    hd = heading(points)
    perp = np.c_[-hd[:, 1], hd[:, 0]]           # left-perpendicular unit vector
    x0 = points.x.to_numpy()
    y0 = points.y.to_numpy()
    t0 = points.t.to_numpy(float)
    rows = []
    for T, sw in scans:
        age = T - t0
        m = (age >= age_min_s) & (age <= age_max_s) & np.isfinite(drift[:, 0])
        if not m.any():
            continue
        idx = np.where(m)[0]
        for off in offsets_m:
            px = x0[idx] + off * perp[idx, 0] + drift[idx, 0] * age[idx]
            py = y0[idx] + off * perp[idx, 1] + drift[idx, 1] * age[idx]
            z = radar.sample(sw, px, py)
            rows.append(pd.DataFrame({"point": idx, "scan_t": T, "age": age[idx], "offset": off,
                                      "z": z}))
    if not rows:
        return pd.DataFrame(columns=["point", "scan_t", "age", "offset", "z"])
    out = pd.concat(rows, ignore_index=True)
    out["z_floor"] = np.where(np.isnan(out.z), ECHO_FLOOR_DBZ, np.maximum(out.z, ECHO_FLOOR_DBZ))
    out["echo"] = (out.z >= 0.0).astype(float)
    return out


def composite(samples: pd.DataFrame, bin_s: float = 300.0) -> pd.DataFrame:
    """Mean floored reflectivity and echo fraction by age bin, for track and control lines."""
    s = samples.copy()
    s["age_bin"] = np.floor(s.age / bin_s) * bin_s + bin_s / 2
    s["line"] = np.where(s.offset == 0, "track", "control")
    g = s.groupby(["line", "age_bin"])
    return g.agg(z=("z_floor", "mean"), echo=("echo", "mean"), n=("z_floor", "size")).reset_index()


def before_after(samples: pd.DataFrame, before=(-1800, -300), after=(900, 2700)) -> dict:
    """Difference-in-differences of mean floored reflectivity: (after - before) on the
    track minus (after - before) on the control lines. Also echo-fraction version."""
    out = {}
    for name, col in (("dz", "z_floor"), ("decho", "echo")):
        def mean(line_mask, lo, hi):
            m = line_mask & (samples.age >= lo) & (samples.age < hi)
            return samples.loc[m, col].mean()
        track = samples.offset == 0
        ctrl = ~track
        t_b, t_a = mean(track, *before), mean(track, *after)
        c_b, c_a = mean(ctrl, *before), mean(ctrl, *after)
        out[name] = (t_a - t_b) - (c_a - c_b)
        out[name + "_track"] = t_a - t_b
        out[name + "_control"] = c_a - c_b
    return out
