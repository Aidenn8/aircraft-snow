"""Inject synthetic aircraft snow streaks into real radar grids (for recall estimates).

A seeding flight's points (x_i, y_i, t_i) in a chosen altitude band leave ice in their air,
which drifts with (u_i + e_u, v_i + e_v): the HRRR drift plus a per-event wind error drawn
from N(0, WIND_ERR_MS) in each component, so the detector does not get the exact wind it
uses. At scan time T, with age a = T - t_i, the streak's reflectivity at distance d from the
drifted point is

    Z_add(d, a) = A * f(a) * exp(-d^2 / (2 s(a)^2))          (linear units, mm^6 m^-3)
    s(a) = S0 + SPREAD * (a - ONSET)                          (slow widening)
    f(a) = 0 before ONSET, rising linearly to 1 at PEAK, then decaying linearly to 0 at END

with A set from a peak reflectivity in dBZ. Streak echo is added to the existing echo in
linear units, so injected snow on top of real snow brightens it as it would in nature.
The shape follows the Denver 25-26 Jan 2025 streaks: detected streak objects there had
mean widths of 1.2-3.5 km (median 2.2 km) at 15-60 min, which for a Gaussian profile cut
6 dB below its peak means s ~ 0.6-0.9 km; they were visible for about an hour. (Hole-punch
holes in cloud widen faster, 1-2 m/s, Iwasaki & Fujii 2025, but the snow falling out of
them stays narrower.)

Physical timing (when the point's height above the beam, dh, and a fall speed v are given):
ice from a point dh metres above the beam centre reaches the beam after a* = dh / v, so
each point has its own envelope, a triangle peaking at a* with half-width
max(FALL_SPREAD * a*, MIN_HALF_S) (a spread of fall speeds), and nothing before GROWTH_S
(the particles need time to grow to radar-visible size). A descending aircraft therefore
leaves a streak whose parts appear at different ages, as in nature. Points at or below the
beam leave nothing.
"""
from __future__ import annotations

import numpy as np

from .event import Grids, FLOOR_DBZ

ONSET_S = 15 * 60
PEAK_S = 35 * 60
END_S = 90 * 60
S0_M = 500.0
SPREAD_MS = 0.2
WIND_ERR_MS = 0.7
GROWTH_S = 10 * 60
FALL_SPREAD = 0.3
MIN_HALF_S = 10 * 60


def envelope(age):
    age = np.asarray(age, float)
    up = np.clip((age - ONSET_S) / (PEAK_S - ONSET_S), 0, 1)
    down = np.clip((END_S - age) / (END_S - PEAK_S), 0, 1)
    return np.where(age < PEAK_S, up, down) * (age >= ONSET_S) * (age <= END_S)


def envelope_fall(age, arrival):
    """Per-point envelope for physical timing: a triangle peaking at the arrival age at the
    beam, zero before GROWTH_S; arrival is NaN (or <= 0) for points that never reach it."""
    age = np.asarray(age, float)
    a = np.asarray(arrival, float)
    half = np.maximum(FALL_SPREAD * a, MIN_HALF_S)
    w = np.clip(1 - np.abs(age - a) / half, 0, 1)
    return np.where(np.isfinite(a) & (a > 0) & (age >= GROWTH_S), w, 0.0)


def inject(grids: Grids, x, y, t, drift, peak_dbz: float, rng=None, wind_err_ms=WIND_ERR_MS,
           dh=None, v_fall=None):
    """Return a copy of grids with a streak along the drifted points (x, y, t, drift).
    Also returns the wind error used. With dh (m above the beam, per point) and v_fall
    (m/s), each point's echo follows the physical timing described above."""
    rng = np.random.default_rng() if rng is None else rng
    err = rng.normal(0, wind_err_ms, 2)
    out = Grids(grids.t.copy(), grids.z.copy(), grids.x0, grids.y0, grids.res)
    A = 10 ** (peak_dbz / 10.0)
    ny, nx = out.z.shape[1:]
    xs = out.x0 + out.res * np.arange(nx)
    ys = out.y0 + out.res * np.arange(ny)
    physical = dh is not None and v_fall is not None
    if physical:
        dh = np.asarray(dh, float)
        arrival = np.where(dh > 0, dh / v_fall, np.nan)
    for k, T in enumerate(out.t):
        age = T - t
        w = envelope_fall(age, arrival) if physical else envelope(age)
        m = w > 0
        if not m.any():
            continue
        px = x[m] + (drift[m, 0] + err[0]) * age[m]
        py = y[m] + (drift[m, 1] + err[1]) * age[m]
        s = S0_M + SPREAD_MS * np.maximum(age[m] - (GROWTH_S if physical else ONSET_S), 0)
        add = np.zeros((ny, nx))
        r = 3.0 * s.max()
        for px_i, py_i, s_i, w_i in zip(px, py, s, w[m]):
            ix0 = max(int((px_i - r - xs[0]) / out.res), 0)
            ix1 = min(int((px_i + r - xs[0]) / out.res) + 1, nx)
            iy0 = max(int((py_i - r - ys[0]) / out.res), 0)
            iy1 = min(int((py_i + r - ys[0]) / out.res) + 1, ny)
            if ix0 >= ix1 or iy0 >= iy1:
                continue
            dx = xs[ix0:ix1][None, :] - px_i
            dy = ys[iy0:iy1][:, None] - py_i
            g = A * w_i * np.exp(-(dx ** 2 + dy ** 2) / (2 * s_i ** 2))
            # a line of points: take the maximum, not the sum, so dense sampling along the
            # track does not inflate the streak's intensity
            add[iy0:iy1, ix0:ix1] = np.maximum(add[iy0:iy1, ix0:ix1], g)
        z = out.z[k]
        valid = np.isfinite(z)
        lin = np.where(valid, 10 ** (np.where(valid, z, FLOOR_DBZ) / 10.0), 0.0)
        # below-threshold background counts as no echo
        lin = np.where(valid & (z > FLOOR_DBZ), lin, 0.0)
        tot = lin + add
        with np.errstate(divide="ignore"):
            zn = np.where(tot > 0, 10 * np.log10(tot), FLOOR_DBZ)
        zn = np.maximum(zn, FLOOR_DBZ).astype(np.float32)
        out.z[k] = np.where(valid, zn, np.nan)
    return out, err
