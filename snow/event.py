"""Event-level test: did aircraft passages at an airport make radar echo appear?

For every track point (position x_i, y_i at time t_i, drift u_i, v_i) and every radar scan
at T with age a = T - t_i in the before or after window, the drifting air parcel is at
    P_i(T) = (x_i + u_i a, y_i + v_i a)
and the reflectivity there is looked up on the scan's Cartesian grid ("on"), and also at
P_i(T) +- CONTROL_OFFSET_M along the perpendicular to the track ("off": parallel lines no
aircraft flew along). For a set of points (an altitude band and a window of passage times)
and a statistic f of the reflectivity, the difference-in-differences is

    S = (f_on,after - f_on,before) - (f_off,after - f_off,before)

with "after" ages 15-60 min and "before" ages -30 to -5 min. Statistics (STATS): "z" the
mean reflectivity floored at FLOOR_DBZ (dB); "e0", "e5", "e10" the fraction of samples with
at least 0, 5, 10 dBZ. Seeding makes S positive; pre-existing precipitation, fixed clutter and
broad intensification cancel. The primary statistic, fixed before the census, is "z".

Sums and counts are accumulated on a fine grid of 1000 ft bands x 1 h passage bins and
summed into analysis cells (e.g. 2000 ft x 3 h), so overlapping cells cost nothing extra.

Null distribution: placebo track sets with every point rotated about the radar by an angle
in PLACEBO_ANGLES (same ranges, track shapes, traffic and wind; wrong places). Each cell's
z-score is (S - placebo mean) / placebo s.d. The search over cells and cell sizes is
accounted for family-wise: the largest real z-score is compared with the largest z-score
of each placebo, computed leave-one-out from the other placebos.
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pandas as pd

FLOOR_DBZ = -10.0
CONTROL_OFFSET_M = 6000.0
BEFORE = (-30 * 60, -5 * 60)
AFTER = (15 * 60, 60 * 60)
FINE_BAND_FT = 1000
BAND_LO_FT = 2000
BAND_HI_FT = 24000
STATS = ("z", "e0", "e5", "e10")
THRESH = {"e0": 0.0, "e5": 5.0, "e10": 10.0}
PRIMARY = "z"
CELL_SIZES = ((2, 1), (2, 3), (4, 3))          # (bands of 1000 ft, hours)
MIN_COUNT = 50
ROTATION_CHUNK = 32
PLACEBO_ANGLES = tuple(np.arange(20.0, 340.1, 2.0))   # 161 rotations, none within 20 deg
AIRCRAFT_MASK_M = 1500.0
RHO_MIN = 0.90            # correlation coefficient below this: not precipitation (Level III N0C)
AIRCRAFT_MATCH_S = 15.0   # positions this close in time to the beam's passage are masked


@dataclasses.dataclass
class Grids:
    """Radar scans resampled to a common Cartesian grid centred on the radar."""
    t: np.ndarray          # (n_scans,) epoch seconds
    z: np.ndarray          # (n_scans, ny, nx) float32, floored dBZ; NaN outside coverage
    x0: float              # coordinate of column 0, metres
    y0: float
    res: float

    @classmethod
    def from_sweeps(cls, scans, half_m: float = 70_000.0, res_m: float = 500.0,
                    min_range_m: float = 8000.0, aircraft: pd.DataFrame | None = None):
        """Resample sweeps to the grid. If `aircraft` is given (columns hex, t = heatmap slot
        start, x, y in the radar plane, gs = ground speed in 0.1 kt), every aircraft is masked
        where the beam actually saw it. For each reported position, f = t - (time the beam
        pointed at its azimuth, Sweep.beam_time). Where f changes sign between two
        consecutive positions of one aircraft (at most 60 s apart), the position is
        interpolated to f = 0 and a disk of AIRCRAFT_MASK_M plus 5 s of flight (the heatmap
        timing uncertainty) is masked. An aircraft seen only on one side of that moment is
        masked at its nearest position if within AIRCRAFT_MATCH_S, with the disk enlarged by
        the time difference. Aircraft echoes therefore never count as snow."""
        from scipy import ndimage

        from . import radar
        xs = np.arange(-half_m, half_m + res_m / 2, res_m)
        X, Y = np.meshgrid(xs, xs)
        R = np.hypot(X, Y)
        bad = (R < min_range_m) | (R > half_m)
        order = np.argsort([t for t, _ in scans])
        if aircraft is not None:
            ac = aircraft.sort_values(["t"]).reset_index(drop=True)
            ac_t = ac.t.to_numpy(float) + 5.0           # slot start + half a slot
            ac_x, ac_y = ac.x.to_numpy(), ac.y.to_numpy()
            ac_h = ac.hex.to_numpy()
            ac_v = np.clip(ac.gs.to_numpy(float), 0, None) * 0.05144   # 0.1 kt -> m/s
            ac_az = np.degrees(np.arctan2(ac_x, ac_y)) % 360.0
            disks = {}
        zs, ts = [], []
        for i in order:
            t, sw = scans[i]
            z = radar.sample(sw, X, Y)
            zf = np.where(np.isnan(z), FLOOR_DBZ, np.maximum(z, FLOOR_DBZ)).astype(np.float32)
            if getattr(sw, "rho", None) is not None:
                # birds, insects, ground clutter, wind turbines and anomalous propagation have
                # a low correlation coefficient; snow and rain are near 1. Such echo counts as
                # no echo (where the coefficient is missing, the reflectivity is kept)
                r = radar.sample(sw.rho, X, Y, "RHO")
                zf[np.isfinite(r) & (r < RHO_MIN)] = FLOOR_DBZ
            zf[bad] = np.nan
            if aircraft is not None:
                dur = sw.duration_s or 60.0
                lo = np.searchsorted(ac_t, t - 70.0)
                hi = np.searchsorted(ac_t, t + dur + 70.0, side="right")
                sel = np.arange(lo, hi)
                sel = sel[np.lexsort((ac_t[sel], ac_h[sel]))]          # by aircraft, then time
                f = ac_t[sel] - (t + sw.beam_time(ac_az[sel]))
                same = (ac_h[sel][1:] == ac_h[sel][:-1]) & (np.diff(ac_t[sel]) <= 60.0)
                cross = same & (f[:-1] <= 0) & (f[1:] > 0)
                j = np.nonzero(cross)[0]
                w = -f[j] / (f[j + 1] - f[j])
                a, b = sel[j], sel[j + 1]
                px = ac_x[a] + w * (ac_x[b] - ac_x[a])
                py = ac_y[a] + w * (ac_y[b] - ac_y[a])
                rad = AIRCRAFT_MASK_M + 5.0 * np.maximum(ac_v[a], ac_v[b])
                # aircraft not bracketed (track starts or ends near the beam's passage)
                bracketed = set(ac_h[a])
                lone = np.array([k for k in range(len(sel)) if ac_h[sel[k]] not in bracketed
                                 and abs(f[k]) <= AIRCRAFT_MATCH_S], int)
                if len(lone):
                    px = np.r_[px, ac_x[sel[lone]]]
                    py = np.r_[py, ac_y[sel[lone]]]
                    rad = np.r_[rad, AIRCRAFT_MASK_M + ac_v[sel[lone]] * (np.abs(f[lone]) + 5.0)]
                radius = np.minimum(rad, 6000.0)
                ix = np.rint((px - xs[0]) / res_m).astype(int)
                iy = np.rint((py - xs[0]) / res_m).astype(int)
                ok = (ix >= 0) & (ix < len(xs)) & (iy >= 0) & (iy < len(xs))
                hit = np.zeros(zf.shape, bool)
                classes = (2000.0, 3000.0, 4500.0, 6000.0)
                which = np.searchsorted(classes, radius)       # smallest disk covering radius
                for k, rmax in enumerate(classes):
                    cls_m = ok & (which == k)
                    if not cls_m.any():
                        continue
                    seed = np.zeros(zf.shape, bool)
                    seed[iy[cls_m], ix[cls_m]] = True
                    if rmax not in disks:
                        rr = int(np.ceil(rmax / res_m))
                        yy, xx = np.mgrid[-rr:rr + 1, -rr:rr + 1]
                        disks[rmax] = np.hypot(yy, xx) * res_m <= rmax
                    hit |= ndimage.binary_dilation(seed, structure=disks[rmax])
                zf[hit] = np.nan
            zs.append(zf)
            ts.append(t)
        return cls(np.array(ts, float), np.stack(zs), float(xs[0]), float(xs[0]), res_m)

    def lookup(self, k, px, py):
        ix = np.rint((px - self.x0) / self.res).astype(np.int64)
        iy = np.rint((py - self.y0) / self.res).astype(np.int64)
        ny, nx = self.z.shape[1:]
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        out = np.full(px.shape, np.nan, np.float32)
        out[ok] = self.z[k, iy[ok], ix[ok]]
        return out


def perpendicular(points: pd.DataFrame) -> np.ndarray:
    from .epoch import heading
    h = heading(points)
    return np.c_[-h[:, 1], h[:, 0]]


def _transforms(z):
    """Stack of per-sample statistics in STATS order."""
    return np.stack([z] + [(z >= THRESH[s]).astype(np.float64) for s in STATS[1:]])


def did_table(grids: Grids, x, y, t, drift, perp, band, hour, n_band, n_hour,
              offset_m=CONTROL_OFFSET_M):
    """Accumulated sums: array (n_band, n_hour, 2 windows, 2 lines, len(STATS) + 1). The
    last slot is the sample count; windows 0=before 1=after; lines 0=on 1=off."""
    ns = len(STATS)
    acc = np.zeros((n_band * n_hour * 4, ns + 1))
    order = np.argsort(t)
    ts = t[order]
    for k, T in enumerate(grids.t):
        lo = np.searchsorted(ts, T - AFTER[1])
        hi = np.searchsorted(ts, T - BEFORE[0], side="right")
        if hi <= lo:
            continue
        idx = order[lo:hi]
        age = T - t[idx]
        win = np.full(len(idx), -1)
        win[(age >= BEFORE[0]) & (age <= BEFORE[1])] = 0
        win[(age >= AFTER[0]) & (age <= AFTER[1])] = 1
        m = win >= 0
        if not m.any():
            continue
        idx, age, win = idx[m], age[m], win[m]
        px = x[idx] + drift[idx, 0] * age
        py = y[idx] + drift[idx, 1] * age
        z_on = grids.lookup(k, px, py)
        ox, oy = offset_m * perp[idx, 0], offset_m * perp[idx, 1]
        z_off = [grids.lookup(k, px + ox, py + oy), grids.lookup(k, px - ox, py - oy)]
        cell = (band[idx] * n_hour + hour[idx]) * 2 + win
        # on-line samples
        ok = np.isfinite(z_on)
        if ok.any():
            flat = cell[ok] * 2 + 0
            T_on = _transforms(z_on[ok].astype(np.float64))
            for s in range(ns):
                acc[:, s] += np.bincount(flat, weights=T_on[s], minlength=len(acc))
            acc[:, ns] += np.bincount(flat, minlength=len(acc))
        # off-line samples: each side counts as a sample
        for zo in z_off:
            ok = np.isfinite(zo)
            if not ok.any():
                continue
            flat = cell[ok] * 2 + 1
            T_off = _transforms(zo[ok].astype(np.float64))
            for s in range(ns):
                acc[:, s] += np.bincount(flat, weights=T_off[s], minlength=len(acc))
            acc[:, ns] += np.bincount(flat, minlength=len(acc))
    return acc.reshape(n_band, n_hour, 2, 2, ns + 1)


NEW_ECHO_DBZ = 5.0       # "new echo": at least this after passage ...
QUIET_DBZ = 0.0          # ... and below this at every scan before passage


def parcel_maxima(grids: Grids, x, y, t, drift, perp, offset_m=CONTROL_OFFSET_M):
    """Per point and line (on, off+, off-): the maximum reflectivity over scans in the
    before window and in the after window, and how many scans were seen in each.
    Returns arrays (n_points, 3) for max_before, max_after, n_before, n_after."""
    n = len(x)
    mb = np.full((n, 3), -np.inf)
    ma = np.full((n, 3), -np.inf)
    nb = np.zeros((n, 3), np.int32)
    na = np.zeros((n, 3), np.int32)
    order = np.argsort(t)
    ts = t[order]
    for k, T in enumerate(grids.t):
        lo = np.searchsorted(ts, T - AFTER[1])
        hi = np.searchsorted(ts, T - BEFORE[0], side="right")
        if hi <= lo:
            continue
        idx = order[lo:hi]
        age = T - t[idx]
        bef = (age >= BEFORE[0]) & (age <= BEFORE[1])
        aft = (age >= AFTER[0]) & (age <= AFTER[1])
        m = bef | aft
        if not m.any():
            continue
        idx, bef, aft, age = idx[m], bef[m], aft[m], age[m]
        px = x[idx] + drift[idx, 0] * age
        py = y[idx] + drift[idx, 1] * age
        ox, oy = offset_m * perp[idx, 0], offset_m * perp[idx, 1]
        for line, (dx, dy) in enumerate(((0, 0), (ox, oy), (-ox, -oy))):
            z = grids.lookup(k, px + dx, py + dy)
            ok = np.isfinite(z)
            for win, wm, mx, cnt in ((0, bef, mb, nb), (1, aft, ma, na)):
                s = ok & wm
                ii = idx[s]
                mx[ii, line] = np.maximum(mx[ii, line], z[s])
                cnt[ii, line] += 1
    return mb, ma, nb, na


def new_echo_table(grids, x, y, t, drift, perp, band, hour, n_band, n_hour,
                   min_before=2, min_after=3):
    """Accumulated counts for the new-echo statistic: array (n_band, n_hour, 3 lines, 2)
    holding (number of parcels with new echo, number of parcels evaluated)."""
    mb, ma, nb, na = parcel_maxima(grids, x, y, t, drift, perp)
    acc = np.zeros((n_band * n_hour, 3, 2))
    cell = band * n_hour + hour
    for line in range(3):
        ev = (nb[:, line] >= min_before) & (na[:, line] >= min_after)
        new = ev & (ma[:, line] >= NEW_ECHO_DBZ) & (mb[:, line] < QUIET_DBZ)
        acc[:, line, 0] = np.bincount(cell[new], minlength=len(acc))
        acc[:, line, 1] = np.bincount(cell[ev], minlength=len(acc))
    return acc.reshape(n_band, n_hour, 3, 2)


def new_echo_did(acc_cells, min_count=MIN_COUNT):
    """S_new = P(new echo on the track) - P(new echo on the side lines) per cell."""
    with np.errstate(invalid="ignore", divide="ignore"):
        p_on = acc_cells[..., 0, 0] / acc_cells[..., 0, 1]
        p_off = (acc_cells[..., 1, 0] + acc_cells[..., 2, 0]) / (acc_cells[..., 1, 1] + acc_cells[..., 2, 1])
    S = p_on - p_off
    S[acc_cells[..., 0, 1] < min_count] = np.nan
    return S


def tables_rotated(grids: Grids, x, y, t, drift, perp, band, hour, n_band, n_hour, angles,
                   offset_m=CONTROL_OFFSET_M, min_before=2, min_after=3):
    """did_table and new_echo_table for many rotations of the same points at once.

    angles: rotation angles in degrees (0 = the real tracks). Returns (acc, new) with shapes
    (n_rot, n_band, n_hour, 2, 2, len(STATS) + 1) and (n_rot, n_band, n_hour, 3, 2), equal
    to calling did_table / new_echo_table on each rotated point set."""
    ns = len(STATS)
    R = len(angles)
    ncell = n_band * n_hour
    a = np.radians(np.asarray(angles, float))
    cs, sn = np.cos(a)[:, None], np.sin(a)[:, None]
    acc = np.zeros((R * ncell * 4, ns + 1))
    n = len(x)
    mb = np.full((R, n, 3), -np.inf, np.float32)
    ma = np.full((R, n, 3), -np.inf, np.float32)
    nb = np.zeros((R, n, 3), np.int16)
    na = np.zeros((R, n, 3), np.int16)
    order = np.argsort(t)
    ts = t[order]
    cell_pt = band * n_hour + hour
    rr = np.arange(R)[:, None]
    for k, T in enumerate(grids.t):
        lo = np.searchsorted(ts, T - AFTER[1])
        hi = np.searchsorted(ts, T - BEFORE[0], side="right")
        if hi <= lo:
            continue
        idx = order[lo:hi]
        age = T - t[idx]
        bef = (age >= BEFORE[0]) & (age <= BEFORE[1])
        aft = (age >= AFTER[0]) & (age <= AFTER[1])
        m = bef | aft
        if not m.any():
            continue
        idx, age, bef, aft = idx[m], age[m], bef[m], aft[m]
        win = aft.astype(np.int64)                                   # 0 before, 1 after
        X = x[idx] * cs - y[idx] * sn                                # (R, n_k)
        Y = x[idx] * sn + y[idx] * cs
        PX = perp[idx, 0] * cs - perp[idx, 1] * sn
        PY = perp[idx, 0] * sn + perp[idx, 1] * cs
        px = X + drift[idx, 0] * age
        py = Y + drift[idx, 1] * age
        base = (rr * ncell + cell_pt[idx]) * 2 + win                 # (R, n_k)
        for line, sgn in enumerate((0.0, 1.0, -1.0)):
            z = grids.lookup(k, px + sgn * offset_m * PX, py + sgn * offset_m * PY)
            ok = np.isfinite(z)
            # difference-in-differences sums (lines 1 and 2 both go to the "off" slot)
            flat = (base * 2 + (0 if line == 0 else 1))[ok]
            zz = z[ok].astype(np.float64)
            acc[:, 0] += np.bincount(flat, weights=zz, minlength=len(acc))
            for s_i, s_name in enumerate(STATS[1:], start=1):
                acc[:, s_i] += np.bincount(flat, weights=(zz >= THRESH[s_name]), minlength=len(acc))
            acc[:, ns] += np.bincount(flat, minlength=len(acc))
            # per-parcel maxima for the new-echo statistic
            zf = np.where(ok, z, -np.inf).astype(np.float32)
            ib = np.nonzero(bef)[0]
            ia = np.nonzero(aft)[0]
            if len(ib):
                g = idx[ib]
                mb[:, g, line] = np.maximum(mb[:, g, line], zf[:, ib])
                nb[:, g, line] += ok[:, ib]
            if len(ia):
                g = idx[ia]
                ma[:, g, line] = np.maximum(ma[:, g, line], zf[:, ia])
                na[:, g, line] += ok[:, ia]
    acc = acc.reshape(R, n_band, n_hour, 2, 2, ns + 1)
    new = np.zeros((R, ncell, 3, 2))
    for line in range(3):
        ev = (nb[:, :, line] >= min_before) & (na[:, :, line] >= min_after)
        hit = ev & (ma[:, :, line] >= NEW_ECHO_DBZ) & (mb[:, :, line] < QUIET_DBZ)
        for r in range(R):
            new[r, :, line, 0] = np.bincount(cell_pt[hit[r]], minlength=ncell)
            new[r, :, line, 1] = np.bincount(cell_pt[ev[r]], minlength=ncell)
    return acc, new.reshape(R, n_band, n_hour, 3, 2)


def cells(acc, bw: int, hw: int):
    """Sum fine accumulators over windows of bw bands and hw hours (sliding, step 1)."""
    c = np.cumsum(np.cumsum(acc, axis=0), axis=1)
    c = np.pad(c, ((1, 0), (1, 0)) + ((0, 0),) * (acc.ndim - 2))
    return c[bw:, hw:] - c[:-bw, hw:] - c[bw:, :-hw] + c[:-bw, :-hw]


def did(acc_cells, min_count=MIN_COUNT):
    """S for every cell and statistic: array (nb, nh, len(STATS)); NaN when the smallest of
    the four (window, line) sample counts is below min_count."""
    cnt = acc_cells[..., -1]
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = acc_cells[..., :-1] / cnt[..., None]
    S = (mean[:, :, 1, 0] - mean[:, :, 0, 0]) - (mean[:, :, 1, 1] - mean[:, :, 0, 1])
    nmin = cnt.reshape(cnt.shape[0], cnt.shape[1], -1).min(axis=-1)
    S[nmin < min_count] = np.nan
    return S


def _fine_index(points: pd.DataFrame, t_start: float, n_hour: int, hour_s: float):
    alt = points.alt_ft.to_numpy()
    edges = np.arange(BAND_LO_FT, BAND_HI_FT + 1, FINE_BAND_FT)
    band = np.digitize(alt, edges) - 1
    hour = np.floor((points.t.to_numpy(float) - t_start) / hour_s).astype(int)
    keep = (band >= 0) & (band < len(edges) - 1) & (hour >= 0) & (hour < n_hour)
    return band, hour, keep, edges


def run(grids: Grids, points: pd.DataFrame, drift: np.ndarray, t_start: float, n_hour: int,
        angles=PLACEBO_ANGLES, hour_s: float = 3600.0) -> dict:
    """Fine accumulators for the real tracks and for each rotated placebo."""
    band, hour, keep, edges = _fine_index(points, t_start, n_hour, hour_s)
    keep &= np.isfinite(drift[:, 0])
    pts = points[keep]
    x, y, t = pts.x.to_numpy(), pts.y.to_numpy(), pts.t.to_numpy(float)
    band, hour, dr = band[keep], hour[keep], drift[keep]
    perp = perpendicular(pts)
    nb = len(edges) - 1
    all_angles = np.r_[0.0, np.asarray(angles, float)]
    accs, news = [], []
    for c0 in range(0, len(all_angles), ROTATION_CHUNK):
        a_c, n_c = tables_rotated(grids, x, y, t, dr, perp, band, hour, nb, n_hour,
                                  all_angles[c0:c0 + ROTATION_CHUNK])
        accs.append(a_c)
        news.append(n_c)
    acc = np.concatenate(accs)
    new = np.concatenate(news)
    return {"acc_real": acc[0], "acc_placebo": acc[1:], "new_real": new[0],
            "new_placebo": new[1:], "band_edges": edges,
            "t_start": t_start, "hour_s": hour_s, "angles": np.asarray(angles)}


CLUSTER_Z = 2.0
CLUSTER_SIZE = (2, 1)          # 2000 ft x 1 h cells, sliding by 1000 ft


def _zmaps(res: dict, stat: str, bw: int, hw: int):
    """Real z-map and leave-one-out placebo z-maps for one statistic and cell size.
    stat is one of STATS or "new" (the per-parcel new-echo statistic)."""
    if stat == "new":
        S_r = new_echo_did(cells(res["new_real"], bw, hw))
        S_p = np.stack([new_echo_did(cells(a, bw, hw)) for a in res["new_placebo"]])
    else:
        si = STATS.index(stat)
        S_r = did(cells(res["acc_real"], bw, hw))[..., si]
        S_p = np.stack([did(cells(a, bw, hw))[..., si] for a in res["acc_placebo"]])
    n_p = S_p.shape[0]
    valid = np.isfinite(S_r) & (np.isfinite(S_p).sum(axis=0) >= max(10, n_p // 2))
    S_p = np.where(valid, S_p, np.nan)
    mu = np.nanmean(S_p, axis=0)
    sd = np.nanstd(S_p, axis=0, ddof=1)
    z_r = np.where(valid, (S_r - mu) / sd, np.nan)
    # leave-one-out: mean/sd of the other n-1 placebos from running sums
    cnt = np.isfinite(S_p).sum(axis=0)
    s1 = np.nansum(S_p, axis=0)
    s2 = np.nansum(S_p ** 2, axis=0)
    z_p = np.full(S_p.shape, np.nan)
    for i in range(n_p):
        xi = S_p[i]
        ok = np.isfinite(xi)
        c = cnt - ok
        m = (s1 - np.where(ok, xi, 0)) / c
        v = (s2 - np.where(ok, xi ** 2, 0) - c * m ** 2) / (c - 1)
        z_p[i] = (xi - m) / np.sqrt(v)
    return z_r, z_p, S_r, sd


def cluster_masses(z: np.ndarray, thr: float = CLUSTER_Z):
    """Connected clusters (4-neighbour, band x hour) of cells with z > thr. Returns the
    label image and a list of (mass, cells) with mass = sum of (z - thr) over the cluster."""
    from scipy import ndimage
    m = np.isfinite(z) & (z > thr)
    lab, n = ndimage.label(m)
    out = []
    for k in range(1, n + 1):
        sel = lab == k
        out.append((float(np.sum(z[sel] - thr)), np.argwhere(sel)))
    return lab, out


def evaluate_clusters(res: dict, stat: str = PRIMARY, size=CLUSTER_SIZE,
                      thr: float = CLUSTER_Z) -> dict:
    """Cluster-mass test: the largest real cluster against the largest cluster of each
    placebo (whose z-maps are leave-one-out). p = (1 + #placebo >= real) / (1 + n)."""
    bw, hw = size
    z_r, z_p, S_r, sd = _zmaps(res, stat, bw, hw)
    _, cl = cluster_masses(z_r, thr)
    best = max(cl, key=lambda c: c[0]) if cl else (0.0, np.empty((0, 2), int))
    pmax = np.array([max([c[0] for c in cluster_masses(zi, thr)[1]] or [0.0]) for zi in z_p])
    edges = res["band_edges"]
    cells_idx = best[1]
    info = {"mass": best[0], "p_cluster": float((1 + np.sum(pmax >= best[0])) / (1 + len(pmax))),
            "placebo_mass_q95": float(np.quantile(pmax, 0.95)), "n_cells": int(len(cells_idx))}
    if len(cells_idx):
        b0, b1 = cells_idx[:, 0].min(), cells_idx[:, 0].max()
        h0, h1 = cells_idx[:, 1].min(), cells_idx[:, 1].max()
        info.update({"band_lo": int(edges[b0]), "band_hi": int(edges[b1] + bw * FINE_BAND_FT),
                     "hour_first": int(h0), "hour_last": int(h1 + hw - 1),
                     "z_peak": float(np.nanmax(z_r[cells_idx[:, 0], cells_idx[:, 1]])),
                     "S_mean": float(np.nanmean(S_r[cells_idx[:, 0], cells_idx[:, 1]]))})
    info["z_map"] = z_r
    info["placebo_masses"] = pmax
    return info


def evaluate(res: dict, stat: str = PRIMARY, sizes=CELL_SIZES) -> dict:
    """z-scores of every cell for each cell size, the best real cell, and the family-wise
    p-value of the best real z against the placebos' own best (leave-one-out)."""
    si = STATS.index(stat)
    out = {"cells": [], "zmax": -np.inf, "best": None}
    n_p = res["acc_placebo"].shape[0]
    zmax_p = np.full(n_p, -np.inf)
    edges = res["band_edges"]
    for bw, hw in sizes:
        S_r = did(cells(res["acc_real"], bw, hw))[..., si]
        S_p = np.stack([did(cells(a, bw, hw))[..., si] for a in res["acc_placebo"]])
        valid = np.isfinite(S_r) & (np.isfinite(S_p).sum(axis=0) >= max(10, n_p // 2))
        mu = np.nanmean(np.where(valid, S_p, np.nan), axis=0)
        sd = np.nanstd(np.where(valid, S_p, np.nan), axis=0, ddof=1)
        z_r = np.where(valid, (S_r - mu) / sd, np.nan)
        for i in range(n_p):
            others = np.delete(S_p, i, axis=0)
            mu_i = np.nanmean(others, axis=0)
            sd_i = np.nanstd(others, axis=0, ddof=1)
            zi = np.where(valid, (S_p[i] - mu_i) / sd_i, np.nan)
            if np.isfinite(zi).any():
                zmax_p[i] = max(zmax_p[i], np.nanmax(zi))
        for b, h in zip(*np.nonzero(np.isfinite(z_r))):
            out["cells"].append({"bw": bw, "hw": hw, "band_lo": int(edges[b]),
                                 "band_hi": int(edges[b] + bw * FINE_BAND_FT),
                                 "hour0": int(h), "S": float(S_r[b, h]), "z": float(z_r[b, h]),
                                 "placebo_sd": float(sd[b, h])})
            if z_r[b, h] > out["zmax"]:
                out["zmax"] = float(z_r[b, h])
                out["best"] = out["cells"][-1]
    out["cells"] = pd.DataFrame(out["cells"])
    out["zmax_placebo"] = zmax_p
    out["p_familywise"] = float((1 + np.sum(zmax_p >= out["zmax"])) / (1 + n_p))
    return out
