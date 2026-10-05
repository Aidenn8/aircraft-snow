"""Segment test: does a line of radar echo appear where an aircraft's air has drifted?

Each flight is cut into short track segments (SEG_LEN_M). For every radar scan, a segment's
points are moved with the drift velocity for the time since the aircraft passed (the
"line"), and the same shape is also sampled at DECOYS: copies displaced by 6-10 km in eight
directions, which see the same kind of background but not the aircraft's air.

Per (segment, scan):
    on      mean reflectivity on the line (dBZ, floored at ECHO_FLOOR_DBZ)
    dec_med median over decoys of their mean reflectivity
    rank    fraction of decoys darker than the line (0..1)
    c       on - dec_med

Per segment, over scans 15-60 min after passage ("after") and 5-30 min before ("before"):
    clear-air hit  (aircraft seeding a non-precipitating cloud, the Denver regime)
        on >= CLEAR_MIN_DBZ, dec_med <= CLEAR_MAX_BG_DBZ and c >= CLEAR_MIN_C, and before
        passage the drifting air had on < CLEAR_BEFORE_MAX_DBZ
    rank hit (any background, including enhancement of existing snow)
        rank >= RANK_HIT and c >= RANK_MIN_C
A segment's persistence is the longest run of consecutive scans that are hits.

Placebos use the same machinery with altered tracks: passage times shifted by hours, or
positions rotated about the radar (same ranges, traffic density and track shapes, but
where no aircraft flew).
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pandas as pd

from . import epoch, radar

SEG_LEN_M = 2000.0
ECHO_FLOOR_DBZ = -10.0
DECOY_RADII_M = (6000.0, 10000.0)
DECOY_DIRS = 8
CLEAR_MIN_DBZ = 3.0
CLEAR_MAX_BG_DBZ = -5.0
CLEAR_MIN_C = 8.0
CLEAR_BEFORE_MAX_DBZ = -3.0
RANK_HIT = 0.95
RANK_MIN_C = 4.0
AGE_AFTER = (10 * 60, 90 * 60)
AGE_BEFORE = (-30 * 60, -5 * 60)


@dataclasses.dataclass
class Scan:
    t: float                  # epoch seconds (volume start)
    sweep: radar.Sweep


def decoy_offsets():
    ang = np.arange(DECOY_DIRS) * 2 * np.pi / DECOY_DIRS
    out = []
    for k, r in enumerate(DECOY_RADII_M):
        a = ang + (np.pi / DECOY_DIRS) * k        # stagger the rings
        out += [(r * np.cos(t), r * np.sin(t)) for t in a]
    return np.array(out)


def make_segments(points: pd.DataFrame, seg_len_m: float = SEG_LEN_M) -> pd.DataFrame:
    """Assign a segment id to every point: consecutive points of a flight are grouped until
    the path length reaches seg_len_m, the altitude changes by more than 1000 ft, or there
    is a gap of more than 90 s."""
    p = points.sort_values(["flight", "t"]).reset_index(drop=True)
    seg = np.zeros(len(p), np.int64)
    sid = -1
    x, y, t, alt, fl = (p.x.to_numpy(), p.y.to_numpy(), p.t.to_numpy(float),
                        p.alt_ft.to_numpy(), p.flight.to_numpy())
    length = 0.0
    alt0 = 0.0
    for i in range(len(p)):
        new = (i == 0 or fl[i] != fl[i - 1] or t[i] - t[i - 1] > 90 or length >= seg_len_m
               or abs(alt[i] - alt0) > 1000)
        if new:
            sid += 1
            length = 0.0
            alt0 = alt[i]
        else:
            length += np.hypot(x[i] - x[i - 1], y[i] - y[i - 1])
        seg[i] = sid
    p["segment"] = seg
    return p


def segment_table(points: pd.DataFrame) -> pd.DataFrame:
    g = points.groupby("segment")
    s = g.agg(flight=("flight", "first"), t=("t", "mean"), x=("x", "mean"), y=("y", "mean"),
              alt_ft=("alt_ft", "mean"), p_hpa=("p_hpa", "mean"), n=("t", "size"))
    for col in ("temp", "RH", "CLMR"):
        if col in points:
            s[col] = g[col].mean()
    return s


def _floor(z):
    return np.where(np.isnan(z), ECHO_FLOOR_DBZ, np.maximum(z, ECHO_FLOOR_DBZ))


def line_samples(points: pd.DataFrame, drift: np.ndarray, scans: list[Scan],
                 t_shift_s: float = 0.0, rotate_deg: float = 0.0,
                 age_range=(AGE_BEFORE[0], AGE_AFTER[1])) -> pd.DataFrame:
    """Per (segment, scan) with age in age_range: line and decoy reflectivity statistics.

    t_shift_s  pretend every passage happened t_shift_s later (time placebo)
    rotate_deg rotate every track point about the radar (position placebo); the drift
               vector is not rotated, so the fake tracks drift with the real wind
    """
    x0 = points.x.to_numpy()
    y0 = points.y.to_numpy()
    if rotate_deg:
        a = np.radians(rotate_deg)
        x0, y0 = x0 * np.cos(a) - y0 * np.sin(a), x0 * np.sin(a) + y0 * np.cos(a)
    t0 = points.t.to_numpy(float) + t_shift_s
    seg = points.segment.to_numpy()
    dec = decoy_offsets()
    out = []
    for sc in scans:
        age = sc.t - t0
        m = (age >= age_range[0]) & (age <= age_range[1]) & np.isfinite(drift[:, 0])
        if not m.any():
            continue
        idx = np.where(m)[0]
        bx = x0[idx] + drift[idx, 0] * age[idx]
        by = y0[idx] + drift[idx, 1] * age[idx]
        cols = {"segment": seg[idx], "age": age[idx], "on": _floor(radar.sample(sc.sweep, bx, by))}
        for k, (dx, dy) in enumerate(dec):
            cols[f"d{k}"] = _floor(radar.sample(sc.sweep, bx + dx, by + dy))
        df = pd.DataFrame(cols)
        g = df.groupby("segment").mean()
        D = g[[f"d{k}" for k in range(len(dec))]].to_numpy()
        res = pd.DataFrame({
            "segment": g.index.to_numpy(), "age": g.age.to_numpy(), "on": g.on.to_numpy(),
            "dec_med": np.median(D, axis=1),
            "rank": (D < g.on.to_numpy()[:, None]).mean(axis=1),
        })
        res["scan_t"] = sc.t
        out.append(res)
    if not out:
        return pd.DataFrame(columns=["segment", "age", "on", "dec_med", "rank", "scan_t", "c"])
    res = pd.concat(out, ignore_index=True)
    res["c"] = res.on - res.dec_med
    return res


def _longest_run(flags: pd.Series, groups: pd.Series) -> pd.Series:
    brk = (~flags).astype(int).groupby(groups).cumsum()
    runs = flags[flags].groupby([groups[flags], brk[flags]]).size()
    return runs.groupby(level=0).max()


def segment_features(samples: pd.DataFrame, segs: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per-segment features from line_samples output (see module docstring)."""
    s = samples.sort_values(["segment", "scan_t"]).reset_index(drop=True)
    bef = s[(s.age >= AGE_BEFORE[0]) & (s.age <= AGE_BEFORE[1])]
    win = s[(s.age >= AGE_AFTER[0]) & (s.age <= AGE_AFTER[1])].reset_index(drop=True)
    clear = (win.on >= CLEAR_MIN_DBZ) & (win.dec_med <= CLEAR_MAX_BG_DBZ) & (win.c >= CLEAR_MIN_C)
    rank = (win["rank"] >= RANK_HIT) & (win.c >= RANK_MIN_C)
    f = pd.DataFrame({
        "n_before": bef.groupby("segment").size(),
        "on_before_max": bef.groupby("segment").on.max(),
        "rank_before": bef.groupby("segment")["rank"].mean(),
        "c_before": bef.groupby("segment").c.mean(),
        "n_after": win.groupby("segment").size(),
        "on_after_max": win.groupby("segment").on.max(),
        "c_after_max": win.groupby("segment").c.max(),
        "persist_clear": _longest_run(clear, win.segment),
        "persist_rank": _longest_run(rank, win.segment),
    })
    for col in ("persist_clear", "persist_rank"):
        f[col] = f[col].fillna(0).astype(int)
    if segs is not None:
        f = f.join(segs, how="left")
    return f


def classify(f: pd.DataFrame, min_persist: int = 3) -> pd.DataFrame:
    """Add boolean columns: clear-air detection (quiet before passage) and rank detection
    (line not already standing out before passage)."""
    f = f.copy()
    ok = (f.n_before >= 2) & (f.n_after >= 4)
    f["hit_clear"] = ok & (f.persist_clear >= min_persist) & (f.on_before_max < CLEAR_BEFORE_MAX_DBZ)
    f["hit_rank"] = ok & (f.persist_rank >= min_persist) & (f.rank_before < 0.8)
    f["hit"] = f.hit_clear | f.hit_rank
    return f
