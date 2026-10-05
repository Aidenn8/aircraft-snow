"""Cheap pre-screens that skip radar-days on which aircraft cannot have made visible snow.

Both are deliberately permissive: a missed event costs more than an extra day processed.

cold_moist(...)   any hourly HRRR analysis over the radar with a level between ALT_LO_FT
                  and ALT_HI_FT (pressure altitude) that is at most MAX_TEMP_C and either
                  has cloud water or relative humidity >= MIN_RH
echo_present(...) any lowest-tilt scan with at least MIN_ECHO_FRAC of the analysis area at
                  >= MIN_ECHO_DBZ (a streak is an echo; a day without echo has none)
"""
from __future__ import annotations

import datetime as dt

import numpy as np
import pandas as pd

from . import hrrr
from .weather import pressure_from_pressure_altitude

ALT_LO_FT = 2000.0
ALT_HI_FT = 24000.0
MAX_TEMP_C = -4.0
MIN_RH = 75.0
MIN_ECHO_DBZ = 5.0
MIN_ECHO_FRAC = 0.005


def cold_moist(lat: float, lon: float, t0: dt.datetime, hours: int = 24) -> dict:
    p_hi = pressure_from_pressure_altitude(ALT_LO_FT)
    p_lo = pressure_from_pressure_altitude(ALT_HI_FT)
    hits, n = [], 0
    for k in range(hours + 1):
        run = (t0 + dt.timedelta(hours=k)).replace(tzinfo=None)
        prof = hrrr.profile(run, lat, lon, box=2)
        n += 1
        lay = prof[(prof.p <= p_hi) & (prof.p >= p_lo)]
        ok = (lay.temp <= MAX_TEMP_C) & ((lay.CLMR > 1e-3) | (lay.RH >= MIN_RH))
        if ok.any():
            hits.append(k)
    return {"pass": bool(hits), "hours": hits, "n_profiles": n}


def echo_present(grids) -> dict:
    valid = np.isfinite(grids.z)
    frac = np.where(valid, grids.z >= MIN_ECHO_DBZ, False).sum(axis=(1, 2)) / np.maximum(valid.sum(axis=(1, 2)), 1)
    return {"pass": bool((frac >= MIN_ECHO_FRAC).any()), "max_frac": float(frac.max()) if len(frac) else 0.0,
            "n_scans_echo": int((frac >= MIN_ECHO_FRAC).sum())}
