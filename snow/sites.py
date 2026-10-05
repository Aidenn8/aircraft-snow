"""NEXRAD station list (NOAA NCEI nexrad-stations.txt, fixed width)."""
from __future__ import annotations

import functools
import pathlib

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
STATIONS = ROOT / "data" / "raw" / "sites" / "nexrad-stations.txt"


@functools.lru_cache(maxsize=1)
def nexrad() -> pd.DataFrame:
    """All NEXRAD sites; column extents are read from the dashed underline of the header."""
    lines = STATIONS.read_text().splitlines()
    hdr, dash = lines[0], lines[1]
    spans, i = [], 0
    while i < len(dash):
        if dash[i] == "-":
            j = i
            while j < len(dash) and dash[j] == "-":
                j += 1
            spans.append((i, j))
            i = j
        else:
            i += 1
    names = [hdr[a:b].strip() for a, b in spans]
    rows = [[ln[a:b].strip() for a, b in spans] for ln in lines[2:] if ln.strip()]
    df = pd.DataFrame(rows, columns=names)
    df = df[df.STNTYPE.str.contains("NEXRAD")].copy()
    df["LAT"] = df.LAT.astype(float)
    df["LON"] = df.LON.astype(float)
    df["ELEV"] = pd.to_numeric(df.ELEV, errors="coerce")
    return df[["ICAO", "NAME", "ST", "LAT", "LON", "ELEV"]].rename(columns=str.lower).set_index("icao")


def location(site: str) -> tuple[float, float]:
    r = nexrad().loc[site]
    return float(r.lat), float(r.lon)
