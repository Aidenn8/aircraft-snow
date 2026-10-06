"""Aircraft identity and engine class from the tar1090-db tables (see scripts/fetch_reference.py).

The ICAO type description is three characters: airframe (L land, S sea, A amphibian,
H helicopter, G gyro, T tiltrotor), engine count, and engine type (J jet, T turboprop or
turboshaft, P piston, E electric). Propeller aircraft and jets are expected to seed
supercooled cloud at different temperatures, so the engine type is the key attribute.
"""
from __future__ import annotations

import functools
import gzip
import json
import pathlib

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
DB = ROOT / "data" / "raw" / "aircraft_db"

ENGINE_CLASS = {"J": "jet", "T": "turboprop", "P": "piston", "E": "electric"}


@functools.lru_cache(maxsize=1)
def types() -> pd.DataFrame:
    """ICAO type designator -> name, description code, wake category, engine class."""
    with gzip.open(DB / "icao_aircraft_types2.js", "rt") as f:
        raw = json.load(f)
    df = pd.DataFrame.from_dict(raw, orient="index", columns=["type_name", "desc", "wtc"])
    df.index.name = "icao_type"
    desc = df["desc"].fillna("")
    df["airframe"] = desc.str[0]
    df["engines"] = pd.to_numeric(desc.str[1], errors="coerce")
    df["engine_class"] = desc.str[2].map(ENGINE_CLASS).fillna("other")
    # helicopters have turboshaft or piston engines but behave nothing like turboprops
    df.loc[df["airframe"] == "H", "engine_class"] = "helicopter"
    return df


@functools.lru_cache(maxsize=1)
def registry() -> pd.DataFrame:
    """24-bit address (int) -> registration, ICAO type, year, owner, engine class."""
    df = pd.read_csv(DB / "aircraft.csv.gz", sep=";", header=None, dtype=str,
                     names=["hex", "reg", "icao_type", "flags", "desc", "year", "owner", "x"],
                     keep_default_na=False)
    df = df[df["hex"].str.fullmatch(r"[0-9A-Fa-f]{6}")]
    df.index = df["hex"].map(lambda h: int(h, 16)).astype("int64")
    df.index.name = "hex"
    df = df[~df.index.duplicated()]
    t = types()
    df = df.join(t[["desc", "wtc", "engine_class", "engines"]].rename(columns={"desc": "type_desc"}),
                 on="icao_type")
    df["engine_class"] = df["engine_class"].fillna("unknown")
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    return df[["reg", "icao_type", "desc", "year", "owner", "type_desc", "wtc", "engines",
               "engine_class", "flags"]]


def describe(hexes) -> pd.DataFrame:
    """Registry rows for the given addresses (missing ones get engine_class 'unknown')."""
    reg = registry()
    idx = pd.Index(pd.unique(pd.Series(hexes, dtype="int64")), name="hex")
    out = reg.reindex(idx)
    out["engine_class"] = out["engine_class"].fillna("unknown")
    return out


# names the type database gets wrong or ambiguous for display
DISPLAY_NAMES = {"B77L": "Boeing 777-200LR / 777F"}
MAKERS = {"MCDONNELL-DOUGLAS": "McDonnell Douglas", "DE HAVILLAND CANADA": "De Havilland Canada"}


def display_name(icao: str) -> str:
    """Readable type name ('Boeing 737-800') for an ICAO type designator ('B738')."""
    if icao in DISPLAY_NAMES:
        return DISPLAY_NAMES[icao]
    t = types()
    if isinstance(icao, str) and icao in t.index and isinstance(t.loc[icao, "type_name"], str):
        maker, _, model = t.loc[icao, "type_name"].partition(" ")
        return f"{MAKERS.get(maker, maker.title())} {model}".strip()
    return icao if isinstance(icao, str) else ""
