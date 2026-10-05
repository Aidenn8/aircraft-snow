"""Stream one day of adsb.lol history and keep low-level aircraft positions over North America.

adsb.lol publishes each day as a tar archive (split into ~2 GB parts) on GitHub, under the
Open Database License. The archive starts with heatmap/NN.bin.ttf: 48 gzip files covering
30 minutes each, followed by gigabytes of per-aircraft trace files. The heatmap files hold
one position per aircraft per 10 s slot, which is all this project needs, so the archive is
streamed and abandoned once the heatmap section is over (about a quarter of the bytes).

Record layout, checked against readsb's globe_index.c (struct heatEntry, 16 bytes, little
endian):
  int32 hex   0x0e7f7c9d marks a slot record. Otherwise bits 0-23 are the 24-bit address,
              bit 24 the readsb "non-ICAO address" flag, bits 27-31 the address type.
  int32 lat   degrees * 1e6. A callsign record has lat = 2**30 + squawk, the squawk written
              as a decimal number (e.g. 7700).
  int32 lon   degrees * 1e6. In a callsign record, bytes 8-15 hold the 8-char callsign.
  int16 alt   altitude / 25 ft: barometric (pressure) altitude when the aircraft reports it,
              otherwise GNSS altitude; -123 = on the ground, -124 = unknown.
  int16 gs    ground speed * 10 kt; -1 = unknown.
A slot record stores the slot start time in ms as (lat << 32) | lon and the slot length in
alt (10 000 ms for adsb.lol). Each aircraft has at most one position per slot: the first
point of its trace at or after the slot start, so positions are timed to within 10 s.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import io
import pathlib
import subprocess
import sys
import tarfile
import time

import numpy as np
import pandas as pd
import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "adsb"
SLOT_MARKER = 0x0E7F7C9D
ALT_GROUND = -123
ALT_UNKNOWN = -124

# Regions kept: contiguous US with southern Canada, and Alaska. (lat0, lat1, lon0, lon1)
BOXES = ((24.0, 50.5, -125.5, -66.0), (55.0, 72.0, -170.0, -129.0))
MAX_ALT_FT = 30_000  # supercooled cloud worth studying lies well below this

INSTANCES = ("prod-0", "prod-0tmp", "staging-0", "staging-0tmp", "prod-1")

_session = requests.Session()


def _token() -> str | None:
    try:
        return subprocess.run(["gh", "auth", "token"], capture_output=True, text=True,
                              check=True).stdout.strip() or None
    except Exception:
        return None


_TOKEN = _token()


def release_tags(year: int, refresh=False) -> set[str]:
    """All release tags of one yearly repository (cached in data/raw/releases_YYYY.txt)."""
    path = ROOT / "data" / "raw" / f"releases_{year}.txt"
    if refresh or not path.exists():
        out = subprocess.run(["gh", "api", "--paginate",
                              f"repos/adsblol/globe_history_{year}/releases?per_page=100",
                              "--jq", ".[].tag_name"], capture_output=True, text=True, check=True)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(out.stdout)
    return set(path.read_text().split())


def release_parts(day: dt.date) -> tuple[str, list[str]]:
    """Return (tag, [part urls in order]) for the preferred instance that published this day.
    A day is normally in its own year's repository, but releases made just after New Year
    (31 December, for one) can sit in the next year's."""
    headers = {"Authorization": f"Bearer {_TOKEN}"} if _TOKEN else {}
    prefix = f"v{day:%Y.%m.%d}-"
    years = [day.year] + ([day.year + 1] if day.month == 12 else [])
    found = {}
    for y in years:
        tags = release_tags(y)
        if not any(t.startswith(prefix) for t in tags):
            tags = release_tags(y, refresh=True)
        for t in tags:
            if t.startswith(prefix):
                found.setdefault(t, f"adsblol/globe_history_{y}")
    for inst in INSTANCES:
        tag = f"v{day:%Y.%m.%d}-planes-readsb-{inst}"
        if tag not in found:
            continue
        repo = found[tag]
        r = _session.get(f"https://api.github.com/repos/{repo}/releases/tags/{tag}",
                         headers=headers, timeout=60)
        if r.status_code != 200:
            continue
        assets = sorted(r.json()["assets"], key=lambda a: a["name"])
        urls = [a["browser_download_url"] for a in assets
                if ".tar" in a["name"] and not a["name"].endswith(".sha256")]
        if urls:
            return tag, urls
    raise FileNotFoundError(f"no release for {day}")


class _Concat(io.RawIOBase):
    """Read several HTTP responses back to back as one stream."""

    def __init__(self, urls):
        self.urls = list(urls)
        self.resp = None
        self.nbytes = 0

    def readable(self):
        return True

    def _next(self):
        if self.resp is not None:
            self.resp.close()
        if not self.urls:
            self.resp = None
            return False
        self.resp = _session.get(self.urls.pop(0), stream=True, timeout=120)
        self.resp.raise_for_status()
        self.resp.raw.decode_content = False
        return True

    def readinto(self, b):
        while True:
            if self.resp is None and not self._next():
                return 0
            n = self.resp.raw.readinto(b)
            if n:
                self.nbytes += n
                return n
            if not self._next():
                return 0

    def close(self):
        if self.resp is not None:
            self.resp.close()
        super().close()


def in_boxes(lat, lon, boxes=BOXES):
    keep = np.zeros(lat.shape, bool)
    for lat0, lat1, lon0, lon1 in boxes:
        keep |= (lat >= lat0) & (lat <= lat1) & (lon >= lon0) & (lon <= lon1)
    return keep


def parse_heatmap(buf: bytes, boxes=BOXES, max_alt_ft=MAX_ALT_FT):
    """Decode one heatmap file; return (positions, callsign records) as DataFrames."""
    a = np.frombuffer(buf, dtype="<i4")
    a = a[: len(a) // 4 * 4].reshape(-1, 4)
    is_slot = a[:, 0] == SLOT_MARKER
    rows = np.arange(len(a))
    last_slot = np.maximum.accumulate(np.where(is_slot, rows, -1))
    t_ms = (a[:, 1].astype(np.int64) << 32) | (a[:, 2].astype(np.int64) & 0xFFFFFFFF)
    t = np.where(last_slot >= 0, t_ms[np.maximum(last_slot, 0)] // 1000, -1)
    # Callsign records have lat = 2**30 + squawk (< 2**30 + 7778). Latitudes are within
    # +-9e7, so a plain bit-30 test would misfire on negative (southern) latitudes.
    is_call = (a[:, 1] >= (1 << 30)) & ~is_slot
    is_pos = ~is_slot & ~is_call & (last_slot >= 0)

    alt = (a[:, 3] & 0xFFFF).astype(np.uint16).view(np.int16)
    lat = a[:, 1] / 1e6
    lon = a[:, 2] / 1e6
    # compare in file units: alt * 25 would overflow int16 above 32 767 ft
    alt_ok = (alt == ALT_GROUND) | ((alt != ALT_UNKNOWN) & (alt <= max_alt_ft // 25))
    keep = is_pos & alt_ok & in_boxes(lat, lon, boxes)
    w = a[keep]
    pos = pd.DataFrame({
        "t": t[keep].astype(np.int64),
        "hex": (w[:, 0] & 0x1FFFFFF).astype(np.int32),
        "atype": ((w[:, 0] >> 27) & 0x1F).astype(np.int8),
        "lat": lat[keep].astype(np.float32),
        "lon": lon[keep].astype(np.float32),
        "alt": alt[keep],
        "gs": ((w[:, 3] >> 16) & 0xFFFF).astype(np.uint16).view(np.int16),
    })

    hexes = np.unique(pos["hex"].to_numpy())
    c = a[is_call]
    ch = c[:, 0] & 0x1FFFFFF
    m = np.isin(ch, hexes)
    raw = np.ascontiguousarray(c[m]).view(np.uint8).reshape(-1, 16)
    calls = pd.DataFrame({
        "t": t[is_call][m].astype(np.int64),
        "hex": ch[m].astype(np.int32),
        "callsign": [bytes(r[8:16]).decode("ascii", "replace").strip() for r in raw],
        "squawk": (c[m][:, 1] & 0x3FFFFFFF).astype(np.int32),
    })
    return pos, calls


def day_path(day: dt.date) -> pathlib.Path:
    return OUT / f"{day:%Y}" / f"{day:%Y-%m-%d}.parquet"


def fetch_day(day: dt.date, force=False) -> pathlib.Path:
    path = day_path(day)
    cpath = path.with_suffix(".calls.parquet")
    if path.exists() and cpath.exists() and not force:
        return path
    tag, urls = release_parts(day)
    t0 = time.time()
    stream = _Concat(urls)
    pos_parts, call_parts, names = [], [], []
    with tarfile.open(fileobj=io.BufferedReader(stream, 1 << 20), mode="r|") as tar:
        for member in tar:
            name = member.name.lstrip("./")
            if name.startswith("heatmap/") and member.isfile():
                buf = gzip.decompress(tar.extractfile(member).read())
                p, c = parse_heatmap(buf)
                pos_parts.append(p)
                call_parts.append(c)
                names.append(name)
            elif names and not name.startswith("heatmap"):
                break  # trace files start; the rest of the archive is not needed
    stream.close()
    if not names:
        raise RuntimeError(f"{tag}: no heatmap files")
    pos = (pd.concat(pos_parts, ignore_index=True)
           .sort_values(["hex", "t"], kind="stable").reset_index(drop=True))
    calls = (pd.concat(call_parts, ignore_index=True)
             .groupby(["hex", "callsign", "squawk"], as_index=False)
             .agg(t0=("t", "min"), t1=("t", "max"), n=("t", "size")))
    path.parent.mkdir(parents=True, exist_ok=True)
    pos.to_parquet(path, index=False, compression="zstd")
    calls.to_parquet(cpath, index=False, compression="zstd")
    meta = {"tag": tag, "heatmap_files": len(names), "positions": len(pos),
            "aircraft": int(pos["hex"].nunique()), "bytes_read": stream.nbytes,
            "seconds": round(time.time() - t0, 1)}
    path.with_suffix(".meta.json").write_text(pd.Series(meta).to_json())
    print(f"{day} {tag}: {len(names)} heatmap files, {len(pos):,} positions, "
          f"{meta['aircraft']:,} aircraft, {stream.nbytes / 1e6:.0f} MB in {meta['seconds']:.0f} s",
          flush=True)
    return path


def load_day(day: dt.date, columns=None) -> pd.DataFrame:
    return pd.read_parquet(day_path(day), columns=columns)


def load_calls(day: dt.date) -> pd.DataFrame:
    return pd.read_parquet(day_path(day).with_suffix(".calls.parquet"))


def parse_dates(specs) -> list[dt.date]:
    days = []
    for d in specs:
        if ":" in d:
            a, b = (dt.date.fromisoformat(x) for x in d.split(":"))
            days += [a + dt.timedelta(n) for n in range((b - a).days + 1)]
        else:
            days.append(dt.date.fromisoformat(d))
    return days


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("dates", nargs="+", help="YYYY-MM-DD or YYYY-MM-DD:YYYY-MM-DD ranges")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    for day in parse_dates(args.dates):
        for attempt in range(3):
            try:
                fetch_day(day, force=args.force)
                break
            except FileNotFoundError as e:
                print(f"{day}: {e}", file=sys.stderr, flush=True)
                break
            except Exception as e:  # network hiccups: retry, then move on
                print(f"{day} attempt {attempt + 1} failed: {e!r}", file=sys.stderr, flush=True)
                time.sleep(10 * (attempt + 1))


if __name__ == "__main__":
    main()
