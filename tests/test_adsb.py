import struct

import numpy as np

from snow import adsb


def slot(t_ms, interval_ms=10_000):
    # readsb stores the low 32 bits in a signed int32 field; pack them unsigned (same bytes)
    return struct.pack("<iiIhh", adsb.SLOT_MARKER, t_ms >> 32, t_ms & 0xFFFFFFFF,
                       interval_ms, 0)


def pos(hex_, lat, lon, alt_ft=None, gs_kt=250.0, addrtype=0, ground=False):
    alt = adsb.ALT_GROUND if ground else (adsb.ALT_UNKNOWN if alt_ft is None else round(alt_ft / 25))
    word0 = hex_ | (addrtype << 27)
    lat_i = round(lat * 1e6)
    lon_i = round(lon * 1e6)
    # struct wants signed int32 for word0; addrtype bits may set the sign bit
    word0 = word0 - (1 << 32) if word0 >= (1 << 31) else word0
    return struct.pack("<iiihh", word0, lat_i, lon_i, alt, round(gs_kt * 10))


def callsign(hex_, cs, squawk):
    return struct.pack("<ii", hex_, (1 << 30) + squawk) + cs.ljust(8).encode()


def test_parse_heatmap_roundtrip():
    t0 = 1_737_853_200_000  # 2025-01-26 01:00:00 UTC in ms
    buf = b"".join([
        slot(t0),
        pos(0xA1B2C3, 39.85, -104.67, alt_ft=9000, gs_kt=180.5),       # Denver approach
        pos(0xA1B2C4, -33.9, 151.2, alt_ft=5000),                      # Sydney: outside boxes
        callsign(0xA1B2C3, "UAL5528", 7700),
        slot(t0 + 10_000),
        pos(0xA1B2C3, 39.86, -104.67, ground=True, gs_kt=20),          # on the ground
        pos(0xA1B2C5, 41.97, -87.90, alt_ft=None),                     # unknown altitude
        pos(0xA1B2C6, 41.97, -87.90, alt_ft=38000),                    # above MAX_ALT_FT
        pos(0xA1B2C7 | (1 << 24), 61.17, -150.0, alt_ft=3000, addrtype=5),  # Alaska, MLAT, non-ICAO
    ])
    p, c = adsb.parse_heatmap(buf)
    assert list(p["hex"]) == [0xA1B2C3, 0xA1B2C3, 0xA1B2C7 | (1 << 24)]
    assert list(p["t"]) == [t0 // 1000, t0 // 1000 + 10, t0 // 1000 + 10]
    assert np.isclose(p["lat"].iloc[0], 39.85, atol=1e-5)
    assert np.isclose(p["lon"].iloc[0], -104.67, atol=1e-5)
    assert p["alt"].iloc[0] * 25 == 9000
    assert p["gs"].iloc[0] == 1805
    assert p["alt"].iloc[1] == adsb.ALT_GROUND
    assert p["atype"].iloc[2] == 5
    # the callsign record is kept because its aircraft has a kept position, and the
    # Sydney position must not be misread as a callsign record (negative latitude)
    assert list(c["callsign"]) == ["UAL5528"]
    assert list(c["squawk"]) == [7700]
    assert list(c["hex"]) == [0xA1B2C3]


def test_negative_latitude_is_not_callsign():
    buf = slot(1_700_000_000_000) + pos(0x123456, -45.0, -70.0, alt_ft=3000)
    p, c = adsb.parse_heatmap(buf, boxes=((-50.0, -40.0, -75.0, -65.0),))
    assert len(p) == 1 and len(c) == 0
