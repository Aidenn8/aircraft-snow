"""Radar reader checks. The Level II/III comparison needs the cached sample files and is
skipped when they are absent (run scripts that fetch the Denver case first)."""
import pathlib

import numpy as np
import pytest

from snow import radar

ROOT = pathlib.Path(__file__).resolve().parents[1]
L2 = ROOT / "data/radar/l2/2025/01/26/KFTG/KFTG20250126_013611_V06"
L3 = ROOT / "data/radar/l3/FTG_N0B_2025_01_26_01_36_11"


def test_beam_height_flat_beam_rises_with_range():
    # at 0 deg elevation the 4/3-earth beam is ~ r^2 / (2 k a) above the radar
    h, s = radar.beam_height_ground_range(np.array([50_000.0, 100_000.0]), 0.0, 0.0)
    ka = radar.K_EFF * radar.EARTH_RADIUS
    assert np.allclose(h, np.array([50e3, 100e3]) ** 2 / (2 * ka), rtol=1e-3)
    assert np.all(s < np.array([50e3, 100e3])) and np.all(s > np.array([49.9e3, 99.9e3]))


def test_beam_height_half_degree_at_40km():
    h, _ = radar.beam_height_ground_range(np.array([40_000.0]), 0.5, 1709.0)
    # 40 km * tan(0.5 deg) = 349 m, plus curvature 40e3^2 / (2 * 8.49e6) = 94 m
    assert abs(h[0] - (1709 + 349 + 94)) < 5


@pytest.mark.skipif(not (L2.exists() and L3.exists()), reason="sample radar files not cached")
def test_level2_and_level3_lowest_sweep_agree():
    sweeps = radar.read_level2(L2, moments=("REF",))
    s2 = sweeps[0]
    s3 = radar.read_level3(L3)
    assert abs(s2.elevation - 0.48) < 0.05 and abs(s3.elevation - 0.5) < 0.05
    assert abs(s2.radar_lat - s3.radar_lat) < 2e-3 and abs(s2.radar_lon - s3.radar_lon) < 2e-3
    # match rays by azimuth and gates by range
    r2, r3 = s2.ranges["REF"], s3.ranges["REF"]
    j3 = np.searchsorted(r3, r2[0] - 1)
    assert abs(r3[j3] - r2[0]) < 1.0
    n = min(len(r2), len(r3) - j3)
    o2, o3 = np.argsort(s2.azimuth), np.argsort(s3.azimuth)
    a2, a3 = s2.azimuth[o2], s3.azimuth[o3]
    assert np.max(np.abs(((a2 - a3) + 180) % 360 - 180)) < 0.3
    v2 = s2.data["REF"][o2, :n]
    v3 = s3.data["REF"][o3, j3:j3 + n]
    both = np.isfinite(v2) & np.isfinite(v3)
    assert both.sum() > 10_000
    assert np.median(np.abs(v2[both] - v3[both])) <= 0.5
    agree_mask = np.mean(np.isfinite(v2) == np.isfinite(v3))
    assert agree_mask > 0.98
