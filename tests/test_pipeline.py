import datetime as dt
import importlib.util
import pathlib

import numpy as np

from snow import pipeline, radar

_spec = importlib.util.spec_from_file_location(
    "census_summary", pathlib.Path(__file__).resolve().parents[1] / "scripts" / "census_summary.py")
census_summary = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(census_summary)


def test_benjamini_hochberg_known_case():
    # n = 8, q = 0.05: step-up thresholds 0.00625, 0.0125, 0.01875, ... ; only the two
    # smallest p (0.001 <= 0.00625, 0.008 <= 0.0125) pass, and no larger k does
    p = np.array([0.039, 0.001, 0.074, 0.041, 0.008, 0.205, 0.042, 0.06])
    m = census_summary.benjamini_hochberg(p, 0.05)
    assert m.tolist() == [False, True, False, False, True, False, False, False]
    # six p-values each just above its threshold: nothing rejected
    assert not census_summary.benjamini_hochberg(np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.5]), 0.05).any()


def test_drop_radar_cache_only_removes_the_day_window(tmp_path, monkeypatch):
    monkeypatch.setattr(radar, "CACHE", tmp_path)
    l3 = tmp_path / "l3"
    l3.mkdir()
    names = [
        "FTG_N0B_2025_01_25_11_00_00",   # before the window (12 UTC - 30 min = 11:30) -> keep
        "FTG_N0B_2025_01_25_11_45_00",   # inside -> delete
        "FTG_N0B_2025_01_25_23_59_00",   # inside -> delete
        "FTG_N0B_2025_01_26_12_30_00",   # inside (window ends 13:00 next day) -> delete
        "FTG_N0B_2025_01_26_13_30_00",   # after -> keep
        "LOT_N0B_2025_01_25_18_00_00",   # other radar -> keep
    ]
    for n in names:
        (l3 / n).write_bytes(b"x")
    removed = pipeline.drop_radar_cache("KFTG", dt.date(2025, 1, 25))
    left = sorted(p.name for p in l3.iterdir())
    assert removed == 3
    assert left == sorted([names[0], names[4], names[5]])
