import datetime as dt

from snow import pipeline, radar


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
