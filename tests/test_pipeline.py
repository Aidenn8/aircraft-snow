import datetime as dt
import importlib.util
import pathlib

import numpy as np

from snow import pipeline, radar

_spec = importlib.util.spec_from_file_location(
    "census_summary", pathlib.Path(__file__).resolve().parents[1] / "scripts" / "census_summary.py")
census_summary = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(census_summary)


def test_empirical_null_recovers_a_shifted_wide_null():
    rng = np.random.default_rng(3)
    z = np.r_[rng.normal(0.5, 1.6, 600), np.full(12, 12.0)]       # 600 null days, 12 strong events
    fit = census_summary.fit_null(z)
    assert abs(fit["centre"] - 0.5) < 0.25 and abs(fit["width"] - 1.6) < 0.25
    p = census_summary.tail_p(z, fit)
    sel = census_summary.benjamini_hochberg(p, 0.05)
    assert sel[-12:].all()                                            # the events are found
    assert sel[:600].sum() <= 3                                       # few null days are
    # the theoretical N(0, 1) null would have picked many null days
    assert census_summary.benjamini_hochberg(
        np.where(np.isfinite(z), census_summary.stats.norm.sf(z), 1.0), 0.05)[:600].sum() > 20
    assert census_summary.tail_p(np.array([np.nan]), fit)[0] == 1.0


def test_event_share_counts_the_excess_of_small_p_values():
    rng = np.random.default_rng(4)
    p = np.r_[rng.uniform(0, 1, 900), rng.uniform(0, 0.002, 100)]     # 10% event days
    es = census_summary.event_share(p, n_boot=200)
    # the excess at p <= 0.01 recovers the strong events sharply (100 events, ~9 null days
    # below 0.01 against 10 expected)
    assert 90 <= es["excess_p01"] <= 110 and es["excess_p01_lo"] <= 100 <= es["excess_p01_hi"]
    # Storey's estimate of all event days is unbiased but noisy at this size (sd ~30 days)
    assert es["n_events_storey_lo"] <= 100 <= es["n_events_storey_hi"]
    null = census_summary.event_share(rng.uniform(0, 1, 1000), n_boot=200)
    assert abs(null["excess_p01"]) < 10
    q = census_summary.q_values(p, es["pi0"])
    assert (q[-100:] < 0.05).mean() > 0.9 and (q[:900] < 0.05).mean() < 0.01


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


def test_screen_passes():
    assert pipeline.screen_passes({"p_day": 0.04, "p_day_gamma": 0.5})
    assert pipeline.screen_passes({"p_day": 0.2, "p_day_gamma": 0.08})
    assert not pipeline.screen_passes({"p_day": 0.08, "p_day_gamma": 0.3})
    assert not pipeline.screen_passes({"p_day": 0.08, "p_day_gamma": float("nan")})
    assert not pipeline.screen_passes({})
