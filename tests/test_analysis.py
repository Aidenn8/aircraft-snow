import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import census_analysis as ca  # noqa: E402


def _table(rows):
    """rows: (stratum, a, b, c, d) -> one record per flight."""
    out = []
    for st, a, b, c, d in rows:
        out += [{"s": st, "e": True, "y": True}] * a + [{"s": st, "e": True, "y": False}] * b
        out += [{"s": st, "e": False, "y": True}] * c + [{"s": st, "e": False, "y": False}] * d
    return pd.DataFrame(out)


def test_mantel_haenszel_common_odds_ratio():
    # both strata have odds ratio 2 (20*10/(10*10), 4*8/(8*2)), with different base rates
    df = _table([(0, 20, 10, 10, 10), (1, 4, 8, 2, 8)])
    r = ca.mantel_haenszel(df, "e", "y", ["s"])
    assert r["or"] == pytest.approx(2.0)
    assert r["lo"] < 2.0 < r["hi"]


def test_mantel_haenszel_removes_confounding():
    # no effect within either stratum (OR 1), but the exposed sit mostly in the low-rate
    # stratum, so pooling without stratifying shows a spurious effect
    df = _table([(0, 20, 10, 40, 20), (1, 2, 18, 1, 9)])
    assert ca.mantel_haenszel(df, "e", "y", ["s"])["or"] == pytest.approx(1.0)
    crude = ca.mantel_haenszel(df.assign(s=0), "e", "y", ["s"])["or"]
    assert crude == pytest.approx((22 * 29) / (28 * 41))
    assert crude < 0.6


def test_wilson_interval_contains_rate():
    lo, hi = ca.wilson(np.array([0, 5, 50]), np.array([10, 10, 100]))
    assert lo[0] == pytest.approx(0.0, abs=1e-12) and 0 < hi[0] < 0.3
    assert lo[1] < 0.5 < hi[1] and lo[2] < 0.5 < hi[2]
    assert hi[2] - lo[2] < hi[1] - lo[1]
