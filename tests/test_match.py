import numpy as np
import pandas as pd

from snow import match


def _table(real, plac):
    return pd.DataFrame({"coverage": real, "placebo_cov": list(plac)})


def test_day_statistic_null_is_calibrated():
    """If the real coverages are just another draw from the placebo distribution, the
    day p-value should be roughly uniform: about 5% of days below 0.05."""
    rng = np.random.default_rng(0)
    ps = []
    for _ in range(300):
        n, R = 20, 60
        allc = np.where(rng.random((n, R + 1)) < 0.4, rng.random((n, R + 1)), 0.0)
        ps.append(match.day_statistic(_table(allc[:, 0], allc[:, 1:]))["p_day"])
    ps = np.array(ps)
    assert 0.01 < np.mean(ps <= 0.05) < 0.11
    assert 0.4 < np.median(ps) < 0.6


def test_day_statistic_detects_signal():
    rng = np.random.default_rng(1)
    n, R = 20, 60
    plac = np.where(rng.random((n, R)) < 0.4, rng.random((n, R)) * 0.8, 0.0)
    real = np.full(n, 0.95)
    s = match.day_statistic(_table(real, plac))
    assert s["p_day"] <= 2 / (R + 1)
    assert s["n_matched"] == n


def test_empirical_fdr_null_and_signal():
    rng = np.random.default_rng(2)
    n, R = 200, 50
    allc = np.where(rng.random((n, R + 1)) < 0.5, rng.random((n, R + 1)), 0.0)
    f = match.empirical_fdr(_table(allc[:, 0], allc[:, 1:]), thresholds=(0.05,))
    assert 0.5 < f.fdr.iloc[0] < 2.0                      # all matches are expected false
    real = allc[:, 0].copy()
    real[:40] = 1.0                                        # 40 genuine streaks beat every placebo
    f = match.empirical_fdr(_table(real, allc[:, 1:]), thresholds=(0.05,))
    assert f.n_real.iloc[0] >= 40 and f.fdr.iloc[0] < 0.4 and f.excess.iloc[0] > 30


def test_newness_quiet_and_echo():
    from snow import event
    res = 500.0
    xs = np.arange(-20e3, 20e3 + res / 2, res)
    times = np.array([0.0, 300.0, 600.0, 900.0, 1200.0])
    z = np.full((5, len(xs), len(xs)), -10.0, np.float32)
    g = event.Grids(times, z, float(xs[0]), float(xs[0]), res)
    sx = np.linspace(-5e3, 5e3, 21)
    sy = np.zeros(21)
    q, n = match.newness(g, sx, sy, T=1200.0, t_pass=1000.0, u=0.0, v=0.0)
    assert n >= 2 and q == 1.0                   # nothing there before passage
    g.z[:] = 10.0
    q, n = match.newness(g, sx, sy, T=1200.0, t_pass=1000.0, u=0.0, v=0.0)
    assert q == 0.0                              # echo already present
