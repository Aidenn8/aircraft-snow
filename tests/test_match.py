import numpy as np
import pandas as pd

from snow import match


def _table(real, plac):
    """Synthetic streak table: every streak matched by a different flight in every set
    (so per-flight aggregation reduces to per-streak)."""
    n, R = plac.shape
    flight = np.where(real > 0, np.arange(n), -1)
    pf = np.where(plac > 0, 1000 * (1 + np.arange(R))[None, :] + np.arange(n)[:, None], -1)
    return pd.DataFrame({"coverage": real, "flight": flight, "placebo_cov": list(plac),
                         "placebo_flight": list(pf)})


def test_day_statistic_counts_each_flight_once():
    """The same flight matched to one echo in five scans must count once, like one match."""
    rng = np.random.default_rng(5)
    R = 40
    plac = np.where(rng.random((6, R)) < 0.3, rng.random((6, R)) * 0.6, 0.0)
    real = np.r_[np.full(5, 0.9), 0.0]
    t = _table(real, plac)
    t["flight"] = [7, 7, 7, 7, 7, -1]                    # one flight, five scans
    single = _table(real[:1], plac[:1])
    single["flight"] = [7]
    assert len(match.per_flight_min_p(match.pseudo_pvalues(t), match.flight_matrix(t), 0)) == 1
    assert np.isclose(match.day_statistic(t)["E"], match.day_statistic(single)["E"])


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


def _exact_coverage(sx, sy, px, py, tol):
    d = np.hypot(sx[:, None] - px[None, :], sy[:, None] - py[None, :])
    return float((d.min(axis=1) <= tol).mean())


def test_candidates_coverage_matches_exact_computation():
    sx = np.arange(0, 21000, 500.0)              # straight skeleton, 500 m pixels
    sy = np.zeros_like(sx)
    xs_trk = np.arange(-2000, 23000, 250.0)      # interpolated track, 250 m spacing
    tracks = {
        1: (xs_trk, np.full_like(xs_trk, 600.0)),             # parallel, 0.6 km off: matches
        # parallel, 5 km off: beyond TOL plus the largest drift correction (1 m/s x 30 min)
        2: (xs_trk, np.full_like(xs_trk, 5000.0)),
        3: (np.linspace(5000, 15000, 81), np.linspace(-8000, 8000, 81)),   # crossing
    }
    x = np.concatenate([v[0] for v in tracks.values()])
    y = np.concatenate([v[1] for v in tracks.values()])
    flight = np.concatenate([np.full(len(v[0]), k) for k, v in tracks.items()])
    t = np.zeros(len(x))
    drift = np.zeros((len(x), 2))
    T = 30 * 60.0
    cands = {c[1]: c for c in match._candidates(sx, sy, x, y, t, T, drift, flight)}
    assert 1 in cands and 2 not in cands and 3 not in cands
    f1, f, du, dv, t_pass, u, v, cov, prec = cands[1]
    exact = _exact_coverage(sx, sy, tracks[1][0], tracks[1][1], match.TOL_M)
    assert abs(cov - exact) < 0.05 and cov > 0.95
    assert prec > 0.8


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


def test_shift_placebos_preserve_orientation_and_exceed_reach():
    shifts = match.shift_placebos()
    assert len(shifts) == 24 and all(isinstance(s, tuple) for s in shifts)
    reach = match.TOL_M + max(match.DRIFT_PERTURB_MS) * match.AGE_MAX_S * np.sqrt(2)
    assert min(np.hypot(*s) for s in shifts) > reach
    x = np.array([0.0, 1000.0])
    y = np.array([0.0, 1000.0])
    sx, sy = match.placebo_transform(x, y, shifts[1])
    assert np.allclose(np.diff(sx), np.diff(x)) and np.allclose(np.diff(sy), np.diff(y))
    rx, ry = match.placebo_transform(x, y, 90.0)
    assert np.allclose([rx[1], ry[1]], [-1000.0, 1000.0])


def test_fall_ok_bounds():
    age = 30 * 60.0
    bw = 200.0
    assert match.fall_ok(1000.0, bw, age)                    # 0.56 m/s: snow
    assert not match.fall_ok(0.0, 0.0, age)                  # at the beam: falls out at once
    assert not match.fall_ok(-500.0, bw, age)                # below the beam
    assert not match.fall_ok(5000.0, bw, age)                # 2.8 m/s: faster than snow
    assert match.fall_ok(1000.0, bw, -age)                   # reversed diagnostic uses |age|
    # the beam width widens the window
    lo = match.V_FALL_MIN_MS * age
    assert match.fall_ok(lo - bw + 1, bw, age) and not match.fall_ok(lo - bw - 1, bw, age)


def _one_streak_case(dh, temp_c):
    """A 15 km streak at y = 0 in a quiet grid, and one flight along it ~25-30 min earlier."""
    from snow import event
    res = 500.0
    xs = np.arange(-30e3, 30e3 + res / 2, res)
    times = np.arange(0.0, 3601.0, 300.0)
    z = np.full((len(times), len(xs), len(xs)), -10.0, np.float32)
    lab = np.zeros((len(xs), len(xs)), np.int32)
    iy = int(np.argmin(np.abs(xs)))
    ix = (xs >= 5e3) & (xs <= 20e3)
    lab[iy, ix] = 1
    z[-1][lab == 1] = 15.0
    g = event.Grids(times, z, float(xs[0]), float(xs[0]), res)
    objects = pd.DataFrame({"label": [1]})
    tx = np.arange(0.0, 25e3, 250.0)
    pts = pd.DataFrame({"flight": 7, "x": tx, "y": 0.0, "t": 1800.0 + tx / 200.0,
                        "alt_ft": 10000.0, "p_hpa": 700.0, "dh": dh, "bw": 200.0, "temp_c": temp_c})
    drift = np.zeros((len(pts), 2))
    return g, lab, objects, xs, pts, drift


def test_match_streaks_physical_requirements():
    T = 3600.0
    angles = (90.0, 180.0, 270.0)
    for dh, temp, tmax, expect in ((1000.0, -15.0, match.T_SEED_MAX_C, True),   # snow-like fall
                                   (0.0, -15.0, match.T_SEED_MAX_C, False),     # no fall: impossible
                                   (1000.0, 2.0, match.T_SEED_MAX_C, False),    # above freezing
                                   (1000.0, 2.0, None, True)):                  # null mode: no T test
        g, lab, objects, xs, pts, drift = _one_streak_case(dh, temp)
        r = match.match_streaks(g, T, lab, objects, xs, xs, pts, drift, angles=angles, max_temp_c=tmax)
        assert (r.flight.iloc[0] == 7) == expect, (dh, temp, tmax)
        if expect:
            assert abs(r.v_fall.iloc[0] - dh / (T - r.t_pass.iloc[0])) < 0.2
