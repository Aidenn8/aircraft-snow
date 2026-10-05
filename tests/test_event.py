import numpy as np
import pandas as pd

from snow import event

SI = event.STATS.index("z")
E5 = event.STATS.index("e5")


def _synthetic(seed_echo=True, drift=(0.0, 3.0)):
    """A flight flies west-east along y = -20 km at t = 0..600 s at 10 000 ft; the wind
    drifts its air north at 3 m/s. If seed_echo, echo of 20 dBZ appears in that air 15 min
    after passage. A second flight along y = +25 km seeds nothing."""
    res = 500.0
    xs = np.arange(-70e3, 70e3 + res / 2, res)
    times = np.arange(-3600, 3 * 3600, 360.0)       # a scan every 6 minutes
    z = np.full((len(times), len(xs), len(xs)), event.FLOOR_DBZ, np.float32)
    x_track = np.linspace(-30e3, 30e3, 61)
    t_track = np.linspace(0, 600, 61)
    for k, T in enumerate(times):
        age = T - t_track
        if seed_echo:
            on = age >= 15 * 60
            for xi, a in zip(x_track[on], age[on]):
                yi = -20e3 + drift[1] * a
                ix = int(round((xi - xs[0]) / res))
                iy = int(round((yi - xs[0]) / res))
                if 0 <= iy < len(xs):
                    z[k, max(iy - 2, 0):iy + 3, max(ix - 2, 0):ix + 3] = 20.0
    grids = event.Grids(times, z, float(xs[0]), float(xs[0]), res)
    pts = pd.DataFrame({
        "flight": np.r_[np.zeros(61, int), np.ones(61, int)],
        "x": np.r_[x_track, x_track], "y": np.r_[np.full(61, -20e3), np.full(61, 25e3)],
        "t": np.r_[t_track, t_track + 120], "alt_ft": np.full(122, 10000.0),
    })
    dr = np.tile(np.array(drift), (len(pts), 1))
    return grids, pts, dr


def _acc(g, pts, dr, hour):
    perp = event.perpendicular(pts)
    return event.did_table(g, pts.x.to_numpy(), pts.y.to_numpy(), pts.t.to_numpy(float), dr, perp,
                           np.zeros(len(pts), int), hour, 1, int(hour.max()) + 1)


def test_did_positive_for_seeded_track_and_zero_for_clean_track():
    g, pts, dr = _synthetic()
    S = event.did(_acc(g, pts, dr, pts.flight.to_numpy()), min_count=10)
    assert S[0, 0, SI] > 20           # seeded: on goes -10 -> ~20 dBZ, off stays at -10
    assert S[0, 0, E5] > 0.8          # echo >= 5 dBZ appears on the line
    assert abs(S[0, 1, SI]) < 1e-9    # clean flight: nothing changes


def test_no_signal_without_seeding():
    g, pts, dr = _synthetic(seed_echo=False)
    S = event.did(_acc(g, pts, dr, np.zeros(len(pts), int)), min_count=10)
    assert abs(S[0, 0, SI]) < 1e-9


def test_wrong_drift_loses_the_signal():
    g, pts, _ = _synthetic()
    wrong = np.tile(np.array([0.0, -3.0]), (len(pts), 1))   # drift south instead of north
    S = event.did(_acc(g, pts, wrong, pts.flight.to_numpy()), min_count=10)
    assert S[0, 0, SI] < 1.0


def test_aircraft_mask_blanks_cells_near_aircraft_at_scan_time():
    from snow import radar
    az = (np.arange(720) * 0.5 + 0.25) % 360
    r = 2125.0 + 250.0 * np.arange(400)
    data = np.full((720, 400), 15.0, np.float32)
    # a 72 s clockwise sweep starting at azimuth 0: the beam points east (90 deg) at 18 s
    sw = radar.Sweep(site="TEST", time=None, elevation=0.5, azimuth=az, ray_time=np.zeros(720),
                     ranges={"REF": r}, data={"REF": data}, radar_lat=40.0, radar_lon=-100.0,
                     radar_height=300.0, duration_s=72.0, az_start=0.0)
    assert np.isclose(sw.beam_time(90.0), 18.0) and np.isclose(sw.beam_time(270.0), 54.0)
    T = 1000.0
    ac = pd.DataFrame({
        # aircraft 1 flies north along x = 20 km at 100 m/s; positions at T+3, T+13, T+23
        # (+5 s = T+8, T+18, T+28). The beam points east of the radar (azimuth ~90 deg) at
        # ~T+18, when the aircraft is at y ~ 0: masked there, not 1 km away at y = +-2 km.
        # aircraft 2: west of the radar, single report 40 s before the beam got there.
        # aircraft 3: north, an hour later.
        "hex": [1, 1, 1, 2, 3],
        "t": [T + 3.0, T + 13.0, T + 23.0, T + 54.0 - 45.0, T + 3600.0],
        "x": [20e3, 20e3, 20e3, -20e3, 0.0], "y": [-1000.0, 0.0, 1000.0, 0.0, 20e3],
        "gs": [1944, 1944, 1944, 0, 0]})          # 194.4 kt = 100 m/s
    g = event.Grids.from_sweeps([(T, sw)], half_m=40e3, aircraft=ac)
    near = g.lookup(0, np.array([20e3, 21.0e3]), np.array([0.0, 0.0]))
    assert np.all(np.isnan(near))          # within ~2 km of where the beam saw aircraft 1
    far = g.lookup(0, np.array([20e3, -20e3, 0.0, 25e3]), np.array([4500.0, 0.0, 20e3, 0.0]))
    assert np.all(far == 15.0)


def test_new_echo_statistic():
    g, pts, dr = _synthetic()
    perp = event.perpendicular(pts)
    acc = event.new_echo_table(g, pts.x.to_numpy(), pts.y.to_numpy(), pts.t.to_numpy(float), dr,
                               perp, np.zeros(len(pts), int), pts.flight.to_numpy(), 1, 2)
    S = event.new_echo_did(acc, min_count=10)
    assert S[0, 0] > 0.9          # nearly every seeded parcel shows new echo, side lines none
    assert S[0, 1] == 0.0         # the clean flight: no new echo anywhere
    g2, pts2, dr2 = _synthetic(seed_echo=False)
    acc2 = event.new_echo_table(g2, pts2.x.to_numpy(), pts2.y.to_numpy(), pts2.t.to_numpy(float),
                                dr2, perp, np.zeros(len(pts2), int), np.zeros(len(pts2), int), 1, 1)
    assert event.new_echo_did(acc2, min_count=10)[0, 0] == 0.0


def test_vectorised_rotations_match_loop():
    g, pts, dr = _synthetic()
    rng = np.random.default_rng(3)
    pts = pts.assign(x=pts.x + rng.normal(0, 3000, len(pts)), y=pts.y + rng.normal(0, 3000, len(pts)))
    perp = event.perpendicular(pts)
    x, y, t = pts.x.to_numpy(), pts.y.to_numpy(), pts.t.to_numpy(float)
    band = np.zeros(len(pts), int)
    hour = pts.flight.to_numpy()
    angles = [0.0, 37.0, 200.0]
    acc, new = event.tables_rotated(g, x, y, t, dr, perp, band, hour, 1, 2, angles)
    for r, a in enumerate(angles):
        th = np.radians(a)
        c, s = np.cos(th), np.sin(th)
        xr, yr = x * c - y * s, x * s + y * c
        pr = np.c_[perp[:, 0] * c - perp[:, 1] * s, perp[:, 0] * s + perp[:, 1] * c]
        ref = event.did_table(g, xr, yr, t, dr, pr, band, hour, 1, 2)
        ref_new = event.new_echo_table(g, xr, yr, t, dr, pr, band, hour, 1, 2)
        assert np.allclose(acc[r], ref)
        assert np.allclose(new[r], ref_new)


def test_cells_sum_matches_direct_accumulation():
    rng = np.random.default_rng(1)
    acc = rng.random((6, 5, 2, 2, len(event.STATS) + 1))
    c = event.cells(acc, 2, 3)
    assert c.shape == (5, 3, 2, 2, len(event.STATS) + 1)
    assert np.allclose(c[1, 2], acc[1:3, 2:5].sum(axis=(0, 1)))
    assert np.allclose(c[0, 0], acc[0:2, 0:3].sum(axis=(0, 1)))
