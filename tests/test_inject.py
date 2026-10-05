import numpy as np

from snow import event, inject


def test_envelope_shape():
    a = np.array([0, inject.ONSET_S, inject.PEAK_S, (inject.PEAK_S + inject.END_S) / 2, inject.END_S + 1])
    e = inject.envelope(a)
    assert e[0] == 0 and e[1] == 0 and e[2] == 1 and np.isclose(e[3], 0.5) and e[4] == 0


def test_injected_streak_peak_and_position():
    res = 500.0
    xs = np.arange(-30e3, 30e3 + res / 2, res)
    times = np.array([0.0, inject.PEAK_S])
    z = np.full((2, len(xs), len(xs)), event.FLOOR_DBZ, np.float32)
    g = event.Grids(times, z, float(xs[0]), float(xs[0]), res)
    x = np.linspace(-10e3, 10e3, 41)
    y = np.zeros(41)
    t = np.zeros(41)
    drift = np.tile([0.0, 2.0], (41, 1))                 # north at 2 m/s
    g2, err = inject.inject(g, x, y, t, drift, peak_dbz=20.0, rng=np.random.default_rng(0), wind_err_ms=0.0)
    assert np.all(g2.z[0] == event.FLOOR_DBZ)            # nothing at age 0
    y_expect = 2.0 * inject.PEAK_S                        # 4.2 km north at peak
    col = g2.z[1][:, len(xs) // 2]
    iy = np.argmax(col)
    assert abs(xs[iy] - y_expect) <= res
    assert abs(col.max() - 20.0) < 0.5                    # peak reflectivity as requested


def test_sigma_for_width_gives_the_threshold_width():
    for peak in (8.0, 15.0, 25.0):
        s = inject.sigma_for_width(2000.0, peak)
        d = np.linspace(-5000, 5000, 20001)
        z = 10 * np.log10(10 ** (peak / 10) * np.exp(-d ** 2 / (2 * s ** 2)))
        above = d[z >= 3.0]
        assert abs((above.max() - above.min()) - 2000.0) < 5.0


def test_envelope_fall_timing():
    arrival = np.array([1800.0, 1800.0, np.nan, -100.0, 300.0])
    age = np.array([1800.0, 1500.0, 1800.0, 1800.0, 400.0])
    w = inject.envelope_fall(age, arrival)
    assert w[0] == 1.0                                     # peak when the ice reaches the beam
    # half-width max(0.3 x 1800 s, 600 s) = 600 s, so 300 s early is half way up
    assert abs(w[1] - 0.5) < 1e-9
    assert w[2] == 0 and w[3] == 0                         # never reaches the beam
    assert w[4] == 0                                       # too young to be visible
