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
