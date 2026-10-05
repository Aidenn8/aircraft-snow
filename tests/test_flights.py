import numpy as np
import pandas as pd

from snow import flights, weather


def test_pressure_altitude_roundtrip_and_known_values():
    # ISA: 10 000 ft -> 696.8 hPa, 5 000 ft -> 843.1 hPa, 0 ft -> 1013.25 hPa
    p = weather.pressure_from_pressure_altitude(np.array([0.0, 5000.0, 10000.0]))
    assert np.allclose(p, [1013.25, 843.07, 696.82], atol=0.1)
    h = weather.pressure_altitude_from_pressure(p)
    assert np.allclose(h, [0.0, 5000.0, 10000.0], atol=0.5)


def test_saturation_over_ice_below_water():
    t = np.array([-5.0, -15.0, -25.0])
    ratio = weather.esat_water(t) / weather.esat_ice(t)
    # known ratios ~1.049, 1.155, 1.274 (e.g. Pruppacher & Klett)
    assert np.allclose(ratio, [1.049, 1.155, 1.274], atol=0.01)
    assert np.isclose(weather.esat_water(0.0), 611.2, atol=1.0)


def test_segment_splits_on_gap_and_hex():
    pos = pd.DataFrame({"hex": [1, 1, 1, 2, 2], "t": [0, 10, 3000, 5, 15]})
    s = flights.segment(pos, gap_s=1200)
    assert list(s.sort_values(["hex", "t"]).flight) == [0, 0, 1, 2, 2]


def test_interp_profile_time_and_log_pressure():
    runs = [pd.Timestamp("2025-01-26 01:00", tz="UTC"), pd.Timestamp("2025-01-26 02:00", tz="UTC")]
    prof = pd.DataFrame({
        "run": [runs[0]] * 2 + [runs[1]] * 2,
        "p": [800, 700, 800, 700],
        "t": [-10.0, -20.0, -12.0, -22.0],
    })
    pts = pd.DataFrame({"t": [runs[0].timestamp() + 1800], "p_hpa": [np.sqrt(800 * 700)]})
    out = flights.interp_profile_to_points(pts, prof, fields=("t",))
    # halfway in time (-11 / -21) and halfway in log-pressure -> -16
    assert np.isclose(out.t.iloc[0], -16.0)
