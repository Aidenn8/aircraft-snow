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
        "temp": [-10.0, -20.0, -12.0, -22.0],
    })
    t_pt = runs[0].timestamp() + 1800
    pts = pd.DataFrame({"t": [t_pt], "p_hpa": [np.sqrt(800 * 700)]})
    out = flights.interp_profile_to_points(pts, prof, fields=("temp",))
    # halfway in time (-11 / -21) and halfway in log-pressure -> -16
    assert np.isclose(out.temp.iloc[0], -16.0)
    assert out.t.iloc[0] == t_pt  # point times untouched


def test_densify_spacing_and_endpoints():
    pts = pd.DataFrame({"flight": [0, 0, 0, 1], "x": [0.0, 1000.0, 3000.0, 0.0],
                        "y": [0.0, 0.0, 0.0, 5.0], "t": [0.0, 10.0, 20.0, 0.0],
                        "alt_ft": [1000.0, 1100.0, 1300.0, 5000.0], "p_hpa": [900.0] * 4})
    d = flights.densify(pts, step_m=250.0)
    f0 = d[d.flight == 0]
    assert np.all(np.diff(f0.x) <= 250.0 + 1e-9)
    assert f0.x.iloc[0] == 0.0 and f0.x.iloc[-1] == 3000.0
    assert np.isclose(f0[f0.x == 500.0].t.iloc[0], 5.0)          # linear in time too
    assert np.isclose(f0[f0.x == 2000.0].alt_ft.iloc[0], 1200.0)
    assert len(d[d.flight == 1]) == 1                            # single point kept as is


def test_densify_does_not_bridge_gaps():
    pts = pd.DataFrame({"flight": [0, 0], "x": [0.0, 10000.0], "y": [0.0, 0.0],
                        "t": [0.0, 300.0], "alt_ft": [1000.0, 1000.0], "p_hpa": [900.0, 900.0]})
    d = flights.densify(pts, step_m=250.0, max_gap_s=60.0)
    assert len(d) == 2


def test_interp_refuses_field_named_t():
    import pytest
    with pytest.raises(ValueError):
        flights.interp_profile_to_points(pd.DataFrame({"t": [0.0], "p_hpa": [700.0]}),
                                         pd.DataFrame({"run": [], "p": [], "t": []}), fields=("t",))
