import numpy as np
import pytest

from snow import hrrr


def _chunk(dtype, fill):
    a = np.linspace(250.0, 280.0, hrrr.CHUNK * hrrr.CHUNK).reshape(hrrr.CHUNK, hrrr.CHUNK)
    a[0, 0] = fill
    return a.astype(dtype)


def test_decode_tile_float16_and_float32():
    # 2023 runs are stored as float16 (fill -10000), 2025 runs as float32 (fill -9999)
    for dtype, fill in (("<f2", -10000.0), ("<f4", -9999.0)):
        a = _chunk(dtype, fill)
        out = hrrr._decode_tile(a.tobytes())
        assert out.dtype == np.float32 and out.shape == (hrrr.CHUNK, hrrr.CHUNK)
        assert np.isnan(out[0, 0])
        assert np.allclose(out[1:, 1:], a[1:, 1:].astype(np.float32))


def test_decode_tile_rejects_other_sizes():
    with pytest.raises(ValueError):
        hrrr._decode_tile(b"\x00" * 1234)
