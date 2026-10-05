import numpy as np

from snow import streaks


def _grid(z):
    n = z.shape[0]
    xs = (np.arange(n) - n // 2) * streaks.GRID_RES_M
    valid = np.ones_like(z, bool)
    return streaks.GridScan(0.0, xs, xs.copy(), z.astype(np.float32), valid)


def test_skeleton_lengths_bar_and_u():
    img = np.zeros((60, 60), bool)
    img[30:33, 5:45] = True                      # 40 px long, 3 px wide bar
    main, total = streaks.skeleton_lengths(img, res_m=1.0)
    assert 36 <= main <= 41 and total / main < 1.2
    u = np.zeros((80, 80), bool)
    u[10:60, 10:13] = True                       # left leg
    u[57:60, 10:50] = True                       # base
    u[10:60, 47:50] = True                       # right leg
    main, total = streaks.skeleton_lengths(u, res_m=1.0)
    assert main > 120 and total / main < 1.3     # the U is one long path


def test_find_streaks_accepts_line_rejects_blob_and_speckle():
    n = 241
    z = np.full((n, n), -10.0)
    # a 20 km long, 1.5 km wide streak of 15 dBZ (500 m grid)
    z[60:63, 40:80] = 15.0
    # a broad 15 km disk of 15 dBZ (precipitation area)
    yy, xx = np.mgrid[:n, :n]
    z[(yy - 170) ** 2 + (xx - 170) ** 2 < 15 ** 2] = 15.0
    # speckle: isolated random pixels
    rng = np.random.default_rng(0)
    sp = rng.random((n, n)) < 0.02
    sp[:, :120] = False
    sp[150:, :] = False
    z[sp] = 20.0
    obj, lab = streaks.find_streaks(_grid(z))
    assert len(obj) == 1
    o = obj.iloc[0]
    assert 17_000 <= o.length_m <= 21_000
    assert o.width_m < 2500


def test_find_streaks_rejects_bright_convective_line():
    n = 241
    z = np.full((n, n), -10.0)
    z[60:63, 40:80] = 15.0                       # light snow streak: kept
    z[160:163, 40:80] = 45.0                     # same shape, convective intensity: rejected
    obj, lab = streaks.find_streaks(_grid(z))
    assert len(obj) == 1 and obj.iloc[0].z_max <= streaks.MAX_PEAK_DBZ
    assert (lab[160:163, 40:80] == 0).all()
