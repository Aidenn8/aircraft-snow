"""Find streak-shaped radar echoes without reference to aircraft.

Each scan is resampled to a GRID_RES_M Cartesian grid (reflectivity floored at
FLOOR_DBZ). A white top-hat with a disk of radius TOPHAT_RADIUS_M (the image minus its
morphological opening) keeps bright features narrower than the disk and removes broad
precipitation. Pixels with top-hat >= MIN_TOPHAT_DB and reflectivity >= MIN_DBZ are
grouped into 8-connected components. The image is first smoothed (Gaussian, SMOOTH_M) and
the mask cleaned (small fragments removed, small holes filled) so that speckle does not
masquerade as thin structure. Shape is measured on the component's skeleton so that bent
streaks (an approach path's downwind, base and final legs make a U) qualify:
    main length   longest path through the skeleton (graph diameter, 8-neighbour steps)
    total length  all skeleton pixels
    mean width    area / total length
    branchiness   total length / main length
A component is a streak candidate when its main length is at least MIN_LENGTH_M, its mean
width at most MAX_MEAN_WIDTH_M, its branchiness at most MAX_BRANCHINESS, its main length at
least MIN_ELONGATION times its mean width, and its area at
least MIN_AREA_KM2.
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pandas as pd
from scipy import ndimage
from skimage import measure, morphology

from . import radar

GRID_HALF_M = 70_000.0
GRID_RES_M = 500.0
FLOOR_DBZ = -10.0
TOPHAT_RADIUS_M = 4000.0
MIN_TOPHAT_DB = 6.0
MIN_DBZ = 3.0
MIN_LENGTH_M = 8000.0
MAX_MEAN_WIDTH_M = 3500.0
MAX_BRANCHINESS = 2.5
MIN_ELONGATION = 4.0      # main length / mean width (Denver streaks: 10-30 km x 2-3 km)
MIN_AREA_KM2 = 4.0
MIN_RANGE_M = 8000.0
SMOOTH_M = 750.0


@dataclasses.dataclass
class GridScan:
    t: float
    xs: np.ndarray
    ys: np.ndarray
    z: np.ndarray          # floored reflectivity on the grid
    valid: np.ndarray      # inside the sweep and beyond MIN_RANGE_M


def grid_scan(t: float, sweep: radar.Sweep, half_m=GRID_HALF_M, res_m=GRID_RES_M) -> GridScan:
    xs, ys, z = radar.to_grid(sweep, half_m, res_m)
    X, Y = np.meshgrid(xs, ys)
    R = np.hypot(X, Y)
    valid = (R >= MIN_RANGE_M) & (R <= half_m)
    zf = np.where(np.isnan(z), FLOOR_DBZ, np.maximum(z, FLOOR_DBZ)).astype(np.float32)
    zf[~valid] = FLOOR_DBZ
    return GridScan(t, xs, ys, zf, valid)


def tophat(z: np.ndarray, radius_m=TOPHAT_RADIUS_M, res_m=GRID_RES_M) -> np.ndarray:
    fp = morphology.disk(int(round(radius_m / res_m)))
    return z - ndimage.grey_opening(z, footprint=fp)


def skeleton_lengths(img: np.ndarray, res_m: float = GRID_RES_M) -> tuple[float, float]:
    """(main length, total length) in metres of the skeleton of a binary image. The main
    length is the graph diameter found by two breadth-first sweeps (exact for trees)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra

    skel = morphology.skeletonize(img)
    yy, xx = np.nonzero(skel)
    n = len(yy)
    if n < 2:
        return float(n * res_m), float(n * res_m)
    index = -np.ones(skel.shape, np.int64)
    index[yy, xx] = np.arange(n)
    rows, cols, w = [], [], []
    for dy, dx in ((0, 1), (1, 0), (1, 1), (1, -1)):
        y2, x2 = yy + dy, xx + dx
        ok = (y2 >= 0) & (y2 < skel.shape[0]) & (x2 >= 0) & (x2 < skel.shape[1])
        j = np.full(n, -1)
        j[ok] = index[y2[ok], x2[ok]]
        m = j >= 0
        rows.append(np.arange(n)[m]); cols.append(j[m])
        w.append(np.full(m.sum(), np.hypot(dy, dx)))
    r, c, ww = np.concatenate(rows), np.concatenate(cols), np.concatenate(w)
    G = coo_matrix((np.r_[ww, ww], (np.r_[r, c], np.r_[c, r])), shape=(n, n)).tocsr()
    d0 = dijkstra(G, indices=0)
    far = int(np.argmax(np.where(np.isfinite(d0), d0, -1)))
    d1 = dijkstra(G, indices=far)
    main = float(np.max(d1[np.isfinite(d1)])) * res_m
    total = float(ww.sum()) * res_m
    return max(main, res_m), max(total, res_m)


def find_streaks(g: GridScan) -> tuple[pd.DataFrame, np.ndarray]:
    """Streak candidates in one scan: a table of properties and the label image."""
    zs = ndimage.gaussian_filter(g.z, SMOOTH_M / GRID_RES_M)
    th = tophat(zs)
    mask = (th >= MIN_TOPHAT_DB) & (zs >= MIN_DBZ) & g.valid
    mask = morphology.remove_small_objects(mask, max_size=7)
    mask = morphology.remove_small_holes(mask, max_size=19)
    lab = measure.label(mask, connectivity=2)
    rows = []
    keep = np.zeros(lab.max() + 1, bool)
    for rp in measure.regionprops(lab, intensity_image=g.z):
        area_km2 = rp.area * (GRID_RES_M / 1000) ** 2
        if area_km2 < MIN_AREA_KM2:
            continue
        L, L_tot = skeleton_lengths(rp.image)
        W = rp.area * GRID_RES_M ** 2 / L_tot
        if (L < MIN_LENGTH_M or W > MAX_MEAN_WIDTH_M or L_tot / L > MAX_BRANCHINESS
                or L / W < MIN_ELONGATION):
            continue
        keep[rp.label] = True
        cy, cx = rp.centroid
        sel = lab == rp.label
        zz = g.z[sel]
        rows.append({"t": g.t, "label": rp.label, "x": np.interp(cx, np.arange(len(g.xs)), g.xs),
                     "y": np.interp(cy, np.arange(len(g.ys)), g.ys), "length_m": L,
                     "skeleton_m": L_tot, "width_m": W,
                     "area_km2": area_km2, "z_mean": float(zz.mean()), "z_max": float(zz.max()),
                     "tophat_mean": float(th[sel].mean())})
    lab = np.where(keep[lab], lab, 0)
    return pd.DataFrame(rows), lab
