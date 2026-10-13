"""Shared figure style for the report (light surface, recessive chrome, validated palette).

apply(dark=True) switches to the web page's look instead: a black surface, white Inter text
(bundled in tools/fonts) and the reflectivity ramp turned round so weak echo is dark and strong
echo bright. Call it before drawing: the colour constants below are rebound in place.

Colours are the reference palette of the data-viz guidance used for this project:
categorical slots in fixed order, a single-hue sequential ramp for magnitude (radar
reflectivity is a magnitude, so no rainbow), solid hairline grids, text in ink tokens.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
BLUE, ORANGE, AQUA, YELLOW = SERIES[:4]

# sequential blue ramp (100 -> 700) for reflectivity; the lightest step sits near the surface
BLUE_RAMP = ["#f4f8fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
REFLECTIVITY = LinearSegmentedColormap.from_list("reflectivity_blue", BLUE_RAMP)
REFLECTIVITY.set_bad(SURFACE)
REFLECTIVITY.set_under(SURFACE)


DARK = {"SURFACE": "#000000", "INK": "#ffffff", "INK_2": "#d4d4d4", "MUTED": "#8e8e8e",
        "GRID": "#1e1e1e", "AXIS": "#3a3a3a"}
DARK_RAMP = ["#0b1d33", "#123058", "#1a4a86", "#2a6fc0", "#4f95e6", "#8fbdf2", "#d3e6fb"]
FONTS = Path(__file__).resolve().parents[1] / "tools" / "fonts"


def _use_dark():
    global SURFACE, INK, INK_2, MUTED, GRID, AXIS, REFLECTIVITY
    SURFACE, INK, INK_2, MUTED, GRID, AXIS = (DARK[k] for k in ("SURFACE", "INK", "INK_2", "MUTED", "GRID", "AXIS"))
    REFLECTIVITY = LinearSegmentedColormap.from_list("reflectivity_blue_dark", DARK_RAMP)
    REFLECTIVITY.set_bad(SURFACE)
    REFLECTIVITY.set_under(SURFACE)
    for f in sorted(FONTS.glob("Inter-*.ttf")):
        font_manager.fontManager.addfont(str(f))


def apply(dark: bool = False):
    if dark:
        _use_dark()
    mpl.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK_2, "axes.titlecolor": INK,
        "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.labelsize": 9.5, "axes.linewidth": 0.8, "axes.grid": True,
        "axes.spines.top": False, "axes.spines.right": False, "axes.prop_cycle": mpl.cycler(color=SERIES),
        "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "xtick.major.size": 0, "ytick.major.size": 0,
        "legend.frameon": False, "legend.fontsize": 8.5, "legend.labelcolor": INK_2,
        "font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"], "font.size": 9.5,
        "text.color": INK, "lines.linewidth": 2.0, "lines.solid_capstyle": "round",
        "lines.solid_joinstyle": "round", "image.cmap": "viridis", "figure.dpi": 100,
        "savefig.dpi": 200, "savefig.bbox": "tight",
    })
    if dark and any(FONTS.glob("Inter-*.ttf")):
        mpl.rcParams["font.family"] = ["Inter", "Helvetica Neue", "Arial", "DejaVu Sans"]


def map_axes(ax, half_km=None):
    """Plan-view axes: equal aspect, light hairline grid, km labels."""
    ax.set_aspect(1)
    ax.grid(True, color=GRID, linewidth=0.5)
    if half_km is not None:
        ax.set_xlim(-half_km, half_km)
        ax.set_ylim(-half_km, half_km)
    ax.set_xlabel("km east of the radar")
    ax.set_ylabel("km north of the radar")
