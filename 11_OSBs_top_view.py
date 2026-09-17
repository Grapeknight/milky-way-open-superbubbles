#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script 11: OSBs top view

Purpose
-------
Draw face-on open-superbubble overview and RGB layered dust map.

Method overview
---------------
1. Load the fitted projected ellipses, raw dust cube and spiral-arm loci.
2. Form the configured vertical dust layers and apply the selected projected-shell
   convention.
3. Render the clean face-on overview and RGB dust-layer map with common annotations and
   legends.

Main inputs
-----------
- ../results/superbubble_final_fit_parameters.csv
- ../results/intermediate_output/3d_dust_map_products/raw_3d_dust_cube.parquet
- ../data/Reid2019_spirals_xy.csv

Main outputs
------------
- ../results/figures/11_topview_rgb_dust_layers.png; source panel for main-text Fig. 1.

Figure/table role
-----------------
../results/figures/11_topview_rgb_dust_layers.png; source panel for main-text Fig. 1.

Runtime and data notes
----------------------
Reference runtime: 21.01 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

External-data caveat
--------------------
Manuscript Fig. 1 is assembled from this output and a JWST image outside this script.

Reading the code
----------------
Start with the path and scientific settings below, then follow main() at
the end of the file. The method overview above describes the order of the
analysis steps; the input/output lists identify the upstream data and
products written by this stage. Relative data paths are resolved from code/.
"""

from __future__ import annotations
import time

from collections import Counter
import os
from pathlib import Path
from typing import Optional

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import pyarrow.parquet as pq
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse
from scipy.ndimage import gaussian_filter

# ----------------------------------------------------------------------------
# Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
# ----------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "11_topview_rgb_dust_layers"

GEOM_CSV = FINAL_DIR / "superbubble_final_fit_parameters.csv"
DUST_PARQUET = Path("..") / "results" / "intermediate_output" / "3d_dust_map_products" / "raw_3d_dust_cube.parquet"
ARM_CSV = DATA_DIR / "Reid2019_spirals_xy.csv"

OUT_ELLIPSE_CSV = OUT_DIR / "11_plotting_ellipse_parameters.csv"
OUT_FIG = FINAL_FIG_DIR / "11_topview_rgb_dust_layers.png"

# ----------------------------------------------------------------------------
# Figure styling and layout configuration for reproducible rendering.
# ----------------------------------------------------------------------------
PLOT_HALF_KPC = 3.5                      # Figure-panel rendering section.
Z_SLICES = [(-0.30, -0.05),
            (-0.05, 0.05),
            (0.05, 0.30)]
RGB_SIGMA = 2.0
RGB_VMIN = 0.0
RGB_VMAX_MAG = 0.05                      # RGB dust-layer rendering section.
RGB_COLORBAR_TICKS = np.arange(0.0, RGB_VMAX_MAG + 0.001, 0.01)
SUN_GC_DISTANCE_KPC = 8.12               # Spiral-arm overlay and coordinate-conversion settings.

# Compromise and maximum-opening ellipse conventions for representative OSB footprints.
# Compromise and maximum-opening ellipse conventions for representative OSB footprints.
# Compromise and maximum-opening ellipse conventions for representative OSB footprints.
ELLIPSOID_CAP_SECTION_FRAC = 0.5

# Geometry and shell-mask conventions used by this analysis stage.
ROW_LABELS = {
    "maxopen": "Maximum-opening method (ellipsoid cap base)",
    "compromise": "Compromise method (mid-cap section)",
}

COLOR_MAP = {0: "black", 1: "blue", 2: "red", 3: "green"}
LABEL_MAP = {0: "Closed", 1: "North Open", 2: "South Open", 3: "Disk-penetrating"}

# Radcliffe Wave reference geometry and slice configuration.
RAD_LINE_ANGLE_DEG = 60.0
RAD_LINE_INTERCEPT_KPC = 0.6
RAD_LINE_HALF_WIDTH_KPC = 0.1

# Spiral-arm overlay and coordinate-conversion settings.
ARM_STYLES = {
    "Perseus": ("r", "-", 2),
    "Local": ("magenta", "-", 2),
    "Sagittarius": ("lime", "-", 2),
    "Carina": ("cyan", "-", 2),
    "LoS": ("black", "--", 1.5),
}
ARM_DISPLAY = {
    "Perseus": "Perseus arm",
    "Local": "Local arm",
    "Sagittarius": "Sagittarius arm",
    "Carina": "Carina arm",
    "LoS": "Local arm spur",
}

# Figure styling and layout configuration for reproducible rendering.
FONT_SIZES = {
    "axis_label": 14,
    "tick_label": 11,
    "panel_tag": 16,
    "legend_text": 9,
    "legend_title": 10,
    "annotation": 11,
    "rgb_subtitle": 13,
}


# ----------------------------------------------------------------------------
# Data-loading helper section.
# ----------------------------------------------------------------------------
def midcap_section_scale(best_fit: pd.DataFrame) -> np.ndarray:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    shape = best_fit["shape"].astype(str).str.lower().to_numpy()
    cz = best_fit["center_z_kpc"].to_numpy(float)
    c = best_fit["c_radius_kpc"].to_numpy(float)
    z_base = best_fit["xy_plane_z_kpc"].to_numpy(float)
    mark = best_fit["mark"].to_numpy(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        k_base = np.sqrt(np.clip(1.0 - ((z_base - cz) / c) ** 2, 0.0, None))
        z_apex = np.where(mark == 1, cz - c, cz + c)
        z_sec = z_base + ELLIPSOID_CAP_SECTION_FRAC * (z_apex - z_base)
        k_sec = np.sqrt(np.clip(1.0 - ((z_sec - cz) / c) ** 2, 0.0, None))
        scale = np.where(k_base > 0, k_sec / k_base, 1.0)
    return np.where(shape == "ellipsoid", scale, 1.0)


def load_xy_plane_ellipses(csv_path: Path) -> Optional[pd.DataFrame]:
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    required = [
        "id", "mark", "shape",
        "xy_plane_center_x_kpc", "xy_plane_center_y_kpc",
        "xy_plane_a_kpc", "xy_plane_b_kpc", "xy_plane_angle_deg",
        "center_z_kpc", "c_radius_kpc", "xy_plane_z_kpc",
    ]
    if not csv_path.exists():
        print(f"Best-fit table not found: {csv_path}")
        return None
    best_fit = pd.read_csv(csv_path)
    missing = [c for c in required if c not in best_fit.columns]
    if missing:
        print(f"Best-fit table missing columns {missing}: {csv_path}")
        return None

    scale = midcap_section_scale(best_fit)
    ellipse = best_fit[[
        "id", "mark", "shape",
        "xy_plane_center_x_kpc", "xy_plane_center_y_kpc",
        "xy_plane_a_kpc", "xy_plane_b_kpc", "xy_plane_angle_deg",
    ]].copy()
    ellipse = ellipse.rename(columns={
        "id": "index",
        "xy_plane_center_x_kpc": "x",
        "xy_plane_center_y_kpc": "y",
        "xy_plane_a_kpc": "a",            # Compromise and maximum-opening ellipse conventions for representative OSB footprints.
        "xy_plane_b_kpc": "b",            # Compromise and maximum-opening ellipse conventions for representative OSB footprints.
        "xy_plane_angle_deg": "angle",
    })
    ellipse["section_scale"] = scale       # Compromise and maximum-opening ellipse conventions for representative OSB footprints.
    ellipse["sort_z"] = best_fit["center_z_kpc"]
    return ellipse


def ellipses_for_method(ellipse_df: Optional[pd.DataFrame], method: str) -> Optional[pd.DataFrame]:
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    if ellipse_df is None:
        return None
    out = ellipse_df.copy()
    if method == "compromise":
        out["a"] = out["a"] * out["section_scale"]
        out["b"] = out["b"] * out["section_scale"]
    return out


def load_arms(csv_path: Path) -> Optional[pd.DataFrame]:
    """Load and validate the input table or array needed by this stage."""
    if not csv_path.exists():
        print(f"Optional spiral-arm table not found: {csv_path}")
        return None
    return pd.read_csv(csv_path)


def load_dust_cube(parquet_path: Path) -> pd.DataFrame:
    """Load and validate the input table or array needed by this stage."""
    half = PLOT_HALF_KPC
    z_lo = min(z[0] for z in Z_SLICES)
    z_hi = max(z[1] for z in Z_SLICES)
    names = set(pq.ParquetFile(str(parquet_path)).schema.names)
    if {"x", "y", "z", "dust"}.issubset(names):
        cols = ["x", "y", "z", "dust"]
        rename = {}
    elif {"X", "Y", "Z", "dust"}.issubset(names):
        cols = ["X", "Y", "Z", "dust"]
        rename = {"X": "x", "Y": "y", "Z": "z"}
    elif {"x", "y", "z", "dust_raw"}.issubset(names):
        cols = ["x", "y", "z", "dust_raw"]
        rename = {"dust_raw": "dust"}
    elif {"X", "Y", "Z", "dust_raw"}.issubset(names):
        cols = ["X", "Y", "Z", "dust_raw"]
        rename = {"X": "x", "Y": "y", "Z": "z", "dust_raw": "dust"}
    else:
        raise ValueError(f"Could not identify XYZ parquet columns: {sorted(names)}")
    x_col, y_col, z_col = cols[0], cols[1], cols[2]
    filters = [
        (z_col, ">=", z_lo), (z_col, "<=", z_hi),
        (x_col, ">=", -half), (x_col, "<=", half),
        (y_col, ">=", -half), (y_col, "<=", half),
    ]
    try:
        return pd.read_parquet(parquet_path, columns=cols, filters=filters).rename(columns=rename)
    except Exception as exc:
        print(f"Parquet filter pushdown failed; falling back to full read: {exc}")
        df = pd.read_parquet(parquet_path, columns=cols).rename(columns=rename)
        return df[
            (df["z"] >= z_lo) & (df["z"] <= z_hi)
            & (df["x"] >= -half) & (df["x"] <= half)
            & (df["y"] >= -half) & (df["y"] <= half)
        ].copy()


def clip_xy_range_to_dust_data(df: pd.DataFrame, x_range, y_range):
    def axis_range(column: str, requested):
        values = np.sort(df[column].dropna().unique().astype(float))
        if values.size == 0:
            raise ValueError('Invalid input or missing required data.')
        diffs = np.diff(values)
        positive = diffs[diffs > 0]
        step = float(np.min(positive)) if positive.size else 0.01
        data_min = float(values[0] - 0.5 * step)
        data_max = float(values[-1] + 0.5 * step)
        clipped = (max(float(requested[0]), data_min), min(float(requested[1]), data_max))
        if not clipped[0] < clipped[1]:
            raise ValueError('Invalid input or missing required data.')
        return clipped

    return axis_range("x", x_range), axis_range("y", y_range)


# ----------------------------------------------------------------------------
# RGB dust-layer rendering section.
# ----------------------------------------------------------------------------
def extract_xy_grid(
    df: pd.DataFrame,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    z_range: tuple[float, float],
    sigma: float,
    vmin: float,
    vmax: float,
) -> Optional[np.ndarray]:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    mask = (
        (df["x"] >= x_range[0]) & (df["x"] <= x_range[1])
        & (df["y"] >= y_range[0]) & (df["y"] <= y_range[1])
        & (df["z"] >= z_range[0]) & (df["z"] <= z_range[1])
    )
    df_cut = df[mask].copy()
    if df_cut.empty:
        return None

    df_cut["x"] = df_cut["x"].round(3)
    df_cut["y"] = df_cut["y"].round(3)

    group = df_cut.groupby(["x", "y"])["dust"]
    mean_val = group.transform("mean")
    std_val = group.transform("std").fillna(1e-9)
    df_cut = df_cut[
        (df_cut["dust"] >= mean_val - 3 * std_val)
        & (df_cut["dust"] <= mean_val + 3 * std_val)
    ]

    g = df_cut.groupby(["x", "y"])["dust"].mean().reset_index()
    if g.empty:
        return None

    grid = g.pivot(index="y", columns="x", values="dust")
    grid = grid.sort_index(ascending=True).sort_index(axis=1, ascending=True)

    z_width_kpc = float(z_range[1] - z_range[0])
    img_array = gaussian_filter(grid.fillna(0).values, sigma=sigma) * z_width_kpc
    img = ma.masked_array(img_array, mask=grid.isnull().values)

    norm = Normalize(vmin=vmin, vmax=vmax)
    normalized = norm(img.data)
    normalized[img.mask] = 0
    return np.clip(normalized, 0, 1)


# ----------------------------------------------------------------------------
# Spiral-arm overlay and coordinate-conversion settings.
# ----------------------------------------------------------------------------
def draw_arms(ax, arm_df: Optional[pd.DataFrame], main_alpha: float = 0.8):
    """Render or save a documented figure product without changing upstream measurements."""
    handles = []
    if arm_df is None:
        return handles
    for name, (color, ls, lw) in ARM_STYLES.items():
        mask = arm_df["arm"] == name
        if not mask.any():
            continue
        line, = ax.plot(
            -(arm_df.loc[mask, "yy"]) + SUN_GC_DISTANCE_KPC, arm_df.loc[mask, "xx"],
            linestyle=ls, color=color, linewidth=lw, alpha=main_alpha,
            label=ARM_DISPLAY[name],
        )
        ax.plot(-(arm_df.loc[mask, "yy0"]) + SUN_GC_DISTANCE_KPC, arm_df.loc[mask, "xx0"],
                linestyle=":", color=color, linewidth=1.2, alpha=0.6)
        ax.plot(-(arm_df.loc[mask, "yy1"]) + SUN_GC_DISTANCE_KPC, arm_df.loc[mask, "xx1"],
                linestyle=":", color=color, linewidth=1.2, alpha=0.6)
        handles.append(line)
    return handles


def ellipse_dims(row) -> tuple[float, float]:
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    return 2.0 * float(row["a"]), 2.0 * float(row["b"])


def configure_matplotlib():
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    mpl.rcParams["font.family"] = "DejaVu Sans"
    mpl.rcParams["axes.unicode_minus"] = False


# ----------------------------------------------------------------------------
# Figure-panel rendering section.
# ----------------------------------------------------------------------------
def draw_clean_map(ax, ellipse_df, arm_df):
    """Render or save a documented figure product without changing upstream measurements."""
    ax.set_facecolor("whitesmoke")


    circle_kwargs = {"edgecolor": "gray", "linewidth": 1, "alpha": 0.6}
    line_kwargs = {"color": "gray", "linestyle": ":", "linewidth": 0.8, "alpha": 0.6}
    for radius in (1, 2, 3):
        ax.add_patch(Circle((0, 0), radius, facecolor="none", **circle_kwargs))
        ax.text(0, radius, f"{radius} kpc", color="gray",
                fontsize=FONT_SIZES["annotation"] - 2, ha="center", va="bottom")


    for angle in np.arange(0, 360, 30):
        theta = np.deg2rad(angle)
        ax.plot([0, 4.5 * np.cos(theta)], [0, 4.5 * np.sin(theta)], **line_kwargs)
    for angle in (0, 180, 270):
        theta = np.deg2rad(angle)
        ax.text(3.1 * np.cos(theta), 3.1 * np.sin(theta), f"{angle}°",
                color="black", fontsize=FONT_SIZES["annotation"], ha="center", va="center")

    # Spiral-arm overlay and coordinate-conversion settings.
    draw_arms(ax, arm_df, main_alpha=0.7)


    if ellipse_df is not None:
        df = ellipse_df.copy()
        if "sort_z" in df.columns:
            df = df.sort_values("sort_z", ascending=True)
        for _, row in df.iterrows():
            if pd.isna(row["mark"]):
                continue
            mark = int(row["mark"])
            c = COLOR_MAP.get(mark, "gray")
            w, h = ellipse_dims(row)
            ax.add_patch(Ellipse((row["x"], row["y"]), w, h, angle=row["angle"],
                                 facecolor=c, edgecolor="black", alpha=0.6))
            tc = "white" if c in ("black", "blue") else "black"
            ax.text(row["x"], row["y"], str(int(row["index"])), color=tc,
                    fontsize=8, ha="center", va="center", fontweight="bold",
                    path_effects=[pe.withStroke(linewidth=1.5,
                                                foreground="black" if tc == "white" else "white")])

    # Radcliffe Wave reference geometry and slice configuration.
    x_l = np.linspace(-4, 4, 100)
    slope = np.tan(np.deg2rad(RAD_LINE_ANGLE_DEG))
    y_l = slope * x_l + RAD_LINE_INTERCEPT_KPC
    offset = RAD_LINE_HALF_WIDTH_KPC * np.sqrt(1 + slope ** 2)
    ax.fill_between(x_l, y_l - offset, y_l + offset, color="gray", alpha=0.3)
    ax.plot(x_l, y_l, color="black", linewidth=2)


    ax.set_xlim(-PLOT_HALF_KPC, PLOT_HALF_KPC)
    ax.set_ylim(-PLOT_HALF_KPC, PLOT_HALF_KPC)
    ax.set_aspect("equal")
    ax.set_xlabel("X [kpc]", fontsize=FONT_SIZES["axis_label"])
    ax.set_ylabel("Y [kpc]", fontsize=FONT_SIZES["axis_label"])
    ax.tick_params(labelsize=FONT_SIZES["tick_label"])
    ax.grid(True, linestyle="--", alpha=0.4)


def ellipse_mark_counts(ellipse_df) -> Counter:
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    counts = Counter()
    if ellipse_df is None:
        return counts
    for mark in ellipse_df["mark"].dropna():
        counts[int(mark)] += 1
    return counts


def add_shared_legend(fig, counts: Counter):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    total = int(sum(counts.values()))
    legend_style = dict(
        fontsize=FONT_SIZES["legend_text"] + 4,
        title_fontsize=FONT_SIZES["legend_title"] + 4,
        frameon=True,
        facecolor="white",
        framealpha=0.95,
        columnspacing=1.1,
        handlelength=2.2,
        handletextpad=0.5,
        borderpad=0.8,
    )

    legends = []
    legends.append(fig.legend(
        [Line2D([0], [0], color="black", linestyle="-", linewidth=2)],
        ["RW line: $y = \\tan(60^\\circ)x + 0.6$"],
        loc="lower center",
        bbox_to_anchor=(0.19, -0.030),
        ncol=1,
        title="Reference",
        **legend_style,
    ))

    osb_handles, osb_labels = [], []
    compact_osb_labels = {
        1: "North",
        2: "South",
        3: "Disk-penetrating",
    }
    for code in (1, 2, 3):
        osb_handles.append(Line2D([0], [0], marker="o", linestyle="None",
                                  markerfacecolor=COLOR_MAP[code],
                                  markeredgecolor="black", markersize=9))
        label = compact_osb_labels.get(code, LABEL_MAP[code])
        osb_labels.append(f"{label} ({int(counts[code]) if code in counts else 0})")
    legends.append(fig.legend(
        osb_handles, osb_labels,
        loc="lower center",
        bbox_to_anchor=(0.49, -0.030),
        ncol=3,
        title=f"OSB types, total: {total}",
        **legend_style,
    ))

    arm_handles = [
        Line2D([0], [0], color=color, linestyle=linestyle, linewidth=2)
        for color, linestyle, _lw in ARM_STYLES.values()
    ]
    arm_labels = list(ARM_DISPLAY.values())
    legends.append(fig.legend(
        arm_handles, arm_labels,
        loc="lower center",
        bbox_to_anchor=(0.84, -0.030),
        ncol=2,
        title="Spiral arms",
        **legend_style,
    ))

    for legend in legends:
        legend.get_frame().set_linewidth(0.8)


# ----------------------------------------------------------------------------
# RGB dust-layer rendering section.
# ----------------------------------------------------------------------------
def draw_rgb_map(fig, ax, dust_df, ellipse_df, arm_df):
    """Render or save a documented figure product without changing upstream measurements."""
    x_range = (-PLOT_HALF_KPC, PLOT_HALF_KPC)
    y_range = (-PLOT_HALF_KPC, PLOT_HALF_KPC)
    x_range, y_range = clip_xy_range_to_dust_data(dust_df, x_range, y_range)
    sigma2 = 2 * RGB_SIGMA

    r_img = extract_xy_grid(dust_df, x_range, y_range, Z_SLICES[0], sigma2, RGB_VMIN, RGB_VMAX_MAG)
    g_img = extract_xy_grid(dust_df, x_range, y_range, Z_SLICES[1], RGB_SIGMA, RGB_VMIN, RGB_VMAX_MAG)
    b_img = extract_xy_grid(dust_df, x_range, y_range, Z_SLICES[2], sigma2, RGB_VMIN, RGB_VMAX_MAG)
    if any(img is None for img in (r_img, g_img, b_img)):
        print("Skipping RGB top-view image because at least one z-slice is empty.")
        return

    rgb = np.stack([r_img, g_img, b_img], axis=-1)
    extent = [x_range[0], x_range[1], y_range[0], y_range[1]]
    ax.imshow(rgb, origin="lower", extent=extent)


    circle_kwargs = {"edgecolor": "white", "facecolor": "none",
                     "linestyle": "--", "linewidth": 0.7, "alpha": 0.5}
    line_kwargs = {"color": "white", "linestyle": "--", "linewidth": 0.6, "alpha": 0.5}
    for radius in (1, 2, 3, 4):
        ax.add_patch(Circle((0, 0), radius, **circle_kwargs))
    radial_length = PLOT_HALF_KPC + 1.0
    for angle in np.arange(0, 360, 30):
        theta = np.deg2rad(angle)
        ax.plot([0, radial_length * np.cos(theta)], [0, radial_length * np.sin(theta)], **line_kwargs)
    for angle in (0, 90, 180, 270):
        theta = np.deg2rad(angle)
        ax.text(3.1 * np.cos(theta), 3.1 * np.sin(theta), f"{angle}°", color="white",
                fontsize=10, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="black", ec="none", alpha=0.6))

    # Spiral-arm overlay and coordinate-conversion settings.
    if arm_df is not None:
        for name, (color, linestyle, lw) in ARM_STYLES.items():
            mask = arm_df["arm"] == name
            if not mask.any():
                continue
            ax.plot(-(arm_df.loc[mask, "yy"]) + SUN_GC_DISTANCE_KPC, arm_df.loc[mask, "xx"],
                    linestyle=linestyle, color=color, linewidth=lw, alpha=0.8)
            ax.plot(-(arm_df.loc[mask, "yy0"]) + SUN_GC_DISTANCE_KPC, arm_df.loc[mask, "xx0"],
                    linestyle=":", color=color, linewidth=1.2, alpha=0.6)
            ax.plot(-(arm_df.loc[mask, "yy1"]) + SUN_GC_DISTANCE_KPC, arm_df.loc[mask, "xx1"],
                    linestyle=":", color=color, linewidth=1.2, alpha=0.6)


    if ellipse_df is not None:
        for _, row in ellipse_df.iterrows():
            if any(pd.isna(row[c]) for c in ("x", "y", "mark", "angle", "index")):
                continue
            mark = int(row["mark"])
            color = COLOR_MAP.get(mark, "black")
            w, h = ellipse_dims(row)
            ax.add_patch(Ellipse((row["x"], row["y"]), width=w, height=h, angle=row["angle"],
                                 facecolor="none", edgecolor="white", linewidth=2.5,
                                 alpha=0.9, zorder=10))
            ax.add_patch(Ellipse((row["x"], row["y"]), width=w, height=h, angle=row["angle"],
                                 facecolor="none", edgecolor=color, linewidth=2,
                                 alpha=0.8, zorder=11))


    z_labels = [f"[{z[0]:.2f}, {z[1]:.2f}]" for z in Z_SLICES]
    ax.set_title(f"R: Z={z_labels[0]}   G: Z={z_labels[1]}   B: Z={z_labels[2]}",
                 fontsize=FONT_SIZES["rgb_subtitle"], pad=14)
    ax.set_xlabel("X [kpc]", fontsize=FONT_SIZES["axis_label"])
    ax.set_ylabel("Y [kpc]", fontsize=FONT_SIZES["axis_label"])
    ax.tick_params(labelsize=FONT_SIZES["tick_label"])
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.set_aspect("equal", adjustable="box")

    _add_rgb_colorbars(fig, ax)


def _add_rgb_colorbars(fig, ax):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    cbar_width = 0.013
    cbar_gap = 0.012
    ax_pos = ax.get_position()
    channels = [
        {"cmap": "Reds", "vmax": RGB_VMAX_MAG},
        {"cmap": "Greens", "vmax": RGB_VMAX_MAG},
        {"cmap": "Blues", "vmax": RGB_VMAX_MAG},
    ]
    first_pos = last_pos = None
    for i, props in enumerate(channels):
        cax_left = ax_pos.x1 + 0.015 + i * (cbar_width + cbar_gap)
        cax = fig.add_axes([cax_left, ax_pos.y0, cbar_width, ax_pos.height])
        norm = mcolors.Normalize(vmin=RGB_VMIN, vmax=props["vmax"])
        sm = cm.ScalarMappable(cmap=props["cmap"], norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=cax, ticks=RGB_COLORBAR_TICKS)
        cb.ax.yaxis.set_ticks_position("right")
        if i == len(channels) - 1:
            cb.set_ticklabels([f"{tick:.2f}" for tick in RGB_COLORBAR_TICKS])
            cb.ax.tick_params(labelsize=11, labelleft=False, labelright=True,
                              left=False, right=True, pad=2)
        else:
            cb.ax.set_yticklabels([])
            cb.ax.tick_params(labelleft=False, labelright=False,
                              left=False, right=True, length=3)
        if i == 0:
            first_pos = cax.get_position()
        if i == len(channels) - 1:
            last_pos = cax.get_position()
    if first_pos and last_pos:
        label_x = (first_pos.x0 + last_pos.x1) / 2
        label_y = min(0.97, last_pos.y1 + 0.04)
        fig.text(label_x, label_y, "ΔE(B-V) [mag]", va="bottom", ha="center", fontsize=12)


# ----------------------------------------------------------------------------
# Main workflow.
# ----------------------------------------------------------------------------
def main():
    """Run Script 11 from validated inputs to the documented outputs."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Data-loading helper section.
    ellipse_df = load_xy_plane_ellipses(GEOM_CSV)
    arm_df = load_arms(ARM_CSV)
    print(f"Reading the dust cube for the plotting region: {DUST_PARQUET}")
    dust_df = load_dust_cube(DUST_PARQUET)
    print(f"  Number of voxels in the plotting region: {len(dust_df):,}")


    ell_maxopen = ellipses_for_method(ellipse_df, "maxopen")      # Compromise and maximum-opening ellipse conventions for representative OSB footprints.
    ell_compromise = ellipses_for_method(ellipse_df, "compromise")  # Compromise and maximum-opening ellipse conventions for representative OSB footprints.


    if ellipse_df is not None:
        out = ellipse_df.rename(columns={"a": "a_maxopen", "b": "b_maxopen"}).copy()
        out["a_compromise"] = out["a_maxopen"] * out["section_scale"]
        out["b_compromise"] = out["b_maxopen"] * out["section_scale"]
        out["area_pct"] = (out["section_scale"] ** 2) * 100.0   # Compromise and maximum-opening ellipse conventions for representative OSB footprints.
        try:
            out[["index", "mark", "shape", "x", "y", "angle",
                 "a_maxopen", "b_maxopen", "a_compromise", "b_compromise",
                 "section_scale", "area_pct"]].to_csv(
                OUT_ELLIPSE_CSV, index=False, encoding="utf-8-sig")
            print(f"Wrote plotting ellipse parameters: {OUT_ELLIPSE_CSV}")
        except PermissionError as exc:
            print(f"Warning: could not overwrite plotting ellipse parameters: {exc}")

        # Compromise and maximum-opening ellipse conventions for representative OSB footprints.
        ell = ellipse_df[ellipse_df["shape"].astype(str).str.lower() == "ellipsoid"]
        if len(ell):
            lin = ell["section_scale"].to_numpy(float)
            area = lin ** 2 * 100.0
            print(
                "Representative ellipsoid footprint scale: "
                f"min={np.min(lin):.3f}, median={np.median(lin):.3f}, max={np.max(lin):.3f}"
            )
            print(
                "Representative compromise area fraction: "
                f"min={np.min(area):.1f}%, median={np.median(area):.1f}%, max={np.max(area):.1f}%"
            )

    configure_matplotlib()

    # Figure styling and layout configuration for reproducible rendering.
    # Figure styling and layout configuration for reproducible rendering.
    fig_w, fig_h = 15.5, 15.0
    panel_h = 0.37
    panel_w = panel_h * fig_h / fig_w
    col1_x = 0.045
    col2_x = col1_x + panel_w + 0.075
    row_top_y = 0.555
    row_bot_y = 0.105
    fig = plt.figure(figsize=(fig_w, fig_h))

    ax_a = fig.add_axes([col1_x, row_top_y, panel_w, panel_h])   # Compromise and maximum-opening ellipse conventions for representative OSB footprints.
    ax_b = fig.add_axes([col2_x, row_top_y, panel_w, panel_h])   # Compromise and maximum-opening ellipse conventions for representative OSB footprints.
    ax_c = fig.add_axes([col1_x, row_bot_y, panel_w, panel_h])   # Compromise and maximum-opening ellipse conventions for representative OSB footprints.
    ax_d = fig.add_axes([col2_x, row_bot_y, panel_w, panel_h])   # Compromise and maximum-opening ellipse conventions for representative OSB footprints.

    draw_clean_map(ax_a, ell_maxopen, arm_df)
    draw_rgb_map(fig, ax_b, dust_df, ell_maxopen, arm_df)
    draw_clean_map(ax_c, ell_compromise, arm_df)
    draw_rgb_map(fig, ax_d, dust_df, ell_compromise, arm_df)

    # Figure-panel rendering section.
    for ax, tag in ((ax_a, "(a)"), (ax_b, "(b)"), (ax_c, "(c)"), (ax_d, "(d)")):
        ax.text(0.02, 0.98, tag, transform=ax.transAxes, fontsize=FONT_SIZES["panel_tag"],
                fontweight="bold", va="top", ha="left",
                bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                          edgecolor="none", alpha=0.85))

    # Figure-panel rendering section.
    for row_y, key, offset in ((row_top_y, "maxopen", 0.052), (row_bot_y, "compromise", 0.020)):
        label = ROW_LABELS[key]
        fig.text(col1_x, row_y + panel_h + offset, label,
                 ha="left", va="bottom", fontsize=18, fontweight="bold")

    add_shared_legend(fig, ellipse_mark_counts(ell_maxopen))

    fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved publication figure: {OUT_FIG}")

def _format_runtime(seconds: float) -> str:
    """Return a compact human-readable wall-clock runtime string."""
    seconds = max(0.0, float(seconds))
    hours, remainder = divmod(seconds, 3600.0)
    minutes, seconds = divmod(remainder, 60.0)
    if hours >= 1.0:
        return f"{int(hours)} h {int(minutes):02d} min {seconds:05.2f} s"
    if minutes >= 1.0:
        return f"{int(minutes)} min {seconds:05.2f} s"
    return f"{seconds:.2f} s"


def _run_with_timing(entrypoint, script_file: str) -> None:
    """Run a script entry point and print its total wall-clock runtime."""
    script_name = Path(script_file).name
    start = time.perf_counter()
    try:
        entrypoint()
    except SystemExit as exc:
        elapsed = time.perf_counter() - start
        if exc.code in (None, 0):
            print(f"[runtime] {script_name} completed in: {_format_runtime(elapsed)} ({elapsed:.2f} s)", flush=True)
        else:
            print(f"[runtime] {script_name} elapsed time before failure: {_format_runtime(elapsed)} ({elapsed:.2f} s)", flush=True)
        raise
    except BaseException:
        elapsed = time.perf_counter() - start
        print(f"[runtime] {script_name} elapsed time before failure: {_format_runtime(elapsed)} ({elapsed:.2f} s)", flush=True)
        raise
    elapsed = time.perf_counter() - start
    print(f"[runtime] {script_name} completed in: {_format_runtime(elapsed)} ({elapsed:.2f} s)", flush=True)


if __name__ == "__main__":
    _run_with_timing(main, __file__)
