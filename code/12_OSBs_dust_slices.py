#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script 12: OSBs dust slices

Purpose
-------
Draw three vertical dust slices with open-superbubble ellipses, bubbles, and HMSFRs.

Method overview
---------------
1. Load the fitted shell geometry, raw dust cube, spiral arms, Grade A/B bubbles and
   HMSFR catalogue.
2. Construct the three configured vertical dust slices and evaluate catalogue contacts
   with the adopted thick-shell geometry.
3. Draw the slice maps with shell intersections and tracer annotations, and save their
   combined figure and contact information.

Main inputs
-----------
- ../results/superbubble_final_fit_parameters.csv
- ../results/intermediate_output/3d_dust_map_products/raw_3d_dust_cube.parquet
- ../data/Reid2019_spirals_xy.csv
- ../data/Bubbles.csv
- ../data/star_cluster_data/Reid2019_HMSFR.csv

Main outputs
------------
- ../results/figures/12_three_dust_slices.png; main-text Fig. 2.

Figure/table role
-----------------
../results/figures/12_three_dust_slices.png; main-text Fig. 2.

Runtime and data notes
----------------------
Reference runtime: 11.47 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

External-data caveat
--------------------
This script uses packaged intermediate or final inputs unless the inputs section explicitly names an external cache or survey product.

Reading the code
----------------
Start with the path and scientific settings below, then follow main() at
the end of the file. The method overview above describes the order of the
analysis steps; the input/output lists identify the upstream data and
products written by this stage. Relative data paths are resolved from code/.
"""

from __future__ import annotations
import time

import argparse
from collections import Counter
import os
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import numpy.ma as ma
import pandas as pd
import pyarrow.parquet as pq
from matplotlib.colors import Normalize
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse
from scipy.ndimage import gaussian_filter


HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "12_three_dust_slices"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"

FINAL_TABLE = FINAL_DIR / "superbubble_final_fit_parameters.csv"
DUST_PARQUET = Path("..") / "results" / "intermediate_output" / "3d_dust_map_products" / "raw_3d_dust_cube.parquet"
ARM_CSV = DATA_DIR / "Reid2019_spirals_xy.csv"
BUBBLE_CSV = DATA_DIR / "Bubbles.csv"
HMSFR_CSV = DATA_DIR / "star_cluster_data" / "Reid2019_HMSFR.csv"

OUT_ELLIPSE_CSV = OUT_DIR / "12_plotting_ellipse_parameters.csv"
OUT_FIG = FINAL_FIG_DIR / "12_three_dust_slices.png"

PLOT_HALF_KPC = 3.5
Z_SLICES = [(-0.20, -0.05), (-0.05, 0.05), (0.05, 0.20)]
SMOOTH_SIGMA = 1.2
DUST_VMIN = 0.0
B_DUST_VMAX = 0.06
CD_DUST_VMAX = 0.04
B_DUST_CBAR_TICKS = np.arange(0.0, 0.061, 0.02)
CD_DUST_CBAR_TICKS = np.arange(0.0, 0.041, 0.01)
SUN_GC_DISTANCE_KPC = 8.12
ELLIPSOID_CAP_SECTION_FRAC = 0.5
SHELL_U_MIN = 0.9
SHELL_U_MAX = 1.2
CYLINDER_HALF_HEIGHT_KPC = 0.025
DEFAULT_HMSFR_DIAMETER_PC = 200.0

COLOR_MAP = {0: "black", 1: "blue", 2: "red", 3: "green"}
LABEL_MAP = {0: "Closed", 1: "North Open", 2: "South Open", 3: "Disk-penetrating"}
ARM_STYLES = {
    "Perseus": ("r", "-", 2.0),
    "Local": ("magenta", "-", 2.0),
    "Sagittarius": ("lime", "-", 2.0),
    "Carina": ("cyan", "-", 2.0),
}
CATALOG_YELLOW = "#FFD500"
HMSFR_COLOR = "darkviolet"
RAD_LINE_ANGLE_DEG = 60.0
RAD_LINE_INTERCEPT_KPC = 0.6
RAD_LINE_HALF_WIDTH_KPC = 0.1
TOPVIEW_FONT_SIZES = {
    "axis_label": 14,
    "tick_label": 11,
    "legend_text": 11,
    "legend_title": 12,
    "annotation": 11,
}
ARM_DISPLAY = {
    "Perseus": "Perseus arm",
    "Local": "Local arm",
    "Sagittarius": "Sagittarius arm",
    "Carina": "Carina arm",
}
TOPVIEW_ARM_STYLES = {
    "Perseus": ("r", "-", 2.0),
    "Local": ("magenta", "-", 2.0),
    "Sagittarius": ("lime", "-", 2.0),
    "Carina": ("cyan", "-", 2.0),
    "LoS": ("black", "--", 1.5),
}
TOPVIEW_ARM_DISPLAY = {
    "Perseus": "Perseus arm",
    "Local": "Local arm",
    "Sagittarius": "Sagittarius arm",
    "Carina": "Carina arm",
    "LoS": "Local arm spur",
}


def ensure_existing_file(path: Path, label: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    return path


def configure_matplotlib():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.unicode_minus": False,
            "figure.dpi": 140,
            "savefig.dpi": 300,
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )


def midcap_section_scale(final_df: pd.DataFrame) -> np.ndarray:
    shape = final_df["shape"].astype(str).str.lower().to_numpy()
    cz = final_df["center_z_kpc"].to_numpy(float)
    c = final_df["c_radius_kpc"].to_numpy(float)
    z_base = final_df["xy_plane_z_kpc"].to_numpy(float)
    mark = final_df["mark"].to_numpy(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        k_base = np.sqrt(np.clip(1.0 - ((z_base - cz) / c) ** 2, 0.0, None))
        z_apex = np.where(mark == 1, cz - c, cz + c)
        z_sec = z_base + ELLIPSOID_CAP_SECTION_FRAC * (z_apex - z_base)
        k_sec = np.sqrt(np.clip(1.0 - ((z_sec - cz) / c) ** 2, 0.0, None))
        scale = np.where(k_base > 0, k_sec / k_base, 1.0)
    return np.where(shape == "ellipsoid", scale, 1.0)


def load_latest_ellipses(final_table: Path, method: str) -> pd.DataFrame:
    df = pd.read_csv(ensure_existing_file(final_table, "final parameter table"), encoding="utf-8-sig")
    required = [
        "id", "mark", "shape",
        "xy_plane_center_x_kpc", "xy_plane_center_y_kpc",
        "xy_plane_a_kpc", "xy_plane_b_kpc", "xy_plane_angle_deg",
        "center_z_kpc", "c_radius_kpc", "xy_plane_z_kpc",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"final parameter table is missing required columns: {', '.join(missing)}")

    scale = midcap_section_scale(df) if method == "compromise" else np.ones(len(df), dtype=float)
    out = pd.DataFrame(
        {
            "index": df["id"].astype(int),
            "mark": df["mark"].astype(int),
            "shape": df["shape"].astype(str),
            "x": df["xy_plane_center_x_kpc"].astype(float),
            "y": df["xy_plane_center_y_kpc"].astype(float),
            "a": df["xy_plane_a_kpc"].astype(float) * scale,
            "b": df["xy_plane_b_kpc"].astype(float) * scale,
            "angle": df["xy_plane_angle_deg"].astype(float),
            "ellipse_method": method,
            "section_scale": scale,
            "sort_z": df["center_z_kpc"].astype(float),
        }
    )
    return out.sort_values("index").reset_index(drop=True)


def rotate_points_to_local(points: np.ndarray, angle_deg: float) -> np.ndarray:
    theta = np.deg2rad(-float(angle_deg))
    c, s = np.cos(theta), np.sin(theta)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)
    return np.asarray(points, dtype=float) @ rot.T


def geometry_for_shell(row: pd.Series):
    shape = str(row["shape"]).lower()
    center = np.array(
        [float(row["center_x_kpc"]), float(row["center_y_kpc"]), float(row["center_z_kpc"])],
        dtype=float,
    )
    if shape == "cylinder":
        axes = np.array(
            [float(row["a_radius_kpc"]), float(row["b_radius_kpc"]), CYLINDER_HALF_HEIGHT_KPC],
            dtype=float,
        )
        return "cylinder", int(row["mark"]), center, axes
    axes = np.array(
        [float(row["a_radius_kpc"]), float(row["b_radius_kpc"]), max(float(row["c_radius_kpc"]), 1e-9)],
        dtype=float,
    )
    return "ellipsoid_cap", int(row["mark"]), center, axes


def extreme_norm_over_ellipsoid(xi_c: np.ndarray, e: np.ndarray, want: str) -> np.ndarray:
    """Metric-coordinate min/max norm over an axis-aligned solid ellipsoid."""
    xi_c = np.asarray(xi_c, dtype=float)
    e = np.maximum(np.asarray(e, dtype=float), 1e-12)
    e2 = e ** 2
    ce = xi_c * e
    maxe2 = np.max(e2, axis=1)
    mine2 = np.min(e2, axis=1)
    uc = np.sqrt(np.sum(xi_c ** 2, axis=1))
    span = (uc + np.max(e, axis=1) + 1.0) ** 2

    def g(mu):
        return np.sum((ce / (mu[:, None] - e2)) ** 2, axis=1)

    def norm_at(mu):
        return np.sqrt(np.sum((mu[:, None] * xi_c / (mu[:, None] - e2)) ** 2, axis=1))

    if want == "max":
        lo = maxe2 + 1e-12 * (maxe2 + 1.0) + 1e-15
        hi = maxe2 + span + 1.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            gt = g(mid) > 1.0
            lo = np.where(gt, mid, lo)
            hi = np.where(gt, hi, mid)
        return norm_at(0.5 * (lo + hi))

    inside = np.sum((xi_c / e) ** 2, axis=1) <= 1.0
    hi = mine2 - 1e-12 * (mine2 + 1.0) - 1e-15
    lo = mine2 - span - 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        gt = g(mid) > 1.0
        hi = np.where(gt, mid, hi)
        lo = np.where(gt, lo, mid)
    return np.where(inside, 0.0, norm_at(0.5 * (lo + hi)))


def thick_shell_intersection_mask(points: np.ndarray, row: pd.Series, radii: np.ndarray | float) -> np.ndarray:
    """True when a tracer sphere intersects the finite-thickness shell u in [0.9, 1.2]."""
    kind, mark, center, axes = geometry_for_shell(row)
    local = rotate_points_to_local(np.asarray(points, dtype=float) - center[None, :], float(row["angle_deg"]))
    X, Y, Z = local[:, 0], local[:, 1], local[:, 2]
    a = max(float(axes[0]), 1e-9)
    b = max(float(axes[1]), 1e-9)
    c = max(float(axes[2]), 1e-9)
    r = np.broadcast_to(np.asarray(radii, dtype=float), X.shape).astype(float)

    if kind == "cylinder":
        h = CYLINDER_HALF_HEIGHT_KPC
        z_star = np.clip(Z, -h, h)
        rho2 = r ** 2 - (z_star - Z) ** 2
        reach = rho2 >= 0.0
        rho = np.sqrt(np.maximum(rho2, 0.0))
        xi_c = np.stack([X / a, Y / b], axis=1)
        e = np.stack([rho / a, rho / b], axis=1)
        u_min = extreme_norm_over_ellipsoid(xi_c, e, "min")
        u_max = extreme_norm_over_ellipsoid(xi_c, e, "max")
        return reach & (u_min <= SHELL_U_MAX) & (u_max >= SHELL_U_MIN)

    xi_c = np.stack([X / a, Y / b, Z / c], axis=1)
    e = np.stack([r / a, r / b, r / c], axis=1)
    u_min = extreme_norm_over_ellipsoid(xi_c, e, "min")
    u_max = extreme_norm_over_ellipsoid(xi_c, e, "max")
    cap_reach = (Z >= -r) if mark == 2 else (Z <= r)
    return cap_reach & (u_min <= SHELL_U_MAX) & (u_max >= SHELL_U_MIN)


def observed_thick_shell_union_mask(
    sample_df: pd.DataFrame,
    catalog_df: pd.DataFrame,
    radii: np.ndarray | float,
) -> tuple[np.ndarray, list[list[int]]]:
    points = catalog_df[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(dtype=float)
    union_mask = np.zeros(len(catalog_df), dtype=bool)
    containing_ids: list[list[int]] = [[] for _ in range(len(catalog_df))]
    for row in sample_df.itertuples(index=False):
        row_series = pd.Series(row._asdict())
        mask = thick_shell_intersection_mask(points, row_series, radii)
        union_mask |= mask
        sb_id = int(row_series["id"])
        for idx in np.flatnonzero(mask):
            containing_ids[int(idx)].append(sb_id)
    return union_mask, containing_ids


def annotate_shell_contacts(
    sample_df: pd.DataFrame,
    bubble_df: pd.DataFrame,
    hmsfr_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    bubble_df = bubble_df.copy()
    hmsfr_df = hmsfr_df.copy()
    bubble_radii = bubble_df["radius_kpc"].to_numpy(dtype=float)
    hmsfr_radius_kpc = DEFAULT_HMSFR_DIAMETER_PC / 2.0 / 1000.0

    bubble_mask, bubble_sb_ids = observed_thick_shell_union_mask(sample_df, bubble_df, bubble_radii)
    hmsfr_mask, hmsfr_sb_ids = observed_thick_shell_union_mask(sample_df, hmsfr_df, hmsfr_radius_kpc)

    bubble_df["in_any_sb_shell"] = bubble_mask
    hmsfr_df["in_any_sb_shell"] = hmsfr_mask
    bubble_df["containing_sb_ids"] = [";".join(f"SB{int(x)}" for x in ids) for ids in bubble_sb_ids]
    hmsfr_df["containing_sb_ids"] = [";".join(f"SB{int(x)}" for x in ids) for ids in hmsfr_sb_ids]
    return bubble_df, hmsfr_df


def load_dust_cube(parquet_path: Path) -> pd.DataFrame:
    ensure_existing_file(parquet_path, "smoothed 3D dust cube parquet")
    half = PLOT_HALF_KPC
    z_lo = min(lo for lo, _ in Z_SLICES)
    z_hi = max(hi for _, hi in Z_SLICES)
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


def spherical_to_cartesian(distance_kpc, lon_deg, lat_deg) -> np.ndarray:
    lon = np.deg2rad(np.asarray(lon_deg, dtype=float))
    lat = np.deg2rad(np.asarray(lat_deg, dtype=float))
    dist = np.asarray(distance_kpc, dtype=float)
    cos_lat = np.cos(lat)
    return np.column_stack(
        [
            dist * cos_lat * np.cos(lon),
            dist * cos_lat * np.sin(lon),
            dist * np.sin(lat),
        ]
    )


def load_ab_bubbles(path: Path) -> pd.DataFrame:
    df = pd.read_csv(ensure_existing_file(path, "small-bubble table"))
    required = ["ID", "l", "b", "physical_size", "best_distance", "Grade"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"{path}  is missing required columns: {', '.join(missing)}")
    df = df[df["Grade"].astype(str).str.strip().str.upper().isin(["A", "B"])].copy()
    df = df.drop_duplicates(subset=["ID", "l", "b", "physical_size", "best_distance"]).copy()
    for col in ["l", "b", "physical_size", "best_distance"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[
        np.isfinite(df["l"])
        & np.isfinite(df["b"])
        & np.isfinite(df["best_distance"])
        & np.isfinite(df["physical_size"])
    ].copy()
    xyz = spherical_to_cartesian(df["best_distance"], df["l"], df["b"])
    df["x_kpc"] = xyz[:, 0]
    df["y_kpc"] = xyz[:, 1]
    df["z_kpc"] = xyz[:, 2]
    df["radius_kpc"] = df["physical_size"] / 2000.0
    return df.reset_index(drop=True)


def load_hmsfr(path: Path) -> pd.DataFrame:
    df = pd.read_csv(ensure_existing_file(path, "HMSFR table")).copy()
    required = ["Name", "x_kpc", "y_kpc", "z_kpc"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"{path}  is missing required columns: {', '.join(missing)}")
    for col in ["x_kpc", "y_kpc", "z_kpc"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[np.isfinite(df["x_kpc"]) & np.isfinite(df["y_kpc"]) & np.isfinite(df["z_kpc"])].copy()
    return df.reset_index(drop=True)


def extract_xy_grid(
    df: pd.DataFrame,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    z_range: tuple[float, float],
    sigma: float,
) -> ma.MaskedArray | None:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    mask = (
        (df["x"] >= x_range[0]) & (df["x"] <= x_range[1])
        & (df["y"] >= y_range[0]) & (df["y"] <= y_range[1])
        & (df["z"] >= z_range[0]) & (df["z"] <= z_range[1])
    )
    df_cut = df.loc[mask].copy()
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
    layer_thickness_kpc = float(z_range[1] - z_range[0])
    g = df_cut.groupby(["x", "y"])["dust"].mean().reset_index()
    if g.empty:
        return None
    g["delta_ebv"] = g["dust"] * layer_thickness_kpc
    grid = g.pivot(index="y", columns="x", values="delta_ebv")
    grid = grid.sort_index(ascending=True).sort_index(axis=1, ascending=True)
    img_array = gaussian_filter(grid.fillna(0).values, sigma=sigma)
    return ma.masked_array(img_array, mask=grid.isnull().values)


def load_arms(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        print(f"Optional spiral-arm file not found: {path}")
        return None
    return pd.read_csv(path)


def draw_arms(ax, arm_df: pd.DataFrame | None):
    if arm_df is None:
        return
    for name, (color, linestyle, linewidth) in ARM_STYLES.items():
        mask = arm_df["arm"] == name
        if not mask.any():
            continue
        ax.plot(
            -(arm_df.loc[mask, "yy"]) + SUN_GC_DISTANCE_KPC,
            arm_df.loc[mask, "xx"],
            linestyle=linestyle, color=color, linewidth=linewidth, alpha=0.8,
        )
        ax.plot(
            -(arm_df.loc[mask, "yy0"]) + SUN_GC_DISTANCE_KPC,
            arm_df.loc[mask, "xx0"],
            linestyle=":", color=color, linewidth=1.2, alpha=0.6,
        )
        ax.plot(
            -(arm_df.loc[mask, "yy1"]) + SUN_GC_DISTANCE_KPC,
            arm_df.loc[mask, "xx1"],
            linestyle=":", color=color, linewidth=1.2, alpha=0.6,
        )


def draw_reference_grid(ax):
    guide_kwargs = {"edgecolor": "k", "linestyle": "--", "linewidth": 1.0, "alpha": 0.7}
    line_kwargs = {"color": "k", "linestyle": "--", "linewidth": 1.0, "alpha": 0.7}
    for radius in (1, 2, 3, 4):
        ax.add_patch(Circle((0, 0), radius, facecolor="none", **guide_kwargs))
    for angle_deg in np.arange(0, 360, 45):
        theta = np.deg2rad(angle_deg)
        ax.plot([0, 5 * np.cos(theta)], [0, 5 * np.sin(theta)], **line_kwargs)
    for angle_deg in (0, 90, 180, 270):
        theta = np.deg2rad(angle_deg)
        ax.text(
            3.2 * np.cos(theta), 3.2 * np.sin(theta), f"{angle_deg}°",
            color="white", fontsize=10, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.2", fc="black", ec="none", alpha=0.6),
        )


def draw_ellipses(ax, ellipse_df: pd.DataFrame, count_marks: bool, label_ids: bool = True) -> Counter:
    counts: Counter = Counter()
    draw_df = ellipse_df.copy()
    if "sort_z" in draw_df.columns:
        draw_df = draw_df.sort_values("sort_z", ascending=True)
    for _, row in draw_df.iterrows():
        mark = int(row["mark"])
        color = COLOR_MAP.get(mark, "gray")
        if mark not in COLOR_MAP:
            continue
        ax.add_patch(
            Ellipse(
                (float(row["x"]), float(row["y"])),
                2.0 * float(row["a"]),
                2.0 * float(row["b"]),
                angle=float(row["angle"]),
                edgecolor=color,
                facecolor="none",
                linewidth=1.8,
                alpha=0.85,
            )
        )
        if label_ids:
            ax.text(
                float(row["x"]), float(row["y"]), str(int(row["index"])),
                color=color, fontsize=9, ha="center", va="center", fontweight="bold",
            )
        if count_marks:
            counts[mark] += 1
    return counts


def draw_catalog_overlays(
    ax,
    z_range: tuple[float, float],
    bubble_df: pd.DataFrame,
    hmsfr_df: pd.DataFrame,
    add_label: bool,
):
    z_lo, z_hi = z_range
    bubble_slice = (bubble_df["z_kpc"] >= z_lo) & (bubble_df["z_kpc"] <= z_hi)
    hmsfr_slice = (hmsfr_df["z_kpc"] >= z_lo) & (hmsfr_df["z_kpc"] <= z_hi)
    bubble_shell = bubble_df["in_any_sb_shell"].astype(bool)
    hmsfr_shell = hmsfr_df["in_any_sb_shell"].astype(bool)
    bubble_in_mask = bubble_slice & bubble_shell
    bubble_out_mask = bubble_slice & ~bubble_shell
    hmsfr_in_mask = hmsfr_slice & hmsfr_shell
    hmsfr_out_mask = hmsfr_slice & ~hmsfr_shell

    hollow_bubble_layers = [
        ("black", 1.3, 6.60, 0.72),
        (CATALOG_YELLOW, 0.8, 6.80, 0.88),
    ]
    for _, row in bubble_df.loc[bubble_out_mask].iterrows():
        center = (float(row["x_kpc"]), float(row["y_kpc"]))
        radius = float(row["radius_kpc"])
        for color, linewidth, zorder, alpha in hollow_bubble_layers:
            ax.add_patch(
                Circle(
                    center,
                    radius,
                    facecolor="none",
                    edgecolor=color,
                    linewidth=linewidth,
                    alpha=alpha,
                    zorder=zorder,
                )
            )

    filled_bubble_layers = [
        (CATALOG_YELLOW, "black", 0.55, 7.80),
    ]
    for _, row in bubble_df.loc[bubble_in_mask].iterrows():
        center = (float(row["x_kpc"]), float(row["y_kpc"]))
        radius = float(row["radius_kpc"])
        for facecolor, edgecolor, linewidth, zorder in filled_bubble_layers:
            ax.add_patch(
                Circle(
                    center,
                    radius,
                    facecolor=facecolor,
                    edgecolor=edgecolor,
                    linewidth=linewidth,
                    alpha=0.78,
                    zorder=zorder,
                )
            )

    if hmsfr_out_mask.any():
        hmsfr_x = hmsfr_df.loc[hmsfr_out_mask, "x_kpc"]
        hmsfr_y = hmsfr_df.loc[hmsfr_out_mask, "y_kpc"]
        for color, size, linewidth, zorder in [
            ("black", 120, 1.15, 6.90),
            (HMSFR_COLOR, 100, 0.9, 7.10),
        ]:
            ax.scatter(
                hmsfr_x,
                hmsfr_y,
                marker="*",
                s=size,
                facecolors="none",
                edgecolors=color,
                linewidths=linewidth,
                alpha=0.90,
                zorder=zorder,
            )

    if hmsfr_in_mask.any():
        hmsfr_x = hmsfr_df.loc[hmsfr_in_mask, "x_kpc"]
        hmsfr_y = hmsfr_df.loc[hmsfr_in_mask, "y_kpc"]
        for color, size, zorder in [
            ("black", 120, 8.80),
            (HMSFR_COLOR, 100, 9.00),
        ]:
            ax.scatter(
                hmsfr_x,
                hmsfr_y,
                marker="*",
                s=size,
                facecolors=color,
                edgecolors=color,
                linewidths=0,
                alpha=0.98,
                zorder=zorder,
            )

    if add_label:
        bubble_in_handle = Line2D(
            [0], [0], marker="o", linestyle="none", markersize=11,
            markerfacecolor=CATALOG_YELLOW, markeredgecolor="black", markeredgewidth=0.65,
        )
        bubble_out_handle = Line2D(
            [0], [0], marker="o", linestyle="none", markersize=11,
            markerfacecolor="none", markeredgecolor=CATALOG_YELLOW, markeredgewidth=1.8,
        )
        bubble_out_handle.set_path_effects([
            pe.Stroke(linewidth=2.5, foreground="black"),
            pe.Normal(),
        ])
        hmsfr_in_handle = Line2D(
            [0], [0], marker="*", linestyle="none", markersize=12,
            markerfacecolor=HMSFR_COLOR, markeredgecolor=HMSFR_COLOR, label="HMSFR",
        )
        hmsfr_in_handle.set_path_effects([
            pe.Stroke(linewidth=2.0, foreground="black"),
            pe.Normal(),
        ])
        hmsfr_out_handle = Line2D(
            [0], [0], marker="*", linestyle="none", markersize=12,
            markerfacecolor="none", markeredgecolor=HMSFR_COLOR, markeredgewidth=1.2,
        )
        hmsfr_out_handle.set_path_effects([
            pe.Stroke(linewidth=2.0, foreground="black"),
            pe.Normal(),
        ])
        handles = [bubble_in_handle, bubble_out_handle, hmsfr_in_handle, hmsfr_out_handle]
        labels = [
            "Shell-contact bubbles",
            "Non-shell bubbles",
            "Shell-contact high-mass star-forming regions",
            "Non-shell high-mass star-forming regions",
        ]
        legend = ax.legend(
            handles=handles,
            labels=labels,
            loc="lower center",
            bbox_to_anchor=(0.52, 0.015),
            ncol=2,
            fontsize=9.5,
            frameon=True,
            framealpha=0.82,
            borderaxespad=0.25,
            columnspacing=0.8,
            handletextpad=0.45,
        )
        legend.set_zorder(30)


def draw_topview_clean_panel(
    ax,
    ellipse_df: pd.DataFrame,
    arm_df: pd.DataFrame | None,
):
    """Render or save a documented figure product without changing upstream measurements."""
    ax.set_facecolor("whitesmoke")

    circle_kwargs = {"edgecolor": "gray", "linewidth": 1.0, "alpha": 0.6}
    line_kwargs = {"color": "gray", "linestyle": ":", "linewidth": 0.8, "alpha": 0.6}
    for radius in (1, 2, 3):
        ax.add_patch(Circle((0, 0), radius, facecolor="none", **circle_kwargs))
        ax.text(0, radius, f"{radius} kpc", color="gray",
                fontsize=TOPVIEW_FONT_SIZES["annotation"] - 2, ha="center", va="bottom")
    for angle in np.arange(0, 360, 30):
        theta = np.deg2rad(angle)
        ax.plot([0, 4.5 * np.cos(theta)], [0, 4.5 * np.sin(theta)], **line_kwargs)
    for angle in (0, 180, 270):
        theta = np.deg2rad(angle)
        ax.text(3.1 * np.cos(theta), 3.1 * np.sin(theta), f"{angle}°",
                color="black", fontsize=TOPVIEW_FONT_SIZES["annotation"],
                ha="center", va="center")

    draw_topview_arms(ax, arm_df)

    counts: Counter = Counter()
    draw_df = ellipse_df.copy()
    if "sort_z" in draw_df.columns:
        draw_df = draw_df.sort_values("sort_z", ascending=True)
    for _, row in draw_df.iterrows():
        mark = int(row["mark"])
        color = COLOR_MAP.get(mark, "gray")
        if mark not in COLOR_MAP:
            continue
        counts[mark] += 1
        ax.add_patch(
            Ellipse(
                (float(row["x"]), float(row["y"])),
                2.0 * float(row["a"]),
                2.0 * float(row["b"]),
                angle=float(row["angle"]),
                facecolor=color,
                edgecolor="black",
                alpha=0.6,
            )
        )
        text_color = "white" if color in ("black", "blue") else "black"
        ax.text(
            float(row["x"]), float(row["y"]), str(int(row["index"])),
            color=text_color, fontsize=8, ha="center", va="center", fontweight="bold",
            path_effects=[pe.withStroke(linewidth=1.5,
                                        foreground="black" if text_color == "white" else "white")],
        )

    x_line = np.linspace(-4, 4, 100)
    slope = np.tan(np.deg2rad(RAD_LINE_ANGLE_DEG))
    y_line = slope * x_line + RAD_LINE_INTERCEPT_KPC
    offset = RAD_LINE_HALF_WIDTH_KPC * np.sqrt(1 + slope ** 2)
    ax.fill_between(x_line, y_line - offset, y_line + offset, color="gray", alpha=0.3)
    ax.plot(x_line, y_line, color="black", linewidth=2.0)

    ax.set_xlim(-PLOT_HALF_KPC, PLOT_HALF_KPC)
    ax.set_ylim(-PLOT_HALF_KPC, PLOT_HALF_KPC)
    ax.set_aspect("equal")
    ax.set_xlabel("X [kpc]", fontsize=TOPVIEW_FONT_SIZES["axis_label"])
    ax.set_ylabel("Y [kpc]", fontsize=TOPVIEW_FONT_SIZES["axis_label"])
    ax.tick_params(labelsize=TOPVIEW_FONT_SIZES["tick_label"])
    ax.grid(True, linestyle="--", alpha=0.4)
    add_topview_legends(ax, counts)


def draw_topview_arms(ax, arm_df: pd.DataFrame | None):
    if arm_df is None:
        return
    for name, (color, linestyle, linewidth) in TOPVIEW_ARM_STYLES.items():
        mask = arm_df["arm"] == name
        if not mask.any():
            continue
        ax.plot(
            -(arm_df.loc[mask, "yy"]) + SUN_GC_DISTANCE_KPC,
            arm_df.loc[mask, "xx"],
            linestyle=linestyle, color=color, linewidth=linewidth, alpha=0.7,
        )
        ax.plot(
            -(arm_df.loc[mask, "yy0"]) + SUN_GC_DISTANCE_KPC,
            arm_df.loc[mask, "xx0"],
            linestyle=":", color=color, linewidth=1.2, alpha=0.6,
        )
        ax.plot(
            -(arm_df.loc[mask, "yy1"]) + SUN_GC_DISTANCE_KPC,
            arm_df.loc[mask, "xx1"],
            linestyle=":", color=color, linewidth=1.2, alpha=0.6,
        )


def add_topview_legends(ax, counts: Counter):
    leg1 = ax.legend(
        [Line2D([0], [0], color="black", linestyle="-", linewidth=2)],
        ["Line of Radcliffe wave\n$y = \\tan(60^\\circ)\\,x + 0.6$"],
        loc="upper left", bbox_to_anchor=(0, 1),
        fontsize=TOPVIEW_FONT_SIZES["legend_text"], frameon=True,
        facecolor="white", framealpha=0.9,
    )
    ax.add_artist(leg1)

    total = int(sum(counts.values()))
    ellipse_handles, ellipse_labels = [], []
    for mark in (1, 2, 3):
        ellipse_handles.append(
            Line2D([0], [0], marker="o", linestyle="None",
                   markerfacecolor=COLOR_MAP[mark], markeredgecolor="black", markersize=10)
        )
        ellipse_labels.append(f"{LABEL_MAP[mark]} ({int(counts[mark]) if mark in counts else 0})")
    leg3 = ax.legend(
        ellipse_handles, ellipse_labels, loc="lower left", bbox_to_anchor=(0, 0),
        title=f"Total: {total}", fontsize=TOPVIEW_FONT_SIZES["legend_text"],
        title_fontsize=TOPVIEW_FONT_SIZES["legend_title"], frameon=True,
        facecolor="white", framealpha=0.9,
    )
    ax.add_artist(leg3)

    arm_handles = [
        Line2D([0], [0], color=color, linestyle=linestyle, linewidth=2)
        for color, linestyle, _ in TOPVIEW_ARM_STYLES.values()
    ]
    arm_labels = [TOPVIEW_ARM_DISPLAY[name] for name in TOPVIEW_ARM_STYLES]
    leg4 = ax.legend(
        arm_handles, arm_labels, loc="lower right", bbox_to_anchor=(1, 0),
        title="Spiral Arms", fontsize=TOPVIEW_FONT_SIZES["legend_text"] - 1,
        title_fontsize=TOPVIEW_FONT_SIZES["legend_title"] - 1, frameon=True,
        facecolor="white", framealpha=0.9,
    )
    ax.add_artist(leg4)


def save_figure(
    dust_df: pd.DataFrame,
    ellipse_df: pd.DataFrame,
    arm_df: pd.DataFrame | None,
    bubble_df: pd.DataFrame,
    hmsfr_df: pd.DataFrame,
    output: Path,
):
    x_range = (-PLOT_HALF_KPC, PLOT_HALF_KPC)
    y_range = (-PLOT_HALF_KPC, PLOT_HALF_KPC)
    x_range, y_range = clip_xy_range_to_dust_data(dust_df, x_range, y_range)
    b_norm = Normalize(vmin=DUST_VMIN, vmax=B_DUST_VMAX)
    cd_norm = Normalize(vmin=DUST_VMIN, vmax=CD_DUST_VMAX)
    cmap = "Greys"

    imgs: list[ma.MaskedArray] = []
    for z_range in Z_SLICES:
        img = extract_xy_grid(dust_df, x_range, y_range, z_range, SMOOTH_SIGMA)
        if img is None:
            raise ValueError('Invalid input or missing required data.')
        imgs.append(img)

    fig = plt.figure(figsize=(14.0, 12.0))
    gs = fig.add_gridspec(
        2, 3,
        width_ratios=[1.0, 1.0, 0.04],
        left=0.055, right=0.935, bottom=0.065, top=0.945,
        wspace=0.08, hspace=0.20,
    )
    ax_topview = fig.add_subplot(gs[0, 0])
    axes = [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    cax_b = fig.add_subplot(gs[0, 2])
    cax_cd = fig.add_subplot(gs[1, 2])

    draw_topview_clean_panel(ax_topview, ellipse_df, arm_df)
    z_panel_order = [1, 0, 2]
    for panel_i, (ax, z_index) in enumerate(zip(axes, z_panel_order)):
        panel_norm = b_norm if z_index == 1 else cd_norm
        ax.imshow(imgs[z_index], origin="lower", cmap=cmap, norm=panel_norm, extent=[*x_range, *y_range])
        ax.set_title(f"Z: [{Z_SLICES[z_index][0]:.2f}, {Z_SLICES[z_index][1]:.2f}] kpc", fontsize=16)
        ax.set_xlabel("X [kpc]", fontsize=14)
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.tick_params(labelsize=12)
        ax.set_aspect("equal", adjustable="box")
        draw_reference_grid(ax)
        draw_ellipses(ax, ellipse_df, count_marks=False, label_ids=False)
        draw_catalog_overlays(ax, Z_SLICES[z_index], bubble_df, hmsfr_df, add_label=(z_index == 1))
        ax.set_xlim(x_range)
        ax.set_ylim(y_range)
    for ax, label in zip([ax_topview, *axes], ["(a)", "(b)", "(c)", "(d)"]):
        ax.text(0.0, 1.02, label, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=15, fontweight="bold")
    axes[1].set_ylabel("Y [kpc]", fontsize=14)

    cbar_b = fig.colorbar(
        plt.cm.ScalarMappable(norm=b_norm, cmap=cmap),
        cax=cax_b, orientation="vertical",
        ticks=B_DUST_CBAR_TICKS,
    )
    cbar_b.set_label(r"$\Delta E(B-V)$ [mag]", fontsize=14)
    cbar_b.ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    cbar_b.ax.tick_params(labelsize=12)

    cbar_cd = fig.colorbar(
        plt.cm.ScalarMappable(norm=cd_norm, cmap=cmap),
        cax=cax_cd, orientation="vertical",
        ticks=CD_DUST_CBAR_TICKS,
    )
    cbar_cd.set_label(r"$\Delta E(B-V)$ [mag]", fontsize=14)
    cbar_cd.ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    cbar_cd.ax.tick_params(labelsize=12)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line interface for Script 12."""
    parser = argparse.ArgumentParser(description='Run Script 12: OSBs dust slices.')
    parser.add_argument(
        "--ellipse-method",
        choices=["maxopen", "compromise"],
        default="compromise",
        help='Command-line option for the documented workflow.',
    )
    return parser


def main():
    """Run Script 12 from validated inputs to the documented outputs."""
    args = build_parser().parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)
    configure_matplotlib()

    sample_df = pd.read_csv(ensure_existing_file(FINAL_TABLE, "final parameter table"), encoding="utf-8-sig")
    sample_df = sample_df.sort_values("id").reset_index(drop=True)
    ellipse_df = load_latest_ellipses(FINAL_TABLE, args.ellipse_method)
    ellipse_df.to_csv(OUT_ELLIPSE_CSV, index=False, encoding="utf-8-sig")
    print(f"Wroteplotting_ellipse_parameters: {OUT_ELLIPSE_CSV}")

    print(f"Reading dust plotting region: {DUST_PARQUET}")
    dust_df = load_dust_cube(DUST_PARQUET)
    print(f"  Dust voxel count: {len(dust_df):,}")
    arm_df = load_arms(ARM_CSV)
    bubble_df = load_ab_bubbles(BUBBLE_CSV)
    hmsfr_df = load_hmsfr(HMSFR_CSV)
    bubble_df, hmsfr_df = annotate_shell_contacts(sample_df, bubble_df, hmsfr_df)
    n_bubble_shell = int(bubble_df["in_any_sb_shell"].sum())
    n_hmsfr_shell = int(hmsfr_df["in_any_sb_shell"].sum())
    print(
        f"  Grade A/B bubbles: {len(bubble_df):,} total, {n_bubble_shell:,} shell-contact; "
        f"HMSFR: {len(hmsfr_df):,} total, {n_hmsfr_shell:,} shell-contact "
        f"(thick shell u=[{SHELL_U_MIN:.1f},{SHELL_U_MAX:.1f}])"
    )
    save_figure(dust_df, ellipse_df, arm_df, bubble_df, hmsfr_df, OUT_FIG)
    print(f"Saved figure: {OUT_FIG}")

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
