# -*- coding: utf-8 -*-
"""Script 4: OSBs plotting

Purpose
-------
Draw the six-panel identification map for each superbubble.

Method overview
---------------
1. Load fitted geometry, cloud constraints, and the per-target publication plotting
   parameters.
2. Build the Cartesian dust slices and Galactic longitude/latitude panels, applying the
   fitted shell geometry to each projection.
3. Draw the six-panel identification maps with the cloud, bubble and HMSFR overlays;
   render the explicitly supplied manual cases through their separate entry point.

Main inputs
-----------
- ../results/intermediate_output/2_automated_fit_results/automated_superbubble_fit_parameters.csv
- ../results/intermediate_output/2_automated_fit_results/json/SB{N}_fit.json
- ../results/superbubble_final_fit_parameters.csv
- ../data/publication_plotting_parameters.csv
- ../data/Bubbles.csv
- ../data/star_cluster_data/Reid2019_HMSFR.csv
- ../results/intermediate_output/3d_dust_map_products/raw_3d_dust_cube.parquet

Main outputs
------------
- ../results/SB_figures/SB{N}.png; used throughout the Supplementary Information and Extended Data SB panels.

Figure/table role
-----------------
../results/SB_figures/SB{N}.png; used throughout the Supplementary Information and Extended Data SB panels.

Runtime and data notes
----------------------
Reference runtime: 1 h 02 min 35.22 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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

import json
import os
import time
from pathlib import Path

import pandas as pd
import types

import astropy.units as u
try:
    import healpy as hp
except ImportError:
    hp = None
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import numpy as np
import numpy.ma as ma
from astropy.coordinates import CartesianRepresentation, SkyCoord
try:
    from dustmaps3d import dustmaps3d
except ImportError:
    dustmaps3d = None
from matplotlib import patheffects as pe
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse
from PIL import Image
from pyarrow.parquet import ParquetFile
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from scipy.spatial import ConvexHull


# Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
HERE = Path(__file__).resolve().parent
os.chdir(HERE)

DATA_DIR = Path("..") / "data"

MID_DIR = Path("..") / "results" / "intermediate_output"
FINAL_DIR = Path("..") / "results"

FIT_DIR = MID_DIR / "2_automated_fit_results"
JSON_DIR = FIT_DIR / "json"
BESTFIT_CSV = FIT_DIR / "automated_superbubble_fit_parameters.csv"
FINAL_TABLE_CSV = FINAL_DIR / "superbubble_final_fit_parameters.csv"
PANEL_DIR = FINAL_DIR / "SB_figures"

PLOTTING_PARAM_CSV = DATA_DIR / "publication_plotting_parameters.csv"
PUBLISH_PARAM_CSV = DATA_DIR / "publication_plotting_parameters.csv"
BUBBLE_CSV = DATA_DIR / "Bubbles.csv"
HMSFR_CSV = DATA_DIR / "star_cluster_data" / "Reid2019_HMSFR.csv"
DUST_PRODUCTS_DIR = Path("..") / "results" / "intermediate_output" / "3d_dust_map_products"
XY_DUST_PATH = DUST_PRODUCTS_DIR / "raw_3d_dust_cube.parquet"
MANUAL_IDS = [23, 26, 27, 37]


DATA_DIR = Path("..") / "data"


def package_path(path_value) -> Path:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    path = Path(path_value)
    return path


# Figure styling and layout configuration for reproducible rendering.
TITLE_SIZE = 16
LABEL_SIZE = 16
TICK_SIZE = 14
LEGEND_SIZE = 18

FIXED_NSIDE = 1024


# ============================================================================
# Data-loading helper section.
# ============================================================================


def normalize_xyz_dataframe(df):
    rename_map = {}
    if "X" in df.columns and "x" not in df.columns:
        rename_map["X"] = "x"
    if "Y" in df.columns and "y" not in df.columns:
        rename_map["Y"] = "y"
    if "Z" in df.columns and "z" not in df.columns:
        rename_map["Z"] = "z"
    if "dust_raw" in df.columns and "dust" not in df.columns:
        rename_map["dust_raw"] = "dust"
    if rename_map:
        df = df.rename(columns=rename_map)
    return df


def estimate_axis_step(values):
    unique_values = np.sort(np.unique(np.asarray(values, dtype=float)))
    diffs = np.diff(unique_values)
    positive = diffs[diffs > 0]
    if positive.size == 0:
        raise ValueError("Could not estimate grid spacing")
    return float(np.min(positive))


def read_local_xyz_cube(parquet_path, x_range, y_range, z_range):
    parquet = ParquetFile(parquet_path)
    schema_names = parquet.schema.names
    if {"x", "y", "z", "dust"}.issubset(schema_names):
        selected_columns = ["x", "y", "z", "dust"]
    elif {"X", "Y", "Z", "dust"}.issubset(schema_names):
        selected_columns = ["X", "Y", "Z", "dust"]
    elif {"x", "y", "z", "dust_raw"}.issubset(schema_names):
        selected_columns = ["x", "y", "z", "dust_raw"]
    elif {"X", "Y", "Z", "dust_raw"}.issubset(schema_names):
        selected_columns = ["X", "Y", "Z", "dust_raw"]
    else:
        raise ValueError(f"Could not identify XYZ parquet columns: {schema_names}")
    chunks = []
    for row_group in range(parquet.num_row_groups):
        table = parquet.read_row_group(row_group, columns=selected_columns)
        df = normalize_xyz_dataframe(table.to_pandas())
        mask = (
            (df["x"] >= x_range[0]) & (df["x"] <= x_range[1]) &
            (df["y"] >= y_range[0]) & (df["y"] <= y_range[1]) &
            (df["z"] >= z_range[0]) & (df["z"] <= z_range[1])
        )
        if mask.any():
            chunks.append(df.loc[mask, ["x", "y", "z", "dust"]].copy())
    if not chunks:
        raise ValueError("No XYZ data exist in the requested range")
    return pd.concat(chunks, ignore_index=True)


def extract_xy_grid(df, x_range, y_range, z_range, smooth_sigma, z_step):
    cut = df[
        (df["x"] >= x_range[0]) & (df["x"] <= x_range[1]) &
        (df["y"] >= y_range[0]) & (df["y"] <= y_range[1]) &
        (df["z"] >= z_range[0]) & (df["z"] <= z_range[1])
    ].copy()
    if cut.empty:
        return None, None
    cut["x"] = cut["x"].round(3)
    cut["y"] = cut["y"].round(3)
    group = cut.groupby(["x", "y"])["dust"]
    mean_val = group.transform("mean")
    std_val = group.transform("std").fillna(1e-9)
    cut = cut[(cut["dust"] >= mean_val - 3 * std_val) & (cut["dust"] <= mean_val + 3 * std_val)]
    if cut.empty:
        return None, None
    grouped = cut.groupby(["x", "y"], as_index=False)["dust"].sum()
    grouped["dust"] = grouped["dust"] * z_step
    grid = grouped.pivot(index="y", columns="x", values="dust").sort_index(ascending=True).sort_index(axis=1, ascending=True)
    grid.index = pd.Index(np.round(grid.index.to_numpy(dtype=float), 3), name=grid.index.name)
    grid.columns = pd.Index(np.round(grid.columns.to_numpy(dtype=float), 3), name=grid.columns.name)

    # Keep image pixels on the same coordinate phase as the source cube even when
    # the requested panel extends beyond the cube boundary. Without these masked
    # padding rows/columns, imshow stretches the available cube subset over the
    # full requested extent and offsets the dust map from the XY overlays.
    def padded_axis(existing_values, requested_range):
        values = np.sort(np.unique(np.asarray(existing_values, dtype=float)))
        if values.size == 0:
            return values
        step = estimate_axis_step(values)
        anchor = float(values[0])
        start = int(np.floor((float(requested_range[0]) - anchor) / step))
        stop = int(np.ceil((float(requested_range[1]) - anchor) / step))
        return np.round(anchor + np.arange(start, stop + 1) * step, 3)

    full_x = padded_axis(grid.columns.to_numpy(dtype=float), x_range)
    full_y = padded_axis(grid.index.to_numpy(dtype=float), y_range)
    grid = grid.reindex(index=full_y, columns=full_x)
    image = ma.masked_array(gaussian_filter(grid.fillna(0).values, sigma=smooth_sigma), mask=grid.isnull().values)
    return image, grid

# ============================================================================
# Geometry and shell-mask conventions used by this analysis stage.
# ============================================================================


def rotate_points_to_global(points, angle_deg):
    rot = rotation_matrix_z(angle_deg)
    return np.asarray(points, dtype=float) @ rot.T


def wrap_longitude_to_center(lon, center_lon):
    lon_array = np.asarray(lon, dtype=float)
    return center_lon + ((lon_array - center_lon + 180.0) % 360.0) - 180.0


def build_pixel_mask(lon, lat, lon_center, delta_lon, lat_center, delta_lat):
    lon_offset = ((lon - lon_center + 180.0) % 360.0) - 180.0
    lat_min = lat_center - delta_lat
    lat_max = lat_center + delta_lat
    mask = np.abs(lon_offset) <= delta_lon
    mask &= (lat >= lat_min) & (lat <= lat_max)
    return mask


def query_shell_delta_ebv(lon, lat, d_low, d_up):
    if d_up <= d_low:
        raise ValueError(f"Invalid distance range: d_low={d_low}, d_up={d_up}")

    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    distance_low = np.full(lon.shape, float(d_low), dtype=float)
    distance_up = np.full(lon.shape, float(d_up), dtype=float)

    ebv_low, _, _, max_distance_low = dustmaps3d(lon, lat, distance_low)
    ebv_up, _, _, max_distance_up = dustmaps3d(lon, lat, distance_up)

    delta_ebv = np.asarray(ebv_up, dtype=float) - np.asarray(ebv_low, dtype=float)
    max_distance = np.minimum(np.asarray(max_distance_low, dtype=float), np.asarray(max_distance_up, dtype=float))
    d_center = 0.5 * (d_low + d_up)
    valid_mask = d_center < max_distance
    return np.where(valid_mask, delta_ebv, hp.UNSEEN)

# ============================================================================
# Figure-panel rendering section.
# ============================================================================


def load_seed_row(param_csv, target_index):
    param_df = pd.read_csv(param_csv)
    matched = param_df[param_df["index"].astype(int) == int(target_index)]
    if matched.empty:
        best_fit_csv = DATA_DIR / "automated_superbubble_fit_parameters.csv"
        if not best_fit_csv.exists():
            raise ValueError('Invalid input or missing required data.')
        best_df = pd.read_csv(best_fit_csv)
        best_match = best_df[best_df["id"].astype(int) == int(target_index)]
        if best_match.empty:
            raise ValueError('Invalid input or missing required data.')
        row = best_match.iloc[0].to_dict()
        center = np.array([row["center_x_kpc"], row["center_y_kpc"], row["center_z_kpc"]], dtype=float)
        axes = np.array([row["a_radius_kpc"], row["b_radius_kpc"], row["c_radius_kpc"]], dtype=float)
        coord = SkyCoord(
            CartesianRepresentation(center[0] * u.kpc, center[1] * u.kpc, center[2] * u.kpc),
            frame="galactic",
        )
        radius_xy = float(np.nanmax(axes[:2]))
        distance = float(coord.distance.to_value(u.kpc))
        angular_size_deg = float(2.0 * np.degrees(np.arctan2(radius_xy, max(distance, 1e-6))))
        return {
            "index": int(target_index),
            "x": float(center[0]),
            "y": float(center[1]),
            "r": radius_xy,
            "ellipticity": 0.0 if axes[0] <= 0 else 1.0 - float(axes[1] / axes[0]),
            "angle": float(row.get("angle_deg", 0.0)),
            "z_low": float(center[2] - axes[2]),
            "z_up": float(center[2] + axes[2]),
            "d_low": max(0.1, distance - float(np.nanmax(axes))),
            "d_up": distance + float(np.nanmax(axes)),
            "delta_lon": np.nan,
            "delta_lat": np.nan,
            "lon_c": np.nan,
            "lat_c": np.nan,
            "l": float(coord.l.deg),
            "b": float(coord.b.deg),
            "d_kpc": distance,
            "diameter_pc": 2.0 * radius_xy * 1000.0,
            "angular_size_deg": angular_size_deg,
            "sigma_deg": 0.2,
            "mark": int(row.get("mark", 0)),
            "major_axis_pc": 2.0 * float(axes[0]) * 1000.0,
            "minor_axis_pc": 2.0 * float(axes[1]) * 1000.0,
            "angle_deg_ellipse": float(row.get("angle_deg", 0.0)),
        }
    return matched.iloc[0].to_dict()


def superbubble_path_effects(linewidth):
    return [pe.Stroke(linewidth=float(linewidth) + 2.8, foreground="white"), pe.Normal()]


def make_superbubble_handle(linewidth):
    handle = Line2D([0], [0], color="black", linestyle="-", linewidth=float(linewidth), label="Superbubble")
    handle.set_path_effects(superbubble_path_effects(linewidth))
    return handle


def make_catalog_handles(publish_params):
    return [
        Line2D([0], [0], color="black", linestyle="--", linewidth=float(publish_params["bubble_linewidth"]), label="Bubbles"),
        Line2D(
            [0],
            [0],
            marker="*",
            markersize=18,
            markerfacecolor="#00d8ff",
            markeredgecolor="black",
            markeredgewidth=1.8,
            linewidth=0,
            label="Star-Forming Region",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            markersize=10,
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=1.6,
            linewidth=0,
            label="Molecular Clouds",
        ),
    ]


def ensure_publish_param_row(csv_path, target_index):
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"not foundpublication plotting parameters CSV: {csv_path}")
    df = pd.read_csv(csv_path)
    if "index" not in df.columns or df.empty or int(target_index) not in df["index"].astype(int).tolist():
        raise ValueError('Invalid input or missing required data.')
    row = df[df["index"].astype(int) == int(target_index)].iloc[0].to_dict()
    return row


def build_fit_data(result_data):
    effective_mark_value = int(result_data.get("mark", 0))
    return {
        "target_index": int(result_data["target_index"]),
        "fit_shape": str(result_data["fit_shape"]),
        "center_kpc": result_data["center_kpc"],
        "axes_fit_kpc": {
            "a": float(result_data["fit_a_kpc"]),
            "b": float(result_data["fit_b_kpc"]),
            "c": float(result_data["fit_c_kpc"]),
        },
        "angle_deg_fixed": float(result_data["angle_deg"]),
        "hemisphere": str(result_data.get("hemisphere", "")),
        "effective_mark_value": effective_mark_value,
        "manual_override": bool(result_data.get("manual_override", False)),
        "outputs": {
            "boundary_points_csv": str(result_data["outputs"]["selected_clouds_csv"]),
        },
    }


def prepare_bubble_catalog(bubble_csv):
    bubble_df = pd.read_csv(bubble_csv).copy()
    for column in ["l", "b", "D", "physical_size", "best_distance"]:
        bubble_df[column] = pd.to_numeric(bubble_df[column], errors="coerce")
    bubble_df["Grade"] = bubble_df["Grade"].astype(str)
    if "calc_failed" in bubble_df.columns:
        bubble_df["calc_failed"] = bubble_df["calc_failed"].astype(str).str.lower().eq("true")
    else:
        bubble_df["calc_failed"] = False
    bubble_df = bubble_df[
        bubble_df["Grade"].isin(["A", "B"])
        & ~bubble_df["calc_failed"]
        & np.isfinite(bubble_df["l"])
        & np.isfinite(bubble_df["b"])
        & np.isfinite(bubble_df["D"])
        & np.isfinite(bubble_df["physical_size"])
        & np.isfinite(bubble_df["best_distance"])
    ].copy()
    coords = SkyCoord(
        l=bubble_df["l"].to_numpy(dtype=float) * u.deg,
        b=bubble_df["b"].to_numpy(dtype=float) * u.deg,
        distance=bubble_df["best_distance"].to_numpy(dtype=float) * u.kpc,
        frame="galactic",
    )
    bubble_df["x_kpc"] = coords.cartesian.x.to_value(u.kpc)
    bubble_df["y_kpc"] = coords.cartesian.y.to_value(u.kpc)
    bubble_df["z_kpc"] = coords.cartesian.z.to_value(u.kpc)
    bubble_df["diameter_kpc"] = bubble_df["physical_size"].to_numpy(dtype=float) / 1000.0
    return bubble_df


def prepare_hmsfr_catalog(hmsfr_csv):
    hmsfr_df = pd.read_csv(hmsfr_csv).copy()
    required = ["l_deg", "b_deg", "distance_kpc", "x_kpc", "y_kpc", "z_kpc"]
    for column in required:
        hmsfr_df[column] = pd.to_numeric(hmsfr_df[column], errors="coerce")
    hmsfr_df = hmsfr_df[np.all([np.isfinite(hmsfr_df[col]) for col in required], axis=0)].copy()
    return hmsfr_df


def section_axes_for_mark(axes_fit, center_z, z_value, mark, fit_shape):
    if fit_shape == "cylinder":
        return float(axes_fit[0]), float(axes_fit[1])
    dz = float(z_value) - float(center_z)
    if int(mark) == 1 and dz > 0:
        return float(axes_fit[0]), float(axes_fit[1])
    if int(mark) == 2 and dz < 0:
        return float(axes_fit[0]), float(axes_fit[1])
    c_axis = float(axes_fit[2])
    if c_axis <= 0:
        return None
    ratio = 1.0 - (dz / c_axis) ** 2
    if ratio <= 0:
        return None
    factor = float(np.sqrt(ratio))
    return float(axes_fit[0]) * factor, float(axes_fit[1]) * factor


def compute_cap_z_from_selected_clouds(selected_clouds, center_z, mark_value):
    if selected_clouds is None or len(selected_clouds) == 0:
        return float("nan")
    z_values = np.asarray(selected_clouds["z"], dtype=float)
    center_z = float(center_z)
    mark_value = int(mark_value)
    if mark_value == 1:
        subset = np.sort(z_values[z_values < center_z])[::-1]
        if subset.size == 0:
            return float("nan")
        rank = max(1, int(np.ceil(subset.size / np.e)))
        return float(subset[min(rank - 1, subset.size - 1)])
    if mark_value == 2:
        subset = np.sort(z_values[z_values > center_z])
        if subset.size == 0:
            return float("nan")
        rank = max(1, int(np.ceil(subset.size / np.e)))
        return float(subset[min(rank - 1, subset.size - 1)])
    return float("nan")


def section_axes_for_cap(axes_fit, center_z, z_value, mark, fit_shape, cap_z_global):
    if fit_shape == "cylinder":
        return float(axes_fit[0]), float(axes_fit[1])
    if int(mark) not in (1, 2) or not np.isfinite(cap_z_global):
        return section_axes_for_mark(axes_fit, center_z, z_value, mark, fit_shape)
    z_value = float(z_value)
    cap_z_global = float(cap_z_global)
    if int(mark) == 1 and z_value > cap_z_global:
        return None
    if int(mark) == 2 and z_value < cap_z_global:
        return None
    c_axis = float(axes_fit[2])
    if c_axis <= 0:
        return None
    dz = z_value - float(center_z)
    ratio = 1.0 - (dz / c_axis) ** 2
    if ratio <= 0:
        return None
    factor = float(np.sqrt(ratio))
    return float(axes_fit[0]) * factor, float(axes_fit[1]) * factor


def cap_plane_section_axes(axes_fit, center_z, fit_shape, cap_z_global):
    if fit_shape == "cylinder":
        return float(axes_fit[0]), float(axes_fit[1])
    if not np.isfinite(cap_z_global):
        return None
    c_axis = float(axes_fit[2])
    if c_axis <= 0:
        return None
    dz = float(cap_z_global) - float(center_z)
    ratio = 1.0 - (dz / c_axis) ** 2
    if ratio <= 0:
        return None
    factor = float(np.sqrt(ratio))
    return float(axes_fit[0]) * factor, float(axes_fit[1]) * factor


def z_slice_is_beyond_cap(z_value, mark, cap_z_global):
    if int(mark) not in (1, 2) or not np.isfinite(cap_z_global):
        return False
    if int(mark) == 1:
        return float(z_value) > float(cap_z_global)
    return float(z_value) < float(cap_z_global)


def require_publish_param_values(publish_params, keys):
    missing = [key for key in keys if pd.isna(publish_params.get(key))]
    if missing:
        target_index = int(publish_params["index"])
        raise ValueError(
            f"Missing publication plotting parameters for target index {target_index}: {missing}. "
            "Fill ../data/publication_plotting_parameters.csv before plotting."
        )
    return [float(publish_params[key]) for key in keys]


def validate_publish_params_for_rendering(publish_params):
    values = require_publish_param_values(
        publish_params,
        [
            "xy_vmin",
            "xy_vmax",
            "lb_vmin",
            "lb_vmax",
            "xy_panel_scale",
            "hmsfr_marker_size",
            "cloud_marker_size",
            "bubble_linewidth",
            "superbubble_linewidth",
            "xy_contour_levels",
            "lb_projection_min_distance_kpc",
        ],
    )
    xy_vmin, xy_vmax, lb_vmin, lb_vmax, xy_panel_scale = values[:5]
    if xy_vmax <= xy_vmin or lb_vmax <= lb_vmin or xy_panel_scale <= 0:
        raise ValueError('Invalid input or missing required data.')
    projection_mode = publish_params.get("lb_projection_mode")
    if pd.isna(projection_mode) or not str(projection_mode).strip():
        raise ValueError('Invalid input or missing required data.')
    xy_panel_shrink_kpc = publish_params.get("xy_panel_shrink_kpc", 0.0)
    if pd.notna(xy_panel_shrink_kpc) and float(xy_panel_shrink_kpc) < 0:
        raise ValueError('Invalid input or missing required data.')
    build_xy_z_slices_from_params(publish_params)
    build_lb_distance_slices(publish_params)
    apply_lb_window_overrides({}, publish_params)


def build_xy_z_slices_from_params(publish_params):
    keys = [f"xy_z{index}_{bound}" for index in range(1, 4) for bound in ("low", "high")]
    values = require_publish_param_values(publish_params, keys)
    z_slices = [[values[index], values[index + 1]] for index in range(0, len(values), 2)]
    for low_value, high_value in z_slices:
        if high_value <= low_value:
            raise ValueError('Invalid input or missing required data.')
    return z_slices


def build_lb_distance_slices(publish_params):
    keys = [f"lb_d{index}_{bound}" for index in range(1, 4) for bound in ("low", "high")]
    values = require_publish_param_values(publish_params, keys)
    slices = [[values[index], values[index + 1]] for index in range(0, len(values), 2)]
    for low_value, high_value in slices:
        if high_value <= low_value:
            raise ValueError('Invalid input or missing required data.')
    return slices


def build_xy_panel_ranges(center, axes_fit, seed_row, publish_params):
    base_half_width = float(publish_params["xy_panel_scale"]) * max(float(axes_fit[0]), float(seed_row.get("r", axes_fit[0])))
    shrink_value = publish_params.get("xy_panel_shrink_kpc", 0.0)
    shrink_kpc = 0.0 if pd.isna(shrink_value) else float(shrink_value)
    panel_half_width = base_half_width - shrink_kpc
    if panel_half_width <= 0:
        raise ValueError(
            "xy_panel_scale and xy_panel_shrink_kpc produce a non-positive "
            "panel half-width. Increase the scale or reduce the shrink value."
        )
    x_range = [float(center[0]) - panel_half_width, float(center[0]) + panel_half_width]
    y_range = [float(center[1]) - panel_half_width, float(center[1]) + panel_half_width]
    return apply_xy_range_overrides(publish_params, x_range, y_range)


def apply_xy_range_overrides(publish_params, x_range, y_range):
    override_map = {
        "xy_x_min_kpc": (x_range, 0),
        "xy_x_max_kpc": (x_range, 1),
        "xy_y_min_kpc": (y_range, 0),
        "xy_y_max_kpc": (y_range, 1),
    }
    for key, (target_range, index) in override_map.items():
        if key in publish_params and pd.notna(publish_params.get(key)):
            target_range[index] = float(publish_params[key])
    if not x_range[0] < x_range[1]:
        raise ValueError('Invalid input or missing required data.')
    if not y_range[0] < y_range[1]:
        raise ValueError('Invalid input or missing required data.')
    return x_range, y_range


def clip_xy_range_to_available_cube(df_xyz, x_range, y_range):
    """Clip XY panel limits to the actual dust cube extent read from disk."""
    def axis_extent(values):
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            raise ValueError('Invalid input or missing required data.')
        if np.unique(values).size > 1:
            step = estimate_axis_step(values)
        else:
            step = 0.01
        return float(np.min(values) - 0.5 * step), float(np.max(values) + 0.5 * step)

    x_min, x_max = axis_extent(df_xyz["x"].to_numpy())
    y_min, y_max = axis_extent(df_xyz["y"].to_numpy())
    clipped_x = [max(float(x_range[0]), x_min), min(float(x_range[1]), x_max)]
    clipped_y = [max(float(y_range[0]), y_min), min(float(y_range[1]), y_max)]
    if not clipped_x[0] < clipped_x[1] or not clipped_y[0] < clipped_y[1]:
        raise ValueError(
            "Requested XY panel does not overlap the available dust cube extent; "
            f"requested X={x_range}, Y={y_range}; cube X={[x_min, x_max]}, Y={[y_min, y_max]}"
        )
    return clipped_x, clipped_y


def apply_lb_window_overrides(seed_row, publish_params):
    lon_c, delta_lon, lat_min, lat_max = require_publish_param_values(
        publish_params,
        ["lb_lon_c", "lb_delta_lon", "lb_lat_min", "lb_lat_max"],
    )
    if delta_lon <= 0 or lat_max <= lat_min:
        raise ValueError('Invalid input or missing required data.')
    seed_row = seed_row.copy()
    seed_row["lon_c"] = lon_c
    seed_row["delta_lon"] = delta_lon
    seed_row["lat_c"] = 0.5 * (lat_min + lat_max)
    seed_row["delta_lat"] = 0.5 * (lat_max - lat_min)
    return seed_row


def longitude_window_to_center_delta(lon_min, lon_max):
    lon_min = float(lon_min) % 360.0
    lon_max = float(lon_max) % 360.0
    lon_max_unwrapped = lon_max
    if lon_max_unwrapped < lon_min:
        lon_max_unwrapped += 360.0
    center = 0.5 * (lon_min + lon_max_unwrapped)
    delta = 0.5 * (lon_max_unwrapped - lon_min)
    return float(center % 360.0), float(delta)


def build_lb_panel_window(publish_params, panel_index, default_lon_c, default_delta_lon, default_lat_min, default_lat_max):
    panel_number = int(panel_index) + 1
    lon_min = publish_params.get(f"lb_p{panel_number}_lon_min")
    lon_max = publish_params.get(f"lb_p{panel_number}_lon_max")
    lat_min = publish_params.get(f"lb_p{panel_number}_lat_min")
    lat_max = publish_params.get(f"lb_p{panel_number}_lat_max")

    if pd.notna(lon_min) != pd.notna(lon_max):
        raise ValueError('Invalid input or missing required data.')
    if pd.notna(lon_min) and pd.notna(lon_max):
        lon_c, delta_lon = longitude_window_to_center_delta(lon_min, lon_max)
    else:
        lon_c, delta_lon = float(default_lon_c), float(default_delta_lon)

    if pd.notna(lat_min) != pd.notna(lat_max):
        raise ValueError('Invalid input or missing required data.')
    if pd.notna(lat_min) and pd.notna(lat_max):
        lat_min_value, lat_max_value = float(lat_min), float(lat_max)
    else:
        lat_min_value, lat_max_value = float(default_lat_min), float(default_lat_max)

    return {
        "lon_c": lon_c,
        "delta_lon": delta_lon,
        "lat_min": lat_min_value,
        "lat_max": lat_max_value,
    }


def draw_solar_xy_grid(ax, x_range, y_range, reference_radii=None):
    radius_max = float(np.hypot(max(abs(x_range[0]), abs(x_range[1])), max(abs(y_range[0]), abs(y_range[1]))))
    for radius in np.arange(0.5, radius_max + 0.25, 0.5):
        ax.add_patch(
            Circle(
                (0.0, 0.0),
                radius=radius,
                edgecolor="black",
                facecolor="none",
                linestyle="--",
                linewidth=0.95,
                alpha=0.34,
                zorder=3,
            )
        )
    for lon_deg in np.arange(0.0, 360.0, 10.0):
        theta = np.deg2rad(lon_deg)
        ax.plot(
            [0.0, radius_max * np.cos(theta)],
            [0.0, radius_max * np.sin(theta)],
            color="black",
            linestyle="--",
            linewidth=0.85,
            alpha=0.30,
            zorder=3,
        )
    if reference_radii is not None:
        for radius in reference_radii:
            radius = float(radius)
            if radius <= 0:
                continue
            ax.add_patch(
                Circle(
                    (0.0, 0.0),
                    radius=radius,
                    edgecolor="white",
                    facecolor="none",
                    linestyle="-",
                    linewidth=1.8,
                    alpha=0.95,
                    zorder=5,
                )
            )


def draw_xy_overlays(ax, z_slice, selected_clouds, bubble_df, hmsfr_df, publish_params):
    z_low, z_high = float(z_slice[0]), float(z_slice[1])

    cloud_mask = (selected_clouds["z"] >= z_low) & (selected_clouds["z"] <= z_high)
    if cloud_mask.any():
        ax.scatter(
            selected_clouds.loc[cloud_mask, "x"],
            selected_clouds.loc[cloud_mask, "y"],
            s=0.5 * float(publish_params["cloud_marker_size"]),
            facecolors="white",
            edgecolors="black",
            linewidths=1.6,
            zorder=7,
        )

    bubble_mask = (bubble_df["z_kpc"] >= z_low) & (bubble_df["z_kpc"] <= z_high)
    for row in bubble_df.loc[bubble_mask].itertuples(index=False):
        ax.add_patch(
            Circle(
                (float(row.x_kpc), float(row.y_kpc)),
                radius=0.5 * float(row.diameter_kpc),
                edgecolor="black",
                facecolor="none",
                linestyle="--",
                linewidth=float(publish_params["bubble_linewidth"]),
                alpha=0.9,
                zorder=5,
            )
        )

    hmsfr_mask = (hmsfr_df["z_kpc"] >= z_low) & (hmsfr_df["z_kpc"] <= z_high)
    if hmsfr_mask.any():
        ax.scatter(
            hmsfr_df.loc[hmsfr_mask, "x_kpc"],
            hmsfr_df.loc[hmsfr_mask, "y_kpc"],
            s=float(publish_params["hmsfr_marker_size"]),
            marker="*",
            c="#00d8ff",
            edgecolors="black",
            linewidths=1.8,
            zorder=8,
        )


def draw_manual_extra_xy_ellipse(ax, z_slice, publish_params):
    keys = [
        "xy_extra_ellipse_center_x_kpc",
        "xy_extra_ellipse_center_y_kpc",
        "xy_extra_ellipse_z_kpc",
        "xy_extra_ellipse_a_pc",
        "xy_extra_ellipse_b_pc",
        "xy_extra_ellipse_pa_deg",
    ]
    if not all(key in publish_params and pd.notna(publish_params.get(key)) for key in keys):
        return

    ax.add_patch(
        Ellipse(
            (
                float(publish_params["xy_extra_ellipse_center_x_kpc"]),
                float(publish_params["xy_extra_ellipse_center_y_kpc"]),
            ),
            width=2.0 * float(publish_params["xy_extra_ellipse_a_pc"]) / 1000.0,
            height=2.0 * float(publish_params["xy_extra_ellipse_b_pc"]) / 1000.0,
            angle=float(publish_params["xy_extra_ellipse_pa_deg"]),
            edgecolor="black",
            facecolor="none",
            linestyle="--",
            linewidth=2.0,
            zorder=7,
        )
    )


def points_to_lb_dataframe(points):
    coords = SkyCoord(
        CartesianRepresentation(
            points[:, 0] * u.kpc,
            points[:, 1] * u.kpc,
            points[:, 2] * u.kpc,
        ),
        frame="galactic",
    )
    return pd.DataFrame({"l": coords.l.deg, "b": coords.b.deg, "d": coords.distance.kpc})


def lb_overlay_x(lon, center_lon):
    lon_values = np.asarray(lon, dtype=float)
    return -(((lon_values - float(center_lon) + 180.0) % 360.0) - 180.0)


def lb_tick_label(x_value, center_lon):
    return f"{int(round((float(center_lon) - float(x_value)) % 360.0))}°"


def split_curve_on_jumps(x_values, y_values, jump_threshold=20.0):
    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)
    if x_values.size <= 1:
        return [(x_values, y_values)]
    split_indices = np.where(np.abs(np.diff(x_values)) > jump_threshold)[0] + 1
    return [
        (x_segment, y_segment)
        for x_segment, y_segment in zip(np.split(x_values, split_indices), np.split(y_values, split_indices))
        if len(x_segment) > 1
    ]


def split_cyclic_valid_indices(valid_mask):
    valid_indices = np.flatnonzero(np.asarray(valid_mask, dtype=bool))
    if valid_indices.size == 0:
        return []
    breaks = np.where(np.diff(valid_indices) > 1)[0] + 1
    groups = [group for group in np.split(valid_indices, breaks) if group.size > 0]
    if len(groups) > 1 and groups[0][0] == 0 and groups[-1][-1] == len(valid_mask) - 1:
        groups[0] = np.concatenate([groups[-1], groups[0]])
        groups.pop()
    return groups


def split_linear_valid_indices(valid_mask):
    valid_indices = np.flatnonzero(np.asarray(valid_mask, dtype=bool))
    if valid_indices.size == 0:
        return []
    breaks = np.where(np.diff(valid_indices) > 1)[0] + 1
    return [group for group in np.split(valid_indices, breaks) if group.size > 0]


def rotation_matrix_z(angle_deg):
    theta = np.deg2rad(float(angle_deg))
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rotate_points_to_local(points, angle_deg):
    rot = rotation_matrix_z(angle_deg)
    return np.asarray(points, dtype=float) @ rot


def spherical_to_galactic_cartesian(distance_kpc, lon_deg, lat_deg):
    lon_rad = np.radians(np.asarray(lon_deg, dtype=float))
    lat_rad = np.radians(np.asarray(lat_deg, dtype=float))
    distance_kpc = np.asarray(distance_kpc, dtype=float)
    cos_lat = np.cos(lat_rad)
    return np.stack(
        [
            distance_kpc * cos_lat * np.cos(lon_rad),
            distance_kpc * cos_lat * np.sin(lon_rad),
            distance_kpc * np.sin(lat_rad),
        ],
        axis=-1,
    )


def ellipsoid_cap_contains_points(local_points, global_points, axes, mark_value, cap_z_global):
    axes = np.asarray(axes, dtype=float)
    local_points = np.asarray(local_points, dtype=float)
    global_points = np.asarray(global_points, dtype=float)
    if np.any(axes <= 0):
        return np.zeros(len(local_points), dtype=bool)
    norm = (
        (local_points[:, 0] / float(axes[0])) ** 2
        + (local_points[:, 1] / float(axes[1])) ** 2
        + (local_points[:, 2] / float(axes[2])) ** 2
    )
    inside = norm <= 1.0 + 1e-9
    if int(mark_value) == 1 and np.isfinite(cap_z_global):
        inside &= global_points[:, 2] <= float(cap_z_global) + 1e-9
    elif int(mark_value) == 2 and np.isfinite(cap_z_global):
        inside &= global_points[:, 2] >= float(cap_z_global) - 1e-9
    elif int(mark_value) == 1:
        inside &= local_points[:, 2] <= 1e-9
    elif int(mark_value) == 2:
        inside &= local_points[:, 2] >= -1e-9
    return inside


def build_cylinder_distance_slice_lb_segments(lbpub, center, axes, angle_deg, d_mid, lat_min, lat_max):
    total_height_deg = 0.5 * max(float(lat_max) - float(lat_min), 0.0)
    half_height_deg = 0.5 * total_height_deg
    half_height_kpc = abs(float(d_mid) * np.tan(np.radians(half_height_deg)))
    if half_height_kpc <= 0:
        return []

    phi = np.linspace(0.0, 2.0 * np.pi, 721)
    local_x = float(axes[0]) * np.cos(phi)
    local_y = float(axes[1]) * np.sin(phi)
    local_xy = np.column_stack([local_x, local_y, np.zeros_like(phi)])
    base_global = lbpub.rotate_points_to_global(local_xy, angle_deg) + center
    radial_sq = np.sum(base_global[:, :2] ** 2, axis=1)
    discriminant = float(d_mid) ** 2 - radial_sq
    if np.nanmax(discriminant) < 0:
        return []

    sqrt_disc = np.sqrt(np.clip(discriminant, 0.0, None))
    z0 = base_global[:, 2]
    lb_segments = []
    for sign in (-1.0, 1.0):
        local_z = -z0 + sign * sqrt_disc
        valid = (discriminant >= 0.0) & (np.abs(local_z) <= half_height_kpc + 1e-9)
        for segment_indices in split_cyclic_valid_indices(valid):
            if segment_indices.size < 2:
                continue
            local_points = np.column_stack(
                [
                    local_x[segment_indices],
                    local_y[segment_indices],
                    local_z[segment_indices],
                ]
            )
            global_points = lbpub.rotate_points_to_global(local_points, angle_deg) + center
            lb_df = points_to_lb_dataframe(global_points)
            if len(lb_df) >= 2:
                lb_segments.append(lb_df)
    return lb_segments


def cylinder_distance_slice_intersects(lbpub, center, axes, angle_deg, d_mid, lat_min, lat_max):
    total_height_deg = 0.5 * max(float(lat_max) - float(lat_min), 0.0)
    half_height_deg = 0.5 * total_height_deg
    half_height_kpc = abs(float(d_mid) * np.tan(np.radians(half_height_deg)))
    if half_height_kpc <= 0:
        return False

    phi = np.linspace(0.0, 2.0 * np.pi, 721)
    local_xy = np.column_stack(
        [
            float(axes[0]) * np.cos(phi),
            float(axes[1]) * np.sin(phi),
            np.zeros_like(phi),
        ]
    )
    base_global = lbpub.rotate_points_to_global(local_xy, angle_deg) + center
    radial_sq = np.sum(base_global[:, :2] ** 2, axis=1)
    discriminant = float(d_mid) ** 2 - radial_sq
    if np.nanmax(discriminant) < 0:
        return False

    sqrt_disc = np.sqrt(np.clip(discriminant, 0.0, None))
    z0 = base_global[:, 2]
    local_z_pos = -z0 + sqrt_disc
    local_z_neg = -z0 - sqrt_disc
    valid = (discriminant >= 0.0) & (
        (np.abs(local_z_pos) <= half_height_kpc + 1e-9)
        | (np.abs(local_z_neg) <= half_height_kpc + 1e-9)
    )
    return bool(np.any(valid))


def build_distance_slice_lb_segments(lbpub, fit_overlay, d_mid, lon_c=None, delta_lon=None, lat_min=None, lat_max=None, include_cap=False, respect_mark=True):
    if fit_overlay is None:
        return []
    fit_shape = str(fit_overlay.get("fit_shape", "")).lower()
    if fit_shape == "cylinder":
        if lat_min is None or lat_max is None:
            return []
        return build_cylinder_distance_slice_lb_segments(
            lbpub,
            np.asarray(fit_overlay["center_kpc"], dtype=float),
            np.asarray(fit_overlay["axes_kpc"], dtype=float),
            float(fit_overlay["angle_deg"]),
            float(d_mid),
            float(lat_min),
            float(lat_max),
        )
    if fit_shape != "ellipsoid":
        return []
    center = np.asarray(fit_overlay["center_kpc"], dtype=float)
    axes = np.asarray(fit_overlay["axes_kpc"], dtype=float)
    angle_deg = float(fit_overlay["angle_deg"])
    mark_value = int(fit_overlay.get("mark_value", 0)) if respect_mark else 0
    if lon_c is None or delta_lon is None or lat_min is None or lat_max is None:
        return []

    lon_c = float(lon_c)
    delta_lon = float(delta_lon)
    lat_min = float(lat_min)
    lat_max = float(lat_max)
    d_mid = float(d_mid)
    if delta_lon <= 0 or lat_max <= lat_min or d_mid <= 0:
        return []

    nx = 721
    ny = max(241, int(round(nx * (lat_max - lat_min) / max(2.0 * delta_lon, 1e-6))))
    x_grid = np.linspace(-delta_lon, delta_lon, nx)
    b_grid = np.linspace(lat_min, lat_max, ny)
    xx, bb = np.meshgrid(x_grid, b_grid)
    ll = (lon_c - xx) % 360.0
    global_points = spherical_to_galactic_cartesian(d_mid, ll, bb).reshape(-1, 3)
    local_points = rotate_points_to_local(global_points - center, angle_deg)
    cap_z_global = float(fit_overlay.get("cap_z_global", np.nan))
    inside_mask = ellipsoid_cap_contains_points(local_points, global_points, axes, mark_value, cap_z_global).reshape(bb.shape)
    if not np.any(inside_mask):
        return []
    if np.all(inside_mask):
        return []

    temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
    try:
        contour = temp_ax.contour(x_grid, b_grid, inside_mask.astype(float), levels=[0.5])
        contour_segments = contour.allsegs[0] if contour.allsegs else []
    finally:
        plt.close(temp_fig)

    lb_segments = []
    for segment in contour_segments:
        if len(segment) < 3:
            continue
        segment_lon = (lon_c - segment[:, 0]) % 360.0
        segment_lat = segment[:, 1]
        segment_df = pd.DataFrame(
            {
                "l": segment_lon,
                "b": segment_lat,
                "d": np.full(len(segment), d_mid, dtype=float),
            }
        )
        lb_segments.append(segment_df)
    return lb_segments


def split_lb_segments_by_cap(lb_segments, d_mid, mark_value, cap_z_global):
    if not lb_segments:
        return {"solid": [], "dashed": []}
    if int(mark_value) not in (1, 2) or not np.isfinite(cap_z_global):
        return {"solid": lb_segments, "dashed": []}
    solid_segments = []
    dashed_segments = []
    for curve_df in lb_segments:
        if curve_df is None or curve_df.empty:
            continue
        global_points = spherical_to_galactic_cartesian(
            float(d_mid),
            curve_df["l"].to_numpy(dtype=float),
            curve_df["b"].to_numpy(dtype=float),
        )
        if int(mark_value) == 1:
            solid_mask = global_points[:, 2] <= float(cap_z_global) + 1e-9
        else:
            solid_mask = global_points[:, 2] >= float(cap_z_global) - 1e-9
        for idx in split_linear_valid_indices(solid_mask):
            if idx.size >= 2:
                solid_segments.append(curve_df.iloc[idx].reset_index(drop=True))
        dashed_mask = ~solid_mask
        for idx in split_linear_valid_indices(dashed_mask):
            if idx.size >= 2:
                dashed_segments.append(curve_df.iloc[idx].reset_index(drop=True))
    return {"solid": solid_segments, "dashed": dashed_segments}


def plot_lb_curve_segments(ax, curve_segments, lon_c, *, color, linestyle, linewidth, alpha, zorder=6, path_effects=None):
    for curve_df in curve_segments:
        if curve_df is None or curve_df.empty:
            continue
        x_values = lb_overlay_x(curve_df["l"].to_numpy(dtype=float), lon_c)
        y_values = curve_df["b"].to_numpy(dtype=float)
        for x_seg, y_seg in split_curve_on_jumps(x_values, y_values):
            if len(x_seg) < 2:
                continue
            (line,) = ax.plot(
                x_seg,
                y_seg,
                color=color,
                linestyle=linestyle,
                linewidth=linewidth,
                alpha=alpha,
                zorder=zorder,
            )
            if path_effects is not None:
                line.set_path_effects(path_effects)


def build_fit_overlay_for_lb(lbpub, seed_row, fit_data, projection_min_distance_kpc=0.0):
    center = np.asarray(fit_data["center_kpc"], dtype=float)
    axes = np.array(
        [
            float(fit_data["axes_fit_kpc"]["a"]),
            float(fit_data["axes_fit_kpc"]["b"]),
            float(fit_data["axes_fit_kpc"]["c"]),
        ],
        dtype=float,
    )
    angle_deg = float(fit_data["angle_deg_fixed"])
    fit_shape = str(fit_data["fit_shape"]).lower()
    mark_value = int(fit_data.get("effective_mark_value", seed_row.get("mark", 0)))
    map_center_lon = float(seed_row.get("lon_c", seed_row.get("l", 0.0)))

    line_lons = None
    sky_outline_df = pd.DataFrame(columns=["l", "b", "d"])
    if fit_shape == "cylinder":
        phi = np.linspace(0.0, 2.0 * np.pi, 360, endpoint=False)
        local = np.stack([axes[0] * np.cos(phi), axes[1] * np.sin(phi), np.zeros_like(phi)], axis=-1)
        global_points = lbpub.rotate_points_to_global(local, angle_deg) + center
        perimeter_df = points_to_lb_dataframe(global_points)
        wrapped_lon = lbpub.wrap_longitude_to_center(perimeter_df["l"].to_numpy(dtype=float), map_center_lon)
        line_lons = [float(np.min(wrapped_lon)), float(np.max(wrapped_lon))]
    else:
        phi = np.linspace(0.0, 2.0 * np.pi, 360, endpoint=False)
        if mark_value == 1:
            theta = np.linspace(-0.5 * np.pi, 0.0, 180)
        elif mark_value == 2:
            theta = np.linspace(0.0, 0.5 * np.pi, 180)
        else:
            theta = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 240)
        uu, vv = np.meshgrid(phi, theta)
        local = np.stack(
            [
                axes[0] * np.cos(vv) * np.cos(uu),
                axes[1] * np.cos(vv) * np.sin(uu),
                axes[2] * np.sin(vv),
            ],
            axis=-1,
        ).reshape(-1, 3)
        projected = points_to_lb_dataframe(lbpub.rotate_points_to_global(local, angle_deg) + center)
        x_proj = lb_overlay_x(projected["l"].to_numpy(dtype=float), map_center_lon)
        y_proj = projected["b"].to_numpy(dtype=float)
        min_distance = max(0.0, float(projection_min_distance_kpc))
        finite = np.isfinite(x_proj) & np.isfinite(y_proj) & (projected["d"].to_numpy(dtype=float) >= min_distance)
        if np.count_nonzero(finite) >= 3:
            points_2d = np.column_stack([x_proj[finite], y_proj[finite]])
            hull = ConvexHull(points_2d)
            hull_indices = np.append(hull.vertices, hull.vertices[0])
            finite_indices = np.flatnonzero(finite)
            sky_outline_df = projected.iloc[finite_indices[hull_indices]].copy()

    boundary_df = None
    raw_boundary = None
    cap_z_global = float("nan")
    boundary_path = package_path(fit_data["outputs"]["boundary_points_csv"])
    if boundary_path.exists():
        raw_boundary = pd.read_csv(boundary_path)
        if {"x", "y", "z"}.issubset(raw_boundary.columns):
            boundary_df = points_to_lb_dataframe(raw_boundary[["x", "y", "z"]].to_numpy(dtype=float))
            cap_z_global = compute_cap_z_from_selected_clouds(raw_boundary, center[2], mark_value)
            # Geometry and shell-mask conventions used by this analysis stage.

            if bool(fit_data.get("manual_override", False)):
                cap_z_global = float(center[2])

    return {
        "fit_shape": fit_shape,
        "line_lons": line_lons,
        "sky_outline_df": sky_outline_df,
        "boundary_df": boundary_df,
        "raw_boundary_df": raw_boundary,
        "cap_z_global": cap_z_global,
        "center_kpc": center.tolist(),
        "axes_kpc": axes.tolist(),
        "angle_deg": angle_deg,
        "mark_value": mark_value,
    }


def plot_lb_curve_in_distance_bin(
    lbpub,
    ax,
    curve_df,
    d_lo,
    d_hi,
    lon_c,
    *,
    color,
    linestyle,
    linewidth,
    alpha,
    label=None,
    zorder=6,
    filter_distance=True,
    path_effects=None,
):
    if curve_df is None or curve_df.empty:
        return
    if filter_distance:
        mask = (curve_df["d"] >= d_lo) & (curve_df["d"] <= d_hi)
    else:
        mask = pd.Series(True, index=curve_df.index)
    if not mask.any():
        return
    indices = np.flatnonzero(mask.to_numpy())
    breaks = np.where(np.diff(indices) > 1)[0] + 1
    segments = np.split(indices, breaks)
    used_label = False
    for segment in segments:
        if len(segment) < 2:
            continue
        subset = curve_df.iloc[segment]
        x_values = lb_overlay_x(subset["l"].to_numpy(dtype=float), lon_c)
        y_values = subset["b"].to_numpy(dtype=float)
        for x_seg, y_seg in split_curve_on_jumps(x_values, y_values):
            (line,) = ax.plot(
                x_seg,
                y_seg,
                color=color,
                linestyle=linestyle,
                linewidth=linewidth,
                alpha=alpha,
                label=label if label and not used_label else None,
                zorder=zorder,
            )
            if path_effects is not None:
                line.set_path_effects(path_effects)
            used_label = True


def draw_lb_bubbles(lbpub, ax, bubble_df, lon_center, publish_params):
    for row in bubble_df.itertuples(index=False):
        center = SkyCoord(l=float(row.l) * u.deg, b=float(row.b) * u.deg, frame="galactic")
        bearings = np.linspace(0.0, 360.0, 361) * u.deg
        circle = center.directional_offset_by(bearings, 0.5 * float(row.D) * u.deg)
        x_values = lb_overlay_x(circle.l.deg, lon_center)
        y_values = circle.b.deg
        for x_seg, y_seg in split_curve_on_jumps(x_values, y_values):
            ax.plot(
                x_seg,
                y_seg,
                color="black",
                linestyle="--",
                linewidth=float(publish_params["bubble_linewidth"]),
                alpha=0.9,
                zorder=5,
            )


def sky_outline_metrics(lbpub, fit_overlay, center_lb, lon_c):
    outline = fit_overlay.get("sky_outline_df")
    if outline is None or outline.empty:
        return 0.0, 0.0
    x_values = lb_overlay_x(outline["l"].to_numpy(dtype=float), lon_c)
    y_values = outline["b"].to_numpy(dtype=float)
    finite = np.isfinite(x_values) & np.isfinite(y_values)
    if not np.any(finite):
        return 0.0, 0.0
    lon_span = float(np.nanmax(x_values[finite]) - np.nanmin(x_values[finite]))
    lat_span = float(np.nanmax(y_values[finite]) - np.nanmin(y_values[finite]))
    return abs(lon_span), abs(lat_span)


def lb_curve_segments_metrics(curve_segments, lon_c):
    x_all = []
    y_all = []
    for curve_df in curve_segments:
        if curve_df is None or curve_df.empty:
            continue
        x_all.append(lb_overlay_x(curve_df["l"].to_numpy(dtype=float), lon_c))
        y_all.append(curve_df["b"].to_numpy(dtype=float))
    if not x_all:
        return 0.0, 0.0
    x_values = np.concatenate(x_all)
    y_values = np.concatenate(y_all)
    finite = np.isfinite(x_values) & np.isfinite(y_values)
    if not np.any(finite):
        return 0.0, 0.0
    lon_span = float(np.nanmax(x_values[finite]) - np.nanmin(x_values[finite]))
    lat_span = float(np.nanmax(y_values[finite]) - np.nanmin(y_values[finite]))
    return abs(lon_span), abs(lat_span)


def build_lb_projected_image(lbpub, npix, pix_ids, lon, lat, d_lo, d_hi, *, lon_c, delta_lon, lat_min, lat_max, smoothing_sigma, vmin, vmax):
    dust = lbpub.query_shell_delta_ebv(lon, lat, d_lo, d_hi)
    dust_map = np.full(npix, hp.UNSEEN)
    dust_map[pix_ids] = dust
    if smoothing_sigma > 0:
        unseen_mask = dust_map == hp.UNSEEN
        dust_map[unseen_mask] = 0.0
        dust_map = hp.smoothing(dust_map, sigma=np.radians(smoothing_sigma))
        dust_map[unseen_mask] = hp.UNSEEN
    dust_map = np.where(dust_map == hp.UNSEEN, np.nan, dust_map)
    xsize = 720
    ysize = max(260, int(round(xsize * (float(lat_max) - float(lat_min)) / (2.0 * float(delta_lon)))))
    projector = hp.projector.CartesianProj(
        lonra=[float(lon_c) - float(delta_lon), float(lon_c) + float(delta_lon)],
        latra=[float(lat_min), float(lat_max)],
        xsize=xsize,
        ysize=ysize,
    )
    image = projector.projmap(dust_map, lambda x, y, z: hp.vec2pix(lbpub.FIXED_NSIDE, x, y, z))
    image = np.asarray(image, dtype=float)
    image[~np.isfinite(image)] = np.nan
    image[(image < vmin - 1e6) | (image > vmax + 1e6)] = np.nan
    return image


def draw_lb_panel_on_axis(lbpub, ax, bubble_df, hmsfr_df, npix, pix_ids, lon, lat, d_lo, d_hi, *, delta_range, lon_c, delta_lon, lat_min, lat_max, smoothing_sigma, fit_overlay, center_lb, publish_params, lb_info_label, show_ylabel, legend_mode):
    image = build_lb_projected_image(
        lbpub,
        npix,
        pix_ids,
        lon,
        lat,
        d_lo,
        d_hi,
        lon_c=lon_c,
        delta_lon=delta_lon,
        lat_min=lat_min,
        lat_max=lat_max,
        smoothing_sigma=smoothing_sigma,
        vmin=float(delta_range[0]),
        vmax=float(delta_range[1]),
    )
    lb_norm_label = str(publish_params.get("lb_norm", "linear")).strip().lower() if publish_params is not None and not pd.isna(publish_params.get("lb_norm", "linear")) else "linear"
    if lb_norm_label == "log":
        log_vmin = float(delta_range[0])
        if log_vmin <= 0:
            log_vmin = max(float(delta_range[1]) * 1e-3, 1e-6)
        norm = LogNorm(vmin=log_vmin, vmax=float(delta_range[1]))
    else:
        norm = Normalize(vmin=float(delta_range[0]), vmax=float(delta_range[1]))
    ax.imshow(
        image,
        origin="lower",
        cmap="Spectral_r",
        norm=norm,
        extent=[-float(delta_lon), float(delta_lon), float(lat_min), float(lat_max)],
        aspect="auto",
    )
    ax.set_title(f"Distance: [{d_lo:.2f}, {d_hi:.2f}] kpc", fontsize=TITLE_SIZE)
    ax.grid(True, color="black", linestyle="--", alpha=0.3)

    draw_lb_bubbles(
        lbpub,
        ax,
        bubble_df[(bubble_df["best_distance"] >= d_lo) & (bubble_df["best_distance"] <= d_hi)],
        lon_c,
        publish_params,
    )

    if fit_overlay is not None:
        fit_shape = str(fit_overlay.get("fit_shape", "ellipsoid"))
        line_lons = fit_overlay.get("line_lons")
        boundary_df = fit_overlay["boundary_df"]
        projection_mode = str(publish_params.get("lb_projection_mode", "")).strip()
        include_cap = projection_mode in {"distance_mid_slice", "distance_mid_slice_with_cap"}
        if projection_mode in {"distance_mid_slice", "distance_mid_slice_with_cap"} and fit_shape == "cylinder" and line_lons is not None:
            intersects = cylinder_distance_slice_intersects(
                lbpub,
                np.asarray(fit_overlay["center_kpc"], dtype=float),
                np.asarray(fit_overlay["axes_kpc"], dtype=float),
                float(fit_overlay["angle_deg"]),
                0.5 * (float(d_lo) + float(d_hi)),
                float(lat_min),
                float(lat_max),
            )
            if intersects:
                x_values = [float(lb_overlay_x(line_lon, lon_c)) for line_lon in line_lons]
                total_height = 0.5 * (float(lat_max) - float(lat_min))
                half_height = 0.5 * total_height
                y_low = max(float(lat_min), float(center_lb[1]) - half_height)
                y_high = min(float(lat_max), float(center_lb[1]) + half_height)
                for x_value in x_values:
                    (line,) = ax.plot([x_value, x_value], [y_low, y_high], color="black", linewidth=float(publish_params["superbubble_linewidth"]), alpha=0.95, zorder=6)
                    line.set_path_effects(superbubble_path_effects(float(publish_params["superbubble_linewidth"])))
        elif projection_mode in {"distance_mid_slice", "distance_mid_slice_with_cap"}:
            d_mid = 0.5 * (float(d_lo) + float(d_hi))
            curve_segments = build_distance_slice_lb_segments(
                lbpub,
                fit_overlay,
                d_mid,
                lon_c=float(lon_c),
                delta_lon=float(delta_lon),
                lat_min=float(lat_min),
                lat_max=float(lat_max),
                include_cap=include_cap,
                respect_mark=False if fit_shape == "ellipsoid" else True,
            )
            if curve_segments and fit_shape == "ellipsoid":
                cap_curve_segments = build_distance_slice_lb_segments(
                    lbpub,
                    fit_overlay,
                    d_mid,
                    lon_c=float(lon_c),
                    delta_lon=float(delta_lon),
                    lat_min=float(lat_min),
                    lat_max=float(lat_max),
                    include_cap=include_cap,
                    respect_mark=True,
                )
                split_segments = split_lb_segments_by_cap(
                    curve_segments,
                    d_mid,
                    int(fit_overlay.get("mark_value", 0)),
                    float(fit_overlay.get("cap_z_global", np.nan)),
                )
                plot_lb_curve_segments(
                    ax,
                    split_segments["dashed"],
                    lon_c,
                    color="black",
                    linestyle="--",
                    linewidth=max(1.2, 0.7 * float(publish_params["superbubble_linewidth"])),
                    alpha=0.95,
                    zorder=6,
                )
                plot_lb_curve_segments(
                    ax,
                    cap_curve_segments,
                    lon_c,
                    color="black",
                    linestyle="-",
                    linewidth=float(publish_params["superbubble_linewidth"]),
                    alpha=0.95,
                    zorder=7,
                    path_effects=superbubble_path_effects(float(publish_params["superbubble_linewidth"])),
                )
            elif curve_segments:
                plot_lb_curve_segments(
                    ax,
                    curve_segments,
                    lon_c,
                    color="black",
                    linestyle="-",
                    linewidth=float(publish_params["superbubble_linewidth"]),
                    alpha=0.95,
                    zorder=6,
                    path_effects=superbubble_path_effects(float(publish_params["superbubble_linewidth"])),
                )
        elif fit_shape == "cylinder" and line_lons is not None:
            x_values = [float(lb_overlay_x(line_lon, lon_c)) for line_lon in line_lons]
            total_height = 0.5 * (float(lat_max) - float(lat_min))
            half_height = 0.5 * total_height
            y_low = max(float(lat_min), float(center_lb[1]) - half_height)
            y_high = min(float(lat_max), float(center_lb[1]) + half_height)
            for x_value in x_values:
                (line,) = ax.plot([x_value, x_value], [y_low, y_high], color="black", linewidth=float(publish_params["superbubble_linewidth"]), alpha=0.95, zorder=6)
                line.set_path_effects(superbubble_path_effects(float(publish_params["superbubble_linewidth"])))
        else:
            plot_lb_curve_in_distance_bin(
                lbpub,
                ax,
                fit_overlay.get("sky_outline_df"),
                d_lo,
                d_hi,
                lon_c,
                color="black",
                linestyle="-",
                linewidth=float(publish_params["superbubble_linewidth"]),
                alpha=0.95,
                zorder=6,
                filter_distance=False,
                path_effects=superbubble_path_effects(float(publish_params["superbubble_linewidth"])),
            )
        if boundary_df is not None:
            boundary_mask = (boundary_df["d"] >= d_lo) & (boundary_df["d"] <= d_hi)
            if boundary_mask.any():
                ax.scatter(
                    lb_overlay_x(boundary_df.loc[boundary_mask, "l"], lon_c),
                    boundary_df.loc[boundary_mask, "b"],
                    s=0.5 * float(publish_params["cloud_marker_size"]),
                    facecolors="white",
                    edgecolors="black",
                    linewidths=1.6,
                    alpha=0.95,
                    zorder=7,
                )

    ax.plot(lb_overlay_x(center_lb[0], lon_c), center_lb[1], marker="x", color="white", markersize=12, mew=4.2, linestyle="None", zorder=9)
    ax.plot(lb_overlay_x(center_lb[0], lon_c), center_lb[1], marker="x", color="black", markersize=12, mew=2.2, linestyle="None", zorder=10)

    hmsfr_mask = (hmsfr_df["distance_kpc"] >= d_lo) & (hmsfr_df["distance_kpc"] <= d_hi)
    if hmsfr_mask.any():
        ax.scatter(
            lb_overlay_x(hmsfr_df.loc[hmsfr_mask, "l_deg"], lon_c),
            hmsfr_df.loc[hmsfr_mask, "b_deg"],
            s=float(publish_params["hmsfr_marker_size"]),
            marker="*",
            c="#00d8ff",
            edgecolors="black",
            linewidths=1.8,
            zorder=8,
        )

    ax.set_xlim(-float(delta_lon), float(delta_lon))
    ax.set_ylim(float(lat_min), float(lat_max))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Galactic Longitude l (°)", fontsize=LABEL_SIZE)
    if show_ylabel:
        ax.set_ylabel("Galactic Latitude b (°)", fontsize=LABEL_SIZE)
    else:
        ax.set_yticklabels([])
    ax.tick_params(labelsize=TICK_SIZE)
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda value, _: lb_tick_label(value, lon_c)))
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda value, _: f"{int(round(value))}°"))

    if legend_mode == "markers":
        handles = make_catalog_handles(publish_params)
        legend = ax.legend(handles=handles, loc="upper right", fontsize=LEGEND_SIZE, framealpha=0.94)
        legend.set_zorder(1000)
    elif legend_mode == "info":
        legend = ax.legend(
            handles=[Line2D([], [], linestyle="None", label=lb_info_label)],
            loc="upper right",
            fontsize=LEGEND_SIZE,
            framealpha=0.94,
            handlelength=0,
            handletextpad=0,
        )
        legend.set_zorder(1000)
    elif legend_mode == "superbubble":
        legend = ax.legend(handles=[make_superbubble_handle(float(publish_params["superbubble_linewidth"]))], loc="upper right", fontsize=LEGEND_SIZE, framealpha=0.94)
        legend.set_zorder(1000)


def render_publication_six_panel(viz, lbpub, seed_row, fit_data, result_data, bubble_df, hmsfr_df, xy_data_path, output_path, publish_params):
    center = np.asarray(fit_data["center_kpc"], dtype=float)
    axes_fit = np.array(
        [fit_data["axes_fit_kpc"]["a"], fit_data["axes_fit_kpc"]["b"], fit_data["axes_fit_kpc"]["c"]],
        dtype=float,
    )
    fit_shape = str(fit_data["fit_shape"]).lower()
    mark = int(fit_data.get("effective_mark_value", result_data.get("mark", seed_row.get("mark", 0))))
    selected_clouds = pd.read_csv(package_path(fit_data["outputs"]["boundary_points_csv"]))
    cap_z_global = compute_cap_z_from_selected_clouds(selected_clouds, center[2], mark)
    # Geometry and shell-mask conventions used by this analysis stage.

    if bool(result_data.get("manual_override", False)):
        cap_z_global = float(center[2])

    requested_x_range, requested_y_range = build_xy_panel_ranges(center, axes_fit, seed_row, publish_params)
    z_slices = build_xy_z_slices_from_params(publish_params)
    z_data_range = [min(z[0] for z in z_slices), max(z[1] for z in z_slices)]
    df_xyz = viz.read_local_xyz_cube(xy_data_path, requested_x_range, requested_y_range, z_data_range)
    x_range, y_range = clip_xy_range_to_available_cube(df_xyz, requested_x_range, requested_y_range)
    z_step = viz.estimate_axis_step(df_xyz["z"])
    xy_images = [viz.extract_xy_grid(df_xyz, x_range, y_range, z_slice, 1.2, z_step)[0] for z_slice in z_slices]

    seed_row = apply_lb_window_overrides(seed_row, publish_params)
    fit_overlay = build_fit_overlay_for_lb(
        lbpub,
        seed_row,
        fit_data,
        projection_min_distance_kpc=float(publish_params.get("lb_projection_min_distance_kpc", 0.0)),
    )
    distance_slices = build_lb_distance_slices(publish_params)
    default_lon_c = float(seed_row["lon_c"])
    default_delta_lon = float(seed_row["delta_lon"])
    default_lat_min = float(publish_params["lb_lat_min"])
    default_lat_max = float(publish_params["lb_lat_max"])
    panel_windows = [
        build_lb_panel_window(
            publish_params,
            index,
            default_lon_c,
            default_delta_lon,
            default_lat_min,
            default_lat_max,
        )
        for index in range(3)
    ]
    smoothing_sigma = float(seed_row.get("sigma_deg", 0.1))
    center_coord = SkyCoord(
        CartesianRepresentation(center[0] * u.kpc, center[1] * u.kpc, center[2] * u.kpc),
        frame="galactic",
    )
    center_lb = (float(center_coord.l.deg), float(center_coord.b.deg), float(center_coord.distance.kpc))
    if len(distance_slices) >= 2:
        mid_d_lo, mid_d_hi = distance_slices[1]
        d_mid = 0.5 * (float(mid_d_lo) + float(mid_d_hi))
        mid_window = panel_windows[1]
        projection_mode = str(publish_params.get("lb_projection_mode", "")).strip()
        include_cap = projection_mode in {"distance_mid_slice", "distance_mid_slice_with_cap"}
        if str(fit_overlay.get("fit_shape", "")).lower() == "cylinder" and fit_overlay.get("line_lons") is not None:
            intersects = cylinder_distance_slice_intersects(
                lbpub,
                np.asarray(fit_overlay["center_kpc"], dtype=float),
                np.asarray(fit_overlay["axes_kpc"], dtype=float),
                float(fit_overlay["angle_deg"]),
                d_mid,
                float(mid_window["lat_min"]),
                float(mid_window["lat_max"]),
            )
            if intersects:
                x_values = [float(lb_overlay_x(line_lon, mid_window["lon_c"])) for line_lon in fit_overlay["line_lons"]]
                lon_span_deg = abs(float(max(x_values) - min(x_values)))
                lat_span_deg = 0.5 * abs(float(mid_window["lat_max"]) - float(mid_window["lat_min"]))
            else:
                lon_span_deg, lat_span_deg = 0.0, 0.0
        else:
            mid_segments = build_distance_slice_lb_segments(
                lbpub,
                fit_overlay,
                d_mid,
                lon_c=float(mid_window["lon_c"]),
                delta_lon=float(mid_window["delta_lon"]),
                lat_min=float(mid_window["lat_min"]),
                lat_max=float(mid_window["lat_max"]),
                include_cap=include_cap,
            )
            lon_span_deg, lat_span_deg = lb_curve_segments_metrics(mid_segments, mid_window["lon_c"])
    else:
        lon_span_deg, lat_span_deg = sky_outline_metrics(lbpub, fit_overlay, center_lb, default_lon_c)
    lb_info_label = (
        f"l, b, d = ({center_lb[0]:.1f}°, {center_lb[1]:.1f}°, {center_lb[2]:.2f} kpc)\n"
        f"longitude span={lon_span_deg:.1f}°\n"
        f"latitude span={lat_span_deg:.1f}°"
    )
    npix = hp.nside2npix(lbpub.FIXED_NSIDE)
    pix_ids_all = np.arange(npix)
    lon_all, lat_all = hp.pix2ang(lbpub.FIXED_NSIDE, pix_ids_all, nest=False, lonlat=True)

    fig = plt.figure(figsize=(18, 12), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        3,
        left=0.06,
        right=0.90,
        bottom=0.07,
        top=0.95,
        wspace=0.16,
        hspace=0.30,
    )
    xy_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    lb_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]

    xy_norm = Normalize(vmin=float(publish_params["xy_vmin"]), vmax=float(publish_params["xy_vmax"]))
    lb_norm_label = str(publish_params.get("lb_norm", "linear")).strip().lower() if not pd.isna(publish_params.get("lb_norm", "linear")) else "linear"
    if lb_norm_label == "log":
        # LogNorm requires strictly positive bounds; clamp vmin if user left it at 0.
        log_vmin = float(publish_params["lb_vmin"])
        if log_vmin <= 0:
            log_vmin = max(float(publish_params["lb_vmax"]) * 1e-3, 1e-6)
        lb_norm = LogNorm(vmin=log_vmin, vmax=float(publish_params["lb_vmax"]))
    else:
        lb_norm = Normalize(vmin=float(publish_params["lb_vmin"]), vmax=float(publish_params["lb_vmax"]))
    xy_handles = make_catalog_handles(publish_params)
    xy_reference_radii = None
    if len(distance_slices) >= 2:
        xy_reference_radii = [float(distance_slices[1][0]), float(distance_slices[1][1])]
    for index, ax in enumerate(xy_axes):
        image = xy_images[index]
        ax.imshow(image, origin="lower", cmap="Spectral_r", norm=xy_norm, extent=[*x_range, *y_range])
        n_contours = int(publish_params["xy_contour_levels"])
        contour_levels = np.linspace(float(publish_params["xy_vmin"]), float(publish_params["xy_vmax"]), n_contours + 2)[1:-1]
        if image is not None and np.ma.count(image) > 0:
            ax.contour(image, levels=contour_levels, origin="lower", extent=[*x_range, *y_range], colors="black", linewidths=0.55, alpha=0.45, zorder=4)
        z_mid = 0.5 * (z_slices[index][0] + z_slices[index][1])
        section_axes = section_axes_for_cap(axes_fit, center[2], z_mid, mark, fit_shape, cap_z_global)
        ellipse_style = "-"
        ellipse_path_effects = superbubble_path_effects(float(publish_params["superbubble_linewidth"]))
        if section_axes is None and z_slice_is_beyond_cap(z_mid, mark, cap_z_global):
            section_axes = cap_plane_section_axes(axes_fit, center[2], fit_shape, cap_z_global)
            ellipse_style = "-"
        if section_axes is not None:
            ellipse = Ellipse(
                (center[0], center[1]),
                width=2.0 * section_axes[0],
                height=2.0 * section_axes[1],
                angle=float(fit_data["angle_deg_fixed"]),
                edgecolor="black",
                facecolor="none",
                linestyle=ellipse_style,
                linewidth=float(publish_params["superbubble_linewidth"]),
                zorder=6,
            )
            if ellipse_path_effects is not None:
                ellipse.set_path_effects(ellipse_path_effects)
            ax.add_patch(ellipse)
        draw_manual_extra_xy_ellipse(ax, z_slices[index], publish_params)
        draw_solar_xy_grid(ax, x_range, y_range, reference_radii=xy_reference_radii)
        draw_xy_overlays(ax, z_slices[index], selected_clouds, bubble_df, hmsfr_df, publish_params)
        ax.plot(center[0], center[1], marker="x", color="black", markersize=8, mew=2.0, zorder=9)
        ax.set_xlim(x_range)
        ax.set_ylim(y_range)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, linestyle="--", alpha=0.22)
        ax.set_title(f"Z: [{z_slices[index][0]:.2f}, {z_slices[index][1]:.2f}] kpc", fontsize=TITLE_SIZE)
        ax.set_xlabel("X [kpc]", fontsize=LABEL_SIZE)
        ax.tick_params(labelsize=TICK_SIZE)
        if index == 0:
            legend = ax.legend(handles=[make_superbubble_handle(float(publish_params["superbubble_linewidth"]))], loc="upper right", fontsize=LEGEND_SIZE, framealpha=0.94)
            legend.set_zorder(1000)
        elif index == 1:
            info_label = (
                f"X, Y, Z = ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}) kpc\n"
                f"a,b,c=({axes_fit[0] * 1000:.0f},{axes_fit[1] * 1000:.0f},{axes_fit[2] * 1000:.0f}) pc\n"
                f"PA={float(fit_data['angle_deg_fixed']):.0f}°"
            )
            legend = ax.legend(handles=[Line2D([], [], linestyle="None", label=info_label)], loc="upper right", fontsize=LEGEND_SIZE, framealpha=0.94, handlelength=0, handletextpad=0)
            legend.set_zorder(1000)
        elif index == 2:
            legend = ax.legend(handles=xy_handles, loc="upper right", fontsize=LEGEND_SIZE, framealpha=0.94)
            legend.set_zorder(1000)
    xy_axes[0].set_ylabel("Y [kpc]", fontsize=LABEL_SIZE)
    for index, ax in enumerate(lb_axes):
        d_lo, d_hi = distance_slices[index]
        panel_window = panel_windows[index]
        mask = lbpub.build_pixel_mask(
            lon_all,
            lat_all,
            float(panel_window["lon_c"]),
            float(panel_window["delta_lon"]),
            0.5 * (float(panel_window["lat_min"]) + float(panel_window["lat_max"])),
            0.5 * (float(panel_window["lat_max"]) - float(panel_window["lat_min"])),
        )
        pix_ids = pix_ids_all[mask]
        lon_masked = lon_all[mask]
        lat_masked = lat_all[mask]
        draw_lb_panel_on_axis(
            lbpub,
            ax,
            bubble_df,
            hmsfr_df,
            npix,
            pix_ids,
            lon_masked,
            lat_masked,
            d_lo,
            d_hi,
            delta_range=(float(publish_params["lb_vmin"]), float(publish_params["lb_vmax"])),
            lon_c=float(panel_window["lon_c"]),
            delta_lon=float(panel_window["delta_lon"]),
            lat_min=float(panel_window["lat_min"]),
            lat_max=float(panel_window["lat_max"]),
            smoothing_sigma=smoothing_sigma,
            fit_overlay=fit_overlay,
            center_lb=center_lb,
            publish_params=publish_params,
            lb_info_label=lb_info_label,
            show_ylabel=(index == 0),
            legend_mode="superbubble" if index == 0 else ("info" if index == 1 else "markers"),
        )

    fig.canvas.draw()
    for axes_row, norm in [(xy_axes, xy_norm), (lb_axes, lb_norm)]:
        reference_ax = axes_row[2]
        reference_box = reference_ax.get_position()
        cax = fig.add_axes([reference_box.x1 + 0.012, reference_box.y0, 0.018, reference_box.height])
        cbar = mpl.colorbar.ColorbarBase(cax, cmap="Spectral_r", norm=norm, orientation="vertical")
        cbar.set_label(r"$\Delta E(B-V)$ [mag]", fontsize=LABEL_SIZE)
        cbar.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
        cbar.ax.tick_params(labelsize=TICK_SIZE)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


# ============================================================================
# ============================================================================
viz = types.SimpleNamespace(
    estimate_axis_step=estimate_axis_step,
    extract_xy_grid=extract_xy_grid,
    normalize_xyz_dataframe=normalize_xyz_dataframe,
    read_local_xyz_cube=read_local_xyz_cube,
)
lbpub = types.SimpleNamespace(
    FIXED_NSIDE=FIXED_NSIDE,
    build_pixel_mask=build_pixel_mask,
    query_shell_delta_ebv=query_shell_delta_ebv,
    rotate_points_to_global=rotate_points_to_global,
    wrap_longitude_to_center=wrap_longitude_to_center,
)


def package_relative(path: Path) -> str:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    path = Path(path)
    return path.relative_to(HERE).as_posix() if path.is_absolute() else path.as_posix()


def package_path(path_value) -> Path:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    path = Path(path_value)
    return path


def build_publication_result_data(row, result_json, fit_table, source_json):
    """Compute an intermediate statistic or table used by this documented workflow."""
    sb = int(row["id"])
    fit_shape = str(row["shape"]).lower()
    c_axis = float(row["least_squares_c_kpc"]) if fit_shape == "cylinder" else float(row["c_radius_kpc"])
    center = [float(row["center_x_kpc"]), float(row["center_y_kpc"]), float(row["center_z_kpc"])]
    axes = [float(row["a_radius_kpc"]), float(row["b_radius_kpc"]), c_axis]
    result_data = dict(source_json)
    outputs = dict(result_data.get("outputs", {}))
    # Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
    outputs["selected_clouds_csv"] = f"../results/intermediate_output/2_automated_fit_results/json/SB{sb}_clouds.csv"
    outputs["json"] = f"../results/intermediate_output/2_automated_fit_results/json/SB{sb}_fit.json"
    result_data.update({
        "target_index": sb,
        "version": str(row.get("version", "peer_review_v4")),
        "mark": int(row["mark"]),
        "fit_shape": fit_shape,
        "center_kpc": center,
        "fit_a_kpc": axes[0],
        "fit_b_kpc": axes[1],
        "fit_c_kpc": axes[2],
        "fit_axes_kpc": axes,
        "angle_deg": float(row["angle_deg"]),
        "publication_geometry_source": package_relative(fit_table) if Path(fit_table).is_absolute() else str(fit_table),
        "outputs": outputs,
    })
    return result_data


def default_json_paths(bestfit_df):
    return {int(row["id"]): f"../results/intermediate_output/2_automated_fit_results/json/SB{int(row['id'])}_fit.json" for _, row in bestfit_df.iterrows()}


def render_six_panels(bestfit_df, json_paths=None, fit_table=BESTFIT_CSV):
    """Render or save a documented figure product without changing upstream measurements."""
    PANEL_DIR.mkdir(parents=True, exist_ok=True)

    if json_paths is None:
        json_paths = default_json_paths(bestfit_df)

    bubble_df = prepare_bubble_catalog(BUBBLE_CSV)
    hmsfr_df = prepare_hmsfr_catalog(HMSFR_CSV)

    for _, row in bestfit_df.sort_values("id").iterrows():
        sb = int(row["id"])
        json_path = package_path(json_paths[sb])
        if not json_path.exists():
            raise FileNotFoundError('Missing input file.')
        source_json = json.loads(json_path.read_text(encoding="utf-8"))
        result_data = build_publication_result_data(row, json_path, fit_table, source_json)
        fit_data = build_fit_data(result_data)
        seed_row = load_seed_row(PLOTTING_PARAM_CSV, sb)
        publish_params = ensure_publish_param_row(PUBLISH_PARAM_CSV, sb)
        validate_publish_params_for_rendering(publish_params)
        out_png = PANEL_DIR / f"SB{sb}.png"
        render_publication_six_panel(
            viz, lbpub, seed_row, fit_data, result_data,
            bubble_df, hmsfr_df, XY_DUST_PATH, out_png, publish_params,
        )
        print(f"  six-panel figure: {out_png}")


def render_manual_six_panels_if_available():
    """Render or save a documented figure product without changing upstream measurements."""
    if not FINAL_TABLE_CSV.exists():
        print(f"Manual-panel rendering skipped: missing final table {FINAL_TABLE_CSV}")
        return 0

    final_df = pd.read_csv(FINAL_TABLE_CSV, encoding="utf-8-sig")
    manual_df = final_df[final_df["id"].astype(int).isin(MANUAL_IDS)].copy()
    if manual_df.empty:
        print("Manual-panel rendering skipped: no manual OSBs are present in the final table.")
        return 0

    json_paths = default_json_paths(manual_df)
    missing = [sb for sb, path in json_paths.items() if not package_path(path).exists()]
    if missing:
        print(f"Manual-panel rendering skipped: missing JSON fits for SB IDs {missing}.")
        return 0

    print(f"Rendering manual six-panel figure for {len(manual_df)} OSBs.")
    render_six_panels(manual_df, json_paths=json_paths, fit_table=FINAL_TABLE_CSV)
    return int(len(manual_df))


def main():
    """Run Script 4 from validated inputs to the documented outputs."""
    if not BESTFIT_CSV.exists():
        raise FileNotFoundError('Missing input file.')

    bestfit = pd.read_csv(BESTFIT_CSV, encoding="utf-8-sig")
    print(f"Rendering publication six-panel figures from {BESTFIT_CSV}")
    render_six_panels(bestfit)
    render_manual_six_panels_if_available()
    print("Script 4 figure rendering complete.")


def format_duration(seconds: float) -> str:
    seconds = max(0.0, float(seconds))
    hours, remainder = divmod(seconds, 3600.0)
    minutes, seconds = divmod(remainder, 60.0)
    if hours >= 1.0:
        return f"{int(hours)} h {int(minutes):02d} min {seconds:05.2f} s"
    if minutes >= 1.0:
        return f"{int(minutes)} min {seconds:05.2f} s"
    return f"{seconds:.2f} s"


def run_timed_main(func, script_file=None):
    """Run this timed entry point and report elapsed wall-clock time."""
    script_name = Path(script_file or __file__).name
    start = time.perf_counter()
    try:
        result = func()
    except SystemExit as exc:
        elapsed = time.perf_counter() - start
        if exc.code in (None, 0):
            print(f"[runtime] {script_name} completed in: {format_duration(elapsed)} ({elapsed:.2f} s)", flush=True)
        else:
            print(f"[runtime] {script_name} elapsed time before failure: {format_duration(elapsed)} ({elapsed:.2f} s)", flush=True)
        raise
    except BaseException:
        elapsed = time.perf_counter() - start
        print(f"[runtime] {script_name} elapsed time before failure: {format_duration(elapsed)} ({elapsed:.2f} s)", flush=True)
        raise
    elapsed = time.perf_counter() - start
    print(f"[runtime] {script_name} completed in: {format_duration(elapsed)} ({elapsed:.2f} s)", flush=True)
    return result


if __name__ == "__main__":
    run_timed_main(main, __file__)
