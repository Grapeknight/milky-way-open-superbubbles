#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script 18: YSO OSBs shell random test

Purpose
-------
Test young-cluster age windows inside open-superbubble volumes and on dust shells with longitude randomization.

Method overview
---------------
1. Build or load the young Hunt-cluster sample and read the final superbubble
   geometries.
2. Count clusters in the configured age windows within the adopted volumes and dust
   shells.
3. Randomize Galactic longitude while retaining the other cluster properties, compare
   the observed counts with the controls, and export per-bubble and overall summaries
   and figures.

Main inputs
-----------
- ../data/star_cluster_data/hunt2024_clusters_full.csv
- ../results/intermediate_output/18_young_cluster_superbubble_shell_random_test/Hunt_young_cluster_sample.csv
- ../results/superbubble_final_fit_parameters.csv

Main outputs
------------
- ../results/figures/18_hunt_young_clusters_shell_relative_offset_trend.png; Extended Data age-stratification figure.

Figure/table role
-----------------
../results/figures/18_hunt_young_clusters_shell_relative_offset_trend.png; Extended Data age-stratification figure.

Runtime and data notes
----------------------
Reference runtime: 25.27 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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
import os
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
RAW_CLUSTER_DIR = DATA_DIR / "star_cluster_data"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "18_young_cluster_superbubble_shell_random_test"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"

HUNT_RAW_CSV = RAW_CLUSTER_DIR / "hunt2024_clusters_full.csv"
HUNT_CLUSTER_CSV = OUT_DIR / "Hunt_young_cluster_sample.csv"
HUNT_FILTER_SUMMARY_CSV = OUT_DIR / "Hunt_young_cluster_filter_summary.csv"
CLUSTER_CSV = HUNT_CLUSTER_CSV
FINAL_SB_CSV = FINAL_DIR / "superbubble_final_fit_parameters.csv"

OUT_SUMMARY_CSV = OUT_DIR / "18_superbubble_volume_shell_random_test_summary.csv"
OUT_RANDOM_CSV = OUT_DIR / "18_superbubble_volume_shell_random_trial_details.csv"
OUT_MEMBERS_CSV = OUT_DIR / "18_superbubble_volume_shell_member_details.csv"
OUT_PLOT_DATA_CSV = OUT_DIR / "18_superbubble_volume_shell_combined_plot_data.csv"
OUT_REPORT_MD = OUT_DIR / "18_superbubble_volume_shell_random_test_notes.md"
OUT_HIST_FIG = FINAL_FIG_DIR / "18_hunt_young_clusters_shell_random_hist_matrix.png"
OUT_LINE_FIG = FINAL_FIG_DIR / "18_hunt_young_clusters_shell_relative_offset_trend.png"

MIN_DIST_PC = 200.0
MAX_DIST_PC = 3000.0
SHELL_U_MIN = 0.9
SHELL_U_MAX = 1.2
DUST_CYLINDER_HALF_HEIGHT_KPC = 0.025
DEFAULT_FULL_CYLINDER_HALF_HEIGHT_KPC = 0.100
HALF_ELLIPSOID_SB_IDS = frozenset({2, 11, 12, 15, 16, 17, 22, 30, 32})
DEFAULT_N_RANDOM = 1000
DEFAULT_SEED = 20260526
DEFAULT_RANDOM_BATCH_SIZE = 50
HUNT_MIN_CST = 4.0
HUNT_MIN_CMDCL50 = 0.3
HUNT_XY_LIMIT_PC = 3000.0
AGE_WINDOW_MARGIN_MYR = 5.0

AGE_BINS = (
    ("0-10", "[0, 10]", 0.0, 10.0, True),
    ("10-20", "(10, 20]", 10.0, 20.0, False),
    ("20-30", "(20, 30]", 20.0, 30.0, False),
    ("30-40", "(30, 40]", 30.0, 40.0, False),
    ("40-50", "(40, 50]", 40.0, 50.0, False),
    ("50-60", "(50, 60]", 50.0, 60.0, False),
)
AGE_LABELS = [item[0] for item in AGE_BINS]
AGE_WINDOWS = {
    label: (max(0.0, lower - AGE_WINDOW_MARGIN_MYR), upper + AGE_WINDOW_MARGIN_MYR)
    for label, _interval, lower, upper, _include_lower in AGE_BINS
}
AGE_MIDPOINTS = {label: (window_lower + window_upper) / 2.0 for label, (window_lower, window_upper) in AGE_WINDOWS.items()}
AGE_TICK_LABELS = [f"{window_lower:g}-{window_upper:g}" for window_lower, window_upper in AGE_WINDOWS.values()]
AGE_TICK_POSITIONS = [AGE_MIDPOINTS[label] for label in AGE_LABELS]
HISTOGRAM_AGE_LABELS = AGE_LABELS

REGIONS = (
    ("cap_interior", "Final superbubble cap/cylinder interiors", "inside_cap_interior", "#377EB8"),
    ("dust_shell", "Final-sample dust shells", "inside_dust_shell", "#984EA3"),
)

PLOT_REGION_LABELS = {
    "cap_interior": "Open superbubble interiors",
    "dust_shell": "Open superbubble shells",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Noto Sans SC", "Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "font.size": 8.8,
        "axes.labelsize": 9.5,
        "axes.titlesize": 8.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "savefig.dpi": 300,
    }
)


def ensure_existing_file(path: Path | str, label: str) -> Path:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{label}  does not exist: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"{label}  is not a file: {path}")
    return path


def require_columns(df: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{label}  is missing required columns: {', '.join(missing)}")


def relative_path_text(path: Path | str) -> str:
    path = Path(path)
    try:
        return str(path.resolve().relative_to(HERE))
    except ValueError:
        return str(path)


def rotation_matrix_z(angle_deg: float) -> np.ndarray:
    theta = np.deg2rad(float(angle_deg))
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rotate_points_to_local(points: np.ndarray, angle_deg: float) -> np.ndarray:
    return np.asarray(points, dtype=float) @ rotation_matrix_z(-float(angle_deg)).T


def build_hunt_cluster_table(raw_path: Path, output_path: Path, summary_path: Path) -> pd.DataFrame:
    df = pd.read_csv(raw_path).copy()
    required = ["recno", "Name", "ID", "CST", "CMDCl50", "GLON", "GLAT", "dist50", "logAge50"]
    require_columns(df, required, str(raw_path))

    numeric = ["recno", "CST", "CMDCl50", "GLON", "GLAT", "dist50", "logAge50"]
    for column in numeric:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    stages: list[tuple[str, int]] = [("raw Hunt catalogue", len(df))]
    current = df.dropna(subset=numeric).copy()
    stages.append(('Selection stage', len(current)))

    lon = np.deg2rad(current["GLON"].to_numpy(dtype=float))
    lat = np.deg2rad(current["GLAT"].to_numpy(dtype=float))
    dist_pc = current["dist50"].to_numpy(dtype=float)
    cos_lat = np.cos(lat)
    current["age_myr"] = np.power(10.0, current["logAge50"].to_numpy(dtype=float) - 6.0)
    current["dist_pc"] = dist_pc
    current["x_helio"] = dist_pc * cos_lat * np.cos(lon)
    current["y_helio"] = dist_pc * cos_lat * np.sin(lon)
    current["z_helio"] = dist_pc * np.sin(lat)

    current = current[current["CST"] > HUNT_MIN_CST].copy()
    stages.append((f"CST > {HUNT_MIN_CST:g}", len(current)))
    current = current[current["CMDCl50"] > HUNT_MIN_CMDCL50].copy()
    stages.append((f"CMDCl50 > {HUNT_MIN_CMDCL50:g}", len(current)))
    current = current[current["x_helio"].abs() < HUNT_XY_LIMIT_PC].copy()
    current = current[current["y_helio"].abs() < HUNT_XY_LIMIT_PC].copy()
    stages.append(('Selection stage', len(current)))

    name = current["Name"].astype(str).str.strip()
    fallback_name = current["ID"].astype(str).str.strip()
    current["name"] = np.where(name.ne("") & name.ne("nan"), name, fallback_name)
    current["row_id"] = np.arange(1, len(current) + 1, dtype=int)
    current["source_catalog"] = "hunt2024_raw_filtered"
    current["hunt_recno"] = current["recno"]
    current["hunt_id"] = current["ID"]

    output_columns = [
        "row_id",
        "name",
        "source_catalog",
        "hunt_recno",
        "hunt_id",
        "age_myr",
        "dist_pc",
        "x_helio",
        "y_helio",
        "z_helio",
        "CST",
        "CMDCl50",
        "GLON",
        "GLAT",
        "dist50",
        "logAge50",
    ]
    output = current.sort_values(["dist_pc", "row_id"], kind="stable")[output_columns].reset_index(drop=True)
    output["row_id"] = np.arange(1, len(output) + 1, dtype=int)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        output.to_csv(output_path, index=False, encoding="utf-8-sig")
    except PermissionError as exc:
        print(f"Warning: cannot write locked cache file {output_path}: {exc}")
    try:
        pd.DataFrame(stages, columns=["stage", "count"]).to_csv(summary_path, index=False, encoding="utf-8-sig")
    except PermissionError as exc:
        print(f"Warning: cannot write locked cache file {summary_path}: {exc}")
    return output


def load_clusters(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path).copy()
    require_columns(df, ["row_id", "name", "age_myr", "dist_pc", "x_helio", "y_helio", "z_helio"], str(path))
    numeric = ["age_myr", "dist_pc", "x_helio", "y_helio", "z_helio"]
    for column in numeric:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    valid = df.dropna(subset=numeric).copy()
    selected = valid[(valid["dist_pc"] > MIN_DIST_PC) & (valid["dist_pc"] <= MAX_DIST_PC)].copy()
    selected["x_kpc"] = selected["x_helio"] / 1000.0
    selected["y_kpc"] = selected["y_helio"] / 1000.0
    selected["z_kpc"] = selected["z_helio"] / 1000.0
    return selected.sort_values(["dist_pc", "row_id"]).reset_index(drop=True)


def age_window_bounds(lower: float, upper: float) -> tuple[float, float]:
    return max(0.0, float(lower) - AGE_WINDOW_MARGIN_MYR), float(upper) + AGE_WINDOW_MARGIN_MYR


def age_window_text(lower: float, upper: float) -> str:
    window_lower, window_upper = age_window_bounds(lower, upper)
    return f"[{window_lower:g}, {window_upper:g}]"


def age_window_mask(df: pd.DataFrame, lower: float, upper: float, include_lower: bool) -> pd.Series:
    window_lower, window_upper = age_window_bounds(lower, upper)
    lower_mask = df["age_myr"] >= window_lower if include_lower else df["age_myr"] > window_lower
    return lower_mask & (df["age_myr"] <= window_upper)


def matching_age_windows(age_myr: float) -> str:
    if not np.isfinite(age_myr):
        return ""
    labels = [
        label
        for label, _interval, lower, upper, _include_lower in AGE_BINS
        if bool(age_window_mask(pd.DataFrame({"age_myr": [float(age_myr)]}), lower, upper, _include_lower).iloc[0])
    ]
    return ";".join(labels)


def assign_age_bins(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["age_bin_myr"] = ""
    out["age_interval"] = ""
    out["age_window_bins_myr"] = [matching_age_windows(value) for value in out["age_myr"].to_numpy(dtype=float)]
    for label, interval, lower, upper, include_lower in AGE_BINS:
        lower_mask = out["age_myr"] >= lower if include_lower else out["age_myr"] > lower
        mask = lower_mask & (out["age_myr"] <= upper)
        out.loc[mask, "age_bin_myr"] = label
        out.loc[mask, "age_interval"] = interval
    return out


def load_final_superbubbles(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path).copy()
    required = [
        "id",
        "mark",
        "shape",
        "center_x_kpc",
        "center_y_kpc",
        "center_z_kpc",
        "a_radius_kpc",
        "b_radius_kpc",
        "c_radius_kpc",
        "least_squares_c_kpc",
        "angle_deg",
        "xy_plane_z_kpc",
        "xy_plane_a_kpc",
        "xy_plane_b_kpc",
        "xy_plane_angle_deg",
    ]
    require_columns(df, required, str(path))
    numeric = [column for column in required if column != "shape"]
    for column in numeric:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df["shape"] = df["shape"].astype(str).str.lower()
    df = df.dropna(subset=["id", "mark", "shape", "center_x_kpc", "center_y_kpc", "center_z_kpc"]).copy()
    df["id"] = df["id"].astype(int)
    df["mark"] = df["mark"].astype(int)

    for column in ["a_radius_kpc", "b_radius_kpc", "angle_deg", "xy_plane_z_kpc", "xy_plane_a_kpc", "xy_plane_b_kpc"]:
        df = df[np.isfinite(df[column])].copy()
    df = df[(df["a_radius_kpc"] > 0.0) & (df["b_radius_kpc"] > 0.0)].copy()

    cylinder = df["shape"] == "cylinder"
    ellipsoid = df["shape"] == "ellipsoid"
    if (~(cylinder | ellipsoid)).any():
        bad = sorted(df.loc[~(cylinder | ellipsoid), "shape"].unique())
        raise ValueError('Invalid input or missing required data.')
    if (ellipsoid & ~(df["c_radius_kpc"] > 0.0)).any():
        bad_ids = df.loc[ellipsoid & ~(df["c_radius_kpc"] > 0.0), "id"].tolist()
        raise ValueError('Invalid input or missing required data.')

    df["volume_c_kpc"] = df["c_radius_kpc"]
    use_ls_c = cylinder & np.isfinite(df["least_squares_c_kpc"]) & (df["least_squares_c_kpc"] > 0.0)
    df.loc[cylinder, "volume_c_kpc"] = DEFAULT_FULL_CYLINDER_HALF_HEIGHT_KPC
    df.loc[use_ls_c, "volume_c_kpc"] = df.loc[use_ls_c, "least_squares_c_kpc"]
    df["xy_plane_angle_deg"] = df["xy_plane_angle_deg"].fillna(df["angle_deg"])
    return df.sort_values("id").reset_index(drop=True)


def ellipsoid_cap_u_and_mask(points: np.ndarray, row: object) -> tuple[np.ndarray, np.ndarray]:
    mark = int(row.mark)
    cap_z = float(row.xy_plane_z_kpc)
    c_full = float(row.c_radius_kpc)
    apex_z = float(row.center_z_kpc - c_full) if mark == 1 else float(row.center_z_kpc + c_full)
    c_cap = abs(cap_z - apex_z)
    if c_cap <= 0.0:
        return np.full(len(points), np.nan, dtype=float), np.zeros(len(points), dtype=bool)

    metric_center = np.array([row.center_x_kpc, row.center_y_kpc, cap_z], dtype=float)
    local = rotate_points_to_local(points - metric_center[None, :], float(row.xy_plane_angle_deg))
    axes = np.array([row.xy_plane_a_kpc, row.xy_plane_b_kpc, c_cap], dtype=float)
    u = np.sqrt(np.sum((local / axes[None, :]) ** 2, axis=1))
    if mark == 1:
        cap_side_mask = points[:, 2] <= cap_z + 1e-12
    elif mark == 2:
        cap_side_mask = points[:, 2] >= cap_z - 1e-12
    else:
        cap_side_mask = np.ones(len(points), dtype=bool)
    return u, cap_side_mask


def ellipsoid_half_u_and_mask(points: np.ndarray, row: object) -> tuple[np.ndarray, np.ndarray]:
    center = np.array([row.center_x_kpc, row.center_y_kpc, row.center_z_kpc], dtype=float)
    local = rotate_points_to_local(points - center[None, :], float(row.angle_deg))
    axes = np.array([row.a_radius_kpc, row.b_radius_kpc, row.c_radius_kpc], dtype=float)
    u = np.sqrt(np.sum((local / axes[None, :]) ** 2, axis=1))
    mark = int(row.mark)
    if mark == 1:
        half_side_mask = local[:, 2] <= 1e-12
    elif mark == 2:
        half_side_mask = local[:, 2] >= -1e-12
    else:
        half_side_mask = np.ones(len(points), dtype=bool)
    return u, half_side_mask


def geometry_model_label(row: object) -> str:
    if row.shape == "ellipsoid" and int(row.id) in HALF_ELLIPSOID_SB_IDS:
        return "half ellipsoid"
    if row.shape == "ellipsoid":
        return "ellipsoid cap"
    return "cylinder"


def cap_interior_membership_matrix(points: np.ndarray, final_df: pd.DataFrame) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    membership = np.zeros((len(points), len(final_df)), dtype=bool)
    for idx, row in enumerate(final_df.itertuples(index=False)):
        center = np.array([row.center_x_kpc, row.center_y_kpc, row.center_z_kpc], dtype=float)
        local = rotate_points_to_local(points - center[None, :], float(row.angle_deg))
        if row.shape == "ellipsoid":
            if int(row.id) in HALF_ELLIPSOID_SB_IDS:
                u, side_mask = ellipsoid_half_u_and_mask(points, row)
            else:
                u, side_mask = ellipsoid_cap_u_and_mask(points, row)
            membership[:, idx] = side_mask & (u <= 1.0)
        elif row.shape == "cylinder":
            q_xy = (local[:, 0] / float(row.a_radius_kpc)) ** 2 + (local[:, 1] / float(row.b_radius_kpc)) ** 2
            membership[:, idx] = (q_xy <= 1.0) & (np.abs(local[:, 2]) <= float(row.volume_c_kpc))
    return membership


def dust_shell_membership_matrix(points: np.ndarray, final_df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    points = np.asarray(points, dtype=float)
    membership = np.zeros((len(points), len(final_df)), dtype=bool)
    u_values = np.full((len(points), len(final_df)), np.nan, dtype=float)

    for idx, row in enumerate(final_df.itertuples(index=False)):
        mark = int(row.mark)
        if row.shape == "cylinder":
            center = np.array([row.center_x_kpc, row.center_y_kpc, row.center_z_kpc], dtype=float)
            local = rotate_points_to_local(points - center[None, :], float(row.angle_deg))
            u = np.sqrt(
                (local[:, 0] / float(row.a_radius_kpc)) ** 2
                + (local[:, 1] / float(row.b_radius_kpc)) ** 2
            )
            height_mask = np.abs(local[:, 2]) <= DUST_CYLINDER_HALF_HEIGHT_KPC
            membership[:, idx] = height_mask & (u >= SHELL_U_MIN) & (u <= SHELL_U_MAX)
            u_values[:, idx] = u
            continue

        u, side_mask = ellipsoid_cap_u_and_mask(points, row)
        membership[:, idx] = side_mask & (u >= SHELL_U_MIN) & (u <= SHELL_U_MAX)
        u_values[:, idx] = u
    return membership, u_values


def add_membership_columns(
    clusters_df: pd.DataFrame,
    final_df: pd.DataFrame,
) -> pd.DataFrame:
    points = clusters_df[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(dtype=float)
    cap_matrix = cap_interior_membership_matrix(points, final_df)
    shell_matrix, shell_u = dust_shell_membership_matrix(points, final_df)

    final_ids = final_df["id"].to_numpy(dtype=int)
    out = assign_age_bins(clusters_df)
    out["inside_cap_interior"] = cap_matrix.any(axis=1)
    out["inside_dust_shell"] = shell_matrix.any(axis=1)
    out["containing_cap_interior_sb_ids"] = [";".join(f"SB{value}" for value in final_ids[row]) for row in cap_matrix]
    out["containing_dust_shell_sb_ids"] = [";".join(f"SB{value}" for value in final_ids[row]) for row in shell_matrix]
    out["containing_dust_shell_u_values"] = [
        ";".join(f"{value:.6f}" for value in shell_u[index, row])
        for index, row in enumerate(shell_matrix)
    ]
    out["inside_any_superbubble_region"] = out["inside_cap_interior"] | out["inside_dust_shell"]
    return out


def random_region_counts_by_sb(
    target: pd.DataFrame,
    final_df: pd.DataFrame,
    n_random: int,
    seed_rng: np.random.Generator,
    batch_size: int,
) -> dict[str, np.ndarray]:
    points = target[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(dtype=float)
    n_sb = len(final_df)
    distance = np.linalg.norm(points, axis=1)
    galactic_lat = np.arcsin(np.divide(points[:, 2], distance, out=np.zeros(len(points), dtype=float), where=distance > 0.0))
    rho_xy = distance * np.cos(galactic_lat)
    z = distance * np.sin(galactic_lat)
    counts = {key: np.zeros((int(n_random), n_sb), dtype=int) for key, _label, _column, _color in REGIONS}
    if len(points) == 0 or n_sb == 0:
        return counts

    for start in range(0, int(n_random), int(batch_size)):
        batch_n = min(int(batch_size), int(n_random) - start)
        random_l_rad = seed_rng.uniform(0.0, 2.0 * np.pi, size=(batch_n, len(points)))
        random_points = np.stack(
            [
                rho_xy[None, :] * np.cos(random_l_rad),
                rho_xy[None, :] * np.sin(random_l_rad),
                np.broadcast_to(z[None, :], (batch_n, len(points))),
            ],
            axis=2,
        ).reshape((-1, 3))

        cap_matrix = cap_interior_membership_matrix(random_points, final_df).reshape((batch_n, len(points), n_sb))
        shell_matrix = dust_shell_membership_matrix(random_points, final_df)[0].reshape((batch_n, len(points), n_sb))

        counts["cap_interior"][start : start + batch_n, :] = np.count_nonzero(cap_matrix, axis=1)
        counts["dust_shell"][start : start + batch_n, :] = np.count_nonzero(shell_matrix, axis=1)
    return counts


def random_region_counts_overall(
    points_df: pd.DataFrame,
    final_df: pd.DataFrame,
    n_random: int,
    rng: np.random.Generator,
    batch_size: int,
) -> dict[str, np.ndarray]:
    points = points_df[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(dtype=float)
    distance = np.linalg.norm(points, axis=1)
    galactic_lat = np.arcsin(np.divide(points[:, 2], distance, out=np.zeros(len(points), dtype=float), where=distance > 0.0))
    rho_xy = distance * np.cos(galactic_lat)
    z = distance * np.sin(galactic_lat)

    counts = {key: np.zeros(int(n_random), dtype=int) for key, _label, _column, _color in REGIONS}
    for start in range(0, int(n_random), int(batch_size)):
        batch_n = min(int(batch_size), int(n_random) - start)
        random_l_rad = rng.uniform(0.0, 2.0 * np.pi, size=(batch_n, len(points)))
        random_points = np.stack(
            [
                rho_xy[None, :] * np.cos(random_l_rad),
                rho_xy[None, :] * np.sin(random_l_rad),
                np.broadcast_to(z[None, :], (batch_n, len(points))),
            ],
            axis=2,
        ).reshape((-1, 3))

        cap_matrix = cap_interior_membership_matrix(random_points, final_df).reshape((batch_n, len(points), len(final_df)))
        shell_matrix = dust_shell_membership_matrix(random_points, final_df)[0].reshape((batch_n, len(points), len(final_df)))
        counts["cap_interior"][start : start + batch_n] = np.count_nonzero(cap_matrix.any(axis=2), axis=1)
        counts["dust_shell"][start : start + batch_n] = np.count_nonzero(shell_matrix.any(axis=2), axis=1)
    return counts


def summarize_one(
    region: str,
    region_label: str,
    age_label: str,
    age_interval: str,
    age_window: str,
    cluster_count: int,
    observed: int,
    random_counts: np.ndarray,
    n_random: int,
    seed: int,
) -> dict:
    total = int(cluster_count)
    counts = np.asarray(random_counts, dtype=float)
    mean = float(np.mean(counts))
    std = float(np.std(counts, ddof=1)) if len(counts) > 1 else 0.0
    if mean > 0.0:
        ratio = observed / mean
        relative_difference = 100.0 * (ratio - 1.0)
    elif observed == 0:
        ratio = np.nan
        relative_difference = np.nan
    else:
        ratio = np.inf
        relative_difference = np.inf
    return {
        "region": region,
        "region_description": region_label,
        "age_bin_myr": age_label,
        "age_interval": age_interval,
        "age_window_myr": age_window,
        "cluster_count": total,
        "observed_count": observed,
        "observed_fraction": observed / float(total) if total else np.nan,
        "random_method": "galactic_longitude_scramble_preserve_distance_and_latitude",
        "random_trials": int(n_random),
        "random_seed": int(seed),
        "random_mean_count": mean,
        "random_std_count": std,
        "random_median_count": float(np.median(counts)),
        "random_p025_count": float(np.percentile(counts, 2.5)),
        "random_p975_count": float(np.percentile(counts, 97.5)),
        "observed_over_random_mean": ratio,
        "relative_difference_percent": relative_difference,
        "relative_random_std_percent": 100.0 * std / mean if mean > 0.0 else np.nan,
    }


def build_results(
    eligible: pd.DataFrame,
    final_df: pd.DataFrame,
    n_random: int,
    seed: int,
    batch_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    members = add_membership_columns(eligible, final_df)
    any_window = np.zeros(len(members), dtype=bool)
    for _age_label, _age_interval, lower, upper, include_lower in AGE_BINS:
        any_window |= age_window_mask(members, lower, upper, include_lower).to_numpy(dtype=bool)
    selected = members[any_window].copy().reset_index(drop=True)
    rng = np.random.default_rng(int(seed))
    summary_rows: list[dict] = []
    trial_frames: list[pd.DataFrame] = []
    final_ids = final_df["id"].to_numpy(dtype=int)
    final_shapes = final_df["shape"].astype(str).to_numpy()
    geometry_models = [geometry_model_label(row) for row in final_df.itertuples(index=False)]

    for age_label, age_interval, lower, upper, include_lower in AGE_BINS:
        age_window = age_window_text(lower, upper)
        target = members[age_window_mask(members, lower, upper, include_lower)].copy()
        target_points = target[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(dtype=float)
        observed_matrices = {
            "cap_interior": cap_interior_membership_matrix(target_points, final_df),
            "dust_shell": dust_shell_membership_matrix(target_points, final_df)[0],
        }
        random_counts = random_region_counts_by_sb(target, final_df, int(n_random), rng, int(batch_size))
        for sb_index, sb_id in enumerate(final_ids):
            for region, region_label, _membership_column, _color in REGIONS:
                counts = random_counts[region][:, sb_index]
                row = summarize_one(
                    region,
                    region_label,
                    age_label,
                    age_interval,
                    age_window,
                    len(target),
                    int(np.count_nonzero(observed_matrices[region][:, sb_index])),
                    counts,
                    int(n_random),
                    int(seed),
                )
                row["sb_id"] = int(sb_id)
                row["sb_label"] = f"SB{int(sb_id)}"
                row["model_shape"] = str(final_shapes[sb_index])
                row["geometry_model"] = geometry_models[sb_index]
                summary_rows.append(row)
                trial_frames.append(
                    pd.DataFrame(
                        {
                            "sb_id": int(sb_id),
                            "sb_label": f"SB{int(sb_id)}",
                            "model_shape": str(final_shapes[sb_index]),
                            "geometry_model": geometry_models[sb_index],
                            "region": region,
                            "age_bin_myr": age_label,
                            "age_interval": age_interval,
                            "age_window_myr": age_window,
                            "trial": np.arange(1, int(n_random) + 1, dtype=int),
                            "random_count": counts,
                        }
                    )
                )

    summary = pd.DataFrame(summary_rows)
    summary["age_midpoint_myr"] = summary["age_bin_myr"].map(AGE_MIDPOINTS)
    summary["shell_u_min_inclusive"] = SHELL_U_MIN
    summary["shell_u_max_inclusive"] = SHELL_U_MAX
    summary["dust_cylinder_half_height_kpc"] = DUST_CYLINDER_HALF_HEIGHT_KPC
    summary["full_cylinder_half_height_default_kpc"] = DEFAULT_FULL_CYLINDER_HALF_HEIGHT_KPC
    return summary, pd.concat(trial_frames, ignore_index=True), selected


def build_overall_results(
    members: pd.DataFrame,
    final_df: pd.DataFrame,
    n_random: int,
    seed: int,
    batch_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(int(seed))
    summary_rows: list[dict] = []
    trial_frames: list[pd.DataFrame] = []

    for age_label, age_interval, lower, upper, include_lower in AGE_BINS:
        age_window = age_window_text(lower, upper)
        target = members[age_window_mask(members, lower, upper, include_lower)].copy()
        target_points = target[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(dtype=float)
        observed_matrices = {
            "cap_interior": cap_interior_membership_matrix(target_points, final_df),
            "dust_shell": dust_shell_membership_matrix(target_points, final_df)[0],
        }
        random_counts = random_region_counts_overall(target, final_df, int(n_random), rng, int(batch_size))
        for region, region_label, _membership_column, _color in REGIONS:
            observed = int(np.count_nonzero(observed_matrices[region].any(axis=1)))
            counts = random_counts[region]
            row = summarize_one(
                region,
                region_label,
                age_label,
                age_interval,
                age_window,
                len(target),
                observed,
                counts,
                int(n_random),
                int(seed),
            )
            row["sb_id"] = 0
            row["sb_label"] = "All SBs"
            row["model_shape"] = "mixed"
            row["geometry_model"] = "mixed"
            summary_rows.append(row)
            trial_frames.append(
                pd.DataFrame(
                    {
                        "sb_id": 0,
                        "sb_label": "All SBs",
                        "model_shape": "mixed",
                        "geometry_model": "mixed",
                        "region": region,
                        "age_bin_myr": age_label,
                        "age_interval": age_interval,
                        "age_window_myr": age_window,
                        "trial": np.arange(1, int(n_random) + 1, dtype=int),
                        "random_count": counts,
                    }
                )
            )

    summary = pd.DataFrame(summary_rows)
    summary["age_midpoint_myr"] = summary["age_bin_myr"].map(AGE_MIDPOINTS)
    summary["shell_u_min_inclusive"] = SHELL_U_MIN
    summary["shell_u_max_inclusive"] = SHELL_U_MAX
    summary["dust_cylinder_half_height_kpc"] = DUST_CYLINDER_HALF_HEIGHT_KPC
    summary["full_cylinder_half_height_default_kpc"] = DEFAULT_FULL_CYLINDER_HALF_HEIGHT_KPC
    return summary, pd.concat(trial_frames, ignore_index=True)


def save_publication_figure(fig: plt.Figure, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")


def plot_histograms(summary: pd.DataFrame, trials: pd.DataFrame, output_path: Path) -> None:
    histogram_bins = [item for item in AGE_BINS if item[0] in HISTOGRAM_AGE_LABELS]
    n_age = len(histogram_bins)
    n_cols = 3
    n_age_blocks = int(np.ceil(n_age / n_cols))
    n_rows = n_age_blocks * len(REGIONS)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.1, 9.4), constrained_layout=True)
    axes = np.asarray(axes)
    panel_index = 0
    for region_index, (region, region_label, _membership_column, color) in enumerate(REGIONS):
        group_start_row = region_index * n_age_blocks
        for age_index, (age_label, _age_interval, _lower, _upper, _include_lower) in enumerate(histogram_bins):
            age_row = age_index // n_cols
            col_index = age_index % n_cols
            ax = axes[group_start_row + age_row, col_index]
            one = summary[(summary["region"] == region) & (summary["age_bin_myr"] == age_label)].iloc[0]
            arr = trials.loc[
                (trials["region"] == region) & (trials["age_bin_myr"] == age_label),
                "random_count",
            ].to_numpy(dtype=float)
            bins = np.arange(np.nanmin(arr) - 0.5, np.nanmax(arr) + 1.5, 1.0)
            ax.hist(arr, bins=bins, color="0.82", edgecolor="white", linewidth=0.45)
            ax.axvline(float(one["observed_count"]), color=color, linewidth=1.75)
            ax.axvline(float(one["random_mean_count"]), color="0.12", linestyle="--", linewidth=1.15)
            ax.set_title(
                f"{one['age_window_myr']} Myr",
                fontsize=11.2,
                pad=3.0,
            )
            ax.text(
                0.03,
                0.94,
                f"{chr(ord('A') + panel_index)}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=12.0,
                fontweight="bold",
            )
            ax.text(
                0.97,
                0.94,
                f"{one['relative_difference_percent']:+.1f}%\n"
                f"$1\\sigma$={one['relative_random_std_percent']:.1f}%",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=11.2,
                bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "none", "alpha": 0.84},
            )
            if col_index == 0:
                region_axis_label = "Interiors" if region == "cap_interior" else "Shells"
                ax.set_ylabel(f"{region_axis_label}\nTrials", fontsize=12.0)
            if age_row == n_age_blocks - 1:
                ax.set_xlabel("Random count", fontsize=12.0)
            ax.grid(axis="y", color="0.90", linewidth=0.7)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(length=3.8, width=0.8, labelsize=11.0)
            panel_index += 1

    for empty_index in range(n_age, n_age_blocks * n_cols):
        age_row = empty_index // n_cols
        col_index = empty_index % n_cols
        for region_index in range(len(REGIONS)):
            group_start_row = region_index * n_age_blocks
            axes[group_start_row + age_row, col_index].set_visible(False)

    handles = [
        Line2D([], [], color="0.82", linewidth=6, label="Random trials"),
        Line2D([], [], color="0.12", linestyle="--", linewidth=1.15, label="Random mean"),
    ]
    for region, _label, _membership_column, color in REGIONS:
        handles.append(Line2D([], [], color=color, linewidth=1.75, label=f"Observed: {PLOT_REGION_LABELS[region]}"))
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, handlelength=2.2, fontsize=11.2)
    save_publication_figure(fig, output_path)
    plt.close(fig)


def plot_percent_trend(summary: pd.DataFrame, output_path: Path, title: str | None = None) -> None:
    fig, ax = plt.subplots(figsize=(4.7, 3.0), constrained_layout=True)
    if title:
        ax.set_title(title, fontsize=10.2, pad=4.0)
    ax.axhline(0.0, color="0.18", linewidth=1.0, linestyle="--", zorder=1)
    x_offsets = {"cap_interior": -0.12, "dust_shell": 0.12}
    markers = {"cap_interior": "o", "dust_shell": "s"}
    linestyles = {"cap_interior": "-", "dust_shell": "--"}
    for region, label, _membership_column, color in REGIONS:
        one = summary[summary["region"] == region].sort_values("age_midpoint_myr")
        x_values = one["age_midpoint_myr"].to_numpy(dtype=float) + x_offsets[region]
        ax.errorbar(
            x_values,
            one["relative_difference_percent"],
            yerr=one["relative_random_std_percent"],
            color=color,
            capsize=2.6,
            elinewidth=1.0,
            fmt=markers[region],
            linestyle=linestyles[region],
            linewidth=1.5,
            markersize=4.6,
            markeredgecolor="white",
            markeredgewidth=0.4,
            label=PLOT_REGION_LABELS[region],
            zorder=3,
        )
    ax.set_xticks(AGE_TICK_POSITIONS, AGE_TICK_LABELS)
    ax.set_xlabel("Star-cluster age window (Myr)", fontsize=10.0)
    ax.set_ylabel("Count difference from random (%)", fontsize=10.0)
    ax.set_xlim(2.5, 57.5)
    ax.tick_params(axis="x", rotation=0, labelsize=9.3)
    ax.tick_params(axis="y", labelsize=9.3)
    ax.grid(axis="y", color="0.90", linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3.2)
    ax.legend(
        loc="upper right",
        ncol=1,
        frameon=True,
        framealpha=0.92,
        edgecolor="0.80",
        facecolor="white",
        borderpad=0.45,
        handlelength=1.7,
        fontsize=8.4,
    )
    save_publication_figure(fig, output_path)
    plt.close(fig)


def write_report(
    path: Path,
    cluster_path: Path,
    eligible: pd.DataFrame,
    selected: pd.DataFrame,
    final_df: pd.DataFrame,
    summary: pd.DataFrame,
    n_random: int,
    seed: int,
    batch_size: int,
) -> None:
    lines = [
        "# Script 18 method notes",
        "",
        "## Input and filtering",
        "",
        f"- Cluster input: `{relative_path_text(cluster_path)}`.",
        f"- Eligible clusters after the young-cluster filters: {len(eligible)}.",
        f"- Clusters assigned to at least one tested region: {len(selected)}.",
        f"- Monte Carlo settings: {n_random} longitude-scramble trials, seed={seed}, batch_size={batch_size}.",
        "",
        "## Two region definitions",
        "",
        "- Interior uses the fitted OSB volume clipped by the cap convention used in the main shell analysis.",
        "- Shell uses the same geometry but tests only the boundary layer around each fitted surface.",
        "- Ellipsoid and cylinder variants are both retained so the diagnostic is transparent to the geometry choice.",
        "",
        "## Random test",
        "",
        "- Each random trial preserves cluster distance and latitude while scrambling Galactic longitude.",
        "- The observed count is compared with the random-trial mean for each SB, age bin, region, and geometry.",
        "- The reported relative offset is `(observed - random_mean) / random_mean` in percent.",
        "",
        "## Result summary",
        "",
        "| SB | Model shape | Geometry model | Region | Age interval | Age window (Myr) | Cluster count | Observed count | Random mean | Relative difference | Random scatter |",
        "| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    shape_order = {"ellipsoid": 0, "cylinder": 1}
    order = {region: i for i, (region, _label, _column, _color) in enumerate(REGIONS)}
    age_order = {label: i for i, label in enumerate(AGE_LABELS)}
    summary = summary.assign(
        _shape_order=summary["model_shape"].map(shape_order),
        _region_order=summary["region"].map(order),
        _age_order=summary["age_bin_myr"].map(age_order),
    )
    for row in summary.sort_values(["sb_id", "_region_order", "_age_order"]).itertuples(index=False):
        lines.append(
            f"| {row.sb_label} | {row.model_shape} | {row.geometry_model} | {row.region} | {row.age_interval} | {row.age_window_myr} | {int(row.cluster_count)} | {int(row.observed_count)} | "
            f"{row.random_mean_count:.2f} | {row.relative_difference_percent:+.1f}% | {row.relative_random_std_percent:.1f}% |"
        )
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    """Build the command-line interface for Script 18."""
    parser = argparse.ArgumentParser(description='Run Script 18: YSO OSBs shell random test.')
    parser.add_argument("--hunt-raw-csv", default=str(HUNT_RAW_CSV), help='Input file used by this stage.')
    parser.add_argument("--hunt-output-csv", default=str(HUNT_CLUSTER_CSV), help='Output path for this product.')
    parser.add_argument("--hunt-filter-summary-csv", default=str(HUNT_FILTER_SUMMARY_CSV), help='Input file used by this stage.')
    parser.add_argument("--cluster-csv", default=str(CLUSTER_CSV), help='Input file used by this stage.')
    parser.add_argument("--skip-build-hunt-table", action="store_true", help='Input file used by this stage.')
    parser.add_argument("--final-sb-csv", default=str(FINAL_SB_CSV), help='Input file used by this stage.')
    parser.add_argument("--n-random", type=int, default=DEFAULT_N_RANDOM, help='Number of random or Monte Carlo realizations.')
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help='Random seed for reproducible controls.')
    parser.add_argument("--random-batch-size", type=int, default=DEFAULT_RANDOM_BATCH_SIZE, help='Command-line option for the documented workflow.')
    return parser.parse_args()


def main() -> None:
    """Run Script 18 from validated inputs to the documented outputs."""
    args = parse_args()
    if int(args.n_random) <= 0:
        raise ValueError('Invalid input or missing required data.')
    if int(args.random_batch_size) <= 0:
        raise ValueError('Invalid input or missing required data.')

    if not args.skip_build_hunt_table:
        hunt_raw_path = ensure_existing_file(args.hunt_raw_csv, "Hunt2024 raw cluster master table")
        hunt_output_path = Path(args.hunt_output_csv)
        hunt_summary_path = Path(args.hunt_filter_summary_csv)
        hunt_table = build_hunt_cluster_table(hunt_raw_path, hunt_output_path, hunt_summary_path)
        print(f"Generated Hunt cluster table: {hunt_output_path} ({len(hunt_table)} rows)")
        print(f"Generated Hunt filter summary: {hunt_summary_path}")
        if args.cluster_csv == str(CLUSTER_CSV):
            args.cluster_csv = str(hunt_output_path)

    cluster_path = ensure_existing_file(args.cluster_csv, "young-cluster input table")
    final_path = ensure_existing_file(args.final_sb_csv, "final OSB fitted-parameter table")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)

    eligible = load_clusters(cluster_path)
    final_df = load_final_superbubbles(final_path)
    summary, trials, selected = build_results(
        eligible,
        final_df,
        int(args.n_random),
        int(args.seed),
        int(args.random_batch_size),
    )
    overall_summary, overall_trials = build_overall_results(
        selected,
        final_df,
        int(args.n_random),
        int(args.seed),
        int(args.random_batch_size),
    )

    for output_df, output_path in [
        (summary, OUT_SUMMARY_CSV),
        (trials, OUT_RANDOM_CSV),
        (selected, OUT_MEMBERS_CSV),
        (overall_summary, OUT_PLOT_DATA_CSV),
    ]:
        try:
            output_df.to_csv(output_path, index=False, encoding="utf-8-sig")
        except PermissionError as exc:
            print(f"Warning: cannot write locked output file {output_path}: {exc}")
    plot_histograms(overall_summary, overall_trials, OUT_HIST_FIG)
    plot_percent_trend(overall_summary, OUT_LINE_FIG)
    try:
        write_report(
            OUT_REPORT_MD,
            cluster_path,
            eligible,
            selected,
            final_df,
            summary,
            int(args.n_random),
            int(args.seed),
            int(args.random_batch_size),
        )
    except PermissionError as exc:
        print(f"Warning: cannot write locked report file {OUT_REPORT_MD}: {exc}")

    selected_unique = selected["row_id"].nunique() if "row_id" in selected.columns else len(selected)
    print(f"Eligible young clusters: {len(eligible)}")
    print(f"Selected cluster-region memberships: {selected_unique}")
    print(f"number of final OSBs: {len(final_df)}")
    print(
        summary[
            [
                "model_shape",
                "geometry_model",
                "sb_label",
                "region",
                "age_interval",
                "cluster_count",
                "observed_count",
                "random_mean_count",
                "relative_difference_percent",
                "relative_random_std_percent",
            ]
        ].to_string(index=False)
    )
    print(f"Saved: {OUT_SUMMARY_CSV}")
    print(f"Saved: {OUT_RANDOM_CSV}")
    print(f"Saved: {OUT_MEMBERS_CSV}")
    print(f"Saved: {OUT_PLOT_DATA_CSV}")
    print(f"Saved: {OUT_REPORT_MD}")
    print(f"Saved: {OUT_HIST_FIG}")
    print(f"Saved: {OUT_LINE_FIG}")

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
