"""Script 5: OSBs parameter statistics

Purpose
-------
Measure cavity/shell/ridge statistics and evidence-grade parameters.

Method overview
---------------
1. Load the final geometry table and extract the raw dust voxels around each fitted
   object.
2. Measure cavity density, shell-density contrast and radial shell-ridge statistics with
   compute_all_dust_metrics().
3. Join the geometry and dust statistics in prepare_table(), then write the per-target
   measurements and parameter-distribution figures.

Main inputs
-----------
- ../results/superbubble_final_fit_parameters.csv
- ../results/intermediate_output/3d_dust_map_products/raw_3d_dust_cube.parquet

Main outputs
------------
- ../results/figures/5_parameter_histogram_matrix.png; Extended Data parameter histogram.

Figure/table role
-----------------
../results/figures/5_parameter_histogram_matrix.png; Extended Data parameter histogram.

Runtime and data notes
----------------------
Reference runtime: 9 min 15.09 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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

import os
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.optimize import curve_fit


# ----------------------------------------------------------------------------
# Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
# ----------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"
GEOM_CANDIDATES = [
    FINAL_DIR / "superbubble_final_fit_parameters.csv",
]
DUST_PARQUET = Path("..") / "results" / "intermediate_output" / "3d_dust_map_products" / "raw_3d_dust_cube.parquet"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "5_dust_parameter_selection"
OUT_DUST_CSV = OUT_DIR / "5_dust_parameter_statistics.csv"
OUT_HIST_CSV = OUT_DIR / "5_eight_parameter_histogram_data.csv"
OUT_FIG = FINAL_FIG_DIR / "5_parameter_histogram_matrix.png"
OUT_FIG_CURV = OUT_DIR / "5_ridge_curvature_smoothness_histogram.png"

# ----------------------------------------------------------------------------
# Geometry and shell-mask conventions used by this analysis stage.
# ----------------------------------------------------------------------------
SHELL_U_MIN = 0.9
SHELL_U_MAX = 1.2
INNER_U_MAX = 0.5
CYLINDER_HALF_HEIGHT_KPC = 0.025

INNER_STD_CLIP_ITERS = 3
INNER_STD_CLIP_SIGMA = 3.0

RIDGE_N_AZIMUTH = 36
RIDGE_N_COSTHETA_FULL = 12
RIDGE_N_COSTHETA_HALF = 6
RIDGE_N_RADIAL = 16
RIDGE_MIN_SECTOR_VOXELS = 10
RIDGE_MIN_RADIAL_VOXELS = 2
RIDGE_PEAK_INNER_RATIO_MIN = 2.0
RIDGE_PEAK_MIN_VOXELS = 30
RIDGE_GMM_MAX_ITER = 200
RIDGE_GMM_TOL = 1e-6
RIDGE_CURV_MIN_CELLS = 5

SEL_FIT_STD_MAX_PC = 100.0
SEL_MEAN_ABS_OFFSET_MAX = 0.10
SEL_CURVATURE_MAX = 0.10
SEL_RIDGE_RATIO_MIN = 5.0
SEL_INNER_MEAN_MAX = 0.2
SEL_INNER_STD_MAX = 0.15

HIGHLIGHT_COLOR = "#F58518"
BLUE = "#4C78A8"

# Histogram figure sizing and typography. The multi-panel validation figure is
# rendered at A4 landscape size so printed text is not down-scaled by page-fit.
A4_LANDSCAPE_FIGSIZE = (11.69, 8.27)
PANEL_TITLE_FONTSIZE = 13.5
AXIS_LABEL_FONTSIZE = 12.5
TICK_LABEL_FONTSIZE = 11.5
PANEL_LETTER_FONTSIZE = 14.0
PANEL_NOTE_FONTSIZE = 11.5
LEGEND_FONTSIZE = 12.5
LEGEND_NOTE_FONTSIZE = 12.5
PANEL_B_FIT_NOTE_FONTSIZE = 14.0
SINGLE_FIG_TITLE_FONTSIZE = 14.0
A_GRADE_EXCLUDE_IDS = {27}

# Figure-panel rendering section.
PANEL_B_BIN_MIN = 170.0
PANEL_B_BIN_MAX = 1400.0
PANEL_B_N_BINS = 8
PANEL_B_FIT_MIN = 350.0
PANEL_B_BOOTSTRAP = 1000
PANEL_B_SEED = 20260427


# ----------------------------------------------------------------------------
# Geometry and shell-mask conventions used by this analysis stage.
# ----------------------------------------------------------------------------
def rotate_points_to_local(points: np.ndarray, angle_deg: float) -> np.ndarray:
    """Rotate coordinates into the local frame used by the shell or projection calculation."""
    theta = np.deg2rad(-float(angle_deg))
    c, s = np.cos(theta), np.sin(theta)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    return np.asarray(points, dtype=float) @ rot.T


def read_local_xyz_cube(parquet_path, x_range, y_range, z_range) -> pd.DataFrame:
    """Load and validate the input table or array needed by this stage."""
    parquet = pq.ParquetFile(str(parquet_path))
    names = parquet.schema.names
    if {"x", "y", "z", "dust"}.issubset(names):
        cols = ["x", "y", "z", "dust"]
    elif {"X", "Y", "Z", "dust"}.issubset(names):
        cols = ["X", "Y", "Z", "dust"]
    elif {"x", "y", "z", "dust_raw"}.issubset(names):
        cols = ["x", "y", "z", "dust_raw"]
    elif {"X", "Y", "Z", "dust_raw"}.issubset(names):
        cols = ["X", "Y", "Z", "dust_raw"]
    else:
        raise ValueError(f"Could not identify XYZ parquet columns: {names}")
    chunks = []
    for rg in range(parquet.num_row_groups):
        df = parquet.read_row_group(rg, columns=cols).to_pandas()
        df = df.rename(columns={"X": "x", "Y": "y", "Z": "z", "dust_raw": "dust"})
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


# ----------------------------------------------------------------------------
# Geometry and shell-mask conventions used by this analysis stage.
# ----------------------------------------------------------------------------
def geometry_kind(row: pd.Series) -> str:
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    if str(row["shape"]) == "cylinder":
        return "cylinder"
    return "ellipsoidal_cap"


def cap_geometry(row: pd.Series):
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    cx = float(row["center_x_kpc"])
    cy = float(row["center_y_kpc"])
    cz = float(row["center_z_kpc"])
    c_full = float(row["c_radius_kpc"])
    xy_z = float(row["xy_plane_z_kpc"])
    A = float(row["xy_plane_a_kpc"])
    B = float(row["xy_plane_b_kpc"])
    mark = int(row["mark"])
    z_apex = cz - c_full if mark == 1 else cz + c_full
    c_cap = abs(xy_z - z_apex)
    cap_sign = -1.0 if mark == 1 else 1.0
    return np.array([cx, cy, xy_z], dtype=float), max(A, 1e-9), max(B, 1e-9), max(c_cap, 1e-9), cap_sign


def cube_ranges(row: pd.Series):
    cx, cy, cz = float(row["center_x_kpc"]), float(row["center_y_kpc"]), float(row["center_z_kpc"])
    a, b, c = float(row["a_radius_kpc"]), float(row["b_radius_kpc"]), float(row["c_radius_kpc"])
    kind = geometry_kind(row)
    if kind == "cylinder":
        xy_margin = max(3.0 * max(a, b), 0.6)
        z_margin = max(3.0 * max(CYLINDER_HALF_HEIGHT_KPC, 1e-6), 0.35)
        return (cx - xy_margin, cx + xy_margin), (cy - xy_margin, cy + xy_margin), (cz - z_margin, cz + z_margin)
    _, A, B, c_cap, _ = cap_geometry(row)
    xy_margin = max(3.0 * max(A, B, a, b), 0.6)
    z_margin = max(3.0 * max(c, c_cap), 0.35)
    return (cx - xy_margin, cx + xy_margin), (cy - xy_margin, cy + xy_margin), (cz - z_margin, cz + z_margin)


def build_local(row: pd.Series, points: np.ndarray):
    """Compute an intermediate statistic or table used by this documented workflow."""
    angle_deg = float(row["angle_deg"])
    kind = geometry_kind(row)

    if kind == "cylinder":
        cx, cy, cz = float(row["center_x_kpc"]), float(row["center_y_kpc"]), float(row["center_z_kpc"])
        a, b = float(row["a_radius_kpc"]), float(row["b_radius_kpc"])
        shifted = np.asarray(points, dtype=float) - np.array([cx, cy, cz])[None, :]
        local = rotate_points_to_local(shifted, angle_deg)
        X, Y, Z = local[:, 0], local[:, 1], local[:, 2]
        q = (X / max(a, 1e-9)) ** 2 + (Y / max(b, 1e-9)) ** 2
        u = np.sqrt(q)
        in_box = np.abs(Z) <= CYLINDER_HALF_HEIGHT_KPC
        shell_mask = in_box & (q >= SHELL_U_MIN ** 2) & (q <= SHELL_U_MAX ** 2)
        inner_mask = in_box & (q <= INNER_U_MAX ** 2)
        return local, kind, u, inner_mask, shell_mask

    base_center, A, B, c_cap, cap_sign = cap_geometry(row)
    shifted = np.asarray(points, dtype=float) - base_center[None, :]
    local = rotate_points_to_local(shifted, angle_deg)
    X, Y, Z = local[:, 0], local[:, 1], local[:, 2]
    q = (X / max(A, 1e-9)) ** 2 + (Y / max(B, 1e-9)) ** 2 + (Z / max(c_cap, 1e-9)) ** 2
    u = np.sqrt(q)
    cap_side = (cap_sign * Z) >= 0.0
    shell_mask = cap_side & (q >= SHELL_U_MIN ** 2) & (q <= SHELL_U_MAX ** 2)
    inner_mask = cap_side & (q <= INNER_U_MAX ** 2)
    return local, kind, u, inner_mask, shell_mask


def cap_mask_local(row: pd.Series, local: np.ndarray, kind: str) -> np.ndarray:
    """Evaluate the membership or shell-selection mask for this analysis stage."""
    if kind == "cylinder":
        return np.ones(local.shape[0], dtype=bool)
    _, _, _, _, cap_sign = cap_geometry(row)
    return (cap_sign * local[:, 2]) >= 0.0


def build_direction_bins(local: np.ndarray, kind: str, cap_mask: np.ndarray):
    """Compute an intermediate statistic or table used by this documented workflow."""
    X, Y, Z = local[:, 0], local[:, 1], local[:, 2]
    phi = (np.arctan2(Y, X) + 2.0 * np.pi) % (2.0 * np.pi)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, RIDGE_N_AZIMUTH + 1)
    phi_idx = np.clip(np.digitize(phi, phi_edges) - 1, 0, RIDGE_N_AZIMUTH - 1)

    if kind == "cylinder":
        return phi_idx.astype(int), RIDGE_N_AZIMUTH

    r = np.sqrt(X * X + Y * Y + Z * Z)
    safe = r > 1e-12
    abs_cos = np.where(safe, np.abs(Z) / np.maximum(r, 1e-12), 0.0)
    n_ct = RIDGE_N_COSTHETA_HALF
    ct_edges = np.linspace(0.0, 1.0, n_ct + 1)
    ct_idx = np.clip(np.digitize(abs_cos, ct_edges) - 1, 0, n_ct - 1)
    sector_idx = np.where(safe & cap_mask, phi_idx * n_ct + ct_idx, -1)
    return sector_idx.astype(int), RIDGE_N_AZIMUTH * n_ct


# ----------------------------------------------------------------------------
# Numerical helper section.
# ----------------------------------------------------------------------------
def fit_two_gaussian_high_peak(values: np.ndarray):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < RIDGE_PEAK_MIN_VOXELS or float(np.std(values)) <= 0:
        return float("nan"), float("nan"), float("nan"), "insufficient_shell_voxels"
    q25, q75 = np.percentile(values, [25, 75])
    means = np.array([q25, q75], dtype=float)
    shared_var = max(float(np.var(values)), 1e-8)
    variances = np.array([shared_var, shared_var], dtype=float)
    weights = np.array([0.5, 0.5], dtype=float)
    eps = 1e-12
    prev = -np.inf
    for _ in range(RIDGE_GMM_MAX_ITER):
        dens = []
        for k in range(2):
            var = max(float(variances[k]), 1e-8)
            coef = 1.0 / np.sqrt(2.0 * np.pi * var)
            dens.append(weights[k] * coef * np.exp(-0.5 * (values - means[k]) ** 2 / var))
        probs = np.vstack(dens).T
        total = np.sum(probs, axis=1) + eps
        resp = probs / total[:, None]
        nk = np.sum(resp, axis=0) + eps
        weights = nk / values.size
        means = np.sum(resp * values[:, None], axis=0) / nk
        variances = np.maximum(np.sum(resp * (values[:, None] - means[None, :]) ** 2, axis=0) / nk, 1e-8)
        loglike = float(np.sum(np.log(total)))
        if abs(loglike - prev) < RIDGE_GMM_TOL:
            break
        prev = loglike
    hi = int(np.argmax(means))
    return float(means[hi]), float(np.sqrt(variances[hi])), float(weights[hi]), "linear_2component_em_high_mean"


def iterative_gaussian_std(values: np.ndarray, n_iter: int = INNER_STD_CLIP_ITERS, n_sigma: float = INNER_STD_CLIP_SIGMA):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    n_total = int(v.size)
    if n_total < 2:
        return float("nan"), float(np.mean(v)) if n_total else float("nan"), n_total, n_total
    keep = np.ones(n_total, dtype=bool)
    for _ in range(int(n_iter)):
        cur = v[keep]
        if cur.size < 2:
            break
        m = float(np.mean(cur))
        s = float(np.std(cur, ddof=1))
        if not np.isfinite(s) or s <= 0:
            break
        new_keep = np.abs(v - m) <= n_sigma * s
        if new_keep.sum() < 2 or np.array_equal(new_keep, keep):
            keep = new_keep if new_keep.sum() >= 2 else keep
            break
        keep = new_keep
    cur = v[keep]
    g_std = float(np.std(cur, ddof=1)) if cur.size > 1 else float("nan")
    g_mean = float(np.mean(cur)) if cur.size else float("nan")
    return g_std, g_mean, int(cur.size), n_total


def ridge_uPeak_curvature_rms(u_peak_grid: np.ndarray):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    grid = np.asarray(u_peak_grid, dtype=float)
    n_phi, n_ct = grid.shape
    residuals: list[float] = []
    for i in range(n_phi):
        ip, im = (i + 1) % n_phi, (i - 1) % n_phi
        for j in range(n_ct):
            c = grid[i, j]
            if not np.isfinite(c):
                continue
            if n_ct == 1:
                neigh = [grid[ip, j], grid[im, j]]
            else:
                if j - 1 < 0 or j + 1 >= n_ct:
                    continue
                neigh = [grid[ip, j], grid[im, j], grid[i, j + 1], grid[i, j - 1]]
            if not all(np.isfinite(v) for v in neigh):
                continue
            residuals.append(float(c - np.mean(neigh)))
    if len(residuals) < RIDGE_CURV_MIN_CELLS:
        return float("nan"), int(len(residuals))
    arr = np.asarray(residuals, dtype=float)
    return float(np.sqrt(np.mean(arr ** 2))), int(arr.size)


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
def compute_metrics(row, dust, local, u, inner_mask, shell_mask, kind):
    dust = np.asarray(dust, dtype=float)
    finite_inner = inner_mask & np.isfinite(dust)
    finite_shell = shell_mask & np.isfinite(dust)
    inner_values = dust[finite_inner]
    shell_values = dust[finite_shell]

    inner_mean = float(np.mean(inner_values)) if inner_values.size else float("nan")
    shell_mean = float(np.mean(shell_values)) if shell_values.size else float("nan")
    inner_std_raw = float(np.std(inner_values, ddof=1)) if inner_values.size > 1 else float("nan")
    inner_std, _, inner_std_kept, _ = iterative_gaussian_std(inner_values)
    shell_inner_ratio = (
        float(shell_mean / inner_mean)
        if np.isfinite(inner_mean) and inner_mean > 0 and np.isfinite(shell_mean) else float("nan")
    )

    ridge_peak_density = ridge_peak_sigma = ridge_peak_weight = ridge_peak_inner_ratio = float("nan")
    ridge_peak_method = "insufficient_shell_voxels"
    shell_pos = shell_values[np.isfinite(shell_values) & (shell_values > 0)]
    if shell_pos.size >= RIDGE_PEAK_MIN_VOXELS and float(np.nanstd(shell_pos)) > 0:
        ridge_peak_density, ridge_peak_sigma, ridge_peak_weight, ridge_peak_method = fit_two_gaussian_high_peak(shell_pos)
        ridge_peak_inner_ratio = (
            float(ridge_peak_density / inner_mean)
            if np.isfinite(inner_mean) and inner_mean > 0 and np.isfinite(ridge_peak_density) else float("nan")
        )

    cap_mask = cap_mask_local(row, local, kind)
    sector_idx, n_sectors_total = build_direction_bins(local, kind, cap_mask)
    u_edges = np.linspace(SHELL_U_MIN, SHELL_U_MAX, RIDGE_N_RADIAL + 1)
    u_centers = 0.5 * (u_edges[:-1] + u_edges[1:])
    u_idx = np.clip(np.digitize(u, u_edges) - 1, 0, RIDGE_N_RADIAL - 1)

    u_peak_by_sector = np.full(n_sectors_total, np.nan, dtype=float)
    for s in range(n_sectors_total):
        sm = finite_shell & (sector_idx == s)
        if int(sm.sum()) < RIDGE_MIN_SECTOR_VOXELS:
            continue
        radial_means = np.full(RIDGE_N_RADIAL, np.nan, dtype=float)
        su, sd = u_idx[sm], dust[sm]
        for j in range(RIDGE_N_RADIAL):
            bv = sd[su == j]
            if bv.size >= RIDGE_MIN_RADIAL_VOXELS:
                radial_means[j] = float(np.nanmean(bv))
        if not np.any(np.isfinite(radial_means)):
            continue
        pi = int(np.nanargmax(radial_means))
        peak_dust = float(radial_means[pi])
        peak_ratio = peak_dust / inner_mean if np.isfinite(inner_mean) and inner_mean > 0 else float("nan")
        if np.isfinite(peak_ratio) and peak_ratio >= RIDGE_PEAK_INNER_RATIO_MIN:
            u_peak_by_sector[s] = float(u_centers[pi])

    valid = u_peak_by_sector[np.isfinite(u_peak_by_sector)]
    mean_signed = float(np.mean(valid - 1.0)) if valid.size else float("nan")
    mean_abs = float(np.mean(np.abs(valid - 1.0))) if valid.size else float("nan")
    n_ct = max(1, n_sectors_total // RIDGE_N_AZIMUTH)
    curv_rms, curv_cells = ridge_uPeak_curvature_rms(u_peak_by_sector.reshape(RIDGE_N_AZIMUTH, n_ct))

    return {
        "dust_inner_voxel_count": int(inner_values.size),
        "dust_shell_voxel_count": int(shell_values.size),
        "dust_inner_mean_mag_kpc": inner_mean,
        "dust_shell_mean_mag_kpc": shell_mean,
        "dust_shell_inner_ratio": shell_inner_ratio,
        "dust_shell_ridge_peak_density_mag_kpc": ridge_peak_density,
        "dust_shell_ridge_peak_sigma_mag_kpc": ridge_peak_sigma,
        "dust_shell_ridge_peak_weight": ridge_peak_weight,
        "dust_shell_ridge_peak_inner_ratio": ridge_peak_inner_ratio,
        "dust_shell_ridge_peak_method": ridge_peak_method,
        "dust_inner_std_mag_kpc": inner_std,
        "dust_inner_std_raw_mag_kpc": inner_std_raw,
        "dust_inner_std_kept_count": int(inner_std_kept),
        "ridge_direction_bin_count": int(n_sectors_total),
        "valid_ridge_sector_count": int(valid.size),
        "shell_ridge_coverage_fraction": float(valid.size / n_sectors_total) if n_sectors_total else float("nan"),
        "shell_radius_mean_signed_offset_over_R_eq": mean_signed,
        "shell_radius_mean_abs_offset_over_R_eq": mean_abs,
        "shell_ridge_uPeak_curvature_rms": curv_rms,
        "shell_ridge_uPeak_curvature_cells": curv_cells,
    }


def compute_all_dust_metrics(geom_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in geom_df.sort_values("id").iterrows():
        kind = geometry_kind(row)
        xr, yr, zr = cube_ranges(row)
        cube = read_local_xyz_cube(DUST_PARQUET, xr, yr, zr)
        pts = cube[["x", "y", "z"]].to_numpy(dtype=float)
        dust = cube["dust"].to_numpy(dtype=float)
        local, _, u, inner_mask, shell_mask = build_local(row, pts)
        metrics = compute_metrics(row, dust, local, u, inner_mask, shell_mask, kind)
        rows.append({
            "id": int(row["id"]), "shape": str(row["shape"]), "mark": int(row["mark"]),
            "final_parameter_estimator": str(row["final_parameter_estimator"]),
            "geometry_kind": kind,
            "center_x_kpc": float(row["center_x_kpc"]), "center_y_kpc": float(row["center_y_kpc"]),
            "center_z_kpc": float(row["center_z_kpc"]),
            "a_radius_kpc": float(row["a_radius_kpc"]), "b_radius_kpc": float(row["b_radius_kpc"]),
            "c_radius_kpc": float(row["c_radius_kpc"]), "angle_deg": float(row["angle_deg"]),
            **metrics,
        })
        print(
            f"SB{int(row['id']):>2} [{kind:>14}]  inner_mean={metrics['dust_inner_mean_mag_kpc']:.3f}  "
            f"inner_std={metrics['dust_inner_std_mag_kpc']:.3f}  ridge/inner={metrics['dust_shell_ridge_peak_inner_ratio']:.3f}  "
            f"mean|uP-1|={metrics['shell_radius_mean_abs_offset_over_R_eq']:.3f}  "
            f"curv={metrics['shell_ridge_uPeak_curvature_rms']:.4f}  "
            f"valid/total={metrics['valid_ridge_sector_count']:>3}/{metrics['ridge_direction_bin_count']}"
        )
    return pd.DataFrame(rows).sort_values("id").reset_index(drop=True)


# ----------------------------------------------------------------------------
# Geometry and shell-mask conventions used by this analysis stage.
# ----------------------------------------------------------------------------
def prepare_table(geom_df: pd.DataFrame, dust_df: pd.DataFrame) -> pd.DataFrame:
    geom_df = geom_df.copy()
    geom_df["id"] = geom_df["id"].astype(int)
    dust_df = dust_df.copy()
    dust_df["id"] = dust_df["id"].astype(int)
    df = geom_df.merge(
        dust_df[[
            "id", "geometry_kind",
            "dust_shell_ridge_peak_inner_ratio", "dust_inner_mean_mag_kpc", "dust_inner_std_mag_kpc",
            "shell_radius_mean_signed_offset_over_R_eq", "shell_radius_mean_abs_offset_over_R_eq",
            "shell_ridge_uPeak_curvature_rms",
        ]],
        on="id", how="inner",
    ).sort_values("id").reset_index(drop=True)

    df["xy_plane_R_eq_pc"] = 1000.0 * np.sqrt(df["xy_plane_a_kpc"] * df["xy_plane_b_kpc"])
    df["ellipsoid_R_eq_pc"] = 1000.0 * np.sqrt(df["a_radius_kpc"] * df["b_radius_kpc"])
    df["fit_std_over_ellipsoid_Req"] = df["fit_std_pc"] / df["ellipsoid_R_eq_pc"]

    df["sel_c_fit_std_lt_100"] = df["fit_std_pc"] < SEL_FIT_STD_MAX_PC
    df["sel_d_mean_abs_offset_lt_0p10"] = df["shell_radius_mean_abs_offset_over_R_eq"] < SEL_MEAN_ABS_OFFSET_MAX
    df["sel_e_uPeak_curvature_lt_0p10"] = df["shell_ridge_uPeak_curvature_rms"] < SEL_CURVATURE_MAX
    df["sel_f_ridge_peak_inner_gt_5"] = df["dust_shell_ridge_peak_inner_ratio"] > SEL_RIDGE_RATIO_MIN
    df["sel_g_inner_mean_lt_0p2"] = df["dust_inner_mean_mag_kpc"] < SEL_INNER_MEAN_MAX
    df["sel_h_inner_std_lt_0p15"] = df["dust_inner_std_mag_kpc"] < SEL_INNER_STD_MAX
    sel_cols = [
        "sel_c_fit_std_lt_100", "sel_d_mean_abs_offset_lt_0p10", "sel_e_uPeak_curvature_lt_0p10",
        "sel_f_ridge_peak_inner_gt_5", "sel_g_inner_mean_lt_0p2", "sel_h_inner_std_lt_0p15",
    ]
    df["pass_selection"] = df[sel_cols].fillna(False).all(axis=1)
    return df


# ----------------------------------------------------------------------------
# Figure styling and layout configuration for reproducible rendering.
# ----------------------------------------------------------------------------
def power_law(x, a, p):
    return a * np.power(x, p)


def plot_panel_b_powerlaw(ax, df, highlight_mask, *, letter="(b)"):
    fit_text = None
    R = pd.to_numeric(df["xy_plane_R_eq_pc"], errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(R)
    hi = np.asarray(highlight_mask, dtype=bool) & finite
    lo = finite & ~hi
    bins = np.logspace(np.log10(PANEL_B_BIN_MIN), np.log10(PANEL_B_BIN_MAX), PANEL_B_N_BINS + 1)
    centers = np.sqrt(bins[:-1] * bins[1:])
    ax.hist(R[finite], bins=bins, color=BLUE, edgecolor="black", linewidth=0.7, alpha=0.72)
    ax.hist(R[hi], bins=bins, color=HIGHLIGHT_COLOR, edgecolor="black", linewidth=0.7, alpha=0.82)
    counts, _ = np.histogram(R[finite], bins=bins)
    fit_mask = (centers > PANEL_B_FIT_MIN) & (counts > 0)
    p_best = None
    if int(fit_mask.sum()) >= 2:
        try:
            p_best, _ = curve_fit(power_law, centers[fit_mask], counts[fit_mask],
                                  p0=(1e6, -2), bounds=((0, -5), (np.inf, 2)), maxfev=20000)
        except RuntimeError:
            p_best = None
    if p_best is not None:
        rng = np.random.default_rng(PANEL_B_SEED)
        Rf = R[finite]
        boot = []
        for _ in range(PANEL_B_BOOTSTRAP):
            res = rng.choice(Rf, size=Rf.size, replace=True)
            c_res, _ = np.histogram(res, bins=bins)
            m = (centers > PANEL_B_FIT_MIN) & (c_res > 0)
            if int(m.sum()) < 2:
                continue
            try:
                po, _ = curve_fit(power_law, centers[m], c_res[m],
                                  p0=(1e6, -2), bounds=((0, -5), (np.inf, 2)), maxfev=5000)
                boot.append(po)
            except RuntimeError:
                continue
        boot = np.array(boot)
        x_model = np.logspace(np.log10(PANEL_B_FIT_MIN), np.log10(PANEL_B_BIN_MAX), 300)
        if boot.size:
            std_p = float(np.std(boot[:, 1]))
            y_lines = np.array([power_law(x_model, *po) for po in boot])
            ax.fill_between(x_model, np.percentile(y_lines, 2.5, axis=0), np.percentile(y_lines, 97.5, axis=0),
                            color="#D62728", alpha=0.18, zorder=2, label="95% CI (bootstrap)")
        else:
            std_p = float("nan")
        fit_text = f"$p={p_best[1]:.1f}\\pm{std_p:.1f}$ (1$\\sigma$)"
        ax.plot(x_model, power_law(x_model, *p_best), color="#D62728", linestyle="--", linewidth=2.2,
                zorder=3, label="Power-law fit")
        ax.axvline(PANEL_B_FIT_MIN, color="gray", linestyle=":", linewidth=1.1, alpha=0.7)
        ax.legend(fontsize=12.5, loc="upper right", frameon=False)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(PANEL_B_BIN_MIN, PANEL_B_BIN_MAX)
    ax.set_ylim(0.8, max(float(counts.max()) * 1.6, 12))
    ax.set_title("XY section equivalent radius", fontsize=PANEL_TITLE_FONTSIZE, pad=8)
    ax.set_xlabel(r"$R_{\rm eq,XY}$ [pc]", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("Number", fontsize=AXIS_LABEL_FONTSIZE)
    xticks = [t for t in (200, 300, 500, 700, 1000) if PANEL_B_BIN_MIN <= t <= PANEL_B_BIN_MAX]
    ax.xaxis.set_major_locator(FixedLocator(xticks))
    ax.xaxis.set_major_formatter(FixedFormatter([str(t) for t in xticks]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    yt = [1, 2, 3, 5, 8, 10]
    ax.yaxis.set_major_locator(FixedLocator(yt))
    ax.yaxis.set_major_formatter(FixedFormatter([str(v) for v in yt]))
    ax.yaxis.set_minor_formatter(NullFormatter())
    if letter:
        ax.text(0.035, 0.93, letter, ha="left", va="top", transform=ax.transAxes, fontsize=PANEL_LETTER_FONTSIZE, fontweight="bold")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE, width=1.0, length=4.5)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
    return fit_text


def nice_bins(values, bins=12, *, zero_floor=False, symmetric=False):
    values = values[np.isfinite(values)]
    if values.size == 0:
        return bins
    if symmetric:
        vmax = max(float(np.nanmax(np.abs(values))), 0.05)
        return np.linspace(-1.08 * vmax, 1.08 * vmax, bins + 1)
    vmin = 0.0 if zero_floor else float(np.nanmin(values))
    vmax = float(np.nanmax(values))
    if np.isclose(vmin, vmax):
        pad = max(abs(vmax) * 0.1, 1.0)
        return np.linspace(vmin - pad, vmax + pad, bins + 1)
    pad = 0.06 * (vmax - vmin)
    vmin = 0.0 if zero_floor else vmin - pad
    vmax += pad
    return np.linspace(vmin, vmax, bins + 1)


def plot_hist(ax, df, column, title, xlabel, *, threshold=None, bins=12, zero_floor=False,
              log_x=False, xlim=None, xticks=None, xticklabels=None, letter=None, highlight_mask=None):
    col = pd.to_numeric(df[column], errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(col)
    values = col[finite]
    if log_x:
        pos = values[values > 0]
        hist_bins = (np.logspace(np.log10(pos.min() * 0.999), np.log10(pos.max() * 1.001), int(bins) + 1)
                     if pos.size else bins)
    else:
        hist_bins = nice_bins(values, bins=bins, zero_floor=zero_floor)
    if highlight_mask is not None:
        hi = np.asarray(highlight_mask, dtype=bool) & finite
        ax.hist(values, bins=hist_bins, color=BLUE, edgecolor="white", linewidth=0.8, alpha=0.72)
        ax.hist(col[hi], bins=hist_bins, color=HIGHLIGHT_COLOR, edgecolor="white", linewidth=0.8, alpha=0.82)
    else:
        ax.hist(values, bins=hist_bins, color=BLUE, edgecolor="white", linewidth=0.8, alpha=0.92)
    if log_x:
        ax.set_xscale("log")
    if values.size:
        med = float(np.nanmedian(values))
        ax.axvline(med, color="#2CA02C", linestyle="-", linewidth=2.3, zorder=5, label=f"median={med:.3g}")
    if threshold is not None:
        ax.axvline(float(threshold), color="#D62728", linestyle="--", linewidth=2.3, zorder=6)
    ax.set_title(title, fontsize=PANEL_TITLE_FONTSIZE, pad=8)
    ax.set_xlabel(xlabel, fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("Number", fontsize=AXIS_LABEL_FONTSIZE)
    ax.text(0.965, 0.90, f"n={values.size}", ha="right", va="top", transform=ax.transAxes, fontsize=PANEL_NOTE_FONTSIZE)
    if letter:
        ax.text(0.035, 0.93, letter, ha="left", va="top", transform=ax.transAxes, fontsize=PANEL_LETTER_FONTSIZE, fontweight="bold")
    if xlim is not None:
        ax.set_xlim(*xlim)
    if xticks is not None:
        ax.xaxis.set_major_locator(FixedLocator(xticks))
        ax.xaxis.set_major_formatter(FixedFormatter(xticklabels or [str(v) for v in xticks]))
        ax.xaxis.set_minor_formatter(NullFormatter())
    ax.grid(axis="y", linestyle="-", linewidth=0.45, alpha=0.22)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE, width=1.0, length=4.5)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def make_figures(df: pd.DataFrame):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "axes.titleweight": "regular",
        "axes.labelweight": "regular", "mathtext.fontset": "dejavusans",
        "font.size": TICK_LABEL_FONTSIZE,
    })
    # Grade A is the strict six-criterion subset excluding manually outlined SB27.
    highlight_mask = (df["pass_selection"] & ~df["id"].isin(A_GRADE_EXCLUDE_IDS)).to_numpy()
    n_hi = int(highlight_mask.sum())
    n_total = len(df)

    panels = [
        ("selected_cloud_count", "Selected molecular clouds", r"$N_{\rm cloud}$",
         {"bins": 9, "zero_floor": True, "log_x": True, "xlim": (4.5, 175),
          "xticks": [5, 10, 20, 50, 100, 150], "xticklabels": ["5", "10", "20", "50", "100", "150"], "letter": "(a)"}),
        ("xy_plane_R_eq_pc", "XY section equivalent radius", r"$R_{\rm eq,XY}$ [pc]", {"letter": "(b)"}),
        ("fit_std_pc", "Shell-fitting residual", r"$\sigma_{\rm fit}$ [pc]",
         {"threshold": SEL_FIT_STD_MAX_PC, "bins": 12, "zero_floor": True, "letter": "(c)"}),
        ("shell_radius_mean_abs_offset_over_R_eq", "Radial offset of shell ridge",
         r"$\langle |u_{\rm peak}-1| \rangle$",
         {"threshold": SEL_MEAN_ABS_OFFSET_MAX, "bins": 12, "letter": "(d)"}),
        ("shell_ridge_uPeak_curvature_rms", "Directional shell-ridge curvature",
         r"$C_u$",
         {"threshold": SEL_CURVATURE_MAX, "bins": 12, "letter": "(e)"}),
        ("dust_shell_ridge_peak_inner_ratio", "Shell-to-cavity density contrast",
         r"$\rho_{\rm ridge}/\rho_{\rm in}$",
         {"threshold": SEL_RIDGE_RATIO_MIN, "bins": 12, "zero_floor": True, "letter": "(f)"}),
        ("dust_inner_mean_mag_kpc", "Mean dust density inside cavity",
         r"$\rho_{\rm in}$ [mag kpc$^{-1}$]",
         {"threshold": SEL_INNER_MEAN_MAX, "bins": 12, "letter": "(g)"}),
        ("dust_inner_std_mag_kpc", "Cavity dust-density dispersion",
         r"$\sigma_{\rm in}$ [mag kpc$^{-1}$]",
         {"threshold": SEL_INNER_STD_MAX, "bins": 12, "letter": "(h)"}),
    ]

    fig, axes = plt.subplots(3, 3, figsize=A4_LANDSCAPE_FIGSIZE, constrained_layout=True)
    plot_axes = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1], axes[1, 2], axes[2, 0], axes[2, 1], axes[2, 2]]
    panel_b_fit_text = None
    for ax, (column, title, xlabel, kwargs) in zip(plot_axes, panels):
        if column == "xy_plane_R_eq_pc":
            panel_b_fit_text = plot_panel_b_powerlaw(ax, df, highlight_mask, letter=kwargs.get("letter", "(b)"))
            continue
        plot_hist(ax, df, column, title, xlabel, highlight_mask=highlight_mask, **kwargs)

    legend_ax = axes[0, 2]
    legend_ax.axis("off")
    handles = [
        plt.Line2D([0], [0], color=BLUE, linewidth=8, label=f"Full sample (n={n_total})"),
        plt.Line2D([0], [0], color=HIGHLIGHT_COLOR, linewidth=8, label=f"Class A targets (n={n_hi})"),
        plt.Line2D([0], [0], color="#2CA02C", linewidth=2.3, label="Median"),
        plt.Line2D([0], [0], color="#D62728", linestyle="--", linewidth=2.1, label="Selection cut"),
    ]
    legend_ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=LEGEND_FONTSIZE, handlelength=2.8)
    if panel_b_fit_text:
        legend_ax.text(
            0.03, 0.25,
            f"Panel (b) power-law fit:\n{panel_b_fit_text}",
            ha="left", va="top", fontsize=PANEL_B_FIT_NOTE_FONTSIZE, linespacing=1.45,
        )
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=300)
    plt.close(fig)


    fig2, ax2 = plt.subplots(1, 1, figsize=(8.4, 6.2), constrained_layout=True)
    plot_hist(ax2, df, "shell_ridge_uPeak_curvature_rms",
              r"Ridge $u_{\rm peak}$ direction-to-direction curvature",
              r"$u_{\rm peak}$ curvature RMS  (smaller = smoother)",
              bins=12, letter="(i)", highlight_mask=highlight_mask)
    fig2.suptitle("Smoothness of shell ridge across adjacent sight-lines", fontsize=SINGLE_FIG_TITLE_FONTSIZE)
    fig2.savefig(OUT_FIG_CURV, dpi=300, bbox_inches="tight")
    plt.close(fig2)


def find_geometry_csv() -> Path:
    for p in GEOM_CANDIDATES:
        if p.exists():
            return p
    raise FileNotFoundError(
        "No final geometry parameter table was found. Checked: "
        + " / ".join(str(p) for p in GEOM_CANDIDATES)
    )


def main():
    """Run Script 5 from validated inputs to the documented outputs."""
    geom_csv = find_geometry_csv()
    if not DUST_PARQUET.exists():
        raise FileNotFoundError('Missing input file.')
    print(f"Using geometry table: {geom_csv}")
    print(f"Using dust cube: {DUST_PARQUET}")

    geom_df = pd.read_csv(geom_csv)
    dust_df = compute_all_dust_metrics(geom_df)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)
    dust_df.to_csv(OUT_DUST_CSV, index=False, encoding="utf-8-sig")

    df = prepare_table(geom_df, dust_df)
    df.to_csv(OUT_HIST_CSV, index=False, encoding="utf-8-sig")
    make_figures(df)

    sel = df[df["pass_selection"]].sort_values("id")
    print(f"\nsaved: {OUT_DUST_CSV}")
    print(f"saved: {OUT_HIST_CSV}")
    print(f"saved: {OUT_FIG}")
    print(f"saved: {OUT_FIG_CURV}")
    print(f"Selected OSBs passing ridge criteria: {len(sel)} / {len(df)}")

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
