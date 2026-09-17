"""Script 17: SB31 RW boundary YSO traceback

Purpose
-------
Produce the SB31/Radcliffe-Wave boundary traceback supplement for G1 group 02 and young RW clusters.

Method overview
---------------
1. Load G1 family labels and orbit tracks, select group 02, and prepare the young
   Konietzka sample used in the comparison.
2. Apply the common median-velocity subtraction for the displayed traceback and project
   the shell overlays into the plotted slices.
3. Write the SB31/Radcliffe-Wave boundary comparison and the accompanying
   velocity/caption information.

Main inputs
-----------
- ../results/intermediate_output/14_g1_traceback_clustering/G1_family_labels.csv
- ../results/intermediate_output/14_g1_traceback_clustering/G1_orbits_55Myr.npz
- ../data/star_cluster_data/Konietzka2023.csv

Main outputs
------------
- ../results/figures/traceback_cluster_figures/17_g1_family2_median_subtracted_traceback_slices.png; SB31/Radcliffe-Wave traceback figure.

Figure/table role
-----------------
../results/figures/traceback_cluster_figures/17_g1_family2_median_subtracted_traceback_slices.png; SB31/Radcliffe-Wave traceback figure.

Runtime and data notes
----------------------
Reference runtime: 1.96 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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
import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from astropy import units as u
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from matplotlib.ticker import MultipleLocator


HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "17_g1_family2_traceback_supplement"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures" / "traceback_cluster_figures"
ELLIPSE_CSV = FINAL_DIR / "superbubble_final_fit_parameters.csv"

G1_DIR = Path("..") / "results" / "intermediate_output" / "14_g1_traceback_clustering"
LABELS_FILE = G1_DIR / "G1_family_labels.csv"
ORBIT_FILE = G1_DIR / "G1_orbits_55Myr.npz"
KONIETZKA_FILE = DATA_DIR / "star_cluster_data" / "Konietzka2023.csv"

OUTPUT_FIG = FINAL_FIG_DIR / "17_g1_family2_median_subtracted_traceback_slices.png"
OUTPUT_SUMMARY = OUT_DIR / "G1_group_02_median_velocity_subtracted_traceback_slice_statistics.csv"
OUTPUT_DETAIL = OUT_DIR / "G1_group_02_median_velocity_subtracted_traceback_slice_details.csv"
OUTPUT_CAPTION = OUT_DIR / "G1_group_02_median_velocity_subtracted_traceback_slice_caption.md"

RO_KPC = 8.122
VO_KMS = 236.0
ZO_KPC = 0.0208
DEFAULT_FAMILY = "group_02"
DEFAULT_STEP_MYR = 5.0
DEFAULT_MAX_LOOKBACK_MYR = 55.0
DEFAULT_XY_LIMIT_PC = 1000.0
DEFAULT_KONIETZKA_X_LIMIT_PC = -250.0
DEFAULT_KONIETZKA_Y_LIMIT_PC = -250.0
VELOCITY_SCALE = 0.045
REFERENCE_VELOCITY_KMS = 10.0
SOLAR_MOTION_LSR_STANDARD = (11.1, 12.24, 7.25)
ELLIPSOID_CAP_SECTION_FRAC = 0.5
DASHED_MAX_OPEN_IDS = {31, 25, 15, 16, 11}


def configure_chinese_plot_style() -> None:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei",
        "SimHei",
        "Noto Sans CJK SC",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    available = {font.name for font in font_manager.fontManager.ttflist}
    for candidate in ("Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "KaiTi", "STKaiti"):
        if candidate in available:
            plt.rcParams["font.family"] = candidate
            break
    plt.rcParams["axes.unicode_minus"] = False


def midcap_section_scale(best_fit: pd.DataFrame) -> np.ndarray:
    shape = best_fit["shape"].astype(str).str.lower().to_numpy()
    cz = best_fit["center_z_kpc"].to_numpy(float)
    c = best_fit["c_radius_kpc"].to_numpy(float)
    z_base = best_fit["xy_plane_z_kpc"].to_numpy(float)
    mark = best_fit["mark"].to_numpy(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        k_base = np.sqrt(np.clip(1.0 - ((z_base - cz) / c) ** 2, 0.0, None))
        z_apex = np.where(mark == 1, cz - c, cz + c)
        z_mid = z_base + ELLIPSOID_CAP_SECTION_FRAC * (z_apex - z_base)
        k_mid = np.sqrt(np.clip(1.0 - ((z_mid - cz) / c) ** 2, 0.0, None))
        scale = np.where(k_base > 0, k_mid / k_base, 1.0)
    return np.where(shape == "ellipsoid", scale, 1.0)


def load_compromise_ellipses(ellipse_csv: Path) -> pd.DataFrame:
    required = [
        "id",
        "mark",
        "shape",
        "xy_plane_center_x_kpc",
        "xy_plane_center_y_kpc",
        "xy_plane_a_kpc",
        "xy_plane_b_kpc",
        "xy_plane_angle_deg",
        "center_z_kpc",
        "c_radius_kpc",
        "xy_plane_z_kpc",
    ]
    best_fit = pd.read_csv(require_file(ellipse_csv))
    missing = [col for col in required if col not in best_fit.columns]
    if missing:
        raise ValueError('Invalid input or missing required data.')

    for col in required:
        if col != "shape":
            best_fit[col] = pd.to_numeric(best_fit[col], errors="coerce")
    scale = midcap_section_scale(best_fit)
    ellipse = pd.DataFrame(
        {
            "index": best_fit["id"].astype("Int64"),
            "mark": best_fit["mark"].astype("Int64"),
            "shape": best_fit["shape"].astype(str),
            "x_kpc": best_fit["xy_plane_center_x_kpc"],
            "y_kpc": best_fit["xy_plane_center_y_kpc"],
            "a_maxopen_kpc": best_fit["xy_plane_a_kpc"],
            "b_maxopen_kpc": best_fit["xy_plane_b_kpc"],
            "angle_deg": best_fit["xy_plane_angle_deg"],
            "section_scale": scale,
        }
    )
    ellipse["a_compromise_kpc"] = ellipse["a_maxopen_kpc"] * ellipse["section_scale"]
    ellipse["b_compromise_kpc"] = ellipse["b_maxopen_kpc"] * ellipse["section_scale"]
    return ellipse.dropna(subset=["x_kpc", "y_kpc", "a_compromise_kpc", "b_compromise_kpc", "angle_deg"]).copy()


def require_file(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path}")
    return path


def get_positions_summary(orbit_file: Path) -> tuple[np.ndarray, np.ndarray]:
    payload = np.load(orbit_file)
    lookback_myr = payload["lookback_myr"].astype(np.float64)
    if "positions_median_pc" in payload:
        positions_pc = payload["positions_median_pc"].astype(np.float64)
    else:
        positions_pc = np.median(payload["positions_pc"].astype(np.float64), axis=0)
    return lookback_myr, positions_pc


def compute_lsr_track(lookback_myr: np.ndarray) -> np.ndarray:
    times = (-lookback_myr) * u.Myr
    orbit = Orbit(vxvv=[1.0, 0.0, 1.0, 0.0, 0.0, 0.0], ro=RO_KPC, vo=VO_KMS, zo=ZO_KPC)
    orbit.integrate(times, MWPotential2014, method="odeint")
    x_pc = -np.asarray(orbit.x(times, use_physical=True)) * 1000.0
    y_pc = np.asarray(orbit.y(times, use_physical=True)) * 1000.0
    z_pc = np.asarray(orbit.z(times, use_physical=True)) * 1000.0
    return np.stack([x_pc, y_pc, z_pc], axis=1)


def cartesian_to_galactic_angles(x_pc: np.ndarray, y_pc: np.ndarray, z_pc: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dist_pc = np.sqrt(x_pc**2 + y_pc**2 + z_pc**2)
    l_deg = np.degrees(np.arctan2(y_pc, x_pc)) % 360.0
    b_deg = np.degrees(np.arcsin(np.divide(z_pc, dist_pc, out=np.zeros_like(z_pc), where=dist_pc > 0.0)))
    return l_deg, b_deg, dist_pc / 1000.0


def integrate_subset_orbits(df: pd.DataFrame, lookback_myr: np.ndarray) -> np.ndarray:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    times = (-lookback_myr) * u.Myr
    positions = np.empty((len(lookback_myr), len(df), 3), dtype=np.float64)
    for cluster_index, row in df.reset_index(drop=True).iterrows():
        x_helio = float(row["x_helio"])
        y_helio = float(row["y_helio"])
        z_helio = float(row["z_helio"])
        u_vel = float(row["U"])
        v_vel = float(row["V"])
        w_vel = float(row["W"])
        l_deg, b_deg, dist_kpc = cartesian_to_galactic_angles(
            np.asarray([x_helio]),
            np.asarray([y_helio]),
            np.asarray([z_helio]),
        )
        orbit = Orbit(
            vxvv=[l_deg[0], b_deg[0], dist_kpc[0], u_vel, v_vel, w_vel],
            lb=True,
            uvw=True,
            ro=RO_KPC,
            vo=VO_KMS,
            zo=ZO_KPC,
            solarmotion=[11.1, 12.24, 7.25],
        )
        orbit.integrate(times, MWPotential2014, method="odeint")
        positions[:, cluster_index, 0] = -np.asarray(orbit.x(times, use_physical=True)) * 1000.0
        positions[:, cluster_index, 1] = np.asarray(orbit.y(times, use_physical=True)) * 1000.0
        positions[:, cluster_index, 2] = np.asarray(orbit.z(times, use_physical=True)) * 1000.0
    return positions


def choose_time_indices(lookback_myr: np.ndarray, step_myr: float, max_lookback_myr: float) -> tuple[list[float], list[int]]:
    requested_times = np.arange(0.0, max_lookback_myr + 0.5 * step_myr, step_myr)
    indices = [int(np.argmin(np.abs(lookback_myr - value))) for value in requested_times]
    actual_times = [float(lookback_myr[index]) for index in indices]
    return actual_times, indices


def prepare_family_detail(labels_df: pd.DataFrame, family: str) -> pd.DataFrame:
    family_col = labels_df["inferred_family"].fillna("none").astype(str)
    selected = labels_df.loc[family_col.eq(family)].copy().reset_index()
    if selected.empty:
        raise ValueError('Invalid input or missing required data.')

    for column in ["x_helio", "y_helio", "z_helio", "U_lsr", "V_lsr", "W_lsr", "age_myr"]:
        selected[column] = pd.to_numeric(selected[column], errors="coerce")
    selected["diagnostic_group"] = family
    selected["plot_group"] = np.where(
        selected["source_catalog"].astype(str).eq("konietzka2023"),
        f"{family} Konietzkamembers",
        f"{family} Huntmembers",
    )
    selected["plot_color"] = "#1f77b4"
    selected["plot_marker"] = np.where(selected["source_catalog"].astype(str).eq("konietzka2023"), "s", "o")
    return selected


def normalize_name(value: object) -> str:
    return "".join(ch for ch in str(value).lower() if ch.isalnum())


def prepare_konietzka_left_lower(
    konietzka_file: Path,
    group_detail: pd.DataFrame,
    x_limit_pc: float,
    y_limit_pc: float,
) -> pd.DataFrame:
    raw = pd.read_csv(require_file(konietzka_file))
    for column in ["x", "y", "z", "vx", "vy", "vz", "age"]:
        raw[column] = pd.to_numeric(raw[column], errors="coerce")
    selected = raw.loc[raw["x"].lt(x_limit_pc) & raw["y"].lt(y_limit_pc)].copy()

    group_names = {normalize_name(value) for value in group_detail["name"]}
    selected = selected.loc[~selected["name"].map(normalize_name).isin(group_names)].copy().reset_index(drop=True)
    result = pd.DataFrame(
        {
            "index": -1,
            "name": selected["name"].astype(str),
            "display_name": selected["name"].astype(str),
            "source_catalog": "konietzka2023_left_lower",
            "reference_literature": selected.get("reference", ""),
            "age_myr": selected["age"],
            "x_helio": selected["x"],
            "y_helio": selected["y"],
            "z_helio": selected["z"],
            "U_lsr": selected["vx"],
            "V_lsr": selected["vy"],
            "W_lsr": selected["vz"],
            "U": selected["vx"] - SOLAR_MOTION_LSR_STANDARD[0],
            "V": selected["vy"] - SOLAR_MOTION_LSR_STANDARD[1],
            "W": selected["vz"] - SOLAR_MOTION_LSR_STANDARD[2],
            "inferred_family": "none",
            "diagnostic_group": "Konietzka lower-left quadrant",
            "plot_group": "Konietzka lower-left quadrant",
            "plot_color": "#d62728",
            "plot_marker": "s",
        }
    )
    return result.dropna(subset=["x_helio", "y_helio", "z_helio", "U_lsr", "V_lsr", "W_lsr"]).reset_index(drop=True)


def apply_common_median_velocity(detail: pd.DataFrame) -> pd.DataFrame:
    result = detail.copy()
    median_u = float(result["U_lsr"].median())
    median_v = float(result["V_lsr"].median())
    median_w = float(result["W_lsr"].median())
    result["common_median_U_lsr_kms"] = median_u
    result["common_median_V_lsr_kms"] = median_v
    result["common_median_W_lsr_kms"] = median_w
    result["dvx_lsr_kms"] = result["U_lsr"] - median_u
    result["dvy_lsr_kms"] = result["V_lsr"] - median_v
    result["dvz_lsr_kms"] = result["W_lsr"] - median_w
    return result


def add_velocity_scale(ax: plt.Axes, xy_limit_pc: float, view_center_xy: np.ndarray) -> None:
    span = 2.0 * xy_limit_pc
    x_min = float(view_center_xy[0]) - xy_limit_pc
    y_min = float(view_center_xy[1]) - xy_limit_pc
    anchor_x = x_min + 0.08 * span
    anchor_y = y_min + 0.08 * span
    ax.quiver(
        anchor_x,
        anchor_y,
        REFERENCE_VELOCITY_KMS,
        0.0,
        angles="xy",
        scale_units="xy",
        scale=VELOCITY_SCALE,
        color="black",
        width=0.0044,
        headwidth=4.1,
        headlength=5.0,
        headaxislength=4.3,
        zorder=9,
    )
    ax.text(anchor_x + 120.0, anchor_y + 35.0, f"{REFERENCE_VELOCITY_KMS:g} km/s", fontsize=7.5, color="black")


def add_subtracted_median_velocity_arrow(
    ax: plt.Axes,
    xy_limit_pc: float,
    view_center_xy: np.ndarray,
    median_u: float,
    median_v: float,
) -> None:
    span = 2.0 * xy_limit_pc
    x_max = float(view_center_xy[0]) + xy_limit_pc
    y_max = float(view_center_xy[1]) + xy_limit_pc
    anchor_x = x_max - 0.23 * span
    anchor_y = y_max - 0.12 * span
    ax.quiver(
        anchor_x,
        anchor_y,
        median_u,
        median_v,
        angles="xy",
        scale_units="xy",
        scale=VELOCITY_SCALE,
        color="#c1121f",
        alpha=0.96,
        width=0.0054,
        headwidth=4.6,
        headlength=5.6,
        headaxislength=4.8,
        zorder=12,
    )


def compute_slice_view_center(group_center_xy: np.ndarray, comparison_center_xy: np.ndarray) -> np.ndarray:
    centers = [group_center_xy]
    if np.all(np.isfinite(comparison_center_xy)):
        centers.append(comparison_center_xy)
    stacked = np.vstack(centers)
    return np.nanmean(stacked, axis=0)


def ellipse_intersects_view(row: pd.Series, a_column: str, b_column: str, view_center_xy: np.ndarray, xy_limit_pc: float) -> bool:
    x0 = float(row["x_kpc"]) * 1000.0
    y0 = float(row["y_kpc"]) * 1000.0
    radius = float(max(row[a_column], row[b_column])) * 1000.0
    x_min = float(view_center_xy[0]) - xy_limit_pc
    x_max = float(view_center_xy[0]) + xy_limit_pc
    y_min = float(view_center_xy[1]) - xy_limit_pc
    y_max = float(view_center_xy[1]) + xy_limit_pc
    return x0 + radius >= x_min and x0 - radius <= x_max and y0 + radius >= y_min and y0 - radius <= y_max


def draw_xy_superbubble_ellipse(
    ax: plt.Axes,
    row: pd.Series,
    a_column: str,
    b_column: str,
    linestyle: str,
    alpha: float,
    zorder: int,
) -> None:
    x0 = float(row["x_kpc"]) * 1000.0
    y0 = float(row["y_kpc"]) * 1000.0
    semi_major = float(row[a_column]) * 1000.0
    semi_minor = float(row[b_column]) * 1000.0
    ax.add_patch(
        Ellipse(
            (x0, y0),
            width=2.0 * semi_major,
            height=2.0 * semi_minor,
            angle=float(row["angle_deg"]),
            facecolor="none",
            edgecolor="black",
            linestyle=linestyle,
            linewidth=0.85,
            alpha=alpha,
            zorder=zorder,
        )
    )


def draw_current_slice_superbubble_ellipses(
    ax: plt.Axes,
    ellipse_df: pd.DataFrame,
    view_center_xy: np.ndarray,
    xy_limit_pc: float,
) -> None:
    x_min = float(view_center_xy[0]) - xy_limit_pc
    x_max = float(view_center_xy[0]) + xy_limit_pc
    y_min = float(view_center_xy[1]) - xy_limit_pc
    y_max = float(view_center_xy[1]) + xy_limit_pc
    for _, row in ellipse_df.iterrows():
        bubble_id = int(row["index"])
        if bubble_id in DASHED_MAX_OPEN_IDS and ellipse_intersects_view(row, "a_maxopen_kpc", "b_maxopen_kpc", view_center_xy, xy_limit_pc):
            draw_xy_superbubble_ellipse(ax, row, "a_maxopen_kpc", "b_maxopen_kpc", (0, (4, 3)), 0.45, 1)
        if ellipse_intersects_view(row, "a_compromise_kpc", "b_compromise_kpc", view_center_xy, xy_limit_pc):
            draw_xy_superbubble_ellipse(ax, row, "a_compromise_kpc", "b_compromise_kpc", "-", 0.68, 2)
            x0 = float(row["x_kpc"]) * 1000.0
            y0 = float(row["y_kpc"]) * 1000.0
            if x_min <= x0 <= x_max and y_min <= y0 <= y_max:
                ax.text(
                    x0,
                    y0,
                    str(bubble_id),
                    color="black",
                    fontsize=6.4,
                    ha="center",
                    va="center",
                    fontweight="bold",
                    clip_on=True,
                    path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
                    zorder=3,
                )


def scatter_source_subset(
    ax: plt.Axes,
    subset: pd.DataFrame,
    coords_xy: np.ndarray,
    velocities_xy: np.ndarray,
    color: str,
    marker: str,
    alpha: float,
    zorder: int,
    draw_vectors: bool = True,
    fill_marker: bool = True,
) -> None:
    if subset.empty:
        return
    row_positions = subset["_row_position"].to_numpy(dtype=int)
    scatter_kwargs = {
        "s": 34,
        "marker": marker,
        "alpha": alpha,
        "zorder": zorder,
    }
    if fill_marker:
        scatter_kwargs.update({"c": color, "edgecolors": "white", "linewidths": 0.25})
    else:
        scatter_kwargs.update({"facecolors": "none", "edgecolors": color, "linewidths": 0.9})
    ax.scatter(coords_xy[row_positions, 0], coords_xy[row_positions, 1], **scatter_kwargs)
    if not draw_vectors:
        return
    ax.quiver(
        coords_xy[row_positions, 0],
        coords_xy[row_positions, 1],
        velocities_xy[row_positions, 0],
        velocities_xy[row_positions, 1],
        angles="xy",
        scale_units="xy",
        scale=VELOCITY_SCALE,
        color=color,
        alpha=min(0.96, alpha + 0.08),
        width=0.0044,
        headwidth=4.1,
        headlength=5.0,
        headaxislength=4.3,
        zorder=zorder + 1,
    )


def plot_family_traceback(
    labels_file: Path,
    orbit_file: Path,
    konietzka_file: Path,
    family: str,
    output_fig: Path,
    output_summary: Path,
    output_detail: Path,
    output_caption: Path,
    step_myr: float,
    max_lookback_myr: float,
    xy_limit_pc: float,
    konietzka_x_limit_pc: float,
    konietzka_y_limit_pc: float,
    apply_age_cutoff: bool,
) -> None:
    configure_chinese_plot_style()
    labels_df = pd.read_csv(require_file(labels_file))
    lookback_myr, positions_pc = get_positions_summary(require_file(orbit_file))
    if positions_pc.shape[1] != len(labels_df):
        raise ValueError('Invalid input or missing required data.')
    ellipse_df = load_compromise_ellipses(ELLIPSE_CSV)

    group_detail = prepare_family_detail(labels_df, family)
    comparison_detail = prepare_konietzka_left_lower(
        konietzka_file=konietzka_file,
        group_detail=group_detail,
        x_limit_pc=konietzka_x_limit_pc,
        y_limit_pc=konietzka_y_limit_pc,
    )
    detail = pd.concat([group_detail, comparison_detail], ignore_index=True)
    detail = apply_common_median_velocity(detail)
    group_n = int(len(group_detail))
    comparison_n = int(len(comparison_detail))
    median_u = float(detail["common_median_U_lsr_kms"].iloc[0])
    median_v = float(detail["common_median_V_lsr_kms"].iloc[0])
    median_w = float(detail["common_median_W_lsr_kms"].iloc[0])

    lsr_track = compute_lsr_track(lookback_myr)
    rel_positions = positions_pc.copy()
    rel_positions[..., 0] -= lsr_track[:, None, 0]
    rel_positions[..., 1] -= lsr_track[:, None, 1]
    rel_positions[..., 2] -= lsr_track[:, None, 2]
    group_indices = group_detail["index"].to_numpy(dtype=int)
    group_positions = rel_positions[:, group_indices, :]
    comparison_positions = integrate_subset_orbits(comparison_detail, lookback_myr)
    comparison_positions[..., 0] -= lsr_track[:, None, 0]
    comparison_positions[..., 1] -= lsr_track[:, None, 1]
    comparison_positions[..., 2] -= lsr_track[:, None, 2]
    plotted_positions = np.concatenate([group_positions, comparison_positions], axis=1)

    detail["_row_position"] = np.arange(len(detail), dtype=int)
    residual_velocities_xy = detail[["dvx_lsr_kms", "dvy_lsr_kms"]].to_numpy(dtype=float)
    ages = detail["age_myr"].to_numpy(dtype=float)

    actual_times, indices = choose_time_indices(lookback_myr, step_myr=step_myr, max_lookback_myr=max_lookback_myr)
    n_panels = len(indices)
    ncols = 4
    nrows = math.ceil(n_panels / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.35 * ncols, 4.0 * nrows))
    axes = np.atleast_1d(axes).ravel()
    summary_rows: list[dict[str, float | int | str]] = []

    for panel_index, (ax, time_value, time_index) in enumerate(zip(axes, actual_times, indices, strict=False)):
        coords = plotted_positions[time_index, :, :]
        if apply_age_cutoff:
            active_mask = np.isfinite(ages) & (time_value <= ages)
        else:
            active_mask = np.ones(len(detail), dtype=bool)
        expired_mask = ~active_mask

        draw_df = detail.copy()
        draw_df["_active"] = active_mask
        for plot_group, source_df in draw_df.groupby("plot_group", sort=True):
            active_df = source_df.loc[source_df["_active"]]
            expired_df = source_df.loc[~source_df["_active"]]
            scatter_source_subset(
                ax,
                expired_df,
                coords[:, :2],
                residual_velocities_xy,
                color="0.62",
                marker=str(source_df["plot_marker"].iloc[0]),
                alpha=0.42,
                zorder=3,
                draw_vectors=True,
                fill_marker=True,
            )
            scatter_source_subset(
                ax,
                active_df,
                coords[:, :2],
                residual_velocities_xy,
                color=str(source_df["plot_color"].iloc[0]),
                marker=str(source_df["plot_marker"].iloc[0]),
                alpha=0.86,
                zorder=5,
                draw_vectors=True,
                fill_marker=True,
            )

        group_coords = coords[:group_n, :]
        comparison_coords = coords[group_n:, :]
        center_xy = np.nanmedian(group_coords[:, :2], axis=0)
        if len(comparison_coords):
            comparison_center_xy = np.nanmedian(comparison_coords[:, :2], axis=0)
            separation = float(np.linalg.norm(center_xy - comparison_center_xy))
        else:
            comparison_center_xy = np.asarray([np.nan, np.nan])
            separation = np.nan
        view_center_xy = compute_slice_view_center(center_xy, comparison_center_xy)

        title = "Present" if abs(time_value) < 1e-6 else f"Lookback {int(round(time_value))} Myr"
        ax.set_title(title, fontsize=11)
        ax.set_xlim(view_center_xy[0] - xy_limit_pc, view_center_xy[0] + xy_limit_pc)
        ax.set_ylim(view_center_xy[1] - xy_limit_pc, view_center_xy[1] + xy_limit_pc)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.18, lw=0.5)
        ax.xaxis.set_major_locator(MultipleLocator(500))
        ax.yaxis.set_major_locator(MultipleLocator(500))
        row_index = panel_index // ncols
        col_index = panel_index % ncols
        ax.tick_params(labelbottom=row_index == nrows - 1, labelleft=col_index == 0)
        if abs(time_value) < 1e-6:
            draw_current_slice_superbubble_ellipses(ax, ellipse_df, view_center_xy, xy_limit_pc)

        summary_rows.append(
            {
                "family": family,
                "lookback_myr": float(time_value),
                "n_group_members": group_n,
                "n_comparison_members": comparison_n,
                "n_active": int(np.sum(active_mask)),
                "n_expired": int(np.sum(expired_mask)),
                "group_median_x_pc": float(center_xy[0]),
                "group_median_y_pc": float(center_xy[1]),
                "comparison_median_x_pc": float(comparison_center_xy[0]),
                "comparison_median_y_pc": float(comparison_center_xy[1]),
                "median_xy_separation_pc": separation,
                "view_center_x_pc": float(view_center_xy[0]),
                "view_center_y_pc": float(view_center_xy[1]),
                "common_median_U_lsr_kms": median_u,
                "common_median_V_lsr_kms": median_v,
                "common_median_W_lsr_kms": median_w,
            }
        )

    for ax in axes[n_panels:]:
        ax.axis("off")

    handles = [
        Line2D([0], [0], marker="o", linestyle="None", color="w", markerfacecolor="#1f77b4", markersize=7, label="Hunt et al. (2024)"),
        Line2D([0], [0], marker="s", linestyle="None", color="w", markerfacecolor="#1f77b4", markersize=7, label="Konietzka et al. (2024)"),
        Line2D([0], [0], marker="s", linestyle="None", color="w", markerfacecolor="#d62728", markersize=7, label="SB31-Radcliffe Wave young clusters"),
        Line2D([0], [0], marker="o", linestyle="None", color="w", markerfacecolor="0.62", markeredgecolor="0.62", markersize=7, label="Older than their traceback age"),
        Line2D([0], [0], color="black", lw=1.0, label="Open superbubbles"),
    ]
    legend_fontsize = 15
    legend_ncol = min(3, len(handles))
    legend_reserved_bottom = 0.140
    fig.supxlabel("x' relative to the LSR (pc)", fontsize=12, y=0.108)
    fig.supylabel("y' relative to the LSR (pc)", fontsize=12, x=0.010)
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.012),
        ncol=legend_ncol,
        fontsize=legend_fontsize,
        frameon=True,
        columnspacing=0.9,
        handletextpad=0.32,
        borderaxespad=0.0,
    )
    fig.tight_layout(rect=(0.0, legend_reserved_bottom, 1.0, 0.985), pad=0.30, w_pad=0.25, h_pad=0.25)

    output_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_fig, dpi=220, bbox_inches="tight")
    plt.close(fig)

    output_summary.parent.mkdir(parents=True, exist_ok=True)
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(output_summary, index=False, encoding="utf-8-sig")

    detail_columns = [
        "name",
        "display_name",
        "source_catalog",
        "age_myr",
        "x_helio",
        "y_helio",
        "z_helio",
        "U_lsr",
        "V_lsr",
        "W_lsr",
        "inferred_family",
        "diagnostic_group",
        "common_median_U_lsr_kms",
        "common_median_V_lsr_kms",
        "common_median_W_lsr_kms",
        "dvx_lsr_kms",
        "dvy_lsr_kms",
        "dvz_lsr_kms",
    ]
    detail.loc[:, [column for column in detail_columns if column in detail.columns]].to_csv(
        output_detail,
        index=False,
        encoding="utf-8-sig",
    )

    caption = (
        "\\includegraphics[width=1\\linewidth]{../results/figures/traceback_cluster_figures/"
        "17_g1_family2_median_subtracted_traceback_slices.png}\n"
        "    \\caption{Traceback time slices for the G1 \\texttt{group\\_02} family and the "
        "SB31-related young-cluster subsample on the Radcliffe Wave. Each panel corresponds to one "
        "traceback time and shows heliocentric \\(X-Y\\) positions after subtracting the LSR trajectory. "
        "Red squares denote the young Radcliffe Wave clusters associated with SB31. Grey symbols indicate "
        "objects shown beyond their cluster ages. Arrows show present-day residual velocities, defined as "
        "each object's current LSR velocity minus the common median velocity of all objects in the figure.}\n"
    )
    output_caption.write_text(caption, encoding="utf-8")

    print(f"family={family}, group_members={group_n}, konietzka_left_lower={comparison_n}, total={len(detail)}")
    print(f"common median LSR velocity: U={median_u:.6g}, V={median_v:.6g}, W={median_w:.6g} km/s")
    print(f"saved figure: {output_fig}")
    print(f"saved summary: {output_summary}")
    print(f"saved detail: {output_detail}")
    print(f"saved caption: {output_caption}")


def main() -> None:
    """Run Script 17 from validated inputs to the documented outputs."""
    parser = argparse.ArgumentParser(description='Run Script 17: SB31 RW boundary YSO traceback.')
    parser.add_argument("--family", default=DEFAULT_FAMILY, help='Command-line option for the documented workflow.')
    parser.add_argument("--labels-file", default=str(LABELS_FILE), help='Input file used by this stage.')
    parser.add_argument("--orbit-file", default=str(ORBIT_FILE), help='Input file used by this stage.')
    parser.add_argument("--konietzka-file", default=str(KONIETZKA_FILE), help='Input file used by this stage.')
    parser.add_argument("--output-fig", default=str(OUTPUT_FIG), help='Output path for this product.')
    parser.add_argument("--output-summary", default=str(OUTPUT_SUMMARY), help='Output path for this product.')
    parser.add_argument("--output-detail", default=str(OUTPUT_DETAIL), help='Output path for this product.')
    parser.add_argument("--output-caption", default=str(OUTPUT_CAPTION), help='Output path for this product.')
    parser.add_argument("--step-myr", type=float, default=DEFAULT_STEP_MYR, help='Command-line option for the documented workflow.')
    parser.add_argument("--max-lookback-myr", type=float, default=DEFAULT_MAX_LOOKBACK_MYR, help='Command-line option for the documented workflow.')
    parser.add_argument("--xy-limit-pc", type=float, default=DEFAULT_XY_LIMIT_PC, help='Command-line option for the documented workflow.')
    parser.add_argument("--konietzka-x-limit-pc", type=float, default=DEFAULT_KONIETZKA_X_LIMIT_PC, help='Command-line option for the documented workflow.')
    parser.add_argument("--konietzka-y-limit-pc", type=float, default=DEFAULT_KONIETZKA_Y_LIMIT_PC, help='Command-line option for the documented workflow.')
    parser.add_argument("--no-age-cutoff", action="store_true", help='Command-line option for the documented workflow.')
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)
    plot_family_traceback(
        labels_file=Path(args.labels_file),
        orbit_file=Path(args.orbit_file),
        konietzka_file=Path(args.konietzka_file),
        family=args.family,
        output_fig=Path(args.output_fig),
        output_summary=Path(args.output_summary),
        output_detail=Path(args.output_detail),
        output_caption=Path(args.output_caption),
        step_myr=args.step_myr,
        max_lookback_myr=args.max_lookback_myr,
        xy_limit_pc=args.xy_limit_pc,
        konietzka_x_limit_pc=args.konietzka_x_limit_pc,
        konietzka_y_limit_pc=args.konietzka_y_limit_pc,
        apply_age_cutoff=not args.no_age_cutoff,
    )

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
