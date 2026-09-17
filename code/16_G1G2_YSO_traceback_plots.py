"""Script 16: G1G2 YSO traceback plots

Purpose
-------
Plot G1/G2 traceback slices, size-age curves, residual velocity fields, and combined family-SB panel.

Method overview
---------------
1. Load the family labels, orbit arrays, property tables and size evolution from Scripts
   14 and 15.
2. Compute the display summaries for backward-time slices, family size versus age, and
   residual velocity maps.
3. Render the G1/G2 figure set and selected superbubble/family comparison, including the
   HBW35 annotation from the bubble catalogue.

Main inputs
-----------
- ../results/intermediate_output/14_g1_traceback_clustering/G1_family_labels.csv
- ../results/intermediate_output/14_g1_traceback_clustering/G1_orbits_55Myr.npz
- ../results/intermediate_output/14_g1_traceback_clustering/G1_family_properties.csv
- ../results/intermediate_output/14_g1_traceback_clustering/G1_family_size_evolution.csv
- ../results/intermediate_output/15_g2_traceback_clustering/G2_family_labels.csv
- ../results/intermediate_output/15_g2_traceback_clustering/G2_orbits_50Myr.npz
- ../results/intermediate_output/15_g2_traceback_clustering/G2_family_properties.csv
- ../results/intermediate_output/15_g2_traceback_clustering/G2_family_size_evolution.csv
- ../results/superbubble_final_fit_parameters.csv
- ../data/Bubbles.csv (HBW35 overlay in the combined family-SB panels)

Main outputs
------------
- ../results/figures/traceback_cluster_figures/; Extended Data traceback panel and young-cluster figures.

Figure/table role
-----------------
../results/figures/traceback_cluster_figures/; Extended Data traceback panel and young-cluster figures.

Runtime and data notes
----------------------
Reference runtime: 11.83 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from matplotlib import colormaps, font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
from matplotlib.patches import Rectangle
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014


HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "16_traceback_cluster_common_plot_data"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures" / "traceback_cluster_figures"
ELLIPSE_CSV = FINAL_DIR / "superbubble_final_fit_parameters.csv"
OUTPUT_ELLIPSE = OUT_DIR / "superbubble_ellipses.csv"

G1_DIR = Path("..") / "results" / "intermediate_output" / "14_g1_traceback_clustering"
G2_DIR = Path("..") / "results" / "intermediate_output" / "15_g2_traceback_clustering"

RO_KPC = 8.122
VO_KMS = 236.0
ZO_KPC = 0.0208


XY_XLIM_PC = (-1500.0, 1500.0)
XY_YLIM_PC = (-1500.0, 1500.0)
SELECTED_XY_HALF_WIDTH_PC = 1000.0
LB_LIM_L_DEG = (0.0, 360.0)
LB_LIM_B_DEG = (-65.0, 35.0)
VELOCITY_SCALE = 0.032
REFERENCE_VELOCITY_KMS = 10.0
W_VELOCITY_SCALE = 0.325
REFERENCE_W_KMS = 1.0
REFERENCE_W_ARROW_POINTS = 34.0
LBW_VERTICAL_TRIM_DEG = 20.0
SELECTED_LEGEND_FONTSIZE = 10.0
SELECTED_SHARED_LEGEND_FONTSIZE = 14.0
SELECTED_SPEED_LABEL_FONTSIZE = 14.0
A4_TEXT_WIDTH_IN = 7.4
PRINT_BODY_FONTSIZE = 12.0
FIGURE_SIZE_SCALE = 1.5


def scaled_figsize(width_in: float, height_in: float) -> tuple[float, float]:
    return width_in * FIGURE_SIZE_SCALE, height_in * FIGURE_SIZE_SCALE
RW_ANGLE_DEG = 60.0
RW_INTERCEPT_KPC = 0.6
RW_HALF_WIDTH_PC = 100.0
CATALOG_BLUE = "#1f77b4"
RW_HIGHLIGHT_RED = "#d62728"
RW_CENTERLINE_COLOR = "#2ca02c"
ELLIPSOID_CAP_SECTION_FRAC = 0.5
HALF_ELLIPSOID_BUFFER_KPC = 0.08
DASHED_MAX_OPEN_IDS = {31, 25, 15, 16, 11}
HBW_BUBBLE_CSV = DATA_DIR / "Bubbles.csv"
HBW_BUBBLE_ID = "hbw35"
ORION_ERIDANUS_LABEL = "Orion-Eridanus superbubble"
ORION_ERIDANUS_XY_FAMILIES = {"G1_group_01", "G1_group_02"}
RW_CLUSTER_LABEL = "Radcliffe Wave-associated clusters"
NON_RW_CLUSTER_LABEL = "Non-Radcliffe Wave clusters"
RADCLIFFE_WAVE_LABEL = "Radcliffe Wave"
SB_LABEL_OFFSETS_PC = {
    15: (0.0, 20.0),
    16: (0.0, -20.0),
}

SOURCE_STYLES = {
    "hunt_stratified": {"color": CATALOG_BLUE, "marker": "o", "label": "Hunt et al. (2024)", "short": "H"},
    "hunt2024": {"color": CATALOG_BLUE, "marker": "o", "label": "Hunt et al. (2024)", "short": "H"},
    "konietzka2023": {"color": CATALOG_BLUE, "marker": "*", "label": "Konietzka et al. (2024)", "short": "K"},
    "li2025_yso_group": {"color": "#2ca02c", "marker": "s", "label": "Li YSO", "short": "Li"},
}
SOURCE_ORDER = ("hunt_stratified", "hunt2024", "konietzka2023", "li2025_yso_group")
DEFAULT_SOURCE_STYLE = {"color": "#9467bd", "marker": "D", "label": "Other", "short": "O"}

DETAIL_COLUMNS = [
    "name",
    "display_name",
    "source_catalog",
    "reference_literature",
    "age_myr",
    "x_pc",
    "y_pc",
    "z_pc",
    "vx_lsr_kms",
    "vy_lsr_kms",
    "vz_lsr_kms",
    "family",
    "final_cluster_label",
    "weight_score",
    "median_vx_lsr_kms",
    "median_vy_lsr_kms",
    "median_vz_lsr_kms",
    "dvx_lsr_kms",
    "dvy_lsr_kms",
    "dvz_lsr_kms",
    "rw_ds_pc",
    "in_rw_band",
]

plt.rcParams["font.sans-serif"] = [
    "Microsoft YaHei",
    "SimHei",
    "Noto Sans CJK SC",
    "Arial Unicode MS",
    "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False


def configure_chinese_plot_style() -> None:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    available = {font.name for font in font_manager.fontManager.ttflist}
    for candidate in ("Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "KaiTi", "STKaiti"):
        if candidate in available:
            plt.rcParams["font.family"] = candidate
            break
    plt.rcParams["axes.unicode_minus"] = False


def compute_family_age_summary(labels_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    family_df = labels_df[labels_df["inferred_family"].fillna("none") != "none"].copy()
    for family, sub in family_df.groupby("inferred_family", sort=True):
        ages = pd.to_numeric(sub["age_myr"], errors="coerce").dropna().to_numpy(dtype=float)
        rows.append(
            {
                "family": family,
                "age_mean_myr": float(np.mean(ages)) if len(ages) else np.nan,
                "age_median_myr": float(np.median(ages)) if len(ages) else np.nan,
                "age_std_myr": float(np.std(ages, ddof=1)) if len(ages) > 1 else np.nan,
                "age_n_valid": int(len(ages)),
            }
        )
    return pd.DataFrame(rows)


def interpolate_size_at_age(curve: pd.DataFrame, age_myr: float) -> float:
    if not np.isfinite(age_myr) or curve.empty:
        return np.nan
    x = curve["lookback_myr"].to_numpy(dtype=float)
    size_col = "size_pc" if "size_pc" in curve.columns else "median_radius_pc"
    y = curve[size_col].to_numpy(dtype=float)
    if age_myr < np.nanmin(x) or age_myr > np.nanmax(x):
        return np.nan
    return float(np.interp(age_myr, x, y))


def family_member_ages(labels_df: pd.DataFrame, family: str, max_lookback: float) -> np.ndarray:
    family_column = "inferred_family" if "inferred_family" in labels_df.columns else "family"
    if family_column not in labels_df.columns or "age_myr" not in labels_df.columns:
        return np.array([], dtype=float)
    mask = labels_df[family_column].astype(str) == str(family)
    ages = pd.to_numeric(labels_df.loc[mask, "age_myr"], errors="coerce").dropna().to_numpy(dtype=float)
    ages = ages[np.isfinite(ages)]
    ages = ages[(ages >= 0.0) & (ages <= max_lookback)]
    return np.sort(ages)


def draw_member_age_ticks(ax: plt.Axes, curve: pd.DataFrame, ages: np.ndarray, color: tuple[float, float, float, float]) -> None:
    if len(ages) == 0 or curve.empty:
        return
    size_col = "size_pc" if "size_pc" in curve.columns else "median_radius_pc"
    y_values = curve[size_col].to_numpy(dtype=float)
    y_span = float(np.nanmax(y_values) - np.nanmin(y_values))
    tick_half_height = max(4.0, 0.055 * y_span)
    for age in ages:
        size_at_age = interpolate_size_at_age(curve, float(age))
        if not np.isfinite(size_at_age):
            continue
        ax.plot(
            [age, age],
            [size_at_age - tick_half_height, size_at_age + tick_half_height],
            color=color,
            lw=1.45,
            alpha=0.86,
            solid_capstyle="butt",
            zorder=7,
        )


def plot_size_evolution_with_mean_age(
    properties_file: Path,
    curve_file: Path,
    labels_file: Path,
    output_file: Path,
) -> Path:
    configure_chinese_plot_style()
    properties_df = pd.read_csv(properties_file).sort_values(["n_members", "family"], ascending=[False, True])
    curve_df = pd.read_csv(curve_file)
    labels_df = pd.read_csv(labels_file)
    properties_df = properties_df.rename(
        columns={
            "compact_lookback_myr": "t_compact_myr",
            "compact_radius_pc": "d_compact_pc",
        }
    )
    if "median_radius_pc" in curve_df.columns and "size_pc" not in curve_df.columns:
        curve_df["size_pc"] = curve_df["median_radius_pc"]
    age_summary = compute_family_age_summary(labels_df)
    properties_df = properties_df.merge(age_summary, on="family", how="left", suffixes=("", "_labels"))
    if "age_mean_myr" not in properties_df.columns and "mean_age_myr" in properties_df.columns:
        properties_df["age_mean_myr"] = properties_df["mean_age_myr"]
    if "age_median_myr" not in properties_df.columns and "median_age_myr" in properties_df.columns:
        properties_df["age_median_myr"] = properties_df["median_age_myr"]

    families = properties_df["family"].tolist()
    n_panels = len(families)
    ncols = 4 if n_panels > 6 else 3
    nrows = math.ceil(n_panels / ncols)
    a4_text_width_in = 10.4 if ncols == 4 else 7.4
    panel_height_in = 3.05
    min_print_fontsize = 10.5
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(a4_text_width_in, 0.65 + panel_height_in * nrows),
        sharex=True,
    )
    axes = np.atleast_1d(axes).ravel()
    palette = colormaps.get_cmap("tab10")
    max_lookback = float(curve_df["lookback_myr"].max())

    for index, family in enumerate(families):
        ax = axes[index]
        color = palette(index % 10)
        family_curve = curve_df[curve_df["family"] == family].sort_values("lookback_myr")
        props = properties_df[properties_df["family"] == family].iloc[0]
        mean_age = float(props["age_mean_myr"]) if pd.notna(props["age_mean_myr"]) else np.nan
        median_age = float(props["age_median_myr"]) if pd.notna(props["age_median_myr"]) else np.nan
        size_at_mean_age = interpolate_size_at_age(family_curve, mean_age)
        member_ages = family_member_ages(labels_df, family, max_lookback)

        ax.plot(family_curve["lookback_myr"], family_curve["size_pc"], color=color, lw=2.0)
        draw_member_age_ticks(ax, family_curve, member_ages, color)
        ax.scatter(
            [props["t_compact_myr"]],
            [props["d_compact_pc"]],
            color=color,
            s=48,
            edgecolors="black",
            linewidths=0.45,
            zorder=12,
        )
        ax.axvline(props["t_compact_myr"], color=color, ls="--", lw=1.0, alpha=0.85, zorder=2)

        if np.isfinite(mean_age) and 0.0 <= mean_age <= max_lookback:
            ax.axvline(mean_age, color="black", ls=":", lw=1.2, alpha=0.75, zorder=3)
            if np.isfinite(size_at_mean_age):
                ax.scatter([mean_age], [size_at_mean_age], marker="D", s=42, c="black", zorder=13)

        ax.set_title(f"{family}  (N={int(props['n_members'])})", fontsize=11.0, pad=5.0)
        ax.grid(alpha=0.22, lw=0.5)
        ax.set_xlim(0, max_lookback)
        ax.tick_params(axis="both", which="major", labelsize=min_print_fontsize, length=4.0, width=0.8)
        annotation_times = np.array([props["t_compact_myr"], mean_age, median_age], dtype=float)
        finite_annotation_times = annotation_times[np.isfinite(annotation_times)]
        annotation_on_left = (
            len(finite_annotation_times) > 0
            and float(np.nanmax(finite_annotation_times)) > 0.68 * max_lookback
        )
        ax.text(
            0.04 if annotation_on_left else 0.96,
            0.96,
            (
                f"tc={props['t_compact_myr']:.1f} Myr\n"
                f"mean={mean_age:.1f} Myr\n"
                f"median={median_age:.1f} Myr"
            ),
            transform=ax.transAxes,
            ha="left" if annotation_on_left else "right",
            va="top",
            fontsize=min_print_fontsize,
            linespacing=1.05,
            bbox={"facecolor": "white", "alpha": 0.72, "edgecolor": "none", "pad": 1.8},
        )

    for ax in axes[n_panels:]:
        ax.axis("off")

    fig.supxlabel("Lookback time (Myr)", fontsize=min_print_fontsize, y=0.024)
    fig.supylabel("Family size D(t) (pc)", fontsize=min_print_fontsize, x=0.040)
    fig.subplots_adjust(left=0.095, right=0.985, top=0.965, bottom=0.100, wspace=0.30, hspace=0.36)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_file


def export_basic_family_catalog(labels_file: Path, output_file: Path) -> Path:
    labels_df = pd.read_csv(labels_file)
    columns = [
        "name",
        "display_name",
        "source_catalog",
        "reference_literature",
        "age_myr",
        "x_helio",
        "y_helio",
        "z_helio",
        "U_lsr",
        "V_lsr",
        "W_lsr",
        "inferred_family",
        "final_cluster_label",
        "weight_score",
    ]
    available = [column for column in columns if column in labels_df.columns]
    result = labels_df[available].copy()
    result = result.rename(
        columns={
            "x_helio": "x_pc",
            "y_helio": "y_pc",
            "z_helio": "z_pc",
            "U_lsr": "vx_lsr_kms",
            "V_lsr": "vy_lsr_kms",
            "W_lsr": "vz_lsr_kms",
            "inferred_family": "family",
        }
    )
    result["family"] = result["family"].fillna("none")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_file, index=False, encoding="utf-8-sig")
    return output_file


def get_positions_summary(orbit_file: Path) -> tuple[np.ndarray, np.ndarray]:
    payload = np.load(orbit_file)
    lookback_myr = payload["lookback_myr"].astype(np.float64)
    if "positions_median_pc" in payload:
        positions = payload["positions_median_pc"].astype(np.float64)
    else:
        positions = np.median(payload["positions_pc"].astype(np.float64), axis=0)
    return lookback_myr, positions


def compute_lsr_track(lookback_myr: np.ndarray) -> np.ndarray:
    times = (-lookback_myr) * u.Myr
    orbit = Orbit(vxvv=[1.0, 0.0, 1.0, 0.0, 0.0, 0.0], ro=RO_KPC, vo=VO_KMS, zo=ZO_KPC)
    orbit.integrate(times, MWPotential2014, method="odeint")
    x_pc = -np.asarray(orbit.x(times, use_physical=True)) * 1000.0
    y_pc = np.asarray(orbit.y(times, use_physical=True)) * 1000.0
    z_pc = np.asarray(orbit.z(times, use_physical=True)) * 1000.0
    return np.stack([x_pc, y_pc, z_pc], axis=1)


def choose_time_indices(lookback_myr: np.ndarray, step_myr: float, max_lookback_myr: float) -> tuple[list[float], list[int]]:
    requested_times = np.arange(0.0, max_lookback_myr + 0.5 * step_myr, step_myr)
    indices = [int(np.argmin(np.abs(lookback_myr - value))) for value in requested_times]
    actual_times = [float(lookback_myr[index]) for index in indices]
    return actual_times, indices


def plot_paperlike_families_across_slices(
    labels_file: str,
    orbit_file: str,
    output_file: str,
    step_myr: float,
    max_lookback_myr: float,
    xy_limit_pc: float,
    apply_age_cutoff: bool = True,
) -> Path:
    configure_chinese_plot_style()
    labels_df = pd.read_csv(labels_file)
    lookback_myr, positions_pc = get_positions_summary(Path(orbit_file))
    actual_times, indices = choose_time_indices(lookback_myr, step_myr=step_myr, max_lookback_myr=max_lookback_myr)

    ages = pd.to_numeric(labels_df.get("age_myr"), errors="coerce").to_numpy(dtype=np.float64)
    families = labels_df["inferred_family"].fillna("none").to_numpy()
    family_names = [name for name in sorted(pd.unique(families)) if name != "none"]
    family_codes = {name: index for index, name in enumerate(family_names)}
    palette = colormaps.get_cmap("tab10")

    lsr_track = compute_lsr_track(lookback_myr)
    rel_positions = positions_pc.copy()
    rel_positions[..., 0] -= lsr_track[:, None, 0]
    rel_positions[..., 1] -= lsr_track[:, None, 1]
    rel_positions[..., 2] -= lsr_track[:, None, 2]

    n_panels = len(indices)
    ncols = 4
    nrows = math.ceil(n_panels / ncols)
    min_print_fontsize = PRINT_BODY_FONTSIZE
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=scaled_figsize(A4_TEXT_WIDTH_IN, 0.35 + 1.78 * nrows),
        sharex=True,
        sharey=True,
    )
    axes = np.atleast_1d(axes).ravel()

    for ax, time_value, index in zip(axes, actual_times, indices, strict=False):
        coords = rel_positions[index]
        if apply_age_cutoff:
            visible_mask = np.isfinite(ages) & (time_value <= ages)
        else:
            visible_mask = np.ones(len(labels_df), dtype=bool)

        ax.scatter(
            coords[visible_mask, 0],
            coords[visible_mask, 1],
            s=10,
            c="black",
            alpha=0.18,
            edgecolors="none",
        )

        for family_name in family_names:
            family_mask = (families == family_name) & visible_mask
            if not np.any(family_mask):
                continue
            ax.scatter(
                coords[family_mask, 0],
                coords[family_mask, 1],
                s=22,
                c=[palette(family_codes[family_name] % 10)],
                alpha=0.86,
                edgecolors="none",
            )

        title = "Present" if abs(time_value) < 1e-6 else f"{int(round(time_value))} Myr ago"
        ax.set_title(title, fontsize=min_print_fontsize, pad=4.0)
        ax.set_xlim(-xy_limit_pc, xy_limit_pc)
        ax.set_ylim(-xy_limit_pc, xy_limit_pc)
        ax.set_xticks([-1000.0, 0.0, 1000.0])
        ax.set_yticks([-1000.0, 0.0, 1000.0])
        ax.set_xticklabels(["-1", "0", "1"])
        ax.set_yticklabels(["-1", "0", "1"])
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.18, lw=0.5)
        ax.tick_params(axis="both", which="major", labelsize=min_print_fontsize, length=4.0, width=0.8)

    for ax in axes[n_panels:]:
        ax.axis("off")

    fig.supxlabel("LSR-relative X' (kpc)", fontsize=min_print_fontsize, y=0.026)
    fig.supylabel("LSR-relative Y' (kpc)", fontsize=min_print_fontsize, x=0.045)
    fig.subplots_adjust(left=0.095, right=0.985, top=0.95, bottom=0.090, wspace=0.08, hspace=0.24)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


@dataclass(frozen=True)
class XYDatasetConfig:
    key: str
    title: str
    source_csv: Path
    output_fig: Path
    output_summary: Path
    output_detail: Path
    source_order: tuple[str, ...]


def safe_name(value: str) -> str:
    return re.sub(r"[^0-9A-Za-z_]+", "_", str(value)).strip("_").lower()


def style_for_source(source_catalog: str) -> dict[str, str]:
    return SOURCE_STYLES.get(str(source_catalog), DEFAULT_SOURCE_STYLE)


def annotate_radcliffe_wave_band(df: pd.DataFrame) -> pd.DataFrame:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    x_column = "x_pc" if "x_pc" in df.columns else "x_helio" if "x_helio" in df.columns else None
    y_column = "y_pc" if "y_pc" in df.columns else "y_helio" if "y_helio" in df.columns else None
    if x_column is None or y_column is None:
        return df
    theta = np.deg2rad(RW_ANGLE_DEG)
    x = pd.to_numeric(df[x_column], errors="coerce").to_numpy(float)
    y = pd.to_numeric(df[y_column], errors="coerce").to_numpy(float)
    intercept_pc = 1000.0 * RW_INTERCEPT_KPC
    ds_pc = -np.sin(theta) * x + np.cos(theta) * y - intercept_pc * np.cos(theta)
    df = df.copy()
    df["rw_ds_pc"] = ds_pc
    df["in_rw_band"] = np.isfinite(ds_pc) & (np.abs(ds_pc) < RW_HALF_WIDTH_PC)
    return df


def point_color_for_rw(subset: pd.DataFrame) -> np.ndarray:
    in_rw = subset.get("in_rw_band", False)
    if isinstance(in_rw, pd.Series):
        mask = in_rw.fillna(False).to_numpy(bool)
    else:
        mask = np.full(len(subset), bool(in_rw))
    return np.where(mask, RW_HIGHLIGHT_RED, CATALOG_BLUE)


def unclustered_mask(subset: pd.DataFrame) -> np.ndarray:
    for column in ("original_family", "family", "inferred_family"):
        if column in subset.columns:
            return subset[column].fillna("none").astype(str).eq("none").to_numpy(bool)
    return np.zeros(len(subset), dtype=bool)


def family_sort_key(value: str) -> int:
    if value.startswith("group_") and value.split("_")[-1].isdigit():
        return int(value.split("_")[-1])
    return 999


def load_catalog(source_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(source_csv)
    if "inferred_family" in df.columns:
        df["family"] = df["inferred_family"].fillna("none").astype(str)
    elif "family" in df.columns:
        df["family"] = df["family"].fillna("none").astype(str)
    else:
        raise ValueError('Invalid input or missing required data.')

    df = df.rename(
        columns={
            "x_helio": "x_pc",
            "y_helio": "y_pc",
            "z_helio": "z_pc",
            "U_lsr": "vx_lsr_kms",
            "V_lsr": "vy_lsr_kms",
            "W_lsr": "vz_lsr_kms",
        }
    )
    required = ["x_pc", "y_pc", "z_pc", "vx_lsr_kms", "vy_lsr_kms", "vz_lsr_kms", "source_catalog"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError('Invalid input or missing required data.')

    for col in ["x_pc", "y_pc", "z_pc", "vx_lsr_kms", "vy_lsr_kms", "vz_lsr_kms", "age_myr", "weight_score"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return annotate_radcliffe_wave_band(df.dropna(subset=["x_pc", "y_pc", "vx_lsr_kms", "vy_lsr_kms"]).copy())


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
    best_fit = pd.read_csv(ellipse_csv)
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
    ellipse["area_pct"] = 100.0 * ellipse["section_scale"] ** 2
    return ellipse.dropna(subset=["x_kpc", "y_kpc", "a_compromise_kpc", "b_compromise_kpc", "angle_deg"]).copy()


def prepare_residuals(catalog_df: pd.DataFrame, source_order: tuple[str, ...]) -> tuple[list[str], pd.DataFrame, pd.DataFrame]:
    family_labels = sorted([family for family in catalog_df["family"].unique() if family != "none"], key=family_sort_key)
    detail_frames = []
    summary_rows = []

    for family in family_labels:
        subset = catalog_df.loc[catalog_df["family"] == family].copy()
        median_vx = float(subset["vx_lsr_kms"].median())
        median_vy = float(subset["vy_lsr_kms"].median())
        median_vz = float(subset["vz_lsr_kms"].median())
        subset["median_vx_lsr_kms"] = median_vx
        subset["median_vy_lsr_kms"] = median_vy
        subset["median_vz_lsr_kms"] = median_vz
        subset["dvx_lsr_kms"] = subset["vx_lsr_kms"] - median_vx
        subset["dvy_lsr_kms"] = subset["vy_lsr_kms"] - median_vy
        subset["dvz_lsr_kms"] = subset["vz_lsr_kms"] - median_vz
        detail_frames.append(subset)

        source_counts = subset["source_catalog"].astype(str).value_counts()
        row = {
            "family": family,
            "n_members": int(len(subset)),
            "median_vx_lsr_kms": median_vx,
            "median_vy_lsr_kms": median_vy,
            "median_vz_lsr_kms": median_vz,
            "mean_age_myr": float(subset["age_myr"].mean()) if "age_myr" in subset and subset["age_myr"].notna().any() else np.nan,
            "median_age_myr": float(subset["age_myr"].median()) if "age_myr" in subset and subset["age_myr"].notna().any() else np.nan,
            "median_weight_score": float(subset["weight_score"].median()) if "weight_score" in subset and subset["weight_score"].notna().any() else np.nan,
        }
        for source in source_order:
            row[f"n_{safe_name(source)}"] = int(source_counts.get(source, 0))
        summary_rows.append(row)

    detail_df = pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame()
    summary_df = pd.DataFrame(summary_rows)
    return family_labels, detail_df, summary_df


def ellipse_may_intersect_view(row: pd.Series, a_column: str = "a_compromise_kpc", b_column: str = "b_compromise_kpc") -> bool:
    x0 = float(row["x_kpc"]) * 1000.0
    y0 = float(row["y_kpc"]) * 1000.0
    radius = float(max(row[a_column], row[b_column])) * 1000.0
    return (
        x0 + radius >= XY_XLIM_PC[0]
        and x0 - radius <= XY_XLIM_PC[1]
        and y0 + radius >= XY_YLIM_PC[0]
        and y0 - radius <= XY_YLIM_PC[1]
    )


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
            linewidth=0.75,
            alpha=alpha,
            zorder=zorder,
        )
    )


def draw_superbubble_ellipses(ax: plt.Axes, ellipse_df: pd.DataFrame) -> None:
    for _, row in ellipse_df.iterrows():
        bubble_id = int(row["index"])
        if bubble_id in DASHED_MAX_OPEN_IDS and ellipse_may_intersect_view(row, "a_maxopen_kpc", "b_maxopen_kpc"):
            draw_xy_superbubble_ellipse(ax, row, "a_maxopen_kpc", "b_maxopen_kpc", (0, (4, 3)), 0.45, 1)
        if ellipse_may_intersect_view(row, "a_compromise_kpc", "b_compromise_kpc"):
            draw_xy_superbubble_ellipse(ax, row, "a_compromise_kpc", "b_compromise_kpc", "-", 0.62, 2)
            x0 = float(row["x_kpc"]) * 1000.0
            y0 = float(row["y_kpc"]) * 1000.0
            if XY_XLIM_PC[0] <= x0 <= XY_XLIM_PC[1] and XY_YLIM_PC[0] <= y0 <= XY_YLIM_PC[1]:
                ax.text(
                    x0,
                    y0,
                    str(bubble_id),
                    color="black",
                    fontsize=6.5,
                    ha="center",
                    va="center",
                    fontweight="bold",
                    clip_on=True,
                    path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
                    zorder=3,
                )


def add_solar_marker(ax: plt.Axes) -> None:
    ax.scatter(0.0, 0.0, marker="*", s=120, color="gold", edgecolors="black", linewidths=0.7, zorder=8)


def add_velocity_scale(ax: plt.Axes) -> None:
    anchor_x = XY_XLIM_PC[0] + 0.08 * (XY_XLIM_PC[1] - XY_XLIM_PC[0])
    anchor_y = XY_YLIM_PC[0] + 0.08 * (XY_YLIM_PC[1] - XY_YLIM_PC[0])
    arrow_length_pc = REFERENCE_VELOCITY_KMS / VELOCITY_SCALE
    box_width = arrow_length_pc + 340.0
    box_height = 170.0
    ax.add_patch(
        Rectangle(
            (anchor_x - 45.0, anchor_y - 70.0),
            box_width,
            box_height,
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.6,
            alpha=0.78,
            zorder=7,
        )
    )
    ax.quiver(
        anchor_x,
        anchor_y,
        REFERENCE_VELOCITY_KMS,
        0.0,
        angles="xy",
        scale_units="xy",
        scale=VELOCITY_SCALE,
        color="black",
        width=0.0032,
        headwidth=4.2,
        headlength=5.0,
        headaxislength=4.4,
        zorder=8,
    )
    ax.text(
        anchor_x + arrow_length_pc + 35.0,
        anchor_y,
        f"{REFERENCE_VELOCITY_KMS:g} km/s",
        fontsize=PRINT_BODY_FONTSIZE,
        color="black",
        ha="left",
        va="center",
        zorder=9,
    )


def add_family_median_velocity_arrow(ax: plt.Axes, stat: pd.Series) -> None:
    vx = float(stat["median_vx_lsr_kms"])
    vy = float(stat["median_vy_lsr_kms"])
    dx_pc = vx / VELOCITY_SCALE
    dy_pc = vy / VELOCITY_SCALE
    x_span = XY_XLIM_PC[1] - XY_XLIM_PC[0]
    y_span = XY_YLIM_PC[1] - XY_YLIM_PC[0]
    box_left = XY_XLIM_PC[1] - 0.36 * x_span
    box_right = XY_XLIM_PC[1] - 0.06 * x_span
    box_bottom = XY_YLIM_PC[1] - 0.26 * y_span
    box_top = XY_YLIM_PC[1] - 0.08 * y_span
    tail_x = np.clip(box_right - max(dx_pc, 0.0), box_left - min(dx_pc, 0.0), box_right - max(dx_pc, 0.0))
    tail_y = np.clip(box_top - max(dy_pc, 0.0), box_bottom - min(dy_pc, 0.0), box_top - max(dy_pc, 0.0))
    ax.quiver(
        tail_x,
        tail_y,
        vx,
        vy,
        angles="xy",
        scale_units="xy",
        scale=VELOCITY_SCALE,
        color="#c1121f",
        width=0.0048,
        headwidth=4.8,
        headlength=5.8,
        headaxislength=5.0,
        zorder=10,
    )


def draw_family_members(
    ax: plt.Axes,
    subset: pd.DataFrame,
    source_order: tuple[str, ...],
    *,
    point_zorder: float = 6,
    arrow_zorder: float = 7,
    arrow_width: float = 0.0032,
) -> None:
    ordered_sources = list(source_order)
    for source in subset["source_catalog"].astype(str).unique():
        if source not in ordered_sources:
            ordered_sources.append(source)

    for source_catalog in ordered_sources:
        source_subset = subset.loc[subset["source_catalog"].astype(str) == source_catalog]
        if source_subset.empty:
            continue
        style = style_for_source(source_catalog)
        colors = point_color_for_rw(source_subset)
        ax.scatter(
            source_subset["x_pc"],
            source_subset["y_pc"],
            s=44 if style["marker"] == "*" else 22,
            marker=style["marker"],
            color=colors,
            edgecolors="white",
            linewidths=0.32,
            alpha=0.96,
            zorder=point_zorder,
        )
        ax.quiver(
            source_subset["x_pc"],
            source_subset["y_pc"],
            source_subset["dvx_lsr_kms"],
            source_subset["dvy_lsr_kms"],
            angles="xy",
            scale_units="xy",
            scale=VELOCITY_SCALE,
            color=colors,
            alpha=0.92,
            width=arrow_width,
            headwidth=4.2,
            headlength=5.0,
            headaxislength=4.4,
            zorder=arrow_zorder,
        )


def draw_background_members(ax: plt.Axes, subset: pd.DataFrame, source_order: tuple[str, ...]) -> None:
    ordered_sources = list(source_order)
    for source in subset["source_catalog"].astype(str).unique():
        if source not in ordered_sources:
            ordered_sources.append(source)

    for source_catalog in ordered_sources:
        source_subset = subset.loc[subset["source_catalog"].astype(str) == source_catalog]
        if source_subset.empty:
            continue
        style = style_for_source(source_catalog)
        ax.scatter(
            source_subset["x_pc"],
            source_subset["y_pc"],
            s=11 if style["marker"] == "*" else 6,
            marker=style["marker"],
            color="0.70",
            edgecolors="none",
            alpha=0.22,
            zorder=3,
        )


def format_source_counts(stat: pd.Series, source_order: tuple[str, ...]) -> str:
    parts = []
    for source in source_order:
        style = style_for_source(source)
        parts.append(f"{style['short']}={int(stat.get(f'n_{safe_name(source)}', 0))}")
    return ", ".join(parts)


def build_legend(catalog_df: pd.DataFrame, source_order: tuple[str, ...]) -> list[Line2D]:
    sources = list(source_order)
    for source in catalog_df["source_catalog"].astype(str).unique():
        if source not in sources:
            sources.append(source)

    handles: list[Line2D] = []
    for source in sources:
        if not (catalog_df["source_catalog"].astype(str) == source).any():
            continue
        style = style_for_source(source)
        handles.append(
            Line2D(
                [0],
                [0],
                marker=style["marker"],
                color="w",
                markerfacecolor=CATALOG_BLUE,
                markersize=10 if style["marker"] == "*" else 7,
                label=style["label"],
            )
        )
    handles.extend(
        [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=RW_HIGHLIGHT_RED, markersize=7, linestyle="None", label="RW band"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="0.70", alpha=0.6, markersize=6, linestyle="None", label="Background"),
            Line2D([0], [0], color="black", lw=1.0, label="Open superbubbles"),
            Line2D([0], [0], marker="*", color="black", markerfacecolor="gold", markersize=10, linestyle="None", label="Sun"),
        ]
    )
    return handles


def make_xy_residual_plot(
    catalog_df: pd.DataFrame,
    ellipse_df: pd.DataFrame,
    family_labels: list[str],
    detail_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    config: XYDatasetConfig,
) -> None:
    n_families = len(family_labels)
    ncols = 3
    nrows = int(np.ceil(n_families / ncols))
    min_print_fontsize = PRINT_BODY_FONTSIZE
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=scaled_figsize(A4_TEXT_WIDTH_IN, 0.65 + 2.65 * nrows),
        sharex=True,
        sharey=True,
    )
    axes = np.asarray(axes).reshape(-1)

    for ax, family in zip(axes, family_labels):
        subset = detail_df.loc[detail_df["family"] == family].copy()
        non_family = catalog_df.loc[catalog_df["family"] != family].copy()
        stat = summary_df.loc[summary_df["family"] == family].iloc[0]

        draw_superbubble_ellipses(ax, ellipse_df)
        draw_background_members(ax, non_family, config.source_order)
        draw_family_members(ax, subset, config.source_order)
        add_solar_marker(ax)
        add_velocity_scale(ax)
        add_family_median_velocity_arrow(ax, stat)

        ax.set_xlim(*XY_XLIM_PC)
        ax.set_ylim(*XY_YLIM_PC)
        ax.set_xticks([-1000.0, 0.0, 1000.0])
        ax.set_yticks([-1000.0, 0.0, 1000.0])
        ax.set_xticklabels(["-1", "0", "1"])
        ax.set_yticklabels(["-1", "0", "1"])
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, linestyle="--", alpha=0.26)
        ax.set_facecolor("whitesmoke")
        ax.tick_params(axis="both", which="major", labelsize=min_print_fontsize, length=4.0, width=0.8)
        ax.set_title(f"{family} (N={int(stat['n_members'])})", fontsize=min_print_fontsize, pad=4.0)

    legend_handles = build_legend(catalog_df, config.source_order)
    empty_axes = axes[n_families:]
    for ax in empty_axes:
        ax.axis("off")

    has_legend_panel = len(empty_axes) > 0
    fig.supxlabel("Heliocentric X (kpc)", fontsize=min_print_fontsize, y=0.055 if has_legend_panel else 0.155)
    fig.supylabel("Heliocentric Y (kpc)", fontsize=min_print_fontsize, x=0.055)
    if has_legend_panel:
        fig.subplots_adjust(left=0.100, right=0.99, top=0.955, bottom=0.115, wspace=0.03, hspace=0.26)
        legend = empty_axes[0].legend(
            handles=legend_handles,
            loc="center left",
            bbox_to_anchor=(-0.03, 0.5),
            ncol=1,
            fontsize=min_print_fontsize,
            frameon=True,
            handlelength=1.5,
            labelspacing=0.65,
            borderpad=0.65,
        )
        legend.set_zorder(30)
    else:
        fig.subplots_adjust(left=0.100, right=0.99, top=0.955, bottom=0.245, wspace=0.03, hspace=0.24)
        fig.legend(
            handles=legend_handles,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.010),
            ncol=3,
            fontsize=min_print_fontsize,
            frameon=True,
        )
    config.output_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(config.output_fig, dpi=235, bbox_inches="tight")
    plt.close(fig)


def write_xy_outputs(
    summary_df: pd.DataFrame,
    detail_df: pd.DataFrame,
    ellipse_df: pd.DataFrame,
    config: XYDatasetConfig,
    output_ellipse: Path | None,
) -> None:
    config.output_summary.parent.mkdir(parents=True, exist_ok=True)
    config.output_detail.parent.mkdir(parents=True, exist_ok=True)
    config.output_fig.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(config.output_summary, index=False, encoding="utf-8-sig")
    detail_df[[col for col in DETAIL_COLUMNS if col in detail_df.columns]].to_csv(config.output_detail, index=False, encoding="utf-8-sig")
    if output_ellipse is not None:
        output_ellipse.parent.mkdir(parents=True, exist_ok=True)
        ellipse_df.to_csv(output_ellipse, index=False, encoding="utf-8-sig")


def render_xy_dataset(config: XYDatasetConfig, ellipse_csv: Path, output_ellipse: Path | None = None) -> tuple[list[str], pd.DataFrame]:
    catalog_df = load_catalog(config.source_csv)
    ellipse_df = load_compromise_ellipses(ellipse_csv)
    family_labels, detail_df, summary_df = prepare_residuals(catalog_df, config.source_order)
    write_xy_outputs(summary_df, detail_df, ellipse_df, config, output_ellipse)
    make_xy_residual_plot(catalog_df, ellipse_df, family_labels, detail_df, summary_df, config)
    return family_labels, summary_df


@dataclass(frozen=True)
class PlotDatasetConfig:
    key: str
    output_prefix: str
    title: str
    labels_file: Path
    orbit_file: Path
    properties_file: Path
    curve_file: Path
    source_order: tuple[str, ...]
    max_lookback_myr: float
    slice_xy_limit_pc: float

    @property
    def xy_config(self) -> XYDatasetConfig:
        return XYDatasetConfig(
            key=self.key,
            title='Residual velocity field',
            source_csv=self.labels_file,
            output_fig=FINAL_FIG_DIR / f"{self.output_prefix}_XY.png",
            output_summary=OUT_DIR / f"{self.output_prefix}_xy_statistics.csv",
            output_detail=OUT_DIR / f"{self.output_prefix}_xy_details.csv",
            source_order=self.source_order,
        )


DATASETS = {
    "first": PlotDatasetConfig(
        key="first",
        output_prefix="G1",
        title='Residual velocity field',
        labels_file=G1_DIR / "G1_family_labels.csv",
        orbit_file=G1_DIR / "G1_orbits_55Myr.npz",
        properties_file=G1_DIR / "G1_family_properties.csv",
        curve_file=G1_DIR / "G1_family_size_evolution.csv",
        source_order=("hunt_stratified", "konietzka2023"),
        max_lookback_myr=55.0,
        slice_xy_limit_pc=2000.0,
    ),
    "second": PlotDatasetConfig(
        key="second",
        output_prefix="G2",
        title='Residual velocity field',
        labels_file=G2_DIR / "G2_family_labels.csv",
        orbit_file=G2_DIR / "G2_orbits_50Myr.npz",
        properties_file=G2_DIR / "G2_family_properties.csv",
        curve_file=G2_DIR / "G2_family_size_evolution.csv",
        source_order=("hunt_stratified", "konietzka2023"),
        max_lookback_myr=50.0,
        slice_xy_limit_pc=2000.0,
    ),
}


SELECTED_XY_FAMILIES = {
    "first": ("group_01", "group_02", "group_04"),
    "second": ("group_02", "group_06", "group_09"),
}
SELECTED_XY_ORDER = (
    "G1_group_01",
    "G1_group_02",
    "G1_group_04",
    "G2_group_02",
    "G2_group_06",
    "G2_group_09",
)
SELECTED_LBW_FAMILIES = SELECTED_XY_FAMILIES
SELECTED_LBW_ORDER = SELECTED_XY_ORDER
LB_PANEL_OVERLAYS = {
    "G1_group_01": (("hbw", "hbw35"),),
    "G1_group_02": (("hbw", "hbw35"), ("sb", 31)),
    "G1_group_04": (("sb", 11),),
    "G2_group_02": (("sb", 16),),
    "G2_group_06": (("sb", 25),),
    "G2_group_09": (("sb", 15),),
}
LBW_PANEL_LIMITS = {
    "G1_group_02": {"xlim": (180.0, 280.0), "ylim": (-30.0, 10.0)},
    "G1_group_04": {"ylim": (-5.0, 15.0)},
    "G2_group_02": {"xlim": (20.0, 120.0)},
    "G2_group_09": {"xlim": (20.0, 140.0)},
}


def require_file(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError('Missing input file.')
    return path


def selected_family_label(config: PlotDatasetConfig, family: str) -> str:
    prefix = "G1" if config.output_prefix == "G1" else "G2"
    return f"{prefix}_{family}"


def build_selected_xy_catalog() -> pd.DataFrame:
    frames = []
    for key, families in SELECTED_XY_FAMILIES.items():
        config = DATASETS[key]
        catalog = load_catalog(require_file(config.labels_file))
        original_family = catalog["family"].astype(str)
        catalog["source_dataset"] = config.output_prefix
        catalog["original_family"] = original_family
        catalog["family"] = "none"
        for family in families:
            catalog.loc[original_family == family, "family"] = selected_family_label(config, family)
        frames.append(catalog)
    combined = pd.concat(frames, ignore_index=True)
    output_file = OUT_DIR / "G1G2_selected_family_xy_residual_input.csv"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_file, index=False, encoding="utf-8-sig")
    return combined


def selected_family_sort_key(value: str) -> int:
    try:
        return SELECTED_XY_ORDER.index(str(value))
    except ValueError:
        return 999


def prepare_selected_xy_residuals(catalog_df: pd.DataFrame) -> tuple[list[str], pd.DataFrame, pd.DataFrame]:
    family_labels = sorted([family for family in catalog_df["family"].unique() if family != "none"], key=selected_family_sort_key)
    selected_catalog = catalog_df.loc[catalog_df["family"].isin(family_labels)].copy()
    _, detail_df, summary_df = prepare_residuals(selected_catalog, SOURCE_ORDER)
    detail_df["family"] = pd.Categorical(detail_df["family"], categories=family_labels, ordered=True)
    detail_df = detail_df.sort_values(["family", "source_catalog", "display_name"]).copy()
    detail_df["family"] = detail_df["family"].astype(str)
    summary_df["family"] = pd.Categorical(summary_df["family"], categories=family_labels, ordered=True)
    summary_df = summary_df.sort_values("family").copy()
    summary_df["family"] = summary_df["family"].astype(str)
    return family_labels, detail_df, summary_df


def selected_family_xy_limits(family: str, subset: pd.DataFrame) -> tuple[tuple[float, float], tuple[float, float]]:
    if family == "G1_group_02":
        center_x = -250.0
        center_y = -500.0
        half_width = 1250.0
    else:
        center_x = float(subset["x_pc"].astype(float).median())
        center_y = float(subset["y_pc"].astype(float).median())
        half_width = SELECTED_XY_HALF_WIDTH_PC
    return (
        (center_x - half_width, center_x + half_width),
        (center_y - half_width, center_y + half_width),
    )


def selected_ellipse_may_intersect_view(
    row: pd.Series,
    a_column: str,
    b_column: str,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> bool:
    x0 = float(row["x_kpc"]) * 1000.0
    y0 = float(row["y_kpc"]) * 1000.0
    radius = float(max(row[a_column], row[b_column])) * 1000.0
    return (
        x0 + radius >= xlim[0]
        and x0 - radius <= xlim[1]
        and y0 + radius >= ylim[0]
        and y0 - radius <= ylim[1]
    )


def draw_selected_ellipse(ax: plt.Axes, row: pd.Series, a_column: str, b_column: str, linestyle: str, alpha: float, zorder: int) -> None:
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


def draw_selected_superbubble_ellipses(
    ax: plt.Axes,
    ellipse_df: pd.DataFrame,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> None:
    for _, row in ellipse_df.iterrows():
        bubble_id = int(row["index"])
        if bubble_id in DASHED_MAX_OPEN_IDS and selected_ellipse_may_intersect_view(row, "a_maxopen_kpc", "b_maxopen_kpc", xlim, ylim):
            draw_selected_ellipse(ax, row, "a_maxopen_kpc", "b_maxopen_kpc", (0, (4, 3)), 0.45, 4)
        if selected_ellipse_may_intersect_view(row, "a_compromise_kpc", "b_compromise_kpc", xlim, ylim):
            draw_selected_ellipse(ax, row, "a_compromise_kpc", "b_compromise_kpc", "-", 0.68, 3)
            x0 = float(row["x_kpc"]) * 1000.0
            y0 = float(row["y_kpc"]) * 1000.0
            if xlim[0] <= x0 <= xlim[1] and ylim[0] <= y0 <= ylim[1]:
                ax.text(
                    x0 + SB_LABEL_OFFSETS_PC.get(bubble_id, (0.0, 0.0))[0],
                    y0 + SB_LABEL_OFFSETS_PC.get(bubble_id, (0.0, 0.0))[1],
                    str(bubble_id),
                    color="black",
                    fontsize=6.4,
                    ha="center",
                    va="center",
                    fontweight="bold",
                    clip_on=True,
                    path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
                    zorder=4.2,
                )


def load_hbw35_bubble() -> pd.Series | None:
    if not HBW_BUBBLE_CSV.exists():
        return None
    bubbles = pd.read_csv(HBW_BUBBLE_CSV)
    required = ["ID", "l", "b", "best_distance", "physical_size"]
    if any(column not in bubbles.columns for column in required):
        return None
    bubbles["ID"] = bubbles["ID"].astype(str)
    row = bubbles.loc[bubbles["ID"].str.lower() == HBW_BUBBLE_ID.lower()].copy()
    if row.empty:
        return None
    for column in ["l", "b", "best_distance", "physical_size"]:
        row[column] = pd.to_numeric(row[column], errors="coerce")
    row = row.dropna(subset=["l", "b", "best_distance", "physical_size"])
    if row.empty:
        return None
    item = row.iloc[0].copy()
    l_rad = np.deg2rad(float(item["l"]))
    b_rad = np.deg2rad(float(item["b"]))
    dist_pc = float(item["best_distance"]) * 1000.0
    item["x_pc"] = dist_pc * np.cos(b_rad) * np.cos(l_rad)
    item["y_pc"] = dist_pc * np.cos(b_rad) * np.sin(l_rad)
    item["radius_pc"] = float(item["physical_size"]) / 2.0
    return item


def draw_hbw35_bubble(ax: plt.Axes, bubble: pd.Series | None) -> None:
    if bubble is None:
        return
    ax.add_patch(
        Circle(
            (float(bubble["x_pc"]), float(bubble["y_pc"])),
            float(bubble["radius_pc"]),
            facecolor="none",
            edgecolor="#00a6d6",
            linewidth=1.45,
            alpha=0.92,
            zorder=2,
        )
    )


def draw_selected_background(ax: plt.Axes, non_family: pd.DataFrame) -> None:
    sources = list(SOURCE_ORDER)
    for source in non_family["source_catalog"].astype(str).unique():
        if source not in sources:
            sources.append(source)
    for source in sources:
        source_subset = non_family.loc[non_family["source_catalog"].astype(str) == source]
        if source_subset.empty:
            continue
        style = style_for_source(source)
        ax.scatter(
            source_subset["x_pc"],
            source_subset["y_pc"],
            s=13 if style["marker"] == "*" else 7,
            marker=style["marker"],
            color="0.70",
            edgecolors="none",
            alpha=0.18,
            zorder=1,
        )


def add_selected_velocity_scale(ax: plt.Axes, xlim: tuple[float, float], ylim: tuple[float, float]) -> None:
    anchor_x = xlim[0] + 0.08 * (xlim[1] - xlim[0])
    anchor_y = ylim[1] - 0.12 * (ylim[1] - ylim[0])
    arrow_dx = REFERENCE_VELOCITY_KMS / VELOCITY_SCALE
    ax.add_patch(
        Rectangle(
            (anchor_x - 35.0, anchor_y - 70.0),
            arrow_dx + 96.0,
            215.0,
            facecolor="0.96",
            edgecolor="0.78",
            linewidth=0.55,
            alpha=0.9,
            zorder=12.5,
        )
    )
    ax.annotate(
        "",
        xy=(anchor_x + arrow_dx, anchor_y),
        xytext=(anchor_x, anchor_y),
        xycoords="data",
        textcoords="data",
        arrowprops={
            "arrowstyle": "-|>",
            "color": "black",
            "lw": 1.7,
            "mutation_scale": 10.5,
            "shrinkA": 0.0,
            "shrinkB": 0.0,
        },
        zorder=13,
    )
    ax.text(
        anchor_x + 0.5 * arrow_dx,
        anchor_y + 62.0,
        f"{REFERENCE_VELOCITY_KMS:g} km/s",
        fontsize=SELECTED_SPEED_LABEL_FONTSIZE,
        color="black",
        ha="center",
        va="bottom",
        zorder=13,
    )


def add_selected_family_median_velocity_arrow(
    ax: plt.Axes,
    stat: pd.Series,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> None:
    vx = float(stat["median_vx_lsr_kms"])
    vy = float(stat["median_vy_lsr_kms"])
    dx_pc = vx / VELOCITY_SCALE
    dy_pc = vy / VELOCITY_SCALE
    x_span = xlim[1] - xlim[0]
    y_span = ylim[1] - ylim[0]
    box_left = xlim[1] - 0.36 * x_span
    box_right = xlim[1] - 0.06 * x_span
    box_bottom = ylim[1] - 0.26 * y_span
    box_top = ylim[1] - 0.08 * y_span
    tail_x = np.clip(box_right - max(dx_pc, 0.0), box_left - min(dx_pc, 0.0), box_right - max(dx_pc, 0.0))
    tail_y = np.clip(box_top - max(dy_pc, 0.0), box_bottom - min(dy_pc, 0.0), box_top - max(dy_pc, 0.0))
    ax.quiver(
        tail_x,
        tail_y,
        vx,
        vy,
        angles="xy",
        scale_units="xy",
        scale=VELOCITY_SCALE,
        color="#c1121f",
        width=0.0062,
        headwidth=4.8,
        headlength=5.8,
        headaxislength=5.0,
        zorder=14,
    )


def selected_panel_legend(subset: pd.DataFrame, non_family: pd.DataFrame, stat: pd.Series, bubble: pd.Series | None) -> list[Line2D]:
    handles: list[Line2D] = []
    for source in SOURCE_ORDER:
        style = style_for_source(source)
        if (subset["source_catalog"].astype(str) == source).any():
            handles.append(
                Line2D(
                    [0],
                    [0],
                    marker=style["marker"],
                    color="w",
                    markerfacecolor=CATALOG_BLUE,
                    markeredgecolor="white",
                    markersize=8.5 if style["marker"] == "*" else 6.5,
                    linestyle="None",
                    label=style["label"],
                )
            )
    handles.extend(
        [
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=RW_HIGHLIGHT_RED,
                markersize=6.5,
                linestyle="None",
                label=RW_CLUSTER_LABEL,
            ),
            Line2D([0], [0], color=RW_CENTERLINE_COLOR, lw=1.35, linestyle="-", label=RADCLIFFE_WAVE_LABEL),
            Line2D([0], [0], color="black", lw=0.9, label="Open superbubbles"),
        ]
    )
    if bubble is not None:
        handles.append(Line2D([0], [0], color="#00a6d6", lw=1.3, label=ORION_ERIDANUS_LABEL))
    return handles


def selected_panel_title(family: str, stat: pd.Series) -> str:
    return f"{family.replace('_group_', ' group ')}  |  N={int(stat['n_members'])}"


def draw_radcliffe_wave_xy_centerline(ax: plt.Axes, xlim: tuple[float, float], ylim: tuple[float, float]) -> None:
    theta = np.deg2rad(RW_ANGLE_DEG)
    x_values = np.linspace(xlim[0], xlim[1], 600)
    y_values = np.tan(theta) * x_values + 1000.0 * RW_INTERCEPT_KPC
    visible = (y_values >= ylim[0]) & (y_values <= ylim[1])
    if not np.any(visible):
        return
    ax.plot(
        x_values[visible],
        y_values[visible],
        color=RW_CENTERLINE_COLOR,
        lw=1.35,
        linestyle="-",
        alpha=0.95,
        zorder=5,
    )


def draw_selected_xy_family_panel(
    ax: plt.Axes,
    family: str,
    catalog_df: pd.DataFrame,
    detail_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    ellipse_df: pd.DataFrame,
    bubble: pd.Series | None,
    *,
    show_legend: bool = True,
    show_orion_eridanus: bool = True,
) -> None:
    subset = detail_df.loc[detail_df["family"] == family].copy()
    non_family = catalog_df.loc[catalog_df["family"] != family].copy()
    stat = summary_df.loc[summary_df["family"] == family].iloc[0]
    xlim, ylim = selected_family_xy_limits(family, subset)

    panel_bubble = bubble if show_orion_eridanus else None
    draw_selected_background(ax, non_family)
    draw_hbw35_bubble(ax, panel_bubble)
    draw_selected_superbubble_ellipses(ax, ellipse_df, xlim, ylim)
    draw_radcliffe_wave_xy_centerline(ax, xlim, ylim)
    draw_family_members(ax, subset, SOURCE_ORDER, point_zorder=11, arrow_zorder=13, arrow_width=0.0046)
    add_solar_marker(ax)
    add_selected_velocity_scale(ax, xlim, ylim)
    add_selected_family_median_velocity_arrow(ax, stat, xlim, ylim)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linestyle=":", alpha=0.22, linewidth=0.55)
    ax.set_facecolor("white")
    ax.set_title(selected_panel_title(family, stat), fontsize=10.2, pad=5)
    ax.tick_params(direction="in", top=True, right=True, length=3.5, width=0.75, labelsize=8.5)
    if show_legend:
        legend = ax.legend(
            handles=selected_panel_legend(subset, non_family, stat, panel_bubble),
            loc="lower right",
            fontsize=SELECTED_LEGEND_FONTSIZE,
            frameon=True,
            framealpha=0.88,
            borderpad=0.32,
            handlelength=1.55,
            handletextpad=0.45,
            labelspacing=0.22,
        )
        legend.set_zorder(30)


def build_selected_lbw_catalog() -> pd.DataFrame:
    frames = []
    for key, families in SELECTED_LBW_FAMILIES.items():
        config = DATASETS[key]
        catalog = load_catalog(require_file(config.labels_file))
        original_family = catalog["family"].astype(str)
        catalog["source_dataset"] = config.output_prefix
        catalog["original_family"] = original_family
        catalog["plot_family"] = "none"
        for family in families:
            catalog.loc[original_family == family, "plot_family"] = selected_family_label(config, family)
        frames.append(catalog)
    combined = pd.concat(frames, ignore_index=True)
    for column in ["x_pc", "y_pc", "z_pc", "vz_lsr_kms"]:
        combined[column] = pd.to_numeric(combined[column], errors="coerce")
    dist = np.sqrt(combined["x_pc"] ** 2 + combined["y_pc"] ** 2 + combined["z_pc"] ** 2)
    combined["distance_pc"] = dist
    combined["l_deg"] = np.degrees(np.arctan2(combined["y_pc"], combined["x_pc"])) % 360.0
    combined["b_deg"] = np.degrees(np.arcsin(np.divide(combined["z_pc"], dist, out=np.zeros(len(combined)), where=dist > 0.0)))
    combined = combined.dropna(subset=["l_deg", "b_deg", "vz_lsr_kms"]).copy()
    combined["plot_order"] = combined["plot_family"].map({family: index for index, family in enumerate(SELECTED_LBW_ORDER)})
    combined = combined.sort_values(["plot_order", "l_deg", "b_deg"], kind="stable").reset_index(drop=True)
    output_file = OUT_DIR / "G1G2_selected_family_lb_w_velocity.csv"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_file, index=False, encoding="utf-8-sig")
    return combined


def load_hbw35_lb_bubble() -> pd.Series | None:
    if not HBW_BUBBLE_CSV.exists():
        return None
    bubbles = pd.read_csv(HBW_BUBBLE_CSV)
    required = ["ID", "l", "b", "D"]
    if any(column not in bubbles.columns for column in required):
        return None
    bubbles["ID"] = bubbles["ID"].astype(str)
    row = bubbles.loc[bubbles["ID"].str.lower() == HBW_BUBBLE_ID.lower()].copy()
    if row.empty:
        return None
    for column in ["l", "b", "D"]:
        row[column] = pd.to_numeric(row[column], errors="coerce")
    row = row.dropna(subset=["l", "b", "D"])
    if row.empty:
        return None
    return row.iloc[0]


def draw_hbw35_lb_circle(
    ax: plt.Axes,
    label_position: tuple[float, float] | None = None,
    label_fontsize: float = 8.2,
) -> tuple[Line2D | None, tuple[float, float, float, float] | None]:
    row = load_hbw35_lb_bubble()
    if row is None:
        return None, None
    l0 = float(row["l"]) % 360.0
    b0 = float(row["b"])
    diameter_b = float(row["D"])
    diameter_l = diameter_b / max(np.cos(np.deg2rad(b0)), 0.1)
    ax.add_patch(
        Ellipse(
            (l0, b0),
            width=diameter_l,
            height=diameter_b,
            angle=0.0,
            facecolor="none",
            edgecolor="#00a6d6",
            linewidth=1.45,
            alpha=0.96,
            zorder=2,
        )
    )
    label_l, label_b = label_position if label_position is not None else (l0, b0)
    ax.text(
        label_l,
        label_b,
        "Orion-Eridanus\nsuperbubble",
        color="#0077aa",
        fontsize=label_fontsize,
        fontweight="bold",
        ha="center",
        va="center",
        bbox={"facecolor": "white", "edgecolor": "none", "boxstyle": "round,pad=0.12", "alpha": 0.72},
        zorder=7,
    )
    handle = Line2D([0], [0], color="#00a6d6", lw=1.45, label=ORION_ERIDANUS_LABEL)
    bounds = (l0 - 0.5 * diameter_l, l0 + 0.5 * diameter_l, b0 - 0.5 * diameter_b, b0 + 0.5 * diameter_b)
    return handle, bounds


def load_superbubble_params(sb_id: int) -> pd.Series:
    bubbles = pd.read_csv(require_file(ELLIPSE_CSV))
    row = bubbles.loc[pd.to_numeric(bubbles["id"], errors="coerce") == int(sb_id)].copy()
    if row.empty:
        raise ValueError(f"{ELLIPSE_CSV} does not contain SB {sb_id}")
    return row.iloc[0]


def family_center_distance_kpc(subset: pd.DataFrame) -> float:
    center_pc = subset[["x_pc", "y_pc", "z_pc"]].astype(float).median(axis=0).to_numpy()
    return float(np.linalg.norm(center_pc) / 1000.0)


def superbubble_center_distance_kpc(row: pd.Series) -> float:
    center = np.asarray([float(row["center_x_kpc"]), float(row["center_y_kpc"]), float(row["center_z_kpc"])])
    return float(np.linalg.norm(center))


def ellipsoid_half_value_on_distance(
    row: pd.Series,
    distance_kpc: float,
    l_grid: np.ndarray,
    b_grid: np.ndarray,
    buffer_kpc: float = 0.0,
) -> np.ma.MaskedArray:
    l_rad = np.deg2rad(l_grid)
    b_rad = np.deg2rad(b_grid)
    x = distance_kpc * np.cos(b_rad) * np.cos(l_rad)
    y = distance_kpc * np.cos(b_rad) * np.sin(l_rad)
    z = distance_kpc * np.sin(b_rad)

    cx = float(row["center_x_kpc"])
    cy = float(row["center_y_kpc"])
    cz = float(row["center_z_kpc"])
    a = float(row["a_radius_kpc"])
    b = float(row["b_radius_kpc"])
    c = float(row["c_radius_kpc"])
    theta = np.deg2rad(float(row["angle_deg"]))

    dx = x - cx
    dy = y - cy
    dz = z - cz
    major = np.cos(theta) * dx + np.sin(theta) * dy
    minor = -np.sin(theta) * dx + np.cos(theta) * dy
    value = (major / a) ** 2 + (minor / b) ** 2 + (dz / c) ** 2

    apex_z = cz - c if int(row["mark"]) == 1 else cz + c
    z_low, z_high = sorted((cz, apex_z))
    half_mask = (z >= z_low - buffer_kpc) & (z <= z_high + buffer_kpc)
    return np.ma.masked_where(~half_mask, value)


def points_inside_true_half(row: pd.Series, distance_kpc: float, points: np.ndarray) -> np.ndarray:
    b_rad = np.deg2rad(points[:, 1])
    z = distance_kpc * np.sin(b_rad)
    cz = float(row["center_z_kpc"])
    c = float(row["c_radius_kpc"])
    apex_z = cz - c if int(row["mark"]) == 1 else cz + c
    z_low, z_high = sorted((cz, apex_z))
    return (z >= z_low) & (z <= z_high)


def contiguous_true_half_segments(row: pd.Series, distance_kpc: float, points: np.ndarray) -> list[np.ndarray]:
    keep = points_inside_true_half(row, distance_kpc, points)
    segments: list[np.ndarray] = []
    start: int | None = None
    for index, is_kept in enumerate(keep):
        if is_kept and start is None:
            start = index
        elif not is_kept and start is not None:
            if index - start >= 2:
                segments.append(points[start:index])
            start = None
    if start is not None and len(points) - start >= 2:
        segments.append(points[start:])
    return segments


def label_superbubble_id(ax: plt.Axes, sb_id: int, segments: list[np.ndarray]) -> None:
    if not segments:
        return
    points = np.vstack(segments)
    point = np.array(
        [
            0.5 * (np.nanmin(points[:, 0]) + np.nanmax(points[:, 0])),
            0.5 * (np.nanmin(points[:, 1]) + np.nanmax(points[:, 1])),
        ]
    )
    ax.text(
        point[0],
        point[1],
        f"SB {sb_id}",
        color="#11823b",
        fontsize=8.2,
        fontweight="bold",
        ha="center",
        va="center",
        bbox={"facecolor": "white", "edgecolor": "none", "boxstyle": "round,pad=0.12", "alpha": 0.75},
        zorder=7,
    )


def draw_superbubble_lb_section(
    ax: plt.Axes,
    row: pd.Series,
    sb_id: int,
    distance_kpc: float,
    linestyle: str,
    label: str,
    *,
    linewidth: float = 1.25,
    alpha: float = 0.88,
    annotate_label: bool = True,
) -> tuple[Line2D, tuple[float, float, float, float] | None]:
    l_values = np.linspace(LB_LIM_L_DEG[0], LB_LIM_L_DEG[1], 1601)
    b_values = np.linspace(LB_LIM_B_DEG[0], LB_LIM_B_DEG[1], 701)
    l_grid, b_grid = np.meshgrid(l_values, b_values)
    values = ellipsoid_half_value_on_distance(row, distance_kpc, l_grid, b_grid, buffer_kpc=HALF_ELLIPSOID_BUFFER_KPC)
    contour = ax.contour(l_grid, b_grid, values, levels=[1.0], colors="#111111", linewidths=0.0, alpha=0.0)
    if hasattr(contour, "collections"):
        for collection in contour.collections:
            collection.remove()

    handle = Line2D([0], [0], color="#111111", lw=linewidth, linestyle=linestyle, label=label)
    raw_segments = [segment for segment in contour.allsegs[0] if len(segment) > 0]
    plotted_segments: list[np.ndarray] = []
    for segment in raw_segments:
        for clipped in contiguous_true_half_segments(row, distance_kpc, segment):
            ax.plot(clipped[:, 0], clipped[:, 1], color="#111111", lw=linewidth, linestyle=linestyle, alpha=alpha, zorder=2)
            plotted_segments.append(clipped)
    if not plotted_segments:
        return handle, None
    if annotate_label:
        label_superbubble_id(ax, sb_id, plotted_segments)
    points = np.vstack(plotted_segments)
    bounds = (float(np.nanmin(points[:, 0])), float(np.nanmax(points[:, 0])), float(np.nanmin(points[:, 1])), float(np.nanmax(points[:, 1])))
    return handle, bounds


def draw_lbw_panel_overlays(ax: plt.Axes, subset: pd.DataFrame, family: str) -> tuple[list[Line2D], list[tuple[float, float, float, float]]]:
    handles: list[Line2D] = []
    bounds: list[tuple[float, float, float, float]] = []
    for kind, value in LB_PANEL_OVERLAYS.get(family, ()):
        if kind == "hbw":
            label_position = (165.0, -15.0) if family == "G1_group_01" else None
            label_fontsize = 9.2 if family == "G1_group_01" else 8.2
            handle, overlay_bounds = draw_hbw35_lb_circle(ax, label_position=label_position, label_fontsize=label_fontsize)
        elif kind == "sb":
            sb_id = int(value)
            row = load_superbubble_params(sb_id)
            handle, overlay_bounds = draw_superbubble_lb_section(
                ax=ax,
                row=row,
                sb_id=sb_id,
                distance_kpc=superbubble_center_distance_kpc(row),
                linestyle="-",
                label=f"SB{sb_id}",
            )
        else:
            continue
        if handle is not None:
            handles.append(handle)
        if overlay_bounds is not None:
            bounds.append(overlay_bounds)
    return handles, bounds


def set_lbw_panel_limits(ax: plt.Axes, subset: pd.DataFrame, overlay_bounds: list[tuple[float, float, float, float]], family: str) -> None:
    l_min = float(subset["l_deg"].min())
    l_max = float(subset["l_deg"].max())
    b_min = float(subset["b_deg"].min())
    b_max = float(subset["b_deg"].max())
    arrow_tips = subset["b_deg"].to_numpy(float) + subset["vz_lsr_kms"].to_numpy(float) / W_VELOCITY_SCALE
    b_min = min(b_min, float(np.nanmin(arrow_tips)))
    b_max = max(b_max, float(np.nanmax(arrow_tips)))
    for left, right, bottom, top in overlay_bounds:
        l_min = min(l_min, left)
        l_max = max(l_max, right)
        b_min = min(b_min, bottom)
        b_max = max(b_max, top)

    fixed_limits = LBW_PANEL_LIMITS.get(family, {})
    fixed_xlim = fixed_limits.get("xlim")
    if fixed_xlim is not None:
        x0, x1 = fixed_xlim
    elif l_max - l_min > 240.0:
        x0, x1 = LB_LIM_L_DEG
    else:
        pad_l = max(8.0, 0.18 * (l_max - l_min))
        x0 = max(LB_LIM_L_DEG[0], l_min - pad_l)
        x1 = min(LB_LIM_L_DEG[1], l_max + pad_l)
    pad_b = max(5.0, 0.18 * (b_max - b_min))
    y0 = max(LB_LIM_B_DEG[0], b_min - pad_b)
    y1 = min(LB_LIM_B_DEG[1], b_max + pad_b)
    if y1 - y0 < 24.0:
        center = 0.5 * (y0 + y1)
        y0 = max(LB_LIM_B_DEG[0], center - 12.0)
        y1 = min(LB_LIM_B_DEG[1], center + 12.0)
    fixed_ylim = fixed_limits.get("ylim")
    if fixed_ylim is not None:
        y0, y1 = fixed_ylim
    else:
        y0 = min(y0 + LBW_VERTICAL_TRIM_DEG, b_min - 2.0)
        y1 = max(y1 - LBW_VERTICAL_TRIM_DEG, b_max + 2.0)
    ax.set_xlim(x1, x0)
    ax.set_ylim(y0, y1)


def draw_lbw_reference_arrow(ax: plt.Axes) -> None:
    x_ref = 0.09
    y_tip = 0.84
    ax.add_patch(
        Rectangle(
            (0.055, 0.62),
            0.25,
            0.28,
            transform=ax.transAxes,
            facecolor="white",
            edgecolor="0.78",
            linewidth=0.55,
            alpha=0.82,
            zorder=7.5,
        )
    )
    ax.annotate(
        "",
        xy=(x_ref, y_tip),
        xytext=(0.0, -REFERENCE_W_ARROW_POINTS),
        xycoords="axes fraction",
        textcoords="offset points",
        arrowprops={
            "arrowstyle": "-|>",
            "color": "black",
            "lw": 1.45,
            "mutation_scale": 10.5,
            "shrinkA": 0.0,
            "shrinkB": 0.0,
        },
        zorder=8,
    )
    ax.annotate(
        f"{REFERENCE_W_KMS:g} km/s",
        xy=(x_ref, y_tip),
        xytext=(7.0, 0.0),
        xycoords="axes fraction",
        textcoords="offset points",
        fontsize=SELECTED_SPEED_LABEL_FONTSIZE,
        va="center",
        ha="left",
        zorder=8,
    )


def lbw_panel_handles(subset: pd.DataFrame, overlay_handles: list[Line2D]) -> list[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=RW_HIGHLIGHT_RED,
            markersize=6.5,
            linestyle="None",
            label=RW_CLUSTER_LABEL,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=CATALOG_BLUE,
            markersize=6.5,
            linestyle="None",
            label=NON_RW_CLUSTER_LABEL,
        ),
    ]


def selected_shared_legend_handles() -> list[Line2D]:
    hunt_style = style_for_source("hunt_stratified")
    konietzka_style = style_for_source("konietzka2023")
    return [
        Line2D(
            [0],
            [0],
            marker=hunt_style["marker"],
            color="w",
            markerfacecolor=CATALOG_BLUE,
            markeredgecolor="white",
            markersize=6.5,
            linestyle="None",
            label=hunt_style["label"],
        ),
        Line2D(
            [0],
            [0],
            marker=konietzka_style["marker"],
            color="w",
            markerfacecolor=CATALOG_BLUE,
            markeredgecolor="white",
            markersize=8.5,
            linestyle="None",
            label=konietzka_style["label"],
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=RW_HIGHLIGHT_RED,
            markersize=6.5,
            linestyle="None",
            label=RW_CLUSTER_LABEL,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=CATALOG_BLUE,
            markersize=6.5,
            linestyle="None",
            label=NON_RW_CLUSTER_LABEL,
        ),
        Line2D([0], [0], color=RW_CENTERLINE_COLOR, lw=1.35, linestyle="-", label=RADCLIFFE_WAVE_LABEL),
        Line2D([0], [0], color="black", lw=0.9, label="Open superbubbles"),
        Line2D([0], [0], color="#00a6d6", lw=1.3, label=ORION_ERIDANUS_LABEL),
    ]


def draw_lbw_family_panel(
    ax: plt.Axes,
    subset: pd.DataFrame,
    family: str,
    *,
    equal_aspect: bool = True,
    aspect_adjustable: str = "box",
    show_title: bool = True,
    box_aspect: float | None = None,
    show_main_legend: bool = True,
) -> None:
    overlay_handles, overlay_bounds = draw_lbw_panel_overlays(ax, subset, family)
    for source in SOURCE_ORDER:
        source_subset = subset.loc[subset["source_catalog"].astype(str) == source]
        if source_subset.empty:
            continue
        style = style_for_source(source)
        colors = point_color_for_rw(source_subset)
        ax.scatter(
            source_subset["l_deg"],
            source_subset["b_deg"],
            s=58 if style["marker"] == "*" else 38,
            marker=style["marker"],
            color=colors,
            edgecolors="white",
            linewidths=0.55,
            alpha=0.94,
            zorder=4,
        )
        ax.quiver(
            source_subset["l_deg"],
            source_subset["b_deg"],
            np.zeros(len(source_subset)),
            source_subset["vz_lsr_kms"],
            angles="xy",
            scale_units="xy",
            scale=W_VELOCITY_SCALE,
            units="inches",
            color=colors,
            width=0.018,
            headwidth=4.2,
            headlength=5.1,
            headaxislength=4.5,
            alpha=0.88,
            zorder=5,
        )

    ax.axhline(0.0, color="0.35", lw=0.7, alpha=0.45)
    set_lbw_panel_limits(ax, subset, overlay_bounds, family)
    if equal_aspect:
        ax.set_aspect("equal", adjustable=aspect_adjustable)
    else:
        ax.set_aspect("auto")
    if box_aspect is not None:
        ax.set_box_aspect(box_aspect)
    if show_title:
        ax.set_title(f"{family.replace('_group_', ' group ')}  |  N={len(subset)}", fontsize=10.2, pad=5)
    draw_lbw_reference_arrow(ax)
    ax.grid(True, linestyle=":", alpha=0.24, linewidth=0.55)
    ax.tick_params(direction="in", top=True, right=True, length=3.5, width=0.75, labelsize=8.5)
    if show_main_legend:
        main_legend = ax.legend(
            handles=lbw_panel_handles(subset, overlay_handles),
            loc="lower left",
            fontsize=SELECTED_LEGEND_FONTSIZE,
            frameon=True,
            framealpha=0.88,
            borderpad=0.32,
            handlelength=1.55,
            handletextpad=0.45,
            labelspacing=0.24,
        )
        main_legend.set_zorder(30)
        ax.add_artist(main_legend)
    median_w = float(subset["vz_lsr_kms"].median())
    ax.text(
        0.98,
        0.96,
        rf"$\bar{{W}}_{{\rm LSR}}={median_w:.1f}$ km/s",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=SELECTED_SPEED_LABEL_FONTSIZE,
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.86, "pad": 1.5},
        zorder=9,
    )


def render_selected_sbs_yso_comparison() -> None:
    xy_catalog = build_selected_xy_catalog()
    ellipse_df = load_compromise_ellipses(require_file(ELLIPSE_CSV))
    xy_family_labels, xy_detail_df, xy_summary_df = prepare_selected_xy_residuals(xy_catalog)
    bubble = load_hbw35_bubble()

    xy_summary_file = OUT_DIR / "G1G2_selected_family_xy_residual_statistics.csv"
    xy_detail_file = OUT_DIR / "G1G2_selected_family_xy_residual_details.csv"
    xy_summary_df.to_csv(xy_summary_file, index=False, encoding="utf-8-sig")
    xy_detail_df[[col for col in DETAIL_COLUMNS if col in xy_detail_df.columns]].to_csv(xy_detail_file, index=False, encoding="utf-8-sig")

    lbw_catalog = build_selected_lbw_catalog()
    lbw_summary_rows = []

    output_file = FINAL_FIG_DIR / "G1G2_SBs_YSO.png"
    fig = plt.figure(figsize=(13.4, 15.8), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.015, h_pad=0.025, wspace=0.005, hspace=0.018)
    grid = fig.add_gridspec(6, 3, height_ratios=[1.0, 0.62, 0.16, 1.0, 0.62, 0.18], hspace=0.018, wspace=0.005)
    upper_lbw_axes: list[plt.Axes] = []
    lower_xy_axes: list[plt.Axes] = []

    for index, family in enumerate(xy_family_labels):
        row_start = 0 if index < 3 else 3
        col = index % 3
        xy_ax = fig.add_subplot(grid[row_start, col])
        lbw_ax = fig.add_subplot(grid[row_start + 1, col])
        if index < 3:
            upper_lbw_axes.append(lbw_ax)
        else:
            lower_xy_axes.append(xy_ax)

        draw_selected_xy_family_panel(
            xy_ax,
            family,
            xy_catalog,
            xy_detail_df,
            xy_summary_df,
            ellipse_df,
            bubble,
            show_legend=False,
            show_orion_eridanus=family in ORION_ERIDANUS_XY_FAMILIES,
        )
        xy_ax.set_xlabel("Heliocentric X (pc)" if col == 1 else "", fontsize=9.2)
        xy_ax.set_ylabel("Heliocentric Y (pc)" if col == 0 else "", fontsize=9.2)

        subset = lbw_catalog.loc[lbw_catalog["plot_family"] == family].copy()
        draw_lbw_family_panel(
            lbw_ax,
            subset,
            family,
            equal_aspect=True,
            aspect_adjustable="datalim",
            show_title=False,
            box_aspect=0.62,
            show_main_legend=False,
        )
        lbw_ax.set_xlabel(r"Galactic longitude $l$ (deg)" if (index >= 3 and col == 1) else "", fontsize=9.2)
        lbw_ax.set_ylabel(r"Galactic latitude $b$ (deg)" if col == 0 else "", fontsize=9.2)
        lbw_summary_rows.append(
            {
                "family": family,
                "n_members": int(len(subset)),
                "median_W_lsr_kms": float(subset["vz_lsr_kms"].median()),
                "mean_age_myr": float(subset["age_myr"].mean()) if "age_myr" in subset and subset["age_myr"].notna().any() else np.nan,
            }
        )

    legend_ax = fig.add_subplot(grid[5, :])
    legend_ax.axis("off")
    shared_legend = legend_ax.legend(
        handles=selected_shared_legend_handles(),
        loc="center",
        ncol=4,
        fontsize=SELECTED_SHARED_LEGEND_FONTSIZE,
        frameon=True,
        framealpha=0.95,
        borderpad=0.36,
        handlelength=1.7,
        handletextpad=0.48,
        columnspacing=1.25,
        labelspacing=0.28,
    )
    shared_legend.set_zorder(30)
    fig.canvas.draw()
    if upper_lbw_axes and lower_xy_axes:
        y_sep = 0.5 * (
            min(ax.get_position().y0 for ax in upper_lbw_axes)
            + max(ax.get_position().y1 for ax in lower_xy_axes)
        )
        x_left = min(ax.get_position().x0 for ax in [*upper_lbw_axes, *lower_xy_axes])
        x_right = max(ax.get_position().x1 for ax in [*upper_lbw_axes, *lower_xy_axes])
        fig.add_artist(Line2D([x_left, x_right], [y_sep, y_sep], transform=fig.transFigure, color="0.25", lw=1.1, alpha=0.85))

    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, dpi=235, bbox_inches="tight")
    plt.close(fig)

    lbw_summary_df = pd.DataFrame(lbw_summary_rows)
    lbw_summary_file = OUT_DIR / "G1G2_selected_family_lb_w_velocity_statistics.csv"
    lbw_summary_df.to_csv(lbw_summary_file, index=False, encoding="utf-8-sig")
    print(f"Saved selected-SB comparison summaries: {xy_summary_file}, {lbw_summary_file}")
    print(xy_summary_df.to_string(index=False))
    print(lbw_summary_df.to_string(index=False))


def render_xy_residual_plot(config: PlotDatasetConfig) -> None:
    families, summary_df = render_xy_dataset(
        config=config.xy_config,
        ellipse_csv=require_file(ELLIPSE_CSV),
        output_ellipse=OUTPUT_ELLIPSE,
    )
    print(f"Saved XY residual plot for {config.output_prefix}.")
    print(summary_df.to_string(index=False))


def render_size_evolution_plot(config: PlotDatasetConfig) -> None:
    output_file = FINAL_FIG_DIR / f"{config.output_prefix}_scale_age.png"
    catalog_file = OUT_DIR / f"{config.output_prefix}_member_summary_table.csv"
    plot_size_evolution_with_mean_age(
        properties_file=require_file(config.properties_file),
        curve_file=require_file(config.curve_file),
        labels_file=require_file(config.labels_file),
        output_file=output_file,
    )
    export_basic_family_catalog(labels_file=config.labels_file, output_file=catalog_file)
    print(f"Saved size-evolution plot and member summary for {config.output_prefix}.")


def render_traceback_slices(config: PlotDatasetConfig) -> None:
    output_file = FINAL_FIG_DIR / f"{config.output_prefix}_traceback_slices.png"
    plot_paperlike_families_across_slices(
        labels_file=str(require_file(config.labels_file)),
        orbit_file=str(require_file(config.orbit_file)),
        output_file=str(output_file),
        step_myr=5.0,
        max_lookback_myr=config.max_lookback_myr,
        xy_limit_pc=config.slice_xy_limit_pc,
        apply_age_cutoff=True,
    )
    print(f"Saved traceback-slice plot for {config.output_prefix}: {output_file}")


def main() -> None:
    """Run Script 16 from validated inputs to the documented outputs."""
    parser = argparse.ArgumentParser(description='Run Script 16: G1G2 YSO traceback plots.')
    parser.add_argument("--dataset", choices=["all", "first", "second"], default="all")
    parser.add_argument(
        "--plot",
        choices=["all", "xy", "scale", "slices", "selected_sbs"],
        default="all",
        help='Command-line option for the documented workflow.',
    )
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)

    dataset_keys = ["first", "second"] if args.dataset == "all" else [args.dataset]
    for key in dataset_keys:
        config = DATASETS[key]
        if args.plot in {"all", "xy"}:
            render_xy_residual_plot(config)
        if args.plot in {"all", "scale"}:
            render_size_evolution_plot(config)
        if args.plot in {"all", "slices"}:
            render_traceback_slices(config)

    if args.plot == "selected_sbs" or (args.plot == "all" and args.dataset == "all"):
        render_selected_sbs_yso_comparison()

    print("Script 16 plotting workflow complete.")

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
