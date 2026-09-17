#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script 24: RW like wave

Purpose
-------
Draw three additional Radcliffe-Wave-like dust/cluster slices.

Method overview
---------------
1. Load the smoothed dust products, young Hunt and G1 cluster samples, and fitted shell
   geometry.
2. Compute the configured oblique-line and constant-Y slices, projecting the dust, shell
   intersections and young-cluster positions into each view.
3. Draw the face-on location overview together with the three slice panels, using the
   configured extinction and vertical-velocity display scales.

Main inputs
-----------
- ../results/intermediate_output/3d_dust_map_products/smoothed_3d_dust_cube.parquet
- ../results/intermediate_output/3d_dust_map_products/xy_mean_dust_map.parquet
- ../results/intermediate_output/18_young_cluster_superbubble_shell_random_test/Hunt_young_cluster_sample.csv
- ../results/intermediate_output/13_traceback_cluster_input_data/G1_traceback_cluster_sample.csv
- ../results/superbubble_final_fit_parameters.csv

Main outputs
------------
- ../results/figures/24_rw_like_three_slices.png; main-text Fig. 4, showing additional Radcliffe-Wave-like structures.

Figure/table role
-----------------
../results/figures/24_rw_like_three_slices.png; main-text Fig. 4, showing additional Radcliffe-Wave-like structures.

Runtime and data notes
----------------------
Reference runtime: 12.18 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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
import importlib.util
import os
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Ellipse, Rectangle
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter


# ----------------------------------------------------------------------------
# Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
# ----------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "24_rw_like_three_slices"
DUST_PRODUCTS_DIR = Path("..") / "results" / "intermediate_output" / "3d_dust_map_products"
TRACEBACK_INPUT_DIR = Path("..") / "results" / "intermediate_output" / "13_traceback_cluster_input_data"
HUNT_OUTPUT_DIR = Path("..") / "results" / "intermediate_output" / "18_young_cluster_superbubble_shell_random_test"

OUT_FIG = FINAL_FIG_DIR / "24_rw_like_three_slices.png"
HUNT_CLUSTER_CSV = HUNT_OUTPUT_DIR / "Hunt_young_cluster_sample.csv"
G1_SAMPLE_CSV = TRACEBACK_INPUT_DIR / "G1_traceback_cluster_sample.csv"
DUST_2D_PARQUET = DUST_PRODUCTS_DIR / "xy_mean_dust_map.parquet"


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


fig3 = _load_module(HERE / "23_RW_cross_sections.py", "script23_radcliffe_fig3")


# ----------------------------------------------------------------------------
# Figure styling and layout configuration for reproducible rendering.
# ----------------------------------------------------------------------------
DEFAULT_SLICES = [
    {"kind": "line", "angle": 50.0, "intercept": -1.70, "extra_cap_section_ids": {7}},
    {"kind": "line", "angle": 127.0, "intercept": 0.90},
    {"kind": "y", "value": -0.90},
]

T_BINS = 301
Z_BINS = 101
TZ_XLIM = (-3.0, 3.0)
TZ_YLIM = (-0.5, 0.5)
LINE_COLORS = ["#d95f02", "#7570b3", "#1b9e77", "#e7298a"]
SB_COLOR_MAP = {0: "#000000", 1: "#0000FF", 2: "#FF0000", 3: "#008000"}
PANEL_DUST_CONFIGS = [
    {"norm": Normalize(vmin=0.0, vmax=1.0), "ticks": [0.0, 0.5, 1.0], "ticklabels": ["0", "0.5", "1"]},
    {"norm": LogNorm(vmin=0.05, vmax=2.0), "ticks": [0.05, 0.1, 0.3, 1.0, 2.0], "ticklabels": ["0.05", "0.1", "0.3", "1", "2"]},
    {"norm": LogNorm(vmin=0.04, vmax=3.0), "ticks": [0.04, 0.1, 0.3, 1.0, 3.0], "ticklabels": ["0.04", "0.1", "0.3", "1", "3"]},
]

FONT_FAMILY = "DejaVu Sans"
LEGEND_FONTSIZE = 11
AXIS_LABEL_FONTSIZE = 11
TICK_FONTSIZE = 9
PANEL_LABEL_FONTSIZE = 12
TEXTBOX_FONTSIZE = 9.5
SB_LABEL_FONTSIZE = 6.4
LEFT_SB_ID_FONTSIZE = SB_LABEL_FONTSIZE + 3.0
RIGHT_SB_ID_FONTSIZE = SB_LABEL_FONTSIZE + 2.0
RIGHT_PANEL_LABEL_ADJUSTMENTS = {
    1: {
        24: {"fontsize_delta": -2.0},
        15: {"dz": -0.05},
        32: {"dz": -0.05},
        16: {"dz": 0.05},
    },
}

VZ_ARROW_SCALE = 0.028
VZ_ARROW_MIN = 0.09
VZ_ARROW_MAX = 0.42
CLUSTER_AGE_CUT_MYR = 30.0
Y_SLICE_SMOOTH_SIGMA = (1.0, 1.2)
ELLIPSOID_CAP_SECTION_FRAC = 0.5
AGE_CMAP_NAME = "twilight_r"
TZ_DUST_CMAP_NAME = "Spectral_r"


def parse_args() -> argparse.Namespace:
    """Build the command-line interface for Script 24."""
    parser = argparse.ArgumentParser(description='Run Script 24: RW like wave.')
    parser.add_argument("--hw", type=float, default=0.05, help='Command-line option for the documented workflow.')
    parser.add_argument("--dpi", type=int, default=300, help='Command-line option for the documented workflow.')
    return parser.parse_args()


def format_delta_ebv_tick(value: float) -> str:
    return f"{value:g}"


def scale_dust_config_to_delta_ebv(config: dict, thickness_kpc: float) -> dict:
    base_norm = config["norm"]
    vmin = base_norm.vmin * thickness_kpc
    vmax = base_norm.vmax * thickness_kpc
    if isinstance(base_norm, LogNorm):
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)
    ticks = [tick * thickness_kpc for tick in config["ticks"]]
    return {
        "norm": norm,
        "ticks": ticks,
        "ticklabels": [format_delta_ebv_tick(tick) for tick in ticks],
    }


def build_y_slice_tz_map(df: pd.DataFrame, y0: float, half_width: float):
    mask = np.abs(df["y"].to_numpy() - y0) <= half_width
    cut = df.loc[mask, ["x", "z", "dust"]]
    t = cut["x"].to_numpy()
    z = cut["z"].to_numpy()
    dust = cut["dust"].to_numpy()
    valid = np.isfinite(dust)

    t_edges = np.linspace(TZ_XLIM[0], TZ_XLIM[1], T_BINS)
    z_edges = np.linspace(TZ_YLIM[0], TZ_YLIM[1], Z_BINS)
    hist, _, _ = np.histogram2d(t[valid], z[valid], bins=[t_edges, z_edges], weights=dust[valid])
    counts, _, _ = np.histogram2d(t[valid], z[valid], bins=[t_edges, z_edges])
    mean_dust = np.divide(
        hist,
        counts,
        out=np.full_like(hist, np.nan, dtype=float),
        where=counts != 0,
    ).T
    mean_dust = smooth_nan_aware(mean_dust, sigma=Y_SLICE_SMOOTH_SIGMA)
    return t_edges, z_edges, mean_dust


def smooth_nan_aware(image: np.ndarray, sigma):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    valid = np.isfinite(image)
    if not np.any(valid):
        return image
    values = np.where(valid, image, 0.0)
    weights = valid.astype(float)
    smooth_values = gaussian_filter(values, sigma=sigma, mode="nearest")
    smooth_weights = gaussian_filter(weights, sigma=sigma, mode="nearest")
    return np.divide(
        smooth_values,
        smooth_weights,
        out=np.full_like(image, np.nan, dtype=float),
        where=smooth_weights > 1e-8,
    )


def prepare_cluster_positions(df, age_cut=CLUSTER_AGE_CUT_MYR):
    needed = ["name", "x_helio", "y_helio", "z_helio", "age_myr"]
    df = df.dropna(subset=needed).copy()
    df = df[df["age_myr"] < age_cut].copy()
    df["x_kpc"] = df["x_helio"] / 1000.0
    df["y_kpc"] = df["y_helio"] / 1000.0
    df["z_kpc"] = df["z_helio"] / 1000.0
    return df


def load_young_hunt_clusters(age_cut=CLUSTER_AGE_CUT_MYR):
    return prepare_cluster_positions(pd.read_csv(HUNT_CLUSTER_CSV), age_cut=age_cut)


def load_young_g1_clusters(age_cut=CLUSTER_AGE_CUT_MYR):
    df = prepare_cluster_positions(pd.read_csv(G1_SAMPLE_CSV), age_cut=age_cut)
    df = df.dropna(subset=["W_lsr"]).copy()
    df["vz_lsr_kms"] = df["W_lsr"]
    return df


def project_clusters_to_line(clusters, theta, point_t0, intercept_kpc, half_width):
    if clusters.empty:
        return clusters.copy()
    x = clusters["x_kpc"].to_numpy()
    y = clusters["y_kpc"].to_numpy()
    signed_dist = -np.sin(theta) * x + np.cos(theta) * y - intercept_kpc * np.cos(theta)
    proj = clusters.copy()
    proj["signed_distance_kpc"] = signed_dist
    proj["T_kpc"] = (x - point_t0[0]) * np.cos(theta) + (y - point_t0[1]) * np.sin(theta)
    proj["Z_kpc"] = proj["z_kpc"]
    return proj[np.abs(signed_dist) <= half_width].copy()


def project_clusters_to_y_slice(clusters, y0, half_width):
    if clusters.empty:
        return clusters.copy()
    proj = clusters.copy()
    proj["signed_distance_kpc"] = proj["y_kpc"] - y0
    proj["T_kpc"] = proj["x_kpc"]
    proj["Z_kpc"] = proj["z_kpc"]
    return proj[np.abs(proj["signed_distance_kpc"]) <= half_width].copy()


def ellipsoid_section_params(cut):
    if cut["cut_kind"] == "rectangle":
        return (
            cut["t_center_kpc"], cut["z_center_kpc"],
            cut["t_radius_kpc"], cut["z_radius_kpc"],
        )
    full_values = [
        cut.get("full_t_center_kpc", np.nan),
        cut.get("full_z_center_kpc", np.nan),
        cut.get("full_t_radius_kpc", np.nan),
        cut.get("full_z_radius_kpc", np.nan),
    ]
    if all(np.isfinite(v) for v in full_values):
        return tuple(full_values)
    return (
        cut["t_center_kpc"], cut["z_center_kpc"],
        cut["t_radius_kpc"], cut["z_radius_kpc"],
    )


def dist_to_arc(t, z, cut):
    tc, zc, tr, zr = ellipsoid_section_params(cut)
    mark = int(cut["mark"])

    if cut["cut_kind"] == "rectangle":
        t_lo, t_hi = tc - tr, tc + tr
        z_lo, z_hi = zc - zr, zc + zr
        n = 60
        ts = np.linspace(t_lo, t_hi, n)
        zs = np.linspace(z_lo, z_hi, n)
        edges_t = np.concatenate([ts, ts, np.full(n, t_lo), np.full(n, t_hi)])
        edges_z = np.concatenate([np.full(n, z_lo), np.full(n, z_hi), zs, zs])
    else:
        if mark == 1:
            angles = np.linspace(np.pi, 2 * np.pi, 120)
        elif mark == 2:
            angles = np.linspace(0, np.pi, 120)
        else:
            angles = np.linspace(0, 2 * np.pi, 180)
        edges_t = tc + tr * np.cos(angles)
        edges_z = zc + zr * np.sin(angles)

    return np.sqrt((t - edges_t) ** 2 + (z - edges_z) ** 2).min()


def assign_clusters_to_nearest_shell(cuts, clusters_on_line):
    records = []
    if cuts.empty or clusters_on_line.empty:
        return records

    shells = []
    for _, cut in cuts.iterrows():
        mark = int(cut["mark"])
        _, _, tr, zr = ellipsoid_section_params(cut)
        if mark in (1, 2) and tr >= 0.01 and zr >= 0.005:
            shells.append(cut)
    if not shells:
        return records

    for _, cl in clusters_on_line.iterrows():
        best_dist = np.inf
        best_shell = None
        for cut in shells:
            dist = dist_to_arc(cl["T_kpc"], cl["Z_kpc"], cut)
            if dist < best_dist:
                best_dist = dist
                best_shell = cut
        if best_shell is None:
            continue
        mark = int(best_shell["mark"])
        vz = cl["vz_lsr_kms"]
        records.append({
            "sb_id": int(best_shell["id"]),
            "mark": mark,
            "consistent": vz < 0 if mark == 1 else vz > 0,
            "dist": best_dist,
        })
    return records


def compute_line_slice(spec, hw, superbubbles, clusters, velocity_clusters, dust_df):
    theta, unit_t, unit_n, point_t0 = fig3.line_geometry(spec["angle"], spec["intercept"])
    cuts = fig3.compute_all_cuts(
        superbubbles, unit_t, point_t0,
        extra_cap_section_ids=spec.get("extra_cap_section_ids"),
    )
    clusters_on_line = project_clusters_to_line(clusters, theta, point_t0, spec["intercept"], hw)
    velocity_clusters_on_line = project_clusters_to_line(
        velocity_clusters, theta, point_t0, spec["intercept"], hw
    )
    t_edges, z_edges, mean_dust = fig3.build_dust_tz_map(
        dust_df, theta, point_t0, spec["intercept"], hw, TZ_XLIM, TZ_YLIM
    )
    records = assign_clusters_to_nearest_shell(cuts, velocity_clusters_on_line)
    n_total = len(records)
    n_consistent = sum(1 for rec in records if rec["consistent"])
    return {
        "kind": "line",
        "angle": spec["angle"],
        "intercept": spec["intercept"],
        "theta": theta,
        "unit_t": unit_t,
        "unit_n": unit_n,
        "point_t0": point_t0,
        "cuts": cuts,
        "clusters_on_line": clusters_on_line,
        "velocity_clusters_on_line": velocity_clusters_on_line,
        "t_edges": t_edges,
        "z_edges": z_edges,
        "mean_dust": mean_dust,
        "n_total": n_total,
        "n_consistent": n_consistent,
        "n_shells": len({rec["sb_id"] for rec in records}),
        "consistency": n_consistent / n_total if n_total else np.nan,
        "n_sb_cut": len(cuts),
    }


def compute_y_slice(spec, hw, superbubbles, clusters, velocity_clusters, dust_df):
    y0 = spec["value"]
    unit_t = np.array([1.0, 0.0])
    unit_n = np.array([0.0, 1.0])
    point_t0 = np.array([0.0, y0])
    cuts = fig3.compute_all_cuts(superbubbles, unit_t, point_t0)
    clusters_on_line = project_clusters_to_y_slice(clusters, y0, hw)
    velocity_clusters_on_line = project_clusters_to_y_slice(velocity_clusters, y0, hw)
    t_edges, z_edges, mean_dust = build_y_slice_tz_map(dust_df, y0, hw)
    records = assign_clusters_to_nearest_shell(cuts, velocity_clusters_on_line)
    n_total = len(records)
    n_consistent = sum(1 for rec in records if rec["consistent"])
    return {
        "kind": "y",
        "y0": y0,
        "theta": 0.0,
        "unit_t": unit_t,
        "unit_n": unit_n,
        "point_t0": point_t0,
        "cuts": cuts,
        "clusters_on_line": clusters_on_line,
        "velocity_clusters_on_line": velocity_clusters_on_line,
        "t_edges": t_edges,
        "z_edges": z_edges,
        "mean_dust": mean_dust,
        "n_total": n_total,
        "n_consistent": n_consistent,
        "n_shells": len({rec["sb_id"] for rec in records}),
        "consistency": n_consistent / n_total if n_total else np.nan,
        "n_sb_cut": len(cuts),
    }


CLUSTER_STAR_SIZE = 82.0


def age_to_size(age_myr):
    return np.full_like(np.asarray(age_myr, dtype=float), CLUSTER_STAR_SIZE)


def draw_young_clusters_and_g1_vz(ax, clusters_on_line, velocity_clusters_on_line, age_norm, age_cmap):
    if clusters_on_line.empty and velocity_clusters_on_line.empty:
        return
    if not clusters_on_line.empty:
        ax.scatter(
            clusters_on_line["T_kpc"], clusters_on_line["Z_kpc"],
            marker="*", s=age_to_size(clusters_on_line["age_myr"].to_numpy()),
            c=clusters_on_line["age_myr"].to_numpy(), cmap=age_cmap, norm=age_norm,
            edgecolor="black", linewidth=0.4, zorder=9,
        )
    for _, row in velocity_clusters_on_line.iterrows():
        vz = row["vz_lsr_kms"]
        sign = 1.0 if vz >= 0 else -1.0
        arrow_dz = sign * np.clip(abs(vz * VZ_ARROW_SCALE), VZ_ARROW_MIN, VZ_ARROW_MAX)
        arrow_color = age_cmap(age_norm(row["age_myr"]))
        ann = ax.annotate(
            "", xy=(row["T_kpc"], row["Z_kpc"] + arrow_dz), xytext=(row["T_kpc"], row["Z_kpc"]),
            arrowprops={"arrowstyle": "-|>", "color": arrow_color, "linewidth": 1.1,
                        "mutation_scale": 10, "shrinkA": 2, "shrinkB": 0},
            zorder=10,
        )
        if ann.arrow_patch is not None:
            ann.arrow_patch.set_path_effects([
                pe.withStroke(linewidth=2.4, foreground="black"),
                pe.Normal(),
            ])
    if not velocity_clusters_on_line.empty:
        ax.scatter(
            velocity_clusters_on_line["T_kpc"], velocity_clusters_on_line["Z_kpc"],
            marker="*", s=age_to_size(velocity_clusters_on_line["age_myr"].to_numpy()),
            c=velocity_clusters_on_line["age_myr"].to_numpy(), cmap=age_cmap, norm=age_norm,
            edgecolor="black", linewidth=0.45, zorder=11,
        )


def load_xy_overview_points():
    df = pd.read_parquet(DUST_2D_PARQUET)
    return df.rename(columns={"X": "x_round", "Y": "y_round"})


def compromise_section_scale(row):
    shape = str(row["shape"]).strip().lower()
    if shape != "ellipsoid":
        return 1.0
    cz = float(row["center_z_kpc"])
    c = float(row["c_radius_kpc"])
    z_base = float(row["xy_plane_z_kpc"])
    mark = int(row["mark"]) if np.isfinite(row["mark"]) else -1
    if not all(np.isfinite(v) for v in [cz, c, z_base]) or c <= 0:
        return 1.0
    z_apex = cz - c if mark == 1 else cz + c
    z_sec = z_base + ELLIPSOID_CAP_SECTION_FRAC * (z_apex - z_base)
    with np.errstate(invalid="ignore", divide="ignore"):
        k_base = np.sqrt(np.clip(1.0 - ((z_base - cz) / c) ** 2, 0.0, None))
        k_sec = np.sqrt(np.clip(1.0 - ((z_sec - cz) / c) ** 2, 0.0, None))
    return float(k_sec / k_base) if k_base > 0 else 1.0


def compromise_xy_ellipse(row):
    cx = row.get("xy_plane_center_x_kpc", row["center_x_kpc"])
    cy = row.get("xy_plane_center_y_kpc", row["center_y_kpc"])
    a = row.get("xy_plane_a_kpc", row["a_radius_kpc"])
    b = row.get("xy_plane_b_kpc", row["b_radius_kpc"])
    angle = row.get("xy_plane_angle_deg", row["angle_deg"])
    scale = compromise_section_scale(row)
    return cx, cy, a * scale, b * scale, angle


def draw_xy_overview(ax, xy_dust, superbubbles, results, hw):
    ax.scatter(
        xy_dust["x_round"], xy_dust["y_round"], c=xy_dust["dust"],
        cmap="Greys", s=9, alpha=0.72, linewidths=0,
        vmin=0, vmax=0.3, rasterized=True, zorder=1,
    )
    for _, row in superbubbles.iterrows():
        sb_id = int(row["id"]) if np.isfinite(row["id"]) else -1
        cx, cy, a, b, angle = compromise_xy_ellipse(row)
        if not all(np.isfinite(v) for v in [cx, cy, a, b, angle]) or a <= 0 or b <= 0:
            continue
        mark = int(row["mark"]) if np.isfinite(row["mark"]) else -1
        color = SB_COLOR_MAP.get(mark, "gray")
        zorder = 2
        ax.add_patch(Ellipse((cx, cy), width=2 * a, height=2 * b, angle=angle,
                             facecolor="none", edgecolor="white",
                             linewidth=1.65, alpha=0.86, zorder=zorder))
        ax.add_patch(Ellipse((cx, cy), width=2 * a, height=2 * b, angle=angle,
                             facecolor="none", edgecolor=color,
                             linewidth=1.15, alpha=1.0, zorder=zorder + 0.1))
        ax.text(cx, cy, f"{sb_id}", color=color, fontsize=LEFT_SB_ID_FONTSIZE,
                fontweight="bold", ha="center", va="center", zorder=zorder + 2,
                path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])

    x_line = np.linspace(-3.4, 3.4, 700)
    t_vals = np.arange(-3, 4)
    for idx, (result, color) in enumerate(zip(results, LINE_COLORS), start=1):
        panel_letter = chr(ord("a") + idx - 1)
        unit_n = result["unit_n"]
        point_t0 = result["point_t0"]
        if result["kind"] == "line":
            theta = result["theta"]
            y_line = np.tan(theta) * x_line + result["intercept"]
            x_t = point_t0[0] + t_vals * np.cos(theta)
            y_t = point_t0[1] + t_vals * np.sin(theta)
            dx = -np.sin(theta) * 0.15
            dy = np.cos(theta) * 0.15
            legend_label = f"{panel_letter}: y=tan({result['angle']:g}°)x{result['intercept']:+.2f}"
        else:
            y_line = np.full_like(x_line, result["y0"])
            x_t = t_vals.astype(float)
            y_t = np.full_like(x_t, result["y0"], dtype=float)
            dx, dy = 0.0, 0.15
            legend_label = f"{panel_letter}: y={result['y0']:+.2f}"

        upper = np.vstack((x_line, y_line)) + unit_n[:, None] * hw
        lower = np.vstack((x_line, y_line)) - unit_n[:, None] * hw
        ax.fill(np.concatenate([upper[0], lower[0][::-1]]),
                np.concatenate([upper[1], lower[1][::-1]]),
                color=color, alpha=0.15, zorder=5)
        ax.plot(x_line, y_line, color=color, linewidth=2.8, zorder=7, label=legend_label)
        ax.plot(x_t, y_t, "o", color=color, markeredgecolor="white",
                markeredgewidth=0.45, markersize=3.4, alpha=0.95, zorder=8)

        for tv, xt, yt in zip(t_vals, x_t, y_t):
            if -3.05 <= xt <= 3.05 and -3.05 <= yt <= 3.05:
                ax.text(xt + dx, yt + dy, f"{tv:g}", fontsize=TEXTBOX_FONTSIZE,
                        ha="center", va="center", color=color, fontweight="bold", zorder=9,
                        path_effects=[pe.withStroke(linewidth=1.6, foreground="white")])

    ax.axhline(0, color="black", linewidth=0.5, alpha=0.55, zorder=3)
    ax.axvline(0, color="black", linewidth=0.5, alpha=0.55, zorder=3)
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    ax.set_aspect("equal", adjustable="box")
    ax.set_anchor("N")
    ax.grid(True, color="0.78", linestyle="--", linewidth=0.45, alpha=0.55)
    ax.set_xlabel("X [kpc]", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("Y [kpc]", fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE, direction="out", length=3, width=0.7)
    legend = ax.legend(loc="lower left", fontsize=LEGEND_FONTSIZE, frameon=True, framealpha=0.82,
                       borderpad=0.35, handlelength=2.0)
    legend.set_zorder(30)


def draw_cut(ax, cut, color):
    tc, zc, tr, zr = ellipsoid_section_params(cut)
    if cut["cut_kind"] == "rectangle":
        ax.add_patch(Rectangle(
            (tc - tr, zc - 2 * zr),
            width=2 * tr, height=4 * zr,
            edgecolor=color, facecolor="none", linewidth=1.0, alpha=0.9, zorder=6,
        ))
        return
    mark = int(cut["mark"])
    if mark == 1:
        angles = np.linspace(np.pi, 2 * np.pi, 180)
    elif mark == 2:
        angles = np.linspace(0, np.pi, 180)
    else:
        angles = np.linspace(0, 2 * np.pi, 240)
    t = tc + tr * np.cos(angles)
    z = zc + zr * np.sin(angles)
    ax.plot(t, z, color=color, linewidth=1.1, alpha=0.95, zorder=6)


def draw_slice_panel(ax, result, norm, age_norm, age_cmap, panel_index, dust_thickness_kpc):
    delta_ebv = result["mean_dust"] * dust_thickness_kpc
    mean_clipped = np.clip(delta_ebv, norm.vmin, norm.vmax)
    mesh = ax.pcolormesh(result["t_edges"], result["z_edges"], mean_clipped,
                         cmap=TZ_DUST_CMAP_NAME, shading="flat", norm=norm, rasterized=True)

    for _, cut in result["cuts"].iterrows():
        color = "black"
        draw_cut(ax, cut, color)
        label_t, label_z, _, label_z_radius = ellipsoid_section_params(cut)
        mark = int(cut["mark"])
        if mark == 1:
            label_z -= 0.35 * label_z_radius
        elif mark == 2:
            label_z += 0.35 * label_z_radius
        sb_id = int(cut["id"])
        label_adjust = RIGHT_PANEL_LABEL_ADJUSTMENTS.get(panel_index, {}).get(sb_id, {})
        label_z += label_adjust.get("dz", 0.0)
        label_fontsize = RIGHT_SB_ID_FONTSIZE + label_adjust.get("fontsize_delta", 0.0)
        ax.text(label_t, label_z, f"SB{sb_id}",
                fontsize=label_fontsize, fontweight="bold", ha="center", va="center",
                color=color, zorder=7,
                path_effects=[pe.withStroke(linewidth=1.4, foreground="white")])

    draw_young_clusters_and_g1_vz(
        ax, result["clusters_on_line"], result["velocity_clusters_on_line"], age_norm, age_cmap
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(*TZ_XLIM)
    ax.set_ylim(*TZ_YLIM)
    ax.set_ylabel("Z [kpc]", fontsize=AXIS_LABEL_FONTSIZE)
    ax.grid(True, color="0.78", linestyle="--", linewidth=0.45, alpha=0.55)
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE, direction="out", length=3, width=0.7)
    return mesh


def main():
    """Run Script 24 from validated inputs to the documented outputs."""
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        "font.family": FONT_FAMILY,
        "font.size": TICK_FONTSIZE,
        "axes.linewidth": 0.75,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })

    print("Loading dust, OSB, and young-cluster inputs.")
    dust_df = fig3.load_dust_subset()
    superbubbles = fig3.load_superbubbles()
    clusters = load_young_hunt_clusters()
    velocity_clusters = load_young_g1_clusters()
    print(
        f"Loaded dust cells={len(dust_df)}, superbubbles={len(superbubbles)}, "
        f"Hunt clusters={len(clusters)}, velocity-selected G1 clusters={len(velocity_clusters)}."
    )

    print("Rendering RW-like line and slab slices.")
    results = []
    for idx, spec in enumerate(DEFAULT_SLICES):
        if spec["kind"] == "line":
            result = compute_line_slice(spec, args.hw, superbubbles, clusters, velocity_clusters, dust_df)
            desc = f"angle={result['angle']:g}, intercept={result['intercept']:+.2f}"
        else:
            result = compute_y_slice(spec, args.hw, superbubbles, clusters, velocity_clusters, dust_df)
            desc = f"y={result['y0']:+.2f}"
        result["line_color"] = LINE_COLORS[idx]
        results.append(result)
        ratio = "NA" if np.isnan(result["consistency"]) else f"{result['consistency']:.0%}"
        print(f"  {desc}: {result['n_consistent']}/{result['n_total']} ({ratio}), "
              f"{result['n_shells']} shells, {result['n_sb_cut']} SB cuts")

    print("Loading XY overview dust points.")
    xy_dust = load_xy_overview_points()

    print("[4/5] plotting...")
    dust_thickness_kpc = 2.0 * args.hw
    panel_dust_configs = [
        scale_dust_config_to_delta_ebv(config, dust_thickness_kpc)
        for config in PANEL_DUST_CONFIGS
    ]
    age_norm = Normalize(vmin=0.0, vmax=CLUSTER_AGE_CUT_MYR)
    age_cmap = mpl.colormaps.get_cmap(AGE_CMAP_NAME)
    fig = plt.figure(figsize=(16.0, 7.4))
    ax_xy = fig.add_axes([0.055, 0.155, 0.350, 0.757])
    right_x, right_w = 0.455, 0.455
    right_h = 0.160
    right_ys = [0.705, 0.455, 0.205]
    axes = [fig.add_axes([right_x, y, right_w, right_h]) for y in right_ys]

    draw_xy_overview(ax_xy, xy_dust, superbubbles, results, args.hw)

    for i, (ax, result, dust_config) in enumerate(zip(axes, results, panel_dust_configs)):
        mesh_ref = draw_slice_panel(
            ax, result, dust_config["norm"], age_norm, age_cmap, i, dust_thickness_kpc
        )
        ax.text(0.006, 0.08, f"({chr(ord('a') + i)})", transform=ax.transAxes,
                fontsize=PANEL_LABEL_FONTSIZE, fontweight="bold", ha="left", va="bottom",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.8},
                zorder=12)
        cax = fig.add_axes([0.925, right_ys[i], 0.012, right_h])
        cbar = fig.colorbar(mesh_ref, cax=cax)
        cbar.set_ticks(dust_config["ticks"])
        cbar.set_ticklabels(dust_config["ticklabels"])
        if i == 1:
            cbar.set_label(r"$\Delta E(B-V)$ [mag]", fontsize=AXIS_LABEL_FONTSIZE)
        cbar.ax.tick_params(labelsize=TICK_FONTSIZE, length=3, width=0.7)
        cbar.outline.set_linewidth(0.7)
    axes[-1].set_xlabel("T [kpc]", fontsize=AXIS_LABEL_FONTSIZE)

    age_cax = fig.add_axes([right_x, 0.095, right_w, 0.018])
    age_sm = mpl.cm.ScalarMappable(norm=age_norm, cmap=age_cmap)
    age_cbar = fig.colorbar(age_sm, cax=age_cax, orientation="horizontal")
    age_cbar.set_label("Cluster age [Myr]", fontsize=AXIS_LABEL_FONTSIZE)
    age_cbar.set_ticks([0, 10, 20, 30])
    age_cbar.ax.tick_params(labelsize=TICK_FONTSIZE, length=3, width=0.7)
    age_cbar.outline.set_linewidth(0.7)

    print(f"Saving RW-like wave figure: {OUT_FIG}")
    fig.savefig(OUT_FIG, dpi=args.dpi)
    plt.close(fig)
    print(f"  saved: {OUT_FIG}")

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
