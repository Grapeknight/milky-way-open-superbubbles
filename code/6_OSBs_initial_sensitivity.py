# -*- coding: utf-8 -*-
"""Script 6: OSBs initial sensitivity

Purpose
-------
Perturb initial geometry and retest first-contact cloud recovery.

Method overview
---------------
1. Construct the reference seed for each target from the final parameter table and
   select its first-contact clouds.
2. Perturb the initial seed, repeat the same cloud-selection procedure, and compare the
   selected clouds with the reference selection.
3. Save per-realization overlap measurements, per-target summaries and the sensitivity
   figure; these trials test initialization sensitivity.

Main inputs
-----------
- ../results/superbubble_final_fit_parameters.csv
- ../data/MCs.csv

Main outputs
------------
- ../results/figures/6_initial_sensitivity_overlap.png; initial-condition sensitivity figure.

Figure/table role
-----------------
../results/figures/6_initial_sensitivity_overlap.png; initial-condition sensitivity figure.

Runtime and data notes
----------------------
Reference runtime: 37 min 30.93 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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

import argparse
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "6_initial_condition_sensitivity"
FINAL_FIG_DIR = Path("..") / "results" / "figures"

# Geometry and shell-mask conventions used by this analysis stage.
FINAL_TABLE_PRIMARY = Path("..") / "results" / "superbubble_final_fit_parameters.csv"
MCS_CSV = DATA_DIR / "MCs.csv"

XY_SEARCH_SCALE = 1.5
N_AZIMUTH = 72
N_ELEVATION = 18
RAY_CONE_APERTURE_DEG = 5.0
RADIUS_TOLERANCE_PC = 0.0
MAX_HITS_PER_RAY = 1
SEARCH_SCALE = 1.5

DEFAULT_N_PERTURB = 1000
DEFAULT_CENTER_PERTURB_FRAC = 0.10
DEFAULT_AXIS_PERTURB_FRAC = 0.10
DEFAULT_RANDOM_SEED = 20260529


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def setup_modules():
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    s1 = _load_module(HERE / "1_OSBs_first_contact_MCs.py", "pr_stage1_core_for_sensitivity")
    s2 = _load_module(HERE / "2_OSBs_3D_fitting.py", "pr_stage2_core_for_sensitivity")
    return s1, s2


# Data-loading helper section.
def resolve_final_table_path() -> Path:
    if FINAL_TABLE_PRIMARY.exists():
        return FINAL_TABLE_PRIMARY
    raise FileNotFoundError(
        f"Final geometry table is missing. Run Script 3 before Script 6: {FINAL_TABLE_PRIMARY}"
    )


def load_final_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [str(c).lstrip("﻿") for c in df.columns]
    df["id"] = df["id"].astype(int)
    df["mark"] = df["mark"].astype(int)
    required = [
        "id", "mark", "center_x_kpc", "center_y_kpc", "center_z_kpc",
        "a_radius_kpc", "b_radius_kpc", "c_radius_kpc", "angle_deg",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"final parameter table is missing required columns: {', '.join(missing)}")
    return df


def parse_targets(text: str) -> list[int]:
    values = []
    for item in str(text).replace(";", ",").split(","):
        item = item.strip()
        if item:
            values.append(int(item))
    return sorted(set(values))


def build_seed_from_row(auto, row):
    """Compute an intermediate statistic or table used by this documented workflow."""
    return auto.BubbleSeed(
        target_index=int(row["id"]),
        center=np.array(
            [row["center_x_kpc"], row["center_y_kpc"], row["center_z_kpc"]],
            dtype=float,
        ),
        axes_init=np.array(
            [row["a_radius_kpc"], row["b_radius_kpc"], row["c_radius_kpc"]],
            dtype=float,
        ),
        angle_deg=float(row["angle_deg"]),
        mark=int(row["mark"]),
        name=f"SB{int(row['id'])}",
    )


def find_first_touch_fast(auto, directions: np.ndarray, seed, candidate: pd.DataFrame, fit_shape: str):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    n_candidate = int(len(candidate))
    if n_candidate == 0:
        return {}, 0, 0
    radius_tol_kpc = RADIUS_TOLERANCE_PC / 1000.0
    cone_half = np.deg2rad(RAY_CONE_APERTURE_DEG / 2.0)

    offsets = candidate[["x", "y", "z"]].to_numpy(dtype=float) - seed.center[None, :]
    radii = candidate["rad_kpc"].to_numpy(dtype=float) + radius_tol_kpc
    dist = np.linalg.norm(offsets, axis=1)
    safe = np.maximum(dist, 1e-12)
    ratio = radii / safe
    env = ratio >= 1.0
    ang_rad = np.where(env, np.pi, np.arcsin(np.clip(ratio, 0.0, 1.0)))
    unit = offsets / safe[:, None]

    # Geometry and shell-mask conventions used by this analysis stage.
    rot = auto.rotation_matrix_z(-float(seed.angle_deg))
    local = directions @ rot.T
    a, b, c = float(seed.axes_init[0]), float(seed.axes_init[1]), float(seed.axes_init[2])
    if fit_shape == "cylinder":
        xy_norm = np.hypot(local[:, 0], local[:, 1])
        denom = (local[:, 0] / a) ** 2 + (local[:, 1] / b) ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            t_guess = 1.0 / np.sqrt(denom)
        t_guess = np.where((xy_norm < 0.2) | (denom <= 0), np.nan, t_guess)
    else:
        denom = (local[:, 0] / a) ** 2 + (local[:, 1] / b) ** 2 + (local[:, 2] / c) ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            t_guess = 1.0 / np.sqrt(denom)
    valid_ray = np.isfinite(t_guess) & (t_guess > 0)
    t_max = SEARCH_SCALE * t_guess

    cos_phi = np.clip(directions @ unit.T, -1.0, 1.0)         # (R, C)
    phi = np.arccos(cos_phi)
    intersect = phi <= (cone_half + ang_rad[None, :])
    phi_minus_alpha = np.maximum(phi - cone_half, 0.0)
    perp = dist[None, :] * np.sin(phi_minus_alpha)
    along = dist[None, :] * np.cos(phi_minus_alpha)
    half_chord = np.sqrt(np.maximum(radii[None, :] ** 2 - perp ** 2, 0.0))
    t_entry = np.maximum(along - half_chord, 0.0)
    t_exit = along + half_chord
    within = (t_exit >= 0.0) & (t_entry <= t_max[:, None])
    hit = intersect & within & valid_ray[:, None]

    t_sel = np.where(hit, t_entry, np.inf)
    best_cloud = np.argmin(t_sel, axis=1)
    ray_idx = np.arange(directions.shape[0])
    has_hit = np.isfinite(t_sel[ray_idx, best_cloud])
    first_touch = best_cloud[has_hit]
    counts = np.bincount(first_touch, minlength=n_candidate)
    seqs = candidate["Seq"].to_numpy()
    weights = {int(seqs[i]): float(counts[i]) for i in np.nonzero(counts)[0]}
    return weights, int(has_hit.sum()), n_candidate


def find_first_touch_canonical(stage1, seed, candidate, rays, fit_shape):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    selected, ray_hits = stage1.find_first_touch_hits_v4(
        candidate, seed, rays,
        RADIUS_TOLERANCE_PC / 1000.0, RAY_CONE_APERTURE_DEG,
        MAX_HITS_PER_RAY, fit_shape, SEARCH_SCALE,
    )
    if not len(selected) or "Seq" not in selected.columns:
        return {}, 0, int(len(candidate))
    counts = selected["ray_hit_count"].to_numpy(dtype=float) if "ray_hit_count" in selected.columns else np.ones(len(selected))
    weights = {int(s): float(w) for s, w in zip(selected["Seq"].to_numpy(), counts)}
    return weights, int(len(ray_hits)), int(len(candidate))


def select_first_touch(stage1, core, filtered_clouds, seed, directions, rays, fit_shape, use_canonical):
    candidate = core.preselect_clouds_v4(filtered_clouds, seed, XY_SEARCH_SCALE)
    if use_canonical:
        return find_first_touch_canonical(stage1, seed, candidate, rays, fit_shape)
    return find_first_touch_fast(core, directions, seed, candidate, fit_shape)


def perturb_seed(auto, seed, rng, center_perturb_frac, axis_perturb_frac):
    axes = np.asarray(seed.axes_init, dtype=float)

    local_delta = rng.uniform(-float(center_perturb_frac), float(center_perturb_frac), size=3) * axes
    global_xy = auto.rotate_points_to_global(
        np.array([[local_delta[0], local_delta[1], 0.0]], dtype=float), seed.angle_deg
    )[0]
    global_delta = np.array([global_xy[0], global_xy[1], local_delta[2]], dtype=float)

    axis_factors = 1.0 + rng.uniform(-float(axis_perturb_frac), float(axis_perturb_frac), size=3)
    pert_seed = auto.BubbleSeed(
        target_index=seed.target_index,
        center=np.asarray(seed.center, dtype=float) + global_delta,
        axes_init=axes * axis_factors,
        angle_deg=float(seed.angle_deg),
        mark=int(seed.mark),
        name=seed.name,
    )
    return pert_seed, local_delta, global_delta, axis_factors


def overlap_metrics(baseline_weights: dict, selected_weights: dict) -> dict:
    baseline_set = set(baseline_weights)
    selected_set = set(selected_weights)
    overlap = baseline_set & selected_set
    union = baseline_set | selected_set
    baseline_total = float(sum(baseline_weights.values()))
    selected_total = float(sum(selected_weights.values()))
    recovered = float(sum(baseline_weights[s] for s in overlap))
    common_min = float(sum(min(baseline_weights.get(s, 0.0), selected_weights.get(s, 0.0)) for s in union))
    union_max = float(sum(max(baseline_weights.get(s, 0.0), selected_weights.get(s, 0.0)) for s in union))
    return {
        "baseline_count": int(len(baseline_set)),
        "perturbed_count": int(len(selected_set)),
        "overlap_count": int(len(overlap)),
        "baseline_ray_hit_weight": baseline_total,
        "perturbed_ray_hit_weight": selected_total,
        "recovered_baseline_ray_hit_weight": recovered,
        "weighted_overlap_fraction": float(recovered / baseline_total) if baseline_total > 0 else 0.0,
        "unweighted_overlap_fraction": float(len(overlap) / max(len(baseline_set), 1)),
        "precision_fraction": float(len(overlap) / max(len(selected_set), 1)) if selected_set else 0.0,
        "weighted_jaccard_fraction": float(common_min / union_max) if union_max > 0 else 0.0,
        "unweighted_jaccard_fraction": float(len(overlap) / max(len(union), 1)),
        "exact_same_set": bool(baseline_set == selected_set),
    }


def summarize_target(detail_df: pd.DataFrame) -> dict:
    out = {
        "n_trials": int(len(detail_df)),
        "baseline_cloud_count": int(detail_df["baseline_count"].iloc[0]),
        "baseline_ray_hit_weight": float(detail_df["baseline_ray_hit_weight"].iloc[0]),
    }
    for col in ["weighted_overlap_fraction", "unweighted_overlap_fraction", "weighted_jaccard_fraction"]:
        v = detail_df[col].to_numpy(dtype=float)
        out[f"{col}_mean"] = float(np.mean(v))
        out[f"{col}_median"] = float(np.median(v))
        out[f"{col}_p10"] = float(np.percentile(v, 10.0))
        out[f"{col}_p90"] = float(np.percentile(v, 90.0))
        out[f"{col}_min"] = float(np.min(v))
    wv = detail_df["weighted_overlap_fraction"].to_numpy(dtype=float)
    out["exact_same_set_fraction"] = float(np.mean(detail_df["exact_same_set"].to_numpy(dtype=bool)))
    out["runs_ge_0p90_weighted_fraction"] = float(np.mean(wv >= 0.90))
    out["runs_ge_0p80_weighted_fraction"] = float(np.mean(wv >= 0.80))
    out["mean_selected_cloud_count"] = float(np.mean(detail_df["perturbed_count"].to_numpy(dtype=float)))
    out["mean_candidate_cloud_count"] = float(np.mean(detail_df["candidate_cloud_count"].to_numpy(dtype=float)))
    return out


def save_summary_figure(summary_df: pd.DataFrame, overall: dict, output_path: Path):
    df = summary_df.sort_values("id").reset_index(drop=True)
    x = np.arange(len(df))
    labels = [f"SB{int(v)}" for v in df["id"]]

    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 9,
        "axes.titlesize": 10, "xtick.labelsize": 7, "ytick.labelsize": 8,
        "legend.fontsize": 7.5, "axes.linewidth": 0.8, "xtick.major.width": 0.8,
        "ytick.major.width": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(8.4, 3.2), constrained_layout=True)
    mean_pct = 100.0 * df["weighted_overlap_fraction_mean"].to_numpy(dtype=float)
    p10_pct = 100.0 * df["weighted_overlap_fraction_p10"].to_numpy(dtype=float)
    p90_pct = 100.0 * df["weighted_overlap_fraction_p90"].to_numpy(dtype=float)
    yerr = np.vstack([mean_pct - p10_pct, p90_pct - mean_pct])
    colors = np.where(mean_pct >= 90.0, "#2A9D8F", np.where(mean_pct >= 80.0, "#E9C46A", "#E76F51"))
    ax.bar(x, mean_pct, width=0.76, color=colors, edgecolor="0.12", linewidth=0.45, zorder=3)
    ax.errorbar(x, mean_pct, yerr=yerr, fmt="none", ecolor="0.12", elinewidth=0.85,
                capsize=1.8, capthick=0.85, zorder=4)
    global_mean = 100.0 * overall["global_ray_weighted_overlap_mean"]
    ax.axhline(global_mean, color="black", linestyle="-", linewidth=1.0,
               label=f"Global weighted mean = {global_mean:.1f}%")
    ax.axhline(90.0, color="0.35", linestyle=(0, (3.0, 2.0)), linewidth=0.8, label="90% reference")
    ax.set_ylabel("Weighted overlap (%)")
    ax.set_ylim(50.0, 102.5)
    ax.set_xticks(x, labels, rotation=45, ha="right")
    ax.set_yticks(np.arange(50, 101, 10))
    ax.grid(True, axis="y", color="0.88", linestyle="-", linewidth=0.55, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", length=2.5, pad=1.5)
    ax.tick_params(axis="y", length=3.0)
    ax.legend(loc="upper right", frameon=False, handlelength=2.8, borderpad=0.2)
    ax.text(0.995, 0.04,
            f"{int(overall['n_targets'])} targets x {int(overall['n_trials_per_target'])} trials; "
            f"center & axes +/-{int(round(100 * overall['center_perturb_frac']))}%; "
            f"error bars: 10th-90th percentiles",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.8, color="0.25")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=240)
    plt.close(fig)


# Main workflow.
def build_parser():
    """Build the command-line interface for Script 6."""
    p = argparse.ArgumentParser(description='Run Script 6: OSBs initial sensitivity.')
    p.add_argument("--targets", default="", help='Comma-separated target IDs; use all for the full sample.')
    p.add_argument("--n-perturb", type=int, default=DEFAULT_N_PERTURB, help='Command-line option for the documented workflow.')
    p.add_argument("--center-perturb-frac", type=float, default=DEFAULT_CENTER_PERTURB_FRAC, help='Command-line option for the documented workflow.')
    p.add_argument("--axis-perturb-frac", type=float, default=DEFAULT_AXIS_PERTURB_FRAC, help='Command-line option for the documented workflow.')
    p.add_argument("--random-seed", type=int, default=DEFAULT_RANDOM_SEED, help='Random seed for reproducible controls.')
    p.add_argument("--use-canonical", action="store_true", help='Command-line option for the documented workflow.')
    return p


def main():
    """Run Script 6 from validated inputs to the documented outputs."""
    args = build_parser().parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)
    stage1, core = setup_modules()

    final_table_path = resolve_final_table_path()
    table = load_final_table(final_table_path)
    if str(args.targets).strip():
        wanted = set(parse_targets(args.targets))
        table = table[table["id"].isin(wanted)].copy()
    table = table.sort_values("id").reset_index(drop=True)
    targets = table["id"].tolist()

    if not MCS_CSV.exists():
        raise FileNotFoundError(f"molecular-cloud catalogue does not exist: {MCS_CSV}")
    all_clouds = core.parse_mcs_catalog(MCS_CSV)
    filtered_clouds = core.filter_local_bubble_clouds(all_clouds)
    rays = stage1.prepare_full_sky_rays_v4(N_AZIMUTH, N_ELEVATION)
    directions = np.asarray([r["direction"] for r in rays], dtype=float)
    rng = np.random.default_rng(int(args.random_seed))

    print(f"Using final geometry table: {final_table_path}")
    print(f"Testing {len(table)} OSBs with {args.n_perturb} perturbations each.")
    print(f"Random seed: {args.random_seed}")

    detail_rows = []
    summary_rows = []
    t0 = time.perf_counter()

    for _, row in table.iterrows():
        target_id = int(row["id"])
        seed = build_seed_from_row(core, row)
        fit_shape = core.fit_shape_for_seed(seed)

        baseline_weights, _, _ = select_first_touch(
            stage1, core, filtered_clouds, seed, directions, rays, fit_shape, args.use_canonical
        )

        target_detail = []
        for trial in range(int(args.n_perturb)):
            pert_seed, local_delta, global_delta, axis_factors = perturb_seed(
                core, seed, rng, args.center_perturb_frac, args.axis_perturb_frac
            )
            selected_weights, n_hit_rays, n_candidate = select_first_touch(
                stage1, core, filtered_clouds, pert_seed, directions, rays, fit_shape, args.use_canonical
            )
            metrics = overlap_metrics(baseline_weights, selected_weights)
            row_out = {
                "id": target_id,
                "trial": int(trial),
                "fit_shape": fit_shape,
                "mark": int(seed.mark),
                "local_dx_kpc": float(local_delta[0]),
                "local_dy_kpc": float(local_delta[1]),
                "local_dz_kpc": float(local_delta[2]),
                "global_dx_kpc": float(global_delta[0]),
                "global_dy_kpc": float(global_delta[1]),
                "global_dz_kpc": float(global_delta[2]),
                "axis_factor_a": float(axis_factors[0]),
                "axis_factor_b": float(axis_factors[1]),
                "axis_factor_c": float(axis_factors[2]),
                "pert_a_kpc": float(pert_seed.axes_init[0]),
                "pert_b_kpc": float(pert_seed.axes_init[1]),
                "pert_c_kpc": float(pert_seed.axes_init[2]),
                "candidate_cloud_count": int(n_candidate),
                "selected_cloud_count": int(len(selected_weights)),
                "ray_hit_count": int(n_hit_rays),
                **metrics,
            }
            detail_rows.append(row_out)
            target_detail.append(row_out)

        target_df = pd.DataFrame(target_detail)
        summary_rows.append({"id": target_id, **summarize_target(target_df)})
        print(
            f"SB{target_id:<2d}: shape={fit_shape:<9s} baseline clouds={len(baseline_weights):<4d} "
            f"weighted overlap mean={100.0 * summary_rows[-1]['weighted_overlap_fraction_mean']:.2f}%"
        )

    detail_df = pd.DataFrame(detail_rows)
    summary_df = pd.DataFrame(summary_rows).sort_values("id").reset_index(drop=True)

    total_recovered = float(detail_df["recovered_baseline_ray_hit_weight"].sum())
    total_baseline = float(detail_df["baseline_ray_hit_weight"].sum())
    wv = detail_df["weighted_overlap_fraction"].to_numpy(dtype=float)
    overall = {
        "global_ray_weighted_overlap_mean": float(total_recovered / max(total_baseline, 1.0)),
        "trial_mean_weighted_overlap": float(np.mean(wv)),
        "trial_median_weighted_overlap": float(np.median(wv)),
        "trial_p10_weighted_overlap": float(np.percentile(wv, 10.0)),
        "trial_p90_weighted_overlap": float(np.percentile(wv, 90.0)),
        "target_mean_weighted_overlap": float(summary_df["weighted_overlap_fraction_mean"].mean()),
        "runs_ge_0p90_weighted_fraction": float(np.mean(wv >= 0.90)),
        "runs_ge_0p80_weighted_fraction": float(np.mean(wv >= 0.80)),
        "exact_same_set_fraction": float(np.mean(detail_df["exact_same_set"].to_numpy(dtype=bool))),
        "n_targets": int(len(targets)),
        "n_trials_per_target": int(args.n_perturb),
        "n_trials_total": int(len(detail_df)),
        "random_seed": int(args.random_seed),
        "center_perturb_frac": float(args.center_perturb_frac),
        "axis_perturb_frac": float(args.axis_perturb_frac),
        "perturbation_distribution": "uniform_local_abc_box_center_plus_multiplicative_axes",
        "first_touch_engine": "canonical_v4" if args.use_canonical else "vectorized_identical_to_v4",
        "cloud_finding_method": "v4_full_sky_xy1p5_local_bubble_removed_cone_ray_radial_upper_limit",
        "baseline_definition": "v4 first-touch clouds from final-fit center & axes (unperturbed)",
        "weight_definition": "baseline selected-cloud ray_hit_count recovery",
        "geometry_source": final_table_path.as_posix(),
        "runtime_sec": float(time.perf_counter() - t0),
    }

    detail_csv = OUT_DIR / "6_initial_condition_sensitivity_trial_details.csv"
    summary_csv = OUT_DIR / "6_initial_condition_sensitivity_by_bubble.csv"
    overall_json = OUT_DIR / "6_initial_condition_sensitivity_summary.json"
    figure_path = FINAL_FIG_DIR / "6_initial_sensitivity_overlap.png"

    detail_df.to_csv(detail_csv, index=False, encoding="utf-8-sig")
    summary_df.to_csv(summary_csv, index=False, encoding="utf-8-sig")
    overall_json.write_text(json.dumps(overall, indent=2, ensure_ascii=False), encoding="utf-8")
    save_summary_figure(summary_df, overall, figure_path)

    print("Initial-condition sensitivity summary:")
    print(f"global_ray_weighted_overlap_mean = {overall['global_ray_weighted_overlap_mean']:.4f}")
    print(f"trial_mean_weighted_overlap      = {overall['trial_mean_weighted_overlap']:.4f}")
    print(f"trial_median_weighted_overlap    = {overall['trial_median_weighted_overlap']:.4f}")
    print(f"trial_p10_weighted_overlap       = {overall['trial_p10_weighted_overlap']:.4f}")
    print(f"runs_ge_0p90_weighted_fraction   = {overall['runs_ge_0p90_weighted_fraction']:.4f}")
    print(f"trial_details: {detail_csv}")
    print(f"by_bubble_summary: {summary_csv}")
    print(f"overall_statistics: {overall_json}")
    print(f"figure: {figure_path}")

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
