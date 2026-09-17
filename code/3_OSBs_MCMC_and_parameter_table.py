# -*- coding: utf-8 -*-
"""Script 3: OSBs MCMC and parameter table

Purpose
-------
Estimate MCMC uncertainties and assemble the final 30-object parameter table.

Method overview
---------------
1. Load the fitted or manually specified OSB geometry and the supporting cloud
   constraints.
2. Sample the local posterior around the adopted geometry to quantify parameter
   uncertainty.
3. Write the final parameter table and diagnostic figures without changing the adopted
   best-fit values.

Main inputs
-----------
- ../results/intermediate_output/1_first_contact_molecular_clouds/
- ../results/intermediate_output/2_automated_fit_results/json/SB{N}_fit.json
- ../data/OSB_initial_pre_identification_parameters.csv

Main outputs
------------
- ../results/superbubble_final_fit_parameters.csv: full fitted-geometry and MCMC table.
- ../results/figures/3_mcmc_analysis.png and corner plots: MCMC diagnostic figure products.

Figure/table role
-----------------
Produces the full fitted-geometry catalogue table and the MCMC diagnostic figure.

Runtime and data notes
----------------------
Reference runtime: 1 min 16.63 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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

import importlib.util
import json
import os
import sys
from pathlib import Path

import corner
import emcee
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks


# Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"

MID_DIR = Path("..") / "results" / "intermediate_output"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"

STAGE1_DIR = MID_DIR / "1_first_contact_molecular_clouds"
JSON_DIR = MID_DIR / "2_automated_fit_results" / "json"
CORNER_DIR = FINAL_FIG_DIR / "mcmc_corner_plots"
MCMC_SAMPLE_DIR = MID_DIR / "3_mcmc_samples"
FINAL_CSV = FINAL_DIR / "superbubble_final_fit_parameters.csv"
MCMC_ANALYSIS_FIG = FINAL_FIG_DIR / "3_mcmc_analysis.png"

PARAM_CSV = DATA_DIR / "OSB_initial_pre_identification_parameters.csv"
MCS_CSV = DATA_DIR / "MCs.csv"

MANUAL_IDS = [23, 26, 27, 37]

# Main workflow.
MCMC_WALKERS = 24
MCMC_STEPS = 320
MCMC_BURN_IN = 100
MCMC_THIN = 4
MCMC_MAX_SAVED = 3000


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def setup():
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    s2 = _load_module(HERE / "2_OSBs_3D_fitting.py", "pr_stage2")
    return s2


def sigma_from_signed_distance(signed_distance):
    values_pc = np.asarray(signed_distance, dtype=float) * 1000.0
    if values_pc.size < 3:
        return {
            "mean_pc": float("nan"),
            "sigma_pc": float("nan"),
            "std_pc": float("nan"),
            "median_pc": float("nan"),
            "mad_pc": float("nan"),
            "iqr_sigma_pc": float("nan"),
            "p68_abs_pc": float("nan"),
        }
    median = float(np.median(values_pc))
    mad = float(np.median(np.abs(values_pc - median)))
    mad_sigma = 1.4826 * mad
    keep = np.abs(values_pc - median) <= max(3.0 * 1.4826 * mad, 1e-6)
    cleaned = values_pc[keep] if keep.sum() >= 3 else values_pc
    q16, q84 = np.percentile(values_pc, [16.0, 84.0])
    q25, q75 = np.percentile(values_pc, [25.0, 75.0])
    return {
        "mean_pc": float(np.mean(cleaned)),
        "sigma_pc": float(mad_sigma),
        "std_pc": float(np.std(cleaned, ddof=1)) if cleaned.size > 1 else 0.0,
        "median_pc": median,
        "mad_pc": mad,
        "iqr_sigma_pc": float((q75 - q25) / 1.349),
        "p68_abs_pc": float(np.percentile(np.abs(values_pc), 68.0)),
        "central_68_half_width_pc": float(0.5 * (q84 - q16)),
        "n_sigma_kept": int(cleaned.size),
    }


def shell_signed_distance_from_params(auto, centered_points, params, fit_shape, c_fixed=0.1, circular_xy=False, fixed_angle=0.0):
    if fit_shape == "cylinder" and circular_xy:
        radius = float(np.exp(params[0]))
        axes = np.array([radius, radius, float(c_fixed)], dtype=float)
        point_radius = np.linalg.norm(centered_points[:, :2], axis=1)
        return point_radius - radius, axes, float(fixed_angle)

    if fit_shape == "cylinder":
        axes = np.array([np.exp(params[0]), np.exp(params[1]), float(c_fixed)], dtype=float)
        angle = float(params[2])
        local = auto.rotate_points_to_local(centered_points, angle)
        point_radius = np.linalg.norm(local[:, :2], axis=1)
        unit_xy = np.divide(
            local[:, :2],
            point_radius[:, None],
            out=np.zeros_like(local[:, :2]),
            where=point_radius[:, None] > 0,
        )
        shell_radius = 1.0 / np.sqrt(
            np.clip((unit_xy[:, 0] / axes[0]) ** 2 + (unit_xy[:, 1] / axes[1]) ** 2, 1e-12, None)
        )
        return point_radius - shell_radius, axes, angle

    if circular_xy:
        radius_xy, radius_z = np.exp(params[:2])
        axes = np.array([radius_xy, radius_xy, radius_z], dtype=float)
        angle = float(fixed_angle)
    else:
        axes = np.exp(params[:3])
        angle = float(params[3])
    local = auto.rotate_points_to_local(centered_points, angle)
    point_radius = np.linalg.norm(local, axis=1)
    unit_vectors = np.divide(local, point_radius[:, None], out=np.zeros_like(local), where=point_radius[:, None] > 0)
    shell_radius = 1.0 / np.sqrt(np.sum((unit_vectors / axes[None, :]) ** 2, axis=1))
    return point_radius - shell_radius, axes, angle


def angle_delta_deg(values, reference):
    return (np.asarray(values, dtype=float) - float(reference) + 90.0) % 180.0 - 90.0


def summarize_linear_samples(values):
    p16, p50, p84 = np.percentile(np.asarray(values, dtype=float), [16.0, 50.0, 84.0])
    return {
        "p16": float(p16),
        "median": float(p50),
        "p84": float(p84),
        "err_minus": float(p50 - p16),
        "err_plus": float(p84 - p50),
    }


def summarize_angle_samples(values, reference):
    deltas = angle_delta_deg(values, reference)
    p16, p50, p84 = np.percentile(deltas, [16.0, 50.0, 84.0])
    return {
        "p16": float((float(reference) + p16) % 180.0),
        "median": float((float(reference) + p50) % 180.0),
        "p84": float((float(reference) + p84) % 180.0),
        "err_minus": float(p50 - p16),
        "err_plus": float(p84 - p50),
    }


def summarize_undefined_angle():
    return {
        "p16": float("nan"),
        "median": float("nan"),
        "p84": float("nan"),
        "err_minus": float("nan"),
        "err_plus": float("nan"),
        "undefined": True,
    }


def detect_histogram_bimodality(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 200:
        return False, {}
    q01, q99 = np.percentile(values, [1.0, 99.0])
    if not np.isfinite(q01) or not np.isfinite(q99) or q99 <= q01:
        return False, {}
    bins = max(24, min(48, int(np.sqrt(values.size))))
    counts, edges = np.histogram(values, bins=bins, range=(q01, q99))
    smooth = gaussian_filter1d(counts.astype(float), sigma=1.0, mode="nearest")
    if np.max(smooth) <= 0:
        return False, {}
    peaks, props = find_peaks(
        smooth,
        prominence=max(1.0, 0.12 * np.max(smooth)),
        distance=max(2, bins // 8),
    )
    if len(peaks) < 2:
        return False, {"n_peaks": int(len(peaks))}
    order = np.argsort(smooth[peaks])[::-1]
    p1, p2 = np.sort(peaks[order[:2]])
    valley = float(np.min(smooth[p1:p2 + 1]))
    peak_low = float(min(smooth[p1], smooth[p2]))
    valley_ratio = valley / max(peak_low, 1e-9)
    center1 = 0.5 * (edges[p1] + edges[p1 + 1])
    center2 = 0.5 * (edges[p2] + edges[p2 + 1])
    separation = abs(center2 - center1) / max(q99 - q01, 1e-9)
    is_bimodal = valley_ratio < 0.72 and separation > 0.18
    return is_bimodal, {
        "n_peaks": int(len(peaks)),
        "valley_ratio": float(valley_ratio),
        "separation_fraction": float(separation),
    }


def detect_quantile_skew(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 100:
        return False, {}
    q05, q16, q50, q84, q95 = np.percentile(values, [5.0, 16.0, 50.0, 84.0, 95.0])
    lower = max(q50 - q16, 1e-9)
    upper = max(q84 - q50, 1e-9)
    asym_ratio = max(upper / lower, lower / upper)
    tail_ratio = max(q95 - q50, q50 - q05) / max(q84 - q16, 1e-9)
    is_skewed = asym_ratio > 2.2 and tail_ratio > 1.25
    return is_skewed, {
        "asym_ratio": float(asym_ratio),
        "tail_ratio": float(tail_ratio),
    }


def analyze_mcmc_posterior_shape(sample_df, fit_shape, circular_xy):
    diagnostics = {}
    problematic_parameters = []
    if fit_shape == "cylinder" and circular_xy:
        param_columns = ["r_xy_kpc"]
    elif circular_xy:
        param_columns = ["r_xy_kpc", "c_kpc"]
    elif fit_shape == "cylinder":
        param_columns = ["a_kpc", "b_kpc", "delta_angle_deg"]
    else:
        param_columns = ["a_kpc", "b_kpc", "c_kpc", "delta_angle_deg"]

    any_bimodal = False
    any_skewed = False
    for column in param_columns:
        if column not in sample_df.columns:
            continue
        values = sample_df[column].to_numpy(dtype=float)
        bimodal, bimodal_info = detect_histogram_bimodality(values)
        skewed, skew_info = detect_quantile_skew(values)
        diagnostics[column] = {
            "bimodal": bool(bimodal),
            "skewed": bool(skewed),
            **bimodal_info,
            **skew_info,
        }
        if bimodal or skewed:
            problematic_parameters.append(column)
        any_bimodal = any_bimodal or bimodal
        any_skewed = any_skewed or skewed

    return {
        "any_bimodal": bool(any_bimodal),
        "any_skewed": bool(any_skewed),
        "use_least_squares": bool(any_bimodal or any_skewed),
        "problematic_parameters": problematic_parameters,
        "diagnostics": diagnostics,
    }


def run_mcmc_shell_fit(
    auto,
    selected_clouds,
    seed,
    center_fit,
    axes_fit,
    angle_fit,
    sigma_info,
    weights,
    fit_shape,
    *,
    circular_xy,
    free_center,
    free_center_xy_fixed_z,
    free_center_z_only,
    center_bounds_scale,
    n_walkers,
    n_steps,
    burn_in,
    thin,
    random_seed,
    max_saved_samples,
    output_csv,
):
    points = selected_clouds[["x", "y", "z"]].to_numpy(dtype=float)
    weights = np.asarray(weights, dtype=float)
    weights = np.where(np.isfinite(weights) & (weights > 0), weights, 1.0)
    sigma_kpc = max(float(sigma_info.get("sigma_pc", np.nan)) / 1000.0, 1e-4)
    center_fit = np.asarray(center_fit, dtype=float)

    center_span = np.maximum(
        np.asarray(seed.axes_init, dtype=float) * float(center_bounds_scale),
        np.array([0.05, 0.05, 0.03], dtype=float),
    )
    xy_axes_init = np.asarray(seed.axes_init[:2], dtype=float)
    xy_max_extent = float(np.max(xy_axes_init)) if xy_axes_init.size else 0.1

    def center_inside_initial_xy_ellipse(center):
        if not (free_center or free_center_xy_fixed_z):
            return True
        offset = np.asarray(center, dtype=float) - np.asarray(seed.center, dtype=float)
        local_offset = auto.rotate_points_to_local(offset[None, :], seed.angle_deg)[0]
        denom_a = max(float(xy_axes_init[0]), 1e-9)
        denom_b = max(float(xy_axes_init[1]), 1e-9)
        q_xy = (float(local_offset[0]) / denom_a) ** 2 + (float(local_offset[1]) / denom_b) ** 2
        z_ok = True if free_center_xy_fixed_z else abs(float(center[2]) - float(seed.center[2])) <= float(center_span[2]) + 1e-9
        return bool(q_xy <= 1.0 + 1e-9 and z_ok)

    def center_z_in_bounds(z_value):
        return bool(abs(float(z_value) - float(seed.center[2])) <= float(center_span[2]) + 1e-9)

    if fit_shape == "cylinder" and circular_xy and free_center_xy_fixed_z:
        ndim = 3
        lower = np.array(
            [
                float(seed.center[0]) - xy_max_extent - 0.05,
                float(seed.center[1]) - xy_max_extent - 0.05,
                np.log(0.01),
            ],
            dtype=float,
        )
        upper = np.array(
            [
                float(seed.center[0]) + xy_max_extent + 0.05,
                float(seed.center[1]) + xy_max_extent + 0.05,
                np.log(5.0),
            ],
            dtype=float,
        )
        x0 = np.array([center_fit[0], center_fit[1], np.log(axes_fit[0])], dtype=float)
        sampled_parameter_names = ["center_x_kpc", "center_y_kpc", "log_r_xy"]
    elif fit_shape == "cylinder" and circular_xy:
        ndim = 1
        lower = np.array([np.log(0.01)], dtype=float)
        upper = np.array([np.log(5.0)], dtype=float)
        x0 = np.array([np.log(axes_fit[0])], dtype=float)
        sampled_parameter_names = ["log_r_xy"]
    elif fit_shape == "cylinder" and free_center_xy_fixed_z:
        ndim = 5
        lower = np.array(
            [
                float(seed.center[0]) - xy_max_extent - 0.05,
                float(seed.center[1]) - xy_max_extent - 0.05,
                np.log(0.01),
                np.log(0.01),
                -90.0,
            ],
            dtype=float,
        )
        upper = np.array(
            [
                float(seed.center[0]) + xy_max_extent + 0.05,
                float(seed.center[1]) + xy_max_extent + 0.05,
                np.log(5.0),
                np.log(5.0),
                90.0,
            ],
            dtype=float,
        )
        x0 = np.array([center_fit[0], center_fit[1], np.log(axes_fit[0]), np.log(axes_fit[1]), 0.0], dtype=float)
        sampled_parameter_names = ["center_x_kpc", "center_y_kpc", "log_a", "log_b", "delta_angle_deg"]
    elif fit_shape == "cylinder":
        ndim = 3
        lower = np.array([np.log(0.01), np.log(0.01), -90.0], dtype=float)
        upper = np.array([np.log(5.0), np.log(5.0), 90.0], dtype=float)
        x0 = np.array([np.log(axes_fit[0]), np.log(axes_fit[1]), 0.0], dtype=float)
        sampled_parameter_names = ["log_a", "log_b", "delta_angle_deg"]
    elif circular_xy and free_center:
        ndim = 5
        lower = np.array(
            [
                float(seed.center[0]) - xy_max_extent - 0.05,
                float(seed.center[1]) - xy_max_extent - 0.05,
                float(seed.center[2]) - float(center_span[2]),
                np.log(0.01),
                np.log(0.005),
            ],
            dtype=float,
        )
        upper = np.array(
            [
                float(seed.center[0]) + xy_max_extent + 0.05,
                float(seed.center[1]) + xy_max_extent + 0.05,
                float(seed.center[2]) + float(center_span[2]),
                np.log(5.0),
                np.log(2.0),
            ],
            dtype=float,
        )
        x0 = np.array([center_fit[0], center_fit[1], center_fit[2], np.log(axes_fit[0]), np.log(axes_fit[2])], dtype=float)
        sampled_parameter_names = ["center_x_kpc", "center_y_kpc", "center_z_kpc", "log_r_xy", "log_c"]
    elif circular_xy:
        ndim = 2
        lower = np.array([np.log(0.01), np.log(0.005)], dtype=float)
        upper = np.array([np.log(5.0), np.log(2.0)], dtype=float)
        x0 = np.array([np.log(axes_fit[0]), np.log(axes_fit[2])], dtype=float)
        sampled_parameter_names = ["log_r_xy", "log_c"]
    elif free_center:
        ndim = 7
        lower = np.array(
            [
                float(seed.center[0]) - xy_max_extent - 0.05,
                float(seed.center[1]) - xy_max_extent - 0.05,
                float(seed.center[2]) - float(center_span[2]),
                np.log(0.01),
                np.log(0.01),
                np.log(0.005),
                -90.0,
            ],
            dtype=float,
        )
        upper = np.array(
            [
                float(seed.center[0]) + xy_max_extent + 0.05,
                float(seed.center[1]) + xy_max_extent + 0.05,
                float(seed.center[2]) + float(center_span[2]),
                np.log(5.0),
                np.log(5.0),
                np.log(2.0),
                90.0,
            ],
            dtype=float,
        )
        x0 = np.array(
            [
                center_fit[0],
                center_fit[1],
                center_fit[2],
                np.log(axes_fit[0]),
                np.log(axes_fit[1]),
                np.log(axes_fit[2]),
                0.0,
            ],
            dtype=float,
        )
        sampled_parameter_names = ["center_x_kpc", "center_y_kpc", "center_z_kpc", "log_a", "log_b", "log_c", "delta_angle_deg"]
    elif free_center_z_only:
        ndim = 5
        lower = np.array(
            [
                float(seed.center[2]) - float(center_span[2]),
                np.log(0.01),
                np.log(0.01),
                np.log(0.005),
                -90.0,
            ],
            dtype=float,
        )
        upper = np.array(
            [
                float(seed.center[2]) + float(center_span[2]),
                np.log(5.0),
                np.log(5.0),
                np.log(2.0),
                90.0,
            ],
            dtype=float,
        )
        x0 = np.array(
            [
                center_fit[2],
                np.log(axes_fit[0]),
                np.log(axes_fit[1]),
                np.log(axes_fit[2]),
                0.0,
            ],
            dtype=float,
        )
        sampled_parameter_names = ["center_z_kpc", "log_a", "log_b", "log_c", "delta_angle_deg"]
    else:
        ndim = 4
        lower = np.array([np.log(0.01), np.log(0.01), np.log(0.005), -90.0], dtype=float)
        upper = np.array([np.log(5.0), np.log(5.0), np.log(2.0), 90.0], dtype=float)
        x0 = np.array([np.log(axes_fit[0]), np.log(axes_fit[1]), np.log(axes_fit[2]), 0.0], dtype=float)
        sampled_parameter_names = ["log_a", "log_b", "log_c", "delta_angle_deg"]

    n_walkers = max(int(n_walkers), 2 * ndim + 2)
    n_steps = int(n_steps)
    burn_in = min(max(int(burn_in), 0), max(n_steps - 1, 0))
    thin = max(int(thin), 1)
    rng = np.random.default_rng(int(random_seed))

    def log_prior(theta):
        theta = np.asarray(theta, dtype=float)
        if np.any(~np.isfinite(theta)) or np.any(theta < lower) or np.any(theta > upper):
            return -np.inf
        if free_center:
            if not center_inside_initial_xy_ellipse(theta[:3]):
                return -np.inf
        elif free_center_xy_fixed_z:
            if not center_inside_initial_xy_ellipse([theta[0], theta[1], float(seed.center[2])]):
                return -np.inf
        elif free_center_z_only:
            if not center_z_in_bounds(theta[0]):
                return -np.inf
        return 0.0

    def log_probability(theta):
        lp = log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        theta = np.asarray(theta, dtype=float)
        if free_center:
            center = np.array(theta[:3], dtype=float)
        elif free_center_xy_fixed_z:
            center = np.array([theta[0], theta[1], float(seed.center[2])], dtype=float)
        elif free_center_z_only:
            center = center_fit.copy()
            center[2] = float(theta[0])
        else:
            center = center_fit.copy()
        centered_points = points - center[None, :]
        model_theta = np.array(
            theta[3:] if free_center else (theta[2:] if free_center_xy_fixed_z else (theta[1:] if free_center_z_only else theta)),
            dtype=float,
            copy=True,
        )
        if not circular_xy:
            model_theta[-1] = (float(angle_fit) + model_theta[-1]) % 180.0
        signed_distance, _, _ = shell_signed_distance_from_params(
            auto,
            centered_points,
            model_theta,
            fit_shape,
            c_fixed=0.1,
            circular_xy=circular_xy,
            fixed_angle=angle_fit,
        )
        if np.any(~np.isfinite(signed_distance)):
            return -np.inf
        chi2 = np.sum(weights * (signed_distance / sigma_kpc) ** 2)
        log_norm = np.sum(np.log(2.0 * np.pi * (sigma_kpc ** 2) / weights))
        return float(lp - 0.5 * (chi2 + log_norm))

    walkers = []
    scale = np.full(ndim, 0.03, dtype=float)
    if not circular_xy:
        scale[-1] = 2.0
    attempts = 0
    while len(walkers) < n_walkers and attempts < n_walkers * 200:
        proposal = x0 + rng.normal(0.0, scale, size=ndim)
        if np.isfinite(log_prior(proposal)):
            walkers.append(proposal)
        attempts += 1
    if len(walkers) < n_walkers:
        raise RuntimeError("Failed to initialize MCMC walkers")

    sampler = emcee.EnsembleSampler(n_walkers, ndim, log_probability)
    sampler.run_mcmc(np.asarray(walkers), n_steps, progress=False)
    flat_params = sampler.get_chain(discard=burn_in, thin=thin, flat=True)
    flat_log_prob = sampler.get_log_prob(discard=burn_in, thin=thin, flat=True)

    if fit_shape == "cylinder" and circular_xy and free_center_xy_fixed_z:
        radius = np.exp(flat_params[:, 2])
        sample_df = pd.DataFrame(
            {
                "center_x_kpc": flat_params[:, 0],
                "center_y_kpc": flat_params[:, 1],
                "center_z_kpc": np.full(len(flat_params), float(seed.center[2]), dtype=float),
                "r_xy_kpc": radius,
                "a_kpc": radius,
                "b_kpc": radius,
                "c_kpc": np.full(len(flat_params), 0.1, dtype=float),
                "angle_deg": np.full(len(flat_params), float(angle_fit), dtype=float),
                "delta_angle_deg": np.full(len(flat_params), np.nan, dtype=float),
                "log_probability": flat_log_prob,
            }
        )
        r_summary = summarize_linear_samples(sample_df["r_xy_kpc"])
        summaries = {
            "center_x_kpc": summarize_linear_samples(sample_df["center_x_kpc"]),
            "center_y_kpc": summarize_linear_samples(sample_df["center_y_kpc"]),
            "center_z_kpc": {
                "p16": float(seed.center[2]),
                "median": float(seed.center[2]),
                "p84": float(seed.center[2]),
                "err_minus": 0.0,
                "err_plus": 0.0,
                "fixed": True,
            },
            "r_xy_kpc": r_summary,
            "a_kpc": r_summary,
            "b_kpc": r_summary,
            "c_kpc": {
                "p16": 0.1,
                "median": 0.1,
                "p84": 0.1,
                "err_minus": 0.0,
                "err_plus": 0.0,
                "fixed": True,
            },
            "angle_deg": summarize_undefined_angle(),
        }
    elif fit_shape == "cylinder" and circular_xy:
        radius = np.exp(flat_params[:, 0])
        sample_df = pd.DataFrame(
            {
                "r_xy_kpc": radius,
                "a_kpc": radius,
                "b_kpc": radius,
                "c_kpc": np.full(len(flat_params), 0.1, dtype=float),
                "angle_deg": np.full(len(flat_params), float(angle_fit), dtype=float),
                "delta_angle_deg": np.full(len(flat_params), np.nan, dtype=float),
                "log_probability": flat_log_prob,
            }
        )
        r_summary = summarize_linear_samples(sample_df["r_xy_kpc"])
        summaries = {
            "r_xy_kpc": r_summary,
            "a_kpc": r_summary,
            "b_kpc": r_summary,
            "c_kpc": {
                "p16": 0.1,
                "median": 0.1,
                "p84": 0.1,
                "err_minus": 0.0,
                "err_plus": 0.0,
                "fixed": True,
            },
            "angle_deg": summarize_undefined_angle(),
        }
    elif fit_shape == "cylinder" and free_center_xy_fixed_z:
        sample_df = pd.DataFrame(
            {
                "center_x_kpc": flat_params[:, 0],
                "center_y_kpc": flat_params[:, 1],
                "center_z_kpc": np.full(len(flat_params), float(seed.center[2]), dtype=float),
                "a_kpc": np.exp(flat_params[:, 2]),
                "b_kpc": np.exp(flat_params[:, 3]),
                "c_kpc": np.full(len(flat_params), 0.1, dtype=float),
                "angle_deg": (float(angle_fit) + flat_params[:, 4]) % 180.0,
                "delta_angle_deg": flat_params[:, 4],
                "log_probability": flat_log_prob,
            }
        )
        summaries = {
            "center_x_kpc": summarize_linear_samples(sample_df["center_x_kpc"]),
            "center_y_kpc": summarize_linear_samples(sample_df["center_y_kpc"]),
            "center_z_kpc": {
                "p16": float(seed.center[2]),
                "median": float(seed.center[2]),
                "p84": float(seed.center[2]),
                "err_minus": 0.0,
                "err_plus": 0.0,
                "fixed": True,
            },
            "a_kpc": summarize_linear_samples(sample_df["a_kpc"]),
            "b_kpc": summarize_linear_samples(sample_df["b_kpc"]),
            "c_kpc": {
                "p16": 0.1,
                "median": 0.1,
                "p84": 0.1,
                "err_minus": 0.0,
                "err_plus": 0.0,
                "fixed": True,
            },
            "angle_deg": summarize_angle_samples(sample_df["angle_deg"], angle_fit),
        }
    elif fit_shape == "cylinder":
        sample_df = pd.DataFrame(
            {
                "a_kpc": np.exp(flat_params[:, 0]),
                "b_kpc": np.exp(flat_params[:, 1]),
                "c_kpc": np.full(len(flat_params), 0.1, dtype=float),
                "angle_deg": (float(angle_fit) + flat_params[:, 2]) % 180.0,
                "delta_angle_deg": flat_params[:, 2],
                "log_probability": flat_log_prob,
            }
        )
        summaries = {
            "a_kpc": summarize_linear_samples(sample_df["a_kpc"]),
            "b_kpc": summarize_linear_samples(sample_df["b_kpc"]),
            "c_kpc": {
                "p16": 0.1,
                "median": 0.1,
                "p84": 0.1,
                "err_minus": 0.0,
                "err_plus": 0.0,
                "fixed": True,
            },
            "angle_deg": summarize_angle_samples(sample_df["angle_deg"], angle_fit),
        }
    elif circular_xy and free_center:
        radius = np.exp(flat_params[:, 3])
        sample_df = pd.DataFrame(
            {
                "center_x_kpc": flat_params[:, 0],
                "center_y_kpc": flat_params[:, 1],
                "center_z_kpc": flat_params[:, 2],
                "r_xy_kpc": radius,
                "a_kpc": radius,
                "b_kpc": radius,
                "c_kpc": np.exp(flat_params[:, 4]),
                "angle_deg": np.full(len(flat_params), float(angle_fit), dtype=float),
                "delta_angle_deg": np.full(len(flat_params), np.nan, dtype=float),
                "log_probability": flat_log_prob,
            }
        )
        r_summary = summarize_linear_samples(sample_df["r_xy_kpc"])
        summaries = {
            "center_x_kpc": summarize_linear_samples(sample_df["center_x_kpc"]),
            "center_y_kpc": summarize_linear_samples(sample_df["center_y_kpc"]),
            "center_z_kpc": summarize_linear_samples(sample_df["center_z_kpc"]),
            "r_xy_kpc": r_summary,
            "a_kpc": r_summary,
            "b_kpc": r_summary,
            "c_kpc": summarize_linear_samples(sample_df["c_kpc"]),
            "angle_deg": summarize_undefined_angle(),
        }
    elif free_center_z_only:
        sample_df = pd.DataFrame(
            {
                "center_x_kpc": np.full(len(flat_params), float(center_fit[0]), dtype=float),
                "center_y_kpc": np.full(len(flat_params), float(center_fit[1]), dtype=float),
                "center_z_kpc": flat_params[:, 0],
                "a_kpc": np.exp(flat_params[:, 1]),
                "b_kpc": np.exp(flat_params[:, 2]),
                "c_kpc": np.exp(flat_params[:, 3]),
                "angle_deg": (float(angle_fit) + flat_params[:, 4]) % 180.0,
                "delta_angle_deg": flat_params[:, 4],
                "log_probability": flat_log_prob,
            }
        )
        summaries = {
            "center_x_kpc": {
                "p16": float(center_fit[0]),
                "median": float(center_fit[0]),
                "p84": float(center_fit[0]),
                "err_minus": 0.0,
                "err_plus": 0.0,
                "fixed": True,
            },
            "center_y_kpc": {
                "p16": float(center_fit[1]),
                "median": float(center_fit[1]),
                "p84": float(center_fit[1]),
                "err_minus": 0.0,
                "err_plus": 0.0,
                "fixed": True,
            },
            "center_z_kpc": summarize_linear_samples(sample_df["center_z_kpc"]),
            "a_kpc": summarize_linear_samples(sample_df["a_kpc"]),
            "b_kpc": summarize_linear_samples(sample_df["b_kpc"]),
            "c_kpc": summarize_linear_samples(sample_df["c_kpc"]),
            "angle_deg": summarize_angle_samples(sample_df["angle_deg"], angle_fit),
        }
    elif circular_xy:
        radius = np.exp(flat_params[:, 0])
        sample_df = pd.DataFrame(
            {
                "r_xy_kpc": radius,
                "a_kpc": radius,
                "b_kpc": radius,
                "c_kpc": np.exp(flat_params[:, 1]),
                "angle_deg": np.full(len(flat_params), float(angle_fit), dtype=float),
                "delta_angle_deg": np.full(len(flat_params), np.nan, dtype=float),
                "log_probability": flat_log_prob,
            }
        )
        r_summary = summarize_linear_samples(sample_df["r_xy_kpc"])
        summaries = {
            "r_xy_kpc": r_summary,
            "a_kpc": r_summary,
            "b_kpc": r_summary,
            "c_kpc": summarize_linear_samples(sample_df["c_kpc"]),
            "angle_deg": summarize_undefined_angle(),
        }
    elif free_center:
        sample_df = pd.DataFrame(
            {
                "center_x_kpc": flat_params[:, 0],
                "center_y_kpc": flat_params[:, 1],
                "center_z_kpc": flat_params[:, 2],
                "a_kpc": np.exp(flat_params[:, 3]),
                "b_kpc": np.exp(flat_params[:, 4]),
                "c_kpc": np.exp(flat_params[:, 5]),
                "angle_deg": (float(angle_fit) + flat_params[:, 6]) % 180.0,
                "delta_angle_deg": flat_params[:, 6],
                "log_probability": flat_log_prob,
            }
        )
        summaries = {
            "center_x_kpc": summarize_linear_samples(sample_df["center_x_kpc"]),
            "center_y_kpc": summarize_linear_samples(sample_df["center_y_kpc"]),
            "center_z_kpc": summarize_linear_samples(sample_df["center_z_kpc"]),
            "a_kpc": summarize_linear_samples(sample_df["a_kpc"]),
            "b_kpc": summarize_linear_samples(sample_df["b_kpc"]),
            "c_kpc": summarize_linear_samples(sample_df["c_kpc"]),
            "angle_deg": summarize_angle_samples(sample_df["angle_deg"], angle_fit),
        }
    else:
        sample_df = pd.DataFrame(
            {
                "a_kpc": np.exp(flat_params[:, 0]),
                "b_kpc": np.exp(flat_params[:, 1]),
                "c_kpc": np.exp(flat_params[:, 2]),
                "angle_deg": (float(angle_fit) + flat_params[:, 3]) % 180.0,
                "delta_angle_deg": flat_params[:, 3],
                "log_probability": flat_log_prob,
            }
        )
        summaries = {
            "a_kpc": summarize_linear_samples(sample_df["a_kpc"]),
            "b_kpc": summarize_linear_samples(sample_df["b_kpc"]),
            "c_kpc": summarize_linear_samples(sample_df["c_kpc"]),
            "angle_deg": summarize_angle_samples(sample_df["angle_deg"], angle_fit),
        }

    saved_df = sample_df
    if max_saved_samples and len(sample_df) > int(max_saved_samples):
        saved_df = sample_df.sample(n=int(max_saved_samples), random_state=int(random_seed)).sort_index()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    saved_df.to_csv(output_csv, index=False)

    try:
        autocorr_time = sampler.get_autocorr_time(tol=0).tolist()
    except Exception:
        autocorr_time = None

    posterior_shape = analyze_mcmc_posterior_shape(sample_df, fit_shape, circular_xy)

    return {
        "enabled": True,
        "sampled_parameter_names": sampled_parameter_names,
        "circular_xy": bool(circular_xy),
        "angle_sampled": bool(not circular_xy),
        "n_dim": int(ndim),
        "n_walkers": int(n_walkers),
        "n_steps": int(n_steps),
        "burn_in": int(burn_in),
        "thin": int(thin),
        "n_flat_samples": int(len(flat_params)),
        "n_saved_samples": int(len(saved_df)),
        "random_seed": int(random_seed),
        "sigma_kpc": float(sigma_kpc),
        "acceptance_fraction_mean": float(np.mean(sampler.acceptance_fraction)),
        "acceptance_fraction_median": float(np.median(sampler.acceptance_fraction)),
        "autocorr_time_steps": autocorr_time,
        "posterior_shape": posterior_shape,
        "summaries": summaries,
        "samples_csv": str(output_csv),
    }


def add_mcmc_columns(result, mcmc_info):
    if not mcmc_info or not mcmc_info.get("enabled"):
        result["mcmc_enabled"] = False
        return
    result["mcmc_enabled"] = True
    result["mcmc_acceptance_fraction_mean"] = mcmc_info.get("acceptance_fraction_mean")
    result["mcmc_n_flat_samples"] = mcmc_info.get("n_flat_samples")
    posterior_shape = mcmc_info.get("posterior_shape", {})
    result["posterior_bimodal"] = bool(posterior_shape.get("any_bimodal", False))
    result["posterior_skewed"] = bool(posterior_shape.get("any_skewed", False))
    result["posterior_shape_problematic"] = bool(posterior_shape.get("use_least_squares", False))
    result["posterior_shape_problematic_parameters"] = ",".join(posterior_shape.get("problematic_parameters", []))
    for label in ["center_x", "center_y", "center_z"]:
        key = f"{label}_kpc"
        summary = mcmc_info["summaries"].get(key)
        if summary is None:
            continue
        prefix = f"mcmc_{label}"
        result[f"{prefix}_median_kpc"] = float(summary["median"])
        result[f"{prefix}_p16_kpc"] = float(summary["p16"])
        result[f"{prefix}_p84_kpc"] = float(summary["p84"])
        result[f"{prefix}_err_minus_kpc"] = float(summary["err_minus"])
        result[f"{prefix}_err_plus_kpc"] = float(summary["err_plus"])
    for label, unit in [("r_xy", "kpc"), ("a", "kpc"), ("b", "kpc"), ("c", "kpc"), ("angle", "deg")]:
        key = f"{label}_{unit}"
        summary = mcmc_info["summaries"].get(key)
        if summary is None:
            continue
        prefix = f"mcmc_{label}"
        result[f"{prefix}_median_{unit}"] = float(summary["median"])
        result[f"{prefix}_p16_{unit}"] = float(summary["p16"])
        result[f"{prefix}_p84_{unit}"] = float(summary["p84"])
        result[f"{prefix}_err_minus_{unit}"] = float(summary["err_minus"])
        result[f"{prefix}_err_plus_{unit}"] = float(summary["err_plus"])


def mcmc_median_geometry(mcmc_info, fallback_axes, fallback_angle, fit_shape, circular_xy):
    axes = np.asarray(fallback_axes, dtype=float).copy()
    summaries = mcmc_info.get("summaries", {}) if mcmc_info else {}
    if circular_xy:
        r_summary = summaries.get("r_xy_kpc") or summaries.get("a_kpc")
        if r_summary is not None:
            axes[0] = float(r_summary["median"])
            axes[1] = float(r_summary["median"])
        c_summary = summaries.get("c_kpc")
        if fit_shape == "cylinder":
            axes[2] = 0.1
        elif c_summary is not None:
            axes[2] = float(c_summary["median"])
        return axes, float(fallback_angle)

    for axis_index, label in enumerate(["a_kpc", "b_kpc", "c_kpc"]):
        summary = summaries.get(label)
        if summary is not None:
            axes[axis_index] = float(summary["median"])
    angle_summary = summaries.get("angle_deg")
    angle = float(fallback_angle)
    if angle_summary is not None and np.isfinite(float(angle_summary.get("median", float("nan")))):
        angle = float(angle_summary["median"])
    return axes, angle


def mcmc_angle_half_width_deg(mcmc_info):
    if not mcmc_info or not mcmc_info.get("enabled"):
        return float("nan")
    summary = mcmc_info.get("summaries", {}).get("angle_deg", {})
    err_minus = float(summary.get("err_minus", float("nan")))
    err_plus = float(summary.get("err_plus", float("nan")))
    if not np.isfinite(err_minus) or not np.isfinite(err_plus):
        return float("nan")
    return 0.5 * (err_minus + err_plus)


def save_mcmc_corner_plot(sample_csv, output_png, result, fit_shape, circular_xy):
    sample_csv = Path(sample_csv)
    if not sample_csv.exists():
        return None
    sample_df = pd.read_csv(sample_csv)
    if sample_df.empty:
        return None

    if fit_shape == "cylinder":
        if circular_xy:
            if {"center_x_kpc", "center_y_kpc"}.issubset(sample_df.columns):
                columns = ["center_x_kpc", "center_y_kpc", "r_xy_kpc"]
                labels = ["X (kpc)", "Y (kpc)", "r_xy (kpc)"]
                truths = [
                    float(result.get("center_kpc", [np.nan, np.nan, np.nan])[0]),
                    float(result.get("center_kpc", [np.nan, np.nan, np.nan])[1]),
                    float(result.get("fit_a_kpc", np.nan)),
                ]
            else:
                columns = ["r_xy_kpc"]
                labels = ["r_xy (kpc)"]
                truths = [float(result.get("fit_a_kpc", np.nan))]
        else:
            if {"center_x_kpc", "center_y_kpc"}.issubset(sample_df.columns):
                columns = ["center_x_kpc", "center_y_kpc", "a_kpc", "b_kpc", "angle_deg"]
                labels = ["X (kpc)", "Y (kpc)", "a (kpc)", "b (kpc)", "PA (deg)"]
                truths = [
                    float(result.get("center_kpc", [np.nan, np.nan, np.nan])[0]),
                    float(result.get("center_kpc", [np.nan, np.nan, np.nan])[1]),
                    float(result.get("fit_a_kpc", np.nan)),
                    float(result.get("fit_b_kpc", np.nan)),
                    float(result.get("angle_deg", np.nan)),
                ]
            else:
                columns = ["a_kpc", "b_kpc", "angle_deg"]
                labels = ["a (kpc)", "b (kpc)", "PA (deg)"]
                truths = [
                    float(result.get("fit_a_kpc", np.nan)),
                    float(result.get("fit_b_kpc", np.nan)),
                    float(result.get("angle_deg", np.nan)),
                ]
    else:
        if circular_xy:
            columns = ["center_x_kpc", "center_y_kpc", "center_z_kpc", "r_xy_kpc", "c_kpc"]
            labels = ["X (kpc)", "Y (kpc)", "Z (kpc)", "r_xy (kpc)", "c (kpc)"]
            truths = [
                float(result.get("center_kpc", [np.nan, np.nan, np.nan])[0]),
                float(result.get("center_kpc", [np.nan, np.nan, np.nan])[1]),
                float(result.get("center_kpc", [np.nan, np.nan, np.nan])[2]),
                float(result.get("fit_a_kpc", np.nan)),
                float(result.get("fit_c_kpc", np.nan)),
            ]
        else:
            columns = ["center_x_kpc", "center_y_kpc", "center_z_kpc", "a_kpc", "b_kpc", "c_kpc", "angle_deg"]
            labels = ["X (kpc)", "Y (kpc)", "Z (kpc)", "a (kpc)", "b (kpc)", "c (kpc)", "PA (deg)"]
            truths = [
                float(result.get("center_kpc", [np.nan, np.nan, np.nan])[0]),
                float(result.get("center_kpc", [np.nan, np.nan, np.nan])[1]),
                float(result.get("center_kpc", [np.nan, np.nan, np.nan])[2]),
                float(result.get("fit_a_kpc", np.nan)),
                float(result.get("fit_b_kpc", np.nan)),
                float(result.get("fit_c_kpc", np.nan)),
                float(result.get("angle_deg", np.nan)),
            ]
    if any(col not in sample_df.columns for col in columns):
        return None
    valid = []
    for col, label, truth in zip(columns, labels, truths):
        values = sample_df[col].to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            continue
        if np.nanmax(finite) - np.nanmin(finite) <= 1e-12:
            continue
        valid.append((col, label, truth))
    if len(valid) < 2:
        return None
    columns = [item[0] for item in valid]
    labels = [item[1] for item in valid]
    truths = [item[2] for item in valid]

    fig = corner.corner(
        sample_df[columns].to_numpy(dtype=float),
        labels=labels,
        truths=truths,
        truth_color="red",
        show_titles=True,
        title_fmt=".3f",
        label_kwargs={"fontsize": 11},
        title_kwargs={"fontsize": 10},
    )
    output_png = Path(output_png)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return str(output_png)



def mcmc_for_auto(sb, fit, core):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    result = fit["result"]
    samples_csv = MCMC_SAMPLE_DIR / f"SB{sb}_mcmc_samples.csv"
    corner_png = CORNER_DIR / f"SB{sb}_corner.png"

    mcmc_info = run_mcmc_shell_fit(
        core,
        fit["selected_clouds"],
        fit["seed"],
        fit["center_fit"],
        fit["axes_fit"],
        fit["angle_fit"],
        fit["sigma_info"],
        fit["weights"],
        fit["fit_shape"],
        circular_xy=fit["circular_xy"],
        free_center=fit["free_center"],
        free_center_xy_fixed_z=fit["free_center_xy_fixed_z"],
        free_center_z_only=False,
        center_bounds_scale=1.5,
        n_walkers=MCMC_WALKERS,
        n_steps=MCMC_STEPS,
        burn_in=MCMC_BURN_IN,
        thin=MCMC_THIN,
        random_seed=20260526 + int(sb),
        max_saved_samples=MCMC_MAX_SAVED,
        output_csv=samples_csv,
    )
    result["mcmc"] = mcmc_info
    result["mcmc_enabled"] = True
    result["outputs"]["mcmc_samples_csv"] = str(samples_csv)
    result["final_estimator_reason"] = "least_squares_with_mcmc_uncertainty"
    add_mcmc_columns(result, mcmc_info)

    saved = save_mcmc_corner_plot(samples_csv, corner_png, result, fit["fit_shape"], fit["circular_xy"])
    if saved:
        result["outputs"]["mcmc_corner_png"] = str(saved)


    Path(fit["json_path"]).write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    nd = mcmc_info.get("n_dim")
    print(f"SB{sb:<2d}: MCMC n_dim={nd} acc={mcmc_info.get('acceptance_fraction_mean', float('nan')):.2f} -> {corner_png.name}")


def _manual_geometry_mcmc(sb, center_fit, axes_fit, angle_fit, circular_xy, mark,
                          override_source, model_class, geometry_note,
                          core, all_clouds):
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    center_fit = np.asarray(center_fit, dtype=float)
    axes_fit = np.asarray(axes_fit, dtype=float)
    angle_fit = float(angle_fit)
    a_init, b_init, c_init = float(axes_fit[0]), float(axes_fit[1]), float(axes_fit[2])

    # Data-loading helper section.
    selected_clouds = pd.read_csv(STAGE1_DIR / f"SB{sb}_first_contact_molecular_clouds.csv")


    seed_for_counts = core.load_seed_v4(PARAM_CSV, sb)
    filtered_clouds = core.filter_local_bubble_clouds(all_clouds)
    candidate_clouds = core.preselect_clouds_v4(filtered_clouds, seed_for_counts, 1.5)
    catalog_cloud_count = int(len(all_clouds))
    local_bubble_removed_cloud_count = int(len(all_clouds) - len(filtered_clouds))
    candidate_cloud_count = int(len(candidate_clouds))

    # Geometry and shell-mask conventions used by this analysis stage.
    centered = selected_clouds[["x", "y", "z"]].to_numpy(dtype=float) - center_fit[None, :]
    if circular_xy:
        params = np.array([np.log(a_init), np.log(c_init)], dtype=float)
        fixed_angle = angle_fit
    else:
        params = np.array([np.log(a_init), np.log(b_init), np.log(c_init), 0.0], dtype=float)
        fixed_angle = 0.0
    signed_distance, _, _ = shell_signed_distance_from_params(
        core, centered, params, "ellipsoid", c_fixed=c_init, circular_xy=circular_xy, fixed_angle=fixed_angle,
    )
    if not circular_xy:
        local = core.rotate_points_to_local(centered, angle_fit)
        point_radius = np.linalg.norm(local, axis=1)
        unit_vectors = np.divide(local, point_radius[:, None], out=np.zeros_like(local), where=point_radius[:, None] > 0)
        shell_radius = 1.0 / np.sqrt(np.sum((unit_vectors / axes_fit[None, :]) ** 2, axis=1))
        signed_distance = point_radius - shell_radius

    sigma_info = sigma_from_signed_distance(signed_distance)
    weights = np.ones(len(selected_clouds), dtype=float)

    # Geometry and shell-mask conventions used by this analysis stage.
    seed = core.load_seed_v4(PARAM_CSV, sb)
    if hasattr(seed, "_replace"):
        seed = seed._replace(center=center_fit, axes_init=axes_fit, angle_deg=angle_fit)
    else:
        seed.center = center_fit
        seed.axes_init = axes_fit
        seed.angle_deg = angle_fit

    samples_csv = MCMC_SAMPLE_DIR / f"SB{sb}_mcmc_samples.csv"
    corner_png = CORNER_DIR / f"SB{sb}_corner.png"
    mcmc_info = run_mcmc_shell_fit(
        core, selected_clouds, seed, center_fit, axes_fit, angle_fit, sigma_info, weights, "ellipsoid",
        circular_xy=circular_xy, free_center=True, free_center_xy_fixed_z=False, free_center_z_only=False,
        center_bounds_scale=1.5, n_walkers=MCMC_WALKERS, n_steps=MCMC_STEPS, burn_in=MCMC_BURN_IN,
        thin=MCMC_THIN, random_seed=20260601 + int(sb), max_saved_samples=MCMC_MAX_SAVED, output_csv=samples_csv,
    )

    clouds_csv = JSON_DIR / f"SB{sb}_clouds.csv"
    json_path = JSON_DIR / f"SB{sb}_fit.json"
    selected_out = selected_clouds.copy()
    selected_out["signed_distance_pc"] = signed_distance * 1000.0
    selected_out["weight"] = weights
    selected_out.to_csv(clouds_csv, index=False)

    result = {
        "target_index": int(sb),
        "version": "peer_review_v4",
        "manual_override": True,
        "manual_override_source": override_source,
        "model_class": model_class,
        "geometry_note": geometry_note,
        "fit_shape": "ellipsoid",
        "xy_model": "circular" if circular_xy else "elliptical",
        "circular_xy": circular_xy,
        "potential_circular_xy_by_axis_ratio": False,
        "final_parameter_estimator": "manual_override",
        "final_estimator_reason": "manual_override_with_mcmc_uncertainty",
        "least_squares_axes_kpc": axes_fit.tolist(),
        "least_squares_angle_deg": float(angle_fit),
        "least_squares_center_kpc": center_fit.tolist(),
        "initial_free_ab_axis_ratio": float(a_init / b_init) if b_init else float("nan"),
        "fit_ab_axis_ratio": float(a_init / b_init) if b_init else float("nan"),
        "center_kpc": center_fit.tolist(),
        "initial_center_kpc": center_fit.tolist(),
        "angle_deg": float(angle_fit),
        "initial_angle_deg": float(angle_fit),
        "mark": int(mark),
        "hemisphere": "full",
        "initial_axes_kpc": axes_fit.tolist(),
        "fit_axes_kpc": axes_fit.tolist(),
        "fit_a_kpc": float(a_init),
        "fit_b_kpc": float(b_init),
        "fit_c_kpc": float(c_init),
        "fit_sigma_pc": float(sigma_info["sigma_pc"]),
        "fit_std_pc": float(sigma_info["std_pc"]),
        "fit_mean_pc": float(sigma_info["mean_pc"]),
        "fit_median_pc": float(sigma_info["median_pc"]),
        "fit_mad_pc": float(sigma_info["mad_pc"]),
        "fit_iqr_sigma_pc": float(sigma_info["iqr_sigma_pc"]),
        "fit_p68_abs_pc": float(sigma_info["p68_abs_pc"]),
        "fit_central_68_half_width_pc": float(sigma_info["central_68_half_width_pc"]),
        "C_local": float("nan"),
        "C_local_pairs": 0,
        "catalog_cloud_count": catalog_cloud_count,
        "local_bubble_filtered_cloud_count": int(len(filtered_clouds)),
        "local_bubble_removed_cloud_count": local_bubble_removed_cloud_count,
        "candidate_cloud_count": candidate_cloud_count,
        "selected_cloud_count": int(len(selected_clouds)),
        "ray_count": 0,
        "ray_hit_count": 0,
        "hit_ray_count": 0,
        "parameters": {
            "catalog_filter": "Flag != 1 (remove Local Bubble boundary clouds)",
            "xy_search_scale": 1.5,
            "fit_center_bounds_scale": 1.5,
            "search_scale": 1.5,
            "max_hits_per_ray": 1,
            "radius_tolerance_pc": 0.0,
            "weight_mode": "uniform",
            "circular_axis_ratio_threshold": 1.30,
        },
        "mcmc": mcmc_info,
        "mcmc_enabled": True,
        "outputs": {
            "selected_clouds_csv": str(clouds_csv),
            "json": str(json_path),
            "mcmc_samples_csv": str(samples_csv),
        },
    }
    add_mcmc_columns(result, mcmc_info)
    saved = save_mcmc_corner_plot(samples_csv, corner_png, result, "ellipsoid", circular_xy)
    if saved:
        result["outputs"]["mcmc_corner_png"] = str(saved)

    json_path.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"SB{sb:<2d}: {model_class} + MCMC (xy={'circular' if circular_xy else 'elliptical'}, "
          f"std={sigma_info['sigma_pc']:.2f} pc, n_dim={mcmc_info.get('n_dim')}) -> {corner_png.name}")
    return {"id": sb, "status": "ok", "json_path": str(json_path)}


def manual_override(sb, manual_row, core, all_clouds):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    cx, cy, cz = float(manual_row["center_x_kpc"]), float(manual_row["center_y_kpc"]), float(manual_row["center_z_kpc"])
    a_init, b_init, c_init = float(manual_row["a_init_kpc"]), float(manual_row["b_init_kpc"]), float(manual_row["c_init_kpc"])
    angle_init = float(manual_row["angle_deg"])
    center_fit = np.array([cx, cy, cz], dtype=float)
    axes_fit = np.array([a_init, b_init, c_init], dtype=float)
    circular_xy = bool(abs(a_init - b_init) <= 1e-9)
    return _manual_geometry_mcmc(
        sb, center_fit, axes_fit, angle_init, circular_xy, int(manual_row.get("mark", 2)),
        override_source=str(PARAM_CSV),
        model_class="half_ellipsoid_manual",
        geometry_note=(
            "user-specified half-ellipsoid (cut at z=center_z); fitting and MCMC use the "
            "full-ellipsoid math anchored at the manual values; cavity/shell defined on the cap."
        ),
        core=core, all_clouds=all_clouds,
    )


# ============================================================================
# Geometry and shell-mask conventions used by this analysis stage.
# ============================================================================
MCMC_FIG_COL_X = "#4c78a8"
MCMC_FIG_COL_Y = "#f28e2b"
MCMC_FIG_COL_Z = "#59a14f"
MCMC_FIG_COL_PA = "#b07aa1"


def _mcmc_fig_sym_err(minus, plus):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    if not np.isfinite(minus) or not np.isfinite(plus):
        return np.nan
    val = 0.5 * (abs(minus) + abs(plus))
    return val if val > 0 else np.nan


def _mcmc_fig_relative(value, minus, plus):
    half = _mcmc_fig_sym_err(minus, plus)
    if not np.isfinite(value) or value == 0 or not np.isfinite(half):
        return np.nan
    return half / abs(value)


def _mcmc_fig_weighted_median(values, weights):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[good], weights[good]
    if values.size == 0:
        return np.nan
    order = np.argsort(values)
    v, w = values[order], weights[order]
    cw = np.cumsum(w)
    return float(v[np.searchsorted(cw, 0.5 * cw[-1])])


def render_mcmc_analysis_figure(final_df, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    from matplotlib.ticker import FuncFormatter, LogLocator

    rcParams["axes.unicode_minus"] = False
    rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial"]
    rcParams.update({
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
    })
    base_fs = 10

    def _plain_log_label(v, _pos=None):
        if v <= 0:
            return ""
        if v >= 1:
            return f"{v:.0f}"
        if v >= 0.1:
            return f"{v:.1f}"
        if v >= 0.01:
            return f"{v:.2f}"
        return f"{v:.3f}"

    def _format_log_axis(ax):
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=15))
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=(2.0, 3.0, 5.0), numticks=15))
        ax.yaxis.set_major_formatter(FuncFormatter(_plain_log_label))
        ax.yaxis.set_minor_formatter(FuncFormatter(_plain_log_label))
        ax.tick_params(axis="y", which="major", labelsize=base_fs)
        ax.tick_params(axis="y", which="minor", labelsize=base_fs)
        ax.grid(True, which="both", alpha=0.3)

    def _plot_skip_nan(ax, x, y, **kwargs):
        y = np.asarray(y, dtype=float)
        mask = np.isfinite(y)
        ax.plot(np.asarray(x)[mask], y[mask], **kwargs)

    def _boxplot_with_points(ax, data_list, labels, colors, ylabel, title):
        clean = [np.asarray([v for v in arr if np.isfinite(v)], dtype=float) for arr in data_list]
        positions = np.arange(1, len(clean) + 1)
        bp = ax.boxplot(clean, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.35)
            patch.set_edgecolor("black")
        for med in bp["medians"]:
            med.set_color("black")
            med.set_linewidth(1.6)
        rng = np.random.default_rng(20260528)
        for pos, arr, color in zip(positions, clean, colors):
            if arr.size == 0:
                continue
            jitter = rng.uniform(-0.16, 0.16, size=arr.size)
            ax.scatter(pos + jitter, arr, s=22, color=color, edgecolor="black", linewidth=0.3, alpha=0.8, zorder=3)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=base_fs)
        ax.set_ylabel(ylabel)
        ax.set_title(title, pad=6)
        _format_log_axis(ax)

    df = final_df.sort_values("id").reset_index(drop=True)
    ids = df["id"].astype(int).tolist()
    x_idx = np.arange(len(ids))
    x_labels = [str(sb) if i % 2 == 0 else "" for i, sb in enumerate(ids)]

    cen_x = np.array([1000.0 * _mcmc_fig_sym_err(r["mcmc_center_x_err_minus_kpc"], r["mcmc_center_x_err_plus_kpc"]) for _, r in df.iterrows()])
    cen_y = np.array([1000.0 * _mcmc_fig_sym_err(r["mcmc_center_y_err_minus_kpc"], r["mcmc_center_y_err_plus_kpc"]) for _, r in df.iterrows()])
    cen_z = np.array([1000.0 * _mcmc_fig_sym_err(r["mcmc_center_z_err_minus_kpc"], r["mcmc_center_z_err_plus_kpc"]) for _, r in df.iterrows()])

    rel_a = np.array([_mcmc_fig_relative(r["a_radius_kpc"], r["mcmc_a_err_minus_kpc"], r["mcmc_a_err_plus_kpc"]) for _, r in df.iterrows()])
    rel_b = np.array([_mcmc_fig_relative(r["b_radius_kpc"], r["mcmc_b_err_minus_kpc"], r["mcmc_b_err_plus_kpc"]) for _, r in df.iterrows()])
    rel_c = np.array([_mcmc_fig_relative(r["c_radius_kpc"], r["mcmc_c_err_minus_kpc"], r["mcmc_c_err_plus_kpc"]) for _, r in df.iterrows()])

    pa_val = df["mcmc_angle_median_deg"].to_numpy(dtype=float)
    pa_lo = df["mcmc_angle_err_minus_deg"].to_numpy(dtype=float)
    pa_hi = df["mcmc_angle_err_plus_deg"].to_numpy(dtype=float)
    pa_sym = np.array([_mcmc_fig_sym_err(lo, hi) for lo, hi in zip(pa_lo, pa_hi)])
    has_pa = np.isfinite(pa_val)
    pa_clean = pa_val[has_pa]
    sigma_pa_clean = pa_sym[has_pa]
    pa_weights = 1.0 / np.maximum(sigma_pa_clean, 1e-6) ** 2
    pa_wmedian = _mcmc_fig_weighted_median(pa_clean, pa_weights)
    pa_plain_median = float(np.nanmedian(pa_clean)) if pa_clean.size else np.nan

    # Direct A4 output avoids later down-scaling that would make text smaller than body text.
    fig, axes = plt.subplots(3, 2, figsize=(8.27, 10.8), constrained_layout=True)

    ax = axes[0, 0]
    _plot_skip_nan(ax, x_idx, cen_x, marker="o", ls="-", color=MCMC_FIG_COL_X, ms=4.5, label=r"$x$")
    _plot_skip_nan(ax, x_idx, cen_y, marker="s", ls="-", color=MCMC_FIG_COL_Y, ms=4.5, label=r"$y$")
    _plot_skip_nan(ax, x_idx, cen_z, marker="^", ls="-", color=MCMC_FIG_COL_Z, ms=4.5, label=r"$z$")
    ax.set_xticks(x_idx)
    ax.set_xticklabels(x_labels, rotation=0, fontsize=base_fs)
    ax.set_xlabel("SB ID")
    ax.set_ylabel("Center error [pc]")
    ax.set_title("Center errors", pad=6)
    _format_log_axis(ax)
    ax.legend(ncol=3, loc="upper left", handlelength=1.4, columnspacing=0.8)
    ax.text(0.02, 0.02, "z fixed for cylinders", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=base_fs)

    _boxplot_with_points(axes[0, 1], [cen_x, cen_y, cen_z], ["x", "y", "z"],
                         [MCMC_FIG_COL_X, MCMC_FIG_COL_Y, MCMC_FIG_COL_Z],
                         "Center error [pc]", "Center-error distribution")

    ax = axes[1, 0]
    _plot_skip_nan(ax, x_idx, rel_a, marker="o", ls="-", color=MCMC_FIG_COL_X, ms=4.5, label=r"$a$")
    _plot_skip_nan(ax, x_idx, rel_b, marker="s", ls="-", color=MCMC_FIG_COL_Y, ms=4.5, label=r"$b$")
    _plot_skip_nan(ax, x_idx, rel_c, marker="^", ls="-", color=MCMC_FIG_COL_Z, ms=4.5, label=r"$c$")
    ax.set_xticks(x_idx)
    ax.set_xticklabels(x_labels, rotation=0, fontsize=base_fs)
    ax.set_xlabel("SB ID")
    ax.set_ylabel(r"Relative axis error")
    ax.set_title("Axis errors", pad=6)
    _format_log_axis(ax)
    ax.legend(ncol=3, loc="upper left", handlelength=1.4, columnspacing=0.8)
    ax.text(0.02, 0.02, "c shown for ellipsoids", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=base_fs)

    _boxplot_with_points(axes[1, 1], [rel_a, rel_b, rel_c], ["a", "b", "c"],
                         [MCMC_FIG_COL_X, MCMC_FIG_COL_Y, MCMC_FIG_COL_Z],
                         r"Relative axis error", "Axis-error distribution")

    ax = axes[2, 0]
    idx_pa = x_idx[has_pa]
    ax.errorbar(idx_pa, pa_val[has_pa], yerr=[np.abs(pa_lo[has_pa]), np.abs(pa_hi[has_pa])],
                fmt="o", color=MCMC_FIG_COL_PA, ecolor="#c9b3c6", ms=4.5, capsize=3)
    ax.set_xticks(x_idx)
    ax.set_xticklabels(x_labels, rotation=90, fontsize=base_fs)
    ax.set_xlabel("SB ID")
    if np.isfinite(pa_wmedian):
        ax.axhline(pa_wmedian, color="#d62728", ls="-", lw=1.8, zorder=1,
                   label=f"Weighted median = {pa_wmedian:.1f} deg")
    ax.set_ylabel("Position angle PA [deg]")
    ax.set_title(f"PA by target ({int(has_pa.sum())} targets)", pad=6)
    ax.set_ylim(-5, 185)
    ax.set_yticks([0, 45, 90, 135, 180])
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="lower right", handlelength=1.4)

    ax = axes[2, 1]
    ax.hist(pa_clean, bins=np.arange(0, 181, 20), color=MCMC_FIG_COL_PA, alpha=0.55, edgecolor="black")
    ax.set_xlabel("Position angle PA [deg]")
    ax.set_ylabel("Number of targets")
    ax.set_title("PA distribution", pad=6)
    ax.set_xlim(0, 180)
    ax.set_xticks([0, 45, 90, 135, 180])
    ax.grid(True, axis="y", alpha=0.3)
    if np.isfinite(pa_wmedian):
        ax.axvline(pa_wmedian, color="#d62728", ls="-", lw=2.0,
                   label=f"Weighted median = {pa_wmedian:.1f} deg")
    if np.isfinite(pa_plain_median):
        ax.axvline(pa_plain_median, color="black", ls="--", lw=1.2,
                   label=f"Median = {pa_plain_median:.1f} deg")
    ax.legend(loc="upper right", handlelength=1.4)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main():
    """Run Script 3 from validated inputs to the documented outputs."""
    CORNER_DIR.mkdir(parents=True, exist_ok=True)
    MCMC_SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_CSV.parent.mkdir(parents=True, exist_ok=True)
    s2 = setup()

    formal = s2.read_formal_parameter_table(PARAM_CSV)
    all_ids = formal["index"].astype(int).tolist()
    auto_ids = [sb for sb in all_ids if sb not in MANUAL_IDS]
    all_clouds = s2.parse_mcs_catalog(MCS_CSV)

    final_rows = []


    print(f"Running MCMC fits for {len(auto_ids)} automatically fitted OSBs.")
    for sb in auto_ids:
        seed = s2.load_seed_v4(PARAM_CSV, sb)
        selected_clouds = pd.read_csv(STAGE1_DIR / f"SB{sb}_first_contact_molecular_clouds.csv")
        rays_csv = STAGE1_DIR / f"SB{sb}_ray_hits.csv"
        ray_hits = pd.read_csv(rays_csv) if rays_csv.exists() else pd.DataFrame(columns=["ray_index", "hit_rank"])
        fit = s2.fit_one_target(sb, seed, selected_clouds, ray_hits, all_clouds)
        mcmc_for_auto(sb, fit, s2)
        final_rows.append({"id": sb, "status": "ok", "json_path": str(fit["json_path"])})

    # Figure-panel rendering section.
    # Geometry and shell-mask conventions used by this analysis stage.
    print(f"Applying manual geometry overrides for {len(MANUAL_IDS)} OSBs.")
    manual = formal.set_index("id")
    for sb in MANUAL_IDS:
        if sb not in manual.index:
            print(f"Manual geometry override skipped for SB{sb}: not present in the formal table.")
            continue
        row = manual_override(sb, manual.loc[sb], s2, all_clouds)
        final_rows.append(row)


    print("Building final fitted-parameter table.")
    final_rows_sorted = sorted(final_rows, key=lambda r: int(r["id"]))
    final = s2.build_best_fit_table(final_rows_sorted)
    final.to_csv(FINAL_CSV, index=False, encoding="utf-8-sig")
    print(f"Final fitted-parameter table: {FINAL_CSV}")
    n_corner = len(list(CORNER_DIR.glob('*.png')))
    print(f"MCMC corner plots written: {n_corner}")


    print("Rendering MCMC analysis summary figure.")
    render_mcmc_analysis_figure(final, MCMC_ANALYSIS_FIG)
    print(f"MCMC analysis figure: {MCMC_ANALYSIS_FIG}")

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
