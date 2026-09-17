# -*- coding: utf-8 -*-
"""Script 2: OSBs 3D fitting

Purpose
-------
Fit ellipsoidal-cap or elliptical-cylinder shell geometry to first-contact clouds.

Method overview
---------------
1. Load the first-contact clouds from Script 1 together with the seed geometry and cloud
   catalogue.
2. Fit the configured shell model for each target with fit_one_target(); the fitting
   helpers handle ellipsoidal caps and elliptical cylinders.
3. Write per-target fit JSON/cloud products and assemble the automated parameter table
   for MCMC and plotting.

Main inputs
-----------
- ../results/intermediate_output/1_first_contact_molecular_clouds/SB{N}_first_contact_molecular_clouds.csv
- ../data/OSB_initial_pre_identification_parameters.csv
- ../data/MCs.csv

Main outputs
------------
- ../results/intermediate_output/2_automated_fit_results/; input to scripts 3 and 4.

Figure/table role
-----------------
../results/intermediate_output/2_automated_fit_results/; input to scripts 3 and 4.

Runtime and data notes
----------------------
Reference runtime: 0.91 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from scipy.optimize import least_squares


# Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"

MID_DIR = Path("..") / "results" / "intermediate_output"

STAGE1_DIR = MID_DIR / "1_first_contact_molecular_clouds"
FIT_DIR = MID_DIR / "2_automated_fit_results"
JSON_DIR = FIT_DIR / "json"
BESTFIT_CSV = FIT_DIR / "automated_superbubble_fit_parameters.csv"

PARAM_CSV = DATA_DIR / "OSB_initial_pre_identification_parameters.csv"
MCS_CSV = DATA_DIR / "MCs.csv"

# Figure styling and layout configuration for reproducible rendering.
MANUAL_IDS = [23, 26, 27, 37]
SKIP_AUTO_IDS = MANUAL_IDS


def package_relative(path: Path) -> str:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    path = Path(path)
    return path.relative_to(HERE).as_posix() if path.is_absolute() else path.as_posix()


def package_path(path_value) -> Path:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    path = Path(path_value)
    return path


WEIGHT_MODE = "uniform"
FIT_CENTER_BOUNDS_SCALE = 1.5
CIRCULAR_AXIS_RATIO_THRESHOLD = 1.30
DISABLE_CIRCULAR_XY_FALLBACK = False
ASYM_OFFSET_OVER_A = 0.45
ASYM_MIN_CLOUDS = 11
ASYM_MAX_AB_RATIO = 3.0
N_AZIMUTH = 72
N_ELEVATION = 18


# ============================================================================
# ----------------------------------------------------------------------------
# Geometry and shell-mask conventions used by this analysis stage.
# ============================================================================
# ----------------------------------------------------------------------------
# Geometry and shell-mask conventions used by this analysis stage.
# ----------------------------------------------------------------------------


def rotation_matrix_z(angle_deg):
    theta = np.deg2rad(float(angle_deg))
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])



def rotate_points_to_local(points, angle_deg):
    rot = rotation_matrix_z(-angle_deg)
    return np.asarray(points) @ rot.T



def rotate_points_to_global(points, angle_deg):
    rot = rotation_matrix_z(angle_deg)
    return np.asarray(points) @ rot.T



@dataclass
class BubbleSeed:
    target_index: int
    center: np.ndarray
    axes_init: np.ndarray
    angle_deg: float
    mark: int
    name: str



def shell_distances(points_global, seed, fit_axes, fit_shape):
    local_points = rotate_points_to_local(points_global - seed.center[None, :], seed.angle_deg)
    if fit_shape == "cylinder":
        point_radius = np.linalg.norm(local_points[:, :2], axis=1)
        unit_xy = np.divide(local_points[:, :2], point_radius[:, None], out=np.zeros_like(local_points[:, :2]), where=point_radius[:, None] > 0)
        shell_radius = 1.0 / np.sqrt(
            np.clip((unit_xy[:, 0] / fit_axes[0]) ** 2 + (unit_xy[:, 1] / fit_axes[1]) ** 2, 1e-12, None)
        )
        signed_distance = point_radius - shell_radius
        shell_points_local = np.column_stack(
            [
                unit_xy[:, 0] * shell_radius,
                unit_xy[:, 1] * shell_radius,
                np.clip(local_points[:, 2], -fit_axes[2], fit_axes[2]),
            ]
        )
        unit_vectors = np.zeros_like(local_points)
        unit_vectors[:, :2] = unit_xy
    else:
        point_radius = np.linalg.norm(local_points, axis=1)
        unit_vectors = np.divide(local_points, point_radius[:, None], out=np.zeros_like(local_points), where=point_radius[:, None] > 0)
        shell_radius = 1.0 / np.sqrt(np.sum((unit_vectors / fit_axes[None, :]) ** 2, axis=1))
        signed_distance = point_radius - shell_radius
        shell_points_local = unit_vectors * shell_radius[:, None]
    shell_points_global = rotate_points_to_global(shell_points_local, seed.angle_deg) + seed.center[None, :]
    return signed_distance, shell_points_global, unit_vectors, shell_points_local


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


def parse_mcs_catalog(path):
    df = pd.read_csv(path)
    rename = {
        "Rad": "Rad_pc",
        "Area": "Area_deg2",
        "Mass": "Mass_Msun",
        "Sigma": "Sigma_Msun_pc2",
        "rho": "rho_mag_kpc",
    }
    df = df.rename(columns=rename)
    required = ["Seq", "GLON", "GLAT", "Dist", "Rad_pc", "Area_deg2", "Mass_Msun", "Sigma_Msun_pc2", "rho_mag_kpc", "Flag"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"molecular-cloud catalogue is missing required columns: {', '.join(missing)}")
    for column in required:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df = df.dropna(subset=required).copy()
    if df.empty:
        raise ValueError('Invalid input or missing required data.')
    df["Seq"] = df["Seq"].astype(int)
    df["Flag"] = df["Flag"].astype(int)
    coords = SkyCoord(l=df["GLON"].to_numpy() * u.deg, b=df["GLAT"].to_numpy() * u.deg, distance=df["Dist"].to_numpy() * u.kpc, frame="galactic")
    cart = coords.cartesian
    df["x"] = cart.x.to_value(u.kpc)
    df["y"] = cart.y.to_value(u.kpc)
    df["z"] = cart.z.to_value(u.kpc)
    df["rad_kpc"] = df["Rad_pc"] / 1000.0
    return df



def fit_shape_for_seed(seed):
    return "cylinder" if int(seed.mark) == 3 else "ellipsoid"



def ab_axis_ratio(axes):
    axes = np.asarray(axes, dtype=float)
    smaller = max(float(np.min(axes[:2])), 1e-9)
    return float(np.max(axes[:2]) / smaller)



def should_use_circular_xy(axes, threshold):
    return ab_axis_ratio(axes) <= float(threshold)



def fit_cloud_shell(
        selected_clouds,
    seed,
    weight_mode,
    fit_shape,
    circular_xy=False,
    free_center=False,
    free_center_xy_fixed_z=False,
    free_center_z_only=False,
    center_bounds_scale=1.5,
):
    points = selected_clouds[["x", "y", "z"]].to_numpy(dtype=float)
    if len(points) < 5:
        raise RuntimeError('Runtime failure in this pipeline stage.')

    if weight_mode == "ray_count":
        weights = np.sqrt(np.maximum(selected_clouds["ray_hit_count"].to_numpy(dtype=float), 1.0))
    elif weight_mode == "rho":
        weights = np.maximum(selected_clouds["rho_mag_kpc"].to_numpy(dtype=float), 1e-3)
    elif weight_mode == "mass_rho":
        weights = np.log1p(np.maximum(selected_clouds["Mass_Msun"].to_numpy(dtype=float), 0.0)) * np.maximum(selected_clouds["rho_mag_kpc"].to_numpy(dtype=float), 1e-3)
    else:
        weights = np.ones(len(points), dtype=float)

    sqrt_w = np.sqrt(np.asarray(weights, dtype=float))
    centered = points - seed.center[None, :]
    initial_angle = float(seed.angle_deg) % 180.0
    center_fit = np.asarray(seed.center, dtype=float).copy()
    xy_axes_init = np.asarray(seed.axes_init[:2], dtype=float)

    def decode_xy_fixed_z_center_from_params(params):
        uv = np.tanh(np.asarray(params[:2], dtype=float))
        norm = float(np.hypot(uv[0], uv[1]))
        if norm > 1.0:
            uv = uv / norm
        local_offset = np.array(
            [
                float(xy_axes_init[0]) * float(uv[0]),
                float(xy_axes_init[1]) * float(uv[1]),
                0.0,
            ],
            dtype=float,
        )
        global_offset = rotate_points_to_global(local_offset[None, :], seed.angle_deg)[0]
        return np.array(
            [
                float(seed.center[0] + global_offset[0]),
                float(seed.center[1] + global_offset[1]),
                float(seed.center[2]),
            ],
            dtype=float,
        )

    if fit_shape == "cylinder" and circular_xy and free_center_xy_fixed_z:
        c_fixed = 0.1

        def residuals(params):
            center = decode_xy_fixed_z_center_from_params(params)
            radius = np.exp(params[2])
            local = points - center[None, :]
            xy_radius = np.linalg.norm(local[:, :2], axis=1)
            model = (xy_radius / radius) ** 2 - 1.0
            return sqrt_w * model

        r0 = np.sqrt(max(seed.axes_init[0] * seed.axes_init[1], 1e-8))
        x0 = np.array([0.0, 0.0, np.log(r0)], dtype=float)
        lower = [-5.0, -5.0, np.log(0.01)]
        upper = [5.0, 5.0, np.log(5.0)]
        result = least_squares(residuals, x0, bounds=(lower, upper), loss="soft_l1")
        center_fit = decode_xy_fixed_z_center_from_params(result.x)
        radius = float(np.exp(result.x[2]))
        axes_fit = np.array([radius, radius, c_fixed], dtype=float)
        angle_fit = initial_angle
    elif fit_shape == "cylinder" and circular_xy:
        c_fixed = 0.1

        def residuals(params):
            radius = np.exp(params[0])
            xy_radius = np.linalg.norm(centered[:, :2], axis=1)
            model = (xy_radius / radius) ** 2 - 1.0
            return sqrt_w * model

        r0 = np.sqrt(max(seed.axes_init[0] * seed.axes_init[1], 1e-8))
        result = least_squares(residuals, np.array([np.log(r0)], dtype=float), bounds=([np.log(0.01)], [np.log(5.0)]), loss="soft_l1")
        radius = float(np.exp(result.x[0]))
        axes_fit = np.array([radius, radius, c_fixed], dtype=float)
        angle_fit = initial_angle
    elif fit_shape == "cylinder" and free_center_xy_fixed_z:
        c_fixed = 0.1  # Fixed cylinder half-height: h = 200 pc.

        def residuals(params):
            center = decode_xy_fixed_z_center_from_params(params)
            axes_xy = np.exp(params[2:4])
            angle = float(params[4])
            local = rotate_points_to_local(points - center[None, :], angle)
            model = np.sum((local[:, :2] / axes_xy[None, :]) ** 2, axis=1) - 1.0
            return sqrt_w * model

        x0 = np.array([0.0, 0.0, np.log(seed.axes_init[0]), np.log(seed.axes_init[1]), initial_angle], dtype=float)
        lower = [-5.0, -5.0, np.log(0.01), np.log(0.01), 0.0]
        upper = [5.0, 5.0, np.log(5.0), np.log(5.0), 180.0]
        result = least_squares(residuals, x0, bounds=(lower, upper), loss="soft_l1")
        center_fit = decode_xy_fixed_z_center_from_params(result.x)
        axes_fit = np.array([np.exp(result.x[2]), np.exp(result.x[3]), c_fixed], dtype=float)
        angle_fit = float(result.x[4]) % 180.0
    elif fit_shape == "cylinder":
        c_fixed = 0.1  # Fixed cylinder half-height: h = 200 pc.

        def residuals(params):
            axes_xy = np.exp(params[:2])
            angle = float(params[2])
            local = rotate_points_to_local(centered, angle)
            model = np.sum((local[:, :2] / axes_xy[None, :]) ** 2, axis=1) - 1.0
            return sqrt_w * model

        x0 = np.array([np.log(seed.axes_init[0]), np.log(seed.axes_init[1]), initial_angle], dtype=float)
        lower = [np.log(0.01), np.log(0.01), 0.0]
        upper = [np.log(5.0), np.log(5.0), 180.0]
        result = least_squares(residuals, x0, bounds=(lower, upper), loss="soft_l1")
        axes_fit = np.array([np.exp(result.x[0]), np.exp(result.x[1]), c_fixed], dtype=float)
        angle_fit = float(result.x[2]) % 180.0
    elif circular_xy and free_center:
        center_span = np.maximum(np.asarray(seed.axes_init, dtype=float) * float(center_bounds_scale), np.array([0.05, 0.05, 0.03], dtype=float))
        def decode_center_from_params(params):
            uv = np.tanh(np.asarray(params[:2], dtype=float))
            norm = float(np.hypot(uv[0], uv[1]))
            if norm > 1.0:
                uv = uv / norm
            local_offset = np.array(
                [
                    float(xy_axes_init[0]) * float(uv[0]),
                    float(xy_axes_init[1]) * float(uv[1]),
                    0.0,
                ],
                dtype=float,
            )
            global_offset = rotate_points_to_global(local_offset[None, :], seed.angle_deg)[0]
            return np.array(
                [
                    float(seed.center[0] + global_offset[0]),
                    float(seed.center[1] + global_offset[1]),
                    float(params[2]),
                ],
                dtype=float,
            )

        def residuals(params):
            center = decode_center_from_params(params)
            radius_xy, radius_z = np.exp(params[3:5])
            local = rotate_points_to_local(points - center[None, :], initial_angle)
            xy_radius = np.linalg.norm(local[:, :2], axis=1)
            model = (xy_radius / radius_xy) ** 2 + (local[:, 2] / radius_z) ** 2 - 1.0
            return sqrt_w * model

        r0 = np.sqrt(max(seed.axes_init[0] * seed.axes_init[1], 1e-8))
        x0 = np.array([0.0, 0.0, seed.center[2], np.log(r0), np.log(seed.axes_init[2])], dtype=float)
        lower = [-5.0, -5.0, seed.center[2] - center_span[2], np.log(0.01), np.log(0.005)]
        upper = [5.0, 5.0, seed.center[2] + center_span[2], np.log(5.0), np.log(2.0)]
        result = least_squares(residuals, x0, bounds=(lower, upper), loss="soft_l1")
        center_fit = decode_center_from_params(result.x)
        radius_xy, radius_z = np.exp(result.x[3:5])
        axes_fit = np.array([radius_xy, radius_xy, radius_z], dtype=float)
        angle_fit = initial_angle
    elif circular_xy:
        def residuals(params):
            radius_xy, radius_z = np.exp(params[:2])
            local = rotate_points_to_local(centered, initial_angle)
            xy_radius = np.linalg.norm(local[:, :2], axis=1)
            model = (xy_radius / radius_xy) ** 2 + (local[:, 2] / radius_z) ** 2 - 1.0
            return sqrt_w * model

        r0 = np.sqrt(max(seed.axes_init[0] * seed.axes_init[1], 1e-8))
        x0 = np.array([np.log(r0), np.log(seed.axes_init[2])], dtype=float)
        lower = [np.log(0.01), np.log(0.005)]
        upper = [np.log(5.0), np.log(2.0)]
        result = least_squares(residuals, x0, bounds=(lower, upper), loss="soft_l1")
        radius_xy, radius_z = np.exp(result.x[:2])
        axes_fit = np.array([radius_xy, radius_xy, radius_z], dtype=float)
        angle_fit = initial_angle
    elif free_center:
        center_span = np.maximum(np.asarray(seed.axes_init, dtype=float) * float(center_bounds_scale), np.array([0.05, 0.05, 0.03], dtype=float))
        xy_axes_init = np.asarray(seed.axes_init[:2], dtype=float)

        def decode_center_from_params(params):
            # Map two free variables into the initial XY ellipse interior.
            uv = np.tanh(np.asarray(params[:2], dtype=float))
            norm = float(np.hypot(uv[0], uv[1]))
            if norm > 1.0:
                uv = uv / norm
            local_offset = np.array(
                [
                    float(xy_axes_init[0]) * float(uv[0]),
                    float(xy_axes_init[1]) * float(uv[1]),
                    0.0,
                ],
                dtype=float,
            )
            global_offset = rotate_points_to_global(local_offset[None, :], seed.angle_deg)[0]
            return np.array(
                [
                    float(seed.center[0] + global_offset[0]),
                    float(seed.center[1] + global_offset[1]),
                    float(params[2]),
                ],
                dtype=float,
            )

        def residuals(params):
            center = decode_center_from_params(params)
            axes = np.exp(params[3:6])
            angle = float(params[6])
            local = rotate_points_to_local(points - center[None, :], angle)
            model = np.sum((local / axes[None, :]) ** 2, axis=1) - 1.0
            return sqrt_w * model

        x0 = np.array(
            [
                0.0,
                0.0,
                seed.center[2],
                np.log(seed.axes_init[0]),
                np.log(seed.axes_init[1]),
                np.log(seed.axes_init[2]),
                initial_angle,
            ],
            dtype=float,
        )
        lower = [
            -5.0,
            -5.0,
            seed.center[2] - center_span[2],
            np.log(0.01),
            np.log(0.01),
            np.log(0.005),
            0.0,
        ]
        upper = [
            5.0,
            5.0,
            seed.center[2] + center_span[2],
            np.log(5.0),
            np.log(5.0),
            np.log(2.0),
            180.0,
        ]
        result = least_squares(residuals, x0, bounds=(lower, upper), loss="soft_l1")
        center_fit = decode_center_from_params(result.x)
        axes_fit = np.exp(result.x[3:6])
        angle_fit = float(result.x[6]) % 180.0
    elif free_center_z_only:
        center_span = np.maximum(np.asarray(seed.axes_init, dtype=float) * float(center_bounds_scale), np.array([0.05, 0.05, 0.03], dtype=float))

        def residuals(params):
            center = np.array([float(seed.center[0]), float(seed.center[1]), float(params[0])], dtype=float)
            axes = np.exp(params[1:4])
            angle = float(params[4])
            local = rotate_points_to_local(points - center[None, :], angle)
            model = np.sum((local / axes[None, :]) ** 2, axis=1) - 1.0
            return sqrt_w * model

        x0 = np.array(
            [
                seed.center[2],
                np.log(seed.axes_init[0]),
                np.log(seed.axes_init[1]),
                np.log(seed.axes_init[2]),
                initial_angle,
            ],
            dtype=float,
        )
        lower = [
            seed.center[2] - center_span[2],
            np.log(0.01),
            np.log(0.01),
            np.log(0.005),
            0.0,
        ]
        upper = [
            seed.center[2] + center_span[2],
            np.log(5.0),
            np.log(5.0),
            np.log(2.0),
            180.0,
        ]
        result = least_squares(residuals, x0, bounds=(lower, upper), loss="soft_l1")
        center_fit = np.array([float(seed.center[0]), float(seed.center[1]), float(result.x[0])], dtype=float)
        axes_fit = np.exp(result.x[1:4])
        angle_fit = float(result.x[4]) % 180.0
    else:
        def residuals(params):
            axes = np.exp(params[:3])
            angle = float(params[3])
            local = rotate_points_to_local(centered, angle)
            model = np.sum((local / axes[None, :]) ** 2, axis=1) - 1.0
            return sqrt_w * model

        x0 = np.array(
            [
                np.log(seed.axes_init[0]),
                np.log(seed.axes_init[1]),
                np.log(seed.axes_init[2]),
                initial_angle,
            ],
            dtype=float,
        )
        lower = [np.log(0.01), np.log(0.01), np.log(0.005), 0.0]
        upper = [np.log(5.0), np.log(5.0), np.log(2.0), 180.0]
        result = least_squares(residuals, x0, bounds=(lower, upper), loss="soft_l1")
        axes_fit = np.exp(result.x[:3])
        angle_fit = float(result.x[3]) % 180.0

    fit_seed = BubbleSeed(
        target_index=seed.target_index,
        center=center_fit,
        axes_init=seed.axes_init,
        angle_deg=angle_fit,
        mark=seed.mark,
        name=seed.name,
    )
    centered_fit = points - center_fit[None, :]
    local_points = rotate_points_to_local(centered_fit, angle_fit)
    if circular_xy and fit_shape == "cylinder" and free_center_xy_fixed_z:
        local_offset_fit = rotate_points_to_local((center_fit - np.asarray(seed.center, dtype=float))[None, :], seed.angle_deg)[0]
        ux = float(np.clip(local_offset_fit[0] / max(float(xy_axes_init[0]), 1e-9), -0.999999, 0.999999))
        uy = float(np.clip(local_offset_fit[1] / max(float(xy_axes_init[1]), 1e-9), -0.999999, 0.999999))
        residual_params = np.array(
            [
                np.arctanh(ux),
                np.arctanh(uy),
                np.log(axes_fit[0]),
            ],
            dtype=float,
        )
    elif circular_xy and fit_shape == "cylinder":
        residual_params = np.array([np.log(axes_fit[0])], dtype=float)
    elif circular_xy and free_center:
        local_offset_fit = rotate_points_to_local((center_fit - np.asarray(seed.center, dtype=float))[None, :], seed.angle_deg)[0]
        ux = float(np.clip(local_offset_fit[0] / max(float(xy_axes_init[0]), 1e-9), -0.999999, 0.999999))
        uy = float(np.clip(local_offset_fit[1] / max(float(xy_axes_init[1]), 1e-9), -0.999999, 0.999999))
        residual_params = np.array(
            [
                np.arctanh(ux),
                np.arctanh(uy),
                float(center_fit[2]),
                np.log(axes_fit[0]),
                np.log(axes_fit[2]),
            ],
            dtype=float,
        )
    elif circular_xy:
        residual_params = np.array([np.log(axes_fit[0]), np.log(axes_fit[2])], dtype=float)
    elif fit_shape == "cylinder" and free_center_xy_fixed_z:
        local_offset_fit = rotate_points_to_local((center_fit - np.asarray(seed.center, dtype=float))[None, :], seed.angle_deg)[0]
        ux = float(np.clip(local_offset_fit[0] / max(float(xy_axes_init[0]), 1e-9), -0.999999, 0.999999))
        uy = float(np.clip(local_offset_fit[1] / max(float(xy_axes_init[1]), 1e-9), -0.999999, 0.999999))
        residual_params = np.array(
            [
                np.arctanh(ux),
                np.arctanh(uy),
                np.log(axes_fit[0]),
                np.log(axes_fit[1]),
                angle_fit,
            ],
            dtype=float,
        )
    elif fit_shape == "cylinder":
        residual_params = np.r_[np.log(axes_fit[:2]), angle_fit]
    elif free_center:
        local_offset_fit = rotate_points_to_local((center_fit - np.asarray(seed.center, dtype=float))[None, :], seed.angle_deg)[0]
        ux = float(np.clip(local_offset_fit[0] / max(float(xy_axes_init[0]), 1e-9), -0.999999, 0.999999))
        uy = float(np.clip(local_offset_fit[1] / max(float(xy_axes_init[1]), 1e-9), -0.999999, 0.999999))
        residual_params = np.array(
            [
                np.arctanh(ux),
                np.arctanh(uy),
                float(center_fit[2]),
                np.log(axes_fit[0]),
                np.log(axes_fit[1]),
                np.log(axes_fit[2]),
                angle_fit,
            ],
            dtype=float,
        )
    elif free_center_z_only:
        residual_params = np.array(
            [
                float(center_fit[2]),
                np.log(axes_fit[0]),
                np.log(axes_fit[1]),
                np.log(axes_fit[2]),
                angle_fit,
            ],
            dtype=float,
        )
    else:
        residual_params = np.r_[np.log(axes_fit), angle_fit]
    raw_residual = residuals(residual_params) / sqrt_w
    signed_distance, shell_points, _, shell_points_local = shell_distances(points, fit_seed, axes_fit, fit_shape)
    return center_fit, axes_fit, angle_fit, raw_residual, signed_distance, shell_points, shell_points_local, weights



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



def compute_ray_clocal(ray_hits_df, n_azimuth, n_elevation, sigma_pc):
    if ray_hits_df.empty or not np.isfinite(sigma_pc) or sigma_pc <= 0:
        return float("nan"), 0
    radius_grid = np.full((n_elevation, n_azimuth), np.nan, dtype=float)
    for row in ray_hits_df.itertuples(index=False):
        ray_index = int(row.ray_index)
        grid_row = ray_index // n_azimuth
        grid_col = ray_index % n_azimuth
        if 0 <= grid_row < n_elevation:
            radius_grid[grid_row, grid_col] = float(row.t_center_kpc)

    tau = float(sigma_pc) / 1000.0
    sims = []
    for r in range(n_elevation):
        for c in range(n_azimuth):
            if not np.isfinite(radius_grid[r, c]):
                continue
            neighbors = [(r, (c + 1) % n_azimuth), (r + 1, c)]
            for rr, cc in neighbors:
                if 0 <= rr < n_elevation and np.isfinite(radius_grid[rr, cc]):
                    delta = radius_grid[r, c] - radius_grid[rr, cc]
                    sims.append(float(np.exp(-((delta / tau) ** 2))))
    if not sims:
        return float("nan"), 0
    return float(np.mean(sims)), int(len(sims))


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


def read_formal_parameter_table(path: str | Path) -> pd.DataFrame:
    required = [
        "index",
        "id",
        "mark",
        "center_x_kpc",
        "center_y_kpc",
        "center_z_kpc",
        "a_init_kpc",
        "b_init_kpc",
        "c_init_kpc",
        "angle_deg",
    ]
    # The current source file has historical rows with an incompatible schema
    # appended below the formal table. They are outside the v4 formal sample.
    table = pd.read_csv(path, on_bad_lines="skip")
    missing = [column for column in required if column not in table.columns]
    if missing:
        raise ValueError(f"initial parameter table is missing required columns: {', '.join(missing)}")
    table = table.dropna(subset=required).copy()
    table["index"] = table["index"].astype(int)
    table["id"] = table["id"].astype(int)
    table["mark"] = table["mark"].astype(int)
    return table



def load_seed_v4(path: str | Path, target_index: int):
    table = read_formal_parameter_table(path)
    matched = table[table["index"] == int(target_index)]
    if matched.empty:
        raise ValueError(f"target index was not found in the formal sample parameter table index={target_index}")
    row = matched.iloc[0]
    return BubbleSeed(
        target_index=int(target_index),
        center=np.array(
            [row["center_x_kpc"], row["center_y_kpc"], row["center_z_kpc"]],
            dtype=float,
        ),
        axes_init=np.array(
            [row["a_init_kpc"], row["b_init_kpc"], row["c_init_kpc"]],
            dtype=float,
        ),
        angle_deg=float(row["angle_deg"]),
        mark=int(row["mark"]),
        name=f"SB{int(row['id'])}",
    )



def filter_local_bubble_clouds(clouds: pd.DataFrame) -> pd.DataFrame:
    return clouds.loc[clouds["Flag"].to_numpy(dtype=int) != 1].copy().reset_index(drop=True)



def preselect_clouds_v4(clouds: pd.DataFrame, seed, xy_search_scale: float) -> pd.DataFrame:
    points = clouds[["x", "y", "z"]].to_numpy(dtype=float)
    local = rotate_points_to_local(points - seed.center[None, :], seed.angle_deg)
    axes_xy = np.asarray(seed.axes_init[:2], dtype=float) * float(xy_search_scale)
    q_xy = np.sum((local[:, :2] / axes_xy[None, :]) ** 2, axis=1)
    mask = q_xy <= 1.0
    selected = clouds.loc[mask].copy()
    selected["local_x"] = local[mask, 0]
    selected["local_y"] = local[mask, 1]
    selected["local_z"] = local[mask, 2]
    selected["expanded_shape_q"] = q_xy[mask]
    return selected.reset_index(drop=True)



def cloud_distribution_metrics_v4(clouds_df: pd.DataFrame, center_xy: np.ndarray) -> dict:
    """Compute XY-plane azimuthal coverage and centroid offset of selected clouds
    relative to the fitted center. Used to detect 'one-sided' elliptical fits."""
    if clouds_df is None or len(clouds_df) == 0:
        return {
            "n_clouds": 0,
            "azimuth_coverage_deg": 0.0,
            "azimuth_max_gap_deg": 360.0,
            "centroid_offset_kpc": 0.0,
        }
    dx = clouds_df["x"].to_numpy(dtype=float) - float(center_xy[0])
    dy = clouds_df["y"].to_numpy(dtype=float) - float(center_xy[1])
    theta = np.degrees(np.arctan2(dy, dx))
    theta_sorted = np.sort(theta)
    # Largest azimuthal gap around the center (with wraparound).
    gaps = np.diff(np.r_[theta_sorted, theta_sorted[0] + 360.0])
    max_gap = float(gaps.max()) if gaps.size else 360.0
    coverage = float(360.0 - max_gap)
    centroid_offset = float(np.hypot(dx.mean(), dy.mean()))
    return {
        "n_clouds": int(len(clouds_df)),
        "azimuth_coverage_deg": coverage,
        "azimuth_max_gap_deg": max_gap,
        "centroid_offset_kpc": centroid_offset,
    }



def detect_asymmetric_one_sided_ellipse_v4(
    clouds_df: pd.DataFrame,
    center_xy: np.ndarray,
    axes_initial: np.ndarray,
    axis_ratio_threshold: float,
    offset_over_a_threshold: float,
    min_clouds_for_robust_ellipse: int,
    max_ab_ratio: float,
    enabled: bool = True,
) -> dict:
    """Decide whether the elliptical fit looks like an asymmetric 'one-sided' fit
    of a structure that would be physically better described by a circle.

    Trigger fires when ALL conditions hold:
      - axis_ratio_threshold < a/b < max_ab_ratio : the fit is only *moderately*
        elongated (the upper cap leaves genuinely elongated structures, e.g. a/b>3,
        as ellipses);
      - either the sample is too small to robustly constrain an ellipse
        (n_clouds < min_clouds_for_robust_ellipse), or the cloud centroid is offset
        from the fitted center by more than offset_over_a_threshold * a.

    NOTE: "clouds concentrated on one side" is captured by the *centroid offset*
    (centroid_offset / a), NOT by the azimuthal coverage. Coverage is reported as a
    diagnostic only and is intentionally NOT used as a gate: a narrow-arc one-sided
    cloud cluster (small coverage, e.g. SB23) and a wide-wrap but lopsided cluster
    (large coverage but high offset, e.g. SB10/15/25/28) are *both* one-sided fits,
    so a single offset criterion handles both correctly.
    """
    metrics = cloud_distribution_metrics_v4(clouds_df, center_xy)
    axes_xy = [float(axes_initial[0]), float(axes_initial[1])]
    major_axis = float(max(axes_xy))
    minor_axis = float(max(min(axes_xy), 1e-12))
    ratio = float(major_axis / minor_axis)
    offset_over_a = float(metrics["centroid_offset_kpc"] / major_axis) if major_axis > 0 else 0.0
    n_clouds = int(metrics["n_clouds"])
    decision = bool(
        enabled
        and ratio > float(axis_ratio_threshold)
        and ratio < float(max_ab_ratio)
        and (
            n_clouds < int(min_clouds_for_robust_ellipse)
            or offset_over_a > float(offset_over_a_threshold)
        )
    )
    return {
        "trigger": decision,
        "ab_ratio_initial": ratio,
        "azimuth_coverage_deg": metrics["azimuth_coverage_deg"],
        "azimuth_max_gap_deg": metrics["azimuth_max_gap_deg"],
        "centroid_offset_kpc": metrics["centroid_offset_kpc"],
        "centroid_offset_over_a": offset_over_a,
        "n_clouds": n_clouds,
        "axis_ratio_threshold": float(axis_ratio_threshold),
        "offset_over_a_threshold": float(offset_over_a_threshold),
        "min_clouds_for_robust_ellipse": int(min_clouds_for_robust_ellipse),
        "max_ab_ratio": float(max_ab_ratio),
        "enabled": bool(enabled),
    }


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


def clean_value(value):
    if pd.isna(value):
        return None
    if isinstance(value, float):
        if not math.isfinite(value):
            return None
        return round(value, 8)
    if hasattr(value, "item"):
        return clean_value(value.item())
    return value



def build_fit_cloud_payload(clouds):
    cloud_columns = [
        "Seq",
        "GLON",
        "GLAT",
        "Dist",
        "Rad_pc",
        "Area_deg2",
        "Mass_Msun",
        "Sigma_Msun_pc2",
        "rho_mag_kpc",
        "Flag",
        "x",
        "y",
        "z",
        "rad_kpc",
        "local_x",
        "local_y",
        "local_z",
        "ray_hit_count",
        "signed_distance_pc",
        "weight",
    ]
    available_columns = [column for column in cloud_columns if column in clouds.columns]
    payload = []
    for record in clouds[available_columns].to_dict(orient="records"):
        payload.append({key: clean_value(value) for key, value in record.items()})
    return payload



def compute_cap_z_from_selected_clouds(clouds, center_z, mark_value):
    if clouds is None or len(clouds) == 0:
        return None
    z_values = np.asarray(clouds["z"], dtype=float)
    center_z = float(center_z)
    mark_value = int(mark_value)
    if mark_value == 1:
        subset = np.sort(z_values[z_values < center_z])[::-1]
        if subset.size == 0:
            return None
        rank = max(1, int(np.ceil(subset.size / np.e)))
        return float(subset[min(rank - 1, subset.size - 1)])
    if mark_value == 2:
        subset = np.sort(z_values[z_values > center_z])
        if subset.size == 0:
            return None
        rank = max(1, int(np.ceil(subset.size / np.e)))
        return float(subset[min(rank - 1, subset.size - 1)])
    return None



def compute_xy_plane_ellipse_params(fit_shape, mark_value, center, axes, angle_deg, clouds):
    center_x = clean_value(center[0])
    center_y = clean_value(center[1])
    if fit_shape == "cylinder" or int(mark_value) == 3:
        return {
            "xy_plane_center_x_kpc": center_x,
            "xy_plane_center_y_kpc": center_y,
            "xy_plane_z_kpc": clean_value(center[2]),
            "xy_plane_a_kpc": clean_value(axes[0]),
            "xy_plane_b_kpc": clean_value(axes[1]),
            "xy_plane_angle_deg": clean_value(angle_deg),
            "xy_plane_definition": "cylinder_base_ellipse",
        }

    cap_z = compute_cap_z_from_selected_clouds(clouds, center[2], mark_value)
    if cap_z is None:
        return {
            "xy_plane_center_x_kpc": center_x,
            "xy_plane_center_y_kpc": center_y,
            "xy_plane_z_kpc": None,
            "xy_plane_a_kpc": None,
            "xy_plane_b_kpc": None,
            "xy_plane_angle_deg": clean_value(angle_deg),
            "xy_plane_definition": "cap_base_ellipse",
        }

    a_axis = float(axes[0])
    b_axis = float(axes[1])
    c_axis = float(axes[2])
    if c_axis <= 0:
        scale = None
    else:
        ratio = 1.0 - ((float(cap_z) - float(center[2])) / c_axis) ** 2
        scale = float(np.sqrt(ratio)) if ratio > 0 else None

    return {
        "xy_plane_center_x_kpc": center_x,
        "xy_plane_center_y_kpc": center_y,
        "xy_plane_z_kpc": clean_value(cap_z),
        "xy_plane_a_kpc": clean_value(a_axis * scale) if scale is not None else None,
        "xy_plane_b_kpc": clean_value(b_axis * scale) if scale is not None else None,
        "xy_plane_angle_deg": clean_value(angle_deg),
        "xy_plane_definition": "cap_base_ellipse",
    }


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


def build_best_fit_table(success_rows: list[dict]) -> pd.DataFrame:
    rows = []
    mcmc_columns = [
        "mcmc_center_x_median_kpc",
        "mcmc_center_x_err_minus_kpc",
        "mcmc_center_x_err_plus_kpc",
        "mcmc_center_y_median_kpc",
        "mcmc_center_y_err_minus_kpc",
        "mcmc_center_y_err_plus_kpc",
        "mcmc_center_z_median_kpc",
        "mcmc_center_z_err_minus_kpc",
        "mcmc_center_z_err_plus_kpc",
        "mcmc_a_median_kpc",
        "mcmc_a_err_minus_kpc",
        "mcmc_a_err_plus_kpc",
        "mcmc_b_median_kpc",
        "mcmc_b_err_minus_kpc",
        "mcmc_b_err_plus_kpc",
        "mcmc_c_median_kpc",
        "mcmc_c_err_minus_kpc",
        "mcmc_c_err_plus_kpc",
        "mcmc_angle_median_deg",
        "mcmc_angle_err_minus_deg",
        "mcmc_angle_err_plus_deg",
        "mcmc_acceptance_fraction_mean",
    ]
    for status_row in sorted(success_rows, key=lambda item: int(item["id"])):
        data = json.loads(package_path(status_row["json_path"]).read_text(encoding="utf-8"))
        clouds = pd.read_csv(package_path(data["outputs"]["selected_clouds_csv"]))
        center = data["center_kpc"]
        axes = data["fit_axes_kpc"]
        fit_shape = str(data["fit_shape"])
        xy_plane = compute_xy_plane_ellipse_params(
            fit_shape,
            int(data["mark"]),
            center,
            axes,
            data["angle_deg"],
            clouds,
        )
        row = {
            "id": int(data["target_index"]),
            "version": "v4_all",
            "mark": int(data["mark"]),
            "shape": fit_shape,
            "xy_model": data.get("xy_model", "elliptical"),
            "circular_xy": bool(data.get("circular_xy", False)),
            "potential_circular_xy_by_axis_ratio": bool(data.get("potential_circular_xy_by_axis_ratio", False)),
            "final_parameter_estimator": data.get("final_parameter_estimator", "least_squares"),
            "final_estimator_reason": data.get("final_estimator_reason"),
            "initial_free_ab_axis_ratio": clean_value(data.get("initial_free_ab_axis_ratio")),
            "fit_ab_axis_ratio": clean_value(data.get("fit_ab_axis_ratio")),
            "least_squares_a_kpc": clean_value(data["least_squares_axes_kpc"][0]),
            "least_squares_b_kpc": clean_value(data["least_squares_axes_kpc"][1]),
            "least_squares_c_kpc": clean_value(data["least_squares_axes_kpc"][2]),
            "least_squares_angle_deg": clean_value(data["least_squares_angle_deg"]),
            "center_x_kpc": clean_value(center[0]),
            "center_y_kpc": clean_value(center[1]),
            "center_z_kpc": clean_value(center[2]),
            "a_radius_kpc": clean_value(axes[0]),
            "b_radius_kpc": clean_value(axes[1]),
            "c_radius_kpc": 0.0 if fit_shape == "cylinder" else clean_value(axes[2]),
            "angle_deg": clean_value(data["angle_deg"]),
            "fit_std_pc": clean_value(data["fit_std_pc"]),
            "C_local": clean_value(data["C_local"]),
            "catalog_cloud_count": int(data["catalog_cloud_count"]),
            "local_bubble_removed_cloud_count": int(data["local_bubble_removed_cloud_count"]),
            "candidate_cloud_count": int(data["candidate_cloud_count"]),
            "selected_cloud_count": int(data["selected_cloud_count"]),
            "ray_count": int(data["ray_count"]),
            "ray_hit_count": int(data["ray_hit_count"]),
            "hit_ray_count": int(data["hit_ray_count"]),
            "catalog_filter": data["parameters"]["catalog_filter"],
            "xy_search_scale": float(data["parameters"]["xy_search_scale"]),
            "fit_center_bounds_scale": float(data["parameters"]["fit_center_bounds_scale"]),
            "search_scale": float(data["parameters"]["search_scale"]),
            "radius_tolerance_pc": float(data["parameters"]["radius_tolerance_pc"]),
            "ray_cone_aperture_deg": float(data["parameters"].get("ray_cone_aperture_deg", 0.0)),
            "max_hits_per_ray": int(data["parameters"]["max_hits_per_ray"]),
            "fit_cloud_seq_list": json.dumps([int(value) for value in clouds["Seq"].tolist()], ensure_ascii=False),
            "fit_clouds_json": json.dumps(
                build_fit_cloud_payload(clouds),
                ensure_ascii=False,
                separators=(",", ":"),
            ),
            **xy_plane,
        }
        for column in mcmc_columns:
            row[column] = clean_value(data.get(column))
        rows.append(row)
    return pd.DataFrame(rows)



def fit_one_target(sb, seed, selected_clouds, ray_hits, all_clouds):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    fit_shape = fit_shape_for_seed(seed)
    free_center = fit_shape == "ellipsoid"
    free_center_xy_fixed_z = fit_shape == "cylinder"


    filtered_clouds = filter_local_bubble_clouds(all_clouds)
    candidate_clouds = preselect_clouds_v4(filtered_clouds, seed, 1.5)


    center_initial, axes_initial, angle_initial, _, signed_distance, _, _, weights = fit_cloud_shell(
        selected_clouds, seed, WEIGHT_MODE, fit_shape,
        circular_xy=False, free_center=free_center,
        free_center_xy_fixed_z=free_center_xy_fixed_z, free_center_z_only=False,
        center_bounds_scale=FIT_CENTER_BOUNDS_SCALE,
    )
    initial_ab_ratio = ab_axis_ratio(axes_initial)


    potential_circular_xy = (
        not DISABLE_CIRCULAR_XY_FALLBACK
        and should_use_circular_xy(axes_initial, CIRCULAR_AXIS_RATIO_THRESHOLD)
    )
    asymmetric_detection = detect_asymmetric_one_sided_ellipse_v4(
        selected_clouds, center_initial[:2], axes_initial,
        axis_ratio_threshold=CIRCULAR_AXIS_RATIO_THRESHOLD,
        offset_over_a_threshold=ASYM_OFFSET_OVER_A,
        min_clouds_for_robust_ellipse=ASYM_MIN_CLOUDS,
        max_ab_ratio=ASYM_MAX_AB_RATIO,
        enabled=True,
    )
    asymmetric_circular = bool(asymmetric_detection["trigger"])
    circular_xy = bool(potential_circular_xy or asymmetric_circular)


    if circular_xy:
        center_fit, axes_fit, angle_fit, _, signed_distance, _, _, weights = fit_cloud_shell(
            selected_clouds, seed, WEIGHT_MODE, fit_shape,
            circular_xy=True, free_center=free_center,
            free_center_xy_fixed_z=free_center_xy_fixed_z, free_center_z_only=False,
            center_bounds_scale=FIT_CENTER_BOUNDS_SCALE,
        )
    else:
        center_fit, axes_fit, angle_fit = center_initial, axes_initial, angle_initial

    sigma_info = sigma_from_signed_distance(signed_distance)
    first_hit_rays = ray_hits[ray_hits["hit_rank"] == 1].copy() if len(ray_hits) else ray_hits
    c_local, n_pairs = compute_ray_clocal(first_hit_rays, N_AZIMUTH, 2 * N_ELEVATION, sigma_info["sigma_pc"])


    selected_rho = selected_clouds["rho_mag_kpc"].to_numpy(dtype=float)
    background = candidate_clouds[~candidate_clouds["Seq"].isin(selected_clouds["Seq"])]["rho_mag_kpc"].to_numpy(dtype=float)
    rho_selected_median = float(np.median(selected_rho)) if selected_rho.size else float("nan")
    rho_background_median = float(np.median(background)) if background.size else float("nan")
    rho_excess = float(rho_selected_median - rho_background_median) if np.isfinite(rho_background_median) else float("nan")

    clouds_csv = JSON_DIR / f"SB{sb}_clouds.csv"
    json_path = JSON_DIR / f"SB{sb}_fit.json"

    result = {
        "target_index": int(seed.target_index),
        "version": "peer_review_v4",
        "method": "molecular_cloud_first_touch_v4_full_sky_xy1p5_local_bubble_removed_cone_ray_radial_upper_limit",
        "fit_shape": fit_shape,
        "xy_model": "circular" if circular_xy else "elliptical",
        "circular_xy": bool(circular_xy),
        "potential_circular_xy_by_axis_ratio": bool(potential_circular_xy),
        "circular_axis_ratio_threshold": float(CIRCULAR_AXIS_RATIO_THRESHOLD),
        "circular_by_asymmetric_distribution": bool(asymmetric_circular),
        "asymmetric_distribution_detection": asymmetric_detection,
        "final_parameter_estimator": "least_squares",
        "final_estimator_reason": "least_squares_initial",
        "least_squares_axes_kpc": axes_fit.tolist(),
        "least_squares_angle_deg": float(angle_fit),
        "least_squares_center_kpc": center_fit.tolist(),
        "initial_free_ab_axis_ratio": float(initial_ab_ratio),
        "fit_ab_axis_ratio": float(ab_axis_ratio(axes_fit)),
        "center_kpc": center_fit.tolist(),
        "initial_center_kpc": seed.center.tolist(),
        "angle_deg": float(angle_fit),
        "initial_angle_deg": float(seed.angle_deg),
        "mark": int(seed.mark),
        "hemisphere": "full",
        "initial_axes_kpc": seed.axes_init.tolist(),
        "fit_axes_kpc": axes_fit.tolist(),
        "fit_a_kpc": float(axes_fit[0]),
        "fit_b_kpc": float(axes_fit[1]),
        "fit_c_kpc": float(axes_fit[2]),
        "fit_sigma_pc": float(sigma_info["sigma_pc"]),
        "fit_std_pc": float(sigma_info["std_pc"]),
        "fit_mean_pc": float(sigma_info["mean_pc"]),
        "fit_median_pc": float(sigma_info["median_pc"]),
        "fit_mad_pc": float(sigma_info["mad_pc"]),
        "fit_iqr_sigma_pc": float(sigma_info["iqr_sigma_pc"]),
        "fit_p68_abs_pc": float(sigma_info["p68_abs_pc"]),
        "fit_central_68_half_width_pc": float(sigma_info["central_68_half_width_pc"]),
        "C_local": float(c_local),
        "C_local_pairs": int(n_pairs),
        "catalog_cloud_count": int(len(all_clouds)),
        "local_bubble_filtered_cloud_count": int(len(filtered_clouds)),
        "local_bubble_removed_cloud_count": int(len(all_clouds) - len(filtered_clouds)),
        "candidate_cloud_count": int(len(candidate_clouds)),
        "selected_cloud_count": int(len(selected_clouds)),
        "ray_count": int(N_AZIMUTH * 2 * N_ELEVATION),
        "ray_hit_count": int(len(ray_hits)),
        "hit_ray_count": int(ray_hits["ray_index"].nunique()) if len(ray_hits) else 0,
        "rho_selected_median_mag_kpc": rho_selected_median,
        "rho_background_median_mag_kpc": rho_background_median,
        "rho_excess_mag_kpc": rho_excess,
        "parameters": {
            "catalog_filter": "Flag != 1 (remove Local Bubble boundary clouds)",
            "xy_search_scale": 1.5,
            "fit_center_bounds_scale": float(FIT_CENTER_BOUNDS_SCALE),
            "ray_model": "circular_cone_with_aperture",
            "ray_cone_aperture_deg": 5.0,
            "search_scale": 1.5,
            "max_hits_per_ray": 1,
            "radius_tolerance_pc": 0.0,
            "n_azimuth": int(N_AZIMUTH),
            "n_elevation_per_hemisphere": int(N_ELEVATION),
            "weight_mode": WEIGHT_MODE,
            "circular_axis_ratio_threshold": float(CIRCULAR_AXIS_RATIO_THRESHOLD),
        },
        "outputs": {
            "selected_clouds_csv": package_relative(clouds_csv),
            "json": package_relative(json_path),
        },
    }

    # Figure styling and layout configuration for reproducible rendering.
    selected_out = selected_clouds.copy()
    selected_out["signed_distance_pc"] = signed_distance * 1000.0
    selected_out["weight"] = weights
    selected_out.to_csv(clouds_csv, index=False)
    json_path.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")

    return {
        "result": result,
        "json_path": json_path,
        "seed": seed,
        "selected_clouds": selected_clouds,
        "center_fit": center_fit,
        "axes_fit": axes_fit,
        "angle_fit": angle_fit,
        "sigma_info": sigma_info,
        "weights": weights,
        "circular_xy": circular_xy,
        "fit_shape": fit_shape,
        "free_center": free_center,
        "free_center_xy_fixed_z": free_center_xy_fixed_z,
    }


def main():
    """Run Script 2 from validated inputs to the documented outputs."""
    JSON_DIR.mkdir(parents=True, exist_ok=True)

    formal = read_formal_parameter_table(PARAM_CSV)
    all_ids = formal["index"].astype(int).tolist()
    auto_ids = [sb for sb in all_ids if sb not in SKIP_AUTO_IDS]
    all_clouds = parse_mcs_catalog(str(MCS_CSV))

    json_paths = {}
    bestfit_rows = []
    print(f"Running 3D geometry fits for {len(auto_ids)} automatically fitted OSBs.")
    for sb in auto_ids:
        seed = load_seed_v4(str(PARAM_CSV), sb)
        # Data-loading helper section.
        clouds_csv = STAGE1_DIR / f"SB{sb}_first_contact_molecular_clouds.csv"
        rays_csv = STAGE1_DIR / f"SB{sb}_ray_hits.csv"
        if not clouds_csv.exists():
            raise FileNotFoundError('Missing input file.')
        selected_clouds = pd.read_csv(clouds_csv)
        ray_hits = pd.read_csv(rays_csv) if rays_csv.exists() else pd.DataFrame(columns=["ray_index", "hit_rank"])

        fit = fit_one_target(sb, seed, selected_clouds, ray_hits, all_clouds)
        json_paths[sb] = fit["json_path"]
        bestfit_rows.append({"id": sb, "status": "ok", "json_path": package_relative(fit["json_path"])})
        r = fit["result"]
        print(f"SB{sb:<2d}: shape={r['fit_shape']:<9s} xy={r['xy_model']:<10s} "
              f"N={r['selected_cloud_count']:<4d} fit_std={r['fit_std_pc']:.2f} pc")


    bestfit = build_best_fit_table(bestfit_rows)
    bestfit.to_csv(BESTFIT_CSV, index=False, encoding="utf-8-sig")
    print(f"Saved best-fit geometry table: {BESTFIT_CSV}")
    print(
        "Stage 2 complete: selected MC subsets, 3D geometry fits, and fit "
        "diagnostics were written."
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
