# -*- coding: utf-8 -*-
"""Script 1: OSBs first contact MCs

Purpose
-------
Identify first-contact molecular clouds along rays for each open-superbubble seed.

Method overview
---------------
1. Load each pre-identified seed and the molecular-cloud catalogue.
2. Use find_first_touch_for_seed() to scan the configured rays and select the
   first-contact clouds for that seed.
3. Write each target's selected clouds and selection metadata for the shell fitting in
   Scripts 2 and 3.

Main inputs
-----------
- ../data/OSB_initial_pre_identification_parameters.csv
- ../data/MCs.csv

Main outputs
------------
- ../results/intermediate_output/1_first_contact_molecular_clouds/; input to scripts 2 and 3.

Figure/table role
-----------------
../results/intermediate_output/1_first_contact_molecular_clouds/; input to scripts 2 and 3.

Runtime and data notes
----------------------
Reference runtime: 3.27 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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
from dataclasses import dataclass
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord


# Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
# Data-loading helper section.
DATA_DIR = Path("..") / "data"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "1_first_contact_molecular_clouds"

PARAM_CSV = DATA_DIR / "OSB_initial_pre_identification_parameters.csv"   # Geometry and shell-mask conventions used by this analysis stage.
MCS_CSV = DATA_DIR / "MCs.csv"


def package_relative(path: Path) -> str:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    path = Path(path)
    return path.relative_to(HERE).as_posix() if path.is_absolute() else path.as_posix()

XY_SEARCH_SCALE = 1.5
N_AZIMUTH = 72
N_ELEVATION = 18
RAY_CONE_APERTURE_DEG = 5.0
RADIUS_TOLERANCE_PC = 0.0
MAX_HITS_PER_RAY = 1
SEARCH_SCALE = 1.5


# ============================================================================
# ----------------------------------------------------------------------------
# Geometry and shell-mask conventions used by this analysis stage.
# ============================================================================


@dataclass
class BubbleSeed:
    """Container for parameters passed between helpers in this pipeline stage."""
    target_index: int
    center: np.ndarray
    axes_init: np.ndarray
    angle_deg: float
    mark: int
    name: str


def rotation_matrix_z(angle_deg):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    theta = np.deg2rad(float(angle_deg))
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rotate_points_to_local(points, angle_deg):
    """Rotate coordinates into the local frame used by the shell or projection calculation."""
    rot = rotation_matrix_z(-angle_deg)
    return np.asarray(points) @ rot.T


def sample_directions(n_azimuth, n_elevation, hemisphere):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    epsilon = np.deg2rad(4.0)
    azimuths = np.linspace(0.0, 2.0 * np.pi, n_azimuth, endpoint=False)
    if hemisphere == "lower":
        elevations = np.linspace(-0.5 * np.pi + epsilon, -epsilon, n_elevation)
    elif hemisphere == "upper":
        elevations = np.linspace(epsilon, 0.5 * np.pi - epsilon, n_elevation)
    else:
        lower = np.linspace(-0.5 * np.pi + epsilon, -epsilon, n_elevation)
        upper = np.linspace(epsilon, 0.5 * np.pi - epsilon, n_elevation)
        elevations = np.concatenate([lower, upper])

    directions = []
    for elevation in elevations:
        cos_el = np.cos(elevation)
        sin_el = np.sin(elevation)
        for azimuth in azimuths:
            directions.append(
                {
                    "direction": np.array(
                        [cos_el * np.cos(azimuth), cos_el * np.sin(azimuth), sin_el],
                        dtype=float,
                    ),
                    "azimuth_deg": float(np.degrees(azimuth)),
                    "elevation_deg": float(np.degrees(elevation)),
                }
            )
    return directions


def ray_radius_to_ellipsoid(direction, axes, angle_deg):
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    local_direction = rotate_points_to_local(np.asarray(direction, dtype=float)[None, :], angle_deg)[0]
    denom = np.sum((local_direction / np.asarray(axes, dtype=float)) ** 2)
    if denom <= 0:
        raise ValueError('Invalid input or missing required data.')
    return 1.0 / np.sqrt(denom)


def ray_radius_to_cylinder(direction, axes, angle_deg, *, min_xy_norm=0.2):
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    local_direction = rotate_points_to_local(np.asarray(direction, dtype=float)[None, :], angle_deg)[0]
    xy_norm = float(np.hypot(local_direction[0], local_direction[1]))
    if xy_norm < min_xy_norm:
        return np.nan
    denom = (local_direction[0] / float(axes[0])) ** 2 + (local_direction[1] / float(axes[1])) ** 2
    if denom <= 0:
        return np.nan
    return 1.0 / np.sqrt(denom)


def parse_mcs_catalog(path):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
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
    coords = SkyCoord(
        l=df["GLON"].to_numpy() * u.deg,
        b=df["GLAT"].to_numpy() * u.deg,
        distance=df["Dist"].to_numpy() * u.kpc,
        frame="galactic",
    )
    cart = coords.cartesian
    df["x"] = cart.x.to_value(u.kpc)
    df["y"] = cart.y.to_value(u.kpc)
    df["z"] = cart.z.to_value(u.kpc)
    df["rad_kpc"] = df["Rad_pc"] / 1000.0
    return df


def fit_shape_for_seed(seed):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    return "cylinder" if int(seed.mark) == 3 else "ellipsoid"


def read_formal_parameter_table(path):
    """Load and validate the input table or array needed by this stage."""
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
    table = pd.read_csv(path, on_bad_lines="skip")
    missing = [column for column in required if column not in table.columns]
    if missing:
        raise ValueError(f"initial parameter table is missing required columns: {', '.join(missing)}")
    table = table.dropna(subset=required).copy()
    table["index"] = table["index"].astype(int)
    table["id"] = table["id"].astype(int)
    table["mark"] = table["mark"].astype(int)
    return table


def load_seed_v4(path, target_index):
    """Load and validate the input table or array needed by this stage."""
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


def filter_local_bubble_clouds(clouds):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    return clouds.loc[clouds["Flag"].to_numpy(dtype=int) != 1].copy().reset_index(drop=True)


def preselect_clouds_v4(clouds, seed, xy_search_scale):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
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


def prepare_full_sky_rays_v4(n_azimuth, n_elevation):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    rays = []
    for ray_index, meta in enumerate(sample_directions(n_azimuth, n_elevation, "full")):
        rays.append(
            {
                "ray_index": int(ray_index),
                "direction": np.asarray(meta["direction"], dtype=float),
                "azimuth_deg": float(meta["azimuth_deg"]),
                "elevation_deg": float(meta["elevation_deg"]),
            }
        )
    return rays


def find_first_touch_hits_v4(
    candidate_clouds,
    seed,
    rays,
    radius_tolerance_kpc,
    ray_cone_aperture_deg,
    max_hits_per_ray,
    fit_shape,
    search_scale,
):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    if candidate_clouds.empty:
        return candidate_clouds.copy(), pd.DataFrame()

    cloud_offsets = candidate_clouds[["x", "y", "z"]].to_numpy(dtype=float) - seed.center[None, :]
    base_radii = candidate_clouds["rad_kpc"].to_numpy(dtype=float)
    cloud_radii = base_radii + float(radius_tolerance_kpc)
    cloud_distances = np.linalg.norm(cloud_offsets, axis=1)
    safe_distances = np.maximum(cloud_distances, 1e-12)
    ratio = cloud_radii / safe_distances
    envelops_origin = ratio >= 1.0
    cloud_angular_radii_rad = np.where(
        envelops_origin, np.pi, np.arcsin(np.clip(ratio, 0.0, 1.0))
    )
    cloud_unit_vectors = cloud_offsets / safe_distances[:, None]
    cone_half_angle_rad = np.deg2rad(float(ray_cone_aperture_deg) / 2.0)
    ray_hits = []

    for ray in rays:
        direction = np.asarray(ray["direction"], dtype=float)
        if fit_shape == "cylinder":
            t_guess = ray_radius_to_cylinder(direction, seed.axes_init, seed.angle_deg, min_xy_norm=1e-6)
        else:
            t_guess = ray_radius_to_ellipsoid(direction, seed.axes_init, seed.angle_deg)
        if not np.isfinite(t_guess) or t_guess <= 0:
            continue
        t_max = float(search_scale) * t_guess

        cos_phi = np.clip(cloud_unit_vectors @ direction, -1.0, 1.0)
        phi = np.arccos(cos_phi)
        intersect_mask = phi <= (cone_half_angle_rad + cloud_angular_radii_rad)
        if not np.any(intersect_mask):
            continue

        phi_minus_alpha = np.maximum(phi - cone_half_angle_rad, 0.0)
        perp_dist = cloud_distances * np.sin(phi_minus_alpha)
        along_axis = cloud_distances * np.cos(phi_minus_alpha)
        half_chord_sq = np.maximum(cloud_radii**2 - perp_dist**2, 0.0)
        half_chord = np.sqrt(half_chord_sq)
        t_entry = np.maximum(along_axis - half_chord, 0.0)
        t_exit = along_axis + half_chord
        within_limit_mask = (t_exit >= 0.0) & (t_entry <= t_max)
        hit_mask = intersect_mask & within_limit_mask
        if not np.any(hit_mask):
            continue

        hit_indices = np.where(hit_mask)[0]
        ordered_positions = np.argsort(t_entry[hit_indices])[: int(max_hits_per_ray)]
        for hit_rank, position in enumerate(ordered_positions, start=1):
            cloud_index = int(hit_indices[position])
            best = candidate_clouds.iloc[cloud_index]
            ray_hits.append(
                {
                    "ray_index": int(ray["ray_index"]),
                    "hit_rank": int(hit_rank),
                    "azimuth_deg": float(ray["azimuth_deg"]),
                    "elevation_deg": float(ray["elevation_deg"]),
                    "direction_x": float(direction[0]),
                    "direction_y": float(direction[1]),
                    "direction_z": float(direction[2]),
                    "t_guess_kpc": float(t_guess),
                    "t_center_kpc": float(along_axis[cloud_index]),
                    "t_hit_kpc": float(t_entry[cloud_index]),
                    "phi_deg": float(np.rad2deg(phi[cloud_index])),
                    "cone_aperture_deg": float(ray_cone_aperture_deg),
                    "cloud_seq": int(best["Seq"]),
                    "cloud_x": float(best["x"]),
                    "cloud_y": float(best["y"]),
                    "cloud_z": float(best["z"]),
                    "cloud_rad_kpc": float(best["rad_kpc"]),
                    "cloud_effective_radius_kpc": float(cloud_radii[cloud_index]),
                    "cloud_rho_mag_kpc": float(best["rho_mag_kpc"]),
                    "cloud_mass_msun": float(best["Mass_Msun"]),
                }
            )

    ray_hits_df = pd.DataFrame(ray_hits)
    if ray_hits_df.empty:
        return candidate_clouds.iloc[0:0].copy(), ray_hits_df

    hit_count = ray_hits_df.groupby("cloud_seq").size().rename("ray_hit_count")
    selected = candidate_clouds[candidate_clouds["Seq"].isin(hit_count.index)].copy()
    selected = selected.merge(hit_count, left_on="Seq", right_index=True, how="left")
    selected["ray_hit_count"] = selected["ray_hit_count"].astype(int)
    return selected.reset_index(drop=True), ray_hits_df


def find_first_touch_for_seed(seed, all_clouds):
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    fit_shape = fit_shape_for_seed(seed)


    filtered_clouds = filter_local_bubble_clouds(all_clouds)

    candidate_clouds = preselect_clouds_v4(filtered_clouds, seed, XY_SEARCH_SCALE)

    rays = prepare_full_sky_rays_v4(N_AZIMUTH, N_ELEVATION)

    selected_clouds, ray_hits = find_first_touch_hits_v4(
        candidate_clouds,
        seed,
        rays,
        RADIUS_TOLERANCE_PC / 1000.0,   # pc -> kpc
        RAY_CONE_APERTURE_DEG,
        MAX_HITS_PER_RAY,
        fit_shape,
        SEARCH_SCALE,
    )
    return {
        "fit_shape": fit_shape,
        "candidate_clouds": candidate_clouds,
        "selected_clouds": selected_clouds,
        "ray_hits": ray_hits,
        "n_rays": len(rays),
        "n_all_clouds": len(all_clouds),
        "n_after_local_bubble": len(filtered_clouds),
    }


def main():
    """Run Script 1 from validated inputs to the documented outputs."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Data-loading helper section.
    formal = read_formal_parameter_table(PARAM_CSV)
    target_ids = formal["index"].astype(int).tolist()


    all_clouds = parse_mcs_catalog(str(MCS_CSV))

    summary_rows = []
    for sb in target_ids:
        seed = load_seed_v4(str(PARAM_CSV), sb)
        res = find_first_touch_for_seed(seed, all_clouds)
        selected = res["selected_clouds"]
        ray_hits = res["ray_hits"]


        clouds_csv = OUT_DIR / f"SB{sb}_first_contact_molecular_clouds.csv"
        rays_csv = OUT_DIR / f"SB{sb}_ray_hits.csv"
        selected.to_csv(clouds_csv, index=False, encoding="utf-8-sig")
        ray_hits.to_csv(rays_csv, index=False, encoding="utf-8-sig")

        n_hit_rays = int(ray_hits["ray_index"].nunique()) if len(ray_hits) else 0
        summary_rows.append(
            {
                "id": sb,
                "fit_shape": res["fit_shape"],
                "catalog_cloud_count": res["n_all_clouds"],
                "after_local_bubble_count": res["n_after_local_bubble"],
                "candidate_cloud_count": int(len(res["candidate_clouds"])),
                "selected_cloud_count": int(len(selected)),
                "ray_count": res["n_rays"],
                "ray_hit_count": int(len(ray_hits)),
                "hit_ray_count": n_hit_rays,
                "clouds_csv": package_relative(clouds_csv),
                "ray_hits_csv": package_relative(rays_csv),
            }
        )
        print(
            f"SB{sb:<2d}: shape={res['fit_shape']:<9s} "
            f"candidates={len(res['candidate_clouds']):<4d} "
            f"first-touch clouds={len(selected):<4d} hit rays={n_hit_rays}/{res['n_rays']}"
        )

    summary = pd.DataFrame(summary_rows)
    summary_csv = OUT_DIR / "first_contact_identification_summary.csv"
    summary.to_csv(summary_csv, index=False, encoding="utf-8-sig")
    print(f"Saved first-contact summary: {summary_csv}")
    print(f"Processed OSBs: {len(summary)}")
    print(f"Per-OSB first-contact products are in: {OUT_DIR}")

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
