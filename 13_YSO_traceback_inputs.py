"""Script 13: YSO traceback inputs

Purpose
-------
Build G1/G2 Hunt and Konietzka young-cluster traceback input tables.

Method overview
---------------
1. Derive the Hunt catalogue kinematics and uncertainties, apply the configured sample
   cuts, and standardize the Hunt and Konietzka columns.
2. Identify cross-catalogue duplicates using the name and astrometric/position matching
   routines, retaining the configured source priorities.
3. Build the combined master table and the G1/G2 young-cluster samples used by the orbit
   and slice analyses.

Main inputs
-----------
- ../data/star_cluster_data/hunt2024_clusters_full.csv
- ../data/star_cluster_data/Konietzka2023.csv
- ../results/intermediate_output/13_traceback_cluster_input_data/G1_traceback_cluster_sample.csv when rebuilding only the G2 stage

Main outputs
------------
- ../results/intermediate_output/13_traceback_cluster_input_data/; input to scripts 14, 15, 16, 23, 24.

Figure/table role
-----------------
../results/intermediate_output/13_traceback_cluster_input_data/; input to scripts 14, 15, 16, 23, 24.

Runtime and data notes
----------------------
Reference runtime: 3.65 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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
import re
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord


HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
RAW_CLUSTER_DIR = DATA_DIR / "star_cluster_data"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "13_traceback_cluster_input_data"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"

SOLAR_MOTION_LSR_STANDARD = (11.1, 12.24, 7.25)
HUNT_VOLUME_LIMITS_PC = {
    "x_min": -3000.0,
    "x_max": 3000.0,
    "y_min": -3000.0,
    "y_max": 3000.0,
    "z_min": -500.0,
    "z_max": 500.0,
}

HUNT_MIN_CST = 4.0
HUNT_MIN_CMD_CLASS50 = 0.3
HUNT_MIN_NDIST = 5
HUNT_MAX_AGE_MYR = 70.0
HUNT_MIN_RV_COUNT = 2
HUNT_MAX_RV_ERROR_KMS = None
HUNT_MAX_UVW_ERROR_KMS = 14.0
HUNT_FAR_DISTANCE_PC = 1000.0
HUNT_MAX_RV_ERROR_FAR_KMS = None
HUNT_MAX_UVW_ERROR_FAR_KMS = 10.0
HUNT_MC_SAMPLES = 256
HUNT_MC_SEED = 42

FIRST_GROUP_SPATIAL_THRESHOLD_PC = 35.0
FIRST_GROUP_TANGENTIAL_VELOCITY_THRESHOLD_KMS = 5.0
SECOND_GROUP_MIN_DISTANCE_PC = 200.0
SECOND_GROUP_MAX_DISTANCE_PC = 2000.0
SECOND_GROUP_MAX_AGE_MYR = 50.0
SECOND_GROUP_DEDUP_SPATIAL_PC = 35.0
SECOND_GROUP_DEDUP_TANGENTIAL_KMS = 5.0

HUNT_NUMERIC_COLUMNS = [
    "ID",
    "CST",
    "N",
    "CSTt",
    "Nt",
    "RA_ICRS",
    "DE_ICRS",
    "GLON",
    "GLAT",
    "r50",
    "rc",
    "rt",
    "rtot",
    "r50pc",
    "rcpc",
    "rtpc",
    "rtotpc",
    "pmRA",
    "s_pmRA",
    "e_pmRA",
    "pmDE",
    "s_pmDE",
    "e_pmDE",
    "Plx",
    "s_Plx",
    "e_Plx",
    "dist16",
    "dist50",
    "dist84",
    "Ndist",
    "globalPlx",
    "X",
    "Y",
    "Z",
    "RV",
    "s_RV",
    "e_RV",
    "n_RV",
    "CMDCl2.5",
    "CMDCl16",
    "CMDCl50",
    "CMDCl84",
    "CMDCl97.5",
    "logAge16",
    "logAge50",
    "logAge84",
    "AV16",
    "AV50",
    "AV84",
    "diffAV16",
    "diffAV50",
    "diffAV84",
    "MOD16",
    "MOD50",
    "MOD84",
    "r50J",
    "rJ",
    "r50Jpc",
    "rJpc",
    "probJ",
    "NJ",
    "MassJ",
    "e_MassJ",
    "MassTot",
    "e_MassTot",
    "minClSize",
    "isMerged",
    "isGMMMemb",
    "NXmatches",
]


def require_file(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path}")
    return path


def add_derived_hunt_columns(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    result["age_myr"] = np.power(10.0, result["logAge50"] - 6.0)
    result["age16_myr"] = np.power(10.0, result["logAge16"] - 6.0)
    result["age84_myr"] = np.power(10.0, result["logAge84"] - 6.0)
    result["distance_sigma_pc"] = 0.5 * (result["dist84"] - result["dist16"])

    glon_rad = np.deg2rad(result["GLON"].to_numpy(dtype=np.float64))
    glat_rad = np.deg2rad(result["GLAT"].to_numpy(dtype=np.float64))
    dist_pc = result["dist50"].to_numpy(dtype=np.float64)

    result["x_helio"] = dist_pc * np.cos(glat_rad) * np.cos(glon_rad)
    result["y_helio"] = dist_pc * np.cos(glat_rad) * np.sin(glon_rad)
    result["z_helio"] = dist_pc * np.sin(glat_rad)
    return result


def load_hunt_catalog(path: Path) -> pd.DataFrame:
    df = pd.read_csv(require_file(path))
    for column in HUNT_NUMERIC_COLUMNS:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return add_derived_hunt_columns(df)


def build_nominal_skycoord(df: pd.DataFrame) -> SkyCoord:
    return SkyCoord(
        ra=df["RA_ICRS"].to_numpy(dtype=np.float64) * u.deg,
        dec=df["DE_ICRS"].to_numpy(dtype=np.float64) * u.deg,
        distance=df["dist50"].to_numpy(dtype=np.float64) * u.pc,
        pm_ra_cosdec=df["pmRA"].to_numpy(dtype=np.float64) * u.mas / u.yr,
        pm_dec=df["pmDE"].to_numpy(dtype=np.float64) * u.mas / u.yr,
        radial_velocity=df["RV"].to_numpy(dtype=np.float64) * u.km / u.s,
        frame="icrs",
    )


def compute_heliocentric_uvw(df: pd.DataFrame) -> pd.DataFrame:
    galactic = build_nominal_skycoord(df).galactic
    result = df.copy()
    result["U"] = galactic.velocity.d_x.to_value(u.km / u.s)
    result["V"] = galactic.velocity.d_y.to_value(u.km / u.s)
    result["W"] = galactic.velocity.d_z.to_value(u.km / u.s)
    result["U_lsr"] = result["U"] + SOLAR_MOTION_LSR_STANDARD[0]
    result["V_lsr"] = result["V"] + SOLAR_MOTION_LSR_STANDARD[1]
    result["W_lsr"] = result["W"] + SOLAR_MOTION_LSR_STANDARD[2]
    return result


def estimate_uvw_uncertainties(df: pd.DataFrame, mc_samples: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n_rows = len(df)

    ra_deg = np.repeat(df["RA_ICRS"].to_numpy(dtype=np.float64), mc_samples)
    dec_deg = np.repeat(df["DE_ICRS"].to_numpy(dtype=np.float64), mc_samples)

    dist_center = df["dist50"].to_numpy(dtype=np.float64)
    dist_sigma = np.clip(df["distance_sigma_pc"].to_numpy(dtype=np.float64), 1.0, None)
    sampled_dist = rng.normal(dist_center[:, None], dist_sigma[:, None], size=(n_rows, mc_samples))
    sampled_dist = np.clip(sampled_dist, 1.0, None)

    pmra_center = df["pmRA"].to_numpy(dtype=np.float64)
    pmra_sigma = np.clip(df["e_pmRA"].to_numpy(dtype=np.float64), 0.0, None)
    sampled_pmra = rng.normal(pmra_center[:, None], pmra_sigma[:, None], size=(n_rows, mc_samples))

    pmde_center = df["pmDE"].to_numpy(dtype=np.float64)
    pmde_sigma = np.clip(df["e_pmDE"].to_numpy(dtype=np.float64), 0.0, None)
    sampled_pmde = rng.normal(pmde_center[:, None], pmde_sigma[:, None], size=(n_rows, mc_samples))

    rv_center = df["RV"].to_numpy(dtype=np.float64)
    rv_sigma = np.clip(df["e_RV"].to_numpy(dtype=np.float64), 0.0, None)
    sampled_rv = rng.normal(rv_center[:, None], rv_sigma[:, None], size=(n_rows, mc_samples))

    coords = SkyCoord(
        ra=ra_deg * u.deg,
        dec=dec_deg * u.deg,
        distance=sampled_dist.reshape(-1) * u.pc,
        pm_ra_cosdec=sampled_pmra.reshape(-1) * u.mas / u.yr,
        pm_dec=sampled_pmde.reshape(-1) * u.mas / u.yr,
        radial_velocity=sampled_rv.reshape(-1) * u.km / u.s,
        frame="icrs",
    ).galactic

    u_samples = coords.velocity.d_x.to_value(u.km / u.s).reshape(n_rows, mc_samples)
    v_samples = coords.velocity.d_y.to_value(u.km / u.s).reshape(n_rows, mc_samples)
    w_samples = coords.velocity.d_z.to_value(u.km / u.s).reshape(n_rows, mc_samples)

    result = df.copy()
    result["U_err"] = np.std(u_samples, axis=1, ddof=1)
    result["V_err"] = np.std(v_samples, axis=1, ddof=1)
    result["W_err"] = np.std(w_samples, axis=1, ddof=1)
    return result


def apply_hunt_preselection(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    stages: list[tuple[str, int]] = []
    current = df.copy()
    stages.append(("raw Hunt catalogue", len(current)))
    stages.append(('Selection stage', len(current)))

    current = current[current["CST"] > HUNT_MIN_CST].copy()
    stages.append((f"CST > {HUNT_MIN_CST}", len(current)))
    current = current[current["CMDCl50"] > HUNT_MIN_CMD_CLASS50].copy()
    stages.append((f"CMDCl50 > {HUNT_MIN_CMD_CLASS50}", len(current)))
    current = current[current["Ndist"] >= HUNT_MIN_NDIST].copy()
    stages.append((f"Ndist >= {HUNT_MIN_NDIST}", len(current)))
    current = current[current["age_myr"] < HUNT_MAX_AGE_MYR].copy()
    stages.append((f"age < {HUNT_MAX_AGE_MYR} Myr", len(current)))

    current = current[current["x_helio"].between(HUNT_VOLUME_LIMITS_PC["x_min"], HUNT_VOLUME_LIMITS_PC["x_max"])].copy()
    current = current[current["y_helio"].between(HUNT_VOLUME_LIMITS_PC["y_min"], HUNT_VOLUME_LIMITS_PC["y_max"])].copy()
    current = current[current["z_helio"].between(HUNT_VOLUME_LIMITS_PC["z_min"], HUNT_VOLUME_LIMITS_PC["z_max"])].copy()
    stages.append(('Selection stage', len(current)))

    current = current[current["n_RV"] >= HUNT_MIN_RV_COUNT].copy()
    stages.append((f"n_RV >= {HUNT_MIN_RV_COUNT}", len(current)))
    stages.append(("no hard cut on Hunt group-level e_RV", len(current)))

    return current.reset_index(drop=True), pd.DataFrame(stages, columns=["stage", "count"])


def finalize_hunt_sample(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    preselected, summary = apply_hunt_preselection(df)
    with_uvw = compute_heliocentric_uvw(preselected)
    with_errors = estimate_uvw_uncertainties(with_uvw, mc_samples=HUNT_MC_SAMPLES, seed=HUNT_MC_SEED)

    base_mask = np.isfinite(with_errors["U"]) & np.isfinite(with_errors["V"]) & np.isfinite(with_errors["W"])
    distances_pc = with_errors["dist50"].to_numpy(dtype=np.float64)
    near_mask = distances_pc <= HUNT_FAR_DISTANCE_PC
    far_mask = ~near_mask
    near_quality = (
        near_mask
        & (with_errors["U_err"] <= HUNT_MAX_UVW_ERROR_KMS)
        & (with_errors["V_err"] <= HUNT_MAX_UVW_ERROR_KMS)
        & (with_errors["W_err"] <= HUNT_MAX_UVW_ERROR_KMS)
    )
    far_quality = (
        far_mask
        & (with_errors["U_err"] <= HUNT_MAX_UVW_ERROR_FAR_KMS)
        & (with_errors["V_err"] <= HUNT_MAX_UVW_ERROR_FAR_KMS)
        & (with_errors["W_err"] <= HUNT_MAX_UVW_ERROR_FAR_KMS)
    )
    final_sample = with_errors.loc[base_mask & (near_quality | far_quality)].copy().reset_index(drop=True)
    uvw_stage = (
        "rows with finite UVW velocities and distance-dependent UVW-error cuts "
        "applied"
    )
    summary = pd.concat(
        [summary, pd.DataFrame([(uvw_stage, len(final_sample))], columns=["stage", "count"])],
        ignore_index=True,
    )
    return final_sample, summary


def normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def extract_aliases(*values: str) -> list[str]:
    aliases: list[str] = []
    for value in values:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            continue
        for token in str(value).split(","):
            token = token.strip()
            if token:
                normalized = normalize_name(token)
                if normalized:
                    aliases.append(normalized)
    return sorted(set(aliases))


def add_traceback_columns(
    result: pd.DataFrame,
    rv_condition_status: str,
    velocity_status: str,
) -> pd.DataFrame:
    positions = result[["x_helio", "y_helio", "z_helio"]].to_numpy(dtype=float)
    velocities = result[["U_lsr", "V_lsr", "W_lsr"]].to_numpy(dtype=float)
    tangential = np.array([tangential_velocity_from_arrays(pos, vel) for pos, vel in zip(positions, velocities)])
    result["U_tan_lsr"] = tangential[:, 0]
    result["V_tan_lsr"] = tangential[:, 1]
    result["W_tan_lsr"] = tangential[:, 2]
    result["has_full_3d_velocity"] = result[["U_lsr", "V_lsr", "W_lsr"]].notna().all(axis=1)
    result["rv_condition_status"] = rv_condition_status
    result["velocity_status"] = np.where(result["has_full_3d_velocity"], velocity_status, "missing_3d_velocity")
    return result


def standardize_konietzka(path: Path) -> pd.DataFrame:
    df = pd.read_csv(require_file(path))
    result = pd.DataFrame(
        {
            "name": df["name"].astype(str),
            "display_name": df["name"].astype(str),
            "source_catalog": "konietzka2023",
            "source_priority": 1,
            "age_myr": pd.to_numeric(df["age"], errors="coerce"),
            "dist_pc": np.sqrt(
                pd.to_numeric(df["x"], errors="coerce") ** 2
                + pd.to_numeric(df["y"], errors="coerce") ** 2
                + pd.to_numeric(df["z"], errors="coerce") ** 2
            ),
            "x_helio": pd.to_numeric(df["x"], errors="coerce"),
            "y_helio": pd.to_numeric(df["y"], errors="coerce"),
            "z_helio": pd.to_numeric(df["z"], errors="coerce"),
            "U_lsr": pd.to_numeric(df["vx"], errors="coerce"),
            "V_lsr": pd.to_numeric(df["vy"], errors="coerce"),
            "W_lsr": pd.to_numeric(df["vz"], errors="coerce"),
            "x_helio_err": pd.to_numeric(df["sigmax"], errors="coerce"),
            "y_helio_err": pd.to_numeric(df["sigmay"], errors="coerce"),
            "z_helio_err": pd.to_numeric(df["sigmaz"], errors="coerce"),
            "U_err": pd.to_numeric(df["sigmavx"], errors="coerce"),
            "V_err": pd.to_numeric(df["sigmavy"], errors="coerce"),
            "W_err": pd.to_numeric(df["sigmavz"], errors="coerce"),
            "rv_kms": pd.NA,
            "rv_err": pd.NA,
            "n_rv": pd.NA,
            "size_pc": pd.NA,
            "quality_note": 'Standardized catalogue input.',
            "reference_literature": df["reference"].astype(str),
            "inlier_flag": pd.to_numeric(df["inlier"], errors="coerce"),
        }
    )
    result = add_traceback_columns(
        result,
        rv_condition_status="not_applicable_konietzka2023",
        velocity_status="full_3d_velocity_available",
    )
    result["alias_keys"] = ["|".join(extract_aliases(name)) for name in result["name"]]
    return result


def standardize_hunt(
    path: Path,
    source_catalog: str = "hunt_stratified",
    source_priority_value: int = 2,
    quality_note: str = 'Standardized catalogue input.',
) -> pd.DataFrame:
    df = pd.read_csv(require_file(path))
    result = pd.DataFrame(
        {
            "name": df["Name"].astype(str),
            "display_name": df["Name"].astype(str),
            "source_catalog": source_catalog,
            "source_priority": source_priority_value,
            "age_myr": pd.to_numeric(df["age_myr"], errors="coerce"),
            "dist_pc": pd.to_numeric(df["dist50"], errors="coerce"),
            "x_helio": pd.to_numeric(df["x_helio"], errors="coerce"),
            "y_helio": pd.to_numeric(df["y_helio"], errors="coerce"),
            "z_helio": pd.to_numeric(df["z_helio"], errors="coerce"),
            "U_lsr": pd.to_numeric(df["U_lsr"], errors="coerce"),
            "V_lsr": pd.to_numeric(df["V_lsr"], errors="coerce"),
            "W_lsr": pd.to_numeric(df["W_lsr"], errors="coerce"),
            "x_helio_err": 0.0,
            "y_helio_err": 0.0,
            "z_helio_err": 0.0,
            "U_err": pd.to_numeric(df["U_err"], errors="coerce"),
            "V_err": pd.to_numeric(df["V_err"], errors="coerce"),
            "W_err": pd.to_numeric(df["W_err"], errors="coerce"),
            "rv_kms": pd.to_numeric(df["RV"], errors="coerce"),
            "rv_err": pd.to_numeric(df["e_RV"], errors="coerce"),
            "n_rv": pd.to_numeric(df["n_RV"], errors="coerce"),
            "size_pc": pd.NA,
            "quality_note": quality_note,
            "reference_literature": "Hunt & Reffert 2024",
            "inlier_flag": 1,
        }
    )
    result = add_traceback_columns(
        result,
        rv_condition_status="pass_n_RV_ge_2",
        velocity_status="full_3d_velocity_available",
    )
    all_names = df["AllNames"] if "AllNames" in df.columns else pd.Series([""] * len(df))
    result["alias_keys"] = [
        "|".join(extract_aliases(name, alias_names))
        for name, alias_names in zip(df["Name"], all_names, strict=False)
    ]
    return result


def read_standardized_li_yso(path: Path) -> pd.DataFrame:
    df = pd.read_csv(require_file(path))
    expected_columns = [
        "name",
        "display_name",
        "source_catalog",
        "reference_literature",
        "age_myr",
        "dist_pc",
        "x_helio",
        "y_helio",
        "z_helio",
        "U_lsr",
        "V_lsr",
        "W_lsr",
        "U_tan_lsr",
        "V_tan_lsr",
        "W_tan_lsr",
        "x_helio_err",
        "y_helio_err",
        "z_helio_err",
        "U_err",
        "V_err",
        "W_err",
        "rv_kms",
        "rv_err",
        "n_rv",
        "has_full_3d_velocity",
        "rv_condition_status",
        "velocity_status",
        "quality_note",
        "inlier_flag",
        "alias_keys",
        "source_priority",
    ]
    result = df.copy()
    for column in expected_columns:
        if column not in result.columns:
            result[column] = pd.NA
    result = result[expected_columns].copy()
    result["source_catalog"] = "li2025_yso_group"
    result["reference_literature"] = result["reference_literature"].fillna("Li Yunqian 2025 YSO catalog")
    result["has_full_3d_velocity"] = True
    result["source_priority"] = 3
    return result


def tangential_velocity_from_arrays(position: np.ndarray, velocity: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(position))
    if not np.isfinite(norm) or norm <= 0.0:
        return np.full(3, np.nan)
    radial_direction = position / norm
    return velocity - float(np.dot(velocity, radial_direction)) * radial_direction


def find_konietzka_hunt_astrometric_duplicate(
    row: pd.Series,
    kept_rows: list[pd.Series],
    max_spatial_sep_pc: float,
    max_tangential_velocity_sep_kms: float,
) -> tuple[str, float, float] | None:
    row_pos = np.array([row["x_helio"], row["y_helio"], row["z_helio"]], dtype=float)
    row_vel = np.array([row["U_lsr"], row["V_lsr"], row["W_lsr"]], dtype=float)
    row_tan = tangential_velocity_from_arrays(row_pos, row_vel)
    if not np.all(np.isfinite(row_tan)):
        return None

    best_match: tuple[str, float, float] | None = None
    best_score = np.inf
    for kept in kept_rows:
        if kept["source_catalog"] != "konietzka2023":
            continue
        kept_pos = np.array([kept["x_helio"], kept["y_helio"], kept["z_helio"]], dtype=float)
        kept_vel = np.array([kept["U_lsr"], kept["V_lsr"], kept["W_lsr"]], dtype=float)
        kept_tan = tangential_velocity_from_arrays(kept_pos, kept_vel)
        if not np.all(np.isfinite(kept_tan)):
            continue
        spatial_sep = float(np.linalg.norm(row_pos - kept_pos))
        tangential_velocity_sep = float(np.linalg.norm(row_tan - kept_tan))
        if spatial_sep <= max_spatial_sep_pc and tangential_velocity_sep <= max_tangential_velocity_sep_kms:
            score = spatial_sep / max_spatial_sep_pc + tangential_velocity_sep / max_tangential_velocity_sep_kms
            if score < best_score:
                best_score = score
                best_match = str(kept["display_name"]), spatial_sep, tangential_velocity_sep
    return best_match


def deduplicate_catalog(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    kept_rows: list[pd.Series] = []
    duplicate_rows: list[dict[str, str]] = []
    seen_aliases: dict[str, str] = {}
    kept_name_to_index: dict[str, int] = {}

    def merge_aliases_into_kept(kept_name: str, aliases: list[str]) -> None:
        kept_index = kept_name_to_index.get(kept_name)
        if kept_index is None:
            return
        kept_aliases = [alias for alias in str(kept_rows[kept_index]["alias_keys"]).split("|") if alias]
        kept_rows[kept_index]["alias_keys"] = "|".join(sorted(set(kept_aliases).union(aliases)))

    for _, row in df.sort_values(["source_priority", "dist_pc", "age_myr"], kind="stable").iterrows():
        aliases = [alias for alias in str(row["alias_keys"]).split("|") if alias]
        duplicate_of = next((seen_aliases[alias] for alias in aliases if alias in seen_aliases), None)
        if duplicate_of is not None:
            merge_aliases_into_kept(duplicate_of, aliases)
            duplicate_rows.append(
                {
                    "dropped_name": str(row["display_name"]),
                    "kept_name": duplicate_of,
                    "source_catalog": str(row["source_catalog"]),
                    "reason": "duplicate name or alias",
                }
            )
            continue

        if row["source_catalog"] == "hunt_stratified":
            astrometric_duplicate = find_konietzka_hunt_astrometric_duplicate(
                row,
                kept_rows,
                max_spatial_sep_pc=30.0,
                max_tangential_velocity_sep_kms=5.0,
            )
            if astrometric_duplicate is not None:
                kept_name, spatial_sep, tangential_velocity_sep = astrometric_duplicate
                merge_aliases_into_kept(kept_name, aliases)
                duplicate_rows.append(
                    {
                        "dropped_name": str(row["display_name"]),
                        "kept_name": kept_name,
                        "source_catalog": str(row["source_catalog"]),
                        "reason": (
                            "alias/name match plus astrometric proximity; "
                            f"Δr={spatial_sep:.1f} pc, Δv_tan={tangential_velocity_sep:.1f} km/s"
                        ),
                    }
                )
                continue

        kept_rows.append(row)
        kept_name_to_index[str(row["display_name"])] = len(kept_rows) - 1
        for alias in aliases:
            seen_aliases[alias] = str(row["display_name"])

    merged = pd.DataFrame(kept_rows).reset_index(drop=True)
    duplicates = pd.DataFrame(duplicate_rows)
    return merged, duplicates


def tangential_velocity(row: pd.Series) -> np.ndarray:
    position = np.array([row["x_helio"], row["y_helio"], row["z_helio"]], dtype=float)
    velocity = np.array([row["U_lsr"], row["V_lsr"], row["W_lsr"]], dtype=float)
    return tangential_velocity_from_arrays(position, velocity)


def find_konietzka_xyz_tangent_duplicate(
    hunt_row: pd.Series,
    konietzka_rows: list[pd.Series],
    spatial_threshold_pc: float,
    tangential_velocity_threshold_kms: float,
) -> tuple[str, float, float] | None:
    hunt_position = np.array([hunt_row["x_helio"], hunt_row["y_helio"], hunt_row["z_helio"]], dtype=float)
    hunt_tangent = tangential_velocity(hunt_row)

    best_match: tuple[str, float, float] | None = None
    for kept in konietzka_rows:
        kept_position = np.array([kept["x_helio"], kept["y_helio"], kept["z_helio"]], dtype=float)
        spatial_sep = float(np.linalg.norm(hunt_position - kept_position))
        if not np.isfinite(spatial_sep) or spatial_sep > spatial_threshold_pc:
            continue
        tangent_sep = float(np.linalg.norm(hunt_tangent - tangential_velocity(kept)))
        if not np.isfinite(tangent_sep) or tangent_sep > tangential_velocity_threshold_kms:
            continue
        if best_match is None or spatial_sep < best_match[1]:
            best_match = (str(kept["display_name"]), spatial_sep, tangent_sep)
    return best_match


def deduplicate_with_konietzka_xyz_tangent(
    combined: pd.DataFrame,
    spatial_threshold_pc: float,
    tangential_velocity_threshold_kms: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    name_merged, name_duplicates = deduplicate_catalog(combined)

    kept_rows: list[pd.Series] = []
    extra_duplicates: list[dict[str, str]] = []
    konietzka_rows: list[pd.Series] = []
    for _, row in name_merged.sort_values(["source_priority", "dist_pc", "age_myr"], kind="stable").iterrows():
        if row["source_catalog"] == "konietzka2023":
            kept_rows.append(row)
            konietzka_rows.append(row)
            continue
        if row["source_catalog"] == "hunt_stratified":
            duplicate = find_konietzka_xyz_tangent_duplicate(
                hunt_row=row,
                konietzka_rows=konietzka_rows,
                spatial_threshold_pc=spatial_threshold_pc,
                tangential_velocity_threshold_kms=tangential_velocity_threshold_kms,
            )
            if duplicate is not None:
                kept_name, spatial_sep, tangent_sep = duplicate
                extra_duplicates.append(
                    {
                        "dropped_name": str(row["display_name"]),
                        "kept_name": kept_name,
                        "source_catalog": str(row["source_catalog"]),
                        "reason": (
                            "Konietzka XYZ/tangential duplicate; "
                            f"Δr={spatial_sep:.1f} pc, Δv_tan={tangent_sep:.1f} km/s"
                        ),
                    }
                )
                continue
        kept_rows.append(row)

    merged = pd.DataFrame(kept_rows).reset_index(drop=True)
    duplicates = pd.concat([name_duplicates, pd.DataFrame(extra_duplicates)], ignore_index=True)
    return merged, duplicates


def deduplicate_by_spatial_tangential_priority(
    df: pd.DataFrame,
    spatial_threshold_pc: float,
    tangential_velocity_threshold_kms: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge all rows connected by XYZ and transverse-velocity proximity."""
    positions = df[["x_helio", "y_helio", "z_helio"]].to_numpy(dtype=float)
    tangential = np.array([tangential_velocity(row) for _, row in df.iterrows()])
    _, find, union = union_find(len(df))
    pair_rows: list[dict[str, object]] = []

    for i in range(len(df)):
        if i + 1 >= len(df):
            break
        dr = np.linalg.norm(positions[i + 1 :] - positions[i], axis=1)
        dtan = np.linalg.norm(tangential[i + 1 :] - tangential[i], axis=1)
        duplicate_indices = np.where((dr <= spatial_threshold_pc) & (dtan <= tangential_velocity_threshold_kms))[0] + i + 1
        for j in duplicate_indices:
            pair_index = int(j) - i - 1
            union(i, int(j))
            pair_rows.append(
                {
                    "left_name": df.loc[i, "display_name"],
                    "right_name": df.loc[int(j), "display_name"],
                    "left_source": df.loc[i, "source_catalog"],
                    "right_source": df.loc[int(j), "source_catalog"],
                    "xyz_sep_pc": float(dr[pair_index]),
                    "tangential_velocity_sep_kms": float(dtan[pair_index]),
                    "reason": (
                        f"XYZ separation<={spatial_threshold_pc:g} pc and "
                        f"tangential-velocity difference<={tangential_velocity_threshold_kms:g} km/s"
                    ),
                }
            )

    components: dict[int, list[int]] = {}
    for index in range(len(df)):
        components.setdefault(find(index), []).append(index)

    keep_indices: list[int] = []
    duplicate_rows: list[dict[str, object]] = []
    for group_id, indices in enumerate(components.values(), start=1):
        ranked = sorted(
            indices,
            key=lambda idx: (
                source_priority(df.loc[idx, "source_catalog"]),
                pd.to_numeric(df.loc[idx, "dist_pc"], errors="coerce"),
                str(df.loc[idx, "display_name"]),
            ),
        )
        keep = ranked[0]
        keep_indices.append(keep)
        for dropped in ranked[1:]:
            spatial_sep = float(np.linalg.norm(positions[dropped] - positions[keep]))
            tangent_sep = float(np.linalg.norm(tangential[dropped] - tangential[keep]))
            duplicate_rows.append(
                {
                    "duplicate_group": group_id,
                    "kept_name": df.loc[keep, "display_name"],
                    "kept_source": df.loc[keep, "source_catalog"],
                    "dropped_name": df.loc[dropped, "display_name"],
                    "dropped_source": df.loc[dropped, "source_catalog"],
                    "xyz_sep_to_kept_pc": spatial_sep,
                    "tangential_velocity_sep_to_kept_kms": tangent_sep,
                    "reason": (
                        f"XYZ separation<={spatial_threshold_pc:g} pc and "
                        f"tangential-velocity difference<={tangential_velocity_threshold_kms:g} km/s; "
                        "retention priority konietzka2023 > hunt2024/hunt_stratified > li2025_yso_group"
                    ),
                }
            )

    output = df.iloc[sorted(keep_indices)].copy().reset_index(drop=True)
    duplicate_log = pd.DataFrame(duplicate_rows)
    pair_log = pd.DataFrame(pair_rows)
    if not pair_log.empty:
        pair_log["is_duplicate_pair"] = True
    if not duplicate_log.empty:
        duplicate_log = duplicate_log.merge(
            pair_log,
            left_on=["kept_name", "dropped_name"],
            right_on=["left_name", "right_name"],
            how="left",
            suffixes=("", "_pair"),
        )
    return output, duplicate_log


def export_minimal_catalog(df: pd.DataFrame) -> pd.DataFrame:
    result = df[
        [
            "name",
            "x_helio",
            "y_helio",
            "z_helio",
            "U_lsr",
            "V_lsr",
            "W_lsr",
            "age_myr",
            "source_catalog",
            "reference_literature",
        ]
    ].copy()
    return result.rename(
        columns={
            "x_helio": "x_pc",
            "y_helio": "y_pc",
            "z_helio": "z_pc",
            "U_lsr": "vx_lsr_kms",
            "V_lsr": "vy_lsr_kms",
            "W_lsr": "vz_lsr_kms",
            "source_catalog": "source_key",
            "reference_literature": "source_literature",
        }
    )


def concat_catalog_frames(frames: list[pd.DataFrame]) -> pd.DataFrame:
    columns: list[str] = []
    trimmed_frames: list[pd.DataFrame] = []
    for frame in frames:
        for col in frame.columns:
            if col not in columns:
                columns.append(col)
        trimmed_frames.append(frame.dropna(axis=1, how="all"))

    combined = pd.concat(trimmed_frames, ignore_index=True)
    for col in columns:
        if col not in combined.columns:
            combined[col] = pd.NA
    return combined[columns]


def prepare_hunt_sample() -> Path:
    hunt_raw = RAW_CLUSTER_DIR / "hunt2024_clusters_full.csv"
    hunt_df = load_hunt_catalog(hunt_raw)
    hunt_sample, hunt_summary = finalize_hunt_sample(hunt_df)

    hunt_sample_file = OUT_DIR / "hunt_3kpc_sample_stratified.csv"
    hunt_summary_file = OUT_DIR / "hunt_3kpc_filter_summary_stratified.csv"
    hunt_sample.to_csv(hunt_sample_file, index=False)
    hunt_summary.to_csv(hunt_summary_file, index=False)
    return hunt_sample_file


def build_hunt_konietzka_table() -> None:
    konietzka_raw = RAW_CLUSTER_DIR / "Konietzka2023.csv"
    hunt_sample_file = prepare_hunt_sample()

    combined = concat_catalog_frames(
        [
            standardize_konietzka(konietzka_raw),
            standardize_hunt(hunt_sample_file),
        ]
    )
    merged, duplicates = deduplicate_with_konietzka_xyz_tangent(
        combined=combined,
        spatial_threshold_pc=FIRST_GROUP_SPATIAL_THRESHOLD_PC,
        tangential_velocity_threshold_kms=FIRST_GROUP_TANGENTIAL_VELOCITY_THRESHOLD_KMS,
    )
    merged = merged.sort_values(["dist_pc", "source_priority", "display_name"], kind="stable").reset_index(drop=True)

    merged_file = OUT_DIR / "G1_traceback_cluster_sample.csv"
    minimal_file = OUT_DIR / "Hunt_Konietzka_master_table_compact_columns.csv"
    duplicate_file = OUT_DIR / "Hunt_Konietzka_processing_log.csv"
    merged.to_csv(merged_file, index=False, encoding="utf-8-sig")
    export_minimal_catalog(merged).to_csv(minimal_file, index=False, encoding="utf-8-sig")
    duplicates.to_csv(duplicate_file, index=False, encoding="utf-8-sig")

    print(f"G1 sample: {merged_file} ({len(merged)} rows)")
    print(merged["source_catalog"].value_counts().to_string())


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def numeric_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    result = df.copy()
    for column in columns:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    return result


def source_priority(value: object) -> int:
    priorities = {
        "konietzka2023": 0,
        "hunt2024": 1,
        "hunt_stratified": 1,
        "li2025_yso_group": 2,
    }
    return priorities.get(str(value), 99)


def union_find(n_items: int) -> tuple[list[int], callable, callable]:
    parent = list(range(n_items))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    return parent, find, union


def build_young_traceback_table() -> None:
    g1_sample_file = OUT_DIR / "G1_traceback_cluster_sample.csv"
    if not g1_sample_file.exists():
        print("not found ../results/intermediate_output/13_traceback_cluster_input_data/G1_traceback_cluster_sample.csv; build the G1 sample first.")
        build_hunt_konietzka_table()

    g1_df = pd.read_csv(require_file(g1_sample_file))
    if "source_priority" not in g1_df.columns:
        g1_df["source_priority"] = g1_df["source_catalog"].map(source_priority)
    df = g1_df.copy()
    initial_rows = len(df)
    need_numeric = ["dist_pc"]
    df = numeric_columns(df, need_numeric)

    too_near = df["dist_pc"] < SECOND_GROUP_MIN_DISTANCE_PC
    too_far = df["dist_pc"] > SECOND_GROUP_MAX_DISTANCE_PC
    in_distance = (~too_near) & (~too_far) & df["dist_pc"].notna()

    output = df.loc[in_distance].copy().reset_index(drop=True)
    output = output.sort_values(["dist_pc", "source_catalog", "display_name"], kind="stable").reset_index(drop=True)
    if "row_id" in output.columns:
        output["row_id"] = np.arange(1, len(output) + 1)

    output_file = OUT_DIR / "G2_traceback_cluster_sample.csv"

    output.to_csv(output_file, index=False, encoding="utf-8-sig")
    for source_name, file_stem in [
        ("hunt_stratified", "G2_Hunt_source_sample"),
        ("konietzka2023", "G2_Konietzka_source_sample"),
    ]:
        source_output = output.loc[output["source_catalog"].astype(str) == source_name].copy()
        source_output.to_csv(OUT_DIR / f"{file_stem}.csv", index=False, encoding="utf-8-sig")

    for stale_name in [
        "G2_LiYSO_source_sample.csv",
        "G2_spatial_tangential_deduplication_log.csv",
        "G2_six_dimensional_deduplication_log.csv",
        "G2_six_dimensional_near_neighbor_log.csv",
    ]:
        stale_path = OUT_DIR / stale_name
        if stale_path.exists():
            stale_path.unlink()

    summary_file = OUT_DIR / "G2_selection_summary.csv"
    pd.DataFrame(
        [
            ("input_total", initial_rows),
            ("remove_dist_lt_200_pc", int(too_near.sum())),
            ("remove_dist_gt_2000_pc", int(too_far.sum())),
            ("output_total", len(output)),
        ],
        columns=["stage", "count"],
    ).to_csv(summary_file, index=False, encoding="utf-8-sig")

    print(f"G2 sample: {output_file} ({len(output)} rows)")
    print(output["source_catalog"].value_counts().to_string())


def main() -> None:
    """Run Script 13 from validated inputs to the documented outputs."""
    parser = argparse.ArgumentParser(description='Run Script 13: YSO traceback inputs.')
    parser.add_argument(
        "--only",
        choices=["all", "first", "second"],
        default="all",
        help='Command-line option for the documented workflow.',
    )
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)

    if args.only in {"all", "first"}:
        build_hunt_konietzka_table()
    if args.only in {"all", "second"}:
        build_young_traceback_table()

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
