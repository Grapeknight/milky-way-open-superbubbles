"""Script 14: G1 YSO traceback

Purpose
-------
Integrate G1 orbits for 55 Myr and run age-weighted HDBSCAN family recovery.

Method overview
---------------
1. Load the G1 sample prepared by Script 13 and integrate or reuse its 55 Myr backward
   orbit tracks.
2. Run the HDBSCAN family inference with the configured age weights, then consolidate
   the family labels.
3. Measure family properties and size evolution and export member tables for Scripts 16
   and 17; the recompute flags control reuse of existing products.

Main inputs
-----------
- ../results/intermediate_output/13_traceback_cluster_input_data/G1_traceback_cluster_sample.csv

Main outputs
------------
- ../results/intermediate_output/14_g1_traceback_clustering/; input to script 16 and 17.

Figure/table role
-----------------
../results/intermediate_output/14_g1_traceback_clustering/; input to script 16 and 17.

Runtime and data notes
----------------------
Reference runtime: 3 min 24.56 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd



HERE = Path(__file__).resolve().parent
os.chdir(HERE)
TRACEBACK_INPUT_DIR = Path("..") / "results" / "intermediate_output" / "13_traceback_cluster_input_data"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "14_g1_traceback_clustering"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"

SAMPLE_FILE = TRACEBACK_INPUT_DIR / "G1_traceback_cluster_sample.csv"
ORBIT_FILE = OUT_DIR / "G1_orbits_55Myr.npz"
LABELS_FILE = OUT_DIR / "G1_family_labels.csv"
REFERENCES_FILE = OUT_DIR / "G1_family_references.csv"
PROPERTIES_FILE = OUT_DIR / "G1_family_properties.csv"
CURVE_FILE = OUT_DIR / "G1_family_size_evolution.csv"
MEMBERS_FILE = OUT_DIR / "G1_family_members.csv"
MEMBER_SUMMARY_FILE = OUT_DIR / "G1_family_member_summary.csv"
PARAMS_FILE = OUT_DIR / "G1_clustering_parameters.json"

MAX_LOOKBACK_MYR = 55.0
STEP_MYR = 0.1
PARAMS = {
    "min_cluster_sizes": "5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30",
    "min_members": 5,
    "max_members": 120,
    "min_density_score": 1.8,
    "max_references": 40,
    "max_jaccard_overlap": 0.55,
    "assignment_overlap_floor": 0.38,
    "min_final_members": 6,
    "sigma_age_myr": 10.0,
    "second_choice_noise_fraction": 0.9,
    "cluster_selection_method": "eom",
    "min_weight_score": 0.25,
    "spatial_scale_pc": 400.0,
    "spatial_power": 0.0,
}



TRACEBACK_RO_KPC = 8.122
TRACEBACK_VO_KMS = 236.0
TRACEBACK_ZO_KPC = 0.0208
TRACEBACK_SOLAR_MOTION_UVW = (-11.1, 12.24, 7.25)
TRACEBACK_SOLAR_MOTION_LSR_STANDARD = (11.1, 12.24, 7.25)


def traceback_load_cluster_sample(path):
    table_path = Path(path)
    df = pd.read_csv(table_path)
    numeric_columns = [
        "id",
        "age_myr",
        "x_helio",
        "y_helio",
        "z_helio",
        "U",
        "V",
        "W",
        "U_lsr",
        "V_lsr",
        "W_lsr",
        "x_helio_err",
        "y_helio_err",
        "z_helio_err",
        "U_err",
        "V_err",
        "W_err",
        "dist_pc",
        "weight_score",
    ]
    for column in numeric_columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    if "id" not in df.columns:
        df["id"] = np.arange(1, len(df) + 1, dtype=np.int64)
    if "name" not in df.columns:
        df["name"] = df["display_name"] if "display_name" in df.columns else df["id"].map(lambda value: f"cluster_{value:04d}")
    if "display_name" not in df.columns:
        df["display_name"] = df["name"]
    if "source_catalog" not in df.columns:
        df["source_catalog"] = "unknown"
    if "family" not in df.columns:
        df["family"] = pd.NA
    if "age_myr" not in df.columns:
        df["age_myr"] = np.nan
    if "dist_pc" not in df.columns:
        df["dist_pc"] = np.sqrt(df["x_helio"] ** 2 + df["y_helio"] ** 2 + df["z_helio"] ** 2)

    if not {"U", "V", "W"}.issubset(df.columns):
        if {"U_lsr", "V_lsr", "W_lsr"}.issubset(df.columns):
            df["U"] = df["U_lsr"] - TRACEBACK_SOLAR_MOTION_LSR_STANDARD[0]
            df["V"] = df["V_lsr"] - TRACEBACK_SOLAR_MOTION_LSR_STANDARD[1]
            df["W"] = df["W_lsr"] - TRACEBACK_SOLAR_MOTION_LSR_STANDARD[2]
        else:
            raise ValueError("The sample table must contain U/V/W or U_lsr/V_lsr/W_lsr.")

    if not {"U_lsr", "V_lsr", "W_lsr"}.issubset(df.columns):
        df["U_lsr"] = df["U"] + TRACEBACK_SOLAR_MOTION_LSR_STANDARD[0]
        df["V_lsr"] = df["V"] + TRACEBACK_SOLAR_MOTION_LSR_STANDARD[1]
        df["W_lsr"] = df["W"] + TRACEBACK_SOLAR_MOTION_LSR_STANDARD[2]

    for column in ("x_helio_err", "y_helio_err", "z_helio_err", "U_err", "V_err", "W_err"):
        if column not in df.columns:
            df[column] = 0.0
        df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0.0)

    required = ["id", "display_name", "source_catalog", "age_myr", "x_helio", "y_helio", "z_helio", "U", "V", "W", "U_lsr", "V_lsr", "W_lsr"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"{table_path}  is missing columns: {missing}")
    return df


def traceback_build_time_grid(max_lookback_myr, step_myr):
    import astropy.units as _u

    lookback_myr = np.arange(0.0, max_lookback_myr + 0.5 * step_myr, step_myr, dtype=np.float64)
    return lookback_myr, (-lookback_myr) * _u.Myr


def traceback_cartesian_to_galactic_angles(x_pc, y_pc, z_pc):
    dist_pc = np.sqrt(x_pc**2 + y_pc**2 + z_pc**2)
    l_deg = np.degrees(np.arctan2(y_pc, x_pc)) % 360.0
    b_deg = np.degrees(np.arcsin(np.divide(z_pc, dist_pc, out=np.zeros_like(z_pc), where=dist_pc > 0.0)))
    return l_deg, b_deg, dist_pc / 1000.0


def traceback_integrate_orbits(sample_file, output_file, max_lookback_myr, step_myr, n_realizations=1, seed=42, sample_errors=False):
    from galpy.orbit import Orbit
    from galpy.potential import MWPotential2014

    df = traceback_load_cluster_sample(sample_file)
    lookback_myr, galpy_times = traceback_build_time_grid(max_lookback_myr, step_myr)
    means = df[["x_helio", "y_helio", "z_helio", "U", "V", "W"]].to_numpy(dtype=np.float64)
    if sample_errors:
        errs = df[["x_helio_err", "y_helio_err", "z_helio_err", "U_err", "V_err", "W_err"]].to_numpy(dtype=np.float64)
        sampled = means[None, :, :] + np.random.default_rng(seed).normal(size=(n_realizations, *means.shape)) * errs[None, :, :]
    else:
        sampled = np.repeat(means[None, :, :], n_realizations, axis=0)

    positions = np.empty((n_realizations, len(lookback_myr), len(df), 3), dtype=np.float32)
    for realization_index in range(n_realizations):
        print(f"integrating realization {realization_index + 1}/{n_realizations}")
        for cluster_index, (x_helio, y_helio, z_helio, u_vel, v_vel, w_vel) in enumerate(sampled[realization_index]):
            l_deg, b_deg, dist_kpc = traceback_cartesian_to_galactic_angles(
                np.asarray([x_helio]), np.asarray([y_helio]), np.asarray([z_helio])
            )
            orbit = Orbit(
                vxvv=[l_deg[0], b_deg[0], dist_kpc[0], u_vel, v_vel, w_vel],
                lb=True,
                uvw=True,
                ro=TRACEBACK_RO_KPC,
                vo=TRACEBACK_VO_KMS,
                zo=TRACEBACK_ZO_KPC,
                solarmotion=[-TRACEBACK_SOLAR_MOTION_UVW[0], TRACEBACK_SOLAR_MOTION_UVW[1], TRACEBACK_SOLAR_MOTION_UVW[2]],
            )
            orbit.integrate(galpy_times, MWPotential2014, method="odeint")
            positions[realization_index, :, cluster_index, 0] = -np.asarray(orbit.x(galpy_times, use_physical=True)) * 1000.0
            positions[realization_index, :, cluster_index, 1] = np.asarray(orbit.y(galpy_times, use_physical=True)) * 1000.0
            positions[realization_index, :, cluster_index, 2] = np.asarray(orbit.z(galpy_times, use_physical=True)) * 1000.0

    payload = {
        "positions_pc": positions,
        "positions_median_pc": np.median(positions, axis=0).astype(np.float32),
        "positions_std_pc": np.std(positions, axis=0).astype(np.float32),
        "lookback_myr": lookback_myr.astype(np.float32),
        "cluster_ids": df["id"].to_numpy(dtype=np.float64),
        "ages_myr": df["age_myr"].to_numpy(dtype=np.float32),
        "n_realizations": np.asarray([n_realizations], dtype=np.int32),
        "sample_errors": np.asarray([int(sample_errors)], dtype=np.int32),
        "random_seed": np.asarray([seed], dtype=np.int32),
    }
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **payload)
    print(f"saved {output_path}")
    return output_path


def traceback_gaussian_age_weights(lookback_myr, ages_myr, sigma_myr):
    weights = np.ones((len(lookback_myr), len(ages_myr)), dtype=np.float64)
    valid = np.isfinite(ages_myr)
    if np.any(valid):
        delta = lookback_myr[:, None] - ages_myr[None, valid]
        weights[:, valid] = np.exp(-0.5 * (delta / sigma_myr) ** 2)
    return weights


def traceback_jaccard(first, second):
    union = set(first) | set(second)
    return 0.0 if not union else len(set(first) & set(second)) / len(union)


def traceback_parse_sizes(value):
    return tuple(int(token.strip()) for token in str(value).split(",") if token.strip())


def traceback_group_name(label):
    return f"group_{label + 1:02d}" if label >= 0 else "none"


def traceback_compactness(coords, members, spatial_scale_pc, spatial_power):
    member_coords = coords[list(members)]
    center = np.median(member_coords, axis=0)
    radius = float(np.median(np.linalg.norm(member_coords - center[None, :], axis=1)))
    factor = (spatial_scale_pc / (spatial_scale_pc + max(radius, 1.0))) ** spatial_power
    return factor, radius


def traceback_collect_hdbscan_runs(
    positions_pc,
    lookback_myr,
    ages_myr,
    min_cluster_sizes,
    min_members,
    max_members,
    min_density_score,
    sigma_age_myr,
    cluster_selection_method,
    spatial_scale_pc,
    spatial_power,
):
    import hdbscan

    age_weights = traceback_gaussian_age_weights(lookback_myr, ages_myr, sigma_age_myr)
    candidates = []
    runs = []
    for time_index, coords in enumerate(positions_pc):
        if time_index % 50 == 0:
            print(f"clustering time step {time_index + 1}/{len(lookback_myr)}")
        for min_cluster_size in min_cluster_sizes:
            labels = hdbscan.HDBSCAN(
                min_cluster_size=int(min_cluster_size),
                cluster_selection_method=cluster_selection_method,
            ).fit_predict(coords)
            runs.append({"time_index": time_index, "min_cluster_size": int(min_cluster_size), "labels": labels.astype(np.int32, copy=False)})
            for label in np.unique(labels):
                if label < 0:
                    continue
                members = tuple(np.flatnonzero(labels == label).tolist())
                n_members = len(members)
                if not (min_members <= n_members <= max_members):
                    continue
                raw_score = float(age_weights[time_index, list(members)].sum())
                factor, median_radius_pc = traceback_compactness(coords, members, spatial_scale_pc, spatial_power)
                density_score = raw_score / np.sqrt(n_members) * factor
                if density_score >= min_density_score:
                    candidates.append(
                        {
                            "time_index": time_index,
                            "min_cluster_size": int(min_cluster_size),
                            "members": members,
                            "n_members": n_members,
                            "raw_score": raw_score,
                            "density_score": density_score,
                            "compactness_factor": factor,
                            "median_radius_pc": median_radius_pc,
                        }
                    )
    candidates.sort(key=lambda item: float(item["density_score"]), reverse=True)
    return candidates, runs, age_weights


def traceback_select_references(candidates, max_references, max_jaccard_overlap):
    references = []
    for candidate in candidates:
        if all(traceback_jaccard(candidate["members"], reference["members"]) < max_jaccard_overlap for reference in references):
            references.append(candidate)
        if len(references) >= max_references:
            break
    return references


def traceback_accumulate_scores(runs, references, age_weights, assignment_overlap_floor, n_clusters, positions_pc, spatial_scale_pc, spatial_power):
    reference_scores = np.zeros((len(references), n_clusters), dtype=np.float64)
    noise_scores = np.zeros(n_clusters, dtype=np.float64)
    noise_counts = np.zeros(n_clusters, dtype=np.int32)
    for run_index, run in enumerate(runs, start=1):
        if run_index % 500 == 0:
            print(f"accumulating run {run_index}/{len(runs)}")
        time_index = int(run["time_index"])
        labels = run["labels"]
        time_weight = age_weights[time_index]
        noise_mask = labels < 0
        noise_scores[noise_mask] += time_weight[noise_mask]
        noise_counts[noise_mask] += 1
        for label in np.unique(labels):
            if label < 0:
                continue
            members = tuple(np.flatnonzero(labels == label).tolist())
            overlaps = [traceback_jaccard(members, reference["members"]) for reference in references]
            if not overlaps:
                continue
            best_index = int(np.argmax(overlaps))
            best_overlap = float(overlaps[best_index])
            if best_overlap < assignment_overlap_floor:
                continue
            factor, _ = traceback_compactness(positions_pc[time_index], members, spatial_scale_pc, spatial_power)
            reference_scores[best_index, list(members)] += time_weight[list(members)] * best_overlap * factor
    return reference_scores, noise_scores, noise_counts.astype(np.float64) / max(len(runs), 1)


def traceback_choose_labels(reference_scores, noise_scores, noise_fraction, second_choice_noise_fraction, min_weight_score):
    n_clusters = noise_scores.shape[0]
    raw_labels = np.full(n_clusters, -1, dtype=np.int32)
    primary_scores = np.zeros(n_clusters, dtype=np.float64)
    secondary_scores = np.zeros(n_clusters, dtype=np.float64)
    weight_scores = np.zeros(n_clusters, dtype=np.float64)
    for cluster_index in range(n_clusters):
        family_scores = reference_scores[:, cluster_index] if reference_scores.size else np.empty(0, dtype=np.float64)
        total_score = float(noise_scores[cluster_index] + family_scores.sum())
        if len(family_scores) == 0 or np.allclose(family_scores, 0.0):
            continue
        best_family = int(np.argmax(family_scores))
        best_score = float(family_scores[best_family])
        sorted_scores = np.sort(family_scores)
        second_score = float(sorted_scores[-2]) if len(sorted_scores) >= 2 else 0.0
        chosen_label = best_family
        chosen_score = best_score
        if noise_scores[cluster_index] >= best_score and noise_fraction[cluster_index] < second_choice_noise_fraction and second_score > 0.0:
            chosen_label = int(np.argsort(family_scores)[-2])
            chosen_score = float(family_scores[chosen_label])
        elif noise_scores[cluster_index] > best_score:
            chosen_label = -1
            chosen_score = 0.0
        chosen_weight = chosen_score / total_score if total_score > 0.0 else 0.0
        if chosen_label >= 0 and chosen_weight < min_weight_score:
            chosen_label = -1
            chosen_score = 0.0
            chosen_weight = 0.0
        raw_labels[cluster_index] = chosen_label
        primary_scores[cluster_index] = chosen_score
        secondary_scores[cluster_index] = second_score
        weight_scores[cluster_index] = chosen_weight
    return raw_labels, primary_scores, secondary_scores, weight_scores


def traceback_enforce_min_members(raw_labels, min_final_members):
    labels = raw_labels.copy()
    for label in [label for label in np.unique(labels) if label >= 0]:
        if int(np.sum(labels == label)) < min_final_members:
            labels[labels == label] = -1
    return labels


def traceback_compact_positive_labels(labels):
    compact = np.full(labels.shape, -1, dtype=np.int32)
    for target, source in enumerate(sorted(int(label) for label in np.unique(labels) if label >= 0)):
        compact[labels == source] = target
    return compact


def traceback_write_reference_summary(references, df, lookback_myr, output_file):
    rows = []
    for rank, reference in enumerate(references, start=1):
        rows.append(
            {
                "reference_group": traceback_group_name(rank - 1),
                "rank": rank,
                "lookback_myr": float(lookback_myr[int(reference["time_index"])]),
                "min_cluster_size": int(reference["min_cluster_size"]),
                "n_members": int(reference["n_members"]),
                "density_score": float(reference["density_score"]),
                "raw_score": float(reference["raw_score"]),
                "median_radius_pc": float(reference.get("median_radius_pc", np.nan)),
                "member_names": json.dumps(df.iloc[list(reference["members"])]["display_name"].tolist(), ensure_ascii=False),
            }
        )
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_path, index=False, encoding="utf-8-sig")
    print(f"saved {output_path}")


def traceback_run_family_inference(sample_file, orbit_file, output_file, references_file, params):
    df = traceback_load_cluster_sample(sample_file)
    payload = np.load(orbit_file)
    positions_pc = payload["positions_median_pc"].astype(np.float64)
    lookback_myr = payload["lookback_myr"].astype(np.float64)
    candidates, runs, age_weights = traceback_collect_hdbscan_runs(
        positions_pc,
        lookback_myr,
        df["age_myr"].to_numpy(dtype=np.float64),
        traceback_parse_sizes(params["min_cluster_sizes"]),
        int(params["min_members"]),
        int(params["max_members"]),
        float(params["min_density_score"]),
        float(params["sigma_age_myr"]),
        str(params["cluster_selection_method"]),
        float(params["spatial_scale_pc"]),
        float(params["spatial_power"]),
    )
    references = traceback_select_references(candidates, int(params["max_references"]), float(params["max_jaccard_overlap"]))
    reference_scores, noise_scores, noise_fraction = traceback_accumulate_scores(
        runs,
        references,
        age_weights,
        float(params["assignment_overlap_floor"]),
        len(df),
        positions_pc,
        float(params["spatial_scale_pc"]),
        float(params["spatial_power"]),
    )
    raw_labels, primary_scores, secondary_scores, weight_scores = traceback_choose_labels(
        reference_scores,
        noise_scores,
        noise_fraction,
        float(params["second_choice_noise_fraction"]),
        float(params["min_weight_score"]),
    )
    final_labels = traceback_compact_positive_labels(traceback_enforce_min_members(raw_labels, int(params["min_final_members"])))
    output_df = df.copy()
    output_df["raw_cluster_label"] = raw_labels
    output_df["final_cluster_label"] = final_labels
    output_df["inferred_family"] = [traceback_group_name(label) for label in final_labels]
    output_df["primary_score"] = primary_scores
    output_df["secondary_score"] = secondary_scores
    output_df["noise_score"] = noise_scores
    output_df["noise_fraction"] = noise_fraction
    output_df["weight_score"] = weight_scores
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_csv(output_path, index=False, encoding="utf-8-sig")
    traceback_write_reference_summary(references, df, lookback_myr, references_file)
    print(output_df["inferred_family"].value_counts(dropna=False).sort_index().to_string())
    return output_df


def traceback_measure_family_properties(sample_file, orbit_file, labels_file, output_summary_file, output_curve_file):
    labels = pd.read_csv(labels_file)
    payload = np.load(orbit_file)
    lookback_myr = payload["lookback_myr"].astype(np.float64)
    positions = payload["positions_median_pc"].astype(np.float64)
    rows = []
    curve_rows = []
    for family in sorted([value for value in labels["inferred_family"].dropna().unique() if value != "none"]):
        idx = np.flatnonzero(labels["inferred_family"].astype(str).to_numpy() == family)
        family_positions = positions[:, idx, :]
        radii = []
        for time_index, time_value in enumerate(lookback_myr):
            center = np.median(family_positions[time_index], axis=0)
            radius = float(np.median(np.linalg.norm(family_positions[time_index] - center[None, :], axis=1)))
            radii.append(radius)
            curve_rows.append({"family": family, "lookback_myr": float(time_value), "median_radius_pc": radius})
        radii = np.asarray(radii)
        min_index = int(np.argmin(radii))
        rows.append(
            {
                "family": family,
                "n_members": int(len(idx)),
                "present_radius_pc": float(radii[0]),
                "compact_radius_pc": float(radii[min_index]),
                "compact_lookback_myr": float(lookback_myr[min_index]),
                "mean_age_myr": float(pd.to_numeric(labels.iloc[idx]["age_myr"], errors="coerce").mean()),
                "median_age_myr": float(pd.to_numeric(labels.iloc[idx]["age_myr"], errors="coerce").median()),
            }
        )
    summary_path = Path(output_summary_file)
    curve_path = Path(output_curve_file)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(summary_path, index=False, encoding="utf-8-sig")
    pd.DataFrame(curve_rows).to_csv(curve_path, index=False, encoding="utf-8-sig")
    print(f"saved {summary_path}")
    print(f"saved {curve_path}")


def traceback_export_family_members(labels_file, members_output_file, summary_output_file):
    labels = pd.read_csv(labels_file)
    members = labels.loc[labels["inferred_family"].astype(str) != "none"].copy()
    members_path = Path(members_output_file)
    summary_path = Path(summary_output_file)
    members_path.parent.mkdir(parents=True, exist_ok=True)
    members.to_csv(members_path, index=False, encoding="utf-8-sig")
    summary = (
        members.groupby("inferred_family", dropna=False)
        .agg(n_members=("display_name", "size"), mean_age_myr=("age_myr", "mean"), median_weight_score=("weight_score", "median"))
        .reset_index()
        .rename(columns={"inferred_family": "family"})
    )
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")
    print(f"saved {members_path}")
    print(f"saved {summary_path}")

def require_file(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path}")
    return path


def main() -> None:
    """Run Script 14 from validated inputs to the documented outputs."""
    parser = argparse.ArgumentParser(description='Run Script 14: G1 YSO traceback.')
    parser.add_argument("--recompute-orbit", action="store_true", help='Recompute the orbit NPZ even if it already exists.')
    parser.add_argument("--recompute-labels", action="store_true", help='Recompute family labels even if they already exist.')
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)
    require_file(SAMPLE_FILE)

    started = time.perf_counter()
    if args.recompute_orbit or not ORBIT_FILE.exists():
        traceback_integrate_orbits(
            sample_file=SAMPLE_FILE,
            output_file=ORBIT_FILE,
            max_lookback_myr=MAX_LOOKBACK_MYR,
            step_myr=STEP_MYR,
        )
    else:
        print(f"Reusing orbit file: {ORBIT_FILE}")

    if args.recompute_labels or not LABELS_FILE.exists():
        traceback_run_family_inference(SAMPLE_FILE, ORBIT_FILE, LABELS_FILE, REFERENCES_FILE, PARAMS)
    else:
        print(f"Reusing family labels: {LABELS_FILE}")

    traceback_measure_family_properties(SAMPLE_FILE, ORBIT_FILE, LABELS_FILE, PROPERTIES_FILE, CURVE_FILE)
    traceback_export_family_members(LABELS_FILE, MEMBERS_FILE, MEMBER_SUMMARY_FILE)
    PARAMS_FILE.write_text(json.dumps(PARAMS, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Saved G1 traceback outputs: {LABELS_FILE}, {MEMBERS_FILE}, {MEMBER_SUMMARY_FILE}")


def format_duration(seconds: float) -> str:
    seconds = max(0.0, float(seconds))
    hours, remainder = divmod(seconds, 3600.0)
    minutes, seconds = divmod(remainder, 60.0)
    if hours >= 1.0:
        return f"{int(hours)} h {int(minutes):02d} min {seconds:05.2f} s"
    if minutes >= 1.0:
        return f"{int(minutes)} min {seconds:05.2f} s"
    return f"{seconds:.2f} s"


def run_timed_main(func, script_file=None):
    """Run this timed entry point and report elapsed wall-clock time."""
    script_name = Path(script_file or __file__).name
    start = time.perf_counter()
    try:
        result = func()
    except SystemExit as exc:
        elapsed = time.perf_counter() - start
        if exc.code in (None, 0):
            print(f"[runtime] {script_name} completed in: {format_duration(elapsed)} ({elapsed:.2f} s)", flush=True)
        else:
            print(f"[runtime] {script_name} elapsed time before failure: {format_duration(elapsed)} ({elapsed:.2f} s)", flush=True)
        raise
    except BaseException:
        elapsed = time.perf_counter() - start
        print(f"[runtime] {script_name} elapsed time before failure: {format_duration(elapsed)} ({elapsed:.2f} s)", flush=True)
        raise
    elapsed = time.perf_counter() - start
    print(f"[runtime] {script_name} completed in: {format_duration(elapsed)} ({elapsed:.2f} s)", flush=True)
    return result


if __name__ == "__main__":
    run_timed_main(main, __file__)
