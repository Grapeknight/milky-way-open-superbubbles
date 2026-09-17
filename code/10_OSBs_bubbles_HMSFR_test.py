"""Script 10: OSBs bubbles HMSFR test

Purpose
-------
Count shell-associated Grade A/B bubbles and HMSFRs; run random geometry controls.

Method overview
---------------
1. Load the fitted OSB shells, Grade A/B bubble catalogue and Reid HMSFR sample.
2. Count observed shell associations and union coverage, then repeat the relevant
   measurements for randomized shell geometries.
3. Export target-level counts, control distributions, coverage summaries and association
   maps, including the configured surface-geometry robustness comparison.

Main inputs
-----------
- ../results/superbubble_final_fit_parameters.csv
- ../data/Bubbles.csv
- ../data/star_cluster_data/Reid2019_HMSFR.csv

Main outputs
------------
- ../results/figures/10_shell_assignment_xy.png, 10_shell_count_random_test_overview.png; feedback-tracer association figures.

Figure/table role
-----------------
../results/figures/10_shell_assignment_xy.png, 10_shell_count_random_test_overview.png; feedback-tracer association figures.

Runtime and data notes
----------------------
Reference runtime: 11 min 09.93 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ----------------------------------------------------------------------------
# Paths and output products. All paths are resolved relative to code/ and point directly into ../data or ../results/.
# ----------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "10_shell_bubbles_hmsfr_random_test"


def ensure_existing_file(path_value, label: str) -> Path:
    """Return an existing input path or raise a clear file-not-found error."""
    path = Path(path_value)
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    return path


SAMPLE_CSV = FINAL_DIR / "superbubble_final_fit_parameters.csv"
BUBBLE_CSV = DATA_DIR / "Bubbles.csv"
HMSFR_4KPC_CSV = DATA_DIR / "star_cluster_data" / "Reid2019_HMSFR.csv"
HMSFR_3KPC_CSV = OUT_DIR / "Reid2019_HMSFR_trigonometric_parallax_sample_within_3kpc.csv"

OUT_SUMMARY_CSV = OUT_DIR / "10_shell_count_random_test_summary.csv"
OUT_RANDOM_CSV = OUT_DIR / "10_shell_count_1000_random_trials.csv"
OUT_BUBBLE_MEMBERS_CSV = OUT_DIR / "10_small_bubble_shell_member_details.csv"
OUT_HMSFR_MEMBERS_CSV = OUT_DIR / "10_HMSFR_3kpc_shell_member_details.csv"
OUT_HMSFR_COVERAGE_CSV = OUT_DIR / "10_HMSFR_3kpc_shell_unique_coverage_random_test.csv"
OUT_BUBBLE_COVERAGE_CSV = OUT_DIR / "10_small_bubble_shell_unique_coverage_random_test.csv"
OUT_HMSFR_SURFACE_COVERAGE_CSV = OUT_DIR / "10_HMSFR_3kpc_shell_unique_coverage_surface_random_test.csv"
OUT_BUBBLE_SURFACE_COVERAGE_CSV = OUT_DIR / "10_small_bubble_shell_unique_coverage_surface_random_test.csv"
OUT_BUBBLE_UNION_MEMBERS_CSV = OUT_DIR / "10_small_bubble_shell_unique_assignment_details.csv"
OUT_BUBBLE_SURFACE_UNION_MEMBERS_CSV = OUT_DIR / "10_small_bubble_shell_unique_assignment_surface_details.csv"
OUT_HMSFR_XY_CSV = OUT_DIR / "10_HMSFR_3kpc_xy_shell_assignment_details.csv"
OUT_FIG = FINAL_FIG_DIR / "10_shell_count_random_test_overview.png"
OUT_HMSFR_XY_FIG = FINAL_FIG_DIR / "10_shell_assignment_xy.png"

# Main shell definition used in the paper. SHELL_MODE can temporarily switch
# shell_mask to the stricter u=1 surface test for robustness checks.
SHELL_U_MIN = 0.9
SHELL_U_MAX = 1.2
SHELL_MODE = "volume"
CYLINDER_HALF_HEIGHT_KPC = 0.025
DEFAULT_N_RANDOM = 1000
DEFAULT_RANDOM_CENTER_HALF_WIDTH_KPC = 3.0
DEFAULT_RANDOM_SEED = 20260523
DEFAULT_HMSFR_DIAMETER_PC = 200.0
ELLIPSOID_CAP_SECTION_FRAC = 0.5
MAX_OPENING_HIGHLIGHT_IDS = {7, 31, 32}
BUBBLE_Z_SLICES = [
    ("Z=[-0.04, 0.04] kpc", -0.04, 0.04),
    ("Z=[-0.20, -0.04] kpc", -0.20, -0.04),
    ("Z=[0.04, 0.20] kpc", 0.04, 0.20),
]


def configure_shell_mode(mode: str) -> None:
    """Set the shell-counting mode used by shell_mask.

    volume: finite shell volume, 0.9 <= u <= 1.2.
    surface: zero-thickness shell surface, u = 1.
    """
    normalized = str(mode).strip().lower()
    if normalized not in {"volume", "surface"}:
        raise ValueError(f"unknown shell mode: {mode}")
    global SHELL_MODE
    SHELL_MODE = normalized


def parse_grade_list(value: str) -> tuple[str, ...]:
    """Parse a comma-separated bubble Grade list from the CLI."""
    grades = tuple(part.strip().upper() for part in str(value).split(",") if part.strip())
    if not grades:
        raise argparse.ArgumentTypeError("at least one bubble grade is required")
    return grades


def coverage_summary_text(
    label: str,
    coverage_df: pd.DataFrame,
    total_col: str,
    count_col: str,
    percent_col: str,
) -> str:
    """Return a one-line observed-vs-random coverage summary for console output."""
    observed = coverage_df.loc[coverage_df["scenario"] == "observed"].iloc[0]
    random_percent = coverage_df.loc[coverage_df["scenario"] == "random", percent_col].to_numpy(dtype=float)
    random_median = float(np.nanmedian(random_percent))
    random_p16 = float(np.nanpercentile(random_percent, 16.0))
    random_p84 = float(np.nanpercentile(random_percent, 84.0))
    random_sigma = (random_p84 - random_p16) / 2.0
    return (
        f"{label}: {int(observed[count_col])}/{int(observed[total_col])} = "
        f"{float(observed[percent_col]):.1f}%; random median = "
        f"{random_median:.1f}+/-{random_sigma:.1f}%"
    )


def spherical_to_cartesian(distance_kpc, lon_deg, lat_deg) -> np.ndarray:
    lon = np.deg2rad(np.asarray(lon_deg, dtype=float))
    lat = np.deg2rad(np.asarray(lat_deg, dtype=float))
    dist = np.asarray(distance_kpc, dtype=float)
    cos_lat = np.cos(lat)
    return np.column_stack(
        [
            dist * cos_lat * np.cos(lon),
            dist * cos_lat * np.sin(lon),
            dist * np.sin(lat),
        ]
    )


def load_small_bubbles(path: Path, allowed_grades: tuple[str, ...] = ("A", "B")) -> pd.DataFrame:
    """Load the bubble catalogue, keep requested grades, remove duplicates, and add XYZ/radius."""
    df = pd.read_csv(path)
    required = ["ID", "l", "b", "physical_size", "best_distance"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"{path}  is missing required columns: {', '.join(missing)}")

    # Default paper sample is Grade A/B. Pass ("A", "B", "C") for comparison runs.
    grade_removed = 0
    if "Grade" in df.columns:
        before_grade = len(df)
        grade_values = tuple(str(grade).strip().upper() for grade in allowed_grades)
        # Original paper-sample filter, kept as a reference for the default:
        # df = df[df["Grade"].astype(str).str.strip().str.upper().isin(["A", "B"])].copy()
        df = df[df["Grade"].astype(str).str.strip().str.upper().isin(grade_values)].copy()
        grade_removed = before_grade - len(df)
        df.attrs["allowed_grades"] = grade_values

    before = len(df)
    df = df.drop_duplicates(subset=["ID", "l", "b", "physical_size", "best_distance"]).copy()
    df["was_duplicate_removed"] = False
    df.attrs["grade_c_removed"] = grade_removed
    if before != len(df):
        df.attrs["duplicate_rows_removed"] = before - len(df)
    else:
        df.attrs["duplicate_rows_removed"] = 0

    for col in ["l", "b", "physical_size", "best_distance"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[np.isfinite(df["l"]) & np.isfinite(df["b"]) & np.isfinite(df["best_distance"])].copy()
    xyz = spherical_to_cartesian(df["best_distance"], df["l"], df["b"])
    df["x_kpc"] = xyz[:, 0]
    df["y_kpc"] = xyz[:, 1]
    df["z_kpc"] = xyz[:, 2]
    df["radius_kpc"] = df["physical_size"] / 2000.0
    return df.reset_index(drop=True)


def save_hmsfr_3kpc_catalog(input_path: Path, output_path: Path) -> pd.DataFrame:
    """Load Reid2019 HMSFRs, keep sources within 3 kpc, and save the filtered catalogue."""
    df = pd.read_csv(input_path).copy()
    required = ["Name", "distance_kpc", "x_kpc", "y_kpc", "z_kpc"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"{input_path}  is missing required columns: {', '.join(missing)}")

    for col in ["distance_kpc", "x_kpc", "y_kpc", "z_kpc"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[np.isfinite(df["distance_kpc"]) & (df["distance_kpc"] <= 3.0)].copy()
    df = df[np.isfinite(df["x_kpc"]) & np.isfinite(df["y_kpc"]) & np.isfinite(df["z_kpc"])].copy()
    df = df.sort_values(["distance_kpc", "Name"]).reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False, encoding="utf-8-sig")
    return df


def rotate_points_to_local(points: np.ndarray, angle_deg: float) -> np.ndarray:
    """Rotate global Cartesian offsets into the local frame of one superbubble."""
    theta = np.deg2rad(-float(angle_deg))
    c, s = np.cos(theta), np.sin(theta)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    return np.asarray(points, dtype=float) @ rot.T


def geometry_kind(row: pd.Series) -> str:
    """Return the simplified geometry family used by the shell-intersection code."""
    return "cylinder" if str(row["shape"]) == "cylinder" else "ellipsoid_cap"


def half_ellipsoid_geometry(row: pd.Series):
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    center = np.array(
        [float(row["center_x_kpc"]), float(row["center_y_kpc"]), float(row["center_z_kpc"])],
        dtype=float,
    )
    axes = np.array(
        [float(row["a_radius_kpc"]), float(row["b_radius_kpc"]), max(float(row["c_radius_kpc"]), 1e-9)],
        dtype=float,
    )
    return center, axes, int(row["mark"])


def geometry_for_row(row: pd.Series):
    """Return (kind, cap-side mark, center, axes, method label) for one fitted superbubble."""
    kind = geometry_kind(row)
    if kind == "cylinder":
        center = np.array(
            [float(row["center_x_kpc"]), float(row["center_y_kpc"]), float(row["center_z_kpc"])],
            dtype=float,
        )
        axes = np.array(
            [float(row["a_radius_kpc"]), float(row["b_radius_kpc"]), CYLINDER_HALF_HEIGHT_KPC],
            dtype=float,
        )
        mark = int(row["mark"])
        mask_method = "cylinder_xy_shell_z25pc"
    else:
        # Geometry and shell-mask conventions used by this analysis stage.
        center, axes, mark = half_ellipsoid_geometry(row)
        mask_method = "half_ellipsoid_shell"
    return kind, mark, center, axes, mask_method


def extreme_norm_over_ellipsoid(xi_c: np.ndarray, e: np.ndarray, want: str) -> np.ndarray:
    """Extremize ||xi|| over axis-aligned ellipsoids in metric space.

    xi_c is the metric-space centre of each tracer sphere, and e is its
    metric-space semi-axis vector. This lets shell_mask decide whether a
    tracer sphere overlaps either the finite shell volume or the u=1 surface
    without sampling points on the sphere.
    """
    'Helper for the documented pipeline stage.'
    xi_c = np.asarray(xi_c, dtype=float)
    e = np.maximum(np.asarray(e, dtype=float), 1e-12)
    e2 = e ** 2
    ce = xi_c * e
    maxe2 = np.max(e2, axis=1)
    mine2 = np.min(e2, axis=1)
    uc = np.sqrt(np.sum(xi_c ** 2, axis=1))
    span = (uc + np.max(e, axis=1) + 1.0) ** 2

    def g(mu):
        return np.sum((ce / (mu[:, None] - e2)) ** 2, axis=1)

    def norm_at(mu):
        return np.sqrt(np.sum((mu[:, None] * xi_c / (mu[:, None] - e2)) ** 2, axis=1))

    if want == "max":
        lo = maxe2 + 1e-12 * (maxe2 + 1.0) + 1e-15
        hi = maxe2 + span + 1.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            gt = g(mid) > 1.0
            lo = np.where(gt, mid, lo)
            hi = np.where(gt, hi, mid)
        return norm_at(0.5 * (lo + hi))

    inside = np.sum((xi_c / e) ** 2, axis=1) <= 1.0
    hi = mine2 - 1e-12 * (mine2 + 1.0) - 1e-15
    lo = mine2 - span - 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        gt = g(mid) > 1.0
        hi = np.where(gt, mid, hi)
        lo = np.where(gt, lo, mid)
    return np.where(inside, 0.0, norm_at(0.5 * (lo + hi)))


def shell_mask(
    points: np.ndarray,
    row: pd.Series,
    center_xy: np.ndarray | None = None,
    angle_deg: float | None = None,
    radii: np.ndarray | None = None,
):
    """Evaluate the membership or shell-selection mask for this analysis stage."""
    kind, mark, center, axes, mask_method = geometry_for_row(row)
    center = np.asarray(center, dtype=float).copy()
    if center_xy is not None:
        center[:2] = np.asarray(center_xy, dtype=float)
    if angle_deg is None:
        angle_deg = float(row["angle_deg"])

    local = rotate_points_to_local(np.asarray(points, dtype=float) - center[None, :], angle_deg)
    X, Y, Z = local[:, 0], local[:, 1], local[:, 2]
    a = max(float(axes[0]), 1e-9)
    b = max(float(axes[1]), 1e-9)
    c = max(float(axes[2]), 1e-9)

    if kind == "cylinder":
        norm_radius = np.sqrt((X / a) ** 2 + (Y / b) ** 2)
        if radii is None:
            half_mask = np.abs(Z) <= CYLINDER_HALF_HEIGHT_KPC
            if SHELL_MODE == "surface":
                mask = half_mask & np.isclose(norm_radius, 1.0, rtol=0.0, atol=1e-8)
            else:
                mask = half_mask & (norm_radius >= SHELL_U_MIN) & (norm_radius <= SHELL_U_MAX)
        else:

            r = np.broadcast_to(np.asarray(radii, dtype=float), X.shape).astype(float)
            h = CYLINDER_HALF_HEIGHT_KPC
            z_star = np.clip(Z, -h, h)
            rho2 = r ** 2 - (z_star - Z) ** 2
            reach = rho2 >= 0.0
            rho = np.sqrt(np.maximum(rho2, 0.0))
            xi_c = np.stack([X / a, Y / b], axis=1)
            e = np.stack([rho / a, rho / b], axis=1)
            u_min = extreme_norm_over_ellipsoid(xi_c, e, "min")
            u_max = extreme_norm_over_ellipsoid(xi_c, e, "max")
            if SHELL_MODE == "surface":
                mask = reach & (u_min <= 1.0) & (u_max >= 1.0)
            else:
                mask = reach & (u_min <= SHELL_U_MAX) & (u_max >= SHELL_U_MIN)
    else:
        norm_radius = np.sqrt((X / a) ** 2 + (Y / b) ** 2 + (Z / c) ** 2)
        if radii is None:
            half_mask = Z >= 0.0 if mark == 2 else Z <= 0.0
            if SHELL_MODE == "surface":
                mask = half_mask & np.isclose(norm_radius, 1.0, rtol=0.0, atol=1e-8)
            else:
                mask = half_mask & (norm_radius >= SHELL_U_MIN) & (norm_radius <= SHELL_U_MAX)
        else:
            # Geometry and shell-mask conventions used by this analysis stage.
            r = np.broadcast_to(np.asarray(radii, dtype=float), X.shape).astype(float)
            xi_c = np.stack([X / a, Y / b, Z / c], axis=1)
            e = np.stack([r / a, r / b, r / c], axis=1)
            u_min = extreme_norm_over_ellipsoid(xi_c, e, "min")
            u_max = extreme_norm_over_ellipsoid(xi_c, e, "max")
            half_mask = (Z >= -r) if mark == 2 else (Z <= r)
            if SHELL_MODE == "surface":
                mask = half_mask & (u_min <= 1.0) & (u_max >= 1.0)
            else:
                mask = half_mask & (u_min <= SHELL_U_MAX) & (u_max >= SHELL_U_MIN)

    return mask, norm_radius, center, axes, mask_method


def sample_random_center_xy(rng: np.random.Generator, center_xy: np.ndarray, half_width_kpc: float) -> np.ndarray:
    """Draw one randomized XY center offset in a square box around the fitted center."""
    return np.asarray(center_xy, dtype=float) + rng.uniform(
        -float(half_width_kpc),
        float(half_width_kpc),
        size=2,
    )


def random_count_stats(prefix: str, observed_count: int, random_counts: np.ndarray) -> dict:
    """Summarize one per-superbubble random-count distribution."""
    random_counts = np.asarray(random_counts, dtype=float)
    mean_count = float(np.mean(random_counts))
    median_count = float(np.median(random_counts))
    return {
        f"observed_{prefix}_count": int(observed_count),
        f"random_{prefix}_mean_count": mean_count,
        f"random_{prefix}_median_count": median_count,
        f"random_{prefix}_std_count": float(np.std(random_counts, ddof=1)) if len(random_counts) > 1 else 0.0,
        f"random_{prefix}_p05_count": float(np.percentile(random_counts, 5)),
        f"random_{prefix}_p95_count": float(np.percentile(random_counts, 95)),
        f"random_{prefix}_max_count": int(np.max(random_counts)),
        f"random_{prefix}_trials_gt_observed": int(np.sum(random_counts > observed_count)),
        f"random_{prefix}_trials_ge_observed": int(np.sum(random_counts >= observed_count)),
        f"p_random_{prefix}_gt_observed": float(np.mean(random_counts > observed_count)),
        f"p_random_{prefix}_ge_observed": float(np.mean(random_counts >= observed_count)),
        f"observed_{prefix}_minus_random_median": float(observed_count - median_count),
        f"observed_{prefix}_over_random_mean": float(observed_count / mean_count) if mean_count > 0 else np.inf,
    }


def observed_members_for_catalog(
    row: pd.Series,
    catalog_df: pd.DataFrame,
    id_column: str,
    source_catalog: str,
    radii: np.ndarray | None = None,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, np.ndarray, str]:
    """Return objects from one catalogue that intersect one superbubble shell."""
    points = catalog_df[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(dtype=float)
    observed_mask, norm_radius, metric_center, axes, mask_method = shell_mask(points, row, radii=radii)
    members = catalog_df.loc[observed_mask].copy()
    members["sb_id"] = int(row["id"])
    members["source_catalog"] = source_catalog
    members["shell_u"] = norm_radius[observed_mask]
    members["shell_mode"] = SHELL_MODE
    members["shell_u_min"] = SHELL_U_MIN
    members["shell_u_max"] = SHELL_U_MAX
    if id_column not in members.columns:
        members[id_column] = ""
    return members, observed_mask, points, metric_center, axes, mask_method


def compute_one_target(
    row: pd.Series,
    bubble_df: pd.DataFrame,
    hmsfr_df: pd.DataFrame,
    rng: np.random.Generator,
    n_random: int,
    center_half_width_kpc: float,
    bubble_radii: np.ndarray | None = None,
    hmsfr_radius_kpc: float | None = None,
):
    """Compute observed and randomized counts for one superbubble."""
    bubble_members, bubble_observed_mask, bubble_points, metric_center, axes, mask_method = observed_members_for_catalog(
        row, bubble_df, "ID", "small_bubbles", radii=bubble_radii
    )
    hmsfr_members, hmsfr_observed_mask, hmsfr_points, _, _, _ = observed_members_for_catalog(
        row, hmsfr_df, "Name", "Reid2019_HMSFR_3kpc", radii=hmsfr_radius_kpc
    )

    req = float(np.sqrt(float(row["a_radius_kpc"]) * float(row["b_radius_kpc"])))
    random_rows = []
    bubble_random_counts = []
    hmsfr_random_counts = []
    total_random_counts = []
    for trial in range(1, n_random + 1):
        random_center_xy = sample_random_center_xy(rng, metric_center[:2], center_half_width_kpc)
        random_angle = float(rng.uniform(0.0, 360.0))
        bubble_random_mask, _, _, _, _ = shell_mask(bubble_points, row, center_xy=random_center_xy, angle_deg=random_angle, radii=bubble_radii)
        hmsfr_random_mask, _, _, _, _ = shell_mask(hmsfr_points, row, center_xy=random_center_xy, angle_deg=random_angle, radii=hmsfr_radius_kpc)
        bubble_count = int(bubble_random_mask.sum())
        hmsfr_count = int(hmsfr_random_mask.sum())
        total_count = bubble_count + hmsfr_count
        bubble_random_counts.append(bubble_count)
        hmsfr_random_counts.append(hmsfr_count)
        total_random_counts.append(total_count)
        random_rows.append(
            {
                "id": int(row["id"]),
                "trial": trial,
                "random_center_x_kpc": float(random_center_xy[0]),
                "random_center_y_kpc": float(random_center_xy[1]),
                "fixed_center_z_kpc": float(metric_center[2]),
                "random_pa_deg": random_angle,
                "random_shell_bubble_count": bubble_count,
                "random_shell_hmsfr_3kpc_count": hmsfr_count,
                "random_shell_total_feedback_count": total_count,
            }
        )

    bubble_observed_count = int(bubble_observed_mask.sum())
    hmsfr_observed_count = int(hmsfr_observed_mask.sum())
    total_observed_count = bubble_observed_count + hmsfr_observed_count
    summary = {
        "id": int(row["id"]),
        "final_parameter_estimator": str(row["final_parameter_estimator"]),
        "mark": int(row["mark"]),
        "shape": str(row["shape"]),
        "mask_method": mask_method,
        "shell_mode": SHELL_MODE,
        "bubble_count_mode": f"sphere_shell_{SHELL_MODE}_intersect" if bubble_radii is not None else f"point_center_in_{SHELL_MODE}",
        "hmsfr_count_mode": f"sphere_shell_{SHELL_MODE}_intersect_r{hmsfr_radius_kpc:.4f}kpc" if hmsfr_radius_kpc else f"point_center_in_{SHELL_MODE}",
        "metric_center_x_kpc": float(metric_center[0]),
        "metric_center_y_kpc": float(metric_center[1]),
        "metric_center_z_kpc": float(metric_center[2]),
        "axis_a_shell_kpc": float(axes[0]),
        "axis_b_shell_kpc": float(axes[1]),
        "axis_c_shell_kpc": float(axes[2]),
        "pa_deg": float(row["angle_deg"]),
        "R_eq_kpc": req,
        "shell_u_min": SHELL_U_MIN,
        "shell_u_max": SHELL_U_MAX,
        "random_n": int(n_random),
        "random_center_half_width_kpc": float(center_half_width_kpc),
        "random_center_distribution": "uniform_xy_box_xplusminus3_yplusminus3_kpc",
        "observed_shell_bubble_ids": ";".join(map(str, bubble_members["ID"].tolist())),
        "observed_shell_hmsfr_3kpc_names": ";".join(map(str, hmsfr_members["Name"].tolist())),
    }
    summary.update(random_count_stats("shell_bubble", bubble_observed_count, np.asarray(bubble_random_counts, dtype=float)))
    summary.update(random_count_stats("shell_hmsfr_3kpc", hmsfr_observed_count, np.asarray(hmsfr_random_counts, dtype=float)))
    summary.update(random_count_stats("shell_total_feedback", total_observed_count, np.asarray(total_random_counts, dtype=float)))
    return summary, pd.DataFrame(random_rows), bubble_members, hmsfr_members


def compute_hmsfr_union_coverage(
    sample_df: pd.DataFrame,
    hmsfr_df: pd.DataFrame,
    random_df: pd.DataFrame,
    n_random: int,
    hmsfr_radius_kpc: float | None = None,
) -> pd.DataFrame:
    """Compute unique HMSFR coverage by any superbubble for observed and randomized geometries."""
    points = hmsfr_df[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(dtype=float)
    total = int(len(hmsfr_df))

    observed_union = np.zeros(total, dtype=bool)
    for row in sample_df.itertuples(index=False):
        row_series = pd.Series(row._asdict())
        observed_union |= shell_mask(points, row_series, radii=hmsfr_radius_kpc)[0]

    rows = [
        {
            "scenario": "observed",
            "trial": 0,
            "hmsfr_total_count": total,
            "hmsfr_unique_shell_count": int(observed_union.sum()),
            "hmsfr_unique_shell_fraction": float(observed_union.mean()) if total else float("nan"),
            "hmsfr_unique_shell_percent": float(observed_union.mean() * 100.0) if total else float("nan"),
        }
    ]

    random_lookup = random_df.set_index(["id", "trial"])
    for trial in range(1, int(n_random) + 1):
        random_union = np.zeros(total, dtype=bool)
        for row in sample_df.itertuples(index=False):
            row_series = pd.Series(row._asdict())
            random_row = random_lookup.loc[(int(row_series["id"]), trial)]
            random_center_xy = np.array(
                [float(random_row["random_center_x_kpc"]), float(random_row["random_center_y_kpc"])],
                dtype=float,
            )
            random_union |= shell_mask(
                points,
                row_series,
                center_xy=random_center_xy,
                angle_deg=float(random_row["random_pa_deg"]),
                radii=hmsfr_radius_kpc,
            )[0]
        rows.append(
            {
                "scenario": "random",
                "trial": trial,
                "hmsfr_total_count": total,
                "hmsfr_unique_shell_count": int(random_union.sum()),
                "hmsfr_unique_shell_fraction": float(random_union.mean()) if total else float("nan"),
                "hmsfr_unique_shell_percent": float(random_union.mean() * 100.0) if total else float("nan"),
            }
        )

    return pd.DataFrame(rows)


def compute_bubble_union_coverage(
    sample_df: pd.DataFrame,
    bubble_df: pd.DataFrame,
    random_df: pd.DataFrame,
    n_random: int,
    bubble_radii: np.ndarray | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute unique bubble coverage by any superbubble and return per-bubble assignments."""
    points = bubble_df[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(dtype=float)
    total = int(len(bubble_df))
    observed_union = np.zeros(total, dtype=bool)
    containing_ids: list[list[int]] = [[] for _ in range(total)]
    for row in sample_df.itertuples(index=False):
        row_series = pd.Series(row._asdict())
        mask = shell_mask(points, row_series, radii=bubble_radii)[0]
        observed_union |= mask
        sb_id = int(row_series["id"])
        for idx in np.flatnonzero(mask):
            containing_ids[int(idx)].append(sb_id)

    member_df = bubble_df.copy()
    member_df["in_any_sb_shell"] = observed_union
    member_df["containing_sb_ids"] = [";".join(f"SB{v}" for v in ids) for ids in containing_ids]

    rows = [
        {
            "scenario": "observed",
            "trial": 0,
            "bubble_total_count": total,
            "bubble_unique_shell_count": int(observed_union.sum()),
            "bubble_unique_shell_fraction": float(observed_union.mean()) if total else float("nan"),
            "bubble_unique_shell_percent": float(observed_union.mean() * 100.0) if total else float("nan"),
        }
    ]
    random_lookup = random_df.set_index(["id", "trial"])
    for trial in range(1, int(n_random) + 1):
        random_union = np.zeros(total, dtype=bool)
        for row in sample_df.itertuples(index=False):
            row_series = pd.Series(row._asdict())
            random_row = random_lookup.loc[(int(row_series["id"]), trial)]
            random_center_xy = np.array(
                [float(random_row["random_center_x_kpc"]), float(random_row["random_center_y_kpc"])],
                dtype=float,
            )
            random_union |= shell_mask(
                points,
                row_series,
                center_xy=random_center_xy,
                angle_deg=float(random_row["random_pa_deg"]),
                radii=bubble_radii,
            )[0]
        rows.append(
            {
                "scenario": "random",
                "trial": trial,
                "bubble_total_count": total,
                "bubble_unique_shell_count": int(random_union.sum()),
                "bubble_unique_shell_fraction": float(random_union.mean()) if total else float("nan"),
                "bubble_unique_shell_percent": float(random_union.mean() * 100.0) if total else float("nan"),
            }
        )
    return pd.DataFrame(rows), member_df


def observed_union_mask(sample_df: pd.DataFrame, catalog_df: pd.DataFrame, radii=None) -> tuple[np.ndarray, list[list[int]]]:
    """Return a union membership mask and containing SB IDs for plotting."""
    points = catalog_df[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(dtype=float)
    union_mask = np.zeros(len(catalog_df), dtype=bool)
    containing_ids: list[list[int]] = [[] for _ in range(len(catalog_df))]
    for row in sample_df.itertuples(index=False):
        row_series = pd.Series(row._asdict())
        mask = shell_mask(points, row_series, radii=radii)[0]
        union_mask |= mask
        sb_id = int(row_series["id"])
        for idx in np.flatnonzero(mask):
            containing_ids[int(idx)].append(sb_id)
    return union_mask, containing_ids


def ellipse_xy_points(center_xy: np.ndarray, axes_xy: np.ndarray, angle_deg: float, scale: float, n: int = 240) -> np.ndarray:
    """Sample points on one rotated XY ellipse."""
    theta = np.linspace(0.0, 2.0 * np.pi, n)
    local = np.column_stack(
        [
            float(axes_xy[0]) * float(scale) * np.cos(theta),
            float(axes_xy[1]) * float(scale) * np.sin(theta),
        ]
    )
    angle = np.deg2rad(float(angle_deg))
    rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]], dtype=float)
    return local @ rot.T + np.asarray(center_xy, dtype=float)


def compromise_section_scale(row: pd.Series) -> float:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    if str(row.get("shape", "")).lower() != "ellipsoid":
        return 1.0
    c = float(row["c_radius_kpc"])
    if not np.isfinite(c) or c <= 0:
        return 1.0
    cz = float(row["center_z_kpc"])
    z_base = float(row["xy_plane_z_kpc"])
    mark = int(row["mark"])
    z_apex = cz - c if mark == 1 else cz + c
    z_sec = z_base + ELLIPSOID_CAP_SECTION_FRAC * (z_apex - z_base)
    k_base = np.sqrt(max(0.0, 1.0 - ((z_base - cz) / c) ** 2))
    k_sec = np.sqrt(max(0.0, 1.0 - ((z_sec - cz) / c) ** 2))
    return float(k_sec / k_base) if k_base > 0 else 1.0


def compromise_xy_ellipse_geometry(row: pd.Series):
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    required = [
        "xy_plane_center_x_kpc", "xy_plane_center_y_kpc",
        "xy_plane_a_kpc", "xy_plane_b_kpc", "xy_plane_angle_deg",
    ]
    if all(col in row.index for col in required):
        scale = compromise_section_scale(row)
        center_xy = np.array([float(row["xy_plane_center_x_kpc"]), float(row["xy_plane_center_y_kpc"])])
        axes_xy = np.array([float(row["xy_plane_a_kpc"]) * scale, float(row["xy_plane_b_kpc"]) * scale])
        angle = float(row["xy_plane_angle_deg"])
        return center_xy, axes_xy, angle
    _, _, center, axes, _ = shell_mask(np.zeros((1, 3), dtype=float), row)
    return center[:2], axes[:2], float(row["angle_deg"])


def max_opening_xy_ellipse_geometry(row: pd.Series):
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    required = [
        "xy_plane_center_x_kpc", "xy_plane_center_y_kpc",
        "xy_plane_a_kpc", "xy_plane_b_kpc", "xy_plane_angle_deg",
    ]
    if not all(col in row.index for col in required):
        return None
    center_xy = np.array([float(row["xy_plane_center_x_kpc"]), float(row["xy_plane_center_y_kpc"])])
    axes_xy = np.array([float(row["xy_plane_a_kpc"]), float(row["xy_plane_b_kpc"])])
    angle = float(row["xy_plane_angle_deg"])
    return center_xy, axes_xy, angle


def _draw_sb_ellipses(
    ax,
    sample_df: pd.DataFrame,
    label_fs: float = 7.0,
    draw_max_opening_highlights: bool = False,
):
    """Construct or evaluate the geometry convention used by this pipeline stage."""
    for row in sample_df.itertuples(index=False):
        row_series = pd.Series(row._asdict())
        sb_id = int(row_series["id"])
        center_xy, axes_xy, angle = compromise_xy_ellipse_geometry(row_series)
        ellipse = ellipse_xy_points(center_xy, axes_xy, angle, 1.0)
        ax.plot(
            ellipse[:, 0],
            ellipse[:, 1],
            color="black",
            linewidth=2.2,
            linestyle="-",
            alpha=1.0,
        )
        if (
            draw_max_opening_highlights
            and
            str(row_series.get("shape", "")).lower() == "ellipsoid"
            and sb_id in MAX_OPENING_HIGHLIGHT_IDS
        ):
            max_opening = max_opening_xy_ellipse_geometry(row_series)
            if max_opening is not None:
                open_center_xy, open_axes_xy, open_angle = max_opening
                open_ellipse = ellipse_xy_points(open_center_xy, open_axes_xy, open_angle, 1.0)
                ax.plot(
                    open_ellipse[:, 0],
                    open_ellipse[:, 1],
                    color="0.45",
                    linewidth=2.4,
                    linestyle="--",
                    alpha=0.95,
                    zorder=4,
                )
        ax.text(
            float(center_xy[0]),
            float(center_xy[1]),
            f"{sb_id}",
            fontsize=label_fs,
            fontweight="bold",
            color="0.1",
            ha="center",
            va="center",
            alpha=0.95,
            bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.70),
            zorder=8,
        )


def save_hmsfr_xy_membership_plot(
    sample_df: pd.DataFrame,
    hmsfr_df: pd.DataFrame,
    output_path: Path,
    csv_path: Path,
    hmsfr_radius_kpc: float | None = None,
    bubble_df: pd.DataFrame | None = None,
    bubble_radii: np.ndarray | None = None,
):
    """Save the four-panel XY assignment figure and the HMSFR membership CSV.

    The caller already filters bubble grades, so this function plots exactly the
    bubble sample passed to it.
    """
    from matplotlib.patches import Circle

    # Geometry and shell-mask conventions used by this analysis stage.
    hmsfr_in_shell, hmsfr_ids = observed_union_mask(sample_df, hmsfr_df, radii=hmsfr_radius_kpc)
    hmsfr_plot = hmsfr_df.copy()
    hmsfr_plot["in_any_sb_shell"] = hmsfr_in_shell
    hmsfr_plot["containing_sb_ids"] = [";".join(f"SB{v}" for v in ids) for ids in hmsfr_ids]
    hmsfr_plot.to_csv(csv_path, index=False, encoding="utf-8-sig")

    # Bubble shell membership for the selected Grade sample.
    bub_in_shell, _ = observed_union_mask(sample_df, bubble_df, radii=bubble_radii)
    bub_x = bubble_df["x_kpc"].to_numpy(dtype=float)
    bub_y = bubble_df["y_kpc"].to_numpy(dtype=float)
    bub_r = (
        np.asarray(bubble_radii, dtype=float)
        if bubble_radii is not None
        else bubble_df["radius_kpc"].to_numpy(dtype=float)
    )
    bub_z = bubble_df["z_kpc"].to_numpy(dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(17.2, 17.6), constrained_layout=False)
    ax_b_mid, ax_h = axes[0]
    ax_b_low, ax_b_up = axes[1]


    xlim = (-3.0, 3.0)
    ylim = (-3.0, 3.0)

    label_fs = 18
    title_fs = 18
    tick_fs = 14
    legend_fs = 18
    in_color = "#D62728"
    out_color = "#7F7F7F"

    def draw_bubble_slice(
        ax,
        title: str,
        z_lo: float,
        z_hi: float,
    ):
        in_slice = (bub_z >= z_lo) & (bub_z <= z_hi)
        _draw_sb_ellipses(ax, sample_df, label_fs=8.5, draw_max_opening_highlights=True)
        for xi, yi, ri, ins in zip(bub_x[in_slice], bub_y[in_slice], bub_r[in_slice], bub_in_shell[in_slice]):
            ax.add_patch(
                Circle(
                    (xi, yi),
                    ri,
                    facecolor="none",
                    edgecolor=in_color if ins else out_color,
                    linewidth=1.8 if ins else 1.25,
                    alpha=0.95 if ins else 0.50,
                    zorder=6 if ins else 5,
                )
            )
        n_in = int(bub_in_shell[in_slice].sum())
        n_all = int(in_slice.sum())
        pct = 100.0 * n_in / n_all if n_all else float("nan")
        ax.set_title(
            f"Bubbles; {title}\n{n_in}/{n_all} = {pct:.1f}% associated with shells",
            fontsize=title_fs,
        )


    draw_bubble_slice(ax_b_mid, *BUBBLE_Z_SLICES[0])
    draw_bubble_slice(ax_b_low, *BUBBLE_Z_SLICES[1])
    draw_bubble_slice(ax_b_up, *BUBBLE_Z_SLICES[2])


    _draw_sb_ellipses(ax_h, sample_df, label_fs=8.5)
    outside = hmsfr_plot[~hmsfr_plot["in_any_sb_shell"]]
    inside = hmsfr_plot[hmsfr_plot["in_any_sb_shell"]]
    ax_h.scatter(outside["x_kpc"], outside["y_kpc"], s=70, facecolors="none", edgecolors=out_color,
                 linewidth=1.3, zorder=5)
    ax_h.scatter(inside["x_kpc"], inside["y_kpc"], s=95, color=in_color, edgecolor="0.15",
                 linewidth=0.7, zorder=6)
    ax_h.scatter([0.0], [0.0], marker="*", s=200, color="gold", edgecolor="0.2", linewidth=0.8, zorder=7)
    ax_h.set_title(
        f"HMSFR; $d \\leq$ 3 kpc\n{len(inside)}/{len(hmsfr_plot)} = {100.0 * hmsfr_in_shell.mean():.1f}% associated with shells",
        fontsize=title_fs,
    )

    for ax in (ax_b_mid, ax_h, ax_b_low, ax_b_up):
        ax.axhline(0.0, color="0.82", linewidth=0.8)
        ax.axvline(0.0, color="0.82", linewidth=0.8)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.tick_params(axis="both", labelsize=tick_fs)
        ax.grid(alpha=0.2, linestyle=":")
        for spine in ax.spines.values():
            spine.set_linewidth(1.1)

    legend_handles = [
        plt.Line2D([0], [0], marker="o", linestyle="none", markersize=13,
                   markerfacecolor="none", markeredgecolor=in_color, markeredgewidth=1.4,
                   label="Bubbles associated with shells"),
        plt.Line2D([0], [0], marker="o", linestyle="none", markersize=13,
                   markerfacecolor="none", markeredgecolor=out_color, markeredgewidth=1.0,
                   label="Bubbles not associated"),
        plt.Line2D([0], [0], marker="o", linestyle="none", markersize=12,
                   markerfacecolor=in_color, markeredgecolor="0.15", markeredgewidth=0.7,
                   label="HMSFR associated with shells"),
        plt.Line2D([0], [0], marker="o", linestyle="none", markersize=12,
                   markerfacecolor="none", markeredgecolor=out_color, markeredgewidth=1.3,
                   label="HMSFR not associated"),
        plt.Line2D([0], [0], marker="*", linestyle="none", markersize=18,
                   markerfacecolor="gold", markeredgecolor="0.2", markeredgewidth=0.8,
                   label="Sun"),
        plt.Line2D([0], [0], color="black", linewidth=2.2, label="Open superbubbles"),
    ]
    fig.subplots_adjust(left=0.075, right=0.985, top=0.940, bottom=0.190, wspace=0.015, hspace=0.19)
    fig.supxlabel("X [kpc]", fontsize=label_fs, y=0.115)
    fig.supylabel("Y [kpc]", fontsize=label_fs, x=0.028)
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=3,
        fontsize=legend_fs,
        framealpha=0.92,
        bbox_to_anchor=(0.5, 0.018),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_total_figure(
    summary_df: pd.DataFrame,
    random_df: pd.DataFrame,
    hmsfr_coverage_df: pd.DataFrame,
    bubble_coverage_df: pd.DataFrame,
    output_path: Path,
    bubble_sample_label: str = "A/B Bubbles",
):
    """Save the six-panel observed-vs-random overview figure."""
    ids = summary_df["id"].to_numpy(dtype=int)
    labels = [str(v) for v in ids]
    staggered_labels = [label if i % 2 == 0 else f"\n{label}" for i, label in enumerate(labels)]
    x = np.arange(len(summary_df))

    label_fs = 23
    title_fs = 24
    tick_fs = 18
    legend_fs = 18

    fig, axes = plt.subplots(2, 3, figsize=(22.0, 12.2), constrained_layout=True)

    def style_count_xaxis(ax):
        ax.set_xticks(x)
        ax.set_xticklabels(staggered_labels, rotation=0, ha="center", fontsize=tick_fs)
        ax.set_xlabel("SB ID", fontsize=label_fs)
        ax.tick_params(axis="y", labelsize=tick_fs)
        ax.tick_params(axis="x", labelsize=tick_fs, pad=2)
        ax.set_xlim(-0.7, len(summary_df) - 0.3)

    ax = axes[0, 0]
    bubble_obs = summary_df["observed_shell_bubble_count"].to_numpy(dtype=float)
    hmsfr_obs = summary_df["observed_shell_hmsfr_3kpc_count"].to_numpy(dtype=float)
    ax.bar(x, bubble_obs, color="#377EB8", edgecolor="0.25", linewidth=0.5, label=bubble_sample_label)
    ax.bar(x, hmsfr_obs, bottom=bubble_obs, color="#E41A1C", edgecolor="0.25", linewidth=0.5, label="HMSFR")
    style_count_xaxis(ax)
    ax.set_ylim(0, float(np.max(bubble_obs + hmsfr_obs)) * 1.15 + 0.5)
    ax.set_ylabel("Count", fontsize=label_fs)
    ax.set_title("Observed counts", fontsize=title_fs)
    ax.legend(fontsize=legend_fs)
    ax.grid(axis="y", alpha=0.25)

    def observed_vs_random_panel(ax, prefix: str, title: str, color: str, ylim: tuple[float, float] | None = None):
        median = summary_df[f"random_{prefix}_median_count"]
        p05 = summary_df[f"random_{prefix}_p05_count"]
        p95 = summary_df[f"random_{prefix}_p95_count"]
        observed = summary_df[f"observed_{prefix}_count"]
        yerr_low = median - p05
        yerr_high = p95 - median
        ax.errorbar(
            x,
            median,
            yerr=[yerr_low, yerr_high],
            fmt="o",
            color="0.25",
            ecolor="0.55",
            elinewidth=1.4,
            capsize=3.0,
            markersize=6,
            label="Random median (5-95%)",
        )
        ax.scatter(x, observed, color=color, s=70, zorder=3, edgecolor="0.15", linewidth=0.5, label="Observed")
        style_count_xaxis(ax)
        top = float(np.nanmax(np.r_[p95.to_numpy(dtype=float), observed.to_numpy(dtype=float)]))
        if ylim is None:
            ax.set_ylim(-0.5, top * 1.12 + 0.5)
        else:
            ax.set_ylim(*ylim)
        ax.set_ylabel("Count", fontsize=label_fs)
        ax.set_title(title, fontsize=title_fs)
        ax.legend(fontsize=legend_fs)
        ax.grid(axis="y", alpha=0.25)

    def probability_panel(ax, prefix: str, title: str):
        p_values = summary_df[f"p_random_{prefix}_gt_observed"].to_numpy(dtype=float)
        colors = np.where(p_values <= 0.05, "#1B9E77", "#7570B3")
        ax.bar(x, p_values, color=colors, edgecolor="0.25", linewidth=0.6)
        ax.axhline(0.05, color="crimson", linestyle="--", linewidth=1.3, label="0.05")
        ax.axhline(0.20, color="0.25", linestyle=":", linewidth=1.3, label="0.20")
        style_count_xaxis(ax)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel("P", fontsize=label_fs)
        ax.set_title(title, fontsize=title_fs)
        ax.legend(fontsize=legend_fs)
        ax.grid(axis="y", alpha=0.25)

    observed_vs_random_panel(axes[0, 1], "shell_bubble", f"{bubble_sample_label}: obs. vs rand.", "#377EB8", ylim=(-0.5, 60))
    observed_vs_random_panel(axes[0, 2], "shell_hmsfr_3kpc", "HMSFR: obs. vs rand.", "#E41A1C")
    probability_panel(axes[1, 0], "shell_bubble", f"{bubble_sample_label}: p(random > obs.)")

    def coverage_hist(ax, observed_percent, random_percent, obs_color, xlabel, title):
        lower = float(np.nanmin(np.r_[random_percent, observed_percent]))
        upper = float(np.nanmax(np.r_[random_percent, observed_percent]))
        margin = max((upper - lower) * 0.12, 0.5)
        bins = np.linspace(max(0.0, lower - margin), upper + margin, 18)
        random_median = float(np.nanmedian(random_percent))
        random_p16 = float(np.nanpercentile(random_percent, 16.0))
        random_p84 = float(np.nanpercentile(random_percent, 84.0))
        random_sigma = (random_p84 - random_p16) / 2.0
        ax.hist(random_percent, bins=bins, color="0.68", edgecolor="white", label="Random")
        ax.axvline(observed_percent, color=obs_color, linewidth=2.4, label=f"Observed={observed_percent:.1f}%")
        ax.axvline(
            random_median,
            color="0.15",
            linestyle="--",
            linewidth=1.6,
            label=f"Random={random_median:.1f}±{random_sigma:.1f}%",
        )
        ax.set_xlim(bins[0], bins[-1])
        ax.tick_params(axis="both", labelsize=tick_fs)
        ax.set_xlabel(xlabel, fontsize=label_fs)
        ax.set_ylabel("Trials", fontsize=label_fs)
        ax.set_title(title, fontsize=title_fs)
        ax.legend(fontsize=legend_fs)
        ax.grid(axis="y", alpha=0.25)

    observed_percent = float(hmsfr_coverage_df.loc[hmsfr_coverage_df["scenario"] == "observed", "hmsfr_unique_shell_percent"].iloc[0])
    random_percent = hmsfr_coverage_df.loc[hmsfr_coverage_df["scenario"] == "random", "hmsfr_unique_shell_percent"].to_numpy(dtype=float)
    coverage_hist(
        axes[1, 2], observed_percent, random_percent, "#E41A1C",
        "Coverage [%]", "HMSFR coverage",
    )

    bubble_observed_percent = float(bubble_coverage_df.loc[bubble_coverage_df["scenario"] == "observed", "bubble_unique_shell_percent"].iloc[0])
    bubble_random_percent = bubble_coverage_df.loc[
        bubble_coverage_df["scenario"] == "random", "bubble_unique_shell_percent"
    ].to_numpy(dtype=float)
    coverage_hist(
        axes[1, 1], bubble_observed_percent, bubble_random_percent, "#377EB8",
        "Coverage [%]", f"{bubble_sample_label} coverage",
    )

    for ax in axes.ravel():
        for spine in ax.spines.values():
            spine.set_linewidth(1.1)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the command-line interface for script 10."""
    parser = argparse.ArgumentParser(
        description=(
            "Count small bubbles and 3 kpc HMSFRs associated with fitted open-superbubble shells, "
            "then compare the observed counts with spatially randomized controls."
        )
    )
    parser.add_argument("--sample-csv", default=str(SAMPLE_CSV))
    parser.add_argument("--bubble-csv", default=str(BUBBLE_CSV))
    parser.add_argument("--hmsfr-4kpc-csv", default=str(HMSFR_4KPC_CSV))
    parser.add_argument("--hmsfr-3kpc-output-csv", default=str(HMSFR_3KPC_CSV))
    parser.add_argument("--n-random", type=int, default=DEFAULT_N_RANDOM)
    parser.add_argument("--random-seed", type=int, default=DEFAULT_RANDOM_SEED)
    parser.add_argument("--center-half-width-kpc", type=float, default=DEFAULT_RANDOM_CENTER_HALF_WIDTH_KPC)
    parser.add_argument(
        "--bubble-grades",
        type=parse_grade_list,
        default=("A", "B"),
        help="Comma-separated bubble grades to include. Default: A,B.",
    )
    parser.add_argument(
        "--shell-mode",
        choices=["volume", "surface"],
        default="volume",
        help="volume: finite shell 0.9<=u<=1.2; surface: zero-thickness u=1 surface.",
    )
    parser.add_argument(
        "--skip-figures",
        action="store_true",
        help="Write CSV outputs without rendering figures.",
    )
    parser.add_argument(
        "--bubble-mode",
        choices=["point", "sphere"],
        default="sphere",
        help="sphere: use radius=physical_size/2 for each bubble; point: use the bubble centre only.",
    )
    parser.add_argument(
        "--hmsfr-diameter-pc",
        type=float,
        default=DEFAULT_HMSFR_DIAMETER_PC,
        help="HMSFR tracer diameter in pc. Use 0 for a centre-point criterion. Default: 200 pc.",
    )
    return parser


def compute_surface_robustness_outputs(
    sample_df: pd.DataFrame,
    bubble_df: pd.DataFrame,
    hmsfr_df: pd.DataFrame,
    random_df: pd.DataFrame,
    n_random: int,
    bubble_radii: np.ndarray | None,
    hmsfr_radius_kpc: float | None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute u=1 surface coverage using the same randomized geometries as the main run."""
    previous_mode = SHELL_MODE
    configure_shell_mode("surface")
    try:
        surface_hmsfr_coverage = compute_hmsfr_union_coverage(
            sample_df,
            hmsfr_df,
            random_df,
            n_random,
            hmsfr_radius_kpc=hmsfr_radius_kpc,
        )
        surface_bubble_coverage, surface_bubble_members = compute_bubble_union_coverage(
            sample_df,
            bubble_df,
            random_df,
            n_random,
            bubble_radii=bubble_radii,
        )
    finally:
        configure_shell_mode(previous_mode)
    return surface_hmsfr_coverage, surface_bubble_coverage, surface_bubble_members


def print_coverage_pair(title: str, bubble_coverage_df: pd.DataFrame, hmsfr_coverage_df: pd.DataFrame) -> None:
    """Print bubble and HMSFR coverage summaries for one shell definition."""
    print(title)
    print(
        "  "
        + coverage_summary_text(
            "Bubbles",
            bubble_coverage_df,
            "bubble_total_count",
            "bubble_unique_shell_count",
            "bubble_unique_shell_percent",
        )
    )
    print(
        "  "
        + coverage_summary_text(
            "HMSFR",
            hmsfr_coverage_df,
            "hmsfr_total_count",
            "hmsfr_unique_shell_count",
            "hmsfr_unique_shell_percent",
        )
    )


def main():
    """Run Script 10 from validated inputs to the documented outputs."""
    parser = build_arg_parser()
    args = parser.parse_args()
    configure_shell_mode(args.shell_mode)

    sample_path = ensure_existing_file(args.sample_csv, "publication OSB sample table")
    bubble_path = ensure_existing_file(args.bubble_csv, "small-bubble table")
    hmsfr_4kpc_path = ensure_existing_file(args.hmsfr_4kpc_csv, "Reid2019 HMSFR 4 kpc table")
    hmsfr_3kpc_output_path = Path(args.hmsfr_3kpc_output_csv)
    sample_df = pd.read_csv(sample_path, encoding="utf-8-sig").sort_values("id").reset_index(drop=True)
    bubble_df = load_small_bubbles(bubble_path, allowed_grades=args.bubble_grades)
    hmsfr_df = save_hmsfr_3kpc_catalog(hmsfr_4kpc_path, hmsfr_3kpc_output_path)

    use_sphere = args.bubble_mode == "sphere"
    bubble_radii = bubble_df["radius_kpc"].to_numpy(dtype=float) if use_sphere else None
    hmsfr_sphere = float(args.hmsfr_diameter_pc) > 0.0
    hmsfr_radius_kpc = (float(args.hmsfr_diameter_pc) / 2.0 / 1000.0) if hmsfr_sphere else None

    bubble_label = (
        "small bubble = sphere (radius = physical_size/2)" if use_sphere else "small bubble = center point"
    )
    hmsfr_label = (
        f"HMSFR = sphere (diameter = {args.hmsfr_diameter_pc:g} pc)" if hmsfr_sphere else "HMSFR = center point"
    )
    shell_label = "0.9<=u<=1.2 shell volume" if SHELL_MODE == "volume" else "u=1 shell surface"
    mode_label = f"{bubble_label}; {hmsfr_label}; bubble grades={','.join(args.bubble_grades)}; counted if it intersects the {shell_label}"
    bubble_sample_label = "/".join(args.bubble_grades) + " Bubbles"

    rng = np.random.default_rng(int(args.random_seed))
    summaries = []
    random_tables = []
    bubble_member_tables = []
    hmsfr_member_tables = []
    for row in sample_df.itertuples(index=False):
        row_series = pd.Series(row._asdict())
        summary, random_detail, bubble_members, hmsfr_members = compute_one_target(
            row_series,
            bubble_df,
            hmsfr_df,
            rng,
            int(args.n_random),
            float(args.center_half_width_kpc),
            bubble_radii=bubble_radii,
            hmsfr_radius_kpc=hmsfr_radius_kpc,
        )
        summaries.append(summary)
        random_tables.append(random_detail)
        bubble_member_tables.append(bubble_members)
        hmsfr_member_tables.append(hmsfr_members)

    summary_df = pd.DataFrame(summaries).sort_values("id").reset_index(drop=True)
    random_df = pd.concat(random_tables, ignore_index=True)
    bubble_members_df = pd.concat(bubble_member_tables, ignore_index=True) if bubble_member_tables else pd.DataFrame()
    hmsfr_members_df = pd.concat(hmsfr_member_tables, ignore_index=True) if hmsfr_member_tables else pd.DataFrame()
    coverage_df = compute_hmsfr_union_coverage(sample_df, hmsfr_df, random_df, int(args.n_random), hmsfr_radius_kpc=hmsfr_radius_kpc)
    bubble_coverage_df, bubble_union_members_df = compute_bubble_union_coverage(
        sample_df, bubble_df, random_df, int(args.n_random), bubble_radii=bubble_radii
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(OUT_SUMMARY_CSV, index=False, encoding="utf-8-sig")
    random_df.to_csv(OUT_RANDOM_CSV, index=False, encoding="utf-8-sig")
    bubble_members_df.to_csv(OUT_BUBBLE_MEMBERS_CSV, index=False, encoding="utf-8-sig")
    hmsfr_members_df.to_csv(OUT_HMSFR_MEMBERS_CSV, index=False, encoding="utf-8-sig")
    coverage_df.to_csv(OUT_HMSFR_COVERAGE_CSV, index=False, encoding="utf-8-sig")
    bubble_coverage_df.to_csv(OUT_BUBBLE_COVERAGE_CSV, index=False, encoding="utf-8-sig")
    bubble_union_members_df.to_csv(OUT_BUBBLE_UNION_MEMBERS_CSV, index=False, encoding="utf-8-sig")
    if not args.skip_figures:
        save_total_figure(
            summary_df,
            random_df,
            coverage_df,
            bubble_coverage_df,
            OUT_FIG,
            bubble_sample_label=bubble_sample_label,
        )
        save_hmsfr_xy_membership_plot(
            sample_df,
            hmsfr_df,
            OUT_HMSFR_XY_FIG,
            OUT_HMSFR_XY_CSV,
            hmsfr_radius_kpc=hmsfr_radius_kpc,
            bubble_df=bubble_df,
            bubble_radii=bubble_radii,
        )

    surface_coverage_df = None
    surface_bubble_coverage_df = None
    surface_bubble_union_members_df = None
    if args.shell_mode == "volume":
        surface_coverage_df, surface_bubble_coverage_df, surface_bubble_union_members_df = compute_surface_robustness_outputs(
            sample_df,
            bubble_df,
            hmsfr_df,
            random_df,
            int(args.n_random),
            bubble_radii,
            hmsfr_radius_kpc=hmsfr_radius_kpc,
        )
        surface_coverage_df.to_csv(OUT_HMSFR_SURFACE_COVERAGE_CSV, index=False, encoding="utf-8-sig")
        surface_bubble_coverage_df.to_csv(OUT_BUBBLE_SURFACE_COVERAGE_CSV, index=False, encoding="utf-8-sig")
        surface_bubble_union_members_df.to_csv(
            OUT_BUBBLE_SURFACE_UNION_MEMBERS_CSV,
            index=False,
            encoding="utf-8-sig",
        )

    duplicate_removed = int(bubble_df.attrs.get("duplicate_rows_removed", 0))
    grade_c_removed = int(bubble_df.attrs.get("grade_c_removed", 0))
    print(f"bubble counting mode: {args.bubble_mode} ({mode_label})")
    print(
        f"small bubbles used (Grade {','.join(args.bubble_grades)}): {len(bubble_df)}; "
        f"Grade C removed: {grade_c_removed}; "
        f"duplicate rows removed: {duplicate_removed}"
    )
    print(f"HMSFR <= 3 kpc used: {len(hmsfr_df)}")
    print(f"saved HMSFR 3kpc catalog: {hmsfr_3kpc_output_path}")
    print(f"saved summary: {OUT_SUMMARY_CSV}")
    print(f"saved random detail: {OUT_RANDOM_CSV}")
    print(f"saved observed bubble members: {OUT_BUBBLE_MEMBERS_CSV}")
    print(f"saved observed HMSFR members: {OUT_HMSFR_MEMBERS_CSV}")
    print(f"saved HMSFR unique coverage: {OUT_HMSFR_COVERAGE_CSV}")
    print(f"saved bubble unique coverage: {OUT_BUBBLE_COVERAGE_CSV}")
    print(f"saved bubble unique membership: {OUT_BUBBLE_UNION_MEMBERS_CSV}")
    print_coverage_pair("finite-thickness shell coverage:", bubble_coverage_df, coverage_df)
    if surface_coverage_df is not None and surface_bubble_coverage_df is not None:
        print(f"saved HMSFR u=1 surface unique coverage: {OUT_HMSFR_SURFACE_COVERAGE_CSV}")
        print(f"saved bubble u=1 surface unique coverage: {OUT_BUBBLE_SURFACE_COVERAGE_CSV}")
        print(f"saved bubble u=1 surface unique membership: {OUT_BUBBLE_SURFACE_UNION_MEMBERS_CSV}")
        print_coverage_pair("zero-thickness u=1 surface coverage:", surface_bubble_coverage_df, surface_coverage_df)
    if args.skip_figures:
        print("skipped figure rendering")
    else:
        print(f"saved HMSFR XY membership: {OUT_HMSFR_XY_CSV}")
        print(f"saved total figure: {OUT_FIG}")
        print(f"saved XY membership figure: {OUT_HMSFR_XY_FIG}")

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
