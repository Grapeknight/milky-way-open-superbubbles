#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Estimate OSB dynamical ages and supernova requirements from fitted geometry.

The input is ../results/superbubble_final_fit_parameters.csv, produced by the
shell-fitting workflow. Each bubble has a fixed effective radius: sqrt(a*b)
for an intermediate horizontal section of its fitted shell. For open
ellipsoids, the default section lies halfway from the tabulated XY plane to
the retained cap apex; cylinders and other marks keep their XY-plane axes.
The geometry is not resampled, so the reported spread describes the adopted
physical parameter ranges rather than uncertainty in the fitted boundaries.

The baseline supernova requirement scales as n_e * R_eff**2 * v_exp**3 /
(E51 * alpha**3), and the age is alpha * R_eff / v_exp. The open-bubble model
increases this requirement by dividing by the retained disk-energy fraction
(1-theta)*(1-f_vent*f_open), where f_open=max(0, 1-H/(alpha*R_eff)). Here alpha
is the expansion-law index, theta is the energy-loss fraction, f_vent controls
additional leakage after opening, H is the adopted gas scale height, and E51
is the energy per supernova in units of 10**51 erg. Radii and H use pc,
velocities use km/s, n_e uses cm^-3, and ages are returned in Myr. Required
supernova numbers are continuous energy-equivalent estimates, not integers.

By default, 10,000 independently uniform parameter sets are generated with
seed 42 and reused across all bubbles. Their ranges are alpha=0.6-1.0,
theta=0.4-0.7, f_vent=0.2-0.5, v_exp=5-10 km/s, and n_e=1-3 cm^-3; H=150 pc
and E51=1 remain fixed. Separately evaluated nominal values (0.75, 0.61,
0.35, 7.5 km/s, 2 cm^-3) are reference cases, not Monte Carlo medians and
not automatically changed when command-line sampling ranges are changed.

Outputs under ../results/intermediate_output/19_superbubble_energy_budget/
are the per-bubble/per-draw sample CSV, a summary CSV with central values and
16/84 and 2.5/97.5 percentiles, and a JSON record of inputs and assumptions.
Script 20 reads the summary, by default using age_myr_median, n_open_median,
and n_base_median plus their 16th/84th percentiles. The figure
../results/figures/19_superbubble_energy_budget_age_sn.png shows medians and
16th-84th percentile intervals versus effective radius. Paths are resolved
from code/; main() writes products without changing the fitted input table."""

from __future__ import annotations
import time

import argparse
import csv
import json
import math
import os
import random
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Any


HERE = Path(__file__).resolve().parent
os.chdir(HERE)
OUT_DIR = Path("..") / "results" / "intermediate_output" / "19_superbubble_energy_budget"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"

DEFAULT_INPUT_TABLE = FINAL_DIR / "superbubble_final_fit_parameters.csv"

OUT_SAMPLES_CSV = OUT_DIR / "19_energy_budget_sample_details.csv"
OUT_SUMMARY_CSV = OUT_DIR / "19_superbubble_energy_budget_summary.csv"
OUT_CONFIG_JSON = OUT_DIR / "19_energy_budget_run_parameters.json"
OUT_FIG_PNG = FINAL_FIG_DIR / "19_superbubble_energy_budget_age_sn.png"

DEFAULT_N_MC = 10_000
DEFAULT_SEED = 42
DEFAULT_H_PC = 150.0
DEFAULT_E51 = 1.0
DEFAULT_COMPROMISE_FRAC = 0.5

# Sampling ranges describe assumed physical variation, not fitted geometry errors.
DEFAULT_ALPHA_RANGE = (0.6, 1.0)
DEFAULT_THETA_RANGE = (0.4, 0.7)
DEFAULT_F_VENT_RANGE = (0.2, 0.5)
DEFAULT_V_EXP_RANGE = (5.0, 10.0)
DEFAULT_NE_RANGE = (1.0, 3.0)

# Fixed illustrative reference values; deliberately distinct from sampled medians.
NOMINAL_ALPHA = 0.75
NOMINAL_THETA = 0.61
NOMINAL_F_VENT = 0.35
NOMINAL_V_EXP_KMS = 7.5
NOMINAL_NE_CM3 = 2.0

PC_PER_KMS_TO_MYR = 0.9777922216731285

REQUIRED_TABLE_COLUMNS = [
    "id",
    "mark",
    "shape",
    "xy_model",
    "center_z_kpc",
    "c_radius_kpc",
    "xy_plane_a_kpc",
    "xy_plane_b_kpc",
    "xy_plane_z_kpc",
]


@dataclass(frozen=True)
class ParameterRanges:
    """Independent uniform sampling intervals; velocity is in km/s and density in cm^-3."""
    alpha: tuple[float, float] = DEFAULT_ALPHA_RANGE
    theta: tuple[float, float] = DEFAULT_THETA_RANGE
    f_vent: tuple[float, float] = DEFAULT_F_VENT_RANGE
    v_exp: tuple[float, float] = DEFAULT_V_EXP_RANGE
    n_e: tuple[float, float] = DEFAULT_NE_RANGE


@dataclass(frozen=True)
class BubbleGeometry:
    """Fixed section geometry in pc, retained for tracing each energy estimate to its fit."""
    sb_id: int
    mark: int
    shape: str
    xy_model: str
    a_compromise_pc: float
    b_compromise_pc: float
    compromise_scale: float
    reff_pc: float


@dataclass(frozen=True)
class ParameterDraw:
    """One shared physical-parameter realization, identified consistently across bubbles."""
    sample_id: int
    alpha: float
    theta: float
    f_vent: float
    v_exp_kms: float
    n_e_cm3: float


@dataclass(frozen=True)
class Sample:
    """One bubble/draw result: pc geometry, Myr age, and dimensionless counts/fractions."""
    sb_id: int
    sample_id: int
    mark: int
    shape: str
    xy_model: str
    a_compromise_pc: float
    b_compromise_pc: float
    compromise_scale: float
    reff_pc: float
    alpha: float
    theta: float
    f_vent: float
    v_exp_kms: float
    n_e_cm3: float
    age_myr: float
    n_base: float
    f_open: float
    epsilon_disk: float
    n_open: float


def baseline_n_sn(
    reff_pc: float,
    alpha: float,
    v_exp_kms: float,
    n_e_cm3: float,
    e51: float = DEFAULT_E51,
) -> float:
    """Return the baseline energy-equivalent supernova count.

    The numerical normalization expects R_eff in pc, v_exp in km/s, n_e in
    cm^-3, and E51 in 10**51 erg per event. Alpha is dimensionless. Keeping
    this expression separate makes the later loss/leakage correction explicit."""
    return (
        (n_e_cm3 / e51)
        * (reff_pc / 97.0) ** 2
        * (v_exp_kms / 5.7) ** 3
        * (0.6 / alpha) ** 3
    )


def dynamical_age_myr(reff_pc: float, alpha: float, v_exp_kms: float) -> float:
    """Convert t=alpha*R_eff/v_exp from pc/(km/s) to Myr.

    This is the instantaneous age implied by R proportional to t**alpha;
    it uses the adopted expansion speed rather than a measured age."""
    return alpha * reff_pc / v_exp_kms * PC_PER_KMS_TO_MYR


def open_fraction(reff_pc: float, alpha: float, h_pc: float) -> float:
    """Return max(0, 1-H/(alpha*R_eff)), the adopted opening factor.

    The model has no opening correction while alpha*R_eff <= H. H and
    R_eff must share pc units; input validation ensures a positive denominator."""
    return max(0.0, 1.0 - h_pc / (alpha * reff_pc))


def open_model_n_sn(
    reff_pc: float,
    alpha: float,
    theta: float,
    f_vent: float,
    v_exp_kms: float,
    n_e_cm3: float,
    h_pc: float,
    e51: float,
) -> tuple[float, float, float]:
    """Return (required SNe, opening factor, retained disk-energy fraction).

    Theta removes a fraction of the baseline energy budget. The remaining
    energy is reduced further by f_vent*f_open, so dividing N_base by
    epsilon_disk increases the requirement when either loss term is active."""
    n_base = baseline_n_sn(
        reff_pc=reff_pc,
        alpha=alpha,
        v_exp_kms=v_exp_kms,
        n_e_cm3=n_e_cm3,
        e51=e51,
    )
    f_open = open_fraction(reff_pc=reff_pc, alpha=alpha, h_pc=h_pc)
    # Both loss factors are dimensionless; validated ranges keep their product positive.
    epsilon_disk = (1.0 - theta) * (1.0 - f_vent * f_open)
    n_open = n_base / epsilon_disk
    return n_open, f_open, epsilon_disk


def percentile(values: list[float], pct: float) -> float:
    """Linearly interpolate sorted samples at rank (N-1)*pct/100; pct is a percentage."""
    if not values:
        raise ValueError('Invalid input or missing required data.')
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]

    rank = (len(ordered) - 1) * pct / 100.0
    lo = int(rank)
    hi = min(lo + 1, len(ordered) - 1)
    weight = rank - lo
    return ordered[lo] * (1.0 - weight) + ordered[hi] * weight


def summarize(values: list[float]) -> dict[str, float]:
    """Summarize the sampled distribution without assuming Gaussian symmetry.

    p16/p84 and p025/p975 bound central 68% and 95% intervals, respectively.
    These intervals reflect the specified input ranges, not a fitted posterior."""
    return {
        "mean": mean(values),
        "median": percentile(values, 50.0),
        "p16": percentile(values, 16.0),
        "p84": percentile(values, 84.0),
        "p025": percentile(values, 2.5),
        "p975": percentile(values, 97.5),
        "min": min(values),
        "max": max(values),
    }


def add_stats(prefix: str, values: list[float], row: dict[str, Any]) -> None:
    """Append prefixed summary fields while retaining the units of the input values."""
    stats = summarize(values)
    for key, value in stats.items():
        row[f"{prefix}_{key}"] = value


def validate_range(name: str, bounds: tuple[float, float], min_allowed: float) -> None:
    low, high = bounds
    if low > high:
        raise ValueError(f"{name} range must satisfy low <= high.")
    if low < min_allowed:
        raise ValueError(f"{name} range lower bound must >= {min_allowed}.")


def validate_inputs(
    n_mc: int,
    h_pc: float,
    e51: float,
    ranges: ParameterRanges,
    compromise_frac: float,
) -> None:
    """Require physical ranges that keep ages and energy-retention denominators valid."""
    if n_mc <= 0:
        raise ValueError("The Monte Carlo sample count must be a positive integer.")
    if h_pc < 0:
        raise ValueError("The gas scale height H must be non-negative.")
    if e51 <= 0:
        raise ValueError("E51 must be positive.")
    if not 0.0 <= compromise_frac <= 1.0:
        raise ValueError("The compromise section fraction must lie within [0, 1].")

    validate_range("alpha", ranges.alpha, 0.0)
    validate_range("theta", ranges.theta, 0.0)
    validate_range("f_vent", ranges.f_vent, 0.0)
    validate_range("v_exp", ranges.v_exp, 0.0)
    validate_range("n_e", ranges.n_e, 0.0)

    if ranges.alpha[0] <= 0.0:
        raise ValueError("alpha range lower bound must > 0.")
    if ranges.v_exp[0] <= 0.0:
        raise ValueError("v_exp range lower bound must > 0.")
    if ranges.n_e[0] <= 0.0:
        raise ValueError("n_e range lower bound must > 0.")
    if ranges.theta[1] >= 1.0:
        raise ValueError("theta range upper bound must < 1.")
    if ranges.f_vent[1] >= 1.0:
        raise ValueError("f_vent range upper bound must < 1.")


def parse_float(row: dict[str, str], column: str) -> float:
    value = row.get(column, "")
    if value is None or str(value).strip() == "":
        raise ValueError('Invalid input or missing required data.')
    return float(value)


def compromise_scale(row: dict[str, str], frac: float = DEFAULT_COMPROMISE_FRAC) -> float:
    """Rescale XY-plane axes to a section between that plane and the cap apex.

    For an ellipsoid, horizontal axes vary as k(z)=sqrt(1-((z-cz)/c)**2).
    The ratio k(z_mid)/k(z_base) moves the tabulated axes to the selected
    section. Mark 1 retains the lower cap (apex cz-c), and mark 2 the upper
    cap (cz+c). The default frac=0.5 is a height midpoint, not an area
    average. Non-ellipsoids, other marks, and degenerate base scales use 1."""
    if row["shape"].strip().lower() != "ellipsoid":
        return 1.0

    mark = int(row["mark"])
    if mark not in (1, 2):
        return 1.0

    cz = parse_float(row, "center_z_kpc")
    c_radius = parse_float(row, "c_radius_kpc")
    z_base = parse_float(row, "xy_plane_z_kpc")
    if c_radius <= 0:
        return 1.0

    z_apex = cz - c_radius if mark == 1 else cz + c_radius
    z_mid = z_base + frac * (z_apex - z_base)
    k_base_sq = max(0.0, 1.0 - ((z_base - cz) / c_radius) ** 2)
    k_mid_sq = max(0.0, 1.0 - ((z_mid - cz) / c_radius) ** 2)
    k_base = math.sqrt(k_base_sq)
    k_mid = math.sqrt(k_mid_sq)
    return k_mid / k_base if k_base > 0 else 1.0


def read_bubble_geometries(table_path: Path, compromise_frac: float) -> list[BubbleGeometry]:
    """Read required fit columns and derive one area-equivalent radius per bubble.

    Convert the rescaled semiaxes from kpc to pc, then use sqrt(a*b), the
    radius of a circle with the same section area. Sort by ID so output
    ordering is stable across input-table row orders."""
    if not table_path.exists():
        raise FileNotFoundError('Missing input file.')

    bubbles: list[BubbleGeometry] = []
    with table_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError('Invalid input or missing required data.')
        missing = [column for column in REQUIRED_TABLE_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise ValueError(f"final parameter table is missing required columns: {', '.join(missing)}")

        for row in reader:
            scale = compromise_scale(row, frac=compromise_frac)
            a_compromise_pc = parse_float(row, "xy_plane_a_kpc") * scale * 1000.0
            b_compromise_pc = parse_float(row, "xy_plane_b_kpc") * scale * 1000.0
            if a_compromise_pc <= 0 or b_compromise_pc <= 0:
                raise ValueError('Invalid input or missing required data.')
            reff_pc = math.sqrt(a_compromise_pc * b_compromise_pc)
            bubbles.append(
                BubbleGeometry(
                    sb_id=int(row["id"]),
                    mark=int(row["mark"]),
                    shape=row["shape"].strip(),
                    xy_model=row["xy_model"].strip(),
                    a_compromise_pc=a_compromise_pc,
                    b_compromise_pc=b_compromise_pc,
                    compromise_scale=scale,
                    reff_pc=reff_pc,
                )
            )

    if not bubbles:
        raise ValueError('Invalid input or missing required data.')
    return sorted(bubbles, key=lambda bubble: bubble.sb_id)


def sample_uniform(rng: random.Random, bounds: tuple[float, float]) -> float:
    low, high = bounds
    return rng.uniform(low, high)


def draw_parameter_sets(n_mc: int, seed: int, ranges: ParameterRanges) -> list[ParameterDraw]:
    """Generate reproducible independent uniform draws from the five ranges.

    This list is built once, then reused for every bubble: the same sample_id
    means the same physical assumptions across different fixed geometries."""
    rng = random.Random(seed)
    return [
        ParameterDraw(
            sample_id=sample_id,
            alpha=sample_uniform(rng, ranges.alpha),
            theta=sample_uniform(rng, ranges.theta),
            f_vent=sample_uniform(rng, ranges.f_vent),
            v_exp_kms=sample_uniform(rng, ranges.v_exp),
            n_e_cm3=sample_uniform(rng, ranges.n_e),
        )
        for sample_id in range(1, n_mc + 1)
    ]


def run_monte_carlo(
    bubbles: list[BubbleGeometry],
    draws: list[ParameterDraw],
    h_pc: float,
    e51: float,
) -> list[Sample]:
    """Evaluate baseline/open counts and age for every fixed geometry and shared draw."""
    samples: list[Sample] = []

    for bubble in bubbles:
        for draw in draws:
            age_myr = dynamical_age_myr(bubble.reff_pc, draw.alpha, draw.v_exp_kms)
            n_base = baseline_n_sn(bubble.reff_pc, draw.alpha, draw.v_exp_kms, draw.n_e_cm3, e51)
            n_open, f_open, epsilon_disk = open_model_n_sn(
                reff_pc=bubble.reff_pc,
                alpha=draw.alpha,
                theta=draw.theta,
                f_vent=draw.f_vent,
                v_exp_kms=draw.v_exp_kms,
                n_e_cm3=draw.n_e_cm3,
                h_pc=h_pc,
                e51=e51,
            )
            samples.append(
                Sample(
                    sb_id=bubble.sb_id,
                    sample_id=draw.sample_id,
                    mark=bubble.mark,
                    shape=bubble.shape,
                    xy_model=bubble.xy_model,
                    a_compromise_pc=bubble.a_compromise_pc,
                    b_compromise_pc=bubble.b_compromise_pc,
                    compromise_scale=bubble.compromise_scale,
                    reff_pc=bubble.reff_pc,
                    alpha=draw.alpha,
                    theta=draw.theta,
                    f_vent=draw.f_vent,
                    v_exp_kms=draw.v_exp_kms,
                    n_e_cm3=draw.n_e_cm3,
                    age_myr=age_myr,
                    n_base=n_base,
                    f_open=f_open,
                    epsilon_disk=epsilon_disk,
                    n_open=n_open,
                )
            )

    return samples


def summarize_by_bubble(
    bubbles: list[BubbleGeometry],
    samples: list[Sample],
    e51: float,
    h_pc: float,
) -> list[dict[str, Any]]:
    """Aggregate draws by ID and evaluate a separate fixed nominal reference.

    The nominal constants are independent of the sampling intervals; nonlinear
    formulas mean their outputs need not match the sampled medians."""
    samples_by_id: dict[int, list[Sample]] = {bubble.sb_id: [] for bubble in bubbles}
    for sample in samples:
        samples_by_id[sample.sb_id].append(sample)

    rows: list[dict[str, Any]] = []
    for bubble in bubbles:
        bubble_samples = samples_by_id[bubble.sb_id]
        row: dict[str, Any] = {
            "sb_id": bubble.sb_id,
            "mark": bubble.mark,
            "shape": bubble.shape,
            "xy_model": bubble.xy_model,
            "a_compromise_pc": bubble.a_compromise_pc,
            "b_compromise_pc": bubble.b_compromise_pc,
            "compromise_scale": bubble.compromise_scale,
            "reff_pc": bubble.reff_pc,
        }
        add_stats("n_base", [sample.n_base for sample in bubble_samples], row)
        add_stats("n_open", [sample.n_open for sample in bubble_samples], row)
        add_stats("age_myr", [sample.age_myr for sample in bubble_samples], row)
        add_stats("f_open", [sample.f_open for sample in bubble_samples], row)
        add_stats("epsilon_disk", [sample.epsilon_disk for sample in bubble_samples], row)

        nominal_base = baseline_n_sn(
            bubble.reff_pc,
            NOMINAL_ALPHA,
            NOMINAL_V_EXP_KMS,
            NOMINAL_NE_CM3,
            e51,
        )
        nominal_open, nominal_f_open, nominal_epsilon = open_model_n_sn(
            reff_pc=bubble.reff_pc,
            alpha=NOMINAL_ALPHA,
            theta=NOMINAL_THETA,
            f_vent=NOMINAL_F_VENT,
            v_exp_kms=NOMINAL_V_EXP_KMS,
            n_e_cm3=NOMINAL_NE_CM3,
            h_pc=h_pc,
            e51=e51,
        )
        row.update(
            {
                "nominal_n_base": nominal_base,
                "nominal_n_open": nominal_open,
                "nominal_age_myr": dynamical_age_myr(
                    bubble.reff_pc,
                    NOMINAL_ALPHA,
                    NOMINAL_V_EXP_KMS,
                ),
                "nominal_f_open": nominal_f_open,
                "nominal_epsilon_disk": nominal_epsilon,
            }
        )
        rows.append(row)

    return rows


def write_samples_csv(samples: list[Sample], path: Path) -> None:
    """Write the full geometry, parameter draw, and model results for each sample."""
    fieldnames = list(Sample.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for sample in samples:
            writer.writerow(sample.__dict__)


def write_summary_csv(rows: list[dict[str, Any]], path: Path) -> None:
    """Write one summary row per bubble, including the columns consumed by Script 20."""
    if not rows:
        raise ValueError('Invalid input or missing required data.')
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_config_json(
    args: argparse.Namespace,
    ranges: ParameterRanges,
    table_path: Path,
    bubbles: list[BubbleGeometry],
    path: Path,
) -> None:
    """Record run settings, fixed nominal assumptions, formulas, and output locations."""
    try:
        table_label = Path(os.path.relpath(table_path.resolve(), HERE)).as_posix()
    except (OSError, ValueError):
        table_label = str(table_path)
    payload = {
        "input_table": table_label,
        "bubble_count": len(bubbles),
        "n_mc_per_bubble": args.n_mc,
        "seed": args.seed,
        "parameter_draws_reused_across_bubbles": True,
        "h_pc": args.h_pc,
        "e51": args.e51,
        "compromise_frac": args.compromise_frac,
        "equivalent_radius_definition": "sqrt(a_compromise_pc*b_compromise_pc)",
        "alpha_range": list(ranges.alpha),
        "theta_range": list(ranges.theta),
        "f_vent_range": list(ranges.f_vent),
        "v_exp_range_kms": list(ranges.v_exp),
        "n_e_range_cm3": list(ranges.n_e),
        "nominal": {
            "alpha": NOMINAL_ALPHA,
            "theta": NOMINAL_THETA,
            "f_vent": NOMINAL_F_VENT,
            "v_exp_kms": NOMINAL_V_EXP_KMS,
            "n_e_cm3": NOMINAL_NE_CM3,
        },
        "formula": {
            "a_compromise": "xy_plane_a_kpc*compromise_scale*1000 pc",
            "b_compromise": "xy_plane_b_kpc*compromise_scale*1000 pc",
            "reff_pc": "sqrt(a_compromise_pc*b_compromise_pc)",
            "n_base": "(n_e/E51)*(R_eff/97 pc)^2*(v_exp/5.7 km/s)^3*(0.6/alpha)^3",
            "age_myr": "alpha*R_eff/v_exp*0.9777922216731285",
            "n_open": "N_base/((1-theta)*(1-f_vent*f_open))",
            "f_open": "max(0, 1 - H/(alpha*R_eff))",
        },
        "outputs": {
            "samples_csv": str(OUT_SAMPLES_CSV),
            "summary_csv": str(OUT_SUMMARY_CSV),
            "figure_png": str(OUT_FIG_PNG),
        },
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def plot_budget_summary(rows: list[dict[str, Any]], path_png: Path) -> None:
    """Plot Monte Carlo medians and asymmetric 16th-84th percentile bars.

    Radius and the supernova-count axis are logarithmic; the independent
    right axis reports dynamical age in Myr. Sorting by radius only controls
    the plotting order and connecting lines."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    rows_sorted = sorted(rows, key=lambda row: row["reff_pc"])
    labels = [f"SB{int(row['sb_id'])}" for row in rows_sorted]
    reff = [row["reff_pc"] for row in rows_sorted]

    def plain_tick(value: float, _: int | None = None) -> str:
        return f"{value:.0f}" if value >= 1.0 else f"{value:g}"

    def stagger_label_y(values: list[float]) -> list[float]:
        label_levels = [0.98, 0.88, 0.78, 0.68, 0.58]
        y_positions: list[float] = []
        previous_log: float | None = None
        close_run_index = 0
        for value in values:
            current_log = math.log10(value)
            if previous_log is not None and current_log - previous_log <= 0.035:
                close_run_index += 1
            else:
                close_run_index = 0
            y_positions.append(label_levels[close_run_index % len(label_levels)])
            previous_log = current_log
        return y_positions

    n_base = [row["n_base_median"] for row in rows_sorted]
    n_base_low = [row["n_base_median"] - row["n_base_p16"] for row in rows_sorted]
    n_base_high = [row["n_base_p84"] - row["n_base_median"] for row in rows_sorted]
    n_open = [row["n_open_median"] for row in rows_sorted]
    n_open_low = [row["n_open_median"] - row["n_open_p16"] for row in rows_sorted]
    n_open_high = [row["n_open_p84"] - row["n_open_median"] for row in rows_sorted]
    age = [row["age_myr_median"] for row in rows_sorted]
    age_low = [row["age_myr_median"] - row["age_myr_p16"] for row in rows_sorted]
    age_high = [row["age_myr_p84"] - row["age_myr_median"] for row in rows_sorted]
    tokens = {
        "surface": "#FCFCFD",
        "panel": "#FFFFFF",
        "ink": "#1F2430",
        "muted": "#6F768A",
        "grid": "#E6E8F0",
        "axis": "#D7DBE7",
        "blue": "#1B6CA8",
        "blue_light": "#9CC3E6",
        "orange": "#D17A22",
        "orange_light": "#EDBD8F",
        "purple": "#8F3F71",
        "purple_light": "#C9A1BA",
    }

    plt.rcParams.update(
        {
            "font.family": ["Segoe UI", "DejaVu Sans", "Arial", "sans-serif"],
            "font.size": 14,
            "axes.labelcolor": tokens["ink"],
            "xtick.color": tokens["muted"],
            "ytick.color": tokens["muted"],
            "axes.edgecolor": tokens["axis"],
            "axes.linewidth": 1.0,
        }
    )

    fig, ax_sn = plt.subplots(figsize=(15.0, 7.2), constrained_layout=False, facecolor=tokens["surface"])
    ax_age = ax_sn.twinx()
    fig.subplots_adjust(left=0.085, right=0.905, top=0.84, bottom=0.17)

    ax_sn.set_facecolor(tokens["panel"])
    ax_sn.grid(axis="y", color=tokens["grid"], linewidth=0.9)
    ax_sn.spines["top"].set_visible(False)
    ax_age.spines["top"].set_visible(False)
    ax_sn.spines["right"].set_visible(False)
    ax_age.spines["left"].set_visible(False)

    for x_value in reff:
        ax_sn.axvline(x_value, color=tokens["grid"], linewidth=0.8, zorder=0)

    base_line = ax_sn.errorbar(
        reff,
        n_base,
        yerr=[n_base_low, n_base_high],
        fmt="^-",
        markersize=6.2,
        linewidth=2.1,
        color=tokens["orange"],
        ecolor=tokens["orange_light"],
        elinewidth=1.6,
        capsize=3.0,
        markeredgecolor=tokens["ink"],
        markeredgewidth=0.35,
        zorder=3,
        label="Baseline required SNe",
    )
    open_line = ax_sn.errorbar(
        reff,
        n_open,
        yerr=[n_open_low, n_open_high],
        fmt="o-",
        markersize=6.4,
        linewidth=2.2,
        color=tokens["purple"],
        ecolor=tokens["purple_light"],
        elinewidth=1.8,
        capsize=3.2,
        markeredgecolor=tokens["ink"],
        markeredgewidth=0.35,
        zorder=4,
        label="Open-superbubble required SNe",
    )
    age_line = ax_age.errorbar(
        reff,
        age,
        yerr=[age_low, age_high],
        fmt="s-",
        markersize=6.0,
        linewidth=2.2,
        color=tokens["blue"],
        ecolor=tokens["blue_light"],
        elinewidth=1.8,
        capsize=3.2,
        markeredgecolor=tokens["ink"],
        markeredgewidth=0.35,
        zorder=5,
        label="Dynamical age",
    )

    ax_sn.set_xscale("log")
    ax_sn.set_yscale("log")
    ax_sn.set_ylim(3.0, 2500.0)
    for x_value, label, label_y in zip(reff, labels, stagger_label_y(reff), strict=True):
        ax_sn.text(
            x_value,
            label_y,
            label,
            ha="center",
            va="top",
            rotation=90,
            fontsize=14,
            color=tokens["ink"],
            zorder=10,
            clip_on=True,
            transform=ax_sn.get_xaxis_transform(),
            bbox={"facecolor": tokens["panel"], "edgecolor": "none", "alpha": 0.72, "pad": 0.5},
        )
    ax_sn.yaxis.set_major_locator(mticker.FixedLocator([5, 10, 20, 50, 100, 200, 500, 1000, 2000]))
    ax_sn.yaxis.set_major_formatter(mticker.FuncFormatter(plain_tick))
    ax_sn.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax_age.yaxis.set_major_locator(mticker.MultipleLocator(20))
    ax_age.yaxis.set_major_formatter(mticker.FuncFormatter(plain_tick))
    ax_age.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax_sn.set_xlabel("Equivalent radius, $R_{\\rm eff}$ (pc)", fontsize=16)
    ax_sn.set_ylabel("Required supernovae", color=tokens["ink"], fontsize=16)
    ax_age.set_ylabel("Dynamical age (Myr)", color=tokens["blue"], fontsize=16)
    ax_sn.tick_params(axis="both", labelsize=14)
    ax_age.tick_params(axis="y", colors=tokens["blue"], labelsize=14)
    ax_sn.tick_params(axis="y", colors=tokens["ink"], labelsize=14)
    ax_sn.set_xlim(min(reff) * 0.92, max(reff) * 1.06)
    ax_sn.xaxis.set_major_locator(mticker.FixedLocator([170, 200, 300, 400, 500, 700, 1000, 1200]))
    ax_sn.xaxis.set_major_formatter(mticker.FuncFormatter(plain_tick))
    ax_sn.xaxis.set_minor_formatter(mticker.NullFormatter())

    handles = [base_line.lines[0], open_line.lines[0], age_line.lines[0]]
    labels_for_legend = [
        "Baseline required SNe",
        "Open-superbubble required SNe",
        "Dynamical age",
    ]
    ax_sn.legend(
        handles,
        labels_for_legend,
        loc="upper left",
        frameon=False,
        ncol=3,
        fontsize=16,
        bbox_to_anchor=(0.0, 1.12),
    )

    fig.savefig(path_png, dpi=300)
    plt.close(fig)


def format_summary_line(label: str, stats: dict[str, float], unit: str = "") -> str:
    unit_text = f" {unit}" if unit else ""
    return (
        f"{label}: median={stats['median']:.2f}{unit_text}, "
        f"68%=[{stats['p16']:.2f}, {stats['p84']:.2f}]{unit_text}, "
        f"95%=[{stats['p025']:.2f}, {stats['p975']:.2f}]{unit_text}"
    )


def parse_args() -> argparse.Namespace:
    """Expose input geometry, sample count/seed, H in pc, E51, and sampling ranges.

    The compromise fraction selects the section height from base toward apex;
    v_exp bounds use km/s and n_e bounds use cm^-3."""
    parser = argparse.ArgumentParser(description='Run Script 19: OSBs energy budget.')
    parser.add_argument(
        "--input-table",
        type=Path,
        default=DEFAULT_INPUT_TABLE,
        help='Input file used by this stage.',
    )
    parser.add_argument("--n-mc", type=int, default=DEFAULT_N_MC, help='Number of random or Monte Carlo realizations.')
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help='Random seed for reproducible controls.')
    parser.add_argument("--h-pc", type=float, default=DEFAULT_H_PC, help='Command-line option for the documented workflow.')
    parser.add_argument("--e51", type=float, default=DEFAULT_E51, help='Command-line option for the documented workflow.')
    parser.add_argument(
        "--compromise-frac",
        type=float,
        default=DEFAULT_COMPROMISE_FRAC,
        help='Command-line option for the documented workflow.',
    )
    parser.add_argument("--alpha-range", nargs=2, type=float, default=DEFAULT_ALPHA_RANGE, metavar=("LOW", "HIGH"))
    parser.add_argument("--theta-range", nargs=2, type=float, default=DEFAULT_THETA_RANGE, metavar=("LOW", "HIGH"))
    parser.add_argument("--f-vent-range", nargs=2, type=float, default=DEFAULT_F_VENT_RANGE, metavar=("LOW", "HIGH"))
    parser.add_argument("--v-exp-range", nargs=2, type=float, default=DEFAULT_V_EXP_RANGE, metavar=("LOW", "HIGH"))
    parser.add_argument("--ne-range", nargs=2, type=float, default=DEFAULT_NE_RANGE, metavar=("LOW", "HIGH"))
    return parser.parse_args()


def main() -> None:
    """Validate assumptions, derive fixed radii, sample the model, and write Script 20 inputs."""
    args = parse_args()
    ranges = ParameterRanges(
        alpha=tuple(args.alpha_range),
        theta=tuple(args.theta_range),
        f_vent=tuple(args.f_vent_range),
        v_exp=tuple(args.v_exp_range),
        n_e=tuple(args.ne_range),
    )
    validate_inputs(args.n_mc, args.h_pc, args.e51, ranges, args.compromise_frac)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)

    table_path = args.input_table.resolve()
    bubbles = read_bubble_geometries(table_path, args.compromise_frac)
    draws = draw_parameter_sets(args.n_mc, args.seed, ranges)
    samples = run_monte_carlo(bubbles=bubbles, draws=draws, h_pc=args.h_pc, e51=args.e51)
    summary_rows = summarize_by_bubble(bubbles=bubbles, samples=samples, e51=args.e51, h_pc=args.h_pc)

    write_samples_csv(samples, OUT_SAMPLES_CSV)
    write_summary_csv(summary_rows, OUT_SUMMARY_CSV)
    write_config_json(args, ranges, table_path, bubbles, OUT_CONFIG_JSON)
    plot_budget_summary(summary_rows, OUT_FIG_PNG)

    reff_values = [bubble.reff_pc for bubble in bubbles]
    n_open_values = [row["n_open_median"] for row in summary_rows]
    age_values = [row["age_myr_median"] for row in summary_rows]

    print(f"Computed energy budget for {len(summary_rows)} OSBs.")
    print(
        f"R_eff range={min(reff_values):.1f}-{max(reff_values):.1f} pc; "
        f"age range={min(age_values):.2f}-{max(age_values):.2f} Myr; "
        f"open-SN median range={min(n_open_values):.2f}-{max(n_open_values):.2f}."
    )
    print(f"Assumed ranges: v_exp={ranges.v_exp} km/s, n_e={ranges.n_e} cm^-3, H={args.h_pc:g} pc")
    print(f"Input geometry table: {Path(os.path.relpath(table_path, HERE)).as_posix()}")
    print(f"Summary table: {OUT_SUMMARY_CSV}")
    print(f"Detail table: {OUT_SAMPLES_CSV}")
    print(f"Run-parameter JSON: {OUT_CONFIG_JSON}")
    print(f"Figure: {OUT_FIG_PNG}")

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
