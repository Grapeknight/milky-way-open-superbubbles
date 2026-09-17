# -*- coding: utf-8 -*-
"""Test dust structure around fitted open superbubbles using 3D IAAFT fields.

Scientific question and null construction
---------------------------------------
At each previously fitted geometry, compare the observed shell contrast, cavity
density and radial ridge alignment with randomized local dust fields. The fitted
centre, axes, orientation and cap/cylinder selection remain fixed throughout;
this calculation does not repeat the upstream search or geometry fitting.

The iterative amplitude-adjusted Fourier transform (IAAFT) target is the entire
local 3D dust box, including the cavity, shell and surroundings. After filling
missing voxels by nearest neighbours and optional block averaging, its sorted
density values and abs(fftn(cube)) define the two constraints. All wavevectors,
including the DC component and resolved large-scale modes, enter the target.
There is no cavity excision, background subtraction, taper, detrending or
high-pass filtering. The target spectrum therefore includes the observed
cavity's contribution together with that of the surrounding dust.

Each realization starts from a random permutation and alternates amplitude
replacement with rank remapping. The returned rank-remapped field preserves
the processed input histogram exactly; its Fourier amplitudes are approximate.
The command-line default allows eight iterations, with earlier stopping when
the change in amplitude error is below 1e-6. Eight iterations do not guarantee
spectral convergence; per-realization errors and first-realization spectra are
saved so the actual agreement can be inspected.

Statistics and finite-sample interpretation
------------------------------------------
The default N = 200 realizations are evaluated on the same processed voxel grid
as the observation. Tail events require a shell-ridge/inner-density ratio at
least as large, an inner mean at most as large, or a mean absolute normalized
ridge offset at most as large as observed. The joint event requires all three
conditions in the same realization. Counts use (1 + count)/(N + 1), giving a
minimum value of 1/201 for the default run. Missing surrogate statistics never
count as tail events but remain in N. An unavailable observed ridge offset gives
NaN offset and joint values; missing observed ratio/inner mean instead prevents
their comparisons from counting. These are conditional fixed-geometry tests,
with the implemented missing-statistic convention, not a repeated detection
experiment or a correction for searching over candidate geometries.
In particular, p_dust_joint3 is a joint tail fraction, without calibration as an
overall test p-value: a threshold of 0.05 does not establish a 5% false-positive
rate. This conditional random-field comparison does not simulate gas dynamics,
dust-map reconstruction uncertainty, or the full discovery/selection procedure,
and does not by itself identify the physical mechanism that formed a cavity.

Inputs, products and units
--------------------------
Paths are resolved from code/ within the data/, code/, results/ package layout.
Default inputs are ../results/superbubble_final_fit_parameters.csv and
../results/intermediate_output/3d_dust_map_products/raw_3d_dust_cube.parquet.
Coordinates and fitted lengths are in kpc; dust density is in mag kpc^-1.
Normal runs write per-target surrogate samples, full summaries and compact
p-value tables under ../results/intermediate_output/8_random_density_fluctuation_test/,
plus ../results/figures/8_condition_test_joint_p_all30.png for the full sample.
Subset runs use target-specific summary filenames.

For each target, the first realization (mc_index = 0; seed + target ID) also
produces a PNG comparison in ../results/figures/random_density_diagnostics/.
Its matched XY/XZ/YZ slices share colour limits and fitted boundaries; only
these displayed slices receive Gaussian smoothing. Unsmoothed processed
cubes (NPZ), radially binned 3D spectra (CSV), and provenance/diagnostics (JSON)
are saved in the intermediate output's diagnostics/ subdirectory. The density
PDF, spectra and test statistics all use the fields before display smoothing.
The realization is selected by order, without scoring or visual selection.
--diagnostics-only generates that one realization per target and omits the
Monte Carlo sample tables, p-value summaries and summary figure.

Code map
--------
read_local_xyz_cube and build_cube_context extract and preprocess the local
field; generate_iaaft_surrogate_3d constructs the random fields. cap_geometry,
build_local and build_geometry_context translate the final fit into fixed
voxel masks and directional/radial bins. compute_dust_metrics measures cavity
and ridge properties, and compute_pvalues counts their empirical tail events.
power_spectrum_diagnostic checks the unsmoothed 3D fields;
save_realization_comparison and save_summary_figure render the diagnostic and
summary figures. main manages target selection, reproducible random streams,
realization loops and output tables.

Run from the package root with:
    python code/08_OSBs_random_dust_test.py --targets all --diagnostics-only
Omit --diagnostics-only to generate the requested Monte Carlo sample and its
tables. The script changes its working directory to code/, so relative paths
passed through --geom-csv, --xy-data-path, --diagnostic-dir and
--comparison-fig-dir are interpreted there. Use --diagnostic-dir and
--comparison-fig-dir to retain separate diagnostic products when comparing
settings. If the iteration budget changes for the
statistical test, regenerate all requested samples and their probability tables
with that budget; a single diagnostic draw cannot establish ensemble fidelity.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.ndimage import distance_transform_edt, gaussian_filter


# Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "8_random_density_fluctuation_test"
FINAL_FIG_DIR = Path("..") / "results" / "figures"
COMPARISON_FIG_DIR = FINAL_FIG_DIR / "random_density_diagnostics"

DEFAULT_GEOM_CSV = Path("..") / "results" / "superbubble_final_fit_parameters.csv"
DEFAULT_XY_DATA_PATH = Path("..") / "results" / "intermediate_output" / "3d_dust_map_products" / "raw_3d_dust_cube.parquet"
BAR_COLOR = "#4C78A8"
COMPARISON_CMAP = "Spectral_r"
COMPARISON_SMOOTH_SIGMA_VOXELS = 1.0  # Gaussian sigma on the actual (possibly rebinned) test grid.


def metadata_path(path: Path) -> str:
    """Record portable paths relative to code/, including external inputs when possible."""
    try:
        return Path(os.path.relpath(Path(path).resolve(), HERE)).as_posix()
    except ValueError:
        # An external input on another Windows drive has no relative path.
        return Path(path).resolve().as_posix()

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Microsoft YaHei", "SimHei"],
        "axes.unicode_minus": False,
        "font.size": 11.0,
        "axes.labelsize": 11.5,
        "axes.titlesize": 12.0,
        "xtick.labelsize": 10.0,
        "ytick.labelsize": 10.5,
        "legend.fontsize": 10.5,
        "axes.linewidth": 1.0,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "savefig.dpi": 300,
    }
)

# Dimensionless elliptical/ellipsoidal radius u defines the cavity and shell;
# the cylindrical slab has a fixed half-height of 0.025 kpc (25 pc).
SHELL_U_MIN = 0.9
SHELL_U_MAX = 1.2
INNER_U_MAX = 0.5
CYLINDER_HALF_HEIGHT_KPC = 0.025

INNER_STD_CLIP_ITERS = 3
INNER_STD_CLIP_SIGMA = 3.0

RIDGE_N_AZIMUTH = 36
RIDGE_N_COSTHETA_FULL = 12
RIDGE_N_COSTHETA_HALF = 6
RIDGE_N_RADIAL = 16
RIDGE_MIN_SECTOR_VOXELS = 10
RIDGE_MIN_RADIAL_VOXELS = 2
RIDGE_PEAK_INNER_RATIO_MIN = 2.0
RIDGE_PEAK_MIN_VOXELS = 30
RIDGE_GMM_MAX_ITER = 200
RIDGE_GMM_TOL = 1e-6


def ensure_existing_file(path, label):
    """Resolve an input as a Path and fail early if it is absent."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{label}  does not exist: {p}")
    return p


def parse_targets(text, geom_df):
    """Use table order for all targets, or parse comma/semicolon-separated IDs."""
    if text is None or str(text).strip().lower() in {"", "all"}:
        return [int(v) for v in geom_df["id"].tolist()]
    values = []
    for item in str(text).replace(";", ",").split(","):
        item = item.strip()
        if item:
            values.append(int(item))
    return values


# ----------------------------------------------------------------------------
# Load the full local box from the packaged dust map, retaining raw densities.
# ----------------------------------------------------------------------------
def read_local_xyz_cube(parquet_path, x_range, y_range, z_range):
    """Read XYZ and dust columns row-group by row-group within inclusive bounds.

    Bounds and coordinates are in kpc. Accepted column aliases are normalized
    to x/y/z/dust; density values are retained in their input mag kpc^-1 units.
    No geometry mask is applied to the selected rectangular volume.
    The default raw product is written by query_and_write_raw_cube in
    00_Data_preparation_1.py: points beyond the queried dmax become NaN and
    finite densities below 1e-4 mag kpc^-1 become zero. That upstream operation
    does not remove structures according to fitted bubble positions.
    """
    parquet = pq.ParquetFile(str(parquet_path))
    names = parquet.schema.names
    if {"x", "y", "z", "dust"}.issubset(names):
        cols = ["x", "y", "z", "dust"]
    elif {"X", "Y", "Z", "dust"}.issubset(names):
        cols = ["X", "Y", "Z", "dust"]
    elif {"x", "y", "z", "dust_raw"}.issubset(names):
        cols = ["x", "y", "z", "dust_raw"]
    elif {"X", "Y", "Z", "dust_raw"}.issubset(names):
        cols = ["X", "Y", "Z", "dust_raw"]
    else:
        raise ValueError(f"Could not identify XYZ parquet columns: {names}")
    chunks = []
    for rg in range(parquet.num_row_groups):
        df = parquet.read_row_group(rg, columns=cols).to_pandas()
        df = df.rename(columns={"X": "x", "Y": "y", "Z": "z", "dust_raw": "dust"})
        mask = (
            (df["x"] >= x_range[0]) & (df["x"] <= x_range[1]) &
            (df["y"] >= y_range[0]) & (df["y"] <= y_range[1]) &
            (df["z"] >= z_range[0]) & (df["z"] <= z_range[1])
        )
        if mask.any():
            chunks.append(df.loc[mask, ["x", "y", "z", "dust"]].copy())
    if not chunks:
        raise ValueError("No XYZ data exist in the requested range")
    return pd.concat(chunks, ignore_index=True)


# ----------------------------------------------------------------------------
# Build the processed test grid and generate fields with IAAFT constraints.
# ----------------------------------------------------------------------------
def choose_rebin_factor(shape, max_voxels):
    """Choose an integer block width from the requested approximate voxel budget."""
    n_voxels = int(np.prod(shape))
    if n_voxels <= int(max_voxels):
        return 1
    factor = int(np.ceil((n_voxels / float(max_voxels)) ** (1.0 / 3.0)))
    return max(factor, 1)


def rebin_axis(values, factor):
    """Return block-centre coordinates, dropping an incomplete high-end block."""
    n = len(values) // factor
    trimmed = values[: n * factor]
    return trimmed.reshape(n, factor).mean(axis=1)


def rebin_cube_mean(cube, x_vals, y_vals, z_vals, factor):
    """Average complete factor-cubed blocks and their coordinate centres.

    Incomplete blocks at the high end of each axis are discarded. This changes
    both the available spatial resolution and the density distribution used by
    the null model; observed and surrogate statistics share the resulting grid.
    """
    if factor <= 1:
        return cube, x_vals, y_vals, z_vals
    nx = (cube.shape[0] // factor) * factor
    ny = (cube.shape[1] // factor) * factor
    nz = (cube.shape[2] // factor) * factor
    cube_trim = cube[:nx, :ny, :nz]
    cube_rebin = cube_trim.reshape(
        nx // factor, factor, ny // factor, factor, nz // factor, factor
    ).mean(axis=(1, 3, 5))
    return (
        cube_rebin,
        rebin_axis(x_vals, factor),
        rebin_axis(y_vals, factor),
        rebin_axis(z_vals, factor),
    )


def build_cube_context(cube_df, iaaft_max_voxels=1500000):
    """Build the *entire local* target field, including cavity and shell.

    Geometry masks are applied later, only for the test statistics. The FFT
    includes the DC term and all resolved long-wavelength modes of this box.
    Missing cells on the Cartesian product of available coordinate axes are
    filled from the nearest finite cell in voxel-index distance. The recorded
    missing fraction refers to this grid before any rebinning. Optional block
    averaging then defines the field whose values and spectrum are preserved;
    neither constraint refers to the unfilled, unreduced input table.
    Filled cells receive the same weight as original finite cells in both the
    target spectrum and statistics. Nearest-neighbour filling can extend
    coherent regions, and block averaging suppresses small-scale fluctuations;
    these preprocessing effects remain part of the conditional null field.

    Each individual 3D Fourier amplitude is constrained, without radial
    averaging or a fitted power-law slope. This is the spectrum of this target's
    local box, including environmental gradients. Conceptually writing
    rho = rho_background + delta_rho gives a total power containing both
    component powers and 2 Re(F_background * conj(F_delta)); the code does not
    estimate or separate these components. Low-frequency power consequently
    includes cavity, shell, environment and cross terms. Phase information also
    controls morphology, so retaining this power does not require a realization
    to retain the original cavity's location, shape or connectivity.

    Coordinates, voxel positions and flattened densities use matching C order
    so the later geometry masks address exactly the same voxel in each field.
    """
    x_vals = np.sort(np.unique(cube_df["x"].to_numpy(dtype=float)))
    y_vals = np.sort(np.unique(cube_df["y"].to_numpy(dtype=float)))
    z_vals = np.sort(np.unique(cube_df["z"].to_numpy(dtype=float)))
    nx, ny, nz = len(x_vals), len(y_vals), len(z_vals)
    x_index = {float(v): i for i, v in enumerate(x_vals)}
    y_index = {float(v): i for i, v in enumerate(y_vals)}
    z_index = {float(v): i for i, v in enumerate(z_vals)}

    cube = np.full((nx, ny, nz), np.nan, dtype=float)
    for row in cube_df.itertuples(index=False):
        cube[x_index[float(row.x)], y_index[float(row.y)], z_index[float(row.z)]] = float(row.dust)
    valid_mask = np.isfinite(cube)
    missing_fraction = float(1.0 - np.mean(valid_mask))
    if not np.isfinite(cube).all():
        if not np.any(valid_mask):
            raise ValueError('Invalid input or missing required data.')
        nearest_idx = distance_transform_edt(~valid_mask, return_distances=False, return_indices=True)
        cube = cube[tuple(nearest_idx)]

    rebin_factor = choose_rebin_factor(cube.shape, iaaft_max_voxels)
    cube, x_vals, y_vals, z_vals = rebin_cube_mean(cube, x_vals, y_vals, z_vals, rebin_factor)

    xx, yy, zz = np.meshgrid(x_vals, y_vals, z_vals, indexing="ij")
    return {
        "array": cube,
        "x_vals": x_vals,
        "y_vals": y_vals,
        "z_vals": z_vals,
        "points_flat": np.column_stack([xx.ravel(order="C"), yy.ravel(order="C"), zz.ravel(order="C")]),
        "values_flat": cube.ravel(order="C"),
        "sorted_values": np.sort(cube.ravel(order="C")),
        "target_amplitude": np.abs(np.fft.fftn(cube)),
        "shape": cube.shape,
        "missing_fraction_filled": missing_fraction,
        "rebin_factor": int(rebin_factor),
        "n_voxels": int(np.prod(cube.shape)),
    }


def generate_iaaft_surrogate_3d(cube_context, rng, max_iter=20, tol=1e-6):
    """Return the density-exact IAAFT stage; its spectrum is only approximate.

    ``tol`` controls the change in amplitude error, not the error itself.
    Starting with a random permutation preserves the empirical density PDF;
    Fourier phases then evolve through the alternating projections.
    Each inverse FFT imposes the target amplitudes; the following stable rank
    remapping restores the exact sorted densities and perturbs those amplitudes.
    The reported error is the mean absolute amplitude difference divided by
    the mean positive target amplitude. The helper default is 20 iterations;
    main explicitly passes the command-line default of eight iterations.
    This whole-spectrum amplitude error is not a relative power error for each
    mode and cannot alone establish fidelity at long wavelengths. Inspect the
    3D mode diagnostics as well as the radially averaged curves; a fixed
    iteration budget and a small change in error are not accuracy guarantees.
    """
    surrogate = rng.permutation(cube_context["sorted_values"]).reshape(cube_context["shape"], order="C")
    prev_error = np.inf
    final_error = np.inf
    iterations = 0
    target_amp = cube_context["target_amplitude"]
    sorted_values = cube_context["sorted_values"]

    for iteration in range(1, int(max_iter) + 1):
        spectrum = np.fft.fftn(surrogate)
        amplitude = np.abs(spectrum)
        phase = np.ones_like(spectrum, dtype=np.complex128)
        np.divide(spectrum, np.where(amplitude > 0, amplitude, 1.0), out=phase, where=amplitude > 0)
        adjusted = np.fft.ifftn(target_amp * phase).real

        ranks = np.argsort(adjusted.ravel(order="C"), kind="mergesort")
        surrogate_flat = np.empty(adjusted.size, dtype=float)
        surrogate_flat[ranks] = sorted_values
        surrogate = surrogate_flat.reshape(cube_context["shape"], order="C")

        current_amp = np.abs(np.fft.fftn(surrogate))
        denom = np.mean(target_amp[target_amp > 0]) if np.any(target_amp > 0) else 1.0
        final_error = float(np.mean(np.abs(current_amp - target_amp)) / max(denom, 1e-12))
        iterations = iteration
        if abs(prev_error - final_error) < tol:
            break
        prev_error = final_error

    return surrogate.ravel(order="C"), {"iaaft_iterations": int(iterations), "iaaft_amp_error": float(final_error)}


# ----------------------------------------------------------------------------
# Convert fitted geometries to fixed masks on the actual processed dust grid.
# ----------------------------------------------------------------------------
def rotate_points_to_local(points, angle_deg):
    """Undo the fitted XY rotation (degrees), leaving the Z coordinate unchanged."""
    theta = np.deg2rad(-float(angle_deg))
    c, s = np.cos(theta), np.sin(theta)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    return np.asarray(points, dtype=float) @ rot.T


def geometry_kind(row):
    """Select the slab-cylinder branch; other shape labels use a half ellipsoid."""
    return "cylinder" if str(row["shape"]) == "cylinder" else "ellipsoid_cap"


def cap_geometry(row):
    """Construct the cap-centred metric from fitted lengths, all in kpc.

    The metric centre lies in the fitted XY cap plane. Its transverse axes are
    the fitted cap-plane axes; the vertical axis is the cap-plane-to-apex
    distance. mark = 1 selects the lower apex, otherwise the upper apex.
    This metric, rather than the parent full ellipsoid, defines normalized u.
    """
    cx, cy, cz = float(row["center_x_kpc"]), float(row["center_y_kpc"]), float(row["center_z_kpc"])
    c_full = float(row["c_radius_kpc"])
    mark = int(row["mark"])
    cap_z = float(row["xy_plane_z_kpc"])
    a_cap = float(row["xy_plane_a_kpc"])
    b_cap = float(row["xy_plane_b_kpc"])
    z_apex = cz - c_full if mark == 1 else cz + c_full
    c_cap = abs(cap_z - z_apex)
    metric_center = np.array([cx, cy, cap_z], dtype=float)
    cap_axes = np.array([a_cap, b_cap, max(c_cap, 1e-9)], dtype=float)
    return metric_center, cap_axes, mark


def cube_ranges(row):
    """Set a surrounding rectangular FFT box from the fitted dimensions in kpc.

    XY half-width is at least 0.6 kpc or three times the largest transverse
    semi-axis. Z half-width is at least 0.35 kpc for cylinders and 0.25 kpc for
    caps, or three times their relevant vertical length. The box includes dust
    outside the fitted boundary and remains unchanged across realizations.
    """
    kind = geometry_kind(row)
    if kind == "cylinder":
        cx, cy, cz = float(row["center_x_kpc"]), float(row["center_y_kpc"]), float(row["center_z_kpc"])
        a, b = float(row["a_radius_kpc"]), float(row["b_radius_kpc"])
        xy_margin = max(3.0 * max(a, b), 0.6)
        z_center, z_margin = cz, max(3.0 * max(CYLINDER_HALF_HEIGHT_KPC, 1e-6), 0.35)
    else:
        metric_center, cap_axes, _ = cap_geometry(row)
        cx, cy, z_center = float(metric_center[0]), float(metric_center[1]), float(metric_center[2])
        xy_margin = max(3.0 * max(cap_axes[0], cap_axes[1]), 0.6)
        z_margin = max(3.0 * max(cap_axes[2], 1e-6), 0.25)
    x_range = (cx - xy_margin, cx + xy_margin)
    y_range = (cy - xy_margin, cy + xy_margin)
    z_range = (z_center - z_margin, z_center + z_margin)
    return x_range, y_range, z_range


def build_local(row, points):
    """Return local coordinates and inner/shell masks at the final fitted geometry.

    u is elliptical radius in XY for cylinders and ellipsoidal radius in XYZ
    for caps. The inner region is u <= 0.5 and the shell is 0.9 <= u <= 1.2.
    Cylinders retain |Z| <= 0.025 kpc; caps retain the upper half for mark = 2
    and the lower half otherwise. Local coordinates remain in kpc.
    """
    angle_deg = float(row["angle_deg"])
    kind = geometry_kind(row)
    points = np.asarray(points, dtype=float)

    if kind == "cylinder":
        cx, cy, cz = float(row["center_x_kpc"]), float(row["center_y_kpc"]), float(row["center_z_kpc"])
        a, b = float(row["a_radius_kpc"]), float(row["b_radius_kpc"])
        mark = int(row["mark"])
        local = rotate_points_to_local(points - np.array([cx, cy, cz])[None, :], angle_deg)
        X, Y, Z = local[:, 0], local[:, 1], local[:, 2]
        q = (X / max(a, 1e-9)) ** 2 + (Y / max(b, 1e-9)) ** 2
        u = np.sqrt(q)
        half_mask = np.abs(Z) <= CYLINDER_HALF_HEIGHT_KPC
        shell_mask = half_mask & (q >= SHELL_U_MIN ** 2) & (q <= SHELL_U_MAX ** 2)
        inner_mask = half_mask & (q <= INNER_U_MAX ** 2)
    else:
        metric_center, cap_axes, mark = cap_geometry(row)
        a_cap, b_cap, c_cap = cap_axes
        local = rotate_points_to_local(points - metric_center[None, :], angle_deg)
        X, Y, Z = local[:, 0], local[:, 1], local[:, 2]
        q = (X / max(a_cap, 1e-9)) ** 2 + (Y / max(b_cap, 1e-9)) ** 2 + (Z / max(c_cap, 1e-9)) ** 2
        u = np.sqrt(q)

        half_mask = Z >= 0.0 if mark == 2 else Z <= 0.0
        shell_mask = half_mask & (q >= SHELL_U_MIN ** 2) & (q <= SHELL_U_MAX ** 2)
        inner_mask = half_mask & (q <= INNER_U_MAX ** 2)
    return local, kind, mark, u, inner_mask, shell_mask, half_mask


def build_direction_bins(local, kind, mark, half_mask):
    """Partition accepted directions for radial ridge profiles.

    Cylinders use 36 azimuth bins. Caps use 36 azimuth by six cos(theta) bins
    on the selected hemisphere, giving equal solid angles in physical local
    coordinates; these are not angles after rescaling by the fitted axes.
    Excluded points receive sector -1; the cap centre has no defined direction.
    """
    X, Y, Z = local[:, 0], local[:, 1], local[:, 2]
    phi = (np.arctan2(Y, X) + 2.0 * np.pi) % (2.0 * np.pi)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, RIDGE_N_AZIMUTH + 1)
    phi_idx = np.clip(np.digitize(phi, phi_edges) - 1, 0, RIDGE_N_AZIMUTH - 1)

    if kind == "cylinder":
        # Within the cylinder slab, azimuth alone identifies a direction.
        sector_idx = np.where(half_mask, phi_idx, -1)
        return sector_idx.astype(int), RIDGE_N_AZIMUTH

    # Uniform cos(theta) intervals avoid oversampling directions near the pole.
    r = np.sqrt(X * X + Y * Y + Z * Z)
    safe = r > 1e-12
    cos_theta = np.where(safe, Z / np.maximum(r, 1e-12), 0.0)
    n_ct = RIDGE_N_COSTHETA_HALF
    ct_edges = np.linspace(0.0, 1.0, n_ct + 1) if mark == 2 else np.linspace(-1.0, 0.0, n_ct + 1)
    ct_idx = np.clip(np.digitize(cos_theta, ct_edges) - 1, 0, n_ct - 1)
    sector_idx = np.where(safe & half_mask, phi_idx * n_ct + ct_idx, -1)
    return sector_idx.astype(int), RIDGE_N_AZIMUTH * n_ct


# ----------------------------------------------------------------------------
# Estimate the shell's high-density component and the clipped cavity scatter.
# ----------------------------------------------------------------------------
def fit_two_gaussian_high_peak(values):
    """Fit two Gaussians in linear density and return the higher-mean component.

    Expectation-maximization starts at the density quartiles with equal weights
    and shared variance. Return its mean and standard deviation in mag kpc^-1,
    mixture weight and method/status label. Fewer than 30 finite samples or
    zero scatter yields NaNs. The higher component mean estimates shell ridge
    density; it is neither the largest voxel value nor a spatial peak location.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < RIDGE_PEAK_MIN_VOXELS or float(np.std(values)) <= 0:
        return float("nan"), float("nan"), float("nan"), "insufficient_shell_voxels"
    q25, q75 = np.percentile(values, [25, 75])
    means = np.array([q25, q75], dtype=float)
    shared_var = max(float(np.var(values)), 1e-8)
    variances = np.array([shared_var, shared_var], dtype=float)
    weights = np.array([0.5, 0.5], dtype=float)
    eps = 1e-12
    prev = -np.inf
    for _ in range(RIDGE_GMM_MAX_ITER):
        dens = []
        for k in range(2):
            var = max(float(variances[k]), 1e-8)
            coef = 1.0 / np.sqrt(2.0 * np.pi * var)
            dens.append(weights[k] * coef * np.exp(-0.5 * (values - means[k]) ** 2 / var))
        probs = np.vstack(dens).T
        total = np.sum(probs, axis=1) + eps
        resp = probs / total[:, None]
        nk = np.sum(resp, axis=0) + eps
        weights = nk / values.size
        means = np.sum(resp * values[:, None], axis=0) / nk
        variances = np.maximum(np.sum(resp * (values[:, None] - means[None, :]) ** 2, axis=0) / nk, 1e-8)
        loglike = float(np.sum(np.log(total)))
        if abs(loglike - prev) < RIDGE_GMM_TOL:
            break
        prev = loglike
    hi = int(np.argmax(means))
    return float(means[hi]), float(np.sqrt(variances[hi])), float(weights[hi]), "linear_2component_em_high_mean"


def iterative_gaussian_std(values, n_iter=INNER_STD_CLIP_ITERS, n_sigma=INNER_STD_CLIP_SIGMA):
    """Estimate inner-density scatter with up to three rounds of 3-sigma clipping.

    Each round compares all original finite values with the current retained
    mean and sample standard deviation (ddof = 1). Return retained scatter,
    mean and counts; units match the input density. This clipped scatter is a
    reported diagnostic, while the cavity mean used in tests remains unclipped.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    n_total = int(v.size)
    if n_total < 2:
        return float("nan"), float(np.mean(v)) if n_total else float("nan"), n_total, n_total
    keep = np.ones(n_total, dtype=bool)
    for _ in range(int(n_iter)):
        cur = v[keep]
        if cur.size < 2:
            break
        m = float(np.mean(cur))
        s = float(np.std(cur, ddof=1))
        if not np.isfinite(s) or s <= 0:
            break
        new_keep = np.abs(v - m) <= n_sigma * s
        if new_keep.sum() < 2 or np.array_equal(new_keep, keep):
            keep = new_keep if new_keep.sum() >= 2 else keep
            break
        keep = new_keep
    cur = v[keep]
    g_std = float(np.std(cur, ddof=1)) if cur.size > 1 else float("nan")
    g_mean = float(np.mean(cur)) if cur.size else float("nan")
    return g_std, g_mean, int(cur.size), n_total


# ----------------------------------------------------------------------------
# Cache geometry once, then measure identical statistics in every dust field.
# ----------------------------------------------------------------------------
def build_geometry_context(row, points):
    """Cache fixed masks, angular sectors and 16 shell-radius bins for one target."""
    local, kind, mark, u, inner_mask, shell_mask, half_mask = build_local(row, points)
    sector_idx, n_sectors_total = build_direction_bins(local, kind, mark, half_mask)
    u_edges = np.linspace(SHELL_U_MIN, SHELL_U_MAX, RIDGE_N_RADIAL + 1)
    u_centers = 0.5 * (u_edges[:-1] + u_edges[1:])
    u_idx = np.clip(np.digitize(u, u_edges) - 1, 0, RIDGE_N_RADIAL - 1)
    return {
        "kind": kind,
        "u": u,
        "inner_mask": inner_mask,
        "shell_mask": shell_mask,
        "sector_idx": sector_idx,
        "n_sectors_total": n_sectors_total,
        "u_idx": u_idx,
        "u_centers": u_centers,
    }


def compute_dust_metrics(geom, dust):
    """Measure cavity density, shell contrast and radial ridge alignment.

    Means use all finite voxels within each fixed mask. The global shell ridge
    density is the higher mean of a two-Gaussian fit to positive shell values;
    divide it by the positive, unclipped inner mean to obtain ridge contrast.
    Densities and their scatter have units mag kpc^-1; ratios are dimensionless.

    For spatial ridge alignment, retain sectors with at least ten shell voxels
    and radial bins with at least two voxels. Select the largest bin mean among
    16 bins over 0.9 <= u <= 1.2, accepting it only when its density is at least
    twice the inner mean. Average u_peak - 1 and |u_peak - 1| equally over valid
    sectors; no valid sector gives NaN. Although legacy output names end in
    ``over_R_eq``, the implemented offsets are dimensionless differences in u,
    not Cartesian distances explicitly divided by an equivalent radius.
    """
    dust = np.asarray(dust, dtype=float)
    inner_mask = geom["inner_mask"]
    shell_mask = geom["shell_mask"]
    sector_idx = geom["sector_idx"]
    n_sectors_total = geom["n_sectors_total"]
    u_idx = geom["u_idx"]
    u_centers = geom["u_centers"]

    finite_inner = inner_mask & np.isfinite(dust)
    finite_shell = shell_mask & np.isfinite(dust)
    inner_values = dust[finite_inner]
    shell_values = dust[finite_shell]

    inner_mean = float(np.mean(inner_values)) if inner_values.size else float("nan")
    inner_std, _, _, _ = iterative_gaussian_std(inner_values)
    shell_mean = float(np.mean(shell_values)) if shell_values.size else float("nan")
    shell_inner_ratio = (
        float(shell_mean / inner_mean)
        if np.isfinite(inner_mean) and inner_mean > 0 and np.isfinite(shell_mean) else float("nan")
    )

    ridge_peak_density = float("nan")
    ridge_peak_inner_ratio = float("nan")
    shell_pos = shell_values[np.isfinite(shell_values) & (shell_values > 0)]
    if shell_pos.size >= RIDGE_PEAK_MIN_VOXELS and float(np.nanstd(shell_pos)) > 0:
        ridge_peak_density, _, _, _ = fit_two_gaussian_high_peak(shell_pos)
        ridge_peak_inner_ratio = (
            float(ridge_peak_density / inner_mean)
            if np.isfinite(inner_mean) and inner_mean > 0 and np.isfinite(ridge_peak_density) else float("nan")
        )

    u_peak_by_sector = np.full(n_sectors_total, np.nan, dtype=float)
    for s in range(n_sectors_total):
        sm = finite_shell & (sector_idx == s)
        if int(sm.sum()) < RIDGE_MIN_SECTOR_VOXELS:
            continue
        radial_means = np.full(RIDGE_N_RADIAL, np.nan, dtype=float)
        su, sd = u_idx[sm], dust[sm]
        for j in range(RIDGE_N_RADIAL):
            bv = sd[su == j]
            if bv.size >= RIDGE_MIN_RADIAL_VOXELS:
                radial_means[j] = float(np.nanmean(bv))
        if not np.any(np.isfinite(radial_means)):
            continue
        pi = int(np.nanargmax(radial_means))
        peak_dust = float(radial_means[pi])
        peak_ratio = peak_dust / inner_mean if np.isfinite(inner_mean) and inner_mean > 0 else float("nan")
        if np.isfinite(peak_ratio) and peak_ratio >= RIDGE_PEAK_INNER_RATIO_MIN:
            u_peak_by_sector[s] = float(u_centers[pi])

    valid = u_peak_by_sector[np.isfinite(u_peak_by_sector)]
    mean_abs = float(np.mean(np.abs(valid - 1.0))) if valid.size else float("nan")
    mean_signed = float(np.mean(valid - 1.0)) if valid.size else float("nan")

    return {
        "dust_inner_voxel_count": int(inner_values.size),
        "dust_shell_voxel_count": int(shell_values.size),
        "dust_inner_mean_mag_kpc": inner_mean,
        "dust_inner_std_mag_kpc": inner_std,
        "dust_shell_mean_mag_kpc": shell_mean,
        "dust_shell_inner_ratio": shell_inner_ratio,
        "dust_shell_ridge_peak_density_mag_kpc": ridge_peak_density,
        "dust_shell_ridge_peak_inner_ratio": ridge_peak_inner_ratio,
        "ridge_direction_bin_count": int(n_sectors_total),
        "valid_ridge_sector_count": int(valid.size),
        "shell_radius_mean_signed_offset_over_R_eq": mean_signed,
        "shell_radius_mean_abs_offset_over_R_eq": mean_abs,
    }


# ----------------------------------------------------------------------------
# Count one-sided tail events, including their same-realization intersection.
# ----------------------------------------------------------------------------
def compute_pvalues(observed, surrogate_rows):
    """Return add-one empirical tail fractions at the fixed observed geometry.

    Count ratio >= observed, inner mean <= observed, and ridge offset <=
    observed. The three-way event is their intersection within a single row,
    not a product of marginal fractions. All supplied rows contribute to N,
    including rows whose missing statistics fail a finite-value comparison.
    A missing ridge offset thus counts as failing the morphology event, rather
    than being removed as unevaluable. With missing rows, retaining the full
    denominator yields a smaller tail fraction than using only evaluable rows.
    The add-one floor is finite Monte Carlo resolution, not zero probability.

    Missing observed offset gives NaN for offset and joint fractions. Missing
    observed ratio or inner mean gives zero events for its comparison and
    therefore the floor 1/(N+1); when the observed offset is finite, such a
    missing comparison also drives the joint fraction to that floor.

    A multivariate orthant tail fraction is not generally uniform under its
    null distribution, even when the component statistics are independent.
    No further calibration of p_dust_joint3 is performed here; its displayed
    threshold therefore does not specify an overall false-positive rate.
    """
    total = len(surrogate_rows)
    obs_ratio = observed["dust_shell_ridge_peak_inner_ratio"]
    obs_inner = observed["dust_inner_mean_mag_kpc"]
    obs_offset = observed["shell_radius_mean_abs_offset_over_R_eq"]

    n_ratio = sum(
        1 for r in surrogate_rows
        if np.isfinite(r["dust_shell_ridge_peak_inner_ratio"]) and np.isfinite(obs_ratio)
        and r["dust_shell_ridge_peak_inner_ratio"] >= obs_ratio
    )
    n_inner = sum(
        1 for r in surrogate_rows
        if np.isfinite(r["dust_inner_mean_mag_kpc"]) and np.isfinite(obs_inner)
        and r["dust_inner_mean_mag_kpc"] <= obs_inner
    )

    p_ratio = (1 + n_ratio) / (total + 1)
    p_inner = (1 + n_inner) / (total + 1)

    if np.isfinite(obs_offset):
        n_offset = sum(
            1 for r in surrogate_rows
            if np.isfinite(r["shell_radius_mean_abs_offset_over_R_eq"])
            and r["shell_radius_mean_abs_offset_over_R_eq"] <= obs_offset
        )
        n_joint3 = sum(
            1 for r in surrogate_rows
            if np.isfinite(r["dust_shell_ridge_peak_inner_ratio"]) and np.isfinite(obs_ratio)
            and r["dust_shell_ridge_peak_inner_ratio"] >= obs_ratio
            and np.isfinite(r["dust_inner_mean_mag_kpc"]) and np.isfinite(obs_inner)
            and r["dust_inner_mean_mag_kpc"] <= obs_inner
            and np.isfinite(r["shell_radius_mean_abs_offset_over_R_eq"])
            and r["shell_radius_mean_abs_offset_over_R_eq"] <= obs_offset
        )
        p_offset = (1 + n_offset) / (total + 1)
        p_joint3 = (1 + n_joint3) / (total + 1)
    else:
        p_offset = float("nan")
        p_joint3 = float("nan")

    return {
        "p_dust_ratio": p_ratio,
        "p_dust_inner_mean": p_inner,
        "p_shell_radius_mean_abs_offset_over_R_eq": p_offset,
        "p_dust_joint3": p_joint3,
    }


# ----------------------------------------------------------------------------
# Inspect full-field spectral fidelity and render matched diagnostic slices.
# ----------------------------------------------------------------------------
def power_spectrum_diagnostic(cube_context, surrogate_flat, long_wave_scale_kpc, n_bins=28):
    """Compare full-box 3D spectra, with a radial average for display only.

    Frequency f is in cycles/kpc (not angular wavenumber); lambda = 1/f.
    P = voxel_volume * abs(FFT(dust))**2 / N. With the DC mode omitted,
    sum(P) / box_volume equals the population variance (Parseval).
    Long modes have 0 < f <= 1 / long_wave_scale_kpc. Their diagnostic uses
    individual 3D modes, not the radially averaged curve.
    Uniform spacing on all three axes is required to assign physical frequency
    and power units: (mag/kpc)^2 kpc^3 = mag^2 kpc. Geometric frequency bins
    summarize mode-averaged power; their centres are geometric means of the
    populated frequencies. The FFT treats opposite box faces as periodic;
    without a taper, boundary mismatches and large-scale gradients contribute
    to the spectrum. Along each axis the box period n * spacing sets the
    longest separately resolved nonzero wavelength; modes larger than the
    local box are not recovered by preserving its low-frequency power.

    These spectra come from the full 3D fields, not the displayed 2D slices.
    Radial averages can conceal directional mode errors. The long-wave total
    power ratio sums individual modes and is not the mean of plotted ratios.
    The accompanying mode-wise relative L1 power error prevents positive and
    negative mode discrepancies from cancelling in that total. Plot bins may
    straddle the long-wave cutoff, so their sums need not reproduce the exact
    mode-based long-wave diagnostics.
    """
    if not np.isfinite(long_wave_scale_kpc) or long_wave_scale_kpc <= 0 or n_bins < 1:
        raise ValueError("Use a positive long-wave scale and at least one spectral bin")
    shape = cube_context["shape"]
    spacings = []
    for key in ("x_vals", "y_vals", "z_vals"):
        axis = cube_context[key]
        delta = np.diff(axis)
        if delta.size == 0 or np.any(delta <= 0) or not np.allclose(delta, np.mean(delta), rtol=1e-4, atol=1e-8):
            raise ValueError("Physical power spectra require at least two uniformly spaced voxels on each axis")
        spacings.append(float(np.mean(delta)))
    frequencies = [np.fft.fftfreq(n, d=d) for n, d in zip(shape, spacings)]
    f = np.sqrt(frequencies[0][:, None, None] ** 2 +
                frequencies[1][None, :, None] ** 2 +
                frequencies[2][None, None, :] ** 2).ravel()
    target_amp = cube_context["target_amplitude"].ravel()
    random_amp = np.abs(np.fft.fftn(np.asarray(surrogate_flat).reshape(shape))).ravel()
    voxel_volume = float(np.prod(spacings))
    normalization = voxel_volume / int(np.prod(shape))
    original_power = target_amp ** 2 * normalization
    random_power = random_amp ** 2 * normalization
    non_dc = f > 0
    f_non_dc = f[non_dc]
    edges = np.geomspace(f_non_dc.min(), f_non_dc.max(), int(n_bins) + 1)
    bin_idx = np.clip(np.searchsorted(edges, f_non_dc, side="right") - 1, 0, n_bins - 1)
    counts = np.bincount(bin_idx, minlength=n_bins)
    populated = counts > 0
    original_binned = np.bincount(bin_idx, weights=original_power[non_dc], minlength=n_bins)[populated] / counts[populated]
    random_binned = np.bincount(bin_idx, weights=random_power[non_dc], minlength=n_bins)[populated] / counts[populated]
    table = pd.DataFrame({
        "frequency_min_cycles_kpc": edges[:-1][populated],
        "frequency_max_cycles_kpc": edges[1:][populated],
        "frequency_cycles_kpc": np.exp(np.bincount(bin_idx, weights=np.log(f_non_dc), minlength=n_bins)[populated] / counts[populated]),
        "mode_count": counts[populated],
        "original_power_mag2_kpc": original_binned,
        "surrogate_power_mag2_kpc": random_binned,
        "surrogate_to_original_power": np.divide(random_binned, original_binned, out=np.full_like(original_binned, np.nan), where=original_binned > 0),
    })
    low = non_dc & (f <= 1.0 / long_wave_scale_kpc)

    def ratio(numerator, denominator):
        return float(numerator / denominator) if denominator > 0 else None

    diagnostics = {
        "spectrum_source": "entire local 3D dust cube: cavity, shell and surroundings; after fill/rebin",
        "power_normalization": "P = voxel_volume * abs(fftn(dust))**2 / N; units (mag/kpc)^2 kpc^3 = mag^2 kpc",
        "frequency_convention": "cycles/kpc; wavelength = 1/f; no taper or detrending; DC omitted only in diagnostics",
        "spacing_kpc": spacings,
        "box_period_kpc": [float(n * d) for n, d in zip(shape, spacings)],
        "long_wave_scale_kpc": float(long_wave_scale_kpc),
        "long_wave_mode_count": int(low.sum()),
        "long_wave_power_ratio": ratio(random_power[low].sum(), original_power[low].sum()),
        "long_wave_mode_power_relative_l1": ratio(np.abs(random_power[low] - original_power[low]).sum(), original_power[low].sum()),
        "non_dc_amplitude_relative_l1": ratio(np.abs(random_amp[non_dc] - target_amp[non_dc]).sum(), target_amp[non_dc].sum()),
        "non_dc_mode_power_relative_l1": ratio(np.abs(random_power[non_dc] - original_power[non_dc]).sum(), original_power[non_dc].sum()),
        "original_variance_from_power": float(original_power[non_dc].sum() / (voxel_volume * np.prod(shape))),
        "surrogate_variance_from_power": float(random_power[non_dc].sum() / (voxel_volume * np.prod(shape))),
        "density_values_match_exactly": bool(np.array_equal(np.sort(surrogate_flat), cube_context["sorted_values"])),
    }
    return table, diagnostics


def save_realization_comparison(row, cube_context, surrogate_flat, meta, seed,
                                output_path, diagnostic_dir, provenance=None):
    """Save matched XY/XZ/YZ slices plus 3D spectrum and density checks.

    Both maps use identical colour limits and a display-only Gaussian with
    sigma = 1 test-grid voxel on each in-plane axis, applied to each 2D slice.
    Spectra, PDFs and saved cubes use the data before display smoothing.
    This is the first realization, without selection by appearance or score.
    The PNG uses coordinate-aligned slices at the nearest grid positions to
    the fitted centre (halfway from cap plane to apex for the cap XY slice).
    Shared colour limits span the original minimum to its 99.5th percentile.
    Values above the colour limit saturate only in the maps; all density values
    remain in the spectrum and histogram. The reflect-mode Gaussian is cut off
    at four sigma, with no upsampling; its physical width follows the processed
    voxel spacing and is recorded separately for each coordinate axis.
    The descriptive long-wave cutoff is twice the largest fitted semi-axis;
    it does not isolate the cavity's contribution to the full-box spectrum.
    NPZ, CSV and JSON products retain the unsmoothed fields, binned spectrum,
    numerical diagnostics, grid properties and provenance for this same draw.

    The bottom panels show mean power by frequency, surrogate/original power
    ratio (unity means matching bin means), and the density PDF on a log axis.
    Low frequency corresponds to large spatial scales. Overlapping PDFs follow
    from exact rank remapping and do not imply identical spatial morphology.
    The shaded long-wave region marks wavelength >= twice the largest fitted
    semi-axis. It is a reproducible scale marker, not a spectral decomposition
    into cavity and environment.
    """
    target_id = int(row["id"])
    original = cube_context["array"]
    surrogate = np.asarray(surrogate_flat).reshape(cube_context["shape"])
    axes = [cube_context[k] for k in ("x_vals", "y_vals", "z_vals")]
    display_sigma_kpc_xyz = [float(COMPARISON_SMOOTH_SIGMA_VOXELS * np.mean(np.diff(axis))) for axis in axes]
    if geometry_kind(row) == "cylinder":
        slice_center = np.array([row["center_x_kpc"], row["center_y_kpc"], row["center_z_kpc"]], dtype=float)
        radii = np.array([row["a_radius_kpc"], row["b_radius_kpc"], CYLINDER_HALF_HEIGHT_KPC], dtype=float)
    else:
        slice_center, radii, mark = cap_geometry(row)
        # XY halfway between the cap plane and apex; XZ/YZ through its centre.
        slice_center[2] += (1 if mark == 2 else -1) * radii[2] / 2.0
    indices = [int(np.argmin(np.abs(axis - coordinate))) for axis, coordinate in zip(axes, slice_center)]
    long_wave_scale = 2.0 * float(np.max(radii))
    spectrum, diagnostics = power_spectrum_diagnostic(cube_context, surrogate_flat, long_wave_scale)
    # Show the actual fitted-volume boundary on both realizations. It is not refit.
    _, _, _, u, _, _, half_mask = build_local(row, cube_context["points_flat"])
    fitted_volume = ((u <= 1.0) & half_mask).reshape(original.shape)
    vmin = float(np.min(original))
    vmax = float(np.percentile(original, 99.5))
    if vmax <= vmin:
        vmax = vmin + 1.0
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    slice_specs = [(0, 1, 2, "XY"), (0, 2, 1, "XZ"), (1, 2, 0, "YZ")]
    spans = [float(len(axis) * np.mean(np.diff(axis))) for axis in axes]
    map_heights = [max(0.25, spans[v] / spans[h]) for h, v, _, _ in slice_specs]
    fig = plt.figure(figsize=(11.8, 5.0 * sum(map_heights) + 3.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[sum(map_heights) + 0.25, 0.65])
    spatial_grid = grid[0].subgridspec(3, 2, height_ratios=map_heights, hspace=0.05)
    diagnostic_grid = grid[1].subgridspec(1, 3, wspace=0.10)
    map_axes = []
    labels = ["X", "Y", "Z"]
    for row_index, (horizontal, vertical, fixed, plane) in enumerate(slice_specs):
        for col, (field, field_title) in enumerate([(original, "Original"), (surrogate, "Random")]):
            ax = fig.add_subplot(spatial_grid[row_index, col])
            map_axes.append(ax)
            image_values = np.take(field, indices[fixed], axis=fixed).T
            boundary = np.take(fitted_volume, indices[fixed], axis=fixed).T
            h, v = axes[horizontal], axes[vertical]
            dh, dv = np.mean(np.diff(h)), np.mean(np.diff(v))
            extent = [h[0] - dh / 2, h[-1] + dh / 2, v[0] - dv / 2, v[-1] + dv / 2]
            # Smooth only the displayed slice, at the test grid's voxel scale.
            image_values = gaussian_filter(
                image_values, sigma=COMPARISON_SMOOTH_SIGMA_VOXELS,
                mode="reflect", truncate=4.0,
            )
            im = ax.imshow(image_values, origin="lower", extent=extent, cmap=COMPARISON_CMAP, norm=norm, interpolation="nearest", aspect="equal")
            if boundary.any() and not boundary.all():
                ax.contour(h, v, boundary.astype(float), levels=[0.5], colors=["black"], linewidths=1.0, linestyles="--")
            panel = chr(ord("a") + row_index * 2 + col)
            ax.set_title(f"({panel})  {field_title}: {plane}", fontsize=11, loc="left", pad=5)
            ax.set_xlabel(f"{labels[horizontal]} (kpc)")
            ax.set_ylabel(f"{labels[vertical]} (kpc)")
            ax.tick_params(labelsize=9)
    fig.colorbar(im, ax=map_axes, shrink=0.75, pad=0.02, fraction=0.025, extend="max", label=r"Dust density (mag kpc$^{-1}$)")
    ax_power, ax_ratio, ax_pdf = [fig.add_subplot(diagnostic_grid[0, i]) for i in range(3)]
    f = spectrum["frequency_cycles_kpc"].to_numpy()
    colors = ["#2D6A9F", "#C86432"]
    for name, column, color, ls in [("Original", "original_power_mag2_kpc", colors[0], "-"), ("Random", "surrogate_power_mag2_kpc", colors[1], "--")]:
        ax_power.loglog(f, spectrum[column], color=color, linestyle=ls, linewidth=1.7, label=name)
    ax_power.set_title("(g)  Power spectrum", fontsize=11, loc="left")
    ax_power.set_ylabel(r"$P(f)$ (mag$^2$ kpc)")
    ax_power.legend(frameon=False)
    ax_ratio.semilogx(f, spectrum["surrogate_to_original_power"], color=colors[1], linewidth=1.6)
    ax_ratio.axhline(1.0, color="0.25", linestyle="--", linewidth=1)
    ax_ratio.set_ylabel(r"$P_{\rm random}/P_{\rm original}$")
    ax_ratio.set_ylim(bottom=0)
    ax_ratio.set_title("(h)  Power ratio", fontsize=11, loc="left")
    for ax in (ax_power, ax_ratio):
        cutoff = 1.0 / long_wave_scale
        if cutoff >= f.min():
            ax.axvspan(f.min(), min(cutoff, f.max()), color="0.88", alpha=0.55, zorder=0)
        ax.set_xlabel(r"$f$ (cycles kpc$^{-1}$)")
        ax.grid(alpha=0.15)
    bins = np.histogram_bin_edges(original.ravel(), bins=60)
    for values, label, color, ls in [(original, "Original", colors[0], "-"), (surrogate, "Random", colors[1], "--")]:
        ax_pdf.hist(values.ravel(), bins=bins, density=True, histtype="step", color=color, linestyle=ls, linewidth=1.6, label=label)
    ax_pdf.set_yscale("log")
    ax_pdf.set_title("(i)  Density distribution", fontsize=11, loc="left")
    ax_pdf.set_xlabel(r"Dust density (mag kpc$^{-1}$)")
    ax_pdf.set_ylabel("Probability density")
    ax_pdf.legend(frameon=False)
    fig.suptitle(f"SB{target_id}", fontsize=14, fontweight="bold")
    output_path = Path(output_path)
    diagnostic_dir = Path(diagnostic_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    diagnostic_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.10)
    plt.close(fig)
    stem = f"8_SB{target_id}_realization_comparison"
    spectrum.to_csv(diagnostic_dir / f"{stem}_power.csv", index=False, encoding="utf-8-sig")
    np.savez_compressed(diagnostic_dir / f"{stem}.npz", original=original, surrogate=surrogate,
                        x=axes[0], y=axes[1], z=axes[2])
    diagnostics.update({
        "id": target_id, "seed": int(seed), "mc_index": 0, **meta,
        "cube_shape": list(original.shape), "rebin_factor": cube_context["rebin_factor"],
        "missing_fraction_filled": cube_context["missing_fraction_filled"],
        "slice_indices_xyz": indices,
        "slice_coordinates_kpc": [float(axis[index]) for axis, index in zip(axes, indices)],
        "map_color_limits_mag_kpc": [vmin, vmax],
        "map_colormap": COMPARISON_CMAP,
        "figure_layout": "rows XY/XZ/YZ; columns original/random; bottom power spectrum, power ratio, density PDF; full caption metadata retained here",
        "comparison_figure_png": metadata_path(output_path),
        "path_base": "code/",
        "map_display_smoothing": {
            "sigma_voxels": COMPARISON_SMOOTH_SIGMA_VOXELS,
            "sigma_kpc_xyz": display_sigma_kpc_xyz,
            "grid": "actual test grid after missing-value filling and optional rebinning",
            "method": "2D gaussian_filter on each native test-grid slice, mode=reflect, truncate=4; no upsampling",
            "scope": "six displayed spatial slices only; spectra, density histograms, saved cubes and test statistics precede display smoothing",
        },
        "long_wave_scale_definition": "twice largest fitted semi-axis (cap axes or cylinder radii/half height); descriptive, not a spectral decomposition of the cavity",
        "realization_selection": "first IAAFT realization from seed + target ID; no selection on appearance or statistics",
        "input_provenance": provenance or {},
    })
    with (diagnostic_dir / f"{stem}.json").open("w", encoding="utf-8") as handle:
        json.dump(diagnostics, handle, indent=2, ensure_ascii=False, allow_nan=False)
    return diagnostics


def style_publication_axis(ax, panel_label=None):
    """Apply consistent frame/tick weights and an optional inset panel letter."""
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.0)
    ax.tick_params(length=3.6, width=0.9, pad=2.2)
    if panel_label is not None:
        ax.text(
            0.015, 0.96, panel_label,
            transform=ax.transAxes,
            ha="left", va="top",
            fontsize=12.0, fontweight="bold",
        )


def save_summary_figure(summary_df, output_path):
    """Plot four tail fractions by target and counts in the joint-p bins.

    Reference lines at 0.05 and 0.12 match the displayed count categories.
    Missing joint values do not enter any category, so counts need not sum to
    the number of requested targets. The output is a 300 dpi PNG in main.
    """
    df = summary_df.sort_values("id").reset_index(drop=True)
    labels = [f"SB{int(v)}" for v in df["id"]]
    x = np.arange(len(df))

    fig = plt.figure(figsize=(10.0, 10.0), constrained_layout=True)
    grid = fig.add_gridspec(5, 2, height_ratios=[1.0, 1.0, 1.0, 1.0, 1.25])
    p_axes = [fig.add_subplot(grid[i, :]) for i in range(4)]
    ax_stat = fig.add_subplot(grid[4, 0])
    legend_ax = fig.add_subplot(grid[4, 1])
    p_panels = [
        ("p_dust_ratio", "Shell ridge / inner dust"),
        ("p_dust_inner_mean", "Inner dust mean"),
        ("p_shell_radius_mean_abs_offset_over_R_eq", "Mean |u_peak - 1|"),
        ("p_dust_joint3", "Joint p-value"),
    ]
    panel_labels = ["a", "b", "c", "d", "e"]
    for ax, (column, title), panel_label in zip(p_axes, p_panels, panel_labels[:4]):
        ax.bar(x, df[column], color=BAR_COLOR, edgecolor="0.20", linewidth=0.35)
        ax.axhline(0.05, color="#B2182B", linestyle="--", linewidth=0.9)
        ax.axhline(0.12, color="#E69F00", linestyle=":", linewidth=0.9)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel("Empirical p-value")
        ax.set_xlabel("")
        ax.set_title(title, pad=5.0)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right", rotation_mode="anchor")
        style_publication_axis(ax, panel_label)

    pj = pd.to_numeric(df["p_dust_joint3"], errors="coerce").to_numpy(dtype=float)
    n_sig = int(np.sum(pj < 0.05))
    n_marg = int(np.sum((pj >= 0.05) & (pj < 0.12)))
    n_ns = int(np.sum(pj >= 0.12))
    ax_stat.bar(["<0.05", "0.05-0.12", ">=0.12"], [n_sig, n_marg, n_ns],
                color=BAR_COLOR, edgecolor="0.20", linewidth=0.35)
    ax_stat.set_title("Joint p-value counts", pad=5.0)
    ax_stat.set_xlabel("Joint p-value bin")
    ax_stat.set_ylabel("Number of SBs")
    ax_stat.set_ylim(0, max(n_sig, n_marg, n_ns) + 5)
    for i, v in enumerate([n_sig, n_marg, n_ns]):
        ax_stat.text(i, v + 0.2, str(v), ha="center", va="bottom", fontsize=11.0)
    style_publication_axis(ax_stat, panel_labels[4])

    legend_ax.axis("off")
    legend_handles = [
        mpl.lines.Line2D([0], [0], color="#B2182B", linestyle="--", linewidth=0.9, label="p = 0.05"),
        mpl.lines.Line2D([0], [0], color="#E69F00", linestyle=":", linewidth=0.9, label="p = 0.12"),
    ]
    legend_ax.legend(
        handles=legend_handles,
        loc="center",
        frameon=False,
        handlelength=2.4,
        borderaxespad=0.0,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Orchestrate per-target random streams and write diagnostic/statistical products.
# ----------------------------------------------------------------------------
def main():
    """Run N = 200 fixed-geometry surrogates per target with an eight-step default.

    seed + target ID gives each object its own reproducible random stream,
    independent of target processing order. Geometry is built once after
    preprocessing and reused for the observation and every realization.
    Diagnostics-only mode generates one draw and exits before sample/summary
    output; a normal run also records every draw's metrics and IAAFT error.
    """
    parser = argparse.ArgumentParser(
        description='Run Script 8: OSBs random dust test.'
    )
    parser.add_argument("--targets", default="all", help='Comma-separated target IDs; use all for the full sample.')
    parser.add_argument("--n-mc", type=int, default=200, help='Number of random or Monte Carlo realizations.')
    parser.add_argument("--seed", type=int, default=20260516, help='Random seed for reproducible controls.')
    parser.add_argument("--geom-csv", default=str(DEFAULT_GEOM_CSV))
    parser.add_argument("--xy-data-path", default=str(DEFAULT_XY_DATA_PATH))
    parser.add_argument("--iaaft-max-iter", type=int, default=8)
    parser.add_argument("--iaaft-max-voxels", type=int, default=1500000)
    parser.add_argument("--diagnostics-only", action="store_true",
                        help="Generate only the first realization and comparison per target; do not write p-value/sample tables.")
    parser.add_argument("--diagnostic-dir", default=str(OUT_DIR / "diagnostics"),
                        help="Directory for comparison cubes, power spectra and metadata.")
    parser.add_argument("--comparison-fig-dir", default=str(COMPARISON_FIG_DIR),
                        help="Directory for the original-versus-random comparison figures.")
    args = parser.parse_args()
    if args.n_mc < 1 or args.iaaft_max_iter < 1 or args.iaaft_max_voxels < 1:
        parser.error("--n-mc, --iaaft-max-iter and --iaaft-max-voxels must be positive")

    geom_csv = ensure_existing_file(args.geom_csv, "final geometry parameter table")
    xy_data_path = ensure_existing_file(args.xy_data_path, "raw 3D XYZ dust cube")
    geom_df = pd.read_csv(geom_csv, encoding="utf-8-sig").sort_values("id").reset_index(drop=True)
    targets = parse_targets(args.targets, geom_df)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)
    diagnostic_dir = Path(args.diagnostic_dir)
    comparison_fig_dir = Path(args.comparison_fig_dir)
    provenance = {
        "geometry_csv": metadata_path(geom_csv),
        "dust_parquet": metadata_path(xy_data_path),
        "iaaft_max_iter": int(args.iaaft_max_iter),
        "iaaft_max_voxels": int(args.iaaft_max_voxels),
        "diagnostics_only": bool(args.diagnostics_only),
    }

    summary_rows = []
    for target_id in targets:
        start_time = time.time()
        rng = np.random.default_rng(int(args.seed) + int(target_id))
        row_df = geom_df[geom_df["id"] == int(target_id)]
        if row_df.empty:
            raise ValueError('Invalid input or missing required data.')
        row = row_df.iloc[0]
        kind = geometry_kind(row)

        xr, yr, zr = cube_ranges(row)
        cube_df = read_local_xyz_cube(xy_data_path, xr, yr, zr)
        cube_context = build_cube_context(cube_df, iaaft_max_voxels=args.iaaft_max_voxels)
        points = cube_context["points_flat"]
        geom = build_geometry_context(row, points)
        observed = compute_dust_metrics(geom, cube_context["values_flat"])

        surrogate_rows = []
        iaaft_iters, iaaft_errors = [], []
        n_realizations = 1 if args.diagnostics_only else int(args.n_mc)
        for mc_index in range(n_realizations):
            dust_flat, meta = generate_iaaft_surrogate_3d(cube_context, rng, max_iter=args.iaaft_max_iter)
            if mc_index == 0:
                comparison_path = comparison_fig_dir / f"8_SB{int(target_id)}_original_vs_random.png"
                diagnostics = save_realization_comparison(
                    row, cube_context, dust_flat, meta, int(args.seed) + int(target_id),
                    comparison_path, diagnostic_dir,
                    provenance={**provenance, "requested_ranges_kpc": [xr, yr, zr]},
                )
                print(f"SB{int(target_id)} comparison: {comparison_path}; "
                      f"long-wave power ratio={diagnostics['long_wave_power_ratio']}", flush=True)
            if args.diagnostics_only:
                continue
            sur = compute_dust_metrics(geom, dust_flat)
            sur["mc_index"] = int(mc_index)
            sur["iaaft_iterations"] = int(meta["iaaft_iterations"])
            sur["iaaft_amp_error"] = float(meta["iaaft_amp_error"])
            surrogate_rows.append(sur)
            iaaft_iters.append(meta["iaaft_iterations"])
            iaaft_errors.append(meta["iaaft_amp_error"])

        if args.diagnostics_only:
            print(f"SB{int(target_id)} diagnostic completed in {time.time() - start_time:.1f}s", flush=True)
            continue

        pvalues = compute_pvalues(observed, surrogate_rows)
        prefix = f"8_SB{int(target_id)}_surrogate_samples"
        pd.DataFrame(surrogate_rows).to_csv(
            OUT_DIR / f"{prefix}.csv", index=False, encoding="utf-8-sig"
        )
        summary_rows.append({
            "id": int(target_id),
            "shape": str(row["shape"]),
            "mark": int(row["mark"]),
            "final_parameter_estimator": str(row["final_parameter_estimator"]),
            "geometry_kind": kind,
            "n_mc": int(args.n_mc),
            "iaaft_max_voxels": int(args.iaaft_max_voxels),
            "cube_rebin_factor": int(cube_context["rebin_factor"]),
            "cube_n_voxels": int(cube_context["n_voxels"]),
            "cube_missing_fraction_filled": float(cube_context["missing_fraction_filled"]),
            "ridge_direction_bin_count": int(observed["ridge_direction_bin_count"]),
            **{f"observed_{k}": v for k, v in observed.items()},
            **pvalues,
            "dust_iaaft_iterations_mean": float(np.mean(iaaft_iters)) if iaaft_iters else float("nan"),
            "dust_iaaft_amp_error_mean": float(np.mean(iaaft_errors)) if iaaft_errors else float("nan"),
            "elapsed_sec": float(time.time() - start_time),
        })
        print(
            f"SB{int(target_id):>2} [{kind:>14}]  "
            f"p_ratio={pvalues['p_dust_ratio']:.4f}  "
            f"p_inner={pvalues['p_dust_inner_mean']:.4f}  "
            f"p_offset={pvalues['p_shell_radius_mean_abs_offset_over_R_eq']!s:>6}  "
            f"p_joint3={pvalues['p_dust_joint3']!s:>6}  "
            f"valid_sec={observed['valid_ridge_sector_count']:>3}/{observed['ridge_direction_bin_count']}  "
            f"elapsed={time.time() - start_time:.1f}s",
            flush=True,
        )

    if args.diagnostics_only:
        return
    summary_df = pd.DataFrame(summary_rows).sort_values("id").reset_index(drop=True)
    target_tag = "ALL30" if len(targets) == len(geom_df) else "SB" + "_".join(str(v) for v in targets)
    summary_csv = OUT_DIR / f"8_condition_test_summary_{target_tag}.csv"
    joint_csv = OUT_DIR / f"8_condition_test_joint_p_values_{target_tag}.csv"
    fig_tag = "all30" if len(targets) == len(geom_df) else "sb" + "_".join(str(v) for v in targets)
    summary_png = FINAL_FIG_DIR / f"8_condition_test_joint_p_{fig_tag}.png"
    summary_df.to_csv(summary_csv, index=False, encoding="utf-8-sig")
    summary_df[[
        "id", "geometry_kind", "p_dust_ratio", "p_dust_inner_mean",
        "p_shell_radius_mean_abs_offset_over_R_eq", "p_dust_joint3",
    ]].to_csv(joint_csv, index=False, encoding="utf-8-sig")
    save_summary_figure(summary_df, summary_png)
    print(f"\nsaved summary csv: {summary_csv}", flush=True)
    print(f"saved joint csv:   {joint_csv}", flush=True)
    print(f"saved summary fig: {summary_png}", flush=True)

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
