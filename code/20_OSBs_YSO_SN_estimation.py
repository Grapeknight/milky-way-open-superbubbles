#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script 20: orbit-gated cluster supernova supply and local-SFR comparison.

All calculation, geometry, orbit, IMF, PARSEC and plotting functions are
contained in this numbered script, grouped by their role in the calculation.
The method below describes the default configuration; the saved run parameters
record the settings and input checksums used by each execution.

Script 20 compares two independent supply estimates with the requirements from
Script 19: supernovae associated with known clusters and an area–time estimate
based on the local star-formation rate. These estimates are not added together.

## Inputs and execution

Run `python 20_OSBs_YSO_SN_estimation.py` from `code/` (the script also resolves
its working directory when launched elsewhere). The normal run requires no
network connection. Python dependencies are NumPy, pandas, SciPy, matplotlib,
Astropy, galpy and Numba. The validated environment used NumPy 1.24.3 and
galpy 1.11.2. IMF and orbit caches are generated under
`../results/intermediate_output/20_cluster_supernova_count_in_superbubble_volume/cache/`.
Allow about 3.5 GB for the IMF libraries, plus the orbit cache and output tables.
The release omits these regenerable caches. The first execution rebuilds them
from the packaged inputs and fixed seeds.

Cluster properties come from `../data/star_cluster_data/hunt2024_clusters_full.csv`.
Member photometry and the PARSEC grid are packaged in
`../data/star_cluster_data/sn_population/`; `INPUTS.md` in that directory gives
their sources, query settings and checksums. Geometry, dynamical ages and SN
requirements are read from the current outputs of Scripts 3 and 19. The compact
catalogue also consumes the existing Script 5 and 8 products.

## Cluster population and explosion times

The valid sample contains 6,956 Hunt & Reffert (2024) clusters with
finite positions, ages and positive catalogue masses. Of these, 5,607 have
finite equatorial coordinates, proper motions and radial velocities; only
these are eligible to supply supernovae. The 1,349 without complete 6D data
are excluded. No extra CST, CMDCl50, radial-velocity error or membership-count
cut is applied. Positive catalogue mass defines the sample but is
not used to normalize the IMF. Stable distance-sorted catalogue row identifiers
determine the cluster-specific random seeds.

Following the public Methods of [Swiggum et al. (2024)](https://doi.org/10.1038/s41586-024-07496-9),
the reconstruction matches the observable stellar population to a stochastic
IMF. The implementation uses the following explicit choices:

- Select catalogue member rows with apparent Gaia G between 12 and 17 and
  positive finite parallax. Retain the original cluster assignments and read
  Gaia identifiers as strings. Shared members can occur across clusters.
- Use PARSEC v1.2S/COLIBRI, Z = 0.0158, Gaia EDR3 passbands and a log-age grid
  spaced by 0.01. Convert each star's magnitude with its parallax and the
  cluster AV50; adopt constant Cardelli Rv = 3.1 coefficients AG/AV = 0.83627
  and ARP/AV = 0.63439.
- Project onto the densely interpolated full colour–magnitude isochrone to
  infer initial mass. Determine the union of all observable PMS/MS initial-mass
  intervals from the G limits and cluster distance. Match the same interval
  union in observations and models, including disconnected intervals caused
  by non-monotonic isochrones.
- Draw Kroupa IMFs from 0.08 to 150 solar masses, with a break at 0.5 and slopes
  1.3/2.3. The main library has 49,900 target masses from 10 to 4,999.9 solar
  masses, step 0.1, with stop-nearest sampling. Only fits outside the supported
  mass range use a combined library with 2,501 additional targets from 5,000
  to 30,000 solar masses, step 10. Preserve every sampled stellar mass.
- For each cluster, repeat the fit 100 times, each time selecting 1,000 distinct
  candidate models and minimizing the observed complete-interval mass
  residual. The base seed is 20260912, extension seed 20260913, and cluster seed
  20260912 + row_id. Relative fit residuals above 10% are flagged in diagnostics.
- Count actual sampled initial masses at least 8 solar masses. Infer explosion
  ages from the maximum surviving initial mass along the PARSEC grid, including
  evolved massive stars. The adopted lifetime range is approximately 2.977 to
  40.997 Myr for 150 to 8 solar masses.

For cluster age A and stellar lifetime tau, lookback time is A − tau. Count an
event only when `0 <= A - tau <= T_dyn` and its trajectory is inside the adopted
bubble boundary then. An old cluster may contribute a recent explosion; it is
not discarded merely because its age exceeds the bubble age.

## Orbits and evolving boundary

Integrate measured cluster astrometry and RV in MWPotential2014 with galpy's
`odeint` method, a 0.05 Myr lookback grid, R0 = 8.122 kpc, V0 = 236 km/s,
z_sun = 0.0208 kpc and `solarmotion = [-11.1, 12.24, 7.25]` km/s. Only the
PARSEC-active age cohort needs integration. Bubble centres start with the local
circular tangential speed and zero radial/vertical speed; axes co-rotate with
the centre's Galactic azimuth. This bubble motion is a model assumption.

In the moving bubble frame, retain the fitted-centre anchor and set

```text
beta = 7.5 * T_dyn / (0.9777922216731285 * R_eff_pc)
q(lookback) = max(0, 1 - lookback/T_dyn)^beta
```

The original geometry contracts by q, reaching zero at T_dyn, with present
effective expansion speed 7.5 km/s. The effective cap-base centre contracts
towards the original fitted centre. For membership, multiply the instantaneous
axes and cylinder half-height by 1.1. For half-ellipsoids and caps, additionally
move the flat opening outward by 0.1 times the instantaneous original cap
height. The expanded ellipsoid supplies the curved boundary; no constant-area
tube is appended. Cylinders do not receive a second axial extension.

Ignoring numerical boundary tolerance, the cap/half-ellipsoid membership volume
is 1.512 times the unbuffered volume; the cylinder factor is 1.331. These are
membership tolerances. They do not change Script 19 radii, ages or requirements,
or the area used for local-SFR supply. Orbit/scale segments are interpolated
linearly and crossings solved analytically, including the exact age endpoint.

## Summary statistics and figure

Sum explosions within each of the 100 IMF realizations before computing the
mean and 16th/50th/84th percentiles. Blue points show the mean; their vertical
intervals show the 16th–84th percentiles. These intervals describe conditional
IMF sampling, not full uncertainties in velocities, ages, extinction, shared
membership, completeness or boundary evolution. Missing candidate population
fits stop the run explicitly rather than being assigned zero yield.

The regional estimates use 922 and 2,000–4,000 solar masses per Myr per square kpc,
scaled by the compromise ellipse area and the bubble age, with 150 solar masses
formed per expected SN. The 3,000 midpoint is a plotting reference; the
2,000–4,000 range is a literature range, not a confidence interval. The sources
are [Quintana et al. (2025)](https://arxiv.org/abs/2503.08286) and
[Kennicutt & Evans (2012)](https://doi.org/10.1146/annurev-astro-081811-125610).

The figure retains distance-ranked bubble IDs, a 1 kpc divider, a linear axis
through 50 SNe and a logarithmic scale above 50. Zero means remain visible.
`FIGURE_SHOW_QUINTANA` in Script 20 controls the Quintana circles and their
connectors to the KE diamonds; it defaults to `False`. Both SFR estimates remain
in the output tables. The KE diamonds, connecting curve and literature-range
bars remain visible.
Different bubbles can share contributions; summed budgets are not a deduplicated
regional total. Excluding incomplete-6D clusters does not imply their true yield
is zero or prove a strict physical lower bound.

## Output fields

The per-bubble summary retains the required and regional-SFR columns and updates
`expected_ccsn_from_clusters`. Its `sne_p16`, `sne_p50` and `sne_p84` are IMF
ensemble quantiles. `expected_ccsn_6d_traceback` equals the total, while
`expected_ccsn_missing_6d_local` is zero by selection. Supply uncertainty
intervals summarize the stochastic IMF ensemble.

Detailed outputs include the valid cluster input table, geometric candidate
pairs, entry/exit intervals, member initial masses, cluster IMF fits, selected
IMF model IDs, per-cluster contributions and per-bubble realization arrays.
`20_run_parameters.json` records input checksums and assumptions. The existing
41-column `Open_superbubbles.csv` keeps its established schema; its SN columns
describe required SNe, not the blue supply estimate.
"""

from __future__ import annotations
import time

import argparse
import json
import hashlib
import math
import os
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from numba import njit
from scipy.spatial import cKDTree


HERE = Path(__file__).resolve().parent
os.chdir(HERE)
os.environ.setdefault("MPLBACKEND", "Agg")

DATA_DIR = Path("..") / "data"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "20_cluster_supernova_count_in_superbubble_volume"

DEFAULT_HUNT_RAW_CSV = DATA_DIR / "star_cluster_data" / "hunt2024_clusters_full.csv"
DEFAULT_POPULATION_DIR = DATA_DIR / "star_cluster_data" / "sn_population"
DEFAULT_MEMBERS_CSV = DEFAULT_POPULATION_DIR / "hunt_members_g12_17.csv"
DEFAULT_PARSEC_GRID = DEFAULT_POPULATION_DIR / "parsec_z0p0158_gaia_logage6_8p4_step0p01.dat"
DEFAULT_FINAL_SB_CSV = FINAL_DIR / "superbubble_final_fit_parameters.csv"
DEFAULT_DUST_STATISTICS_CSV = Path("..") / "results" / "intermediate_output" / "5_dust_parameter_selection" / "5_dust_parameter_statistics.csv"
DEFAULT_DUST_JOINT_P_CSV = Path("..") / "results" / "intermediate_output" / "8_random_density_fluctuation_test" / "8_condition_test_joint_p_values_ALL30.csv"
DEFAULT_STEP19_SUMMARY_CSV = Path("..") / "results" / "intermediate_output" / "19_superbubble_energy_budget" / "19_superbubble_energy_budget_summary.csv"

OUT_FILTER_CSV = OUT_DIR / "20_Hunt_valid_sample_statistics.csv"
OUT_DETAIL_CSV = OUT_DIR / "20_cluster_supernova_contribution_details.csv"
OUT_SUMMARY_CSV = OUT_DIR / "20_cluster_supernova_count_by_superbubble.csv"
OUT_CONFIG_JSON = OUT_DIR / "20_run_parameters.json"
OUT_METHOD_MD = OUT_DIR / "20_method_notes.md"
OUT_OPEN_SUPERBUBBLES_CSV = FINAL_DIR / "Open_superbubbles.csv"
OUT_FIG_PNG = FINAL_FIG_DIR / "20_cluster_sn_estimate_in_superbubble_volume.png"
OUT_CLUSTER_INPUTS_CSV = OUT_DIR / "20_cluster_supply_inputs.csv"

# Adopt a 25 pc half-height for planar cylinders when testing 3D membership.
DEFAULT_FULL_CYLINDER_HALF_HEIGHT_KPC = 0.025
HALF_ELLIPSOID_SB_IDS = frozenset({2, 11, 12, 15, 16, 17, 22, 30, 32})

DEFAULT_MASS_COLUMN = "MassTot"
DEFAULT_MASS_ERROR_COLUMN = "e_MassTot"
DEFAULT_AGE_COLUMN = "age_myr_median"
DEFAULT_REQUIRED_SN_COLUMN = "n_open_median"
DEFAULT_BASE_SN_COLUMN = "n_base_median"

DEFAULT_COMPROMISE_FRAC = 0.5
# Local star-formation-rate surface densities in Msun Myr^-1 kpc^-2.
# Quintana et al. (2025), Table 2: https://arxiv.org/html/2503.08286v1
# Kennicutt & Evans (2012): https://doi.org/10.1146/annurev-astro-081811-125610
DEFAULT_REGIONAL_SFR_QUINTANA = 922.0
DEFAULT_REGIONAL_SFR_KE_RANGE = (2000.0, 4000.0)
DEFAULT_SN_MASS_PER_EVENT_MSUN = 150.0
FIGURE_LINEAR_THRESHOLD_SN = 50.0
FIGURE_LINEAR_SCALE = 1.5
# Display-only switch: False hides Quintana circles and their connectors.
# Quintana expectations are still calculated and retained in the summary CSV.
FIGURE_SHOW_QUINTANA = False

OPEN_SUPERBUBBLES_COLUMNS = [
    "id",
    "x_c_kpc",
    "x_c_err_minus_kpc",
    "x_c_err_plus_kpc",
    "y_c_kpc",
    "y_c_err_minus_kpc",
    "y_c_err_plus_kpc",
    "z_c_kpc",
    "z_c_err_minus_kpc",
    "z_c_err_plus_kpc",
    "a_kpc",
    "a_err_minus_kpc",
    "a_err_plus_kpc",
    "b_kpc",
    "b_err_minus_kpc",
    "b_err_plus_kpc",
    "c_kpc",
    "c_err_minus_kpc",
    "c_err_plus_kpc",
    "pa_deg",
    "pa_err_minus_deg",
    "pa_err_plus_deg",
    "class",
    "type",
    "a_compromise_kpc",
    "b_compromise_kpc",
    "r_eff_compromise_kpc",
    "pa_compromise_deg",
    "fit_std_pc",
    "shell_radius_mean_abs_offset_over_R_eq",
    "shell_ridge_uPeak_curvature_rms",
    "dust_shell_ridge_peak_inner_ratio",
    "dust_inner_mean_mag_kpc",
    "dust_inner_std_mag_kpc",
    "p_dust_joint3",
    "age_dyn_myr",
    "age_myr_err_minus",
    "age_myr_err_plus",
    "n_sn_required",
    "n_sn_required_err_minus",
    "n_sn_required_err_plus",
]

OPEN_TYPE_BY_MARK = {1: "N", 2: "S", 3: "P"}

# Compact evidence classes used in the public catalogue. Quantitative columns in
# Open_superbubbles.csv are regenerated from the full fit table and downstream
# analysis products; these labels preserve the final publication grouping.
OPEN_CLASS_BY_ID = {
    1: "A",
    2: "B",
    3: "B",
    6: "B",
    7: "A",
    8: "B",
    9: "B",
    10: "C",
    11: "A",
    12: "A",
    13: "B",
    15: "A",
    16: "A",
    17: "A",
    18: "B",
    19: "B",
    20: "B",
    21: "B",
    22: "C",
    23: "B",
    24: "C",
    25: "A",
    26: "C",
    27: "B",
    28: "B",
    29: "B",
    30: "A",
    31: "A",
    32: "B",
    37: "B",
}


@dataclass(frozen=True)
class RegionalSfrConfig:
    """Regional comparison assumptions, independent of individual cluster fits.

    Surface rates use Msun Myr^-1 kpc^-2; mass per event uses Msun. With
    timescale_myr=None, each bubble uses its Script 19 age. compromise_frac
    sets the slice between the opening plane (0) and the ellipsoid apex (1).
    """
    quintana_msun_myr_kpc2: float = DEFAULT_REGIONAL_SFR_QUINTANA
    ke_low_msun_myr_kpc2: float = DEFAULT_REGIONAL_SFR_KE_RANGE[0]
    ke_high_msun_myr_kpc2: float = DEFAULT_REGIONAL_SFR_KE_RANGE[1]
    sn_mass_per_event_msun: float = DEFAULT_SN_MASS_PER_EVENT_MSUN
    compromise_frac: float = DEFAULT_COMPROMISE_FRAC
    timescale_myr: float | None = None

    @property
    def ke_mid_msun_myr_kpc2(self) -> float:
        """Plotting reference: the arithmetic midpoint, not a separate measurement."""
        return 0.5 * (self.ke_low_msun_myr_kpc2 + self.ke_high_msun_myr_kpc2)


# ============================================================================
# Three-dimensional bubble geometry
# ============================================================================

_SIDE_TOLERANCE_KPC = 1e-12  # Match Script 20's half/cap membership convention.


def _value(row: Any, name: str) -> Any:
    """Read pandas Series fields before attributes (Series.shape is not a field)."""
    if hasattr(row, "__getitem__"):
        try:
            return row[name]
        except (KeyError, TypeError, IndexError):
            pass
    try:
        return getattr(row, name)
    except AttributeError as exc:
        raise ValueError(f"Bubble geometry is missing {name!r}.") from exc


def _number(row: Any, name: str, *, positive: bool = False) -> float:
    try:
        value = float(_value(row, name))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Bubble geometry {name!r} must be numeric.") from exc
    if not np.isfinite(value) or (positive and value <= 0.0):
        qualifier = "finite and positive" if positive else "finite"
        raise ValueError(f"Bubble geometry {name!r} must be {qualifier}.")
    return value


def _geometry(row: Any) -> tuple[np.ndarray, np.ndarray, float, int, list[tuple[float, float]]] | None:
    """Return metric center, axes, angle, quadric dimension and z halfspaces.

    Coordinates, axes and halfspace bounds are in kpc; the angle is in degrees.
    Each halfspace is sign * local_z <= bound. Listed half-ellipsoid IDs use
    their fitted axes. Other ellipsoids use the opening-plane axes and the
    plane-to-apex height as an effective cap. Marks 1 and 2 select opposite
    sides of the opening. Cylinders use an elliptical cross-section and two
    height bounds. A zero-height cap returns None and contributes no volume.
    """
    shape = str(_value(row, "shape")).strip().lower()
    center = np.array([_number(row, f"center_{axis}_kpc") for axis in "xyz"])
    if shape == "cylinder":
        axes = np.array([_number(row, "a_radius_kpc", positive=True),
                         _number(row, "b_radius_kpc", positive=True)])
        height = _number(row, "volume_c_kpc", positive=True)
        return center, axes, _number(row, "angle_deg"), 2, [(1.0, height), (-1.0, height)]
    if shape != "ellipsoid":
        raise ValueError(f"Unsupported bubble shape: {shape!r}.")

    sb_id = int(_number(row, "id"))
    mark = int(_number(row, "mark"))
    if sb_id in HALF_ELLIPSOID_SB_IDS:
        axes = np.array([_number(row, f"{axis}_radius_kpc", positive=True) for axis in "abc"])
        angle = _number(row, "angle_deg")
    else:
        cap_z = _number(row, "xy_plane_z_kpc")
        full_c = _number(row, "c_radius_kpc", positive=True)
        apex_z = center[2] - full_c if mark == 1 else center[2] + full_c
        cap_c = abs(cap_z - apex_z)
        if cap_c <= 0.0:
            return None
        center[2] = cap_z
        axes = np.array([_number(row, "xy_plane_a_kpc", positive=True),
                         _number(row, "xy_plane_b_kpc", positive=True), cap_c])
        angle = _number(row, "xy_plane_angle_deg")
    sides = [(1.0, _SIDE_TOLERANCE_KPC)] if mark == 1 else (
        [(-1.0, _SIDE_TOLERANCE_KPC)] if mark == 2 else [])
    return center, axes, angle, 3, sides


# ============================================================================
# Measured cluster orbits and circular bubble reference
# ============================================================================

@dataclass(frozen=True)
class OrbitTracebackConfig:
    """Galpy reference scales and nominal orbit integration settings.

    ro/zo are the solar radius/height in kpc, vo is the velocity scale in
    km s^-1, and solar_motion is galpy's [-U, V, W] convention in km s^-1.
    The reusable helper defaults to 0.1 Myr; the production calculation
    explicitly passes STEP_MYR=0.05 Myr below.
    """
    step_myr: float = 0.1
    ro_kpc: float = 8.122
    vo_kms: float = 236.0
    zo_kpc: float = 0.0208
    solar_motion: tuple[float, float, float] = (-11.1, 12.24, 7.25)
    bubble_motion: str = "circular"

    def __post_init__(self) -> None:
        for name in ("step_myr", "ro_kpc", "vo_kms"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if not np.isfinite(self.zo_kpc):
            raise ValueError("zo_kpc must be finite.")
        motion = np.asarray(self.solar_motion, dtype=float)
        if motion.shape != (3,) or not np.isfinite(motion).all():
            raise ValueError("solar_motion must contain three finite [-U, V, W] values.")
        if self.bubble_motion not in {"circular", "lsr"}:
            raise ValueError("bubble_motion must be 'circular' or 'lsr'.")


def _validate_lookback(lookback_myr: np.ndarray) -> np.ndarray:
    lookback = np.asarray(lookback_myr, dtype=np.float64)
    if lookback.ndim != 1 or len(lookback) == 0 or not np.isfinite(lookback).all():
        raise ValueError("lookback_myr must be a nonempty finite one-dimensional array.")
    if lookback[0] != 0.0 or np.any(np.diff(lookback) <= 0.0):
        raise ValueError("lookback_myr must start at zero and increase strictly.")
    return lookback


def heliocentric_to_galcen(positions_kpc: np.ndarray, config: OrbitTracebackConfig) -> np.ndarray:
    """Exact galpy position conversion, including solar-height tilt."""
    from galpy.util import coords

    positions = np.asarray(positions_kpc, dtype=np.float64)
    if positions.shape[-1:] != (3,) or not np.isfinite(positions).all():
        raise ValueError("Positions must have a final dimension of three finite coordinates.")
    if positions.size == 0:
        return positions.copy()
    flat = positions.reshape(-1, 3)
    converted = coords.XYZ_to_galcenrect(
        *flat.T, Xsun=config.ro_kpc, Zsun=config.zo_kpc
    )
    return np.asarray(converted).reshape(positions.shape)


def _integrated_positions(orbit, lookback_myr: np.ndarray) -> np.ndarray:
    """Return (time, orbit, XYZ); never launch galpy's default process pool."""
    import astropy.units as u
    from galpy.potential import MWPotential2014

    times = -lookback_myr * u.Myr
    if len(times) > 1:
        orbit.integrate(
            times, MWPotential2014, method="odeint", numcores=1, progressbar=False
        )
        xyz = [np.asarray(getattr(orbit, axis)(times, use_physical=True)) for axis in "xyz"]
    else:
        xyz = [np.asarray(getattr(orbit, axis)(0.0, use_physical=True))[..., None] for axis in "xyz"]
    # A scalar reference orbit returns (time,), whereas a batch returns (orbit,time).
    return np.stack([np.atleast_2d(values).T for values in xyz], axis=-1)


def integrate_cluster_orbits(
    clusters: pd.DataFrame,
    lookback_myr: np.ndarray,
    config: OrbitTracebackConfig,
) -> np.ndarray:
    """Integrate finite nominal 6D Hunt inputs without any age/spatial selection.

    Positions are defined by GLON/GLAT/dist50 (pc); heliocentric U/V/W are
    calculated from ICRS astrometry as in script 13. A caller must report any
    unavailable 6D inputs instead of imputing radial velocities as zero.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from galpy.orbit import Orbit

    lookback = _validate_lookback(lookback_myr)
    positions = np.empty((len(lookback), len(clusters), 3), dtype=np.float64)
    if len(clusters) == 0:
        return positions
    required = ["GLON", "GLAT", "dist50", "RA_ICRS", "DE_ICRS", "pmRA", "pmDE", "RV"]
    values = clusters[required].to_numpy(dtype=np.float64)
    if not np.isfinite(values).all() or np.any(clusters["dist50"].to_numpy(dtype=float) <= 0.0):
        raise ValueError("Cluster orbit inputs must be finite and distances positive.")

    sky = SkyCoord(
        ra=clusters["RA_ICRS"].to_numpy(dtype=float) * u.deg,
        dec=clusters["DE_ICRS"].to_numpy(dtype=float) * u.deg,
        distance=clusters["dist50"].to_numpy(dtype=float) * u.pc,
        pm_ra_cosdec=clusters["pmRA"].to_numpy(dtype=float) * u.mas / u.yr,
        pm_dec=clusters["pmDE"].to_numpy(dtype=float) * u.mas / u.yr,
        radial_velocity=clusters["RV"].to_numpy(dtype=float) * u.km / u.s,
        frame="icrs",
    ).galactic
    uvw = sky.velocity.d_xyz.to_value(u.km / u.s).T
    initial = np.column_stack(
        [values[:, 0], values[:, 1], values[:, 2] / 1000.0, uvw]
    )
    batch_size = 128
    for start in range(0, len(clusters), batch_size):
        stop = min(start + batch_size, len(clusters))
        orbit = Orbit(
            initial[start:stop], lb=True, uvw=True,
            ro=config.ro_kpc, vo=config.vo_kms, zo=config.zo_kpc,
            solarmotion=list(config.solar_motion),
        )
        positions[:, start:stop, :] = _integrated_positions(orbit, lookback)
        if len(clusters) > batch_size:
            print(f"Traceback orbits: {stop}/{len(clusters)} clusters", flush=True)
    return positions


def bubble_reference_track(
    lookback_myr: np.ndarray,
    center_helio_kpc: np.ndarray,
    config: OrbitTracebackConfig,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the reference center (GC kpc) and shape rotation (radians).

    For 'circular', vcirc is evaluated at the present cylindrical radius in
    the MWPotential2014 midplane, and the center is integrated from its actual
    height with initial vR=vZ=0. Thus its height is not artificially held fixed.
    For 'lsr', the entire center translation is the solar-radius LSR track,
    anchored so its present center exactly equals the supplied fitted center.
    """
    from galpy.orbit import Orbit
    from galpy.potential import MWPotential2014, vcirc

    lookback = _validate_lookback(lookback_myr)
    center = np.asarray(center_helio_kpc, dtype=float)
    if center.shape != (3,) or not np.isfinite(center).all():
        raise ValueError("center_helio_kpc must contain three finite coordinates.")
    center_gc = heliocentric_to_galcen(center, config)
    if config.bubble_motion == "circular":
        radius = np.hypot(center_gc[0], center_gc[1])
        if radius <= 0.0:
            raise ValueError("A circular bubble reference cannot start on the Galactic axis.")
        phi0 = np.arctan2(center_gc[1], center_gc[0])
        orbit = Orbit(
            [radius / config.ro_kpc, 0.0,
             float(vcirc(MWPotential2014, radius / config.ro_kpc)),
             center_gc[2] / config.ro_kpc, 0.0, phi0],
            ro=config.ro_kpc, vo=config.vo_kms, zo=config.zo_kpc,
        )
        track = _integrated_positions(orbit, lookback)[:, 0, :]
        phi = np.unwrap(np.arctan2(track[:, 1], track[:, 0]))
        return track, phi - phi[0]

    orbit = Orbit(
        [1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        ro=config.ro_kpc, vo=config.vo_kms, zo=config.zo_kpc,
    )
    lsr_track = _integrated_positions(orbit, lookback)[:, 0, :]
    return center_gc + lsr_track - lsr_track[0], np.zeros(len(lookback))


def bubble_frame_tracks(
    positions_gc_kpc: np.ndarray,
    lookback_myr: np.ndarray,
    center_helio_kpc: np.ndarray,
    config: OrbitTracebackConfig,
) -> np.ndarray:
    """Map orbit positions to coordinates of the transported present-day shape.

    Shape centers and orientation are transported with the chosen reference.
    Each object's relative motion is retained; moving both the shape and the
    object does not force the relative trajectory to be stationary.
    """
    lookback = _validate_lookback(lookback_myr)
    positions = np.asarray(positions_gc_kpc, dtype=np.float64)
    if positions.ndim != 3 or positions.shape[0] != len(lookback) or positions.shape[2] != 3:
        raise ValueError("positions_gc_kpc must have shape (time, cluster, 3).")
    if not np.isfinite(positions).all():
        raise ValueError("Orbit positions must be finite.")
    center = np.asarray(center_helio_kpc, dtype=float)
    center_track, angles = bubble_reference_track(lookback, center, config)
    if positions.shape[1] == 0:
        return positions.copy()
    relative = positions - center_track[:, None, :]
    cosine, sine = np.cos(angles)[:, None], np.sin(angles)[:, None]
    # Inverse of the volume's Galactic-Z rotation, using unreflected GC X/Y.
    unrotated = relative.copy()
    unrotated[..., 0] = cosine * relative[..., 0] + sine * relative[..., 1]
    unrotated[..., 1] = -sine * relative[..., 0] + cosine * relative[..., 1]
    origin_gc = heliocentric_to_galcen(np.zeros(3), config)
    helio_to_gc_basis = (heliocentric_to_galcen(np.eye(3), config) - origin_gc).T
    return unrotated @ helio_to_gc_basis + center


# ============================================================================
# Shrinking geometry and opening-plane envelope
# ============================================================================

PC_PER_KMS_TO_MYR = 0.9777922216731285


def scale_at_lookback(lookback_myr, age_myr, reff_pc, velocity_kms=7.5, law="powerlaw"):
    """Return the dimensionless radius scale at positive lookback times in Myr.

    For age T and present effective radius R (pc), q=(1-t/T)^beta with
    beta=v*T/(0.977792...*R); this matches the adopted present expansion speed
    v (km s^-1). The default is 7.5 km s^-1. All axes shrink by q, reaching
    zero at birth. The optional rigid law keeps the present size at all times.
    """
    values = np.asarray(lookback_myr, dtype=float)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("Lookback times must be finite and nonnegative.")
    if not all(np.isfinite(x) and x > 0 for x in (age_myr, reff_pc, velocity_kms)):
        raise ValueError("Age, effective radius and expansion velocity must be positive.")
    if law == "rigid":
        return np.ones_like(values)
    if law != "powerlaw":
        raise ValueError("law must be powerlaw or rigid.")
    beta = velocity_kms * age_myr / (PC_PER_KMS_TO_MYR * reff_pc)
    return np.maximum(0.0, 1.0 - values / age_myr) ** beta


def _rotation(angle_deg):
    angle = np.deg2rad(angle_deg)
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _shell_multiplier(fraction):
    if not np.isfinite(fraction) or fraction < 0:
        raise ValueError("Shell fraction must be finite and nonnegative.")
    return 1.0 + fraction


def _expanded_sides(axes, dim, sides, multiplier, opening_fraction):
    """Return buffered halfspace bounds in the unscaled, present-day frame.

    Axis dilation multiplies each existing bound. An additional opening slice
    has thickness opening_fraction times the original vertical cap axis;
    the evolving-track caller multiplies this by q as the bubble shrinks.
    Cylinder ends receive only axis dilation, so their height is not padded twice.
    """
    if not np.isfinite(opening_fraction) or opening_fraction < 0:
        raise ValueError("Opening fraction must be finite and nonnegative.")
    padding = opening_fraction * axes[2] if dim == 3 and sides else 0.0
    return [(sign, multiplier * bound + padding) for sign, bound in sides]


def contains_current(points_kpc, sb_row, shell_fraction=0.1, opening_fraction=0.0):
    """Test heliocentric XYZ positions in kpc against the present buffered volume.

    shell_fraction=0.1 expands the axes by 10 percent; it is an inclusion
    envelope, not a hollow shell. opening_fraction separately extends the cap
    plane while retaining the expanded quadric. This helper defaults to no
    opening padding; the production call supplies OPENING_FRACTION=0.1.
    """
    points = np.asarray(points_kpc, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("Positions must be a finite (N, 3) array.")
    geometry = _geometry(sb_row)
    if geometry is None:
        return np.zeros(len(points), dtype=bool)
    center, axes, angle, dim, sides = geometry
    multiplier = _shell_multiplier(shell_fraction)
    local = (points - center) @ _rotation(angle)
    inside = np.sum((local[:, :dim] / (multiplier * axes)) ** 2, axis=1) <= 1.0
    for sign, bound in _expanded_sides(axes, dim, sides, multiplier, opening_fraction):
        inside &= sign * local[:, 2] <= bound
    return inside


def _quadratic_bounds(a, b, c, magnitude_a, magnitude_b, magnitude_c):
    """Solve a*u^2+b*u+c <= 0 on [0,1], allowing zero/negative a.

    Up to two pieces are retained per segment. Linear halfspace and positive
    radius constraints are applied by the caller before chronological merging.
    """
    count = len(a)
    lower = np.zeros((count, 2))
    upper = np.zeros((count, 2))
    upper[:, 0] = 1.0
    eps = 64.0 * np.finfo(float).eps
    linear = np.abs(a) <= eps * magnitude_a
    constant = linear & (np.abs(b) <= eps * magnitude_b)
    upper[constant & (c > eps * magnitude_c), 0] = 0.0
    sloped = linear & ~constant
    rising = sloped & (b > 0)
    falling = sloped & (b < 0)
    upper[rising, 0] = np.minimum(1.0, -c[rising] / b[rising])
    lower[falling, 0] = np.maximum(0.0, -c[falling] / b[falling])
    curved = np.flatnonzero(~linear)
    if not len(curved):
        return lower, upper
    aa, bb, cc = a[curved], b[curved], c[curved]
    disc = bb * bb - 4 * aa * cc
    tolerance = 64 * np.finfo(float).eps * (bb * bb + np.abs(4 * aa * cc))
    disc = np.where((disc < 0) & (disc >= -tolerance), 0.0, disc)
    no_roots = disc < 0
    upper[curved[no_roots & (aa > 0)], 0] = 0.0
    hit = ~no_roots
    indices = curved[hit]
    if len(indices):
        aa, bb, cc = aa[hit], bb[hit], cc[hit]
        q = -0.5 * (bb + np.copysign(np.sqrt(disc[hit]), bb))
        first = q / aa
        second = np.divide(cc, q, out=-bb / (2 * aa), where=q != 0)
        left, right = np.minimum(first, second), np.maximum(first, second)
        upward = aa > 0
        idx = indices[upward]
        lower[idx, 0] = np.maximum(0.0, left[upward])
        upper[idx, 0] = np.minimum(1.0, right[upward])
        idx = indices[~upward]
        upper[idx, 0] = np.minimum(1.0, left[~upward])
        lower[idx, 1] = np.maximum(0.0, right[~upward])
        upper[idx, 1] = 1.0
    return lower, upper


def inside_intervals_for_evolving_track(
    lookback_myr, positions_kpc, sb_row, max_lookback_myr, reff_pc,
    velocity_kms=7.5, shell_fraction=0.1, law="powerlaw", opening_fraction=0.0,
):
    """Return finite-duration entry/exit intervals in the moving, scaled shape.

    Input positions are in kpc in the transported bubble frame, and returned
    entry/exit bounds are positive lookback times in Myr. Within each sampled
    orbit segment, both position and q are interpolated linearly. Quadratic
    volume inequalities and linear opening-plane bounds give the accepted
    subsegments, which are then merged chronologically.

    The q=0 endpoint has zero volume and cannot itself add any supernovae.
    A final partial orbit segment is interpolated to the exact bubble age.
    No division by q or artificial nonzero birth radius is used.
    """
    times = np.asarray(lookback_myr, dtype=float)
    points = np.asarray(positions_kpc, dtype=float)
    if times.ndim != 1 or not len(times) or points.shape != (len(times), 3):
        raise ValueError("Expected one-dimensional times and (N, 3) positions.")
    if not np.isfinite(times).all() or not np.isfinite(points).all() or np.any(np.diff(times) <= 0):
        raise ValueError("Track must be finite with strictly increasing times.")
    age = float(max_lookback_myr)
    if not np.isfinite(age) or age < 0:
        raise ValueError("Bubble age must be finite and nonnegative.")
    multiplier = _shell_multiplier(shell_fraction)
    geometry = _geometry(sb_row)
    empty = np.empty((0, 2), dtype=float)
    if len(times) < 2 or age == 0 or geometry is None:
        return empty
    begin, end = max(0.0, times[0]), min(age, times[-1])
    if end <= begin:
        return empty
    selected = (times > begin) & (times < end)
    clipped_times = np.concatenate(([begin], times[selected], [end]))
    clipped_points = np.column_stack([
        np.interp(clipped_times, times, points[:, k]) for k in range(3)
    ])
    q = scale_at_lookback(clipped_times, age, reff_pc, velocity_kms, law)
    metric_center, axes, angle, dim, sides = geometry
    anchor = np.array([_number(sb_row, f"center_{axis}_kpc") for axis in "xyz"])
    rotation = _rotation(angle)
    offset = (metric_center - anchor) @ rotation
    local = (clipped_points - anchor) @ rotation - q[:, None] * offset
    scaled = local[:, :dim] / (multiplier * axes)
    # A necessary coordinate-box condition cheaply rejects distant trajectories.
    maximum_scale = np.max(q)
    if np.any(np.min(scaled, axis=0) > maximum_scale) or np.any(np.max(scaled, axis=0) < -maximum_scale):
        return empty
    start = scaled[:-1]
    direction = np.diff(scaled, axis=0)
    q0, dq = q[:-1], np.diff(q)
    dd = np.einsum("ij,ij->i", direction, direction)
    sd = np.einsum("ij,ij->i", start, direction)
    ss = np.einsum("ij,ij->i", start, start)
    a, b, c = dd - dq*dq, 2*(sd - q0*dq), ss - q0*q0
    lower, upper = _quadratic_bounds(
        a, b, c, dd + dq*dq, 2*(np.abs(sd) + np.abs(q0*dq)), ss + q0*q0,
    )
    constraints = [(-q0, -dq)]
    for sign, bound in _expanded_sides(axes, dim, sides, multiplier, opening_fraction):
        margin = sign * local[:, 2] - q * bound
        constraints.append((margin[:-1], np.diff(margin)))
    for offset0, slope in constraints:
        constant = slope == 0.0
        upper[constant & (offset0 > 0), :] = 0.0
        rising, falling = slope > 0, slope < 0
        upper[rising, :] = np.minimum(upper[rising, :], (-offset0[rising] / slope[rising])[:, None])
        lower[falling, :] = np.maximum(lower[falling, :], (-offset0[falling] / slope[falling])[:, None])
    valid = upper > lower
    if not np.any(valid):
        return empty
    duration = np.diff(clipped_times)[:, None]
    entries = (clipped_times[:-1, None] + duration * lower)[valid]
    exits = (clipped_times[:-1, None] + duration * upper)[valid]
    order = np.argsort(entries, kind="stable")
    merged = []
    for entry, exit_time in zip(entries[order], exits[order]):
        tolerance = 64*np.finfo(float).eps * max(1.0, abs(entry), abs(exit_time))
        if exit_time <= entry:
            continue
        if merged and entry <= merged[-1][1] + tolerance:
            merged[-1][1] = max(merged[-1][1], float(exit_time))
        else:
            merged.append([float(entry), float(exit_time)])
    return np.asarray(merged).reshape(-1, 2)


# ============================================================================
# Stochastic Kroupa initial mass functions
# ============================================================================

# Stellar initial masses are in Msun; dN/dM is proportional to M^(-alpha).
# HIGH_COEFFICIENT makes the two branches continuous at 0.5 Msun.
SCHEMA = 1
M_MIN = 0.08
M_BREAK = 0.5
M_MAX = 150.0
ALPHA_LOW = 1.3
ALPHA_HIGH = 2.3
HIGH_COEFFICIENT = M_BREAK ** (ALPHA_HIGH - ALPHA_LOW)
MODEL_COUNT = 49_900


def _integral(lo: float, hi: float, exponent: float) -> float:
    power = exponent + 1.0
    return (hi ** power - lo ** power) / power if power != 0.0 else math.log(hi / lo)


NUMBER_LOW = _integral(M_MIN, M_BREAK, -ALPHA_LOW)
NUMBER_HIGH = HIGH_COEFFICIENT * _integral(M_BREAK, M_MAX, -ALPHA_HIGH)
NUMBER_TOTAL = NUMBER_LOW + NUMBER_HIGH
P_LOW = NUMBER_LOW / NUMBER_TOTAL
MEAN_MASS = (
    _integral(M_MIN, M_BREAK, 1.0 - ALPHA_LOW)
    + HIGH_COEFFICIENT * _integral(M_BREAK, M_MAX, 1.0 - ALPHA_HIGH)
) / NUMBER_TOTAL


def sample_kroupa(rng: np.random.Generator, size: int) -> np.ndarray:
    """Independent initial masses, using the exact inverse number CDF."""
    uniforms = rng.random(size)
    low = uniforms < P_LOW
    masses = np.empty(size, dtype=np.float64)
    q = 1.0 - ALPHA_LOW
    masses[low] = (
        M_MIN ** q + (uniforms[low] / P_LOW) * (M_BREAK ** q - M_MIN ** q)
    ) ** (1.0 / q)
    q = 1.0 - ALPHA_HIGH
    masses[~low] = (
        M_BREAK ** q
        + ((uniforms[~low] - P_LOW) / (1.0 - P_LOW)) * (M_MAX ** q - M_BREAK ** q)
    ) ** (1.0 / q)
    return masses


def _sample_cluster(target: float, rng: np.random.Generator) -> np.ndarray:
    """Stop at the nearer of the totals immediately before/after crossing."""
    parts = []
    accumulated = 0.0
    while accumulated < target:
        size = max(64, int(math.ceil((target - accumulated) / MEAN_MASS * 1.15)) + 32)
        part = sample_kroupa(rng, size)
        cumulative = np.cumsum(part, dtype=np.float64)
        crossing = int(np.searchsorted(cumulative, target - accumulated, side="left"))
        if crossing == size:
            parts.append(part)
            accumulated += float(cumulative[-1])
            continue
        before = accumulated + (float(cumulative[crossing - 1]) if crossing else 0.0)
        after = accumulated + float(cumulative[crossing])
        # In an exact tie retain the crossing star; no final-mass rescaling.
        keep = crossing + int((after - target) <= (target - before))
        if keep:
            parts.append(part[:keep])
        break
    return np.concatenate(parts) if parts else np.empty(0, dtype=np.float64)


@njit(cache=False)
def _range_sums(masses, prefix, offsets, lo, hi):
    output = np.empty(len(offsets) - 1, dtype=np.float64)
    for model in range(len(output)):
        start, end = offsets[model], offsets[model + 1]
        left, right = start, end
        while left < right:
            middle = (left + right) // 2
            if masses[middle] < lo:
                left = middle + 1
            else:
                right = middle
        first = left
        left, right = first, end
        while left < right:
            middle = (left + right) // 2
            if masses[middle] <= hi:
                left = middle + 1
            else:
                right = middle
        last = left
        if last <= first:
            output[model] = 0.0
        else:
            value = prefix[last - 1]
            if first > start:
                value -= prefix[first - 1]
            output[model] = value
    return output


class ImfBankRangeError(ValueError):
    """The observed complete-interval mass lies outside the available library."""

    def __init__(self, observed_mass, lo, hi, available_min, available_max):
        self.observed_mass = float(observed_mass)
        self.lo = float(lo)
        self.hi = float(hi)
        self.available_min = float(available_min)
        self.available_max = float(available_max)
        super().__init__(
            f"Observed complete-interval mass {observed_mass:.9g} Msun is outside "
            f"the bank range [{available_min:.9g}, {available_max:.9g}] Msun "
            f"for stellar initial masses [{lo:.9g}, {hi:.9g}] Msun. "
            "Do not treat the nearest library boundary as an adequate fit."
        )


@dataclass
class ImfBank:
    """Sampled birth populations, stored by model in sorted mass arrays.

    targets are requested total birth masses (Msun); total_mass records the
    realized stop-nearest totals. offsets delimit models in the memory-mapped
    arrays, and prefix_mass enables mass sums without resampling any stars.
    """
    targets: np.ndarray
    total_mass: np.ndarray
    sorted_masses: np.memmap
    offsets: np.ndarray
    prefix_mass: np.memmap
    manifest: dict

    def mass_in_range(self, lo: float, hi: float) -> np.ndarray:
        """Actual sampled mass in closed initial-mass interval [lo, hi]."""
        lo, hi = float(lo), float(hi)
        if not np.isfinite(lo) or np.isnan(hi) or lo < 0.0 or hi <= lo:
            raise ValueError("Require finite 0 <= lo < hi; positive infinity is allowed for hi.")
        return _range_sums(self.sorted_masses, self.prefix_mass, self.offsets, lo, hi)

    def massive_stars(self, model_id: int) -> np.ndarray:
        """Read-only sorted sampled stellar initial masses >= 8 Msun."""
        model_id = int(model_id)
        if not 0 <= model_id < len(self.targets):
            raise IndexError("IMF model index is outside the bank.")
        start, end = self.offsets[model_id:model_id + 2]
        model = self.sorted_masses[start:end]
        return model[np.searchsorted(model, 8.0, side="left"):]

    def fit(self, observed_mass: float, lo: float, hi: float,
            rng: np.random.Generator, repeats: int = 100, candidates: int = 1000):
        """Select one best candidate per repeat; return IDs and interval masses.

        Candidates are sampled uniformly without replacement within a repeat;
        repeats are independent. The objective is the absolute difference from
        observed total stellar mass in the specified complete interval.
        """
        observed_mass = float(observed_mass)
        if not np.isfinite(observed_mass) or observed_mass <= 0.0:
            raise ValueError("observed_mass must be finite and positive; an empty sample is not a zero-mass fit.")
        if int(repeats) != repeats or repeats < 1 or int(candidates) != candidates or not 1 <= candidates <= len(self.targets):
            raise ValueError("Require repeats >= 1 and 1 <= candidates <= library size, both integers.")
        interval_mass = self.mass_in_range(lo, hi)
        minimum, maximum = float(interval_mass.min()), float(interval_mass.max())
        tolerance = 1e-10 * max(1.0, abs(maximum), abs(observed_mass))
        if maximum <= 0.0 or observed_mass < minimum - tolerance or observed_mass > maximum + tolerance:
            raise ImfBankRangeError(observed_mass, lo, hi, minimum, maximum)
        ids = np.empty(int(repeats), dtype=np.int64)
        matched_mass = np.empty(int(repeats), dtype=np.float64)
        for repeat in range(int(repeats)):
            choices = rng.choice(len(self.targets), size=int(candidates), replace=False)
            winner = int(choices[np.argmin(np.abs(interval_mass[choices] - observed_mass))])
            ids[repeat] = winner
            matched_mass[repeat] = interval_mass[winner]
        return ids, matched_mass


def _paths_from_stem(cache_dir: Path, stem: str) -> dict:
    return {name: cache_dir / f"{stem}.{suffix}" for name, suffix in {
        "masses": "masses.f32", "prefix": "prefix.f64", "metadata": "metadata.npz",
        "manifest": "manifest.json"}.items()}


def _paths(cache_dir: Path, seed: int) -> dict:
    return _paths_from_stem(cache_dir, f"kroupa_{MODEL_COUNT}_v{SCHEMA}_seed{seed}")


def _load(paths: dict) -> ImfBank:
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    with np.load(paths["metadata"], allow_pickle=False) as metadata:
        targets = metadata["targets"]
        total_mass = metadata["total_mass"]
        offsets = metadata["offsets"]
    count = int(offsets[-1])
    model_count = int(manifest["model_count"])
    if (model_count <= 0 or targets.shape != (model_count,)
            or total_mass.shape != (model_count,)
            or offsets.shape != (model_count + 1,)
            or offsets[0] != 0 or np.any(np.diff(offsets) < 0)):
        raise ValueError("Corrupt IMF cache metadata.")
    if paths["masses"].stat().st_size != count * 4 or paths["prefix"].stat().st_size != count * 8:
        raise ValueError("IMF cache file lengths do not match model offsets.")
    masses = np.memmap(paths["masses"], dtype="<f4", mode="r", shape=(count,))
    prefix = np.memmap(paths["prefix"], dtype="<f8", mode="r", shape=(count,))
    return ImfBank(targets, total_mass, masses, offsets, prefix, manifest)


def _build_target_bank(cache_dir, targets, seed, paths, library_kind) -> ImfBank:
    """Common sampler/cache writer; published and extended files stay separate."""
    cache_dir = Path(cache_dir).resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    targets = np.asarray(targets, dtype=np.float64)
    model_count = len(targets)
    if paths["manifest"].exists():
        bank = _load(paths)
        if bank.manifest.get("schema") != SCHEMA or bank.manifest.get("seed") != int(seed):
            raise ValueError("IMF bank cache configuration mismatch.")
        if not np.array_equal(bank.targets, targets):
            raise ValueError("IMF cache targets differ from the requested mass grid.")
        print(f"Loaded IMF bank: {model_count} models, {len(bank.sorted_masses):,} stars", flush=True)
        return bank

    started = time.perf_counter()
    total_mass = np.empty(model_count, dtype=np.float64)
    offsets = np.zeros(model_count + 1, dtype=np.int64)
    rng = np.random.default_rng(int(seed))
    mass_hash, prefix_hash = hashlib.sha256(), hashlib.sha256()
    mass_partial = paths["masses"].with_suffix(paths["masses"].suffix + ".building")
    prefix_partial = paths["prefix"].with_suffix(paths["prefix"].suffix + ".building")
    max_storage_mass_error = 0.0
    max_target_error = 0.0

    with mass_partial.open("wb") as mass_file, prefix_partial.open("wb") as prefix_file:
        for model_id, target in enumerate(targets):
            sampled = _sample_cluster(float(target), rng)
            original_total = float(sampled.sum(dtype=np.float64))
            sampled.sort()
            stored = sampled.astype("<f4")
            prefix = np.cumsum(stored, dtype="<f8")
            total_mass[model_id] = float(prefix[-1]) if len(prefix) else 0.0
            offsets[model_id + 1] = offsets[model_id] + len(stored)
            max_storage_mass_error = max(max_storage_mass_error, abs(total_mass[model_id] - original_total))
            max_target_error = max(max_target_error, abs(original_total - target))
            mass_bytes, prefix_bytes = stored.tobytes(), prefix.tobytes()
            mass_file.write(mass_bytes); prefix_file.write(prefix_bytes)
            mass_hash.update(mass_bytes); prefix_hash.update(prefix_bytes)
            if (model_id + 1) % 2500 == 0 or model_id + 1 == model_count:
                print(f"IMF bank: {model_id + 1}/{model_count} models; "
                      f"{offsets[model_id + 1]:,} stars; {time.perf_counter() - started:.1f} s", flush=True)

    if max_storage_mass_error >= 1e-3:
        raise AssertionError(f"Storing float32 masses changed total mass by {max_storage_mass_error} Msun.")
    os.replace(mass_partial, paths["masses"])
    os.replace(prefix_partial, paths["prefix"])
    with paths["metadata"].open("wb") as handle:
        np.savez(handle, targets=targets, total_mass=total_mass, offsets=offsets)
    count = int(offsets[-1])
    manifest = {
        "schema": SCHEMA, "seed": int(seed), "model_count": model_count,
        "library_kind": library_kind,
        "target_formula": ("10.0 + 0.1 * arange(49900); last target 4999.9 Msun"
                           if library_kind == "published_grid"
                           else "Explicit target grid stored in metadata; independent mass-range extension"),
        "targets_sha256": hashlib.sha256(targets.tobytes()).hexdigest(),
        "target_min_msun": float(targets.min()), "target_max_msun": float(targets.max()),
        "imf": {"min": M_MIN, "break": M_BREAK, "max": M_MAX,
                "alpha_low": ALPHA_LOW, "alpha_high": ALPHA_HIGH,
                "high_coefficient": HIGH_COEFFICIENT},
        "sampling": "NumPy default_rng PCG64; inverse-CDF stellar number sampling; stop-nearest, ties retain crossing star; no mass rescaling",
        "fit": "100 independent selections of 1000 distinct models; minimize absolute complete-interval mass residual",
        "stored_stars": count, "masses_dtype": "<f4", "prefix_dtype": "<f8",
        "mass_interval": "Closed [lo,hi]; exact binary search in each sorted model",
        "files": {name: str(path) for name, path in paths.items()},
        "masses_sha256": mass_hash.hexdigest(), "prefix_sha256": prefix_hash.hexdigest(),
        "max_float32_model_mass_error_msun": max_storage_mass_error,
        "max_stop_nearest_target_error_msun": max_target_error,
        "elapsed_seconds": time.perf_counter() - started,
        "numpy_version": np.__version__,
    }
    # Manifest is written last only after numerical range queries are verified.
    bank = ImfBank(targets, total_mass,
                   np.memmap(paths["masses"], dtype="<f4", mode="r", shape=(count,)),
                   offsets, np.memmap(paths["prefix"], dtype="<f8", mode="r", shape=(count,)), manifest)
    full_range = bank.mass_in_range(0.0, np.inf)
    max_query_error = float(np.max(np.abs(full_range - total_mass)))
    if max_query_error > 1e-9:
        raise AssertionError("Full mass-range queries did not recover model total masses.")
    paths["manifest"].write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Built IMF bank in {time.perf_counter() - started:.1f} s; "
          f"{count:,} stars, {count * 12 / 1e9:.3f} GB of mass/prefix data", flush=True)
    return bank


def build_or_load_bank(cache_dir, seed: int = 20260912) -> ImfBank:
    """Original 49,900-model grid, retaining its original cache paths and seed."""
    cache_dir = Path(cache_dir).resolve()
    targets = 10.0 + np.arange(MODEL_COUNT, dtype=np.float64) * 0.1
    return _build_target_bank(cache_dir, targets, int(seed),
                              _paths(cache_dir, int(seed)), "published_grid")


def build_or_load_extension(cache_dir, targets=None, seed: int = 20260913) -> ImfBank:
    """Optional independent models above the published library's mass range.

    The default grid is arange(5000, 30001, 10) Msun (2,501 models). Calling
    this function builds the extension; importing the module does not. The
    target grid and seed are hashed into distinct file names. No original
    model, seed, target or cache file is changed or rescaled.
    """
    cache_dir = Path(cache_dir).resolve()
    if targets is None:
        targets = np.arange(5000.0, 30001.0, 10.0)
    targets = np.asarray(targets, dtype=np.float64)
    if (targets.ndim != 1 or targets.size == 0 or not np.isfinite(targets).all()
            or np.any(targets < 5000.0) or np.any(np.diff(targets) <= 0.0)):
        raise ValueError("Extension targets must be finite, strictly increasing, and >= 5000 Msun.")
    digest = hashlib.sha256(targets.tobytes()).hexdigest()
    stem = f"kroupa_extension_n{len(targets)}_v{SCHEMA}_seed{int(seed)}_{digest[:16]}"
    return _build_target_bank(cache_dir, targets, int(seed),
                              _paths_from_stem(cache_dir, stem), "extended_mass_grid")


class CompositeImfBank:
    """A virtual common candidate pool, without copying the large star arrays.

    Global model IDs index concatenated targets/total_mass. massive_stars()
    dispatches to the corresponding component library. Uniform candidate
    selection weights each MODEL equally: the original 0.1-Msun grid therefore
    has 100 times the model density per unit target mass of a 10-Msun extension.
    This sampling choice is explicit and is not a uniform-in-mass prior.
    """

    def __init__(self, banks):
        self.banks = tuple(banks)
        if not self.banks:
            raise ValueError("At least one IMF bank is required.")
        if not all(isinstance(bank, ImfBank) for bank in self.banks):
            raise TypeError("Composite components must be ImfBank objects.")
        self.model_offsets = np.r_[0, np.cumsum([len(bank.targets) for bank in self.banks])].astype(np.int64)
        self.targets = np.concatenate([bank.targets for bank in self.banks])
        self.total_mass = np.concatenate([bank.total_mass for bank in self.banks])
        self.manifest = {
            "library_kind": "composite_candidate_pool",
            "model_count": int(len(self.targets)),
            "component_model_offsets": self.model_offsets.tolist(),
            "candidate_sampling": "Uniform over all component models; no mass-grid-density reweighting",
            "mass_grid_density_caveat": "A 0.1-Msun grid contributes 100 times as many models per target-mass interval as a 10-Msun grid.",
            "components": [bank.manifest for bank in self.banks],
        }

    def mass_in_range(self, lo: float, hi: float) -> np.ndarray:
        return np.concatenate([bank.mass_in_range(lo, hi) for bank in self.banks])

    def model_origin(self, model_id: int) -> dict:
        model_id = int(model_id)
        if not 0 <= model_id < len(self.targets):
            raise IndexError("IMF model index is outside the composite bank.")
        component = int(np.searchsorted(self.model_offsets, model_id, side="right") - 1)
        local_id = model_id - int(self.model_offsets[component])
        return {"component": component, "local_model_id": local_id,
                "library_kind": self.banks[component].manifest.get("library_kind", "published_grid")}

    def massive_stars(self, model_id: int) -> np.ndarray:
        origin = self.model_origin(model_id)
        return self.banks[origin["component"]].massive_stars(origin["local_model_id"])

    def fit(self, observed_mass: float, lo: float, hi: float,
            rng: np.random.Generator, repeats: int = 100, candidates: int = 1000):
        return ImfBank.fit(self, observed_mass, lo, hi, rng,
                           repeats=repeats, candidates=candidates)


def combine_banks(*banks) -> CompositeImfBank:
    """Combine original/extension pools; retain model IDs via component offsets."""
    components = []
    for bank in banks:
        components.extend(bank.banks if isinstance(bank, CompositeImfBank) else [bank])
    return CompositeImfBank(components)


# ============================================================================
# PARSEC photometry, initial masses and lifetimes
# ============================================================================

# Fixed Gaia extinction ratios A_G/A_V and A_RP/A_V for the adopted law.
AG_AV, ARP_AV = 0.83627, 0.63439
class ParsecGrid:
    """Use PARSEC isochrones for photometric initial masses and CCSN lifetimes.

    Mini is initial mass in Msun and logAge is log10(age/yr). The maximum
    living initial mass at each age supplies an explosion-time proxy for
    8-150 Msun progenitors, including evolved stars in the surviving envelope.
    This prescription assigns an event to each sampled star in that domain;
    it does not model mass-dependent failed explosions or remnant outcomes.
    """
    def __init__(self,path):
        path=Path(path)
        with path.open(encoding='utf-8') as f:
            names=next(line[2:].split() for line in f if line.startswith('# Zini'))
        self.table=pd.read_csv(path,sep=r'\s+',comment='#',names=names)
        self.groups={float(a):g.reset_index(drop=True) for a,g in self.table.groupby('logAge')}
        self.logages=np.array(sorted(self.groups))
        # Use the final living isochrone, including evolved massive stars,
        # to construct the lifetime envelope; the grid has no remnant extension.
        self.max_masses=np.array([self.groups[a].Mini.max() for a in self.logages])
        rises=np.diff(self.max_masses)
        if np.any(rises > 0.01):raise ValueError('Non-monotonic surviving-mass envelope: inspect tracks')
        self.max_masses=np.minimum.accumulate(self.max_masses)
        if self.max_masses[0]<150 or self.max_masses[-1]>8:raise ValueError('Grid does not cover 8-150 Msun lifetimes')
        self.last_sn=float(self.lifetime(np.array([8.]))[0])
        self.first_sn=float(self.lifetime(np.array([150.]))[0])
    def lifetime(self,mass):
        """Interpolate the surviving-mass crossing in log mass/log age; return Myr."""
        mass=np.asarray(mass,dtype=float)
        if np.any((mass<8)|(mass>150)):raise ValueError('Progenitor mass outside IMF CCSN domain')
        # Choose the actual crossing segment, including the end of any
        # plateau in maximum surviving mass.
        index=np.searchsorted(-self.max_masses,-mass,side='right')
        index=np.clip(index,1,len(self.max_masses)-1)
        left,right=index-1,index
        fraction=np.log(mass/self.max_masses[left])/np.log(self.max_masses[right]/self.max_masses[left])
        return 10**(self.logages[left]+fraction*(self.logages[right]-self.logages[left])-6)
    def surviving_mass(self,age_myr):
        a=np.asarray(age_myr,dtype=float)
        return np.exp(np.interp(np.log10(np.maximum(a,1e-12)*1e6),self.logages,np.log(self.max_masses),left=np.log(350),right=np.log(self.max_masses[-1])))
    def nearest(self,age_myr):
        index=np.argmin(abs(self.logages-np.log10(age_myr*1e6)))
        return self.groups[self.logages[index]].copy(),float(self.logages[index])
    def mass_from_cmd(self,iso,g_abs,colour):
        """Project onto the densified G-RP versus absolute-G isochrone in magnitudes.

        Both coordinates have equal weight in the Euclidean nearest-neighbour
        distance. Return initial masses in Msun and CMD residuals in magnitudes.
        """
        # Densify each segment, retaining post-main-sequence structure. This
        # avoids treating a sparse tabulation as a set of discrete mass bins.
        xy=np.column_stack([iso.Gmag-iso.G_RPmag,iso.Gmag]); mass=iso.Mini.to_numpy()
        points=[]; masses=[]
        for j in range(len(xy)-1):
            n=max(2,min(300,int(np.linalg.norm(xy[j+1]-xy[j])/0.008)+2))
            w=np.linspace(0,1,n,endpoint=False)
            points.append(xy[j]+w[:,None]*(xy[j+1]-xy[j]))
            masses.append(np.exp(np.log(mass[j])+w*np.log(mass[j+1]/mass[j])))
        points=np.concatenate(points); masses=np.concatenate(masses)
        distance,index=cKDTree(points).query(np.column_stack([colour,g_abs]))
        return masses[index],distance
    def fit_input(self,cluster,members):
        """Build the observed mass normalization and matching completeness domain.

        Use the nearest age isochrone, cluster AV50 and valid member parallaxes
        (mas). Select apparent 12 <= G <= 17 before dereddening; member parallaxes
        give individual absolute G, while the cluster distance sets the model
        magnitude window. Summed initial mass in the observable PMS/MS interval
        union normalizes the IMF. Returned per-star flags retain this selection
        for inspection; CMD residuals are diagnostics, not an extra rejection cut.
        """
        iso,logage=self.nearest(float(cluster.age_myr));av=float(cluster.AV50)
        dm=5*np.log10(float(cluster.dist_pc))-5
        good=np.isfinite(members.Plx)&(members.Plx>0)&np.isfinite(members.Gmag)&np.isfinite(members.RPmag)
        stars=members.loc[good].copy()
        # Observed apparent G selection precedes extinction and distance shifts.
        stars=stars.loc[stars.Gmag.between(12,17)].copy()
        g_abs=stars.Gmag.to_numpy()+5*np.log10(stars.Plx.to_numpy())-10-AG_AV*av
        colour=(stars.Gmag-stars.RPmag).to_numpy()-(AG_AV-ARP_AV)*av
        values,residual=self.mass_from_cmd(iso,g_abs,colour)
        # The completeness domain is the observable PMS/MS mass-interval union.
        # Evolved stars are projected using the full isochrone but included only
        # when their initial mass belongs to this observable mass interval.
        ms=iso.loc[iso.label<=1].sort_values('Mini').drop_duplicates('Mini')
        mx=ms.Mini.to_numpy(); gx=ms.Gmag.to_numpy()+dm+AG_AV*av
        # Intersect every PMS/MS segment with the apparent-G window. A
        # massive-star hook can make G(Mini) non-monotonic; do not invert it
        # with an unchecked np.interp. Keep the union of observable intervals.
        intervals=[]
        for j in range(len(mx)-1):
            delta=gx[j+1]-gx[j]
            if abs(delta)<1e-12:
                if 12<=gx[j]<=17:intervals.append((mx[j],mx[j+1]))
                continue
            f0,f1=sorted(((12-gx[j])/delta,(17-gx[j])/delta))
            f0=max(0.,f0);f1=min(1.,f1)
            if f1>f0:
                intervals.append(tuple(np.exp(np.log(mx[j])+np.array([f0,f1])*np.log(mx[j+1]/mx[j]))))
        merged=[]
        for low_i,high_i in intervals:
            if merged and low_i<=merged[-1][1]*(1+1e-8):merged[-1][1]=max(high_i,merged[-1][1])
            else:merged.append([low_i,high_i])
        if not merged:raise ValueError('Empty observable initial-mass interval')
        merged=[[max(.08,a),min(150.,b)] for a,b in merged if min(150.,b)>max(.08,a)]
        low,high=merged[0][0],merged[-1][1]
        if high<=low: raise ValueError('Empty observable initial-mass interval')
        included=np.zeros(len(values),dtype=bool)
        for lower,upper in merged:included|=(values>=lower)&(values<=upper)
        if not included.any():raise ValueError('No members in complete mass interval')
        stars['initial_mass_msun']=values;stars['cmd_residual_mag']=residual
        stars['absolute_g_dereddened']=g_abs;stars['g_minus_rp_dereddened']=colour
        stars['used_for_imf']=included
        metadata=dict(isochrone_logage=logage,complete_mass_low=low,complete_mass_high=high,
            complete_mass_components=len(merged),complete_mass_intervals_json=json.dumps(merged),observed_complete_mass=float(values[included].sum()),n_valid_photometry=len(stars),
            n_complete_members=int(included.sum()),median_cmd_residual=float(np.median(residual[included])),
            AV50=av,highest_surviving_initial_mass=float(self.surviving_mass(float(cluster.age_myr))))
        return metadata,stars,iso


def fit_visible_imf(bank,meta,rng,repeats=100,candidates=1000):
    """Match observed and model mass over the same observable interval union.

    Gaps between components contribute no normalization mass. Each of 100
    default repeats chooses the closest mass sum from 1000 randomly selected
    library models, retaining a discrete population and its massive stars.
    """
    intervals=json.loads(meta['complete_mass_intervals_json'])
    class VisibleBank:
        targets=bank.targets
        def mass_in_range(self,lo,hi):
            return np.sum([bank.mass_in_range(a,b) for a,b in intervals],axis=0)
    return ImfBank.fit(VisibleBank(),meta['observed_complete_mass'],meta['complete_mass_low'],
        meta['complete_mass_high'],rng,repeats=repeats,candidates=candidates)


# ============================================================================
# Orbit-gated supernova population calculation
# ============================================================================

# Production settings: fixed seeds make IMF selections reproducible. The orbit
# grid is in Myr, expansion speed in km s^-1, and both buffers are fractions
# of the corresponding present-day axes (then scaled with bubble size).
SEED = 20260912
REPEATS = 100
STEP_MYR = 0.05
VELOCITY_KMS = 7.5
SHELL_FRACTION = 0.1
OPENING_FRACTION = 0.1


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")


def load_or_integrate_traceback(clusters, times, config, cache_dir):
    """Key nominal orbit caches by the actual rows, time grid and galpy version."""
    import galpy

    columns = ["hunt_recno", "GLON", "GLAT", "dist50", "RA_ICRS", "DE_ICRS", "pmRA", "pmDE", "RV"]
    parameters = {k: v for k, v in asdict(config).items() if k != "bubble_motion"}
    configuration = {"schema": 1, "galpy_version": galpy.__version__, **parameters}
    digest = hashlib.sha256(json.dumps(configuration, sort_keys=True).encode())
    digest.update(clusters[columns].to_numpy(dtype=np.float64).tobytes())
    digest.update(np.asarray(times, dtype=np.float64).tobytes())
    fingerprint = digest.hexdigest()
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"20_hunt_nominal_orbits_{fingerprint[:16]}.npz"
    if path.exists():
        with np.load(path, allow_pickle=False) as cache:
            if str(cache["fingerprint"].item()) != fingerprint:
                raise ValueError("Orbit-cache fingerprint mismatch.")
            np.testing.assert_array_equal(cache["hunt_recno"], clusters.hunt_recno.to_numpy())
            np.testing.assert_array_equal(cache["lookback_myr"], times)
            positions = cache["positions_galpy_gc_kpc"]
        print(f"Reused nominal orbits: {len(clusters)} clusters", flush=True)
    else:
        print(f"Integrating nominal orbits: {len(clusters)} clusters", flush=True)
        positions = integrate_cluster_orbits(clusters, times, config)
        np.savez_compressed(path, positions_galpy_gc_kpc=positions, lookback_myr=times,
                            hunt_recno=clusters.hunt_recno.to_numpy(),
                            fingerprint=np.asarray(fingerprint))
    if positions.shape != (len(times), len(clusters), 3) or not np.isfinite(positions).all():
        raise ValueError("Invalid nominal traceback positions.")
    return positions, {**configuration, "sha256": fingerprint,
                       "cache_file": path.name, "integrated_cluster_count": len(clusters)}


def geometry_intervals(clusters, bubbles, ages, raw_ages, grid, out):
    """Intersect bubble membership with the possible PARSEC explosion window.

    First retain clusters whose ages allow any 8-150 Msun explosion within
    the largest bubble age. Integrate only those with complete nominal 6D
    data, then test each bubble's moving, shrinking volume. Cluster age minus
    the longest/shortest progenitor lifetimes supplies the pair-specific
    lookback limits. Current membership is recorded only as a diagnostic:
    a cluster can contribute after leaving the present volume.
    """
    eligible = ((clusters.age_myr > grid.first_sn)
                & (clusters.age_myr - grid.last_sn < ages.bubble_age_myr.max()))
    clusters = clusters.copy()
    clusters["sn_time_window_eligible"] = eligible
    clusters["orbit_integrated"] = eligible & clusters.has_6d_kinematics
    orbit = clusters.loc[clusters.orbit_integrated].reset_index(drop=True)
    config = OrbitTracebackConfig(step_myr=STEP_MYR, bubble_motion="circular")
    times = np.arange(int(np.ceil(ages.bubble_age_myr.max() / STEP_MYR)) + 1) * STEP_MYR
    positions, metadata = load_or_integrate_traceback(orbit, times, config, out / "cache/orbits")
    intervals, pairs = [], []
    age_by_id = ages.set_index("sb_id").bubble_age_myr
    radius_by_id = raw_ages.set_index("sb_id").reff_pc
    for b in bubbles.itertuples(index=False):
        sb = int(b.id)
        age, radius = float(age_by_id.loc[sb]), float(radius_by_id.loc[sb])
        near = np.maximum(0, orbit.age_myr.to_numpy() - grid.last_sn)
        far = np.minimum(age, orbit.age_myr.to_numpy() - grid.first_sn)
        candidates = np.flatnonzero(far > near)
        stop = min(len(times), int(np.searchsorted(times, age, side="left")) + 1)
        bt = times[:stop]
        tracks = bubble_frame_tracks(positions[:stop, candidates, :], bt,
                                     np.array([b.center_x_kpc, b.center_y_kpc, b.center_z_kpc]), config)
        current = contains_current(orbit[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(),
                                   b, SHELL_FRACTION, OPENING_FRACTION)
        count = 0
        for j, k in enumerate(candidates):
            c = orbit.iloc[k]
            iv = inside_intervals_for_evolving_track(
                bt, tracks[:, j, :], b, age, radius, velocity_kms=VELOCITY_KMS,
                shell_fraction=SHELL_FRACTION, opening_fraction=OPENING_FRACTION,
                law="powerlaw")
            if not len(iv):
                continue
            iv[:, 0] = np.maximum(iv[:, 0], near[k])
            iv[:, 1] = np.minimum(iv[:, 1], far[k])
            iv = iv[iv[:, 1] > iv[:, 0]]
            if not len(iv):
                continue
            count += 1
            rid = int(c.row_id)
            pairs.append(dict(sb_id=sb, cluster_row_id=rid, cluster_name=c["name"],
                              has_6d=True, currently_inside=bool(current[k]), bubble_age_myr=age))
            for start, end in iv:
                intervals.append(dict(sb_id=sb, cluster_row_id=rid,
                                      start_myr=float(start), end_myr=float(end)))
        print(f"SB{sb}: {count} orbit candidates", flush=True)
    iv = pd.DataFrame(intervals, columns=["sb_id", "cluster_row_id", "start_myr", "end_myr"])
    pair = pd.DataFrame(pairs, columns=["sb_id", "cluster_row_id", "cluster_name", "has_6d", "currently_inside", "bubble_age_myr"])
    iv.to_csv(out / "20_traceback_sn_intervals.csv", index=False)
    pair.to_csv(out / "20_traceback_candidate_pairs.csv", index=False)
    return clusters, iv, pair, metadata


def compute_counts(masses, cluster_age, bubble_age, intervals, grid):
    """Return historical, bubble-age-window and orbit-gated explosion counts.

    Explosion lookback time is cluster_age minus PARSEC lifetime (all Myr).
    Historical events have nonnegative lookback; the age-window subset also
    lies within the bubble age; the final subset must lie in an accepted
    geometric interval. Bounds are inclusive and the same stars supply all
    three counts, so inside <= window <= historical for every realization.
    """
    lookback = cluster_age - grid.lifetime(masses)
    historical = lookback >= 0
    window = historical & (lookback <= bubble_age)
    inside = np.zeros(len(masses), dtype=bool)
    for start, end in intervals:
        inside |= (lookback >= start) & (lookback <= end)
    inside &= window
    counts = int(historical.sum()), int(window.sum()), int(inside.sum())
    assert 0 <= counts[2] <= counts[1] <= counts[0]
    return counts


def calculate_cluster_supply(clusters, bubbles, ages, raw_ages, members_path, parsec_path, out):
    """Fit candidate populations and sum explosions within each IMF realization.

    Only clusters with nonempty orbit/time intersections need photometric
    fits. The original IMF bank is tried first; an out-of-range normalization
    triggers the larger-mass extension. Catalogue mass is retained as a
    diagnostic and sample-selection field, not used to normalize these fits.

    Bubble totals sum the per-cluster counts at each realization index before
    computing p16, p50 and p84; individual cluster percentiles are not added.
    These intervals describe stochastic IMF fitting with fixed measured inputs.
    Astrometry, ages, extinction, membership and boundary assumptions are not
    varied. Original catalogue membership assignments can overlap.

    Missing population fits produce an explicit error before final catalogue
    and figure replacement, rather than silently turning an unknown yield into zero.
    """
    grid = ParsecGrid(parsec_path)
    clusters, intervals, pairs, orbit_metadata = geometry_intervals(
        clusters, bubbles, ages, raw_ages, grid, out)
    members = pd.read_csv(members_path, dtype={"GaiaDR3": str})
    members["Name"] = members.Name.str.strip()
    groups = dict(tuple(members.groupby("Name")))
    base_bank = build_or_load_bank(out / "cache/imf", seed=SEED)
    bank, extended_bank = base_bank, None
    selected = clusters.loc[clusters.row_id.isin(pairs.cluster_row_id)].sort_values("row_id")
    fits, draws, star_tables, failures = [], {}, [], []
    for number, c in enumerate(selected.itertuples(index=False), 1):
        rid = int(c.row_id)
        row = dict(cluster_row_id=rid, cluster_name=c.name, hunt_mass_msun=c.mass_msun,
                   cluster_age_myr=c.age_myr)
        try:
            data = groups.get(c.name)
            if data is None:
                raise ValueError("No member photometry")
            meta, stars, _ = grid.fit_input(c, data)
            row.update(meta)
            try:
                ids, matched = fit_visible_imf(base_bank, meta, np.random.default_rng(SEED + rid),
                                               repeats=REPEATS, candidates=1000)
                row["imf_library"] = "paper_mass_range"
            except ImfBankRangeError:
                if extended_bank is None:
                    extension = build_or_load_extension(out / "cache/imf",
                        targets=np.arange(5000., 30001., 10.), seed=SEED + 1)
                    extended_bank = combine_banks(base_bank, extension)
                    bank = extended_bank
                ids, matched = fit_visible_imf(extended_bank, meta, np.random.default_rng(SEED + rid),
                                               repeats=REPEATS, candidates=1000)
                row["imf_library"] = "extended_mass_range_5000_30000"
            draws[rid] = np.asarray(ids, dtype=int)
            masses = bank.total_mass[ids]
            historical = [int((grid.lifetime(bank.massive_stars(int(mid))) <= c.age_myr).sum()) for mid in ids]
            row.update(status="ok", birth_mass_mean=float(masses.mean()),
                       birth_mass_p16=float(np.percentile(masses, 16)),
                       birth_mass_p50=float(np.median(masses)), birth_mass_p84=float(np.percentile(masses, 84)),
                       historical_sne_mean=float(np.mean(historical)),
                       median_relative_mass_residual=float(np.median(abs(matched / meta["observed_complete_mass"] - 1))))
            stars["cluster_row_id"] = rid
            star_tables.append(stars)
        except ValueError as exc:
            row.update(status="unavailable", reason=str(exc))
            failures.append(row)
        fits.append(row)
        if number % 100 == 0:
            print(f"Fitted {number}/{len(selected)} cluster populations", flush=True)
    fit = pd.DataFrame(fits)
    fit.to_csv(out / "20_cluster_imf_fits.csv", index=False)
    if failures:
        write_json(out / "20_population_fit_failures.json", {"failures": failures})
        raise ValueError(f"{len(failures)} candidate clusters lack a valid IMF fit; inspect fit diagnostics.")
    if star_tables:
        pd.concat(star_tables, ignore_index=True).to_csv(out / "20_member_initial_masses.csv", index=False)
    np.savez_compressed(out / "20_selected_imf_models.npz", **{str(k): v for k, v in draws.items()})
    by_id = clusters.set_index("row_id")
    interval_groups = {key: g[["start_myr", "end_myr"]].to_numpy()
                       for key, g in intervals.groupby(["sb_id", "cluster_row_id"])}
    # A realization is shared across every bubble tested for the same cluster.
    # Empty bubble ensembles remain zero when no eligible orbit contributes.
    ensemble = {int(sb): np.zeros(REPEATS) for sb in ages.sb_id}
    details = []
    for pair in pairs.itertuples(index=False):
        sb, rid = int(pair.sb_id), int(pair.cluster_row_id)
        c = by_id.loc[rid]
        iv = interval_groups[(sb, rid)]
        counts = np.asarray([compute_counts(bank.massive_stars(int(mid)), float(c.age_myr),
                                             pair.bubble_age_myr, iv, grid) for mid in draws[rid]])
        sample = counts[:, 2]
        ensemble[sb] += sample
        # full_window_mean ignores position; expected_ccsn applies the orbit gate.
        # Duration is the union of allowed lookback intervals, not a star count.
        details.append(dict(sb_id=sb, cluster_row_id=rid, cluster_name=pair.cluster_name,
                            has_6d=True, currently_inside=pair.currently_inside,
                            cluster_age_myr=c.age_myr, bubble_age_myr=pair.bubble_age_myr,
                            cluster_older_than_bubble=bool(c.age_myr > pair.bubble_age_myr),
                            cluster_catalogue_mass_msun=c.mass_msun,
                            full_window_mean=float(counts[:, 1].mean()), expected_ccsn=float(sample.mean()),
                            sne_p16=float(np.percentile(sample, 16)), sne_p50=float(np.median(sample)),
                            sne_p84=float(np.percentile(sample, 84)),
                            sn_active_inside_duration_myr=float(np.diff(iv, axis=1).sum())))
    detail = pd.DataFrame(details, columns=[
        "sb_id", "cluster_row_id", "cluster_name", "has_6d", "currently_inside", "cluster_age_myr",
        "bubble_age_myr", "cluster_older_than_bubble", "cluster_catalogue_mass_msun", "full_window_mean",
        "expected_ccsn", "sne_p16", "sne_p50", "sne_p84", "sn_active_inside_duration_myr"])
    rows = []
    for sb, sample in ensemble.items():
        rows.append(dict(sb_id=sb, expected_ccsn_from_clusters=float(sample.mean()),
                         sne_p16=float(np.percentile(sample, 16)), sne_p50=float(np.median(sample)),
                         sne_p84=float(np.percentile(sample, 84))))
    np.savez_compressed(out / "20_bubble_sne_realizations.npz", **{str(k): v for k, v in ensemble.items()})
    metadata = dict(
        method="member_imf_parsec_six_d_traceback", citation="https://doi.org/10.1038/s41586-024-07496-9",
        selection="Finite position/age/positive catalogue mass; complete nominal astrometry and RV required for supply; no additional velocity-precision cut",
        normalization="Member photometry and complete initial-mass interval union, not catalogue total mass",
        excluded_missing_6d=int((~clusters.has_6d_kinematics).sum()),
        with_6d=int(clusters.has_6d_kinematics.sum()), candidate_clusters=len(selected),
        orbits=orbit_metadata, bubble_motion="circular, vR=vz=0 initially, axes co-rotate",
        powerlaw="q=(1-lookback/T)^beta; beta=7.5*T/(0.9777922216731285*R_eff)",
        shell_fraction=SHELL_FRACTION, opening_fraction=OPENING_FRACTION,
        shell="Axes times 1.1; cap plane displaced outward by 0.1*q*original height; expanded quadric retained; cylinder ends not expanded twice",
        sn_lifetimes="PARSEC maximum surviving initial-mass envelope; 8 to 150 Msun",
        first_sn_myr=grid.first_sn, last_sn_myr=grid.last_sn,
        isochrone="PARSEC v1.2S/COLIBRI Z=.0158 Gaia EDR3; logage step .01; no remnant extension",
        extinction="Constant Cardelli Rv3.1 coefficients AG/AV=.83627, ARP/AV=.63439",
        imf="Kroupa .08-.5-150 Msun; slopes 1.3/2.3; stop-nearest actual stellar draws",
        library="49900 targets 10..4999.9 step .1; only out-of-range fits use additional 2501 targets 5000..30000 step 10",
        random_seed=SEED, cluster_seed="20260912 + catalogue row_id", repeats=REPEATS, candidates_per_fit=1000,
        uncertainty="p16-p84 from summed IMF realizations; excludes astrometric, age, extinction, membership and boundary-model uncertainty",
        shared_members="Original catalogue assignments retained; cluster membership can overlap",
        poor_fit_relative_mass_residual_threshold=.1,
        poor_fit_cluster_count=int((fit.median_relative_mass_residual > .1).sum()) if len(fit) else 0,
        reproduction="Independent implementation of public Methods, not the authors' private code",
    )
    return detail, pd.DataFrame(rows), clusters, metadata


# ============================================================================
# Input tables, regional supply, catalogue and figure
# ============================================================================

def require_columns(df: pd.DataFrame, required: list[str], source: str) -> None:
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"{source}  is missing required columns: {missing}")


def relative_path_text(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(HERE))
    except ValueError:
        return str(path)


def load_final_superbubbles(path: Path) -> pd.DataFrame:
    """Validate fitted geometry in kpc/degrees and set the cylinder's 25 pc half-height."""
    df = pd.read_csv(path).copy()
    required = [
        "id",
        "mark",
        "shape",
        "center_x_kpc",
        "center_y_kpc",
        "center_z_kpc",
        "a_radius_kpc",
        "b_radius_kpc",
        "c_radius_kpc",
        "least_squares_c_kpc",
        "angle_deg",
        "xy_plane_z_kpc",
        "xy_plane_a_kpc",
        "xy_plane_b_kpc",
        "xy_plane_angle_deg",
    ]
    require_columns(df, required, str(path))
    numeric = [column for column in required if column != "shape"]
    for column in numeric:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df["shape"] = df["shape"].astype(str).str.lower()
    df = df.dropna(subset=["id", "mark", "shape", "center_x_kpc", "center_y_kpc", "center_z_kpc"]).copy()
    df["id"] = df["id"].astype(int)
    df["mark"] = df["mark"].astype(int)

    for column in ["a_radius_kpc", "b_radius_kpc", "angle_deg", "xy_plane_z_kpc", "xy_plane_a_kpc", "xy_plane_b_kpc"]:
        df = df[np.isfinite(df[column])].copy()
    df = df[(df["a_radius_kpc"] > 0.0) & (df["b_radius_kpc"] > 0.0)].copy()

    cylinder = df["shape"] == "cylinder"
    ellipsoid = df["shape"] == "ellipsoid"
    if (~(cylinder | ellipsoid)).any():
        bad = sorted(df.loc[~(cylinder | ellipsoid), "shape"].unique())
        raise ValueError('Invalid input or missing required data.')
    if (ellipsoid & ~(df["c_radius_kpc"] > 0.0)).any():
        bad_ids = df.loc[ellipsoid & ~(df["c_radius_kpc"] > 0.0), "id"].tolist()
        raise ValueError('Invalid input or missing required data.')

    df["volume_c_kpc"] = df["c_radius_kpc"]
    df.loc[cylinder, "volume_c_kpc"] = DEFAULT_FULL_CYLINDER_HALF_HEIGHT_KPC
    df["xy_plane_angle_deg"] = df["xy_plane_angle_deg"].fillna(df["angle_deg"])
    return df.sort_values("id").reset_index(drop=True)


def build_hunt_cluster_table(
    raw_path: Path,
    mass_column: str,
    mass_error_column: str,
) -> tuple[pd.DataFrame, list[tuple[str, int]]]:
    """Retain Hunt source IDs and flag incomplete kinematics before orbit selection.

    Convert log10(age/yr) to Myr and Galactic positions from pc to kpc.
    hunt_recno preserves the source catalogue identifier; row_id is a new
    sequential index assigned after distance sorting and also seeds IMF fits.
    Positive catalogue mass defines the inherited valid sample, not the IMF
    amplitude. No CST, CMDCl50 or velocity-precision threshold is imposed.
    """
    df = pd.read_csv(raw_path).copy()
    required = [
        "AV50",
        "recno",
        "Name",
        "ID",
        "GLON",
        "GLAT",
        "dist50",
        "logAge50",
        mass_column,
        "RA_ICRS", "DE_ICRS", "pmRA", "pmDE", "RV",
    ]
    require_columns(df, required, str(raw_path))

    optional_mass_error = mass_error_column in df.columns
    optional_quality_columns = [column for column in ["CST", "CMDCl50"] if column in df.columns]
    kinematic_columns = ["RA_ICRS", "DE_ICRS", "pmRA", "pmDE", "RV"]
    optional_kinematic_columns = [column for column in ["e_RV", "n_RV"] if column in df.columns]
    numeric = ["recno", "GLON", "GLAT", "dist50", "logAge50", mass_column] + optional_quality_columns + kinematic_columns + optional_kinematic_columns
    if optional_mass_error:
        numeric.append(mass_error_column)
    for column in numeric:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    stages: list[tuple[str, int]] = [("raw Hunt catalogue", len(df))]
    required_numeric = ["recno", "GLON", "GLAT", "dist50", "logAge50", mass_column]
    current = df.loc[np.isfinite(df[required_numeric]).all(axis=1) & (df["dist50"] > 0.0)].copy()
    stages.append(("finite position, age, mass; positive distance", len(current)))

    lon = np.deg2rad(current["GLON"].to_numpy(dtype=float))
    lat = np.deg2rad(current["GLAT"].to_numpy(dtype=float))
    dist_pc = current["dist50"].to_numpy(dtype=float)
    cos_lat = np.cos(lat)
    current["age_myr"] = np.power(10.0, current["logAge50"].to_numpy(dtype=float) - 6.0)
    current["dist_pc"] = dist_pc
    current["x_helio"] = dist_pc * cos_lat * np.cos(lon)
    current["y_helio"] = dist_pc * cos_lat * np.sin(lon)
    current["z_helio"] = dist_pc * np.sin(lat)

    current = current[current[mass_column] > 0.0].copy()
    stages.append((f"{mass_column} > 0", len(current)))
    current["has_6d_kinematics"] = np.isfinite(current[kinematic_columns]).all(axis=1)
    stages.append(("with finite proper motions and radial velocity", int(current["has_6d_kinematics"].sum())))

    name = current["Name"].astype(str).str.strip()
    fallback_name = current["ID"].astype(str).str.strip()
    current["name"] = np.where(name.ne("") & name.ne("nan"), name, fallback_name)
    current["row_id"] = np.arange(1, len(current) + 1, dtype=int)
    current["source_catalog"] = "hunt2024_raw_all_valid"
    current["hunt_recno"] = current["recno"]
    current["hunt_id"] = current["ID"]
    current["mass_msun"] = current[mass_column]
    current["mass_error_msun"] = current[mass_error_column] if optional_mass_error else np.nan
    current["x_kpc"] = current["x_helio"] / 1000.0
    current["y_kpc"] = current["y_helio"] / 1000.0
    current["z_kpc"] = current["z_helio"] / 1000.0
    for column in ["CST", "CMDCl50"]:
        if column not in current.columns:
            current[column] = np.nan

    current["AV50"] = pd.to_numeric(current["AV50"], errors="coerce")
    output_columns = [
        "row_id",
        "name",
        "source_catalog",
        "hunt_recno",
        "hunt_id",
        "age_myr",
        "dist_pc",
        "x_helio",
        "y_helio",
        "z_helio",
        "x_kpc",
        "y_kpc",
        "z_kpc",
        "mass_msun",
        "mass_error_msun",
        "CST",
        "CMDCl50",
        "GLON",
        "GLAT",
        "dist50",
        "logAge50",
        *kinematic_columns,
        *optional_kinematic_columns,
        "has_6d_kinematics", "AV50",
    ]
    output = current.sort_values(["dist_pc", "row_id"], kind="stable")[output_columns].reset_index(drop=True)
    output["row_id"] = np.arange(1, len(output) + 1, dtype=int)
    return output, stages


def load_step19_summary(path: Path, age_column: str, required_sn_column: str, base_sn_column: str) -> pd.DataFrame:
    df = pd.read_csv(path).copy()
    uncertainty_columns = ["age_myr_p16", "age_myr_p84", "n_open_p16", "n_open_p84", "n_base_p16", "n_base_p84"]
    required = ["sb_id", age_column, required_sn_column, base_sn_column] + uncertainty_columns
    require_columns(df, required, str(path))
    for column in required:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df = df.dropna(subset=required).copy()
    df["sb_id"] = df["sb_id"].astype(int)
    return df[["sb_id", age_column, required_sn_column, base_sn_column] + uncertainty_columns].rename(
        columns={
            age_column: "bubble_age_myr",
            required_sn_column: "required_n_open",
            base_sn_column: "required_n_base",
            "age_myr_p16": "bubble_age_p16_myr",
            "age_myr_p84": "bubble_age_p84_myr",
            "n_open_p16": "required_n_open_p16",
            "n_open_p84": "required_n_open_p84",
            "n_base_p16": "required_n_base_p16",
            "n_base_p84": "required_n_base_p84",
        }
    )


def load_dust_statistics(path: Path) -> pd.DataFrame:
    """Load shell and dust metrics measured by script 5."""
    df = pd.read_csv(path, encoding="utf-8-sig").copy()
    required = [
        "id",
        "shell_radius_mean_abs_offset_over_R_eq",
        "shell_ridge_uPeak_curvature_rms",
        "dust_shell_ridge_peak_inner_ratio",
        "dust_inner_mean_mag_kpc",
        "dust_inner_std_mag_kpc",
    ]
    require_columns(df, required, str(path))
    for column in required:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df = df.dropna(subset=["id"]).copy()
    df["id"] = df["id"].astype(int)
    return df[required]


def load_dust_joint_p_values(path: Path) -> pd.DataFrame:
    """Load script-8 joint dust p-values without rerunning the expensive random test."""
    df = pd.read_csv(path, encoding="utf-8-sig").copy()
    required = ["id", "p_dust_joint3"]
    require_columns(df, required, str(path))
    for column in required:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df = df.dropna(subset=["id"]).copy()
    df["id"] = df["id"].astype(int)
    return df[required]


def _open_round(values: Any, decimals: int) -> pd.Series:
    """Round numeric values for the compact science-facing catalogue."""
    return pd.to_numeric(pd.Series(values), errors="coerce").round(decimals)


def build_open_superbubbles_table(
    superbubbles: pd.DataFrame,
    dust_statistics: pd.DataFrame,
    dust_joint_p: pd.DataFrame,
    step19_summary: pd.DataFrame,
    sn_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Assemble the compact catalogue by bubble ID from completed analysis products.

    Geometry errors are inherited from the fit table. Age and required-SN
    errors are central values minus p16 and p84 minus central values from
    Script 19. Cylinder vertical dimensions and circular position angles
    remain blank because those fitted quantities are not applicable. Output
    precision is column-specific; the detailed source tables retain precision.
    """
    required_geometry = [
        "id",
        "mark",
        "shape",
        "circular_xy",
        "center_x_kpc",
        "center_y_kpc",
        "center_z_kpc",
        "a_radius_kpc",
        "b_radius_kpc",
        "c_radius_kpc",
        "angle_deg",
        "fit_std_pc",
        "mcmc_center_x_err_minus_kpc",
        "mcmc_center_x_err_plus_kpc",
        "mcmc_center_y_err_minus_kpc",
        "mcmc_center_y_err_plus_kpc",
        "mcmc_center_z_err_minus_kpc",
        "mcmc_center_z_err_plus_kpc",
        "mcmc_a_err_minus_kpc",
        "mcmc_a_err_plus_kpc",
        "mcmc_b_err_minus_kpc",
        "mcmc_b_err_plus_kpc",
        "mcmc_c_err_minus_kpc",
        "mcmc_c_err_plus_kpc",
        "mcmc_angle_err_minus_deg",
        "mcmc_angle_err_plus_deg",
        "xy_plane_angle_deg",
    ]
    require_columns(superbubbles, required_geometry, "final superbubble geometry table")
    require_columns(
        sn_summary,
        ["sb_id", "regional_a_compromise_pc", "regional_b_compromise_pc"],
        "script-20 SN summary table",
    )

    geom = superbubbles.copy().sort_values("id").reset_index(drop=True)
    geom["id"] = geom["id"].astype(int)
    ids = geom["id"].astype(int)

    dust = dust_statistics.set_index("id")
    pvals = dust_joint_p.set_index("id")
    age_required = step19_summary.set_index("sb_id")
    sn = sn_summary.copy()
    sn["sb_id"] = pd.to_numeric(sn["sb_id"], errors="coerce")
    sn = sn.dropna(subset=["sb_id"]).copy()
    sn["sb_id"] = sn["sb_id"].astype(int)
    sn = sn.set_index("sb_id")

    missing_sources = {
        "script 5 dust statistics": sorted(set(ids) - set(dust.index.astype(int))),
        "script 8 joint p-values": sorted(set(ids) - set(pvals.index.astype(int))),
        "script 19 energy budget": sorted(set(ids) - set(age_required.index.astype(int))),
        "script 20 SN summary": sorted(set(ids) - set(sn.index.astype(int))),
        "compact class map": sorted(set(ids) - set(OPEN_CLASS_BY_ID)),
    }
    missing_text = [f"{name}: {values}" for name, values in missing_sources.items() if values]
    if missing_text:
        raise ValueError("Open_superbubbles.csv cannot be assembled; missing IDs in " + "; ".join(missing_text))

    dust_rows = dust.reindex(ids)
    pval_rows = pvals.reindex(ids)
    age_rows = age_required.reindex(ids)
    sn_rows = sn.reindex(ids)

    is_cylinder = geom["shape"].astype(str).str.lower().eq("cylinder")
    circular_xy = geom["circular_xy"].astype(bool)

    out = pd.DataFrame({"id": ids.to_numpy(dtype=int)})
    out["x_c_kpc"] = _open_round(geom["center_x_kpc"], 2).to_numpy()
    out["x_c_err_minus_kpc"] = _open_round(geom["mcmc_center_x_err_minus_kpc"], 2).to_numpy()
    out["x_c_err_plus_kpc"] = _open_round(geom["mcmc_center_x_err_plus_kpc"], 2).to_numpy()
    out["y_c_kpc"] = _open_round(geom["center_y_kpc"], 2).to_numpy()
    out["y_c_err_minus_kpc"] = _open_round(geom["mcmc_center_y_err_minus_kpc"], 2).to_numpy()
    out["y_c_err_plus_kpc"] = _open_round(geom["mcmc_center_y_err_plus_kpc"], 2).to_numpy()
    out["z_c_kpc"] = _open_round(geom["center_z_kpc"].where(~is_cylinder), 2).to_numpy()
    out["z_c_err_minus_kpc"] = _open_round(geom["mcmc_center_z_err_minus_kpc"].where(~is_cylinder), 2).to_numpy()
    out["z_c_err_plus_kpc"] = _open_round(geom["mcmc_center_z_err_plus_kpc"].where(~is_cylinder), 2).to_numpy()
    out["a_kpc"] = _open_round(geom["a_radius_kpc"], 2).to_numpy()
    out["a_err_minus_kpc"] = _open_round(geom["mcmc_a_err_minus_kpc"], 2).to_numpy()
    out["a_err_plus_kpc"] = _open_round(geom["mcmc_a_err_plus_kpc"], 2).to_numpy()
    out["b_kpc"] = _open_round(geom["b_radius_kpc"], 2).to_numpy()
    out["b_err_minus_kpc"] = _open_round(geom["mcmc_b_err_minus_kpc"], 2).to_numpy()
    out["b_err_plus_kpc"] = _open_round(geom["mcmc_b_err_plus_kpc"], 2).to_numpy()
    out["c_kpc"] = _open_round(geom["c_radius_kpc"].where(~is_cylinder), 2).to_numpy()
    out["c_err_minus_kpc"] = _open_round(geom["mcmc_c_err_minus_kpc"].where(~is_cylinder), 2).to_numpy()
    out["c_err_plus_kpc"] = _open_round(geom["mcmc_c_err_plus_kpc"].where(~is_cylinder), 2).to_numpy()
    out["pa_deg"] = _open_round(geom["angle_deg"].where(~circular_xy), 0).to_numpy()
    out["pa_err_minus_deg"] = _open_round(geom["mcmc_angle_err_minus_deg"].where(~circular_xy), 0).to_numpy()
    out["pa_err_plus_deg"] = _open_round(geom["mcmc_angle_err_plus_deg"].where(~circular_xy), 0).to_numpy()
    out["class"] = ids.map(OPEN_CLASS_BY_ID).to_numpy()
    out["type"] = geom["mark"].astype(int).map(OPEN_TYPE_BY_MARK).to_numpy()
    out["a_compromise_kpc"] = _open_round(sn_rows["regional_a_compromise_pc"].to_numpy(dtype=float) / 1000.0, 3).to_numpy()
    out["b_compromise_kpc"] = _open_round(sn_rows["regional_b_compromise_pc"].to_numpy(dtype=float) / 1000.0, 3).to_numpy()
    out["r_eff_compromise_kpc"] = _open_round(
        np.sqrt(sn_rows["regional_a_compromise_pc"].to_numpy(dtype=float) * sn_rows["regional_b_compromise_pc"].to_numpy(dtype=float)) / 1000.0,
        3,
    ).to_numpy()
    out["pa_compromise_deg"] = _open_round(geom["xy_plane_angle_deg"].where(~circular_xy), 1).to_numpy()
    out["fit_std_pc"] = _open_round(geom["fit_std_pc"], 1).to_numpy()
    out["shell_radius_mean_abs_offset_over_R_eq"] = _open_round(dust_rows["shell_radius_mean_abs_offset_over_R_eq"], 3).to_numpy()
    out["shell_ridge_uPeak_curvature_rms"] = _open_round(dust_rows["shell_ridge_uPeak_curvature_rms"], 3).to_numpy()
    out["dust_shell_ridge_peak_inner_ratio"] = _open_round(dust_rows["dust_shell_ridge_peak_inner_ratio"], 2).to_numpy()
    out["dust_inner_mean_mag_kpc"] = _open_round(dust_rows["dust_inner_mean_mag_kpc"], 3).to_numpy()
    out["dust_inner_std_mag_kpc"] = _open_round(dust_rows["dust_inner_std_mag_kpc"], 3).to_numpy()
    out["p_dust_joint3"] = _open_round(pval_rows["p_dust_joint3"], 3).to_numpy()
    out["age_dyn_myr"] = _open_round(age_rows["bubble_age_myr"], 1).to_numpy()
    out["age_myr_err_minus"] = _open_round(age_rows["bubble_age_myr"] - age_rows["bubble_age_p16_myr"], 1).to_numpy()
    out["age_myr_err_plus"] = _open_round(age_rows["bubble_age_p84_myr"] - age_rows["bubble_age_myr"], 1).to_numpy()
    out["n_sn_required"] = _open_round(age_rows["required_n_open"], 1).to_numpy()
    out["n_sn_required_err_minus"] = _open_round(age_rows["required_n_open"] - age_rows["required_n_open_p16"], 1).to_numpy()
    out["n_sn_required_err_plus"] = _open_round(age_rows["required_n_open_p84"] - age_rows["required_n_open"], 1).to_numpy()

    return out[OPEN_SUPERBUBBLES_COLUMNS]


def write_open_superbubbles_table(
    path: Path,
    superbubbles: pd.DataFrame,
    dust_statistics: pd.DataFrame,
    dust_joint_p: pd.DataFrame,
    step19_summary: pd.DataFrame,
    sn_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Write the compact catalogue produced by the completed analysis chain."""
    compact = build_open_superbubbles_table(
        superbubbles=superbubbles,
        dust_statistics=dust_statistics,
        dust_joint_p=dust_joint_p,
        step19_summary=step19_summary,
        sn_summary=sn_summary,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    compact.to_csv(path, index=False, encoding="utf-8-sig")
    return compact


def compromise_scale_from_row(row: Any, frac: float) -> float:
    """Scale the opening ellipse to a slice a fraction of the way to the apex.

    frac=0.5 uses the halfway height for marks 1/2; the ellipse scale follows
    the fitted ellipsoid cross-section. Other shapes retain their base scale.
    """
    if str(row.shape).strip().lower() != "ellipsoid":
        return 1.0
    mark = int(row.mark)
    if mark not in (1, 2):
        return 1.0
    c_radius = float(row.c_radius_kpc)
    if c_radius <= 0.0:
        return 1.0

    center_z = float(row.center_z_kpc)
    base_z = float(row.xy_plane_z_kpc)
    apex_z = center_z - c_radius if mark == 1 else center_z + c_radius
    mid_z = base_z + float(frac) * (apex_z - base_z)
    base_scale_sq = max(0.0, 1.0 - ((base_z - center_z) / c_radius) ** 2)
    mid_scale_sq = max(0.0, 1.0 - ((mid_z - center_z) / c_radius) ** 2)
    base_scale = math.sqrt(base_scale_sq)
    mid_scale = math.sqrt(mid_scale_sq)
    return mid_scale / base_scale if base_scale > 0.0 else 1.0


def compromise_ellipse_geometry(row: Any, frac: float) -> tuple[float, float, float, float, float]:
    scale = compromise_scale_from_row(row, frac=frac)
    a_pc = float(row.xy_plane_a_kpc) * scale * 1000.0
    b_pc = float(row.xy_plane_b_kpc) * scale * 1000.0
    if a_pc <= 0.0 or b_pc <= 0.0:
        return scale, np.nan, np.nan, np.nan, np.nan
    area_pc2 = math.pi * a_pc * b_pc
    area_kpc2 = area_pc2 / 1_000_000.0
    return scale, a_pc, b_pc, area_pc2, area_kpc2


def regional_expected_supernovae(
    area_kpc2: float,
    timescale_myr: float,
    sigma_sfr_msun_myr_kpc2: float,
    config: RegionalSfrConfig,
) -> tuple[float, float]:
    """Return formed mass (Msun) and expected events from a uniform regional rate.

    Multiply surface SFR (Msun Myr^-1 kpc^-2), area (kpc^2) and duration (Myr),
    then divide by the adopted 150 Msun per event by default. This comparison
    uses no individual stellar lifetimes or cluster-orbit membership gate.
    """
    if not np.isfinite(area_kpc2) or not np.isfinite(timescale_myr) or timescale_myr <= 0.0:
        return np.nan, np.nan
    stellar_mass_msun = sigma_sfr_msun_myr_kpc2 * area_kpc2 * timescale_myr
    expected_sn = stellar_mass_msun / config.sn_mass_per_event_msun
    return stellar_mass_msun, expected_sn


def build_summary(
    detail: pd.DataFrame,
    step19_summary: pd.DataFrame,
    superbubbles: pd.DataFrame,
    regional_sfr: RegionalSfrConfig,
    population_summary: pd.DataFrame,
    raw_ages: pd.DataFrame,
) -> pd.DataFrame:
    """Combine orbit-gated population statistics with regional SFR comparisons.

    expected_ccsn_from_clusters is the ensemble mean; sne_p50 is its median.
    Required-event columns and their p16/p84 values come from Script 19.
    Regional low/mid/high columns reflect the adopted literature SFR range,
    not a probability interval. A missing-6D contribution of zero records
    exclusion from this calculation, not a measured absence of supernovae.
    """
    rows: list[dict[str, Any]] = []
    population = population_summary.set_index("sb_id")
    radii = raw_ages.set_index("sb_id").reff_pc
    required_lookup = step19_summary.set_index("sb_id").to_dict(orient="index")
    for sb in superbubbles.itertuples(index=False):
        sb_id = int(sb.id)
        members = detail[detail["sb_id"] == sb_id]
        contributors = members[members["expected_ccsn"] > 0.0]
        required = required_lookup[sb_id]
        sn_sum = float(population.loc[sb_id, "expected_ccsn_from_clusters"])
        np.testing.assert_allclose(contributors.expected_ccsn.sum(), sn_sum, atol=1e-10)
        required_n_open = float(required.get("required_n_open", np.nan))
        required_n_base = float(required.get("required_n_base", np.nan))
        required_n_open_p16 = float(required.get("required_n_open_p16", np.nan))
        required_n_open_p84 = float(required.get("required_n_open_p84", np.nan))
        required_n_base_p16 = float(required.get("required_n_base_p16", np.nan))
        required_n_base_p84 = float(required.get("required_n_base_p84", np.nan))
        ratio_open = sn_sum / required_n_open if np.isfinite(required_n_open) and required_n_open > 0.0 else np.nan
        ratio_base = sn_sum / required_n_base if np.isfinite(required_n_base) and required_n_base > 0.0 else np.nan
        bubble_age_myr = float(required.get("bubble_age_myr", np.nan))
        regional_timescale_myr = bubble_age_myr if regional_sfr.timescale_myr is None else float(regional_sfr.timescale_myr)
        regional_scale, regional_a_pc, regional_b_pc, regional_area_pc2, regional_area_kpc2 = compromise_ellipse_geometry(
            sb,
            frac=regional_sfr.compromise_frac,
        )
        regional_values: dict[str, float] = {}
        for scenario, sfr in (
            ("quintana", regional_sfr.quintana_msun_myr_kpc2),
            ("ke_low", regional_sfr.ke_low_msun_myr_kpc2),
            ("ke_mid", regional_sfr.ke_mid_msun_myr_kpc2),
            ("ke_high", regional_sfr.ke_high_msun_myr_kpc2),
        ):
            stellar_mass, expected_sn = regional_expected_supernovae(
                area_kpc2=regional_area_kpc2,
                timescale_myr=regional_timescale_myr,
                sigma_sfr_msun_myr_kpc2=sfr,
                config=regional_sfr,
            )
            regional_values[f"regional_sfr_{scenario}_msun_myr_kpc2"] = sfr
            regional_values[f"regional_stellar_mass_{scenario}_msun"] = stellar_mass
            regional_values[f"regional_expected_sn_{scenario}"] = expected_sn
            regional_values[f"regional_sn_to_required_base_ratio_{scenario}"] = (
                expected_sn / required_n_base if np.isfinite(expected_sn) and required_n_base > 0.0 else np.nan
            )
            regional_values[f"regional_sn_to_required_open_ratio_{scenario}"] = (
                expected_sn / required_n_open if np.isfinite(expected_sn) and required_n_open > 0.0 else np.nan
            )
        rows.append({
                "sb_id": sb_id,
                "mark": int(sb.mark),
                "shape": str(sb.shape),
                "bubble_age_myr": bubble_age_myr,
                "heliocentric_distance_kpc": float(np.linalg.norm([
                    sb.center_x_kpc, sb.center_y_kpc, sb.center_z_kpc,
                ])),
                "orbit_candidate_cluster_count": int(len(members)),
                "currently_inside_candidate_count": int(members.currently_inside.sum()),
                "contributing_cluster_count": int(len(contributors)),
                "contributors_older_than_bubble_count": int(contributors["cluster_older_than_bubble"].sum()),
                "expected_ccsn_from_clusters": sn_sum,
                "sne_p16": float(population.loc[sb_id, "sne_p16"]),
                "sne_p50": float(population.loc[sb_id, "sne_p50"]),
                "sne_p84": float(population.loc[sb_id, "sne_p84"]),
                "expected_ccsn_6d_traceback": sn_sum,
                "expected_ccsn_missing_6d_local": 0.0,
                "cluster_selection": "six_d_only",
                "shell_fraction": SHELL_FRACTION,
                "opening_fraction": OPENING_FRACTION,
                "boundary_law": "powerlaw",
                "expansion_reff_pc": float(radii.loc[sb_id]),
                "expansion_beta": VELOCITY_KMS * bubble_age_myr / (PC_PER_KMS_TO_MYR * radii.loc[sb_id]),
                "present_expansion_velocity_kms": VELOCITY_KMS,
                "required_n_base_from_step19": required_n_base,
                "required_n_base_p16_from_step19": required_n_base_p16,
                "required_n_base_p84_from_step19": required_n_base_p84,
                "required_n_open_from_step19": required_n_open,
                "required_n_open_p16_from_step19": required_n_open_p16,
                "required_n_open_p84_from_step19": required_n_open_p84,
                "cluster_sn_to_required_base_ratio": ratio_base,
                "cluster_sn_to_required_open_ratio": ratio_open,
                "regional_compromise_scale": regional_scale,
                "regional_a_compromise_pc": regional_a_pc,
                "regional_b_compromise_pc": regional_b_pc,
                "regional_area_pc2": regional_area_pc2,
                "regional_area_kpc2": regional_area_kpc2,
                "regional_timescale_myr": regional_timescale_myr,
                "regional_timescale_source": "fixed" if regional_sfr.timescale_myr is not None else "bubble_age_median_from_step19",
                "regional_sn_mass_per_event_msun": regional_sfr.sn_mass_per_event_msun,
                **regional_values,
            }
        )
    return pd.DataFrame(rows)

def write_method_note(path: Path, config: dict, regional_sfr: RegionalSfrConfig) -> None:
    """Save the method documented in this script with the actual run settings.

    The module docstring is the single source of method documentation, so a
    normal execution does not depend on a separate explanatory Markdown file.
    """
    text = "# " + (__doc__ or "").strip() + "\n"
    text += "\n## Run parameters\n\n```json\n" + json.dumps(config, indent=2, default=str) + "\n```\n"
    path.write_text(text, encoding="utf-8")


def plot_summary(summary: pd.DataFrame, output_png: Path) -> None:
    """Compare supply and required counts in increasing heliocentric distance.

    Blue markers show IMF ensemble means with p16-p84 spans, purple bars
    inherit Script 19 uncertainties, and green diamonds span the adopted
    Kennicutt-Evans SFR range. Quintana values remain in tables but their
    markers/connectors are hidden by default. The count axis is linear up to
    50 events and logarithmic above; zero-count markers remain visible.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

    if not np.isfinite(summary["heliocentric_distance_kpc"].to_numpy(dtype=float)).all():
        raise ValueError("Figure 20 requires finite heliocentric bubble-center distances.")
    plot_df = summary.sort_values(["heliocentric_distance_kpc", "sb_id"], kind="stable").copy()
    labels = [f"SB{int(value)}" for value in plot_df["sb_id"]]
    x = np.arange(len(plot_df), dtype=float)

    def count_values(column: str) -> np.ndarray:
        values = plot_df[column].to_numpy(dtype=float, copy=True)
        if np.any(values < 0.0):
            raise ValueError(f"Figure 20 requires nonnegative supernova counts: {column}")
        return values

    def count_yerr(median_column: str, low_column: str, high_column: str) -> list[np.ndarray]:
        median = count_values(median_column)
        low = plot_df[low_column].to_numpy(dtype=float)
        high = plot_df[high_column].to_numpy(dtype=float)
        lower_err = median - low
        upper_err = high - median
        lower_err[~np.isfinite(median) | (lower_err < 0.0)] = np.nan
        upper_err[~np.isfinite(median) | (upper_err < 0.0)] = np.nan
        return [lower_err, upper_err]

    fig, ax = plt.subplots(figsize=(11.5, 5.8), dpi=180)
    ax.plot(
        x,
        count_values("expected_ccsn_from_clusters"),
        marker="o",
        linewidth=2.0,
        label="Cluster historical SNe",
        color="#1b6ca8",
    )
    # Blue intervals are p16-p84 across stochastic IMF realizations only.
    ax.vlines(x, count_values("sne_p16"), count_values("sne_p84"),
              color="#1b6ca8", alpha=0.55, linewidth=1.0, zorder=3)
    ax.plot(
        x,
        count_values("required_n_base_from_step19"),
        marker="s",
        linewidth=1.8,
        label="Baseline required SNe",
        color="#d17a22",
    )
    ax.errorbar(
        x,
        count_values("required_n_open_from_step19"),
        yerr=count_yerr(
            "required_n_open_from_step19",
            "required_n_open_p16_from_step19",
            "required_n_open_p84_from_step19",
        ),
        marker="^",
        markersize=5.5,
        linewidth=1.8,
        label="Open-superbubble required SNe",
        color="#8f3f71",
        ecolor="#c9a1ba",
        elinewidth=1.1,
        capsize=2.4,
        zorder=4,
    )
    # Offset the regional estimates slightly so their bars remain distinct from
    # the Script 19 uncertainty bars at the same categorical position.
    x_regional = x + 0.13
    regional_color = "#16a08a"
    ke_mid_sn = count_values("regional_expected_sn_ke_mid")
    if FIGURE_SHOW_QUINTANA:
        quintana_sn = count_values("regional_expected_sn_quintana")
        ax.vlines(
            x_regional,
            quintana_sn,
            ke_mid_sn,
            color=regional_color,
            linewidth=0.85,
            alpha=0.55,
            zorder=3,
        )
        ax.plot(
            x_regional,
            quintana_sn,
            linestyle="none",
            marker="o",
            markersize=4.4,
            markerfacecolor="white",
            markeredgecolor=regional_color,
            markeredgewidth=1.1,
            zorder=5,
        )
    ax.errorbar(
        x_regional,
        ke_mid_sn,
        yerr=count_yerr(
            "regional_expected_sn_ke_mid",
            "regional_expected_sn_ke_low",
            "regional_expected_sn_ke_high",
        ),
        fmt="D",
        linestyle="-",
        linewidth=1.5,
        markersize=4.3,
        label="Local-SFR supplied SNe",
        color=regional_color,
        ecolor=regional_color,
        elinewidth=1.5,
        capsize=3.0,
        capthick=1.2,
        zorder=5,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Number of supernovae", fontsize=16)
    threshold = FIGURE_LINEAR_THRESHOLD_SN
    ax.set_yscale("symlog", linthresh=threshold, linscale=FIGURE_LINEAR_SCALE, base=10)
    upper_columns = [
        "expected_ccsn_from_clusters", "sne_p84", "required_n_base_from_step19",
        "required_n_open_from_step19", "required_n_open_p84_from_step19",
        "regional_expected_sn_ke_high",
    ]
    if FIGURE_SHOW_QUINTANA:
        upper_columns.append("regional_expected_sn_quintana")
    upper_count = float(np.nanmax(plot_df[upper_columns].to_numpy(dtype=float)))
    upper_limit = max(threshold * 1.25, upper_count * 1.4)
    # A small display margin below zero keeps the zero-count markers intact.
    # Tick labels and all measured values remain nonnegative.
    ax.set_ylim(-0.03 * threshold, upper_limit)
    major_ticks = list(np.linspace(0.0, threshold, 6))
    for decade in range(int(np.floor(np.log10(threshold))), int(np.ceil(np.log10(upper_limit))) + 1):
        major_ticks.extend(value for factor in (1.0, 3.0)
                           if threshold < (value := factor * 10.0**decade) <= upper_limit)
    ax.yaxis.set_major_locator(FixedLocator(major_ticks))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f"{value:g}"))
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(axis="y", alpha=0.25, linewidth=0.8)
    ax.axhline(threshold, color="#858585", linestyle="--", linewidth=1.0, zorder=1)
    for label, offset, alignment in [("log", 4, "bottom"), ("linear", -4, "top")]:
        ax.annotate(label, (0.008, threshold), xycoords=("axes fraction", "data"),
                    xytext=(0, offset), textcoords="offset points", ha="left", va=alignment,
                    color="#6b6b6b", fontsize=9)
    distances = plot_df["heliocentric_distance_kpc"].to_numpy(dtype=float)
    nearby_count = int(np.count_nonzero(distances < 1.0))
    if 0 < nearby_count < len(plot_df):
        # This is a distance-ranked categorical axis. Put the threshold between
        # the last nearby bubble and the first distant bubble, not at x=1.
        ax.axvline(nearby_count - 0.5, color="#737373", linestyle="--", linewidth=1.3, zorder=1)
        ax.text(
            (nearby_count - 1) / 2.0, 0.98, r"$d_\odot < 1\,\mathrm{kpc}$",
            transform=ax.get_xaxis_transform(), ha="center", va="top",
            fontsize=12, color="#4a4a4a",
        )
        far_label = r"$d_\odot \geq 1\,\mathrm{kpc}$" if np.any(distances == 1.0) else r"$d_\odot > 1\,\mathrm{kpc}$"
        ax.text(
            (nearby_count + len(plot_df) - 1) / 2.0, 0.98, far_label,
            transform=ax.get_xaxis_transform(), ha="center", va="top",
            fontsize=12, color="#4a4a4a",
        )
    handles, legend_labels = ax.get_legend_handles_labels()
    legend_by_label = dict(zip(legend_labels, handles))
    legend_order = ["Baseline required SNe", "Open-superbubble required SNe",
                    "Cluster historical SNe", "Local-SFR supplied SNe"]
    ax.legend(
        [legend_by_label[label] for label in legend_order], legend_order,
        frameon=False,
        ncol=2,
        loc="lower left",
        bbox_to_anchor=(0.0, 1.02),
        borderaxespad=0.0,
        fontsize=13.2,
    )
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    """Build the command-line interface for Script 20."""
    parser = argparse.ArgumentParser(description='Run Script 20: OSBs YSO SN estimation.')
    parser.add_argument("--hunt-raw-csv", type=Path, default=DEFAULT_HUNT_RAW_CSV)
    parser.add_argument("--members-csv", type=Path, default=DEFAULT_MEMBERS_CSV)
    parser.add_argument("--parsec-grid", type=Path, default=DEFAULT_PARSEC_GRID)
    parser.add_argument("--final-sb-csv", type=Path, default=DEFAULT_FINAL_SB_CSV)
    parser.add_argument("--dust-statistics-csv", type=Path, default=DEFAULT_DUST_STATISTICS_CSV)
    parser.add_argument("--dust-joint-p-csv", type=Path, default=DEFAULT_DUST_JOINT_P_CSV)
    parser.add_argument("--step19-summary-csv", type=Path, default=DEFAULT_STEP19_SUMMARY_CSV)
    parser.add_argument("--mass-column", default=DEFAULT_MASS_COLUMN)
    parser.add_argument("--mass-error-column", default=DEFAULT_MASS_ERROR_COLUMN)
    parser.add_argument("--age-column", default=DEFAULT_AGE_COLUMN)
    parser.add_argument("--required-sn-column", default=DEFAULT_REQUIRED_SN_COLUMN)
    parser.add_argument("--base-sn-column", default=DEFAULT_BASE_SN_COLUMN)
    parser.add_argument("--compromise-frac", type=float, default=DEFAULT_COMPROMISE_FRAC)
    parser.add_argument(
        "--regional-sfr-quintana",
        type=float,
        default=DEFAULT_REGIONAL_SFR_QUINTANA,
        help="Quintana local SFR surface density in Msun Myr^-1 kpc^-2 (default: 922).",
    )
    parser.add_argument(
        "--regional-sfr-ke-range",
        type=float,
        nargs=2,
        default=DEFAULT_REGIONAL_SFR_KE_RANGE,
        metavar=("LOW", "HIGH"),
        help="Kennicutt & Evans SFR range in Msun Myr^-1 kpc^-2 (default: 2000 4000).",
    )
    parser.add_argument(
        "--regional-sn-mass-per-event",
        type=float,
        default=DEFAULT_SN_MASS_PER_EVENT_MSUN,
        help="Formed stellar mass per expected SN in Msun (default: 150).",
    )
    parser.add_argument(
        "--regional-timescale-myr",
        type=float,
        default=None,
        help="Optional fixed regional integration time; default: each Script 19 median bubble age.",
    )
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not 0.0 <= args.compromise_frac <= 1.0:
        raise ValueError('Invalid input or missing required data.')
    regional_rates = [args.regional_sfr_quintana, *args.regional_sfr_ke_range]
    if not all(np.isfinite(rate) and rate > 0.0 for rate in regional_rates):
        raise ValueError("Regional SFR surface densities must be finite and positive.")
    if args.regional_sfr_ke_range[0] > args.regional_sfr_ke_range[1]:
        raise ValueError("--regional-sfr-ke-range requires LOW <= HIGH.")
    if not np.isfinite(args.regional_sn_mass_per_event) or args.regional_sn_mass_per_event <= 0.0:
        raise ValueError("--regional-sn-mass-per-event must be finite and positive.")
    if args.regional_timescale_myr is not None and (
        not np.isfinite(args.regional_timescale_myr) or args.regional_timescale_myr <= 0.0
    ):
        raise ValueError("--regional-timescale-myr must be finite and positive.")


def main() -> None:
    """Run catalogue preparation, orbital gates, population fitting and comparisons.

    Scripts 5, 8 and 19 must already supply dust metrics, joint p-values and
    bubble ages/energy requirements. Detailed diagnostics and reusable caches
    go to this script's intermediate directory; final catalogue and figure
    go to results/. Recorded input hashes and settings identify the run.
    """
    import hashlib
    args = parse_args()
    validate_args(args)
    regional_sfr = RegionalSfrConfig(
        quintana_msun_myr_kpc2=args.regional_sfr_quintana,
        ke_low_msun_myr_kpc2=args.regional_sfr_ke_range[0],
        ke_high_msun_myr_kpc2=args.regional_sfr_ke_range[1],
        sn_mass_per_event_msun=args.regional_sn_mass_per_event,
        compromise_frac=args.compromise_frac, timescale_myr=args.regional_timescale_myr,
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)
    clusters, filter_stages = build_hunt_cluster_table(args.hunt_raw_csv, args.mass_column, args.mass_error_column)
    superbubbles = load_final_superbubbles(args.final_sb_csv)
    dust_statistics = load_dust_statistics(args.dust_statistics_csv)
    dust_joint_p = load_dust_joint_p_values(args.dust_joint_p_csv)
    ages = load_step19_summary(args.step19_summary_csv, args.age_column, args.required_sn_column, args.base_sn_column)
    missing = sorted(set(superbubbles.id.astype(int)) - set(ages.sb_id.astype(int)))
    if missing or not (np.isfinite(ages.bubble_age_myr).all() and (ages.bubble_age_myr > 0).all()):
        raise ValueError(f"Missing or invalid Script 19 ages: {missing}")
    raw_ages = pd.read_csv(args.step19_summary_csv)
    detail, population, clusters, supply_config = calculate_cluster_supply(
        clusters, superbubbles, ages, raw_ages, args.members_csv, args.parsec_grid, OUT_DIR)
    summary = build_summary(detail, ages, superbubbles, regional_sfr, population, raw_ages)
    filter_stages += [
        ("excluded missing complete 6D", int((~clusters.has_6d_kinematics).sum())),
        ("6D clusters with PARSEC time-window overlap, integrated", int(clusters.orbit_integrated.sum())),
        ("unique orbit candidates in at least one bubble", detail.cluster_row_id.nunique()),
        ("unique positive contributors", detail.loc[detail.expected_ccsn > 0, "cluster_row_id"].nunique()),
    ]
    pd.DataFrame(filter_stages, columns=["stage", "count"]).to_csv(OUT_FILTER_CSV, index=False, encoding="utf-8-sig")
    clusters.to_csv(OUT_CLUSTER_INPUTS_CSV, index=False, encoding="utf-8-sig")
    detail.to_csv(OUT_DETAIL_CSV, index=False, encoding="utf-8-sig")
    summary.to_csv(OUT_SUMMARY_CSV, index=False, encoding="utf-8-sig")
    compact = write_open_superbubbles_table(OUT_OPEN_SUPERBUBBLES_CSV, superbubbles,
        dust_statistics, dust_joint_p, ages, summary)
    sources = [args.hunt_raw_csv, args.members_csv, args.parsec_grid, args.final_sb_csv,
               args.step19_summary_csv, args.dust_statistics_csv, args.dust_joint_p_csv]
    config = dict(
        hunt_raw_csv=args.hunt_raw_csv, members_csv=args.members_csv, parsec_grid=args.parsec_grid,
        final_sb_csv=args.final_sb_csv, step19_summary_csv=args.step19_summary_csv,
        age_column=args.age_column, mass_column=args.mass_column,
        required_sn_column=args.required_sn_column, base_sn_column=args.base_sn_column,
        cluster_supply=supply_config,
        regional_sfr={**regional_sfr.__dict__, "ke_mid_msun_myr_kpc2": regional_sfr.ke_mid_msun_myr_kpc2,
                      "ke_interval_meaning": "Literature range, not a confidence interval",
                      "formula": "Sigma_SFR * area_kpc2 * time_Myr / sn_mass_per_event_msun"},
        figure=dict(order_by="heliocentric_distance_kpc, sb_id", distance_divider_kpc=1.,
                    y_scale="symlog", linear_threshold_sn=50., linear_scale=1.5,
                    show_quintana=FIGURE_SHOW_QUINTANA,
                    blue_interval="p16-p84 of 100 IMF realizations"),
        input_sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources})
    OUT_CONFIG_JSON.write_text(json.dumps(config, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    write_method_note(OUT_METHOD_MD, config, regional_sfr)
    if not args.no_plot:
        plot_summary(summary, OUT_FIG_PNG)
    print(f"Complete-6D catalogue clusters: {clusters.has_6d_kinematics.sum()}")
    print(f"Bubbles summarized: {len(summary)}")
    print(f"Per-bubble summary: {OUT_SUMMARY_CSV}")
    print(f"Compact catalogue: {OUT_OPEN_SUPERBUBBLES_CSV} ({len(compact)} rows)")
    if not args.no_plot:
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
