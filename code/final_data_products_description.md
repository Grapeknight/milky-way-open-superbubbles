# Final Data Products Description

This documentation-only note describes the two catalogue-level final data
products written to `../results/` by the reproducible pipeline:

- `superbubble_final_fit_parameters.csv`: the full 30-object open-superbubble
  geometry and provenance table used by the analysis scripts.
- `Open_superbubbles.csv`: a simplified, rounded catalogue table intended for
  manuscript-facing use and convenient reuse. In a reproducible run, script 20
  assembles this compact table after the upstream geometry, dust, random-density,
  energy-budget, and SN-summary products are available.

`superbubble_final_fit_parameters.csv` is the authoritative machine-readable
table for reproducing the pipeline. `Open_superbubbles.csv` is the compact
science-facing summary.

## Coordinate System and Units

- Coordinates are heliocentric Galactic Cartesian coordinates in kpc.
- The convention is:
  - `x = d cos(b) cos(l)`
  - `y = d cos(b) sin(l)`
  - `z = d sin(b)`
- Columns ending in `_kpc` are in kpc; columns ending in `_pc` are in pc; columns
  ending in `_deg` are in degrees.
- MCMC uncertainty columns store one-sided 16th/84th percentile offsets around
  the median. The adopted values used by the scripts remain the non-MCMC final
  geometry columns unless a script explicitly states otherwise.

## Object and Geometry Classes

### `mark` in `superbubble_final_fit_parameters.csv`

| Value | Meaning | Geometry convention |
|---:|---|---|
| `1` | North-open OSB | Open cap points toward negative local/global Z in the cap convention used by the scripts. |
| `2` | South-open OSB | Open cap points toward positive local/global Z in the cap convention used by the scripts. |
| `3` | Disk-penetrating OSB | Usually represented by an elliptical-cylinder model. |

### `shape` in `superbubble_final_fit_parameters.csv`

| Value | Meaning |
|---|---|
| `ellipsoid` | Ellipsoidal target. The table gives the fitted ellipsoid and the dust-plane cap-base parameters in `xy_plane_*`. |
| `cylinder` | Elliptical-cylinder target. `c_radius_kpc` is set to `0` in the full table to mark the cylinder convention; individual scripts define the vertical half-height they need. |

### `type` in `Open_superbubbles.csv`

| Value | Meaning |
|---|---|
| `N` | North-open OSB, corresponding to `mark = 1`. |
| `S` | South-open OSB, corresponding to `mark = 2`. |
| `P` | Disk-penetrating OSB, corresponding to `mark = 3`. |

### `class` in `Open_superbubbles.csv`

`class` is the compact evidence class used for the simplified catalogue. The
full table retains the numerical diagnostics from which the simplified
catalogue was assembled.

## Full Table: `superbubble_final_fit_parameters.csv`

Rows: 30. Columns: 69.

| Column | Description |
|---|---|
| `id` | OSB identifier used throughout the package. |
| `version` | Final-parameter table version label. |
| `mark` | Opening or morphology class, defined above. |
| `shape` | Three-dimensional geometry model, defined above. |
| `xy_model` | Whether the XY cross-section is elliptical or circularized. |
| `circular_xy` | Boolean flag for the circularized-XY final fit. |
| `potential_circular_xy_by_axis_ratio` | Whether the initial axis-ratio rule suggested circularization. |
| `final_parameter_estimator` | Source of the adopted final parameters, such as `least_squares` or `manual_override`. |
| `final_estimator_reason` | Short provenance note for why the final parameter source was adopted. |
| `initial_free_ab_axis_ratio` | Axis ratio from the unconstrained initial/free XY solution. |
| `fit_ab_axis_ratio` | Adopted final XY axis ratio `a/b`; circularized targets have value `1`. |
| `least_squares_a_kpc` | Least-squares semi-axis `a` before any reporting simplification. |
| `least_squares_b_kpc` | Least-squares semi-axis `b` before any reporting simplification. |
| `least_squares_c_kpc` | Least-squares semi-axis `c`; cylinder targets use the cylinder convention described above. |
| `least_squares_angle_deg` | Least-squares XY position angle. |
| `center_x_kpc` | Adopted final center X coordinate. |
| `center_y_kpc` | Adopted final center Y coordinate. |
| `center_z_kpc` | Adopted final center Z coordinate. |
| `a_radius_kpc` | Adopted final local `a` semi-axis. |
| `b_radius_kpc` | Adopted final local `b` semi-axis. |
| `c_radius_kpc` | Adopted final local `c` semi-axis for ellipsoids; `0` marks the cylinder convention. |
| `angle_deg` | Position angle of the local `a` axis, measured counter-clockwise from global `+X`. |
| `fit_std_pc` | Standard deviation of selected first-contact molecular-cloud residuals from the fitted shell. |
| `C_local` | Local consistency statistic; may be blank for manually handled targets. |
| `catalog_cloud_count` | Number of molecular clouds in the source catalogue. |
| `local_bubble_removed_cloud_count` | Number after removing Local Bubble boundary clouds. |
| `candidate_cloud_count` | Number of preselected candidate clouds. |
| `selected_cloud_count` | Number of first-contact clouds used in the fit. |
| `ray_count` | Number of rays used by the first-contact search. |
| `ray_hit_count` | Number of retained ray-hit records. |
| `hit_ray_count` | Number of rays with at least one hit. |
| `catalog_filter` | Molecular-cloud catalogue filtering rule. |
| `xy_search_scale` | Dimensionless XY search-scale control. |
| `fit_center_bounds_scale` | Dimensionless center-bound scale used during fitting. |
| `search_scale` | Dimensionless search-scale control used by the first-contact stage. |
| `radius_tolerance_pc` | First-contact radial tolerance. |
| `ray_cone_aperture_deg` | Ray-cone aperture used to associate clouds with a ray. |
| `max_hits_per_ray` | Maximum number of retained hits per ray. |
| `fit_cloud_seq_list` | Molecular-cloud sequence IDs used by the fit. |
| `fit_clouds_json` | JSON-encoded selected cloud records and local-frame fit diagnostics. |
| `xy_plane_center_x_kpc` | X coordinate of the representative cap-base or cylinder-base ellipse. |
| `xy_plane_center_y_kpc` | Y coordinate of the representative cap-base or cylinder-base ellipse. |
| `xy_plane_z_kpc` | Z height of the representative dust-plane ellipse. |
| `xy_plane_a_kpc` | Semi-axis `a` of the representative dust-plane ellipse. |
| `xy_plane_b_kpc` | Semi-axis `b` of the representative dust-plane ellipse. |
| `xy_plane_angle_deg` | Position angle of the representative dust-plane ellipse. |
| `xy_plane_definition` | Provenance of the representative ellipse, such as `cap_base_ellipse` or `cylinder_base_ellipse`. |
| `mcmc_center_x_median_kpc` | MCMC median for center X. |
| `mcmc_center_x_err_minus_kpc` | Lower one-sided uncertainty for center X. |
| `mcmc_center_x_err_plus_kpc` | Upper one-sided uncertainty for center X. |
| `mcmc_center_y_median_kpc` | MCMC median for center Y. |
| `mcmc_center_y_err_minus_kpc` | Lower one-sided uncertainty for center Y. |
| `mcmc_center_y_err_plus_kpc` | Upper one-sided uncertainty for center Y. |
| `mcmc_center_z_median_kpc` | MCMC median for center Z. |
| `mcmc_center_z_err_minus_kpc` | Lower one-sided uncertainty for center Z. |
| `mcmc_center_z_err_plus_kpc` | Upper one-sided uncertainty for center Z. |
| `mcmc_a_median_kpc` | MCMC median for semi-axis `a`. |
| `mcmc_a_err_minus_kpc` | Lower one-sided uncertainty for semi-axis `a`. |
| `mcmc_a_err_plus_kpc` | Upper one-sided uncertainty for semi-axis `a`. |
| `mcmc_b_median_kpc` | MCMC median for semi-axis `b`. |
| `mcmc_b_err_minus_kpc` | Lower one-sided uncertainty for semi-axis `b`. |
| `mcmc_b_err_plus_kpc` | Upper one-sided uncertainty for semi-axis `b`. |
| `mcmc_c_median_kpc` | MCMC median for semi-axis `c` or the reporting vertical scale. |
| `mcmc_c_err_minus_kpc` | Lower one-sided uncertainty for `c`. |
| `mcmc_c_err_plus_kpc` | Upper one-sided uncertainty for `c`. |
| `mcmc_angle_median_deg` | MCMC median for the position angle. |
| `mcmc_angle_err_minus_deg` | Lower one-sided uncertainty for the position angle. |
| `mcmc_angle_err_plus_deg` | Upper one-sided uncertainty for the position angle. |
| `mcmc_acceptance_fraction_mean` | Mean MCMC acceptance fraction for the target. |

## Simplified Table: `Open_superbubbles.csv`

Rows: 30. Columns: 41.

This table is a compact reporting catalogue. Values are rounded for readability.
Script 20 regenerates it from the current analysis products: geometry and MCMC
uncertainties from script 3, shell and dust metrics from script 5, joint dust
p-values from script 8, dynamical ages and required supernova counts from script
19, and the compact evidence-class labels from the class map encoded in script
20. No static `Open_superbubbles.csv` input is required.
The table omits internal ray/cloud provenance and MCMC bookkeeping.

| Column | Description |
|---|---|
| `id` | OSB identifier. |
| `x_c_kpc` | Reported center X coordinate. |
| `x_c_err_minus_kpc` | Lower uncertainty for center X. |
| `x_c_err_plus_kpc` | Upper uncertainty for center X. |
| `y_c_kpc` | Reported center Y coordinate. |
| `y_c_err_minus_kpc` | Lower uncertainty for center Y. |
| `y_c_err_plus_kpc` | Upper uncertainty for center Y. |
| `z_c_kpc` | Reported center Z coordinate. For disk-penetrating targets, this follows the simplified reporting convention rather than the cylinder marker in `c_radius_kpc`. |
| `z_c_err_minus_kpc` | Lower uncertainty for center Z. |
| `z_c_err_plus_kpc` | Upper uncertainty for center Z. |
| `a_kpc` | Reported semi-axis `a`. |
| `a_err_minus_kpc` | Lower uncertainty for `a`. |
| `a_err_plus_kpc` | Upper uncertainty for `a`. |
| `b_kpc` | Reported semi-axis `b`. |
| `b_err_minus_kpc` | Lower uncertainty for `b`. |
| `b_err_plus_kpc` | Upper uncertainty for `b`. |
| `c_kpc` | Reported vertical semi-axis or reporting vertical scale. |
| `c_err_minus_kpc` | Lower uncertainty for `c`. |
| `c_err_plus_kpc` | Upper uncertainty for `c`. |
| `pa_deg` | Reported position angle. |
| `pa_err_minus_deg` | Lower uncertainty for position angle. |
| `pa_err_plus_deg` | Upper uncertainty for position angle. |
| `class` | Compact evidence class. |
| `type` | Compact morphology type: `N`, `S`, or `P`, as defined above. |
| `a_compromise_kpc` | Projected compromise semi-axis `a` used by reporting/derived calculations. |
| `b_compromise_kpc` | Projected compromise semi-axis `b` used by reporting/derived calculations. |
| `r_eff_compromise_kpc` | Effective projected compromise radius. |
| `pa_compromise_deg` | Position angle used with the compromise projected ellipse. |
| `fit_std_pc` | First-contact shell fit residual scatter. |
| `shell_radius_mean_abs_offset_over_R_eq` | Mean absolute shell-radius offset normalized by the equivalent radius. |
| `shell_ridge_uPeak_curvature_rms` | RMS curvature statistic of the shell-ridge peak in normalized radius. |
| `dust_shell_ridge_peak_inner_ratio` | Ratio of dust-shell ridge peak to inner-cavity dust level. |
| `dust_inner_mean_mag_kpc` | Mean inner-cavity dust density in mag/kpc. |
| `dust_inner_std_mag_kpc` | Standard deviation of inner-cavity dust density in mag/kpc. |
| `p_dust_joint3` | Joint dust-significance p-value used for the compact evidence summary. |
| `age_dyn_myr` | Dynamical age estimate in Myr. |
| `age_myr_err_minus` | Lower uncertainty for the dynamical age. |
| `age_myr_err_plus` | Upper uncertainty for the dynamical age. |
| `n_sn_required` | Estimated number of required supernovae. |
| `n_sn_required_err_minus` | Lower uncertainty for the required supernova count. |
| `n_sn_required_err_plus` | Upper uncertainty for the required supernova count. |

## Geometry Reconstruction Notes

For an ellipsoid, first rotate a point into the local frame:

```python
import numpy as np

def rotate_to_local(points, center, angle_deg):
    q = np.asarray(points, dtype=float) - np.asarray(center, dtype=float)
    th = np.deg2rad(-angle_deg)
    c, s = np.cos(th), np.sin(th)
    rot = np.array([[c, -s, 0.0],
                    [s,  c, 0.0],
                    [0.0, 0.0, 1.0]])
    return q @ rot.T
```

The local `x'` axis corresponds to `a`, the local `y'` axis to `b`, and the
local `z'` axis is parallel to global Z. The full ellipsoid surface is:

```text
u_full = sqrt((x'/a)^2 + (y'/b)^2 + (z'/c)^2) = 1
```

For dust-shell statistics, ellipsoidal OSBs use the self-similar cap geometry
defined by `xy_plane_*`. Typical shell/interior masks are documented in the
individual scripts that compute the corresponding statistics.
