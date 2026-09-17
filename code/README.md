# The Milky Way as a "Phantom Galaxy": A Bubble-Dominated Disk Sculpted by Stellar Feedback and the Origin of the Radcliffe Wave

This repository is the end-to-end code package for the paper *The Milky Way as
a "Phantom Galaxy": A Bubble-Dominated Disk Sculpted by Stellar Feedback and
the Origin of the Radcliffe Wave*. It reproduces the data analysis, validation
analyses, figures, data tables, and interactive visualization used in the study.

AI assistance disclosure: AI tools were used to help organize, document, and
consistency-check this code package. The scientific analysis, data
interpretation, and final package content remain the responsibility of the
authors.

The static runtime overview figure `pipeline_runtime_summary.png` in this
`code/` directory summarizes reference wall-clock runtimes of the numbered
scripts. Script 0b is an external raw-data export estimate.

This package is organized for Code Ocean's standard capsule layout. Scripts and
documentation are in `code/`, input data are in `data/` at the
capsule root (`../data/` from `code/`), and reproducible-run outputs are written
to `../results/` from the `code/` working directory.

Scripts run from `code/` and use direct relative paths in the current file tree:
input files are read from `../data/`, intermediate outputs are written to
`../results/intermediate_output/`, and final outputs are written to
`../results/`. 

## Download the complete publication package

The [GitHub repository](https://github.com/Grapeknight/milky-way-open-superbubbles)
contains the scripts and documentation in `code/`. Input data and saved results
are distributed as assets of [release v1.0.0](https://github.com/Grapeknight/milky-way-open-superbubbles/releases/tag/v1.0.0).
GitHub's automatically generated **Source code** archives contain the code only.

For the complete package, download these three ZIP files from the release:

- `milky-way-open-superbubbles-v1.0.0-code-data.zip`
- `milky-way-open-superbubbles-v1.0.0-results-1.zip`
- `milky-way-open-superbubbles-v1.0.0-results-2.zip`

Extract all three archives into the same parent directory, merging the common
`milky-way-open-superbubbles-v1.0.0/` folder. The resulting capsule contains
`code/`, `data/`, and `results/`; the main instructions remain in `code/README.md`.
The archives are ordinary independent ZIP files and do not require a split-archive utility.

`SHA256SUMS.txt` verifies the downloaded assets, and `file-manifest.csv` records
the size and SHA-256 checksum of every file inside the complete capsule.
On Linux/macOS/WSL, verify downloads with `sha256sum -c SHA256SUMS.txt` (or
`shasum -a 256 -c SHA256SUMS.txt` on macOS). On PowerShell, compare
`Get-FileHash -Algorithm SHA256 <filename>` with the corresponding checksum.

## Input data and sources

The main input catalogues and data products are:

- Wang, T., Yuan, H., et al. 2025, ApJS, 280, 15: all-sky 3D dust map,
  accessed through `dustmaps3d`.
- Wang, T., Yuan, H., et al. 2025, ApJS, 280, 16: molecular-cloud catalogue
  used for shell-tracing clouds.
- Wang et al. 2026, submitted: Galactic dust-bubble catalogue
  (`../data/Bubbles.csv`) from a separate companion work.
  The direct consumers are scripts 4, 10, 12, 16, 22, and 25; rerun all six
  after replacing this catalogue so their figures, tables, and annotations
  remain synchronized.
- Reid et al. 2019, ApJ, 885, 131: spiral-arm loci and HMSFR sample.
- Hunt, E. L., and Reffert, S. 2024, Astron. Astrophys., 686, A42:
  open-cluster catalogue in `../data/star_cluster_data/hunt2024_clusters_full.csv`
  from the script working directory.
- Konietzka, R., et al. 2024, Nature, 628, 62-65: Radcliffe-Wave young-cluster
  sample in `../data/star_cluster_data/Konietzka2023.csv`.
- MWISP / Galactic Plane Survey `12CO`, CfA `12CO`, and HI4PI HI products:
  compact PV caches under `../data/COdata/` and `../data/HIdata/`.

## Program guide and runtime

All numbered scripts are part of the complete analysis workflow. Script 0b is
the raw CO/HI export step and requires external raw survey cubes; the compact
PV outputs from that step are already included for downstream reproduction.
Shorter descriptions are given here; each script contains a more detailed header
at the top of the file.

In the table, paths are direct relative paths from `code/`: generated intermediate products are under `../results/intermediate_output/...`, generated final products are under `../results/...`, and inputs are under `../data/...`. Runtimes are representative recorded wall-clock measurements and depend on the execution environment. Script 0b is an external raw-data export estimate because the raw survey cubes are not packaged.

| Order | Script | Runtime | Main task | Main outputs and manuscript use |
|---:|---|---:|---|---|
| 0a | `0_Data_preparation_1.py` | 20 min 58.49 s | Build raw/smoothed 3D dust cubes, XY mean dust map, and Radcliffe-Wave dust slice from `dustmaps3d` using memory-bounded chunked output. | Reproducible intermediate inputs in `../results/intermediate_output/3d_dust_map_products/` for scripts 5, 8, 11, 12, 23, 24, 25. |
| 0b | `0_Data_preparation_2.py` | > 6 h (external raw-data export) | Export compact MWISP/CfA/HI4PI PV caches from external original FITS cubes. This step cannot be run from the package alone because the several-hundred-GB raw cubes are not included. | `../data/COdata/`, `../data/HIdata/`; these compact caches are included and used by scripts 21 and 22. |
| 1 | `1_OSBs_first_contact_MCs.py` | 3.27 s | Identify first-contact molecular clouds along rays for each open-superbubble seed. | `../results/intermediate_output/1_first_contact_molecular_clouds/`; input to scripts 2 and 3. |
| 2 | `2_OSBs_3D_fitting.py` | 0.91 s | Fit ellipsoidal-cap or elliptical-cylinder shell geometry to first-contact clouds. | `../results/intermediate_output/2_automated_fit_results/`; input to scripts 3 and 4. |
| 3 | `3_OSBs_MCMC_and_parameter_table.py` | 1 min 16.63 s | Estimate MCMC uncertainties and assemble the final 30-object parameter table. | `../results/superbubble_final_fit_parameters.csv`, `../results/figures/3_mcmc_analysis.png`, corner plots; source of the full catalogue table and MCMC diagnostics. |
| 4 | `4_OSBs_plotting.py` | 1 h 02 min 35.22 s | Draw the six-panel identification map for each superbubble. | `../results/SB_figures/SB{N}.png`; used in the Supplementary Information and Extended Data SB panels. |
| 5 | `5_OSBs_parameter_statistics.py` | 9 min 15.09 s | Measure cavity/shell/ridge statistics and evidence-grade parameters. | `../results/figures/5_parameter_histogram_matrix.png`; Extended Data parameter histogram. |
| 6 | `6_OSBs_initial_sensitivity.py` | 37 min 30.93 s | Perturb initial geometry and retest first-contact cloud recovery. | `../results/figures/6_initial_sensitivity_overlap.png`; initial-condition sensitivity figure. |
| 7 | `7_OSBs_recall_analysis.py` | 0.95 s | Analyze independent repeat-identification and observer agreement. | `../results/figures/7_independent_repeat_recall.png`; independent-identification recall figure. |
| 8 | `8_OSBs_random_dust_test.py` | 3 h 12 min 21.40 s | Run 3D IAAFT random-density-fluctuation tests at fixed fitted geometry. | `../results/figures/8_condition_test_joint_p_all30.png`; random-density test figure. |
| 9 | `9_OSBs_shell_amplitude_test.py` | 35.22 s | Compare shell jump amplitudes with extinction-map and stellar-sample uncertainties; automatically detects an installed `dustmaps3d` cache when needed. | `../results/figures/9_shell_amplitude_vs_extinction_error_3sigma.png`; shell-significance figure. |
| 10 | `10_OSBs_bubbles_HMSFR_test.py` | 11 min 09.93 s | Count shell-associated Grade A/B bubbles and HMSFRs; run random geometry controls. | `../results/figures/10_shell_assignment_xy.png`, `10_shell_count_random_test_overview.png`; feedback-tracer association figures. |
| 11 | `11_OSBs_top_view.py` | 21.01 s | Draw face-on open-superbubble overview and RGB layered dust map. | `../results/figures/11_topview_rgb_dust_layers.png`; source panel for main-text Fig. 1. |
| 12 | `12_OSBs_dust_slices.py` | 11.47 s | Draw three vertical dust slices with open-superbubble ellipses, bubbles, and HMSFRs. | `../results/figures/12_three_dust_slices.png`; main-text Fig. 2. |
| 13 | `13_YSO_traceback_inputs.py` | 3.65 s | Build G1/G2 Hunt and Konietzka young-cluster traceback input tables. | `../results/intermediate_output/13_traceback_cluster_input_data/`; input to scripts 14, 15, 16, 23, 24. |
| 14 | `14_G1_YSO_traceback.py` | 3 min 24.56 s | Integrate G1 orbits for 55 Myr and run age-weighted HDBSCAN family recovery. | `../results/intermediate_output/14_g1_traceback_clustering/`; input to script 16 and 17. |
| 15 | `15_G2_YSO_traceback.py` | 1 min 26.16 s | Integrate G2 orbits for 50 Myr and run no-age-weighted HDBSCAN family recovery. | `../results/intermediate_output/15_g2_traceback_clustering/`; input to script 16. |
| 16 | `16_G1G2_YSO_traceback_plots.py` | 11.83 s | Plot G1/G2 traceback slices, size-age curves, residual velocity fields, and the combined family-SB panel, including the HBW35 overlay read from `../data/Bubbles.csv`. | `../results/figures/traceback_cluster_figures/`; Extended Data traceback panel and young-cluster figures. |
| 17 | `17_SB31_RW_boundary_YSO_traceback.py` | 1.96 s | Produce the SB31/Radcliffe-Wave boundary traceback supplement for G1 group 02 and young RW clusters. | `../results/figures/traceback_cluster_figures/17_g1_family2_median_subtracted_traceback_slices.png`; SB31/RW traceback figure. |
| 18 | `18_YSO_OSBs_shell_random_test.py` | 25.27 s | Test young-cluster age windows inside open-superbubble volumes and on dust shells with longitude randomization. | `../results/figures/18_hunt_young_clusters_shell_relative_offset_trend.png`; Extended Data age-stratification figure. |
| 19 | `19_OSBs_energy_budget.py` | 10.28 s | Estimate dynamical ages and required supernova counts by Monte Carlo energy-budget modelling. | `../results/figures/19_superbubble_energy_budget_age_sn.png`; energy-budget figure and input to script 20. |
| 20 | `20_OSBs_YSO_SN_estimation.py` | 1 min 37.11 s (first run, including cache build) | Estimate historical cluster SNe by matching member-photometry constraints to 100 stochastic IMF populations per cluster, then gating sampled progenitors by PARSEC isochrone ages and complete-6D orbits through each bubble's evolving, contraction-calibrated volume; compare with local-SFR supply using 922 and 2000–4000 Msun Myr^-1 kpc^-2. | `../results/Open_superbubbles.csv`, `../results/figures/20_cluster_sn_estimate_in_superbubble_volume.png`; member-level contributions, per-bubble SN supply, selection statistics, and method parameters. |
| 21 | `21_OSBs_gas_PV_diagrams.py` | 27.98 s | Draw CO/HI longitude-velocity comparisons for three shell velocity prescriptions. | `../results/figures/21_PV_*.png`; CO/HI PV figures. |
| 22 | `22_bubbles_gas_PV_diagram.py` | 5.61 s | Overlay low-latitude Grade A/B bubbles on MWISP CO PV diagrams. | `../results/figures/22_gradeAB_low_lat_bubbles_PV_overlay.png`; closed-bubble PV comparison. |
| 23 | `23_RW_cross_sections.py` | 17.46 s | Draw Radcliffe-Wave dust and young-cluster cross sections with open-superbubble shell overlays. | `../results/figures/23_fig3.png`; combined dust and young-cluster panels for main-text Fig. 3. |
| 24 | `24_RW_like_wave.py` | 12.18 s | Draw three additional Radcliffe-Wave-like dust/cluster slices. | `../results/figures/24_rw_like_three_slices.png`; main-text Fig. 4. |
| 25 | `25_interactive_three_volume_viewer.py` | 7.05 s | Build the Three.js 3D dust-volume viewer with embedded scientific data. | `../results/interactive_3d_viewer.html`; local source for the OSB 3D visualization. |

## Figure provenance by document figure number

### Manuscript figures

| Manuscript figure | Package source product | Corresponding script |
|---|---|---|
| Manuscript Fig. 1 | Composite figure. The Milky Way OSB/RGB panel comes from `../results/figures/11_topview_rgb_dust_layers.png`; the JWST/NGC 628 panel is external. | Script 11 plus external JWST-image assembly |
| Manuscript Fig. 2 | `../results/figures/12_three_dust_slices.png` | Script 12 |
| Manuscript Fig. 3 | `../results/figures/23_fig3.png` | Script 23 |
| Manuscript Fig. 4 | `../results/figures/24_rw_like_three_slices.png` | Script 24 |

### Extended Data figures

In the current manuscript PDF, the Extended Data figures are printed as
Extended Data Fig. 5 through Extended Data Fig. 12 because the figure counter
continues after the four main-text figures.

| Extended Data figure | Package source product | Corresponding script |
|---|---|---|
| Extended Data Fig. 5 | `../results/figures/5_parameter_histogram_matrix.png` | Script 5 |
| Extended Data Fig. 6 | `../results/figures/18_hunt_young_clusters_shell_relative_offset_trend.png` | Script 18 |
| Extended Data Fig. 7 | `../results/figures/20_cluster_sn_estimate_in_superbubble_volume.png`; this figure uses the energy-budget summary generated by script 19. | Script 20, using script 19 outputs |
| Extended Data Fig. 8 | `../results/figures/traceback_cluster_figures/G1G2_SBs_YSO.png` | Script 16 |
| Extended Data Fig. 9 | `../results/SB_figures/SB2.png` | Script 4 |
| Extended Data Fig. 10 | `../results/SB_figures/SB15.png`, `../results/SB_figures/SB16.png` | Script 4 |
| Extended Data Fig. 11 | `../results/SB_figures/SB25.png` | Script 4 |
| Extended Data Fig. 12 | `../results/SB_figures/SB31.png` | Script 4 |

The Extended Data catalogue table is based on
`../results/superbubble_final_fit_parameters.csv` and the compact
`../results/Open_superbubbles.csv`, generated by scripts 3 and 20,
respectively.

### Supplementary Information

The Supplementary Information figure provenance is not itemized here. All
Supplementary Information open-superbubble identification figures are generated
by `4_OSBs_plotting.py` and stored in `../results/SB_figures/`.

## Repository contents

```text
peer_reviewed_code/
|-- code/
|   |-- 0_Data_preparation_1.py
|   |-- 0_Data_preparation_2.py
|   |-- 1_OSBs_first_contact_MCs.py ... 25_interactive_three_volume_viewer.py
|   |-- final_data_products_description.md
|   |-- pipeline_runtime_summary.png
|   |-- README.md
|   |-- run
|   `-- run.sh
|-- data/
|   |-- COdata/
|   |-- HIdata/
|   |-- Cantilever/
|   |-- star_cluster_data/
|   |-- Bubbles.csv
|   |-- MCs.csv
|   |-- OSB_initial_pre_identification_parameters.csv
|   |-- publication_plotting_parameters.csv
|   |-- Reid2019_spirals_xy.csv
|   |-- extinction_error.fits
|   `-- third_party_recall_test.csv
`-- results/
    |-- intermediate_output/       # reproducible intermediate products
    |-- figures/                   # generated manuscript and diagnostic figures
    |-- SB_figures/                # generated six-panel OSB figures
    |-- Open_superbubbles.csv
    |-- superbubble_final_fit_parameters.csv
    `-- interactive_3d_viewer.html
```

Approximate package sizes in the current local copy:

| Directory | Approximate size | Role |
|---|---:|---|
| `code/` | 1.5 MB | Scripts, documentation, runtime summary, and Code Ocean run entry point |
| `data/` | 0.64 GB | Input catalogues and compact survey products |
| `results/` | 2.58 GB | Final outputs and scientific intermediate products; Script 20 rebuilds its IMF/orbit caches when run |

## System requirements

- Operating systems: Windows 11, macOS, or Linux.
- Python: 3.11 recommended to match the documented reference environment.
- Disk space: packaged inputs are about 0.64 GB; at least 10 GB free space is
  recommended for the full workflow and generated caches.
- Internet: required for dependency installation, the first `dustmaps3d` cache
  download if needed, and loading the interactive viewer's Three.js modules.
- Memory: at least 16 GB recommended for the full workflow; large raw FITS
  preprocessing in `0_Data_preparation_2.py` can require more, depending on the
  input cubes.

Windows note: `4_OSBs_plotting.py` uses `healpy` for all-sky panels. In the
tested native Windows Python environment, run the workflow through WSL. On
Code Ocean, Linux, or macOS, it can be run directly if `healpy` is installed.

## Tested environment

The current local environment checked for this README used:

| Component | Version or note |
|---|---|
| Python | 3.11.3 |
| NumPy | 1.24.3 |
| pandas | 2.3.1 |
| Matplotlib | 3.7.2 |
| SciPy | 1.10.1 |
| Astropy | 7.1.0 |
| dustmaps3d | 2.2.13 |
| PyArrow | 20.0.0 |
| emcee | 3.1.4 |
| corner | 2.2.2 |
| galpy | 1.11.2 |
| hdbscan | 0.8.40 |
| astropy-healpix | 1.0.3 |
| Pillow | 9.4.0 |
| healpy | Not installed in the native Windows Python; use WSL/Linux for script 4 |

The environment above records the reference versions; other environments can
produce small numerical or rendering differences. The install commands below
list the required Python packages directly.

## Installation

Create and activate a virtual environment from the repository root. For the
complete workflow on Linux, macOS, or WSL:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy pandas scipy matplotlib emcee corner pyarrow astropy galpy numba hdbscan healpy astropy-healpix dustmaps3d platformdirs Pillow tqdm
```

The dependency list includes `numba` for Script 20 and `tqdm` for the optional
raw-survey export progress display. On Windows, create this environment inside
WSL so the Linux Python and `healpy` installation are used by `run.sh`.

First use of `dustmaps3d` may download the Wang et al. 3D dust-map cache. The
processed dust products needed by most scripts are generated by
`0_Data_preparation_1.py` under `../results/intermediate_output/3d_dust_map_products/`
inside the run output area.

## Quick start

On Code Ocean, set `/code/run` as the run file. The run script creates the
output directories and executes the 26 reproducible scripts in order.

For local Linux or WSL use from the capsule root:

```bash
cd code
bash run.sh
```

On Windows, open WSL, change to the repository root, activate the WSL virtual
environment created above, and run the same commands. A native Windows Python
virtual environment is separate from the WSL environment.

The exception is `0_Data_preparation_2.py`: it cannot be executed from this
package alone because it requires the original HI4PI, CfA CO, and MWISP FITS
cubes. Those raw data are several hundred GB and are not included here. This
package includes the compact PV cache products and the code used to export those
PV products from the raw data.

For manual, individual-script execution, first prepare the output directories:

```bash
cd code
mkdir -p ../results ../results/intermediate_output
```

Then run the scripts in order:

```powershell
# Builds the reproducible 3D dust products under ../results/intermediate_output/3d_dust_map_products/.
python 0_Data_preparation_1.py

# Optional from-raw-data step only.
# Requires external raw HI4PI, CfA CO, and Galactic Plane Survey / MWISP data.
# The resulting compact PV caches are already included in ../data/COdata/ and ../data/HIdata/.
# python 0_Data_preparation_2.py --mwisp-fits <path> --cfa-fits <path> --hi4pi-fits <path>

python 1_OSBs_first_contact_MCs.py
python 2_OSBs_3D_fitting.py
python 3_OSBs_MCMC_and_parameter_table.py
python 4_OSBs_plotting.py
python 5_OSBs_parameter_statistics.py
python 6_OSBs_initial_sensitivity.py
python 7_OSBs_recall_analysis.py
python 8_OSBs_random_dust_test.py
python 9_OSBs_shell_amplitude_test.py
python 10_OSBs_bubbles_HMSFR_test.py
python 11_OSBs_top_view.py
python 12_OSBs_dust_slices.py
python 13_YSO_traceback_inputs.py
python 14_G1_YSO_traceback.py
python 15_G2_YSO_traceback.py
python 16_G1G2_YSO_traceback_plots.py
python 17_SB31_RW_boundary_YSO_traceback.py
python 18_YSO_OSBs_shell_random_test.py
python 19_OSBs_energy_budget.py
python 20_OSBs_YSO_SN_estimation.py
python 21_OSBs_gas_PV_diagrams.py
python 22_bubbles_gas_PV_diagram.py
python 23_RW_cross_sections.py
python 24_RW_like_wave.py
python 25_interactive_three_volume_viewer.py
```

Script 4 is not expected to run in the tested native Windows Python
environment. On Windows, run `bash run.sh` through WSL; on Code Ocean, Linux, or
macOS, run `python 4_OSBs_plotting.py` directly if `healpy` is installed.

## Data preparation notes

The two `0_` scripts generate shared data products used by later scripts.

`0_Data_preparation_1.py` builds:

- `../results/intermediate_output/3d_dust_map_products/raw_3d_dust_cube.parquet`
- `../results/intermediate_output/3d_dust_map_products/smoothed_3d_dust_cube.parquet`
- `../results/intermediate_output/3d_dust_map_products/xy_mean_dust_map.parquet`
- `../results/intermediate_output/3d_dust_map_products/radcliffe_wave_dust_slice.csv`

`0_Data_preparation_2.py` is provided to document how the compact CO/HI PV
caches were exported from the original raw survey cubes:

- code-relative `../data/COdata/MWISP_12CO_PV_cache.npz`
- code-relative `../data/COdata/CfA_12CO_PV_cache.npz`
- code-relative `../data/HIdata/HI4PI_HI_PV_cache.npy`
- code-relative `../data/HIdata/HI4PI_HI_PV_axes.npz`
- code-relative `../data/HIdata/HI4PI_HI_PV_cache_summary.json`

These files are physically stored under `../data/COdata/` and
`../data/HIdata/` in the capsule root.

This script cannot be run directly from the distributed package unless the user
also supplies the original HI4PI, CfA CO, and Galactic Plane Survey / MWISP raw
FITS cubes. These raw files are several hundred GB and are not included. The
compact PV caches needed by scripts 21 and 22 are included instead. To rebuild
the caches from raw data, pass the raw data paths explicitly:

```powershell
python 0_Data_preparation_2.py `
  --mwisp-fits <path-to-mosaic_12CO.fits> `
  --cfa-fits <path-to-CfA-FITS-cube> `
  --hi4pi-fits <path-to-HI4PI-CAR.fits>
```

The `MWISP_DATA_PATH` environment variable can also be used for the MWISP FITS
path.

## Local interactive 3D visualization

The local source file is built by:

```bash
cd code
python3 25_interactive_three_volume_viewer.py
```

Main inputs:

- `../results/intermediate_output/3d_dust_map_products/smoothed_3d_dust_cube.parquet`
- `../data/Bubbles.csv`
- `../results/superbubble_final_fit_parameters.csv`
- `../results/Open_superbubbles.csv` (optional simplified display table; the viewer falls back to the full table if this file is absent)

Outputs:

- `../results/interactive_3d_viewer.html`

The HTML embeds the uint8 dust-volume texture, open-superbubble geometry, and
bubble annotations. It loads Three.js 0.160.0 and its add-ons from `unpkg.com`,
so viewing requires network access. An optional `--output-bin <path>` exports a
separate copy of the volume payload; the HTML does not require that file.
In Code Ocean, the HTML is available after a run at
`/results/interactive_3d_viewer.html`.

## Output conventions

- `../data/` contains original catalogues and reusable non-generated caches.
- `../results/intermediate_output/` contains reproducible products, including
  generated 3D dust-map products, per-script CSV/JSON summaries, random trials,
  and diagnostic tables.
- `../results/` contains final manuscript-level products after a run.
- `../results/superbubble_final_fit_parameters.csv` is the
  authoritative 30-object geometry table generated by script 3.
- `../results/Open_superbubbles.csv` is the simplified catalogue table
  generated by script 20 after all upstream ingredients are available. Geometry
  and MCMC uncertainty columns come from script 3, shell and dust metrics from
  script 5, random-density p-values from script 8, dynamical ages and required
  SN counts from script 19, and the published evidence-class labels from the
  compact-catalogue class map encoded in script 20.
- `final_data_products_description.md` in this `code/` directory is a
  documentation file only; it explains the final data products and their
  columns and is not an analysis output.
- `pipeline_runtime_summary.png` in this `code/` directory is a static runtime
  summary figure for planning pipeline execution; it is
  not regenerated by the pipeline.
- `../results/SB_figures/` contains the six-panel maps for all 30
  objects after script 4 runs.
- `../results/figures/` contains manuscript, Extended Data, and
  diagnostic plots generated by scripts 3 and 5-24.

Most scripts set their working directory to `code/` and use code-relative paths.
The scripts use these direct relative paths both in Code Ocean and in the local capsule layout.


## Reproducibility caveats

Script 8 uses the full local raw 3D dust cube, including the cavity, shell,
and surroundings, to set the IAAFT Fourier amplitudes after missing-value
filling and optional block averaging. It preserves the density distribution
exactly and approximates the spectrum. The default is 8 iterations;
iteration count alone is not a convergence guarantee. See
[the method and diagnostic documentation in Script 8](8_OSBs_random_dust_test.py).

Each of the 30 targets also saves `8_SB{id}_original_vs_random.png` (PNG only)
in `../results/figures/random_density_diagnostics/`, using
the first realization and the standard 8-iteration setting. The comparison
has three rows of matched spatial slices (XY, XZ, YZ), with the original
on the left and the random realization on the right. Three panels beneath
them show the full-cube 3D power spectrum, power ratio and density distribution.
Figures retain concise titles, axes, units and legends; detailed method
descriptions and diagnostic numbers are recorded in the notes and metadata
for use in publication captions.
The six spatial slices use a shared `Spectral_r` colour scale and a display-only
Gaussian with sigma 1 voxel per in-plane axis on the actual test grid
(after optional rebinning), applied directly to each 2D slice.
Spectra, density distributions and saved arrays use
the cubes before display smoothing; details are in the method notes.
The arrays, spectral tables and metadata go to
`../results/intermediate_output/8_random_density_fluctuation_test/diagnostics/`.
Paths recorded in the diagnostic metadata are relative to `code/` for
portability; this path convention does not alter the scientific statistics.
From `code/`, run
`python 8_OSBs_random_dust_test.py --targets all --diagnostics-only --iaaft-max-iter 8`
to generate one comparison per target without rewriting Monte Carlo sample or
probability tables.
Use `--diagnostic-dir` and `--comparison-fig-dir` for separate diagnostic runs.

- `0_Data_preparation_2.py` processes the original HI4PI, CfA CO, and Galactic
  Plane Survey / MWISP FITS cubes. Those raw inputs are several hundred GB and
  are not included; the compact downstream PV caches are included.
- Script 8 is the longest standard analysis script because it runs 3D IAAFT
  surrogate-field tests.
- Script 9 automatically checks whether `dustmaps3d` is installed and, when it
  is available, locates the local 3D dust-map data package automatically. For
  normal execution, the standard `dustmaps3d` package cache is sufficient and
  no extra data-path argument is required.
- Scripts 14 and 15 can reuse cached orbit products by default. For a full orbit
  recomputation, run them with `--recompute-orbit --recompute-labels`; this can
  increase runtime.
- The composite final manuscript Fig. 1 is assembled from generated panels,
  as documented above. Script 23 writes manuscript Fig. 3 directly as one figure.
- If `../results/intermediate_output/` and `../results/` are deleted, the
  ordered script sequence can regenerate the reproducible products, provided the
  external `dustmaps3d` cache is available and the capsule `../data/`
  directory is kept.

## Companion online 3D dust-map platform

The companion online 3D dust-map platform for the open-superbubble
visualization is available at:

- https://nadc.china-vo.org/data/dustmaps/OSBs

## License

MIT.

## Contact

Prof. Haibo Yuan: yuanhb@bnu.edu.cn

Tao Wang: wt@mail.bnu.edu.cn
