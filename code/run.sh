#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAPSULE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

DATA_ROOT="${CAPSULE_DIR}/data"
RESULTS_ROOT="${CAPSULE_DIR}/results"

DATA_INPUT_DIR="${DATA_ROOT}"
INTERMEDIATE_OUTPUT_DIR="${RESULTS_ROOT}/intermediate_output"

if [[ ! -d "${DATA_INPUT_DIR}" ]]; then
  echo "Missing input data directory: ${DATA_INPUT_DIR}" >&2
  exit 1
fi

mkdir -p "${RESULTS_ROOT}" "${INTERMEDIATE_OUTPUT_DIR}"

cd "${SCRIPT_DIR}"

python3 0_Data_preparation_1.py

# 0_Data_preparation_2.py is intentionally not executed in this capsule.
# It documents PV-cache export from the original MWISP/CfA/HI4PI survey data,
# which are too large to redistribute with the publication package.

python3 1_OSBs_first_contact_MCs.py
python3 2_OSBs_3D_fitting.py
python3 3_OSBs_MCMC_and_parameter_table.py
python3 4_OSBs_plotting.py
python3 5_OSBs_parameter_statistics.py
python3 6_OSBs_initial_sensitivity.py
python3 7_OSBs_recall_analysis.py
python3 8_OSBs_random_dust_test.py
python3 9_OSBs_shell_amplitude_test.py
python3 10_OSBs_bubbles_HMSFR_test.py
python3 11_OSBs_top_view.py
python3 12_OSBs_dust_slices.py
python3 13_YSO_traceback_inputs.py
python3 14_G1_YSO_traceback.py
python3 15_G2_YSO_traceback.py
python3 16_G1G2_YSO_traceback_plots.py
python3 17_SB31_RW_boundary_YSO_traceback.py
python3 18_YSO_OSBs_shell_random_test.py
python3 19_OSBs_energy_budget.py
python3 20_OSBs_YSO_SN_estimation.py
python3 21_OSBs_gas_PV_diagrams.py
python3 22_bubbles_gas_PV_diagram.py
python3 23_RW_cross_sections.py
python3 24_RW_like_wave.py
python3 25_interactive_three_volume_viewer.py
