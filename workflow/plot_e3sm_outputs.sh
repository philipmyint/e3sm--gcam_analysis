#!/bin/bash

# set some sbatch directives for optional job submission (sbatch plot_e3sm_outputs.sh)
#SBATCH --job-name=plot_e3sm_outputs
#SBATCH --output=plot_e3sm_outputs.out
#SBATCH --error=plot_e3sm_outputs.err
#SBATCH --time=08:00:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -C cpu
#SBATCH -A e3sm

set -eo pipefail

# This script runs the Python scripts that extract and plot E3SM outputs (ELM and EAM).
# These JSON config files are configured to make simple plots for two cases for a default set of E3SM outputs.

# There are many plot options available once the desired data are extracted.
# See documentation and the JSON file in the scripts directory for examples of other plot options.

# All files will be written by default to an 'output' directory in the current working directory (./output).
# To change this default, change the text in the JSON files in this workflow directory
#   e3sm_extract...json: output_files(s)
#   e3sm_plot...json: plot_directory, netcdf_files
# Note that the output directory specifed in the JSON files will be created if it does not already exist.

# It is recommended to set the name of the 'output' directory to a unique name for each set of
#    E3SM-GCAM outputs that are to be analyzed by extracting data and plotting.
# This will avoid overwriting previous outputs and plots.

# Check that the correct grid scrip file is set in e3sm_plot_spatial_data.json for EAM spatial plots.
#    The default grid file is set to the ne30pg2 grid file in this repository.

MACHINE="perlmutter"
OVERWRITE_ARGS=()
E3SM_ENV_SCRIPT=""
LOAD_ENV=true

# Parse optional command-line arguments. The --machine flag selects a known
# E3SM unified environment to source; with no flag, default to Perlmutter.
# The --overwrite flag is forwarded to extraction scripts so existing extracted
# output files are replaced before plotting.
while [[ $# -gt 0 ]]; do
    case "$1" in
        --machine)
            MACHINE="${2:-}"
            if [[ -z "${MACHINE}" ]]; then
                echo "ERROR: --machine requires one of: perlmutter, compy, chrysalis" >&2
                exit 1
            fi
            shift 2
            ;;
        --overwrite)
            OVERWRITE_ARGS=(--overwrite)
            shift
            ;;
        --env-script)
            E3SM_ENV_SCRIPT="${2:-}"
            if [[ -z "${E3SM_ENV_SCRIPT}" ]]; then
                echo "ERROR: --env-script requires a path to an environment loader script" >&2
                exit 1
            fi
            shift 2
            ;;
        --skip-env-load)
            LOAD_ENV=false
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [--machine perlmutter|compy|chrysalis] [--overwrite] [--env-script PATH] [--skip-env-load]"
            echo "Default: --machine perlmutter"
            echo "--overwrite replaces existing extracted output files before plotting"
            echo "--env-script PATH sources a specific E3SM Unified loader script"
            echo "--skip-env-load uses the currently active Python environment"
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            echo "Usage: $0 [--machine perlmutter|compy|chrysalis] [--overwrite] [--env-script PATH] [--skip-env-load]" >&2
            exit 1
            ;;
    esac
done

echo "Starting job at: $(date)"

# Load the E3SM unified environment for the selected machine.
if [[ "${LOAD_ENV}" == true && -z "${E3SM_ENV_SCRIPT}" ]]; then
    case "${CONDA_DEFAULT_ENV:-}" in
        e3sm_unified*|e3sm-unified*)
            echo "Using active E3SM Unified environment: ${CONDA_DEFAULT_ENV}"
            LOAD_ENV=false
            ;;
        "")
            ;;
        *)
            echo "Active Conda environment is '${CONDA_DEFAULT_ENV}', not an E3SM Unified environment; loading ${MACHINE} default."
            ;;
    esac
fi
if [[ "${LOAD_ENV}" == true && -z "${E3SM_ENV_SCRIPT}" ]]; then
    case "${MACHINE}" in
        perlmutter)
            E3SM_ENV_SCRIPT="/global/common/software/e3sm/anaconda_envs/load_e3sm_unified_1.12.0_pm-cpu.sh"
            ;;
        compy)
            E3SM_ENV_SCRIPT="/share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh"
            ;;
        chrysalis)
            E3SM_ENV_SCRIPT="/lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh"
            ;;
        *)
            echo "ERROR: unknown machine: ${MACHINE}" >&2
            echo "Valid machines: perlmutter, compy, chrysalis" >&2
            exit 1
        ;;
    esac
fi
if [[ "${LOAD_ENV}" == true ]]; then
    if [[ ! -f "${E3SM_ENV_SCRIPT}" ]]; then
        echo "ERROR: E3SM environment loader not found: ${E3SM_ENV_SCRIPT}" >&2
        echo "Use --env-script PATH to point to a valid E3SM Unified loader script." >&2
        exit 1
    fi
    if ! source "${E3SM_ENV_SCRIPT}"; then
        echo "ERROR: failed to source E3SM environment loader: ${E3SM_ENV_SCRIPT}" >&2
        echo "Submit from an active E3SM Unified environment or use --env-script PATH / --skip-env-load." >&2
        exit 1
    fi
else
    echo "Skipping E3SM environment load; using active Python: $(command -v python3)"
fi

# Resolve the directory containing this script so JSON files are found
# regardless of the working directory from which this script is called.
# Under sbatch, the script is copied to a spool file, so BASH_SOURCE no longer
# points here; use SLURM_SUBMIT_DIR (the sbatch submission directory) instead.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# extract time series of specified H0 variables from E3SM outputs and write to text files
python3 "${SCRIPT_DIR}/../scripts/e3sm_extract_time_series_h0.py" \
    "${OVERWRITE_ARGS[@]}" \
    "${SCRIPT_DIR}/e3sm_extract_time_series_h0.json"

# extract time series of specified surfdata variables from E3SM outputs and write to text files
python3 "${SCRIPT_DIR}/../scripts/e3sm_extract_time_series_landuse.py" \
    "${OVERWRITE_ARGS[@]}" \
    "${SCRIPT_DIR}/e3sm_extract_time_series_landuse.json"

# extract spatial data of specified H0 variables from E3SM outputs and write to netcdf files
python3 "${SCRIPT_DIR}/../scripts/e3sm_extract_spatial_data_h0.py" \
    "${OVERWRITE_ARGS[@]}" \
    "${SCRIPT_DIR}/e3sm_extract_spatial_data_h0.json"

# plot time series of specified H0 and surfdata variables from E3SM outputs
python3 "${SCRIPT_DIR}/../scripts/e3sm_plot_time_series.py" \
    "${SCRIPT_DIR}/e3sm_plot_time_series.json"

# plot spatial data of specified H0 variables from E3SM outputs
python3 "${SCRIPT_DIR}/../scripts/e3sm_plot_spatial_data.py" \
    "${SCRIPT_DIR}/e3sm_plot_spatial_data.json"

echo "Finished job at: $(date)"
