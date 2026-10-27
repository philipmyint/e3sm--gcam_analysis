#!/bin/bash

set -eo pipefail
trap 'echo "ERROR: script failed at line $LINENO. Exiting." >&2' ERR

MAX_PROCESSES="${MAX_PROCESSES:-2}"
LOAD_R_MODULE=true

while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-processes)
            MAX_PROCESSES="${2:-}"
            if [[ -z "${MAX_PROCESSES}" ]]; then
                echo "ERROR: --max-processes requires a positive integer" >&2
                exit 1
            fi
            shift 2
            ;;
        --skip-r-module-load)
            LOAD_R_MODULE=false
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [--max-processes N] [--skip-r-module-load]"
            echo "Runs GCAM R extraction, processing, synthetic-data generation, and plotting."
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            echo "Usage: $0 [--max-processes N] [--skip-r-module-load]" >&2
            exit 1
            ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/${USER:-e3sm_gcam}_matplotlib}"
mkdir -p "${MPLCONFIGDIR}"

PRE_R_PROJ_LIB_WAS_SET="${PROJ_LIB+x}"
PRE_R_PROJ_LIB="${PROJ_LIB:-}"
PRE_R_PROJ_DATA_WAS_SET="${PROJ_DATA+x}"
PRE_R_PROJ_DATA="${PROJ_DATA:-}"

echo "Starting GCAM synthetic plotting run at: $(date)"
echo "Using Python: $(command -v python3)"
echo "Synthetic-data max processes: ${MAX_PROCESSES}"

python3 - <<'PY'
import importlib.util
import sys

missing = [name for name in ("geopandas", "matplotlib", "numpy", "pandas", "seaborn") if importlib.util.find_spec(name) is None]
if missing:
    print("ERROR: missing required Python package(s): " + ", ".join(missing), file=sys.stderr)
    print("Run this script from an E3SM Unified Python environment.", file=sys.stderr)
    sys.exit(1)
PY

if [[ "${LOAD_R_MODULE}" == true && -z "${__LMOD_REF_COUNT_R:-}" ]]; then
    if command -v module >/dev/null 2>&1; then
        module load R/4.3.3
    fi
fi

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript was not found. Load an R module or rerun with an R environment active." >&2
    exit 1
fi

# Extract GCAM query results from the existing project files.
Rscript gcam_extract_csv_from_xml_or_project_files.R gcam_extract_csv_from_xml_or_project_files.json

# Keep Python GIS libraries from inheriting an incompatible PROJ database path
# from an R module.
if [[ -n "${PRE_R_PROJ_LIB_WAS_SET}" ]]; then
    export PROJ_LIB="${PRE_R_PROJ_LIB}"
else
    unset PROJ_LIB
fi
if [[ -n "${PRE_R_PROJ_DATA_WAS_SET}" ]]; then
    export PROJ_DATA="${PRE_R_PROJ_DATA}"
else
    unset PROJ_DATA
fi

# Prepare processed CSVs and EHC scalar products used by the plotting configs.
python3 gcam_process_extracted_data.py gcam_process_extracted_data.json
python3 gcam_compile_ehc_feedbacks.py gcam_compile_ehc_feedbacks.json

# Generate synthetic ensemble rows used by the ensemble plotting examples.
python3 gcam_produce_synthetic_data.py --max-processes "${MAX_PROCESSES}"

# Produce all GCAM plot families.
python3 gcam_plot_time_series.py gcam_plot_time_series.json
python3 gcam_plot_spatial_data.py gcam_plot_spatial_data.json
python3 gcam_plot_box_and_whiskers.py gcam_plot_box_and_whiskers.json

echo "Finished GCAM synthetic plotting run at: $(date)"
