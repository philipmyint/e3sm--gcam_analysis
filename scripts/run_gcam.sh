#!/bin/bash

set -eo pipefail
trap 'echo "ERROR: script failed at line $LINENO. Exiting." >&2' ERR

MAX_PROCESSES="${MAX_PROCESSES:-}"
LOAD_R_MODULE=true
R_MODULE="${R_MODULE:-R/4.3.3}"
USE_LATEX=true
TMP_JSON_DIR=""

cleanup() {
    if [[ -n "${TMP_JSON_DIR}" && -d "${TMP_JSON_DIR}" ]]; then
        rm -rf "${TMP_JSON_DIR}"
    fi
}
trap cleanup EXIT

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
        --r-module)
            R_MODULE="${2:-}"
            if [[ -z "${R_MODULE}" ]]; then
                echo "ERROR: --r-module requires a module name" >&2
                exit 1
            fi
            shift 2
            ;;
        --no-latex)
            USE_LATEX=false
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [--max-processes N] [--r-module NAME] [--skip-r-module-load] [--no-latex]"
            echo "Runs GCAM R extraction, processing, synthetic-data generation, and plotting."
            echo "By default, plot max_processes and use_latex come from each JSON file."
            echo "By default, loads R/4.3.3 only if Rscript is not already available."
            echo "Temporary plotting JSON copies are used only with --max-processes or --no-latex."
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            echo "Usage: $0 [--max-processes N] [--r-module NAME] [--skip-r-module-load] [--no-latex]" >&2
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
if [[ -n "${MAX_PROCESSES}" ]]; then
    echo "Synthetic-data and plot max processes override: ${MAX_PROCESSES}"
else
    echo "Synthetic-data max processes: script default"
    echo "Plot max processes: JSON settings"
fi
echo "Use LaTeX for plots: ${USE_LATEX}"

python3 - <<'PY'
import importlib.util
import sys

missing = [name for name in ("geopandas", "matplotlib", "numpy", "pandas", "seaborn") if importlib.util.find_spec(name) is None]
if missing:
    print("ERROR: missing required Python package(s): " + ", ".join(missing), file=sys.stderr)
    print("Run this script from an E3SM Unified Python environment.", file=sys.stderr)
    sys.exit(1)
PY

if [[ "${LOAD_R_MODULE}" == true ]] && ! command -v Rscript >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        if module is-avail "${R_MODULE}" >/dev/null 2>&1; then
            module load "${R_MODULE}"
        else
            echo "Warning: R module '${R_MODULE}' is not available; continuing with current environment." >&2
            echo "Use --r-module NAME for this machine, or --skip-r-module-load if Rscript is already active." >&2
        fi
    fi
fi

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript was not found. Load an R module or rerun with an R environment active." >&2
    echo "You can also pass --r-module NAME if this machine uses a different R module name." >&2
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
SYNTHETIC_MAX_PROCESS_ARGS=()
if [[ -n "${MAX_PROCESSES}" ]]; then
    SYNTHETIC_MAX_PROCESS_ARGS=(--max-processes "${MAX_PROCESSES}")
fi
python3 gcam_produce_synthetic_data.py "${SYNTHETIC_MAX_PROCESS_ARGS[@]}"

TIME_SERIES_JSON="gcam_plot_time_series.json"
SPATIAL_JSON="gcam_plot_spatial_data.json"
BOX_JSON="gcam_plot_box_and_whiskers.json"
if [[ -n "${MAX_PROCESSES}" || "${USE_LATEX}" != true ]]; then
    TMP_JSON_DIR="$(mktemp -d "${TMPDIR:-/tmp}/gcam-run-json.XXXXXX")"
    python3 - "${TMP_JSON_DIR}" "${MAX_PROCESSES}" "${USE_LATEX}" "${TIME_SERIES_JSON}" "${SPATIAL_JSON}" "${BOX_JSON}" <<'PY'
import json
import os
import sys

output_dir = sys.argv[1]
max_processes = int(sys.argv[2]) if sys.argv[2] else None
use_latex = sys.argv[3] == "true"
input_files = sys.argv[4:]

def update_plot_options(value):
    if isinstance(value, dict):
        updated = {key: update_plot_options(item) for key, item in value.items()}
        if max_processes is not None:
            updated["max_processes"] = max_processes
        if not use_latex:
            updated["use_latex"] = False
        return updated
    if isinstance(value, list):
        return [update_plot_options(item) for item in value]
    return value

for input_file in input_files:
    with open(input_file, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    output_file = os.path.join(output_dir, os.path.basename(input_file))
    with open(output_file, "w", encoding="utf-8") as handle:
        json.dump(update_plot_options(data), handle, indent=4)
        handle.write("\n")
PY
    TIME_SERIES_JSON="${TMP_JSON_DIR}/${TIME_SERIES_JSON}"
    SPATIAL_JSON="${TMP_JSON_DIR}/${SPATIAL_JSON}"
    BOX_JSON="${TMP_JSON_DIR}/${BOX_JSON}"
fi

# Produce all GCAM plot families.
python3 gcam_plot_time_series.py "${TIME_SERIES_JSON}"
python3 gcam_plot_spatial_data.py "${SPATIAL_JSON}"
python3 gcam_plot_box_and_whiskers.py "${BOX_JSON}"

echo "Finished GCAM synthetic plotting run at: $(date)"
