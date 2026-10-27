#!/bin/bash
#SBATCH --time=01:00:00

set -eo pipefail
trap 'echo "ERROR: script failed at line $LINENO. Exiting." >&2' ERR

MAX_PROCESSES="${MAX_PROCESSES:-}"
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
        --no-latex)
            USE_LATEX=false
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [--max-processes N] [--no-latex]"
            echo "Runs E3SM synthetic-data generation and plotting from the scripts directory."
            echo "Does not run E3SM extraction scripts."
            echo "By default, plot max_processes and use_latex come from each JSON file."
            echo "Temporary plotting JSON copies are used only with --max-processes or --no-latex."
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            echo "Usage: $0 [--max-processes N] [--no-latex]" >&2
            exit 1
            ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/${USER:-e3sm_gcam}_matplotlib}"
mkdir -p "${MPLCONFIGDIR}"

echo "Starting E3SM synthetic plotting run at: $(date)"
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

missing = [name for name in ("xarray", "matplotlib", "numpy", "pandas") if importlib.util.find_spec(name) is None]
if missing:
    print("ERROR: missing required Python package(s): " + ", ".join(missing), file=sys.stderr)
    print("Run this script from an E3SM Unified Python environment.", file=sys.stderr)
    sys.exit(1)
PY

# Generate synthetic ensemble members used by the ensemble plotting examples.
SYNTHETIC_MAX_PROCESS_ARGS=()
if [[ -n "${MAX_PROCESSES}" ]]; then
    SYNTHETIC_MAX_PROCESS_ARGS=(--max-processes "${MAX_PROCESSES}")
fi
python3 e3sm_produce_synthetic_time_series.py "${SYNTHETIC_MAX_PROCESS_ARGS[@]}"
python3 e3sm_produce_synthetic_spatial_data.py "${SYNTHETIC_MAX_PROCESS_ARGS[@]}"

TIME_SERIES_JSON="e3sm_plot_time_series.json"
SPATIAL_ELM_JSON="e3sm_plot_spatial_data_elm.json"
SPATIAL_EAM_JSON="e3sm_plot_spatial_data_eam.json"
if [[ -n "${MAX_PROCESSES}" || "${USE_LATEX}" != true ]]; then
    TMP_JSON_DIR="$(mktemp -d "${TMPDIR:-/tmp}/e3sm-run-json.XXXXXX")"
    python3 - "${TMP_JSON_DIR}" "${MAX_PROCESSES}" "${USE_LATEX}" "${TIME_SERIES_JSON}" "${SPATIAL_ELM_JSON}" "${SPATIAL_EAM_JSON}" <<'PY'
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
    SPATIAL_ELM_JSON="${TMP_JSON_DIR}/${SPATIAL_ELM_JSON}"
    SPATIAL_EAM_JSON="${TMP_JSON_DIR}/${SPATIAL_EAM_JSON}"
fi

# Plot from existing E3SM time-series and spatial NetCDF products. The extraction
# scripts are intentionally skipped because their source simulation paths are not
# expected to exist in this repository checkout.
python3 e3sm_plot_time_series.py "${TIME_SERIES_JSON}"
python3 e3sm_plot_spatial_data.py "${SPATIAL_ELM_JSON}"
python3 e3sm_plot_spatial_data.py "${SPATIAL_EAM_JSON}"

echo "Finished E3SM synthetic plotting run at: $(date)"
