#!/bin/bash

set -eo pipefail
trap 'echo "ERROR: script failed at line $LINENO. Exiting." >&2' ERR

MAX_PROCESSES="${MAX_PROCESSES:-2}"

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
        -h|--help)
            echo "Usage: $0 [--max-processes N]"
            echo "Runs E3SM synthetic-data generation and plotting from the scripts directory."
            echo "Does not run E3SM extraction scripts."
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            echo "Usage: $0 [--max-processes N]" >&2
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
echo "Synthetic-data max processes: ${MAX_PROCESSES}"

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
python3 e3sm_produce_synthetic_time_series.py --max-processes "${MAX_PROCESSES}"
python3 e3sm_produce_synthetic_spatial_data.py --max-processes "${MAX_PROCESSES}"

# Plot from existing E3SM time-series and spatial NetCDF products. The extraction
# scripts are intentionally skipped because their source simulation paths are not
# expected to exist in this repository checkout.
python3 e3sm_plot_time_series.py e3sm_plot_time_series.json
python3 e3sm_plot_spatial_data.py e3sm_plot_spatial_data_elm.json
python3 e3sm_plot_spatial_data.py e3sm_plot_spatial_data_eam.json

echo "Finished E3SM synthetic plotting run at: $(date)"
