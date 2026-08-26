#!/bin/bash

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

# First load the E3SM unified environment

# Perlmutter (login or cpu nodes)
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

# compy
#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

# chrysalis
#source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh

# Resolve the directory containing this script so JSON files are found
# regardless of the working directory from which this script is called.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# extract time series of specified H0 variables from E3SM outputs and write to text files
python3 "${SCRIPT_DIR}/../scripts/e3sm_extract_time_series_h0.py" \
    "${SCRIPT_DIR}/e3sm_extract_time_series_h0.json"

# extract time series of specified surfdata variables from E3SM outputs and write to text files
python3 "${SCRIPT_DIR}/../scripts/e3sm_extract_time_series_surfdata_iesm_dyn.py" \
    "${SCRIPT_DIR}/e3sm_extract_time_series_surfdata_iesm_dyn.json"

# extract spatial data of specified H0 variables from E3SM outputs and write to netcdf files
python3 "${SCRIPT_DIR}/../scripts/e3sm_extract_spatial_data_h0.py" \
    "${SCRIPT_DIR}/e3sm_extract_spatial_data_h0.json"

# plot time series of specified H0 and surfdata variables from E3SM outputs
python3 "${SCRIPT_DIR}/../scripts/e3sm_plot_time_series.py" \
    "${SCRIPT_DIR}/e3sm_plot_time_series.json"

# plot spatial data of specified H0 variables from E3SM outputs
python3 "${SCRIPT_DIR}/../scripts/e3sm_plot_spatial_data.py" \
    "${SCRIPT_DIR}/e3sm_plot_spatial_data.json"
