#!/bin/bash
set -e
trap 'echo "ERROR: script failed at line $LINENO. Exiting." >&2' ERR

# This script runs the Python scripts that extract and plot GCAM outputs.
# These JSON config files are configured to make simple plots for two cases for a default set of GCAM outputs.

# There are many plot options available once the desired data are extracted.
# See documentation and the JSON file in the scripts directory for examples of other plot options.

# All files will be written by default to an 'output_test' directory in the current working directory (./output_test).
# To change this default, change the text in the JSON files in this workflow directory
#   gcam_extract...json: output_files(s)
#   gcam_process...json: output_files(s)
#   gcam_plot...json: plot_directory, netcdf_files
# Note that the output directory specifed in the JSON files will be created if it does not already exist.

# It is recommended to set the name of the 'output' directory to a unique name for each set of
#    E3SM-GCAM outputs that are to be analyzed by extracting data and plotting.
# This will avoid overwriting previous outputs and plots.

# Resolve the directory containing this script so JSON files are found
# regardless of the working directory from which this script is called.
export SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# First load the R module and make sure the needed libraries are installed.
module load R/4.3.3

# perlmutter (and compy?) has the required R packages, but chrysalis does not
# chrysalis also does not have the required system libraries to install the required R packages
# so cannot run the R script on chrysalis.

# Ensure environment includes all required paths for compiling, and also the custom user library path for R packages.
#Rscript -e '.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.3", .libPaths()))'

# Perform local user install of required R libraries if they are not already installed.
# Change lib="" to a different directory if you want to install the libraries somewhere else,
#    but note that the one listed below is the current default on Chrysalis (R v4.3.3) as of 23 july 2026 (for R versions 4.3.#).

# devtools, not necessary to run rgcam, only needed for installing rgcam, but remotes will also work
#Rscript -e 'if(!requireNamespace("curl", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, quietly=TRUE)) install.packages("curl", lib="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, repos="https://cloud.r-project.org")'
#Rscript -e 'if(!requireNamespace("httr2", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, quietly=TRUE)) install.packages("httr2", lib="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, repos="https://cloud.r-project.org")'
#Rscript -e 'if(!requireNamespace("xml2", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, quietly=TRUE)) install.packages("xml2", lib="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, repos="https://cloud.r-project.org")'
#Rscript -e 'if(!requireNamespace("gert", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, quietly=TRUE)) install.packages("gert", lib="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, repos="https://cloud.r-project.org")'
#Rscript -e 'if(!requireNamespace("gh", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, quietly=TRUE)) install.packages("gh", lib="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, repos="https://cloud.r-project.org")'
#Rscript -e 'if(!requireNamespace("usethis", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, quietly=TRUE)) install.packages("usethis", lib="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, repos="https://cloud.r-project.org")'
#Rscript -e 'if(!requireNamespace("devtools", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, quietly=TRUE)) install.packages("devtools", lib="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, repos="https://cloud.r-project.org")'

# rgcam and remotes
#Rscript -e 'target_dir <- "~/R/x86_64-pc-linux-gnu-library/4.3"; if (!requireNamespace("rgcam", lib.loc = target_dir, quietly = TRUE)) { if (!requireNamespace("remotes", lib.loc = target_dir, quietly = TRUE)) install.packages("remotes", lib = target_dir, repos = "https://cloud.r-project.org"); remotes::install_github("JGCRI/rgcam", lib = target_dir, dependencies = TRUE) } else { message("rgcam is already installed.") }'
# dplyr - it installed on chrysalis, but with errors, so it may not work correctly
#Rscript -e 'if(!requireNamespace("dplyr", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, quietly=TRUE)) install.packages("dplyr", lib="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, repos="https://cloud.r-project.org")'
# rjson - this installed successfully on chrysalis
#Rscript -e 'if(!requireNamespace("rjson", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, quietly=TRUE)) install.packages("rjson", lib="~/R/x86_64-pc-linux-gnu-library/4.3", dependencies=TRUE, repos="https://cloud.r-project.org")'

# extract variables from GCAM outputs and write each variable to a csv file
# the first json block extracts data from the xml databases into .dat files and writes the first csv files,
#    and the rest of the blocks extract data from the .dat files into csv files
Rscript "${SCRIPT_DIR}/../scripts/gcam_extract_csv_from_xml_or_project_files.R" \
    "${SCRIPT_DIR}/gcam_extract_csv_from_xml_or_project_files.json"

# Next load the E3SM unified environment

# Perlmutter (login or cpu nodes)
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

# compy
#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

# chrysalis
#source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh

# process the extracted csv files to create 'processed' csv files for plotting
python3 "${SCRIPT_DIR}/../scripts/gcam_process_extracted_data.py" \
    "${SCRIPT_DIR}/gcam_process_extracted_data.json"

# process the ehc diagnostics into single variable csv files for plotting
python3 "${SCRIPT_DIR}/../scripts/gcam_compile_ehc_scalars.py" \
    "${SCRIPT_DIR}/gcam_compile_ehc_scalars.json"

# plot time series of specified variables from GCAM outputs
python3 "${SCRIPT_DIR}/../scripts/gcam_plot_time_series.py" \
    "${SCRIPT_DIR}/gcam_plot_time_series.json"
