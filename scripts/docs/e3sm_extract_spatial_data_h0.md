# E3SM Spatial Data Extraction Script Documentation

## Overview

**Script Name:** `e3sm_extract_spatial_data_h0.py`

**Purpose:** Extracts spatially-gridded annual-mean data from E3SM monthly NetCDF h0 output files (ELM and/or EAM) and creates smaller, focused NetCDF files containing only user-specified variables for a specific time period. Unlike time series extraction which aggregates spatially, this script preserves full latitude/longitude resolution while aggregating temporally (annual means from monthly data).

**Key Use Case:** Creating manageable NetCDF files for spatial analysis, mapping, and visualization of E3SM output variables over specific time periods without the overhead of processing entire simulation archives.

**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Requirements](#installation--requirements)
3. [Basic Usage](#basic-usage)
4. [Key Concepts](#key-concepts)
5. [Complete Parameter Reference Table](#complete-parameter-reference-table)
6. [Detailed Parameter Descriptions](#detailed-parameter-descriptions)
7. [Variable Processing](#variable-processing)
8. [JSON Configuration Examples](#json-configuration-examples)
9. [Output Files](#output-files)
10. [Troubleshooting](#troubleshooting)
11. [Best Practices](#best-practices)
12. [Comparison with Time Series Extraction](#comparison-with-time-series-extraction)

---

## Installation & Requirements

### Required Python Libraries

```bash
pip install pandas numpy xarray netCDF4 dask multiprocessing
```

### Required Utility Modules

The script imports several utility modules that must be in the same directory or Python path:
- `utility_constants` - Physical constants (conversion factors)
- `utility_functions` - General utility functions
- `utility_e3sm_netcdf` - E3SM-specific NetCDF handling functions

### System Requirements

- Python 3.7+
- Multi-core processor (script uses parallel processing)
- Sufficient RAM (8+ GB recommended for large datasets)
- Sufficient disk space for output NetCDF files
- Access to E3SM simulation output directories

---

## Basic Usage

### Command Line Execution

```bash
python e3sm_extract_spatial_data_h0.py path/to/config.json
```

**Multiple Configuration Files:**
```bash
python e3sm_extract_spatial_data_h0.py config1.json config2.json config3.json
```

### What the Script Does

For each configuration block in the JSON file, the script:
1. **Locates** monthly E3SM NetCDF h0 files (elm.h0.*, eam.h0.*)
2. **Filters** files by year range (start_years to end_years)
3. **Opens** multiple NetCDF files efficiently using xarray's parallel processing
4. **Extracts** only user-specified variables
5. **Computes** annual means from monthly data (temporal aggregation)
6. **Processes** variables (adds total precipitation, CO₂ mole fraction)
7. **Writes** compact NetCDF file with spatial data
8. **Uses** parallel processing for multiple output files

---

## Key Concepts

### Spatial vs Time Series Extraction

**This Script (Spatial Extraction):**
- **Preserves:** Full latitude/longitude resolution
- **Aggregates:** Temporal (monthly → annual mean)
- **Output:** NetCDF file with gridded data
- **Use:** Mapping, spatial analysis, visualization

**Time Series Extraction Script:**
- **Preserves:** Temporal resolution (monthly or annual)
- **Aggregates:** Spatial (lat/lon → global or regional mean/sum)
- **Output:** Text file with time series
- **Use:** Trend analysis, statistical analysis

### Annual Mean Calculation

The script converts monthly E3SM output to annual means:

```
Monthly files for 2085:
elm.h0.2085-01.nc, elm.h0.2085-02.nc, ..., elm.h0.2085-12.nc

↓ [Annual mean calculated]

Annual mean for 2085 at each lat/lon gridcell
```

**Formula:**
```
Annual_mean(lat, lon) = mean(Jan, Feb, Mar, ..., Dec)
```

### Multi-File Processing

The script efficiently handles multiple monthly files using xarray's `open_mfdataset`:

```python
# Opens all monthly files in parallel
ds = xr.open_mfdataset(
    netcdf_files,              # List of monthly files
    parallel=True,             # Parallel I/O
    combine='nested',          # Combine along time
    concat_dim='time'          # Concatenate on time dimension
)
```

---

## Complete Parameter Reference Table

### Required Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `simulation_path` | string | **Yes** | Path to E3SM simulation directory |
| `output_files` | string or list | **Yes** | Path(s) to output NetCDF file(s) |
| `netcdf_substrings` | nested list | **Yes** | File type identifiers for NetCDF files |
| `variables` | nested list | **Yes** | Variables to extract from each file type |
| `start_years` | integer or list | **Yes** | First year(s) to extract |
| `end_years` | integer or list | **Yes** | Last year(s) to extract |

### No Optional Parameters

This script requires all parameters to be specified. Unlike other scripts in the suite, there are no default values for optional parameters.

---

## Detailed Parameter Descriptions

### Core Required Parameters

#### `simulation_path`
**Type:** String (directory path)  
**Required:** Yes  
**Description:** Complete path to directory containing E3SM NetCDF h0 output files.

**Examples:**
```json
"simulation_path": "/lcrc/group/e3sm/ac.eva.sinha/20240730_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3_without_feedbacks/run"
```

**Requirements:**
- Must be valid directory path
- Must contain monthly NetCDF (.nc) files
- Files must follow E3SM naming convention

**Typical Structure:**
```
simulation_path/
├── elm.h0.2015-01.nc
├── elm.h0.2015-02.nc
├── ...
├── eam.h0.2015-01.nc
├── eam.h0.2015-02.nc
└── ...
```

---

#### `output_files`
**Type:** String OR List of strings  
**Required:** Yes  
**Description:** Path(s) to output NetCDF file(s). Can be a single file or list of files.

**Format 1 - Single Output File:**
```json
"output_files": "./output/spatial_data.nc"
```

**Format 2 - Multiple Output Files:**
```json
"output_files": [
    "./output/spatial_elm.nc",
    "./output/spatial_eam.nc"
]
```

**Automatic Conversion:**
If you provide a single string, the script converts it to a list:
```json
// You write:
"output_files": "file.nc"

// Script converts to:
"output_files": ["file.nc"]
```

**File Format:**
- Must be NetCDF format (.nc extension)
- Creates standard NetCDF4 file
- Contains spatial gridded data

---

#### `netcdf_substrings`
**Type:** Nested list of strings OR list of strings  
**Required:** Yes  
**Description:** Identifies different types of NetCDF files by substrings in filenames.

**Format:** List of lists, where each inner list contains substrings for one file type.

**Examples:**

**Single file type (ELM only):**
```json
"netcdf_substrings": [["elm.h0"]]
```

**Two file types (ELM and EAM):**
```json
"netcdf_substrings": [["elm.h0"], ["eam.h0"]]
```

**Automatic List Conversion:**
```json
// If you write:
"netcdf_substrings": ["elm.h0"]

// Script converts to:
"netcdf_substrings": [["elm.h0"]]
```

**How It Works:**
- Script searches `simulation_path` for .nc files
- Files must contain ALL substrings in their corresponding list
- Each file type creates separate output file

---

#### `variables`
**Type:** Nested list of strings OR list of strings  
**Required:** Yes  
**Description:** Variables to extract from each NetCDF file type.

**Format:** List of lists, where each inner list corresponds to one file type.

**Structure:**
```
variables[i] = variables to extract from netcdf_substrings[i] files
```

**Example:**
```json
{
    "output_files": ["elm_spatial.nc", "eam_spatial.nc"],
    "netcdf_substrings": [["elm.h0"], ["eam.h0"]],
    "variables": [
        ["GPP", "NPP", "NBP", "TOTECOSYSC"],
        ["PRECC", "PRECL", "TREFHT"]
    ]
}
```

**Breakdown:**
- `["GPP", "NPP", "NBP", "TOTECOSYSC"]` → elm_spatial.nc (from elm.h0 files)
- `["PRECC", "PRECL", "TREFHT"]` → eam_spatial.nc (from eam.h0 files)

**Automatic List Conversion:**
```json
// If you write:
"variables": ["GPP", "NPP"]

// And output_files has length 1:
// Script keeps as:
"variables": [["GPP", "NPP"]]
```

**Common ELM Variables:**
```json
"variables": [
    ["GPP", "NPP", "NBP", "ER", "HR",
     "TOTECOSYSC", "TOTVEGC", "TOTSOMC",
     "LAND_UPTAKE", "LAND_USE_FLUX"]
]
```

**Common EAM Variables:**
```json
"variables": [
    ["PRECC", "PRECL", "TREFHT", "PBOT",
     "SFCO2", "SFCO2_LND", "SFCO2_OCN"]
]
```

---

#### `start_years`
**Type:** Integer OR List of integers  
**Required:** Yes  
**Description:** First year(s) for time period to extract.

**Format 1 - Single Value (Applied to All):**
```json
"start_years": 2085
```

**Format 2 - Multiple Values (One Per Output File):**
```json
"start_years": [2085, 2085]
```

**Automatic Value Replication:**
```json
// If you write:
"start_years": 2085

// And output_files has length 2:
// Script converts to:
"start_years": [2085, 2085]
```

**Usage:**
- Defines start of time window
- Script extracts monthly files from start_years to end_years
- Computes annual means for each year in range

---

#### `end_years`
**Type:** Integer OR List of integers  
**Required:** Yes  
**Description:** Last year(s) for time period to extract.

**Format 1 - Single Value (Applied to All):**
```json
"end_years": 2090
```

**Format 2 - Multiple Values (One Per Output File):**
```json
"end_years": [2090, 2090]
```

**Automatic Value Replication:**
```json
// If you write:
"end_years": 2090

// And output_files has length 2:
// Script converts to:
"end_years": [2090, 2090]
```

**Time Period:**
```
start_years = 2085, end_years = 2090

→ Extracts: 2085, 2086, 2087, 2088, 2089, 2090 (6 years)
→ Processes: 72 monthly files (6 years × 12 months)
→ Output: 6 annual means at each gridcell
```

---

## Variable Processing

The script automatically processes certain variables when they are present:

### 1. Total Precipitation Calculation

**Requires:** `PRECC`, `PRECL`, `PRECSC`, `PRECSL` (all must be present)

**Processing:**
1. Convert units: m/s → mm/year
2. Create new variable: `PRECIP` = sum of all precipitation components

**Formula:**
```python
# Convert each component
PRECC (mm/year) = PRECC (m/s) × seconds_per_year × 1000 mm/m
PRECL (mm/year) = PRECL (m/s) × seconds_per_year × 1000 mm/m
PRECSC (mm/year) = PRECSC (m/s) × seconds_per_year × 1000 mm/m
PRECSL (mm/year) = PRECSL (m/s) × seconds_per_year × 1000 mm/m

# Calculate total
PRECIP (mm/year) = PRECC + PRECL + PRECSC + PRECSL
```

**Output Variables:**
- `PRECC` - Convective precipitation (mm/year)
- `PRECL` - Large-scale precipitation (mm/year)
- `PRECSC` - Convective snow (mm/year)
- `PRECSL` - Large-scale snow (mm/year)
- `PRECIP` - **New:** Total precipitation (mm/year)

**Attributes:**
```
PRECIP:
  units: "mm/year"
  description: "Total precipitation rate"
```

---

### 2. CO₂ Mole Fraction Calculation

**Requires:** `PBOT`, `PCO2`, `QBOT` (all must be present)

**Processing:**
Calculate atmospheric CO₂ mole fraction in dry air (ppm)

**Formula:**
```python
# Partial pressure of water vapor (Pa)
# See: https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
pH2O = PBOT × QBOT / (0.622 + 0.378 × QBOT)

# CO₂ mole fraction in dry air (ppm)
ZCO2 (ppm) = 1e6 × PCO2 / (PBOT - pH2O)
```

**Output Variable:**
- `ZCO2` - **New:** CO₂ mole fraction in dry air (ppm)

**Attributes:**
```
ZCO2:
  units: "ppm"
  description: "CO2 mole fraction in dry air"
```

**Scientific Background:**
- CO₂ concentrations reported as mole fraction in *dry* air
- Must remove water vapor contribution
- Uses standard psychrometric formula

---

## JSON Configuration Examples

### Example 1: Basic Configuration (ELM + EAM)

**Purpose:** Extract land and atmosphere variables for 2085-2090

```json
{
    "simulation_path": "/lcrc/group/e3sm/ac.eva.sinha/20240730_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3_without_feedbacks/run",
    "output_files": [
        "./../2025_DiVittorio_et_al_e3sm/control_spatial_data_elm.nc",
        "./../2025_DiVittorio_et_al_e3sm/control_spatial_data_eam.nc"
    ],
    "netcdf_substrings": [["elm.h0"], ["eam.h0"]],
    "variables": [
        ["DWT_CONV_CFLUX_GRC", "ER", "GPP", "HR", "LAND_UPTAKE", "LAND_USE_FLUX",
         "NBP", "NEE", "NEP", "NPP", "PFT_FIRE_CLOSS", "WOOD_HARVESTC",
         "WOODC_ALLOC", "WOODC_LOSS", "PBOT", "PCO2", "QBOT", "TBOT",
         "TOTECOSYSC", "TOTSOMC", "TOTLITC", "TOTVEGC", "TOTVEGC_ABG", "WOODC"],
        ["PRECC", "PRECL", "PRECSC", "PRECSL", "SFCO2", "SFCO2_FFF",
         "SFCO2_LND", "SFCO2_OCN", "TMCO2", "TMCO2_FFF", "TMCO2_LND",
         "TMCO2_OCN", "TREFHT"]
    ],
    "start_years": 2085,
    "end_years": 2090
}
```

**What it does:**
- Creates 2 output files (ELM and EAM)
- ELM file: 24 land variables
- EAM file: 13 atmospheric variables + PRECIP (added automatically)
- Time period: 2085-2090 (6 years)
- Output includes ZCO2 (derived from PBOT, PCO2, QBOT)

---

### Example 2: Single File Type

**Purpose:** Extract only land model outputs

```json
{
    "simulation_path": "/lcrc/group/e3sm/simulation/run",
    "output_files": "./output/elm_spatial_data.nc",
    "netcdf_substrings": [["elm.h0"]],
    "variables": [["GPP", "NPP", "NBP", "TOTECOSYSC"]],
    "start_years": 2085,
    "end_years": 2090
}
```

---

### Example 3: Multiple Simulations

**Purpose:** Process multiple simulations in batch

```json
[
    {
        "simulation_path": "/lcrc/group/e3sm/control/run",
        "output_files": ["./control_elm.nc", "./control_eam.nc"],
        "netcdf_substrings": [["elm.h0"], ["eam.h0"]],
        "variables": [["GPP", "NPP", "NBP"], ["PRECC", "PRECL", "TREFHT"]],
        "start_years": 2085,
        "end_years": 2090
    },
    {
        "simulation_path": "/lcrc/group/e3sm/full_feedback/run",
        "output_files": ["./full_feedback_elm.nc", "./full_feedback_eam.nc"],
        "netcdf_substrings": [["elm.h0"], ["eam.h0"]],
        "variables": [["GPP", "NPP", "NBP"], ["PRECC", "PRECL", "TREFHT"]],
        "start_years": 2085,
        "end_years": 2090
    }
]
```

---

## Output Files

### NetCDF File Structure

**Format:** NetCDF4

**Typical Dimensions:**
```
Dimensions:
  lat: 192 (grid-dependent)
  lon: 288 (grid-dependent)
  year: 6 (2085-2090)
```

**Variables:**
```
lat(lat): Latitude coordinates
lon(lon): Longitude coordinates
year(year): Year values
[User variables](year, lat, lon): Data arrays
```

### Using Output Files

**Load in Python:**
```python
import xarray as xr
ds = xr.open_dataset('spatial_elm.nc')
print(ds)

# Plot
ds['GPP'].sel(year=2085).plot()
```

---

## Troubleshooting

### Issue 1: Memory Error

**Solutions:**
- Reduce variables
- Shorten time period
- Process in batches
- Use HPC with more RAM

### Issue 2: File Not Found

**Solutions:**
- Verify simulation_path
- Check netcdf_substrings match files
- Verify year range has files

### Issue 3: Variable Not Found

**Solutions:**
- Check variable exists in files
- Verify correct file type (ELM vs EAM)
- Check for typos (case-sensitive)

---

## Best Practices

1. **Test first** with short periods and few variables
2. **Name files descriptively** (include simulation, years, component)
3. **Process in parallel** by using multiple output_files
4. **Document** your configurations
5. **Verify outputs** after extraction

---

## Comparison with Time Series Extraction

| Feature | Spatial Script | Time Series Script |
|---------|----------------|-------------------|
| Output | NetCDF (gridded) | Text (time series) |
| Spatial | Full resolution | Aggregated |
| Temporal | Annual means | Monthly/annual |
| Use | Mapping | Trend analysis |

---

## References

- E3SM Documentation: [https://e3sm.org/](https://e3sm.org/)
- xarray Documentation: [https://xarray.pydata.org/](https://xarray.pydata.org/)
- DiVittorio et al. (2025). [DOI: 10.1029/2024MS004806](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024MS004806)
- GitHub Repository: [https://github.com/philipmyint/e3sm--gcam_analysis](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Contact

- **Philip Myint**: myint1@llnl.gov
- **Dalei Hao**: dalei.hao@pnnl.gov

---

*Documentation for `e3sm_extract_spatial_data_h0.py` - Last Updated: January 2026*
