# E3SM Time Series Extraction Script Documentation

## Overview

**Script Name:** `e3sm_extract_time_series_h0.py`

**Purpose:** Extracts time series data from E3SM (Energy Exascale Earth System Model) NetCDF h0 output files, aggregates data spatially (global or regional), processes variables (unit conversions, derived quantities), and outputs to formatted text files. The script handles both EAM (atmosphere) and ELM (land) model outputs.

**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Requirements](#installation--requirements)
3. [Basic Usage](#basic-usage)
4. [E3SM Model Background](#e3sm-model-background)
5. [Complete Parameter Reference Table](#complete-parameter-reference-table)
6. [Detailed Parameter Descriptions](#detailed-parameter-descriptions)
7. [Spatial Aggregation Types](#spatial-aggregation-types)
8. [Variable Processing](#variable-processing)
9. [Regional Analysis](#regional-analysis)
10. [JSON Configuration Examples](#json-configuration-examples)
11. [Output Files](#output-files)
12. [Troubleshooting](#troubleshooting)
13. [Best Practices](#best-practices)

---

## Installation & Requirements

### Required Python Libraries

```bash
pip install pandas numpy xarray netCDF4 multiprocessing
```

### Required Utility Modules

The script imports several utility modules that must be in the same directory or Python path:
- `utility_constants` - Physical constants (conversion factors, molecular masses)
- `utility_dataframes` - DataFrame manipulation functions
- `utility_functions` - General utility functions
- `utility_e3sm_netcdf` - E3SM-specific NetCDF handling functions

### System Requirements

- Python 3.7+
- Multi-core processor (script uses parallel processing)
- Sufficient RAM for NetCDF file processing (8+ GB recommended)
- Access to E3SM simulation output directories

---

## Basic Usage

### Command Line Execution

```bash
python e3sm_extract_time_series_h0.py path/to/config.json
```

**Multiple Configuration Files:**
```bash
python e3sm_extract_time_series_h0.py config1.json config2.json config3.json
```

### What the Script Does

For each configuration block in the JSON file, the script:
1. **Locates** E3SM NetCDF h0 files in simulation directory
2. **Filters** files by year range and file type (ELM/EAM)
3. **Extracts** specified variables from each monthly NetCDF file
4. **Aggregates** spatially (area-weighted mean/sum over lat/lon)
5. **Filters** by region (optional - global, Amazon, CONUS, etc.)
6. **Processes** variables (unit conversions, derived quantities)
7. **Combines** data across months and years
8. **Outputs** time series to formatted text file
9. **Uses** parallel processing for fast execution

---

## E3SM Model Background

### E3SM Components

**E3SM (Energy Exascale Earth System Model)** consists of multiple component models:

- **EAM (E3SM Atmosphere Model)** - Atmospheric processes, radiation, clouds
- **ELM (E3SM Land Model)** - Land surface, vegetation, carbon cycle
- **MPAS-Ocean** - Ocean circulation
- **MPAS-Seaice** - Sea ice dynamics
- **MOSART** - River transport

This script processes outputs from **EAM** and **ELM** components.

### NetCDF h0 Files

**File Naming Convention:**
```
[component].h0.YYYY-MM.nc

Examples:
elm.h0.2015-01.nc    # ELM output for January 2015
eam.h0.2015-01.nc    # EAM output for January 2015
```

**File Contents:**
- One file per month per component
- Gridded data (latitude × longitude)
- Multiple variables in each file
- Metadata (units, descriptions)

### Common Variables

**ELM Variables (Land):**
- `GPP` - Gross Primary Production (gC/m²/s)
- `NPP` - Net Primary Production
- `NBP` - Net Biome Production
- `TOTECOSYSC` - Total ecosystem carbon (gC/m²)
- `TOTVEGC` - Total vegetation carbon
- `TOTSOMC` - Total soil carbon
- `PCT_LANDUNIT` - Land unit percentages
- `PCT_NAT_PFT` - Plant functional type percentages

**EAM Variables (Atmosphere):**
- `PRECC`, `PRECL` - Convective and large-scale precipitation (m/s)
- `SFCO2` - Surface CO₂ flux (kg/m²/s)
- `TMCO2` - Atmospheric CO₂ mass (kg)
- `TREFHT` - Reference height temperature (K)

---

## Complete Parameter Reference Table

### Required Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `simulation_path` | string | **Yes** | Path to E3SM simulation output directory |
| `output_file` | string | **Yes** | Path/name for output time series file |
| `netcdf_substrings` | nested list | **Yes** | File type identifiers for NetCDF files |
| `variables` | nested list | **Yes** | Variables to extract from each file type |

### Optional Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `lat_lon_aggregation_types` | list | No | `["area_weighted_mean_or_sum"]` | `"area_weighted_mean_or_sum"`, `"mean"`, `"sum"` | Spatial aggregation method |
| `regions` | list | No | `[None]` (global) | Region names or `None` | Geographic regions to extract |
| `process_variables` | boolean | No | `true` | `true`, `false` | Apply unit conversions and derived variables |
| `start_year` | integer | No | `2015` | Any year | First year to extract |
| `end_year` | integer | No | `2100` | Any year | Last year to extract |
| `write_to_csv` | boolean | No | `false` | `true`, `false` | Output as CSV instead of fixed-width |

---

## Detailed Parameter Descriptions

### Core Required Parameters

#### `simulation_path`
**Type:** String (directory path)  
**Required:** Yes  
**Description:** Complete path to directory containing E3SM NetCDF output files.

**Examples:**
```json
"simulation_path": "/lcrc/group/e3sm/ac.eva.sinha/20240730_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3_without_feedbacks/run"
```

**Requirements:**
- Must be valid directory path
- Must contain NetCDF (.nc) files
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

#### `output_file`
**Type:** String (file path)  
**Required:** Yes  
**Description:** Path and filename for output time series data.

**Examples:**
```json
"output_file": "./../2025_DiVittorio_et_al_e3sm/control_time_series.dat"
"output_file": "./output/amazon_time_series.csv"
```

**File Formats:**
- **Fixed-width (.dat)** - Default, human-readable, aligned columns
- **CSV (.csv)** - If filename ends in .csv or `write_to_csv: true`

**Output Structure:**
```
Year  Month  GPP (PgC/month)  NPP (PgC/month)  TREFHT (K)  ...
2015      1        10.234           5.123      288.45
2015      2         9.876           4.987      289.12
...
```

---

#### `netcdf_substrings`
**Type:** Nested list of strings  
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

**Multiple substrings (must all be present):**
```json
"netcdf_substrings": [["elm", "h0", "clm2"]]
```

**How It Works:**
- Script searches `simulation_path` for .nc files
- Each file must contain ALL substrings in its corresponding list
- Files grouped by type based on which substring list they match

---

#### `variables`
**Type:** Nested list of strings  
**Required:** Yes  
**Description:** Variables to extract from each NetCDF file type.

**Format:** List of lists, where each inner list corresponds to one file type in `netcdf_substrings`.

**Structure:**
```
variables[i] = variables to extract from netcdf_substrings[i] files
```

**Example:**
```json
{
    "netcdf_substrings": [["elm.h0"], ["eam.h0"]],
    "variables": [
        ["GPP", "NPP", "NBP", "TOTECOSYSC"],
        ["PRECC", "PRECL", "TREFHT"]
    ]
}
```

**Breakdown:**
- `["GPP", "NPP", "NBP", "TOTECOSYSC"]` extracted from ELM files
- `["PRECC", "PRECL", "TREFHT"]` extracted from EAM files
- Final output has columns for ALL variables (merged by Year/Month)

**Common ELM Variables:**
```json
"variables": [
    ["GPP", "NPP", "NBP", "ER", "HR", "NEE", "NEP",
     "TOTECOSYSC", "TOTVEGC", "TOTSOMC", "TOTLITC",
     "WOODC", "LAND_UPTAKE", "LAND_USE_FLUX",
     "PCT_LANDUNIT", "PCT_NAT_PFT"]
]
```

**Common EAM Variables:**
```json
"variables": [
    ["PRECC", "PRECL", "PRECSC", "PRECSL",
     "SFCO2", "SFCO2_LND", "SFCO2_OCN",
     "TMCO2", "TMCO2_LND", "TMCO2_OCN",
     "TREFHT", "PBOT", "PCO2", "QBOT"]
]
```

---

### Optional Parameters

#### `lat_lon_aggregation_types`
**Type:** List of strings  
**Required:** No  
**Default:** `["area_weighted_mean_or_sum"]` for all file types  
**Possible Values:** `"area_weighted_mean_or_sum"`, `"mean"`, `"sum"`

**Description:** How to aggregate data spatially across lat/lon grid.

**Format:** One aggregation type per file type in `netcdf_substrings`.

**Example:**
```json
{
    "netcdf_substrings": [["elm.h0"], ["eam.h0"]],
    "lat_lon_aggregation_types": ["area_weighted_mean_or_sum", "area_weighted_mean_or_sum"]
}
```

**Options:**

1. **`"area_weighted_mean_or_sum"`** (recommended, default)
   - Intelligent aggregation based on variable type
   - **Fluxes/Stocks** (per area units): Area-weighted sum → global total
   - **Intensive variables** (temperature, pressure): Area-weighted mean
   - Accounts for gridcell area variations
   - Accounts for land/ocean fractions

2. **`"mean"`**
   - Simple arithmetic mean across all gridcells
   - Ignores gridcell area differences
   - Use with caution (biases towards high-latitude cells)

3. **`"sum"`**
   - Simple sum across all gridcells
   - Appropriate only for specific use cases
   - Generally not recommended

---

#### `regions`
**Type:** List of strings or `None`  
**Required:** No  
**Default:** `[None]` (global - no filtering)  
**Possible Values:** Region names (see below) or `None`

**Description:** Geographic regions to filter data extraction.

**Format:** One region per file type in `netcdf_substrings`.

**Example:**
```json
{
    "netcdf_substrings": [["elm.h0"], ["eam.h0"]],
    "regions": ["amazon", "amazon"]
}
```

**Available Regions:**

**North America:**
- `"noam"` - North America
- `"bona"` - Boreal North America
- `"tena"` - Temperate North America
- `"conus"` - Continental United States
- `"columbia"` - Columbia River watershed
- `"ceam"` - Central America

**South America:**
- `"soam"` - South America
- `"nhsa"` - Northern Hemisphere South America
- `"shsa"` - Southern Hemisphere South America
- `"amazon"` - Amazon Basin

**Europe/Middle East:**
- `"euro"` - Europe
- `"mide"` - Middle East

**Africa:**
- `"afrc"` - Africa
- `"nhaf"` - Northern Hemisphere Africa
- `"shaf"` - Southern Hemisphere Africa

**Asia:**
- `"asia"` - Asia
- `"boas"` - Boreal Asia
- `"ceas"` - Central Asia
- `"seas"` - Southeast Asia
- `"eqas"` - Equatorial Asia

**Australia:**
- `"aust"` - Australia

**Regional Boundaries:**
Defined in `utility_e3sm_netcdf.py` as lat/lon bounding boxes.

**Example:**
```python
amazon = [-85°, -35°] longitude, [-25°, 15°] latitude
conus = [-125.25°, -66.25°] longitude, [23.25°, 54.75°] latitude
```

---

#### `process_variables`
**Type:** Boolean  
**Required:** No  
**Default:** `true`  
**Possible Values:** `true`, `false`

**Description:** Whether to apply post-processing to extracted variables.

**When `true` (recommended):**
- Unit conversions (g → Pg, m/s → mm/month, /s → /month)
- Derived variables (total precipitation, CO₂ mole fraction)
- Improved readability and scientific usability

**When `false`:**
- Raw values from NetCDF files
- Original units
- Minimal processing
- Use for debugging or custom processing

**Example:**
```json
"process_variables": true
```

**Transformations Applied:**
See [Variable Processing](#variable-processing) section.

---

#### `start_year` and `end_year`
**Type:** Integer  
**Required:** No  
**Defaults:** `start_year: 2015`, `end_year: 2100`

**Description:** Year range for time series extraction.

**Example:**
```json
"start_year": 2020,
"end_year": 2050
```

**Usage:**
- Extracts only NetCDF files within this range
- Reduces processing time for subset analysis
- Useful for specific time periods (e.g., 21st century only)

**File Selection:**
```
elm.h0.2019-12.nc  → Excluded (year 2019 < start_year 2020)
elm.h0.2020-01.nc  → Included
elm.h0.2050-12.nc  → Included
elm.h0.2051-01.nc  → Excluded (year 2051 > end_year 2050)
```

---

#### `write_to_csv`
**Type:** Boolean  
**Required:** No  
**Default:** `false`

**Description:** Force output to CSV format.

**Behavior:**
- `false` → Fixed-width format (.dat)
- `true` → CSV format
- Automatically `true` if `output_file` ends in `.csv`

**Example:**
```json
"write_to_csv": true
```

**Comparison:**

**Fixed-width (.dat):**
```
Year  Month  GPP (PgC/month)  NPP (PgC/month)
2015      1        10.234           5.123
2015      2         9.876           4.987
```

**CSV (.csv):**
```
Year,Month,GPP (PgC/month),NPP (PgC/month)
2015,1,10.234,5.123
2015,2,9.876,4.987
```

---

## Spatial Aggregation Types

### Area-Weighted Mean or Sum (Recommended)

**Name:** `"area_weighted_mean_or_sum"`

**How It Works:**

1. **Multiply by gridcell area:**
   - Each gridcell has different area (latitude-dependent)
   - Multiply all values by their gridcell area

2. **Account for land/ocean fractions:**
   - EAM variables with `_LND` suffix → multiply by land fraction
   - EAM variables with `_OCN` suffix → multiply by ocean fraction

3. **Sum across all gridcells:**
   - Create global totals (area-weighted sums)

4. **Decide mean vs sum based on units:**
   - **Fluxes/Stocks** (with `/m²` or `/m2`): Keep as sum → global total
   - **Intensive variables** (temperature, pressure): Divide by total area → mean

**Examples:**

**Flux Variable (GPP in gC/m²/s):**
```
Original units: gC/m²/s
After aggregation: gC/s (global total)
After processing: PgC/month (global total)
```

**Intensive Variable (TREFHT in K):**
```
Original units: K
After aggregation: K (global mean temperature)
Stays as: K
```

**Formula for area-weighted mean:**
```
mean = Σ(value × area) / Σ(area)
```

---

### Simple Mean

**Name:** `"mean"`

**How It Works:**
- Arithmetic mean across all gridcells
- Ignores gridcell area differences
- `mean = Σ(values) / n`

**Use Cases:**
- Quick estimates (non-rigorous)
- Debugging

**Limitations:**
- Biases towards high-latitude gridcells (smaller area but equal weight)
- Not appropriate for fluxes/stocks
- Not scientifically rigorous for global estimates

---

### Simple Sum

**Name:** `"sum"`

**How It Works:**
- Simple sum across all gridcells
- `sum = Σ(values)`
- Ignores area

**Use Cases:**
- Very specific applications
- Generally not recommended

**Limitations:**
- Meaningless for intensive variables
- Incorrect for fluxes (doesn't account for area)

---

## Variable Processing

When `process_variables: true`, the script applies multiple transformations:

### 1. Precipitation Unit Conversion

**Variables:** `PRECC`, `PRECL`, `PRECSC`, `PRECSL`

**Transformation:**
```
Original: m/s
Convert to: mm/month
Factor: seconds_in_month × 1000 mm/m
```

**Derived Variable:**
```
PRECIP (mm/month) = PRECC + PRECL + PRECSC + PRECSL
```

**Result:**
- Individual precipitation components in mm/month
- Total precipitation as new column

---

### 2. CO₂ Mole Fraction Calculation

**Variables Required:** `PBOT`, `PCO2`, `QBOT`

**Calculation:**
```python
# Partial pressure of water vapor
pH2O = PBOT × QBOT / (0.622 + 0.378 × QBOT)

# Dry air mole fraction (ppm)
ZCO2 (ppm) = 1e6 × PCO2 / (PBOT - pH2O)
```

**Derived Variable:**
```
ZCO2 (ppm) - Atmospheric CO₂ mole fraction in dry air
```

**Scientific Background:**
- CO₂ concentration typically reported as mole fraction in *dry* air
- Must remove water vapor contribution
- Uses humidity formula from psychrometric calculations

---

### 3. Carbon Mass Unit Conversion

**Variables:** All with units containing `gC`, `g/`, or `kg`

**Transformation:**
```
g → Pg (1 Pg = 10¹⁵ g)
kg → Pg (1 Pg = 10¹² kg)
```

**Examples:**
```
Original: gC/m²/s
After: PgC/m²/s

Original: gC/m²
After: PgC/m²

Original: kg/m²/s
After: Pg/m²/s
```

**Rationale:**
- Pg (Petagrams) = 10¹⁵ grams
- Appropriate scale for global carbon budgets
- Standard in climate science literature

---

### 4. Temporal Unit Conversion

**Variables:** All with units containing `/s`

**Transformation:**
```
per second → per month
Factor: seconds_in_month (varies by month)
```

**Account for varying month lengths:**
```python
seconds_in_month = days_in_month × 86400
# January: 31 days = 2,678,400 seconds
# February: 28/29 days = 2,419,200 / 2,505,600 seconds
# etc.
```

**Examples:**
```
Original: PgC/m²/s
After: PgC/m²/month

Original: kg/m²/s
After: Pg/m²/month
```

---

### 5. CO₂ to C Conversion

**Variables:** `SFCO2`, `TMCO2` (and variants)

**Transformation:**
```
PgCO₂ → PgC
Factor: MM_C / MM_CO₂ = 12 / 44 = 0.2727
```

**Examples:**
```
SFCO2: Surface CO₂ flux
Original: Pg/m²/month (as CO₂)
After: PgC/m²/month (as C)

TMCO2: Total atmospheric CO₂ mass
Original: Pg (as CO₂)
After: PgC (as C)
```

**Rationale:**
- Carbon budgets reported in C, not CO₂
- Molecular mass conversion: CO₂ (44 g/mol) → C (12 g/mol)

---

### 6. Plant Functional Type (PFT) Analysis

**Variables:** `PCT_LANDUNIT`, `PCT_NAT_PFT`

**Processing:**

**Land Unit Extraction:**
```
PCT_LANDUNIT → 9 land units
Select: Land unit 0 (vegetation)
Convert: Percent → Fraction (/100)
Output: FRAC_VEG
```

**PFT Areas:**
```
17 PFTs in E3SM:
0: Bare soil
1-8: Trees (8 types)
9-11: Shrubs (3 types)
12-14: Grasses (3 types)
15: Crop
16: Empty (ignored)

Aggregate into:
- Individual PFTs: PFT_1_AREA through PFT_16_AREA (km²)
- Bare soil: BARE_AREA (km²)
- Forest: FOREST_AREA (km²) [sum of PFTs 1-8]
- Shrub: SHRUB_AREA (km²) [sum of PFTs 9-11]
- Grass: GRASS_AREA (km²) [sum of PFTs 12-14]
- Crop: CROP_AREA (km²) [PFT 15]
```

**Calculation:**
```
PFT_AREA (km²) = AREA (km²) × FRAC_VEG × PFT_fraction
```

**Output Columns:**
```
FRAC_VEG
AREA (km^2)
PFT_1_AREA (km^2)
...
PFT_16_AREA (km^2)
BARE_AREA (km^2)
FOREST_AREA (km^2)
SHRUB_AREA (km^2)
GRASS_AREA (km^2)
CROP_AREA (km^2)
```

---

## Regional Analysis

### Defining Regions

Regions defined by lat/lon bounding boxes in `utility_e3sm_netcdf.py`:

```python
def get_regional_bounds(region):
    if region == 'amazon':
        bounds = [-85., -35., -25., 15.]  # [min_lon, max_lon, min_lat, max_lat]
    elif region == 'conus':
        bounds = [-125.25, -66.25, 23.25, 54.75]
    # ... etc
```

### How Regional Filtering Works

1. **Before spatial aggregation:**
   - Load NetCDF file
   - Filter gridcells outside bounding box
   - Keep only gridcells within region

2. **Apply aggregation:**
   - Area-weighted mean/sum over regional gridcells only
   - Land/ocean fractions applied within region

3. **Result:**
   - Regional totals or means instead of global

### Example: Amazon Carbon Fluxes

```json
{
    "simulation_path": "/path/to/simulation",
    "output_file": "./amazon_carbon.dat",
    "netcdf_substrings": [["elm.h0"]],
    "variables": [["GPP", "NPP", "NBP", "ER"]],
    "regions": ["amazon"],
    "start_year": 2015,
    "end_year": 2100
}
```

**Output:**
- GPP, NPP, NBP, ER for Amazon Basin only
- Units: PgC/month (regional totals)

---

## JSON Configuration Examples

### Example 1: Global Carbon Cycle (ELM + EAM)

**Purpose:** Extract comprehensive carbon cycle variables globally

```json
{
    "simulation_path": "/lcrc/group/e3sm/simulation_directory/run",
    "output_file": "./global_carbon_time_series.dat",
    "netcdf_substrings": [["elm.h0"], ["eam.h0"]],
    "variables": [
        ["GPP", "NPP", "NBP", "ER", "HR", "TOTECOSYSC", "TOTVEGC", "TOTSOMC"],
        ["SFCO2", "SFCO2_LND", "SFCO2_OCN", "TMCO2"]
    ],
    "lat_lon_aggregation_types": ["area_weighted_mean_or_sum", "area_weighted_mean_or_sum"],
    "process_variables": true,
    "start_year": 2015,
    "end_year": 2100
}
```

**What it does:**
- Extracts ELM land carbon fluxes and stocks
- Extracts EAM atmospheric CO₂ fluxes and mass
- Global aggregation with area weighting
- Processes units (g→Pg, /s→/month)
- Output: 86 years × 12 months = 1,032 rows

---

### Example 2: Regional Analysis (Amazon)

**Purpose:** Focus on Amazon Basin carbon dynamics

```json
{
    "simulation_path": "/lcrc/group/e3sm/simulation_directory/run",
    "output_file": "./amazon_time_series.dat",
    "netcdf_substrings": [["elm.h0"], ["eam.h0"]],
    "variables": [
        ["GPP", "NPP", "NBP", "NEE", "TOTECOSYSC"],
        ["PRECC", "PRECL", "TREFHT"]
    ],
    "lat_lon_aggregation_types": ["area_weighted_mean_or_sum", "area_weighted_mean_or_sum"],
    "regions": ["amazon", "amazon"],
    "process_variables": true,
    "start_year": 2015,
    "end_year": 2100
}
```

**What it does:**
- Filters to Amazon Basin coordinates
- Extracts carbon fluxes and climate variables
- Regional totals/means (not global)
- Total precipitation derived from PRECC + PRECL

---

### Example 3: Climate Variables Only (EAM)

**Purpose:** Extract atmospheric/climate data without land carbon

```json
{
    "simulation_path": "/lcrc/group/e3sm/simulation_directory/run",
    "output_file": "./climate_time_series.dat",
    "netcdf_substrings": [["eam.h0"]],
    "variables": [
        ["PRECC", "PRECL", "PRECSC", "PRECSL", "TREFHT", "PBOT", "PCO2", "QBOT"]
    ],
    "lat_lon_aggregation_types": ["area_weighted_mean_or_sum"],
    "process_variables": true,
    "start_year": 2015,
    "end_year": 2100
}
```

**What it does:**
- EAM variables only
- Calculates PRECIP (total precipitation in mm/month)
- Calculates ZCO2 (CO₂ mole fraction in ppm)
- Global mean temperature (TREFHT)

---

### Example 4: Land Use and Vegetation

**Purpose:** Track land use changes and vegetation types

```json
{
    "simulation_path": "/lcrc/group/e3sm/simulation_directory/run",
    "output_file": "./landuse_time_series.dat",
    "netcdf_substrings": [["elm.h0"]],
    "variables": [
        ["LAND_USE_FLUX", "WOOD_HARVESTC", "DWT_CONV_CFLUX_GRC",
         "PCT_LANDUNIT", "PCT_NAT_PFT"]
    ],
    "lat_lon_aggregation_types": ["area_weighted_mean_or_sum"],
    "process_variables": true,
    "start_year": 2015,
    "end_year": 2100
}
```

**What it does:**
- Land use change fluxes
- Wood harvest carbon
- Vegetation fractions
- PFT areas (individual and aggregated)
- Output includes: FRAC_VEG, FOREST_AREA, CROP_AREA, etc.

---

### Example 5: Short Time Period for Testing

**Purpose:** Quick test with subset of data

```json
{
    "simulation_path": "/lcrc/group/e3sm/simulation_directory/run",
    "output_file": "./test_time_series.dat",
    "netcdf_substrings": [["elm.h0"]],
    "variables": [["GPP", "NPP", "NBP"]],
    "lat_lon_aggregation_types": ["area_weighted_mean_or_sum"],
    "process_variables": true,
    "start_year": 2015,
    "end_year": 2020
}
```

**What it does:**
- Only 6 years of data
- Fast execution for testing
- Verify configuration before full run

---

### Example 6: CSV Output

**Purpose:** Output in CSV format for easy import to other tools

```json
{
    "simulation_path": "/lcrc/group/e3sm/simulation_directory/run",
    "output_file": "./carbon_time_series.csv",
    "netcdf_substrings": [["elm.h0"]],
    "variables": [["GPP", "NPP", "NBP", "ER"]],
    "lat_lon_aggregation_types": ["area_weighted_mean_or_sum"],
    "process_variables": true,
    "start_year": 2015,
    "end_year": 2100,
    "write_to_csv": true
}
```

**What it does:**
- Same as Example 1 but CSV output
- Comma-separated values
- Easy import to Excel, R, Python

---

### Example 7: Raw Variables (No Processing)

**Purpose:** Extract raw NetCDF values without unit conversions

```json
{
    "simulation_path": "/lcrc/group/e3sm/simulation_directory/run",
    "output_file": "./raw_data.dat",
    "netcdf_substrings": [["elm.h0"]],
    "variables": [["GPP", "NPP", "ER"]],
    "lat_lon_aggregation_types": ["area_weighted_mean_or_sum"],
    "process_variables": false,
    "start_year": 2015,
    "end_year": 2100
}
```

**What it does:**
- No unit conversions
- Original units: gC/m²/s, etc.
- Use for custom processing or debugging

---

## Output Files

### Fixed-Width Format (.dat)

**Default format** - Human-readable with aligned columns

**Example:**
```
Year  Month  GPP (PgC/month)  NPP (PgC/month)  NBP (PgC/month)  TREFHT (K)
2015      1       10.234           5.123           -0.456      288.45
2015      2        9.876           4.987           -0.421      289.12
2015      3       11.345           5.678           -0.534      290.23
```

**Characteristics:**
- Spaces used for padding
- Columns aligned
- Easy to read
- Standard for E3SM analysis

---

### CSV Format (.csv)

**Optional format** - Machine-readable

**Example:**
```
Year,Month,GPP (PgC/month),NPP (PgC/month),NBP (PgC/month),TREFHT (K)
2015,1,10.234,5.123,-0.456,288.45
2015,2,9.876,4.987,-0.421,289.12
2015,3,11.345,5.678,-0.534,290.23
```

**Characteristics:**
- Comma-separated
- Smaller file size
- Easy import to Excel, pandas, R

---

### Output Structure

**Temporal Organization:**
- One row per month
- Sorted by Year, then Month
- Continuous monthly time series

**Column Organization:**
```
[Year] [Month] [ELM variables...] [EAM variables...]
```

**Column Headers:**
```
Variable_name (units)

Examples:
GPP (PgC/month)
TREFHT (K)
PRECIP (mm/month)
ZCO2 (ppm)
```

---

## Troubleshooting

### Issue 1: No NetCDF Files Found

**Error Message:**
```
No files found matching criteria
```

**Causes:**
- Incorrect `simulation_path`
- Wrong `netcdf_substrings`
- Files don't match year range

**Solutions:**
```bash
# Check directory exists
ls -la /path/to/simulation/

# Check for NetCDF files
ls /path/to/simulation/*.nc | head

# Verify filenames match substrings
ls /path/to/simulation/elm.h0.*.nc
ls /path/to/simulation/eam.h0.*.nc

# Check year range
ls /path/to/simulation/elm.h0.2015-*.nc
```

---

### Issue 2: Variable Not Found in NetCDF

**Error Message:**
```
KeyError: 'VARIABLE_NAME'
```

**Causes:**
- Variable not in NetCDF file
- Typo in variable name
- Wrong file type (ELM vs EAM)

**Solutions:**
```python
# Check what variables are in file
import xarray as xr
ds = xr.open_dataset('elm.h0.2015-01.nc')
print(list(ds.data_vars))

# Common mistakes:
# "GPP" is in ELM, not EAM
# "TREFHT" is in EAM, not ELM
# Variable names are case-sensitive
```

---

### Issue 3: Memory Error

**Error Message:**
```
MemoryError
```

**Causes:**
- Large NetCDF files
- Many variables
- Long time series
- Insufficient RAM

**Solutions:**
1. **Reduce year range:**
```json
"start_year": 2015,
"end_year": 2050  // Instead of 2100
```

2. **Process fewer variables:**
```json
"variables": [["GPP", "NPP"]]  // Instead of 20+ variables
```

3. **Run on larger memory machine**

4. **Process in chunks:**
```json
// Config 1
{"start_year": 2015, "end_year": 2050}

// Config 2
{"start_year": 2051, "end_year": 2100}
```

---

### Issue 4: Slow Execution

**Symptom:** Script taking hours to complete

**Causes:**
- Many NetCDF files (86 years × 12 months × 2 file types = 2,064 files)
- Large files
- Limited CPU cores

**Solutions:**

1. **Use multiprocessing (already built-in):**
Script automatically uses all CPU cores

2. **Reduce year range:**
```json
"start_year": 2080,
"end_year": 2100  // Only 21 years instead of 86
```

3. **Run on high-performance computing (HPC)**

4. **Monitor progress:**
Script prints completion time for each output file

---

### Issue 5: Regional Bounds Error

**Error Message:**
```
Did not recognize the selected region!
```

**Cause:**
- Region name typo
- Region not in predefined list

**Solutions:**
```json
// Check spelling
"regions": ["amazon"]  // Correct
"regions": ["Amazon"]  // Wrong (case-sensitive)

// Use None for global
"regions": [null]

// Check available regions in utility_e3sm_netcdf.py
```

---

### Issue 6: Mismatched List Lengths

**Error Message:**
```
IndexError: list index out of range
```

**Cause:**
- `netcdf_substrings`, `variables`, `lat_lon_aggregation_types`, `regions` have different lengths

**Solution:**
```json
// All lists must have same length
{
    "netcdf_substrings": [["elm.h0"], ["eam.h0"]],      // Length 2
    "variables": [["GPP"], ["TREFHT"]],                 // Length 2
    "lat_lon_aggregation_types": ["area_weighted_mean_or_sum", "area_weighted_mean_or_sum"],  // Length 2
    "regions": ["amazon", "amazon"]                     // Length 2
}

// Or use defaults (don't specify)
{
    "netcdf_substrings": [["elm.h0"], ["eam.h0"]],
    "variables": [["GPP"], ["TREFHT"]]
    // lat_lon_aggregation_types and regions will default to length 2
}
```

---

## Best Practices

### 1. Directory Organization

**Recommended Structure:**
```
project/
├── simulations/
│   ├── control/
│   │   └── run/          # simulation_path
│   ├── scenario1/
│   │   └── run/
│   └── scenario2/
│       └── run/
├── output/
│   ├── time_series/      # output files
│   └── processed/
├── configs/
│   └── extraction_configs.json
└── analysis/
    └── scripts/
```

---

### 2. Configuration File Organization

**Use descriptive names:**
```
configs/
├── global_carbon_extraction.json
├── regional_amazon_extraction.json
├── climate_variables_extraction.json
└── testing_2015_2020.json
```

**Group related extractions:**
```json
[
    {"output_file": "control_global.dat", ...},
    {"output_file": "control_amazon.dat", ...},
    {"output_file": "scenario_global.dat", ...},
    {"output_file": "scenario_amazon.dat", ...}
]
```

---

### 3. Testing Strategy

**Start small:**
```json
{
    "start_year": 2015,
    "end_year": 2016,  // Just 2 years
    "variables": [["GPP", "NPP"]]  // Few variables
}
```

**Then scale up:**
```json
{
    "start_year": 2015,
    "end_year": 2100,  // Full range
    "variables": [["GPP", "NPP", "NBP", ...]]  // All variables
}
```

---

### 4. Variable Selection

**Consider your analysis needs:**

**Carbon budget analysis:**
```json
"variables": [["GPP", "NPP", "NBP", "ER", "HR", "NEE"]]
```

**Climate impact:**
```json
"variables": [["TREFHT", "PRECC", "PRECL"]]
```

**Land use:**
```json
"variables": [["LAND_USE_FLUX", "WOOD_HARVESTC", "PCT_LANDUNIT", "PCT_NAT_PFT"]]
```

**Comprehensive:**
```json
"variables": [[
    "GPP", "NPP", "NBP", "ER", "HR", "NEE",
    "TOTECOSYSC", "TOTVEGC", "TOTSOMC",
    "LAND_USE_FLUX", "PCT_LANDUNIT", "PCT_NAT_PFT"
]]
```

---

### 5. Regional vs Global

**When to use global:**
- Earth system scale analysis
- Comparing total carbon budgets
- Climate projections

**When to use regional:**
- Regional impacts (Amazon, CONUS, etc.)
- Hotspot analysis
- Policy-relevant scales

**Multiple regions:**
```json
[
    {"output_file": "global.dat", "regions": [null]},
    {"output_file": "amazon.dat", "regions": ["amazon"]},
    {"output_file": "conus.dat", "regions": ["conus"]}
]
```

---

### 6. Performance Optimization

**Parallel processing:**
- Script uses all CPU cores automatically
- Process multiple output files in one run

**Efficient configuration:**
```json
// Good: One JSON with multiple outputs
[
    {"output_file": "control.dat", ...},
    {"output_file": "scenario.dat", ...}
]

// Less efficient: Multiple script calls
python script.py config1.json
python script.py config2.json
```

---

### 7. Data Validation

**After extraction, verify:**

```python
import pandas as pd

df = pd.read_csv('output.csv')  # or read_fwf for .dat

# Check dimensions
print(f"Shape: {df.shape}")
# Expected: 1032 rows (86 years × 12 months) for 2015-2100

# Check for missing values
print(df.isnull().sum())

# Check value ranges
print(df.describe())

# Verify temporal continuity
print(df[['Year', 'Month']].head(15))
```

---

### 8. Reproducibility

**Document your extraction:**

```json
{
    "comment": "Global carbon cycle for control simulation",
    "date_extracted": "2026-01-24",
    "simulation": "SSP245_ZATM_BGC_without_feedbacks",
    "simulation_path": "/lcrc/group/e3sm/.../run",
    ...
}
```

**Version control configurations:**
```bash
git add configs/
git commit -m "Add extraction configs for control simulation"
```

---

## Integration with Other Scripts

### Typical E3SM-GCAM Analysis Workflow

```
1. Run E3SM simulation
   ↓ (Produces NetCDF h0 files)
   
2. e3sm_extract_time_series_h0.py  ← THIS SCRIPT
   ↓ (Extracts time series from NetCDF)
   
3. Further analysis scripts:
   - Time series plotting
   - Statistical analysis
   - Model comparison
   - GCAM coupling validation
```

### Downstream Analysis

**Time Series Plotting:**
```python
# Use extracted data in plotting scripts
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('carbon_time_series.csv')
plt.plot(df['Year'], df['GPP (PgC/month)'])
plt.xlabel('Year')
plt.ylabel('GPP (PgC/month)')
plt.show()
```

**Statistical Analysis:**
```python
# Analyze trends
from scipy import stats

df = pd.read_csv('carbon_time_series.csv')
slope, intercept, r_value, p_value, std_err = stats.linregress(
    df['Year'], df['NBP (PgC/month)']
)
print(f"NBP trend: {slope:.4f} PgC/month/year, p={p_value:.4f}")
```

---

## Appendix: Common E3SM Variables

### ELM Carbon Cycle Variables

| Variable | Units | Description |
|----------|-------|-------------|
| `GPP` | gC/m²/s | Gross Primary Production |
| `NPP` | gC/m²/s | Net Primary Production |
| `NBP` | gC/m²/s | Net Biome Production |
| `NEE` | gC/m²/s | Net Ecosystem Exchange |
| `NEP` | gC/m²/s | Net Ecosystem Production |
| `ER` | gC/m²/s | Ecosystem Respiration |
| `HR` | gC/m²/s | Heterotrophic Respiration |
| `LAND_UPTAKE` | gC/m²/s | Net land carbon uptake |
| `LAND_USE_FLUX` | gC/m²/s | Land use change emissions |
| `WOOD_HARVESTC` | gC/m²/s | Wood harvest carbon |
| `PFT_FIRE_CLOSS` | gC/m²/s | Fire carbon loss |

### ELM Carbon Stock Variables

| Variable | Units | Description |
|----------|-------|-------------|
| `TOTECOSYSC` | gC/m² | Total ecosystem carbon |
| `TOTVEGC` | gC/m² | Total vegetation carbon |
| `TOTSOMC` | gC/m² | Total soil organic carbon |
| `TOTLITC` | gC/m² | Total litter carbon |
| `WOODC` | gC/m² | Wood carbon |

### EAM Atmospheric Variables

| Variable | Units | Description |
|----------|-------|-------------|
| `TREFHT` | K | Reference height temperature (2m) |
| `PBOT` | Pa | Surface pressure |
| `PCO2` | Pa | CO₂ partial pressure |
| `QBOT` | kg/kg | Specific humidity at surface |
| `PRECC` | m/s | Convective precipitation rate |
| `PRECL` | m/s | Large-scale precipitation rate |
| `PRECSC` | m/s | Convective snow rate |
| `PRECSL` | m/s | Large-scale snow rate |

### EAM CO₂ Flux Variables

| Variable | Units | Description |
|----------|-------|-------------|
| `SFCO2` | kg/m²/s | Total surface CO₂ flux |
| `SFCO2_LND` | kg/m²/s | Land surface CO₂ flux |
| `SFCO2_OCN` | kg/m²/s | Ocean surface CO₂ flux |
| `SFCO2_FFF` | kg/m²/s | Fossil fuel CO₂ flux |
| `TMCO2` | kg | Total atmospheric CO₂ mass |
| `TMCO2_LND` | kg | CO₂ mass from land |
| `TMCO2_OCN` | kg | CO₂ mass from ocean |
| `TMCO2_FFF` | kg | CO₂ mass from fossil fuels |

---

## References

- E3SM Documentation: [https://e3sm.org/](https://e3sm.org/)
- E3SM Model Description: [https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2018MS001350](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2018MS001350)
- DiVittorio et al. (2025). "E3SM-GCAM coupling methodology and applications." *Journal of Advances in Modeling Earth Systems*. [DOI: 10.1029/2024MS004806](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024MS004806)
- GitHub Repository: [https://github.com/philipmyint/e3sm--gcam_analysis](https://github.com/philipmyint/e3sm--gcam_analysis)
- xarray Documentation: [https://xarray.pydata.org/](https://xarray.pydata.org/)

---

## Contact

For questions or issues:
- **Philip Myint**: myint1@llnl.gov
- **Dalei Hao**: dalei.hao@pnnl.gov

---

## Version Information

**Script:** e3sm_extract_time_series_h0.py  
**Utility Modules:** utility_e3sm_netcdf.py, utility_constants.py, utility_dataframes.py, utility_functions.py  
**Documentation Version:** 1.0  
**Last Updated:** January 2026  
**Python Version:** 3.7+  
**xarray Version:** 0.16+

---

*This documentation provides comprehensive guidance for using the `e3sm_extract_time_series_h0.py` script to extract, aggregate, and process time series data from E3SM model output NetCDF files for climate and carbon cycle analysis.*
