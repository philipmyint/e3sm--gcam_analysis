# E3SM Surface Data Time Series Extraction Script Documentation

## Overview

**Script Name:** `e3sm_extract_time_series_surfdata_iesm_dyn.py`

**Purpose:** Extracts annual time series data from the E3SM Human Component (EHC) dynamically-generated land surface data file (`surfdata_iESM_dyn.nc`). This script specifically processes human-related land management activities including grazing and wood harvesting (various harvest types), aggregates data spatially, and outputs annual time series for global or regional analysis.

**Key Difference from h0 Script:** Unlike `e3sm_extract_time_series_h0.py` which processes monthly ELM/EAM output files, this script processes a **single NetCDF file** containing time-varying land use data that is updated dynamically during E3SM-GCAM coupled simulations.

**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Requirements](#installation--requirements)
3. [Basic Usage](#basic-usage)
4. [E3SM Human Component Background](#e3sm-human-component-background)
5. [Complete Parameter Reference Table](#complete-parameter-reference-table)
6. [Detailed Parameter Descriptions](#detailed-parameter-descriptions)
7. [Available Variables](#available-variables)
8. [Variable Processing](#variable-processing)
9. [Regional Analysis](#regional-analysis)
10. [JSON Configuration Examples](#json-configuration-examples)
11. [Output Files](#output-files)
12. [Troubleshooting](#troubleshooting)
13. [Best Practices](#best-practices)
14. [Comparison with h0 Extraction Script](#comparison-with-h0-extraction-script)

---

## Installation & Requirements

### Required Python Libraries

```bash
pip install pandas numpy xarray netCDF4 multiprocessing
```

### Required Utility Modules

The script imports several utility modules that must be in the same directory or Python path:
- `utility_constants` - Physical constants (conversion factors)
- `utility_dataframes` - DataFrame manipulation functions
- `utility_functions` - General utility functions
- `utility_e3sm_netcdf` - E3SM-specific NetCDF handling functions

### System Requirements

- Python 3.7+
- Multi-core processor (script uses parallel processing)
- Sufficient RAM for NetCDF file processing (4+ GB recommended)
- Access to E3SM simulation output directories

---

## Basic Usage

### Command Line Execution

```bash
python e3sm_extract_time_series_surfdata_iesm_dyn.py path/to/config.json
```

**Multiple Configuration Files:**
```bash
python e3sm_extract_time_series_surfdata_iesm_dyn.py config1.json config2.json config3.json
```

### What the Script Does

For each configuration block in the JSON file, the script:
1. **Locates** the `surfdata_iESM_dyn.nc` file in simulation directory
2. **Extracts** specified variables for each year in the range
3. **Processes** PFT (Plant Functional Type) data automatically
4. **Converts** grazing/harvest fractions to areas (km²)
5. **Aggregates** spatially (area-weighted sum over lat/lon)
6. **Filters** by region (optional - global, Amazon, CONUS, etc.)
7. **Calculates** total harvest from individual harvest types
8. **Outputs** annual time series to formatted text file
9. **Uses** parallel processing (one process per year)

---

## E3SM Human Component Background

### The E3SM Human Component (EHC)

The **E3SM Human Component (EHC)** is a module within E3SM that represents human activities affecting the land surface, particularly in E3SM-GCAM coupled simulations. The EHC:

- **Bridges** E3SM (Earth system model) and GCAM (integrated assessment model)
- **Updates** land use dynamically during simulation
- **Tracks** human land management activities
- **Generates** `surfdata_iESM_dyn.nc` file at runtime

### The surfdata_iESM_dyn.nc File

**Purpose:** Land surface data file dynamically updated during E3SM-GCAM coupled runs

**Key Characteristics:**
- **Single file** (not monthly like h0 files)
- **Multi-year** (contains data for all simulation years)
- **Annual data** (one timestep per year)
- **Human activities** (grazing, harvest)
- **PFT distributions** (plant functional types)

**File Location:**
```
simulation_path/
└── surfdata_iESM_dyn.nc  (single file)
```

**Typical Size:** 1-5 GB (depending on simulation length)

**Time Dimension:**
```
time = year values (e.g., 2015, 2016, ..., 2100)
```

### Human Land Management Variables

**Grazing:**
- `GRAZING` - Fraction of grassland grazed by livestock
- Units: Dimensionless fraction (0-1)
- Converted to area (km²) by script

**Wood Harvesting:**
- `HARVEST_VH1` - Very heavy harvest intensity 1
- `HARVEST_VH2` - Very heavy harvest intensity 2
- `HARVEST_SH1` - Secondary heavy harvest intensity 1
- `HARVEST_SH2` - Secondary heavy harvest intensity 2
- `HARVEST_SH3` - Secondary heavy harvest intensity 3
- Units: Dimensionless fractions (0-1)
- Converted to forest area (km²) by script

**Land Cover:**
- `PCT_NATVEG` - Percent natural vegetation
- `PCT_NAT_PFT` - Percent of each plant functional type
- Used to calculate PFT areas

---

## Complete Parameter Reference Table

### Required Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `simulation_path` | string | **Yes** | Path to E3SM simulation directory containing surfdata_iESM_dyn.nc |
| `output_file` | string | **Yes** | Path/name for output time series file |
| `variables` | list | **Yes** | Variables to extract from surfdata file |

### Optional Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `region` | string | No | `None` (global) | Region name or `None` | Geographic region to extract |
| `start_year` | integer | No | `2015` | Any year | First year to extract |
| `end_year` | integer | No | `2100` | Any year | Last year to extract |
| `write_to_csv` | boolean | No | `false` | `true`, `false` | Output as CSV instead of fixed-width |

---

## Detailed Parameter Descriptions

### Core Required Parameters

#### `simulation_path`
**Type:** String (directory path)  
**Required:** Yes  
**Description:** Complete path to directory containing the `surfdata_iESM_dyn.nc` file.

**Examples:**
```json
"simulation_path": "/lcrc/group/e3sm/ac.eva.sinha/20240730_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3_without_feedbacks/run"
```

**Requirements:**
- Must be valid directory path
- Must contain `surfdata_iESM_dyn.nc` file
- File must follow E3SM naming convention

**File Location:**
```
simulation_path/
└── surfdata_iESM_dyn.nc  ← Script looks for this exact filename
```

---

#### `output_file`
**Type:** String (file path)  
**Required:** Yes  
**Description:** Path and filename for output time series data.

**Examples:**
```json
"output_file": "./../2025_DiVittorio_et_al_e3sm/control_time_series_surfdata_iESM_dyn_20240730.dat"
"output_file": "./output/amazon_harvest_grazing.csv"
```

**File Formats:**
- **Fixed-width (.dat)** - Default, human-readable, aligned columns
- **CSV (.csv)** - If filename ends in .csv or `write_to_csv: true`

**Output Structure:**
```
Year  FRAC_VEG  GRAZING_AREA (km^2)  HARVEST_VH1_AREA (km^2)  HARVEST_AREA (km^2)  ...
2015    0.6234        1234.56                567.89                2345.67
2016    0.6198        1289.34                578.12                2401.23
...
```

---

#### `variables`
**Type:** List of strings  
**Required:** Yes  
**Description:** Variables to extract from surfdata_iESM_dyn.nc file.

**Format:** Simple list (not nested like in h0 script)

**Common Usage:**
```json
"variables": ["GRAZING", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3", "HARVEST_VH1", "HARVEST_VH2"]
```

**Available Variables:**

**Human Management:**
- `GRAZING` - Grazing fraction on grasslands
- `HARVEST_VH1` - Very heavy harvest 1
- `HARVEST_VH2` - Very heavy harvest 2
- `HARVEST_SH1` - Secondary heavy harvest 1
- `HARVEST_SH2` - Secondary heavy harvest 2
- `HARVEST_SH3` - Secondary heavy harvest 3

**Land Cover (Automatically included if needed):**
- `PCT_NATVEG` - Natural vegetation percent
- `PCT_NAT_PFT` - PFT percentages

**Note:** If you request grazing or harvest variables, the script automatically includes `PCT_NATVEG` and `PCT_NAT_PFT` to calculate areas.

---

### Optional Parameters

#### `region`
**Type:** String or `None`  
**Required:** No  
**Default:** `None` (global - no filtering)  
**Possible Values:** Region names or `None`

**Description:** Geographic region to filter data extraction.

**Example:**
```json
"region": "amazon"
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

**Usage:**
```json
// Global (default)
"region": null

// Amazon only
"region": "amazon"
```

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
- Extracts data only for years within this range
- Reduces processing time for subset analysis
- Useful for specific time periods

**Data Selection:**
```
Year 2019 → Excluded (< start_year 2020)
Year 2020 → Included
Year 2050 → Included
Year 2051 → Excluded (> end_year 2050)
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
Year  FRAC_VEG  GRAZING_AREA (km^2)
2015    0.6234        1234.56
2016    0.6198        1289.34
```

**CSV (.csv):**
```
Year,FRAC_VEG,GRAZING_AREA (km^2)
2015,0.6234,1234.56
2016,0.6198,1289.34
```

---

## Available Variables

### Direct Extraction Variables

These are variables that exist directly in the surfdata_iESM_dyn.nc file:

#### Grazing

**Variable:** `GRAZING`  
**Original Units:** Dimensionless fraction (0-1)  
**Output Units:** km² (after processing)  
**Description:** Fraction of grassland area used for livestock grazing  
**Output Column:** `GRAZING_AREA (km^2)`

**Calculation:**
```
GRAZING_AREA (km²) = GRAZING × GRASS_AREA (km²)
```

#### Wood Harvesting

**Variable:** `HARVEST_VH1`  
**Original Units:** Dimensionless fraction (0-1)  
**Output Units:** km² (after processing)  
**Description:** Very heavy harvest intensity 1 - fraction of forest harvested  
**Output Column:** `HARVEST_VH1_AREA (km^2)`

**Variable:** `HARVEST_VH2`  
**Original Units:** Dimensionless fraction (0-1)  
**Output Units:** km² (after processing)  
**Description:** Very heavy harvest intensity 2 - fraction of forest harvested  
**Output Column:** `HARVEST_VH2_AREA (km^2)`

**Variable:** `HARVEST_SH1`  
**Original Units:** Dimensionless fraction (0-1)  
**Output Units:** km² (after processing)  
**Description:** Secondary heavy harvest 1 - fraction of forest harvested  
**Output Column:** `HARVEST_SH1_AREA (km^2)`

**Variable:** `HARVEST_SH2`  
**Original Units:** Dimensionless fraction (0-1)  
**Output Units:** km² (after processing)  
**Description:** Secondary heavy harvest 2 - fraction of forest harvested  
**Output Column:** `HARVEST_SH2_AREA (km^2)`

**Variable:** `HARVEST_SH3`  
**Original Units:** Dimensionless fraction (0-1)  
**Output Units:** km² (after processing)  
**Description:** Secondary heavy harvest 3 - fraction of forest harvested  
**Output Column:** `HARVEST_SH3_AREA (km^2)`

**Harvest Calculation:**
```
HARVEST_*_AREA (km²) = HARVEST_* × FOREST_AREA (km²)
```

#### Vegetation Cover

**Variable:** `PCT_NATVEG`  
**Original Units:** Percent (0-100)  
**Output Units:** Fraction (0-1)  
**Description:** Percentage of grid cell covered by natural vegetation  
**Output Column:** `FRAC_VEG`

**Variable:** `PCT_NAT_PFT`  
**Original Units:** Percent (0-100)  
**Output Units:** km² (individual PFT areas)  
**Description:** Percentage of each plant functional type  
**Output Columns:** Individual and aggregated PFT areas

---

### Derived Variables

These are calculated by the script from other variables:

#### Total Harvest Area

**Variable:** Automatically created  
**Output Column:** `HARVEST_AREA (km^2)`  
**Description:** Total forest area harvested (sum of all harvest types)

**Calculation:**
```
HARVEST_AREA (km²) = HARVEST_VH1_AREA + HARVEST_VH2_AREA + 
                      HARVEST_SH1_AREA + HARVEST_SH2_AREA + HARVEST_SH3_AREA
```

#### Grid Cell Area

**Variable:** Automatically included  
**Output Column:** `AREA (km^2)`  
**Description:** Total grid cell area

#### PFT Areas

**Variables:** Automatically calculated when PCT_NAT_PFT requested  
**Output Columns:**
- Individual PFTs: `PFT_1_AREA (km^2)` through `PFT_16_AREA (km^2)`
- Aggregated groups:
  - `BARE_AREA (km^2)` - Bare soil
  - `FOREST_AREA (km^2)` - All tree PFTs (1-8)
  - `SHRUB_AREA (km^2)` - All shrub PFTs (9-11)
  - `GRASS_AREA (km^2)` - All grass PFTs (12-14)
  - `CROP_AREA (km^2)` - Crop PFT (15)

**17 PFTs in E3SM:**
```
0:      Bare soil
1-8:    Trees (8 types)
9-11:   Shrubs (3 types)
12-14:  Grasses (3 types)
15:     Crop
16:     Empty (ignored)
```

**Calculation:**
```
PFT_AREA (km²) = AREA (km²) × FRAC_VEG × PFT_fraction
```

---

## Variable Processing

### Automatic Processing Steps

The script performs several automatic processing steps:

#### 1. PFT Data Inclusion

**Automatic Action:** If any grazing or harvest variable is requested, script automatically includes `PCT_NAT_PFT` and `PCT_NATVEG`

**Reason:** Needed to convert fractions to areas

**Example:**
```json
// You specify:
"variables": ["GRAZING", "HARVEST_VH1"]

// Script automatically adds:
// - PCT_NATVEG
// - PCT_NAT_PFT
```

---

#### 2. Duplicate Time Value Removal

**Processing:** Removes duplicate time coordinate values, keeping last occurrence

**Reason:** Simulation restarts can create duplicate time entries

**Example:**
```
Before: time = [2015, 2016, 2016*, 2017]  (* from restart)
After:  time = [2015, 2016, 2017]
```

---

#### 3. Vegetation Fraction Calculation

**Input:** `PCT_NATVEG` (percent, 0-100)

**Processing:**
```python
FRAC_VEG = PCT_NATVEG / 100
```

**Output:** `FRAC_VEG` (fraction, 0-1)

---

#### 4. PFT Area Calculation

**For Each PFT:**
```python
PFT_AREA (km²) = AREA (km²) × FRAC_VEG × (PCT_NAT_PFT / 100)
```

**Aggregated Groups:**
```python
FOREST_AREA (km²) = sum of PFTs 1-8
SHRUB_AREA (km²) = sum of PFTs 9-11
GRASS_AREA (km²) = sum of PFTs 12-14
```

---

#### 5. Grazing Area Conversion

**Input:** `GRAZING` (fraction, 0-1)

**Processing:**
```python
GRAZING_AREA (km²) = GRAZING × GRASS_AREA (km²)
```

**Reason:** Grazing only occurs on grasslands

**Output:** `GRAZING_AREA (km^2)` column

---

#### 6. Harvest Area Conversion

**Input:** `HARVEST_*` variables (fractions, 0-1)

**Processing:**
```python
HARVEST_VH1_AREA (km²) = HARVEST_VH1 × FOREST_AREA (km²)
HARVEST_VH2_AREA (km²) = HARVEST_VH2 × FOREST_AREA (km²)
HARVEST_SH1_AREA (km²) = HARVEST_SH1 × FOREST_AREA (km²)
HARVEST_SH2_AREA (km²) = HARVEST_SH2 × FOREST_AREA (km²)
HARVEST_SH3_AREA (km²) = HARVEST_SH3 × FOREST_AREA (km²)
```

**Reason:** Harvesting only occurs on forests

**Output:** Individual `HARVEST_*_AREA (km^2)` columns

---

#### 7. Total Harvest Calculation

**Processing:**
```python
HARVEST_AREA (km²) = sum of all HARVEST_*_AREA columns present
```

**Output:** `HARVEST_AREA (km^2)` column (placed after individual harvest columns)

---

#### 8. Spatial Aggregation

**For All Variables:**
- Sum over all latitude/longitude coordinates
- Creates global or regional totals

**For FRAC_VEG:**
- Calculate area-weighted mean
- Formula: `Σ(FRAC_VEG × AREA) / Σ(AREA)`

---

## Regional Analysis

### How Regional Filtering Works

1. **Load NetCDF file**
2. **Apply regional bounds** (if region specified)
   - Filter gridcells outside lat/lon bounding box
3. **Process data** for regional gridcells only
4. **Aggregate** over regional gridcells

### Regional Bounds

Defined by lat/lon bounding boxes in `utility_e3sm_netcdf.py`:

**Example Regions:**
```python
amazon = [-85°, -35°] longitude, [-25°, 15°] latitude
conus = [-125.25°, -66.25°] longitude, [23.25°, 54.75°] latitude
```

### Example: Amazon Harvest Analysis

```json
{
    "simulation_path": "/path/to/simulation",
    "output_file": "./amazon_harvest.dat",
    "variables": ["HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3", "HARVEST_VH1", "HARVEST_VH2"],
    "region": "amazon",
    "start_year": 2015,
    "end_year": 2100
}
```

**Output:**
- Harvest areas for Amazon Basin only
- Regional totals (not global)

---

## JSON Configuration Examples

### Example 1: Global Grazing and Harvest

**Purpose:** Extract all grazing and harvest activities globally

```json
{
    "simulation_path": "/lcrc/group/e3sm/ac.eva.sinha/20240730_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3_without_feedbacks/run",
    "output_file": "./../2025_DiVittorio_et_al_e3sm/control_time_series_surfdata_iESM_dyn_20240730.dat",
    "variables": ["GRAZING", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3", "HARVEST_VH1", "HARVEST_VH2"],
    "start_year": 2015,
    "end_year": 2100
}
```

**What it does:**
- Extracts grazing on grasslands
- Extracts all 5 harvest types on forests
- Calculates total harvest area automatically
- Global aggregation (no regional filtering)
- 86 years of annual data (2015-2100)

**Output Columns:**
```
Year
FRAC_VEG
AREA (km^2)
PFT_1_AREA (km^2) ... PFT_16_AREA (km^2)
BARE_AREA (km^2)
FOREST_AREA (km^2)
SHRUB_AREA (km^2)
GRASS_AREA (km^2)
CROP_AREA (km^2)
GRAZING_AREA (km^2)
HARVEST_SH1_AREA (km^2)
HARVEST_SH2_AREA (km^2)
HARVEST_SH3_AREA (km^2)
HARVEST_VH1_AREA (km^2)
HARVEST_VH2_AREA (km^2)
HARVEST_AREA (km^2)
```

---

### Example 2: Regional Analysis (Amazon)

**Purpose:** Focus on Amazon Basin grazing and harvest

```json
{
    "simulation_path": "/lcrc/group/e3sm/ac.eva.sinha/20240730_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3_without_feedbacks/run",
    "output_file": "./../2025_DiVittorio_et_al_e3sm/control_time_series_surfdata_iESM_dyn_20240730_amazon.dat",
    "variables": ["GRAZING", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3", "HARVEST_VH1", "HARVEST_VH2"],
    "region": "amazon",
    "start_year": 2015,
    "end_year": 2100
}
```

**What it does:**
- Filters to Amazon Basin coordinates
- Extracts grazing and harvest for Amazon only
- Regional totals (not global)
- Same variables as Example 1, but for Amazon

**Use Case:**
- Analyzing deforestation in Amazon
- Tracking agricultural expansion
- Regional land use impacts

---

### Example 3: Grazing Only

**Purpose:** Extract only grazing data (no harvest)

```json
{
    "simulation_path": "/lcrc/group/e3sm/simulation/run",
    "output_file": "./grazing_only.dat",
    "variables": ["GRAZING"],
    "start_year": 2015,
    "end_year": 2100
}
```

**What it does:**
- Extracts grazing fraction
- Automatically includes PCT_NAT_PFT and PCT_NATVEG
- Calculates grass area
- Converts grazing to area (km²)

**Output Focus:**
- `GRAZING_AREA (km^2)`
- `GRASS_AREA (km^2)`
- PFT areas included for context

---

### Example 4: Harvest Only

**Purpose:** Extract only harvest activities (no grazing)

```json
{
    "simulation_path": "/lcrc/group/e3sm/simulation/run",
    "output_file": "./harvest_only.dat",
    "variables": ["HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3", "HARVEST_VH1", "HARVEST_VH2"],
    "start_year": 2015,
    "end_year": 2100
}
```

**What it does:**
- Extracts all 5 harvest types
- Calculates total harvest automatically
- Converts to forest area harvested (km²)

**Output Focus:**
- Individual harvest areas
- `HARVEST_AREA (km^2)` (total)
- `FOREST_AREA (km^2)` for context

---

### Example 5: Short Time Period for Testing

**Purpose:** Quick test with subset of data

```json
{
    "simulation_path": "/lcrc/group/e3sm/simulation/run",
    "output_file": "./test_surfdata.dat",
    "variables": ["GRAZING", "HARVEST_VH1"],
    "start_year": 2015,
    "end_year": 2020
}
```

**What it does:**
- Only 6 years of data (2015-2020)
- Fast execution for testing
- Verify configuration before full run
- Two representative variables

---

### Example 6: CSV Output

**Purpose:** Output in CSV format for easy import to other tools

```json
{
    "simulation_path": "/lcrc/group/e3sm/simulation/run",
    "output_file": "./harvest_grazing.csv",
    "variables": ["GRAZING", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3", "HARVEST_VH1", "HARVEST_VH2"],
    "start_year": 2015,
    "end_year": 2100,
    "write_to_csv": true
}
```

**What it does:**
- Same as Example 1 but CSV output
- Comma-separated values
- Easy import to Excel, R, Python pandas

---

### Example 7: Multiple Regions in One JSON

**Purpose:** Process multiple regions with one command

```json
[
    {
        "simulation_path": "/lcrc/group/e3sm/simulation/run",
        "output_file": "./global_harvest_grazing.dat",
        "variables": ["GRAZING", "HARVEST_VH1", "HARVEST_VH2"],
        "start_year": 2015,
        "end_year": 2100
    },
    {
        "simulation_path": "/lcrc/group/e3sm/simulation/run",
        "output_file": "./amazon_harvest_grazing.dat",
        "variables": ["GRAZING", "HARVEST_VH1", "HARVEST_VH2"],
        "region": "amazon",
        "start_year": 2015,
        "end_year": 2100
    },
    {
        "simulation_path": "/lcrc/group/e3sm/simulation/run",
        "output_file": "./conus_harvest_grazing.dat",
        "variables": ["GRAZING", "HARVEST_VH1", "HARVEST_VH2"],
        "region": "conus",
        "start_year": 2015,
        "end_year": 2100
    }
]
```

**What it does:**
- Processes three configurations in one run
- Global, Amazon, and CONUS analyses
- Efficient batch processing

---

## Output Files

### Fixed-Width Format (.dat)

**Default format** - Human-readable with aligned columns

**Example:**
```
Year  FRAC_VEG  GRAZING_AREA (km^2)  HARVEST_VH1_AREA (km^2)  HARVEST_AREA (km^2)
2015    0.6234        1234.56                567.89                2345.67
2016    0.6198        1289.34                578.12                2401.23
2017    0.6175        1301.45                589.34                2456.78
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
Year,FRAC_VEG,GRAZING_AREA (km^2),HARVEST_VH1_AREA (km^2),HARVEST_AREA (km^2)
2015,0.6234,1234.56,567.89,2345.67
2016,0.6198,1289.34,578.12,2401.23
2017,0.6175,1301.45,589.34,2456.78
```

**Characteristics:**
- Comma-separated
- Smaller file size
- Easy import to Excel, pandas, R

---

### Output Structure

**Temporal Organization:**
- One row per year
- Sorted by Year
- Annual time series (not monthly)

**Column Organization:**
```
[Year] [FRAC_VEG] [AREA] [PFT areas...] [Grazing area] [Harvest areas...] [Total harvest]
```

**Column Headers:**
```
Variable_name (units)

Examples:
Year
FRAC_VEG
GRAZING_AREA (km^2)
HARVEST_VH1_AREA (km^2)
HARVEST_AREA (km^2)
FOREST_AREA (km^2)
```

---

## Troubleshooting

### Issue 1: File Not Found

**Error Message:**
```
FileNotFoundError: surfdata_iESM_dyn.nc
```

**Causes:**
- Incorrect `simulation_path`
- File doesn't exist
- Wrong simulation directory

**Solutions:**
```bash
# Check directory exists
ls -la /path/to/simulation/

# Check for surfdata file
ls /path/to/simulation/surfdata_iESM_dyn.nc

# Verify exact filename
# Must be: surfdata_iESM_dyn.nc
```

---

### Issue 2: Variable Not Found

**Error Message:**
```
KeyError: 'VARIABLE_NAME'
```

**Causes:**
- Variable not in surfdata file
- Typo in variable name
- File doesn't contain requested variable

**Solutions:**
```python
# Check what variables are in file
import xarray as xr
ds = xr.open_dataset('surfdata_iESM_dyn.nc')
print(list(ds.data_vars))

# Common variables:
# GRAZING, HARVEST_VH1, HARVEST_VH2, HARVEST_SH1, HARVEST_SH2, HARVEST_SH3
# PCT_NATVEG, PCT_NAT_PFT
```

---

### Issue 3: Memory Error

**Error Message:**
```
MemoryError
```

**Causes:**
- Large surfdata file
- Long time series
- Many variables
- Insufficient RAM

**Solutions:**
1. **Reduce year range:**
```json
"start_year": 2015,
"end_year": 2050  // Instead of 2100
```

2. **Process fewer variables:**
```json
"variables": ["GRAZING", "HARVEST_VH1"]  // Instead of all harvest types
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

### Issue 4: Duplicate Time Values

**Error Message:**
```
Warning: Duplicate time values found
```

**Cause:**
- Simulation restarts create duplicate time entries
- Script handles this automatically

**Behavior:**
- Script removes duplicates, keeping last occurrence
- No user action needed
- Informational warning only

---

### Issue 5: Slow Execution

**Symptom:** Script taking longer than expected

**Causes:**
- Large surfdata file
- Many years to process
- Limited CPU cores

**Solutions:**

1. **Parallel processing (already built-in):**
Script automatically uses all CPU cores

2. **Reduce year range:**
```json
"start_year": 2080,
"end_year": 2100  // Only 21 years
```

3. **Run on HPC:**
Use high-performance computing cluster

4. **Monitor progress:**
Script prints completion time

**Typical Timing:**
- 86 years, 6 variables: 30-120 seconds
- Depends on file size and CPU count

---

### Issue 6: Regional Bounds Error

**Error Message:**
```
Did not recognize the selected region!
```

**Cause:**
- Region name typo
- Region not in predefined list

**Solutions:**
```json
// Check spelling (case-sensitive)
"region": "amazon"  // Correct
"region": "Amazon"  // Wrong

// Use null for global
"region": null

// Check available regions in utility_e3sm_netcdf.py
```

---

## Best Practices

### 1. Directory Organization

**Recommended Structure:**
```
project/
├── simulations/
│   ├── control/
│   │   └── run/
│   │       └── surfdata_iESM_dyn.nc
│   ├── scenario1/
│   │   └── run/
│   │       └── surfdata_iESM_dyn.nc
│   └── scenario2/
│       └── run/
│           └── surfdata_iESM_dyn.nc
├── output/
│   ├── surfdata_time_series/
│   └── processed/
├── configs/
│   └── surfdata_extraction_configs.json
└── analysis/
    └── scripts/
```

---

### 2. Configuration File Organization

**Use descriptive names:**
```
configs/
├── global_harvest_grazing_extraction.json
├── regional_amazon_extraction.json
└── testing_2015_2020.json
```

**Group related extractions:**
```json
[
    {"output_file": "control_global.dat", ...},
    {"output_file": "control_amazon.dat", "region": "amazon", ...},
    {"output_file": "scenario_global.dat", ...},
    {"output_file": "scenario_amazon.dat", "region": "amazon", ...}
]
```

---

### 3. Testing Strategy

**Start small:**
```json
{
    "start_year": 2015,
    "end_year": 2016,  // Just 2 years
    "variables": ["GRAZING"]  // One variable
}
```

**Then scale up:**
```json
{
    "start_year": 2015,
    "end_year": 2100,  // Full range
    "variables": ["GRAZING", "HARVEST_VH1", "HARVEST_VH2", ...]
}
```

---

### 4. Variable Selection

**All human activities (comprehensive):**
```json
"variables": ["GRAZING", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3", "HARVEST_VH1", "HARVEST_VH2"]
```

**Grazing focus:**
```json
"variables": ["GRAZING"]
```

**Harvest focus:**
```json
"variables": ["HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3", "HARVEST_VH1", "HARVEST_VH2"]
```

**Simplified (only major harvest):**
```json
"variables": ["HARVEST_VH1", "HARVEST_VH2"]
```

---

### 5. Regional vs Global

**When to use global:**
- Earth system scale analysis
- Total land use impacts
- Global harvest/grazing trends

**When to use regional:**
- Hotspot analysis (Amazon deforestation)
- Regional impacts (CONUS agriculture)
- Policy-relevant scales

**Multiple regions:**
```json
[
    {"output_file": "global.dat"},
    {"output_file": "amazon.dat", "region": "amazon"},
    {"output_file": "conus.dat", "region": "conus"}
]
```

---

### 6. Performance Optimization

**Parallel processing:**
- Script uses all CPU cores automatically
- One process per year
- Scales well with more cores

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
# Expected: 86 rows for 2015-2100

# Check for missing values
print(df.isnull().sum())

# Check value ranges
print(df.describe())

# Verify temporal continuity
print(df['Year'].head(10))
# Should be consecutive years
```

**Check totals:**
```python
# Total harvest should equal sum of individual harvests
harvest_cols = [col for col in df.columns if 'HARVEST_' in col and col != 'HARVEST_AREA (km^2)']
calculated_total = df[harvest_cols].sum(axis=1)
reported_total = df['HARVEST_AREA (km^2)']
print(f"Match: {np.allclose(calculated_total, reported_total)}")
```

---

### 8. Reproducibility

**Document your extraction:**

```json
{
    "comment": "Harvest and grazing for control simulation",
    "date_extracted": "2026-01-24",
    "simulation": "SSP245_ZATM_BGC_without_feedbacks",
    "simulation_path": "/lcrc/group/e3sm/.../run",
    ...
}
```

**Version control configurations:**
```bash
git add configs/
git commit -m "Add surfdata extraction configs for control simulation"
```

---

## Comparison with h0 Extraction Script

### Key Differences

| Feature | e3sm_extract_time_series_h0.py | e3sm_extract_time_series_surfdata_iesm_dyn.py |
|---------|-------------------------------|-----------------------------------------------|
| **Input Files** | Multiple monthly NetCDF files (elm.h0.*, eam.h0.*) | Single NetCDF file (surfdata_iESM_dyn.nc) |
| **Temporal Resolution** | Monthly | Annual |
| **Variables** | Climate, carbon cycle (100+ available) | Human land management (6 main variables) |
| **File Type Config** | `netcdf_substrings` (nested list) | Not needed (always surfdata file) |
| **Variables Config** | `variables` (nested list, one per file type) | `variables` (simple list) |
| **Aggregation Config** | `lat_lon_aggregation_types` parameter | Always area-weighted sum |
| **Processing** | Optional (`process_variables`) | Automatic (always converts fractions to areas) |
| **Parallel Processing** | One process per file (per month per file type) | One process per year |
| **Output Frequency** | Monthly rows | Annual rows |

### When to Use Each Script

**Use h0 script when:**
- Need climate variables (temperature, precipitation)
- Need carbon cycle fluxes (GPP, NPP, NBP)
- Need monthly resolution
- Analyzing atmosphere or land model outputs

**Use surfdata script when:**
- Need human land management data
- Need grazing information
- Need wood harvest information
- Analyzing E3SM-GCAM coupled simulations
- Annual resolution sufficient

### Example Workflow Using Both Scripts

```
1. Extract climate and carbon cycle (monthly):
   python e3sm_extract_time_series_h0.py h0_config.json
   → Output: Monthly GPP, NPP, temperature, precipitation

2. Extract human land management (annual):
   python e3sm_extract_time_series_surfdata_iesm_dyn.py surfdata_config.json
   → Output: Annual grazing, harvest areas

3. Analyze together:
   - Monthly climate and carbon dynamics
   - Annual human activities
   - Correlate human activities with ecosystem responses
```

---

## Integration with Other Scripts

### Typical E3SM-GCAM Analysis Workflow

```
1. Run E3SM-GCAM coupled simulation
   ↓ (Produces surfdata_iESM_dyn.nc + monthly h0 files)
   
2a. e3sm_extract_time_series_h0.py
    ↓ (Extract climate and carbon cycle - monthly)
    
2b. e3sm_extract_time_series_surfdata_iesm_dyn.py  ← THIS SCRIPT
    ↓ (Extract harvest and grazing - annual)
    
3. Combine analyses:
   - Plot harvest trends over time
   - Correlate with carbon cycle changes
   - Compare regional impacts
   - Validate coupling between E3SM and GCAM
```

### Downstream Analysis

**Time Series Plotting:**
```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('harvest_grazing.csv')

# Plot total harvest over time
plt.plot(df['Year'], df['HARVEST_AREA (km^2)'])
plt.xlabel('Year')
plt.ylabel('Total Harvest Area (km²)')
plt.title('Global Forest Harvest')
plt.show()
```

**Regional Comparison:**
```python
# Load global and regional data
df_global = pd.read_csv('global_harvest.csv')
df_amazon = pd.read_csv('amazon_harvest.csv')

# Calculate Amazon fraction
amazon_fraction = df_amazon['HARVEST_AREA (km^2)'] / df_global['HARVEST_AREA (km^2)']

plt.plot(df_global['Year'], amazon_fraction * 100)
plt.xlabel('Year')
plt.ylabel('Amazon % of Global Harvest')
plt.show()
```

**Statistical Analysis:**
```python
from scipy import stats

df = pd.read_csv('harvest_grazing.csv')

# Test for trend
slope, intercept, r_value, p_value, std_err = stats.linregress(
    df['Year'], df['HARVEST_AREA (km^2)']
)
print(f"Harvest trend: {slope:.2f} km²/year, p={p_value:.4f}")
```

---

## Appendix: Complete Variable Reference

### Available in surfdata_iESM_dyn.nc

| Variable | Units | Description | Output After Processing |
|----------|-------|-------------|-------------------------|
| `GRAZING` | Fraction (0-1) | Grazing intensity on grasslands | `GRAZING_AREA (km^2)` |
| `HARVEST_VH1` | Fraction (0-1) | Very heavy harvest intensity 1 | `HARVEST_VH1_AREA (km^2)` |
| `HARVEST_VH2` | Fraction (0-1) | Very heavy harvest intensity 2 | `HARVEST_VH2_AREA (km^2)` |
| `HARVEST_SH1` | Fraction (0-1) | Secondary heavy harvest 1 | `HARVEST_SH1_AREA (km^2)` |
| `HARVEST_SH2` | Fraction (0-1) | Secondary heavy harvest 2 | `HARVEST_SH2_AREA (km^2)` |
| `HARVEST_SH3` | Fraction (0-1) | Secondary heavy harvest 3 | `HARVEST_SH3_AREA (km^2)` |
| `PCT_NATVEG` | Percent (0-100) | Natural vegetation cover | `FRAC_VEG` (fraction, 0-1) |
| `PCT_NAT_PFT` | Percent (0-100) | PFT percentages | Individual and aggregated PFT areas |

### Automatically Generated Output Columns

| Column | Units | Description |
|--------|-------|-------------|
| `Year` | Year | Year of data |
| `FRAC_VEG` | Fraction (0-1) | Vegetation fraction (area-weighted mean) |
| `AREA (km^2)` | km² | Total gridcell area (global or regional sum) |
| `PFT_1_AREA (km^2)` ... `PFT_16_AREA (km^2)` | km² | Individual PFT areas |
| `BARE_AREA (km^2)` | km² | Bare soil area |
| `FOREST_AREA (km^2)` | km² | Total forest area (PFTs 1-8) |
| `SHRUB_AREA (km^2)` | km² | Total shrub area (PFTs 9-11) |
| `GRASS_AREA (km^2)` | km² | Total grass area (PFTs 12-14) |
| `CROP_AREA (km^2)` | km² | Crop area (PFT 15) |
| `GRAZING_AREA (km^2)` | km² | Grazed grassland area |
| `HARVEST_*_AREA (km^2)` | km² | Harvested forest area by type |
| `HARVEST_AREA (km^2)` | km² | Total harvested area (sum of all types) |

---

## References

- E3SM Documentation: [https://e3sm.org/](https://e3sm.org/)
- E3SM-GCAM Coupling: DiVittorio et al. (2025). "E3SM-GCAM coupling methodology and applications." *Journal of Advances in Modeling Earth Systems*. [DOI: 10.1029/2024MS004806](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024MS004806)
- GitHub Repository: [https://github.com/philipmyint/e3sm--gcam_analysis](https://github.com/philipmyint/e3sm--gcam_analysis)
- xarray Documentation: [https://xarray.pydata.org/](https://xarray.pydata.org/)

---

## Contact

For questions or issues:
- **Philip Myint**: myint1@llnl.gov
- **Dalei Hao**: dalei.hao@pnnl.gov

---

## Version Information

**Script:** e3sm_extract_time_series_surfdata_iesm_dyn.py  
**Utility Modules:** utility_e3sm_netcdf.py, utility_constants.py, utility_dataframes.py, utility_functions.py  
**Documentation Version:** 1.0  
**Last Updated:** January 2026  
**Python Version:** 3.7+  
**xarray Version:** 0.16+

---

*This documentation provides comprehensive guidance for using the `e3sm_extract_time_series_surfdata_iesm_dyn.py` script to extract annual time series data of human land management activities (grazing and harvest) from E3SM-GCAM coupled simulation surface data files.*
