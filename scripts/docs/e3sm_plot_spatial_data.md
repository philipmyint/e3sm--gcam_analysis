# E3SM Spatial Data Plotting Script Documentation

## Overview

**Script Name:** `e3sm_plot_spatial_data.py`

**Purpose:** Creates publication-quality spatial (map) plots from E3SM NetCDF data files with support for absolute differences, percent differences, means, sums, ensemble analysis, statistical significance testing (stippling), and separate visualizations for different scenarios.

**Key Capabilities:**
- Global and regional spatial maps
- Absolute and percent difference plots
- Mean or sum across multiple individual files
- Ensemble mean with statistical significance stippling
- Separate plots for individual scenarios or ensemble members
- Support for both ELM (structured grid) and EAM (unstructured grid)
- Customizable colormaps and projections
- Statistical significance testing (aggregated and per-gridcell)
- Parallel processing for efficiency
- Flexible ensemble file organization

**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Requirements](#installation--requirements)
3. [Basic Usage](#basic-usage)
4. [Complete Parameter Reference](#complete-parameter-reference)
5. [Plot Types Explained](#plot-types-explained)
6. [Individual vs Ensemble Plots](#individual-vs-ensemble-plots)
7. [Statistical Testing Explained](#statistical-testing-explained)
8. [Ensemble File Organization](#ensemble-file-organization)
9. [ELM vs EAM Differences](#elm-vs-eam-differences)
10. [Comprehensive JSON Examples](#comprehensive-json-examples)
11. [Output Files](#output-files)
12. [Best Practices](#best-practices)
13. [Troubleshooting](#troubleshooting)

---

## Installation & Requirements

### Required Python Libraries

```bash
pip install xarray uxarray cartopy matplotlib scipy numpy pandas multiprocessing
```

### Required Utility Modules

- `utility_constants` - Physical constants
- `utility_dataframes` - DataFrame operations (t-tests)
- `utility_functions` - General utilities (includes `transpose_scenarios_if_needed()`)
- `utility_plots` - Plotting defaults and functions
- `utility_xarray` - xarray/uxarray operations

### Additional Requirements

- **For EAM plots:** Grid file (unstructured mesh, e.g., `ne30pg2_scrip_c20191218.nc`)
- **For ELM plots:** No additional files needed (structured lat/lon grid)

---

## Basic Usage

### Command Line Execution

```bash
python e3sm_plot_spatial_data.py config.json
```

**Multiple Configuration Files:**
```bash
python e3sm_plot_spatial_data.py elm_config.json eam_config.json
```

### What the Script Does

1. Reads JSON configuration file(s)
2. Loads NetCDF spatial data files (from `e3sm_extract_spatial_data_h0.py`)
3. Calculates temporal means/sums/differences over specified years
4. Performs statistical testing (aggregated and per-gridcell)
5. Creates spatial maps with:
   - Coastlines
   - Colorbars
   - Statistical annotations (min, mean, median, max, std)
   - Stippling for significance (ELM only)
6. Saves plots as PDF or PNG
7. Outputs p-value files

---

## Complete Parameter Reference

### Required Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `netcdf_files` | string/list/nested list | **Yes** | - | NetCDF file path(s) |
| `plot_directory` | string | No | `'./'` | Output directory for plots |

### Core Optional Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `variables` | list or `'all'` | No | `'all'` | Variables to plot |
| `plot_type` | string | No | `'absolute_difference'` | Plot type (see below) |
| `start_year` | integer | No | `2071` | First year for temporal averaging |
| `end_year` | integer | No | `2090` | Last year for temporal averaging |
| `time_calculation` | string | No | `'mean'` | `'mean'` or `'sum'` for temporal aggregation |
| `grid_file` | string | No | `None` | **Required for EAM**, grid file path |
| `netcdf_file_sets` | list | No | Auto-generated | Labels for ensemble file sets |

**Plot Type Options:**
- `'absolute_difference'` - Scenario minus control (default)
- `'percent_difference'` - (Scenario - Control) / Control × 100
- `'separate_plots'` - Individual plots for each scenario/ensemble
- `'mean'` - Mean across multiple individual files
- `'sum'` - Sum across multiple individual files

### Visual Customization

| Parameter | Type | Required | Default | Description |
|-----------|------|---------|---------|-------------|
| `width` | float | No | `10` | Figure width (inches) |
| `height` | float | No | `8` | Figure height (inches) |
| `cmap` | string | No | `'bwr'` | Colormap (blue-white-red) |
| `projection` | cartopy projection | No | `Robinson` | Map projection |
| `use_latex` | boolean | No | `false` | Use LaTeX fonts |
| `produce_png` | boolean | No | `false` | Output PNG instead of PDF |
| `title_size` | integer | No | `24` | Title font size |
| `statistics_panel_size` | integer | No | `14` | Statistics text size |

### Colorbar Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|---------|---------|-------------|
| `cbar_on` | boolean | No | `true` | Show colorbar |
| `cbar_limits` | list `[min, max]` | No | `None` | Colorbar limits (auto if None) |
| `cbar_label_size` | integer | No | `20` | Colorbar tick label size |
| `cbar_x_offset` | float | No | `0.06` | Colorbar horizontal offset |

### Stippling Parameters (Statistical Significance)

| Parameter | Type | Required | Default | Description |
|-----------|------|---------|---------|-------------|
| `stippling_on` | boolean | No | `false` | Enable stippling |
| `stippling_hatches` | string | No | `'xxxx'` | Hatch pattern (e.g., `'///'`, `'xxx'`) |
| `stippling_std_multiple` | float | No | `2` | Std deviation multiple for stippling |

**Stippling Behavior:**
- **Ensemble data (≥2 members per set):** Stipples where p < threshold (ELM only)
- **Single dataset:** Stipples where |value| > mean ± N×std

### Statistical Testing Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|---------|---------|-------------|
| `p_value_threshold` | float | No | `0.05` | Significance threshold |
| `p_value_file` | string | No | `'p_values.dat'` | Output file for p-values |
| `p_value_file_print_only_if_below_threshold` | boolean | No | `true` | Only print significant p-values |

### Advanced Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|---------|---------|-------------|
| `multiplier` | float/dict | No | `1` | Scale data values |
| `plot_name` | dict | No | Auto-generated | Custom plot filenames per variable |
| `title` | dict | No | Auto-generated | Custom titles per variable |

---

## Plot Types Explained

### 1. Absolute Difference (`plot_type: "absolute_difference"`)

**Description:** Shows scenario minus control in original units

**When to Use:**
- Comparing two scenarios/ensembles
- Want actual magnitude of changes
- Variables where absolute change is meaningful

**Requirements:** Exactly 2 files/ensembles

**Formula:**
```
Difference = Scenario - Control
```

---

### 2. Percent Difference (`plot_type: "percent_difference"`)

**Description:** Shows percent change relative to control

**When to Use:**
- Emphasizing relative changes
- Variables with widely varying magnitudes
- Comparing fractional impacts

**Requirements:** Exactly 2 files/ensembles

**Formula:**
```
% Difference = (Scenario - Control) / |Control| × 100
```

**Auto-limits:** Sets colorbar to [-100, 100]% if max exceeds 100%

---

### 3. Separate Plots (`plot_type: "separate_plots"`)

**Description:** Creates individual maps for each scenario or ensemble

**When to Use:**
- Want to see absolute values for each scenario
- Side-by-side comparison in publications
- Showing spatial patterns individually

**Works with:** Any number of files/ensembles

**Output:**
- Individual plots: One map per file
- Ensemble plots: One ensemble mean map per set

---

### 4. Mean (`plot_type: "mean"`)

**Description:** Averages across multiple individual files

**When to Use:**
- Multiple files from same scenario (e.g., ensemble members)
- Want to see average spatial pattern
- Reducing noise across realizations

**Requirements:** Simple list of files (not nested)

**Processing:**
```
Mean = (File1 + File2 + ... + FileN) / N
```

**Example:**
```json
{
    "netcdf_files": ["file1.nc", "file2.nc", "file3.nc", "file4.nc", "file5.nc"],
    "plot_type": "mean"
}
```

**Output:** Single map showing average across all 5 files

---

### 5. Sum (`plot_type: "sum"`)

**Description:** Sums across multiple individual files

**When to Use:**
- Variables that should be aggregated (e.g., precipitation components)
- Combining related variables
- Total fluxes across ensemble

**Requirements:** Simple list of files (not nested)

**Processing:**
```
Sum = File1 + File2 + ... + FileN
```

**Example:**
```json
{
    "netcdf_files": ["precip_conv.nc", "precip_large.nc"],
    "plot_type": "sum",
    "variables": ["PRECIP"]
}
```

**Output:** Single map showing total precipitation (convective + large-scale)

---

## Individual vs Ensemble Plots

### Individual Plots (Simple List)

**Structure:** Flat list of files
```json
{
    "netcdf_files": ["file1.nc", "file2.nc", "file3.nc"]
}
```

**Processing Options:**
- **2 files + `absolute_difference`:** Shows difference
- **2 files + `percent_difference`:** Shows % difference
- **N files + `separate_plots`:** Shows N individual maps
- **N files + `mean`:** Shows average across N files
- **N files + `sum`:** Shows sum across N files

**Statistical Testing:** Aggregated only (when differencing)

**No per-gridcell testing** (no ensemble structure)

---

### Ensemble Plots (Nested List)

**Structure:** List of lists (rows = file sets)
```json
{
    "netcdf_files": [
        ["ctrl.nc", "ctrl_2.nc", "ctrl_3.nc"],
        ["scen.nc", "scen_2.nc", "scen_3.nc"]
    ],
    "netcdf_file_sets": ["Control", "Scenario"]
}
```

**Limitation:** Currently limited to **exactly 2 ensembles**

**Processing:**
- Calculates ensemble mean for each set
- Then applies plot_type (absolute_difference, percent_difference, separate_plots)

**Statistical Testing:**
- Aggregated across all gridcells
- Per-gridcell (ELM only) with stippling

**Requirements for stippling:**
- Nested list structure
- ≥2 members per ensemble
- ELM data only
- `stippling_on: true`

---

## Statistical Testing Explained

The script performs **two types of statistical tests**:

### 1. Aggregated t-test (Always Performed When Comparing)

**What:** Tests if means differ across all gridcells combined

**When:** Both individual and ensemble configurations with 2 files/sets

**How:**
- Takes all gridcell values for first file/ensemble
- Takes all gridcell values for second file/ensemble
- Performs two-sample t-test
- Records p-value in file

**Example Output:**
```
p_values.dat:
GPP in Scenario: 1.234e-04
```

**Interpretation:** Scenario is significantly different from Control globally

---

### 2. Per-Gridcell t-test (Ensemble Only, ELM Only)

**What:** Tests if means differ at each individual gridcell

**When:** 
- Nested list configuration (ensemble)
- ≥2 members per ensemble
- **ELM data only** (too slow for EAM unstructured grid)
- `stippling_on: true`

**How:**
- At gridcell (45°N, 90°W): Compare Control ensemble vs Scenario ensemble
- Performs t-test for that gridcell
- If p < threshold, adds stippling to that gridcell
- Repeats for all gridcells

**Visual Result:** Stippling (hash marks) appears where differences are significant

**Example:** Amazon basin shows dense stippling → significant GPP changes there

**EAM Limitation:** Per-gridcell tests disabled for EAM due to computational cost with unstructured grids

---

### Statistical Testing Summary

**Individual Files (simple list):**
```
2 files + difference:
  ✓ Aggregated (all gridcells) - saved to file
  ✗ Per-gridcell - not applicable

N files + mean/sum:
  ✗ No statistical testing
```

**Ensemble Files (nested list, 2 ensembles only):**
```
✓ Aggregated (all gridcells) - saved to file
✓ Per-gridcell (ELM only) - stippling on map
✗ Per-gridcell (EAM) - disabled (too slow)
```

---

## Ensemble File Organization

When plotting ensembles, you can organize your file lists in **two equivalent ways**:

### Option 1: Organized by File Set (Intuitive)

Each row groups all ensemble members for one scenario:

```json
{
    "netcdf_files": [
        ["control.nc", "control_2.nc", "control_3.nc", "control_4.nc", "control_5.nc"],
        ["feedback.nc", "feedback_2.nc", "feedback_3.nc", "feedback_4.nc", "feedback_5.nc"]
    ],
    "netcdf_file_sets": ["Control", "Full feedback"]
}
```

---

### Option 2: Organized by Ensemble Member

Each row contains one member from each scenario:

```json
{
    "netcdf_files": [
        ["control.nc", "feedback.nc"],
        ["control_2.nc", "feedback_2.nc"],
        ["control_3.nc", "feedback_3.nc"],
        ["control_4.nc", "feedback_4.nc"],
        ["control_5.nc", "feedback_5.nc"]
    ],
    "netcdf_file_sets": ["Control", "Full feedback"]
}
```

---

### Automatic Detection

The script automatically detects which organization you're using.

**Both produce identical results!** Use whichever is more intuitive.

**Tip:** Always provide `netcdf_file_sets` to label your scenarios clearly.

**Note:** Ensemble plots currently limited to 2 file sets.

---

## ELM vs EAM Differences

### ELM (Land Model)

**Grid:** Structured latitude-longitude grid

**Configuration:**
```json
{
    "netcdf_files": "./elm_spatial_data.nc",
    "plot_directory": "./plots/elm/"
}
```

**No grid_file needed** - uses native lat/lon coordinates

**Variables:** GPP, NPP, NBP, ER, TOTECOSYSC, etc.

**Stippling:** ✓ **Available** for per-gridcell significance testing

---

### EAM (Atmosphere Model)

**Grid:** Unstructured spectral element mesh (ne30, ne120, etc.)

**Configuration:**
```json
{
    "netcdf_files": "./eam_spatial_data.nc",
    "plot_directory": "./plots/eam/",
    "grid_file": "./ne30pg2_scrip_c20191218.nc"
}
```

**Requires grid_file** - SCRIP format grid description

**Variables:** PRECC, PRECL, TREFHT, ZCO2, etc.

**Stippling:** ✗ **Not available** for per-gridcell tests (too slow for unstructured grids)

**Grid Files:**
- ne30: `ne30pg2_scrip_c20191218.nc`
- ne120: `ne120pg2_scrip_c20200803.nc`

---

## Comprehensive JSON Examples

### ELM Examples

#### Example 1: Single Control Map (ELM Config #1)

**From JSON file:**
```json
{
    "netcdf_files": "./control_spatial_data_elm.nc",
    "plot_directory": "./spatial_plots/elm_control",
    "use_latex": true
}
```

**What This Creates:**

**Individual spatial maps for all variables:**
- `spatial_GPP.pdf` - Control GPP map (2071-2090 mean)
- `spatial_NPP.pdf` - Control NPP map
- `spatial_NBP.pdf` - Control NBP map
- ... (all variables in NetCDF file)

**Statistics panel on each map:**
```
Max: 1.23e+01
Mean: 4.56e+00
Median: 3.89e+00
Min: -2.34e-01
Std: 2.67e+00
```

**No statistical testing** (only one file)

---

#### Example 2: Single Scenario Map (ELM Config #2)

**From JSON file:**
```json
{
    "netcdf_files": "./full_feedback_spatial_data_elm.nc",
    "plot_directory": "./spatial_plots/elm_full_feedback",
    "use_latex": true
}
```

**What This Creates:**

Same as Example 1, but for Full feedback scenario

---

#### Example 3: Individual Absolute Difference (ELM Config #3)

**From JSON file:**
```json
{
    "netcdf_files": [
        "./control_spatial_data_elm.nc", 
        "./full_feedback_spatial_data_elm.nc"
    ],
    "plot_directory": "./spatial_plots/elm_full_feedback_minus_control_individual_absolute_difference",
    "use_latex": true
}
```

**What This Creates:**

**Difference maps (Full feedback - Control):**
- Positive values (red): GPP increased in Full feedback
- Negative values (blue): GPP decreased in Full feedback
- White/near-zero: Little change

**Example for GPP:**
- Amazon: +5 to +10 gC/m²/month (strong increase)
- Sahara: ±0.1 gC/m²/month (minimal change)
- Northern Canada: +2 to +4 gC/m²/month (moderate increase)

**Statistical testing:** Aggregated t-test in `p_values.dat`

**No stippling** (simple list, not ensemble)

---

#### Example 4: Individual Percent Difference (ELM Config #4)

**From JSON file:**
```json
{
    "netcdf_files": [
        "./control_spatial_data_elm.nc", 
        "./full_feedback_spatial_data_elm.nc"
    ],
    "plot_directory": "./spatial_plots/elm_full_feedback_minus_control_individual_percent_difference",
    "use_latex": true,
    "plot_type": "percent_difference"
}
```

**What This Creates:**

**Percent change maps:**
- Amazon: +25% to +40% (large relative increase)
- Boreal forests: +10% to +15% (moderate relative increase)

**Statistical testing:** Aggregated only

---

#### Example 5: Mean Across Multiple Files (ELM Config #5)

**From JSON file:**
```json
{
    "netcdf_files": [
        "./full_feedback_spatial_data_elm.nc",
        "./full_feedback_spatial_data_elm_2.nc",
        "./full_feedback_spatial_data_elm_3.nc",
        "./full_feedback_spatial_data_elm_4.nc",
        "./full_feedback_spatial_data_elm_5.nc"
    ],
    "plot_directory": "./spatial_plots/elm_full_feedback_only_individual_mean",
    "variables": ["HR", "NPP", "ZCO2"],
    "use_latex": true,
    "plot_type": "mean"
}
```

**What This Creates:**

**Average maps across 5 ensemble members:**
- `spatial_HR.pdf` - Mean HR across 5 files
- `spatial_NPP.pdf` - Mean NPP across 5 files
- `spatial_ZCO2.pdf` - Mean ZCO2 across 5 files

**Processing:**
```
For each gridcell:
  Mean_GPP = (GPP_file1 + GPP_file2 + GPP_file3 + GPP_file4 + GPP_file5) / 5
```

**No statistical testing** (no comparison being made)

**Use case:** 
- Visualizing ensemble mean without needing control comparison
- Showing average pattern across realizations
- Publication plot of ensemble mean alone

---

#### Example 6: Ensemble Absolute Difference (ELM Config #6)

**From JSON file:**
```json
{
    "netcdf_files": [
        ["./control.nc", "./control_2.nc", "./control_3.nc", "./control_4.nc", "./control_5.nc"],
        ["./feedback.nc", "./feedback_2.nc", "./feedback_3.nc", "./feedback_4.nc", "./feedback_5.nc"]
    ],
    "netcdf_file_sets": ["Control", "Full feedback"],
    "plot_directory": "./spatial_plots/elm_full_feedback_minus_control_ensemble_absolute_difference",
    "plot_type": "absolute_difference",
    "use_latex": true
}
```

**What This Creates:**

**Ensemble mean difference maps:**
- Control ensemble mean: 120 ± 5 gC/m²/month (varies by location)
- Full feedback ensemble mean: 128 ± 6 gC/m²/month
- Difference map shows: +8 gC/m²/month

**Statistical testing:**
1. **Aggregated:** Global t-test in `p_values.dat`
2. **Per-gridcell:** NO stippling (stippling_on not set)

---

#### Example 7: Ensemble Percent Difference with Stippling (ELM Config #7)

**From JSON file:**
```json
{
    "netcdf_files": [
        ["./control.nc", "./control_2.nc", "./control_3.nc", "./control_4.nc", "./control_5.nc"],
        ["./feedback.nc", "./feedback_2.nc", "./feedback_3.nc", "./feedback_4.nc", "./feedback_5.nc"]
    ],
    "netcdf_file_sets": ["Control", "Full feedback"],
    "plot_directory": "./spatial_plots/elm_full_feedback_minus_control_ensemble_percent_difference",
    "plot_type": "percent_difference",
    "use_latex": true,
    "stippling_on": true
}
```

**What This Creates:**

**Ensemble mean percent difference with significance:**
- Shows % change with **stippling (xxxx pattern)** where p < 0.05

**Example for GPP:**
- Amazon basin: +35% change, **dense stippling** (highly significant)
- Sahara: +5% change, no stippling (not significant)
- Boreal: +15% change, **moderate stippling** (some gridcells significant)

**Interpretation:**
- Stippled regions: Changes are statistically robust at that location
- Non-stippled regions: Changes may be due to natural variability

---

#### Example 8: Ensemble Separate Plots with Stippling (ELM Config #8)

**From JSON file:**
```json
{
    "netcdf_files": [
        ["./control.nc", "./control_2.nc", "./control_3.nc", "./control_4.nc", "./control_5.nc"],
        ["./feedback.nc", "./feedback_2.nc", "./feedback_3.nc", "./feedback_4.nc", "./feedback_5.nc"]
    ],
    "netcdf_file_sets": ["Control", "Full feedback"],
    "plot_directory": "./spatial_plots/elm_full_feedback_minus_control_ensemble_separate_plots",
    "plot_type": "separate_plots",
    "use_latex": true,
    "stippling_on": true
}
```

**What This Creates:**

**Two separate ensemble mean maps:**
- `spatial_GPP_set_1.pdf` - Control ensemble mean (100-150 gC/m²/month)
- `spatial_GPP_set_2.pdf` - Full feedback ensemble mean (105-160 gC/m²/month)

**Stippling meaning:** Shows regions with high ensemble confidence

---

### EAM Examples

#### Example 9: Single Control Map (EAM Config #1)

**From JSON file:**
```json
{
    "netcdf_files": "./control_spatial_data_eam.nc",
    "plot_directory": "./spatial_plots/eam_control",
    "use_latex": true,
    "grid_file": "./eam_grid_file/ne30pg2_scrip_c20191218.nc"
}
```

**What This Creates:**

**Maps for all EAM variables:**
- `spatial_TREFHT.pdf` - Surface temperature
- `spatial_PRECC.pdf` - Convective precipitation
- `spatial_PRECL.pdf` - Large-scale precipitation
- `spatial_ZCO2.pdf` - CO₂ concentration

**Key:** `grid_file` parameter **required** for EAM

---

#### Example 10: Single Scenario Map (EAM Config #2)

Same as Example 9, but for Full feedback scenario

---

#### Example 11: Individual Absolute Difference (EAM Config #3)

**From JSON file:**
```json
{
    "netcdf_files": [
        "./control_spatial_data_eam.nc", 
        "./full_feedback_spatial_data_eam.nc"
    ],
    "plot_directory": "./spatial_plots/eam_full_feedback_minus_control_individual_absolute_difference",
    "use_latex": true,
    "grid_file": "./eam_grid_file/ne30pg2_scrip_c20191218.nc"
}
```

**What This Creates:**

**Difference maps for EAM variables:**
- Temperature differences
- Precipitation changes
- CO₂ concentration differences

---

#### Example 12: Individual Percent Difference (EAM Config #4)

**From JSON file:**
```json
{
    "netcdf_files": [
        "./control_spatial_data_eam.nc", 
        "./full_feedback_spatial_data_eam.nc"
    ],
    "plot_directory": "./spatial_plots/eam_full_feedback_minus_control_individual_percent_difference",
    "use_latex": true,
    "plot_type": "percent_difference",
    "grid_file": "./eam_grid_file/ne30pg2_scrip_c20191218.nc"
}
```

**What This Creates:**

**Percent change maps:**
- TREFHT: +1% to +5% (temperature)
- PRECC: +10% to +30% (convective precipitation intensification)

---

#### Example 13: Sum Across Multiple Files (EAM Config #5)

**From JSON file:**
```json
{
    "netcdf_files": [
        "./full_feedback_spatial_data_eam.nc",
        "./full_feedback_spatial_data_eam_2.nc",
        "./full_feedback_spatial_data_eam_3.nc",
        "./full_feedback_spatial_data_eam_4.nc",
        "./full_feedback_spatial_data_eam_5.nc"
    ],
    "plot_directory": "./spatial_plots/eam_full_feedback_only_individual_sum",
    "variables": ["PRECIP", "PRECSL"],
    "use_latex": true,
    "plot_type": "sum",
    "grid_file": "./eam_grid_file/ne30pg2_scrip_c20191218.nc"
}
```

**What This Creates:**

**Summed maps across 5 files:**
- `spatial_PRECIP.pdf` - Sum of PRECIP across 5 files
- `spatial_PRECSL.pdf` - Sum of PRECSL across 5 files

**Processing:**
```
For each gridcell:
  Sum_PRECIP = PRECIP_file1 + PRECIP_file2 + PRECIP_file3 + PRECIP_file4 + PRECIP_file5
```

**Use case:**
- Combining precipitation components
- Total fluxes across ensemble members
- Aggregating related variables

**No statistical testing** (no comparison)

---

#### Example 14: Ensemble Absolute Difference (EAM Config #6)

**From JSON file:**
```json
{
    "netcdf_files": [
        ["./control_eam.nc", "./control_eam_2.nc", "./control_eam_3.nc", "./control_eam_4.nc", "./control_eam_5.nc"],
        ["./feedback_eam.nc", "./feedback_eam_2.nc", "./feedback_eam_3.nc", "./feedback_eam_4.nc", "./feedback_eam_5.nc"]
    ],
    "netcdf_file_sets": ["Control", "Full feedback"],
    "plot_directory": "./spatial_plots/eam_full_feedback_minus_control_ensemble_absolute_difference",
    "plot_type": "absolute_difference",
    "use_latex": true,
    "grid_file": "./eam_grid_file/ne30pg2_scrip_c20191218.nc"
}
```

**What This Creates:**

**Ensemble mean difference maps for EAM:**
- TREFHT difference: +1.7 K average
  - Arctic: +3 to +5 K (strong amplification)
  - Tropics: +1 to +2 K

**Statistical testing:** Aggregated only (no stippling for EAM)

---

#### Example 15: Ensemble Percent Difference (EAM Config #7)

**From JSON file:**
```json
{
    "netcdf_files": [
        ["./control_eam.nc", "./control_eam_2.nc", "./control_eam_3.nc", "./control_eam_4.nc", "./control_eam_5.nc"],
        ["./feedback_eam.nc", "./feedback_eam_2.nc", "./feedback_eam_3.nc", "./feedback_eam_4.nc", "./feedback_eam_5.nc"]
    ],
    "netcdf_file_sets": ["Control", "Full feedback"],
    "plot_directory": "./spatial_plots/eam_full_feedback_minus_control_ensemble_percent_difference",
    "plot_type": "percent_difference",
    "use_latex": true,
    "stippling_on": false,
    "grid_file": "./eam_grid_file/ne30pg2_scrip_c20191218.nc"
}
```

**Key:** `stippling_on: false` (EAM limitation)

**What This Creates:**

**Percent difference without stippling:**
- PRECC: +40% to +60% (tropical convection)
- Mediterranean: -20% to -30% (drying)

---

#### Example 16: Ensemble Separate Plots (EAM Config #8)

**From JSON file:**
```json
{
    "netcdf_files": [
        ["./control_eam.nc", "./control_eam_2.nc", "./control_eam_3.nc", "./control_eam_4.nc", "./control_eam_5.nc"],
        ["./feedback_eam.nc", "./feedback_eam_2.nc", "./feedback_eam_3.nc", "./feedback_eam_4.nc", "./feedback_eam_5.nc"]
    ],
    "netcdf_file_sets": ["Control", "Full feedback"],
    "plot_directory": "./spatial_plots/eam_full_feedback_minus_control_ensemble_separate_plots",
    "plot_type": "separate_plots",
    "use_latex": true,
    "stippling_on": false,
    "grid_file": "./eam_grid_file/ne30pg2_scrip_c20191218.nc"
}
```

**What This Creates:**

**Two separate ensemble mean maps:**
- `spatial_TREFHT_set_1.pdf` - Control temperature (288-290 K)
- `spatial_TREFHT_set_2.pdf` - Full feedback temperature (289-291 K)

---

## Output Files

### Plot Files

**Naming Convention:**
```
{plot_directory}/spatial_{variable}.pdf
{plot_directory}/spatial_{variable}_set_1.pdf  (separate_plots)
{plot_directory}/spatial_{variable}_set_2.pdf  (separate_plots)
```

### P-Value Files

**Purpose:** Records aggregated statistical significance across all gridcells

**Format:**
```
GPP in Full feedback: 1.234e-05
NPP in Full feedback: 3.456e-03
```

**Location:** `{plot_directory}/p_values.dat`

**Note:** Contains **aggregated** p-values (all gridcells combined). Per-gridcell p-values shown as stippling on maps (ELM only).

---

## Best Practices

### 1. Choose Appropriate Plot Type

**Use `absolute_difference` when:**
- Comparing two scenarios/ensembles
- Want actual magnitude of changes
- Publishing quantitative results

**Use `percent_difference` when:**
- Emphasizing relative changes
- Variables with varying baseline magnitudes

**Use `separate_plots` when:**
- Need to show absolute values for each scenario
- Side-by-side comparison

**Use `mean` when:**
- Multiple files from same scenario
- Want average spatial pattern
- Reducing noise

**Use `sum` when:**
- Aggregating related variables (e.g., precipitation components)
- Total fluxes across files

---

### 2. Individual vs Ensemble

**Use individual (simple list) when:**
- Comparing 2 files: absolute or percent difference
- Averaging N files: mean or sum
- N separate maps: separate_plots

**Use ensemble (nested list) when:**
- Need statistical significance testing
- Per-gridcell confidence via stippling (ELM)
- Robust ensemble mean comparison

**Note:** Ensembles currently limited to 2 file sets

---

### 3. Temporal Averaging Period

**Standard:** 20-year average (2071-2090)
```json
"start_year": 2071,
"end_year": 2090
```

**Why average?** Reduces interannual variability

---

### 4. Stippling Strategy

**For ELM ensemble plots:**
```json
"stippling_on": true
```

**Benefits:** Shows where changes are statistically significant

**For EAM or individual plots:**
```json
"stippling_on": false
```

**Why:** EAM too slow, individuals lack ensemble structure

---

### 5. Variable Selection

**Plot all:**
```json
"variables": "all"
```

**Plot specific:**
```json
"variables": ["GPP", "NPP", "NBP"]
```

---

## Troubleshooting

### Issue 1: Grid File Missing (EAM)

**Error:**
```
TypeError: 'NoneType' object is not callable
```

**Cause:** Forgot `grid_file` parameter for EAM data

**Solution:**
```json
{
    "netcdf_files": "./eam_data.nc",
    "grid_file": "./ne30pg2_scrip_c20191218.nc"
}
```

---

### Issue 2: Stippling Not Appearing (EAM)

**Cause:** Stippling disabled for EAM (computational limitation)

**Solution:** Accept limitation or use ELM data

---

### Issue 3: Stippling Not Appearing (ELM)

**Possible Causes:**
1. `stippling_on: false` (default)
2. Simple list (not ensemble)
3. All gridcells non-significant

**Solutions:**
```json
"stippling_on": true  // Must enable
```

Verify nested list structure

---

### Issue 4: Wrong Plot Type for Configuration

**Symptom:** Error or unexpected output

**Cause:** Using incompatible plot_type

**Rules:**
- `absolute_difference`, `percent_difference`: Need exactly 2 files/ensembles
- `mean`, `sum`: Simple list only (not nested)
- `separate_plots`: Works with any configuration

---

## References

- E3SM Documentation: [https://e3sm.org/](https://e3sm.org/)
- Cartopy Documentation: [https://scitools.org.uk/cartopy/](https://scitools.org.uk/cartopy/)
- xarray Documentation: [https://xarray.pydata.org/](https://xarray.pydata.org/)
- uxarray Documentation: [https://uxarray.readthedocs.io/](https://uxarray.readthedocs.io/)
- GitHub Repository: [https://github.com/philipmyint/e3sm--gcam_analysis](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Contact

For questions or issues:
- **Philip Myint**: myint1@llnl.gov
- **Dalei Hao**: dalei.hao@pnnl.gov

---

## Version Information

**Script:** e3sm_plot_spatial_data.py  
**Documentation Version:** 2.0  
**Last Updated:** January 2026  
**Python Version:** 3.7+  
**Dependencies:** xarray, uxarray, cartopy, matplotlib, scipy, numpy, pandas

---

*This documentation provides comprehensive guidance for using the `e3sm_plot_spatial_data.py` script to create publication-quality spatial maps from E3SM model output with ensemble analysis and statistical significance testing.*
