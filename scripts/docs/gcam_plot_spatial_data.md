# GCAM Spatial Data Plotting Script Documentation


**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Overview

**Script Name:** `gcam_plot_spatial_data.py`

**Purpose:** Creates spatial (geographic) plots from GCAM (Global Change Analysis Model) output files, displaying data on maps using GCAM region/basin boundaries. The script supports individual scenario plots, ensemble means, absolute differences, and percent differences between scenarios, with optional statistical significance stippling.

**Author:** Philip Myint (myint1@llnl.gov) and Dalei Hao (dalei.hao@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Requirements](#installation--requirements)
3. [Basic Usage](#basic-usage)
4. [Plot Types](#plot-types)
5. [Complete Parameter Reference Table](#complete-parameter-reference-table)
6. [Detailed Parameter Descriptions](#detailed-parameter-descriptions)
7. [JSON Configuration Examples](#json-configuration-examples)
8. [Statistical Significance Stippling](#statistical-significance-stippling)
9. [Colormap and Colorbar](#colormap-and-colorbar)
10. [Shape Files and Geographic Data](#shape-files-and-geographic-data)
11. [Output Files](#output-files)
12. [Troubleshooting](#troubleshooting)
13. [Best Practices](#best-practices)

---

## Installation & Requirements

### Required Python Libraries

```bash
pip install pandas numpy scipy matplotlib geopandas
```

### Required Utility Modules

The script imports several utility modules that must be in the same directory or Python path:
- `utility_constants` - Default constants for plotting
- `utility_dataframes` - Functions for reading files and performing t-tests
- `utility_functions` - General utility functions
- `utility_gcam` - GCAM-specific functions for landtype grouping
- `utility_plots` - Plotting utility functions (provides default values)

### System Requirements

- Python 3.7+
- GeoPandas for handling shapefiles
- GCAM region/basin boundary shapefiles
- LaTeX installation (optional, for publication-quality typography)

---

## Basic Usage

### Command Line Execution

```bash
python gcam_plot_spatial_data.py path/to/config.json
```

**Multiple Configuration Files:**
```bash
python gcam_plot_spatial_data.py config1.json config2.json config3.json
```

### What the Script Does

For each configuration block in the JSON file, the script:
1. **Reads** processed GCAM output CSV files
2. **Reads** GCAM region/basin boundary shapefiles
3. **Aggregates** data across time (mean or sum over years)
4. **Aggregates** data spatially (for landtype groups, regions/basins)
5. **Joins** data with geographic boundaries
6. **Creates** spatial maps with color-coded values
7. **Performs** statistical tests (t-tests) when comparing scenarios
8. **Adds** stippling to indicate statistical significance
9. **Generates** publication-quality PDF (and optionally PNG) figures

---

## Plot Types

The script supports four main plot types:

### 1. Single Scenario (`plot_type` not specified or single scenario)

Shows absolute values for one scenario.

**Use for:**
- Displaying baseline conditions
- Showing spatial patterns of a single scenario
- Initial data exploration

**Example:** Control scenario commodity prices across GCAM regions

---

### 2. Mean (`plot_type: "mean"`)

Shows the mean across multiple scenarios (ensemble mean).

**Use for:**
- Displaying ensemble average
- Reducing noise from individual ensemble members
- Showing central tendency of ensemble

**Example:** Mean forest area across 5 ensemble members

**Requirements:**
- Multiple scenarios in simple list format
- All scenarios averaged together

---

### 3. Absolute Difference (`plot_type: "absolute_difference"` - default)

Shows the difference between scenarios: Scenario 2 - Scenario 1 (or Ensemble 2 mean - Ensemble 1 mean).

**Use for:**
- Comparing two scenarios directly
- Showing magnitude of changes
- Identifying regions with largest differences

**Example:** Full feedback emissions - Control emissions

**Color Interpretation:**
- Positive values (red/orange): Scenario 2 > Scenario 1
- Negative values (blue): Scenario 2 < Scenario 1
- Zero (white/yellow): No difference

---

### 4. Percent Difference (`plot_type: "percent_difference"`)

Shows the percent change: (Scenario 2 - Scenario 1) / Scenario 1 × 100

**Use for:**
- Comparing relative changes
- Normalizing differences across regions with different baseline values
- Highlighting proportional impacts

**Example:** Percent change in land area from Control to Full feedback

**Color Interpretation:**
- Positive values: Percent increase
- Negative values: Percent decrease
- Zero: No change

**Formula:**
```
percent_difference = ((value_2 - value_1) / value_1) × 100
```

---

## Complete Parameter Reference Table

### Required Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `shape_file` | string | **Yes** | Path to GCAM region/basin boundary shapefile |
| `output_file` | string | **Yes** | Path to processed CSV file containing GCAM data |
| `scenarios` | list or nested list | **Yes** | Scenario names (list) or ensemble members (nested list) |

### Data Selection Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `categories` | list | No | All categories | List of category names | Specific categories to plot |
| `categories_to_exclude` | list | No | `None` | List of category names | Categories to exclude |
| `category_label` | string | No | `"sector"` | Column name | Column identifying categories |
| `scenario_label` | string | No | `"scenario"` | Column name | Column identifying scenarios |
| `scenario_sets` | list | No | `None` | List of strings | Names for ensemble groups |
| `notify_scenarios_transposed` | boolean | No | `false` | `true`, `false` | Print console message when scenarios are automatically transposed |
| `year_label` | string | No | `"year"` | Column name | Column identifying years |
| `value_label` | string | No | `"value"` | Column name | Column with numeric data |
| `region_label` | string | No | `"region"` | Column name | Column identifying regions |
| `basin_label` | string | No | `"basin"` | Column name | Column identifying basins |
| `start_year` | integer | No | `2070` | Any year | Start year for time aggregation |
| `end_year` | integer | No | `2090` | Any year | End year for time aggregation |
| `key_columns` | list | No | `None` | List of column names | Columns for grouping |

### Geographic Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `shape_file` | string | **Yes** | - | Valid shapefile path | GCAM boundary shapefile (.shp) |
| `shape_file_region_label` | string | No | `"reg_nm"` | Column name | Region name column in shapefile |
| `shape_file_basin_label` | string | No | `"glu_nm"` | Column name | Basin name column in shapefile |

### Aggregation Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `mean_or_sum_for_time_aggregation` | string | No | `"mean"` | `"mean"`, `"sum"` | How to aggregate across years |
| `mean_or_sum_if_more_than_one_row_in_same_landtype_group` | string | No | `"area_weighted_mean"` | `"mean"`, `"area_weighted_mean"`, `"sum"` | Aggregation for landtype groups |
| `mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin` | string | No | `"mean"` | `"mean"`, `"sum"` | Aggregation within region/basin |
| `landtype_groups` | string | No | `"modified"` | `"modified"`, `"original"` | Landtype grouping dictionary |
| `multiplier` | number | No | `1` | Any number | Multiply values (unit conversion) |

### Plot Type and Display Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `plot_type` | string | No | `"absolute_difference"` | `"mean"`, `"absolute_difference"`, `"percent_difference"` | Type of spatial plot |
| `plot_directory` | string | No | `"./"` | Valid directory path | Output directory |
| `plot_name` | string | No | Auto-generated | Filename with .pdf | Plot filename |
| `title` | string | No | Auto-generated | Any string | Plot title |

### Visual Styling Parameters (from utility_plots.py)

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `width` | number | No | `10` inches | Positive number | Figure width |
| `height` | number | No | `8` inches | Positive number | Figure height |
| `use_latex` | boolean | No | `false` | `true`, `false` | LaTeX text rendering |
| `produce_png` | boolean | No | `false` | `true`, `false` | Also create PNG |
| `title_size` | number | No | `24` points | Positive number | Plot title font size |
| `x_tick_label_size` | number | No | `20` points | Positive number | X-axis tick font size |
| `y_tick_label_size` | number | No | `20` points | Positive number | Y-axis tick font size |

### Colormap and Colorbar Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `cmap` | string | No | `"viridis"` | Matplotlib colormap name | Color scheme for map |
| `cbar_on` | boolean | No | `true` | `true`, `false` | Show colorbar |
| `cbar_limits` | list | No | `None` | `[min, max]` | Colorbar value range |
| `linewidth` | number | No | `0.5` | Positive number | Width of boundary lines |

### Statistical Significance Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `stippling_on` | boolean | No | `true` | `true`, `false` | Show statistical significance |
| `stippling_hatches` | string | No | `"xxxx"` | Matplotlib hatch pattern | Stippling pattern |
| `p_value_threshold` | number | No | `0.05` | 0 to 1 | Significance threshold |
| `p_value_file` | string | No | `"p_values.dat"` | Filename | P-value output file |
| `p_value_file_print_only_if_below_threshold` | boolean | No | `true` | `true`, `false` | Only print significant values |

---

## Detailed Parameter Descriptions

### Core Required Parameters

#### `shape_file`
**Type:** String (file path)  
**Required:** Yes  
**Description:** Path to GCAM region/basin boundary shapefile.

**Examples:**
```json
"shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp"
```

**Requirements:**
- Must be a valid shapefile (.shp) with associated files (.shx, .dbf, .prj)
- Must contain GCAM region and/or basin boundaries
- Must have columns matching `shape_file_region_label` and `shape_file_basin_label`

**Common GCAM Shapefiles:**
- Combined region + basin boundaries: `reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp`
- Region-only boundaries: `reg_boundaries_moirai_3p1_0p5arcmin.shp`
- Basin-only boundaries: `glu_boundaries_moirai_3p1_0p5arcmin.shp`

---

#### `output_file`
**Type:** String (file path)  
**Required:** Yes  
**Description:** Path to processed CSV file containing GCAM output data.

**Examples:**
```json
"output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv"
"output_file": "./data/land_allocation_processed.csv"
```

**Requirements:**
- CSV format with headers
- Must have columns for: scenario, region/basin, year, category (if applicable), value
- Processed using `gcam_process_extracted_data.py`

---

#### `scenarios`
**Type:** List of strings OR nested list  
**Required:** Yes  
**Description:** Scenarios to plot or compare.

**Format 1 - Single or Multiple Scenarios (Simple List):**
```json
"scenarios": ["Control"]
"scenarios": ["Control", "Full feedback"]
"scenarios": ["Control", "Control_2", "Control_3"]
```

**Format 2 - Ensemble Comparison (Nested List):**

Users can specify ensemble scenarios in **either of two formats**. The script automatically detects and handles both formats using the `transpose_scenarios_if_needed()` function in `utility_functions.py`.

**Format 2A - Organized by Scenario Set (Recommended):**
```json
"scenarios": [
    ["Control", "Control_2", "Control_3"],
    ["Full feedback", "Full feedback_2", "Full feedback_3"]
],
"scenario_sets": ["Control", "Full feedback"]
```
Each inner list contains all ensemble members for one scenario set. This format is more intuitive as it groups related scenarios together.

**Format 2B - Organized by Ensemble Member (Alternative):**
```json
"scenarios": [
    ["Control", "Full feedback"],
    ["Control_2", "Full feedback_2"],
    ["Control_3", "Full feedback_3"]
],
"scenario_sets": ["Control", "Full feedback"]
```
Each inner list is one ensemble member.

**Behavior by Plot Type:**
- **Single scenario:** Shows absolute values
- **Multiple scenarios (list) + `plot_type: "mean"`:** Shows ensemble mean
- **Multiple scenarios (list):** Shows absolute or percent difference (2nd - 1st)
- **Nested list:** Compares ensemble means (Set 2 mean - Set 1 mean)

**Important Notes:**
- The `scenario_sets` parameter helps the script detect which format you're using
- Scenario names do not need to follow any specific naming convention

---

### Year Range Parameters

#### `start_year` and `end_year`
**Type:** Integer  
**Required:** No  
**Defaults:** `start_year: 2070`, `end_year: 2090`

**Description:** Define temporal range for aggregation. Data between these years is aggregated (mean or sum) to create the spatial plot.

**Example:**
```json
"start_year": 2080,
"end_year": 2100
```

**Usage:**
- Shorter range: More specific time period, less temporal averaging
- Longer range: More stable estimates, reduces year-to-year variability
- Single year: Set both to same value

---

### Aggregation Parameters

#### `mean_or_sum_for_time_aggregation`
**Type:** String  
**Required:** No  
**Default:** `"mean"`  
**Possible Values:** `"mean"`, `"sum"`

**Description:** How to aggregate data across the year range.

**Use Cases:**
- `"mean"` - Average across years (prices, scalars, rates)
- `"sum"` - Total across years (cumulative emissions, production)

**Example:**
```json
"mean_or_sum_for_time_aggregation": "mean"
```

---

#### `mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin`
**Type:** String  
**Required:** No  
**Default:** `"mean"`  
**Possible Values:** `"mean"`, `"sum"`

**Description:** How to aggregate when multiple data rows exist for the same geographic unit.

**Use Cases:**
- `"mean"` - Average values (e.g., average price across crops in a region)
- `"sum"` - Total values (e.g., total land area across all landtypes)

---

#### `mean_or_sum_if_more_than_one_row_in_same_landtype_group`
**Type:** String  
**Required:** No  
**Default:** `"area_weighted_mean"`  
**Possible Values:** `"mean"`, `"area_weighted_mean"`, `"sum"`

**Description:** How to aggregate when grouping landtypes (e.g., all crops into "crop" group).

**Use Cases:**
- `"area_weighted_mean"` - Weighted average by land area
- `"mean"` - Simple average
- `"sum"` - Total (for areas, production)

---

### Plot Type Parameter

#### `plot_type`
**Type:** String  
**Required:** No  
**Default:** `"absolute_difference"`  
**Possible Values:** `"mean"`, `"absolute_difference"`, `"percent_difference"`

**Description:** Type of spatial plot to create.

**`"mean"`** - Ensemble mean
- Requires: Multiple scenarios in simple list
- Shows: Average across all scenarios
- Use for: Central tendency of ensemble

**`"absolute_difference"`** - Absolute difference (default)
- Requires: 2 scenarios or 2 ensemble sets
- Shows: Scenario2 - Scenario1
- Use for: Magnitude of change

**`"percent_difference"`** - Percent difference
- Requires: 2 scenarios or 2 ensemble sets
- Shows: (Scenario2 - Scenario1) / Scenario1 × 100
- Use for: Relative/normalized change

**Example:**
```json
"plot_type": "percent_difference"
```

---

### Colormap Parameters

#### `cmap`
**Type:** String  
**Required:** No  
**Default:** `"viridis"`

**Description:** Matplotlib colormap for coloring the map.

**Common Colormaps:**

**Sequential (for single-valued data):**
- `"viridis"` - Blue to yellow (perceptually uniform, colorblind-friendly)
- `"plasma"` - Purple to yellow
- `"cividis"` - Blue to yellow (optimized for colorblind)
- `"Blues"`, `"Reds"`, `"Greens"` - Single hue gradients

**Diverging (for difference plots):**
- `"RdBu"` - Red (positive) to Blue (negative)
- `"RdBu_r"` - Reversed: Blue (positive) to Red (negative)
- `"RdYlBu"` - Red-Yellow-Blue
- `"BrBG"` - Brown (positive) to Blue-Green (negative)
- `"PiYG"` - Pink (positive) to Green (negative)

**Example:**
```json
"cmap": "RdBu_r"
```

**Best Practices:**
- Use diverging colormaps for difference/percent difference plots
- Use sequential colormaps for absolute value plots
- Consider colorblind accessibility

---

#### `cbar_limits`
**Type:** List of two numbers  
**Required:** No  
**Default:** `None` (auto-scaled)

**Description:** Manually set colorbar value range.

**Example:**
```json
"cbar_limits": [-20, 60]
```

**Usage:**
- Ensures consistent scaling across multiple plots
- Prevents outliers from dominating color scale
- Useful for comparison across categories
- For symmetric diverging: use [-X, X]

**Auto-scaling behavior:**
- Scales to data min/max
- May differ between plots
- Outliers can compress color range

---

#### `cbar_on`
**Type:** Boolean  
**Required:** No  
**Default:** `true`

**Description:** Whether to display colorbar.

**Example:**
```json
"cbar_on": false
```

**Use Cases:**
- `false` - For panel plots with shared colorbar
- `false` - When colorbar added separately
- `true` - Standard single plots

---

### Stippling Parameters

#### `stippling_on`
**Type:** Boolean  
**Required:** No  
**Default:** `true`

**Description:** Whether to show statistical significance stippling on difference plots.

**Example:**
```json
"stippling_on": false
```

**Behavior:**
- Automatically applied to absolute_difference and percent_difference plots
- Performs t-test for each region/basin
- Adds hatching pattern where p < threshold
- Only relevant when comparing scenarios/ensembles

---

#### `stippling_hatches`
**Type:** String  
**Required:** No  
**Default:** `"xxxx"`

**Description:** Matplotlib hatch pattern for stippling.

**Common Patterns:**
- `"xxxx"` - Dense X pattern (default)
- `"////"` - Diagonal lines
- `"\\\\\\\\"`- Reverse diagonal lines
- `"...."` - Dots
- `"++++"` - Plus signs
- `"||||"` - Vertical lines

**Example:**
```json
"stippling_hatches": "////"
```

---

### Shape File Label Parameters

#### `shape_file_region_label` and `shape_file_basin_label`
**Type:** String  
**Required:** No  
**Defaults:** `"reg_nm"` (region), `"glu_nm"` (basin)

**Description:** Column names in shapefile containing region and basin names.

**Example:**
```json
"shape_file_region_label": "region_name",
"shape_file_basin_label": "basin_name"
```

**Requirements:**
- Must match actual column names in shapefile
- Values must match region/basin names in CSV data
- Case-sensitive

**Verification:**
```python
import geopandas as gpd
gdf = gpd.read_file('your_shapefile.shp')
print(gdf.columns)
```

---

## JSON Configuration Examples

### Example 1: Single Scenario Map

**Purpose:** Display agricultural commodity prices for Control scenario

```json
{
    "shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
    "scenarios": ["Control"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/spatial_plots",
    "plot_name": "spatial_ag_commodity_control_only.pdf",
    "use_latex": true
}
```

**What it does:**
- Shows absolute values for Control scenario
- Creates separate maps for 5 crops
- Averages over default years (2070-2090)
- Uses viridis colormap (default)

---

### Example 2: Vegetation Scalars with Custom Labels

**Purpose:** Display EHC vegetation scalars with custom column names

```json
{
    "shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback.csv",
    "scenarios": ["Control"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/spatial_plots",
    "plot_name": "spatial_vegetation_scalars_control_only.pdf",
    "category_label": "landtype",
    "value_label": "vegetation",
    "title": "vegetation scalars",
    "use_latex": true
}
```

**What it does:**
- Uses "landtype" column instead of default "sector"
- Uses "vegetation" column instead of default "value"
- Custom plot title

---

### Example 3: Land Allocation with Landtype Groups

**Purpose:** Show forest area with landtype aggregation and summation

```json
{
    "shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed.csv",
    "scenarios": ["Control"],
    "categories": ["forest"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/spatial_plots",
    "plot_name": "spatial_land_allocation_forest_control_only.pdf",
    "category_label": "landtype",
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum",
    "mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin": "sum",
    "key_columns": ["scenario", "region", "basin", "year"],
    "title": "Forest area (thousands km$^2$)",
    "use_latex": true
}
```

**What it does:**
- Aggregates all forest types into "forest" group
- Sums areas within landtype group
- Sums across regions/basins
- LaTeX title with superscript

---

### Example 4: Ensemble Mean

**Purpose:** Show mean forest area across ensemble members

```json
{
    "shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_ensemble.csv",
    "scenarios": ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
    "categories": ["forest"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/spatial_plots",
    "plot_name": "spatial_land_allocation_forest_mean.pdf",
    "plot_type": "mean",
    "category_label": "landtype",
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum",
    "mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin": "sum",
    "key_columns": ["scenario", "region", "basin", "year"],
    "title": "Forest area (thousands km$^2$)",
    "use_latex": true
}
```

**What it does:**
- Averages forest area across 5 ensemble members
- Reduces ensemble variability
- Shows central tendency

---

### Example 5: Absolute Difference (Default)

**Purpose:** Show absolute difference in forest area between scenarios

```json
{
    "shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_ensemble.csv",
    "scenarios": ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
    "categories": ["forest"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/spatial_plots",
    "plot_name": "spatial_land_allocation_forest_absolute_difference.pdf",
    "category_label": "landtype",
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum",
    "mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin": "sum",
    "key_columns": ["scenario", "region", "basin", "year"],
    "title": "Forest area (thousands km$^2$)",
    "use_latex": true
}
```

**What it does:**
- Compares 2nd scenario minus 1st scenario (Control_2 - Control)
- Note: plot_type defaults to "absolute_difference"
- Shows where and how much forest area changes

---

### Example 6: Percent Difference

**Purpose:** Show percent change in forest area

```json
{
    "shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_ensemble.csv",
    "scenarios": ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
    "categories": ["forest"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/spatial_plots",
    "plot_name": "spatial_land_allocation_forest_percent_difference.pdf",
    "plot_type": "percent_difference",
    "category_label": "landtype",
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum",
    "mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin": "sum",
    "key_columns": ["scenario", "region", "basin", "year"],
    "title": "Forest area",
    "use_latex": true
}
```

**What it does:**
- Shows percent change: (Control_2 - Control) / Control × 100
- Title automatically gets "% difference" appended
- Normalizes by baseline values

---

### Example 7: Ensemble Comparison

**Purpose:** Compare ensemble means between two scenario groups

```json
{
    "shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_ensemble.csv",
    "scenarios": [
        ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
        ["Full feedback", "Full feedback_2", "Full feedback_3", "Full feedback_4", "Full feedback_5"]
    ],
    "scenario_sets": ["Control", "Full feedback"],
    "categories": ["forest"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/spatial_plots",
    "plot_name": "spatial_land_allocation_forest_ensemble_absolute_difference.pdf",
    "category_label": "landtype",
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum",
    "mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin": "sum",
    "key_columns": ["scenario", "region", "basin", "year"],
    "title": "Forest area (thousands km$^2$)",
    "use_latex": true
}
```

**What it does:**
- Calculates mean of Control ensemble (Control, Control_2, ...)
- Calculates mean of Full feedback ensemble (Full feedback, Full feedback_2, ...)
- Shows difference: Full feedback mean - Control mean
- Performs t-tests and adds stippling

---

### Example 8: Percent Difference with Custom Colorbar

**Purpose:** Show percent change in vegetation scalars with fixed color scale

```json
{
    "shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback_ensemble.csv",
    "scenarios": [
        ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
        ["Full feedback", "Full feedback_2", "Full feedback_3", "Full feedback_4", "Full feedback_5"]
    ],
    "scenario_sets": ["Control", "Full feedback"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/spatial_plots",
    "plot_name": "spatial_vegetation_scalars_ensemble_percent_difference.pdf",
    "plot_type": "percent_difference",
    "cbar_limits": [-20, 60],
    "category_label": "landtype",
    "value_label": "vegetation",
    "title": "vegetation scalars",
    "use_latex": true
}
```

**What it does:**
- Ensemble comparison with percent difference
- Fixed colorbar from -20% to +60%
- Ensures consistent scaling across crops
- Custom value column

---

### Example 9: No Stippling

**Purpose:** Show differences without statistical significance markers

```json
{
    "shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback_ensemble.csv",
    "scenarios": [
        ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
        ["Full feedback", "Full feedback_2", "Full feedback_3", "Full feedback_4", "Full feedback_5"]
    ],
    "scenario_sets": ["Control", "Full feedback"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/spatial_plots",
    "plot_name": "spatial_vegetation_scalars_ensemble_absolute_difference_no_stippling.pdf",
    "category_label": "landtype",
    "value_label": "vegetation",
    "title": "vegetation scalars",
    "stippling_on": false,
    "use_latex": true
}
```

**What it does:**
- Disables statistical significance stippling
- Cleaner visual appearance
- Useful when significance not needed or clutters plot

---

## Statistical Significance Stippling

### How It Works

**For Difference Plots:**
When comparing scenarios or ensembles, the script:
1. Performs independent t-test for each region/basin
2. Compares data from scenario 1 vs scenario 2 (or ensemble 1 vs ensemble 2)
3. Calculates p-value for each region/basin
4. Adds hatching pattern where p < threshold

**Sample Sizes:**
- Individual scenarios: n = number of years in range
- Ensemble comparison: n = number of ensemble members × number of years

**Statistical Test:**
```python
# For each region/basin
data_1 = values from scenario/ensemble 1
data_2 = values from scenario/ensemble 2
t_statistic, p_value = scipy.stats.ttest_ind(data_1, data_2)

if p_value < p_value_threshold:
    # Add stippling to this region
```

### Interpreting Stippling

**Stippled regions:** Difference is statistically significant (p < 0.05 default)
- Unlikely difference occurred by chance
- Strong evidence of real difference between scenarios

**Non-stippled regions:** Difference is not statistically significant
- Could be due to chance/variability
- Insufficient evidence of real difference
- Doesn't mean no difference exists

### Controlling Stippling

**Enable/Disable:**
```json
"stippling_on": true,  // Show stippling (default)
"stippling_on": false  // No stippling
```

**Change Threshold:**
```json
"p_value_threshold": 0.05,  // Standard
"p_value_threshold": 0.01,  // More stringent
"p_value_threshold": 0.1    // More lenient
```

**Change Pattern:**
```json
"stippling_hatches": "xxxx",  // Dense X (default)
"stippling_hatches": "////"   // Diagonal lines
```

### P-Value Output File

**Location:** `plot_directory/p_values.dat`

**Content:**
```
============================================================
Plot: ./plots/spatial_forest.pdf
------------------------------------------------------------
category=forest, region=USA: p=0.0234 *
category=forest, region=China: p=0.3421
category=forest, region=Brazil: p=0.0012 **
============================================================
```

---

## Colormap and Colorbar

### Choosing Colormaps

**For Single Scenario (Absolute Values):**
- Use sequential colormaps
- Examples: `"viridis"`, `"plasma"`, `"cividis"`, `"Blues"`
- Low to high progression

**For Difference Plots:**
- Use diverging colormaps
- Examples: `"RdBu_r"`, `"RdYlBu"`, `"BrBG"`, `"PiYG"`
- Center at zero
- Clearly distinguish positive vs negative

**Colorblind-Friendly:**
- `"viridis"`, `"plasma"`, `"cividis"` - Sequential
- `"RdBu"` - Diverging (but challenging)
- Avoid `"jet"` rainbow colormap

### Setting Colorbar Limits

**Auto-scaling (default):**
```json
"cbar_limits": null
```
- Scales to data min/max
- Different for each plot
- Outliers can compress range

**Manual limits:**
```json
"cbar_limits": [-100, 100]
```
- Consistent across plots
- Good for comparisons
- May clip outliers

**For symmetric diverging:**
```json
"cbar_limits": [-50, 50]  // Symmetric around zero
```

**For percent difference:**
```json
"cbar_limits": [-20, 60]  // Asymmetric if needed
```

---

## Shape Files and Geographic Data

### GCAM Boundary Shapefiles

**Structure:**
GCAM uses a hierarchical geographic structure:
- **Regions:** 32 political/economic regions (e.g., USA, China, EU-15)
- **Basins:** ~235 water basins (GLUs - Geographic Land Units)
- **Region-Basin combinations:** Intersections of regions and basins

**Required Files:**
A shapefile consists of multiple files:
- `.shp` - Feature geometry
- `.shx` - Shape index
- `.dbf` - Attribute data
- `.prj` - Projection information

All must be present in the same directory.

### Shapefile Attributes

**Typical Columns:**
- `reg_nm` - Region name (default for `shape_file_region_label`)
- `glu_nm` - Basin name (default for `shape_file_basin_label`)
- `reg_id` - Region ID number
- `glu_id` - Basin ID number
- `geometry` - Geographic boundaries

**Verify Attributes:**
```python
import geopandas as gpd
gdf = gpd.read_file('your_shapefile.shp')
print(gdf.columns)
print(gdf[['reg_nm', 'glu_nm']].head())
```

### Joining Data with Shapefile

The script automatically:
1. Reads shapefile
2. Aggregates CSV data by region/basin
3. Joins data to shapefile using region/basin names
4. Handles missing data (regions with no data shown in gray/white)

**Requirements:**
- Region/basin names in CSV must match shapefile
- Case-sensitive matching
- Whitespace matters

---

## Output Files

### Plot Files

**PDF Format (default):**
- Vector graphics (scalable)
- High resolution
- Small file size
- Publication quality

**PNG Format (optional via `produce_png: true`):**
- Raster graphics
- Fixed resolution
- Larger file size
- Web/presentation use

**File Naming:**
- Default: `spatial_[output_file_name].pdf`
- Custom: Specified by `plot_name`

### P-Value Files

**Location:** `plot_directory/p_values.dat`

**Generated for:** Difference and percent difference plots with multiple scenarios/ensembles

**Content:**
- Plot filename
- P-values for each region/basin
- Significance markers (*)

---

## Troubleshooting

### Issue 1: Shapefile Not Found

**Error:** `FileNotFoundError` or `DriverError`

**Solutions:**
- Verify shapefile path is correct
- Ensure all shapefile components (.shp, .shx, .dbf, .prj) are present
- Use absolute paths if unsure
- Check file permissions

---

### Issue 2: No Data on Map

**Symptom:** Map shows boundaries but no colors/all white

**Causes:**
- Region/basin names don't match between CSV and shapefile
- No data after aggregation
- All values are NaN

**Solutions:**
```python
# Check CSV region names
import pandas as pd
df = pd.read_csv('your_file.csv')
print(df['region'].unique())

# Check shapefile region names
import geopandas as gpd
gdf = gpd.read_file('your_shapefile.shp')
print(gdf['reg_nm'].unique())

# Compare - must match exactly
```

---

### Issue 3: Colorbar Range Issues

**Symptom:** All regions same color or very compressed range

**Cause:** Outliers or inappropriate auto-scaling

**Solutions:**
```json
// Set manual limits
"cbar_limits": [-100, 100]

// Or remove outliers first in data processing
```

---

### Issue 4: Stippling Too Dense

**Symptom:** Hatching pattern overwhelming or obscuring colors

**Solutions:**
```json
// Use lighter pattern
"stippling_hatches": ".."  // Sparse dots

// Or disable
"stippling_on": false

// Or increase threshold (fewer significant regions)
"p_value_threshold": 0.01
```

---

### Issue 5: LaTeX Errors

**Error:** `RuntimeError: Failed to process string with tex`

**Solutions:**
- Install LaTeX or set `"use_latex": false`
- Check title for special characters
- Escape LaTeX special characters

**Example:**
```json
// Wrong
"title": "Area (km^2)"

// Correct
"title": "Area (thousands km$^2$)"
```

---

### Issue 6: Memory Issues with Large Shapefiles

**Symptom:** Slow performance or out of memory errors

**Solutions:**
- Use simplified shapefiles (lower resolution)
- Process fewer categories at once
- Reduce figure size
- Close other applications

---

## Best Practices

### 1. Colormap Selection

**Guidelines:**
- **Sequential for absolute values:** Single scenario, ensemble mean
- **Diverging for differences:** Absolute difference, percent difference
- **Consistent colors:** Use same colormap across related figures
- **Accessibility:** Choose colorblind-friendly options

**Recommendations:**
```json
// Absolute values
"cmap": "viridis"  // Default, good choice

// Positive differences expected
"cmap": "YlOrRd"  // Yellow to red

// Differences (positive and negative)
"cmap": "RdBu_r"  // Red (positive) to Blue (negative)

// Percent differences
"cmap": "RdYlGn_r"  // Red (negative) to Green (positive)
```

---

### 2. Time Range Selection

**Consider:**
- **Longer range (20-30 years):** More stable, averages out variability
- **Shorter range (5-10 years):** More specific, but may be noisier
- **End-of-century (2080-2100):** Common for long-term projections
- **Mid-century (2040-2060):** Transition period

**Example:**
```json
// Long-term average
"start_year": 2070,
"end_year": 2100

// Specific decade
"start_year": 2090,
"end_year": 2100

// Single year
"start_year": 2100,
"end_year": 2100
```

---

### 3. Statistical Rigor

**For Ensemble Comparisons:**
- Use adequate ensemble size (5+ members preferred)
- Be aware of multiple testing issues
- Consider adjusted thresholds

**Bonferroni Correction:**
```json
// If testing 20 regions
"p_value_threshold": 0.0025  // 0.05 / 20
```

**Reporting:**
- Always report threshold used
- Report all comparisons, not just significant
- Consider effect size alongside significance

---

### 4. Figure Quality

**For Publications:**
```json
{
    "width": 10,
    "height": 8,
    "use_latex": true,
    "produce_png": false,
    "linewidth": 0.3
}
```

**For Presentations:**
```json
{
    "width": 12,
    "height": 8,
    "use_latex": false,
    "produce_png": true,
    "linewidth": 0.5,
    "stippling_on": false  // Cleaner for projector
}
```

**For Posters:**
```json
{
    "width": 14,
    "height": 10,
    "use_latex": false,
    "produce_png": true,
    "linewidth": 0.8
}
```

---

### 5. Organization

**Directory Structure:**
```
project/
├── data/
│   ├── processed/
│   └── shapefiles/
├── plots/
│   ├── spatial/
│   │   ├── absolute/
│   │   ├── percent/
│   │   └── means/
│   └── statistics/
└── configs/
    └── spatial_configs.json
```

**File Naming:**
```
spatial_[variable]_[scenario]_[plottype].pdf

Examples:
spatial_forest_control_only.pdf
spatial_emissions_ensemble_percent_difference.pdf
spatial_scalars_mean.pdf
```

---

### 6. Reproducibility

**Document Settings:**
- Save configuration files
- Note shapefile version/source
- Record data processing steps
- Version control configs

**Configuration Comments:**
Use external documentation since JSON doesn't support comments:
```
config_notes.txt:
- spatial_forest.json: Forest area comparison for paper Figure 3
- Used Moirai 3.1 boundaries at 0.5 arcmin resolution
- Years 2070-2090 average for end-of-century projection
- RdBu_r colormap with symmetric ±100 limits
```

---

## Integration with Other Scripts

### Typical Workflow

```
1. gcam_extract_csv_from_project_files.R
   ↓ (Extract from GCAM project files)
   
2. gcam_process_extracted_data.py
   ↓ (Clean and organize)
   
3. gcam_add_areas_to_files.py (optional)
   ↓ (Add land area data)
   
4. gcam_plot_spatial_data.py  ← THIS SCRIPT
   ↓ (Create spatial maps)
   
5. Analysis and publication
```

### Complementary Visualizations

**Together with:**
- **Time series plots** (`gcam_plot_time_series.py`) - Temporal evolution
- **Box plots** (`gcam_plot_box_and_whiskers.py`) - Distributions
- **Spatial plots** (this script) - Geographic patterns

**Provides:**
- Temporal: How does it change over time?
- Statistical: What's the distribution/variability?
- Spatial: Where do changes occur?

---

## Appendix: Complete Default Values

### From gcam_plot_spatial_data.py

```python
{
    'basin_label': 'basin',
    'category_label': 'sector',
    'cbar_limits': None,
    'cbar_on': True,
    'cmap': 'viridis',
    'end_year': 2090,
    'height': 8,  # inches (from utility_plots.py)
    'key_columns': None,
    'landtype_groups': 'modified',
    'linewidth': 0.5,
    'multiplier': 1,
    'mean_or_sum_for_time_aggregation': 'mean',
    'mean_or_sum_if_more_than_one_row_in_same_landtype_group': 'area_weighted_mean',
    'mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin': 'mean',
    'p_value_file': 'p_values.dat',
    'p_value_file_print_only_if_below_threshold': True,
    'p_value_threshold': 0.05,
    'plot_directory': './',
    'plot_type': 'absolute_difference',
    'produce_png': False,  # from utility_plots.py
    'region_label': 'region',
    'scenario_label': 'scenario',
    'scenario_sets': None,
    'shape_file_basin_label': 'glu_nm',
    'shape_file_region_label': 'reg_nm',
    'start_year': 2070,
    'stippling_hatches': 'xxxx',
    'stippling_on': True,
    'title_size': 24,  # points (from utility_plots.py)
    'use_latex': False,  # from utility_plots.py
    'value_label': 'value',
    'width': 10,  # inches (from utility_plots.py)
    'x_tick_label_size': 20,  # points (from utility_plots.py)
    'y_tick_label_size': 20,  # points (from utility_plots.py)
    'year_label': 'year'
}
```

---

## References

- DiVittorio et al. (2025). "E3SM-GCAM coupling methodology and applications." *Journal of Advances in Modeling Earth Systems*. [DOI: 10.1029/2024MS004806](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024MS004806)
- GitHub Repository: [https://github.com/philipmyint/e3sm--gcam_analysis](https://github.com/philipmyint/e3sm--gcam_analysis)
- GeoPandas Documentation: [https://geopandas.org/](https://geopandas.org/)
- Matplotlib Colormaps: [https://matplotlib.org/stable/gallery/color/colormap_reference.html](https://matplotlib.org/stable/gallery/color/colormap_reference.html)

---

## Contact

For questions or issues:
- **Philip Myint**: myint1@llnl.gov
- **Dalei Hao**: dalei.hao@pnnl.gov

---

## Version Information

**Script:** gcam_plot_spatial_data.py  
**Utility Module:** utility_plots.py  
**Documentation Version:** 1.0  
**Last Updated:** January 2026  
**Python Version:** 3.7+  
**GeoPandas Version:** 0.8+

---

*This documentation provides comprehensive guidance for using the `gcam_plot_spatial_data.py` script to create geographic visualizations of GCAM output data with statistical analysis and publication-quality cartographic displays.*
