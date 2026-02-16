# GCAM Box-and-Whisker Plotting Script Documentation

## Overview

**Script Name:** `gcam_plot_box_and_whiskers.py`

**Purpose:** Creates box-and-whisker plots (box plots) from GCAM (Global Change Analysis Model) output files to visualize data distributions across different scenarios, categories, and regions. Box plots provide a statistical summary showing median, quartiles, and outliers, making them ideal for comparing distributions and identifying variability in GCAM outputs.

**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Requirements](#installation--requirements)
3. [Basic Usage](#basic-usage)
4. [Understanding Box Plots](#understanding-box-plots)
5. [Plot Types](#plot-types)
6. [Configuration Parameters](#configuration-parameters)
7. [Complete Parameter Reference Table](#complete-parameter-reference-table)
8. [Detailed Parameter Descriptions](#detailed-parameter-descriptions)
9. [JSON Configuration Examples](#json-configuration-examples)
10. [The x_variable and hue System](#the-x_variable-and-hue-system)
11. [Output Files](#output-files)
12. [Advanced Features](#advanced-features)
13. [Troubleshooting](#troubleshooting)
14. [Best Practices](#best-practices)

---

## Installation & Requirements

### Required Python Libraries

```bash
pip install pandas numpy scipy matplotlib seaborn
```

### Required Utility Modules

The script imports several utility modules that must be in the same directory or Python path:
- `utility_constants` - Default constants for plotting
- `utility_dataframes` - Functions for reading files
- `utility_functions` - General utility functions
- `utility_gcam` - GCAM-specific functions for landtype grouping
- `utility_plots` - Plotting utility functions

### System Requirements

- Python 3.7+
- Multi-core processor (script uses parallel processing for multiple plots)
- LaTeX installation (optional, for publication-quality typography)
- Seaborn library (for enhanced box plot styling)

---

## Basic Usage

### Command Line Execution

```bash
python gcam_plot_box_and_whiskers.py path/to/config.json
```

**Multiple Configuration Files:**
```bash
python gcam_plot_box_and_whiskers.py config1.json config2.json config3.json
```

### What the Script Does

For each configuration block in the JSON file, the script:
1. **Reads** processed GCAM output CSV files
2. **Filters** data by scenarios, categories, regions, and years
3. **Aggregates** data across time to create distributions
4. **Creates** box plots showing statistical summaries
5. **Applies** styling and customization options
6. **Generates** publication-quality PDF (and optionally PNG) figures

---

## Understanding Box Plots

### Box Plot Components

A box-and-whisker plot displays the distribution of data through quartiles:

```
        │                          outlier ○
        │                          
    ────┼────  ← maximum (excluding outliers)
        │    │
        │    │
    ────┼────  ← 75th percentile (Q3)
        │ ▓▓ │
        │ ▓▓ │ ← median (50th percentile / Q2)
        │ ▓▓ │
    ────┼────  ← 25th percentile (Q1)
        │    │
        │    │
    ────┼────  ← minimum (excluding outliers)
        │
        │      outlier ○
```

**Key Statistics Shown:**
- **Box:** Interquartile range (IQR) from Q1 to Q3 (contains middle 50% of data)
- **Line in box:** Median value
- **Whiskers:** Extend to min/max values within 1.5 × IQR
- **Points:** Outliers beyond whiskers

---

## Plot Types

The script supports two main plot types:

### 1. Individual Plots

Each scenario plotted as a separate box. Best for:
- Comparing specific scenarios directly
- Showing distribution of values across time for each scenario
- Analyzing individual simulation runs

**Example:** Comparing "Control" vs "Full feedback" scenarios

### 2. Ensemble Plots

Scenarios are grouped into sets, boxes show combined distribution from ensemble members. Best for:
- Analyzing ensemble spread and uncertainty
- Comparing groups of scenarios
- Understanding variability across ensemble members

**Example:** 5 ensemble members for "Control" vs 5 ensemble members for "Full feedback"

---

## Configuration Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `output_file` | string | Path to the processed CSV file containing GCAM data |
| `scenarios` | list or nested list | Scenario names to plot (list for individual, nested list for ensemble) |

### Commonly Used Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `plot_directory` | string | `"./"` | Directory where plots will be saved |
| `categories` | list | All categories | Specific categories to plot (e.g., crop types, sectors) |
| `regions` | list | `["Global"]` | Regions to include in plots |
| `hue` | string | `None` | Variable to use for color-coding (e.g., "scenario", "region") |
| `x_variable` | string | `category_label` | Variable to plot on x-axis |
| `y_label` | string | Auto-generated | Y-axis label |
| `y_scale` | string | `"linear"` | Y-axis scale ("linear" or "log") |
| `use_latex` | boolean | `false` | Use LaTeX for text rendering |

---

## Complete Parameter Reference Table

### Data Selection Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `output_file` | string | **Yes** | - | Valid file path | Path to processed CSV file |
| `scenarios` | list/nested list | **Yes** | - | List of scenario names or nested lists | Scenarios to plot |
| `scenario_sets` | list | No | `None` | List of strings | Names for ensemble groups (when using nested scenarios) |
| `notify_scenarios_transposed` | boolean | No | `false` | `true`, `false` | Print console message when scenarios are automatically transposed |
| `categories` | list | No | All categories | List of category names | Specific categories to plot |
| `categories_to_exclude` | list | No | `None` | List of category names | Categories to exclude from plotting |
| `category_label` | string | No | `"sector"` | Column name | Column identifying categories |
| `regions` | list | No | `["Global"]` | Region names | Regions to include in plots |
| `region_label` | string | No | `"region"` | Column name | Column identifying regions |
| `basins` | list | No | `None` | Basin names | Water basins to include (alternative to regions) |
| `basin_label` | string | No | `"basin"` | Column name | Column identifying basins |
| `scenario_label` | string | No | `"scenario"` | Column name | Column identifying scenarios |
| `year_label` | string | No | `"year"` | Column name | Column identifying years |
| `value_label` | string | No | `"value"` | Column name | Column containing numeric data |
| `start_year` | integer | No | `2015` | Any year | Start year for data inclusion |
| `end_year` | integer | No | `2100` | Any year | End year for data inclusion |
| `key_columns` | list | No | `None` | List of column names | Columns for grouping (used with landtype groups) |

### Plot Configuration Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `x_variable` | string | No | Same as `category_label` | Column name | Variable to plot on x-axis |
| `hue` | string | No | `None` | Column name or `None` | Variable for color-coding boxes |
| `plot_directory` | string | No | `"./"` | Valid directory path | Directory for output plots |
| `plot_name` | string | No | Auto-generated | Filename with .pdf extension | Name of output plot file |
| `plot_type` | string | No | `"ensemble"` | `"ensemble"`, `"individual"` | Type of plot to create |
| `plot_percent_difference` | boolean | No | `false` | `true`, `false` | Plot percent difference from first scenario |

### Aggregation & Processing Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `mean_or_sum_if_more_than_one_row_in_same_landtype_group` | string | No | `"area_weighted_mean"` | `"area_weighted_mean"`, `"sum"`, `"mean"` | How to aggregate when grouping landtypes |
| `landtype_groups` | string | No | `"modified"` | `"modified"`, `"original"` | Which landtype grouping dictionary to use |
| `multiplier` | number | No | `1` | Any number | Multiply all values by this factor |

### Visual Styling Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `y_label` | string | No | Auto-generated | Any string | Y-axis label |
| `x_label` | string | No | `None` | Any string | X-axis label |
| `y_scale` | string | No | `"linear"` | `"linear"`, `"log"` | Y-axis scale |
| `x_scale` | string | No | `None` | `"linear"`, `"log"`, `None` | X-axis scale |
| `y_limits` | list | No | `None` | `[min, max]` | Y-axis limits |
| `width` | number | No | Default value | Positive number | Figure width in inches |
| `height` | number | No | Default value | Positive number | Figure height in inches |
| `fill_boxes` | boolean | No | `true` | `true`, `false` | Fill boxes with color |
| `linewidth` | number | No | `1` | Positive number | Width of box outlines |
| `marker_size` | number | No | `6` | Positive number | Size of outlier markers |
| `plot_colors` | list | No | Default colors | List of color codes | Colors for boxes |
| `use_latex` | boolean | No | `false` | `true`, `false` | Use LaTeX for text rendering |
| `produce_png` | boolean | No | `false` | `true`, `false` | Also produce PNG version of plot |

### Legend Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `legend_on` | boolean | No | `true` | `true`, `false` | Show legend on plot |
| `legend_place_outside` | boolean | No | `false` | `true`, `false` | Place legend outside plot area |
| `legend_x_offset` | number | No | `None` | Positive number | Horizontal offset for legend placement |
| `legend_num_columns` | integer | No | `1` | Positive integer | Number of columns in legend |
| `legend_label_size` | number | No | Default value | Positive number | Font size for legend labels |

### Font Size Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `x_label_size` | number | No | Default value | Positive number | X-axis label font size |
| `y_label_size` | number | No | Default value | Positive number | Y-axis label font size |
| `x_tick_label_size` | number | No | Default value | Positive number | X-axis tick label font size |
| `y_tick_label_size` | number | No | Default value | Positive number | Y-axis tick label font size |

---

## Detailed Parameter Descriptions

### Core Data Parameters

#### `output_file`
**Type:** String (file path)  
**Required:** Yes  
**Description:** Path to the processed CSV file containing GCAM output data.

**Examples:**
```json
"output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv"
"output_file": "./data/scalars_control+full_feedback.csv"
```

**Requirements:**
- File must exist and be readable
- Must be CSV format with headers
- Should have been processed using `gcam_process_extracted_data.py`

---

#### `scenarios`
**Type:** List of strings OR nested list (list of lists)  
**Required:** Yes  
**Description:** Specifies which scenarios to include in the plot.

**Format 1 - Individual Plots (Simple List):**
```json
"scenarios": ["Control", "Full feedback"]
```
Each scenario creates separate box(es) in the plot.

**Format 2 - Ensemble Plots (Nested List):**

Users can specify ensemble scenarios in **either of two formats**. The script automatically detects and handles both formats using the `transpose_scenarios_if_needed()` function in `utility_functions.py`.

**Format 2A - Organized by Scenario Set (Recommended):**
```json
"scenarios": [
    ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
    ["Full feedback", "Full feedback_2", "Full feedback_3", "Full feedback_4", "Full feedback_5"]
],
"scenario_sets": ["Control", "Full feedback"]
```
Each inner list contains all ensemble members for one scenario set. This format is more intuitive as it groups related scenarios together.

**Format 2B - Organized by Ensemble Member (Alternative):**
```json
"scenarios": [
    ["Control", "Full feedback"],
    ["Control_2", "Full feedback_2"],
    ["Control_3", "Full feedback_3"],
    ["Control_4", "Full feedback_4"],
    ["Control_5", "Full feedback_5"]
],
"scenario_sets": ["Control", "Full feedback"]
```
Each inner list is one ensemble member. Data from all members at the same position are combined.

**Important Notes:**
- For ensemble plots, all inner lists must have the same length
- Scenario names must match exactly with values in the CSV file's scenario column
- Order matters for percent difference calculations (first scenario/group is the reference)
- The `scenario_sets` parameter helps the script detect which format you're using
- Scenario names do not need to follow any specific naming convention

---

#### `scenario_sets`
**Type:** List of strings  
**Required:** No (required for ensemble plots)  
**Default:** `None`  
**Description:** Names for each ensemble group when using nested scenarios.

**Example:**
```json
"scenarios": [
    ["Control", "Control_2"],
    ["Full feedback", "Full feedback_2"]
],
"scenario_sets": ["Control", "Full feedback"]
```

**Usage:**
- Assigns readable names to each position across inner lists
- Used in legend labels and for organizing data

---

#### `categories`
**Type:** List of strings  
**Required:** No  
**Default:** All categories in the data (minus excluded ones)  
**Description:** Specific categories to plot (e.g., crop types, economic sectors, landtypes).

**Examples:**
```json
"categories": ["Rice", "Wheat", "Corn", "Soybean"]
"categories": ["crop", "forest", "shrub", "pasture", "grass"]
"categories": ["agricultural energy use", "construction energy use"]
```

**Special Values:**
- If not specified, uses all unique values from the `category_label` column
- Can specify landtype groups like "crop", "forest" which aggregate multiple landtypes

---

#### `regions`
**Type:** List of strings  
**Required:** No  
**Default:** `["Global"]`  
**Description:** Regions to include in the plot.

**Example:**
```json
"regions": ["Global", "USA", "Brazil", "Canada"]
"regions": ["UsaPacNW", "UsaCstNE", "Caribbean"]
```

**Special Value:**
- `"Global"` - Aggregates data across all regions
- Can be used with basins for basin-level analysis

---

### Plot Configuration Parameters

#### `x_variable`
**Type:** String  
**Required:** No  
**Default:** Same as `category_label`  
**Description:** Determines what variable is plotted on the x-axis.

**Possible Values:**
- `"sector"` or `"landtype"` (category-based, default)
- `"scenario"` (scenario-based)
- `"region"` (region-based)
- `"basin"` (basin-based)
- Any valid column name

**Example:**
```json
"x_variable": "scenario"
```

**Behavior:**
- Default: Categories appear on x-axis
- When set to "scenario": Scenarios appear on x-axis instead
- Works in combination with `hue` parameter

---

#### `hue`
**Type:** String  
**Required:** No  
**Default:** `None`  
**Description:** Variable used to color-code different boxes within the same x-axis position.

**Possible Values:**
- `None` - No color coding, all boxes same color
- `"scenario"` - Color by scenario
- `"region"` - Color by region
- `"sector"` or `"landtype"` - Color by category
- Any valid column name

**Examples:**

**No hue (default):**
```json
"hue": null
```
Result: All boxes same color, categories on x-axis.

**Hue by scenario:**
```json
"hue": "scenario",
"scenarios": ["Control", "Full feedback"]
```
Result: For each category on x-axis, two boxes (one per scenario) with different colors.

**Hue by region:**
```json
"hue": "region",
"regions": ["USA", "Brazil", "Canada"]
```
Result: For each category on x-axis, three boxes (one per region) with different colors.

**Hue by sector (with x_variable as scenario):**
```json
"x_variable": "scenario",
"hue": "sector",
"scenarios": ["Control", "Full feedback"],
"categories": ["Rice", "Wheat", "Corn"]
```
Result: Scenarios on x-axis, each with three colored boxes (one per crop).

---

#### `plot_percent_difference`
**Type:** Boolean  
**Required:** No  
**Default:** `false`  
**Description:** Plot percent difference relative to the first scenario instead of absolute values.

**Example:**
```json
"plot_percent_difference": true
```

**Behavior:**
- First scenario (or first scenario set) becomes the reference (0% difference)
- All other scenarios shown as percent change from reference
- Y-axis label automatically updated to include "% difference"
- Useful for highlighting relative changes

**Formula:**
```
percent_difference = (value - reference_value) / reference_value × 100
```

---

### Aggregation Parameters

#### `mean_or_sum_if_more_than_one_row_in_same_landtype_group`
**Type:** String  
**Required:** No  
**Default:** `"area_weighted_mean"`  
**Possible Values:** `"area_weighted_mean"`, `"sum"`, `"mean"`

**Description:** How to aggregate when multiple landtypes are grouped together (e.g., grouping Forest, ProtectedForest, UnmanagedForest into "forest").

**Example:**
```json
"mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum"
```

**Usage:**
- Only applies when using landtype groups like "crop", "forest", "shrub", etc.
- Choose based on whether the quantity is intensive (use mean/area_weighted_mean) or extensive (use sum)

---

### Visual Styling Parameters

#### `y_scale`
**Type:** String  
**Required:** No  
**Default:** `"linear"`  
**Possible Values:** `"linear"`, `"log"`

**Description:** Scale for y-axis.

**Example:**
```json
"y_scale": "log"
```

**Use Cases:**
- `"log"` - For data spanning multiple orders of magnitude (e.g., commodity prices)
- `"linear"` - For most standard plots

---

#### `y_limits`
**Type:** List of two numbers  
**Required:** No  
**Default:** `None` (auto-scaled)  
**Description:** Manually set y-axis range.

**Example:**
```json
"y_limits": [0, 100]
"y_limits": [-20, 100]
```

**Usage:**
- Useful for consistent scaling across multiple plots
- First value is minimum, second is maximum
- Particularly useful with percent difference plots

---

#### `fill_boxes`
**Type:** Boolean  
**Required:** No  
**Default:** `true`  
**Description:** Whether to fill boxes with color or just show outlines.

**Example:**
```json
"fill_boxes": false
```

---

#### `use_latex`
**Type:** Boolean  
**Required:** No  
**Default:** `false`  
**Description:** Use LaTeX for rendering text in plots.

**Example:**
```json
"use_latex": true
```

**Benefits:**
- Professional typography
- Proper math symbols and formatting (e.g., km$^2$)
- Better for publications

**Requirements:**
- LaTeX must be installed on the system

---

## JSON Configuration Examples

### Example 1: Simple Box Plot - No Hue

**Purpose:** Show distribution of commodity prices across crops for a single scenario

```json
{
    "output_file": "./data/ag_commodity_prices_processed.csv",
    "scenarios": ["Control"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "y_scale": "log",
    "plot_directory": "./plots/box_plots/",
    "plot_name": "box_plot_ag_commodity_no_hue.pdf",
    "use_latex": true
}
```

**Result:**
- X-axis: Categories (BioenergyCrop, Rice, etc.)
- Each category has one box showing distribution across all years (2015-2100)
- All boxes same color
- Logarithmic y-axis

---

### Example 2: Box Plot with Hue by Region

**Purpose:** Compare regional distributions for each crop

```json
{
    "output_file": "./data/ag_commodity_prices_processed.csv",
    "hue": "region",
    "scenarios": ["Control"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "regions": ["Global", "USA", "Brazil", "Canada"],
    "y_scale": "log",
    "plot_directory": "./plots/box_plots/",
    "plot_name": "box_plot_ag_commodity_hue_region.pdf",
    "use_latex": true
}
```

**Result:**
- X-axis: Categories (crops)
- For each crop: 4 colored boxes (Global, USA, Brazil, Canada)
- Colors distinguish regions
- Shows how distributions vary by region for each crop

---

### Example 3: Box Plot with Hue by Scenario

**Purpose:** Compare scenario distributions for each crop

```json
{
    "output_file": "./data/ag_commodity_prices_processed.csv",
    "hue": "scenario",
    "scenarios": ["Control", "Full feedback"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "regions": ["USA", "Brazil", "Canada"],
    "y_scale": "log",
    "plot_directory": "./plots/box_plots/",
    "plot_name": "box_plot_ag_commodity_hue_scenario.pdf",
    "use_latex": true
}
```

**Result:**
- X-axis: Categories (crops)
- For each crop: 2 colored boxes (Control and Full feedback)
- Colors distinguish scenarios
- Data aggregated across specified regions

---

### Example 4: Scenario on X-axis, Category as Hue

**Purpose:** Flip the typical layout - scenarios on x-axis, categories colored

```json
{
    "output_file": "./data/ag_commodity_prices_processed.csv",
    "x_variable": "scenario",
    "hue": "sector",
    "scenarios": ["Control", "Full feedback"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "regions": ["USA", "Brazil", "Canada"],
    "y_scale": "log",
    "plot_directory": "./plots/box_plots/",
    "plot_name": "box_plot_ag_commodity_hue_sector.pdf",
    "use_latex": true
}
```

**Result:**
- X-axis: Scenarios (Control, Full feedback)
- For each scenario: 5 colored boxes (one per crop)
- Colors distinguish crops
- Useful for emphasizing scenario differences

---

### Example 5: Percent Difference Plot

**Purpose:** Show relative changes from baseline scenario

```json
{
    "output_file": "./data/ag_commodity_prices_processed.csv",
    "hue": "scenario",
    "scenarios": ["Control", "Full feedback"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "regions": ["USA", "Brazil", "Canada"],
    "plot_percent_difference": true,
    "plot_directory": "./plots/box_plots/",
    "plot_name": "box_plot_ag_commodity_percent_difference.pdf",
    "use_latex": true
}
```

**Result:**
- X-axis: Categories
- Only Full feedback shown (as percent difference from Control)
- Y-axis shows % difference
- Distribution shows variability in relative change

---

### Example 6: Landtype Groups

**Purpose:** Use GCAM landtype aggregations (crop, forest, etc.)

```json
{
    "output_file": "./data/ag_commodity_prices_processed.csv",
    "hue": "scenario",
    "scenarios": ["Control", "Full feedback"],
    "categories": ["crop", "BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "regions": ["USA", "Brazil", "Canada"],
    "key_columns": ["scenario", "region", "year"],
    "plot_directory": "./plots/box_plots/",
    "plot_name": "box_plot_ag_commodity_landtype_groups.pdf",
    "y_scale": "log",
    "use_latex": true
}
```

**Result:**
- X-axis: Both landtype group ("crop") and individual crops
- "crop" box aggregates all crop types
- Allows comparison of aggregated vs individual categories

---

### Example 7: Vegetation Scalars

**Purpose:** Plot EHC vegetation scalars with specific value column

```json
{
    "output_file": "./data/scalars_control+full_feedback.csv",
    "hue": "scenario",
    "scenarios": ["Control", "Full feedback"],
    "category_label": "landtype",
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "regions": ["USA", "Brazil", "Canada"],
    "value_label": "vegetation",
    "y_label": "Vegetation scalars",
    "plot_directory": "./plots/box_plots/",
    "plot_name": "box_plot_scalars_vegetation.pdf",
    "use_latex": true
}
```

**Result:**
- Uses "vegetation" column instead of default "value"
- Custom y-axis label
- Shows distribution of vegetation scalar values

---

### Example 8: Basin-Level Analysis

**Purpose:** Analyze data at basin (watershed) level instead of regional

```json
{
    "output_file": "./data/scalars_control+full_feedback.csv",
    "hue": "scenario",
    "scenarios": ["Control", "Full feedback"],
    "category_label": "landtype",
    "categories": ["crop", "forest", "shrub", "pasture", "grass"],
    "key_columns": ["scenario", "region", "basin", "year"],
    "regions": ["UsaPacNW", "UsaCstNE", "Caribbean"],
    "region_label": "basin",
    "value_label": "vegetation",
    "y_label": "Vegetation scalars",
    "plot_directory": "./plots/box_plots/",
    "plot_name": "box_plot_scalars_vegetation_basin.pdf",
    "use_latex": true
}
```

**Result:**
- Uses basin-level geographic units instead of regions
- Shows landtype groups
- Data aggregated at finer spatial resolution

---

### Example 9: Ensemble Plot

**Purpose:** Show ensemble spread across multiple simulation members

```json
{
    "output_file": "./data/scalars_control+full_feedback.csv",
    "hue": "scenario",
    "scenarios": [
        ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
        ["Full feedback", "Full feedback_2", "Full feedback_3", "Full feedback_4", "Full feedback_5"]
    ],
    "scenario_sets": ["Control", "Full feedback"],
    "category_label": "landtype",
    "categories": ["crop", "forest", "shrub", "pasture", "grass"],
    "key_columns": ["scenario", "region", "basin", "year"],
    "value_label": "vegetation",
    "y_label": "Vegetation scalars",
    "plot_directory": "./plots/box_plots/ensemble/",
    "plot_name": "box_plot_scalars_vegetation.pdf",
    "use_latex": true
}
```

**Result:**
- Combines data from 5 ensemble members
- Two boxes per category (Control and Full feedback groups)
- Distribution includes variability across ensemble members AND across time

---

### Example 10: Land Allocation with Custom Y-limits

**Purpose:** Plot land areas with consistent y-axis scaling

```json
{
    "output_file": "./data/land_allocation_processed_ensemble.csv",
    "hue": "scenario",
    "scenarios": [
        ["Control", "Control_2", "Control_3"],
        ["Full feedback", "Full feedback_2", "Full feedback_3"]
    ],
    "scenario_sets": ["Control", "Full feedback"],
    "category_label": "landtype",
    "categories": ["crop", "forest", "shrub", "pasture", "grass"],
    "key_columns": ["scenario", "region", "basin", "year"],
    "region_label": "region",
    "y_label": "Global land area (thousands km$^2$)",
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum",
    "y_limits": [0.1, 3000],
    "plot_directory": "./plots/box_plots/ensemble/",
    "plot_name": "box_plot_land_allocation.pdf",
    "use_latex": true
}
```

**Result:**
- Sums land areas within landtype groups
- Fixed y-axis range for consistency
- LaTeX formatting in y-axis label (superscript for km²)

---

## The x_variable and hue System

### Understanding the Layout

Box plots organize data along two dimensions:
1. **X-axis position** - Controlled by `x_variable`
2. **Color grouping** - Controlled by `hue`

### Default Behavior

**Without specifying x_variable or hue:**
```json
{
    "categories": ["Rice", "Wheat", "Corn"],
    "scenarios": ["Control"]
}
```
- X-axis: Categories (Rice, Wheat, Corn)
- All boxes same color
- One box per category

### Common Patterns

#### Pattern 1: Categories on X-axis, Scenarios as Hue
```json
{
    "categories": ["Rice", "Wheat", "Corn"],
    "scenarios": ["Control", "Full feedback"],
    "hue": "scenario"
}
```
**Result:**
```
X-axis:    Rice         Wheat        Corn
Boxes:   [C][FF]      [C][FF]      [C][FF]
Colors:  Blue Red    Blue Red    Blue Red
```

#### Pattern 2: Scenarios on X-axis, Categories as Hue
```json
{
    "x_variable": "scenario",
    "categories": ["Rice", "Wheat", "Corn"],
    "scenarios": ["Control", "Full feedback"],
    "hue": "sector"
}
```
**Result:**
```
X-axis:      Control              Full feedback
Boxes:   [R][W][C]                [R][W][C]
Colors:  Blu Grn Red             Blu Grn Red
```

#### Pattern 3: Categories on X-axis, Regions as Hue
```json
{
    "categories": ["Rice", "Wheat"],
    "regions": ["USA", "China", "India"],
    "hue": "region"
}
```
**Result:**
```
X-axis:        Rice                    Wheat
Boxes:   [USA][CHN][IND]         [USA][CHN][IND]
Colors:   Blu  Grn   Red          Blu  Grn   Red
```

### Valid Combinations

| x_variable | hue | Result |
|------------|-----|--------|
| category (default) | None | One box per category |
| category | scenario | Multiple scenarios per category |
| category | region | Multiple regions per category |
| scenario | None | One box per scenario |
| scenario | category | Multiple categories per scenario |
| scenario | region | Multiple regions per scenario |
| region | None | One box per region |
| region | scenario | Multiple scenarios per region |
| region | category | Multiple categories per region |

---

## Output Files

### Plot Files

**PDF Format (default):**
- High-resolution vector graphics
- Suitable for publications
- Filename: `plot_name` parameter or auto-generated

**PNG Format (optional):**
- Raster graphics
- Enabled with `produce_png: true`
- Same filename as PDF with .png extension

### Auto-generated Filenames

**Default naming pattern:**
```
box_plot_[output_file_name].pdf
```

**Example:**
```
Input: ag_commodity_prices_processed.csv
Output: box_plot_ag_commodity_prices_processed.pdf
```

---

## Advanced Features

### Landtype Grouping

The script supports two levels of landtype aggregation:

**Modified Groups (default):**
- Aggregates GCAM landtypes into major categories
- Example: "crop" includes all crop types
- Set with `"landtype_groups": "modified"`

**Original Groups:**
- Uses original GCAM landtype names with broader groupings
- Set with `"landtype_groups": "original"`

**Landtype Group Examples:**
```python
crop = [BioenergyCrop, Corn, FodderGrass, Rice, Wheat, ...]
forest = [Forest, ProtectedUnmanagedForest, UnmanagedForest]
shrub = [Shrubland, ProtectedUnmanagedShrubland, UnmanagedShrubland]
pasture = [Pasture, ProtectedUnmanagedPasture]
grass = [Grassland, ProtectedUnmanagedGrassland, UnmanagedGrassland]
```

**Usage:**
```json
"categories": ["crop", "forest", "shrub"],
"key_columns": ["scenario", "region", "basin", "year"]
```

---

### Time Aggregation

Unlike time series plots, box plots aggregate all data within the specified year range into a single distribution for each box.

**What this means:**
- Start year: 2015, End year: 2100 → Each box shows distribution of 86 data points (one per year)
- Useful for seeing overall variability across time
- Box components represent temporal variability, not spatial or ensemble variability

**Example:**
For "Rice" in "Control" scenario from 2015-2100:
- Box shows quartiles of 86 annual values
- Median line = median price across all years
- IQR = middle 50% of annual values

---

### Combining Temporal and Ensemble Variability

**Ensemble Plots:**
```json
"scenarios": [
    ["Control", "Control_2", "Control_3"],
    ["Full feedback", "Full feedback_2", "Full feedback_3"]
]
```

**Distribution includes:**
- Temporal variability (2015-2100, 86 years)
- Ensemble variability (3 members)
- Total data points per box: 86 years × 3 members = 258 points

**Useful for:**
- Understanding total uncertainty (time + ensemble)
- Comparing ensemble spread between scenarios

---

### Regional and Basin Aggregation

**Global Aggregation:**
```json
"regions": ["Global"]
```
Creates a copy of entire dataset labeled "Global"

**Multiple Regions:**
```json
"regions": ["USA", "China", "Brazil"]
```
Filters to only specified regions

**Basin-Level:**
```json
"regions": ["UsaPacNW", "UsaCstNE"],
"region_label": "basin"
```
Uses basin column instead of region column

---

### Parallel Processing

The script processes multiple plot configurations in parallel:

**Command:**
```bash
python gcam_plot_box_and_whiskers.py config1.json config2.json config3.json
```

**Behavior:**
- Each JSON configuration file processed independently
- All plots within a file created sequentially
- Utilizes parallel processing with limited cores to reduce memory pressure
- Significantly faster for many plots

---

## Troubleshooting

### Common Issues

#### Issue 1: File Not Found Error
**Error:** `FileNotFoundError: output_file not found`

**Solutions:**
- Verify file path is correct
- Use absolute paths if unsure
- Check file has been processed by `gcam_process_extracted_data.py`

---

#### Issue 2: Empty Boxes or Missing Data
**Symptom:** Boxes don't appear or have no height

**Causes:**
- Year range doesn't overlap with data
- Category names don't match CSV file
- All values identical (no variation)

**Solutions:**
- Check `start_year` and `end_year` overlap with data
- Verify category names match exactly
- Print unique values: `df['category_label'].unique()`

---

#### Issue 3: Hue Not Working
**Symptom:** All boxes same color despite specifying hue

**Causes:**
- Hue variable not in filtered data
- Typo in hue parameter value
- Only one unique value in hue variable

**Solutions:**
- Verify hue column exists in CSV
- Check spelling and case
- Ensure multiple values present after filtering

---

#### Issue 4: X-axis Labels Overlapping
**Symptom:** Category names overlap on x-axis

**Solutions:**
- Reduce number of categories
- Increase figure width: `"width": 12`
- Rotate labels (requires manual matplotlib adjustment)
- Use abbreviations in category names

---

#### Issue 5: LaTeX Errors
**Error:** `RuntimeError: Failed to process string with tex`

**Solutions:**
- Install LaTeX or set `"use_latex": false`
- Escape special characters in labels (e.g., use `\$` for `$`)
- Check y_label for proper LaTeX formatting

---

#### Issue 6: Ensemble Structure Mismatch
**Error:** `ValueError: All inner lists must have same length`

**Cause:** Nested scenario lists have different lengths

**Solution:**
```json
// Incorrect
"scenarios": [
    ["Control", "Control_2"],
    ["Full feedback"]  // Missing Full feedback_2
]

// Correct
"scenarios": [
    ["Control", "Control_2"],
    ["Full feedback", "Full feedback_2"]
]
```

---

## Best Practices

### 1. Organization

**Directory Structure:**
```
project/
├── data/
│   ├── raw/
│   └── processed/
├── plots/
│   ├── box_plots/
│   │   ├── individual/
│   │   ├── ensemble/
│   │   └── percent_difference/
└── configs/
    ├── box_plot_config.json
    └── time_series_config.json
```

---

### 2. Choosing x_variable and hue

**Question: What do I want to compare?**

**Answer: Compare scenarios for each category**
```json
{
    "x_variable": "sector",  // or use default
    "hue": "scenario"
}
```

**Answer: Compare categories for each scenario**
```json
{
    "x_variable": "scenario",
    "hue": "sector"
}
```

**Answer: Compare regions for each category**
```json
{
    "x_variable": "sector",  // or use default
    "hue": "region"
}
```

---

### 3. When to Use Box Plots vs Time Series

**Use Box Plots When:**
- Want to see overall distribution across time
- Interested in quartiles, outliers, and spread
- Comparing statistical summaries
- Space-constrained publications (one plot vs many time points)

**Use Time Series When:**
- Temporal trends are important
- Want to see trajectory over time
- Interested in specific years or periods
- Need to show year-by-year changes

**Use Both:**
- Box plot for overview
- Time series for detailed temporal analysis

---

### 4. Styling for Publications

**High Quality Settings:**
```json
{
    "use_latex": true,
    "width": 7,
    "height": 5,
    "linewidth": 1.5,
    "marker_size": 6,
    "legend_label_size": 10,
    "x_label_size": 12,
    "y_label_size": 12,
    "produce_png": false,
    "fill_boxes": true
}
```

**Presentation Settings:**
```json
{
    "use_latex": false,
    "width": 10,
    "height": 6,
    "linewidth": 2,
    "marker_size": 8,
    "legend_label_size": 14,
    "x_label_size": 16,
    "y_label_size": 16,
    "produce_png": true,
    "fill_boxes": true
}
```

---

### 5. Data Preparation

**Before Plotting:**
- Process with `gcam_process_extracted_data.py`
- Add areas with `gcam_add_areas_to_files.py` (if needed for area-weighted means)
- Verify data structure matches expectations
- Check for missing values or outliers

**Verify Year Range:**
```python
import pandas as pd
df = pd.read_csv('your_file.csv')
print(f"Years: {df['year'].min()} to {df['year'].max()}")
```

---

### 6. Color Schemes

**Default Colors (Matplotlib Tableau palette):**
- Automatically applied when `hue` is specified
- Up to 10 colors, then extends to XKCD colors
- Defined in `utility_plots.py`

**Custom Colors:**
```json
"plot_colors": ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
```

---

### 7. Iterative Workflow

**Recommended Process:**
1. Start with simple configuration (one scenario, few categories)
2. Verify plot looks reasonable
3. Experiment with x_variable and hue combinations
4. Add complexity (more scenarios, categories, regions)
5. Adjust styling parameters
6. Run full ensemble when configuration is finalized

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
   
4. gcam_plot_box_and_whiskers.py  ← THIS SCRIPT
   ↓ (Create box plots)
   
5. Analysis and publication
```

### Complementary Visualizations

**Also create:**
- Time series plots (`gcam_plot_time_series.py`) - Show temporal trends
- Spatial plots (`gcam_plot_spatial_data.py`) - Show geographic patterns

**Together they provide:**
- Box plots: Statistical summary across time
- Time series: Temporal evolution
- Spatial plots: Geographic distribution

---

## Performance Considerations

### Execution Time

**Factors:**
- Number of configurations
- Number of scenarios/ensemble members
- Number of categories, regions, and years
- File size
- LaTeX rendering (slower)

**Typical Times:**
- Simple box plot: 1-2 seconds
- Ensemble plot with 5 members: 2-5 seconds
- Complex multi-category, multi-region: 5-15 seconds

**Optimization:**
- Use PNG instead of PDF for faster rendering
- Disable LaTeX for draft plots
- Reduce year range for testing
- Process multiple configs in parallel

---

## Appendix: Complete Default Values

### All Default Parameters

```python
{
    'basin_label': 'basin',
    'basins': None,
    'category_label': 'sector',
    'end_year': 2100,
    'fill_boxes': True,
    'height': 8,  # Default from utility_plots
    'hue': None,
    'key_columns': None,
    'landtype_groups': 'modified',
    'legend_label_size': 14,  # Default from utility_plots
    'legend_num_columns': 1,
    'legend_on': True,
    'legend_place_outside': False, 
    'legend_x_offset': None,
    'linewidth': 1,
    'marker_size': 6,
    'multiplier': 1,
    'mean_or_sum_if_more_than_one_row_in_same_landtype_group': 'area_weighted_mean',
    'plot_colors': Matplotlib Tableau colors (see utility_plots.py),
    'plot_directory': './',
    'plot_percent_difference': False,
    'plot_type': 'ensemble',
    'produce_png': False, 
    'region_label': 'region', 
    'regions': ['Global'],
    'scenario_label': 'scenario', 
    'scenario_sets': None,
    'start_year': 2015,
    'use_latex': False, 
    'value_label': 'value',
    'width': 10,  # Default from utility_plots
    'x_label': None,
    'x_label_size': 24,  # Default from utility_plots
    'x_scale': None,
    'x_tick_label_size': 20,  # Default from utility_plots
    'y_label_size': 24,  # Default from utility_plots
    'y_limits': None,
    'y_scale': 'linear',
    'y_tick_label_size': 20,  # Default from utility_plots
    'year_label': 'year'
}
```

---

## Quick Reference Guide

### Minimal Configuration

```json
{
    "output_file": "./data/file.csv",
    "scenarios": ["Control", "Full feedback"]
}
```

### Recommended Configuration

```json
{
    "output_file": "./data/file.csv",
    "scenarios": ["Control", "Full feedback"],
    "categories": ["Rice", "Wheat", "Corn"],
    "hue": "scenario",
    "plot_directory": "./plots/",
    "y_scale": "log",
    "use_latex": true
}
```

### Full Featured Configuration

```json
{
    "output_file": "./data/file.csv",
    "scenarios": [
        ["Control", "Full feedback"],
        ["Control_2", "Full feedback_2"],
        ["Control_3", "Full feedback_3"]
    ],
    "scenario_sets": ["Control", "Full feedback"],
    "categories": ["Rice", "Wheat", "Corn"],
    "regions": ["USA", "China", "Brazil"],
    "hue": "scenario",
    "x_variable": "sector",
    "plot_directory": "./plots/ensemble/",
    "plot_name": "box_plot_crops.pdf",
    "y_label": "Commodity price ($/kg)",
    "y_scale": "log",
    "y_limits": [1, 1000],
    "start_year": 2020,
    "end_year": 2100,
    "use_latex": true,
    "width": 10,
    "height": 6,
    "legend_num_columns": 1,
    "legend_place_outside": true,
    "legend_x_offset": 1.3
}
```

---

## References

- DiVittorio et al. (2025). "E3SM-GCAM coupling methodology and applications." *Journal of Advances in Modeling Earth Systems*. [DOI: 10.1029/2024MS004806](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024MS004806)
- GitHub Repository: [https://github.com/philipmyint/e3sm--gcam_analysis](https://github.com/philipmyint/e3sm--gcam_analysis)
- Seaborn Box Plot Documentation: [https://seaborn.pydata.org/generated/seaborn.boxplot.html](https://seaborn.pydata.org/generated/seaborn.boxplot.html)

---

## Contact

For questions or issues:
- **Philip Myint**: myint1@llnl.gov
- **Dalei Hao**: dalei.hao@pnnl.gov

---

## Version Information

**Script:** gcam_plot_box_and_whiskers.py  
**Documentation Version:** 1.0  
**Last Updated:** January 2026  
**Python Version:** 3.7+

---

*This documentation provides comprehensive guidance for using the `gcam_plot_box_and_whiskers.py` script to create statistical visualizations of GCAM output data distributions.*