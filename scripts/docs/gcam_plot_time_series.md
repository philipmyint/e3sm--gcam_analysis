# GCAM Time Series Plotting Script Documentation


**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Overview

**Script Name:** `gcam_plot_time_series.py`

**Purpose:** Creates time series plots from GCAM (Global Change Analysis Model) output files with years on the x-axis and user-specified quantities on the y-axis. The script supports both individual scenario plots and ensemble plots (where scenarios are grouped), performs statistical analyses including t-tests, and provides extensive customization options for publication-quality figures.

**Author:** Philip Myint (myint1@llnl.gov) and Dalei Hao (dalei.hao@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Requirements](#installation--requirements)
3. [Basic Usage](#basic-usage)
4. [Default Values from utility_plots.py](#default-values-from-utility_plotspy)
5. [Plot Types](#plot-types)
6. [Configuration Parameters](#configuration-parameters)
7. [Complete Parameter Reference Table](#complete-parameter-reference-table)
8. [Detailed Parameter Descriptions](#detailed-parameter-descriptions)
9. [JSON Configuration Examples](#json-configuration-examples)
10. [Statistical Analysis](#statistical-analysis)
11. [Plotting Customization Details](#plotting-customization-details)
12. [Output Files](#output-files)
13. [Advanced Features](#advanced-features)
14. [Troubleshooting](#troubleshooting)
15. [Best Practices](#best-practices)

---

## Installation & Requirements

### Required Python Libraries

```bash
pip install pandas numpy scipy matplotlib
```

### Required Utility Modules

The script imports several utility modules that must be in the same directory or Python path:
- `utility_constants` - Default constants for plotting
- `utility_dataframes` - Functions for reading files and performing t-tests
- `utility_functions` - General utility functions
- `utility_gcam` - GCAM-specific functions for landtype grouping
- **`utility_plots`** - **Critical plotting utility functions and default values**

### System Requirements

- Python 3.7+
- Multi-core processor (script uses parallel processing for multiple plots)
- LaTeX installation (optional, for publication-quality typography)

---

## Default Values from utility_plots.py

The `utility_plots.py` module defines all default plotting parameters. Understanding these defaults helps you know what happens when you don't specify certain parameters.

### Figure Dimensions (from utility_plots.py)

| Parameter | Default Value | Unit | Description |
|-----------|---------------|------|-------------|
| `width_default` | `10` | inches | Default figure width |
| `height_default` | `8` | inches | Default figure height |

**Usage:** Figures are 10" × 8" by default, suitable for presentations and papers.

### Scale and Limits (from utility_plots.py)

| Parameter | Default Value | Options | Description |
|-----------|---------------|---------|-------------|
| `scale_default` | `'linear'` | `'linear'`, `'log'` | Default for both x and y axes |
| `axis_limits_default` | `None` | `None` or `[min, max]` | Auto-scaled if not specified |

### Font Sizes (from utility_plots.py)

| Parameter | Default Value | Unit | Applied To |
|-----------|---------------|------|------------|
| `axis_label_size_default` | `24` | points | X and Y axis labels |
| `tick_label_size_default` | `20` | points | X and Y tick labels (numbers on axes) |
| `legend_label_size_default` | `14` | points | Legend text |

**Note:** These are publication-quality sizes. For presentations, consider increasing these values.

### Legend Configuration (from utility_plots.py)

| Parameter | Default Value | Description |
|-----------|---------------|-------------|
| `legend_num_columns_default` | `1` | Number of columns in legend |
| `legend_on_default` | `True` | Whether legend is displayed |
| `legend_place_outside_default` | `False` | Legend inside plot area by default |

**Auto-calculated Legend Offset:**
- If `legend_place_outside = True` and `legend_num_columns = 1`: `x_offset = 1.3`
- If `legend_place_outside = True` and `legend_num_columns > 1`: `x_offset = 1 + 0.45 * num_columns`

### Line and Marker Styles (from utility_plots.py)

| Parameter | Default Value | Description |
|-----------|---------------|-------------|
| `linewidth_default` | `2` | Width of plot lines in points |

**Default Colors (plot_colors_default):**
```python
# Matplotlib Tableau palette (10 colors)
['#1f77b4',  # Blue
 '#ff7f0e',  # Orange
 '#2ca02c',  # Green
 '#d62728',  # Red
 '#9467bd',  # Purple
 '#8c564b',  # Brown
 '#e377c2',  # Pink
 '#7f7f7f',  # Gray
 '#bcbd22',  # Olive
 '#17becf']  # Cyan

# Extended with XKCD colors (11 more)
['#029386',  # Teal
 '#c20078',  # Magenta
 '#53fca1',  # Sea green
 '#fe01b1',  # Bright pink
 '#c65102',  # Dark orange
 '#fac205',  # Goldenrod
 '#0b5509',  # Forest
 '#8a6e45',  # Dirt
 '#fc5a50',  # Coral
 '#a2cffe',  # Baby blue
 '#ffb07c']  # Peach
```

**Default Line Styles (linestyle_tuples_default):**
```python
[
    ('solid',                 (0, ())),
    ('dashed',                (0, (5, 5))),
    ('dotted',                (0, (1, 5))),
    ('dashdot',               (0, (3, 5, 1, 5))),
    ('loosely dotted',        (0, (1, 10))),
    ('loosely dashed',        (0, (5, 10))),
    ('densely dashed',        (0, (5, 1))),
    ('loosely dashdotted',    (0, (3, 10, 1, 10))),
    ('densely dashdotted',    (0, (3, 1, 1, 1))),
    ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
    ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
    ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))
]
```

**Default Markers (markers_default):**
```python
['o',   # Circle
 'v',   # Triangle down
 '^',   # Triangle up
 '<',   # Triangle left
 '>',   # Triangle right
 '8',   # Octagon
 's',   # Square
 'p',   # Pentagon
 '*',   # Star
 'h',   # Hexagon 1
 'H',   # Hexagon 2
 'D',   # Diamond
 'd',   # Thin diamond
 'P',   # Plus (filled)
 'X',   # X (filled)
 '1', '2', '3', '4',  # Tri down, up, left, right
 '+',   # Plus
 'x',   # X
 '|']   # Vertical line
```

### Output Options (from utility_plots.py)

| Parameter | Default Value | Description |
|-----------|---------------|-------------|
| `produce_png_default` | `False` | Only PDF generated by default |
| `use_latex_default` | `False` | Plain text rendering by default |
| `bbox_inches_default` | `'tight'` | Tight bounding box (crops whitespace) |

---

## Basic Usage

### Command Line Execution

```bash
python gcam_plot_time_series.py path/to/config.json
```

**Multiple Configuration Files:**
```bash
python gcam_plot_time_series.py config1.json config2.json config3.json
```

### What the Script Does

For each configuration block in the JSON file, the script:
1. **Reads** processed GCAM output CSV files
2. **Filters** data by scenarios, categories, regions, and years
3. **Aggregates** data (mean, area-weighted mean, or sum) for each year
4. **Creates** time series plots with customizable styling
5. **Performs** statistical tests (t-tests) comparing scenarios
6. **Generates** publication-quality PDF (and optionally PNG) figures
7. **Outputs** p-value statistics to text files

---

## Plot Types

The script supports two main plot types:

### 1. Individual Plots

Each scenario is plotted as a separate line. Best for:
- Comparing a small number of scenarios directly
- Showing all individual simulation runs
- Regional comparisons within a single scenario

**Example:** Comparing "Control" vs "Full feedback" scenarios

**Visual Example:**
```
Value
  ^
  |     ___
  |    /   \___    Legend:
  |   /        \__ — Control
  | _/            — Full feedback
  |
  +-------------------> Year
  2015            2100
```

### 2. Ensemble Plots

Scenarios are grouped into sets, and the mean of each set is plotted with error bars representing ±N standard deviations. Best for:
- Analyzing uncertainty across multiple ensemble members
- Comparing groups of scenarios with different conditions
- Showing ensemble spread and statistical significance

**Example:** 5 ensemble members for "Control" vs 5 ensemble members for "Full feedback"

**Visual Example:**
```
Value
  ^
  |     ___
  |    /█░█\___    Legend:
  |   /█░░░█   \__ — Control (mean ± 1σ)
  | _/░░░░░░█  ░__ — Full feedback (mean ± 1σ)
  |
  +-------------------> Year
  2015            2100
  
  █ = Mean line
  ░ = ± 1 standard deviation envelope
```

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
| `regions` | list or dict | `["Global"]` | Regions to include in plots |
| `y_label` | string | Auto-generated | Y-axis label |
| `y_scale` | string | `"linear"` | Y-axis scale ("linear" or "log") |
| `width` | number | `10` inches | Figure width (from utility_plots.py) |
| `height` | number | `8` inches | Figure height (from utility_plots.py) |
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
| `regions` | list or dict | No | `["Global"]` | Region names or dict mapping | Regions to plot for each category |
| `region_label` | string | No | `"region"` | Column name | Column identifying regions |
| `scenario_label` | string | No | `"scenario"` | Column name | Column identifying scenarios |
| `year_label` | string | No | `"year"` | Column name | Column identifying years |
| `value_label` | string | No | `"value"` | Column name | Column containing numeric data |
| `start_year` | integer | No | `2015` | Any year | Start year for time series |
| `end_year` | integer | No | `2100` | Any year | End year for time series |
| `key_columns` | list | No | `None` | List of column names | Columns for grouping (used with landtype groups) |

### Aggregation & Processing Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `aggregation_type_in_each_year` | string | No | `"area_weighted_mean"` | `"mean"`, `"area_weighted_mean"`, `"sum"` | How to aggregate data within each year |
| `mean_or_sum_if_more_than_one_row_in_same_landtype_group` | string | No | `"area_weighted_mean"` | `"area_weighted_mean"`, `"sum"`, `"mean"` | How to aggregate when grouping landtypes |
| `landtype_groups` | string | No | `"modified"` | `"modified"`, `"original"` | Which landtype grouping dictionary to use |
| `multiplier` | number | No | `1` | Any number | Multiply all values by this factor |
| `set_nan_to_zero` | boolean | No | `false` | `true`, `false` | Replace NaN values with zero |

### Figure Dimension Parameters (from utility_plots.py)

| Parameter | Type | Required | Default (from utility_plots.py) | Possible Values | Description |
|-----------|------|----------|--------------------------------|-----------------|-------------|
| `width` | number | No | `10` inches | Positive number | Figure width in inches |
| `height` | number | No | `8` inches | Positive number | Figure height in inches |

### Plot Styling Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `plot_directory` | string | No | `"./"` | Valid directory path | Directory for output plots |
| `plot_name` | string | No | Auto-generated | Filename with .pdf extension | Name of output plot file |
| `plot_type` | string | No | `"ensemble"` | `"ensemble"`, `"individual"` | Type of plot to create |
| `plot_percent_difference` | boolean | No | `false` | `true`, `false` | Plot percent difference from first scenario |
| `y_label` | string | No | Auto-generated | Any string | Y-axis label |
| `y_scale` | string | No | `"linear"` (from utility_plots.py) | `"linear"`, `"log"` | Y-axis scale |
| `y_limits` | list | No | `None` (from utility_plots.py) | `[min, max]` | Y-axis limits |
| `x_scale` | string | No | `"linear"` (from utility_plots.py) | `"linear"`, `"log"` | X-axis scale |
| `x_limits` | list | No | `None` (from utility_plots.py) | `[min, max]` | X-axis limits |
| `use_latex` | boolean | No | `false` (from utility_plots.py) | `true`, `false` | Use LaTeX for text rendering |
| `produce_png` | boolean | No | `false` (from utility_plots.py) | `true`, `false` | Also produce PNG version of plot |

### Legend Parameters

| Parameter | Type | Required | Default (from utility_plots.py) | Possible Values | Description |
|-----------|------|----------|--------------------------------|-----------------|-------------|
| `legend_on` | boolean | No | `true` | `true`, `false` | Show legend on plot |
| `legend_place_outside` | boolean | No | `true` | `true`, `false` | Place legend outside plot area |
| `legend_x_offset` | number | No | Auto-calculated (see below) | Positive number | Horizontal offset for legend placement |
| `legend_num_columns` | integer | No | Auto-calculated | Positive integer | Number of columns in legend |
| `legend_label_size` | number | No | `14` points | Positive number | Font size for legend labels |

**Auto-calculated legend_x_offset (from utility_plots.py):**
- If `legend_num_columns = 1`: Default offset = `1.3`
- If `legend_num_columns > 1`: Default offset = `1 + 0.45 * legend_num_columns`

### Line & Marker Styling Parameters

| Parameter | Type | Required | Default (from utility_plots.py) | Possible Values | Description |
|-----------|------|----------|--------------------------------|-----------------|-------------|
| `plot_colors` | list | No | 21-color palette (see above) | List of color codes | Colors for lines |
| `markers` | list | No | 22 marker styles (see above) | List of marker styles | Marker shapes for data points |
| `linestyle_tuples` | list | No | 12 line styles (see above) | List of linestyle tuples | Line styles (solid, dashed, etc.) |
| `linewidth` | number | No | `2` points | Positive number | Width of plot lines |
| `marker_size` | number | No | `6` | Positive number | Size of data point markers |

### Font Size Parameters (from utility_plots.py)

| Parameter | Type | Required | Default (from utility_plots.py) | Possible Values | Description |
|-----------|------|----------|--------------------------------|-----------------|-------------|
| `x_label_size` | number | No | `24` points | Positive number | X-axis label font size |
| `y_label_size` | number | No | `24` points | Positive number | Y-axis label font size |
| `x_tick_label_size` | number | No | `20` points | Positive number | X-axis tick label font size |
| `y_tick_label_size` | number | No | `20` points | Positive number | Y-axis tick label font size |

**Font Size Recommendations:**
- **Publications:** Use defaults (24pt labels, 20pt ticks, 14pt legend)
- **Presentations:** Increase by 20-30% (30-32pt labels, 24-26pt ticks, 16-18pt legend)
- **Posters:** Increase by 50-70% (36-40pt labels, 30-34pt ticks, 20-24pt legend)

### Statistical Analysis Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `p_value_threshold` | number | No | `0.05` | 0 to 1 | Significance threshold for t-tests |
| `p_value_file` | string | No | `"p_values.dat"` | Filename | File to store p-value results |
| `p_value_file_print_only_if_below_threshold` | boolean | No | `true` | `true`, `false` | Only print significant p-values |
| `p_value_marker_size` | number | No | `10` | Positive number | Size of markers showing significance |
| `std_multiplier` | number | No | `1` | Positive number | Multiplier for error bars (±N std dev) |
| `std_mean_across_all_data_multiplier` | number | No | `1` | Positive number | Multiplier for overall mean error bars |
| `include_mean_across_all_data` | boolean | No | `false` | `true`, `false` | Include overall mean across all time series |

### Error Bar Parameters

| Parameter | Type | Required | Default | Possible Values | Description |
|-----------|------|----------|---------|-----------------|-------------|
| `error_bars_alpha` | number | No | `0.2` | 0 to 1 | Transparency of error bar shading |

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
"output_file": "./data/co2_emissions_regions_processed.csv"
```

**Requirements:**
- File must exist and be readable
- Must be CSV format with headers
- Should have been processed using `gcam_process_extracted_data.py`

---

#### `scenarios`
**Type:** List of strings OR nested list (list of lists)  
**Required:** Yes  
**Description:** Specifies which scenarios to plot.

**Format 1 - Individual Plots (Simple List):**
```json
"scenarios": ["Control", "Full feedback"]
```
Each scenario plotted as a separate line.

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
Each inner list is one ensemble member. Scenarios at the same position across inner lists are grouped together.

**Important Notes:**
- For ensemble plots, all inner lists must have the same length
- Scenario names must match exactly with values in the CSV file's scenario column
- Order matters for statistical comparisons (first scenario/group is the reference)
- The `scenario_sets` parameter helps the script detect which format you're using
- Scenario names do not need to follow any specific naming convention (e.g., they don't need `_2`, `_3` suffixes)

---

#### `scenario_sets`
**Type:** List of strings  
**Required:** No (required for ensemble plots)  
**Default:** `None`  
**Description:** Names for each ensemble group when using nested scenarios.

**Example:**
```json
"scenarios": [
    ["Control", "Full feedback"],
    ["Control_2", "Full feedback_2"]
],
"scenario_sets": ["Control", "Full feedback"]
```

**Usage:**
- First inner list of scenarios → labeled "Control"
- Second inner list → labeled "Full feedback"
- Used in legend labels and statistical output

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
- Can specify landtype groups like "crop", "forest" which will aggregate multiple landtypes

---

#### `categories_to_exclude`
**Type:** List of strings  
**Required:** No  
**Default:** `None`  
**Description:** Categories to exclude from automatic category selection.

**Example:**
```json
"categories_to_exclude": ["UnmanagedLand", "Forest"]
```

**Usage:**
- Only applies when `categories` is not explicitly specified
- Useful for excluding irrelevant or problematic categories

---

#### `regions`
**Type:** List of strings OR dictionary  
**Required:** No  
**Default:** `["Global"]`  
**Description:** Regions to include in the plot.

**Format 1 - Simple List (same regions for all categories):**
```json
"regions": ["Global", "USA", "China", "EU-15"]
```

**Format 2 - Dictionary (different regions per category):**
```json
"regions": {
    "Rice": ["Global", "China", "India"],
    "Wheat": ["Global", "USA", "Russia"],
    "Corn": ["Global", "USA", "Brazil"]
}
```

**Special Value:**
- `"Global"` - Aggregates data across all regions

---

### Aggregation Parameters

#### `aggregation_type_in_each_year`
**Type:** String  
**Required:** No  
**Default:** `"area_weighted_mean"`  
**Possible Values:** `"mean"`, `"area_weighted_mean"`, `"sum"`

**Description:** How to aggregate data across multiple rows within each year.

**Use Cases:**
- `"mean"` - Simple average (e.g., average price across regions)
- `"area_weighted_mean"` - Weighted by area (e.g., average yield weighted by cropland area)
- `"sum"` - Total (e.g., total CO2 emissions, total land area)

**Example:**
```json
"aggregation_type_in_each_year": "sum"
```

**Important:** 
- Area-weighted mean requires an `area` column in the CSV file
- Use `sum` for extensive quantities (totals)
- Use `mean` or `area_weighted_mean` for intensive quantities (rates, prices, scalars)

---

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

### Plot Styling Parameters

#### `plot_directory`
**Type:** String (directory path)  
**Required:** No  
**Default:** `"./"`  
**Description:** Directory where output plots will be saved.

**Examples:**
```json
"plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/individual"
"plot_directory": "./output/ensemble"
```

**Behavior:**
- Directory is created automatically if it doesn't exist
- Can use relative or absolute paths

---

#### `plot_name`
**Type:** String (filename)  
**Required:** No  
**Default:** Auto-generated from output_file name  
**Description:** Name of the output plot file.

**Examples:**
```json
"plot_name": "time_series_vegetation_scalars.pdf"
"plot_name": "emissions_comparison.pdf"
```

**Behavior:**
- If only filename provided (no path), placed in `plot_directory`
- If full path provided, saved to that location
- Typically `.pdf` extension (PNG also generated if `produce_png` is true)

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

#### `y_scale` and `x_scale`
**Type:** String  
**Required:** No  
**Default:** `"linear"` (from utility_plots.py)  
**Possible Values:** `"linear"`, `"log"`

**Description:** Scale for y-axis or x-axis.

**Example:**
```json
"y_scale": "log"
```

**Use Cases:**
- `"log"` - For data spanning multiple orders of magnitude
- `"linear"` - For most standard plots

---

#### `width` and `height`
**Type:** Number  
**Required:** No  
**Default:** `width = 10` inches, `height = 8` inches (from utility_plots.py)

**Description:** Figure dimensions in inches.

**Examples:**
```json
"width": 12,
"height": 8
```

**Aspect Ratio Guide:**
- Standard (4:3): `width: 12, height: 9`
- Widescreen (16:9): `width: 16, height: 9`
- Publication (varies): `width: 7, height: 5`
- Square: `width: 10, height: 10`

**Resolution Note:**
- PDF output is vector (infinite resolution)
- PNG output (if `produce_png: true`): Resolution determined by DPI (typically 300 DPI for publications)

---

#### `linewidth`
**Type:** Number  
**Required:** No  
**Default:** `2` points (from utility_plots.py)

**Description:** Width of plot lines in points.

**Recommendations:**
- Thin lines (1-1.5): Dense plots with many lines
- Medium lines (2-2.5): Standard use (default)
- Thick lines (3-4): Presentations and posters
- Very thick (5+): Only for extreme cases (large posters)

**Example:**
```json
"linewidth": 2.5
```

---

#### `use_latex`
**Type:** Boolean  
**Required:** No  
**Default:** `false` (from utility_plots.py)

**Description:** Use LaTeX for rendering text in plots.

**When True:**
- Enables professional typography
- Math symbols render correctly (e.g., km$^2$ → km²)
- Requires LaTeX installation on system
- May slow down plot generation

**When False:**
- Uses matplotlib's default text rendering
- Faster generation
- No LaTeX installation required
- Math symbols may not display properly

**Example with LaTeX formatting:**
```json
{
    "y_label": "Global land area (thousands km$^2$)",
    "use_latex": true
}
```
Result: "Global land area (thousands km²)"

---

#### `produce_png`
**Type:** Boolean  
**Required:** No  
**Default:** `false` (from utility_plots.py)

**Description:** Additionally produce PNG version of plot.

**When False (default):**
- Only PDF generated
- Vector graphics (scalable, small file size)
- Best for publications

**When True:**
- Both PDF and PNG generated
- PNG is raster graphics (fixed resolution)
- Useful for presentations, websites, quick previews

**Example:**
```json
"produce_png": true
```

**Output:**
- `plot_name.pdf` (always)
- `plot_name.png` (only if `produce_png: true`)

---

### Legend Configuration

#### `legend_place_outside`
**Type:** Boolean  
**Required:** No  
**Default:** `false` (from utility_plots.py), but script changes to `true`

**Description:** Whether to place legend outside the plot area.

**When False:**
- Legend placed inside plot area at best location
- Matplotlib chooses location to minimize data overlap
- More compact figure

**When True:**
- Legend placed to the right of plot
- No data obstruction
- Figure wider (consider adjusting `width`)

**Example:**
```json
"legend_place_outside": true,
"legend_x_offset": 1.4
```

---

#### `legend_x_offset`
**Type:** Number  
**Required:** No  
**Default:** Auto-calculated (from utility_plots.py)
- If `legend_num_columns = 1`: `1.3`
- If `legend_num_columns > 1`: `1 + 0.45 * legend_num_columns`

**Description:** Horizontal position of legend when placed outside.

**How It Works:**
- Value of `1.0` = right edge of plot area
- Value of `1.3` = 30% of plot width to the right
- Larger values move legend further right

**Examples:**
```json
"legend_x_offset": 1.3   // Close to plot
"legend_x_offset": 1.5   // Further from plot
"legend_x_offset": 1.1   // Very close to plot
```

**Adjustment Guide:**
- Long legend labels: Increase offset (1.4-1.6)
- Short legend labels: Use default (1.3)
- Multiple columns: Auto-calculation usually sufficient

---

### Font Size Configuration

All font sizes specified in **points** (from utility_plots.py).

#### Default Font Sizes

| Element | Default Size | Typical Range |
|---------|-------------|---------------|
| Axis labels (`x_label_size`, `y_label_size`) | `24` points | 12-40 points |
| Tick labels (`x_tick_label_size`, `y_tick_label_size`) | `20` points | 10-34 points |
| Legend labels (`legend_label_size`) | `14` points | 9-24 points |

#### Font Size Scaling Guidelines

**Rule of Thumb:** All text should be readable at the figure's final display size.

**For Different Media:**

**Academic Papers:**
- Figures typically 3.5" (single column) or 7" (double column) wide
- Use smaller fonts: labels 12-14pt, ticks 10-12pt, legend 8-10pt

**Presentations:**
- Viewed from distance
- Use larger fonts: labels 28-36pt, ticks 24-30pt, legend 18-24pt

**Posters:**
- Viewed from several feet away
- Use largest fonts: labels 36-48pt, ticks 30-40pt, legend 22-28pt

**Example Configurations:**

**Publication Quality:**
```json
{
    "x_label_size": 12,
    "y_label_size": 12,
    "x_tick_label_size": 10,
    "y_tick_label_size": 10,
    "legend_label_size": 9
}
```

**Presentation:**
```json
{
    "x_label_size": 32,
    "y_label_size": 32,
    "x_tick_label_size": 26,
    "y_tick_label_size": 26,
    "legend_label_size": 20
}
```

---

## JSON Configuration Examples

### Example 1: Simple Individual Plot

**Purpose:** Compare two scenarios for agricultural commodity prices

```json
{
    "output_file": "./data/ag_commodity_prices_processed.csv",
    "scenarios": ["Control", "Full feedback"],
    "categories": ["Rice", "Wheat", "Corn"],
    "y_scale": "log",
    "plot_directory": "./plots/",
    "use_latex": true
}
```

**What it does:**
- Plots Control and Full feedback as separate lines
- Shows Rice, Wheat, and Corn on the same plot
- Uses logarithmic y-axis scale
- Saves to `./plots/` directory
- Uses LaTeX for text rendering

---

### Example 2: Excluding Specific Categories

**Purpose:** Plot all categories except unwanted ones

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
    "scenarios": ["Control", "Full feedback"],
    "categories_to_exclude": ["UnmanagedLand", "Forest"],
    "y_scale": "log",
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/individual",
    "use_latex": true
}
```

**What it does:**
- Automatically includes all categories in the file
- Excludes "UnmanagedLand" and "Forest"
- Useful when you want most categories but need to remove a few

---

### Example 3: Percent Difference Plot

**Purpose:** Show relative changes from baseline scenario

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
    "scenarios": ["Control", "Full feedback"],
    "categories_to_exclude": ["UnmanagedLand", "Forest"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/individual_percent_difference",
    "plot_percent_difference": true,
    "use_latex": true
}
```

**What it does:**
- Calculates percent difference of Full feedback relative to Control
- Only shows Full feedback line (Control is the 0% reference)
- Y-axis automatically labeled with "% difference"
- Useful for highlighting relative changes

---

### Example 4: Regional Comparison

**Purpose:** Compare specific regions for selected categories

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
    "scenarios": ["Control"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "regions": ["Global", "USA", "Brazil", "Canada"],
    "legend_x_offset": 1.4,
    "legend_num_columns": 1,
    "y_scale": "log",
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/individual_regions",
    "use_latex": true
}
```

**What it does:**
- Single scenario with multiple regions
- Different line styles for each region
- Each category gets its own marker shape
- Legend positioned outside plot area

---

### Example 5: Ensemble Plot with Error Bars

**Purpose:** Compare ensemble sets with uncertainty quantification

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed_ensemble.csv",
    "scenarios": [
        ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
        ["Full feedback", "Full feedback_2", "Full feedback_3", "Full feedback_4", "Full feedback_5"]
    ],
    "scenario_sets": ["Control", "Full feedback"],
    "categories": ["Rice", "Soybean", "Wheat"],
    "y_scale": "log",
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/ensemble",
    "legend_num_columns": 1,
    "legend_x_offset": 1.4,
    "use_latex": true
}
```

**What it does:**
- Groups 5 ensemble members into 2 scenario sets
- Plots mean of each set with ±1 std deviation error bars
- Shows uncertainty across ensemble members
- Custom legend positioning

---

### Example 6: Ensemble Percent Difference

**Purpose:** Show ensemble uncertainty in relative changes

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed_ensemble.csv",
    "scenarios": [
        ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
        ["Full feedback", "Full feedback_2", "Full feedback_3", "Full feedback_4", "Full feedback_5"]
    ],
    "scenario_sets": ["Control", "Full feedback"],
    "categories": ["Rice", "Soybean", "Wheat"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/ensemble_percent_difference",
    "legend_num_columns": 1,
    "legend_x_offset": 1.2,
    "plot_percent_difference": true,
    "use_latex": true
}
```

**What it does:**
- Combines ensemble and percent difference features
- Shows Full feedback ensemble as % change from Control ensemble
- Error bars represent uncertainty in the percent difference

---

### Example 7: CO2 Emissions (Sum Aggregation)

**Purpose:** Plot total CO2 emissions by region

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/co2_emissions_regions_processed.csv",
    "scenarios": ["Control", "Full feedback"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/individual",
    "aggregation_type_in_each_year": "sum",
    "legend_num_columns": 1,
    "use_latex": true
}
```

**What it does:**
- Uses "sum" aggregation (appropriate for emissions totals)
- Plots all regions in the file
- Compares two scenarios

---

### Example 8: Sectoral Emissions with Regional Detail

**Purpose:** Compare specific sectors across multiple regions

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/co2_emissions_sectors_processed.csv",
    "scenarios": ["Control"],
    "categories": ["agricultural energy use", "construction energy use", "trn_aviation_intl", "refining"],
    "regions": ["Global", "EU-15", "Japan", "India"],
    "legend_x_offset": 1.45,
    "legend_num_columns": 1,
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/individual_regions",
    "aggregation_type_in_each_year": "sum",
    "use_latex": true
}
```

**What it does:**
- Specific sectors and regions
- Sum aggregation for emissions
- Different combinations of line colors/styles/markers for visual distinction

---

### Example 9: Vegetation Scalars

**Purpose:** Plot EHC vegetation scalars for landtype groups

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback.csv",
    "scenarios": ["Full feedback"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/individual_full_feedback_only",
    "plot_name": "time_series_vegetation_scalars.pdf",
    "categories": ["crop", "forest", "shrub", "pasture", "grass"],
    "category_label": "landtype",
    "key_columns": ["scenario", "region", "basin", "year"],
    "value_label": "vegetation",
    "y_label": "vegetation scalars",
    "legend_num_columns": 1,
    "legend_x_offset": 1.3,
    "use_latex": true
}
```

**What it does:**
- Plots vegetation scalar values (not the default "value" column)
- Uses landtype groups (crop, forest, etc.)
- Custom y-axis label
- Specific plot filename

---

### Example 10: Soil Scalars

**Purpose:** Plot EHC soil scalars for landtype groups

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback.csv",
    "scenarios": ["Full feedback"],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/individual_full_feedback_only",
    "plot_name": "time_series_soil_scalars.pdf",
    "categories": ["crop", "forest", "shrub", "pasture", "grass"],
    "category_label": "landtype",
    "key_columns": ["scenario", "region", "basin", "year"],
    "value_label": "soil",
    "y_label": "soil scalars",
    "legend_num_columns": 1,
    "legend_x_offset": 1.3,
    "use_latex": true
}
```

**What it does:**
- Similar to Example 9 but for soil scalars
- Different value column ("soil" instead of "vegetation")
- Different y-axis label

---

### Example 11: Single Ensemble Set (Control Only)

**Purpose:** Show uncertainty within Control scenario ensemble

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback_ensemble.csv",
    "scenarios": [
        ["Control", "Control_2", "Control_3", "Control_4", "Control_5"]
    ],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/ensemble_control_only",
    "plot_name": "time_series_vegetation_scalars.pdf",
    "scenario_sets": ["Control"],
    "categories": ["crop", "forest", "shrub", "pasture", "grass"],
    "category_label": "landtype",
    "key_columns": ["scenario", "region", "basin", "year"],
    "value_label": "vegetation",
    "y_label": "vegetation scalars",
    "legend_num_columns": 1,
    "legend_x_offset": 1.2,
    "use_latex": true
}
```

**What it does:**
- Ensemble plot with single scenario set
- Shows internal ensemble variability
- Useful for understanding baseline uncertainty

---

### Example 12: Land Allocation with Area Aggregation

**Purpose:** Show global land area changes over time with original crop names

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_original_crop_names_ensemble.csv",
    "scenarios": [
        ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
        ["Full feedback", "Full feedback_2", "Full feedback_3", "Full feedback_4", "Full feedback_5"]
    ],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/ensemble",
    "scenario_sets": ["Control", "Full feedback"],
    "categories": ["crop", "forest", "shrub", "pasture", "grass"],
    "category_label": "landtype",
    "landtype_groups": "original",
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum",
    "aggregation_type_in_each_year": "sum",
    "y_label": "Global land area (thousands km$^2$)",
    "key_columns": ["scenario", "region", "basin", "year"],
    "legend_num_columns": 1,
    "legend_x_offset": 1.35,
    "use_latex": true
}
```

**What it does:**
- Ensemble plot of land allocation
- Uses "original" landtype groupings
- Sums areas within landtype groups
- Sums across all regions for global total
- LaTeX formatting in y-axis label ($^2$ for superscript)

---

### Example 13: Land Allocation with Modified Crop Names

**Purpose:** Same as Example 12 but with modified (standardized) crop names

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_ensemble.csv",
    "scenarios": [
        ["Control", "Full feedback"],
        ["Control_2", "Full feedback_2"],
        ["Control_3", "Full feedback_3"],
        ["Control_4", "Full feedback_4"],
        ["Control_5", "Full feedback_5"]
    ],
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots/ensemble",
    "scenario_sets": ["Control", "Full feedback"],
    "categories": ["crop", "forest", "shrub", "pasture", "grass"],
    "category_label": "landtype",
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum",
    "aggregation_type_in_each_year": "sum",
    "y_label": "Global land area (thousands km$^2$)",
    "key_columns": ["scenario", "region", "basin", "year"],
    "legend_num_columns": 1,
    "legend_x_offset": 1.35,
    "use_latex": true
}
```

**What it does:**
- Uses "modified" landtype groups (default)
- Standardized crop names for consistency
- Otherwise identical to Example 12

---

## Statistical Analysis

### T-Tests Performed

The script automatically performs t-tests to compare:

#### 1. Scenario Comparisons (Individual Plots)
- Compares first scenario vs. each other scenario
- Performed for each category and region combination
- Tests entire time series (all years combined)

**Example:**
```
Comparing: Control (reference) vs. Full feedback
For: Rice, Global region
Data: All years from start_year to end_year
Test: Independent samples t-test
```

#### 2. Ensemble Comparisons (Ensemble Plots)
- Compares first ensemble set vs. each other set
- Year-by-year t-tests
- Overall time series comparison
- Statistical significance markers on plots

**Example:**
```
Comparing: Control ensemble (reference) vs. Full feedback ensemble
For each year: t-test of 5 Control values vs. 5 Full feedback values
Overall: t-test of all Control data vs. all Full feedback data
```

#### 3. Regional Comparisons
- When single scenario with multiple regions
- Compares first region vs. each other region

**Example:**
```
Comparing: Global (reference) vs. USA, Brazil, Canada
For: Each category separately
Test: Independent samples t-test
```

### Interpretation of Results

**P-Value Output File:**
```
============================================================
Plot: ./plots/time_series_emissions.pdf
------------------------------------------------------------
scenario=Full feedback, category=Rice, region=Global: p=0.0234 *
scenario=Full feedback, category=Wheat, region=Global: p=0.3421
scenario=Full feedback, category=Corn, region=Global: p=0.0012 **
============================================================
```

**Significance Markers:**
- `*` indicates p-value below threshold (default 0.05)
- `**` could indicate p < 0.01 (very significant)
- Only significant results printed if `p_value_file_print_only_if_below_threshold` is true

**On Ensemble Plots:**
- Markers appear on plot lines at years where difference is significant
- Size controlled by `p_value_marker_size` (default 10)
- Color typically matches the line color

### Controlling Statistical Analysis

**P-Value Threshold:**
```json
{
    "p_value_threshold": 0.05,  // Standard significance level
    "p_value_threshold": 0.01,  // More stringent
    "p_value_threshold": 0.1    // More lenient
}
```

**Output Control:**
```json
{
    "p_value_file_print_only_if_below_threshold": true,  // Only significant
    "p_value_file_print_only_if_below_threshold": false  // All p-values
}
```

**P-Value File Location:**
```json
{
    "p_value_file": "./plots/statistics/p_values.dat"  // Custom location
}
```

### Statistical Considerations

**Multiple Testing:**
- Be aware of multiple comparison problem
- With many categories/regions/scenarios, some significant results expected by chance
- Consider Bonferroni correction: Divide threshold by number of tests
- Example: 20 tests, use threshold = 0.05/20 = 0.0025

**Sample Size:**
- For individual plots: Sample size = number of years in range
- For ensemble plots: Sample size = number of ensemble members × number of years
- Larger sample sizes → more statistical power

**Assumptions:**
- T-test assumes normally distributed data
- Check data distributions if concerned
- For highly skewed data, consider transformations (e.g., log)

**Interpretation:**
- Significant p-value: Difference unlikely due to chance alone
- Non-significant: Cannot conclude there is a difference (doesn't prove no difference)
- Consider effect size along with statistical significance

---

## Plotting Customization Details

### Understanding Color Assignment

**For Individual Plots:**
- If 1-2 scenarios: Colors assigned by **category**
- If 3+ scenarios: Colors assigned by **scenario**

**For Ensemble Plots:**
- Colors always assigned by **scenario set**

**Color Palette (from utility_plots.py):**
The script uses a 21-color palette combining:
1. Matplotlib Tableau colors (10 colors, professional)
2. XKCD colors (11 colors, distinctive)

This ensures plots with up to 21 different scenarios/categories have unique colors.

**Example:**
```
Individual plot with 2 scenarios:
- Rice: Blue (first color)
- Wheat: Orange (second color)
- Corn: Green (third color)

Individual plot with 5 scenarios:
- Control: Blue (first color)
- Scenario_2: Orange (second color)
- Scenario_3: Green (third color)
...
```

### Understanding Line Styles

**Assignment:**
- Line styles (solid, dashed, dotted, etc.) are assigned by **region**
- First region: Solid line
- Second region: Dashed line
- Third region: Dotted line
- And so on...

**Available Styles (from utility_plots.py):**
12 distinct line styles ensure plots with up to 12 regions are clearly distinguishable.

**Example:**
```
Plot with regions Global, USA, China:
- Global: Solid line (──────)
- USA: Dashed line (– – – –)
- China: Dotted line (· · · ·)
```

### Understanding Markers

**Assignment:**
- Markers (shapes) are assigned by **category**
- First category: Circle (o)
- Second category: Triangle down (v)
- Third category: Triangle up (^)
- And so on...

**Available Markers (from utility_plots.py):**
22 distinct markers ensure plots with up to 22 categories have unique symbols.

**Example:**
```
Plot with categories Rice, Wheat, Corn:
- Rice: Circle markers (○)
- Wheat: Triangle down markers (▽)
- Corn: Triangle up markers (△)
```

### 3-Dimensional Visual Encoding

The script uses a sophisticated system to ensure visual distinction:

**3-Dimensional Visual Encoding:**
1. **Color** → Assigned by scenario (or category if few scenarios)
2. **Line Style** → Assigned by region
3. **Marker** → Assigned by category

**This allows clear visualization of:**
- Up to 21 scenarios/categories (via color from utility_plots.py)
- Up to 12 regions (via line style from utility_plots.py)
- Up to 22 categories (via markers from utility_plots.py)

**Total Theoretical Combinations:** 21 × 12 × 22 = 5,544 unique visual combinations!

**Example Visual Encoding:**
```
Scenario A (Blue) + Region USA (Solid) + Category Rice (Circle)
= Blue solid line with circle markers

Scenario B (Orange) + Region China (Dashed) + Category Wheat (Triangle)
= Orange dashed line with triangle markers
```

### Customizing for Different Media

**For Papers (7" × 5" typical journal figure):**
```json
{
    "width": 7,
    "height": 5,
    "linewidth": 1.5,
    "marker_size": 5,
    "x_label_size": 12,
    "y_label_size": 12,
    "x_tick_label_size": 10,
    "y_tick_label_size": 10,
    "legend_label_size": 9,
    "use_latex": true
}
```

**For Presentations (widescreen):**
```json
{
    "width": 12,
    "height": 6.75,
    "linewidth": 3,
    "marker_size": 10,
    "x_label_size": 32,
    "y_label_size": 32,
    "x_tick_label_size": 26,
    "y_tick_label_size": 26,
    "legend_label_size": 20,
    "use_latex": false,
    "produce_png": true
}
```

**For Posters:**
```json
{
    "width": 14,
    "height": 10,
    "linewidth": 4,
    "marker_size": 12,
    "x_label_size": 40,
    "y_label_size": 40,
    "x_tick_label_size": 34,
    "y_tick_label_size": 34,
    "legend_label_size": 24,
    "use_latex": false,
    "produce_png": true
}
```

---

## Output Files

### Plot Files

**PDF Format (default from utility_plots.py):**
- Vector graphics (scalable without quality loss)
- Small file size
- Professional quality
- Suitable for publications
- Default bounding box: `'tight'` (crops whitespace)

**PNG Format (optional via `produce_png: true`):**
- Raster graphics
- Resolution-dependent
- Typical DPI: 300 for publications, 150 for web, 72 for screen
- Larger file size than PDF
- Embedded in presentations

**File Naming:**
- PDF: `plot_name.pdf` or auto-generated name
- PNG: `plot_name.png` (same base name as PDF)

### P-Value Files

**Location:** `plot_directory/p_values.dat` (or custom location)

**Content:**
```
============================================================
Plot: ./plots/time_series_emissions.pdf
------------------------------------------------------------
scenario=Full feedback, category=agricultural energy use, region=Global: p=0.0123 *
scenario=Full feedback, category=construction energy use, region=Global: p=0.4567
scenario=Full feedback, category=trn_aviation_intl, region=Global: p=0.0089 **
scenario=Full feedback, category=refining, region=Global: p=0.1234
============================================================
```

**Interpretation:**
- Organized by plot
- One line per statistical comparison
- Asterisk (*) marks significant results (p < threshold)
- Double asterisk (**) for highly significant (p < 0.01, if used)

---

## Advanced Features

### Landtype Grouping

The script supports two levels of landtype aggregation:

**Modified Groups (default):**
- Aggregates GCAM landtypes into major categories
- Example: "crop" includes all crop types
- Defined in `utility_functions.py`

**Original Groups:**
- Uses original GCAM landtype names
- Set with `"landtype_groups": "original"`

**Landtype Group Examples:**
```python
crop = [BioenergyCrop, Corn, FodderGrass, Rice, Wheat, Soybean, ...]
forest = [Forest, ProtectedUnmanagedForest, UnmanagedForest]
shrub = [Shrubland, ProtectedUnmanagedShrubland, UnmanagedShrubland]
pasture = [Pasture, ProtectedUnmanagedPasture]
grass = [Grassland, ProtectedUnmanagedGrassland, UnmanagedGrassland]
```

---

### Area-Weighted Calculations

When data includes an `area` column:

**Area-Weighted Mean:**
```python
weighted_mean = sum(value × area) / sum(area)
```

**Use Cases:**
- Agricultural yields weighted by cropland area
- Vegetation scalars weighted by land cover
- Any intensive property where area matters

**Requirements:**
- CSV file must have `area` column
- Added using `gcam_add_areas_to_files.py`

**Example:**
```
Year 2020:
Region A: yield=5 t/ha, area=100 km²
Region B: yield=7 t/ha, area=50 km²

Simple mean: (5+7)/2 = 6 t/ha
Area-weighted mean: (5×100 + 7×50)/(100+50) = 5.67 t/ha
```

---

### Parallel Processing

The script processes multiple plot configurations in parallel:

**Command:**
```bash
python gcam_plot_time_series.py config1.json config2.json config3.json
```

**Behavior:**
- Each JSON configuration file processed independently
- All plots within a file created sequentially
- Utilizes parallel processing with limited cores to reduce memory pressure
- Significantly faster for many plots

---

### Error Handling

**NaN Values in Area-Weighted Means:**

When total area is zero, area-weighted mean becomes NaN.

**Option 1 - Set to Zero:**
```json
"set_nan_to_zero": true
```

**Option 2 - Skip (default):**
```json
"set_nan_to_zero": false
```
Years with NaN are not plotted (gaps in line).

**Console Warning:**
```
Area-weighted mean resulted in NaN because total area is zero for some years.
No fix applied as set_nan_to_zero is False.
Years with NaN values will not appear on the plot.
```

---

## Troubleshooting

### Common Issues

#### Issue 1: File Not Found Error
**Error:** `FileNotFoundError: output_file not found`

**Solutions:**
- Verify file path is correct
- Use absolute paths if unsure: `/full/path/to/file.csv`
- Check file has been processed by `gcam_process_extracted_data.py`
- Ensure CSV extension is included in filename

---

#### Issue 2: Column Not Found
**Error:** `KeyError: 'landtype'`

**Causes:**
- `category_label` doesn't match actual column name
- File structure different than expected

**Solutions:**
- Check column names in Python:
```python
import pandas as pd
df = pd.read_csv('your_file.csv')
print(df.columns.tolist())
```
- Verify `category_label` spelling and case
- Ensure file was processed correctly

---

#### Issue 3: LaTeX Errors
**Error:** `RuntimeError: Failed to process string with tex`

**Causes:**
- LaTeX not installed
- Invalid LaTeX syntax in labels

**Solutions:**
- Install LaTeX or set `"use_latex": false`
- Escape special characters in labels (use `\$` for `$`)
- Check y_label for proper LaTeX formatting
- Test without LaTeX first, then add it back

**Common LaTeX Issues:**
```json
// Wrong
"y_label": "Area (km^2)"

// Correct
"y_label": "Area (km$^2$)"
```

---

#### Issue 4: Empty or Missing Categories
**Error:** No data appears on plot

**Causes:**
- Category names don't match CSV file
- All categories excluded
- Year range doesn't overlap with data

**Solutions:**
- Print unique categories:
```python
import pandas as pd
df = pd.read_csv('your_file.csv')
print(df['sector'].unique())  # or whatever category_label is
```
- Verify category names match exactly (case-sensitive)
- Check `start_year` and `end_year` overlap with data
- Ensure `categories_to_exclude` isn't too broad

---

#### Issue 5: Legend Overlap
**Symptom:** Legend covers plot data

**Solutions:**
```json
"legend_place_outside": true,
"legend_x_offset": 1.4,
"legend_num_columns": 1
```

**Or increase figure width:**
```json
"width": 12,  // Was 10
"legend_place_outside": true
```

**Or reduce legend font size:**
```json
"legend_label_size": 10  // Was 14
```

---

#### Issue 6: Ensemble Structure Mismatch
**Error:** `ValueError: All inner lists must have same length`

**Cause:** Nested scenario lists have different lengths

**Solution:**
```json
// Incorrect
"scenarios": [
    ["Control", "Full feedback"],
    ["Control_2"]  // Missing Full feedback_2
]

// Correct
"scenarios": [
    ["Control", "Full feedback"],
    ["Control_2", "Full feedback_2"]
]
```

---

#### Issue 7: Plot Text Too Small/Large

**Symptom:** Text unreadable or overwhelming

**Solution - Increase all font sizes proportionally:**
```json
{
    "x_label_size": 32,      // Was 24
    "y_label_size": 32,      // Was 24
    "x_tick_label_size": 26, // Was 20
    "y_tick_label_size": 26, // Was 20
    "legend_label_size": 18  // Was 14
}
```

**Or use multiplier approach:**
- Publication × 1.0 (use defaults from utility_plots.py)
- Presentation × 1.3-1.4
- Poster × 1.5-2.0

---

#### Issue 8: Colors Look Similar

**Cause:** Using more than 21 scenarios/categories (exceeds color palette from utility_plots.py)

**Solutions:**
1. Reduce number of scenarios on single plot
2. Create multiple plots
3. Define custom color palette:
```json
"plot_colors": [
    "#FF0000", "#00FF00", "#0000FF",  // RGB primaries
    "#FFFF00", "#FF00FF", "#00FFFF",  // CMY secondaries
    "#FF8000", "#8000FF", "#00FF80",  // Additional colors
    // Add more as needed
]
```

**Color Selection Tools:**
- [Coolors.co](https://coolors.co) - Generate color palettes
- [Color Oracle](https://colororacle.org) - Test for color blindness
- [ColorBrewer](https://colorbrewer2.org) - Scientific color schemes

---

#### Issue 9: Lines Too Thin/Thick

**Solution:**
Adjust `linewidth` (default 2 from utility_plots.py):

```json
"linewidth": 1.5  // Thinner
"linewidth": 3    // Thicker
```

**Guidance:**
- Multi-line plots (many scenarios): 1-1.5
- Standard plots: 2-2.5 (default)
- Presentations: 3-4
- Posters: 4-5

---

#### Issue 10: Area-Weighted Mean Gives NaN

**Error Message:**
```
Area-weighted mean resulted in NaN because total area is zero for some years.
```

**Cause:** No area data for certain year/category/region combinations

**Solutions:**
1. Check if area column exists and has values
2. Verify `gcam_add_areas_to_files.py` was run successfully
3. Use simple mean instead:
```json
"aggregation_type_in_each_year": "mean"
```
4. Or set NaN to zero:
```json
"set_nan_to_zero": true
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
│   ├── time_series/
│   │   ├── individual_plots/
│   │   ├── ensemble_plots/
│   │   └── percent_difference/
│   └── box_plots/
├── configs/
│   ├── time_series_config.json
│   └── ensemble_config.json
└── statistics/
    └── p_values/
```

---

### 2. Naming Conventions

**Scenario Names:**
- Use descriptive names: "Control", "Full feedback", not "S1", "S2"
- Be consistent across ensemble members: "Control", "Control_2", not "Control", "Control_v2"

**File Names:**
- Include scenario info: `emissions_control_vs_feedback.pdf`
- Include date for versions: `emissions_2026-01-20.pdf`
- Use underscores not spaces: `time_series_vegetation_scalars.pdf`

---

### 3. Configuration Management

**Separate Configurations:**
- One config per analysis type
- Individual plots vs. ensemble plots
- Different time periods
- Different categories/regions

**Documentation:**
- Comment JSON files using external notes
- Keep track of which configs produced which figures
- Version control configuration files

**Example Organization:**
```
configs/
├── publication/
│   ├── fig1_emissions.json
│   ├── fig2_land_use.json
│   └── fig3_scalars.json
├── exploratory/
│   ├── test_all_scenarios.json
│   └── regional_detail.json
└── presentation/
    └── talk_2026-01-20.json
```

---

### 4. Choosing Figure Dimensions

**Start with defaults (from utility_plots.py):**
- Width: 10 inches
- Height: 8 inches
- Aspect ratio: 5:4

**Adjust based on:**
1. **Publication requirements:** Check journal guidelines
2. **Number of panels:** Multi-panel figures may need different ratios
3. **Legend size:** More legend items → wider figure
4. **Data density:** Dense data → larger figure

**Common Dimensions:**
```json
// Journal single column
{"width": 3.5, "height": 3}

// Journal double column
{"width": 7, "height": 5}

// Presentation slide (16:9)
{"width": 12, "height": 6.75}

// Presentation slide (4:3)
{"width": 12, "height": 9}

// Poster
{"width": 14, "height": 10}
```

---

### 5. Color Scheme Best Practices

**Use Default Palette (from utility_plots.py):**
The 21-color palette is carefully chosen for:
- Distinction (colors easily told apart)
- Accessibility (some consideration for color blindness)
- Professional appearance

**When to Customize Colors:**
1. Specific color requirements (e.g., blue for water, green for vegetation)
2. More than 21 scenarios needed
3. Matching colors across multiple figures
4. Institutional/publication requirements

**Custom Color Tips:**
- Use color picker tools (e.g., coolors.co)
- Test for color blindness accessibility (e.g., Color Oracle)
- Maintain sufficient contrast
- Use consistent colors across related figures

**Example Custom Palette:**
```json
"plot_colors": [
    "#0066CC",  // Control (blue)
    "#FF6600",  // Feedback (orange)
    "#00CC66",  // Scenario 3 (green)
    "#CC0066"   // Scenario 4 (magenta)
]
```

---

### 6. Statistical Rigor

**Multiple Testing:**
- Be aware of multiple comparison problem
- With many categories/regions/scenarios, some significant results expected by chance (5% with p=0.05)
- Consider Bonferroni correction: threshold = 0.05 / number_of_tests
- Example: 20 tests, use threshold = 0.05/20 = 0.0025

**Reporting:**
```json
{
    "p_value_threshold": 0.0025,  // Bonferroni corrected
    "p_value_file_print_only_if_below_threshold": false  // Report all
}
```

**Error Bars:**
- Use appropriate `std_multiplier` (1 for 1σ, 2 for ~95% CI)
- Understand what error bars represent:
  - Ensemble plots: Spread across ensemble members
  - Not measurement uncertainty
  - Not confidence intervals (unless specifically calculated)

**Sample Size:**
- Larger sample sizes → more statistical power
- Individual plots: n = number of years
- Ensemble plots: n = number of members × number of years

---

### 7. Data Preparation

**Before Plotting:**
1. Process with `gcam_process_extracted_data.py`
2. Add areas with `gcam_add_areas_to_files.py` (if needed for area-weighted means)
3. Verify data structure matches expectations
4. Check for missing values or outliers

**Verification Steps:**
```python
import pandas as pd

# Load data
df = pd.read_csv('your_file.csv')

# Check structure
print(df.columns)
print(df.head())

# Check years
print(f"Years: {df['year'].min()} to {df['year'].max()}")

# Check scenarios
print(f"Scenarios: {df['scenario'].unique()}")

# Check categories
print(f"Categories: {df['sector'].unique()}")

# Check for missing values
print(f"Missing values:\n{df.isnull().sum()}")
```

---

### 8. Iterative Workflow

**Recommended Process:**
1. **Start simple:** One scenario, one category, short time range
2. **Verify plot:** Check it looks reasonable
3. **Add complexity:** More scenarios, categories, regions
4. **Adjust styling:** Font sizes, colors, legend position
5. **Test different formats:** Try log scale, percent difference
6. **Run full analysis:** Complete time range, all categories
7. **Generate final versions:** Publication-quality settings

**Example Progression:**
```json
// Step 1: Test
{
    "scenarios": ["Control"],
    "categories": ["Rice"],
    "start_year": 2015,
    "end_year": 2030,
    "use_latex": false
}

// Step 2: Add complexity
{
    "scenarios": ["Control", "Full feedback"],
    "categories": ["Rice", "Wheat", "Corn"],
    "start_year": 2015,
    "end_year": 2050
}

// Step 3: Final version
{
    "scenarios": ["Control", "Full feedback"],
    "categories": ["Rice", "Wheat", "Corn", "Soybean"],
    "start_year": 2015,
    "end_year": 2100,
    "use_latex": true,
    "width": 7,
    "height": 5
}
```

---

### 9. Performance Optimization

**For Faster Execution:**
1. Disable LaTeX for drafts: `"use_latex": false`
2. Use PDF only (no PNG): `"produce_png": false`
3. Reduce year range for testing: `"start_year": 2015, "end_year": 2050`
4. Process multiple configs in parallel:
```bash
python gcam_plot_time_series.py config1.json config2.json config3.json
```

**Memory Considerations:**
- Large ensemble plots (many members, categories, regions) use more memory
- If memory issues occur, reduce number of elements or split into multiple plots

---

### 10. Quality Control Checklist

Before finalizing plots:

**Visual:**
- [ ] All text legible at final size
- [ ] Colors distinguishable
- [ ] Legend complete and positioned well
- [ ] Axis labels clear and properly formatted
- [ ] Line styles/markers distinct
- [ ] No data obscured

**Technical:**
- [ ] Correct scenarios included
- [ ] Appropriate aggregation method used
- [ ] Correct time range
- [ ] Units specified in labels
- [ ] Statistical significance properly marked

**Files:**
- [ ] Saved to correct directory
- [ ] Named descriptively
- [ ] PDF format for publications
- [ ] PNG format for presentations (if needed)
- [ ] Configuration file saved for reproducibility

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
   
4. gcam_plot_time_series.py  ← THIS SCRIPT
   ↓ (Create time series plots)
   
5. gcam_plot_box_and_whiskers.py (optional)
   ↓ (Create distribution plots)
   
6. Analysis and publication
```

### Related Scripts

**`gcam_plot_box_and_whiskers.py`**
- Complementary visualization
- Shows distributions instead of time series
- Uses similar configuration structure
- Good for showing variability across time

**`gcam_plot_spatial_data.py`**
- Spatial maps of GCAM data
- Geographic visualization
- Shows regional patterns

**Together they provide:**
- Time series: Temporal evolution
- Box plots: Statistical distribution
- Spatial plots: Geographic patterns

---

## Performance Considerations

### Execution Time

**Factors:**
- Number of configurations
- Number of scenarios/ensemble members
- Number of categories and regions
- File size (number of rows)
- Year range
- LaTeX rendering (slower than plain text)
- PDF vs PNG generation

**Typical Times:**
- Simple individual plot: 1-3 seconds
- Ensemble plot with 5 members: 3-8 seconds
- Complex multi-category, multi-region: 10-30 seconds
- With LaTeX: Add 20-50% to times above

**Optimization:**
- Use PNG instead of PDF for faster rendering (but lower quality)
- Disable LaTeX for draft plots (`"use_latex": false`)
- Reduce year range for testing (`"end_year": 2050` instead of 2100)
- Process multiple configs in parallel

**Memory Usage:**
- Typical plot: 50-200 MB RAM
- Large ensemble plots: Up to 500 MB RAM
- Python process overhead: ~100 MB

---

## Appendix: Complete Default Values

### From gcam_plot_time_series.py

```python
{
    'aggregation_type_in_each_year': 'area_weighted_mean',
    'category_label': 'sector',
    'end_year': 2100,
    'error_bars_alpha': 0.2,
    'include_mean_across_all_data': False,
    'key_columns': None,
    'landtype_groups': 'modified',
    'legend_on': True,
    'legend_place_outside': True,  # Script default (overrides utility_plots)
    'legend_x_offset': None,  # Auto-calculated by utility_plots.py
    'marker_size': 6,
    'mean_or_sum_if_more_than_one_row_in_same_landtype_group': 'area_weighted_mean',
    'multiplier': 1,
    'p_value_file': 'p_values.dat',
    'p_value_file_print_only_if_below_threshold': True,
    'p_value_marker_size': 10,
    'p_value_threshold': 0.05,
    'plot_directory': './',
    'plot_percent_difference': False,
    'plot_type': 'ensemble',
    'region_label': 'region',
    'regions': ['Global'],
    'scenario_label': 'scenario',
    'scenario_sets': None,
    'set_nan_to_zero': False,
    'start_year': 2015,
    'std_mean_across_all_data_multiplier': 1,
    'std_multiplier': 1,
    'value_label': 'value',
    'year_label': 'year'
}
```

### From utility_plots.py

```python
{
    'width': 10,  # inches
    'height': 8,  # inches
    'scale_default': 'linear',  # 'linear' or 'log'
    'axis_limits_default': None,  # Auto-scaled
    'axis_label_size': 24,  # points
    'tick_label_size': 20,  # points
    'legend_label_size': 14,  # points
    'legend_num_columns': 1,
    'legend_on': True,
    'legend_place_outside': False,  # Inside by default
    'linewidth': 2,  # points
    'produce_png': False,  # Only PDF by default
    'use_latex': False,  # Plain text by default
    'bbox_inches': 'tight',  # Crop whitespace
    
    # 21-color palette (10 Tableau + 11 XKCD)
    'plot_colors': [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        '#029386', '#c20078', '#53fca1', '#fe01b1', '#c65102',
        '#fac205', '#0b5509', '#8a6e45', '#fc5a50', '#a2cffe', '#ffb07c'
    ],
    
    # 12 line styles
    'linestyle_tuples': [
        ('solid', (0, ())),
        ('dashed', (0, (5, 5))),
        ('dotted', (0, (1, 5))),
        ('dashdot', (0, (3, 5, 1, 5))),
        ('loosely dotted', (0, (1, 10))),
        ('loosely dashed', (0, (5, 10))),
        ('densely dashed', (0, (5, 1))),
        ('loosely dashdotted', (0, (3, 10, 1, 10))),
        ('densely dashdotted', (0, (3, 1, 1, 1))),
        ('dashdotdotted', (0, (3, 5, 1, 5, 1, 5))),
        ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
        ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))
    ],
    
    # 22 markers
    'markers': [
        'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H',
        'D', 'd', 'P', 'X', '1', '2', '3', '4', '+', 'x', '|'
    ]
}
```

---

## Quick Reference: Common Configurations

### Minimal Configuration

```json
{
    "output_file": "./data/file.csv",
    "scenarios": ["Control", "Full feedback"]
}
```

### Standard Configuration

```json
{
    "output_file": "./data/file.csv",
    "scenarios": ["Control", "Full feedback"],
    "categories": ["Rice", "Wheat", "Corn"],
    "plot_directory": "./plots/",
    "y_scale": "log",
    "use_latex": true
}
```

### Publication Configuration

```json
{
    "output_file": "./data/file.csv",
    "scenarios": ["Control", "Full feedback"],
    "categories": ["Rice", "Wheat", "Corn"],
    "plot_directory": "./plots/publication/",
    "width": 7,
    "height": 5,
    "linewidth": 1.5,
    "marker_size": 5,
    "x_label_size": 12,
    "y_label_size": 12,
    "x_tick_label_size": 10,
    "y_tick_label_size": 10,
    "legend_label_size": 9,
    "use_latex": true,
    "produce_png": false
}
```

### Presentation Configuration

```json
{
    "output_file": "./data/file.csv",
    "scenarios": ["Control", "Full feedback"],
    "categories": ["Rice", "Wheat", "Corn"],
    "plot_directory": "./plots/presentation/",
    "width": 12,
    "height": 6.75,
    "linewidth": 3,
    "marker_size": 10,
    "x_label_size": 32,
    "y_label_size": 32,
    "x_tick_label_size": 26,
    "y_tick_label_size": 26,
    "legend_label_size": 20,
    "use_latex": false,
    "produce_png": true
}
```

### Full Ensemble Configuration

```json
{
    "output_file": "./data/file_ensemble.csv",
    "scenarios": [
        ["Control", "Control_2", "Control_3"],
        ["Full feedback", "Full feedback_2", "Full feedback_3"]
    ],
    "scenario_sets": ["Control", "Full feedback"],
    "categories": ["Rice", "Wheat", "Corn"],
    "regions": ["USA", "China", "Brazil"],
    "plot_directory": "./plots/ensemble/",
    "y_label": "Commodity price ($/kg)",
    "y_scale": "log",
    "start_year": 2020,
    "end_year": 2100,
    "std_multiplier": 1,
    "p_value_threshold": 0.05,
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
- Matplotlib Documentation: [https://matplotlib.org/](https://matplotlib.org/)
- Matplotlib Color Reference: [https://matplotlib.org/stable/gallery/color/named_colors.html](https://matplotlib.org/stable/gallery/color/named_colors.html)
- Matplotlib Markers: [https://matplotlib.org/stable/gallery/lines_bars_and_markers/marker_reference.html](https://matplotlib.org/stable/gallery/lines_bars_and_markers/marker_reference.html)
- Matplotlib Line Styles: [https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html](https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html)

---

## Contact

For questions or issues:
- **Philip Myint**: myint1@llnl.gov
- **Dalei Hao**: dalei.hao@pnnl.gov

---

## Version Information

**Script:** gcam_plot_time_series.py  
**Utility Module:** utility_plots.py  
**Documentation Version:** 2.0  
**Last Updated:** January 2026  
**Python Version:** 3.7+

---

*This documentation provides comprehensive guidance for using the `gcam_plot_time_series.py` script, including detailed default values from `utility_plots.py` for precise control over plot appearance.*
