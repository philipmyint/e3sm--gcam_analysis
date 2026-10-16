# GCAM Time-Series Plotting

**Script:** `gcam_plot_time_series.py`

## Overview

`gcam_plot_time_series.py` creates line plots from processed GCAM CSV files. It supports individual scenario plots, ensemble plots, percent differences, landtype group aggregation, multi-column categories, p-value output, and optional uncertainty bands.

The script expects at least:

```text
scenario, region, year, value
```

Category columns such as `sector` or `landtype` are optional. Regional scalar files such as `hdd_scalars.csv` and `cdd_scalars.csv` can be plotted without category settings.

## Usage

```bash
python gcam_plot_time_series.py path/to/config.json
```

Multiple JSON files can be supplied:

```bash
python gcam_plot_time_series.py config1.json config2.json
```

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `output_file` | Yes | | Processed CSV file to plot. Missing or empty files are skipped. |
| `scenarios` | Yes | | Scenario list, or nested list for ensemble plots. |
| `scenario_sets` | No | `None` | Labels for scenario groups when `scenarios` is nested. |
| `plot_type` | No | `"ensemble"` | `"ensemble"` groups nested scenarios by scenario set; `"individual"` plots every scenario separately. |
| `plot_percent_difference` | No | `false` | Plot percent difference from the first scenario or scenario set. |
| `start_year` | No | `2015` | First year to plot. |
| `end_year` | No | `2050` | Last year to plot. |
| `category_label` | No | `"sector"` | Category column, or list of category columns. |
| `categories` | No | Auto-detected or `"All"` | Categories to plot. |
| `regions` | No | `["Global"]` | Regions to include. |
| `aggregation_type_in_each_year` | No | `"sum"` | `"sum"`, `"mean"`, or `"area_weighted_mean"` across rows in a year. |
| `mean_or_sum_if_more_than_one_row_in_same_landtype_group` | No | `"sum"` | Landtype-group aggregation; still the active option name for this script. |
| `key_columns` | No | `None` | Grouping columns used when aggregating GCAM landtype groups. |
| `value_label` | No | `"value"` | Numeric column to plot. |
| `units_label` | No | `"units"` | Units column used for labels. |
| `max_processes` | No | `MAX_PROCESSES` | Maximum parallel workers. Defaults to all CPUs minus one, with a minimum of one. |

Styling options include `plot_directory`, `plot_name`, `plot_colors`, `markers`, `linestyle_tuples`, `linewidth`, `marker_size`, `y_label`, `y_label_abbreviated`, `x_limits`, `y_limits`, `legend_on`, `legend_place_outside`, `use_latex`, `produce_png`, `width`, and `height`.

## Scenario Formats

Flat scenario lists are plotted as individual scenario curves:

```json
"scenarios": ["Control", "Full feedback"]
```

Nested scenario lists are ensemble inputs:

```json
"scenarios": [
    ["Control", "Control_2", "Control_3"],
    ["Full feedback", "Full feedback_2", "Full feedback_3"]
],
"scenario_sets": ["Control", "Full feedback"]
```

For percent-difference plots, the first scenario or scenario set is the reference.

## Categories

`category_label` can be a single column:

```json
"category_label": "landtype",
"categories": ["crop", "forest", "pasture"]
```

It can also be a list of columns for multi-column categories:

```json
"category_label": ["sector", "subsector", "technology"],
"categories": {
    "sector": ["electricity"],
    "subsector": ["coal", "gas"],
    "technology": ["generation"]
}
```

For multi-column categories, the legend uses color, marker, and linewidth to encode up to three category columns.

## Aggregation

`aggregation_type_in_each_year` controls how rows are collapsed within each year:

| Value | Behavior | Typical use |
|-------|----------|-------------|
| `sum` | Sum values. | Emissions, land area, production. |
| `mean` | Simple average. | Equal-row averages. |
| `area_weighted_mean` | `sum(value * area) / sum(area)`. | Prices, scalars, rates. |

`area_weighted_mean` requires an `area` column, normally added by `gcam_process_extracted_data.py` or `gcam_compile_ehc_feedbacks.py`.

When requesting landtype groups such as `crop` or `forest`, use `mean_or_sum_if_more_than_one_row_in_same_landtype_group` to control the group aggregation. This option name is intentionally unchanged for the time-series script.

## Examples

### Vegetation scalars

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/vegetation_scalars.csv",
    "scenarios": ["Control", "Full feedback"],
    "category_label": "landtype",
    "categories": ["crop", "forest", "shrub", "pasture", "grass"],
    "regions": ["Global"],
    "aggregation_type_in_each_year": "area_weighted_mean",
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "area_weighted_mean",
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/time_series_plots",
    "plot_name": "time_series_vegetation_scalars.pdf",
    "y_label": "vegetation scalars"
}
```

### Soil scalar percent difference

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/soil_scalars.csv",
    "scenarios": ["Control", "Full feedback"],
    "category_label": "landtype",
    "categories": ["crop", "forest", "shrub", "pasture", "grass"],
    "plot_percent_difference": true,
    "aggregation_type_in_each_year": "area_weighted_mean",
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "area_weighted_mean"
}
```

### HDD scalar time series

```json
{
    "output_file": "./output/gcam_csv/hdd_scalars.csv",
    "scenarios": ["control", "case1"],
    "plot_directory": "./output/time_series_plots_gcam/individual",
    "aggregation_type_in_each_year": "area_weighted_mean",
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "area_weighted_mean"
}
```

### Ensemble percent difference

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/vegetation_scalars_ensemble.csv",
    "scenarios": [
        ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
        ["Full feedback", "Full feedback_2", "Full feedback_3", "Full feedback_4", "Full feedback_5"]
    ],
    "scenario_sets": ["Control", "Full feedback"],
    "category_label": "landtype",
    "categories": ["crop", "forest"],
    "plot_percent_difference": true,
    "aggregation_type_in_each_year": "area_weighted_mean",
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "area_weighted_mean"
}
```

## P-Values

For ensemble comparisons with multiple scenario sets, the script performs t-tests and writes results to `p_value_file`. Use `p_value_file_print_only_if_below_threshold` and `p_value_threshold` to control printed output.

## Workflow Notes

Run `gcam_process_extracted_data.py` first for ordinary GCAM query outputs. Run `gcam_compile_ehc_feedbacks.py` for scalar and diagnostic outputs. Area enrichment now happens inside those processing steps.
