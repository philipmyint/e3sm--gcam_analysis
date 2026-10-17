# GCAM Spatial Plotting

**Script:** `gcam_plot_spatial_data.py`

## Overview

`gcam_plot_spatial_data.py` creates map plots from processed GCAM CSV files and GCAM region/basin shapefiles. It supports single-scenario maps, means across scenarios, absolute differences, percent differences, ensemble comparisons, year panels, category panels, and significance stippling.

The script expects processed data with at least:

```text
scenario, region, year, value
```

Files with basin-level values can include `basin`. Files with categories can include a column such as `sector`, `landtype`, or a multi-column category such as `["sector", "subsector", "technology"]`.

## Usage

```bash
python gcam_plot_spatial_data.py path/to/config.json
```

Multiple JSON files can be supplied:

```bash
python gcam_plot_spatial_data.py config1.json config2.json
```

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `shape_file` | Yes | | GCAM boundary shapefile. |
| `output_file` | Yes | | Processed CSV file to plot. Missing or empty files are skipped. |
| `scenarios` | Yes | | Scenario list, or nested list for ensemble comparisons. |
| `scenario_sets` | No | `None` | Labels for scenario groups when `scenarios` is nested. |
| `plot_type` | No | `"mean"` | `"mean"`, `"absolute_difference"`, or `"percent_difference"`. A single scenario plots that scenario directly. |
| `plot_years` | No | `[2015, 2025, 2035, 2045]` | Years to plot. Replaces the older `start_year`/`end_year` behavior. |
| `plot_years_ncols` | No | `2` | Number of columns for year-panel plots. |
| `category_label` | No | `"sector"` | Category column, or list of category columns. |
| `categories` | No | Auto-detected or `"All"` | Data categories to include before plotting. |
| `plot_categories` | No | `None` | Categories to display as separate panels. Can be a list or a dict for multi-column categories. |
| `plot_categories_ncols` | No | `2` | Number of columns for category-panel plots. |
| `multipanel` | No | `true` | If `true`, `plot_years` and `plot_categories` produce panel figures. If `false`, each requested year/category is written as a separate single-map figure. |
| `aggregation_function` | No | `"sum"` | `"sum"`, `"mean"`, or `"area_weighted_mean"` when rows must collapse to one map value. |
| `aggregate_to_region` | No | `false` | Drop basin granularity and dissolve the shapefile to region polygons. |
| `key_columns` | No | `None` | Grouping columns used when aggregating GCAM landtype groups. |
| `landtype_groups` | No | `"modified"` | `"modified"` or `"original"` GCAM landtype grouping dictionary. |
| `value_label` | No | `"value"` | Numeric column to plot. |
| `units` | No | From `units` column when possible | Unit label for titles/colorbars. |
| `max_processes` | No | `MAX_PROCESSES` | Maximum parallel workers. Defaults to all CPUs minus one, with a minimum of one. |

Styling options include `plot_directory`, `plot_name`, `title`, `cmap`, `cbar_on`, `cbar_limits`, `linewidth`, `stippling_on`, `stippling_hatches`, `p_value_threshold`, `use_latex`, `produce_png`, `width`, and `height`.

## Scenario Formats

Individual or direct comparison plots use a flat list:

```json
"scenarios": ["Control", "Full feedback"]
```

Ensemble plots use a nested list. The script accepts either organization by scenario set or by ensemble member and transposes when needed:

```json
"scenarios": [
    ["Control", "Control_2", "Control_3"],
    ["Full feedback", "Full feedback_2", "Full feedback_3"]
],
"scenario_sets": ["Control", "Full feedback"]
```

For differences, the first scenario or scenario set is the reference.

## Plot Layouts

By default, `multipanel` is `true`. In this mode, if `plot_years` is set and `plot_categories` is omitted, the script creates one multi-panel figure with one map per year.

If `plot_categories` is set, the script creates one output file per year, with one map panel per requested category. When `plot_years` has one year, the configured `plot_name` is used directly.

Set `"multipanel": false` to split those same requests into separate single-map figures. The script appends the category and/or year to the configured `plot_name` so the individual outputs do not overwrite one another.

For regional scalar files such as `hdd_scalars.csv` and `cdd_scalars.csv`, omit `category_label`, `categories`, and `plot_categories` unless a real category column exists.

## Aggregation

`aggregation_function` controls how rows are collapsed into one value per scenario and map polygon:

| Value | Behavior | Typical use |
|-------|----------|-------------|
| `sum` | Sum rows. | Land area, emissions, production. |
| `mean` | Average rows. | Simple averages where all rows should count equally. |
| `area_weighted_mean` | Compute `sum(value * area) / sum(area)`. | Prices, scalars, rates, temperatures. |

`area_weighted_mean` requires an `area` column. Add it upstream with `gcam_process_extracted_data.py` or `gcam_compile_ehc_feedbacks.py`.

## Examples

### Single-scenario commodity price panels

```json
{
    "shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
    "scenarios": ["Control"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "plot_categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean"],
    "plot_years": [2045],
    "aggregation_function": "area_weighted_mean",
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/spatial_plots",
    "plot_name": "spatial_ag_commodity_prices_control_2045.pdf"
}
```

### Land allocation ensemble percent difference

```json
{
    "shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_ensemble.csv",
    "scenarios": [
        ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
        ["Full feedback", "Full feedback_2", "Full feedback_3", "Full feedback_4", "Full feedback_5"]
    ],
    "scenario_sets": ["Control", "Full feedback"],
    "category_label": "landtype",
    "plot_categories": ["crop", "forest", "pasture", "grass"],
    "plot_years": [2045],
    "plot_type": "percent_difference",
    "aggregation_function": "sum",
    "key_columns": ["scenario", "region", "basin", "year"]
}
```

### Vegetation scalar percent difference

```json
{
    "shape_file": "./../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./../2025_DiVittorio_et_al_gcam/vegetation_scalars.csv",
    "scenarios": ["Control", "Full feedback"],
    "category_label": "landtype",
    "plot_categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean"],
    "plot_years": [2045],
    "plot_type": "percent_difference",
    "aggregation_function": "area_weighted_mean",
    "cbar_limits": [-20, 60]
}
```

### Regional HDD scalar maps

```json
{
    "shape_file": "../2025_DiVittorio_et_al_gcam/gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
    "output_file": "./output/gcam_csv/hdd_scalars.csv",
    "scenarios": ["control", "case1"],
    "plot_years": [2025, 2035, 2045],
    "plot_type": "percent_difference",
    "aggregation_function": "area_weighted_mean"
}
```

## Significance Stippling

For ensemble comparisons with exactly two scenario sets, the script performs t-tests at each mapped polygon and can stipple significant regions. Global p-value summaries are written to `p_value_file`.

## Notes

Percent-difference maps use `(test - control) / control * 100`, with the first scenario or scenario set as the control. If no `cbar_limits` are supplied and values exceed 100 percent, the script caps the panel color scale at `[-100, 100]` to keep maps readable.
