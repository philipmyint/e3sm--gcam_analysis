# GCAM Box-And-Whisker Plotting

**Script:** `gcam_plot_box_and_whiskers.py`

## Overview

`gcam_plot_box_and_whiskers.py` creates box plots from processed GCAM CSV files. It can compare individual scenarios or scenario-set ensembles across categories, regions, basins, and years, either as multi-panel year figures or as separate single-year figures.

The script expects at least:

```text
scenario, region, year, value
```

Category-rich files usually include `sector` or `landtype`. Regional scalar files such as `hdd_scalars.csv` and `cdd_scalars.csv` do not have those category columns, so their box plot JSON should set `x_variable` explicitly, usually to `"scenario"` or `"region"`.

## Usage

```bash
python gcam_plot_box_and_whiskers.py path/to/config.json
```

Multiple JSON files can be supplied:

```bash
python gcam_plot_box_and_whiskers.py config1.json config2.json
```

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `output_file` | Yes | | Processed CSV file to plot. Missing or empty files are skipped. |
| `scenarios` | Yes | | Scenario list, or nested list for ensemble plots. |
| `scenario_sets` | No | `None` | Labels for scenario groups when `scenarios` is nested. |
| `plot_type` | No | `"ensemble"` | `"ensemble"` groups nested scenarios by scenario set; `"individual"` plots every scenario separately. |
| `plot_percent_difference` | No | `false` | Plot percent difference from the first scenario or scenario set. |
| `plot_years` | No | `[2015, 2025, 2035, 2045]` | Years to show as panels. |
| `plot_years_ncols` | No | `2` | Number of columns for year panels. |
| `multipanel` | No | `true` | If `true`, requested years are shown as panels in one figure. If `false`, each requested year is written as a separate single-year figure. |
| `category_label` | No | `"sector"` | Category column used when `x_variable` is not supplied. |
| `categories` | No | Auto-detected | Categories to include. |
| `regions` | No | `["Global"]` | Regions to include. `"Global"` pools all regions. |
| `basins` | No | `None` | Basin abbreviations to include. |
| `x_variable` | No | `category_label` | Column used on the x-axis. |
| `hue` | No | `None` | Column used to color boxes. Cannot be `"scenario"` when `plot_percent_difference` is true. |
| `aggregation_function` | No | `"sum"` | `"sum"`, `"mean"`, or `"area_weighted_mean"` for landtype group aggregation. |
| `key_columns` | No | `None` | Grouping columns used when aggregating GCAM landtype groups. |
| `landtype_groups` | No | `"modified"` | `"modified"` or `"original"` GCAM landtype grouping dictionary. |
| `value_label` | No | `"value"` | Numeric column to plot. |
| `max_processes` | No | `MAX_PROCESSES` | Maximum parallel workers. Defaults to all CPUs minus one, with a minimum of one. |

Styling options include `plot_directory`, `plot_name`, `plot_colors`, `fill_boxes`, `linewidth`, `marker_size`, `x_tick_labels`, `y_label`, `y_scale`, `y_limits`, `legend_on`, `legend_place_outside`, `use_latex`, `produce_png`, `width`, and `height`.

## Scenario Formats

Individual comparisons use a flat list:

```json
"scenarios": ["Control", "Full feedback"]
```

Ensemble plots use a nested list plus `scenario_sets`:

```json
"scenarios": [
    ["Control", "Control_2", "Control_3"],
    ["Full feedback", "Full feedback_2", "Full feedback_3"]
],
"scenario_sets": ["Control", "Full feedback"]
```

For percent-difference plots, the first scenario or scenario set is the reference and is not plotted.

## Percent Difference And Hue

Do not combine:

```json
"hue": "scenario",
"plot_percent_difference": true
```

The script intentionally skips that combination because percent-difference plots subtract the first scenario or scenario set and then drop it from the plotted data. Use the scenario order to define the comparison instead:

```json
{
    "scenarios": ["Control", "Full feedback"],
    "plot_percent_difference": true
}
```

## Aggregation

`aggregation_function` is used when a requested category is a GCAM landtype group such as `crop`, `forest`, `pasture`, `grass`, or `shrub`.

| Value | Behavior | Typical use |
|-------|----------|-------------|
| `sum` | Sum member landtypes. | Land area and totals. |
| `mean` | Average member landtypes equally. | Simple category averages. |
| `area_weighted_mean` | Compute `sum(value * area) / sum(area)`. | Prices, scalars, rates. |

`area_weighted_mean` requires an `area` column, normally added by `gcam_process_extracted_data.py` or `gcam_compile_ehc_feedbacks.py`.

## Examples

### Commodity prices by scenario

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
    "hue": "scenario",
    "scenarios": ["Control", "Full feedback"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "regions": ["USA", "Brazil", "Canada"],
    "aggregation_function": "area_weighted_mean",
    "y_scale": "log",
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/box_plots/individual"
}
```

### Commodity price percent difference

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
    "scenarios": ["Control", "Full feedback"],
    "categories": ["BioenergyCrop", "Rice", "SugarCrop", "Soybean", "Wheat"],
    "regions": ["USA", "Brazil", "Canada"],
    "aggregation_function": "area_weighted_mean",
    "plot_percent_difference": true,
    "plot_directory": "./../2025_DiVittorio_et_al_gcam/box_plots/individual"
}
```

### Vegetation scalar ensemble plot

```json
{
    "output_file": "./../2025_DiVittorio_et_al_gcam/vegetation_scalars_ensemble.csv",
    "hue": "scenario",
    "scenarios": [
        ["Control", "Control_2", "Control_3", "Control_4", "Control_5"],
        ["Full feedback", "Full feedback_2", "Full feedback_3", "Full feedback_4", "Full feedback_5"]
    ],
    "scenario_sets": ["Control", "Full feedback"],
    "category_label": "landtype",
    "categories": ["crop", "forest", "shrub", "pasture", "grass"],
    "key_columns": ["scenario", "region", "basin", "year"],
    "aggregation_function": "area_weighted_mean",
    "y_label": "Vegetation scalars"
}
```

### Regional HDD scalar box plot

```json
{
    "output_file": "./output/gcam_csv/hdd_scalars.csv",
    "scenarios": ["control", "case1"],
    "x_variable": "scenario",
    "plot_years": [2025, 2035, 2045],
    "regions": ["Global"],
    "aggregation_function": "area_weighted_mean",
    "plot_directory": "./output/box_plots_gcam"
}
```

Use `"x_variable": "region"` instead when the goal is to compare regional distributions.

## Notes

By default, `multipanel` is `true`, so each requested year is rendered as a panel in the same output figure. If `plot_years` is omitted, the default years are used. If `plot_years` is set to `null`, all years in the file are plotted.

For box-and-whisker plots, `multipanel` controls only the year layout:

| JSON options | Output |
|--------------|--------|
| `multipanel: true`, `plot_years: [2025, 2035, 2045]` | One output file with three panels: 2025, 2035, and 2045. |
| `multipanel: false`, `plot_years: [2025, 2035, 2045]` | Three separate single-year output files. |

For example, the `scripts/gcam_plot_box_and_whiskers.json` commodity-price multipanel block uses:

```json
"plot_years": [2025, 2035, 2045],
"plot_name": "box_plot_ag_commodity_prices_control_years_multipanel.pdf",
"multipanel": true
```

and produces one PDF with year panels. The nearby non-multipanel commodity-price examples use the same years with `"multipanel": false`, so the script writes separate files such as `box_plot_ag_commodity_prices_control_2025.pdf`, `box_plot_ag_commodity_prices_control_2035.pdf`, and `box_plot_ag_commodity_prices_control_2045.pdf`.

The example JSON in `scripts/` includes both multipanel and non-multipanel blocks.
