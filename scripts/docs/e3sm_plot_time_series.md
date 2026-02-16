# E3SM Time Series Plotting Script Documentation

## Overview

**Script Name:** `e3sm_plot_time_series.py`

**Purpose:** Creates comprehensive time series plots from E3SM extracted data files with support for annual means, seasonal averages, monthly climatologies, ensemble analysis, statistical significance testing, and percent difference comparisons.

**Key Capabilities:**
- Annual and seasonal time series plots
- Monthly climatology subplots
- Ensemble mean with uncertainty bands
- Individual ensemble member plots
- Statistical significance testing (multiple types)
- Percent difference plots
- Multiple scenarios comparison
- Flexible ensemble file organization

**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Requirements](#installation--requirements)
3. [Basic Usage](#basic-usage)
4. [Complete Parameter Reference](#complete-parameter-reference)
5. [Ensemble vs Individual Plots](#ensemble-vs-individual-plots)
6. [Statistical Testing Explained](#statistical-testing-explained)
7. [Ensemble File Organization](#ensemble-file-organization)
8. [Comprehensive JSON Examples](#comprehensive-json-examples)
9. [Output Files](#output-files)
10. [Best Practices](#best-practices)
11. [Troubleshooting](#troubleshooting)

---

## Installation & Requirements

### Required Python Libraries

```bash
pip install pandas numpy matplotlib scipy json multiprocessing
```

### Required Utility Modules

- `utility_constants` - Physical constants
- `utility_dataframes` - DataFrame operations
- `utility_functions` - General utilities
- `utility_plots` - Plotting defaults and functions

---

## Basic Usage

### Command Line Execution

```bash
python e3sm_plot_time_series.py config.json
```

**Multiple Configuration Files:**
```bash
python e3sm_plot_time_series.py config1.json config2.json
```

### What the Script Does

1. Reads JSON configuration file(s)
2. Loads time series data files (.dat or .csv)
3. Processes data (temporal aggregation, seasonal means)
4. Creates publication-quality plots
5. Performs statistical testing (multiple types available)
6. Saves plots as PDF or PNG
7. Outputs p-value files

---

## Complete Parameter Reference

### Required Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|---------|---------|-------------|
| `output_files` | list/nested list | **Yes** | - | Data file paths |
| `output_labels` | list | **Yes** | - | Labels for legend |
| `plot_directory` | string | No | `'./'` | Output directory for plots |

### Core Optional Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|---------|---------|-------------|
| `variables` | list or `'all'` | No | `'all'` | Variables to plot |
| `plot_type` | string | No | `'ensemble'` | `'ensemble'` or `'individual'` |
| `plot_percent_difference` | boolean | No | `false` | Plot % difference vs first file |
| `start_year` | integer | No | `2015` | First year to plot |
| `end_year` | integer | No | `2100` | Last year to plot |
| `use_latex` | boolean | No | `false` | Use LaTeX fonts |

### Visual Customization

| Parameter | Type | Required | Default | Description |
|-----------|------|---------|---------|-------------|
| `width` | float | No | `10` | Figure width (inches) |
| `height` | float | No | `8` | Figure height (inches) |
| `x_limits` | list `[min, max]` | No | `None` | X-axis limits |
| `y_limits` | list `[min, max]` | No | `None` | Y-axis limits |
| `y_label` | dict | No | Auto | Y-axis labels per variable |
| `linewidth` | float | No | `2` | Line width |
| `plot_colors` | list | No | Matplotlib default | Line colors |
| `legend_on` | boolean | No | `true` | Show legend |
| `legend_label_size` | integer | No | `14` | Legend font size |
| `x_label_size` | integer | No | `24` | X-axis label font size |
| `y_label_size` | integer | No | `24` | Y-axis label font size |
| `x_tick_label_size` | integer | No | `20` | X-axis tick label font size |
| `y_tick_label_size` | integer | No | `20` | Y-axis tick label font size |

### Seasonal Analysis

| Parameter | Type | Required | Default | Description |
|-----------|------|---------|---------|-------------|
| `include_seasons` | dict | No | All `false` | Include seasons on main plot |
| `seasons_to_plot_separately` | dict | No | All `false` | Create separate seasonal plots |
| `monthly_aggregation_type` | string/dict | No | `'mean'` | `'mean'` or `'sum'` |

**Seasons:**
- **MAM:** March-April-May (spring)
- **JJA:** June-July-August (summer)
- **SON:** September-October-November (fall)
- **DJF:** December-January-February (winter)

### Monthly Time Series

| Parameter | Type | Required | Default | Description |
|-----------|------|---------|---------|-------------|
| `monthly_time_series_plot` | boolean/dict | No | `false` | Create monthly climatology subplot |
| `monthly_time_series_start_year` | integer | No | `2071` | Start year for monthly plot |
| `monthly_time_series_end_year` | integer | No | `2090` | End year for monthly plot |
| `monthly_time_series_x_limits` | list | No | `None` | X-axis limits for monthly plot |
| `monthly_time_series_y_limits` | list | No | `None` | Y-axis limits for monthly plot |

### Statistical Testing

| Parameter | Type | Required | Default | Description |
|-----------|------|---------|---------|-------------|
| `p_value_threshold` | float | No | `0.05` | Significance threshold |
| `p_value_file` | string | No | `'p_values.dat'` | Output file for p-values |
| `p_value_file_print_only_if_below_threshold` | boolean | No | `true` | Only print significant p-values |
| `p_value_marker` | string | No | `'o'` | Marker for significant points |
| `p_value_marker_size` | integer | No | `6` | Size of significance markers |
| `include_annual_mean_across_all_sets` | boolean | No | `false` | Plot ensemble mean across all scenarios |
| `include_monthly_mean_across_all_sets` | boolean | No | `false` | Plot monthly ensemble mean |
| `std_annual_multiplier` | float | No | `1` | Error bar multiplier for annual data |
| `std_monthly_multiplier` | float | No | `1` | Error bar multiplier for monthly data |
| `std_annual_mean_across_all_sets_multiplier` | float | No | `1` | Error bar multiplier for ensemble mean |
| `std_monthly_mean_across_all_sets_multiplier` | float | No | `1` | Error bar multiplier for monthly ensemble mean |
| `std_seasons_multiplier` | float | No | `None` | Error bar multiplier for seasons |
| `error_bars_alpha` | float | No | `0.2` | Transparency of error bands |

### Advanced Options

| Parameter | Type | Required | Default | Description |
|-----------|------|---------|---------|-------------|
| `multiplier` | float/dict | No | `1` | Scale data values |
| `areas_in_thousands_km2` | boolean | No | `true` | Convert area units |
| `produce_png` | boolean | No | `false` | Output PNG instead of PDF |
| `notify_output_files_transposed` | boolean | No | `false` | Print message if ensemble format auto-converted |

---

## Ensemble vs Individual Plots

### Plot Type: `'ensemble'` (Default)

**Purpose:** Calculate and display ensemble mean with uncertainty quantification

**When to Use:**
- Multiple realizations available (ensemble members)
- Need statistical confidence in results
- Want to show uncertainty bands
- Publication-quality significance testing

**What It Does:**
1. Averages across ensemble members for each scenario
2. Calculates standard deviation/error bars
3. Performs statistical tests (t-tests)
4. Displays ensemble mean lines with shaded uncertainty
5. Marks statistically significant points

**Example Configuration:**
```json
{
    "output_files": [
        ["control.dat", "control_2.dat", "control_3.dat"],
        ["scenario.dat", "scenario_2.dat", "scenario_3.dat"]
    ],
    "output_labels": ["Control", "Scenario"],
    "plot_type": "ensemble"
}
```

**Result:**
- 2 lines (Control ensemble mean, Scenario ensemble mean)
- Shaded bands showing ±1σ uncertainty
- Markers at points where scenarios are significantly different (p < 0.05)

---

### Plot Type: `'individual'`

**Purpose:** Display each ensemble member as a separate line

**When to Use:**
- Want to see range of individual realizations
- Exploring ensemble spread
- Comparing specific ensemble members
- No need for statistical testing

**What It Does:**
1. Plots each file as a separate line
2. No averaging or statistics
3. Shows full range of variability
4. Good for understanding ensemble behavior

**Example Configuration:**
```json
{
    "output_files": [
        ["control.dat", "control_2.dat", "control_3.dat"],
        ["scenario.dat", "scenario_2.dat", "scenario_3.dat"]
    ],
    "output_labels": ["Control", "Scenario"],
    "plot_type": "individual"
}
```

**Result:**
- 6 lines total (3 Control + 3 Scenario)
- Each ensemble member visible
- No averaging or error bars
- No statistical testing

---

### Comparison Table

| Feature | Ensemble | Individual |
|---------|----------|------------|
| **Lines per scenario** | 1 (mean) | N (all members) |
| **Uncertainty bands** | Yes (±σ) | No |
| **Statistical testing** | Yes (multiple types) | Limited (aggregated only) |
| **Markers for significance** | Yes | No (ensemble) or Yes (aggregated) |
| **Use case** | Publication, significance | Exploration, range |
| **Recommended for** | Final analysis | Initial investigation |

---

## Statistical Testing Explained

The script performs **three types of statistical tests** depending on configuration:

### 1. Aggregated t-test (Always Available)

**What:** Tests if the overall means differ across entire time period

**When:** Both ensemble and individual plots

**How:**
- Takes mean of all values for Control across all time
- Takes mean of all values for Scenario across all time
- Performs two-sample t-test
- Records p-value in file

**Example Output:**
```
p_values.dat:
GPP in Scenario: 1.234e-04
```

**Interpretation:** Scenario is significantly different from Control over entire period

---

### 2. Per-Year t-test (Ensemble Only)

**What:** Tests if means differ at each individual year

**When:** `plot_type = 'ensemble'` with ≥2 members per scenario

**How:**
- At year 2050: Compare Control ensemble [ctrl, ctrl_2, ctrl_3] vs Scenario ensemble [scen, scen_2, scen_3]
- Performs t-test for that year
- If p < threshold, places marker at that year
- Repeats for each year

**Example:** Markers appear at years 2075, 2080, 2085, 2090 where p < 0.05

**Interpretation:** These specific years show significant differences

---

### 3. Per-Month t-test (Ensemble + Monthly Plot)

**What:** Tests if means differ for each month (averaged across years)

**When:** `plot_type = 'ensemble'`, `monthly_time_series_plot = true`, ≥2 members

**How:**
- For January: Average all Januaries (2071-2090), compare Control vs Scenario
- Performs t-test for January
- If p < threshold, places marker at January
- Repeats for each month

**Example:** Markers at June, July, August (summer months significantly different)

**Interpretation:** Seasonal differences identified

---

### Statistical Testing Summary

**Configuration Controls:**
- `p_value_threshold`: Significance level (default: 0.05)
- `p_value_marker`: Shape of significance markers (default: 'o')
- `p_value_file`: Where to save aggregated p-values

**What Gets Tested:**
```
Ensemble plot:
  ✓ Aggregated (whole period) - saved to file
  ✓ Per-year - markers on plot
  ✓ Per-month (if monthly plot enabled) - markers on monthly plot

Individual plot:
  ✓ Aggregated (whole period) - saved to file
  ✗ Per-year - not applicable
  ✗ Per-month - not applicable
```

**P-Value File Contents:**
```
GPP in Full feedback: 2.345e-05
NPP in Full feedback: 1.234e-03
NBP in Ag scaling: 8.765e-02
```

**Interpretation:**
- Values < 0.05: Statistically significant
- Values ≥ 0.05: Not statistically significant at 95% confidence

---

## Ensemble File Organization

When plotting ensembles, you can organize your file lists in **two equivalent ways**:

### Option 1: Organized by Scenario (Intuitive)

Each row groups all ensemble members for one scenario together:

```json
{
    "output_files": [
        ["control.dat", "control_2.dat", "control_3.dat", "control_4.dat", "control_5.dat"],
        ["feedback.dat", "feedback_2.dat", "feedback_3.dat", "feedback_4.dat", "feedback_5.dat"],
        ["scenario.dat", "scenario_2.dat", "scenario_3.dat", "scenario_4.dat", "scenario_5.dat"]
    ],
    "output_labels": ["Control", "Full feedback", "Scenario"]
}
```

**Structure:**
- Row 1: All Control files
- Row 2: All Full feedback files  
- Row 3: All Scenario files
- Easy to see which files belong to each scenario

---

### Option 2: Organized by Ensemble Member

Each row contains one member from each scenario:

```json
{
    "output_files": [
        ["control.dat", "feedback.dat", "scenario.dat"],
        ["control_2.dat", "feedback_2.dat", "scenario_2.dat"],
        ["control_3.dat", "feedback_3.dat", "scenario_3.dat"],
        ["control_4.dat", "feedback_4.dat", "scenario_4.dat"],
        ["control_5.dat", "feedback_5.dat", "scenario_5.dat"]
    ],
    "output_labels": ["Control", "Full feedback", "Scenario"]
}
```

**Structure:**
- Row 1: Member 1 for all scenarios
- Row 2: Member 2 for all scenarios
- Row 3: Member 3 for all scenarios
- Easy to add/remove ensemble members

---

### Automatic Detection

The script automatically detects which organization you're using by comparing:
1. Number of items in `output_labels`
2. Dimensions of `output_files` matrix

**Both produce identical results!** Use whichever is more intuitive for your workflow.

**Tip:** Always provide `output_labels` for unambiguous detection.

---

## Comprehensive JSON Examples

### Example 1: Simple Individual Scenarios

**From JSON file (Config #1):**
```json
{
    "output_files": [
        "./control_time_series.dat", 
        "./full_feedback_time_series.dat",
        "./ag_scaling_time_series.dat", 
        "./carbon_scaling_time_series.dat"
    ],
    "output_labels": ["Control", "Full feedback", "Ag scaling", "Carbon scaling"],
    "plot_directory": "./time_series_plots/individual",
    "y_label": {"ZCO2": "CO$_2$ concentration (ppm)"},
    "monthly_aggregation_type": {"PRECIP": "sum", "NPP": "sum"},
    "x_limits": {"ZCO2": [2015, 2100]},
    "y_limits": {"ZCO2": [400, 700]},
    "use_latex": true
}
```

**What This Creates:**

**Annual time series plots for all variables:**
- `time_series_GPP.pdf` - 4 lines (Control, Full feedback, Ag scaling, Carbon scaling)
- `time_series_NPP.pdf` - 4 lines
- `time_series_NBP.pdf` - 4 lines
- `time_series_ZCO2.pdf` - 4 lines, custom y-label, limits [400, 700] ppm
- `time_series_PRECIP.pdf` - 4 lines (monthly data summed to annual)
- ... (all other variables)

**Statistical Testing:**
- Aggregated t-tests comparing each scenario to Control
- P-values saved to `./time_series_plots/individual/p_values.dat`
- Example: `ZCO2 in Full feedback: 1.234e-05`

**No per-year markers** (individual plots don't show these)

---

### Example 2: Regional Comparison

**From JSON file (Config #2):**
```json
{
    "output_files": [
        "./control_time_series.dat", 
        "./control_time_series_amazon.dat"
    ],
    "output_labels": ["Control (Global)", "Control (Amazon)"],
    "plot_directory": "./time_series_plots/individual_amazon",
    "include_seasons": {"ZCO2": {"JJA": true, "DJF": true}},
    "seasons_to_plot_separately": {"ZCO2": {"MAM": true, "JJA": true, "SON": true, "DJF": true}},
    "monthly_time_series_plot": {"ZCO2": true}
}
```

**What This Creates:**

**For ZCO2 variable:**
- `time_series_ZCO2.pdf` - Main plot with:
  - Annual mean (2 lines: Global, Amazon)
  - JJA seasonal line (summer, 2 lines)
  - DJF seasonal line (winter, 2 lines)
  - Total: 6 lines on one plot

**Separate seasonal plots:**
- `time_series_ZCO2_MAM.pdf` - Spring only (2 lines)
- `time_series_ZCO2_JJA.pdf` - Summer only (2 lines)
- `time_series_ZCO2_SON.pdf` - Fall only (2 lines)
- `time_series_ZCO2_DJF.pdf` - Winter only (2 lines)

**Monthly climatology:**
- `time_series_ZCO2_monthly.pdf` - Months Jan-Dec on x-axis
  - Shows average seasonal cycle
  - 2 lines (Global, Amazon)

**For all other variables:**
- Standard annual plots (no seasons, no monthly)

---

### Example 3: Percent Difference

**From JSON file (Config #3):**
```json
{
    "output_files": [
        "./control_time_series.dat", 
        "./full_feedback_time_series.dat",
        "./ag_scaling_time_series.dat", 
        "./carbon_scaling_time_series.dat"
    ],
    "output_labels": ["Control", "Full feedback", "Ag scaling", "Carbon scaling"],
    "plot_directory": "./time_series_plots/individual_percent_difference",
    "plot_percent_difference": true,
    "y_limits": {"ZCO2": [-10, 10]}
}
```

**What This Creates:**

**Percent difference plots:**
- Control: Not shown (reference = 0%)
- Full feedback: Shows % change relative to Control
- Ag scaling: Shows % change relative to Control
- Carbon scaling: Shows % change relative to Control

**Example values:**
- Year 2050: Full feedback = +2.5% (2.5% higher than Control)
- Year 2075: Ag scaling = -1.8% (1.8% lower than Control)
- Year 2100: Carbon scaling = +5.2% (5.2% higher than Control)

**Y-axis for ZCO2:** [-10, 10]% (custom limits)

**Useful for:** Highlighting relative impacts rather than absolute values

---

### Example 4: Ensemble Analysis

**From JSON file (Config #4):**
```json
{
    "output_files": [
        ["./control.dat", "./control_2.dat", "./control_3.dat", "./control_4.dat", "./control_5.dat"],
        ["./feedback.dat", "./feedback_2.dat", "./feedback_3.dat", "./feedback_4.dat", "./feedback_5.dat"],
        ["./ag_scaling.dat", "./ag_scaling_2.dat", "./ag_scaling_3.dat", "./ag_scaling_4.dat", "./ag_scaling_5.dat"],
        ["./carbon.dat", "./carbon_2.dat", "./carbon_3.dat", "./carbon_4.dat", "./carbon_5.dat"]
    ],
    "output_labels": ["Control", "Full feedback", "Ag scaling", "Carbon scaling"],
    "plot_directory": "./time_series_plots/ensemble"
}
```

**Note:** This uses Option 1 format (organized by scenario)

**What This Creates:**

**Ensemble mean plots:**
- `time_series_GPP.pdf` - 4 lines (ensemble means)
  - Control: Mean of 5 members ± shaded uncertainty band
  - Full feedback: Mean of 5 members ± shaded band
  - Ag scaling: Mean of 5 members ± shaded band
  - Carbon scaling: Mean of 5 members ± shaded band
  - Markers (circles) at years where p < 0.05

**Example interpretation:**
- Years 2060-2070: No markers → scenarios not significantly different
- Years 2075-2100: Markers on Full feedback line → significantly different from Control

**Statistical Testing:**
- **Per-year t-tests:** Markers on plot at significant years
- **Aggregated t-tests:** In p_values.dat file
- **Example p-value file:**
  ```
  GPP in Full feedback: 1.23e-06
  GPP in Ag scaling: 4.56e-02
  GPP in Carbon scaling: 7.89e-01
  ```

**Uncertainty bands:** Default ±1σ (can adjust with `std_annual_multiplier`)

---

### Example 5: Ensemble Percent Difference

**From JSON file (Config #5):**
```json
{
    "output_files": [
        ["./control.dat", "./control_2.dat", "./control_3.dat", "./control_4.dat", "./control_5.dat"],
        ["./feedback.dat", "./feedback_2.dat", "./feedback_3.dat", "./feedback_4.dat", "./feedback_5.dat"],
        ["./ag_scaling.dat", "./ag_scaling_2.dat", "./ag_scaling_3.dat", "./ag_scaling_4.dat", "./ag_scaling_5.dat"],
        ["./carbon.dat", "./carbon_2.dat", "./carbon_3.dat", "./carbon_4.dat", "./carbon_5.dat"]
    ],
    "output_labels": ["Control", "Full feedback", "Ag scaling", "Carbon scaling"],
    "plot_directory": "./time_series_plots/ensemble_percent_difference",
    "plot_percent_difference": true,
    "y_limits": {"ZCO2": [-10, 10]}
}
```

**What This Creates:**

**Ensemble mean percent difference plots:**
- Control: Not shown (reference = 0%)
- Full feedback: Ensemble mean % change ± uncertainty
- Ag scaling: Ensemble mean % change ± uncertainty
- Carbon scaling: Ensemble mean % change ± uncertainty
- Markers where ensemble means significantly different (p < 0.05)

**Example values for ZCO2 in year 2100:**
- Full feedback: +3.2% ± 0.8% (shaded band from +2.4% to +4.0%)
- Ag scaling: +1.5% ± 0.5%
- Carbon scaling: +2.1% ± 0.6%

**Combines:** Ensemble uncertainty + percent difference + significance testing

---

### Example 6: Individual Members Despite Ensemble Format

**From JSON file (Config #6):**
```json
{
    "output_files": [
        ["./control.dat", "./control_2.dat", "./control_3.dat", "./control_4.dat", "./control_5.dat"],
        ["./feedback.dat", "./feedback_2.dat", "./feedback_3.dat", "./feedback_4.dat", "./feedback_5.dat"],
        ["./ag_scaling.dat", "./ag_scaling_2.dat", "./ag_scaling_3.dat", "./ag_scaling_4.dat", "./ag_scaling_5.dat"],
        ["./carbon.dat", "./carbon_2.dat", "./carbon_3.dat", "./carbon_4.dat", "./carbon_5.dat"]
    ],
    "output_labels": ["Control", "Full feedback", "Ag scaling", "Carbon scaling"],
    "plot_type": "individual",
    "plot_directory": "./time_series_plots/individual_override_ensemble_default"
}
```

**Key:** `"plot_type": "individual"` overrides default ensemble behavior

**What This Creates:**

**Individual member plots:**
- `time_series_GPP.pdf` - 20 lines total!
  - 5 Control lines (control, control_2, control_3, control_4, control_5)
  - 5 Full feedback lines
  - 5 Ag scaling lines
  - 5 Carbon scaling lines

**No averaging, no statistics, no uncertainty bands**

**Use case:** Want to see full range of ensemble spread without averaging

---

### Example 7: Simple Surfdata Comparison

**From JSON file (Config #7):**
```json
{
    "output_files": [
        "./control_time_series_surfdata_iESM_dyn_20240730.dat", 
        "./full_feedback_time_series_surfdata_iESM_dyn_20240730.dat"
    ],
    "output_labels": ["Control", "Full feedback"],
    "plot_directory": "./time_series_plots/surfdata/",
    "use_latex": true
}
```

**What This Creates:**

**Simple comparison plots:**
- All variables in surfdata files
- 2 lines per plot (Control, Full feedback)
- LaTeX formatting for publication quality
- Aggregated t-tests in p_values.dat

**Typical surfdata variables:**
- PCT_CROP, PCT_NATVEG, AREA_CROP, etc.
- Land use change variables

---

### Example 8: Version Comparison

**From JSON file (Config #9):**
```json
{
    "output_files": [
        "./control_surfdata_20240730.dat", 
        "./feedback_surfdata_20240730.dat",
        "./control_surfdata_20240820.dat", 
        "./feedback_surfdata_20240820.dat"
    ],
    "output_labels": [
        "Control (20240730)", 
        "Full feedback (20240730)", 
        "Control (20240820)", 
        "Full feedback (20240820)"
    ],
    "plot_directory": "./time_series_plots/surfdata_0730_vs_0820/"
}
```

**What This Creates:**

**Comparison across different data versions:**
- 4 lines per plot
- Control from two dates
- Full feedback from two dates
- Useful for tracking how data/code changes affect results

**Example interpretation:**
- If 20240730 and 20240820 lines identical → code changes didn't affect results
- If different → can quantify impact of updates

---

### Example 9: Surfdata Ensemble

**From JSON file (Config #10):**
```json
{
    "output_files": [
        ["./control_surfdata.dat", "./control_surfdata_2.dat", "./control_surfdata_3.dat", 
         "./control_surfdata_4.dat", "./control_surfdata_5.dat"],
        ["./feedback_surfdata.dat", "./feedback_surfdata_2.dat", "./feedback_surfdata_3.dat", 
         "./feedback_surfdata_4.dat", "./feedback_surfdata_5.dat"]
    ],
    "output_labels": ["Control", "Full feedback"],
    "plot_directory": "./time_series_plots/surfdata_ensemble/"
}
```

**What This Creates:**

**Ensemble analysis for surfdata:**
- 2 ensemble mean lines (Control, Full feedback)
- Shaded uncertainty bands (±1σ)
- Per-year significance markers
- Aggregated p-values in file

**Example for PCT_CROP:**
- Control ensemble mean: 15.2% ± 0.3% (2015) → 18.5% ± 0.5% (2100)
- Full feedback ensemble mean: 15.2% ± 0.3% (2015) → 19.8% ± 0.6% (2100)
- Markers appear starting around year 2050 where differences emerge

---

### Example 10: Surfdata Ensemble Percent Difference

**From JSON file (Config #11):**
```json
{
    "output_files": [
        ["./control_surfdata.dat", "./control_surfdata_2.dat", "./control_surfdata_3.dat", 
         "./control_surfdata_4.dat", "./control_surfdata_5.dat"],
        ["./feedback_surfdata.dat", "./feedback_surfdata_2.dat", "./feedback_surfdata_3.dat", 
         "./feedback_surfdata_4.dat", "./feedback_surfdata_5.dat"]
    ],
    "output_labels": ["Control", "Full feedback"],
    "plot_directory": "./time_series_plots/surfdata_ensemble_percent_difference/",
    "plot_percent_difference": true
}
```

**What This Creates:**

**Percent difference with ensemble uncertainty:**
- Control: Not shown (reference = 0%)
- Full feedback: Ensemble mean % change with uncertainty bands
- Markers at significant years

**Example for PCT_CROP in year 2100:**
- Full feedback: +7.0% ± 1.2% relative to Control
- Means: If Control = 18.5%, Full feedback = 19.8%
- Calculation: (19.8 - 18.5) / 18.5 × 100 = +7.0%
- Uncertainty from ensemble spread

---

## Output Files

### Plot Files

**Naming Convention:**
```
{plot_directory}/time_series_{variable}.pdf
{plot_directory}/time_series_{variable}_MAM.pdf  (seasonal)
{plot_directory}/time_series_{variable}_JJA.pdf
{plot_directory}/time_series_{variable}_SON.pdf
{plot_directory}/time_series_{variable}_DJF.pdf
{plot_directory}/time_series_{variable}_monthly.pdf
```

**Example Output Directory:**
```
./time_series_plots/ensemble/
├── time_series_GPP.pdf
├── time_series_GPP_JJA.pdf
├── time_series_GPP_monthly.pdf
├── time_series_NPP.pdf
├── time_series_NBP.pdf
├── time_series_ZCO2.pdf
├── time_series_ZCO2_MAM.pdf
├── time_series_ZCO2_JJA.pdf
├── time_series_ZCO2_SON.pdf
├── time_series_ZCO2_DJF.pdf
├── time_series_ZCO2_monthly.pdf
└── p_values.dat
```

### P-Value Files

**Purpose:** Records aggregated statistical test results

**Format:**
```
GPP in Full feedback: 1.234e-05
GPP in Ag scaling: 3.456e-03
NPP in Full feedback: 8.901e-07
NBP in Full feedback: 5.678e-02
```

**Location:** `{plot_directory}/p_values.dat`

**Interpretation:**
- < 0.05: Statistically significant difference from Control
- ≥ 0.05: Not statistically significant at 95% confidence level

**Note:** This file contains **aggregated** p-values (whole period). Per-year and per-month p-values shown as markers on plots.

---

## Best Practices

### 1. Choose Appropriate Plot Type

**Use `plot_type: "ensemble"` when:**
- Multiple ensemble members available
- Need statistical confidence
- Publishing results
- Want uncertainty quantification

**Use `plot_type: "individual"` when:**
- Exploring ensemble spread
- Diagnosing outliers
- Comparing specific members
- Initial investigation phase

---

### 2. Variable Selection

**Plot all variables:**
```json
"variables": "all"
```

**Plot specific variables:**
```json
"variables": ["GPP", "NPP", "NBP", "ER"]
```

**Recommendation:** Start specific for testing, expand to "all" for comprehensive analysis

---

### 3. Monthly Aggregation

**Use `"sum"` for:**
- Fluxes (GPP, NPP, NBP, ER)
- Precipitation (PRECIP)
- Any rate that should be integrated over time

**Use `"mean"` for:**
- State variables (temperature, CO₂ concentration)
- Fractions (land cover percentages)
- Dimensionless quantities

**Example:**
```json
"monthly_aggregation_type": {
    "GPP": "sum",
    "NPP": "sum",
    "PRECIP": "sum",
    "TREFHT": "mean",
    "ZCO2": "mean"
}
```

---

### 4. Ensemble Size

**Recommended:** 5-10 members

**Small (n=3-5):** Quick analysis, limited statistical power
**Medium (n=5-10):** Standard analysis, good balance (recommended)
**Large (n=10-20):** Publication-quality, robust statistics, high computational cost

---

### 5. Seasonal Analysis Strategy

**Step 1:** Include key seasons on main plot
```json
"include_seasons": {"GPP": {"JJA": true, "DJF": true}}
```

**Step 2:** Create separate plots for detailed analysis
```json
"seasons_to_plot_separately": {
    "GPP": {"MAM": true, "JJA": true, "SON": true, "DJF": true}
}
```

**Avoid:** Including all 4 seasons on main plot (too cluttered)

---

### 6. Ensemble File Organization

**Recommendation:** Use Option 1 (organized by scenario) for clarity

```json
"output_files": [
    ["ctrl.dat", "ctrl_2.dat", "ctrl_3.dat"],
    ["scen.dat", "scen_2.dat", "scen_3.dat"]
]
```

**Always provide `output_labels`** for unambiguous detection:
```json
"output_labels": ["Control", "Scenario"]
```

---

### 7. File Organization

```
project/
├── data/
│   ├── control_time_series.dat
│   ├── control_time_series_2.dat
│   ├── scenario_time_series.dat
│   └── scenario_time_series_2.dat
├── configs/
│   ├── ensemble_config.json
│   ├── individual_config.json
│   └── percent_diff_config.json
└── plots/
    ├── ensemble/
    ├── individual/
    ├── percent_diff/
    └── seasonal/
```

---

## Troubleshooting

### Issue 1: Wrong Number of Lines

**Symptom:** Plot has unexpected number of lines

**Possible Causes:**
1. `plot_type` set incorrectly
2. Ensemble format misdetected

**Solutions:**

**Check plot_type:**
```json
"plot_type": "ensemble"  // 1 line per scenario (mean)
"plot_type": "individual"  // N lines per scenario (all members)
```

**Verify ensemble detection:**
```json
"notify_output_files_transposed": true
```

**Check file count:**
```python
import json
with open('config.json') as f:
    config = json.load(f)[0]
    files = config['output_files']
    print(f"Structure: {len(files)} x {len(files[0])}")
```

---

### Issue 2: No Significance Markers

**Symptom:** Expected markers but none appear

**Possible Causes:**
1. Using `plot_type: "individual"` (no per-year tests)
2. All p-values > threshold
3. Ensemble size too small (n=1)

**Solutions:**

**Check plot type:**
```json
"plot_type": "ensemble"  // Required for per-year markers
```

**Lower threshold (exploratory only):**
```json
"p_value_threshold": 0.10  // More lenient
```

**Check p-value file:**
```bash
cat plots/p_values.dat
```
If aggregated p-values high, per-year likely also high

---

### Issue 3: File Not Found

**Error:**
```
FileNotFoundError: control_time_series.dat
```

**Solutions:**
```bash
# Verify files exist
ls -la control_time_series*.dat

# Use absolute paths
"output_files": ["/absolute/path/to/file.dat"]

# Check working directory
pwd
```

---

### Issue 4: Variable Not Found

**Error:**
```
KeyError: 'GPP'
```

**Solution:**
```python
import pandas as pd
df = pd.read_fwf('time_series.dat')
print(df.columns)
```

Check exact variable names (case-sensitive)

---

### Issue 5: Memory Error

**Error:**
```
MemoryError: Unable to allocate array
```

**Causes:** Too many variables × files processed simultaneously

**Solutions:**

**Process fewer variables:**
```json
"variables": ["GPP", "NPP"]  // Instead of "all"
```

**Run configurations separately:**
```bash
python script.py config1.json
python script.py config2.json
```

---

### Issue 6: LaTeX Errors

**Error:**
```
RuntimeError: Failed to process string with tex
```

**Solutions:**

**Disable LaTeX:**
```json
"use_latex": false
```

**Or install LaTeX:**
```bash
# Ubuntu/Debian
sudo apt-get install texlive texlive-latex-extra

# macOS
brew install mactex
```

---

## Advanced Features

### Per-Variable Configuration

Most parameters accept dictionaries for per-variable customization:

```json
{
    "variables": ["GPP", "NBP", "TREFHT"],
    "y_limits": {
        "GPP": [100, 180],
        "NBP": [-10, 20],
        "TREFHT": [285, 295]
    },
    "monthly_time_series_plot": {
        "GPP": true,
        "NBP": false,
        "TREFHT": true
    },
    "monthly_aggregation_type": {
        "GPP": "sum",
        "NBP": "sum",
        "TREFHT": "mean"
    }
}
```

---

### Multiple Configurations in One File

```json
[
    {
        "output_files": [...],
        "plot_directory": "./plots/set1/"
    },
    {
        "output_files": [...],
        "plot_directory": "./plots/set2/"
    }
]
```

**Result:** All plots generated in one execution

---

## References

- E3SM Documentation: [https://e3sm.org/](https://e3sm.org/)
- Matplotlib Documentation: [https://matplotlib.org/](https://matplotlib.org/)
- SciPy Stats: [https://docs.scipy.org/doc/scipy/reference/stats.html](https://docs.scipy.org/doc/scipy/reference/stats.html)
- GitHub Repository: [https://github.com/philipmyint/e3sm--gcam_analysis](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Contact

For questions or issues:
- **Philip Myint**: myint1@llnl.gov
- **Dalei Hao**: dalei.hao@pnnl.gov

---

## Version Information

**Script:** e3sm_plot_time_series.py  
**Documentation Version:** 2.0  
**Last Updated:** January 2026  
**Python Version:** 3.7+  
**Dependencies:** pandas, numpy, matplotlib, scipy

---

*This documentation provides comprehensive guidance for using the `e3sm_plot_time_series.py` script to create publication-quality time series plots from E3SM model output with ensemble analysis and statistical significance testing.*
