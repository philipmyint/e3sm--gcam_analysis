# GCAM Extracted Data Processing

**Script:** `gcam_process_extracted_data.py`

## Overview

`gcam_process_extracted_data.py` converts raw CSV files extracted from GCAM XML or project-file queries into the processed CSV shape used by the plotting scripts.

It can:

- drop unneeded columns,
- lowercase column names,
- split composite columns,
- split columns at a named substring,
- aggregate duplicate key rows,
- add land-area metadata from a processed land allocation file,
- standardize GCAM crop names,
- write columns in a plotting-friendly order.

The previous standalone area-enrichment workflow is no longer used. Area columns are added here through `land_allocation_file`.

## Usage

```bash
python gcam_process_extracted_data.py path/to/config.json
```

Multiple JSON files can be supplied:

```bash
python gcam_process_extracted_data.py config1.json config2.json
```

Land allocation blocks are run first and serially because other blocks may depend on their output for `land_allocation_file`. Remaining blocks are processed in parallel.

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `input_file` | Yes | | Raw GCAM CSV file to process. |
| `output_file` | Yes | | Processed CSV file to write. |
| `columns_to_drop` | No | `None` | Columns to remove before processing. |
| `columns_to_split` | No | `None` | Dictionary mapping one source column to a list of new columns. Splits on underscores. |
| `columns_to_split_on_substring` | No | `None` | Dictionary describing columns split at the first occurrence of a substring. |
| `key_columns` | No | `None` | Columns used to sort, select, and aggregate rows. Names are lowercased internally. |
| `aggregate_function` | No | `"sum"` | Aggregation for duplicate `key_columns`: `"sum"` or `"area_weighted_mean"`. |
| `call_modify_crop_names` | No | `false` | Standardize original GCAM crop names to the common modified names. |
| `land_allocation_file` | Required for non-land-allocation files that need area metadata | `None` | Processed detailed land allocation file used to add `area` and `area_units`. |
| `max_processes` | No | `MAX_PROCESSES` | Maximum parallel workers. Defaults to all CPUs minus one, with a minimum of one. |

## Output Shape

When `key_columns` is supplied, the output is reordered so metadata keys come first and these columns come last when present:

```text
value, units, area, area_units
```

Most downstream plotting JSON assumes the numeric column is called `value`.

## Aggregation

`aggregate_function` controls what happens when multiple rows share the same `key_columns`.

| Value | Behavior | Typical use |
|-------|----------|-------------|
| `sum` | Sums numeric columns within each key group. | Land area, emissions, production, totals. |
| `area_weighted_mean` | Computes `sum(value * area) / sum(area)` and carries summed area. | Prices, scalars, rates, intensive quantities. |

String metadata columns such as `units` and `area_units` are preserved by taking the first value within each key group.

## Area Metadata

For non-land-allocation files, the script uses `land_allocation_file` to add `area` and `area_units`. This is required for later area-weighted plotting.

The land allocation file should usually be a processed detailed land allocation file that keeps original crop names when it will be used by `gcam_compile_ehc_feedbacks.py`, for example:

```json
"land_allocation_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_original_crop_names.csv"
```

## Crop Name Standardization

Set `call_modify_crop_names: true` to map original GCAM crop names to the common modified names used by most examples:

```text
biomass, biomassGrass, biomassTree -> BioenergyCrop
CornC4 -> Corn
OilPalmTree -> OilPalm
SugarCropC4 -> SugarCrop
```

When both area metadata and crop-name standardization are needed, area is added first and crop names are standardized afterward.

## Examples

### Process land allocation with standardized crop names

```json
{
    "input_file": "./../2025_DiVittorio_et_al_gcam/land_allocation.csv",
    "output_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed.csv",
    "columns_to_drop": ["Unnamed: 0", "Units"],
    "columns_to_split": {"landleaf": ["landtype", "basin"]},
    "key_columns": ["scenario", "region", "basin", "landtype", "year"],
    "aggregate_function": "sum",
    "call_modify_crop_names": true
}
```

### Process land allocation while preserving original crop names

```json
{
    "input_file": "./../2025_DiVittorio_et_al_gcam/land_allocation.csv",
    "output_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_original_crop_names.csv",
    "columns_to_drop": ["Unnamed: 0", "Units"],
    "columns_to_split": {"landleaf": ["landtype", "basin"]},
    "key_columns": ["scenario", "region", "basin", "landtype", "year"],
    "aggregate_function": "sum",
    "call_modify_crop_names": false
}
```

### Process commodity prices with area-weighted aggregation

```json
{
    "input_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices.csv",
    "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
    "columns_to_drop": ["Unnamed: 0", "Units"],
    "key_columns": ["scenario", "region", "sector", "year"],
    "aggregate_function": "area_weighted_mean",
    "call_modify_crop_names": true,
    "land_allocation_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_original_crop_names.csv"
}
```

## Troubleshooting

If `area_weighted_mean` fails with a missing `area` column, make sure the JSON block includes `land_allocation_file` and that the referenced processed land allocation file exists.

If crop-name standardization produces unexpected groups, check whether the input file uses original or already modified crop names. Use the original crop-name land allocation file when merging area into EHC scalar outputs.
