# GCAM Synthetic Ensemble Data

**Script:** `gcam_produce_synthetic_data.py`

## Overview

`gcam_produce_synthetic_data.py` creates synthetic ensemble versions of selected processed GCAM CSV files. It keeps the original scenarios and appends perturbed copies named with numeric suffixes such as `Control_2`, `Control_3`, and so on.

This script is intended for testing ensemble plotting and workflow examples. It is not a replacement for real ensemble simulations.

## Usage

```bash
python gcam_produce_synthetic_data.py
```

Limit parallel workers:

```bash
python gcam_produce_synthetic_data.py --max-processes 4
```

`--max-processes` overrides `MAX_PROCESSES`. If it is omitted, the script uses `MAX_PROCESSES`, which is all available CPUs minus one, with a minimum of one.

## Current Built-In Inputs

The script currently creates ensemble versions of:

```text
ag_commodity_prices_processed.csv
co2_emissions_regions_processed.csv
co2_emissions_sectors_processed.csv
land_allocation_processed.csv
land_allocation_processed_original_crop_names.csv
vegetation_scalars.csv
soil_scalars.csv
```

For each CSV input, the output is named by inserting `_ensemble` before `.csv`, for example:

```text
vegetation_scalars.csv -> vegetation_scalars_ensemble.csv
soil_scalars.csv -> soil_scalars_ensemble.csv
```

The previous combined scalar file is no longer used; vegetation and soil scalars are handled as separate standard `value` tables.

## Perturbation Method

For each scenario in each file, the script:

1. keeps the original rows,
2. creates synthetic copies,
3. multiplies selected numeric columns by a base multiplier between `1.02` and `1.05`,
4. adds a random perturbation between `-0.02` and `0.02`,
5. writes the combined original and synthetic rows to an `_ensemble` file.

The default number of total time series per scenario is five: the original plus four synthetic variants.

## Columns Modified

The current built-in files all perturb the standard `value` column:

```python
columns_to_modify = [['value'], ['value'], ['value'], ['value'], ['value'], ['value'], ['value']]
```

This matches the current scalar compiler, which writes `vegetation_scalars.csv` and `soil_scalars.csv` with a common `value` column.

## Editing The Built-In File List

This script does not read a JSON file. To synthesize different files, edit the `files`, `columns_to_modify`, and `scenario_labels` lists in the main execution block so all three lists have the same length.

Each input file must include the scenario column named by the corresponding `scenario_labels` entry.

## Workflow Position

Run this after:

1. `gcam_process_extracted_data.py` has created processed GCAM CSV files.
2. `gcam_compile_ehc_feedbacks.py` has created scalar outputs such as `vegetation_scalars.csv` and `soil_scalars.csv`.

Then run the plotting scripts against the generated `_ensemble.csv` files.
