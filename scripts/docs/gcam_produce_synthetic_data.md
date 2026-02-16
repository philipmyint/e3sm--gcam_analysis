# GCAM Synthetic Data Generation Script Documentation

## Overview

**Script Name:** `gcam_produce_synthetic_data.py`

**Purpose:** Generates synthetic ensemble sets of time series data by applying random perturbations to base GCAM (Global Change Analysis Model) simulation outputs. This script is designed for testing and validation purposes, allowing users to create multiple variations of existing time series data to simulate ensemble runs.

**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Script Description

The `gcam_produce_synthetic_data.py` script reads base time series data from CSV or fixed-width format (.dat) files and creates synthetic variations by applying controlled random perturbations. Each synthetic time series maintains the general pattern of the original data while introducing realistic variability through multiplicative factors.

### Key Features

- **Parallel Processing:** Utilizes Python's multiprocessing capabilities to generate synthetic data for multiple files simultaneously
- **Flexible Input Formats:** Supports both CSV and fixed-width format (.dat) files
- **Controlled Perturbations:** Applies base multipliers combined with random variations to maintain data consistency
- **Scenario Management:** Processes multiple scenarios within each file independently
- **Automatic Output Naming:** Generates ensemble files with intuitive naming conventions

### Workflow

1. Reads base time series data from input files
2. For each scenario in the file:
   - Creates multiple synthetic variations
   - Applies base multipliers (ranging from 1.02 to 1.05)
   - Adds random perturbations (±2%)
   - Multiplies specified data columns by these combined factors
3. Concatenates all variations with the original data
4. Outputs results to a new ensemble file

---

## Function Documentation

### `produce_synthetic_time_series(inputs)`

**Purpose:** Core function that generates synthetic ensemble time series from a base dataset.

**Parameters:**
- `inputs` (list): A list containing five elements:
  1. `file` (str): Path to the input file containing base time series data
  2. `scenarios` (list): List of scenario names present in the file
  3. `scenario_label` (str): Column name identifying scenarios in the dataframe
  4. `columns_to_modify` (list): Column names containing numerical data to perturb
  5. `num_synthetic_sets_in_ensemble` (int): Total number of time series including the base

**Returns:** None (writes output to file)

**Process:**
1. Loads data from the specified file
2. Creates a copy for the ensemble dataset
3. Generates base multipliers linearly spaced between 1.02 and 1.05
4. For each scenario and each synthetic variation:
   - Applies base multiplier plus random variation (±2%)
   - Creates new scenario label (e.g., `scenario_2`, `scenario_3`, etc.)
   - Multiplies specified columns by the combined multiplier
5. Outputs to CSV or fixed-width format file with `_ensemble` suffix

**Output Naming Convention:**
- CSV input: `original_name.csv` → `original_name_ensemble.csv`
- DAT input: `original_name.dat` → `original_name_ensemble.dat`
- Other formats: `original_name` → `original_name_ensemble`

---

## Main Execution Block

The script's main execution block (runs when `if __name__ == '__main__':`) processes multiple GCAM output files in parallel.

### Default Configuration

The script is pre-configured to process six different GCAM output files from the DiVittorio et al. 2025 study:

1. **Agricultural commodity prices** (`ag_commodity_prices_processed.csv`)
2. **CO2 emissions by region** (`co2_emissions_regions_processed.csv`)
3. **CO2 emissions by sector** (`co2_emissions_sectors_processed.csv`)
4. **Land allocation** (`land_allocation_processed.csv`)
5. **Land allocation with original crop names** (`land_allocation_processed_original_crop_names.csv`)
6. **Vegetation and soil scalars** (`scalars_control+full_feedback.csv`)

---

## Input Parameters Table

| Parameter | Default Value | Possible Values | Required? | Description |
|-----------|--------------|-----------------|-----------|-------------|
| `files` | List of 6 CSV file paths | Any valid file path(s) to CSV or .dat files | Yes | List of input files containing base time series data |
| `columns_to_modify` | `[['value'], ['value'], ['value'], ['value'], ['value'], ['vegetation', 'soil']]` | Any valid column name(s) from the input files | Yes | Columns containing numerical data to be perturbed for each file |
| `scenario_labels` | `['scenario', 'scenario', 'scenario', 'scenario', 'scenario', 'scenario']` | Any valid column name identifying scenarios | Yes | Column name that identifies different scenarios in each file |
| `num_variations_for_each_scenario` | `[5, 5, 5, 5, 5, 5]` | Any positive integer ≥ 1 | Yes | Total number of time series (including base) for each file |
| `base_multipliers` | `np.linspace(1.02, 1.05, num_synthetic_sets_in_ensemble)` | Any array of positive float values | No | Base multiplicative factors applied to data (hardcoded in function) |
| `random_multipliers` | `np.random.uniform(low=-0.02, high=0.02, size=len(df_this_scenario))` | Any range of float values | No | Random perturbations added to base multipliers (hardcoded in function) |
| `processes` | `MAX_PROCESSES` | Any positive integer | No | Number of parallel processes to use |

---

## Detailed Parameter Descriptions

### Files Configuration

**Default Files:**
```python
files = [
    "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
    "./../2025_DiVittorio_et_al_gcam/co2_emissions_regions_processed.csv",
    "./../2025_DiVittorio_et_al_gcam/co2_emissions_sectors_processed.csv",
    "./../2025_DiVittorio_et_al_gcam/land_allocation_processed.csv",
    "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_original_crop_names.csv",
    "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback.csv"
]
```

**Customization:** Users can modify this list to include any CSV or .dat files containing GCAM time series data. All files must contain a column identifying scenarios.

### Columns to Modify

**Default Configuration:**
```python
columns_to_modify = [
    ['value'],      # For agricultural commodity prices
    ['value'],      # For regional CO2 emissions
    ['value'],      # For sectoral CO2 emissions
    ['value'],      # For land allocation
    ['value'],      # For land allocation (original crop names)
    ['vegetation', 'soil']  # For vegetation and soil scalars
]
```

**Customization:** Each sub-list corresponds to one file in the `files` list. Users must specify which columns contain numerical data that should be perturbed. Column names must exactly match those in the input files.

### Scenario Labels

**Default Configuration:**
```python
scenario_labels = ['scenario', 'scenario', 'scenario', 'scenario', 'scenario', 'scenario']
```

**Customization:** Each element corresponds to one file. The value should be the exact column name that identifies different scenarios in that file (e.g., 'scenario', 'case_name', 'simulation_id').

### Number of Variations

**Default Configuration:**
```python
num_variations_for_each_scenario = [5, 5, 5, 5, 5, 5]
```

**Interpretation:** For each file, this creates 5 total time series (1 original + 4 synthetic variations). The synthetic variations are numbered 2, 3, 4, 5, etc.

**Customization:** Users can specify different numbers for each file. For example, `[10, 5, 3, 5, 5, 8]` would create:
- 10 total variations for file 1 (1 base + 9 synthetic)
- 5 total variations for file 2 (1 base + 4 synthetic)
- 3 total variations for file 3 (1 base + 2 synthetic)
- And so on...

---

## Perturbation Methodology

### Base Multipliers

The script generates base multipliers using:
```python
base_multipliers = np.linspace(1.02, 1.05, num_synthetic_sets_in_ensemble)
```

**Effect:** Creates systematic variation with multipliers ranging from 1.02 (2% increase) to 1.05 (5% increase), evenly distributed across the number of synthetic sets.

**Example:** For 4 synthetic variations:
- Variation 1: base multiplier = 1.02
- Variation 2: base multiplier = 1.03
- Variation 3: base multiplier = 1.04
- Variation 4: base multiplier = 1.05

### Random Perturbations

Additional random variation is added:
```python
random_multipliers = np.random.uniform(low=-0.02, high=0.02, size=len(df_this_scenario))
```

**Effect:** Adds ±2% random variation to each data point, creating realistic noise and variability within each synthetic time series.

### Combined Effect

Final multiplier for each data point:
```
final_multiplier = base_multiplier + random_perturbation
```

**Range:** Approximately 1.00 to 1.07 (with variation between -2% and +7% relative to original values)

---

## Usage Instructions

### Running the Script

Since this script is designed for testing purposes and does not use JSON configuration files, users must modify the script directly to customize behavior.

**Command Line Execution:**
```bash
python gcam_produce_synthetic_data.py
```

**Note:** Unlike other scripts in the repository, this does not accept JSON file arguments.

### Modifying Input Files

To process different files, edit the `files` list in the main execution block:

```python
files = [
    "path/to/your/file1.csv",
    "path/to/your/file2.dat",
    # Add more files as needed
]
```

### Adjusting Perturbation Parameters

To modify the perturbation behavior, users can edit the function directly:

**Change base multiplier range:**
```python
# Original:
base_multipliers = np.linspace(1.02, 1.05, num_synthetic_sets_in_ensemble)

# Example modification (1% to 3% increase):
base_multipliers = np.linspace(1.01, 1.03, num_synthetic_sets_in_ensemble)
```

**Change random perturbation range:**
```python
# Original:
random_multipliers = np.random.uniform(low=-0.02, high=0.02, size=len(df_this_scenario))

# Example modification (±5% random variation):
random_multipliers = np.random.uniform(low=-0.05, high=0.05, size=len(df_this_scenario))
```

### Adjusting Number of Variations

Modify the `num_variations_for_each_scenario` list:

```python
# Original (5 variations for each file):
num_variations_for_each_scenario = [5]*len(files)

# Custom (different numbers for each file):
num_variations_for_each_scenario = [10, 8, 6, 5, 5, 5]

# Same number for all files:
num_variations_for_each_scenario = [15]*len(files)
```

---

## Output Files

### Output Format

The script preserves the input file format:
- **CSV input** → CSV output with comma separation
- **DAT input** → Fixed-width format output
- **Other formats** → Fixed-width format output

### Output Structure

Each output file contains:
1. **Original data:** All rows from the base time series
2. **Synthetic variations:** New rows for each synthetic scenario

**Example scenario naming:**
- Original scenario: `control`
- Synthetic variations: `control_2`, `control_3`, `control_4`, `control_5`

### File Size Considerations

Output files will be approximately N times larger than input files, where N is the number of variations specified in `num_variations_for_each_scenario`.

**Example:** If the input file has 1,000 rows and 5 variations are requested, the output will have approximately 5,000 rows (1,000 original + 4,000 synthetic).

---

## Performance Considerations

### Multiprocessing

The script uses parallel processing with limited cores to reduce memory pressure:
```python
with multiprocessing.Pool(processes=MAX_PROCESSES) as pool:
    pool.map(produce_synthetic_time_series, inputs)
```

**Advantages:**
- Processes multiple files simultaneously
- Significantly reduces execution time for large ensembles
- Scales well with available hardware

**Timing Information:**
The script prints execution time for each file and total execution time:
```
Elapsed time for producing the synthetic data for file1.csv: 12.45 seconds
Elapsed time for producing the synthetic data for file2.csv: 8.32 seconds
...
Elapsed time for producing all the synthetic data: 45.67 seconds
```

### Memory Requirements

Memory usage depends on:
- Size of input files
- Number of scenarios per file
- Number of variations requested

**Estimation:** Peak memory ≈ (file size) × (num_variations) × (number of simultaneous processes)

---

## Dependencies

### Required Python Modules

```python
import multiprocessing  # Standard library
import numpy as np      # Install: pip install numpy
import pandas as pd     # Install: pip install pandas
import time            # Standard library
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_fwf  # Custom module
```

### Custom Utility Functions

The script requires two custom functions from `utility_dataframes.py`:

1. **`read_file_into_dataframe(file)`**
   - Reads CSV or fixed-width format files into pandas DataFrame
   - Automatically detects file format

2. **`write_dataframe_to_fwf(new_file, df_ensemble)`**
   - Writes DataFrame to fixed-width format file
   - Maintains proper column alignment

**Note:** These utility functions must be available in the Python path.

---

## Important Notes and Warnings

### 1. Testing Purpose Only

This script is designed for **testing and validation purposes only**. It should not be used when actual ensemble simulation data are available, as it creates synthetic data rather than processing real GCAM outputs.

### 2. Data Integrity

The perturbations applied are multiplicative, meaning:
- Relative relationships between data points are partially preserved
- Negative values will become more negative (multiplied by factors >1)
- Zero values remain zero

### 3. Scenario Numbering

Synthetic scenarios start from `_2` to avoid confusion with the base time series:
- Base scenario: `scenario_name`
- Synthetic variations: `scenario_name_2`, `scenario_name_3`, etc.

### 4. No JSON Configuration

Unlike most other scripts in the repository, this script **does not use JSON configuration files**. All modifications must be made directly in the Python script.

### 5. File Overwriting

The script creates new files with `_ensemble` suffix. Original files are not modified. However, if ensemble files already exist, they will be overwritten without warning.

### 6. Random Seed

The script does not set a random seed, meaning each execution will produce different random perturbations. For reproducible results, add at the beginning of the script:
```python
np.random.seed(42)  # Or any integer
```

---

## Example Workflow

### Typical Usage Scenario

1. **Prepare base time series files** from GCAM extraction and processing scripts
2. **Edit the script** to specify desired files and parameters
3. **Run the script** to generate synthetic ensemble
4. **Use ensemble files** for testing plotting and analysis scripts

### Example Modification

```python
# Customize for 3 files with 10 variations each
files = [
    "my_emissions_data.csv",
    "my_land_use_data.csv",
    "my_prices_data.csv"
]

columns_to_modify = [
    ['emissions'],
    ['area', 'yield'],
    ['price']
]

scenario_labels = [
    'scenario',
    'case_id',
    'run_name'
]

num_variations_for_each_scenario = [10, 10, 10]
```

---

## Troubleshooting

### Common Issues

**Issue:** `FileNotFoundError`
- **Cause:** Incorrect file path
- **Solution:** Verify file paths are correct and files exist

**Issue:** `KeyError` when modifying columns
- **Cause:** Column name doesn't exist in the file
- **Solution:** Check column names match exactly (case-sensitive)

**Issue:** Memory error
- **Cause:** Too many variations or files too large
- **Solution:** Reduce number of variations or process fewer files at once

**Issue:** Slow execution
- **Cause:** Large files or limited available cores
- **Solution:** Reduce number of variations or process files sequentially

---

## Related Scripts

This script is part of the E3SM-GCAM analysis toolkit. Related scripts include:

- **`gcam_extract_csv_from_project_files.R`**: Extracts data from GCAM project files
- **`gcam_process_extracted_data.py`**: Processes extracted CSV files
- **`gcam_plot_time_series.py`**: Plots time series from GCAM data
- **`e3sm_produce_synthetic_spatial_data.py`**: E3SM equivalent for spatial data
- **`e3sm_produce_synthetic_time_series.py`**: E3SM equivalent for time series

---

## References

- DiVittorio et al. (2025). "E3SM-GCAM coupling methodology and applications." *Journal of Advances in Modeling Earth Systems*. [DOI: 10.1029/2024MS004806](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024MS004806)
- GitHub Repository: [https://github.com/philipmyint/e3sm--gcam_analysis](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Contact

For questions or issues:
- **Philip Myint**: myint1@llnl.gov
- **Dalei Hao**: dalei.hao@pnnl.gov

---

## Version History

- **Current Version:** As documented in repository (January 2025)
- **Last Updated:** January 2026

---

*This documentation is intended to help users understand and modify the `gcam_produce_synthetic_data.py` script for creating synthetic GCAM ensemble datasets for testing purposes.*
