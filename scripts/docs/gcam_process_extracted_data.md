# GCAM Data Processing Script Documentation


**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Overview

The `gcam_process_extracted_data.py` script is a Python tool for processing raw data extracted from GCAM (Global Change Analysis Model) project files. It transforms, cleans, aggregates, and standardizes GCAM output data, preparing it for analysis and visualization. The script is particularly useful for handling land allocation data, commodity prices, emissions data, and other GCAM outputs.

## Table of Contents

1. [Installation & Requirements](#installation--requirements)
2. [Basic Usage](#basic-usage)
3. [Configuration File Structure](#configuration-file-structure)
4. [Configuration Parameters](#configuration-parameters)
5. [Data Processing Operations](#data-processing-operations)
6. [Examples](#examples)
7. [Advanced Features](#advanced-features)
8. [Understanding GCAM Land Types](#understanding-gcam-land-types)
9. [Troubleshooting](#troubleshooting)

---

## Installation & Requirements

### Required Python Libraries

```bash
pip install pandas numpy
```

### Required Utility Modules

The script imports several utility modules that must be in the same directory or Python path:
- `utility_dataframes` - Functions for reading/writing data files
- `utility_gcam` - GCAM-specific functions for crop name standardization and land type grouping

---

## Basic Usage

### Command Line Execution

```bash
python gcam_process_extracted_data.py path/to/config.json
```

You can also specify multiple JSON configuration files:

```bash
python gcam_process_extracted_data.py config1.json config2.json config3.json
```

The script will process all files in parallel using all available CPU cores for maximum efficiency.

### What the Script Does

The script performs the following operations on each input file:

1. **Reads** the raw CSV data extracted from GCAM
2. **Drops** unnecessary columns (e.g., index columns, unit columns)
3. **Splits** composite columns into separate columns
4. **Sorts** data by year and other key columns
5. **Aggregates** data when multiple rows share the same keys
6. **Standardizes** crop names to a common naming convention
7. **Writes** the processed data to an output file

---

## Configuration File Structure

The configuration file is a JSON array where each object specifies the processing parameters for one input file. Here's a minimal example:

```json
[
    {
        "input_file": "./data/land_allocation.csv",
        "output_file": "./data/land_allocation_processed.csv",
        "key_columns": ["scenario", "region", "basin", "landtype", "year"]
    }
]
```

---

## Configuration Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `input_file` | string | Path to the raw CSV file extracted from GCAM |
| `output_file` | string | Path where the processed CSV file will be saved |

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `columns_to_drop` | list | `None` | Column names to remove from the dataset |
| `columns_to_split` | dictionary | `None` | Columns to split (key: column name, value: new column names) |
| `key_columns` | list | `None` | Columns that uniquely identify rows for aggregation |
| `mean_or_sum_if_more_than_one_row_in_same_landtype_group` | string | `None` | Aggregation method: `"mean"` or `"sum"` |
| `call_modify_crop_names` | boolean | `false` | Whether to modify/standardize crop names |

---

## Data Processing Operations

### 1. Column Dropping

Remove unnecessary columns that were included in the GCAM extraction but aren't needed for analysis.

**Example:**
```json
{
    "columns_to_drop": ["Unnamed: 0", "Units"]
}
```

**Common columns to drop:**
- `Unnamed: 0` - Auto-generated index column
- `Units` - Unit information (often redundant if units are consistent)

### 2. Column Splitting

Split composite columns (where multiple pieces of information are concatenated with underscores) into separate columns.

**Example:**
```json
{
    "columns_to_split": {
        "landleaf": ["landtype", "basin"]
    }
}
```

**What this does:**
- Takes a column like `landleaf` with values like `"Forest_Amazon"`
- Splits it into `landtype` = `"Forest"` and `basin` = `"Amazon"`

### 3. Column Sorting and Selection

Automatically sorts data and selects relevant columns based on `key_columns`.

**Example:**
```json
{
    "key_columns": ["scenario", "region", "basin", "landtype", "year"]
}
```

**Result:**
- Data is sorted by these columns in order
- Only these columns plus `value` are retained in output

### 4. Data Aggregation

When multiple rows have the same key column values, they can be aggregated.

**Aggregation Methods:**

| Method | Description | Use Case |
|--------|-------------|----------|
| `"sum"` | Add all values together | Land areas, total production |
| `"mean"` | Calculate average | Prices, scalars, ratios |
| `None` | No aggregation | When each row is already unique |

**Example:**
```json
{
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum"
}
```

### 5. Crop Name Standardization

GCAM uses various crop naming conventions. This feature standardizes them to a common set.

**Example:**
```json
{
    "call_modify_crop_names": true
}
```

**Crop modification/standardization examples:**
- `biomass`, `biomassGrass`, `biomassTree` → `BioenergyCrop`
- `CornC4` → `Corn`
- `OilPalmTree` → `OilPalm`
- `SugarCropC4` → `SugarCrop`

**Complete mapping:**
```
Original Name       →  Standard Name
─────────────────────────────────────
biomass             →  BioenergyCrop
biomassGrass        →  BioenergyCrop
biomassTree         →  BioenergyCrop
CornC4              →  Corn
FodderHerbC4        →  FodderHerb
FruitsTree          →  Fruits
MiscCropC4          →  MiscCrop
MiscCropTree        →  MiscCrop
NutsSeedsTree       →  NutsSeeds
OilCropTree         →  OilCrop
OilPalmTree         →  OilPalm
OtherGrainC4        →  OtherGrain
SugarCropC4         →  SugarCrop
```

---

## Examples

### Example 1: Processing Land Allocation Data

**Purpose:** Convert raw GCAM land allocation data into a clean format with standardized crop names.

```json
{
    "input_file": "./data/land_allocation.csv",
    "output_file": "./data/land_allocation_processed.csv",
    "columns_to_drop": ["Unnamed: 0", "Units"],
    "columns_to_split": {
        "landleaf": ["landtype", "basin"]
    },
    "key_columns": ["scenario", "region", "basin", "landtype", "year"],
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum",
    "call_modify_crop_names": true
}
```

**What happens:**
1. Drops the index column and units column
2. Splits `landleaf` (e.g., `"Corn_Amazon"`) into `landtype` and `basin`
3. Sorts by scenario, region, basin, landtype, and year
4. Sums areas when multiple crops map to the same standardized name
5. Standardizes crop names (e.g., `CornC4` → `Corn`)
6. Saves processed data

**Input data structure:**
```
scenario | region | landleaf        | year | value | Units | Unnamed: 0
---------|--------|-----------------|------|-------|-------|------------
Control  | USA    | Corn_Amazon     | 2020 | 150.5 | km²   | 0
Control  | USA    | CornC4_Amazon   | 2020 | 75.3  | km²   | 1
```

**Output data structure:**
```
scenario | region | basin  | landtype | year | value
---------|--------|--------|----------|------|-------
Control  | USA    | Amazon | Corn     | 2020 | 225.8
```

### Example 2: Processing Land Allocation with Original Crop Names

**Purpose:** Process the same data but keep original GCAM crop names.

```json
{
    "input_file": "./data/land_allocation.csv",
    "output_file": "./data/land_allocation_processed_original_crop_names.csv",
    "columns_to_drop": ["Unnamed: 0", "Units"],
    "columns_to_split": {
        "landleaf": ["landtype", "basin"]
    },
    "key_columns": ["scenario", "region", "basin", "landtype", "year"],
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum",
    "call_modify_crop_names": false
}
```

**Difference from Example 1:**
- `call_modify_crop_names` is `false`
- Crop names remain as they appear in GCAM output
- Useful when you need to track specific crop variants (e.g., C3 vs C4)

### Example 3: Processing Agricultural Commodity Prices

**Purpose:** Clean and aggregate commodity price data.

```json
{
    "input_file": "./data/ag_commodity_prices.csv",
    "output_file": "./data/ag_commodity_prices_processed.csv",
    "columns_to_drop": ["Unnamed: 0", "Units"],
    "key_columns": ["scenario", "region", "sector", "year"],
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "mean",
    "call_modify_crop_names": true
}
```

**Key differences:**
- No column splitting (price data doesn't need it)
- Uses `"mean"` aggregation (prices should be averaged, not summed)
- `sector` instead of `landtype` and `basin`
- Standardizes crop names for consistency with other analyses

**Why use "mean" for prices:**
- If multiple price quotes exist for the same commodity/region/year, averaging gives a representative price
- Summing prices would give meaningless values

### Example 4: Processing CO₂ Emissions by Sector

**Purpose:** Process emissions data with no crop name standardization needed.

```json
{
    "input_file": "./data/co2_emissions_sectors.csv",
    "output_file": "./data/co2_emissions_sectors_processed.csv",
    "columns_to_drop": ["Unnamed: 0", "Units"],
    "key_columns": ["scenario", "region", "sector", "year"]
}
```

**What happens:**
1. Drops unnecessary columns
2. Sorts by scenario, region, sector, and year
3. No aggregation (each row should already be unique)
4. No crop name modification (not applicable to emissions data)

### Example 5: Processing CO₂ Emissions by Region

**Purpose:** Process regional aggregate emissions data.

```json
{
    "input_file": "./data/co2_emissions_regions.csv",
    "output_file": "./data/co2_emissions_regions_processed.csv",
    "columns_to_drop": ["Unnamed: 0", "Units"],
    "key_columns": ["scenario", "region", "year"]
}
```

**Simplest case:**
- Minimal processing needed
- No sector breakdown (regional totals only)
- Straightforward cleaning and sorting

---

## Advanced Features

### Understanding Aggregation with Crop Name Standardization

When you enable crop name standardization, multiple original crop types may map to the same standardized name. The script handles this intelligently:

**Example scenario:**
```
Original Data:
scenario | region | landtype      | year | value
---------|--------|---------------|------|-------
Control  | USA    | biomass       | 2030 | 50
Control  | USA    | biomassGrass  | 2030 | 30
Control  | USA    | biomassTree   | 2030 | 20
```

**With `"mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum"`:**
```
Processed Data:
scenario | region | landtype       | year | value
---------|--------|----------------|------|-------
Control  | USA    | BioenergyCrop  | 2030 | 100
```

**With `"mean_or_sum_if_more_than_one_row_in_same_landtype_group": "mean"`:**
```
Processed Data:
scenario | region | landtype       | year | value
---------|--------|----------------|------|-------
Control  | USA    | BioenergyCrop  | 2030 | 33.33
```

### Column Name Normalization

All column names are automatically converted to lowercase. This ensures consistency:

**Input columns:** `Scenario`, `Region`, `YEAR`, `Value`  
**Output columns:** `scenario`, `region`, `year`, `value`

### Automatic Sorting

Data is always sorted by year first, then by key columns if specified. This ensures:
- Chronological ordering for time series analysis
- Consistent row ordering across files
- Easier debugging and manual inspection

---

## Understanding GCAM Land Types

### Standard Land Type Groups

GCAM organizes land into several major categories. The `utility_gcam.py` module defines these groups:

| Group | Land Types Included |
|-------|---------------------|
| **forest** | Forest, ProtectedUnmanagedForest, UnmanagedForest |
| **pasture** | Pasture, UnmanagedPasture, ProtectedUnmanagedPasture |
| **grass** | Grassland, ProtectedGrassland |
| **crop** | 18 standardized crop types (see below) |
| **shrub** | Shrubland, ProtectedShrubland |
| **urban** | UrbanLand |
| **other** | RockIceDesert, Tundra |

### Standardized Crop Types

After name standardization, crops are consolidated into 18 types:

1. BioenergyCrop
2. Corn
3. FiberCrop
4. FodderGrass
5. FodderHerb
6. Fruits
7. Legumes
8. MiscCrop
9. NutsSeeds
10. OilCrop
11. OilPalm
12. OtherGrain
13. Rice
14. RootTuber
15. Soybean
16. SugarCrop
17. Vegetables
18. Wheat

### Original Crop Types (Before Standardization)

GCAM's original crop list includes 29 types, many with C3/C4 variants or tree/grass variants:

- biomass, biomassGrass, biomassTree
- Corn, CornC4
- FiberCrop
- FodderGrass, FodderHerb, FodderHerbC4
- Fruits, FruitsTree
- Legumes
- MiscCrop, MiscCropC4, MiscCropTree
- NutsSeeds, NutsSeedsTree
- OilCrop, OilCropTree, OilPalmTree
- OtherGrain, OtherGrainC4
- Rice
- RootTuber
- Soybean
- SugarCrop, SugarCropC4
- Vegetables
- Wheat

---

## Troubleshooting

### Common Issues

**Problem:** `FileNotFoundError` when running the script
```
Solution: Check that the input_file path is correct. Use absolute paths 
if running from a different directory.
```

**Problem:** Missing columns after processing
```
Solution: Verify that key_columns exactly match the column names in your 
input file (case-insensitive after processing). Check for typos.
```

**Problem:** Unexpected aggregation results
```
Solution: Check your choice of "mean" vs "sum". Prices should use "mean", 
areas should use "sum". If you're getting zeros or very large values, you 
may have chosen the wrong method.
```

**Problem:** Crop names not being standardized
```
Solution: Ensure call_modify_crop_names is set to true (not "true" as a 
string). Check that your data actually contains crop types that need 
standardization.
```

**Problem:** Output file is empty or has very few rows
```
Solution: Check if your columns_to_drop list accidentally includes critical 
columns. Verify that key_columns doesn't over-specify uniqueness, causing 
rows to be dropped.
```

### Input File Requirements

Your CSV file should have:
- A `year` column for temporal data
- A `value` column containing the numeric data
- Other columns for categorization (scenario, region, sector, landtype, etc.)
- Consistent naming conventions

**Example of well-formed input:**
```csv
scenario,region,landleaf,year,value,Units,Unnamed: 0
Control,USA,Corn_Amazon,2020,150.5,km²,0
Control,USA,Wheat_Mississippi,2020,200.3,km²,1
```

### Performance Tips

1. **Multiple Files:** Process multiple files at once by listing them all in one JSON config
2. **Parallel Processing:** The script automatically uses all CPU cores
3. **Large Files:** For very large files (>1GB), ensure sufficient RAM is available
4. **Debugging:** Process one file at a time first to verify configuration before batch processing

---

## Best Practices

### 1. Keep Original Data

Always preserve original GCAM extraction files. Create processed versions with different names:

```json
{
    "input_file": "./data/raw/land_allocation.csv",
    "output_file": "./data/processed/land_allocation_processed.csv"
}
```

### 2. Maintain Two Versions When Needed

For crop data, consider creating both standardized and original versions:

```json
[
    {
        "input_file": "./data/land_allocation.csv",
        "output_file": "./data/land_allocation_standardized.csv",
        "call_modify_crop_names": true
    },
    {
        "input_file": "./data/land_allocation.csv",
        "output_file": "./data/land_allocation_original.csv",
        "call_modify_crop_names": false
    }
]
```

### 3. Document Your Aggregation Choices

Keep notes about why you chose "mean" vs "sum" for different datasets:

- **Sum:** Land areas, production quantities, emissions totals
- **Mean:** Prices, rates, scalars, percentages

### 4. Use Consistent Key Columns

Maintain consistent key column ordering across related datasets for easier merging later:

```json
"key_columns": ["scenario", "region", "basin", "landtype", "year"]
```

### 5. Verify Output

After processing, check the output file:
- Row counts make sense
- Value ranges are reasonable
- No unexpected NaN or zero values
- Column names are as expected

---

## Integration with Other GCAM Scripts

This data processing script works in conjunction with other GCAM analysis tools:

### Workflow Example

```
1. Extract data from GCAM → raw CSV files
2. Process with gcam_process_extracted_data.py → cleaned CSV files
3. Visualize with gcam_plot_box_and_whiskers.py → box plots
4. Further analysis with time series scripts → trend plots
```

### Data Flow

```
GCAM Model Output
       ↓
[Raw CSV Files]
       ↓
gcam_process_extracted_data.py
       ↓
[Processed CSV Files]
       ↓
┌──────┴──────┬──────────────┬────────────┐
│             │              │            │
Box Plots   Time Series   Statistics   Custom
                                        Analysis
```

---

## Output File Characteristics

### Structure

Processed files have a consistent structure:
- **Columns:** Only key columns plus `value`
- **Sorting:** By year, then by other key columns
- **Format:** CSV with headers
- **Encoding:** UTF-8

### Example Output

```csv
scenario,region,basin,landtype,year,value
Control,USA,Amazon,BioenergyCrop,2015,45.2
Control,USA,Amazon,BioenergyCrop,2020,52.1
Control,USA,Amazon,BioenergyCrop,2025,58.7
Control,USA,Amazon,Corn,2015,145.8
Control,USA,Amazon,Corn,2020,152.3
```

---

## Execution Time

Processing times vary based on:
- **File size:** 10MB file ≈ 1-2 seconds
- **Number of files:** Parallel processing scales well
- **Operations:** Column splitting and crop name standardization add minimal overhead
- **System:** More CPU cores = faster parallel processing

**Example output:**
```
Elapsed time for producing land_allocation_processed.csv: 1.23 seconds
Elapsed time for producing ag_commodity_prices_processed.csv: 0.87 seconds
Elapsed time for producing all files: 1.45 seconds
```

---

## Version Information

**Script:** gcam_process_extracted_data.py  
**Documentation Version:** 1.0  
**Last Updated:** January 2026  
**Python Version:** 3.7+

---

## Appendix: Complete Configuration Reference

### Full Parameter List

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `input_file` | string | ✓ | - | Path to input CSV file |
| `output_file` | string | ✓ | - | Path to output CSV file |
| `columns_to_drop` | list | ✗ | `None` | Columns to remove |
| `columns_to_split` | dict | ✗ | `None` | Columns to split into new columns |
| `key_columns` | list | ✗ | `None` | Columns for sorting and grouping |
| `mean_or_sum_if_more_than_one_row_in_same_landtype_group` | string | ✗ | `None` | Aggregation method: `"mean"` or `"sum"` |
| `call_modify_crop_names` | boolean | ✗ | `false` | Enable crop name standardization |

### Template Configuration

```json
[
    {
        "input_file": "./path/to/input.csv",
        "output_file": "./path/to/output.csv",
        "columns_to_drop": ["Unnamed: 0", "Units"],
        "columns_to_split": {
            "composite_column": ["new_col1", "new_col2"]
        },
        "key_columns": ["scenario", "region", "sector", "year"],
        "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "mean",
        "call_modify_crop_names": false
    }
]
```

---

## Support and Feedback

For issues or questions:
1. Verify your JSON configuration is valid (use a JSON validator)
2. Check that utility modules are in the Python path
3. Review console output for specific error messages
4. Ensure input file format matches expected structure
5. Test with a small subset of data first

This documentation provides comprehensive guidance for using the GCAM data processing script to prepare GCAM model outputs for analysis and visualization.