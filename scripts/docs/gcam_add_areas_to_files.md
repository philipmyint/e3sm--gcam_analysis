# GCAM Add Areas to Files Script Documentation


**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Overview

**Script Name:** `gcam_add_areas_to_files.py`

**Purpose:** Adds land allocation area data as an additional column to GCAM (Global Change Analysis Model) output files. This script matches area information from detailed land allocation files with other GCAM datasets (such as agricultural commodity prices, CO2 emissions, or vegetation/soil scalars) based on scenario, geographical unit, category, and year.

**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Script Description

The `gcam_add_areas_to_files.py` script enriches GCAM output files by adding land area information, enabling area-weighted calculations and analyses. The script matches entries between a target data file and a land allocation reference file based on common identifiers (scenario, region/basin, category, year).

### Key Features

- **Parallel Processing:** Utilizes multiprocessing to add areas to different data subsets simultaneously
- **Flexible Matching:** Supports matching at both regional and basin (watershed) levels
- **JSON Configuration:** All parameters specified through JSON files for easy customization
- **Multiple File Processing:** Can process multiple files in sequence from a single JSON configuration
- **Crop Name Standardization:** Optional conversion to standardized crop naming conventions
- **Area Aggregation:** Automatically sums areas across basins when matching at regional level

### Workflow

1. Reads target data file and land allocation reference file
2. Extracts unique scenarios, geographical units, and categories
3. Creates Cartesian product of all combinations
4. For each combination:
   - Matches entries in both files
   - Extracts corresponding area values
   - Handles missing matches (sets area to 0)
5. Concatenates all subsets into complete dataset
6. Optionally standardizes crop names
7. Outputs enhanced file with area column

---

## Function Documentation

### `add_areas_to_subset_of_file(df, df_land, geographical_label, category_label, scenario, geography, category)`

**Purpose:** Core function that adds land allocation areas to a specific subset of data matching given criteria.

**Parameters:**
- `df` (DataFrame): Target data containing the quantity of interest
- `df_land` (DataFrame): Reference data containing land allocation areas
- `geographical_label` (str): Column name for geographical unit ('region' or 'basin')
- `category_label` (str): Column name for category ('sector' or 'landtype')
- `scenario` (str): Scenario/simulation name to match
- `geography` (str): Specific geographical unit name (e.g., 'USA', 'Amazon_Basin')
- `category` (str): Specific category name (e.g., 'Corn', 'Forest')

**Returns:** DataFrame with added 'area' column

**Process:**
1. Filters both DataFrames to matching scenario, geography, and category
2. Identifies year range from target data
3. If no matching land allocation data found, sets areas to 0
4. For regional matching: Sums areas across all basins for each year
5. For basin matching: Aggregates areas across all regions containing that basin
6. Returns enhanced DataFrame with area column

---

### `add_areas_to_file(inputs)`

**Purpose:** Orchestrates the process of adding areas to an entire file based on JSON configuration.

**Parameters:**
- `inputs` (dict): Dictionary containing all configuration parameters (see Input Parameters Table below)

**Returns:** None (writes output to file)

**Process:**
1. Unpacks input parameters from dictionary
2. Reads target file and land allocation file into DataFrames
3. Extracts unique values for scenarios, geographies, and categories
4. Creates Cartesian product of all combinations
5. Processes each combination in parallel using multiprocessing
6. Concatenates results and sorts by key columns
7. Optionally applies crop name standardization
8. Writes enhanced data to output file
9. Reports execution time

---

## Input Parameters Table

| Parameter | Default Value | Possible Values | Required? | Description |
|-----------|--------------|-----------------|-----------|-------------|
| `input_file` | None | Valid file path to CSV or .dat file | **Yes** | Path to the input file that needs area data added |
| `output_file` | None | Valid file path to CSV or .dat file | **Yes** | Path where the enhanced file with areas will be saved |
| `key_columns` | None | List of valid column names | **Yes** | Columns used for sorting the final output (e.g., `["scenario", "region", "sector", "year"]`) |
| `geographical_label` | None | `"region"` or `"basin"` | **Yes** | Column name identifying geographical units in the data |
| `category_label` | None | `"sector"`, `"landtype"`, or any valid column name | **Yes** | Column name identifying categories (e.g., crop types, sectors) |
| `land_allocation_file` | None | Valid file path to land allocation CSV or .dat file | **Yes** | Path to reference file containing land allocation area data |
| `call_modify_crop_names` | `False` | `true` or `false` | No | Whether to standardize crop names using GCAM mapping conventions |
| `mean_or_sum_if_more_than_one_row_in_same_landtype_group` | `None` | `"area_weighted_mean"`, `"sum"`, or `None` | No | How to aggregate when multiple rows map to same standardized landtype |

---

## Detailed Parameter Descriptions

### Required Parameters

#### `input_file`
**Type:** String (file path)

**Description:** Path to the GCAM output file that requires area information.

**Examples:**
```json
"input_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv"
"input_file": "/path/to/your/data/emissions_by_sector.csv"
```

**Supported Formats:** CSV (.csv) or fixed-width format (.dat) files

**Expected Structure:** Must contain columns for scenario, geographical unit, category (if applicable), and year

---

#### `output_file`
**Type:** String (file path)

**Description:** Path where the enhanced file with area column will be saved.

**Examples:**
```json
"output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv"
"output_file": "/path/to/output/emissions_with_areas.csv"
```

**Note:** Can be the same as `input_file` to overwrite the original file, or a different path to preserve the original.

**Output Format:** Matches input format (CSV or .dat)

---

#### `key_columns`
**Type:** List of strings

**Description:** Column names used for sorting the final output DataFrame. Order matters - data will be sorted by the first column, then second, etc.

**Examples:**
```json
"key_columns": ["scenario", "region", "sector", "year"]
"key_columns": ["scenario", "region", "basin", "landtype", "year"]
```

**Common Patterns:**
- Regional data: `["scenario", "region", "category", "year"]`
- Basin data: `["scenario", "region", "basin", "landtype", "year"]`
- Sector data: `["scenario", "region", "sector", "year"]`

**Important:** All columns listed must exist in the input file.

---

#### `geographical_label`
**Type:** String

**Description:** Column name that identifies geographical units in the data.

**Possible Values:**
- `"region"` - For GCAM regions (e.g., USA, China, Brazil)
- `"basin"` - For watersheds/catchments (e.g., Amazon_Basin, Mississippi_Basin)

**Examples:**
```json
"geographical_label": "region"
"geographical_label": "basin"
```

**Behavior:**
- **Region level:** Areas are summed across all basins within each region
- **Basin level:** Areas are collected from all regions that contain portions of that basin

---

#### `category_label`
**Type:** String

**Description:** Column name that identifies the category of data (e.g., crop types, sectors, land uses).

**Possible Values:** Any valid column name from the input file
- `"sector"` - For economic sectors
- `"landtype"` - For land use types or crop types
- `"commodity"` - For agricultural commodities
- Custom column names as needed

**Examples:**
```json
"category_label": "sector"
"category_label": "landtype"
```

**Usage Context:**
- Agricultural commodity prices: `"sector"`
- Land allocation data: `"landtype"`
- Emissions data: `"sector"`
- Vegetation/soil scalars: `"landtype"`

---

#### `land_allocation_file`
**Type:** String (file path)

**Description:** Path to the reference file containing detailed land allocation area data. This file must have been previously extracted and processed from GCAM outputs.

**Examples:**
```json
"land_allocation_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed.csv"
"land_allocation_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_original_crop_names.csv"
```

**Required Structure:**
- Must contain columns: `scenario`, `region`, `basin`, `landtype`, `year`, `value`
- The `value` column contains area measurements
- Should cover all scenarios, regions, basins, and years present in input file

**Important Notes:**
- Choose the correct land allocation file based on crop naming convention
- For standardized crop names: use `land_allocation_processed.csv`
- For original GCAM crop names: use `land_allocation_processed_original_crop_names.csv`

---

### Optional Parameters

#### `call_modify_crop_names`
**Type:** Boolean

**Default:** `false`

**Description:** Determines whether to apply crop name standardization after adding areas. When enabled, converts original GCAM crop names to standardized names defined in the `gcam_crop_mappings` dictionary.

**Possible Values:**
- `true` - Apply crop name standardization
- `false` - Keep original crop names

**Examples:**
```json
"call_modify_crop_names": true
"call_modify_crop_names": false
```

**When to Use:**
- Set to `true` when working with vegetation/soil scalar files
- Set to `true` when harmonizing crop names across different data sources
- Set to `false` for regional/sectoral data that doesn't involve crop-specific categories

**Standardization Examples:**
- `biomass`, `biomassGrass`, `biomassTree` → `BioenergyCrop`
- `Forest`, `ProtectedUnmanagedForest`, `UnmanagedForest` → `forest` (aggregate)
- See `utility_gcam.py` for complete mapping definitions

---

#### `mean_or_sum_if_more_than_one_row_in_same_landtype_group`
**Type:** String or None

**Default:** `None`

**Description:** Specifies how to aggregate data when crop name standardization causes multiple rows to map to the same standardized landtype.

**Possible Values:**
- `"area_weighted_mean"` - Calculate weighted average using area as weight
- `"sum"` - Sum all values
- `None` - No aggregation (may result in duplicate rows)

**Examples:**
```json
"mean_or_sum_if_more_than_one_row_in_same_landtype_group": "area_weighted_mean"
"mean_or_sum_if_more_than_one_row_in_same_landtype_group": "sum"
```

**When to Use:**
- `"area_weighted_mean"` - For intensive properties (e.g., vegetation scalars, soil carbon density)
- `"sum"` - For extensive properties (e.g., total emissions, total production)
- `None` - When no standardization is applied or aggregation is not needed

**Usage Context:**
- Required when `call_modify_crop_names` is `true`
- Particularly important for scalar files where values represent rates or densities

---

## JSON Configuration File Structure

### File Format

The script accepts JSON configuration files containing an array of configuration objects. Each object represents one file to process.

**Basic Structure:**
```json
[
    {
        "input_file": "path/to/input.csv",
        "output_file": "path/to/output.csv",
        "key_columns": ["scenario", "region", "sector", "year"],
        "geographical_label": "region",
        "category_label": "sector",
        "land_allocation_file": "path/to/land_allocation.csv"
    },
    {
        // Additional file configurations...
    }
]
```

### Example Configurations

#### Example 1: Regional Agricultural Commodity Prices
```json
{
    "input_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
    "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
    "key_columns": ["scenario", "region", "sector", "year"],
    "geographical_label": "region",
    "category_label": "sector",
    "land_allocation_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed.csv"
}
```

**Use Case:** Adding areas to agricultural commodity price data at the regional level.

**Key Features:**
- Matches at regional level (areas summed across basins)
- Uses "sector" as category label
- No crop name standardization needed
- Output overwrites input file

---

#### Example 2: Basin-Level Vegetation and Soil Scalars with Crop Standardization
```json
{
    "input_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback.csv",
    "output_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback.csv",
    "key_columns": ["scenario", "region", "basin", "landtype", "year"],
    "geographical_label": "basin",
    "category_label": "landtype",
    "land_allocation_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_original_crop_names.csv",
    "call_modify_crop_names": true,
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "area_weighted_mean"
}
```

**Use Case:** Adding areas to vegetation/soil scalar files at the basin level with crop name harmonization.

**Key Features:**
- Matches at basin level (more granular than regional)
- Uses "landtype" as category label
- Applies crop name standardization
- Uses area-weighted mean for aggregation (appropriate for intensive properties like scalars)
- Uses land allocation file with original crop names for proper matching

---

#### Example 3: Ensemble Data Processing
```json
{
    "input_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed_ensemble.csv",
    "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed_ensemble.csv",
    "key_columns": ["scenario", "region", "sector", "year"],
    "geographical_label": "region",
    "category_label": "sector",
    "land_allocation_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_ensemble.csv"
}
```

**Use Case:** Processing ensemble simulation data with multiple scenario variations.

**Key Features:**
- Works with ensemble files containing multiple scenario variations
- Uses corresponding ensemble land allocation file
- Same configuration structure as non-ensemble files

---

## Usage Instructions

### Command Line Execution

**Basic Syntax:**
```bash
python gcam_add_areas_to_files.py path/to/config.json
```

**Multiple Configuration Files:**
```bash
python gcam_add_areas_to_files.py config1.json config2.json config3.json
```

**With Full Paths:**
```bash
python gcam_add_areas_to_files.py /full/path/to/gcam_add_areas_to_files.json
```

### Example Workflow

#### Step 1: Prepare Land Allocation File
Ensure you have extracted and processed land allocation data using `gcam_extract_csv_from_project_files.R` and `gcam_process_extracted_data.py`.

```bash
# Extract data from GCAM project files
Rscript gcam_extract_csv_from_project_files.R extraction_config.json

# Process extracted data
python gcam_process_extracted_data.py processing_config.json
```

#### Step 2: Create JSON Configuration
Create a JSON file specifying which files need areas added:

**my_config.json:**
```json
[
    {
        "input_file": "./data/emissions_by_region.csv",
        "output_file": "./data/emissions_by_region_with_areas.csv",
        "key_columns": ["scenario", "region", "sector", "year"],
        "geographical_label": "region",
        "category_label": "sector",
        "land_allocation_file": "./data/land_allocation_processed.csv"
    }
]
```

#### Step 3: Run the Script
```bash
python gcam_add_areas_to_files.py my_config.json
```

#### Step 4: Verify Output
Check the output file to ensure the 'area' column has been added:
```python
import pandas as pd
df = pd.read_csv('./data/emissions_by_region_with_areas.csv')
print(df.head())
print(df.columns)  # Should include 'area' column
```

---

## Understanding Area Matching Logic

### Regional Matching (geographical_label: "region")

**Process:**
1. Filter land allocation data to matching scenario, region, and landtype
2. For each year, sum areas across all basins within that region
3. Assign summed area to corresponding entry in target data

**Example:**
```
Target Data:
scenario | region | sector | year | value
---------|--------|--------|------|-------
base     | USA    | Corn   | 2020 | 150.5

Land Allocation Data:
scenario | region      | basin           | landtype | year | value (area)
---------|-------------|-----------------|----------|------|-------------
base     | USA         | Mississippi     | Corn     | 2020 | 45.2
base     | USA         | Colorado        | Corn     | 2020 | 32.1
base     | USA         | Columbia        | Corn     | 2020 | 28.8

Result:
scenario | region | sector | year | value | area
---------|--------|--------|------|-------|-------
base     | USA    | Corn   | 2020 | 150.5 | 106.1  (sum of all basins)
```

---

### Basin Matching (geographical_label: "basin")

**Process:**
1. Filter land allocation data to matching scenario, basin, and landtype
2. Some basins span multiple regions - collect areas from all relevant regions
3. Assign area to corresponding entry for each region-basin combination

**Example:**
```
Target Data:
scenario | region    | basin       | landtype | year | value
---------|-----------|-------------|----------|------|-------
base     | USA       | Colorado    | Forest   | 2020 | 82.3
base     | Mexico    | Colorado    | Forest   | 2020 | 18.5

Land Allocation Data:
scenario | region    | basin       | landtype | year | value (area)
---------|-----------|-------------|----------|------|-------------
base     | USA       | Colorado    | Forest   | 2020 | 125.4
base     | Mexico    | Colorado    | Forest   | 2020 | 34.2

Result:
scenario | region    | basin       | landtype | year | value | area
---------|-----------|-------------|----------|------|-------|-------
base     | USA       | Colorado    | Forest   | 2020 | 82.3  | 125.4
base     | Mexico    | Colorado    | Forest   | 2020 | 18.5  | 34.2
```

---

### Handling Missing Matches

When no corresponding land allocation data is found for a particular combination:
- Area is set to **0**
- No error is raised
- Processing continues normally

**Common Reasons for Missing Matches:**
- Category exists in target data but not in land allocation
- Geographical unit has no land allocation for that category
- Year range mismatch between files
- Scenario name mismatch

---

## Crop Name Standardization

### Overview

When `call_modify_crop_names` is set to `true`, the script applies standardization mappings defined in the `gcam_crop_mappings` dictionary from `utility_gcam.py`.

### Standard Mappings

**Bioenergy Crops:**
```
biomass, biomassGrass, biomassTree → BioenergyCrop
```

**Forest Aggregation:**
```
Forest, ProtectedUnmanagedForest, UnmanagedForest → forest
```

**Grassland Aggregation:**
```
Grassland, ProtectedUnmanagedGrassland, UnmanagedGrassland → grassland
```

**Shrubland Aggregation:**
```
Shrubland, ProtectedUnmanagedShrubland, UnmanagedShrubland → shrubland
```

**Other Examples:**
```
Corn → Corn (unchanged)
FodderGrass → FodderGrass (unchanged)
Rice → Rice (unchanged)
```

### Aggregation After Standardization

When multiple original crop types map to the same standardized name, aggregation is controlled by `mean_or_sum_if_more_than_one_row_in_same_landtype_group`:

**Area-Weighted Mean (for intensive properties):**
```python
# Example: Vegetation scalars
biomass:      area=50, vegetation=1.2
biomassGrass: area=30, vegetation=1.5
biomassTree:  area=20, vegetation=1.1

# After standardization to BioenergyCrop:
BioenergyCrop: area=100, vegetation=(50*1.2 + 30*1.5 + 20*1.1)/100 = 1.27
```

**Sum (for extensive properties):**
```python
# Example: Total emissions
biomass:      area=50, emissions=120
biomassGrass: area=30, emissions=85
biomassTree:  area=20, emissions=45

# After standardization to BioenergyCrop:
BioenergyCrop: area=100, emissions=250 (120+85+45)
```

---

## Performance Considerations

### Multiprocessing

The script uses parallel processing with limited cores to reduce memory pressure for parallel processing:
```python
with multiprocessing.Pool(processes=MAX_PROCESSES) as pool:
    dataframes_for_each_subset = list(pool.starmap(add_areas_to_subset_of_file, cartesian_product))
```

**Performance Factors:**
- Number of scenarios
- Number of geographical units (regions/basins)
- Number of categories (sectors/landtypes)
- Number of years
- CPU core count

**Scalability:**
- Processing time scales roughly linearly with Cartesian product size
- More CPU cores = faster processing for large datasets

### Execution Time

**Typical Performance:**
- Small files (< 10,000 rows): 1-5 seconds
- Medium files (10,000-100,000 rows): 5-30 seconds
- Large files (> 100,000 rows): 30-120 seconds
- Ensemble files: Multiply by number of scenarios

**Timing Output:**
```
Elapsed time for adding areas to ag_commodity_prices_processed.csv: 12.34 seconds
Elapsed time for adding areas to scalars_control+full_feedback.csv: 45.67 seconds
Elapsed time for adding areas to all files: 58.01 seconds
```

### Memory Requirements

**Memory Usage Estimation:**
```
Peak Memory ≈ (Input File Size + Land Allocation File Size) × Number of CPU Cores × 1.5
```

**Example:**
- Input file: 50 MB
- Land allocation file: 200 MB
- CPU cores: 8
- Estimated peak memory: (50 + 200) × 8 × 1.5 = 3 GB

**Optimization Tips:**
- Process files with fewer scenarios/categories first
- Close unnecessary applications to free memory
- For very large ensembles, process subset of files at a time

---

## Dependencies

### Required Python Modules

```python
import itertools        # Standard library
import json            # Standard library
import multiprocessing # Standard library
import pandas as pd    # Install: pip install pandas
import sys            # Standard library
import time           # Standard library
```

### Custom Utility Modules

**From `utility_constants.py`:**
- Various constants used throughout the analysis pipeline

**From `utility_dataframes.py`:**
- `read_file_into_dataframe(file)` - Reads CSV or .dat files
- `write_dataframe_to_file(df, file)` - Writes to CSV or .dat files

**From `utility_gcam.py`:**
- `modify_crop_names(df, key_columns, aggregation_method)` - Standardizes crop names

**Installation:**
Ensure all utility modules are available in the Python path or the same directory as the script.

---

## Output File Structure

### Added Column

The script adds a single new column to the input file:

**Column Name:** `area`

**Data Type:** Float

**Units:** Typically thousands of square kilometers (km²), depending on GCAM output units

**Position:** Appended as the last column (rightmost)

### Example Output Structure

**Before (Input):**
```csv
scenario,region,sector,year,value
base,USA,Corn,2020,150.5
base,USA,Corn,2025,158.2
base,USA,Wheat,2020,95.3
```

**After (Output):**
```csv
scenario,region,sector,year,value,area
base,USA,Corn,2020,150.5,106.1
base,USA,Corn,2025,158.2,108.7
base,USA,Wheat,2020,95.3,82.4
```

### Data Sorting

Output is sorted according to `key_columns` specification:
- Primary sort: First column in key_columns
- Secondary sort: Second column in key_columns
- And so on...

---

## Common Use Cases

### Use Case 1: Area-Weighted Regional Analysis

**Goal:** Calculate area-weighted average agricultural commodity prices by region

**Configuration:**
```json
{
    "input_file": "./data/ag_prices.csv",
    "output_file": "./data/ag_prices_with_areas.csv",
    "key_columns": ["scenario", "region", "sector", "year"],
    "geographical_label": "region",
    "category_label": "sector",
    "land_allocation_file": "./data/land_allocation.csv"
}
```

**Post-Processing:**
```python
import pandas as pd
df = pd.read_csv('./data/ag_prices_with_areas.csv')

# Calculate area-weighted average price
weighted_avg = (df['value'] * df['area']).sum() / df['area'].sum()
```

---

### Use Case 2: Basin-Level Ecosystem Analysis

**Goal:** Analyze vegetation scalars at the basin level with standardized crop names

**Configuration:**
```json
{
    "input_file": "./data/veg_scalars.csv",
    "output_file": "./data/veg_scalars_with_areas.csv",
    "key_columns": ["scenario", "region", "basin", "landtype", "year"],
    "geographical_label": "basin",
    "category_label": "landtype",
    "land_allocation_file": "./data/land_allocation_original_names.csv",
    "call_modify_crop_names": true,
    "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "area_weighted_mean"
}
```

**Benefits:**
- Harmonized crop names across datasets
- Area-weighted aggregation for physically meaningful averages
- Basin-level spatial resolution for watershed analysis

---

### Use Case 3: Ensemble Time Series Comparison

**Goal:** Compare land use changes across ensemble members with area weighting

**Configuration:**
```json
{
    "input_file": "./data/land_use_ensemble.csv",
    "output_file": "./data/land_use_ensemble_with_areas.csv",
    "key_columns": ["scenario", "region", "landtype", "year"],
    "geographical_label": "region",
    "category_label": "landtype",
    "land_allocation_file": "./data/land_allocation_ensemble.csv"
}
```

**Analysis:**
```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('./data/land_use_ensemble_with_areas.csv')

# Group by scenario and year, calculate area-weighted mean
grouped = df.groupby(['scenario', 'year']).apply(
    lambda x: (x['value'] * x['area']).sum() / x['area'].sum()
)

# Plot ensemble spread
grouped.unstack(level=0).plot()
plt.ylabel('Area-weighted mean land use change')
plt.show()
```

---

## Troubleshooting

### Common Issues and Solutions

#### Issue 1: FileNotFoundError
**Error Message:**
```
FileNotFoundError: [Errno 2] No such file or directory: './data/input.csv'
```

**Causes:**
- Incorrect file path in JSON configuration
- File hasn't been created yet
- Working directory is different than expected

**Solutions:**
- Verify file paths are correct (use absolute paths if uncertain)
- Check that input files exist before running
- Use `pwd` command to check current working directory

---

#### Issue 2: KeyError on Column Names
**Error Message:**
```
KeyError: 'sector'
```

**Causes:**
- Column specified in JSON doesn't exist in input file
- Typo in column name (case-sensitive)
- File structure different than expected

**Solutions:**
- Verify column names using: `df.columns.tolist()`
- Check for exact spelling and capitalization
- Inspect a few rows of the input file: `df.head()`

---

#### Issue 3: No Areas Added (All Zeros)
**Symptom:** Area column added but all values are 0

**Causes:**
- Mismatch between category names in input and land allocation files
- Scenario names don't match
- Geographical units don't match
- Year ranges don't overlap

**Solutions:**
- Check unique values: `df['category_label'].unique()`
- Verify scenario names match exactly
- Ensure land allocation file covers the same time period
- If using crop name standardization, ensure land allocation file uses corresponding naming convention

---

#### Issue 4: Memory Error
**Error Message:**
```
MemoryError: Unable to allocate array
```

**Causes:**
- Files too large for available RAM
- Too many simultaneous parallel processes
- Cartesian product creates enormous intermediate dataset

**Solutions:**
- Process fewer files at once
- Reduce number of CPU cores (won't help much, but worth trying)
- Split large files into smaller subsets
- Close other applications to free memory
- Use machine with more RAM

---

#### Issue 5: Slow Performance
**Symptom:** Script takes much longer than expected

**Causes:**
- Very large Cartesian product (many scenarios × regions × categories)
- Large land allocation file requiring repeated reads
- Insufficient CPU cores

**Solutions:**
- Check size of Cartesian product: scenarios × geographies × categories
- Consider processing subsets of data separately
- Ensure multiprocessing is working (should use multiple cores)
- Verify no other intensive processes running

---

#### Issue 6: Incorrect Aggregation After Crop Standardization
**Symptom:** Values don't make sense after aggregation

**Causes:**
- Using "sum" for intensive properties (should use area-weighted mean)
- Using area-weighted mean for extensive properties (should use sum)

**Solutions:**
- Use `"area_weighted_mean"` for: scalars, densities, rates, concentrations
- Use `"sum"` for: totals, absolute quantities, production volumes
- Verify units of the quantity being aggregated

---

## Validation and Quality Checks

### Recommended Checks After Running

#### 1. Verify Area Column Added
```python
import pandas as pd
df = pd.read_csv('output_file.csv')
assert 'area' in df.columns, "Area column not added!"
print(f"✓ Area column successfully added")
```

#### 2. Check for Missing Areas
```python
zero_areas = df[df['area'] == 0]
if len(zero_areas) > 0:
    print(f"⚠ Warning: {len(zero_areas)} rows have zero area")
    print(zero_areas[['scenario', 'region', 'category', 'year']].head())
else:
    print("✓ No zero areas found")
```

#### 3. Verify Area Distributions
```python
import matplotlib.pyplot as plt

df.boxplot(column='area', by='category')
plt.title('Area Distribution by Category')
plt.ylabel('Area (1000 km²)')
plt.show()
```

#### 4. Check Total Areas by Scenario
```python
total_areas = df.groupby('scenario')['area'].sum()
print("Total areas by scenario:")
print(total_areas)
```

#### 5. Validate Sorting
```python
# Verify data is sorted according to key_columns
key_cols = ['scenario', 'region', 'sector', 'year']
is_sorted = (df[key_cols] == df[key_cols].sort_values(key_cols)).all().all()
print(f"✓ Data properly sorted: {is_sorted}")
```

---

## Best Practices

### 1. File Organization
- Keep land allocation files in a central reference directory
- Use descriptive output file names
- Maintain separate directories for processed vs. raw data
- Document which land allocation file corresponds to which datasets

### 2. Configuration Management
- Create separate JSON files for different analysis workflows
- Comment JSON files (use external documentation file)
- Version control JSON configurations
- Use meaningful file and parameter names

### 3. Data Quality
- Always validate land allocation files before use
- Check for completeness (all scenarios, regions, years covered)
- Verify units are consistent across files
- Inspect sample outputs before processing entire ensemble

### 4. Processing Strategy
- Process smaller test files first
- Start with single files before batch processing
- Monitor memory usage for large datasets
- Keep backups of original files when overwriting

### 5. Documentation
- Document which land allocation file was used for each analysis
- Record any crop name standardization applied
- Note aggregation methods used
- Keep log of execution times for performance tracking

---

## Integration with Other Scripts

### Typical Workflow Position

```
1. gcam_extract_csv_from_project_files.R
   ↓ (Extract raw data from GCAM project files)
   
2. gcam_process_extracted_data.py
   ↓ (Clean and organize extracted data)
   
3. gcam_add_areas_to_files.py  ← THIS SCRIPT
   ↓ (Enrich with land allocation areas)
   
4. gcam_plot_time_series.py / gcam_plot_spatial_data.py
   ↓ (Visualize and analyze with area weighting)
```

### Upstream Dependencies

**Required predecessor scripts:**
- `gcam_extract_csv_from_project_files.R` - Must run first to extract land allocation data
- `gcam_process_extracted_data.py` - Should run to clean and organize data

### Downstream Usage

**Scripts that benefit from added area data:**
- `gcam_plot_time_series.py` - Can create area-weighted plots
- `gcam_plot_box_and_whiskers.py` - Can display area-weighted distributions
- `gcam_plot_spatial_data.py` - Can show area-normalized spatial patterns

---

## Related Scripts

### Complementary GCAM Scripts

**`gcam_extract_csv_from_project_files.R`**
- Extracts data from GCAM project files (including land allocations)
- Creates the land allocation reference files used by this script

**`gcam_process_extracted_data.py`**
- Cleans and organizes extracted CSV files
- Prepares files for area addition

**`gcam_compile_ehc_scalars.py`**
- Compiles vegetation and soil scalar files
- Output files often need areas added for weighted analysis

**`gcam_plot_time_series.py`**
- Plots GCAM time series with optional area weighting
- Directly uses output from this script

### Related E3SM Scripts

**`e3sm_extract_time_series_surfdata_iesm_dyn.py`**
- Extracts land surface data from E3SM
- Creates analogous area data for E3SM analysis

---

## Example JSON Configuration File

Here is the complete example configuration file from the repository:

```json
[
    {
        "input_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
        "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed.csv",
        "key_columns": ["scenario", "region", "sector", "year"],
        "geographical_label": "region",
        "category_label": "sector",
        "land_allocation_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed.csv"
    },
    {
        "input_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed_ensemble.csv",
        "output_file": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices_processed_ensemble.csv",
        "key_columns": ["scenario", "region", "sector", "year"],
        "geographical_label": "region",
        "category_label": "sector",
        "land_allocation_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_ensemble.csv"
    },
    {
        "input_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback.csv",
        "output_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback.csv",
        "key_columns": ["scenario", "region", "basin", "landtype", "year"],
        "geographical_label": "basin",
        "category_label": "landtype",
        "land_allocation_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_original_crop_names.csv",
        "call_modify_crop_names": true,
        "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "area_weighted_mean"
    },
    {
        "input_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback_ensemble.csv",
        "output_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback_ensemble.csv",
        "key_columns": ["scenario", "region", "basin", "landtype", "year"],
        "geographical_label": "basin",
        "category_label": "landtype",
        "land_allocation_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_original_crop_names_ensemble.csv",
        "call_modify_crop_names": true,
        "mean_or_sum_if_more_than_one_row_in_same_landtype_group": "area_weighted_mean"
    }
]
```

This configuration processes four files:
1. Base agricultural commodity prices
2. Ensemble agricultural commodity prices
3. Base vegetation/soil scalars with crop standardization
4. Ensemble vegetation/soil scalars with crop standardization

---

## References

- DiVittorio et al. (2025). "E3SM-GCAM coupling methodology and applications." *Journal of Advances in Modeling Earth Systems*. [DOI: 10.1029/2024MS004806](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024MS004806)
- GitHub Repository: [https://github.com/philipmyint/e3sm--gcam_analysis](https://github.com/philipmyint/e3sm--gcam_analysis)
- GCAM Documentation: [https://gcims.pnnl.gov/modeling/gcam-global-change-analysis-model](https://gcims.pnnl.gov/modeling/gcam-global-change-analysis-model)

---

## Contact

For questions or issues:
- **Philip Myint**: myint1@llnl.gov
- **Dalei Hao**: dalei.hao@pnnl.gov

For bug reports or feature requests, please create an issue on the [GitHub repository](https://github.com/philipmyint/e3sm--gcam_analysis/issues).

---

## Version History

- **Current Version:** As documented in repository (January 2025)
- **Last Updated:** January 2026

---

## Appendix: Quick Reference

### Quick Start Checklist

- [ ] Install required Python packages (pandas)
- [ ] Ensure utility modules are accessible
- [ ] Prepare land allocation reference file
- [ ] Create JSON configuration file
- [ ] Verify all file paths are correct
- [ ] Run script with JSON file
- [ ] Validate output has 'area' column
- [ ] Check for zero or missing areas
- [ ] Proceed with downstream analysis

### Command Quick Reference

```bash
# Basic execution
python gcam_add_areas_to_files.py config.json

# Multiple configs
python gcam_add_areas_to_files.py config1.json config2.json

# View help (add to script if desired)
python gcam_add_areas_to_files.py --help
```

### JSON Template

```json
[
    {
        "input_file": "path/to/input.csv",
        "output_file": "path/to/output.csv",
        "key_columns": ["scenario", "region", "category", "year"],
        "geographical_label": "region or basin",
        "category_label": "sector or landtype",
        "land_allocation_file": "path/to/land_allocation.csv",
        "call_modify_crop_names": false,
        "mean_or_sum_if_more_than_one_row_in_same_landtype_group": null
    }
]
```

---

*This documentation provides comprehensive guidance for using the `gcam_add_areas_to_files.py` script to enrich GCAM output files with land allocation area data for area-weighted analyses.*
