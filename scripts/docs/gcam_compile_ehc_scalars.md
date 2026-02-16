# GCAM EHC Scalar Compilation Script Documentation

**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Overview

The `gcam_compile_ehc_scalars.py` script is a Python tool for compiling and processing E3SM Human Component (EHC) scalar files generated during coupled E3SM–GCAM simulations. During runtime, the EHC dynamically generates CSV files containing soil and vegetation multipliers (referred to as "scalars") that are passed to GCAM. This script aggregates multiple scalar files from different directories (representing different scenarios) into a single consolidated CSV file for analysis and visualization.

**Related Publication:** The scalars compiled by this script are described in detail in [DiVittorio et al. (2025)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024MS004806), published in the Journal of Advances in Modeling Earth Systems.

## Table of Contents

1. [Installation & Requirements](#installation--requirements)
2. [Background: What are EHC Scalars?](#background-what-are-ehc-scalars)
3. [Basic Usage](#basic-usage)
4. [Configuration File Structure](#configuration-file-structure)
5. [Configuration Parameters](#configuration-parameters)
6. [Data Processing Operations](#data-processing-operations)
7. [Examples](#examples)
8. [Understanding EHC Scalar Files](#understanding-ehc-scalar-files)
9. [Integration with E3SM–GCAM Workflow](#integration-with-e3smgcam-workflow)
10. [Troubleshooting](#troubleshooting)

---

## Installation & Requirements

### Required Python Libraries

```bash
pip install pandas
```

### Required Utility Modules

The script imports several utility modules that must be in the same directory or Python path:
- `utility_dataframes` - Functions for reading/writing data files
- `utility_functions` - General utility functions including file path operations
- `utility_gcam` - GCAM-specific functions for crop name standardization

### System Requirements

- Python 3.7+
- Multi-core processor (script uses parallel processing with limited cores to reduce memory pressure for parallel processing)
- Sufficient memory to handle multiple CSV files simultaneously

---

## Background: What are EHC Scalars?

### E3SM Human Component (EHC)

In coupled E3SM–GCAM simulations, the E3SM Human Component (EHC) serves as the intermediary between:
- **E3SM** (Energy Exascale Earth System Model) - Earth system model
- **GCAM** (Global Change Analysis Model) - Human-Earth system model

### Scalars Explained

**Scalars** are multipliers that represent environmental conditions affecting agricultural productivity:

| Scalar Type | Description | Use in GCAM |
|-------------|-------------|-------------|
| **Vegetation Scalars** | Multipliers representing above-ground biomass productivity | Adjust crop yields based on climate conditions |
| **Soil Scalars** | Multipliers representing soil carbon and nutrient availability | Adjust land productivity based on soil health |

### Runtime Generation

During a coupled simulation:
1. E3SM calculates land surface conditions (temperature, precipitation, soil moisture, etc.)
2. EHC translates these conditions into vegetation and soil scalars
3. EHC outputs these scalars to CSV files for each time step
4. These scalars are passed to GCAM to inform land use decisions
5. GCAM makes land allocation decisions based on these productivity signals

### File Organization

EHC generates scalar files organized by:
- **Scenario** - Different simulation scenarios (e.g., Control, Full feedback)
- **Time step** - Typically annual or 5-year intervals
- **Region** - GCAM geographic regions
- **Basin** - Water basins/catchments within regions
- **Land type** - Different land cover types (crop types, forest, grassland, etc.)

---

## Basic Usage

### Command Line Execution

```bash
python gcam_compile_ehc_scalars.py path/to/config.json
```

You can specify multiple JSON configuration files:

```bash
python gcam_compile_ehc_scalars.py config1.json config2.json config3.json
```

The script processes all configurations in parallel using parallel processing with limited cores to reduce memory pressure.

### What the Script Does

For each configuration in the JSON file, the script:

1. **Locates** all CSV files in the specified input directories
2. **Reads** all scalar files from each directory
3. **Concatenates** files within each scenario
4. **Splits** the `landtype_basin` column into separate `landtype` and `basin` columns
5. **Assigns** scenario labels to distinguish different simulation runs
6. **Sorts** data by key columns (scenario, region, basin, landtype, year)
7. **Optionally standardizes** crop names to a common naming convention
8. **Writes** the compiled data to a single output CSV file

---

## Configuration File Structure

The configuration file is a JSON array where each object specifies one compilation task. Here's a minimal example:

```json
[
    {
        "input_directories": [
            "./scalars_output/control",
            "./scalars_output/full_feedback"
        ],
        "scenarios": ["Control", "Full feedback"],
        "output_file": "./compiled_scalars.csv"
    }
]
```

---

## Configuration Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `input_directories` | list of strings | Paths to directories containing EHC scalar CSV files |
| `scenarios` | list of strings | Names to assign to each scenario (must match order of `input_directories`) |
| `output_file` | string | Path where the compiled CSV file will be saved |

**Critical Requirement:** The number of entries in `scenarios` must exactly match the number of entries in `input_directories`, and they must be in the same order.

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `call_modify_crop_names` | boolean | `false` | Whether to standardize crop names to common naming convention |

---

## Data Processing Operations

### 1. Directory Scanning and File Reading

The script automatically finds all files in each specified directory:

```
./scalars_output/control/
├── scalars_2015.csv
├── scalars_2020.csv
├── scalars_2025.csv
└── scalars_2030.csv
```

All CSV files are read and concatenated together for each scenario.

### 2. Column Name Normalization

All column names are automatically converted to lowercase for consistency:

**Before:** `Landtype_Basin`, `Region`, `YEAR`, `Vegetation`, `Soil`  
**After:** `landtype_basin`, `region`, `year`, `vegetation`, `soil`

### 3. Column Splitting

The `landtype_basin` column contains composite information that gets split:

**Input format:**
```
landtype_basin
--------------
Corn_Amazon
Wheat_Mississippi
Forest_Nile
```

**Output format:**
```
landtype    basin
--------    -----
Corn        Amazon
Wheat       Mississippi
Forest      Nile
```

The original `landtype_basin` column is removed after splitting.

### 4. Scenario Labeling

Each row is tagged with its corresponding scenario name:

```
scenario        region  basin   landtype  year  vegetation  soil
-----------     ------  -----   --------  ----  ----------  ----
Control         USA     Amazon  Corn      2020  1.05        0.98
Full feedback   USA     Amazon  Corn      2020  1.08        0.95
```

### 5. Data Sorting

Data is sorted by key columns in this order:
1. scenario
2. region
3. basin
4. landtype
5. year

This ensures consistent organization across all output files.

### 6. Crop Name Standardization (Optional)

When `call_modify_crop_names` is `true`, crop names are standardized to match a common naming convention:

**Example transformations:**
- `biomass`, `biomassGrass`, `biomassTree` → `BioenergyCrop`
- `CornC4` → `Corn`
- `OilPalmTree` → `OilPalm`
- `SugarCropC4` → `SugarCrop`

See the `utility_gcam.py` module for the complete crop name mapping.

**Why standardize?** Different GCAM output files use different crop naming conventions. Standardization ensures consistency when combining data from multiple sources.

---

## Examples

### Example 1: Compiling Four Scenarios

**Purpose:** Compile scalars from four different simulation scenarios into a single file.

```json
{
    "input_directories": [
        "./../2025_DiVittorio_et_al_gcam/scalars_ehc_output/control",
        "./../2025_DiVittorio_et_al_gcam/scalars_ehc_output/full_feedback",
        "./../2025_DiVittorio_et_al_gcam/scalars_ehc_output/ag_scaling",
        "./../2025_DiVittorio_et_al_gcam/scalars_ehc_output/carbon_scaling"
    ],
    "scenarios": ["Control", "Full feedback", "Ag scaling", "Carbon scaling"],
    "output_file": "./../2025_DiVittorio_et_al_gcam/scalars_all_scenarios.csv"
}
```

**What happens:**
1. Reads all CSV files from each of the four directories
2. Labels each dataset with its scenario name
3. Combines all four scenarios into one file
4. Splits `landtype_basin` column
5. Sorts by scenario, region, basin, landtype, and year
6. Outputs to `scalars_all_scenarios.csv`

**Input structure (files in each directory):**
```
control/
├── scalars_2015_region1.csv
├── scalars_2015_region2.csv
├── scalars_2020_region1.csv
└── scalars_2020_region2.csv

full_feedback/
├── scalars_2015_region1.csv
├── scalars_2015_region2.csv
├── scalars_2020_region1.csv
└── scalars_2020_region2.csv
...
```

**Output structure:**
```csv
scenario,region,basin,landtype,year,vegetation,soil
Control,USA,Amazon,BioenergyCrop,2015,1.05,0.98
Control,USA,Amazon,BioenergyCrop,2020,1.08,0.97
Control,USA,Amazon,Corn,2015,1.02,1.01
Full feedback,USA,Amazon,BioenergyCrop,2015,1.12,0.95
Full feedback,USA,Amazon,BioenergyCrop,2020,1.15,0.93
...
```

### Example 2: Compiling Two Scenarios for Comparison

**Purpose:** Compile only Control and Full feedback scenarios for focused analysis.

```json
{
    "input_directories": [
        "./../2025_DiVittorio_et_al_gcam/scalars_ehc_output/control",
        "./../2025_DiVittorio_et_al_gcam/scalars_ehc_output/full_feedback"
    ],
    "scenarios": ["Control", "Full feedback"],
    "output_file": "./../2025_DiVittorio_et_al_gcam/scalars_control+full_feedback.csv"
}
```

**Use case:** When you only need to compare a baseline scenario against one alternative scenario.

### Example 3: With Crop Name Standardization

**Purpose:** Compile scenarios while ensuring consistent crop naming.

```json
{
    "input_directories": [
        "./scalars_output/scenario_a",
        "./scalars_output/scenario_b"
    ],
    "scenarios": ["Scenario A", "Scenario B"],
    "output_file": "./scalars_standardized.csv",
    "call_modify_crop_names": true
}
```

**Effect of standardization:**

**Before:**
```
landtype
-----------
biomass
biomassGrass
biomassTree
CornC4
Corn
```

**After:**
```
landtype
-----------
BioenergyCrop
BioenergyCrop
BioenergyCrop
Corn
Corn
```

### Example 4: Multiple Compilations in One Run

**Purpose:** Create both a comprehensive file (all scenarios) and a focused file (two scenarios only).

```json
[
    {
        "input_directories": [
            "./scalars/control",
            "./scalars/full_feedback",
            "./scalars/ag_scaling",
            "./scalars/carbon_scaling"
        ],
        "scenarios": ["Control", "Full feedback", "Ag scaling", "Carbon scaling"],
        "output_file": "./scalars_all.csv"
    },
    {
        "input_directories": [
            "./scalars/control",
            "./scalars/full_feedback"
        ],
        "scenarios": ["Control", "Full feedback"],
        "output_file": "./scalars_baseline_comparison.csv"
    }
]
```

Both compilations are processed in parallel for efficiency.

---

## Understanding EHC Scalar Files

### Typical File Structure

Individual scalar files generated by EHC typically have this structure:

```csv
Region,Landtype_Basin,Year,Vegetation,Soil
USA,Corn_Mississippi,2015,1.05,0.98
USA,Corn_Mississippi,2020,1.08,0.97
USA,Wheat_Amazon,2015,1.02,1.01
Brazil,Soybean_Amazon,2015,1.15,0.95
```

### Column Descriptions

| Column | Description | Typical Values |
|--------|-------------|----------------|
| `Region` | GCAM geographic region | USA, Brazil, China, India, etc. |
| `Landtype_Basin` | Combined land type and water basin | Corn_Amazon, Forest_Nile, etc. |
| `Year` | Simulation year | 2015, 2020, 2025, etc. |
| `Vegetation` | Vegetation scalar (multiplier) | 0.5 to 2.0 (typically near 1.0) |
| `Soil` | Soil scalar (multiplier) | 0.5 to 2.0 (typically near 1.0) |

### Scalar Value Interpretation

**Vegetation and Soil Scalars:**
- **= 1.0** - Normal/baseline productivity
- **> 1.0** - Enhanced productivity (favorable conditions)
- **< 1.0** - Reduced productivity (unfavorable conditions)

**Example:**
- Vegetation scalar of 1.15 = 15% increase in above-ground productivity
- Soil scalar of 0.85 = 15% decrease in soil-related productivity

### File Naming Conventions

EHC scalar files may follow various naming patterns:
- `scalars_YYYY.csv` (by year)
- `scalars_YYYY_regionX.csv` (by year and region)
- `vegetation_soil_scalars_YYYY.csv`

The script handles all files in the directory regardless of naming convention.

---

## Integration with E3SM–GCAM Workflow

### Complete Workflow

```
1. E3SM Simulation Run
   ↓
2. EHC generates scalar files
   ├── control/scalars_*.csv
   ├── full_feedback/scalars_*.csv
   └── other_scenarios/scalars_*.csv
   ↓
3. gcam_compile_ehc_scalars.py ← YOU ARE HERE
   ↓
4. Compiled scalars CSV
   ↓
5. Analysis & Visualization
   ├── gcam_plot_time_series.py
   ├── gcam_plot_box_and_whiskers.py
   └── gcam_plot_spatial_data.py
```

### Relationship to Other Scripts

**Preceding Steps:**
- E3SM–GCAM coupled simulation generates EHC scalar files
- Files are organized by scenario in separate directories

**Following Steps:**
1. **gcam_plot_time_series.py** - Create time series plots of vegetation/soil scalars
2. **gcam_plot_box_and_whiskers.py** - Create box plots showing scalar distributions
3. **gcam_plot_spatial_data.py** - Create spatial maps of scalar values
4. **Statistical analysis** - Compare scenarios, perform hypothesis tests

### Data Flow Example

**Input directories:**
```
scalars_ehc_output/
├── control/
│   ├── scalars_2015.csv (500 KB)
│   ├── scalars_2020.csv (500 KB)
│   └── scalars_2025.csv (500 KB)
└── full_feedback/
    ├── scalars_2015.csv (500 KB)
    ├── scalars_2020.csv (500 KB)
    └── scalars_2025.csv (500 KB)
```

**Output file:**
```
scalars_control+full_feedback.csv (3 MB)
```

All six input files compiled into one organized dataset ready for analysis.

---

## Troubleshooting

### Common Issues

**Problem:** `FileNotFoundError` - Cannot find input directory
```
Error: [Errno 2] No such file or directory: './scalars/control'

Solution: Verify that input_directories paths are correct. Use absolute 
paths if running from a different directory. Check for typos in directory 
names.
```

**Problem:** Scenario/directory count mismatch
```
Error: Length of scenarios does not match length of input_directories

Solution: Ensure the scenarios list has exactly the same number of entries 
as input_directories, and they are in matching order.
```

**Problem:** Empty output file
```
Solution: Check that the input directories actually contain CSV files. 
Verify file permissions allow reading. Ensure CSV files have the expected 
structure with required columns.
```

**Problem:** Missing columns after processing
```
Error: KeyError: 'landtype_basin'

Solution: Verify that input CSV files contain a 'landtype_basin' column 
(case-insensitive). If EHC output has changed format, you may need to 
modify the column splitting logic.
```

**Problem:** Crop name standardization not working
```
Solution: Ensure call_modify_crop_names is set to true (boolean, not 
string "true"). Verify that crop names in your files match those defined 
in utility_gcam.py's gcam_crop_mappings dictionary.
```

### Input File Requirements

Individual EHC scalar files must contain:
- **Region** column - GCAM geographic regions
- **Landtype_Basin** column - Combined land type and basin information
- **Year** column - Simulation year
- **Vegetation** column - Vegetation scalar values
- **Soil** column - Soil scalar values

**Column names are case-insensitive** (they're converted to lowercase during processing).

### Data Validation Tips

After running the script, verify the output:

```python
import pandas as pd

# Load the compiled file
df = pd.read_csv('scalars_all_scenarios.csv')

# Check structure
print(df.head())
print(df.columns)
print(f"Total rows: {len(df)}")
print(f"Scenarios: {df['scenario'].unique()}")
print(f"Years: {sorted(df['year'].unique())}")
print(f"Land types: {df['landtype'].unique()}")

# Check for missing values
print(df.isnull().sum())

# Check value ranges
print(df[['vegetation', 'soil']].describe())
```

### Performance Considerations

**Processing Speed:**
- Small datasets (< 100 files, < 10 MB total): < 5 seconds
- Medium datasets (100-500 files, 10-100 MB): 10-30 seconds
- Large datasets (> 500 files, > 100 MB): 1-5 minutes

**Memory Usage:**
- Scales with total size of input files
- Multiple configurations processed in parallel
- Peak memory ≈ 2-3× total input file size

**Optimization Tips:**
1. **Fewer scenarios per run** - Process subsets if memory is limited
2. **Pre-filter data** - Remove unnecessary regions/years from source files before compilation
3. **Use SSD storage** - Significantly faster file I/O for large datasets

---

## Best Practices

### 1. Organize Source Files

Keep EHC scalar files organized by scenario:
```
scalars_ehc_output/
├── control/
├── full_feedback/
├── ag_scaling/
└── carbon_scaling/
```

### 2. Use Descriptive Scenario Names

Choose clear, meaningful names:
```json
"scenarios": ["Control", "Full feedback", "Ag scaling", "Carbon scaling"]
```

Not:
```json
"scenarios": ["S1", "S2", "S3", "S4"]
```

### 3. Document Your Configuration

Add comments to your JSON (though note that standard JSON doesn't support comments, you can maintain a separate documentation file):

```
scalars_compilation_notes.txt:
- control: Baseline simulation with no feedbacks
- full_feedback: Complete two-way coupling between E3SM and GCAM
- ag_scaling: Agricultural scaling scenario
- carbon_scaling: Carbon scaling scenario
```

### 4. Version Your Output Files

Include dates or version numbers in output filenames:
```json
"output_file": "./scalars_all_scenarios_2025-01-15.csv"
```

### 5. Verify Before Downstream Analysis

Always check compiled files before using them for plotting or analysis:
- Row counts match expectations
- All scenarios present
- Year coverage is complete
- No unexpected missing values
- Scalar values are in reasonable ranges

---

## Advanced Usage

### Compiling Ensemble Members

For ensemble simulations (same scenario with different initial conditions):

```json
[
    {
        "input_directories": [
            "./scalars/control_member1",
            "./scalars/control_member2",
            "./scalars/control_member3",
            "./scalars/control_member4",
            "./scalars/control_member5"
        ],
        "scenarios": [
            "Control_1", "Control_2", "Control_3", 
            "Control_4", "Control_5"
        ],
        "output_file": "./scalars_control_ensemble.csv"
    }
]
```

This creates individual identifiers for each ensemble member while grouping them under the same base scenario name.

### Processing Time Slices

Compile only specific time periods from a long simulation:

**Approach:** Pre-filter files in input directories to include only desired years before running the script.

```bash
# Example: Move only 2020-2050 files to a temporary directory
mkdir temp_2020_2050
cp scalars/*_202[0-9].csv temp_2020_2050/
cp scalars/*_203[0-9].csv temp_2020_2050/
cp scalars/*_204[0-9].csv temp_2020_2050/
cp scalars/*_2050.csv temp_2020_2050/
```

Then reference `temp_2020_2050` in your JSON configuration.

### Parallel Processing of Multiple Compilations

The script automatically processes multiple JSON configuration blocks in parallel:

```json
[
    {
        "input_directories": [...],
        "scenarios": [...],
        "output_file": "./output1.csv"
    },
    {
        "input_directories": [...],
        "scenarios": [...],
        "output_file": "./output2.csv"
    },
    {
        "input_directories": [...],
        "scenarios": [...],
        "output_file": "./output3.csv"
    }
]
```

All three compilations run simultaneously on different CPU cores.

---

## Output File Characteristics

### Structure

Compiled files have a consistent structure:

```csv
scenario,region,basin,landtype,year,vegetation,soil
Control,USA,Amazon,BioenergyCrop,2015,1.05,0.98
Control,USA,Amazon,BioenergyCrop,2020,1.08,0.97
...
```

### Column Order

Fixed column order (left to right):
1. scenario
2. region
3. basin
4. landtype
5. year
6. vegetation
7. soil

### Sorting

Data is sorted hierarchically:
1. Primary sort: scenario (alphabetical)
2. Secondary sort: region (alphabetical)
3. Tertiary sort: basin (alphabetical)
4. Quaternary sort: landtype (alphabetical)
5. Quinary sort: year (numerical, ascending)

### File Format

- **Format:** CSV (Comma-Separated Values)
- **Encoding:** UTF-8
- **Line endings:** Platform-dependent (Unix/Windows)
- **Header:** Always included as first row

---

## Execution Time Examples

Typical execution times based on dataset size:

**Small dataset:**
```
Input: 4 scenarios × 10 files × 50 KB each = 2 MB total
Output: 1 file × 1.5 MB
Elapsed time: 1.23 seconds
```

**Medium dataset:**
```
Input: 4 scenarios × 50 files × 200 KB each = 40 MB total
Output: 1 file × 30 MB
Elapsed time: 8.45 seconds
```

**Large dataset:**
```
Input: 4 scenarios × 200 files × 500 KB each = 400 MB total
Output: 1 file × 300 MB
Elapsed time: 45.67 seconds
```

Console output example:
```
Elapsed time processing/compiling the data for scalars_all_scenarios.csv: 8.45 seconds
Elapsed time processing/compiling the data for scalars_control+full_feedback.csv: 3.21 seconds
Elapsed time processing/compiling the data for all files: 8.89 seconds
```

---

## Related Documentation

For comprehensive understanding of the E3SM–GCAM analysis workflow, refer to:

1. **GCAM Data Processing Script** (`gcam_process_extracted_data.py`)
   - Processes other GCAM CSV outputs
   - Similar workflow but for different data types

2. **GCAM Plotting Scripts**
   - `gcam_plot_time_series.py` - Visualize scalar changes over time
   - `gcam_plot_box_and_whiskers.py` - Show scalar distributions
   - `gcam_plot_spatial_data.py` - Map scalars geographically

3. **GitHub Repository**
   - [E3SM–GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)
   - Additional context, examples, and updates

4. **Scientific Background**
   - [DiVittorio et al. (2025)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024MS004806)
   - Detailed explanation of EHC scalars and their role in coupled modeling

---

## Version Information

**Script:** gcam_compile_ehc_scalars.py  
**Documentation Version:** 1.0  
**Last Updated:** January 2026  
**Python Version:** 3.7+  
**Related Project:** E3SM–GCAM Coupled Modeling Framework

---

## Support and Feedback

For issues, questions, or feedback:

1. **Check this documentation** for configuration examples and troubleshooting
2. **Verify input file formats** match expected EHC scalar structure
3. **Review console output** for specific error messages
4. **Test with small dataset** before processing large compilations
5. **Consult utility modules** (`utility_gcam.py`) for crop name mappings and other details

### Contact Information

**Authors:**
- Philip Myint (myint1@llnl.gov) - Lead developer
- Dalei Hao (dalei.hao@pnnl.gov) - Co-developer

**Contributors:**
- Alan DiVittorio (LBNL)
- Sha Feng (PNNL)
- Eva Sinha (PNNL)

---

## Appendix: Complete Parameter Reference

### Quick Reference Table

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `input_directories` | list of strings | ✓ | - | Paths to directories with scalar CSV files |
| `scenarios` | list of strings | ✓ | - | Names for each scenario (must match directory order) |
| `output_file` | string | ✓ | - | Path to output compiled CSV file |
| `call_modify_crop_names` | boolean | ✗ | `false` | Enable crop name standardization |

### Template Configuration

```json
[
    {
        "input_directories": [
            "./path/to/scenario1/scalars",
            "./path/to/scenario2/scalars",
            "./path/to/scenario3/scalars"
        ],
        "scenarios": [
            "Scenario 1 Name",
            "Scenario 2 Name",
            "Scenario 3 Name"
        ],
        "output_file": "./output/compiled_scalars.csv",
        "call_modify_crop_names": false
    }
]
```

---

This documentation provides comprehensive guidance for using the GCAM EHC scalar compilation script to aggregate and prepare EHC-generated scalar files for analysis in coupled E3SM–GCAM simulations.