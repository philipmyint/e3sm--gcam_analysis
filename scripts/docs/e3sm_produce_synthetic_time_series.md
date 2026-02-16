# E3SM Synthetic Time Series Generation Script Documentation


**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Overview

**Script Name:** `e3sm_produce_synthetic_time_series.py`

**Purpose:** Generates synthetic ensemble members from existing E3SM time series data files by applying systematic and random perturbations. This creates multiple realizations for statistical testing and uncertainty quantification when only single simulation outputs are available.

**Key Difference from Spatial Script:** Uses both systematic (progressive base multipliers) and random (varying across time) perturbations, whereas the spatial script uses only uniform random multipliers.

**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Requirements](#installation--requirements)
3. [Basic Usage](#basic-usage)
4. [How It Works](#how-it-works)
5. [Configuration Parameters](#configuration-parameters)
6. [Perturbation Method](#perturbation-method)
7. [Usage Examples](#usage-examples)
8. [Output Files](#output-files)
9. [Integration with Other Scripts](#integration-with-other-scripts)
10. [Best Practices](#best-practices)
11. [Troubleshooting](#troubleshooting)
12. [Comparison with Spatial Script](#comparison-with-spatial-script)

---

## Installation & Requirements

### Required Python Libraries

```bash
pip install numpy pandas multiprocessing
```

### Required Utility Modules

The script imports:
- `utility_dataframes` - **Required** for reading/writing data files
  - `read_file_into_dataframe()`
  - `write_dataframe_to_fwf()`

### System Requirements

- Python 3.7+
- Multi-core processor (script uses parallel processing)
- Sufficient RAM for DataFrame operations (2+ GB recommended)
- Sufficient disk space for synthetic ensemble files

---

## Basic Usage

### Direct Execution (Hardcoded Configuration)

```bash
python e3sm_produce_synthetic_time_series.py
```

**Note:** This script does NOT use JSON configuration files. All parameters are hardcoded and must be edited directly in the script.

### Editing the Script

Open the script and modify the configuration section:

```python
# Edit these lists in the script:
files = [
    "./control_time_series.dat",
    "./scenario_time_series.dat"
]
num_files_in_each_set = [5, 5]  # Total ensemble members per base file
```

### What the Script Does

1. **Reads** base time series file(s) (.dat or .csv)
2. **Generates** systematic base multipliers (linearly spaced)
3. **Adds** random temporal perturbations
4. **Creates** synthetic copies with combined perturbations
5. **Preserves** Year and Month columns (unperturbed)
6. **Saves** synthetic files with numbered suffixes
7. **Processes** multiple base files in parallel

---

## How It Works

### Dual Perturbation Method

Unlike the spatial script, this uses **two types of perturbations**:

**1. Systematic Base Multipliers (Progressive)**
```python
base_multipliers = np.linspace(1.02, 1.05, num_synthetic_members)
# Example for 4 synthetic members:
# [1.02, 1.03, 1.04, 1.05]
```

**2. Random Temporal Multipliers (Time-Varying)**
```python
random_multipliers = np.random.uniform(low=-0.02, high=0.02, size=num_timesteps)
# Example for each timestep:
# [-0.015, 0.008, -0.003, 0.018, ...]
```

**Combined Perturbation:**
```python
final_multiplier[timestep] = base_multiplier + random_multiplier[timestep]

# Example:
timestep 1: 1.02 + (-0.015) = 1.005 → 0.5% increase
timestep 2: 1.02 + (0.008)  = 1.028 → 2.8% increase
timestep 3: 1.02 + (-0.003) = 1.017 → 1.7% increase
```

### Example Calculation

**Base time series:**
```
Year  Month  GPP (PgC/month)
2015    1         10.0
2015    2         12.0
2015    3         11.5
```

**Synthetic member #2 (base_multiplier = 1.02):**
```python
random_multipliers = [-0.01, 0.015, -0.005]

# Month 1: 10.0 × (1.02 + (-0.01))  = 10.0 × 1.01  = 10.1
# Month 2: 12.0 × (1.02 + 0.015)    = 12.0 × 1.035 = 12.42
# Month 3: 11.5 × (1.02 + (-0.005)) = 11.5 × 1.015 = 11.67
```

**Result:**
```
Year  Month  GPP (PgC/month)
2015    1         10.1
2015    2         12.42
2015    3         11.67
```

---

## Configuration Parameters

### Hardcoded Parameters (Edit in Script)

#### 1. `files`
**Type:** List of strings  
**Description:** Paths to base time series files (.dat or .csv).

**Location in Script:**
```python
files = [
    "./../2025_DiVittorio_et_al_e3sm/control_time_series.dat",
    "./../2025_DiVittorio_et_al_e3sm/full_feedback_time_series.dat",
    "./../2025_DiVittorio_et_al_e3sm/ag_scaling_time_series.dat",
    "./../2025_DiVittorio_et_al_e3sm/carbon_scaling_time_series.dat",
    "./../2025_DiVittorio_et_al_e3sm/control_time_series_amazon.dat"
]
files.extend([
    "./../2025_DiVittorio_et_al_e3sm/control_time_series_surfdata_iESM_dyn_20240730.dat",
    "./../2025_DiVittorio_et_al_e3sm/control_time_series_surfdata_iESM_dyn_20240730_amazon.dat",
    ...
])
```

**File Formats Supported:**
- `.dat` - Fixed-width format (default E3SM output)
- `.csv` - Comma-separated values

**Example Modification:**
```python
# Single file
files = ["./my_time_series.dat"]

# Multiple files
files = [
    "./control_time_series.dat",
    "./scenario_time_series.dat",
    "./control_amazon.csv"
]
```

---

#### 2. `num_files_in_each_set`
**Type:** List of integers  
**Description:** Total number of ensemble members (including base file) for each base file.

**Location in Script:**
```python
num_files_in_each_set = [5]*len(files)  # 5 members for each base file
```

**Interpretation:**
- `num_files_in_each_set = [5]` means 5 total members:
  - 1 original file (base_file.dat)
  - 4 synthetic files (base_file_2.dat, _3.dat, _4.dat, _5.dat)

**Example Modifications:**
```python
# 10 members for each file
num_files_in_each_set = [10]*len(files)

# Different ensemble sizes
num_files_in_each_set = [5, 5, 10, 10]  # First 2 get 5, next 2 get 10

# Single file with 20 members
files = ["./time_series.dat"]
num_files_in_each_set = [20]
```

---

#### 3. Base Multiplier Range (Inside Function)
**Type:** Float range  
**Description:** Range for systematic base multipliers.

**Location in Script:**
```python
def produce_synthetic_time_series(inputs):
    ...
    base_multipliers = np.linspace(1.02, 1.05, num_synthetic_sets_in_ensemble)
```

**Default:** 1.02 to 1.05 (2% to 5% systematic increase)

**Modification Examples:**

**Conservative (1% to 3%):**
```python
base_multipliers = np.linspace(1.01, 1.03, num_synthetic_sets_in_ensemble)
```

**Aggressive (3% to 10%):**
```python
base_multipliers = np.linspace(1.03, 1.10, num_synthetic_sets_in_ensemble)
```

**Symmetric around 1.0:**
```python
base_multipliers = np.linspace(0.98, 1.02, num_synthetic_sets_in_ensemble)
```

---

#### 4. Random Multiplier Range (Inside Function)
**Type:** Float range  
**Description:** Range for random temporal perturbations.

**Location in Script:**
```python
def produce_synthetic_time_series(inputs):
    ...
    random_multipliers = np.random.uniform(low=-0.02, high=0.02, size=len(df))
```

**Default:** -0.02 to +0.02 (±2% random variation)

**Modification Examples:**

**Smaller variation (±1%):**
```python
random_multipliers = np.random.uniform(low=-0.01, high=0.01, size=len(df))
```

**Larger variation (±5%):**
```python
random_multipliers = np.random.uniform(low=-0.05, high=0.05, size=len(df))
```

**Asymmetric (0% to +4%):**
```python
random_multipliers = np.random.uniform(low=0.0, high=0.04, size=len(df))
```

---

## Perturbation Method

### Step-by-Step Process

**For each base file:**

1. **Read base time series**
```python
df = read_file_into_dataframe('control_time_series.dat')
# Contains: Year, Month, GPP, NPP, NBP, etc.
```

2. **Identify data columns**
```python
# Columns to perturb (excludes Year, Month)
columns = [column for column in df.columns if column not in ['Year', 'Month']]
# Example: ['GPP (PgC/month)', 'NPP (PgC/month)', 'NBP (PgC/month)', ...]
```

3. **Generate base multipliers**
```python
# For num_files_in_each_set = 5 (4 synthetic members):
base_multipliers = [1.02, 1.03, 1.04, 1.05]
```

4. **Generate random multipliers**
```python
# One per timestep (e.g., 1032 months for 2015-2100):
random_multipliers = [-0.015, 0.008, -0.003, ..., 0.011]  # 1032 values
```

5. **Apply combined perturbations**
```python
for index in range(4):  # 4 synthetic members
    # Combined multiplier for each timestep
    multipliers = base_multipliers[index] + random_multipliers
    # Example: [1.005, 1.028, 1.017, ..., 1.061]
    
    # Apply to all data columns
    df_new[columns] = df[columns].multiply(multipliers, axis='index')
```

6. **Save synthetic file**
```python
# First synthetic: control_time_series_2.dat
# Second synthetic: control_time_series_3.dat
# ...
```

---

### Mathematical Details

**Full Formula:**
```
synthetic_value[t] = base_value[t] × (base_multiplier + random_multiplier[t])

where:
  t = timestep (month or year)
  base_multiplier ~ linearly spaced in [1.02, 1.05]
  random_multiplier[t] ~ Uniform(-0.02, 0.02)
```

**Properties:**
- **Progressive systematic increase:** Each synthetic member has progressively higher base multiplier
- **Temporal variability:** Random component varies at each timestep
- **Column-wise uniform:** Same multiplier applied to all variables at each timestep
- **Year/Month preserved:** Temporal indices unchanged

**Example with 4 synthetic members:**
```
Synthetic #2: base = 1.02, random varies → overall ~1.00 to 1.04 per timestep
Synthetic #3: base = 1.03, random varies → overall ~1.01 to 1.05 per timestep
Synthetic #4: base = 1.04, random varies → overall ~1.02 to 1.06 per timestep
Synthetic #5: base = 1.05, random varies → overall ~1.03 to 1.07 per timestep
```

---

## Usage Examples

### Example 1: Default Configuration

**As provided in script:**
```python
files = [
    "./../2025_DiVittorio_et_al_e3sm/control_time_series.dat",
    "./../2025_DiVittorio_et_al_e3sm/full_feedback_time_series.dat",
    "./../2025_DiVittorio_et_al_e3sm/ag_scaling_time_series.dat",
    "./../2025_DiVittorio_et_al_e3sm/carbon_scaling_time_series.dat",
    "./../2025_DiVittorio_et_al_e3sm/control_time_series_amazon.dat",
    ...
]
num_files_in_each_set = [5]*len(files)
```

**Execution:**
```bash
python e3sm_produce_synthetic_time_series.py
```

**Output:**
```
Creates 4 synthetic files per base file × 10 base files = 40 synthetic files:
control_time_series_2.dat
control_time_series_3.dat
control_time_series_4.dat
control_time_series_5.dat
full_feedback_time_series_2.dat
...

Total: 10 base + 40 synthetic = 50 files
```

---

### Example 2: Single File, Large Ensemble

**Edit script:**
```python
files = ["./my_time_series.dat"]
num_files_in_each_set = [20]  # 20 total members
```

**Output:**
```
Creates 19 synthetic files:
my_time_series_2.dat
my_time_series_3.dat
...
my_time_series_20.dat

Total ensemble: 1 base + 19 synthetic = 20 members
```

**Base multipliers:** Linearly spaced from 1.02 to 1.05
```
Member 2:  1.02
Member 3:  1.0216
Member 4:  1.0232
...
Member 20: 1.05
```

---

### Example 3: CSV Output

**Edit script:**
```python
files = ["./time_series.csv"]  # Note: .csv extension
num_files_in_each_set = [5]
```

**Output:**
```
Script automatically detects .csv and outputs:
time_series.csv       (base, unchanged)
time_series_2.csv     (synthetic)
time_series_3.csv     (synthetic)
time_series_4.csv     (synthetic)
time_series_5.csv     (synthetic)
```

---

### Example 4: Conservative Perturbations

**Edit script function:**
```python
def produce_synthetic_time_series(inputs):
    start_time = time.time()
    file = inputs[0]
    num_synthetic_sets_in_ensemble = inputs[1] - 1
    df = read_file_into_dataframe(file)
    columns = [column for column in df.columns if column not in ['Year', 'Month']]
    
    # Change these lines:
    base_multipliers = np.linspace(1.005, 1.015, num_synthetic_sets_in_ensemble)  # 0.5% to 1.5%
    random_multipliers = np.random.uniform(low=-0.005, high=0.005, size=len(df))  # ±0.5%
    
    for index in range(len(base_multipliers)):
        ...
```

**Use case:** More conservative uncertainty estimates

---

### Example 5: Symmetric Perturbations

**Edit script function:**
```python
# Symmetric around 1.0 (some decrease, some increase)
base_multipliers = np.linspace(0.98, 1.02, num_synthetic_sets_in_ensemble)  # -2% to +2%
random_multipliers = np.random.uniform(low=-0.01, high=0.01, size=len(df))  # ±1%
```

**Result:** Some synthetic members have lower values than base

---

## Output Files

### File Naming Convention

**Pattern:**
```
base_filename.dat           → Original (untouched)
base_filename_2.dat         → First synthetic member
base_filename_3.dat         → Second synthetic member
...
base_filename_N.dat         → (N-1)th synthetic member
```

**For CSV files:**
```
base_filename.csv           → Original (untouched)
base_filename_2.csv         → First synthetic member
...
```

**Example:**
```
Input: control_time_series.dat
Output:
  control_time_series.dat     (original, unchanged)
  control_time_series_2.dat   (synthetic, base × 1.02 + random)
  control_time_series_3.dat   (synthetic, base × 1.03 + random)
  control_time_series_4.dat   (synthetic, base × 1.04 + random)
  control_time_series_5.dat   (synthetic, base × 1.05 + random)
```

---

### File Structure

Each synthetic file has identical structure to base file:

**Fixed-Width Format (.dat):**
```
Year  Month  GPP (PgC/month)  NPP (PgC/month)  NBP (PgC/month)
2015      1        10.1234          8.5678          1.2345
2015      2        12.3456         10.1234          2.1234
2015      3        11.5678          9.8765          1.6789
...
```

**CSV Format (.csv):**
```
Year,Month,GPP (PgC/month),NPP (PgC/month),NBP (PgC/month)
2015,1,10.1234,8.5678,1.2345
2015,2,12.3456,10.1234,2.1234
2015,3,11.5678,9.8765,1.6789
...
```

**Preserved:**
- Year and Month columns (unchanged)
- Column headers
- File format

**Modified:**
- All data columns (perturbed by combined multipliers)

---

### File Sizes

**Typical sizes:**
- Base file: 100 KB (86 years × 12 months × 10 variables)
- Each synthetic file: 100 KB
- 5-member ensemble: ~500 KB (5 × 100 KB)
- 10-member ensemble: ~1 MB (10 × 100 KB)

**Storage requirement:**
```
Total storage = base_file_size × num_files_in_each_set × num_base_files
```

**Example:**
```
10 base files × 100 KB each = 1 MB
5 members each × 10 files = 50 total files
Total storage: 50 × 100 KB = 5 MB
```

---

## Integration with Other Scripts

### Typical Workflow

#### 1. Extract Time Series from E3SM

```bash
# Extract h0 time series
python e3sm_extract_time_series_h0.py config.json
# Output: control_time_series.dat

# Or extract surfdata time series
python e3sm_extract_time_series_surfdata_iesm_dyn.py config.json
# Output: control_time_series_surfdata.dat
```

#### 2. Create Synthetic Ensemble

```bash
# Generate synthetic members
python e3sm_produce_synthetic_time_series.py

# Output:
#   control_time_series.dat (base)
#   control_time_series_2.dat (synthetic)
#   control_time_series_3.dat (synthetic)
#   ...
```

#### 3. Plot with Statistical Analysis

```bash
# Use ensemble for statistical testing
python gcam_plot_time_series.py plot_config.json

# plot_config.json specifies ensemble members for comparison
```

---

### Use with Time Series Plotting Script

**Enable ensemble mean and spread:**

```json
{
    "files": [
        "control_time_series.dat",
        "control_time_series_2.dat",
        "control_time_series_3.dat",
        "control_time_series_4.dat",
        "control_time_series_5.dat"
    ],
    "variable": "GPP (PgC/month)",
    "plot_ensemble_mean": true,
    "plot_ensemble_spread": true
}
```

**Result:** Plot shows ensemble mean with shaded spread (min-max or std deviation)

---

## Best Practices

### 1. Ensemble Size Selection

**Small ensembles (n=3-5):**
- **Pros:** Fast, low storage, good for initial exploration
- **Cons:** Limited statistical power
- **Use:** Quick analysis, method development

**Medium ensembles (n=5-10):**
- **Pros:** Good statistical power, reasonable storage
- **Cons:** Moderate computation
- **Use:** Standard analysis (recommended)

**Large ensembles (n=10-20):**
- **Pros:** High statistical power, robust estimates
- **Cons:** Storage requirements, processing time
- **Use:** Publication-quality analysis

**Recommendation:** Start with n=5

---

### 2. Perturbation Range Selection

**Base Multiplier Range:**

**Conservative (0.5% to 1.5%):**
```python
base_multipliers = np.linspace(1.005, 1.015, ...)
```
- Tight systematic increase
- Use when model output highly certain

**Standard (2% to 5%):**
```python
base_multipliers = np.linspace(1.02, 1.05, ...)  # Default
```
- Moderate systematic increase
- General use case

**Wide (5% to 10%):**
```python
base_multipliers = np.linspace(1.05, 1.10, ...)
```
- Large systematic increase
- Sensitivity testing

**Random Multiplier Range:**

**Conservative (±0.5%):**
```python
random_multipliers = np.random.uniform(low=-0.005, high=0.005, ...)
```

**Standard (±2%):**
```python
random_multipliers = np.random.uniform(low=-0.02, high=0.02, ...)  # Default
```

**Wide (±5%):**
```python
random_multipliers = np.random.uniform(low=-0.05, high=0.05, ...)
```

---

### 3. File Organization

**Recommended structure:**
```
project/
├── base_time_series/
│   ├── control_time_series.dat
│   └── scenario_time_series.dat
├── synthetic_ensembles/
│   ├── control_time_series_2.dat
│   ├── control_time_series_3.dat
│   └── ...
└── plots/
```

**Or keep together:**
```
project/
├── time_series/
│   ├── control_time_series.dat       (base)
│   ├── control_time_series_2.dat     (synthetic)
│   ├── control_time_series_3.dat     (synthetic)
│   └── ...
```

---

### 4. Reproducibility

**Set random seed:**
```python
def produce_synthetic_time_series(inputs):
    # Add at beginning:
    np.random.seed(42)  # Fixed seed
    
    start_time = time.time()
    file = inputs[0]
    ...
```

**Caution:** Same seed = same synthetic ensemble every run

**Alternative:** Save random values
```python
# After generating:
np.save('random_multipliers.npy', random_multipliers)
```

---

### 5. Validation

**After generation, verify:**

```python
import pandas as pd
import numpy as np

# Load base and synthetic
base = pd.read_fwf('control_time_series.dat')
synthetic = pd.read_fwf('control_time_series_2.dat')

# Calculate ratio for data columns
data_cols = [col for col in base.columns if col not in ['Year', 'Month']]
ratio = synthetic[data_cols] / base[data_cols]

print(f"Mean ratio: {ratio.mean().mean():.4f}")
print(f"Min ratio: {ratio.min().min():.4f}")
print(f"Max ratio: {ratio.max().max():.4f}")

# For synthetic #2 with base=1.02, random±0.02:
# Expected mean: ~1.02
# Expected min: ~1.00
# Expected max: ~1.04
```

---

## Troubleshooting

### Issue 1: Memory Error

**Error:**
```
MemoryError: Unable to allocate array
```

**Cause:**
- Large time series files
- Many files processed simultaneously

**Solutions:**

1. **Reduce parallel processing:**
```python
# Instead of all CPUs:
with multiprocessing.Pool(processes=4) as pool:  # Use only 4
    pool.map(produce_synthetic_time_series, inputs)
```

2. **Process sequentially:**
```python
# Replace parallel loop:
for inp in inputs:
    produce_synthetic_time_series(inp)
```

3. **Reduce ensemble size:**
```python
num_files_in_each_set = [3]*len(files)  # 3 instead of 5
```

---

### Issue 2: File Not Found

**Error:**
```
FileNotFoundError: [Errno 2] No such file or directory
```

**Cause:**
- Incorrect file path in `files` list
- Base file doesn't exist

**Solutions:**
```bash
# Verify file exists
ls -la control_time_series.dat

# Use absolute paths
files = ["/absolute/path/to/control_time_series.dat"]

# Or correct relative path
files = ["./output/control_time_series.dat"]
```

---

### Issue 3: Year/Month Column Issues

**Error:**
```
KeyError: 'Year' or KeyError: 'Month'
```

**Cause:**
- Time series file doesn't have Year/Month columns
- Column names different (e.g., "year" vs "Year")

**Solutions:**

**Check columns:**
```python
import pandas as pd
df = pd.read_fwf('time_series.dat')
print(df.columns)
```

**Modify script if needed:**
```python
# If columns are 'year' and 'month' (lowercase):
columns = [column for column in df.columns if column not in ['year', 'month']]
```

---

### Issue 4: All Values Same

**Symptom:** All synthetic files identical

**Cause:**
- Random seed set elsewhere in code
- numpy random not seeded properly

**Solution:**
```python
# Ensure numpy random is working:
import numpy as np
print(np.random.uniform(low=-0.02, high=0.02, size=5))
# Should show different values each time (unless seed set)
```

---

### Issue 5: Negative Values

**Symptom:** Some values become negative

**Cause:**
- Original values close to zero
- Perturbations too large
- Combined perturbation < 0

**Example:**
```
Original NBP = 0.1 PgC/month
Multiplier = 1.02 + (-0.03) = 0.99
Result = 0.1 × 0.99 = 0.099  (OK)

But if:
Original NBP = 0.05 PgC/month
Multiplier = 1.02 + (-0.04) = 0.98  (random can go to -0.02, base 1.02 → 1.00)
Result still positive

However, if values very close to zero, issues can arise
```

**Solutions:**

1. **Reduce random range:**
```python
random_multipliers = np.random.uniform(low=-0.01, high=0.01, ...)  # ±1% instead of ±2%
```

2. **Ensure multipliers stay positive:**
```python
# Add clipping:
multipliers = np.maximum(base_multipliers[index] + random_multipliers, 0.5)  # Min 50%
```

3. **Check original data:**
```python
# Identify columns with near-zero values
print(df[columns].min())
```

---

## Comparison with Spatial Script

### Key Differences

| Feature | Time Series Script | Spatial Script |
|---------|-------------------|----------------|
| **Perturbation Type** | Systematic + Random | Random only |
| **Base Multipliers** | Linearly spaced (1.02-1.05) | Uniform random (0.99-1.02) |
| **Random Component** | Temporal (-0.02 to +0.02) | Spatially uniform |
| **Variation** | Across time | Across space |
| **File Format** | Text (.dat, .csv) | NetCDF (.nc) |
| **Data Type** | Time series (1D/2D) | Spatial grids (3D) |
| **Preserved Columns** | Year, Month | Coordinates (lat, lon, year) |

### Perturbation Comparison

**Time Series Script:**
```python
# Progressive + temporal variability
base_multipliers = [1.02, 1.03, 1.04, 1.05]  # Systematic
random_multipliers = [-0.015, 0.008, ...]     # Per timestep

Synthetic #2: 1.02 + random[t]  → varies 1.00-1.04 over time
Synthetic #3: 1.03 + random[t]  → varies 1.01-1.05 over time
```

**Spatial Script:**
```python
# Uniform random only
random_multipliers = [1.015, 0.992, 1.008, 0.997]  # Per member

Synthetic #2: 1.015 × base  → uniform 1.5% increase everywhere
Synthetic #3: 0.992 × base  → uniform 0.8% decrease everywhere
```

---

### When to Use Each Script

**Use Time Series Script when:**
- Working with time series data (.dat, .csv files)
- Need temporal variability in synthetic members
- Want progressive systematic increases across ensemble
- Analyzing trends, time evolution
- Input from time series extraction scripts

**Use Spatial Script when:**
- Working with spatial data (NetCDF files)
- Need spatial variability representation
- Want uniform perturbations across space
- Analyzing spatial patterns, maps
- Input from spatial data extraction script

---

## Scientific Considerations

### Interpretation of Perturbations

**Systematic Component (Base Multipliers):**
- Represents progressive uncertainty levels
- Member #2 (1.02) = conservative estimate
- Member #5 (1.05) = less conservative estimate
- Allows exploring range of possible outcomes

**Random Component:**
- Represents temporal variability/uncertainty
- Different from interannual variability (which is in base data)
- Adds "noise" to time series
- Not physically based

**Combined:**
- Creates ensemble with both systematic spread and temporal variability
- Useful for statistical testing
- Not substitute for true model ensemble

---

### Limitations

**1. Not True Ensemble Members:**
- Don't represent different initial conditions
- Don't capture true model physics uncertainty
- Synthetic, not physically based

**2. Artificial Perturbations:**
- Linear systematic increase not physically motivated
- Random temporal component arbitrary
- Doesn't represent real climate variability

**3. Statistical Limitations:**
- May not represent true uncertainty distribution
- P-values approximate
- Ensemble spread artificial

---

### Appropriate Use Cases

**Good for:**
- Exploratory statistical testing
- Method development
- Sensitivity analysis
- When true ensemble unavailable

**Not appropriate for:**
- Representing true climate variability
- Publishing as actual ensemble
- Quantifying model uncertainty (use true ensemble)
- Attribution studies

---

## Advanced Modifications

### Normal Distribution Random Perturbations

**Current:** Uniform distribution  
**Modification:** Normal (Gaussian) distribution

```python
# In function, replace:
random_multipliers = np.random.uniform(low=-0.02, high=0.02, size=len(df))

# With:
random_multipliers = np.random.normal(loc=0.0, scale=0.01, size=len(df))
# loc = mean (0.0 = centered on base)
# scale = standard deviation (0.01 = 1% std)
```

---

### Variable-Specific Perturbations

**Current:** Same perturbation for all variables  
**Modification:** Different perturbation per variable

```python
def produce_synthetic_time_series(inputs):
    start_time = time.time()
    file = inputs[0]
    num_synthetic_sets_in_ensemble = inputs[1] - 1
    df = read_file_into_dataframe(file)
    columns = [column for column in df.columns if column not in ['Year', 'Month']]
    base_multipliers = np.linspace(1.02, 1.05, num_synthetic_sets_in_ensemble)
    
    for index in range(len(base_multipliers)):
        df_new = df.copy()
        
        # Different random perturbation for each column
        for column in columns:
            random_multipliers = np.random.uniform(low=-0.02, high=0.02, size=len(df))
            multipliers = base_multipliers[index] + random_multipliers
            df_new[column] = df[column] * multipliers
        
        # Save file...
```

---

### Autocorrelated Random Perturbations

**Current:** Independent random at each timestep  
**Modification:** Temporally correlated (more realistic)

```python
# Add autocorrelation:
def generate_autocorrelated_noise(n, alpha=0.5):
    """Generate autocorrelated noise with AR(1) process"""
    noise = np.zeros(n)
    noise[0] = np.random.normal(0, 0.01)
    for i in range(1, n):
        noise[i] = alpha * noise[i-1] + np.random.normal(0, 0.01)
    return noise

# In function:
random_multipliers = generate_autocorrelated_noise(len(df), alpha=0.7)
```

---

## Appendix: Complete Modified Example

### Custom Script for 3 Files, Different Settings

```python
import multiprocessing
import numpy as np
import time
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_fwf

def produce_synthetic_time_series(inputs):
    start_time = time.time()
    file = inputs[0]
    num_synthetic_sets_in_ensemble = inputs[1] - 1
    df = read_file_into_dataframe(file)
    columns = [column for column in df.columns if column not in ['Year', 'Month']]
    
    # Custom: Conservative perturbations
    base_multipliers = np.linspace(1.005, 1.015, num_synthetic_sets_in_ensemble)
    random_multipliers = np.random.uniform(low=-0.005, high=0.005, size=len(df))
    
    for index in range(len(base_multipliers)):
        df_new = df
        multipliers = base_multipliers[index] + random_multipliers
        df_new[columns] = df[columns].multiply(multipliers, axis='index')
        
        if file.endswith('.csv'):
            new_file = file.replace('.csv', f'_{index+2}.csv')
            df_new.to_csv(new_file, index=False)
        else:
            new_file = file.replace('.dat', f'_{index+2}.dat')
            write_dataframe_to_fwf(new_file, df_new)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for {file}: {elapsed_time:.2f} seconds")

if __name__ == '__main__':
    start_time = time.time()
    
    # Configuration
    files = [
        "./control_time_series.dat",
        "./scenario_time_series.dat",
        "./regional_time_series.csv"
    ]
    num_files_in_each_set = [10, 10, 5]  # Different sizes
    
    inputs = list(zip(files, num_files_in_each_set))
    
    with multiprocessing.Pool(processes=MAX_PROCESSES) as pool:
        pool.map(produce_synthetic_time_series, inputs)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total time: {elapsed_time:.2f} seconds")
```

---

## References

- E3SM Documentation: [https://e3sm.org/](https://e3sm.org/)
- NumPy Documentation: [https://numpy.org/doc/stable/](https://numpy.org/doc/stable/)
- Pandas Documentation: [https://pandas.pydata.org/](https://pandas.pydata.org/)
- GitHub Repository: [https://github.com/philipmyint/e3sm--gcam_analysis](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Contact

For questions or issues:
- **Philip Myint**: myint1@llnl.gov
- **Dalei Hao**: dalei.hao@pnnl.gov

---

## Version Information

**Script:** e3sm_produce_synthetic_time_series.py  
**Documentation Version:** 1.0  
**Last Updated:** January 2026  
**Python Version:** 3.7+  
**Dependencies:** numpy, pandas, multiprocessing, utility_dataframes

---

*This documentation provides comprehensive guidance for using the `e3sm_produce_synthetic_time_series.py` script to generate synthetic ensemble members from E3SM time series data for statistical testing and uncertainty quantification.*
