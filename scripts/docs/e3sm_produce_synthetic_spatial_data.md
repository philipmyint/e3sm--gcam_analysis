# E3SM Synthetic Spatial Data Generation Script Documentation


**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Overview

**Script Name:** `e3sm_produce_synthetic_spatial_data.py`

**Purpose:** Generates synthetic ensemble members from existing E3SM spatial NetCDF files by applying random perturbations to create multiple realizations for uncertainty quantification and ensemble analysis. This is particularly useful for statistical significance testing, sensitivity analysis, and creating pseudo-ensembles when limited simulation data is available.

**Key Use Case:** Creating synthetic ensemble members to enable statistical testing (e.g., t-tests, ensemble means) when only single simulation outputs exist. Common in climate model analysis for assessing spatial variability and statistical significance of differences between scenarios.

**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Requirements](#installation--requirements)
3. [Basic Usage](#basic-usage)
4. [How It Works](#how-it-works)
5. [Configuration Parameters](#configuration-parameters)
6. [Synthetic Data Generation Method](#synthetic-data-generation-method)
7. [Usage Examples](#usage-examples)
8. [Output Files](#output-files)
9. [Integration with Other Scripts](#integration-with-other-scripts)
10. [Best Practices](#best-practices)
11. [Troubleshooting](#troubleshooting)
12. [Scientific Considerations](#scientific-considerations)

---

## Installation & Requirements

### Required Python Libraries

```bash
pip install numpy xarray netCDF4 multiprocessing
```

### Required Utility Modules

The script imports utility modules (though none are strictly required for this script):
- `utility_e3sm_netcdf` (not used in current implementation)
- `utility_functions` (not used in current implementation)
- `utility_dataframes` (not used in current implementation)
- `utility_xarray` (not used in current implementation)

**Note:** These imports can be removed if not needed.

### System Requirements

- Python 3.7+
- Multi-core processor (script uses parallel processing)
- Sufficient RAM to load NetCDF files (4+ GB recommended)
- Sufficient disk space for synthetic ensemble files

---

## Basic Usage

### Direct Execution (Hardcoded Configuration)

```bash
python e3sm_produce_synthetic_spatial_data.py
```

**Note:** This script does NOT use JSON configuration files. All parameters are hardcoded in the script and must be edited directly.

### Editing the Script

Open the script and modify the configuration section:

```python
# Edit these lists in the script:
files = [
    "./path/to/spatial_data_elm.nc",
    "./path/to/spatial_data_eam.nc"
]
num_files_in_each_set = [5, 5]  # Total ensemble members per base file
```

### What the Script Does

1. **Reads** base spatial NetCDF file(s)
2. **Generates** random multipliers (default: uniform 0.99-1.02)
3. **Creates** synthetic copies by multiplying all variables by random factor
4. **Saves** synthetic files with numbered suffixes (_2.nc, _3.nc, etc.)
5. **Processes** multiple base files in parallel

---

## How It Works

### Perturbation Method

The script creates synthetic ensemble members by applying random multiplicative perturbations:

```python
# For each base file:
base_file.nc → Read original data

# Create N-1 synthetic members (N = total ensemble size):
synthetic_1 = base_data × random_multiplier_1  # e.g., × 1.015
synthetic_2 = base_data × random_multiplier_2  # e.g., × 0.997
synthetic_3 = base_data × random_multiplier_3  # e.g., × 1.008
...

# Save as:
base_file_2.nc
base_file_3.nc
base_file_4.nc
...
```

### Random Multiplier Generation

**Distribution:** Uniform random distribution

**Range:** 0.99 to 1.02 (default)

**Interpretation:**
- Multiplier = 0.99 → 1% decrease
- Multiplier = 1.00 → No change
- Multiplier = 1.02 → 2% increase

**Example:**
```python
# If base GPP = 10.0 gC/m²/s at a gridcell:
multiplier = 1.015  # Randomly drawn
synthetic_GPP = 10.0 × 1.015 = 10.15 gC/m²/s
```

### Parallel Processing

The script uses Python's multiprocessing to create synthetic ensembles in parallel:

```python
# Process multiple base files simultaneously
with multiprocessing.Pool(processes=MAX_PROCESSES) as pool:
    pool.map(produce_synthetic_spatial_data, inputs)
```

**Efficiency:**
- 6 base files × 4 synthetic members each = 24 files
- On 8-core machine: Processes ~8 files simultaneously
- Typical speedup: 5-8x compared to sequential

---

## Configuration Parameters

### Hardcoded Parameters (Edit in Script)

#### 1. `files`
**Type:** List of strings  
**Description:** Paths to base spatial NetCDF files from which to generate synthetic ensembles.

**Location in Script:**
```python
files = [
    "./../2025_DiVittorio_et_al_e3sm/control_spatial_data_elm.nc",
    "./../2025_DiVittorio_et_al_e3sm/full_feedback_spatial_data_elm.nc",
    "./../2025_DiVittorio_et_al_e3sm/carbon_scaling_spatial_data_elm.nc"
]
files.extend([
    "./../2025_DiVittorio_et_al_e3sm/control_spatial_data_eam.nc",
    "./../2025_DiVittorio_et_al_e3sm/full_feedback_spatial_data_eam.nc",
    "./../2025_DiVittorio_et_al_e3sm/carbon_scaling_spatial_data_eam.nc"
])
```

**Example Modification:**
```python
# Single file
files = ["./my_spatial_data.nc"]

# Multiple files
files = [
    "./control_elm.nc",
    "./control_eam.nc",
    "./scenario_elm.nc",
    "./scenario_eam.nc"
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
  - 1 original file (base_file.nc)
  - 4 synthetic files (base_file_2.nc, _3.nc, _4.nc, _5.nc)

**Example Modifications:**
```python
# 10 members for each file
num_files_in_each_set = [10]*len(files)

# Different ensemble sizes for different files
num_files_in_each_set = [5, 5, 10, 10]  # First 2 get 5, next 2 get 10

# Single file with 20 members
files = ["./spatial_data.nc"]
num_files_in_each_set = [20]
```

---

#### 3. Random Multiplier Range (Inside Function)
**Type:** Float range  
**Description:** Range for random perturbations.

**Location in Script:**
```python
def produce_synthetic_spatial_data(inputs):
    ...
    random_multipliers = np.random.uniform(low=0.99, high=1.02, size=num_synthetic_sets_in_ensemble)
```

**Default:** 0.99 to 1.02 (±1-2% variation)

**Modification Examples:**

**Smaller perturbations (±0.5%):**
```python
random_multipliers = np.random.uniform(low=0.995, high=1.005, size=num_synthetic_sets_in_ensemble)
```

**Larger perturbations (±5%):**
```python
random_multipliers = np.random.uniform(low=0.95, high=1.05, size=num_synthetic_sets_in_ensemble)
```

**Conservative (±1%):**
```python
random_multipliers = np.random.uniform(low=0.99, high=1.01, size=num_synthetic_sets_in_ensemble)
```

---

## Synthetic Data Generation Method

### Step-by-Step Process

**For each base file:**

1. **Generate random multipliers**
```python
# If num_files_in_each_set = 5, create 4 multipliers
multipliers = [1.015, 0.992, 1.008, 0.997]  # Example random values
```

2. **Load base NetCDF file**
```python
ds = xr.open_dataset('control_spatial_data_elm.nc')
# Contains: GPP, NPP, NBP, TOTECOSYSC, etc.
```

3. **Apply perturbation to ALL variables**
```python
# For each synthetic member:
for variable in ['GPP', 'NPP', 'NBP', 'TOTECOSYSC', ...]:
    ds[variable] = ds[variable] × multiplier
```

4. **Save synthetic file**
```python
# First synthetic: control_spatial_data_elm_2.nc
# Second synthetic: control_spatial_data_elm_3.nc
# ...
```

5. **Repeat** for all ensemble members

---

### Mathematical Details

**Perturbation Formula:**
```
synthetic_value(lat, lon, year) = base_value(lat, lon, year) × multiplier

where:
  multiplier ~ Uniform(0.99, 1.02)
```

**Properties:**
- **Spatially uniform:** Same multiplier applied to all gridcells
- **Variable-specific:** Different multiplier for each ensemble member
- **Temporally uniform:** Same multiplier for all years
- **Preserves structure:** Spatial patterns maintained, only magnitudes change

**Example:**
```
Base GPP at gridcell (45°N, 120°W) = 12.5 gC/m²/s
Multiplier = 1.015
Synthetic GPP = 12.5 × 1.015 = 12.6875 gC/m²/s

Base GPP at gridcell (30°N, 90°W) = 8.3 gC/m²/s
Same multiplier = 1.015
Synthetic GPP = 8.3 × 1.015 = 8.4245 gC/m²/s
```

---

## Usage Examples

### Example 1: Default Configuration

**As provided in script:**
```python
files = [
    "./../2025_DiVittorio_et_al_e3sm/control_spatial_data_elm.nc",
    "./../2025_DiVittorio_et_al_e3sm/full_feedback_spatial_data_elm.nc",
    "./../2025_DiVittorio_et_al_e3sm/carbon_scaling_spatial_data_elm.nc",
    "./../2025_DiVittorio_et_al_e3sm/control_spatial_data_eam.nc",
    "./../2025_DiVittorio_et_al_e3sm/full_feedback_spatial_data_eam.nc",
    "./../2025_DiVittorio_et_al_e3sm/carbon_scaling_spatial_data_eam.nc"
]
num_files_in_each_set = [5]*len(files)
```

**Execution:**
```bash
python e3sm_produce_synthetic_spatial_data.py
```

**Output:**
```
Creates 24 synthetic files (4 per base file × 6 base files):
control_spatial_data_elm_2.nc
control_spatial_data_elm_3.nc
control_spatial_data_elm_4.nc
control_spatial_data_elm_5.nc
full_feedback_spatial_data_elm_2.nc
full_feedback_spatial_data_elm_3.nc
...
(Total: 6 base + 24 synthetic = 30 files)
```

---

### Example 2: Single File, Large Ensemble

**Edit script:**
```python
files = ["./my_spatial_data.nc"]
num_files_in_each_set = [20]  # 20 total members
```

**Output:**
```
Creates 19 synthetic files:
my_spatial_data_2.nc
my_spatial_data_3.nc
...
my_spatial_data_20.nc

Total ensemble: 1 base + 19 synthetic = 20 members
```

**Use case:** Large ensemble for robust statistical testing

---

### Example 3: Different Ensemble Sizes

**Edit script:**
```python
files = [
    "./control_elm.nc",
    "./control_eam.nc",
    "./scenario_elm.nc",
    "./scenario_eam.nc"
]
num_files_in_each_set = [10, 10, 5, 5]
```

**Output:**
```
control_elm: 9 synthetic files (_2 through _10)
control_eam: 9 synthetic files (_2 through _10)
scenario_elm: 4 synthetic files (_2 through _5)
scenario_eam: 4 synthetic files (_2 through _5)

Total: 4 base + 26 synthetic = 30 files
```

---

### Example 4: Conservative Perturbations

**Edit script function:**
```python
def produce_synthetic_spatial_data(inputs):
    start_time = time.time()
    file = inputs[0]
    num_synthetic_sets_in_ensemble = inputs[1] - 1
    
    # Change this line:
    random_multipliers = np.random.uniform(low=0.995, high=1.005, size=num_synthetic_sets_in_ensemble)
    # Now ±0.5% instead of ±2%
    
    for index in range(len(random_multipliers)):
        ...
```

**Use case:** More conservative uncertainty estimates

---

## Output Files

### File Naming Convention

**Pattern:**
```
base_filename.nc           → Original (untouched)
base_filename_2.nc         → First synthetic member
base_filename_3.nc         → Second synthetic member
...
base_filename_N.nc         → (N-1)th synthetic member
```

**Example:**
```
Input: control_spatial_data_elm.nc
Output:
  control_spatial_data_elm.nc     (original, unchanged)
  control_spatial_data_elm_2.nc   (synthetic)
  control_spatial_data_elm_3.nc   (synthetic)
  control_spatial_data_elm_4.nc   (synthetic)
  control_spatial_data_elm_5.nc   (synthetic)
```

---

### File Structure

Each synthetic file is identical in structure to the base file, only values differ:

**Dimensions:** Same as base file
```
year: 6
lat: 192
lon: 288
```

**Variables:** Same as base file
```
GPP(year, lat, lon)
NPP(year, lat, lon)
NBP(year, lat, lon)
...
```

**Coordinates:** Identical to base file
```
lat, lon, year values unchanged
```

**Metadata:** Copied from base file

**Data Values:** Multiplied by random factor

---

### File Sizes

**Typical sizes:**
- Base file: 100 MB
- Each synthetic file: 100 MB
- 5-member ensemble: ~500 MB (5 × 100 MB)
- 10-member ensemble: ~1 GB (10 × 100 MB)

**Storage requirement:**
```
Total storage = base_file_size × num_files_in_each_set × num_base_files
```

**Example:**
```
6 base files × 100 MB each = 600 MB
5 members each × 6 files = 30 total files
Total storage: 30 × 100 MB = 3 GB
```

---

## Integration with Other Scripts

### Typical Workflow

#### 1. Generate Base Spatial Data

```bash
# Extract spatial data from E3SM simulation
python e3sm_extract_spatial_data_h0.py config.json

# Output: control_spatial_data_elm.nc
```

#### 2. Create Synthetic Ensemble

```bash
# Generate synthetic members
python e3sm_produce_synthetic_spatial_data.py

# Output: 
#   control_spatial_data_elm.nc (base)
#   control_spatial_data_elm_2.nc (synthetic)
#   control_spatial_data_elm_3.nc (synthetic)
#   ...
```

#### 3. Plot with Statistical Significance

```bash
# Use ensemble for t-tests in spatial plotting
python gcam_plot_spatial_data.py plot_config.json

# plot_config.json specifies ensemble members:
{
    "scenarios": [
        ["control_spatial_data_elm.nc",
         "control_spatial_data_elm_2.nc",
         "control_spatial_data_elm_3.nc",
         "control_spatial_data_elm_4.nc",
         "control_spatial_data_elm_5.nc"]
    ]
}
```

---

### Use with Spatial Plotting Script

**Enable statistical significance testing:**

```json
{
    "scenarios": [
        ["control_elm.nc", "control_elm_2.nc", "control_elm_3.nc", 
         "control_elm_4.nc", "control_elm_5.nc"],
        ["scenario_elm.nc", "scenario_elm_2.nc", "scenario_elm_3.nc",
         "scenario_elm_4.nc", "scenario_elm_5.nc"]
    ],
    "plot_type": "absolute_difference",
    "stippling_on": true,
    "p_value_threshold": 0.05
}
```

**Result:** Spatial map showing difference with stippling where statistically significant (p < 0.05)

---

## Best Practices

### 1. Ensemble Size Selection

**Small ensembles (n=3-5):**
- **Pros:** Fast, low storage
- **Cons:** Limited statistical power
- **Use:** Quick exploratory analysis

**Medium ensembles (n=5-10):**
- **Pros:** Good statistical power, reasonable storage
- **Cons:** Moderate computation time
- **Use:** Standard analysis (recommended)

**Large ensembles (n=10-20):**
- **Pros:** High statistical power, robust estimates
- **Cons:** High storage, longer processing
- **Use:** Publication-quality statistical testing

**Recommendation:** Start with n=5, increase if needed

---

### 2. Perturbation Range Selection

**Conservative (±0.5%):**
```python
random_multipliers = np.random.uniform(low=0.995, high=1.005, ...)
```
- Tight uncertainty bounds
- Use when model output is highly certain

**Standard (±2%):**
```python
random_multipliers = np.random.uniform(low=0.99, high=1.02, ...)  # Default
```
- Moderate uncertainty
- General use case

**Wide (±5%):**
```python
random_multipliers = np.random.uniform(low=0.95, high=1.05, ...)
```
- Large uncertainty bounds
- Use for sensitivity testing

**Guidance:** Match to expected model uncertainty

---

### 3. File Organization

**Recommended structure:**
```
project/
├── base_spatial_data/
│   ├── control_spatial_data_elm.nc
│   └── control_spatial_data_eam.nc
├── synthetic_ensembles/
│   ├── control_spatial_data_elm_2.nc
│   ├── control_spatial_data_elm_3.nc
│   ├── ...
│   ├── control_spatial_data_eam_2.nc
│   └── ...
└── plots/
```

**Or keep together:**
```
project/
├── spatial_data/
│   ├── control_spatial_data_elm.nc       (base)
│   ├── control_spatial_data_elm_2.nc     (synthetic)
│   ├── control_spatial_data_elm_3.nc     (synthetic)
│   └── ...
```

---

### 4. Version Control

**Document your configuration:**

```python
# Create README.txt
"""
Synthetic ensemble created: 2026-01-24
Base files: control_spatial_data_elm.nc, control_spatial_data_eam.nc
Ensemble size: 5 members per file
Perturbation range: 0.99-1.02 (±2%)
Random seed: Not set (different each run)
Purpose: Statistical significance testing for Figure 3
"""
```

---

### 5. Reproducibility

**Set random seed for reproducibility:**

```python
def produce_synthetic_spatial_data(inputs):
    # Add this at the beginning:
    np.random.seed(42)  # Fixed seed for reproducibility
    
    start_time = time.time()
    file = inputs[0]
    ...
```

**Caution:** Setting seed means same synthetic ensemble every time

**Alternative:** Save random multipliers to file

```python
# After generating multipliers:
np.save('random_multipliers.npy', random_multipliers)

# Document in README:
# Random multipliers saved to random_multipliers.npy
```

---

### 6. Validation

**After generation, verify:**

```python
import xarray as xr
import numpy as np

# Load base and synthetic
base = xr.open_dataset('control_elm.nc')
synthetic = xr.open_dataset('control_elm_2.nc')

# Calculate ratio (should be close to 1.0)
ratio = synthetic['GPP'] / base['GPP']

print(f"Min ratio: {ratio.min().item():.4f}")
print(f"Max ratio: {ratio.max().item():.4f}")
print(f"Mean ratio: {ratio.mean().item():.4f}")

# Should be approximately:
# Min: ~0.99
# Max: ~1.02
# Mean: ~1.005
```

---

## Troubleshooting

### Issue 1: Memory Error

**Error:**
```
MemoryError: Unable to allocate array
```

**Cause:**
- Large NetCDF files
- Many files processed simultaneously
- Insufficient RAM

**Solutions:**

1. **Reduce parallel processing:**
```python
# Instead of all CPUs:
with multiprocessing.Pool(processes=4) as pool:  # Use only 4 cores
    pool.map(produce_synthetic_spatial_data, inputs)
```

2. **Process files sequentially:**
```python
# Replace parallel processing with loop:
for inp in inputs:
    produce_synthetic_spatial_data(inp)
```

3. **Reduce ensemble size:**
```python
num_files_in_each_set = [3]*len(files)  # 3 instead of 5
```

---

### Issue 2: File Not Found

**Error:**
```
FileNotFoundError: [Errno 2] No such file or directory: 'control_spatial_data_elm.nc'
```

**Cause:**
- Incorrect file path in `files` list
- Base file doesn't exist

**Solutions:**

```bash
# Verify file exists
ls -la control_spatial_data_elm.nc

# Use absolute paths
files = ["/absolute/path/to/control_spatial_data_elm.nc"]

# Or correct relative path
files = ["./output/control_spatial_data_elm.nc"]
```

---

### Issue 3: Disk Space

**Error:**
```
OSError: [Errno 28] No space left on device
```

**Cause:**
- Insufficient disk space for synthetic files
- Large ensembles × large base files

**Solutions:**

1. **Check available space:**
```bash
df -h .  # Check disk space
```

2. **Estimate required space:**
```python
# If base file is 100 MB and you want 5 members:
# Required space ≈ 100 MB × 5 members × num_files
```

3. **Reduce ensemble size:**
```python
num_files_in_each_set = [3]*len(files)  # Smaller ensemble
```

4. **Use different output location:**
```python
# Save to larger disk
files = ["/large_disk/path/to/file.nc"]
```

---

### Issue 4: Overwriting Files

**Symptom:** Script reruns create same files again

**Behavior:**
```python
ds.to_netcdf(new_file, mode='w')  # mode='w' overwrites
```

**If files exist, they are overwritten without warning**

**Solutions:**

1. **Check before running:**
```bash
# List existing synthetic files
ls *_2.nc *_3.nc *_4.nc *_5.nc
```

2. **Add file existence check:**
```python
import os

new_file = file.replace('.nc', f'_{index+2}.nc')
if os.path.exists(new_file):
    print(f"Warning: {new_file} already exists. Skipping.")
    continue
ds.to_netcdf(new_file, mode='w')
```

3. **Backup existing files:**
```bash
mkdir backup
mv *_2.nc *_3.nc *_4.nc *_5.nc backup/
```

---

### Issue 5: Variable-Specific Perturbations

**Current behavior:** All variables perturbed by same multiplier

**If you want different perturbations per variable:**

```python
def produce_synthetic_spatial_data(inputs):
    start_time = time.time()
    file = inputs[0]
    num_synthetic_sets_in_ensemble = inputs[1] - 1
    
    for index in range(num_synthetic_sets_in_ensemble):
        ds = xr.open_dataset(file)
        variables = ds.data_vars
        
        # Different multiplier for each variable
        for variable in variables:
            random_multiplier = np.random.uniform(low=0.99, high=1.02)
            ds[variable] *= random_multiplier
        
        new_file = file.replace('.nc', f'_{index+2}.nc')
        ds.to_netcdf(new_file, mode='w')
```

---

## Scientific Considerations

### When to Use Synthetic Ensembles

**Appropriate uses:**
1. **Statistical significance testing** when single simulations available
2. **Sensitivity analysis** to perturbations
3. **Uncertainty quantification** when true ensemble unavailable
4. **Method development** before full ensemble runs complete

**Inappropriate uses:**
1. **Replacing true model ensembles** (different initial conditions, physics)
2. **Representing true climate variability**
3. **Publishing as actual ensemble members**

---

### Limitations

**1. Not True Ensemble Members:**
- Synthetic members don't represent true model uncertainty
- Don't capture different initial conditions
- Don't reflect model physics variability

**2. Artificial Perturbations:**
- Random multipliers don't represent physical processes
- All gridcells perturbed by same factor (spatially uniform)
- Doesn't represent realistic spatial covariance

**3. Statistical Limitations:**
- May underestimate or overestimate true variability
- P-values may not reflect true statistical significance
- Conservative for exploratory analysis only

---

### Best Practices for Publication

**If using synthetic ensembles in publications:**

1. **Clearly state synthetic nature:**
   - "Synthetic ensemble members created by random perturbations"
   - "Not true model ensemble members"

2. **Document method:**
   - Perturbation range (e.g., ±2%)
   - Ensemble size
   - Purpose (e.g., statistical testing)

3. **Appropriate interpretation:**
   - "Regions with statistically significant differences in synthetic ensemble"
   - Not "Regions with robust differences across ensemble members"

4. **Compare with true ensemble if available:**
   - Validate synthetic approach against real ensemble
   - Document any differences

---

### Alternative Approaches

**Instead of synthetic ensembles, consider:**

1. **Bootstrap resampling** of time dimension
2. **Block bootstrap** for spatial data
3. **True ensemble simulations** (preferred)
4. **Multi-model ensembles** from different E3SM configurations

---

## Advanced Modifications

### Spatially-Varying Perturbations

**Current:** Same multiplier for all gridcells  
**Modification:** Different multiplier per gridcell

```python
def produce_synthetic_spatial_data(inputs):
    start_time = time.time()
    file = inputs[0]
    num_synthetic_sets_in_ensemble = inputs[1] - 1
    
    ds = xr.open_dataset(file)
    
    for index in range(num_synthetic_sets_in_ensemble):
        ds_synthetic = ds.copy(deep=True)
        
        for variable in ds.data_vars:
            # Create spatial field of random multipliers
            shape = ds[variable].shape
            random_field = np.random.uniform(low=0.99, high=1.02, size=shape)
            ds_synthetic[variable] = ds[variable] * random_field
        
        new_file = file.replace('.nc', f'_{index+2}.nc')
        ds_synthetic.to_netcdf(new_file, mode='w')
```

---

### Normal Distribution Perturbations

**Current:** Uniform distribution  
**Modification:** Normal (Gaussian) distribution

```python
# In function:
# Instead of:
random_multipliers = np.random.uniform(low=0.99, high=1.02, size=...)

# Use:
random_multipliers = np.random.normal(loc=1.0, scale=0.01, size=...)
# loc = mean (1.0 = no change on average)
# scale = standard deviation (0.01 = 1% std)
```

---

## Appendix: Complete Modified Examples

### Example A: Custom Script for 3 Files, 10 Members Each

```python
import multiprocessing
import numpy as np
import time
import xarray as xr

def produce_synthetic_spatial_data(inputs):
    start_time = time.time()
    file = inputs[0]
    num_synthetic_sets_in_ensemble = inputs[1] - 1
    random_multipliers = np.random.uniform(low=0.99, high=1.02, size=num_synthetic_sets_in_ensemble)
    
    for index in range(len(random_multipliers)):
        ds = xr.open_dataset(file)
        variables = ds.data_vars
        for variable in variables:
            ds[variable] *= random_multipliers[index]
        new_file = file.replace('.nc', f'_{index+2}.nc')
        ds.to_netcdf(new_file, mode='w')
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for {file}: {elapsed_time:.2f} seconds")

if __name__ == '__main__':
    start_time = time.time()
    
    # Configuration
    files = [
        "./control_elm.nc",
        "./control_eam.nc",
        "./scenario_elm.nc"
    ]
    num_files_in_each_set = [10, 10, 10]  # 10 members each
    
    inputs = list(zip(files, num_files_in_each_set))
    
    with multiprocessing.Pool(processes=MAX_PROCESSES) as pool:
        pool.map(produce_synthetic_spatial_data, inputs)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total time: {elapsed_time:.2f} seconds")
```

---

## References

- E3SM Documentation: [https://e3sm.org/](https://e3sm.org/)
- NumPy Random Documentation: [https://numpy.org/doc/stable/reference/random/](https://numpy.org/doc/stable/reference/random/)
- xarray Documentation: [https://xarray.pydata.org/](https://xarray.pydata.org/)
- GitHub Repository: [https://github.com/philipmyint/e3sm--gcam_analysis](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Contact

For questions or issues:
- **Philip Myint**: myint1@llnl.gov
- **Dalei Hao**: dalei.hao@pnnl.gov

---

## Version Information

**Script:** e3sm_produce_synthetic_spatial_data.py  
**Documentation Version:** 1.0  
**Last Updated:** January 2026  
**Python Version:** 3.7+  
**Dependencies:** numpy, xarray, multiprocessing

---

*This documentation provides comprehensive guidance for using the `e3sm_produce_synthetic_spatial_data.py` script to generate synthetic ensemble members from E3SM spatial NetCDF files for statistical testing and uncertainty quantification.*
