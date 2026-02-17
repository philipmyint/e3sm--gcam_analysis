# E3SM Analysis Debugging Tools

Quick local testing tools to debug E3SM extraction configurations without waiting in SLURM queues.

## Problem

SLURM queues can take hours to start jobs, making configuration debugging very slow:
```bash
# This could take 4+ hours to even start running
sbatch submit_e3sm_time_series.sh  

# Then fail immediately due to a simple config error
# ✗ IndexError: list index out of range
# ✗ KeyError: 'VARIABLE_NAME'  
# ✗ Did not recognize region...
```

## Solution

Three local testing tools that run in **seconds to minutes** instead of hours:

| Tool | Speed | Purpose | Dependencies |
|------|-------|---------|--------------|
| `test_config.py` | 30 sec | JSON structure validation | None (basic Python) |
| `dry_run_test.py` | 1 min | Logic validation | None |
| `small_data_test.py` | 5-10 min | Full pipeline test | E3SM environment |

## Tools Overview

### 1. `test_config.py` - JSON Validation
**Purpose**: Validate JSON configuration structure and file paths  
**Speed**: ~30 seconds  
**Dependencies**: Basic Python only

**What it checks**:
- ✅ JSON file loads correctly
- ✅ Required fields present (`variables`, `netcdf_substrings`, etc.)
- ✅ List lengths match between parameters
- ✅ Simulation paths exist
- ✅ NetCDF files found in directories
- ✅ Output directories accessible

**Usage**:
```bash
python test_config.py e3sm_extract_time_series_h0_regs.json
```

**Example output**:
```
✓ JSON file loaded successfully
✓ Found 14 configuration entries
✓ Variables groups: 1, NetCDF substring groups: 1 ✓ Lengths match!
✓ Simulation path exists: /global/cfs/cdirs/e3sm/lvroekel/archive/lnd/hist
  Found 5098 .nc files
✓ Configuration appears valid for SLURM submission!
```

### 2. `dry_run_test.py` - Validation Logic Testing
**Purpose**: Test the E3SM script's validation logic without processing files  
**Speed**: ~1 minute  
**Dependencies**: Basic Python only

**What it checks**:
- ✅ List length validation (variables vs netcdf_substrings vs regions)
- ✅ Default parameter creation logic
- ✅ Parameter consistency across all entries
- ✅ File path accessibility

**Usage**:
```bash
python dry_run_test.py e3sm_extract_time_series_h0_regs.json
```

**Example output**:
```
✓ Loaded 14 entries from JSON
✓ Length validation: PASS
✓ All parameter lengths consistent
✓ Simulation path accessible
  Group 1: 3 variables, substring ['elm.h0']
✓ All validations passed! Ready for SLURM submission.
```

### 3. `small_data_test.py` - Full Pipeline Test
**Purpose**: Run the actual E3SM extraction code on a small subset of data  
**Speed**: 5-10 minutes  
**Dependencies**: E3SM environment 

**What it does**:
- ✅ Loads E3SM extraction modules
- ✅ Creates test config with limited years (2-3 years instead of 392)
- ✅ Runs full extraction pipeline on real data
- ✅ Tests variable processing, region recognition, file I/O
- ✅ Outputs test files prefixed with `test_`

**Usage**:
```bash
<source E3SM-unified environment>   # Load E3SM environment first
python small_data_test.py e3sm_extract_time_series_h0_regs.json
```

**Example output**:
```
Created test configuration: test_config_small.json
- Limited to 3 years per entry
- Output files prefixed with 'test_'
✓ Entry 1 completed successfully
✓ Entry 2 completed successfully
✓ Small data test completed successfully!
✓ Ready for full SLURM job submission
```

## Recommended Workflow

**Step 1: Quick JSON Check** (30 seconds)
```bash
python test_config.py your_config.json
```
🎯 **Catches**: JSON syntax errors, missing fields, file path issues

**Step 2: Validation Logic Test** (1 minute)
```bash  
python dry_run_test.py your_config.json
```
🎯 **Catches**: List length mismatches, parameter inconsistencies

**Step 3: Small Scale Test** (5-10 minutes) 
```bash
source_e3sm && python small_data_test.py your_config.json
```
🎯 **Catches**: Variable name errors, region recognition issues, processing errors

**Step 4: Submit Full Job** 
```bash
sbatch submit_e3sm_time_series.sh
```

## Common Issues These Tools Catch

### ✗ List Length Mismatches
**Error**: `IndexError: list index out of range`
```json
{
  "variables": [["GPP", "NPP"]],           // Length 1
  "netcdf_substrings": [["elm.h0"], ["eam.h0"]]  // Length 2 ❌
}
```
**Caught by**: `dry_run_test.py`

### ✗ Region Name Issues  
**Error**: `Did not recognize the selected region`
```json
{
  "regions": [["amazon", "boam", "tena"]]  // ❌ Should be separate entries
}
```
**Should be**:
```json
[
  {"regions": ["amazon"], "output_file": "amazon.dat"},
  {"regions": ["boam"], "output_file": "boam.dat"}
]
```
**Caught by**: `small_data_test.py`

### ✗ Variable Name Errors
**Error**: `KeyError: 'VARIABLE_NAME'`  
```json
{
  "variables": [["GPP", "TREFHT"]]  // ❌ TREFHT is EAM, not ELM
}
```
**Caught by**: `small_data_test.py`

### ✗ Path Issues
**Error**: `FileNotFoundError`  
```json
{
  "simulation_path": "/wrong/path/to/simulation"  // ❌
}
```
**Caught by**: `test_config.py`

## Time Savings

| Method | Time to Results | Success Rate |
|--------|----------------|--------------|
| **Direct SLURM submission** | 4-8 hours (queue + runtime) | ~30% (many config errors) |
| **Local debugging → SLURM** | 10 minutes + 2-3 hours | ~95% (pre-validated) |

**Result**: Save 4-6 hours per configuration iteration!

## Advanced Usage

### Test Multiple Configurations
```bash
# Test all configs in a directory
for config in configs/*.json; do
    echo "Testing $config..."
    python test_config.py "$config" || echo "❌ $config has issues"
done
```

### Create Test Configs for Development
```bash  
# The small_data_test.py automatically creates test_config_small.json
# Use this for rapid iteration:
python small_data_test.py myconfig.json  # Creates test_config_small.json
python e3sm_extract_time_series_h0.py test_config_small.json  # Quick runs
```

### Batch Testing
```bash
# Test validation on all entries
python dry_run_test.py large_config_with_many_entries.json

# Test just first few entries with real data  
head -50 large_config.json > test_subset.json  # Edit to valid JSON
python small_data_test.py test_subset.json
```

## Best Practices

1. **Always run all 3 tools before SLURM submission**
2. **Fix issues in order**: JSON → validation → small data → full job
3. **Use descriptive output filenames** to track which test produced what
4. **Keep test configs around** for future quick iterations
5. **Run tests after any configuration changes**

## Output Files

The tools create minimal output:
- `test_config.py`: No files (just validation output)
- `dry_run_test.py`: No files (just validation output)  
- `small_data_test.py`: Creates `test_config_small.json` and `test_*.dat` files

Clean up with:
```bash
rm -f test_config_small.json output/test_*.dat
```

## Troubleshooting

**"ModuleNotFoundError: No module named 'numpy'"**
- Solution: Run <E3SM unified environemnt> first (only needed for `small_data_test.py`)

**"No files found matching criteria"**  
- Check simulation paths in your JSON
- Run `ls /your/simulation/path/*.nc | head` to verify files exist

**"Did not recognize the selected region"**
- Each region needs its own configuration entry
- Valid regions: `amazon`, `boam`, `tena`, `boas`, `eqas`, `noam`, `soam`, `conus`, etc.

**Tests pass but SLURM job still fails**
- Check SLURM resource requirements (memory, time limits)
- Verify the SLURM script uses the same JSON file as your tests
- Check for environment differences between login nodes and compute nodes

---

## Summary

These tools turn configuration debugging from a **multi-hour ordeal** into a **10-minute process**. Always test locally before submitting to SLURM!

**Quick check sequence**:
```bash
python test_config.py config.json && \
python dry_run_test.py config.json && \
source_e3sm && python small_data_test.py config.json && \
echo "✅ Ready for SLURM!" || echo "❌ Fix issues first"
```