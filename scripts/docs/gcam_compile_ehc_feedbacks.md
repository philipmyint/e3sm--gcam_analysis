# GCAM EHC Feedback Compilation

**Script:** `gcam_compile_ehc_feedbacks.py`

## Overview

`gcam_compile_ehc_feedbacks.py` compiles CSV files written by the E3SM Human Component (EHC) into the standard long-table format used by the GCAM plotting scripts.

The script can compile three EHC output families:

| `file_type` | Input files | Output files |
|-------------|-------------|--------------|
| `scalars` | `scalars_<year>.csv` | `vegetation_scalars.csv`, `soil_scalars.csv` |
| `hdd_cdd_scalars` | `hdd_cdd_scalars_<year>.csv` | `hdd_scalars.csv`, `cdd_scalars.csv` |
| `scaled_hdd_cdd_diag` | `scaled_hdd_cdd_diag<year>.csv` | `comm_heating_diag.csv`, `comm_cooling_diag.csv`, `resid_heating_diag.csv`, `resid_cooling_diag.csv` |

All output files use a common plotting shape:

```text
scenario, region, year, value, units
```

Vegetation and soil scalar outputs also include `basin` and `landtype`. Files can include `area` and `area_units` when a processed land allocation file is supplied.

## Usage

```bash
python gcam_compile_ehc_feedbacks.py path/to/config.json
```

Multiple JSON files can be supplied:

```bash
python gcam_compile_ehc_feedbacks.py config1.json config2.json
```

Each configuration block is processed according to its `file_type`. If `file_type` is omitted or set to `all`, the script attempts all supported EHC feedback families.

## Parameters

| Parameter | Required | Description |
|-----------|----------|-------------|
| `input_directories` | Yes | One directory per scenario. Each directory is scanned for CSV files matching the selected file root. |
| `scenarios` | Yes | Scenario labels, in the same order as `input_directories`. |
| `output_dir` | Yes | Directory where compiled CSV files are written. |
| `file_type` | No | One of `all`, `scalars`, `hdd_cdd_scalars`, `scaled_hdd_cdd_diag`, or a list of those values. Default is `all`. |
| `file_root_names` | No | Optional dictionary mapping scenario names to custom file prefixes. |
| `call_modify_crop_names` | No | If true, standardize vegetation/soil landtype names after area metadata is added. |
| `land_allocation_file` | Required for `scalars`; optional for HDD/CDD files | Processed detailed land allocation file used to add `area` and `area_units`. |
| `max_processes` | No | Maximum parallel workers. Defaults to `MAX_PROCESSES`, which is all available CPUs minus one, with a minimum of one. |

## Scalars

For `file_type: "scalars"`, each input row is expected to include:

```text
year, region, landtype_basin, vegetation, soil
```

The script splits `landtype_basin` into `landtype` and `basin`, adds the scenario label, adds area metadata from `land_allocation_file`, optionally standardizes crop names, and writes two output files:

```text
vegetation_scalars.csv
soil_scalars.csv
```

Each output has a single numeric `value` column, so plotting JSON can use the default `value_label: "value"`.

## HDD/CDD Scalars

For `file_type: "hdd_cdd_scalars"`, input files are expected to include:

```text
year, region, hdd, cdd
```

The outputs are regional scalar tables:

```text
hdd_scalars.csv
cdd_scalars.csv
```

These files do not have `landtype`, `basin`, or `sector` columns. Spatial and time-series plots can use them directly. Box plots need an explicit `x_variable`, such as `"scenario"` or `"region"`.

## Scaled HDD/CDD Diagnostics

For `file_type: "scaled_hdd_cdd_diag"`, the year is read from the filename, such as `scaled_hdd_cdd_diag2025.csv`. Each input file is expected to include:

```text
region, comm_heating, comm_cooling, resid_heating, resid_cooling
```

The outputs are:

```text
comm_heating_diag.csv
comm_cooling_diag.csv
resid_heating_diag.csv
resid_cooling_diag.csv
```

## Example

```json
[
    {
        "input_directories": [
            "./../2025_DiVittorio_et_al_gcam/scalars_ehc_output/control",
            "./../2025_DiVittorio_et_al_gcam/scalars_ehc_output/full_feedback",
            "./../2025_DiVittorio_et_al_gcam/scalars_ehc_output/ag_scaling",
            "./../2025_DiVittorio_et_al_gcam/scalars_ehc_output/carbon_scaling"
        ],
        "scenarios": [
            "Control",
            "Full feedback",
            "Ag scaling",
            "Carbon scaling"
        ],
        "output_dir": "./../2025_DiVittorio_et_al_gcam",
        "file_type": "scalars",
        "call_modify_crop_names": true,
        "land_allocation_file": "./../2025_DiVittorio_et_al_gcam/land_allocation_processed_original_crop_names.csv"
    }
]
```

## Workflow Notes

Run `gcam_process_extracted_data.py` first to create the processed land allocation file. The scalar compiler needs the original crop-name version of that file because EHC scalar inputs still use original GCAM crop names when area metadata is attached.

Run `gcam_produce_synthetic_data.py` after this script if synthetic ensemble versions are needed. For scalar outputs, the synthetic script now reads `vegetation_scalars.csv` and `soil_scalars.csv`, producing `vegetation_scalars_ensemble.csv` and `soil_scalars_ensemble.csv`.
