# GCAM CSV Extraction Script Documentation

**Authors:** Philip Myint (myint1@llnl.gov), Dalei Hao (dalei.hao@pnnl.gov), Sha Feng (sha.feng@pnnl.gov), and Eva Sinha (eva.sinha@pnnl.gov)

**Repository:** [E3SM-GCAM Analysis Scripts](https://github.com/philipmyint/e3sm--gcam_analysis)

---

## Overview

The `gcam_extract_csv_from_xml_or_project_files.R` script extracts variables from GCAM output and writes them to CSV files, one file per variable. It supports two modes of operation: reading from pre-existing project files, or creating new project files from raw GCAM XML output before extracting.

## Purpose

This script automates the extraction of GCAM output data by:
- Reading from existing GCAM project files (.dat format), or
- Creating project files from raw GCAM XML output (GCAMDBOutput.xml) using a query file
- Extracting specified variables from the project files
- Combining data from multiple scenarios into single CSV files
- Caching project files in memory to avoid redundant loading across multiple variables

## Requirements

### R Packages
- `dplyr` — Data manipulation and transformation
- `rgcam` — GCAM-specific data reading and processing functions
- `rjson` — JSON file parsing

Install missing packages using:
```r
install.packages(c("dplyr", "rjson"))
devtools::install_github("JGCRI/rgcam")
```

---

## Input Files

### 1. JSON Configuration File

Each block in the JSON file specifies one variable to extract and write to a CSV file.

#### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `variable` | string | Exact name of the GCAM query to extract; must match a query title in the project file |
| `outputFile` | string | Path and filename for the CSV output |
| `scenarios` | array | Descriptive scenario names, one per project file |
| `projectFiles` | array | Paths to GCAM project (.dat) files, in the same order as `scenarios` |

#### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `queryFile` | string | `Main_queries.xml` | Path to the XML query file. Only required when `createProjectFiles` is `true`; ignored when reading from existing project files |
| `createProjectFiles` | boolean | `false` | If `true`, creates project files from raw XML output before extracting. Requires `xmlFiles` and `queryFile` |
| `xmlFiles` | array | — | Paths to raw GCAM XML output files (one per scenario). Required when `createProjectFiles` is `true` |
| `maxMemory` | string | `"8g"` | Java heap memory limit for the `addScenario` step. Increase for large XML files (e.g., `"32g"`) |

### 2. GCAM Project Files (.dat)

Binary project files containing pre-computed query results, created by `addScenario` from raw GCAM XML output. These may already exist (e.g., transferred from an HPC cluster) or can be created by the script itself using `createProjectFiles: true`.

### 3. Query File (optional)

An XML file containing query definitions used to create project files. Only needed when `createProjectFiles: true`. Two query files are included in the repository:

- **`energy_and_water_queries.xml`** — A small targeted file containing 4 queries (irrigation water withdrawals, biophysical water demand, total final energy, building total final energy). Use this for fast extraction of specific variables since `addScenario` only runs the queries defined in the file.
- **`Main_queries.xml`** — A comprehensive file with 268 queries covering all major GCAM output variables. Note that 10 of these use complex XPath syntax incompatible with older BaseX versions (see Troubleshooting).

The `variable` name in the JSON must exactly match a query title in the query file used to create the project file.

---

## Modes of Operation

### Mode 1: Reading from Existing Project Files

Set `createProjectFiles` to `false` (or omit it). The script loads the project file and calls `getQuery` to retrieve the pre-computed result. No `queryFile` or `xmlFiles` are needed.

```json
{
    "variable": "ag commodity prices",
    "outputFile": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices.csv",
    "scenarios": ["Control", "Full feedback"],
    "projectFiles": [
        "./../2025_DiVittorio_et_al_gcam/20240730_SSP245_ZATM_control.dat",
        "./../2025_DiVittorio_et_al_gcam/20240730_SSP245_ZATM_full_feedback.dat"
    ]
}
```

### Mode 2: Creating Project Files from Raw XML

Set `createProjectFiles: true` and provide `xmlFiles` and `queryFile`. The script runs `addScenario` to create the project file (running all queries in `queryFile` against the XML database), then immediately extracts the requested variable.

**Important:** `addScenario` runs all queries in the query file, not just the one requested. Use a small targeted query file (e.g., `energy_and_water_queries.xml`) rather than `Main_queries.xml` when only a few variables are needed, to avoid waiting for hundreds of queries to complete.

Only the first block that creates a given project file needs `createProjectFiles: true`. Subsequent blocks that reference the same project file can set `createProjectFiles: false` — the script caches the loaded project file in memory so it is not reloaded.

```json
{
    "variable": "irrigation water withdrawals by region",
    "queryFile": "energy_and_water_queries.xml",
    "outputFile": "./../2025_DiVittorio_et_al_gcam/irrigation_water_withdrawals_by_region.csv",
    "scenarios": ["Control"],
    "projectFiles": ["./../2025_DiVittorio_et_al_gcam/extraction_from_gcamdboutput.dat"],
    "xmlFiles": ["./../2025_DiVittorio_et_al_gcam/GCAMDBOutput.xml"],
    "createProjectFiles": true,
    "maxMemory": "32g"
}
```

---

## Example JSON Configuration

```json
[
    {
        "variable": "ag commodity prices",
        "outputFile": "./../2025_DiVittorio_et_al_gcam/ag_commodity_prices.csv",
        "scenarios": ["Control", "Full feedback"],
        "projectFiles": [
            "./../2025_DiVittorio_et_al_gcam/20240730_SSP245_ZATM_control.dat",
            "./../2025_DiVittorio_et_al_gcam/20240730_SSP245_ZATM_full_feedback.dat"
        ]
    },
    {
        "variable": "CO2 emissions by sector",
        "outputFile": "./../2025_DiVittorio_et_al_gcam/co2_emissions_sectors.csv",
        "scenarios": ["Control", "Full feedback"],
        "projectFiles": [
            "./../2025_DiVittorio_et_al_gcam/20240730_SSP245_ZATM_control.dat",
            "./../2025_DiVittorio_et_al_gcam/20240730_SSP245_ZATM_full_feedback.dat"
        ]
    },
    {
        "variable": "irrigation water withdrawals by region",
        "queryFile": "energy_and_water_queries.xml",
        "outputFile": "./../2025_DiVittorio_et_al_gcam/irrigation_water_withdrawals_by_region.csv",
        "scenarios": ["Control"],
        "projectFiles": ["./../2025_DiVittorio_et_al_gcam/extraction_from_gcamdboutput.dat"],
        "xmlFiles": ["./../2025_DiVittorio_et_al_gcam/GCAMDBOutput.xml"],
        "createProjectFiles": true,
        "maxMemory": "32g"
    },
    {
        "variable": "total final energy by region",
        "outputFile": "./../2025_DiVittorio_et_al_gcam/total_final_energy_by_region.csv",
        "scenarios": ["Control"],
        "projectFiles": ["./../2025_DiVittorio_et_al_gcam/extraction_from_gcamdboutput.dat"]
    }
]
```

In this example, blocks 1 and 2 read from existing project files. Block 3 creates a new project file from `GCAMDBOutput.xml` using `energy_and_water_queries.xml`, then extracts irrigation water withdrawals. Block 4 reads from the same project file created by block 3 — no `createProjectFiles` needed since it is already cached.

---

## Running the Script

```bash
Rscript gcam_extract_csv_from_xml_or_project_files.R path/to/config.json
```

Multiple JSON files can be specified:
```bash
Rscript gcam_extract_csv_from_xml_or_project_files.R config1.json config2.json
```

---

## Output Files

Each output CSV file contains all columns from the GCAM query result plus a `scenario` column identifying which scenario each row belongs to. Data from all specified scenarios is row-stacked into a single file.

---

## Script Workflow

1. Parse command line arguments and validate that at least one JSON file was provided
2. For each JSON block:
   - Validate that `scenarios` and `projectFiles` are the same length
   - If `createProjectFiles: true`, validate `xmlFiles` and `queryFile`; run `addScenario` to create the project file
   - For each scenario: load the project file (using the cache if already loaded), call `getQuery` to retrieve the variable, add the scenario label
   - If a query returns no data, print a warning and skip to the next block
   - Combine all scenario dataframes and write to the output CSV

---

## Querying Available Variables

To see which variables are stored in an existing project file:
```r
library(rgcam)
project = loadProject("path/to/project.dat")
print(listQueries(project))
print(listScenarios(project))
```

The 268 queries available in `Main_queries.xml` are listed in `Main_queries_list.md` in the repository. Of these, 258 are compatible with BaseX 9.5.0 (bundled with rgcam 1.2.0); the 10 incompatible ones use complex XPath expressions. Use `Main_queries_compatible.xml` (with the 10 incompatible queries removed) if running `addScenario` on a system with an older BaseX version.

---

## Error Handling

The script stops with an error if:
- No JSON file is provided on the command line
- A JSON file does not exist
- `scenarios` and `projectFiles` have different lengths
- `createProjectFiles: true` but `xmlFiles` is not specified
- A project file does not exist (and `createProjectFiles` is `false`)
- The output directory does not exist
- The query file does not exist (when `createProjectFiles: true`)

If a query returns no data (e.g., the variable is not tracked in that simulation), the script prints a warning and skips to the next block rather than halting.

---

## Troubleshooting

**"Query not found" warning**: The variable name does not match any query stored in the project file. This can mean the query was not included in the query file when the project file was created, or the variable is simply not tracked in this simulation. Check `listQueries(project)` to see what is available.

**`addScenario` fails with `[XPST0003]` error**: The query file contains XPath expressions incompatible with your BaseX version. Use `Main_queries_compatible.xml` (which removes the 10 incompatible queries) or a smaller targeted query file like `energy_and_water_queries.xml`.

**`addScenario` is very slow**: This step runs all queries in the query file against the full XML database. Use a small targeted query file rather than `Main_queries.xml` when only a few variables are needed.

**Scenario name mismatch**: The scenario name stored in the project file (visible via `listScenarios(project)`) must match the name specified in `scenarios` in the JSON. For project files created from `GCAMDBOutput.xml`, the scenario name is typically `gcam_test`.
