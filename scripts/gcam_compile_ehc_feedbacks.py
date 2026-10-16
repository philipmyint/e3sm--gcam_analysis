import json
import multiprocessing
import os
import pandas as pd
import sys
import time
from utility_constants import MAX_PROCESSES
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_file
from utility_functions import get_all_files_in_path, get_max_processes_from_inputs
from utility_gcam import modify_crop_names, add_areas_to_dataframe

"""
Compile EHC feedback files into the same long-table format used by the GCAM
plotting scripts.

The E3SM human component writes several small CSV files while it runs. Those
files are organized by scenario and year, and some of them contain multiple
physical quantities in separate columns. This script combines the per-year files
across scenarios, adds scenario labels, optionally attaches land area metadata,
and writes one plottable CSV per quantity with a common column named 'value'.
"""

def compile_ehc_scalars(inputs):
    """ 
    Processes and compiles the scalars files generated at run time by the E3SM human component (EHC) and subsequently gets passed into GCAM. 
    Writes two output files: one for vegetation scalars and one for soil scalars, each with a 'value' column and
    optional 'area' and 'area_units' columns added from a processed land allocation file.

    Parameters:
        inputs: Dictionary with user-specified inputs for the directory paths, the name of the output file, whether we want to standardize the
                crop names, a land_allocation_file path, and a list of scenario names, with each scenario in the list
                corresponding to one directory.

    Returns:
        N/A.
    """
    start_time = time.time()
    input_directories = inputs['input_directories']
    output_dir = inputs['output_dir']
    call_modify_crop_names = inputs.get('call_modify_crop_names', False)
    scenarios = inputs['scenarios']
    # file_root_names is an optional dict mapping scenario names to file root names.
    # Scenarios not present in the dict fall back to the default root name 'scalars'.
    file_root_names = inputs.get('file_root_names', {})
    land_allocation_file = inputs.get('land_allocation_file', None)

    # Scalar input files contain both vegetation and soil columns. Split them
    # into separate output files so downstream plotting JSON can treat each
    # quantity the same way as other processed GCAM tables: one 'value' column
    # plus metadata columns.
    vegetation_output_file = os.path.join(output_dir, 'vegetation_scalars.csv')
    soil_output_file = os.path.join(output_dir, 'soil_scalars.csv')

    # Check that input_directories and scenarios are the same length.
    if len(input_directories) != len(scenarios):
        print(f"Error: 'input_directories' (length {len(input_directories)}) and 'scenarios' "
              f"(length {len(scenarios)}) must be the same length.")
        return

    # Check that both output files are writable before doing any processing.
    for out_file in [vegetation_output_file, soil_output_file]:
        try:
            open(out_file, 'a').close()
        except OSError as e:
            print(f"Error: cannot write to output file '{out_file}': {e}")
            return

    # Each input directory corresponds to one scenario. Read all matching
    # per-year scalar files in that directory, combine them, and tag the rows
    # with the scenario label from the JSON.
    dataframes_for_all_scenarios = []
    for index, scenario in enumerate(scenarios):
        input_directory = input_directories[index]

        # Check that the input directory exists before scanning for files.
        if not os.path.exists(input_directory):
            print(f"Warning: input directory not found: '{input_directory}' (scenario '{scenario}'). Skipping.")
            continue

        # Look up the root name for this scenario. Most EHC scalar files are
        # named scalars_<year>.csv, but file_root_names and the dispatcher's
        # _file_root_name override allow the same reader to handle alternate
        # root names without duplicating the loop below.
        file_root_name = file_root_names.get(scenario, inputs.get('_file_root_name', 'scalars'))

        # Filter to CSV files whose names start with file_root_name, to avoid reading non-scalars
        # files (e.g. land allocation, commodity prices) that may be in the same directory.
        files = get_all_files_in_path(input_directory, file_extension='.csv')
        files = [f for f in files if os.path.basename(f).startswith(file_root_name)]

        # Check that the directory contains at least one matching CSV file.
        if len(files) == 0:
            print(f"Warning: no CSV files starting with '{file_root_name}' found in '{input_directory}' "
                  f"(scenario '{scenario}'). Skipping.")
            continue

        print(f"  [{index+1}/{len(scenarios)}] Reading {len(files)} file(s) for scenario '{scenario}'...")
        dataframes_for_all_files_for_this_scenario = [read_file_into_dataframe(file) for file in files]
        df = pd.concat(dataframes_for_all_files_for_this_scenario, ignore_index=True)

        # Convert all column names to lowercase.
        df.columns = df.columns.str.lower()

        # Check that the landtype_basin column exists after lowercasing before splitting.
        if 'landtype_basin' not in df.columns:
            print(f"Error: 'landtype_basin' column not found in files in '{input_directory}'. "
                  f"Available columns: {list(df.columns)}")
            return

        df = df.sort_values(by='year')
        # The EHC scalar files encode two plotting keys together, for example
        # Corn_UsaPacNW. Split this into separate landtype and basin columns so
        # the output has the same key structure as the processed land allocation
        # tables.
        df[['landtype', 'basin']] = df['landtype_basin'].str.split('_', expand=True)
        df.drop('landtype_basin', axis=1, inplace=True)
        df['scenario'] = scenario
        dataframes_for_all_scenarios.append(df)

    if not dataframes_for_all_scenarios:
        print("Warning: no data collected for any scenario. Skipping output.")
        return
    df = pd.concat(dataframes_for_all_scenarios, ignore_index=True)
    key_columns = ['scenario', 'region', 'basin', 'landtype', 'year']

    # Check that vegetation and soil columns exist before selecting them.
    for col in ['vegetation', 'soil']:
        if col not in df.columns:
            print(f"Error: '{col}' column not found in compiled data. Available columns: {list(df.columns)}")
            return

    df = df.sort_values(by=key_columns)
    df = df[key_columns + ['vegetation', 'soil']]

    # Add areas before standardizing crop names. The scalar files still use the
    # original GCAM crop names here, so the land allocation file should be the
    # processed land allocation table that also keeps original crop names
    # (call_modify_crop_names = false in gcam_process_extracted_data.json).
    if land_allocation_file is None:
        print(f"Error: A 'land_allocation_file' must be specified in the JSON input")
        print(f"Error: Use the processed detailed land allocation file that keeps original crop names, "
              f"for example 'land_allocation_processed_original_crop_names.csv'.")
        return
    if not os.path.exists(land_allocation_file):
        print(f"Error: land_allocation_file not found: '{land_allocation_file}'")
        return
    df_land = read_file_into_dataframe(land_allocation_file)
    df = add_areas_to_dataframe(df, df_land)

    # Standardize crop names after adding areas so the merge above uses original names.
    if call_modify_crop_names:
        df = modify_crop_names(df, key_columns, 'area_weighted_mean')

    # Keep area metadata at the end of the output tables if it was successfully
    # added. Plotting scripts can use the area column for area-weighted means.
    trailing = [c for c in ['area', 'area_units'] if c in df.columns]

    # Write vegetation scalars as a standard processed GCAM table. Renaming the
    # data column to 'value' lets plotting JSON omit value_label or set it to the
    # default 'value'.
    veg_df = df[key_columns + ['vegetation'] + trailing].rename(columns={'vegetation': 'value'}).copy()
    veg_df['units'] = 'Unitless'
    veg_df = veg_df[key_columns + ['value', 'units'] + trailing]
    write_dataframe_to_file(veg_df, vegetation_output_file)

    # Write soil scalars with the same structure as vegetation_scalars.csv.
    soil_df = df[key_columns + ['soil'] + trailing].rename(columns={'soil': 'value'}).copy()
    soil_df['units'] = 'Unitless'
    soil_df = soil_df[key_columns + ['value', 'units'] + trailing]
    write_dataframe_to_file(soil_df, soil_output_file)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time processing/compiling the data for {vegetation_output_file} and {soil_output_file}: {elapsed_time:.2f} seconds")


def compile_hdd_cdd_scalars(inputs):
    """
    Processes and compiles hdd_cdd_scalars files generated at run time by the EHC.
    Writes two output files: one for HDD scalars and one for CDD scalars, each with a 'value' column.

    Parameters:
        inputs: Dictionary with 'input_directories', 'output_dir', 'scenarios', and optionally
                'file_root_names' (dict mapping scenario names to file root names; default 'hdd_cdd_scalars').

    Returns:
        N/A.
    """
    start_time = time.time()
    input_directories = inputs['input_directories']
    output_dir = inputs['output_dir']
    scenarios = inputs['scenarios']
    file_root_names = inputs.get('file_root_names', {})
    land_allocation_file = inputs.get('land_allocation_file', None)

    # The HDD/CDD scalar input stores two quantities in one row. As above, write
    # one output file per quantity so each downstream plot reads a single 'value'
    # column.
    hdd_output_file = os.path.join(output_dir, 'hdd_scalars.csv')
    cdd_output_file = os.path.join(output_dir, 'cdd_scalars.csv')

    if len(input_directories) != len(scenarios):
        print(f"Error: 'input_directories' (length {len(input_directories)}) and 'scenarios' "
              f"(length {len(scenarios)}) must be the same length.")
        return

    for out_file in [hdd_output_file, cdd_output_file]:
        try:
            open(out_file, 'a').close()
        except OSError as e:
            print(f"Error: cannot write to output file '{out_file}': {e}")
            return

    # Combine all per-year HDD/CDD scalar files across the requested scenarios.
    dataframes_for_all_scenarios = []
    for index, scenario in enumerate(scenarios):
        input_directory = input_directories[index]

        if not os.path.exists(input_directory):
            print(f"Warning: input directory not found: '{input_directory}' (scenario '{scenario}'). Skipping.")
            continue

        # Default to hdd_cdd_scalars_<year>.csv, unless file_type or
        # file_root_names asks this input block to use a different prefix.
        file_root_name = file_root_names.get(scenario, inputs.get('_file_root_name', 'hdd_cdd_scalars'))

        files = get_all_files_in_path(input_directory, file_extension='.csv')
        files = [f for f in files if os.path.basename(f).startswith(file_root_name)]

        if len(files) == 0:
            print(f"Warning: no CSV files starting with '{file_root_name}' found in '{input_directory}' "
                  f"(scenario '{scenario}'). Skipping.")
            continue

        print(f"  [{index+1}/{len(scenarios)}] Reading {len(files)} file(s) for scenario '{scenario}'...")
        dataframes_for_all_files_for_this_scenario = [read_file_into_dataframe(file) for file in files]
        df = pd.concat(dataframes_for_all_files_for_this_scenario, ignore_index=True)

        df.columns = df.columns.str.lower()

        for col in ['year', 'region', 'hdd', 'cdd']:
            if col not in df.columns:
                print(f"Error: '{col}' column not found in files in '{input_directory}'. "
                      f"Available columns: {list(df.columns)}")
                return

        df = df.sort_values(by='year')
        df['scenario'] = scenario
        dataframes_for_all_scenarios.append(df)

    if not dataframes_for_all_scenarios:
        print("Warning: no data collected for any scenario. Skipping output.")
        return
    df = pd.concat(dataframes_for_all_scenarios, ignore_index=True)
    key_columns = ['scenario', 'region', 'year']

    for col in ['hdd', 'cdd']:
        if col not in df.columns:
            print(f"Error: '{col}' column not found in compiled data. Available columns: {list(df.columns)}")
            return

    df = df.sort_values(by=key_columns)
    df = df[key_columns + ['hdd', 'cdd']]

    # HDD/CDD scalar files are regional rather than landtype/basin-specific.
    # If a land allocation file is provided, add aggregate regional areas for
    # use in area-weighted downstream calculations; otherwise leave the data as
    # pure regional time series.
    if land_allocation_file is not None:
        if not os.path.exists(land_allocation_file):
            print(f"Error: land_allocation_file not found: '{land_allocation_file}'")
            return
        df_land = read_file_into_dataframe(land_allocation_file)
        df = add_areas_to_dataframe(df, df_land)
    trailing = [c for c in ['area', 'area_units'] if c in df.columns]

    # Write one standardized output table per quantity.
    hdd_df = df[key_columns + ['hdd'] + trailing].rename(columns={'hdd': 'value'}).copy()
    hdd_df['units'] = 'Unitless'
    hdd_df = hdd_df[key_columns + ['value', 'units'] + trailing]
    write_dataframe_to_file(hdd_df, hdd_output_file)

    cdd_df = df[key_columns + ['cdd'] + trailing].rename(columns={'cdd': 'value'}).copy()
    cdd_df['units'] = 'Unitless'
    cdd_df = cdd_df[key_columns + ['value', 'units'] + trailing]
    write_dataframe_to_file(cdd_df, cdd_output_file)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time processing/compiling the data for {hdd_output_file} and {cdd_output_file}: {elapsed_time:.2f} seconds")


def compile_scaled_hdd_cdd_diag(inputs):
    """
    Processes and compiles scaled_hdd_cdd_diag files generated at run time by the EHC.
    The year is embedded in each filename (e.g. 'scaled_hdd_cdd_diag2025.csv') rather than in
    the data, so it is extracted from the filename and added as a column.
    Writes four output files, one per value column: comm_heating, comm_cooling, resid_heating,
    and resid_cooling.

    Parameters:
        inputs: Dictionary with 'input_directories', 'output_dir', 'scenarios', and optionally
                'file_root_names' (dict mapping scenario names to file root names;
                default 'scaled_hdd_cdd_diag').

    Returns:
        N/A.
    """
    start_time = time.time()
    input_directories = inputs['input_directories']
    output_dir = inputs['output_dir']
    scenarios = inputs['scenarios']
    file_root_names = inputs.get('file_root_names', {})
    land_allocation_file = inputs.get('land_allocation_file', None)

    # The diagnostic file contains four end-use columns. Create four output CSVs
    # so each file has one plottable quantity in the standard 'value' column.
    value_columns = ['comm_heating', 'comm_cooling', 'resid_heating', 'resid_cooling']
    output_files = {col: os.path.join(output_dir, f'{col}_diag.csv') for col in value_columns}

    if len(input_directories) != len(scenarios):
        print(f"Error: 'input_directories' (length {len(input_directories)}) and 'scenarios' "
              f"(length {len(scenarios)}) must be the same length.")
        return

    for out_file in output_files.values():
        try:
            open(out_file, 'a').close()
        except OSError as e:
            print(f"Error: cannot write to output file '{out_file}': {e}")
            return

    # Combine all matching diagnostic files across scenarios. Unlike the scalar
    # files above, the year is not stored in the file contents; it is encoded in
    # the filename and is added explicitly below.
    dataframes_for_all_scenarios = []
    for index, scenario in enumerate(scenarios):
        input_directory = input_directories[index]

        if not os.path.exists(input_directory):
            print(f"Warning: input directory not found: '{input_directory}' (scenario '{scenario}'). Skipping.")
            continue

        # Default to scaled_hdd_cdd_diag<year>.csv. The same override machinery
        # lets a JSON block use another prefix if a run produced differently
        # named diagnostic files.
        file_root_name = file_root_names.get(scenario, inputs.get('_file_root_name', 'scaled_hdd_cdd_diag'))

        files = get_all_files_in_path(input_directory, file_extension='.csv')
        files = [f for f in files if os.path.basename(f).startswith(file_root_name)]

        if len(files) == 0:
            print(f"Warning: no CSV files starting with '{file_root_name}' found in '{input_directory}' "
                  f"(scenario '{scenario}'). Skipping.")
            continue

        print(f"  [{index+1}/{len(scenarios)}] Reading {len(files)} file(s) for scenario '{scenario}'...")
        dfs = []
        for file in files:
            basename = os.path.basename(file)
            # File names are expected to look like
            # scaled_hdd_cdd_diag2025.csv. Slice off the prefix and suffix to
            # recover the year.
            year_str = basename[len(file_root_name):-len('.csv')]
            try:
                year = int(year_str)
            except ValueError:
                print(f"Error: could not extract year from filename '{basename}' "
                      f"(expected '{file_root_name}<yyyy>.csv').")
                return
            df_file = read_file_into_dataframe(file)
            df_file['year'] = year
            dfs.append(df_file)

        df = pd.concat(dfs, ignore_index=True)
        df.columns = df.columns.str.lower()

        for col in ['region'] + value_columns:
            if col not in df.columns:
                print(f"Error: '{col}' column not found in files in '{input_directory}'. "
                      f"Available columns: {list(df.columns)}")
                return

        df = df.sort_values(by='year')
        df['scenario'] = scenario
        dataframes_for_all_scenarios.append(df)

    if not dataframes_for_all_scenarios:
        print("Warning: no data collected for any scenario. Skipping output.")
        return
    df = pd.concat(dataframes_for_all_scenarios, ignore_index=True)
    key_columns = ['scenario', 'region', 'year']

    for col in value_columns:
        if col not in df.columns:
            print(f"Error: '{col}' column not found in compiled data. Available columns: {list(df.columns)}")
            return

    df = df.sort_values(by=key_columns)
    df = df[key_columns + value_columns]

    # Area metadata is optional for these regional diagnostic outputs. When it
    # is present, the plotting scripts can use it for area-weighted aggregation.
    if land_allocation_file is not None:
        if not os.path.exists(land_allocation_file):
            print(f"Error: land_allocation_file not found: '{land_allocation_file}'")
            return
        df_land = read_file_into_dataframe(land_allocation_file)
        df = add_areas_to_dataframe(df, df_land)
    trailing = [c for c in ['area', 'area_units'] if c in df.columns]

    # Write one file per diagnostic column, normalizing the selected column name
    # to 'value' in each output.
    for col in value_columns:
        out_df = df[key_columns + [col] + trailing].rename(columns={col: 'value'}).copy()
        out_df['units'] = 'degF/yr'
        out_df = out_df[key_columns + ['value', 'units'] + trailing]
        write_dataframe_to_file(out_df, output_files[col])

    end_time = time.time()
    elapsed_time = end_time - start_time
    output_names = ', '.join(output_files.values())
    print(f"Elapsed time processing/compiling the data for {output_names}: {elapsed_time:.2f} seconds")


def compile_ehc_feedbacks(inputs):
    """Dispatcher that calls compile functions based on the 'file_type' key.
    'file_type' can be:
      'all' (default) — runs all three types with their default root names
      a string        — runs that type, using the string as both type and root name
      a list          — runs those types, each string as both type and root name
      a dict          — keys are type identifiers, values are the root names to search for
    Valid type identifiers: 'scalars', 'hdd_cdd_scalars', 'scaled_hdd_cdd_diag'."""
    # Map each file_type identifier to the function that knows that input
    # layout. The JSON file_type field controls which of these readers runs for
    # a given input block.
    all_types = ['scalars', 'hdd_cdd_scalars', 'scaled_hdd_cdd_diag']
    type_to_func = {
        'scalars': compile_ehc_scalars,
        'hdd_cdd_scalars': compile_hdd_cdd_scalars,
        'scaled_hdd_cdd_diag': compile_scaled_hdd_cdd_diag,
    }
    # Normalize file_type into (type identifier, file prefix) pairs:
    #   "scalars" runs the scalar reader on files starting with "scalars";
    #   ["scalars", "hdd_cdd_scalars"] runs both default prefixes;
    #   {"scalars": "my_scalars"} runs the scalar reader with a custom prefix.
    raw = inputs.get('file_type', 'all')
    if raw == 'all':
        type_root_pairs = [(t, t) for t in all_types]
    elif isinstance(raw, dict):
        type_root_pairs = list(raw.items())
    elif isinstance(raw, str):
        type_root_pairs = [(raw, raw)]
    else:
        type_root_pairs = [(t, t) for t in raw]
    for file_type, root_name in type_root_pairs:
        if file_type not in type_to_func:
            print(f"Error: unknown file_type '{file_type}'. Valid types: {all_types}.")
            return
        modified_inputs = dict(inputs)
        # Store the chosen root name in a private key so the lower-level compile
        # function can use the same input dictionary without changing its public
        # JSON interface.
        modified_inputs['_file_root_name'] = root_name
        type_to_func[file_type](modified_inputs)


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file(s) on the command line.
    start_time = time.time()
    if len(sys.argv) < 2:
        print('Usage: python gcam_compile_ehc_feedbacks.py `path/to/json/input/file(s)\'')
        sys.exit()

    # Read and load the JSON file(s) into a list of dictionaries.
    list_of_inputs = []
    for i in range(1, len(sys.argv)):
        input_file = sys.argv[i]
        with open(input_file) as f:
            list_of_inputs.extend(json.load(f))

    # Use any max_processes value requested in the inputs, otherwise use the default.
    # Do not start more worker processes than there are independent input blocks.
    max_processes = get_max_processes_from_inputs(list_of_inputs, MAX_PROCESSES)
    num_processes = min(max_processes, max(1, len(list_of_inputs)))

    # Each top-level JSON block is independent, so multiple blocks can be
    # compiled in parallel. Inside each block, the files for that block are read
    # serially to keep the per-output ordering and error messages simple.
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(compile_ehc_feedbacks, list_of_inputs)
    
    # Print the total execution time needed to process/compile all requested
    # feedback files.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time processing/compiling the data for all files: {elapsed_time:.2f} seconds")
