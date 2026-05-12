import json
import multiprocessing
import os
import pandas as pd
import sys
import time
from utility_constants import MAX_PROCESSES
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_file
from utility_functions import get_all_files_in_path
from utility_gcam import modify_crop_names

def compile_ehc_scalars(inputs):
    """ 
    Processes and compiles the scalars files generated at run time by the E3SM human component (EHC) and subsequently gets passed into GCAM. 
    The files are located in multiple directories into a single output file that gets written to a .dat or .csv file.

    Parameters:
        inputs: Dictionary with user-specified inputs for the directory paths, the name of the output file, whether we want to standardize the
                crop names, and a list of scenario names, with each scenario in the list corresponding to one directory.

    Returns:
        N/A.
    """
    start_time = time.time()
    input_directories = inputs['input_directories']
    output_file = inputs['output_file']
    call_modify_crop_names = inputs.get('call_modify_crop_names', False)
    scenarios = inputs['scenarios']
    # file_root_names is an optional dict mapping scenario names to file root names.
    # Scenarios not present in the dict fall back to the default root name 'scalars'.
    file_root_names = inputs.get('file_root_names', {})

    # Check that input_directories and scenarios are the same length.
    if len(input_directories) != len(scenarios):
        print(f"Error: 'input_directories' (length {len(input_directories)}) and 'scenarios' "
              f"(length {len(scenarios)}) must be the same length.")
        return

    # Check that the output file is writable before doing any processing.
    try:
        open(output_file, 'a').close()
    except OSError as e:
        print(f"Error: cannot write to output file '{output_file}': {e}")
        return

    dataframes_for_all_scenarios = []
    for index, scenario in enumerate(scenarios):
        input_directory = input_directories[index]

        # Check that the input directory exists before scanning for files.
        if not os.path.exists(input_directory):
            print(f"Error: input directory not found: '{input_directory}' (scenario '{scenario}').")
            return

        # Look up the root name for this scenario from the file_root_names dict; default to 'scalars'.
        file_root_name = file_root_names.get(scenario, 'scalars')

        # Filter to CSV files whose names start with file_root_name, to avoid reading non-scalars
        # files (e.g. land allocation, commodity prices) that may be in the same directory.
        files = get_all_files_in_path(input_directory, file_extension='.csv')
        files = [f for f in files if os.path.basename(f).startswith(file_root_name)]

        # Check that the directory contains at least one matching CSV file.
        if len(files) == 0:
            print(f"Error: no CSV files starting with '{file_root_name}' found in '{input_directory}' "
                  f"(scenario '{scenario}').")
            return

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
        # Split the landtype_basin column into 'landtype' and 'basin' columns. Delete the original 'landtype_basin' column.
        df[['landtype', 'basin']] = df['landtype_basin'].str.split('_', expand=True)
        df.drop('landtype_basin', axis=1, inplace=True)
        df['scenario'] = scenario
        dataframes_for_all_scenarios.append(df)

    df = pd.concat(dataframes_for_all_scenarios, ignore_index=True)
    key_columns = ['scenario', 'region', 'basin', 'landtype', 'year']

    # Check that vegetation and soil columns exist before selecting them.
    for col in ['vegetation', 'soil']:
        if col not in df.columns:
            print(f"Error: '{col}' column not found in compiled data. Available columns: {list(df.columns)}")
            return

    df = df.sort_values(by=key_columns)
    df = df[key_columns + ['vegetation', 'soil']]

    # Update original crop names to a common, standardized set of names. 
    if call_modify_crop_names:
        df = modify_crop_names(df, key_columns)

    write_dataframe_to_file(df, output_file)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time processing/compiling the data for {output_file}: {elapsed_time:.2f} seconds")


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file(s) on the command line.
    start_time = time.time()
    if len(sys.argv) < 2:
        print('Usage: python gcam_compile_ehc_scalars.py `path/to/json/input/file(s)\'')
        sys.exit()

    # Read and load the JSON file(s) into a list of dictionaries.
    list_of_inputs = []
    for i in range(1, len(sys.argv)):
        input_file = sys.argv[i]
        with open(input_file) as f:
            list_of_inputs.extend(json.load(f))

    # Produce data for each file in parallel.
    with multiprocessing.Pool(processes=MAX_PROCESSES) as pool:
        pool.map(compile_ehc_scalars, list_of_inputs)
    
    # Print the total execution time needed to process/compile the scalars for all files.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time processing/compiling the data for all files: {elapsed_time:.2f} seconds")