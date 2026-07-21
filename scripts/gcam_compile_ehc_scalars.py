import json
import multiprocessing
import os
import pandas as pd
import sys
import time
from utility_constants import MAX_PROCESSES
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_file
from utility_functions import get_all_files_in_path
from utility_gcam import modify_crop_names, add_areas_to_dataframe

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

    # Derive output paths for vegetation and soil files from output_dir.
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

    # Add areas from the detailed land allocation file before standardizing crop names,
    # since the detailed land allocation file uses the original crop names.
    if land_allocation_file is None:
        print(f"Error: A 'land_allocation_file' must be specified in the JSON input")
        print(f"Error: The detailed land allocation file ('detailed_land_allocation.csv') must be provided "
              f"here (not the processed version), since it uses the original crop names.")
        print(f"Error: The default 'land_allocation_file' is called 'detailed_land_allocation_processed.csv' "
              f"and should be located in the 'output/gcam_csv/' directory.")
        return
    if not os.path.exists(land_allocation_file):
        print(f"Error: land_allocation_file not found: '{land_allocation_file}'")
        return
    df_land = read_file_into_dataframe(land_allocation_file)
    df = add_areas_to_dataframe(df, df_land)

    # Standardize crop names after adding areas so the merge above uses original names.
    if call_modify_crop_names:
        df = modify_crop_names(df, key_columns, 'area_weighted_mean')

    # Trailing columns present after add_areas (area and area_units if available).
    trailing = [c for c in ['area', 'area_units'] if c in df.columns]

    # Write vegetation scalars: rename vegetation -> value, add units column.
    veg_df = df[key_columns + ['vegetation'] + trailing].rename(columns={'vegetation': 'value'}).copy()
    veg_df['units'] = 'Unitless'
    veg_df = veg_df[key_columns + ['value', 'units'] + trailing]
    write_dataframe_to_file(veg_df, vegetation_output_file)

    # Write soil scalars: rename soil -> value, add units column.
    soil_df = df[key_columns + ['soil'] + trailing].rename(columns={'soil': 'value'}).copy()
    soil_df['units'] = 'Unitless'
    soil_df = soil_df[key_columns + ['value', 'units'] + trailing]
    write_dataframe_to_file(soil_df, soil_output_file)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time processing/compiling the data for {vegetation_output_file} and {soil_output_file}: {elapsed_time:.2f} seconds")


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