import json
import multiprocessing
import pandas as pd
import sys
import time
from utility_dataframes import move_columns_next_to_each_other_in_dataframe, read_file_into_dataframe, write_dataframe_to_fwf
from utility_functions import get_all_files_in_path

def compile_time_series(inputs):
    """ 
    Compiles the scalars time series files located in multiple directories into a single output file that gets written to a .dat or .csv file.

    Parameters:
        inputs: Dictionary with user-specified inputs for the directory paths, the name of the output file, and a list of scenario names,
                with each scenario in the list corresponding to one directory.

    Returns:
        N/A.
    """
    start_time = time.time()
    input_directories = inputs['input_directories']
    output_file = inputs['output_file']
    scenarios = inputs['scenarios']
    dataframes_for_all_scenarios = []
    for index, scenario in enumerate(scenarios):
        input_directory = input_directories[index]
        files = get_all_files_in_path(input_directory)
        dataframes_for_all_files_for_this_scenario = [read_file_into_dataframe(file) for file in files]
        df = pd.concat(dataframes_for_all_files_for_this_scenario, ignore_index=True)
        df = df.sort_values(by='Year')
        df['Scenario'] = scenario
        dataframes_for_all_scenarios.append(df)
    df = pd.concat(dataframes_for_all_scenarios, ignore_index=True)
    df = move_columns_next_to_each_other_in_dataframe(df, 'Year', 'Scenario')
    if output_file.endswith('.csv'):
        df.to_csv(output_file, index=False)
    else:
        write_dataframe_to_fwf(output_file, df)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for compiling the data for {output_file}: {elapsed_time:.2f} seconds")


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file(s) on the command line.
    start_time = time.time()
    if len(sys.argv) < 2:
        print('Usage: python plot_time_series.py `path/to/json/input/file(s)\'')
        sys.exit()

    # Read and load the JSON file(s) into a list of dictionaries.
    list_of_inputs = []
    for i in range(1, len(sys.argv)):
        input_file = sys.argv[i]
        with open(input_file) as f:
            list_of_inputs.extend(json.load(f))

    # Produce data for each file in parallel.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(compile_time_series, list_of_inputs)
    
    # Print the total execution time needed to compile the scalars for all files.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for compiling the data for all files: {elapsed_time:.2f} seconds")