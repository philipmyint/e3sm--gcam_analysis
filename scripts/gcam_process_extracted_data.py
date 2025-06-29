import json
import multiprocessing
import pandas as pd
import sys
import time
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_fwf
from utility_gcam import standardize_crop_names

def process_extracted_data(inputs):
    """ 
    Processes a file containing data extracted from GCAM project files using XML queries. The data are rearranged, split, aggregated in various ways,
    including potentially aggregating landtype groups, and the resulting processed data then gets written as an output file in .dat or .csv format.

    Parameters:
        inputs: Dictionary with user-specified inputs, like the names of the input file (the one extracted from GCAM) and the output (processed) file.

    Returns:
        N/A.
    """
    start_time = time.time()
    input_file = inputs['input_file']
    output_file = inputs['output_file']
    columns_to_drop = inputs.get('columns_to_drop', None) 
    columns_to_split = inputs.get('columns_to_split', None) 
    key_columns = inputs.get('key_columns', None) 
    mean_or_sum_if_more_than_one_row_in_same_landtype_group = inputs.get('mean_or_sum_if_more_than_one_row_in_same_landtype_group', None) 
    call_standardize_crop_names = inputs.get('call_standardize_crop_names', False)
    
    df = read_file_into_dataframe(input_file)
    if columns_to_drop:
        df = df.drop(columns_to_drop, axis=1)
    df.columns = df.columns.str.lower()
    df = df.sort_values(by='year')

    # Split each specified column (if any) separated by underscores into a set of new columns. 
    if columns_to_split:
        for column_to_split, new_columns in columns_to_split.items():
            df[new_columns] = df[column_to_split].str.split('_', expand=True).iloc[:, :len(new_columns)]

    # Sort columns by a set of keys if specified to do so.
    if key_columns:
        df = df.sort_values(by=key_columns)
        df = df[key_columns + ['value']]

    # Perform either a sum or a mean if there happens to be more than one row that agrees on all the table keys (i.e., belongs to the same group).
    if mean_or_sum_if_more_than_one_row_in_same_landtype_group == 'sum' and key_columns:
        df = df.groupby(key_columns).sum().reset_index()
    elif mean_or_sum_if_more_than_one_row_in_same_landtype_group == 'mean' and key_columns:
        df = df.groupby(key_columns).mean().reset_index()

    # Update any non-standard crop names to belong to the standard set. No action is performed by this function is there are no non-standard names.
    if call_standardize_crop_names:
        df = standardize_crop_names(df, key_columns, mean_or_sum_if_more_than_one_row_in_same_landtype_group)

    if output_file.endswith('.csv'):
        df.to_csv(output_file, index=False)
    else:
        write_dataframe_to_fwf(output_file, df)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing {output_file}: {elapsed_time:.2f} seconds")


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file(s) on the command line.
    start_time = time.time()
    if len(sys.argv) < 2:
        print('Usage: python gcam_process_extracted_data.py `path/to/json/input/file(s)\'')
        sys.exit()

    # Read and load the JSON file(s) into a list of dictionaries.
    list_of_inputs = []
    for i in range(1, len(sys.argv)):
        input_file = sys.argv[i]
        with open(input_file) as f:
            list_of_inputs.extend(json.load(f))

    # Produce data for each file in parallel.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(process_extracted_data, list_of_inputs)
    
    # Print the total execution time needed to process/compile the scalars for all files.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all files: {elapsed_time:.2f} seconds")