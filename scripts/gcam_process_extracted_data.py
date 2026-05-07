import json
import multiprocessing
import os
import pandas as pd
import sys
import time
from utility_constants import MAX_PROCESSES
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_file
from utility_gcam import modify_crop_names

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
    call_modify_crop_names = inputs.get('call_modify_crop_names', False)
    
    # Verify input file exists and output file is writable before doing any processing.
    if not os.path.exists(input_file):
        print(f"Error: input file '{input_file}' does not exist.")
        return
    try:
        open(output_file, 'a').close()
    except OSError as e:
        print(f"Error: cannot write to output file '{output_file}': {e}")
        return

    df = read_file_into_dataframe(input_file)

    # Check that columns_to_drop columns actually exist before dropping, and warn about any that are missing.
    if columns_to_drop:
        missing_drop = [c for c in columns_to_drop if c not in df.columns]
        if missing_drop:
            print(f"Warning: the following columns_to_drop were not found in '{input_file}' and will be skipped: {missing_drop}")
        columns_to_drop_existing = [c for c in columns_to_drop if c in df.columns]
        df = df.drop(columns_to_drop_existing, axis=1)
    df.columns = df.columns.str.lower()

    # Check that the 'year' column exists after lowercasing before sorting by it.
    if 'year' not in df.columns:
        print(f"Error: no 'year' column found in '{input_file}' after lowercasing. Available columns: {list(df.columns)}")
        return
    df = df.sort_values(by='year')

    # Split each specified column (if any) separated by underscores into a set of new columns.
    if columns_to_split:
        for column_to_split, new_columns in columns_to_split.items():
            # Check that the column to split exists after lowercasing.
            if column_to_split not in df.columns:
                print(f"Error: columns_to_split column '{column_to_split}' not found in '{input_file}' after lowercasing. "
                      f"Available columns: {list(df.columns)}")
                return
            # Check that the split produces enough columns for the requested new column names.
            split_result = df[column_to_split].str.split('_', expand=True)
            if split_result.shape[1] < len(new_columns):
                print(f"Error: splitting '{column_to_split}' on '_' produced {split_result.shape[1]} column(s) "
                      f"but {len(new_columns)} new column(s) were requested: {new_columns}.")
                return
            df[new_columns] = split_result.iloc[:, :len(new_columns)]

    # Sort columns by a set of keys if specified to do so.
    if key_columns:
        # Check that all key_columns exist (including any newly created columns from columns_to_split).
        missing_keys = [c for c in key_columns if c not in df.columns]
        if missing_keys:
            print(f"Error: the following key_columns were not found in '{input_file}' after all processing steps: {missing_keys}. "
                  f"Available columns: {list(df.columns)}")
            return
        # Check that a 'value' column exists before selecting key_columns + ['value'].
        if 'value' not in df.columns:
            print(f"Error: no 'value' column found in '{input_file}'. Available columns: {list(df.columns)}")
            return
        df = df.sort_values(by=key_columns)
        df = df[key_columns + ['value']]

    # Perform either a sum or a mean if there happens to be more than one row that agrees on all the table keys (i.e., belongs to the same group).
    if mean_or_sum_if_more_than_one_row_in_same_landtype_group == 'sum' and key_columns:
        df = df.groupby(key_columns).sum().reset_index()
    elif mean_or_sum_if_more_than_one_row_in_same_landtype_group == 'mean' and key_columns:
        df = df.groupby(key_columns).mean().reset_index()

    # Update original crop names to a common, standardized set of names. 
    if call_modify_crop_names:
        df = modify_crop_names(df, key_columns, mean_or_sum_if_more_than_one_row_in_same_landtype_group)

    write_dataframe_to_file(df, output_file)
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
    with multiprocessing.Pool(processes=MAX_PROCESSES) as pool:
        pool.map(process_extracted_data, list_of_inputs)
    
    # Print the total execution time needed to process/compile the scalars for all files.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all files: {elapsed_time:.2f} seconds")