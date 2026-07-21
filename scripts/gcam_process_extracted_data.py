import json
import multiprocessing
import os
import pandas as pd
import sys
import time
from utility_constants import MAX_PROCESSES
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_file
from utility_gcam import modify_crop_names, add_areas_to_dataframe

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
    key_columns = [c.lower() for c in inputs['key_columns']] if inputs.get('key_columns') else None
    aggregate_function = inputs.get('aggregate_function', "sum") 
    call_modify_crop_names = inputs.get('call_modify_crop_names', False)
    land_allocation_file = inputs.get('land_allocation_file', None)
    
    # Verify input file exists and output file is writable before doing any processing.
    if not os.path.exists(input_file):
        print(f"Error: input file '{input_file}' does not exist.")
        return
    if os.path.getsize(input_file) == 0:
        print(f"Warning: input file '{input_file}' is empty. Skipping.")
        return
    try:
        open(output_file, 'a').close()
    except OSError as e:
        print(f"Error: cannot write to output file '{output_file}': {e}")
        return

    df = read_file_into_dataframe(input_file)
    if df.empty or len(df.columns) == 0:
        print(f"Warning: input file '{input_file}' contains no data. Skipping.")
        return

    # Check that columns_to_drop columns actually exist before dropping, and warn about any that are missing.
    if columns_to_drop:
        missing_drop = [c for c in columns_to_drop if c not in df.columns]
        if missing_drop:
            print(f"Warning: the following columns_to_drop were not found in '{input_file}' and will be skipped: {missing_drop}")
        columns_to_drop_existing = [c for c in columns_to_drop if c in df.columns]
        df = df.drop(columns_to_drop_existing, axis=1)
    
    # Lowercase all column names to avoid case sensitivity issues.
    df.columns = df.columns.str.lower()

    # Check that the 'year' column exists after lowercasing before sorting by it.
    if 'year' not in df.columns:
        print(f"Error: no 'year' column found in '{input_file}' after lowercasing. Available columns: {list(df.columns)}")
        return
    df = df.sort_values(by='year')

    # For carbon-density files, replace 'None Specified' in the 'units' column with 'MgC/ha'.
    if 'carbon-density' in input_file and 'units' in df.columns:
        df.loc[df['units'] == 'None Specified', 'units'] = 'MgC/ha'
    
    # For scaler files, replace NA in the 'units' column with 'Unitless'.
    if 'scaler' in input_file and 'units' in df.columns:
        df['units'] = df['units'].fillna('Unitless')

    # Split each specified column (if any) separated by underscores into a set of new columns.
    # Rows with no '_' have the full value placed in the first new column and remaining new
    # columns set to "". Rows with '_' are split across the new columns; extra new columns
    # beyond the number of split elements are set to "". It is an error if any row produces
    # more split elements than there are new columns.
    if columns_to_split:
        for column_to_split, new_columns in columns_to_split.items():
            # Check that the column to split exists after lowercasing.
            if column_to_split not in df.columns:
                print(f"Error: columns_to_split column '{column_to_split}' not found in '{input_file}' after lowercasing. "
                      f"Available columns: {list(df.columns)}")
                return
            mask = df[column_to_split].str.contains('_', na=False)
            # Initialize all new columns to "".
            for col_name in new_columns:
                df[col_name] = ""
            if mask.any():
                split_result = df.loc[mask, column_to_split].str.split('_', expand=True)
                if split_result.shape[1] > len(new_columns):
                    print(f"Error: columns_to_split column '{column_to_split}' in '{input_file}' produces "
                          f"{split_result.shape[1]} split elements but only {len(new_columns)} new columns were specified.")
                    return
                for i, col_name in enumerate(new_columns):
                    if i < split_result.shape[1]:
                        df.loc[mask, col_name] = split_result.iloc[:, i].fillna('')
            # Rows without '_': the full value goes into the first new column.
            df.loc[~mask, new_columns[0]] = df.loc[~mask, column_to_split]

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

        # Add areas from the processed land allocation file.
        # Skip this for land allocation data itself, since those files ARE the area data.
        if 'land_allocation' not in input_file:
            if land_allocation_file is None:
                print(f"Error: A 'land_allocation_file' file must be specified in the JSON input")
                print(f"Error: So first the 'detailed_land_allocation_processed.csv' must be created by running "
                      f"the 'gcam_process_extracted_data.py' script on the 'detailed_land_allocation.csv' file."
                      f"without modifying the crop names (call_modify_crop_names = false).")
                print(f"Error: The default 'land_allocation_file' is called 'detailed_land_allocation_processed.csv' "
                      f"and is should be located in the 'output/gcam_csv/' directory.")
                return
            if not os.path.exists(land_allocation_file):
                print(f"Error: land_allocation_file not found: '{land_allocation_file}'")
                return
            df_land = read_file_into_dataframe(land_allocation_file)
            df = add_areas_to_dataframe(df, df_land)

    # Perform either a sum (default) or a weighted mean if there happens to be more than one row that agrees on all the table keys
    #    (i.e., belongs to the same group).
    # This also works for means where the areas are the same, e.g., groups that are in the same area unit,
    #  in which case the weighted mean is equivalent to a simple mean.
    if aggregate_function == 'sum' and key_columns:
        # Save string non-key columns (e.g. area_units) to rejoin after aggregation,
        # since groupby().sum() cannot aggregate string columns.
        str_non_key = [c for c in df.columns if c not in key_columns and df[c].dtype == object]
        str_first = df.groupby(key_columns)[str_non_key].first().reset_index() if str_non_key else None
        df = df.drop(columns=str_non_key).groupby(key_columns).sum().reset_index()
        if str_first is not None:
            df = df.merge(str_first, on=key_columns, how='left')
    elif aggregate_function == 'area_weighted_mean' and key_columns:
        # Save string non-key columns (e.g. area_units) to rejoin after aggregation.
        str_non_key = [c for c in df.columns if c not in key_columns and df[c].dtype == object]
        str_first = df.groupby(key_columns)[str_non_key].first().reset_index() if str_non_key else None
        df_numeric = df.drop(columns=str_non_key)
        # Compute area-weighted mean: for each group, sum(value * area) / sum(area).
        df_numeric['weighted_value'] = df_numeric['value'] * df_numeric['area']
        grouped = df_numeric.groupby(key_columns)[['weighted_value', 'area']].sum().reset_index()
        grouped['value'] = grouped['weighted_value'] / grouped['area'].replace(0, float('nan'))
        df = grouped.drop('weighted_value', axis=1)
        if str_first is not None:
            df = df.merge(str_first, on=key_columns, how='left')

    # Update original crop names to a common, standardized set of names. 
    if call_modify_crop_names:
        df = modify_crop_names(df, key_columns, aggregate_function)

    # Reorder columns so the last four are: value, units, area, area_units.
    trailing = [c for c in ['value', 'units', 'area', 'area_units'] if c in df.columns]
    other_cols = [c for c in df.columns if c not in trailing]
    df = df[other_cols + trailing]

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

    # Split inputs: land allocation blocks must run first (serially) because other blocks
    # depend on their output as the land_allocation_file for area-weighted calculations.
    land_alloc_inputs = [i for i in list_of_inputs if 'land_allocation' in i['input_file']]
    other_inputs = [i for i in list_of_inputs if 'land_allocation' not in i['input_file']]

    for inputs in land_alloc_inputs:
        process_extracted_data(inputs)

    # Produce the remaining data in parallel.
    with multiprocessing.Pool(processes=MAX_PROCESSES) as pool:
        pool.map(process_extracted_data, other_inputs)
    
    # Print the total execution time needed to process/compile the scalars for all files.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all files: {elapsed_time:.2f} seconds")