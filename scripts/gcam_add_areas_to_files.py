import json
import multiprocessing
import os
import pandas as pd
import sys
import time
from utility_constants import *
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_file
from utility_gcam import modify_crop_names

def add_areas_to_file(inputs):
    """ 
    Adds areas from the land allocation file as an extra column to the file specified in the inputs dictionary.
    Areas are matched to rows in the input file using a vectorized merge on scenario, geographical unit,
    category, and year — replacing a previous approach that used a multiprocessing pool over a Cartesian
    product of all (scenario, geography, category) combinations.

    Parameters:
        inputs: Dictionary with user-specified inputs for the name of the input and output files, the name of the land allocation file, etc.

    Returns:
        N/A.
    """
    # Unpack the inputs.
    start_time = time.time()
    input_file = inputs['input_file']
    output_file = inputs['output_file']
    key_columns = inputs['key_columns']
    geographical_label = inputs['geographical_label']
    category_label = inputs.get('category_label', None)
    land_allocation_file = inputs['land_allocation_file']
    mean_or_sum_if_more_than_one_row_in_same_landtype_group = inputs.get('mean_or_sum_if_more_than_one_row_in_same_landtype_group', None)
    call_modify_crop_names = inputs.get('call_modify_crop_names', False)

    # Check that the input and land allocation files exist, and that the output file is writable.
    if not os.path.exists(input_file):
        print(f"Error: input file not found: '{input_file}'")
        return
    if not os.path.exists(land_allocation_file):
        print(f"Error: land allocation file not found: '{land_allocation_file}'")
        return
    try:
        open(output_file, 'a').close()
    except OSError as e:
        print(f"Error: cannot write to output file '{output_file}': {e}")
        return

    # Read both the input file and the land allocation file into Pandas DataFrames.
    df = read_file_into_dataframe(input_file)
    df_land = read_file_into_dataframe(land_allocation_file)

    # Drop any existing 'area' column from df. Since input_file and output_file are often the same path,
    # the file may already contain an 'area' column from a previous run. If not dropped, the merge would
    # produce 'area_x' and 'area_y' columns instead of 'area', causing a KeyError downstream.
    if 'area' in df.columns:
        df = df.drop('area', axis=1)

    # Check that the required columns exist in both DataFrames before attempting the merge.
    for col in ['scenario', geographical_label, 'year']:
        if col not in df.columns:
            print(f"Error: '{col}' column not found in input file '{input_file}'. Available columns: {list(df.columns)}")
            return
    for col in ['scenario', geographical_label, 'year', 'value', 'landtype']:
        if col not in df_land.columns:
            print(f"Error: '{col}' column not found in land allocation file '{land_allocation_file}'. Available columns: {list(df_land.columns)}")
            return
    if category_label and category_label not in df.columns:
        print(f"Error: category_label '{category_label}' not found in input file '{input_file}'. Available columns: {list(df.columns)}")
        return

    # Add areas using a vectorized merge. This replaces the previous approach of creating a Cartesian product
    # of all (scenario, geography, category) combinations and processing each in a separate worker process.
    # The merge achieves the same result in a single pass without the overhead of pickling full DataFrames
    # to many workers.
    if geographical_label == 'region':
        # Sum land areas across all basins for each (scenario, region, landtype, year) to get the total
        # area per region — equivalent to what the previous code computed inside each worker.
        area_df = df_land.groupby(['scenario', 'region', 'landtype', 'year'])['value'].sum().reset_index()
        area_df = area_df.rename(columns={'value': 'area'})
        if category_label:
            area_df = area_df.rename(columns={'landtype': category_label})
            merge_cols = ['scenario', 'region', category_label, 'year']
        else:
            area_df = area_df.groupby(['scenario', 'region', 'year'])['area'].sum().reset_index()
            merge_cols = ['scenario', 'region', 'year']
        df = df.merge(area_df, on=merge_cols, how='left')

    elif geographical_label == 'basin':
        # For basin-level matching, areas are already per (scenario, region, basin, landtype, year),
        # so no groupby is needed — just rename and merge directly.
        area_df = df_land[['scenario', 'region', 'basin', 'landtype', 'year', 'value']].copy()
        area_df = area_df.rename(columns={'value': 'area'})
        if category_label:
            area_df = area_df.rename(columns={'landtype': category_label})
            merge_cols = ['scenario', 'region', 'basin', category_label, 'year']
        else:
            area_df = area_df.groupby(['scenario', 'region', 'basin', 'year'])['area'].sum().reset_index()
            merge_cols = ['scenario', 'region', 'basin', 'year']
        df = df.merge(area_df, on=merge_cols, how='left')

    # Fill any unmatched rows with 0, consistent with the previous approach which set area to 0
    # when no matching land allocation data was found.
    df['area'] = df['area'].fillna(0)
    df.sort_values(key_columns, inplace=True)

    # Update original crop names to a common, standardized set of names.
    if call_modify_crop_names:
        df = modify_crop_names(df, key_columns, mean_or_sum_if_more_than_one_row_in_same_landtype_group)

    write_dataframe_to_file(df, output_file)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for adding areas to {output_file}: {elapsed_time:.2f} seconds")


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file(s) on the command line.
    start_time = time.time()
    if len(sys.argv) < 2:
        print('Usage: python gcam_add_areas_to_files.py `path/to/json/input/file(s)\'')
        sys.exit()

    # Read and load the JSON file(s) into a list of dictionaries, where each dictionary corresponds to a particular file.
    list_of_inputs = []
    for i in range(1, len(sys.argv)):
        input_file = sys.argv[i]
        with open(input_file) as f:
            list_of_inputs.extend(json.load(f))

    # Add areas to each file in parallel.
    with multiprocessing.Pool(processes=MAX_PROCESSES) as pool:
        pool.map(add_areas_to_file, list_of_inputs)

    # Print the total execution time to add areas to all the files.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for adding areas to all files: {elapsed_time:.2f} seconds")
