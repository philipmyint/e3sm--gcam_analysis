import itertools
import json
import multiprocessing
import pandas as pd
import sys
import time
from utility_constants import *
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_fwf

def add_areas_to_subset_of_file(df, df_land, geographical_label, category_label, scenario, geography, category):
    """ 
    Adds areas from land allocation data contained in a Pandas DataFrame to another DataFrame containing data for the quantity of interest.
    This function adds areas to a subset of the DataFrame that matches the given scenario, geographical unit, and category.

    Parameters:
        df: DataFrame containing the data of interest.
        df_land: DataFrame for the land allocation areas.
        geographical_label: String specifying the label for the geographical unit (e.g., 'region' or 'basin').
        category_label: String specifying the label for the appropriate category (e.g., 'sector' or 'landtype').
        scenario: String for the scenario (simulation name) of interest.
        geography: String specifying the name of the geographical unit for the subset of interest (e.g., the name of the region or the basin).
        category: String specifying the name of the category for the subset of interest (e.g., the name of the sector or the landtype).

    Returns:
        DataFrame that is the same as the input df, but with an extra column for the corresponding land allocation areas.
    """
    # Read in the DataFrames for the data and the land allocations, get the appropriate subset for them, and put these subset DataFrames into a list.
    dataframes = [df, df_land]
    columns = {0: ['scenario', geographical_label, category_label], 1: ['scenario', geographical_label, 'landtype']}
    matches = [scenario, geography, category]
    for df_index, dataframe in enumerate(dataframes):
        for index, column in enumerate(columns[df_index]):
            dataframe = dataframe[dataframe[column] == matches[index]]
        dataframes[df_index] = dataframe

    df = dataframes[0]
    start_year, end_year = min(df['year'].unique()), max(df['year'].unique())
    df_land = dataframes[1]
    if df_land.empty:
        # If there are no corresponding land allocations, set the areas to 0.
        df['area'] = 0
    else:
        df_land = df_land[(df_land['year'] >= start_year) & (df_land['year'] <= end_year)]
        if geographical_label == 'region':
            # For each year, add up the areas of all rows (e.g., from all basins) that correspond to the given scenario, category, and region.
            df['area'] = df_land.groupby('year').sum()['value'].to_numpy()
        if geographical_label == 'basin':
            # If matching on the basin, there could be multiple regions that contain parts of this basin. As a result, for each year, 
            # add up the areas of all matching regions that correspond to the given scenario, category, and basin.
            regions = df['region'].unique()
            dataframes_for_this_region = []
            for region in regions:
                df_this_region = df[df['region'] == region].copy()
                df_land_this_region = df_land[df_land['region'] == region]
                if df_land_this_region.empty:
                    df_this_region['area'] = 0
                else:
                    df_this_region['area'] = df_land_this_region['value'].to_numpy()
                dataframes_for_this_region.append(df_this_region)
            # Concatenate the DataFrame from all regions into a single DataFrame for the given scenario, category, and basin.
            df = pd.concat(dataframes_for_this_region)
    return df

def add_areas_to_file(inputs):
    """ 
    Adds areas from the land allocation file as an extra column to the file specified in the inputs dictionary.

    Parameters:
        inputs: Dictionary with user-specified inputs for the name of the input and output files, the name of the land allocation file, etc.

    Returns:
        N/A.
    """
    # Unpack the inputs and read both the input file and the land allocation file into Pandas DataFrames.
    start_time = time.time()
    input_file = inputs['input_file']
    output_file = inputs['output_file']
    key_columns = inputs['key_columns']
    geographical_label = inputs['geographical_label']
    category_label = inputs.get('category_label', None)
    land_allocation_file = inputs['land_allocation_file']
    df = read_file_into_dataframe(input_file)
    df_land = read_file_into_dataframe(land_allocation_file)
    
    # Form a list of tuples that represents the Cartesian product of all scenarios, geographies, and categories, along with the DataFrames and labels.
    scenarios = df['scenario'].unique()
    geographies = df[geographical_label].unique()
    if category_label:
        categories = df[category_label].unique()
    else:
        categories = [None]
    cartesian_product = list(itertools.product([df], [df_land], [geographical_label], [category_label], 
                                               scenarios, geographies, categories))

    # Add areas for each subset of the data (each tuple in the Cartesian product) in parallel. Store the subsets in a list of DataFrames.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        dataframes_for_each_subset = list(pool.starmap(add_areas_to_subset_of_file, cartesian_product))

    # Concatenate all DataFrames in the list together to form a single DataFrame for this file. Sort by all the given key columns.
    df = pd.concat(dataframes_for_each_subset)
    df.sort_values(key_columns, inplace=True)

    if output_file.endswith('.csv'):
        df.to_csv(output_file, index=False)
    else:
        write_dataframe_to_fwf(output_file, df)
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
    
    # Add areas to each file sequentially.
    for inputs in list_of_inputs:
        add_areas_to_file(inputs)

    # Print the total execution time to add areas to all the files.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for adding areas to all files: {elapsed_time:.2f} seconds")
    