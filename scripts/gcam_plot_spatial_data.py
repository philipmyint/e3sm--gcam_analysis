import geopandas as gpd
import json
from matplotlib import pyplot as plt
import multiprocessing
import os
import pandas as pd
from scipy import stats
import sys
import time
from utility_constants import *
from utility_dataframes import perform_ttest, read_file_into_dataframe
from utility_functions import check_is_list_of_lists, print_p_values, sort_file
from utility_gcam import *
from utility_plots import *

""" Dictionary of default input values for spatial plots. """
default_inputs_time_series = {
    'basin_label': 'basin',
    'category_label': 'sector',
    'cbar_limits': None,
    'cbar_on': True,
    'cmap': 'viridis',
    'end_year': 2090,
    'height': height_default,
    'key_columns': None,
    'landtype_groups': 'standard',
    'linewidth': 0.5,
    'multiplier': 1,
    'mean_or_sum_if_more_than_one_row_in_same_landtype_group': 'area_weighted_mean',
    'mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin': 'mean',
    'p_value_file': 'p_values.dat',
    'p_value_file_print_only_if_below_threshold': True,
    'p_value_threshold': 0.05,
    'plot_directory': './',
    'plot_type': 'absolute_difference',
    'produce_png': produce_png_default, 
    'region_label': 'region', 
    'scenario_label': 'scenario', 
    'scenario_sets': None,
    'shape_file_basin_label': 'glu_nm',
    'shape_file_region_label': 'reg_nm',
    'start_year': 2070,
    'stippling_hatches': 'xxxx',
    'stippling_on': True,
    'time_calculation': 'mean',
    'use_latex': use_latex_default, 
    'value_label': 'value',
    'width': width_default,   
    'x_tick_label_size': tick_label_size_default,   
    'y_tick_label_size': tick_label_size_default,
    'year_label': 'year'            
}

def process_inputs(inputs):
    """ 
    Processes a dictionary of inputs (keys are options, values are choices for those options) for creating spatial plots.

    Parameters:
        inputs: Dictionary containing the user plotting choice inputs for different options. This dictionary may be incomplete or have invalid values.

    Returns:
        Dictionary that completely specifies all plotting options. 
        If the user did not select a plotting option fora particular category, the default choice for that plotting option will be selected.
    """
    df = read_file_into_dataframe(inputs['output_file'])

    if 'category_label' not in inputs:
        category_label = default_inputs_time_series['category_label']
        inputs['category_label'] = category_label
    else:
        category_label = inputs['category_label']

    # If the list of categories (e.g., sectors or landtypes) have not been specified, populate the list with all categories except the excluded ones.
    if 'categories' not in inputs:
        if category_label in df.columns:
            inputs['categories'] = df[category_label].unique()
            if 'categories_to_exclude' in inputs:
                inputs['categories'] = [column for column in inputs['categories'] if column not in inputs['categories_to_exclude']]
        else:
            inputs['categories'] = 'All'
    categories = inputs['categories']
    # If the user entered a string indicating a single category, put that string in a list.
    if isinstance(categories, str):
        categories = [categories]
        inputs['categories'] = categories

    # Create the plot directory if it does not already exist. By default, put the name of the file containing p-value results in this directory.
    if 'plot_directory' not in inputs:
        inputs['plot_directory'] = default_inputs_time_series['plot_directory']
    if not os.path.exists(inputs['plot_directory']):
        os.makedirs(inputs['plot_directory'])
    if 'p_value_file' not in inputs:
        inputs['p_value_file'] = os.path.join(inputs['plot_directory'], default_inputs_time_series['p_value_file'])

    # Use the name of the output file itself (without its path) to set defaults for the title and the name of the plot.
    index_of_last_backslash = inputs['output_file'].rfind('/')
    index_of_dot_csv = inputs['output_file'].find('.csv')
    if index_of_last_backslash == -1:
        output_file_name = inputs['output_file'][:index_of_dot_csv]
    else:
        output_file_name = inputs['output_file'][index_of_last_backslash+1:index_of_dot_csv]
    if 'title' not in inputs:
        inputs['title'] = output_file_name
    if 'plot_name' not in inputs:
        inputs['plot_name'] = os.path.join(inputs['plot_directory'], 'spatial_' + output_file_name)

    # Add keys for plotting options that have not been specified in the inputs dictionary and use default values for them.
    for key in default_inputs_time_series.keys():
        if key not in inputs:
            inputs[key] = default_inputs_time_series[key]

    # If the scenarios have not been specified, use all the scenarios in the Pandas DataFrame.
    if 'scenarios' not in inputs:
        scenario_label = inputs['scenario_label']
        inputs['scenarios'] = df[scenario_label].unique()
    
    return inputs

def plot_spatial_data(inputs):
    """ 
    Creates a spatial plot and perform statistical analysis for a single output file. The data in the file are organized
    into scenarios or scenario sets, categories, and regions.

    Parameters:
        input: Dictionary containing user inputs for different plotting options, where the keys are options and values are choices for those options.
               This dictionary is assumed to be complete (pre-processed).

    Returns:
        N/A.
    """
    # This function creates a spatial plot where the data (scenarios or scenario sets, categories, regions) all come from a single output file.
    start_time = time.time()
    output_file = inputs['output_file']
    
    # Extract all other plotting options.
    basin_label = inputs['basin_label']
    categories = inputs['categories']
    category_label = inputs['category_label']
    cbar_limits = inputs['cbar_limits']
    cbar_on = inputs['cbar_on']
    cmap_color = inputs['cmap']
    end_year = inputs['end_year']
    height = inputs['height']
    key_columns = inputs['key_columns']
    landtype_groups = inputs['landtype_groups']
    linewidth = inputs['linewidth']
    multiplier = inputs['multiplier']
    mean_or_sum_if_more_than_one_row_in_same_landtype_group = inputs['mean_or_sum_if_more_than_one_row_in_same_landtype_group']
    mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin = inputs['mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin']
    p_value_file = inputs['p_value_file']
    p_value_file_print_only_if_below_threshold = inputs['p_value_file_print_only_if_below_threshold']
    p_value_threshold = inputs['p_value_threshold']
    plot_name = inputs['plot_name']
    plot_type = inputs['plot_type']
    produce_png = inputs['produce_png']
    region_label = inputs['region_label']
    scenario_label = inputs['scenario_label']
    scenarios = inputs['scenarios']
    shape_file = inputs['shape_file']
    shape_file_basin_label = inputs['shape_file_basin_label']
    shape_file_region_label = inputs['shape_file_region_label']
    start_year = inputs['start_year']
    stippling_hatches = inputs['stippling_hatches']
    stippling_on = inputs['stippling_on']
    time_calculation = inputs['time_calculation']
    title = inputs['title']
    use_latex = inputs['use_latex']
    value_label = inputs['value_label']
    width = inputs['width']
    x_tick_label_size = inputs['x_tick_label_size']
    y_tick_label_size = inputs['y_tick_label_size']
    year_label = inputs['year_label']

    # Set the plotting options.
    if plot_type == 'percent_difference' and title and (f'%' not in title or 'percent' not in title):
        title += rf' (\% difference)'
    plot_options = dict(width=width, height=height, name=plot_name, produce_png=produce_png)
    plot_options.update(zip(['x_tick_label_size', 'y_tick_label_size', 'use_latex'], [x_tick_label_size, y_tick_label_size, use_latex]))

    # Use LaTeX fonts for figures and set font size of tick labels.
    if use_latex:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', weight='bold')

    # Read the data file into a Pandas DataFrame and select rows between the start and end years.
    df = read_file_into_dataframe(output_file)
    df = df[(df[year_label] >= start_year) & (df[year_label] <= end_year)]
    df[value_label] *= multiplier

    # Set the appropriate dictionary for the landtype group.
    if landtype_groups == 'standard':
        landtype_groups = gcam_landtype_groups
    elif landtype_groups == 'nonstandard':
        landtype_groups = gcam_landtype_groups_nonstandard

    # Create separate DataFrames for the following cases: 1) all categories (create a copy of the entire DataFrame), 2) categories that correspond to
    # a group of landtypes (e.g., forest, crop, grass, shrub, pasture), and 3) a set of individual categories. Concatenate into a single DataFrame.
    dataframes = []
    if 'All' in categories:
        df_all = df.copy()
        df_all[category_label] = 'All'
        dataframes.append(df_all) 
        categories.remove('All')
    if any(item in landtype_groups for item in categories):
        for landtype_group in landtype_groups:
            if landtype_group in categories:
                df_this_landtype_group = produce_dataframe_for_landtype_group(df, landtype_group, category_label, 
                    value_label, landtype_groups, mean_or_sum_if_more_than_one_row_in_same_landtype_group, key_columns)
                dataframes.append(df_this_landtype_group)
                categories.remove(landtype_group)
    if categories:
        df = df[df[category_label].isin(categories)]
        dataframes.append(df)
    df = pd.concat(dataframes).reset_index()

    # Read the shape file into a GeoDataFrame from the GeoPandas library and get all regions in the file.
    gdf = gpd.read_file(shape_file)
    regions = gdf[shape_file_region_label].unique()

    # Create a common variable called scenario_list that will contain a 1D list of all scenarios for both the direct plots and ensemble plots.
    if not check_is_list_of_lists(scenarios):
        # Direct plots, in which each such spatial plot can include one or more individual (not grouped) data sets and the scenarios is already a 1D list.
        scenario_list = scenarios
        # Column labels for these scenarios.
        scenario_columns = [f'scen={i}' for i in range(len(scenarios))]
    else:
        # Ensemble plots, in which the ensemble is further subdivided into groups, where each group represents a set of scenarios.
        # The scenarios in this case are contained in a list of lists (i.e., a 2D list).
        scenario_list = []
        scenario_columns = []
        num_scenario_sets = len(scenarios[0])
        num_scenarios_in_each_set = len(scenarios)
        # Group all scenarios that belong to the same set (they share the same column in the 2D list) and put them together in the 1D scenario_list.
        for index in range(num_scenario_sets):
            scenarios_in_set = [scenarios[i][index] for i in range(num_scenarios_in_each_set)]
            scenario_list.extend(scenarios_in_set)
            scenarios_indices_this_set = [f'scen={index}_{i}' for i in range(num_scenarios_in_each_set)]
            scenario_columns.extend(scenarios_indices_this_set)

    # For each scenario, region, and basin, calculate the mean or sum of all relevant categories over all years and put into the GeoDataFrame.
    for index, scenario in enumerate(scenario_list):
        df_this_scenario = df[df[scenario_label] == scenario]
        for region in regions:
            region_filter = gdf[shape_file_region_label] == region
            if basin_label in df_this_scenario.columns:
                # If the DataFrame containing the data of interest are organized into basins, add basin-level information.
                basins = gdf[region_filter][shape_file_basin_label].unique()
                for basin in basins:
                    basin_filter = gdf[shape_file_basin_label] == basin
                    if basin not in gcam_basin_names_and_abbrevations:
                        # Set the value to 0 if there is no data for the current region-and-basin combination.
                        gdf.loc[region_filter & basin_filter, scenario_columns[index]] = 0
                    else:
                        # Basin names are fully written out in the GeoDataFrame (shape file), while they are abbreviated in the DataFrame (data file).
                        basin_abbrv = gcam_basin_names_and_abbrevations[basin]
                        # Calculate the mean or sum of all relevant categories over all years for the current region in the for-loop iteration.
                        if mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin == 'mean':
                            df_grouped_by_year = \
                            df_this_scenario[(df_this_scenario[region_label] == region) & \
                                        (df_this_scenario[basin_label] == basin_abbrv)].groupby('year')[value_label].mean()
                        elif mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin == 'sum':
                            df_grouped_by_year = \
                            df_this_scenario[(df_this_scenario[region_label] == region) & \
                                        (df_this_scenario[basin_label] == basin_abbrv)].groupby('year')[value_label].sum()
                        if time_calculation == 'mean':
                            gdf.loc[region_filter & basin_filter, scenario_columns[index]] = df_grouped_by_year.mean()
                        elif time_calculation == 'sum':
                            gdf.loc[region_filter & basin_filter, scenario_columns[index]] = df_grouped_by_year.sum()
            else:
                # If the DataFrame containing the data of interest does not provide basin-level information (only regions), calculate the mean or sum
                # of all relevant categories over all years for the current region in the for-loop iteration. On the spatial plot, the displayed
                # quantity will take on a uniform value over the entire region (over all basins that encompass that region).
                if mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin == 'mean':
                    df_grouped_by_year = df_this_scenario[df_this_scenario[region_label] == region].groupby('year')[value_label].mean()
                elif mean_or_sum_if_more_than_one_row_in_same_region_and_or_basin == 'sum':
                    df_grouped_by_year = df_this_scenario[df_this_scenario[region_label] == region].groupby('year')[value_label].sum()
                if time_calculation == 'mean':
                    gdf.loc[region_filter, scenario_columns[index]] = df_grouped_by_year.mean()
                elif time_calculation == 'sum':
                    gdf.loc[region_filter, scenario_columns[index]] = df_grouped_by_year.sum()
    gdf.fillna(0, inplace=True)

    if len(scenario_list) == 1:
        # If there is only one scenario in the list, simply plot that column.
        gdf['plot'] = gdf[scenario_columns[0]]
    elif plot_type == 'mean':
        # If there are multiple data columns, one option is to plot the mean across all columns.
        gdf['plot'] = gdf[scenario_columns].mean(axis=1)
    elif plot_type == 'absolute_difference' or plot_type == 'percent_difference':
        if not check_is_list_of_lists(scenarios):
            # If a direct plot, assume the first column is the reference for the absolute difference or percent difference calculation.
            columns_control = [scenario_columns[0]]
            columns_test = [x for x in scenario_columns[1:]] 
        else:
            # If an ensemble plot, group the scenarios into sets and calculate set means. Assume the first scenario set is the reference set.
            columns_control = [x for x in scenario_columns if x.startswith('scen=0')]
            columns_test = [x for x in scenario_columns if not x.startswith('scen=0')]
            if check_is_list_of_lists(scenarios) and num_scenario_sets == 2:
                # If there are two scenario sets, do a t-test at each individual region-and-basin combination and put results into the GeoDataFrame.
                df = pd.DataFrame()
                df[columns_control], df[columns_test] = gdf[columns_control], gdf[columns_test]
                gdf['p_value'] = df.apply(perform_ttest, columns_set_1=columns_control, \
                                                               columns_set_2=columns_test, axis=1).fillna(1)
        control_data = gdf.loc[:, columns_control].mean(axis=1)
        test_data = gdf.loc[:, columns_test].mean(axis=1)
        if plot_type == 'percent_difference':
            gdf['plot'] = (test_data - control_data)/(control_data + EPSILON)*100
        else:
            # Absolute difference.
            gdf['plot'] = test_data - control_data
        # Do a t-test for the test and control data sets as a whole (global comparison). Print the output from this into a data file.
        ttest = stats.ttest_ind(control_data, test_data)
        print_p_values(ttest, columns_test, p_value_threshold, p_value_file, plot_name, p_value_file_print_only_if_below_threshold)

    # Create figure and axis objects for the plot, set the title and colorbar limits.
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.set_title(title)
    if cbar_limits:
        vmin, vmax = cbar_limits[0], cbar_limits[1]
    else:
        vmin, vmax = None, None
    # Generate the plot and optionally add stippling to indicate statistically significant differences at individual regions and/or basins.
    gdf.plot('plot', ax=ax, legend=cbar_on, cmap=cmap_color, vmin=vmin, vmax=vmax, legend_kwds={"shrink": .5}, edgecolor='k', linewidth=linewidth)
    if 'p_value' in gdf.columns and stippling_on:
        mask = gdf['p_value'] <= p_value_threshold
        gdf[mask].plot(ax=ax, facecolor='none', color='none', hatch=stippling_hatches, linewidth=0)
    
    # Finalize the spatial plot now that all data sets have been processed.
    plot_options['name'] = plot_name
    fig.set_size_inches(width, height)
    save_figure(plot_name, fig, plot_options)

    # Close the plot now that we are done with it. Record the elapsed time.
    plt.close(fig)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing plot {plot_name}: {elapsed_time:.2f} seconds")


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file(s) on the command line.
    start_time = time.time()
    if len(sys.argv) < 2:
        print('Usage: python gcam_plot_spatial_data.py `path/to/json/input/file(s)\'')
        sys.exit()

    # Read and load the JSON file(s) into a list of dictionaries.
    inputs = []
    for i in range(1, len(sys.argv)):
        input_file = sys.argv[i]
        with open(input_file) as f:
            inputs.extend(json.load(f))

    # Process each dictionary so that each of them specifies a complete set of options (e.g., by adding default values) for a single plot.
    start_time = time.time()
    list_of_inputs = []
    for index in range(len(inputs)):
        list_of_inputs.append(process_inputs(inputs[index]))

    # Delete all the p-value files before we do any calculations to start a fresh run.
    for inputs in list_of_inputs:
        file = inputs['p_value_file']
        if os.path.exists(file): 
            os.remove(file)

    # Create all of the bpx plots in parallel.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(plot_spatial_data, list_of_inputs)
    
    # Sort all the p-value files alphabetically.
    for inputs in list_of_inputs:
        file = inputs['p_value_file']
        if os.path.exists(file): 
            sort_file(file)
    
    # Print the total execution time to produce all the plots.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all plots: {elapsed_time:.2f} seconds")