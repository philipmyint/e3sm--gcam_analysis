import itertools
import json
from matplotlib import pyplot as plt
import multiprocessing
import os
import pandas as pd
import seaborn as sns
import sys
import time
from utility_constants import *
from utility_dataframes import read_file_into_dataframe
from utility_functions import check_is_list_of_lists
from utility_gcam import gcam_landtype_groups, gcam_landtype_groups_nonstandard, produce_dataframe_for_landtype_group
from utility_plots import *

""" Dictionary of default input values for box plots (i.e., box-and-whisker plots). """
default_inputs_time_series = {
    'basin_label': 'basin',
    'basins': None,
    'category_label': 'sector',
    'end_year': 2100,
    'fill_boxes': True,
    'height': height_default,
    'hue': None,
    'key_columns': None,
    'landtype_groups': 'standard',
    'legend_x_offset': None,
    'legend_label_size': legend_label_size_default,
    'legend_num_columns': 1,
    'legend_on': legend_on_default,
    'legend_place_outside': False, 
    'linewidth': 1,
    'marker_size': 6,
    'multiplier': 1,
    'mean_or_sum_if_more_than_one_row_in_same_landtype_group': 'area_weighted_mean',
    'plot_colors': plot_colors_default,
    'plot_directory': './',
    'plot_percent_difference': False,
    'plot_type': 'ensemble_averages',
    'produce_png': produce_png_default, 
    'region_label': 'region', 
    'regions': ['Global'],
    'scenario_label': 'scenario', 
    'scenario_sets': None,
    'start_year': 2015,
    'use_latex': use_latex_default, 
    'value_label': 'value',
    'width': width_default,   
    'x_label': None,
    'x_label_size': axis_label_size_default,
    'x_scale': None,           
    'x_tick_label_size': tick_label_size_default,   
    'y_label_size': axis_label_size_default,
    'y_limits': axis_limits_default,
    'y_scale': scale_default,
    'y_tick_label_size': tick_label_size_default,
    'year_label': 'year'            
}

def process_inputs(inputs):
    """ 
    Processes a dictionary of inputs (keys are options, values are choices for those options) for creating box (or box-and-whisker) plots.

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

    # Create the plot directory if it does not already exist.
    if 'plot_directory' not in inputs:
        inputs['plot_directory'] = default_inputs_time_series['plot_directory']
    if not os.path.exists(inputs['plot_directory']):
        os.makedirs(inputs['plot_directory'])

    # Use the name of the output file itself (without its path) to set defaults for the y-axis label and the name of the plot.
    index_of_last_backslash = inputs['output_file'].rfind('/')
    index_of_dot_csv = inputs['output_file'].find('.csv')
    if index_of_last_backslash == -1:
        output_file_name = inputs['output_file'][:index_of_dot_csv]
    else:
        output_file_name = inputs['output_file'][index_of_last_backslash+1:index_of_dot_csv]
    if 'y_label' not in inputs:
        inputs['y_label'] = output_file_name
    if 'plot_name' not in inputs:
        inputs['plot_name'] = os.path.join(inputs['plot_directory'], 'box_plot_' + output_file_name)

    # Add keys for plotting options that have not been specified in the inputs dictionary and use default values for them.
    for key in default_inputs_time_series.keys():
        if key not in inputs:
            inputs[key] = default_inputs_time_series[key]

    # If the scenarios have not been specified, use all the scenarios in the Pandas DataFrame.
    if 'scenarios' not in inputs:
        scenario_label = inputs['scenario_label']
        inputs['scenarios'] = df[scenario_label].unique()
    
    return inputs

def plot_box_and_whiskers(inputs):
    """ 
    Creates a box plot (box-and-whisker plot) and perform statistical analysis for a single output file. The data in the file are organized
    into scenarios or scenario sets, categories, and regions.
    Types of plots: 1) Direct plots in which each scenario is treated as an individual data set in the collection. 
    2) Ensemble plots in which the scenarios are grouped into sets and results from of each set are plotted. The scenarios must be specified as a 
    list of lists for ensemble plots.

    Parameters:
        input: Dictionary containing user inputs for different plotting options, where the keys are options and values are choices for those options.
               This dictionary is assumed to be complete (pre-processed).

    Returns:
        N/A.
    """
    # This function creates a box plot where the data (scenarios or scenario sets, categories, regions) all come from a single output file.
    start_time = time.time()
    output_file = inputs['output_file']
    
    # Extract all other plotting options.
    basin_label = inputs['basin_label']
    basins = inputs['basins']
    categories = inputs['categories']
    category_label = inputs['category_label']
    end_year = inputs['end_year']
    fill_boxes = inputs['fill_boxes']
    height = inputs['height']
    hue = inputs['hue']
    key_columns = inputs['key_columns']
    landtype_groups = inputs['landtype_groups']
    legend_x_offset = inputs['legend_x_offset']
    legend_label_size = inputs['legend_label_size']
    legend_num_columns = inputs['legend_num_columns']
    legend_on = inputs['legend_on']
    legend_place_outside = inputs['legend_place_outside']
    linewidth = inputs['linewidth']
    marker_size = inputs['marker_size']
    multiplier = inputs['multiplier']
    mean_or_sum_if_more_than_one_row_in_same_landtype_group = inputs['mean_or_sum_if_more_than_one_row_in_same_landtype_group']
    plot_colors = inputs['plot_colors']
    plot_name = inputs['plot_name']
    plot_percent_difference = inputs['plot_percent_difference']
    plot_type = inputs['plot_type']
    produce_png = inputs['produce_png']
    region_label = inputs['region_label']
    regions = inputs['regions']
    scenario_label = inputs['scenario_label']
    scenario_sets = inputs['scenario_sets']
    scenarios = inputs['scenarios']
    start_year = inputs['start_year']
    use_latex = inputs['use_latex']
    value_label = inputs['value_label']
    width = inputs['width']
    x_label = inputs['x_label']
    x_label_size = inputs['x_label_size']
    x_scale = inputs['x_scale']
    x_tick_label_size = inputs['x_tick_label_size']
    x_variable = inputs.get('x_variable', category_label)
    y_label = inputs['y_label']
    y_label_size = inputs['y_label_size']
    y_limits = inputs['y_limits']
    y_scale = inputs['y_scale']
    y_tick_label_size = inputs['y_tick_label_size']
    year_label = inputs['year_label']

    # Set the plotting options.
    if plot_percent_difference and (f'%' not in y_label or 'percent' not in y_label):
        y_label += rf' (\% difference)'
    plot_options = dict(width=width, height=height, name=plot_name, x_label=x_label, y_label=fr'{y_label}')
    plot_options.update(zip(['x_scale', 'y_scale', 'y_limits', 'use_latex'], [x_scale, y_scale, y_limits, use_latex]))
    plot_options.update(zip(['x_tick_label_size', 'y_tick_label_size', 'legend_on'], [x_tick_label_size, y_tick_label_size, legend_on]))
    plot_options.update(zip(['x_label_size', 'y_label_size', 'legend_label_size'], [x_label_size, y_label_size, legend_label_size]))
    plot_options.update(zip(['produce_png', 'legend_x_offset', 'legend_place_outside'], [produce_png, legend_x_offset, legend_place_outside]))
    plot_options.update(zip(['legend_num_columns'], [legend_num_columns]))

    # Use LaTeX fonts for figures and set font size of tick labels.
    if use_latex:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', weight='bold')
    plt.rcParams['xtick.labelsize'] = x_tick_label_size
    plt.rcParams['ytick.labelsize'] = y_tick_label_size

    # Read the file, select rows between the start and end years, create the figure and axis objects for the plot.
    df = read_file_into_dataframe(output_file)
    df = df[(df[year_label] >= start_year) & (df[year_label] <= end_year)]
    df[value_label] *= multiplier
    fig, ax = plt.subplots(nrows=1, ncols=1)

    # Create and concatenate Pandas DataFrames for the regions and basins; 'Global' involving creating a copy of the entire DataFrame.
    dataframes = []
    if 'Global' in regions:
        df_global = df.copy()
        df_global[region_label] = 'Global'
        dataframes.append(df_global) 
        regions.remove('Global')
    if regions:
        df_these_regions = df[df[region_label].isin(regions)]
        dataframes.append(df_these_regions) 
    if basins:
        df_these_basins = df[df[basin_label].isin(basins)]
        dataframes.append(df_these_basins) 
    df = pd.concat(dataframes)

    # Direct plots, in which each such box plot can include one or more individual (not grouped) data sets.
    if not check_is_list_of_lists(scenarios) or plot_type == 'direct':
        
        # If scenarios is a list of lists, but we want a direct plot, unpack scenarios into a list of strings, with one string for each output file.
        if check_is_list_of_lists(scenarios) and plot_type == 'direct':
            scenarios = list(itertools.chain.from_iterable(scenarios))
            # Turn the corresponding plot_colors to list of lists with enough rows to match, then unpack them.
            plot_colors = list(itertools.chain.from_iterable(plot_colors))

        # Select the scenarios of interest.
        df = df[df[scenario_label].isin(scenarios)]
        scenario_list = scenarios
    
    # Ensemble plots, in which the ensemble is further subdivided into groups, where each group represents a set of scenarios.
    elif check_is_list_of_lists(scenarios) and plot_type == 'ensemble_averages':
        # For each scenario set, collect all scenarios that belong to that set and put in a DataFrame, then concatenate the DataFrames for all sets.
        dataframes = []
        for index, scenario_set_name in enumerate(scenario_sets):
            scenarios_in_set = [scenarios[i][index] for i in range(len(scenarios))]
            df_this_set = df[df[scenario_label].isin(scenarios_in_set)].copy()
            df_this_set[scenario_label] = scenario_set_name
            dataframes.append(df_this_set)
        df = pd.concat(dataframes)
        scenario_list = scenario_sets

    # Set the appropriate dictionary for the landtype group.
    if landtype_groups == 'standard':
        landtype_groups = gcam_landtype_groups
    elif landtype_groups == 'nonstandard':
        landtype_groups = gcam_landtype_groups_nonstandard

    # Create DataFrames for 1) all categories (create a copy of the entire DataFrame), 2) categories that correspond to
    # a group of landtypes (e.g., forest, crop, grass, shrub, pasture), and 3) a set of individual categories.
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

    # Calculate a percent difference with respect to the first scenario or first scenario set as the reference.
    if plot_percent_difference:
        reference = df[df[scenario_label] == scenario_list[0]][value_label].to_numpy()
        for scenario in scenario_list:
            df.loc[df[scenario_label] == scenario, value_label] -= reference
            df.loc[df[scenario_label] == scenario, value_label] *= 100/(reference + EPSILON)
        # Delete the reference from the DataFrame.
        df = df[df[scenario_label] != scenario_list[0]]
        
    # Depending on if the plot involves hues or not, set the number of colors in the palette to match the number of x-axis quantities in the box plot.
    if hue:
        num_x_variables = len(df[hue].unique())
    else:
        num_x_variables = len(df[x_variable].unique())

    if not plot_percent_difference:
        palette = plot_colors[:num_x_variables]
    else:
        # When plotting a percent difference, offset the index by 1 because index 0 is for the now-deleted first scenario or scenario set (reference).
        palette = plot_colors[1:num_x_variables+1]
    
    # If all boxes are of the same color (the 'hue' variable is None since the default color will be used for all boxes), we do not need a legend.
    if hue:
        sns.boxplot(df, x=x_variable, y=value_label, hue=hue, linewidth=linewidth, palette=palette, fliersize=marker_size, fill=fill_boxes)
    else:
        plot_options['legend_on'] = False
        sns.boxplot(df, x=x_variable, y=value_label, legend=None, linewidth=linewidth, fliersize=marker_size, fill=fill_boxes)

    # Finalize the box plot now that all x-axis variables have been processed.
    plot_options['name'] = plot_name
    set_figure_options(fig, ax, plot_options)

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
        print('Usage: python gcam_plot_box_and_whiskers.py `path/to/json/input/file(s)\'')
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

    # Create all of the bpx plots in parallel.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(plot_box_and_whiskers, list_of_inputs)
    
    # Print the total execution time to produce all the plots.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all plots: {elapsed_time:.2f} seconds")