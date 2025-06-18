import itertools
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
from utility_plots import *

""" Dictionary of default input values for time series plots. """
default_inputs_time_series = {
    'areas_in_thousands_km2': True,
    'category_label': 'sector',
    'end_year': 2100,
    'error_bars_alpha': 0.2,
    'height': height_default,
    'include_mean_across_all_data': True,
    'legend_bbox_x': None,
    'legend_label_size': legend_label_size_default,
    'legend_num_columns': None,
    'legend_on': legend_on_default,
    'legend_place_outside': True, 
    'linestyle_tuples': linestyle_tuples_default,
    'linewidth': linewidth_default,
    'marker_size': 6,
    'markers': markers_default,
    'multiplier': 1,
    'p_value_file': 'p_values.dat',
    'p_value_file_print_only_if_below_threshold': True,
    'p_value_marker_size': 10, 
    'p_value_threshold': 0.05,
    'plot_colors': plot_colors_default,
    'plot_directory': './',
    'plot_percent_difference': False,
    'plot_type': 'ensemble_averages',
    'produce_png': produce_png_default, 
    'region_label': 'region', 
    'regions': ['global'],
    'scenario_label': 'scenario', 
    'start_year': 2015,
    'std_mean_across_all_data_multiplier': 1,
    'std_multiplier': 1, 
    'use_latex': use_latex_default, 
    'value_label': 'value',
    'width': width_default,   
    'x_label_size': axis_label_size_default,
    'x_limits': axis_limits_default,             
    'x_scale': scale_default,   
    'x_tick_label_size': tick_label_size_default,   
    'y_label_size': axis_label_size_default,
    'y_limits': axis_limits_default,
    'y_scale': scale_default,
    'y_tick_label_size': tick_label_size_default,
    'year_label': 'year'            
}

def process_inputs(inputs):
    """ 
    Processes a dictionary of inputs (keys are options, values are choices for those options) for creating time series plots.

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

    if 'categories' not in inputs:
        if category_label in df.columns:
            inputs['categories'] = df[category_label].unique()
            if 'categories_to_exclude' in inputs:
                inputs['categories'] = [column for column in inputs['categories'] if column not in inputs['categories_to_exclude']]
        else:
            inputs['categories'] = 'all_categories'
    categories = inputs['categories']
    # If the user entered a string indicating a single category, put that string in a list.
    if isinstance(categories, str):
        categories = [categories]
        inputs['categories'] = categories

    if 'plot_directory' not in inputs:
        inputs['plot_directory'] = default_inputs_time_series['plot_directory']
    if not os.path.exists(inputs['plot_directory']):
        os.makedirs(inputs['plot_directory'])
    if 'p_value_file' not in inputs:
        inputs['p_value_file'] = os.path.join(inputs['plot_directory'], default_inputs_time_series['p_value_file'])

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
        inputs['plot_name'] = os.path.join(inputs['plot_directory'], 'time_series_' + output_file_name)

    # Add keys for plotting options that have not been specified in the inputs dictionary and use default values for them.
    for key in default_inputs_time_series.keys():
        if key not in inputs:
            inputs[key] = default_inputs_time_series[key]

    # If the user specified anything other than a dictionary for the 'regions', create a dictionary with the keys given by the categories. 
    if 'regions' not in inputs or not isinstance(inputs['regions'], dict):
        inputs['regions'] = dict.fromkeys(categories, inputs['regions'])
    
    return inputs

def plot_time_series(inputs):
    """ 
    Creates a time series plot with years on the x-axis and perform statistical analysis for a single output file. The data in the file are organized
    into scenarios or scenario sets, categories, and regions.
    Types of plots: 1) Direct plots in which each scenario is treated as an individual curve in the time series collection. 
    2) Ensemble plots in which the scenarios are grouped into sets and averages of each set are plotted. The scenarios must be specified as a 
    list of lists for ensemble plots.

    Parameters:
        input: Dictionary containing user inputs for different plotting options, where the keys are options and values are choices for those options.
               This dictionary is assumed to be complete (pre-processed).

    Returns:
        N/A.
    """
    # This function creates a time series plot where the data (scenarios or scenario sets, categories, regions) all come from a single output file.
    start_time = time.time()
    output_file = inputs['output_file']
    
    # Extract all other plotting options.
    categories = inputs['categories']
    category_label = inputs['category_label']
    end_year = inputs['end_year']
    error_bars_alpha = inputs['error_bars_alpha']
    height = inputs['height']
    include_mean_across_all_data = inputs['include_mean_across_all_data']
    legend_bbox_x = inputs['legend_bbox_x']
    legend_label_size = inputs['legend_label_size']
    legend_num_columns = inputs['legend_num_columns']
    legend_on = inputs['legend_on']
    legend_place_outside = inputs['legend_place_outside']
    linestyle_tuples = inputs['linestyle_tuples']
    linewidth = inputs['linewidth']
    marker_size = inputs['marker_size']
    markers = inputs['markers']
    multiplier = inputs['multiplier']
    p_value_file = inputs['p_value_file']
    p_value_file_print_only_if_below_threshold = inputs['p_value_file_print_only_if_below_threshold']
    p_value_marker_size = inputs['p_value_marker_size']
    p_value_threshold = inputs['p_value_threshold']
    plot_colors = inputs['plot_colors']
    plot_directory = inputs['plot_directory']
    plot_name = inputs['plot_name']
    plot_percent_difference = inputs['plot_percent_difference']
    plot_type = inputs['plot_type']
    produce_png = inputs['produce_png']
    region_label = inputs['region_label']
    scenario_label = inputs['scenario_label']
    scenarios = inputs['scenarios']
    start_year = inputs['start_year']
    std_mean_across_all_data_multiplier = inputs['std_mean_across_all_data_multiplier']
    std_multiplier = inputs['std_multiplier']
    use_latex = inputs['use_latex']
    value_label = inputs['value_label']
    width = inputs['width']
    x_label_size = inputs['x_label_size']
    x_limits = inputs['x_limits']
    x_scale = inputs['x_scale']
    x_tick_label_size = inputs['x_tick_label_size']
    y_label = inputs['y_label']
    y_label_size = inputs['y_label_size']
    y_limits = inputs['y_limits']
    y_scale = inputs['y_scale']
    y_tick_label_size = inputs['y_tick_label_size']
    year_label = inputs['year_label']

    # Initialize dictionary to store data from the first scenario or first scenario set, to be the reference for percent difference calculations. 
    reference_data = {}

    # Set the plotting options.
    if plot_percent_difference and (f'%' not in y_label or 'percent' not in y_label):
        y_label += rf' (\% difference)'
    plot_options = dict(width=width, height=height, name=plot_name, x_label='Year', y_label=fr'{y_label}')
    plot_options.update(zip(['x_scale', 'y_scale', 'x_limits', 'y_limits', 'use_latex'], [x_scale, y_scale, x_limits, y_limits, use_latex]))
    plot_options.update(zip(['x_tick_label_size', 'y_tick_label_size', 'legend_on'], [x_tick_label_size, y_tick_label_size, legend_on]))
    plot_options.update(zip(['x_label_size', 'y_label_size', 'legend_label_size'], [x_label_size, y_label_size, legend_label_size]))
    plot_options.update(zip(['produce_png', 'legend_bbox_x', 'legend_place_outside'], [produce_png, legend_bbox_x, legend_place_outside]))

    # Use LaTeX fonts for figures and set font size of tick labels.
    if use_latex:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', weight='bold')
    plt.rcParams['xtick.labelsize'] = x_tick_label_size
    plt.rcParams['ytick.labelsize'] = y_tick_label_size

    # Create an empty DataFrame that will later be used to store all of the time series data and to calculate their statistics.
    df_all_time_series = pd.DataFrame()
    have_not_stored_time_in_df = True

    df = read_file_into_dataframe(output_file)
    df = df[(df[year_label] >= start_year) & (df[year_label] <= end_year)]

    fig, ax = plt.subplots(nrows=1, ncols=1)

    # Direct plots, in which each such time series plot can include one or more sets of individual (not grouped) curves.
    if not check_is_list_of_lists(scenarios) or plot_type == 'direct':
        
        # If scenarios is a list of lists, but we want a direct plot, unpack scenarios into a list of strings, one string for each output file.
        if check_is_list_of_lists(scenarios) and plot_type == 'direct':
            scenarios = list(itertools.chain.from_iterable(scenarios))
            # Turn the corresponding plot_colors to list of lists with enough rows to match, then unpack them.
            plot_colors = list(itertools.chain.from_iterable(plot_colors))

        num_scenarios = len(scenarios)
        for scenario_index, scenario in enumerate(scenarios):
            scenario = scenarios[scenario_index]
            df_this_scenario = df[df[scenario_label] == scenario]
            num_categories = len(categories)
            for category_index, category in enumerate(categories):

                # If we do not want to consider all categories, focus on just the subset pertaining to the current category.
                if category != 'all_categories':
                     df_this_category = df_this_scenario[df_this_scenario[category_label].str.contains(category)]
                else:
                    df_this_category = df_this_scenario
                # The markers (symbols) on the plot will vary by category.
                marker = markers[category_index]
                
                if num_scenarios == 1 or (num_scenarios == 2 and plot_percent_difference):
                    # If we have 2 scenarios but plotting a percent difference, 
                    # the reference scenario is not included, so only one set of curves will appear on the plot.
                    line_color = plot_colors[category_index]
                    if not legend_num_columns:
                        plot_options.update(zip(['legend_num_columns'], [1]))
                else:
                    # Base the line colors on the scenario index if there will be more than one set of curves on the plot.
                    line_color = plot_colors[scenario_index]
                    if not legend_num_columns:
                        if plot_percent_difference:
                            plot_options.update(zip(['legend_num_columns'], [num_scenarios-1]))
                        else:
                            plot_options.update(zip(['legend_num_columns'], [num_scenarios]))

                num_regions = len(inputs['regions'][category])
                for region_index, region in enumerate(inputs['regions'][category]):
                    
                    if region != 'global':
                        df_this_region = df_this_category[df_this_category[region_label] == region]
                    else:
                        df_this_region = df_this_category
                    # The linestyles (e.g., solid, dotted, dashed) on the plot will vary by region.
                    linestyle = linestyle_tuples[region_index][1]
                    
                    # For each year, take the mean over all categories and regions that are subsets of current category and region.
                    x = df_this_region[year_label].unique()
                    y = df_this_region.groupby(year_label)[value_label].mean()*multiplier

                    if num_categories > 1:
                        label = category
                        if num_regions > 1:
                            label += f' ({region})'
                    elif num_regions > 1:
                        label = region
                    else:
                        label = scenario
                    
                    # If plotting a percent difference, store data for the first scenario and calculate the percent difference with respect to that.
                    if plot_percent_difference:
                        if scenario_index == 0:
                            reference_data[f'{category}_{region}'] = y
                        y = (y - reference_data[f'{category}_{region}'])/(reference_data[f'{category}_{region}'] + EPSILON)*100
                    # Do not include the first scenario if plotting a percent difference. Include the data in the plot otherwise.
                    if not plot_percent_difference or scenario_index != 0:
                        ax.plot(x, y, label=label, color=line_color, linestyle=linestyle, linewidth=linewidth, marker=marker, markersize=marker_size)
                    # Join the time series data from the current scenario into the larger DataFrame.
                    if have_not_stored_time_in_df:
                        x_series = pd.Series(x, name='Year')
                        df_all_time_series = pd.concat([df_all_time_series, x_series], axis=1)
                        have_not_stored_time_in_df = False
                    y_series = pd.Series(list(y), name=f'scen={scenario_index}_cat={category_index}_reg={region_index}')
                    df_all_time_series = pd.concat([df_all_time_series, y_series], axis=1)

        # Calculate the overall mean and standard deviation (for error bars) across all time series.
        df = df_all_time_series
        if include_mean_across_all_data:
            if plot_percent_difference:
                # If plotting a percent difference with respect to the first scenario, do not include that scenario in calculating the overall mean.
                columns = [column for column in df.columns if column != 'Year' and 'scen=0' not in column]
            else:
                columns = [column for column in df.columns if column != 'Year']
            x = df['Year']
            y = df[columns].mean(axis=1)
            ax.plot(x, y, label='Mean', color='k', linestyle=linestyle_tuples[0][1], linewidth=linewidth)
            if std_mean_across_all_data_multiplier:
                error = df[columns].std(axis=1)*std_mean_across_all_data_multiplier
                ax.fill_between(x, y-error, y+error, color='k', alpha=error_bars_alpha)

        if num_scenarios > 1:
            # Perform t-tests to compare the entire time series in the first scenario (assumed to be control) vs. other scenarios.
            for scenario_index in range(1, num_scenarios):
                scenario = scenarios[scenario_index]
                for category_index, category in enumerate(categories):
                    for region_index, region in enumerate(inputs['regions'][category]):
                        # Categories and regions between the current scenario and the control scenario are matched up before performing the t-test.
                        control_data = df[f'scen=0_cat={category_index}_reg={region_index}']
                        test_data = df[f'scen={scenario_index}_cat={category_index}_reg={region_index}']
                        ttest = stats.ttest_ind(control_data, test_data)
                        label = f'scenario={scenario}, category={category}, region={region}'
                        print_p_values(ttest, label, p_value_threshold, p_value_file, plot_directory, p_value_file_print_only_if_below_threshold)
        elif num_scenarios == 1 and num_regions > 1:
            # If there is only one scenario, but multiple regions for that scenario, perform inter-regional t-tests for each category.
            for category_index, category in enumerate(categories):
                for region_index in range(1, len(inputs['regions'][category])):
                    control_data = df[f'scen=0_cat={category_index}_reg=0']
                    test_data = df[f'scen=0_cat={category_index}_reg={region_index}']
                    ttest = stats.ttest_ind(control_data, test_data)
                    label = f'scenario={scenario}, category={category}, region={region}'
                    print_p_values(ttest, label, p_value_threshold, p_value_file, plot_directory, p_value_file_print_only_if_below_threshold)
    
    # Ensemble plots, in which the ensemble is further subdivided into groups of curves, where each group represents a set of scenarios.
    elif check_is_list_of_lists(scenarios) and plot_type == 'ensemble_averages':

        # Get the total number of scenario sets (groups) in the ensemble and the number of scenarios in each set.
        num_scenario_sets = len(scenarios[0])
        num_scenario_in_each_set = len(scenarios)
        # Unlike in the case of direct plots, the plot colors in ensemble plots are determined only by the scenario sets.
        if not legend_num_columns:
            if plot_percent_difference:
                plot_options.update(zip(['legend_num_columns'], [num_scenario_sets-1]))
            else:
                plot_options.update(zip(['legend_num_columns'], [num_scenario_sets]))
        
        for scenario_set_index in range(num_scenario_sets):
            scenarios_in_set = [scenarios[i][scenario_set_index] for i in range(num_scenario_in_each_set)]
            for scenario_index, scenario in enumerate(scenarios_in_set):
                scenario = scenarios_in_set[scenario_index]
                df_this_scenario = df[df[scenario_label] == scenario]
                for category_index, category in enumerate(categories):
                    if category != 'all_categories':
                        # Like in the direct plots, if we do not want to consider all categories, focus on just the subset of the current category.
                        df_this_category = df_this_scenario[df_this_scenario[category_label].str.contains(category)]
                    else:
                        df_this_category = df_this_scenario
                    for region_index, region in enumerate(inputs['regions'][category]):

                        if region != 'global':
                            df_this_region = df_this_category[df_this_category[region_label] == region]
                        else:
                            df_this_region = df_this_category

                        x = df_this_region[year_label].unique()
                        y = df_this_region.groupby(year_label)[value_label].mean()*multiplier

                        # Join the time series data for the current region and category in the current scenario set into the larger DataFrame.
                        if have_not_stored_time_in_df:
                            x_series = pd.Series(x, name='Year')
                            df_all_time_series = pd.concat([df_all_time_series, x_series], axis=1)
                            have_not_stored_time_in_df = False
                        y_series = pd.Series(list(y), name=f'set={scenario_set_index}_scen={scenario_index}_cat={category_index}_reg={region_index}')
                        df_all_time_series = pd.concat([df_all_time_series, y_series], axis=1)
        
        # Now that each individual time series has been stored, group them into scenario sets and perform analysis on the group/set means.
        x = x_series
        # Initialize DataFrame to store all of the set means.
        df_all_set_means = pd.DataFrame()
        have_not_stored_time_in_df = True
        num_categories = len(categories)
        for category_index, category in enumerate(categories):
            # Like in the case of direct plots, markers are determined by the category.
            marker = markers[category_index]
            columns = [column for column in df_all_time_series.columns if f'cat={category_index}' in column]
            df_this_category = df_all_time_series[columns]
            num_regions = len(inputs['regions'][category])
            for region_index, region in enumerate(inputs['regions'][category]):
                # Like in the case of direct plots, markers are determined by the region.
                linestyle = linestyle_tuples[region_index][1]
                columns = [column for column in  df_this_category.columns if f'reg={region_index}' in column]
                df_this_region = df_this_category[columns]
                for scenario_set_index in range(num_scenario_sets):

                    line_color = plot_colors[scenario_set_index]
                    columns = [column for column in df_this_region.columns if f'set={scenario_set_index}' in column]
                    df_this_set = df_this_region[columns]

                    if num_categories > 1:
                        label = category
                        if num_regions > 1:
                            label += f' ({region})'
                    elif num_regions > 1:
                        label = region
                    else:
                        label = scenario_set_index

                    # Calculate the mean and standard deviation over all scenarios in the current scenario set.
                    y = df_this_set[columns].mean(axis=1)
                    y_std = df_this_set[columns].std(axis=1)

                    # If plotting a percent difference, store data for the first set and calculate the percent difference with respect to that.
                    if plot_percent_difference:
                        if scenario_set_index == 0:
                            reference_data[f'{category}_{region}'] = y
                        y = (y - reference_data[f'{category}_{region}'])/(reference_data[f'{category}_{region}'] + EPSILON)*100
                        y_std = y_std/(reference_data[f'{category}_{region}'] + EPSILON)*100

                    # Join the time series data from the current set into the larger DataFrame that stores the set means.
                    if have_not_stored_time_in_df:
                        df_all_set_means = pd.concat([df_all_set_means, x_series], axis=1)
                        have_not_stored_time_in_df = False
                    y_series = pd.Series(list(y), name=f'set={scenario_set_index}_cat={category_index}_reg={region_index}')
                    df_all_set_means = pd.concat([df_all_set_means, y_series], axis=1)

                    # Plot the annual time series, including possibly the error bars. 
                    # Do not include the first file set if plotting a percent difference. Include the data in the plot otherwise.
                    if not plot_percent_difference or scenario_set_index != 0:
                        ax.plot(x, y, label=label, color=line_color, linestyle=linestyle, linewidth=linewidth, marker=marker, markersize=marker_size)
                        if std_multiplier:
                            error = y_std*std_multiplier
                            ax.fill_between(x, y-error, y+error, color=line_color, alpha=error_bars_alpha)

                    if num_scenario_sets > 1 and scenario_set_index > 0:
                        # Perform t-test to compare the entire time series in the first scenario set (assumed to be control) vs. other sets.
                        control_data = df_all_set_means[f'set=0_cat={category_index}_reg={region_index}']
                        test_data = df_all_set_means[f'set={scenario_set_index}_cat={category_index}_reg={region_index}']
                        ttest = stats.ttest_ind(control_data, test_data)
                        label = f'set={scenario_set_index}, category={category}, region={region}'
                        print_p_values(ttest, label, p_value_threshold, p_value_file, plot_directory, p_value_file_print_only_if_below_threshold)
                    elif num_scenario_sets == 1 and num_regions > 1:
                        # If there is only one scenario set, but multiple regions for that set, perform inter-regional t-tests for each category.
                        control_data = df_all_set_means[f'set=0_cat={category_index}_reg=0']
                        test_data = df_all_set_means[f'set=0_cat={category_index}_reg={region_index}']
                        ttest = stats.ttest_ind(control_data, test_data)
                        label = f'category={category}, region={region}'
                        print_p_values(ttest, label, p_value_threshold, p_value_file, plot_directory, p_value_file_print_only_if_below_threshold)

                    # Perform a t-test to compare the control against the current data set at each time period (each year).
                    if num_scenario_sets > 1 and scenario_set_index > 0:
                        columns_control_set = [column for column in df_all_time_series.columns if 
                                        (f'set=0' in column and f'cat={category_index}_reg={region_index}' in column)]
                        columns_this_set = [column for column in df_all_time_series.columns if 
                                        (f'set={scenario_set_index}' in column and f'cat={category_index}_reg={region_index}' in column)]
                        p_values = df_all_time_series.apply(perform_ttest, columns_set_1=columns_control_set, \
                                                               columns_set_2=columns_this_set, axis=1).fillna(1)
                        mask = p_values <= p_value_threshold
                        ax.plot(df_all_time_series['Year'][mask], y_series[mask], color=line_color, linestyle='None', linewidth=linewidth, \
                                  marker=marker, markersize=p_value_marker_size)

                # Calculate the overall mean across all scenario sets and a standard deviation for group/set means.
                if include_mean_across_all_data:
                    if plot_percent_difference:
                    # If plotting a percent difference, do not include the first scenario set in calculating the overall mean.
                        columns = [column for column in df_all_set_means.columns if column != 'Year' \
                                   and ('set=0' not in column and f'cat={category_index}_reg={region_index}' in column)]
                    else:
                        columns = [column for column in df_all_set_means.columns if column != 'Year' \
                                   and f'cat={category_index}_reg={region_index}' in column]
                    x = df_all_set_means['Year']
                    y = df_all_set_means[columns].mean(axis=1)
                    label = 'Mean'
                    if num_categories > 1 and num_regions > 1:
                        label += f' ({category}_{region})'
                    elif num_categories > 1:
                        label += f' ({category})'
                    elif num_categories > 1:
                        label += f' ({region})'
                    ax.plot(x, y, label=label, color='k', linestyle=linestyle, linewidth=linewidth, marker=marker, markersize=marker_size)
                    if std_mean_across_all_data_multiplier:
                        error = df_all_set_means[columns].std(axis=1)*std_mean_across_all_data_multiplier
                        ax.fill_between(x, y-error, y+error, color='k', alpha=error_bars_alpha)

    # Finalize the time series plot now that all curves for the output file have been processed.
    plot_options['name'] = plot_name
    set_figure_options(fig, ax, plot_options)

    # Close the plot now that we are done with it. Record the elapsed time.
    plt.close(fig)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing plots in {plot_directory}: {elapsed_time:.2f} seconds")


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file(s) on the command line.
    start_time = time.time()
    if len(sys.argv) < 2:
        print('Usage: python plot_time_series.py `path/to/json/input/file(s)\'')
        sys.exit()

    # Read and load the JSON file(s) into a list of dictionaries.
    inputs = []
    for i in range(1, len(sys.argv)):
        input_file = sys.argv[i]
        with open(input_file) as f:
            inputs.extend(json.load(f))

    # Process each dictionary to produce a list of smaller dictionaries, where each smaller dictionary specifies options for a single plot.
    start_time = time.time()
    list_of_inputs = []
    for index in range(len(inputs)):
        # Process the inputs to fill in missing plotting input choices with default values, etc., and add to the list of dictionaries.
        list_of_inputs.append(process_inputs(inputs[index]))

    # Delete all the p-value files before we do any calculations to start a fresh run.
    for inputs in list_of_inputs:
        file = inputs['p_value_file']
        if os.path.exists(file): 
            os.remove(file)

    # Create all of the times series plots in parallel.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(plot_time_series, list_of_inputs)

    # Sort all the p-value files alphabetically.
    for inputs in list_of_inputs:
        file = inputs['p_value_file']
        if os.path.exists(file): 
            sort_file(file)
    
    # Print the total execution time to produce all the plots.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all plots: {elapsed_time:.2f} seconds")