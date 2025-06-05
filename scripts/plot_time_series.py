import itertools
import json
from matplotlib import pyplot as plt
import multiprocessing
import numpy as np
import os
import pandas as pd
from scipy import stats
import sys
import time
from utility_constants import *
from utility_dataframes import get_columns_without_units_in_dataframe, get_matching_column_in_dataframe, perform_ttest, read_file_into_dataframe
from utility_functions import check_is_list_of_lists, check_substrings_in_list, print_p_values, replace_inside_parentheses, sort_file
from utility_plots import *

""" Dictionary of default input values for time series plots. """
default_inputs_time_series = {'plot_directory': './',
                    'calculation_type': 'mean',
                    'plot_type': 'ensemble_averages',
                    'plot_percent_difference': False,
                    'multiplier': 1,
                    'std_annual_multiplier': 1, 
                    'std_seasons_multiplier': None,
                    'std_monthly_multiplier': 1,
                    'error_bars_alpha': 0.2,
                    'areas_in_thousands_km2': True,
                    'start_year': 2015,
                    'end_year': 2100,
                    'width': width_default,
                    'height': height_default,
                    'x_scale': scale_default,
                    'y_scale': scale_default,
                    'x_limits': axis_limits_default,
                    'y_limits': axis_limits_default,
                    'x_tick_label_size': tick_label_size_default,
                    'y_tick_label_size': tick_label_size_default,
                    'x_label_size': axis_label_size_default,
                    'y_label_size': axis_label_size_default,
                    'legend_label_size': legend_label_size_default,
                    'legend_on': legend_on_default,
                    'linewidth': linewidth_default,
                    'plot_colors': plot_colors_default,
                    'linestyle_tuples': linestyle_tuples_default,
                    'use_latex': use_latex_default, 
                    'produce_png': produce_png_default, 
                    'include_annual_mean_across_all_sets': True,
                    'include_monthly_mean_across_all_sets': True,
                    'std_annual_mean_across_all_sets_multiplier': 1,
                    'std_monthly_mean_across_all_sets_multiplier': 1,
                    'p_value_file': None,
                    'p_value_threshold': 0.05,
                    'p_value_file_print_only_if_below_threshold': True,
                    'p_value_marker': 'o',
                    'p_value_marker_size': 6, 
                    'include_seasons': {'MAM': False, 'JJA': False, 'SON': False, 'DJF': False},
                    'seasons_to_plot_separately': {'MAM': False, 'JJA': False, 'SON': False, 'DJF': False},
                    'monthly_time_series_plot': False,
                    'monthly_time_series_start_year': 2071,
                    'monthly_time_series_end_year': 2090,
                    'monthly_time_series_x_limits': axis_limits_default,
                    'monthly_time_series_y_limits': axis_limits_default
}

def process_inputs(inputs):
    """ 
    Processes a dictionary of inputs (keys are options, values are choices for those options) for creating time series plots.

    Parameters:
        inputs: Dictionary containing the user plotting choice inputs for different options. This dictionary may be incomplete or have invalid values.

    Returns:
        List of dictionaries, each of which specifies all plotting options for a single variable. If the user did not select a plotting option for
        a particular variable, the default choice for that plotting option will be selected.
    """
    # The user can specify output_files to be a string (indicating a single output file), a list (or a list of lists) of different output files, 
    # or a dictionary of output files. Read the output file (any one will do if there is more than one such file) and put into a Pandas DataFrame.
    if isinstance(inputs['output_files'], str):
        file = inputs['output_files']
    elif check_is_list_of_lists(inputs['output_files']):
        file = inputs['output_files'][0][0]
    elif isinstance(inputs['output_files'], list):
        file = inputs['output_files'][0]
    elif isinstance(inputs['output_files'], dict):
        file = inputs['output_files'][list(inputs['output_files'].keys())[0]]
    df = read_file_into_dataframe(file)

    # If the user entered the string 'all' for the variables or no input at all for the variables, assume that they want to make plots for
    # all columns (all variables) that are in the DataFrame, except for the Year and Month columns.
    if 'variables' not in inputs:
        inputs['variables'] = None
    variables = inputs['variables']
    if not variables or variables == 'all':
        variables = get_columns_without_units_in_dataframe(df)
        variables.remove('Year')
        variables.remove('Month')
        inputs['variables'] = variables
    # If the user entered a string indicating a single variable, put that string in a list.
    if isinstance(variables, str):
        variables = [variables]
        inputs['variables'] = variables

    # Add keys for plotting options that have not been specified in the inputs dictionary, and set to an empty dictionary or use default values.
    for key in default_inputs_time_series.keys():
        if key not in inputs:
            if key in ['include_seasons', 'seasons_to_plot_separately']:
                # For these two plotting options, set them as empty dictionaries if they are not keys in the inputs dictionary.
                inputs[key] = {}
            else:
                # Use the default values for the other plotting options.
                inputs[key] = default_inputs_time_series[key]

    # If the user specified anything other than a dictionary (e.g., a single value [string, integer, float] or a list) for the other plotting options, 
    # assume that they want to use that value/list for all the variables. Enable this by creating dictionaries with the keys given by the variables. 
    # This also covers the case above, where default values were added for all plotting options that are not specified in the inputs dictionary.
    for key, value in inputs.items():
        if key not in ['variables', 'include_seasons', 'seasons_to_plot_separately']:
            if not isinstance(value, dict):
                # If we have a single output file and corresponding label, and these are specified as strings, turn them into a list with this string.
                if (key == 'output_files' or key == 'output_labels') and isinstance(value, str):
                    value = [value]
                inputs[key] = dict.fromkeys(variables, value)
        elif key in ['include_seasons', 'seasons_to_plot_separately']:
            # If the user entered a dictionary like {'MAM': True, 'JJA': True, 'SON': True, 'DJF': True}, use that for all of the variables.
            if check_substrings_in_list(inputs[key].keys(), ['MAM', 'JJA', 'SON', 'DJF'], all_or_any='all'):
                inputs[key] = dict.fromkeys(variables, value)

    # For each variable, if a plotting option is missing (has not been specified), fill in with the default corresponding to that plotting option.
    if 'plot_name' not in inputs:
        inputs['plot_name'] = {}
    if 'y_label' not in inputs:
        inputs['y_label'] = {}
    for variable in variables:
        # Use the default for the plot directory and create the directory if it does not already exist.
        if not any(key == variable for key in inputs['plot_directory'].keys()):
            inputs['plot_directory'][variable] = default_inputs_time_series['plot_directory']
        if not os.path.exists(inputs['plot_directory'][variable]):
            os.makedirs(inputs['plot_directory'][variable])
        # Default for the plot names is to call it 'time_series_[var_name]', where '[var_name]' is the name of the variable.
        if not any(key == variable for key in inputs['plot_name'].keys()):
            inputs['plot_name'][variable] = os.path.join(inputs['plot_directory'][variable], 'time_series_' + variable)
        # Default for the y_label of a variable is to use the column header for that variable from the DataFrame.
        if not any(key == variable for key in inputs['y_label'].keys()):
            inputs['y_label'][variable] = get_matching_column_in_dataframe(df, variable)
        # Default for the other plotting options are specified in the default_inputs_time_series dictionary.
        for key, value in inputs.items():
            if key not in ['variables', 'plot_directory', 'plot_name', 'y_label']:
                if not any(value_key == variable for value_key in value.keys()):
                    inputs[key][variable] = default_inputs_time_series[key]

    # Now that the dictionary has been populated with complete plotting options for each variable, separate it into a list of dictionaries,
    # where each of these smaller dictionaries contain the plotting options for a single variable. Return this list of dictionaries.
    list_of_inputs = []
    for variable in inputs['variables']:
        inputs_for_this_variable = {'variable': variable}
        for input_option in inputs.keys():
            if input_option != 'variables':
                inputs_for_this_variable[input_option] = inputs[input_option][variable]
        list_of_inputs.append(inputs_for_this_variable)
    return list_of_inputs

def plot_seasons(include_seasons, ax, df, x, output_label, plot_colors, linestyle_tuples, linewidth, columns, 
                 reference_data=None, file_or_file_set_index=None, std_multiplier=None, error_bars_alpha=None):
    """ 
    Adds seasons (as specified in the include_seasons dictionary) to a given time series plot. 

    Parameters:
        include_seasons: Dictionary specifying the seasons to include in the plot.
        ax: Axis object for the plot.
        df: Pandas DataFrame containing the time series data.
        x: Years to put on the x-axis.
        output_label: Base output label.
        plot_colors: List of plot line colors.
        linestyle_tuples: Tuples for the linestyles.
        linewidth: Thickness of the lines.
        columns: String (if a single column) or list (multiple columns) for the datasets in the DataFrame that we want to plot. 
                 If a list, plot the mean across all columns.
        reference_data: Dictionary to store the data from the first file or first file set, which will serve as the reference set for 
                 percent difference calculations.
        file_or_file_set_index: Index for the particular file or file set, so that we can flag whether the current data pertain to the reference set.
        std_multiplier: Multiplier on the standard deviation for plotting error bars.
        error_bars_alpha: Opacity value for the error bars.

    Returns:
        The reference_data, which will be updated with entries for the reference data for each season if this dictionary was not initially empty.
    """
    seasons = ['MAM', 'JJA', 'SON', 'DJF']
    months_geq = [3, 6, 9, 12]
    months_leq = [5, 8, 11, 2]
    linestyle_tuples = [linestyle_tuples[2][1], linestyle_tuples[2][1], linestyle_tuples[1][1], linestyle_tuples[1][1]]
    linewidths = [linewidth+2, linewidth, linewidth+2, linewidth]
    # Iterate over all seasons and include only the ones that are specified in the include_seasons dictionary.
    for season_index in range(len(seasons)):
        season = seasons[season_index]
        if include_seasons[season]:
            month_geq = months_geq[season_index]
            month_leq = months_leq[season_index]
            linestyle=linestyle_tuples[season_index]
            linewidth=linewidths[season_index]
            # December--January--February needs to be handled differently than the other seasons.
            if season == 'DJF':
                condition = (df['Month'] >= month_geq) | (df['Month'] <= month_leq)
            else:
                condition = (df['Month'] >= month_geq) & (df['Month'] <= month_leq)
            # If we have multiple datasets to plot, calculate the standard deviation across the datasets for each year for each of the seasons.
            if isinstance(columns, list):
                y = df[condition].groupby('Year', as_index=False).mean()[columns].mean(axis=1)
                y_std = df[condition].groupby('Year', as_index=False).mean()[columns].std(axis=1)
            else:
                y = df[condition].groupby('Year', as_index=False).mean()[columns]
            # If the reference data is not None or not empty, we assume that the user wants to perform a percent difference calculation.
            if reference_data:
                if file_or_file_set_index == 0:
                    reference_data[season] = y
                y = (y - reference_data[season])/(reference_data[season] + EPSILON)*100
                if isinstance(columns, list):
                    y_std = y_std/(reference_data[season] + EPSILON)*100
            # Add the season to the plot, including its error bar if the multiplier is not None.
            # Do not include the reference data if plotting a percent difference. Include the data in the plot otherwise.
            if not reference_data or file_or_file_set_index != 0:
                ax.plot(x, y, label=output_label + f' ({season})', color=plot_colors, linestyle=linestyle, linewidth=linewidth)
                if std_multiplier and isinstance(columns, list):
                    error = y_std*std_multiplier
                    ax.fill_between(x, y-error, y+error, color=plot_colors, alpha=error_bars_alpha)
    return reference_data

def read_file_into_single_variable_dataframe(file, variable, start_year, end_year, multiplier, return_column=False):
    """ 
    Reads a file and puts the time series data for one specific variable from that file into a Pandas DataFrame. 

    Parameters:
        file: Directory and name for the file we want to read.
        variable: The variable of interest.
        start_year: First year in the extracted time series data.
        end_year: Last year in the extracted time series data.
        multiplier: Multiplier on the time series (in case the user wants to change the units).
        return_column: Boolean that specifies if we want to additionally return the column in the DataFrame corresponding to the variable.

    Returns:
        DataFrame containing the data for the variable and also the column in the DataFrame for that variable if return_column is True.
    """
    df = read_file_into_dataframe(file)
    df = df[(df['Year'] >= start_year) & (df['Year'] <= end_year)]
    # Find the column label for the variable, reduce the DataFrame to that column (plus columns for Year and Month).
    column = get_matching_column_in_dataframe(df, variable)
    df = df[['Year', 'Month', column]]
    df[column] *= multiplier
    if return_column:
        return df, column
    else:
        return df
    
def read_file_set_into_single_variable_dataframe(output_files, file_set_index, variable, start_year, end_year, multiplier):
    """ 
    Reads a list of output files and the time series data for one specific variable in those files into a Pandas DataFrame. 

    Parameters:
        output_files: List of lists (2D matrix) for the files. We want to read in one column of this matrix, which corresponds to a particular group
                      of files that we want to collect together into a set and later plot the various set means of a larger ensemble of sets.
        file_set_index: The column number for the output_files matrix. This indicates the index for the group/set of files that are of interest.
        variable: The variable of interest.
        start_year: First year in the extracted time series data.
        end_year: Last year in the extracted time series data.
        multiplier: Multiplier on the time series (in case the user wants to change the units).

    Returns:
        DataFrame containing data from all files for the variable and also all DataFrame columns for the variable, one column from each of the files.
    """
    file = output_files[0][file_set_index]
    num_files_in_set = len(output_files)
    df, column = read_file_into_single_variable_dataframe(file, variable, start_year, end_year, multiplier, return_column=True)
    # Rename the columns before doing the join operation to avoid naming conflicts for the various columns that correspond to the variable.
    df = df.rename(columns={column: column+f'_0_{file_set_index}'})
    for file_index in range(1, num_files_in_set):
        file = output_files[file_index][file_set_index]
        df_for_this_file = read_file_into_single_variable_dataframe(file, variable, start_year, end_year, multiplier)
        df_for_this_file = df_for_this_file.rename(columns={column: column+f'_{file_index}_{file_set_index}'})
        df = pd.merge(df, df_for_this_file, on=['Year', 'Month'], how='inner')
    # After the join operation, there will be multiple columns for the variable. Return the list of all of these columns along with the DataFrame.
    columns = get_matching_column_in_dataframe(df, variable, all_matches=True)
    return df, columns

def plot_time_series(inputs):
    """ 
    Parses a dictionary of inputs (keys are options, values are choices for those options) to create time series plots for a single variable. 
    Time series plots can be annual, with time presented in years (including seasonal variants that can be plotted separately), or 
    monthly, in which the values that are plotted indicate averages between specified start and end years for each month. 
    Types of plots: 1) Direct plots in which each output file is treated as an individual curve in the time series collection. 2) Ensemble plots in 
    which the files are grouped into sets and averages of each set are plotted. Output files must be specified as a list of lists for ensemble plots.

    Parameters:
        input: Dictionary containing the user plotting choice inputs for different options. This dictionary is assumed to be complete (pre-processed).

    Returns:
        N/A.
    """
    # This function creates time series plots for a single variable, and so it assumes that there is only one variable in the inputs dictionary.
    start_time = time.time()
    variable = inputs['variable']

    # Extract all plotting options for both annual and monthly time series plots.
    output_files = inputs['output_files']
    output_labels = inputs['output_labels']
    plot_directory = inputs['plot_directory']
    plot_name = inputs['plot_name']
    y_label = inputs['y_label']
    calculation_type = inputs['calculation_type']
    plot_type = inputs['plot_type']
    plot_percent_difference = inputs['plot_percent_difference']
    multiplier = inputs['multiplier']
    std_annual_multiplier = inputs['std_annual_multiplier']
    std_seasons_multiplier = inputs['std_seasons_multiplier']
    std_monthly_multiplier = inputs['std_monthly_multiplier']
    error_bars_alpha = inputs['error_bars_alpha']
    areas_in_thousands_km2 = inputs['areas_in_thousands_km2']
    start_year = inputs['start_year']
    end_year = inputs['end_year']
    width = inputs['width']
    height = inputs['height']
    x_scale = inputs['x_scale']
    y_scale = inputs['y_scale']
    x_limits = inputs['x_limits']
    y_limits = inputs['y_limits']
    x_tick_label_size = inputs['x_tick_label_size']
    y_tick_label_size = inputs['y_tick_label_size']
    x_label_size = inputs['x_label_size']
    y_label_size = inputs['y_label_size']
    legend_on = inputs['legend_on']
    legend_label_size = inputs['legend_label_size']
    linewidth = inputs['linewidth']
    plot_colors = inputs['plot_colors']
    linestyle_tuples = inputs['linestyle_tuples']
    use_latex = inputs['use_latex']
    produce_png = inputs['produce_png']
    include_seasons = inputs['include_seasons']
    seasons_to_plot_separately = inputs['seasons_to_plot_separately']
    include_annual_mean_across_all_sets = inputs['include_annual_mean_across_all_sets']
    std_annual_mean_across_all_sets_multiplier = inputs['std_annual_mean_across_all_sets_multiplier']
    include_monthly_mean_across_all_sets = inputs['include_monthly_mean_across_all_sets']
    std_monthly_mean_across_all_sets_multiplier = inputs['std_monthly_mean_across_all_sets_multiplier']
    p_value_threshold = inputs['p_value_threshold']
    p_value_file = inputs['p_value_file']
    p_value_file_print_only_if_below_threshold = inputs['p_value_file_print_only_if_below_threshold']
    p_value_marker = inputs['p_value_marker']
    p_value_marker_size = inputs['p_value_marker_size']

    # For area variables that are in units of km^2, plot them in units of thousands of km^2 if specified to do so.
    if areas_in_thousands_km2 and 'AREA' in variable and 'km^2' in y_label:
        multiplier = 1/1000
        y_label = y_label.replace('km^2', rf'thousands km$^2$')

    # Fix the labels involving m^2 or km^2 so that they get rendered properly.
    if 'm^2' in y_label:
        y_label = y_label.replace('m^2', rf'm$^2$')

    # Initialize dictionary to store the data from the first file or first file set, which will serve as the reference for percent difference 
    # calculations. Also update the units in the y-axis label if plotting a percent difference.
    reference_data = {}
    if plot_percent_difference:
        y_label = replace_inside_parentheses(y_label, rf'($\%$ difference)')

    # Set the plotting options
    plot_options = dict(width=width, height=height, name=plot_name, x_label='Year', y_label=fr'{y_label}')
    plot_options.update(zip(['x_scale', 'y_scale', 'x_limits', 'y_limits', 'use_latex'], [x_scale, y_scale, x_limits, y_limits, use_latex]))
    plot_options.update(zip(['x_tick_label_size', 'y_tick_label_size', 'legend_on'], [x_tick_label_size, y_tick_label_size, legend_on]))
    plot_options.update(zip(['x_label_size', 'y_label_size', 'legend_label_size'], [x_label_size, y_label_size, legend_label_size]))
    plot_options.update(zip(['produce_png'], [produce_png]))

    # Get plotting options relevant only to the monthly time series.
    monthly_time_series_plot = inputs['monthly_time_series_plot']
    monthly_time_series_start_year = inputs['monthly_time_series_start_year']
    monthly_time_series_end_year = inputs['monthly_time_series_end_year']
    monthly_time_series_x_limits = inputs['monthly_time_series_x_limits']
    monthly_time_series_y_limits = inputs['monthly_time_series_y_limits']

    # Use LaTeX fonts for figures and set font size of tick labels.
    if use_latex:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', weight='bold')
    plt.rcParams['xtick.labelsize'] = x_tick_label_size
    plt.rcParams['ytick.labelsize'] = y_tick_label_size

    # Figure and axis objects for annual time series (+ possibly including seasons), seasons separately, and monthly time series plots.
    fig, ax = plt.subplots(nrows=1, ncols=1)
    fig_seasons, ax_seasons = plt.subplots(nrows=1, ncols=1)
    fig_monthly, ax_monthly = plt.subplots(nrows=1, ncols=1)

    # Create an empty DataFrame that will later be used to store all of the annual and monthly time series data and to calculate their statistics.
    df_all_annual_time_series = pd.DataFrame()
    df_all_monthly_time_series = pd.DataFrame()

    # Direct plots, in which each such time series plot can include one or more sets of individual (not grouped) curves.
    if not check_is_list_of_lists(output_files) or plot_type == 'direct':
        
        # If output_files is a list of lists, but we want a direct plot, unpack output_files into a list of strings, one string for each output file.
        if check_is_list_of_lists(output_files) and plot_type == 'direct':
            # Turn the corresponding output_labels and plot_colors to list of lists with enough rows to match output_labels, then unpack them.
            output_labels = [output_labels]*len(output_files)
            output_labels = list(itertools.chain.from_iterable(output_labels))
            plot_colors = [plot_colors]*len(output_files)
            plot_colors = list(itertools.chain.from_iterable(plot_colors))
            output_files = list(itertools.chain.from_iterable(output_files))

        # Get the total number of files and iterate over each file.
        num_files = len(output_files)
        for file_index in range(num_files):

            # Get the label and line color for this output file.
            file = output_files[file_index]
            output_label = output_labels[file_index]
            line_color = plot_colors[file_index]

            # Time series data can be either averaged or summed over each month of the year.
            df, column = read_file_into_single_variable_dataframe(file, variable, start_year, end_year, multiplier, return_column=True)
            if calculation_type == 'mean':
                y = df.groupby('Year', as_index=False).mean()[column]
            elif calculation_type == 'sum':
                # Update the labels accordingly if a sum.
                y = df.groupby('Year', as_index=False).sum()[column]
                plot_options['y_label'] = y_label.replace('/month', '/year')
            
            # If plotting a percent difference with respect to the first file, store the data for that file and calculate the percent difference.
            if plot_percent_difference:
                if file_index == 0:
                    reference_data['annual'] = y
                y = (y - reference_data['annual'])/(reference_data['annual'] + EPSILON)*100
                    
            # Plot the annual time series, including possibly the error bars. 
            x = df.groupby('Year', as_index=False).mean()['Year']
            # Do not include the first file if plotting a percent difference. Include the data in the plot otherwise.
            if not plot_percent_difference or file_index != 0:
                ax.plot(x, y, label=output_label, color=line_color, linestyle=linestyle_tuples[0][1], linewidth=linewidth)
            # Join the annual time series data from the current file into the larger DataFrame.
            if file_index == 0:
                x_series = pd.Series(x, name='Year')
                df_all_annual_time_series = pd.concat([df_all_annual_time_series, x_series], axis=1)
            y_series = pd.Series(y, name=f'{column}_{file_index}')
            df_all_annual_time_series = pd.concat([df_all_annual_time_series, y_series], axis=1)

            # Add time series for one or more seasonal averages if specified to do so.
            if calculation_type == 'mean' and any(include_seasons.values()):
                reference_data = plot_seasons(include_seasons, ax, df, x, output_label, line_color, linestyle_tuples, linewidth, columns=column, 
                                              reference_data=reference_data, file_or_file_set_index=file_index)

            # Create a separate time series plot for the seasonal averages if specified to do so. Set name of this seasons-only plot accordingly.
            if calculation_type == 'mean' and any(seasons_to_plot_separately.values()):
                reference_data = plot_seasons(seasons_to_plot_separately, ax_seasons, df, x, output_label, line_color, linestyle_tuples, linewidth,
                                               columns=column, reference_data=reference_data, file_or_file_set_index=file_index)
                plot_options['name'] = plot_name + '_seasons'   
                set_figure_options(fig_seasons, ax_seasons, plot_options)

            # Create the monthly time series plot for the variable if specified to do so.
            if monthly_time_series_plot:
                df = read_file_into_single_variable_dataframe(file, variable, monthly_time_series_start_year, 
                                                                monthly_time_series_end_year, multiplier, return_column=False)
                y = df.groupby('Month', as_index=False).mean()[column]
                x = (df.groupby('Month', as_index=False).mean()['Month']).to_numpy()
                if plot_percent_difference:
                    if file_index == 0:
                        reference_data['monthly'] = y
                    y = (y - reference_data['monthly'])/(reference_data['monthly'] + EPSILON)*100
                #x = np.vectorize(MONTH_NUM_TO_NAME.get)(x)
                if not plot_percent_difference or file_index != 0:
                    ax_monthly.plot(x, y, label=output_label, color=line_color, linestyle=linestyle_tuples[0][1], linewidth=linewidth)
                #ax_monthly.tick_params(axis='x', labelrotation=45)
                # Join the monthly time series data from the current file into the larger DataFrame.
                if file_index == 0:
                    x_series = pd.Series(x, name='Month')
                    df_all_monthly_time_series = pd.concat([df_all_monthly_time_series, x_series], axis=1)
                y_series = pd.Series(y, name=f'{column}_{file_index}')
                df_all_monthly_time_series = pd.concat([df_all_monthly_time_series, y_series], axis=1)

        # Calculate the mean and standard deviation (for error bars) across all annual and monthly time series. Perform t-tests.
        conditions = [include_annual_mean_across_all_sets, monthly_time_series_plot and include_monthly_mean_across_all_sets]
        dataframes = [df_all_annual_time_series, df_all_monthly_time_series]
        time_columns = ['Year', 'Month']
        axes = [ax, ax_monthly]
        std_multipliers = [std_annual_mean_across_all_sets_multiplier, std_monthly_mean_across_all_sets_multiplier]
        for time_index, condition in enumerate(conditions):
            df = dataframes[time_index]
            time_column = time_columns[time_index]
            axis = axes[time_index]
            std_multiplier = std_multipliers[time_index]
            if condition:
                if plot_percent_difference:
                    # If plotting a percent difference with respect to the first file, do not include that file in calculating the overall mean.
                    columns = [column for column in df.columns if column != time_column and not column.endswith('_0')]
                else:
                    columns = [column for column in df.columns if column != time_column]
                x = df[time_column]
                y = df[columns].mean(axis=1)
                axis.plot(x, y, label='Mean', color='k', linestyle=linestyle_tuples[0][1], linewidth=linewidth)
                if std_multiplier:
                    error = df[columns].std(axis=1)*std_multiplier
                    axis.fill_between(x, y-error, y+error, color='k', alpha=error_bars_alpha)
            # Perform t-test to compare the entire annual or monthly time series in the first file (assumed to be control) vs. other files.
            if not df.empty:
                control_data = df[column + '_0']
                for file_index in range(1, num_files):
                    output_file = output_files[file_index]
                    test_data = df[column + f'_{file_index}']
                    ttest = stats.ttest_ind(control_data, test_data)
                    print_p_values(ttest, variable, p_value_threshold, p_value_file, output_file, p_value_file_print_only_if_below_threshold)
    
    # Ensemble plots, in which the ensemble is further subdivided into sets/groups of curves.
    elif check_is_list_of_lists(output_files) and plot_type == 'ensemble_averages':

        # Get the total number of file sets (groups) in the ensemble and iterate over each set.
        num_file_sets = len(output_files[0])
        for file_set_index in range(num_file_sets):

            # Get the label and line color for this file set.
            output_label = output_labels[file_set_index]
            line_color = plot_colors[file_set_index]

            # Calculate the mean over all files in this set. Time series data can be either averaged or summed over each month of the year.
            df, columns = read_file_set_into_single_variable_dataframe(output_files, file_set_index, variable, start_year, end_year, multiplier)
            if calculation_type == 'mean':
                y = df.groupby('Year', as_index=False).mean()[columns].mean(axis=1)
                y_std = df.groupby('Year', as_index=False).mean()[columns].std(axis=1)
                # Join the annual time series data from the current set into the larger DataFrame.
                df_all_annual_time_series = pd.concat([df_all_annual_time_series, df.groupby('Year', as_index=False).mean()], axis=1)
            elif calculation_type == 'sum':
                y = df.groupby('Year', as_index=False).sum()[columns].mean(axis=1)
                y_std = df.groupby('Year', as_index=False).sum()[columns].std(axis=1)
                # Update the labels accordingly if a sum.
                plot_options['y_label'] = y_label.replace('/month', '/year')
                # Join the annual time series data from the current set into the larger DataFrame.
                df_all_annual_time_series = pd.concat([df_all_annual_time_series, df.groupby('Year', as_index=False).sum()], axis=1)

            # If plotting a percent difference with respect to the first file set, store the data for that set and calculate the percent difference.
            if plot_percent_difference:
                if file_set_index == 0:
                    reference_data['annual'] = y
                # Add a tiny number to avoid divide-by-zero errors.
                y = (y - reference_data['annual'])/(reference_data['annual'] + EPSILON)*100
                y_std = y_std/(reference_data['annual'] + EPSILON)*100

            # Plot the annual time series, including possibly the error bars. 
            x = df.groupby('Year', as_index=False).mean()['Year']
            # Do not include the first file set if plotting a percent difference. Include the data in the plot otherwise.
            if not plot_percent_difference or file_set_index != 0:
                ax.plot(x, y, label=output_label, color=line_color, linestyle=linestyle_tuples[0][1], linewidth=linewidth)
                if std_annual_multiplier:
                    error = y_std*std_annual_multiplier
                    ax.fill_between(x, y-error, y+error, color=plot_colors[file_set_index], alpha=error_bars_alpha)
            
            # Add time series for one or more seasonal averages if specified to do so.
            if calculation_type == 'mean' and any(include_seasons.values()):
                plot_seasons(include_seasons, ax, df, x, output_label, line_color, linestyle_tuples, linewidth, columns=columns, 
                             reference_data=reference_data, file_or_file_set_index=file_set_index, std_multiplier=std_seasons_multiplier, 
                             error_bars_alpha=error_bars_alpha)

            # Create a separate time series plot for the seasonal averages if specified to do so. Set name of this seasons-only plot accordingly.
            if calculation_type == 'mean' and any(seasons_to_plot_separately.values()):
                plot_seasons(seasons_to_plot_separately, ax_seasons, df, x, output_label, line_color, linestyle_tuples, linewidth, 
                            columns=columns, reference_data=reference_data, file_or_file_set_index=file_set_index, 
                            std_multiplier=std_seasons_multiplier, error_bars_alpha=error_bars_alpha)
                plot_options['name'] = plot_name + '_seasons'   
                set_figure_options(fig_seasons, ax_seasons, plot_options)

            # Create the monthly time series plot for the variable if specified to do so.
            if monthly_time_series_plot:
                df, columns = read_file_set_into_single_variable_dataframe(output_files, file_set_index, variable, start_year, end_year, multiplier)
                y = df.groupby('Month', as_index=False).mean()[columns].mean(axis=1)
                y_std = df.groupby('Month', as_index=False).mean()[columns].std(axis=1)
                x = (df.groupby('Month', as_index=False).mean()['Month']).to_numpy()
                if plot_percent_difference:
                    if file_set_index == 0:
                        reference_data['monthly'] = y
                    y = (y - reference_data['monthly'])/(reference_data['monthly'] + EPSILON)*100
                    y_std = y_std/(reference_data['monthly'] + EPSILON)*100
                if not plot_percent_difference or file_set_index != 0:
                    ax_monthly.plot(x, y, label=output_label, color=line_color, linestyle=linestyle_tuples[0][1], linewidth=linewidth)
                    if std_monthly_multiplier:
                        error = y_std*std_monthly_multiplier
                        ax_monthly.fill_between(x, y-error, y+error, color=line_color, alpha=error_bars_alpha)
                # Join the monthly time series data from the current file into the larger DataFrame.
                df_all_monthly_time_series = pd.concat([df_all_monthly_time_series, df.groupby('Month', as_index=False).mean()], axis=1)

        # Calculate the mean and error bars for each of the annual and monthly time series data set groups. Perform multiple different t-tests.
        conditions = [include_annual_mean_across_all_sets, monthly_time_series_plot and include_monthly_mean_across_all_sets]
        dataframes = [df_all_annual_time_series, df_all_monthly_time_series]
        time_columns = ['Year', 'Month']
        axes = [ax, ax_monthly]
        std_multipliers = [std_annual_mean_across_all_sets_multiplier, std_monthly_mean_across_all_sets_multiplier]
        for time_index, condition in enumerate(conditions):
            df = dataframes[time_index]
            if not df.empty:
                time_column = time_columns[time_index]
                axis = axes[time_index]

                # Remove the duplicate 'Year' or 'Month' columns.
                df = df.loc[:,~df.columns.duplicated()]
                # Create DataFrame to store the means of each of the data set groups.
                df_all_data_set_means = pd.DataFrame()

                # Get all columns (except for the time) and find the columns for the first (assumed to be control) set. Calculate control set mean.
                columns = get_matching_column_in_dataframe(df, variable, all_matches=True)
                columns_control_set = [column for column in columns if column.endswith(f'_0')]
                df_control_mean = df[columns_control_set].mean(axis=1)

                # For each of the non-control set, perform a t-test at each individual time period against the control, as well as an overall t-test.
                for file_set_index in range(1, num_file_sets):

                    # Calculate the mean of the current set.
                    columns_this_set = [column for column in columns if column.endswith(f'_{file_set_index}')]
                    df_this_set_mean = df[columns_this_set].mean(axis=1)

                    # Perform a t-test to compare the means of the control against the current data set taken as a whole over the entire time period.
                    ttest = stats.ttest_ind(df_control_mean, df_this_set_mean)
                    output_label = output_labels[file_set_index]
                    print_p_values(ttest, variable, p_value_threshold, p_value_file, output_label, p_value_file_print_only_if_below_threshold)

                    # Perform a t-test to compare the control against the current data set at each time period (each year or each month).
                    p_values = df.apply(perform_ttest, columns_set_1=columns_control_set, columns_set_2=columns_this_set, axis=1).fillna(1)
                    mask = p_values <= p_value_threshold
                    color = plot_colors[file_set_index]
                    if plot_percent_difference:
                        df_this_set_mean = (df_this_set_mean - df_control_mean)/df_control_mean*100
                    df_all_data_set_means = pd.concat([df_all_data_set_means, df_this_set_mean], axis=1)
                    axis.plot(df[time_column][mask], df_this_set_mean[mask], marker=p_value_marker, markersize=p_value_marker_size, 
                              linestyle='None', color=color)
                    
                if plot_percent_difference:
                    # If plotting a percent difference with respect to the first file (control) set, the data for the control should be set to zeros.
                    df_control_mean -= df_control_mean
                else:
                    # If not plotting a percent difference, include the control set in the ensemble mean over all data sets.
                    df_all_data_set_means = pd.concat([df_all_data_set_means, df_control_mean], axis=1)

                # Calculate the ensemble mean, which is the mean of each of the data set means, and the error bars on this ensemble mean.
                if condition:
                    x = df[time_column]
                    y = df_all_data_set_means.mean(axis=1)
                    axis.plot(x, y, label='Mean', color='k', linestyle=linestyle_tuples[0][1], linewidth=linewidth)
                    std_multiplier = std_multipliers[time_index]
                    if std_multiplier:
                        error = df_all_data_set_means.std(axis=1)*std_multiplier
                        axis.fill_between(x, y-error, y+error, color='k', alpha=error_bars_alpha)

    # Finalize the annual time series plot for the variable now that all output files for the variable have been processed.
    plot_options['name'] = plot_name
    set_figure_options(fig, ax, plot_options)
    
    # Finalize the monthly time series plot for the variable now that all output files for the variable have been processed.
    if monthly_time_series_plot:
        plot_options['name'] = plot_name + '_monthly'
        plot_options['x_label'] = 'Month'
        plot_options['x_limits'] = monthly_time_series_x_limits
        plot_options['y_limits'] = monthly_time_series_y_limits
        set_figure_options(fig_monthly, ax_monthly, plot_options)

    # Close all the plots for the variable now that we are done with it. Record the elapsed time.
    plt.close(fig)
    plt.close(fig_seasons)
    plt.close(fig_monthly)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing plots for {variable} in {plot_directory}: {elapsed_time:.2f} seconds")


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file(s) on the command line.
    start_time = time.time()
    if len(sys.argv) < 2:
        print('Usage: plot_spatial_data.py `path/to/json/input/file(s)\'')
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
        list_of_inputs.extend(process_inputs(inputs[index]))

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