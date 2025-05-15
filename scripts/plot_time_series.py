import itertools
import json
from matplotlib import pyplot as plt
import multiprocessing
import numpy as np
import os
import pandas as pd
import sys
import time
from utility_constants import *
from utility_dataframes import get_columns_without_units_in_dataframe, get_matching_column_in_dataframe, read_file_into_dataframe
from utility_functions import check_is_list_of_lists
from utility_plots import *

""" Dictionary of default input values for time series plots. """
default_inputs_time_series = {'plot_directory': './',
                    'calculation_type': 'mean',
                    'plot_type': 'ensemble_means',
                    'multiplier': 1,
                    'std_annual_multiplier': None, 
                    'std_seasons_multiplier': None,
                    'std_monthly_multiplier': None,
                    'error_bars_alpha': 0.2,
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
    # If the user entered a string to indicate a single variable, put that string in a list.
    if isinstance(variables, str):
        variables = [variables]
        inputs['variables'] = variables

    # For plotting options that have not been specified in the inputs dictionary, add keys for them and use the default values.
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
        if key != 'variables':
            if not isinstance(value, dict):
                # If we have a single output file and corresponding label, and these are specified as strings, turn them into a list with this string.
                if (key == 'output_files' or key == 'output_labels') and isinstance(value, str):
                    value = [value]
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
                 std_multiplier=None, error_bars_alpha=None):
    """ 
    Adds seasons (as specified in the include_seasons dictionary) to a given time series plot. 

    Parameters:
        include_seasons: Dictionary specifying the seasons to include in the plot.
        ax: Axis object for the plot.
        df: Pandas DataFrame containing the time series data.
        columns: String (if a single column) or list (multiple columns) for the variable in the DataFrame that we want to plot.
        output_label: Base output label.
        plot_colors: di

    Returns:
        N/A.
    """
    seasons = ['MAM', 'JJA', 'SON', 'DJF']
    months_geq = [3, 6, 9, 12]
    months_leq = [5, 8, 11, 2]
    linestyle_tuples = [linestyle_tuples[2][1], linestyle_tuples[2][1], linestyle_tuples[1][1], linestyle_tuples[1][1]]
    linewidths = [linewidth+2, linewidth, linewidth+2, linewidth]
    for index in range(len(seasons)):
        season = seasons[index]
        if include_seasons[season]:

            month_geq = months_geq[index]
            month_leq = months_leq[index]
            linestyle=linestyle_tuples[index]
            linewidth=linewidths[index]
            
            if season == 'DJF':
                condition = (df['Month'] >= month_geq) | (df['Month'] <= month_leq)
            else:
                condition = (df['Month'] >= month_geq) & (df['Month'] <= month_leq)

            if isinstance(columns, list):
                y = df[condition].groupby('Year', as_index=False).mean()[columns].mean(axis=1)
                y_std = df[condition].groupby('Year', as_index=False).mean()[columns].std(axis=1)
            else:
                y = df[condition].groupby('Year', as_index=False).mean()[columns]

            ax.plot(x, y, label=output_label + f' ({season})', color=plot_colors, linestyle=linestyle, linewidth=linewidth)
            if std_multiplier and isinstance(columns, list):
                error = y_std*std_multiplier
                ax.fill_between(x, y-error, y+error, color=plot_colors, alpha=error_bars_alpha)

def read_file_into_single_variable_dataframe(file, variable, start_year, end_year, multiplier, return_column=False):
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
    
def read_file_into_single_variable_dataframe_ensemble(output_files, file_index, variable, start_year, end_year, multiplier):
    file = output_files[0][file_index]
    num_file_sets = len(output_files)
    df, column = read_file_into_single_variable_dataframe(file, variable, start_year, end_year, multiplier, return_column=True)
    for file_set_index in range(num_file_sets):
        file = output_files[file_set_index][file_index]
        df_for_this_file = read_file_into_single_variable_dataframe(file, variable, start_year, end_year, multiplier)
        df_for_this_file = df_for_this_file.rename(columns={column: column+f'_{file_set_index+1}'})
        df = pd.merge(df, df_for_this_file, on=['Year', 'Month'], how='inner')
    columns = get_matching_column_in_dataframe(df, variable, all_matches=True)
    return df, columns

def plot_time_series(inputs):
    """ 
    Parses a dictionary of inputs (keys are options, values are choices for those options) to create time series plots for a single variable. 
    Time series plots can be annual, with time presented in years (including seasonal variants that can be plotted separately), or 
    monthly, in which the values that are plotted indicate averages between specified start and end years for each month.
    If the output files are specified by a list of strings, the plots will show the direct values of the variables.
    If the output files are specified by a list of lists, representing an ensemble collection, the plots will show differences between ensemble averages.

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
    multiplier = inputs['multiplier']
    std_annual_multiplier = inputs['std_annual_multiplier']
    std_seasons_multiplier = inputs['std_seasons_multiplier']
    std_monthly_multiplier = inputs['std_monthly_multiplier']
    error_bars_alpha = inputs['error_bars_alpha']
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

    # Direct value plots in which each such time series plot can include one or more sets of curves, with one set for each output file.
    if not check_is_list_of_lists(output_files) or plot_type == 'direct' or (check_is_list_of_lists(output_files) and plot_type == 'ensemble_means'):
        
        # If output_files is a list of lists, unpack it into a list of strings, where each string specifies an output file.
        if check_is_list_of_lists(output_files) and plot_type == 'direct':
            # Turn the corresponding output_labels and plot_colors to list of lists with enough rows to match output_labels, then unpack.
            output_labels = [output_labels]*len(output_files)
            output_labels = list(itertools.chain.from_iterable(output_labels))
            plot_colors = [plot_colors]*len(output_files)
            plot_colors = list(itertools.chain.from_iterable(plot_colors))
            output_files = list(itertools.chain.from_iterable(output_files))

        if check_is_list_of_lists(output_files) and plot_type == 'ensemble_means':
            num_files_in_each_set = len(output_files[0])
        else:
            num_files_in_each_set = len(output_files)

        df_all_annual_time_series = pd.DataFrame()
        df_all_monthly_time_series = pd.DataFrame()

        for file_index in range(num_files_in_each_set):

            # Get the label for this output file and restrict between user-specified start and end years.
            output_label = output_labels[file_index]
            line_color = plot_colors[file_index]

            if not check_is_list_of_lists(output_files) or plot_type == 'direct':
                file = output_files[file_index]
                df, column = read_file_into_single_variable_dataframe(file, variable, start_year, end_year, multiplier, return_column=True)
                if calculation_type == 'mean':
                    y = df.groupby('Year', as_index=False).mean()[column]
                elif calculation_type == 'sum':
                    y = df.groupby('Year', as_index=False).sum()[column]
                    plot_options['y_label'] = y_label.replace('/month', '/year')
                
            elif check_is_list_of_lists(output_files) and plot_type == 'ensemble_means':
                df, columns = read_file_into_single_variable_dataframe_ensemble(output_files, file_index, variable, start_year, end_year, multiplier)
                if calculation_type == 'mean':
                    y = df.groupby('Year', as_index=False).mean()[columns].mean(axis=1)
                    y_std = df.groupby('Year', as_index=False).mean()[columns].std(axis=1)
                elif calculation_type == 'sum':
                    y = df.groupby('Year', as_index=False).sum()[columns].mean(axis=1)
                    y_std = df.groupby('Year', as_index=False).sum()[columns].std(axis=1)
                    plot_options['y_label'] = y_label.replace('/month', '/year')
            
            # Collect and plot the annual time series. Time series data can be either averaged or summed over each month of the year.
            x = df.groupby('Year', as_index=False).mean()['Year']
            ax.plot(x, y, label=output_label, color=line_color, linestyle=linestyle_tuples[0][1], linewidth=linewidth)
            if file_index == 0:
                x_series = pd.Series(x, name='Year')
                df_all_annual_time_series = pd.concat([df_all_annual_time_series, x_series], axis=1)
            y_series = pd.Series(y, name=f'{variable}_{file_index}')
            df_all_annual_time_series = pd.concat([df_all_annual_time_series, y_series], axis=1)
            if std_annual_multiplier and (check_is_list_of_lists(output_files) and plot_type == 'ensemble_means'):
                error = y_std*std_annual_multiplier
                ax.fill_between(x, y-error, y+error, color=plot_colors[file_index], alpha=error_bars_alpha)

            # Add time series for one or more seasonal averages if specified to do so.
            if calculation_type == 'mean' and any(include_seasons.values()):
                if not check_is_list_of_lists(output_files) or plot_type == 'direct':
                    plot_seasons(include_seasons, ax, df, x, output_label, line_color, linestyle_tuples, linewidth, columns=column)
                elif check_is_list_of_lists(output_files) and plot_type == 'ensemble_means':
                    plot_seasons(include_seasons, ax, df, x, output_label, line_color, linestyle_tuples, linewidth, columns=columns, 
                                 std_multiplier=std_seasons_multiplier, error_bars_alpha=error_bars_alpha)

            # Create a separate time series plot for the seasonal averages if specified to do so. Set name of this seasons-only plot accordingly.
            if calculation_type == 'mean' and any(seasons_to_plot_separately.values()):
                #print(f'{variable} {seasons_to_plot_separately}')
                if not check_is_list_of_lists(output_files) or plot_type == 'direct':
                    plot_seasons(seasons_to_plot_separately, ax_seasons, df, x, output_label, line_color, linestyle_tuples, linewidth, columns=column)
                elif check_is_list_of_lists(output_files) and plot_type == 'ensemble_means':
                    plot_seasons(seasons_to_plot_separately, ax_seasons, df, x, output_label, line_color, linestyle_tuples, linewidth, 
                                 columns=columns, std_multiplier=std_seasons_multiplier, error_bars_alpha=error_bars_alpha)
                plot_options['name'] = plot_name + '_seasons'   
                set_figure_options(fig_seasons, ax_seasons, plot_options)

            # Create the monthly time series plot for the variable if specified to do so.
            if monthly_time_series_plot:
                if not check_is_list_of_lists(output_files) or plot_type == 'direct':
                    df = read_file_into_single_variable_dataframe(file, variable, monthly_time_series_start_year, 
                                                                    monthly_time_series_end_year, multiplier, return_column=False)
                    y = df.groupby('Month', as_index=False).mean()[column]
                    x = (df.groupby('Month', as_index=False).mean()['Month']).to_numpy()
                    #x = np.vectorize(MONTH_NUM_TO_NAME.get)(x)
                    ax_monthly.plot(x, y, label=output_label, color=line_color, linestyle=linestyle_tuples[0][1], linewidth=linewidth)
                    #ax_monthly.tick_params(axis='x', labelrotation=45)
                elif check_is_list_of_lists(output_files) and plot_type == 'ensemble_means':
                    df, columns = read_file_into_single_variable_dataframe_ensemble(output_files, file_index, variable, start_year, end_year, multiplier)
                    y = df.groupby('Month', as_index=False).mean()[columns].mean(axis=1)
                    y_std = df.groupby('Month', as_index=False).mean()[columns].std(axis=1)
                    x = (df.groupby('Month', as_index=False).mean()['Month']).to_numpy()
                    ax_monthly.plot(x, y, label=output_label, color=line_color, linestyle=linestyle_tuples[0][1], linewidth=linewidth)
                    if std_monthly_multiplier:
                        error = y_std*std_monthly_multiplier
                        ax_monthly.fill_between(x, y-error, y+error, color=line_color, alpha=error_bars_alpha)
                
                if file_index == 0:
                    x_series = pd.Series(x, name='Month')
                    df_all_monthly_time_series = pd.concat([df_all_monthly_time_series, x_series], axis=1)
                y_series = pd.Series(y, name=f'{variable}_{file_index}')
                df_all_monthly_time_series = pd.concat([df_all_monthly_time_series, y_series], axis=1)

        if include_annual_mean_across_all_sets:
            columns = [column for column in df_all_annual_time_series.columns if column != 'Year']
            x = df_all_annual_time_series['Year']
            y = df_all_annual_time_series[columns].mean(axis=1)
            ax.plot(x, y, label='Mean', color='k', linestyle=linestyle_tuples[0][1], linewidth=linewidth)
            if std_annual_mean_across_all_sets_multiplier:
                error = df_all_annual_time_series[columns].std(axis=1)*std_annual_mean_across_all_sets_multiplier
                ax.fill_between(x, y-error, y+error, color='k', alpha=error_bars_alpha)

        if monthly_time_series_plot and include_monthly_mean_across_all_sets:
            columns = [column for column in df_all_monthly_time_series.columns if column != 'Month']
            x = df_all_monthly_time_series['Month']
            y = df_all_monthly_time_series[columns].mean(axis=1)
            ax_monthly.plot(x, y, label='Mean', color='k', linestyle=linestyle_tuples[0][1], linewidth=linewidth)
            if std_monthly_mean_across_all_sets_multiplier:
                error = df_all_monthly_time_series[columns].std(axis=1)*std_monthly_mean_across_all_sets_multiplier
                ax_monthly.fill_between(x, y-error, y+error, color='k', alpha=error_bars_alpha)

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

    # Ensemble average plots where the file sets are contained in a list of lists, and each inner list represents one set.

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

    # Process each dictionary to produce two lists of smaller dictionaries, where each smaller dictionary specifies options for a single plot.
    start_time = time.time()
    list_of_inputs = []
    for index in range(len(inputs)):
        # Process the inputs to fill in missing plotting input choices with default values, etc., and add to the list of dictionaries.
        list_of_inputs.extend(process_inputs(inputs[index]))

    # Create all of the times series plots in parallel.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(plot_time_series, list_of_inputs)
    
    # Print the total execution time to produce all the plots.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all plots: {elapsed_time:.2f} seconds")