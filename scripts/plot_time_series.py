import json
from matplotlib import pyplot as plt
import multiprocessing
import numpy as np
import pandas as pd
import os
import sys
import time
from utility_constants import MONTH_NUM_TO_NAME
from utility_dataframes import get_columns_without_units_in_dataframe, get_matching_column_in_dataframe, read_file_into_dataframe
from utility_plots import default_inputs_time_series, set_figure_options

def process_inputs(inputs):
    """ 
    Processes a dictionary of inputs (keys are options, values are choices for those options) for creating time series plots.

    Parameters:
        inputs: Dictionary containing the user plotting choice inputs for different options. This dictionary may be incomplete or have invalid values.

    Returns:
        List of dictionaries, each of which specifies all plotting options for a single variable. If the user did not select a plotting option for
        a particular variable, the default choice for that plotting option will be selected.
    """
    # The user can specify output_files to be a string (indicating a single output file), a list of different output files, 
    # or a dictionary of output files. Read the output file (any one will do if there is more than one such file) and put into a Pandas DataFrame.
    if isinstance(inputs['output_files'], str):
        file = inputs['output_files']
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

    # For the plotting options that have not been specified in the inputs dictionary, add keys for them if necessary and use the default values.
    for key in default_inputs_time_series.keys():
        if key not in inputs:
            inputs[key] = default_inputs_time_series[key]

    # If the user specified only a single value (string, integer, float, or a single-element list) for the other plotting options, assume that they
    # want to use that value for all of the variables. Enable this by creating dictionaries with the keys given by the variables. This also covers the
    # case above, where default values are added for all plotting options that are not specified in the inputs dictionary.
    for key, value in inputs.items():
        if key != 'variables':
            if not isinstance(value, dict):
                # If we have a single output file and corresponding label, and these are specified as strings, turn them into a list with this string.
                if (key == 'output_files' or key == 'output_labels') and isinstance(value, str):
                    value = [value]
                inputs[key] = dict.fromkeys(variables, value)

    # For each variable, if a plotting option is missing (has not been specified), fill in with the default corresponding to that plotting option.
    if 'plot_names' not in inputs:
        inputs['plot_names'] = {}
    if 'y_labels' not in inputs:
        inputs['y_labels'] = {}
    for variable in variables:
        # Default for the plot names is to call it 'time_series_[var_name]', where '[var_name]' is the name of the variable.
        if not any(key == variable for key in inputs['plot_names'].keys()):
            inputs['plot_names'][variable] = os.path.join(inputs['plot_directories'][variable], 'time_series_' + variable)
        # Default for the y_label of a variable is to use the column header for that variable from the DataFrame.
        if not any(key == variable for key in inputs['y_labels'].keys()):
            inputs['y_labels'][variable] = get_matching_column_in_dataframe(df, variable)
        # Default for the other plotting options are specified in the default_inputs_time_series dictionary.
        for key, value in inputs.items():
            if key not in ['variables', 'plot_names', 'y_labels']:
                if not any(value_key == variable for value_key in value.keys()):
                    inputs[key][variable] = default_inputs_time_series[key]

    # Now that the dictionary has been populated with a complete plotting options for each variable, separate it into a list of dictionaries,
    # where each of these smaller dictionaries contain the plotting options for a single variable. Return this list of dictionaries.
    list_of_inputs_for_each_output_file = []
    for variable in inputs['variables']:
        inputs_for_this_variable = {'variables': variable}
        for input_option in inputs.keys():
            if input_option != 'variables':
                inputs_for_this_variable[input_option] = {variable: inputs[input_option][variable]}
        list_of_inputs_for_each_output_file.append(inputs_for_this_variable)
    return list_of_inputs_for_each_output_file

def plot_time_series(inputs):
    """ 
    Parses a dictionary of inputs (keys are options, values are choices for those options) to create time series plots for a single variable. 
    Time series plots can be annual, with time presented in years (including seasonal variants that can be plotted separately), or 
    monthly, in which the values that are plotted indicate averages between specified start and end years for each month.

    Parameters:
        input: Dictionary containing the user plotting choice inputs for different options. This dictionary is assumed to be complete (pre-processed).

    Returns:
        N/A.
    """
    # This function creates time series plots for a single variable, and so it assumes that there is only one variable in the inputs dictionary.
    start_time = time.time()
    variable = inputs['variables']

    # Get all plotting options for both annual and monthly time series plots.
    output_files = inputs['output_files'][variable]
    output_labels = inputs['output_labels'][variable]
    plot_directory = inputs['plot_directories'][variable]
    plot_name = os.path.join(plot_directory, inputs['plot_names'][variable])
    y_label = inputs['y_labels'][variable]
    calculation_type = inputs['calculation_types'][variable]
    multiplier = inputs['multipliers'][variable]
    start_year = inputs['start_years'][variable]
    end_year = inputs['end_years'][variable]
    width = inputs['widths'][variable]
    height = inputs['heights'][variable]
    x_scale = inputs['x_scales'][variable]
    y_scale = inputs['y_scales'][variable]
    x_limits = inputs['x_limits'][variable]
    y_limits = inputs['y_limits'][variable]
    x_tick_label_size = inputs['x_tick_label_sizes'][variable]
    y_tick_label_size = inputs['y_tick_label_sizes'][variable]
    x_label_size = inputs['x_label_sizes'][variable]
    y_label_size = inputs['y_label_sizes'][variable]
    legend_on = inputs['legends_on'][variable]
    legend_label_size = inputs['legend_label_sizes'][variable]
    linewidth = inputs['linewidths'][variable]
    plot_colors = inputs['plot_colors'][variable]
    linestyle_tuples = inputs['linestyle_tuples'][variable]
    use_latex = inputs['use_latex'][variable]
    include_seasons = inputs['include_seasons'][variable]
    seasons_to_plot_separately = inputs['seasons_to_plot_separately'][variable]

    # Set the plotting options
    plot_options = dict(width=width, height=height, name=plot_name, x_label='Year', y_label=fr'{y_label}')
    plot_options.update(zip(['x_scale', 'y_scale', 'x_limits', 'y_limits', 'use_latex'], [x_scale, y_scale, x_limits, y_limits, use_latex]))
    plot_options.update(zip(['x_tick_label_size', 'y_tick_label_size', 'legend_on'], [x_tick_label_size, y_tick_label_size, legend_on]))
    plot_options.update(zip(['x_label_size', 'y_label_size', 'legend_label_size'], [x_label_size, y_label_size, legend_label_size]))

    # Get plotting options relevant only to the monthly time series.
    monthly_time_series_plot = inputs['monthly_time_series_plot'][variable]
    monthly_time_series_start_year = inputs['monthly_time_series_start_year'][variable]
    monthly_time_series_end_year = inputs['monthly_time_series_end_year'][variable]
    monthly_time_series_x_limits = inputs['monthly_time_series_x_limits'][variable]
    monthly_time_series_y_limits = inputs['monthly_time_series_y_limits'][variable]

    # Figure and axis objects for annual time series (+ possibly including seasons), seasons separately, and monthly time series plots.
    fig, ax = plt.subplots(nrows=1, ncols=1)
    fig_seasons, ax_seasons = plt.subplots(nrows=1, ncols=1)
    fig_monthly, ax_monthly = plt.subplots(nrows=1, ncols=1)

    # Each time series plot can include one or more sets of curves, with one set for each output file.
    for file_index, file in enumerate(output_files):

        # Get the label for this output file, restrict between user-specified start and end years, find the column label for the variable.
        output_label = output_labels[file_index]
        df = read_file_into_dataframe(file)
        df = df[(df['Year'] >= start_year) & (df['Year'] <= end_year)]
        for column in df.columns:
            if variable in column:
                break
        
        # Collect and plot the annual time series. Time series data can be either averaged or summed over each month of the year.
        x = df.groupby('Year', as_index=False).mean()['Year']
        if calculation_type == 'mean':
            y = df.groupby('Year', as_index=False).mean()[column]*multiplier
            ax.plot(x, y, label=output_label, color=plot_colors[file_index], linestyle=linestyle_tuples[0][1], linewidth=linewidth)
            # Add time series for one or more seasonal averages if specified to do so.
            if any(include_seasons.values()):
                if include_seasons['MAM']:
                    y_season = df[(df['Month'] >= 3) & (df['Month'] <= 5)].groupby('Year', as_index=False).mean()[column]*multiplier
                    ax.plot(x, y_season, label=output_label + ' (MAM)', color=plot_colors[file_index], 
                            linestyle=linestyle_tuples[2][1], linewidth=linewidth+2)
                if include_seasons['JJA']:     
                    y_season = df[(df['Month'] >= 6) & (df['Month'] <= 8)].groupby('Year', as_index=False).mean()[column]*multiplier
                    ax.plot(x, y_season, label=output_label + ' (JJA)', color=plot_colors[file_index], 
                            linestyle=linestyle_tuples[2][1], linewidth=linewidth)
                if include_seasons['SON']:
                    y_season = df[(df['Month'] >= 9) & (df['Month'] <= 11)].groupby('Year', as_index=False).mean()[column]*multiplier
                    ax.plot(x, y_season, label=output_label + ' (SON)', color=plot_colors[file_index], 
                            linestyle=linestyle_tuples[1][1], linewidth=linewidth+2)
                if include_seasons['DJF']:
                    y_season = df[(df['Month'] == 12) | (df['Month'] <= 2)].groupby('Year', as_index=False).mean()[column]*multiplier
                    ax.plot(x, y_season, label=output_label + ' (DJF)', color=plot_colors[file_index], 
                            linestyle=linestyle_tuples[1][1], linewidth=linewidth)
            # Create a separate time series plot for the seasonal averages if specified to do so.
            if any(seasons_to_plot_separately.values()):
                if seasons_to_plot_separately['MAM']:
                    y_season = df[(df['Month'] >= 3) & (df['Month'] <= 5)].groupby('Year', as_index=False).mean()[column]*multiplier
                    ax_seasons.plot(x, y_season, label=output_label + ' (MAM)', color=plot_colors[file_index], 
                            linestyle=linestyle_tuples[2][1], linewidth=linewidth+2)
                if seasons_to_plot_separately['JJA']:     
                    y_season = df[(df['Month'] >= 6) & (df['Month'] <= 8)].groupby('Year', as_index=False).mean()[column]*multiplier
                    ax_seasons.plot(x, y_season, label=output_label + ' (JJA)', color=plot_colors[file_index], 
                            linestyle=linestyle_tuples[2][1], linewidth=linewidth)
                if seasons_to_plot_separately['SON']:
                    y_season = df[(df['Month'] >= 9) & (df['Month'] <= 11)].groupby('Year', as_index=False).mean()[column]*multiplier
                    ax_seasons.plot(x, y_season, label=output_label + ' (SON)', color=plot_colors[file_index], 
                            linestyle=linestyle_tuples[1][1], linewidth=linewidth+2)
                if seasons_to_plot_separately['DJF']:
                    y_season = df[(df['Month'] == 12) | (df['Month'] <= 2)].groupby('Year', as_index=False).mean()[column]*multiplier
                    ax_seasons.plot(x, y_season, label=output_label + ' (DJF)', color=plot_colors[file_index], 
                            linestyle=linestyle_tuples[1][1], linewidth=linewidth)
                # Plot the seasons separately, making sure to set the name of this seasons-only plot accordingly.
                plot_options['name'] = plot_name + '_seasons'   
                set_figure_options(fig_seasons, ax_seasons, plot_options)
        elif calculation_type == 'sum':
            # If the variable is a sum over the entire year, update the '/month' in the label to '/year' if it is a per-month quantity.
            y = df.groupby('Year', as_index=False).sum()[column]*multiplier
            plot_options['y_label'] = y_label.replace('/month', '/year')
            ax.plot(x, y, label=output_label, color=plot_colors[file_index], linestyle='solid', linewidth=linewidth)

        # Create the monthly time series plot for the variable if specified to do so.
        if monthly_time_series_plot:
            df = read_file_into_dataframe(file)
            df = df[(df['Year'] >= monthly_time_series_start_year) & (df['Year'] <= monthly_time_series_end_year)]
            y = df.groupby('Month', as_index=False).mean()[column]*multiplier
            x = (df.groupby('Month', as_index=False).mean()['Month']).to_numpy()
            #x = np.vectorize(MONTH_NUM_TO_NAME.get)(x)
            ax_monthly.plot(x, y, label=output_label, color=plot_colors[file_index], linestyle=linestyle_tuples[0][1], linewidth=linewidth)
            #ax_monthly.tick_params(axis='x', labelrotation=45)
            
    # Finalize the annual time series plot for tthe variable now that all output files for the variable have been processed.
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

    # Run this script together with the input JSON file on the command line.
    if len(sys.argv) != 2:
        print('Usage: plot_time_series.py `path/to/json/input/file\'')
        sys.exit()

    # Read and load the JSON file.
    input_file = sys.argv[1]
    with open(input_file) as f:
        inputs = json.load(f)

    # Process each block in the JSON file to produce a list of dictionaries, where each specifies time series plotting options for a single variable.
    start_time = time.time()
    list_of_inputs_for_each_plot = []
    for index in range(len(inputs)):
        # Process the inputs to fill in missing plotting input choices with default values, etc., and add to the list of dictionaries.
        list_of_inputs_for_each_plot.extend(process_inputs(inputs[index]))
    
    # Create all of the times series plots in parallel.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(plot_time_series, list_of_inputs_for_each_plot)
    
    # Print the total execution time.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all plots: {elapsed_time:.2f} seconds")