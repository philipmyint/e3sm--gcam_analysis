import cartopy.crs as ccrs
import json
from matplotlib import pyplot as plt
import multiprocessing
import numpy as np
import os
import pandas as pd
from scipy import stats
import sys
import time
import xarray as xr
from utility_constants import *
from utility_plots import *
from utility_xarray import calculate_mean_and_std_of_da_list, calculate_min_mean_max_std_of_da

""" Dictionary of default input values for spatial plots. """
default_inputs_spatial_data = {'plot_directory': './',
                    'calculation_type': 'mean',
                    'multiplier': 1,
                    'start_year': 2071,
                    'end_year': 2090,
                    'projection': ccrs.Robinson,
                    'cmap': 'bwr',
                    'width': width_default,
                    'height': height_default,
                    'cbar_on': True,
                    'cbar_label_size': tick_label_size_default,
                    'cbar_limits': None,
                    'cbar_x_offset': 0.04,
                    'title_size': axis_label_size_default, 
                    'use_latex': False,
                    'produce_png': False,
                    'bbox_inches': 'tight',
                    'min_mean_max_size': legend_label_size_default,
                    'p_value_threshold': 0.05,
                    'p_value_file': None,
                    'p_value_file_print_only_if_below_threshold': True,
                    'stippling_on': False,
                    'stippling_std_multiple': 2,
                    'stippling_hatches': 'xxxx'
}

def process_inputs(inputs):
    """ 
    Processes a dictionary of inputs (keys are options, values are choices for those options) for creating spatial plots from data in NetCDF files.

    Parameters:
        inputs: Dictionary containing the user plotting choice inputs for different options. This dictionary may be incomplete or have invalid values.

    Returns:
        List of dictionaries, each of which specifies all plotting options for a single variable. If the user did not select a plotting option for
        a particular variable, the default choice for that plotting option will be selected.
    """
    # Turn the list of NetCDF files and their corresponding labels into a list of lists, if they are not already in that form.
    for input_type in ['netcdf_files']:
        if isinstance(inputs[input_type], str):
            inputs[input_type] = [[inputs[input_type]]]
        elif isinstance(inputs[input_type], list) and isinstance(inputs[input_type][0], str):
            inputs[input_type] = [inputs[input_type]]

    # Read one of the NetCDF files into an xarray Dataset so that we can later get the variables contained in them and the units of these variables.
    ds = xr.open_dataset(inputs['netcdf_files'][0][0])

    # If the user entered the string 'all' for the variables or no input at all for the variables, assume that they want to make plots for
    # all variables that are in the Dataset.
    if 'variables' not in inputs:
        inputs['variables'] = None
    variables = inputs['variables']
    if not variables or variables == 'all':
        variables = list(ds.keys())
        inputs['variables'] = variables
    # If the user entered a string to indicate a single variable, put that string in a list.
    if isinstance(variables, str):
        variables = [variables]
        inputs['variables'] = variables

    # For the plotting options that have not been specified in the inputs dictionary, add keys for them if necessary and use the default values.
    for key in default_inputs_spatial_data.keys():
        if key not in inputs:
            inputs[key] = default_inputs_spatial_data[key]

    # If the user specified anything other than a dictionary (e.g., a single value [string, integer, float] or a list) for the other plotting options, 
    # assume that they want to use that value/list for all the variables. Enable this by creating dictionaries with the keys given by the variables. 
    # This also covers the case above, where default values were added for all plotting options that are not specified in the inputs dictionary.
    for key, value in inputs.items():
        if key != 'variables':
            if not isinstance(value, dict):
                inputs[key] = dict.fromkeys(variables, value)
    
    # For each variable, if a plotting option is missing (has not been specified), fill in with the default corresponding to that plotting option.
    if 'plot_name' not in inputs:
        inputs['plot_name'] = {}
    if 'title' not in inputs:
        inputs['title'] = {}
    for variable in variables:
        # Use the default for the plot directory and create the directory if it does not already exist.
        if not any(key == variable for key in inputs['plot_directory'].keys()):
            inputs['plot_directory'][variable] = default_inputs_spatial_data['plot_directory']
        if not os.path.exists(inputs['plot_directory'][variable]):
            os.makedirs(inputs['plot_directory'][variable])
        # Default for the plot names is to call it 'spatial_[var_name]', where '[var_name]' is the name of the variable.
        if not any(key == variable for key in inputs['plot_name'].keys()):
            inputs['plot_name'][variable] = os.path.join(inputs['plot_directory'][variable], 'spatial_' + variable)
        # Default for the title of a variable is to use the column header for that variable from the DataFrame.
        if not any(key == variable for key in inputs['title'].keys()):
            # Replace ^2 with $^2$ in the units for the variable, so that the intended exponent (e.g., m^2) gets rendered correctly.
            units = ds[variable].attrs['units'].replace('^2', '$^2$')
            inputs['title'][variable] = rf"{variable} ({units})"
        # Default for the other plotting options are specified in the default_inputs_spatial_data dictionary.
        for key, value in inputs.items():
            if key not in ['variables', 'plot_directory', 'plot_name', 'title']:
                if not any(value_key == variable for value_key in value.keys()):
                    inputs[key][variable] = default_inputs_spatial_data[key]

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

def plot_spatial_data_from_netcdf_files(inputs):
    """ 
    Parses a dictionary of inputs (keys are options, values are choices for those options) to create spatial plots for a single variable. 
    The data for these spatial plots are stored in NetCDF files specified by the inputs dictionary.

    Parameters:
        input: Dictionary containing the user plotting choice inputs for different options. This dictionary is assumed to be complete (pre-processed).

    Returns:
        N/A.
    """
    # This function creates spatial plots for a single variable, and so it assumes that there is only one variable in the inputs dictionary.
    start_time = time.time()
    variable = inputs['variable']

    # Extract all other plotting options.
    plot_directory = inputs['plot_directory']
    multiplier = inputs['multiplier']
    netcdf_files = inputs['netcdf_files']
    cmap_color = inputs['cmap']
    projection = inputs['projection']
    title = inputs['title'] 
    title_size = inputs['title_size']
    width = inputs['width'] 
    height = inputs['height']
    plot_name = inputs['plot_name']
    cbar_on = inputs['cbar_on']
    cbar_label_size = inputs['cbar_label_size']
    cbar_limits = inputs['cbar_limits']
    cbar_x_offset = inputs['cbar_x_offset']
    min_mean_max_size = inputs['min_mean_max_size']
    calculation_type = inputs['calculation_type']
    use_latex = inputs['use_latex']
    start_year = inputs['start_year']
    end_year = inputs['end_year']
    p_value_threshold = inputs['p_value_threshold']
    p_value_file = inputs['p_value_file']
    p_value_file_print_only_if_below_threshold = inputs['p_value_file_print_only_if_below_threshold']
    stippling_on = inputs['stippling_on']
    stippling_std_multiple = inputs['stippling_std_multiple']
    stippling_hatches = inputs['stippling_hatches']
    
    # Calculate either the mean or sum between the start and end years for each lat/lon coordinate, and display this mean or sum on the spatial plot.
    # Store these means or sums of the variable in dataArrays. The dataArrays correspond to the NetCDF files, which are arranged in a list of lists, 
    # where each of the inner nested lists can specify either one or two files. If they specify only one file, the spatial plots will show the value
    # for that file; if there are two files, the spatial plots will show values for the difference between the two files (e.g., a control and a test).
    # If the file set (outer list) contains more than one inner list, the plot will show the mean over all the inner lists (an ensemble average).
    num_file_sets = len(netcdf_files)
    num_files_in_each_set = len(netcdf_files[0])
    dataArrays_to_include_in_plot = []
    all_dataArrays = [[0 for _ in range(num_files_in_each_set)] for _ in range(num_file_sets)]
    for file_set_index in range(num_file_sets):
        for file_index in range(num_files_in_each_set):
            file = netcdf_files[file_set_index][file_index]
            if calculation_type == 'mean':
                da = xr.open_dataset(file).sel(year=slice(start_year, end_year)).mean(dim='year')[variable]*multiplier
            elif calculation_type == 'sum':
                da = xr.open_dataset(file).sel(year=slice(start_year, end_year)).sum(dim='year')[variable]*multiplier
                # If calculating the sum, we have to change the per-time quantities and their units accordingly.
                per_time_labels = ['/year', '/month', '/day', '/hour', '/min', '/s']
                time_multipliers = np.array([1, years_TO_months, years_TO_days, years_TO_hours, years_TO_mins, years_TO_s])
                for index, per_time_label in enumerate(per_time_labels):
                    if per_time_label in title:
                        title = title.replace(per_time_label, '')
                        da *= time_multipliers[index]
                        break
            all_dataArrays[file_set_index][file_index] = da
        if num_files_in_each_set > 1:
            da = all_dataArrays[file_set_index][1] - all_dataArrays[file_set_index][0]
        else:
            da = all_dataArrays[file_set_index][0]
        dataArrays_to_include_in_plot.append(da)

    # Calculate the average over all file sets to produce a mean dataArray, and calculate the statistics (e.g., min, max) of this mean dataArray.
    da = calculate_mean_and_std_of_da_list(dataArrays_to_include_in_plot, calculate_std=False)
    min, mean, max, std = calculate_min_mean_max_std_of_da(da, calculate_std=True)

    # If the inner list(s) contain 2 files, concatenate the dataArrays over each set of inner lists and perform a t-test on the concatenated dataArray.
    # The t-test is used to analyze whether there is a statistically significant difference between the 2 sets of files (e.g., control and test).
    if num_files_in_each_set > 1:
        df1 = all_dataArrays[0][0].to_dataframe().dropna()[variable]
        df2 = all_dataArrays[0][1].to_dataframe().dropna()[variable]
        for file_set_index in range(1, num_file_sets):
            df1 = pd.concat([df1, all_dataArrays[file_set_index][0].to_dataframe().dropna()[variable]])
            df2 = pd.concat([df2, all_dataArrays[file_set_index][1].to_dataframe().dropna()[variable]])
        ttest = stats.ttest_ind(df1, df2)
        # Print the p-values to the console and optionally print to an output file.
        if ttest.pvalue < p_value_threshold:
            print(f'p-value of {variable}: {ttest.pvalue:.4e}, which is less than {p_value_threshold}')
            if p_value_file:
                with open(p_value_file, 'a+') as f:
                    f.write(f'{plot_directory}: ({variable}, {ttest.pvalue:.4e})\n')
        else:
            print(f'p-value of {variable}: {ttest.pvalue:.4e}')
            if p_value_file and not p_value_file_print_only_if_below_threshold:
                with open(p_value_file, 'a+') as f:
                    f.write(f'{plot_directory}: ({variable}, {ttest.pvalue:.4e})\n')

    # Use LaTeX for the labels if specified to do so.
    if use_latex:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', weight='bold') 

    # Plot the DataArray and optionally add the title and colorbar.
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=projection())
    da_fig = da.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=plt.get_cmap(cmap_color), extend='both', add_colorbar=False)
    if title:
        plt.rcParams['axes.titlesize'] = title_size
        ax.set_title(title)
    if cbar_on:
        cbar_ax = fig.add_axes([ax.get_position().x1+cbar_x_offset, ax.get_position().y0, 0.02, ax.get_position().height])
        cbar = plt.colorbar(da_fig, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=cbar_label_size, length=0)
        if cbar_limits:
            cbar.ax.set_ylim(cbar_limits[0], cbar_limits[1])

    # Add stippling to indicate regions where the value is +/- some multiple of the standard deviation (default is 2*std) away from the mean.
    if stippling_on:
        mask = np.abs(da) >= mean + stippling_std_multiple*std
        ax.contourf(da.lon, da.lat, mask, levels=1, hatches=['', stippling_hatches], alpha=0, transform=ccrs.PlateCarree())
    
    # Display min, mean, max.
    ax.text(x=0.88, y=0.9, s=f'Max:{max:.2e}\nMean:{mean:.2e}\nMin:{min:.2e}', ha='left', fontsize=min_mean_max_size, transform=ax.transAxes)

    # Add additional features like coastline and oceans.
    ax.coastlines(lw=0.6)

    # Save the figure and then close it. Record the elapsed time.
    save_figure(plot_name, fig, inputs)
    plt.close(fig)
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
    for index in range(1, len(sys.argv)):
        input_file = sys.argv[index]
        with open(input_file) as f:
            inputs.extend(json.load(f))
    
    # Process each dictionary to produce a list of smaller dictionaries, where each smaller dictionary specifies options for a single plot.
    start_time = time.time()
    list_of_inputs_for_each_plot = []
    for index in range(len(inputs)):
        # Process the inputs to fill in missing plotting input choices with default values, etc., and add to the list of dictionaries.
        list_of_inputs_for_each_plot.extend(process_inputs(inputs[index]))

    # Create all of the spatial plots in parallel.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(plot_spatial_data_from_netcdf_files, list_of_inputs_for_each_plot)

    # Print the total execution time to produce all the plots.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all plots: {elapsed_time:.2f} seconds")