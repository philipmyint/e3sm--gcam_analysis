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
import uxarray as ux
import xarray as xr
from utility_constants import *
from utility_dataframes import perform_ttest
from utility_functions import check_is_list_of_lists, print_p_values, replace_inside_parentheses, sort_file, transpose_scenarios_if_needed
from utility_plots import *
from utility_xarray import calculate_statistics_of_xarray, convert_xarray_to_uxarray

""" Dictionary of default input values for spatial plots. """
default_inputs = {
    'cbar_label_size': tick_label_size_default,
    'cbar_limits': None,
    'cbar_on': True,
    'cbar_x_offset': 0.06,
    'cmap': 'bwr',
    'end_year': 2090,
    'grid_file': None,
    'height': height_default,
    'multiplier': 1,
    'p_value_file': "p_values.dat",
    'p_value_file_print_only_if_below_threshold': True,
    'p_value_threshold': 0.05,
    'plot_directory': './',
    'plot_type': 'absolute_difference',
    'produce_png': False,
    'projection': ccrs.Robinson,
    'start_year': 2071,
    'statistics_panel_size': legend_label_size_default,
    'stippling_hatches': 'xxxx',
    'stippling_on': False,
    'stippling_std_multiple': 2,
    'time_calculation': 'mean',
    'title_size': axis_label_size_default, 
    'use_latex': False,
    'width': width_default               
}

def process_inputs(inputs):
    """ 
    Processes a dictionary of inputs (keys are options, values are choices for those options) for creating spatial plots from data in NetCDF files.

    Parameters:
        inputs: Dictionary containing the user plotting choice inputs for different options. This dictionary may be incomplete or have invalid values.

    Returns:
        List of dictionaries, each of which specifies completed plotting options for a single variable. If the user did not select a plotting option 
        for a particular variable, the default choice for that plotting option will be selected.
    """
    # If the NetCDF files are specified as a list of lists (for ensemble plots), check if they need to be transposed.
    # Users can now specify output files in two formats:
    #   Format A (organized by ensemble member - original format):
    #       [["Control", "Full feedback"], ["Control_2", "Full feedback_2"], ...]
    #   Format B (organized by file set - new user-friendly format):
    #       [["Control", "Control_2", ...], ["Full feedback", "Full feedback_2", ...]]
    # The plotting functions expect Format A internally, so Format B will be automatically transposed.
    if check_is_list_of_lists(inputs['netcdf_files']):
        inputs['netcdf_files'], was_transposed = transpose_scenarios_if_needed(inputs['netcdf_files'])
        notify_netcdf_files_transposed = inputs.get('notify_netcdf_files_transposed', False)
        if was_transposed and notify_netcdf_files_transposed:
            print(f"Note: NetCDF files were automatically transposed from 'organized by file set' format to 'organized by ensemble member' format.")

    # Turn the list of NetCDF files and their corresponding labels into a list of lists, if they are not already in that form.
    if not check_is_list_of_lists(inputs['netcdf_files']):
        if isinstance(inputs['netcdf_files'], str):
            inputs['netcdf_files'] = [[inputs['netcdf_files']]]
        elif isinstance(inputs['netcdf_files'], list):
            inputs['netcdf_files'] = [inputs['netcdf_files']]

    # Collect all NetCDF file paths from the (potentially nested) list and check that each one exists.
    netcdf_files_flat = []
    raw_files = inputs['netcdf_files']
    if check_is_list_of_lists(raw_files):
        for sublist in raw_files:
            netcdf_files_flat.extend(sublist)
    elif isinstance(raw_files, list):
        netcdf_files_flat.extend(raw_files)
    else:
        netcdf_files_flat.append(raw_files)
    missing_files = [f for f in netcdf_files_flat if not os.path.exists(f)]
    if missing_files:
        raise FileNotFoundError(f"Error: the following NetCDF files were not found: {missing_files}")

    # Read one of the NetCDF files into an xarray Dataset so that we can later get the variables contained in them and the units of these variables.
    ds = xr.open_dataset(inputs['netcdf_files'][0][0] if check_is_list_of_lists(inputs['netcdf_files']) else
                         inputs['netcdf_files'][0] if isinstance(inputs['netcdf_files'], list) else
                         inputs['netcdf_files'])

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

    # Check that all requested variables exist in the first NetCDF file.
    missing_vars = [v for v in variables if v not in ds]
    if missing_vars:
        available = list(ds.keys())
        ds.close()
        raise KeyError(f"Error: the following variables were not found in '{netcdf_files_flat[0]}': {missing_vars}. "
                       f"Available variables: {available}")

    # For the plotting options that have not been specified in the inputs dictionary, add keys for them if necessary and use the default values.
    for key in default_inputs.keys():
        if key not in inputs:
            inputs[key] = default_inputs[key]

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
            inputs['plot_directory'][variable] = default_inputs['plot_directory']
        if not os.path.exists(inputs['plot_directory'][variable]):
            os.makedirs(inputs['plot_directory'][variable])
        # Default for the plot names is to call it 'spatial_[var_name]', where '[var_name]' is the name of the variable.
        if not any(key == variable for key in inputs['plot_name'].keys()):
            inputs['plot_name'][variable] = os.path.join(inputs['plot_directory'][variable], 'spatial_' + variable)
        # Append the plot directory with the name of the p-value file.
            inputs['p_value_file'][variable] = os.path.join(inputs['plot_directory'][variable], inputs['p_value_file'][variable])
        # Default for the title of a variable is to use the column header for that variable from the Dataset.
        if not any(key == variable for key in inputs['title'].keys()):
            # Replace ^2 with $^2$ in the units for the variable, so that the intended exponent (e.g., m^2) gets rendered correctly.
            units = ds[variable].attrs['units'].replace('^2', '$^2$')
            # Replace /m2 with /m$^2$ in the units for the variable, so that it gets rendered correctly.
            units = units.replace('/m2', '/m$^2$')
            inputs['title'][variable] = rf'{variable} ({units})'
        # Default for the other plotting options are specified in the default_inputs dictionary.
        for key, value in inputs.items():
            if key not in ['variables', 'plot_directory', 'plot_name', 'title']:
                if not any(value_key == variable for value_key in value.keys()):
                    inputs[key][variable] = default_inputs[key]

    # Close the dataset now that all variable and unit information has been extracted.
    ds.close()

    # Check that the grid file exists if one has been specified (required for EAM plots).
    grid_file = inputs['grid_file'][variables[0]] if isinstance(inputs['grid_file'], dict) else inputs['grid_file']
    if grid_file and not os.path.exists(grid_file):
        raise FileNotFoundError(f"Error: grid file not found: '{grid_file}'")

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

def plot_spatial_data_eam(inputs, grid_file):
    """ 
    Creates spatial plots and perform statistical analysis for a single variable from E3SM EAM outputs. 
    The data for these spatial plots are stored in NetCDF files specified by the inputs dictionary.

    Parameters:
        input: Dictionary containing user inputs for different plotting options, where the keys are options and values are choices for those options.
               This dictionary is assumed to be complete (pre-processed).
        grid_file: Path and name of the grid file for the EAM unstructured mesh.

    Returns:
        N/A.
    """
    # This function creates spatial plots for a single variable, and so it assumes that there is only one variable in the inputs dictionary.
    start_time = time.time()
    variable = inputs['variable']

    # Extract all other plotting options.
    cbar_label_size = inputs['cbar_label_size']
    cbar_limits = inputs['cbar_limits']
    cbar_on = inputs['cbar_on']
    cbar_x_offset = inputs['cbar_x_offset']
    cmap_color = inputs['cmap']
    end_year = inputs['end_year']
    height = inputs['height']
    multiplier = inputs['multiplier']
    netcdf_files = inputs['netcdf_files']
    p_value_file = inputs['p_value_file']
    p_value_file_print_only_if_below_threshold = inputs['p_value_file_print_only_if_below_threshold']
    p_value_threshold = inputs['p_value_threshold']
    plot_directory = inputs['plot_directory']
    plot_name = inputs['plot_name']
    plot_type = inputs['plot_type']
    projection = inputs['projection']
    statistics_panel_size = inputs['statistics_panel_size']
    start_year = inputs['start_year']
    stippling_hatches = inputs['stippling_hatches']
    stippling_on = inputs['stippling_on']
    stippling_std_multiple = inputs['stippling_std_multiple']
    time_calculation = inputs['time_calculation']
    title = inputs['title'] 
    title_size = inputs['title_size']
    use_latex = inputs['use_latex']
    width = inputs['width'] 
 
    # Store the grid file in an uxarray Dataset.
    grid = ux.open_grid(grid_file)

    # We either have individual plots, in which case there could be multiple files arranged like [[file1, file2, file3, ...]],
    # or we could have ensemble plots, in which case there are at most two data sets, but potentially multiple files in each of the sets.
    # The preprocessing done earlier has formatted the files into a list of lists so that they are in the form:
    #   For individual plots: [[file1, file2, file3, ...]]    
    #   For ensemble plots: [["Control", "Full feedback"], ["Control_2", "Full feedback_2"], ...]
    num_files_in_each_set = len(netcdf_files)
    num_file_sets = len(netcdf_files[0])
   
    # Note that unlike time series, it is not meaningful to have a spatial plot with more than two ensembles, 
    # which is why we are restricting the total number of data sets to 2 when producing ensemble plots.
    # Note that num_file_sets and num_files_in_each_set represent different things depending on whether there
    #    are one or two ensembles.
    # If there is one ensemble (of any length, including 1 file):
    #    num_file_sets = the number of files in the ensemble
    #    num_files_in_each_set = 1 (the number of ensembles)
    #    only mean or sum can be calculated across the files and plotted (no comparison)
    # If there are two ensembles (or more, but the code exits if there are more than two):
    #    num_file_sets = 2 (the number of ensembles)
    #    num_files_in_each_set = the number of files in each ensemble
    #    the ensemble means are calculated and can be plotted separately or as one of two difference comparisons

    if (num_files_in_each_set == 1) or (num_file_sets == 2):
        pass
    else:
        error_message = "Error: For spatial plots, there can be only one or two sets of netcdf_files, " \
            + "and if there are two sets there must be at least two files in each set (for ensemble plots)."
        raise ValueError(error_message)

    # Read each of the NetCDF output files, which are arranged in a list of lists (2D matrix), into an uxarray DataArray and then add each of these 
    # DataArrays to a single uxarray Dataset that will store the data from all of the files. To form the DataArrays, calculate either the mean or sum 
    # between the start and end years for each lat/lon coordinate. We will later display some function of this mean or sum on the spatial plot.
    # Use a plain xr.Dataset as the intermediate container to avoid uxarray constructor
    # incompatibility with xarray >= 2025.7.1; convert to uxarray only at the plotting step.
    uxds = xr.Dataset()
    for file_set_index in range(num_file_sets):
        for file_index in range(num_files_in_each_set):
            file = netcdf_files[file_index][file_set_index]
            # Check that the requested start and end years fall within the year range of the file before slicing,
            # so the error message can report the available range rather than checking an already-filtered dataset.
            ds_full = xr.open_dataset(file)
            file_start = int(ds_full['year'].min())
            file_end = int(ds_full['year'].max())
            if start_year < file_start or start_year > file_end:
                ds_full.close()
                raise ValueError(f"Plot start_year {start_year} not in '{file}' (available: {file_start}–{file_end}); update input json file")
            if end_year < file_start or end_year > file_end:
                ds_full.close()
                raise ValueError(f"Plot end_year {end_year} not in '{file}' (available: {file_start}–{file_end}); update input json file")

            # Create a temporary NetCDF file with data between only the start and end years.
            ds = ds_full.sel(year=slice(start_year, end_year))[variable]
            ds_full.close()
            wfile = os.path.join(plot_directory, f'temp_{variable}_{file_index}_{file_set_index}.nc')
            ds.to_netcdf(wfile, 'w')
            ds.close()
            # Load the temporal aggregate as a plain xr.DataArray rather than a ux.UxDataArray.
            # ux.UxDataset.from_xarray and ux.open_dataset both use a UxDataset constructor
            # that is incompatible with xarray >= 2025.7.1. Since uxarray is only needed at
            # the final plotting step (where convert_xarray_to_uxarray is called), we use
            # plain xarray DataArrays here and store them in an xr.Dataset (uxds).
            if time_calculation == 'mean':
                xr_ds = xr.open_dataset(wfile).mean(dim='year')
                uxda = xr_ds[variable]*multiplier
                xr_ds.close()
            elif time_calculation == 'sum':
                xr_ds = xr.open_dataset(wfile).sum(dim='year')
                uxda = xr_ds[variable]*multiplier
                xr_ds.close()
                # If calculating the sum, change the per-time quantities and their units accordingly.
                per_time_labels = ['/year', '/month', '/day', '/hour', '/min', '/s']
                time_multipliers = np.array([1, years_TO_months, years_TO_days, years_TO_hours, years_TO_mins, years_TO_s])
                for index, per_time_label in enumerate(per_time_labels):
                    if per_time_label in title:
                        title = title.replace(per_time_label, '')
                        uxda.attrs['units'] = uxda.attrs['units'].replace(per_time_label, '')
                        uxda *= time_multipliers[index]
                        break
            # If there is more than one data set, modify the labels so that we know which data set corresponds to which variable of the DataFrame.
            if num_file_sets >= 2:
                uxda = uxda.rename(f'{variable}_{file_index}_{file_set_index}')
            uxds[uxda.name] = uxda
            # Delete the temporary NetCDF file now that the data have been read.
            os.system(f'rm {wfile}')
    
    # Initialize list that will store all the uxarray DataArrays that we will want to plot for the variable.
    uxDataArrays_to_plot = []
    df = uxds.to_dataframe()

    # If we have a single data set (but potentially multiple files in this set), take the mean over all files
    # for each spatial coordinate. For a single file this is equivalent to just plotting that file's values directly.
    if num_files_in_each_set == 1:
        df = df.mean(axis=1)
        uxDataArrays_to_plot.append(convert_xarray_to_uxarray(df.to_xarray(), grid, variable=variable))
        print(f"Plotting single file set of {num_file_sets} file(s) for {variable}.")
    
    # If we have two data sets, we can plot either the absolute difference, percent difference, or the two data sets separately.
    elif num_file_sets == 2:
        # Take the mean over all files for each lat/lon coordinate in each data set.
        columns_control_set = [column for column in df.columns if column.endswith(f'_0')]
        columns_test_set = [column for column in df.columns if column.endswith(f'_1')]
        df_control_set = df[columns_control_set].mean(axis=1)
        df_test_set = df[columns_test_set].mean(axis=1)

        # Perform a t-test to compare the two spatial data sets as whole over all coordinates. Print the results to the console and to an output file.
        ttest = stats.ttest_ind(df_control_set, df_test_set)
        print_p_values(ttest, variable, p_value_threshold, p_value_file, plot_name, p_value_file_print_only_if_below_threshold)

        # Plot either an absolute difference, percent difference, or the two data sets separately.
        # If we have more than one file per data set (meaning that we have an ensemble), we will be examining the ensemble means in all cases.
        if plot_type == 'absolute_difference':
            # Plot absolute differences between the two data sets. 
            df = df_test_set - df_control_set
            uxDataArrays_to_plot.append(convert_xarray_to_uxarray(df.to_xarray(), grid, variable=variable))
            print(f"Calculating {plot_type} between {num_file_sets} ensembles with {num_files_in_each_set} files in each ensemble.")
        elif plot_type == 'percent_difference':
            # Plot percent differences between the two data sets. Add a tiny number to avoid a divide-by-zero error. Take the absolute value
            # so that if the control is negative, while the test set is positive, we get a positive value for the percent difference.
            df = ((df_test_set - df_control_set)/(df_control_set.abs() + EPSILON))*100
            uxDataArrays_to_plot.append(convert_xarray_to_uxarray(df.to_xarray(), grid, variable=variable))
            title = replace_inside_parentheses(title, rf'($\%$ difference)')
            print(f"Calculating {plot_type} between {num_file_sets} ensembles with {num_files_in_each_set} files in each ensemble.")
        elif plot_type == 'separate_plots':
            # Plot the two data sets individually in their own separate plots. 
            uxDataArrays_to_plot.append(convert_xarray_to_uxarray(df_control_set.to_xarray(), grid, variable=variable))
            uxDataArrays_to_plot.append(convert_xarray_to_uxarray(df_test_set.to_xarray(), grid, variable=variable))
            print(f"Calculating {plot_type} of {num_file_sets} ensembles with {num_files_in_each_set} files in each ensemble.")

    # This is when there are more than three files listed in one input set, so differences don't make sense.
    else:
        error_message = (f"Error: variable '{variable}', num_files_in_each_set={num_files_in_each_set}, "
                         f"num_file_sets={num_file_sets}, plot_type='{plot_type}': "
                         "If there are more than two netcdf_files in one single input set "
                         "plot_type must be mean or sum.")
        raise ValueError(error_message)

    # If stippling_on is True, meaning we want to add markers on plot to indicate potential regions of statistical significance, 
    if stippling_on:
        uxds['lat'] = grid.face_lat
        uxds['lon'] = grid.face_lon
        # If there is more than one file per data set (meaning that we have an ensemble) and we do not want separate plots,
        # we can compare the two data sets by performing a t-test at each individual lat/lon coordinate and later adding stippling.
        if num_files_in_each_set >= 2 and plot_type != 'separate_plots':
            df = uxds.to_dataframe()
            df = df.groupby(['lat', 'lon']).mean()
            da_pvalues = df.apply(perform_ttest, columns_set_1=columns_control_set, \
                                columns_set_2=columns_test_set, axis=1).fillna(1).to_xarray().fillna(1)
            #uxda_pvalues = convert_xarray_to_uxarray(da_pvalues, grid, variable=variable, fillna=1)

    # Iterate over all uxDataArrays in the list to create a plot for each one.
    for uxda_index, uxda in enumerate(uxDataArrays_to_plot):

        # Calculate some basic statistics of the current uxDataArray.
        min, mean, median, max, std = calculate_statistics_of_xarray(uxda, variable)

        # Use LaTeX for the labels if specified to do so.
        if use_latex:
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif', weight='bold') 
        
        # If plotting a percent difference and no colorbar limits are specified, set them to be -100 and 100 percent if the max exceeds 100 percent.
        if plot_type == 'percent_difference' and not cbar_limits and max > 100:
            cbar_limits = [-100, 100]

        # Plot the uxDataArray and optionally add the title and colorbar.
        fig = plt.figure(figsize=(width, height))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=projection())
        uxda = uxda[variable]
        uxda_fig = uxda.to_polycollection(cache=True)
        uxda_fig.set_transform(ccrs.PlateCarree())
        uxda_fig.set_cmap(cmap_color)
        ax.set_title(title)
        ax.add_collection(uxda_fig)
        ax.set_global()
        if title:
            plt.rcParams['axes.titlesize'] = title_size
            ax.set_title(title)
        if cbar_on:
            # Extra for the colorbar padding in case using the default (not LaTeX) font.
            if not use_latex and cbar_x_offset < 0.09:
                cbar_x_offset = 0.09
            cbar_ax = fig.add_axes([ax.get_position().x1+cbar_x_offset, ax.get_position().y0, 0.02, ax.get_position().height])
            if cbar_limits:
                uxda_fig.set_clim(vmin=cbar_limits[0], vmax=cbar_limits[1])
            cbar = plt.colorbar(uxda_fig, cax=cbar_ax)
            cbar.ax.tick_params(labelsize=cbar_label_size, length=0)

        # Add stippling to indicate regions of potential statistical significance (this is very slow, so not really practical at the moment). 
        if stippling_on:
            plt.rcParams['hatch.linewidth'] = 0.5
            plt.rcParams['hatch.color'] = 'gray'
            if num_file_sets == 2 and num_files_in_each_set >= 2 and plot_type != 'separate_plots':
                # If there are two data sets and at least two files in each data set so that we will have previously calculated p-values at each
                # lat/lon coordinate, the stippling will indicate regions where the p-value is less than the designated threshold.
                mask = da_pvalues <= p_value_threshold
                #print(mask)
                ax.contourf(da_pvalues.lon, da_pvalues.lat, mask, levels=1, hatches=['', stippling_hatches], 
                            alpha=0, transform=ccrs.PlateCarree())
            else:
                # For all other cases, the stippling will indicate regions where the value is +/- some multiple of the standard deviation 
                # (default is 2*std) away from the mean.
                mask = np.abs(uxda) >= mean + stippling_std_multiple*std
                #print(uxds)
                ax.contourf(uxds['lon'], uxds['lat'], mask, levels=1, hatches=['', stippling_hatches], alpha=0, transform=ccrs.PlateCarree())
        
        # Display statistics.
        ax.text(x=0.88, y=0.9, s=f'Max:{max:.2e}\nMean:{mean:.2e}\nMedian:{median:.2e}\nMin:{min:.2e}', ha='left', 
                fontsize=statistics_panel_size, transform=ax.transAxes)
        ax.text(x=0.88, y=0.05, s=f'Std:{std:.2e}', ha='left', fontsize=statistics_panel_size, transform=ax.transAxes)

        # Add additional features like coastline and oceans.
        ax.coastlines(lw=0.6)

        # Save the figure and then close it. Record the elapsed time.
        end_time = time.time()
        elapsed_time = end_time - start_time

        if len(uxDataArrays_to_plot) > 1:
            # If we want to plot two data sets separately, add a set number to the figure name for each set.
            save_figure(plot_name + f'_set_{uxda_index+1}', fig, inputs)
            print(f"Elapsed time for producing plots for {variable} (set {uxda_index+1}) in {plot_directory}: {elapsed_time:.2f} seconds") 
        else:
            save_figure(plot_name, fig, inputs)
            print(f"Elapsed time for producing plots for {variable} in {plot_directory}: {elapsed_time:.2f} seconds") 
        plt.close(fig)

def plot_spatial_data_elm(inputs):
    """ 
    Creates spatial plots and perform statistical analysis for a single variable from E3SM ELM outputs. 
    The data for these spatial plots are stored in NetCDF files specified by the inputs dictionary.

    Parameters:
        input: Dictionary containing user inputs for different plotting options, where the keys are options and values are choices for those options.
               This dictionary is assumed to be complete (pre-processed).

    Returns:
        N/A.
    """
    # This function creates spatial plots for a single variable, and so it assumes that there is only one variable in the inputs dictionary.
    start_time = time.time()
    variable = inputs['variable']

    # Extract all other plotting options.
    cbar_label_size = inputs['cbar_label_size']
    cbar_limits = inputs['cbar_limits']
    cbar_on = inputs['cbar_on']
    cbar_x_offset = inputs['cbar_x_offset']
    cmap_color = inputs['cmap']
    end_year = inputs['end_year']
    height = inputs['height']
    multiplier = inputs['multiplier']
    netcdf_files = inputs['netcdf_files']
    p_value_file = inputs['p_value_file']
    p_value_file_print_only_if_below_threshold = inputs['p_value_file_print_only_if_below_threshold']
    p_value_threshold = inputs['p_value_threshold']
    plot_directory = inputs['plot_directory']
    plot_name = inputs['plot_name']
    plot_type = inputs['plot_type']
    projection = inputs['projection']
    statistics_panel_size = inputs['statistics_panel_size']
    start_year = inputs['start_year']
    stippling_hatches = inputs['stippling_hatches']
    stippling_on = inputs['stippling_on']
    stippling_std_multiple = inputs['stippling_std_multiple']
    time_calculation = inputs['time_calculation']
    title = inputs['title'] 
    title_size = inputs['title_size']
    use_latex = inputs['use_latex']
    width = inputs['width'] 

    # We either have individual plots, in which case there could be multiple files arranged like [[file1, file2, file3, ...]],
    # or we could have ensemble plots, in which case there are at most two data sets, but potentially multiple files in each of the sets.
    # The preprocessing done earlier has formatted the files into a list of lists so that they are in the form:
    #   For individual plots: [[file1, file2, file3, ...]]    
    #   For ensemble plots: [["Control", "Full feedback"], ["Control_2", "Full feedback_2"], ...]
    num_files_in_each_set = len(netcdf_files)
    num_file_sets = len(netcdf_files[0])
    
    # Note that unlike time series, it is not meaningful to have a spatial plot with more than two ensembles, 
    # which is why we are restricting the total number of data sets to 2 when producing ensemble plots.
    # Note that num_file_sets and num_files_in_each_set represent different things depending on whether there
    #    are one or two ensembles.
    # If there is one ensemble (of any length, including 1 file):
    #    num_file_sets = the number of files in the ensemble
    #    num_files_in_each_set = 1 (the number of ensembles)
    #    only mean or sum can be calculated across the files and plotted (no comparison)
    # If there are two ensembles (or more, but the code exits if there are more than two):
    #    num_file_sets = 2 (the number of ensembles)
    #    num_files_in_each_set = the number of files in each ensemble
    #    the ensemble means are calculated and can be plotted separately or as one of two difference comparisons

    if (num_files_in_each_set == 1) or (num_file_sets == 2):
        pass
    else:
        error_message = "Error: For spatial plots, there can be only one or two sets of netcdf_files, " \
            + "and if there are two sets there must be at least two files in each set (for ensemble plots)."
        raise ValueError(error_message)

    # Read each of the NetCDF output files, which are arranged in a list of lists (2D matrix), into an xarray DataArray and then add each of these 
    # DataArrays to a single Pandas DataFrame that will store the data from all of the files. To form the DataArrays, calculate either the mean or sum 
    # between the start and end years for each lat/lon coordinate. We will later display some function of this mean or sum on the spatial plot.
    num_files_in_each_set = len(netcdf_files)
    num_file_sets = len(netcdf_files[0])
    df = pd.DataFrame()
    for file_set_index in range(num_file_sets):
        for file_index in range(num_files_in_each_set):
            file = netcdf_files[file_index][file_set_index]

            # Check that the requested start and end years fall within the year range of the file before slicing,
            # so the error message can report the available range rather than checking an already-filtered dataset.
            ds_full = xr.open_dataset(file)
            file_start = int(ds_full['year'].min())
            file_end = int(ds_full['year'].max())
            if start_year < file_start or start_year > file_end:
                ds_full.close()
                raise ValueError(f"Plot start_year {start_year} not in '{file}' (available: {file_start}–{file_end}); update input json file")
            if end_year < file_start or end_year > file_end:
                ds_full.close()
                raise ValueError(f"Plot end_year {end_year} not in '{file}' (available: {file_start}–{file_end}); update input json file")
            ds_full.close()

            if time_calculation == 'mean':
                da = xr.open_dataset(file).sel(year=slice(start_year, end_year)).mean(dim='year')[variable]*multiplier
            elif time_calculation == 'sum':
                da = xr.open_dataset(file).sel(year=slice(start_year, end_year)).sum(dim='year')[variable]*multiplier
                # If calculating the sum, change the per-time quantities and their units accordingly.
                per_time_labels = ['/year', '/month', '/day', '/hour', '/min', '/s']
                time_multipliers = np.array([1, years_TO_months, years_TO_days, years_TO_hours, years_TO_mins, years_TO_s])
                for index, per_time_label in enumerate(per_time_labels):
                    if per_time_label in title:
                        title = title.replace(per_time_label, '')
                        da.attrs['units'] = da.attrs['units'].replace(per_time_label, '')
                        da *= time_multipliers[index]
                        break
            # If there is more than one data set, modify the labels so that we know which data set corresponds to which column of the DataFrame.
            if num_file_sets >= 2:
                da = da.rename(f'{variable}_{file_index}_{file_set_index}')
            df = pd.concat([df, da.to_dataframe().dropna()], axis=1)
            da.close()
    
    # Initialize list that will store all the DataArrays that we will want to plot for the variable.
    dataArrays_to_plot = []

    # If we have a single data set (but potentially multiple files in this set), take the mean over all files
    # for each spatial coordinate. For a single file this is equivalent to just plotting that file's values directly.
    if num_files_in_each_set == 1:
        df = df.mean(axis=1)
        da = df.to_xarray()
        dataArrays_to_plot.append(da)
        print(f"Plotting single file set of {num_file_sets} file(s) for {variable}.")

    # If we have two data sets, we can plot either the absolute difference, percent difference, or the two data sets separately.
    elif num_file_sets == 2:
        # Take the mean over all files for each lat/lon coordinate in each data set.
        columns_control_set = [column for column in df.columns if column.endswith(f'_0')]
        columns_test_set = [column for column in df.columns if column.endswith(f'_1')]
        df_control_set = df[columns_control_set].mean(axis=1)
        df_test_set = df[columns_test_set].mean(axis=1)

        # If there is more than one file per data set (meaning that we have an ensemble) and we do not want separate plots,
        # we can compare the two data sets by performing a t-test at each individual lat/lon coordinate and later adding stippling.
        if num_files_in_each_set >= 2 and stippling_on and plot_type != 'separate_plots':
            # Perform this per-pixel t-test only if we do not want separate plots and if stippling_on is True (want to add p-value markers on plot).
            da_pvalues = df.apply(perform_ttest, columns_set_1=columns_control_set, columns_set_2=columns_test_set, axis=1).fillna(1).to_xarray()

        # Perform a t-test to compare the two spatial data sets as whole over all coordinates. Print the results to the console and to an output file.
        ttest = stats.ttest_ind(df_control_set, df_test_set)
        print_p_values(ttest, variable, p_value_threshold, p_value_file, plot_name, p_value_file_print_only_if_below_threshold)

        # Plot either an absolute difference, percent difference, or the two data sets separately.
        # If we have more than one file per data set (meaning that we have an ensemble), we will be examining the ensemble means in all cases.
        if plot_type == 'absolute_difference':
            # Plot absolute differences between the two data sets. 
            df = df_test_set - df_control_set
            da = df.to_xarray()
            dataArrays_to_plot.append(da)
            print(f"Calculating {plot_type} between {num_file_sets} ensembles with {num_files_in_each_set} files in each ensemble.")
        elif plot_type == 'percent_difference':
            # Plot percent differences between the two data sets. Add a tiny number to avoid a divide-by-zero error. Take the absolute value
            # so that if the control is negative, while the test set is positive, we get a positive value for the percent difference.
            df = ((df_test_set - df_control_set)/(df_control_set.abs() + EPSILON))*100
            #mask = df_test_set > 1.0e4*df_control_set
            #print('test', df_test_set[mask].head(10))
            #print('control', df_control_set[mask].head(10))
            da = df.to_xarray()
            dataArrays_to_plot.append(da)
            title = replace_inside_parentheses(title, rf'($\%$ difference)')
            print(f"Calculating {plot_type} between {num_file_sets} ensembles with {num_files_in_each_set} files in each ensemble.")
        elif plot_type == 'separate_plots':
            # Plot the two data sets individually in their own separate plots. 
            dataArrays_to_plot.append(df_control_set.to_xarray())
            dataArrays_to_plot.append(df_test_set.to_xarray())
            print(f"Calculating {plot_type} of {num_file_sets} ensembles with {num_files_in_each_set} files in each ensemble.")

    # this is when there are more than three files listed in one input set, so differences don't make sense
    else:
        error_message = (f"Error: variable '{variable}', num_files_in_each_set={num_files_in_each_set}, "
                         f"num_file_sets={num_file_sets}, plot_type='{plot_type}': "
                         "If there are more than two netcdf_files in one single input set "
                         "plot_type must be mean or sum.")
        raise ValueError(error_message)

    # Iterate over all dataArrays in the list to create a plot for each one.
    for da_index, da in enumerate(dataArrays_to_plot):

        # Calculate some basic statistics of the current dataArray.
        min, mean, median, max, std = calculate_statistics_of_xarray(da)

        # Use LaTeX for the labels if specified to do so.
        if use_latex:
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif', weight='bold') 
        
        # If plotting a percent difference and no colorbar limits are specified, set them to be -100 and 100 percent if the max exceeds 100 percent.
        if plot_type == 'percent_difference' and not cbar_limits and max > 100:
            cbar_limits = [-100, 100]

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
                cbar.mappable.set_clim(cbar_limits[0], cbar_limits[1])

        # Add stippling to indicate regions of potential statistical significance. 
        if stippling_on:
            plt.rcParams['hatch.linewidth'] = 0.5
            plt.rcParams['hatch.color'] = 'gray'
            if num_file_sets == 2 and num_files_in_each_set >= 2 and plot_type != 'separate_plots':
                # If there are two data sets and at least two files in each data set so that we will have previously calculated p-values at each
                # lat/lon coordinate, the stippling will indicate regions where the p-value is less than the designated threshold.
                mask = da_pvalues <= p_value_threshold
                ax.contourf(da_pvalues.lon, da_pvalues.lat, mask, levels=1, hatches=['', stippling_hatches], alpha=0, transform=ccrs.PlateCarree())
            else:
                # For all other cases, the stippling will indicate regions where the value is +/- some multiple of the standard deviation 
                # (default is 2*std) away from the mean.
                mask = np.abs(da) >= mean + stippling_std_multiple*std
                ax.contourf(da.lon, da.lat, mask, levels=1, hatches=['', stippling_hatches], alpha=0, transform=ccrs.PlateCarree())
        
        # Display statistics.
        ax.text(x=0.88, y=0.9, s=f'Max:{max:.2e}\nMean:{mean:.2e}\nMedian:{median:.2e}\nMin:{min:.2e}', ha='left', 
                fontsize=statistics_panel_size, transform=ax.transAxes)
        ax.text(x=0.88, y=0.05, s=f'Std:{std:.2e}', ha='left', fontsize=statistics_panel_size, transform=ax.transAxes)

        # Add additional features like coastline and oceans.
        ax.coastlines(lw=0.6)

        # Save the figure and then close it. Record the elapsed time.
        end_time = time.time()
        elapsed_time = end_time - start_time
        if len(dataArrays_to_plot) > 1:
            # If we want to plot two data sets separately, add a set number to the figure name for each set.
            save_figure(plot_name + f'_set_{da_index+1}', fig, inputs)
            print(f"Elapsed time for producing plots for {variable} (set {da_index+1}) in {plot_directory}: {elapsed_time:.2f} seconds") 
        else:
            save_figure(plot_name, fig, inputs)
            print(f"Elapsed time for producing plots for {variable} in {plot_directory}: {elapsed_time:.2f} seconds") 
        plt.close(fig)

def plot_spatial_data_from_netcdf_files(inputs):
    """ 
    Parses a dictionary of inputs (keys are options, values are choices for those options) to create spatial plots for a single variable from either 
    E3SM atmosphere (EAM) or land (ELM) outputs. The data for these spatial plots are stored in NetCDF files specified by the inputs dictionary.

    Parameters:
        input: Dictionary containing the user plotting choice inputs for different options. This dictionary is assumed to be complete (pre-processed).

    Returns:
        N/A.
    """
    # EAM requires a separate grid file since the simulations are run on an unstructured grid.
    grid_file = inputs['grid_file']
    if grid_file:
        plot_spatial_data_eam(inputs, grid_file)
    else:
        plot_spatial_data_elm(inputs)

###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file(s) on the command line.
    start_time = time.time()
    if len(sys.argv) < 2:
        print('Usage: python e3sm_plot_spatial_data.py `path/to/json/input/file(s)\'')
        sys.exit()

    # Read and load the JSON file(s) into a list of dictionaries.
    inputs = []
    for index in range(1, len(sys.argv)):
        input_file = sys.argv[index]
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

    # Create all of the spatial plots in parallel.
    try:
        with multiprocessing.Pool(processes=MAX_PROCESSES) as pool:
            pool.map(plot_spatial_data_from_netcdf_files, list_of_inputs)
    except ValueError as e:
        print(e)
        sys.exit(1)

    # Sort all the p-value files alphabetically.
    for inputs in list_of_inputs:
        file = inputs['p_value_file']
        if os.path.exists(file): 
            sort_file(file)

    # Print the total execution time to produce all the plots.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all the plots: {elapsed_time:.2f} seconds")
