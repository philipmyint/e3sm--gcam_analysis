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
from utility_functions import check_is_list_of_lists, print_p_values, replace_inside_parentheses, sort_file
from utility_plots import *
from utility_xarray import calculate_statistics_of_xarray, convert_xarray_to_uxarray

""" Dictionary of default input values for spatial plots. """
default_inputs_spatial_data = {'plot_directory': './',
                    'calculation_type': 'mean',
                    'plot_type_for_2_sets': 'absolute_difference',
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
                    'cbar_x_offset': 0.06,
                    'title_size': axis_label_size_default, 
                    'use_latex': False,
                    'produce_png': False,
                    'bbox_inches': 'tight',
                    'statistics_panel_size': legend_label_size_default,
                    'p_value_threshold': 0.05,
                    'p_value_file': None,
                    'p_value_file_print_only_if_below_threshold': True,
                    'stippling_on': False,
                    'stippling_std_multiple': 2,
                    'stippling_hatches': 'xxxx',
                    'grid_file': None
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
    # Turn the list of NetCDF files and their corresponding labels into a list of lists, if they are not already in that form.
    for input_type in ['netcdf_files']:
        if isinstance(inputs[input_type], str):
            inputs[input_type] = [[inputs[input_type]]]
        elif isinstance(inputs[input_type], list) and not check_is_list_of_lists(inputs[input_type]):
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
        # Default for the title of a variable is to use the column header for that variable from the Dataset.
        if not any(key == variable for key in inputs['title'].keys()):
            # Replace ^2 with $^2$ in the units for the variable, so that the intended exponent (e.g., m^2) gets rendered correctly.
            units = ds[variable].attrs['units'].replace('^2', '$^2$')
            # Replace /m2 with /m$^2$ in the units for the variable, so that it gets rendered correctly.
            units = units.replace('/m2', '/m$^2$')
            inputs['title'][variable] = rf'{variable} ({units})'
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

def plot_spatial_data_eam(inputs, grid_file):
    """ 
    Parses a dictionary of inputs (keys are options, values are choices for those options) to create spatial plots for a single variable from E3SM EAM
    outputs. The data for these spatial plots are stored in NetCDF files specified by the inputs dictionary.

    Parameters:
        input: Dictionary containing the user plotting choice inputs for different options. This dictionary is assumed to be complete (pre-processed).
        grid_file: Path and name of the grid file for the EAM unstructured mesh.

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
    statistics_panel_size = inputs['statistics_panel_size']
    calculation_type = inputs['calculation_type']
    plot_type_for_2_sets = inputs['plot_type_for_2_sets']
    use_latex = inputs['use_latex']
    start_year = inputs['start_year']
    end_year = inputs['end_year']
    p_value_threshold = inputs['p_value_threshold']
    p_value_file = inputs['p_value_file']
    p_value_file_print_only_if_below_threshold = inputs['p_value_file_print_only_if_below_threshold']
    stippling_on = inputs['stippling_on']
    stippling_std_multiple = inputs['stippling_std_multiple']
    stippling_hatches = inputs['stippling_hatches']

    # Store the grid file in an uxarray Dataset.
    grid = ux.open_grid(grid_file)

    # Read each of the NetCDF output files, which are arranged in a list of lists (2D matrix), into an uxarray DataArray and then add each of these 
    # DataArrays to a single uxarray Dataset that will store the data from all of the files. To form the DataArrays, calculate either the mean or sum 
    # between the start and end years for each lat/lon coordinate. We will later display some function of this mean or sum on the spatial plot.
    num_files_in_each_set = len(netcdf_files)
    num_file_sets = len(netcdf_files[0])
    uxds = ux.UxDataset()
    for file_set_index in range(num_file_sets):
        for file_index in range(num_files_in_each_set):
            file = netcdf_files[file_index][file_set_index]
            # Create a temporary NetCDF file with data between only the start and end years. 
            ds = xr.open_dataset(file).sel(year=slice(start_year, end_year))[variable]
            file = os.path.join(plot_directory, f'temp_{variable}.nc')
            ds.to_netcdf(file, 'w')
            if calculation_type == 'mean':
                uxda = ux.open_dataset(grid_file, file).mean(dim='year')[variable]*multiplier
            elif calculation_type == 'sum':
                uxda = ux.open_dataset(grid_file, file).sum(dim='year')[variable]*multiplier
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
            if num_file_sets == 2:
                uxda = uxda.rename(f'{variable}_{file_index}_{file_set_index}')
            uxds[uxda.name] = uxda
            # Delete the temporary NetCDF file now that the data have been read.
            os.system(f'rm {file}')
    
    # Initialize list that will store all the uxarray DataArrays that we will want to plot for the variable.
    uxDataArrays_to_plot = []
    df = uxds.to_dataframe()
    if num_file_sets == 2:
        columns_control_set = [column for column in df.columns if column.endswith(f'_0')]
        columns_test_set = [column for column in df.columns if column.endswith(f'_1')]
        # Perform a t-test to compare the two spatial data sets as whole over all coordinates. Print the results to the console and to an output file.
        df_control_set = df[columns_control_set].mean(axis=1)
        df_test_set = df[columns_test_set].mean(axis=1)
        ttest = stats.ttest_ind(df_control_set, df_test_set)
        print_p_values(ttest, variable, p_value_threshold, p_value_file, plot_directory, p_value_file_print_only_if_below_threshold)
        if plot_type_for_2_sets == 'absolute_difference':
            # Plot absolute differences between the two data sets. 
            df = df_test_set - df_control_set
            uxDataArrays_to_plot.append(convert_xarray_to_uxarray(df.to_xarray(), grid, variable=variable))
        elif plot_type_for_2_sets == 'percent_difference':
            # Plot percent differences between the two data sets. Add a tiny number to avoid a divide-by-zero error. Take the absolute value
            # so that if the control is negative, while the test set is positive, we get a positive value for the percent difference.
            df = ((df_test_set - df_control_set)/(df_control_set.abs() + EPSILON))*100
            uxDataArrays_to_plot.append(convert_xarray_to_uxarray(df.to_xarray(), grid, variable=variable))
            title = replace_inside_parentheses(title, rf'($\%$ difference)')
        elif plot_type_for_2_sets == 'separate_plots':
            # Plot the two data sets individually in their own separate plots. 
            uxDataArrays_to_plot.append(convert_xarray_to_uxarray(df_control_set.to_xarray(), grid, variable=variable))
            uxDataArrays_to_plot.append(convert_xarray_to_uxarray(df_test_set.to_xarray(), grid, variable=variable))
    elif num_file_sets == 1:
        # If we have a single data set (but potentially multiple files in this data set), take the mean over all files for each lat/lon coordinate.
        df = df.mean(axis=1)
        uxDataArrays_to_plot.append(convert_xarray_to_uxarray(df.to_xarray(), grid, variable=variable))

    # If stippling_on is True, meaning we want to add markers on plot to indicate potential regions of statistical significance, 
    if stippling_on:
        uxds['lat'] = grid.face_lat
        uxds['lon'] = grid.face_lon
        # If there is more than one file per data set and we do not want separate plots, perform a t-test at each individual lat/lon coordinate.
        if num_files_in_each_set >= 2 and plot_type_for_2_sets != 'separate_plots':
            df = uxds.to_dataframe()
            df = df.groupby(['lat', 'lon']).mean()
            da_pvalues = df.apply(perform_ttest, columns_set_1=columns_control_set, columns_set_2=columns_test_set, axis=1).fillna(1).to_xarray().fillna(1)
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
        if plot_type_for_2_sets == 'percent_difference' and not cbar_limits and max > 100:
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

        # Add stippling to indicate regions of potential statistical significance. 
        if stippling_on:
            plt.rcParams['hatch.linewidth'] = 0.5
            plt.rcParams['hatch.color'] = 'gray'
            if num_file_sets == 2 and num_files_in_each_set >= 2 and plot_type_for_2_sets != 'separate_plots':
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
    Parses a dictionary of inputs (keys are options, values are choices for those options) to create spatial plots for a single variable from E3SM ELM
    outputs. The data for these spatial plots are stored in NetCDF files specified by the inputs dictionary.

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
    statistics_panel_size = inputs['statistics_panel_size']
    calculation_type = inputs['calculation_type']
    plot_type_for_2_sets = inputs['plot_type_for_2_sets']
    use_latex = inputs['use_latex']
    start_year = inputs['start_year']
    end_year = inputs['end_year']
    p_value_threshold = inputs['p_value_threshold']
    p_value_file = inputs['p_value_file']
    p_value_file_print_only_if_below_threshold = inputs['p_value_file_print_only_if_below_threshold']
    stippling_on = inputs['stippling_on']
    stippling_std_multiple = inputs['stippling_std_multiple']
    stippling_hatches = inputs['stippling_hatches']

    # Read each of the NetCDF output files, which are arranged in a list of lists (2D matrix), into an xarray DataArray and then add each of these 
    # DataArrays to a single Pandas DataFrame that will store the data from all of the files. To form the DataArrays, calculate either the mean or sum 
    # between the start and end years for each lat/lon coordinate. We will later display some function of this mean or sum on the spatial plot.
    num_files_in_each_set = len(netcdf_files)
    num_file_sets = len(netcdf_files[0])
    df = pd.DataFrame()
    for file_set_index in range(num_file_sets):
        for file_index in range(num_files_in_each_set):
            file = netcdf_files[file_index][file_set_index]
            if calculation_type == 'mean':
                da = xr.open_dataset(file).sel(year=slice(start_year, end_year)).mean(dim='year')[variable]*multiplier
            elif calculation_type == 'sum':
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
            if num_file_sets == 2:
                da = da.rename(f'{variable}_{file_index}_{file_set_index}')
            df = pd.concat([df, da.to_dataframe().dropna()], axis=1)

    # Initialize list that will store all the DataArrays that we will want to plot for the variable.
    dataArrays_to_plot = []
    if num_file_sets == 2:
        columns_control_set = [column for column in df.columns if column.endswith(f'_0')]
        columns_test_set = [column for column in df.columns if column.endswith(f'_1')]
        # If there is more than one file per data set, we can compare the two data sets by performing a t-test at each individual lat/lon coordinate.
        if num_files_in_each_set >= 2 and stippling_on and plot_type_for_2_sets != 'separate_plots':
            # Perform this per-pixel t-test only if we do not want separate plots and if stippling_on is True (want to add p-value markers on plot).
            da_pvalues = df.apply(perform_ttest, columns_set_1=columns_control_set, columns_set_2=columns_test_set, axis=1).fillna(1).to_xarray()
        # Perform a t-test to compare the two spatial data sets as whole over all coordinates. Print the results to the console and to an output file.
        df_control_set = df[columns_control_set].mean(axis=1)
        df_test_set = df[columns_test_set].mean(axis=1)
        ttest = stats.ttest_ind(df_control_set, df_test_set)
        print_p_values(ttest, variable, p_value_threshold, p_value_file, plot_directory, p_value_file_print_only_if_below_threshold)
        if plot_type_for_2_sets == 'absolute_difference':
            # Plot absolute differences between the two data sets. 
            df = df_test_set - df_control_set
            da = df.to_xarray()
            dataArrays_to_plot.append(da)
        elif plot_type_for_2_sets == 'percent_difference':
            # Plot percent differences between the two data sets. Add a tiny number to avoid a divide-by-zero error. Take the absolute value
            # so that if the control is negative, while the test set is positive, we get a positive value for the percent difference.
            df = ((df_test_set - df_control_set)/(df_control_set.abs() + EPSILON))*100
            #mask = df_test_set > 1.0e4*df_control_set
            #print('test', df_test_set[mask].head(10))
            #print('control', df_control_set[mask].head(10))
            da = df.to_xarray()
            dataArrays_to_plot.append(da)
            title = replace_inside_parentheses(title, rf'($\%$ difference)')
        elif plot_type_for_2_sets == 'separate_plots':
            # Plot the two data sets individually in their own separate plots. 
            dataArrays_to_plot.append(df_control_set.to_xarray())
            dataArrays_to_plot.append(df_test_set.to_xarray())
    elif num_file_sets == 1:
        # If we have a single data set (but potentially multiple files in this data set), take the mean over all files for each lat/lon coordinate.
        df = df.mean(axis=1)
        da = df.to_xarray()
        dataArrays_to_plot.append(da)

    # Iterate over all dataArrays in the list to create a plot for each one.
    for da_index, da in enumerate(dataArrays_to_plot):

        # Calculate some basic statistics of the current dataArray.
        min, mean, median, max, std = calculate_statistics_of_xarray(da)

        # Use LaTeX for the labels if specified to do so.
        if use_latex:
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif', weight='bold') 
        
        # If plotting a percent difference and no colorbar limits are specified, set them to be -100 and 100 percent if the max exceeds 100 percent.
        if plot_type_for_2_sets == 'percent_difference' and not cbar_limits and max > 100:
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
            if num_file_sets == 2 and num_files_in_each_set >= 2 and plot_type_for_2_sets != 'separate_plots':
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
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(plot_spatial_data_from_netcdf_files, list_of_inputs)
    
    # Sort all the p-value files alphabetically.
    for inputs in list_of_inputs:
        file = inputs['p_value_file']
        if os.path.exists(file): 
            sort_file(file)

    # Print the total execution time to produce all the plots.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all plots: {elapsed_time:.2f} seconds")