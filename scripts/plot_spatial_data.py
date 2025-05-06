import cartopy.crs as ccrs
import json
from matplotlib import pyplot as plt
import multiprocessing
import numpy as np
import pandas as pd
import os
import sys
import time
import xarray as xr
from utility_constants import MONTH_NUM_TO_NAME
from utility_dataframes import get_columns_without_units_in_dataframe, get_matching_column_in_dataframe
from utility_plots import *

""" Dictionary of default input values for spatial plots. """
default_inputs_spatial_data = {'plot_directories': './',
                    'calculation_types': 'mean',
                    'multipliers': 1,
                    'start_years': 2071,
                    'end_years': 2090,
                    'projections': ccrs.Robinson,
                    'cbars_on': True,
                    'cmap': 'bwr',
                    'widths': width_default,
                    'heights': height_default,
                    'x_scales': scale_default,
                    'y_scales': scale_default,
                    'x_limits': axis_limits_default,
                    'y_limits': axis_limits_default,
                    'x_tick_label_sizes': tick_label_size_default,
                    'y_tick_label_sizes': tick_label_size_default,
                    'x_label_sizes': axis_label_size_default,
                    'y_label_sizes': axis_label_size_default,
                    'legend_label_sizes': legend_label_size_default,
                    'legends_on': legend_on_default,
                    'linewidths': linewidth_default,
                    'plot_colors': plot_colors_default,
                    'linestyle_tuples': linestyle_tuples_default,
                    'use_latex': False
}

def xr_plot_global(inputs):

    start_time = time.time()
    da = inputs['da']
    cmap_col = inputs['cmap_col']
    projection = inputs['projections']
    title = inputs['title'] 
    fig_wt = inputs['widths'] 
    fig_ht = inputs['heights']
    levels = inputs['levels'] 
    fname = inputs['fname']
    cbar_on = inputs['cbar_on'] 

    fig = plt.figure(figsize=(fig_wt, fig_ht))
    ax  = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=projection())

    cmap = plt.get_cmap(cmap_col)

    fg = da.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap, extend='both', add_colorbar= False)

    # Add title
    ax.set_title(title)

    # Add colorbar
    cax = fig.add_axes([ax.get_position().x1+0.01, ax.get_position().y0, 0.02, ax.get_position().height])
    #cbar = plt.colorbar(fg, cax=cax)

    #if levels is None:
    #    cbar.ax.tick_params(labelsize=15, length=0)

    # Min, Mean, Max
    ax.text(x=0.88, y=0.9, s='Max\nMean\nMin', ha='left', fontsize=15, transform=ax.transAxes)
    #ax.text(x=1.0, y=0.9, s='%.2f\n%.2f\n%.2f'%stats[0:3], ha='right', fontsize=15, transform=ax.transAxes)

    # Add additional features like coastline and oceans
    ax.coastlines(lw = 0.6)

    plt.savefig(fname, bbox_inches='tight')

    plt.close(fig=None)         

###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file on the command line.
    if len(sys.argv) != 2:
        print('Usage: plot_spatial_data.py `path/to/json/input/file\'')
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

    list2 = []
    variables = ['TOTECOSYSC']
    start_time = time.time()
    
    '''options = []
    for variable in variables:
        files = ['./../2025_DiVittorio_et_al/control_spatial_data_elm.nc', './../2025_DiVittorio_et_al/full_feedback_spatial_data_elm.nc']
        for index, file in enumerate(files):
            ds = xr.open_dataset(files[index]).mean(dim='year')[variable]
            list2.append(ds)
        da = list2[0] - list2[1]
        cmap_col = 'bwr'
        data = {'da': da, 'cmap_col': cmap_col, 'title': 'title', 'width': 10, 'height': 8, 'fname':'plot_spatial', 'levels':None, 'projection':ccrs.Robinson}
        options.append(data)
    
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(xr_plot_global, options)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for: {elapsed_time:.2f} seconds")'''