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
from utility_plots import default_inputs_time_series, set_figure_options
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs

def xr_plot_global(data):

    da_plot = data['da']
    cmap_col = data['cmap_col']
    title = data['title'] 
    fig_wt = data['width'] 
    fig_ht = data['height']
    levels = data['levels'] 
    fname = data['fname'] 
    stipple_data = None

    fig = plt.figure(figsize=(fig_wt, fig_ht))
    ax  = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection = ccrs.Robinson(central_longitude=0))

    cmap = plt.get_cmap(cmap_col)
    if(cmap_col == 'bwr'):
        cmap.set_extremes(under='darkblue', over='darkred')
    else:
        cmap.set_extremes(under='white', over='darkred')

    # spatial plot using xarray
    #fg = da_plot.plot.contourf(ax          = ax,
    fg = da_plot.plot(ax          = ax,
                      transform   = ccrs.PlateCarree(), # coordinate system of data
                      levels      = levels, 
                      #norm        = norm,
                      cmap        = cmap,
                      extend      = 'both',
                      add_colorbar= False)

    # Add stippling
    if stipple_data is not None:
        mask = stipple_data <= 0.05
        tmp = mask.values
        tmp = tmp[tmp == True]
        print(tmp.size)
        ax.contourf(stipple_data.lon, stipple_data.lat, mask,  1, hatches=['', 'xx'], alpha=0,
                    transform=ccrs.PlateCarree())

    # Add title
    ax.set_title(title)

    # Add colorbar
    cax = fig.add_axes([ax.get_position().x1+0.01, ax.get_position().y0, 0.02, ax.get_position().height])
    cbar = plt.colorbar(fg, cax=cax)

    if levels is None:
        cbar.ax.tick_params(labelsize=15, length=0)
    '''else:
        maxval = np.amax(np.absolute(levels[1:-1]))
        if maxval < 0.2:
            fmt = "%5.3f"
            pad = 28
        elif maxval < 10.0:
            fmt = "%5.2f"
            pad = 40
        elif maxval < 100.0:
            fmt = "%5.0f"
            pad = 25
        elif maxval < 10000.0:
            fmt = "%5.0f"
            pad = 30
        elif maxval > 9999.0:
            fmt = "%.0f"
            pad = 50
        else:
            fmt = "%6.1f"
            pad = 40

        cbar.set_ticks(levels[1:-1])
        labels = [fmt % level for level in levels[1:-1]]
        cbar.ax.set_yticklabels(labels, ha='right')
        cbar.ax.tick_params(labelsize=15, pad=pad, length=0)'''

    # Min, Mean, Max
    ax.text(x=0.88, y=0.9, s='Max\nMean\nMin', ha='left', fontsize=15, transform=ax.transAxes)
    #ax.text(x=1.0, y=0.9, s='%.2f\n%.2f\n%.2f'%stats[0:3], ha='right', fontsize=15, transform=ax.transAxes)

    # Add additional features like coastline and oceans
    ax.coastlines(lw = 0.6)

    plt.savefig(fname, bbox_inches='tight')

    plt.close(fig=None)         

###---------------Begin execution---------------###
if __name__ == '__main__':

    list = []
    variables = ['NPP', 'TOTECOSYSC']
    start_time = time.time()
    
    options = []
    for variable in variables:
        files = ['./../2025_DiVittorio_et_al/control_spatial_data_elm.nc', './../2025_DiVittorio_et_al/full_feedback_spatial_data_elm.nc']
        for index, file in enumerate(files):
            ds = xr.open_dataset(files[index])[variable].mean(dim='year')
            list.append(ds)
        da = list[0] - list[1]
        cmap_col = 'bwr'
        data = {'da': da, 'cmap_col': cmap_col, 'title': 'title', 'width': 10, 'height': 8, 'fname':'plot_spatial', 'levels':None}
        options.append(data)
    
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(xr_plot_global, options)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for: {elapsed_time:.2f} seconds")