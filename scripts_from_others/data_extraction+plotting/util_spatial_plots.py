"""
Python modules for plotting time series
"""

import os
import matplotlib as mpl
mpl.use('Agg')
import pandas as pd
import xarray as xr
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib.colors as colors

__author__ = 'Eva Sinha'
__email__  = 'eva.sinha@pnnl.gov'

from util_myDict import *

plt.rc('figure', titlesize=20)
plt.rc('legend', fontsize=20)
plt.rc('axes',   labelsize=20, titlesize=20)
plt.rc('xtick',  labelsize=20)
plt.rc('ytick',  labelsize=20)
plt.rc('figure', figsize=(11, 8.5))

# Colormap for sptial plots
colormap = 'WhiteBlueGreenYellowRed.rgb'
rgb_arr = np.loadtxt(colormap)
rgb_arr = rgb_arr / 255.0
cmap = LinearSegmentedColormap.from_list(name=colormap, colors=rgb_arr)

diff_cmap = 'BrBG'

# -----------------------------------------------------------
def mean(variable, axis='xy'):
    return cdutil.averager(variable, axis=axis, weights='generate')

# -----------------------------------------------------------
# Source: https://github.com/E3SM-Project/e3sm_diags/blob/master/e3sm_diags/driver/lat_lon_driver.py
def create_metrics(da):
    """Creates the mean, max, min, rmse, corr in a dictionary"""
    # For input None, metrics are instantiated to 999.999.
    # Apply float() to make sure the elements in metrics_dict are JSON serializable, i.e. np.float64 type is JSON serializable, but not np.float32.
    missing_value = 999.999
    metrics_dict = {}

    metrics_dict = {
        'min': float(da.min()) if da is not None else missing_value,
        'max': float(da.max()) if da is not None else missing_value,
        #'mean': float(mean(da)) if da is not None else missing_value,# Was giving error - You have specified an invalid axis
        'mean': float(da.mean()) if da is not None else missing_value,
    }

    return(metrics_dict)

# -----------------------------------------------------------
def xr_plot_global(da_plot, stats, cmap_col, title, fig_wt, fig_ht, levels, norm, fname, stipple_data=None):

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
    else:
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
       cbar.ax.tick_params(labelsize=15, pad=pad, length=0)

    # Min, Mean, Max
    ax.text(x=0.88, y=0.9, s='Max\nMean\nMin', ha='left', fontsize=15, transform=ax.transAxes)
    ax.text(x=1.0, y=0.9, s='%.2f\n%.2f\n%.2f'%stats[0:3], ha='right', fontsize=15, transform=ax.transAxes)

    # Add additional features like coastline and oceans
    ax.coastlines(lw = 0.6)

    plt.savefig(fname, bbox_inches='tight')

    plt.close(fig=None)

# -----------------------------------------------------------
