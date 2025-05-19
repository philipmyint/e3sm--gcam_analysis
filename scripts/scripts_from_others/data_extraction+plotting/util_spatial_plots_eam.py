"""
Python modules for plotting time series
"""

import os
import matplotlib as mpl
mpl.use('Agg')
import pandas as pd
import geopandas
import xarray as xr
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import uxarray as ux
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
def uxds_plot_global_polycollection(uxds, title, levels, fig_wt, fig_ht, fname):

    # Change directory    
    os.chdir('../../figures/e3sm_figures/')

    # Estimate temporal mean
    uxds = uxds.mean(dim='time')
    print(uxds)

    pc = uxds.to_polycollection(cache=False)

    print(pc)
    norm = colors.CenteredNorm(vcenter = 0)

    pc.set_transform(ccrs.PlateCarree())
    plt.figure(figsize=(fig_wt, fig_ht))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.coastlines()
    pc.set_cmap('bwr')
    pc.set_norm(norm)
    ax.add_collection(pc, levels=levels)
    ax.set_title(title)
    ax.set_global()
    plt.colorbar(pc, shrink= 0.5, pad=0.02),

    plt.savefig(fname, bbox_inches='tight')

    plt.close(fig=None)

    # Change directory    
    os.chdir('../../workflow/e3sm_analysis/')

# -----------------------------------------------------------
