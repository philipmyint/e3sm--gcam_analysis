'''
Python modules for making spatial plots of EAM outputs
'''
import os
import numpy as np
import uxarray as ux
import holoviews as hv
hv.extension('matplotlib')

import matplotlib
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

import cartopy.crs as ccrs
import matplotlib.colors as colors

__author__ = 'Eva Sinha'
__email__  = 'eva.sinha@pnnl.gov'

from util_myDict import *
from util_spatial_plots import *

plt.rc('figure', titlesize=20)
plt.rc('legend', fontsize=20)
plt.rc('axes',   labelsize=20, titlesize=20)
plt.rc('xtick',  labelsize=20)
plt.rc('ytick',  labelsize=20)
plt.rc('figure', figsize=(11, 8.5))

# -----------------------------------------------------------

fpath     = '/lcrc/group/e3sm/data/inputdata/iac/giac/atm/cam/ggas/'
fname     = 'ne30pg2_CSEM_hist_ssp370_ocean_flux_1850-2100_c20231108.nc'

#fpath     = '/lcrc/group/acme/public_html/e3sm_support/compset_generation/ssp245/ocn/co2_flux/merged_regridded_data/'
#fname     = 'fgco2_CESM2_SSP245_ne30pg2_2015-2100.nc'
var       = 'fgco2'
yr_start  = 2071
yr_end    = 2100
time_period  = 'Annual'

grid_path = '/lcrc/group/e3sm/data/inputdata/share/meshes/homme/ne30pg2_scrip_c20191218.nc'

uxds = ux.open_mfdataset(grid_path, fpath+fname)
uxds_yr = uxds.time.dt.year.values
sel_index = list(np.where((uxds_yr >= yr_start) * (uxds_yr <=yr_end)))
# Subset for select years and scale by conversion factor for time
uxds = uxds[var]
uxds = uxds.isel(time=sel_index[0])
print(uxds)
uxds = uxds/conv_factor[time_period][var]

out_fname = 'ocean_CO2_flux_SSP370.png'
title = myDict_units[var] + '\nSSP370 - Average from '+ str(yr_start) + ' to ' + str(yr_end)
#out_fname = 'ocean_CO2_flux_SSP245.png'
#title = myDict_units[var] + '\nSSP245 - Average from '+ str(yr_start) + ' to ' + str(yr_end)

uxds_plot_global_polycollection(uxds, title=title, fig_wt=8*2, fig_ht=14, 
                                  fname=out_fname + '_' + str(yr_start) + '_' + str(yr_end) + '.png')
