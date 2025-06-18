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

plt.rc('figure', titlesize=20)
plt.rc('legend', fontsize=20)
plt.rc('axes',   labelsize=20, titlesize=20)
plt.rc('xtick',  labelsize=20)
plt.rc('ytick',  labelsize=20)
plt.rc('figure', figsize=(11, 8.5))

myDict_clims = {'SFCO2': {'min_clim': -10, 'max_clim':10},
                'SFCO2_FFF': {'min_clim': -10, 'max_clim':10},
                'SFCO2_LND': {'min_clim': -2, 'max_clim':2},
                'SFCO2_OCN': {'min_clim': -0.3, 'max_clim':0.3},
                'CO2': {'min_clim': 400, 'max_clim':800}}

# -----------------------------------------------------------
def uxds_plot_global_polycollection(uxds, var, title, cmap, fig_wt, fig_ht, fname):

    # Change directory    
    os.chdir('../../figures/e3sm_figures/')

    uxds = uxds.mean(dim='time')
    print(uxds)
    print(np.min(uxds.values), np.max(uxds.values))

    if(var == 'CO2'):
      #min_clim = np.min(uxds.values)
      #max_clim = np.max(uxds.values)
      norm = None
    else:
      norm = colors.CenteredNorm(vcenter = 0)

    min_clim = myDict_clims[var]['min_clim']
    max_clim = myDict_clims[var]['max_clim']

    pc = uxds.to_polycollection(cache=False)

    pc.set_transform(ccrs.PlateCarree())
    plt.figure(figsize=(fig_wt, fig_ht))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.coastlines()
    pc.set_cmap(cmap)
    pc.set_norm(norm)
    pc.set_clim(min_clim, max_clim)
    ax.add_collection(pc)
    ax.set_title(title)
    ax.set_global()
    plt.colorbar(pc, shrink= 0.5, pad=0.02, extend='both'),

    plt.savefig(fname, bbox_inches='tight')

    plt.close(fig=None)

    # Change directory    
    os.chdir('../../workflow/e3sm_analysis/')

# -----------------------------------------------------------
# List of variable names that we want to keep
varnames = ['SFCO2', 'SFCO2_FFF', 'SFCO2_LND', 'SFCO2_OCN']
varnames = ['CO2']

#fpath_scratch = '/lcrc/group/e3sm/ac.dalei.hao/E3SM_GCAM_simulations/'
#run   = 'SSP245_FULL_FDBK'

fpath_scratch = '/lcrc/group/e3sm/ac.eva.sinha/'
run    = 'Historical'
#run   = 'SSP370_FULL_FDBK'

caseid        = myDict_caseid[run]
e3sm_comp   = 'eam'
yr_start    = 1995
yr_end      = 2014
#yr_start    = 2071
#yr_end      = 2090
time_period = 'Annual'
cmap = 'bwr'

#fpath = fpath_scratch + caseid + '/run/'
fpath = fpath_scratch + caseid + '/archive/atm/hist/'
fnames = []
for yr in range(int(yr_start), int(yr_end)+1, yr_step):
  for mon in range(1,13):
    fnames.append(fpath + '/' + caseid + '.' + e3sm_comp + '.h0.' + str(yr).zfill(4) + '-' + str(mon).zfill(2) + '.nc')

grid_path = '/lcrc/group/e3sm/data/inputdata/share/meshes/homme/ne30pg2_scrip_c20191218.nc'

uxds = ux.open_mfdataset(grid_path, fnames)

for ind, var in enumerate(varnames):

  # Scale by conversion factor for time
  uxds[var] = uxds[var]/conv_factor[time_period][var]

  if( var == 'CO2'):
    cmap = 'jet'

    # Use concentration at last level
    uxds[var] = uxds[var].isel(lev = uxds.dims['lev'] - 1)

    # Convert CO_2 from units of kg/kg to ppm
    uxds[var] = uxds[var] * ppmfac
    print(np.min(uxds[var].values), np.max(uxds[var].values))

  out_fname = run + '_' + var + '_emissions'
  title=run + '; ' + myDict_units[var] + '\nAverage from '+ str(yr_start) + ' to ' + str(yr_end)

  uxds_plot_global_polycollection(uxds[var], var, title=title, cmap=cmap,
                                  fig_wt=8*2, fig_ht=14, 
                                  fname=out_fname + '_' + str(yr_start) + '_' + str(yr_end) + '.png')

