'''
Python modules for making spatial plots of ELM outputs
'''
import os
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from scipy import stats
#import cdutil
import xarray as xr
from datetime import timedelta

__author__ = 'Eva Sinha'
__email__  = 'eva.sinha@pnnl.gov'

from util_spatial_plots import *
from util_read_data import *
from util_myDict import *
from util_contour_levels import *

#----------------------------------------------------------
# Calculate the T-test for the means of two independent samples of scores
def apply_ttest_ind(df, var, sel_set):

   a = df[df['Set'] == 'CONTROL'][var]
   b = df[df['Set'] == sel_set][var]

   ttest = stats.ttest_ind(a,b)

   return(ttest.pvalue)

# -----------------------------------------------------------

# List of variable names that we want to keep
varnames = ['TBOT','NPP','NBP','TOTECOSYSC','TOTVEGC', 'LAND_USE_FLUX', 'WOOD_HARVESTC']
varnames = ['TOTECOSYSC','TOTVEGC', 'LAND_USE_FLUX', 'WOOD_HARVESTC']
varnames = ['TOTECOSYSC','TOTVEGC','TBOT']
fpath_scratch = '/lcrc/group/e3sm/ac.eva.sinha/'
runs = ['FULL_FDBK', 'CONTROL']
#runs = ['AG_FDBK', 'CONTROL']
#runs = ['C_FDBK', 'CONTROL']

conv_unit   = 'gCm2_gCm2'
time_period = 'Daily'
multi_area  = False
#conv_unit   = 'gC_PgC'
#time_period = 'Annual'
#multi_area  = True
e3sm_comp   = 'elm'

if(conv_unit == 'gCm2_gCm2'):
   CONV_C       = 1
   myDict_units = myDict_g_m2_day_units[time_period]
if(conv_unit == 'gC_PgC'):
   CONV_C       = CONV_gC_PgC
   myDict_units = myDict_PgC_units[time_period]

yr_step  = 1
yr_start = 2071
yr_end   = 2090

# Make spatial plot of individual variables
for ind, var in enumerate(varnames):

  # iterate through runs
  for i, run in enumerate(runs):

    caseid   = myDict_caseid[run]

    # Read ELM model outputs for the transient period
    ds_model = read_transient_outputs(fpath_scratch, caseid, e3sm_comp, varnames, yr_start, yr_end, yr_step, CONV_C, time_period, multi_area)

    # Average across selected years
    ds_model = ds_model[var].groupby('time.year').mean() # Mean for each year
    #ds_model = ds_model.mean(dim='year')                # Mean across years

    ds_model.to_netcdf(path = var + '_' + run + '_' + str(yr_start) + '_' + str(yr_end) + '.nc', mode='w')
    print('netcdf file generated')

    # Merge all datasets into a single dataset
    ds_model = ds_model.expand_dims(Set = [run])

    if (i == 0):
      da_plot_merge = ds_model
    else:
      da_plot_merge = xr.merge([da_plot_merge, ds_model])

  # Convert datset to dataarray
  da_plot_merge = da_plot_merge.to_array()

  # Convert to pandas dataframe
  # and grouby lat and lon so that T-test can be applied on each grid
  df = da_plot_merge.to_dataframe(name=var).reset_index()
  df = df.dropna()
  df = df.groupby(['lat','lon'])

  # Calculate the T-test for the means of two independent samples of scores
  p_values = df.apply(apply_ttest_ind, var, sel_set='FULL_FDBK')
  #p_values = df.apply(apply_ttest_ind, var, sel_set='AG_FDBK')
  #p_values = df.apply(apply_ttest_ind, var, sel_set='C_FDBK')
  # Convert to xarray
  p_values = p_values.to_xarray()

  # Mean across years
  da_plot_merge = da_plot_merge.mean(dim='year')

  # Compute difference between carbon scaling and no carbon scaling
  da_plot_merge.loc[:,'FULL_FDBK',:,:] = da_plot_merge.loc[:,'FULL_FDBK',:,:] - da_plot_merge.loc[:,'CONTROL',:,:]
  #da_plot_merge.loc[:,'AG_FDBK',:,:] = da_plot_merge.loc[:,'AG_FDBK',:,:] - da_plot_merge.loc[:,'CONTROL',:,:]
  #da_plot_merge.loc[:,'C_FDBK',:,:] = da_plot_merge.loc[:,'C_FDBK',:,:] - da_plot_merge.loc[:,'CONTROL',:,:]

  # Plot after dropping no carbon scaling sets
  da_plot_merge = da_plot_merge.drop_sel(Set = 'CONTROL')
  # Now drop cordinate with only one dimension
  da_plot_merge = da_plot_merge.isel(Set=0)
  da_plot_merge = da_plot_merge.isel(variable=0)

  cmap_col = 'bwr' #'PiYG'
  #facet_plot(da_plot_merge, colplot='Set', colwrap=len(da_plot_merge.Set)-1,
  #                  cmap_col=cmap_col, cbar_label=myDict_units[var], \
  #                  fig_wt=10*(len(da_plot_merge.Set)-1)+0.3, fig_ht=10+0.2, \
  #                  fname=var+'_diff.png')

  # Estimate metrics for the data
  metrics_dict = create_metrics(da_plot_merge)
  plot_stats=(metrics_dict['max'], metrics_dict['mean'], metrics_dict['min'])

  levels = None
  norm   = None
  clevels = myDict_diff_clevels[var]
  if len(clevels) > 0:
    levels = clevels
    norm = colors.BoundaryNorm(boundaries=levels, ncolors=len(levels))

  title = myDict_units[var] + '(' + str(yr_start) + '-' + str(yr_end) + ')\nFULL_FDBK - CONTROL'
  #title = myDict_units[var] + '(' + str(yr_start) + '-' + str(yr_end) + ')\nAG_FDBK - CONTROL'
  #title = myDict_units[var] + '(' + str(yr_start) + '-' + str(yr_end) + ')\nC_FDBK - CONTROL'

  fname = var+'_FULL_FDBK_diff.png'
  #fname = var+'_AG_FDBK_diff.png'
  #fname = var+'_C_FDBK_diff.png'
  xr_plot_global(da_plot_merge, stats=plot_stats, \
                  cmap_col=cmap_col, title=title, fig_wt=8*2, fig_ht=12, \
                  levels=levels, norm=norm, \
                  fname=fname,
                  stipple_data = p_values)
