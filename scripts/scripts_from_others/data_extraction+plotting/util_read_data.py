import xarray as xr
import os
import glob

from util_myDict import *

# -----------------------------------------------------------
def read_model_output(yr_start, yr_end, yr_step, fpath, caseid, e3sm_comp, mon_day_str, varnames, decode_times=True):
    """Read ELM model output for select variables
        param: yr_start      start year for reading model output
        param: yr_end        end year for reading model output
        param: fpath         directory path
        param: caseid        model run case id
        param: varnames      variable name for subsetting
        :return:             data array with model output
        """

    # Read names of all NetCDF files within the given year range
    fnames = []
    for yr in range(int(yr_start), int(yr_end)+1, yr_step):
        fnames.append(fpath + '/' + caseid + '.' + e3sm_comp + '.h0.' + str(yr).zfill(4) + mon_day_str + '-00000.nc')

    # Open a multiple netCDF data file and load the data into xarrays
    with xr.open_mfdataset(fnames, decode_times=decode_times, combine='nested', concat_dim='time', data_vars='minimal') as ds:

        # Only keep select variables in the data array
        ds = ds[varnames]

    return(ds)

# -----------------------------------------------------------
def read_monthly_model_output(yr_start, yr_end, yr_step, fpath, caseid, e3sm_comp, varnames, decode_times=True):
    """Read ELM model output for select variables
        param: yr_start      start year for reading model output
        param: yr_end        end year for reading model output
        param: fpath         directory path
        param: caseid        model run case id
        param: varnames      variable name for subsetting
        :return:             data array with model output
        """

    # Read names of all NetCDF files within the given year range
    fnames = []
    for yr in range(int(yr_start), int(yr_end)+1, yr_step):
       for mon in range(1,13):
          fnames.append(fpath + '/' + caseid + '.' + e3sm_comp + '.h0.' + str(yr).zfill(4) + '-' + str(mon).zfill(2) + '.nc')
       #fnames = glob.glob(fpath + '/' + caseid + '.' + e3sm_comp + '.h0.' + str(yr).zfill(4) + '*.nc')

    # Open a multiple netCDF data file and load the data into xarrays
    with xr.open_mfdataset(fnames, decode_times=decode_times, combine='nested', concat_dim='time', data_vars='minimal') as ds:

       # Only keep select variables in the data array
       ds = ds[varnames]

    return(ds)

# -----------------------------------------------------------
def read_monthly_model_output_large(yr_start, yr_end, yr_step, fpath, caseid, e3sm_comp, varnames, decode_times=True):
    """Read ELM model output for select variables
        param: yr_start      start year for reading model output
        param: yr_end        end year for reading model output
        param: fpath         directory path
        param: caseid        model run case id
        param: varnames      variable name for subsetting
        :return:             data array with model output
        """

    # Read names of all NetCDF files within the given year range
    for yr in range(int(yr_start), int(yr_end)+1, yr_step):
       fnames = []
       for mon in range(1,13):
          fnames.append(fpath + '/' + caseid + '.' + e3sm_comp + '.h0.' + str(yr).zfill(4) + '-' + str(mon).zfill(2) + '.nc')
       #fnames = glob.glob(fpath + '/' + caseid + '.' + e3sm_comp + '.h0.' + str(yr).zfill(4) + '*.nc')

       # Open a multiple netCDF data file and load the data into xarrays
       with xr.open_mfdataset(fnames, combine='nested', concat_dim='time', data_vars='minimal') as ds:

          # Only keep select variables in the data array
          ds = ds[varnames]

       if (yr == yr_start):
          ds_merge = ds
       else:
          ds_merge = xr.concat([ds_merge, ds], dim='time')

    return(ds_merge)

# -----------------------------------------------------------
def read_col_lev_model_output(yr_start, yr_end, fpath, caseid, e3sm_comp):
    """Read ELM model output for select variables
    :param: yr_start     start year for reading model output
    :param: yr_end       end year for reading model output
    :param: fpath        directory path
    :param: caseid       model run case id
    :return:             data array with model output
    """

    # Read names of all NetCDF files within the given year range
    fnames = []
    for yr in range(int(yr_start), int(yr_end)+1):
       for mon in range(1,13):
          fnames.append(fpath + '/' + caseid + '.' + e3sm_comp + '.h1.' + str(yr).zfill(4) + '-' + str(mon).zfill(2) + '.nc')

    # Open a multiple netCDF data file and load the data into xarray
    ds = xr.open_mfdataset(fnames, combine='nested', concat_dim='time')

    return(ds)

# -----------------------------------------------------------
# Read ELM model outputs for the transient period
def read_transient_outputs(fpath_scratch, caseid, e3sm_comp, varnames, yr_start, yr_end, yr_step, CONV_C, time_period, multi_area=True):

   # Read grid level ELM model output for transient period
   if(e3sm_comp == 'elm'):
     fpath = fpath_scratch + caseid + '/archive/lnd/hist/'
     var_landfrac = 'landfrac'
   elif(e3sm_comp == 'eam'):
     fpath = fpath_scratch + caseid + '/archive/atm/hist/'
     var_landfrac = 'LANDFRAC'
   if(os.path.exists(fpath) == False):
     fpath = fpath_scratch + caseid + '/run/'

   ds_model    = read_monthly_model_output_large(yr_start, yr_end, yr_step, fpath, caseid, e3sm_comp, varnames)
   ds_area     = read_monthly_model_output(yr_start, yr_start, yr_step, fpath, caseid, e3sm_comp, ['area'])
   ds_landfrac = read_monthly_model_output(yr_start, yr_end, yr_step, fpath, caseid, e3sm_comp, [var_landfrac])

   # Multiply by area and convert gC to specified  unit
   if(multi_area):
      if(e3sm_comp == 'elm'):
        ds_model = ds_model * ds_landfrac[var_landfrac] * ds_area['area'] * CONV_KM2_M2 * CONV_C
      if(e3sm_comp == 'eam'):
        re              = SHR_CONST_REARTH * 0.001   # radius of earth (km)
        ds_area['area'] = ds_area['area'] * (re*re) # Convert EAM area from steradians to km2
        for ind, var in enumerate(varnames):
          if(var in ['SFCO2_LND']):
            ds_model[var] = ds_model[var] * ds_landfrac[var_landfrac] * ds_area['area'] * CONV_KM2_M2 * CONV_C * CONC_CO2_C
          elif(var in ['SFCO2_OCN']):
            ds_model[var] = ds_model[var] * (1 - ds_landfrac[var_landfrac]) * ds_area['area'] * CONV_KM2_M2 * CONV_C * CONC_CO2_C
          elif(var in ['SFCO2_FFF','SFCO2']):
            ds_model[var] = ds_model[var] * ds_area['area'] * CONV_KM2_M2 * CONV_C * CONC_CO2_C
   else:
      ds_model = ds_model * CONV_C

   # Shift output by one month
   ds_model['time'] = xr.CFTimeIndex(ds_model.get_index('time').shift(-1, 'M'))

   # Convert time to datetime
   ds_model['time'] = ds_model.indexes['time'].to_datetimeindex()

   # Scale by conversion factor for time
   for ind, var in enumerate(varnames):
      ds_model[var]   = ds_model[var]/conv_factor[time_period][var]

   return(ds_model)

# -----------------------------------------------------------
# Read ELM model outputs for the ad and final spinup period
def read_spinup_outputs(fpath_scratch, caseid, e3sm_comp, varnames, CONV_C, time_period, multi_area=True):

   yr_step  = 1

   # Read ELM model output for ad spinup
   run = 'ad_spinup'
   caseid = myDict_caseid[run]
   fpath = fpath_scratch + caseid + '/run/'
   ds_model_ad_spinup    = read_monthly_model_output_large(yr_start[run], yr_end[run], yr_step, fpath, caseid, e3sm_comp, varnames)
   ds_area_ad_spinup     = read_monthly_model_output(yr_start[run], yr_start[run], yr_step, fpath, caseid, e3sm_comp, ['area'])
   ds_landfrac_ad_spinup = read_monthly_model_output(yr_start[run], yr_start[run], yr_step, fpath, caseid, e3sm_comp, ['landfrac'])

   # Read ELM model output for final spinup
   run = 'final_spinup'
   caseid = myDict_caseid[run]
   fpath = fpath_scratch + caseid + '/run/'
   ds_model_final_spinup    = read_monthly_model_output_large(yr_start[run], yr_end[run], yr_step, fpath, caseid, e3sm_comp, varnames)
   ds_area_final_spinup     = read_monthly_model_output(yr_start[run], yr_start[run], yr_step, fpath, caseid, e3sm_comp, ['area'])
   ds_landfrac_final_spinup = read_monthly_model_output(yr_start[run], yr_start[run], yr_step, fpath, caseid, e3sm_comp, ['landfrac'])

   # Shift final spinup time
   ds_model_final_spinup['time'] = xr.CFTimeIndex(ds_model_final_spinup.get_index('time').shift(200*365, 'D'))
   print(ds_model_final_spinup)

   # Multiply by area and convert gC to specified unit
   if(multi_area):
      ds_model_ad_spinup    = ds_model_ad_spinup    * ds_area_ad_spinup['area']    * ds_landfrac_ad_spinup['landfrac']    * CONV_KM2_M2 * CONV_C
      ds_model_final_spinup = ds_model_final_spinup * ds_area_final_spinup['area'] * ds_landfrac_final_spinup['landfrac'] * CONV_KM2_M2 * CONV_C
   else:
      ds_model_ad_spinup    = ds_model_ad_spinup    * CONV_C
      ds_model_final_spinup = ds_model_final_spinup * CONV_C

   # Shift output by one month
   ds_model_ad_spinup['time']    = xr.CFTimeIndex(ds_model_ad_spinup.get_index('time').shift(-1, 'M'))
   ds_model_final_spinup['time'] = xr.CFTimeIndex(ds_model_final_spinup.get_index('time').shift(-1, 'M'))

   # Estimate total annual flux
   ds_model_ad_spinup    = create_mean_annual_da(ds_model_ad_spinup, est_mon_total=False)
   ds_model_final_spinup = create_mean_annual_da(ds_model_final_spinup, est_mon_total=False)

   # Merge ad and final spinups
   ds_model = ds_model_ad_spinup
   ds_model = xr.merge([ds_model_ad_spinup, ds_model_final_spinup])

   # Scale by conversion factor
   for ind, var in enumerate(varnames):
      ds_model[var] = ds_model[var]/conv_factor[time_period][var]

   return(ds_model)

#----------------------------------------------------------
