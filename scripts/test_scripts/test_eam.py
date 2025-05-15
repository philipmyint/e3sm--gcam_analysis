import json
import multiprocessing
import sys
sys.path.append('./../')
import os
import time
import xarray as xr
from utility_constants import *
from utility_functions import check_substrings_in_list, get_all_files_in_path
from utility_e3sm_netcdf import get_netcdf_files_between_start_and_end_years

###---------------Begin execution---------------###
if __name__ == '__main__':

    start_time = time.time()

    directory = '/lcrc/group/e3sm/ac.myint1/e3sm_scratch/20250401_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3/run'
    file = '20250401_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3.eam.h0.2019-09.nc'
    file = os.path.join(directory, file)
    ds = xr.open_dataset(file)[['PRECC', 'lat', 'lon']]
    ds['PRECC'].expand_dims(dim=['lat', 'lon'])

    print(ds['lat'].to_dataframe().unique())
    print(ds['PRECC'].min().item(), ds['PRECC'].mean().item())
    '''df = ds.to_dataframe().groupby(['lat', 'lon']).mean()
    print(df)
    ds = df.to_xarray()
    print(ds)'''

    # Print the total execution time needed to complete all data extraction operations.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time to extract all spatial data outputs: {elapsed_time:.2f} seconds")