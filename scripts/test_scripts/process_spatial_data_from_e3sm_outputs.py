import json
import multiprocessing
import sys
sys.path.append('./../')
import time
import xarray as xr
from utility_constants import *
from utility_functions import check_substrings_in_list, get_all_files_in_path

def process_netcdf_files(file):    
    """ 
    Processes an E3SM-generated NetCDF file to add a total precipitation variable (including change units from m/s to mm/year) and 
    adding a variable for the atmospheric mole fraction of CO2 in units of ppm.

    Parameters:
        file: NetCDF to process.

    Returns:
        The processed NetCDF file.
    """
    ds = xr.open_dataset(file)
    variables = list(ds.keys())

    # Add a new variable for the total precipitation, and update the variable labels.
    precip_variables = ['PRECC', 'PRECL', 'PRECSC', 'PRECSL']
    if check_substrings_in_list(precip_variables, variables, all_or_any='all'):
        ds[precip_variables] *= years_TO_s*m_TO_mm
        for precip_variable in precip_variables:
            ds[precip_variable].attrs['units'] = 'mm/year'
        ds['PRECIP'] = ds[precip_variables].to_array(dim='variable').sum(dim='variable')
        ds['PRECIP'].attrs = {'units':'mm/year', 'description':'Total precipitation rate'}

    # Calculate the atmospheric mole fraction of CO2 in units of ppm, which is defined as a mole fraction in dry air (with humidity subtracted out).
    pressure_variables = ['PBOT', 'PCO2', 'QBOT']
    if check_substrings_in_list(pressure_variables, variables, all_or_any='all'):
        # See derivation of H2O partial pressure formula at https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html.
        partial_pressure_H2O = ds['PBOT']*ds['QBOT']/(0.622 + (0.378*ds['QBOT']))
        # Add this CO2 variable to 
        ds['XCO2'] = mole_fraction_TO_ppm*ds['PCO2']/(ds['PBOT'] - partial_pressure_H2O)
        ds['XCO2'].attrs = {'units':'ppm', 'description':'CO2 mole fraction in dry air'}

    print(ds)

###---------------Begin execution---------------###
if __name__ == '__main__':

    start_time = time.time()

    directory = './../../2025_DiVittorio_et_al/'
    netcdf_files = get_all_files_in_path(directory, file_extension='.nc')
    print(netcdf_files)

    # Create all output NetCDF files in parallel.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(process_netcdf_files, netcdf_files)

    # Print the total execution time needed to complete all data extraction operations.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time to process all spatial data outputs in {directory}: {elapsed_time:.2f} seconds")