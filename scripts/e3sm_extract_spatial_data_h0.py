import json
import multiprocessing
import sys
import time
import xarray as xr
from utility_constants import *
from utility_functions import check_substrings_in_list, get_all_files_in_path
from utility_e3sm_netcdf import get_netcdf_files_between_start_and_end_years

def process_inputs(inputs):    
    """ 
    Processes a dictionary of inputs (keys are options, values are choices for those options) for extracting spatial data from E3SM-generated NetCDF 
    files and printing that data to smaller, focused NetCDF files that contain spatial data only for the variables specified in the dictionary.

    Parameters:
        inputs: Dictionary containing the user choice inputs for different options. This dictionary may be incomplete or have invalid values.

    Returns:
        List of dictionaries, where each dictionary has been processed so that it is complete in all options for a single focused NetCDF file.
    """
    # If the user specified output_files to be a string (indicating a single output file), turn that string into a list containing that string.
    if isinstance(inputs['output_files'], (str, int, float)):
        inputs['output_files'] = [inputs['output_files']]

    # If the user specified only a single value (string, integer, float, or a single-element list) for the other plotting options, assume that they
    # want to use that value for all of the variables.
    input_types_to_modify = ['netcdf_substrings', 'variables', 'start_years', 'end_years']
    for input_type in input_types_to_modify:
        if isinstance(inputs[input_type], (str, int, float)):
            inputs[input_type] = [inputs[input_type]]
        if len(inputs[input_type]) == 1:
            inputs[input_type] = inputs[input_type]*len(inputs['output_files'])
        # This processes 'netcdf_substrings' and 'variables' so that they are lists of lists.
        if input_type in ['netcdf_substrings', 'variables'] and not isinstance(inputs[input_type][0], list):
            inputs[input_type] = [inputs[input_type] for i in range(len(inputs[input_type]))]

    # Now that the dictionary has been populated with complete data extraction options for each output file, separate it into a list of dictionaries,
    # where each of these smaller dictionaries contain the data extraction options for a single output file. Return this list of dictionaries.
    list_of_inputs = []
    for file_index in range(len(inputs['output_files'])):
        output_file = inputs['output_files'][file_index]
        netcdf_substrings = inputs['netcdf_substrings'][file_index]
        variables = inputs['variables'][file_index]
        start_years = inputs['start_years'][file_index]
        end_years = inputs['end_years'][file_index]
        inputs_for_this_output_file = {'simulation_path': inputs['simulation_path'], 'output_files': output_file, 'netcdf_substrings': netcdf_substrings}
        inputs_for_this_output_file.update({'variables': variables, 'start_years': start_years, 'end_years': end_years})
        list_of_inputs.append(inputs_for_this_output_file)
    return list_of_inputs

def process_dataset(ds):    
    """ 
    Processes an xarray Dataset to add a total precipitation variable (including change units from m/s to mm/year) and 
    add a variable for the atmospheric mole fraction of CO2 in units of ppm.

    Parameters:
        ds: Dataset to process.

    Returns:
        The processed Dataset.
    """
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
        ds['ZCO2'] = mole_fraction_TO_ppm*ds['PCO2']/(ds['PBOT'] - partial_pressure_H2O)
        ds['ZCO2'].attrs = {'units':'ppm', 'description':'CO2 mole fraction in dry air'}
    
    return ds

def extract_spatial_data_from_netcdf_files(inputs):
    """ 
    Extracts time series data from E3SM-generated h0 NetCDF files of a particular type that are located in a simulation directory and puts this data 
    into a smaller, focused NetCDF file that contains annual-mean spatial data only for the variables specified by the user in the inputs dictionary.

    Parameters:
        input: Dictionary containing the user data-extraction inputs for different options. This dictionary is assumed to be complete (pre-processed).

    Returns:
        N/A.
    """
    # Record the start time and extract all the user-chosen options from the inputs dictionary.
    start_time = time.time()
    simulation_path = inputs['simulation_path']
    output_file = inputs['output_files'] 
    netcdf_substrings = inputs['netcdf_substrings'] 
    variables = inputs['variables']  
    start_year = inputs['start_years']   
    end_year = inputs['end_years']   

    # Get all NetCDF files for this particular type that fall within the specified start and end years.
    netcdf_files = get_all_files_in_path(simulation_path, file_name_substrings=netcdf_substrings, file_extension='.nc')
    netcdf_files = get_netcdf_files_between_start_and_end_years(netcdf_files, start_year, end_year)
    
    # Collect the NetCDF files (one for each month between the start and end years) in an xarray Dataset and store only the specified variables.
    ds = xr.open_mfdataset(netcdf_files, decode_times=True, combine='nested', concat_dim='time', data_vars='minimal', parallel=True)[variables]
    
    # Shift output back by one month to get rid of the extra month (January in the next year after end_year) that somehow gets added.
    ds['time'] = xr.CFTimeIndex(ds.get_index('time').shift(-1, 'ME'))

    # Take the mean over all months in each year so that the Dataset records only an annual mean value for each variable at each lat/lon coordinate.
    ds = ds.groupby('time.year').mean()

    # Add total precipitation (in units of mm/year) and CO2 concentration variables to the Dataset.
    ds = process_dataset(ds)

    # Write the Dataset to a NetCDF file.
    ds.to_netcdf(output_file, mode='w')

    # Print the time needed to create the smaller NetCDF file.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for {output_file}: {elapsed_time:.2f} seconds")


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file(s) on the command line.
    start_time = time.time()
    if len(sys.argv) < 2:
        print('Usage: python e3sm_extract_spatial_data_h0.py `path/to/json/input/file(s)\'')
        sys.exit()

    # Read and load the JSON file(s) into a list of dictionaries.
    inputs = []
    for index in range(1, len(sys.argv)):
        input_file = sys.argv[index]
        with open(input_file) as f:
            inputs.extend(json.load(f))

    # Process each dictionary to produce a list of smaller dictionaries, where each specifies data extraction options for a single NetCDF output file.
    list_of_inputs_for_each_output_file = []
    for index in range(len(inputs)):
        list_of_inputs_for_each_output_file.extend(process_inputs(inputs[index]))

    # Create all output NetCDF files in parallel.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(extract_spatial_data_from_netcdf_files, list_of_inputs_for_each_output_file)

    # Print the total execution time needed to complete all data extraction operations.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time to extract all spatial data outputs: {elapsed_time:.2f} seconds")