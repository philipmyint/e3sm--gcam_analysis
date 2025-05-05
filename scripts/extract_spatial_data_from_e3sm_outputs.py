import json
import multiprocessing
import sys
import time
import xarray as xr
from utility_constants import *
from utility_functions import get_all_files_in_path
from utility_e3sm_netcdf import get_netcdf_files_between_start_and_end_years

def process_inputs(inputs):    
    """ 
    Processes a dictionary of inputs (keys are options, values are choices for those options) for extracting spatial data from E3SM-generated NetCDF 
    files and printing the output to smaller, focused NetCDF files that contain spatial data only for the variables specified in the dictionary.

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

    # Now that the dictionary has been populated with a complete data extraction options for each variable, separate it into a list of dictionaries,
    # where each of these smaller dictionaries contain the data extraction options for a single variable. Return this list of dictionaries.
    list_of_inputs_for_each_output_file = []
    for file_index in range(len(inputs['output_files'])):
        output_file = inputs['output_files'][file_index]
        netcdf_substrings = inputs['netcdf_substrings'][file_index]
        variables = inputs['variables'][file_index]
        start_years = inputs['start_years'][file_index]
        end_years = inputs['end_years'][file_index]
        inputs_for_this_output_file = {'simulation_path': inputs['simulation_path'], 'output_files': output_file, 'netcdf_substrings': netcdf_substrings}
        inputs_for_this_output_file.update({'variables': variables, 'start_years': start_years, 'end_years': end_years})
        list_of_inputs_for_each_output_file.append(inputs_for_this_output_file)
    return list_of_inputs_for_each_output_file

def extract_spatial_data_from_netcdf_files(inputs):
    """ 
    Extracts time series data from E3SM-generated NetCDF files of a particular type that are located in a simulation directory and puts this data into
    a smaller, focused NetCDF file that contains annual-mean spatial data only for the variables specified by the user in the inputs dictionary.

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
    ds = xr.open_mfdataset(netcdf_files, decode_times=True, combine='nested', concat_dim='time', data_vars='minimal')[variables]
    
    # Shift output back by one month to get rid of the extra month (January in the next year after end_year) that somehow gets added.
    ds['time'] = xr.CFTimeIndex(ds.get_index('time').shift(-1, 'ME'))

    # Take the mean over all months in each year so that the Dataset records only an annual mean value for each variable at each lat/lon coordinate.
    ds = ds.groupby('time.year').mean()

    # Write the Dataset to a NetCDF file.
    ds.to_netcdf(output_file, mode='w')

    # Print the time needed to create the smaller NetCDF file.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for {output_file}: {elapsed_time:.2f} seconds")


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file on the command line.
    if len(sys.argv) != 2:
        print('Usage: extract_spatial_data_from_netcdf_files.py `path/to/json/input/file\'')
        sys.exit()

    # Read and load the JSON file.
    start_time = time.time()
    input_file = sys.argv[1]
    with open(input_file) as f:
        inputs = json.load(f)

    # Process each block in the JSON file to produce a list of dictionaries, where each specifies data extraction options for a NetCDF output file.
    list_of_inputs_for_each_output_file = []
    for index in range(len(inputs)):
        list_of_inputs_for_each_output_file.extend(process_inputs(inputs[index]))

    # Create all output NetCDF files in parallel.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(extract_spatial_data_from_netcdf_files, list_of_inputs_for_each_output_file)

    # Print the total execution time needed to complete all data extraction operations.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time to extract all spatial data outputs specified in {sys.argv[1]}: {elapsed_time:.2f} seconds")