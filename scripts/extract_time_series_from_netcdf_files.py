import json
import multiprocessing
import numpy as np
import pandas as pd
import sys
import time
import xarray as xr
from utility_constants import *
from utility_functions import *
from utility_dataframes import write_dataframe_to_fwf, move_columns_next_to_each_other_in_dataframe

def extract_year_and_month_from_name_of_netcdf_file(file):

    # Files are named *YYYY-MM.nc, so the year starts 7 characters before the .nc extension.
    index = file.find('.nc')
    index -= len('YYYY-MM')
    year = int(file[index:index+len('YYYY')])
    month = int(file[index+len('YYYY-'):index+len('YYYY-MM')])
    return year, month

def get_netcdf_files_between_start_and_end_years(file_paths, start_year, end_year):
    new_file_paths = []
    for file in file_paths:
        year, _ = extract_year_and_month_from_name_of_netcdf_file(file)
        if year >= start_year and year <= end_year:
            new_file_paths.append(file)
    return new_file_paths

def find_gridcell_areas_in_netcdf_file(file):
    xrds = xr.open_dataset(file)
    variables = ['area']
    areas = create_numpy_array_from_xrds(xrds, variables, [0]*len(variables))
    if 'elm.h0' in file:
        variables = ['landfrac', 'landmask', 'pftmask']
        landfrac, landmask, pftmask = create_numpy_array_from_xrds(xrds, variables, [0]*len(variables))
        areas *= landfrac*landmask*pftmask*km2_TO_m2
    elif 'eam.h0' in file:
        areas *= SURF_AREA
    return areas, xrds

def process_dataframe(df):

    start_year = df['Year'].min()
    end_year = df['Year'].max()
    num_years = end_year - start_year

    # Convert precipitation variables from m/s to mm/month, add a new column for the total precipitation, and update the column labels.
    substrings = ['PRECC', 'PRECL', 'PRECSC', 'PRECSL']
    if check_substrings_in_list(substrings, df.columns, all_or_any='all'):
        columns_to_modify = [label for label in df.columns for substring in substrings if substring in label]
        seconds_in_months_tiled = np.tile(NUM_SECONDS_IN_MONTHS, num_years+1).reshape(-1,1)
        df[columns_to_modify] *= seconds_in_months_tiled*m_TO_mm
        df['PRECIP (mm/month)'] = df[columns_to_modify].sum(axis=1)
        condition = lambda x: any(substring in x for substring in ['PRECC', 'PRECL', 'PRECSC', 'PRECSL'])
        new_column_label_function = lambda x: x.replace('(m/s)', '(mm/month)')
        df.columns = modify_list_based_on_condition(df.columns, condition, new_column_label_function)

    # Calculate the atmospheric mole fraction of CO2 in units of ppm, which is defined as a mole fraction in dry air (with humidity subtracted out).
    substrings = ['PBOT', 'PCO2', 'QBOT']
    if check_substrings_in_list(substrings, df.columns, all_or_any='all'):
        columns_to_modify = [label for label in df.columns for substring in substrings if substring in label]
        # See derivation of H2O partial pressure formula at https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html.
        partial_pressure_H2O = df['PBOT (Pa)']*df['QBOT (kg/kg)']/(0.622 + (0.378*df['QBOT (kg/kg)']))
        df['XCO2 (ppm)'] = mole_fraction_TO_ppm*df['PCO2 (Pa)']/(df['PBOT (Pa)'] - partial_pressure_H2O)
        df = move_columns_next_to_each_other_in_dataframe(df, 'PCO2 (Pa)', 'XCO2 (ppm)')

    # Convert fluxes and stocks that have units of g or kg to Pg.
    old_labels = ['(gC', '(g/', '(kg']
    new_labels = ['(PgC', '(Pg/', '(Pg']
    multipliers = [1/Pg_TO_g, 1/Pg_TO_g, 1/Pg_TO_kg]
    for index, old_label in enumerate(old_labels):
        # Make sure that the quantity is actually a flux or stock and not a mass fraction (i.e., does not have units of kg/kg).
        columns_to_modify = [label for label in df.columns if (old_label in label and '(kg/kg)' not in label)]
        condition = lambda x: old_label in x and '(kg/kg)' not in x
        new_column_label_function = lambda x: x.replace(old_label, new_labels[index]) 
        df[columns_to_modify] *= multipliers[index] 
        df.columns = modify_list_based_on_condition(df.columns, condition, new_column_label_function)    
    
    # Convert flux units from per second to per month.
    old_label = '/s)'
    columns_to_modify = [label for label in df.columns if old_label in label]
    seconds_in_months_tiled = np.tile(NUM_SECONDS_IN_MONTHS, num_years+1).reshape(-1,1)
    df[columns_to_modify] *= seconds_in_months_tiled
    condition = lambda x: old_label in x
    new_column_label_function = lambda x: x.replace(old_label, '/month)')
    df.columns = modify_list_based_on_condition(df.columns, condition, new_column_label_function)

    # Convert CO2 fluxes and stocks from PgCO2 to PgC.
    substrings = ['SFCO2', 'TMCO2']
    columns_to_modify = [label for label in df.columns if check_substrings_in_string(substrings, label, all_or_any='any')]
    df[columns_to_modify] *= MM_C/MM_CO2
    old_label = 'Pg'
    condition = lambda x: old_label in x and check_substrings_in_string(substrings, x, all_or_any='any')
    new_column_label_function = lambda x: x.replace(old_label, 'PgC')
    df.columns = modify_list_based_on_condition(df.columns, condition, new_column_label_function)

    return df

def extract_netcdf_file_into_dataframe(file, variables, calculation_type='mean'):
    """ 
    Extracts time series data from NetCDF files located in a simulation directory into a Pandas DataFrame that then gets written to an output file. 
    The NetCDF files can be of more than one type (e.g., one set generated from the ELM model and another set from the EAM model in E3SM).
    Each NetCDF file corresponds to simulation results for one month in a particular year.

    Parameters:
        outputs_path: Complete path of the directory containing the NetCDF output files from running a simulation.
        extracted_outputs_file: Name of the file to which the DataFrame will be printed.
        output_substrings: List where each element is itself a list of substrings. Each list corresponds to a particular type of NetCDF file and
        indicates the substrings that must be in the names of that NetCDF file type.
        variables: List where each element is itself a list of variables we want to extract for each type of NetCDF file. 
        The aggregate of all variables contained these lists together will form the columns of the DataFrame.
        calculation_types: List where each element is itself a list of strings that indicate what type of calculation we want to perform to produce 
        time series data for each type of NetCDF file. Current options include 'area_weighted_mean_or_sum', 'mean', and 'sum'.
        process_variables: Boolean that indicates if we want to perform further processing on the outputs in the DataFrame. This includes adding
        new variables not in the NetCDF files (e.g., total precipitation, mole fraction CO2) or changing the units (e.g., Pg instead of g).

    Returns:
        N/A
    """
    areas, xrds = find_gridcell_areas_in_netcdf_file(file)
    xrds = xrds[variables]
    df = xrds.to_dataframe()

    units = [xrds[var].attrs['units'] for var in xrds.data_vars]
    df.columns = add_lists_elementwise(df.columns, units, list2_are_units=True)

    if calculation_type == 'area_weighted_mean_or_sum':
        df *= areas
        df = df.sum().to_frame().T

        # Variables that are not fluxes and stocks (so that they are not per-area quantities), should be global area-weighted means
        # rather than area-weighted sums, and therefore we need to divide them by the total area.
        total_area = np.sum(areas)
        columns_to_modify = [label for label in df.columns if not check_substrings_in_string(['/m^2', '/m2'], label, all_or_any='any')]
        df[columns_to_modify] /= total_area

        # Fluxes and stocks are global area-weighted sums and have been multiplied by areas, so we must update the labels to remove the '/m2' part.
        for old_label in ['/m^2', '/m2']:
            condition = lambda x: old_label in x
            new_column_label_function = lambda x: x.replace(old_label, '')
            df.columns = modify_list_based_on_condition(df.columns, condition, new_column_label_function)    

    elif calculation_type == 'sum':
        df = df.sum().to_frame().T
    elif calculation_type == 'mean':
        df = df.mean().to_frame().T

    year, month = extract_year_and_month_from_name_of_netcdf_file(file)
    column_names_with_year_and_month_first = ['Year', 'Month']
    column_names_with_year_and_month_first.extend(df.columns)
    df['Year'] = year
    df['Month'] = month
    return df[column_names_with_year_and_month_first]

def extract_time_series_data_from_netcdf_files(outputs_path, extracted_outputs_file, outputs_substrings, 
            variables, calculation_types, process_variables=True, start_year=2015, end_year=2100):
    """ 
    Extracts time series data from NetCDF files located in a simulation directory into a Pandas DataFrame that then gets written to an output file. 
    The NetCDF files can be of more than one type (e.g., one set generated from the ELM model and another set from the EAM model in E3SM).
    Each NetCDF file corresponds to simulation results for one month in a particular year.

    Parameters:
        outputs_path: Complete path of the directory containing the NetCDF output files from running a simulation.
        extracted_outputs_file: Name of the file to which the DataFrame will be printed.
        output_substrings: List where each element is itself a list of substrings. Each list corresponds to a particular type of NetCDF file and
        indicates the substrings that must be in the names of that NetCDF file type.
        variables: List where each element is itself a list of variables we want to extract for each type of NetCDF file. 
        The aggregate of all variables contained these lists together will form the columns of the DataFrame.
        calculation_types: List where each element is itself a list of strings that indicate what type of calculation we want to perform to produce 
        time series data for each type of NetCDF file. Current options include 'area_weighted_mean_or_sum', 'mean', and 'sum'.
        process_variables: Boolean that indicates if we want to perform further processing on the outputs in the DataFrame. This includes adding
        new variables not in the NetCDF files (e.g., total precipitation, mole fraction CO2) or changing the units (e.g., Pg instead of g).

    Returns:
        N/A
    """
    # This list will store a DataFrame for each type of NetCDF file, and all elements of this list will later be merged into a single DataFrame.
    dataframes = []

    # Iterate over all NetCDF file types.
    for index in range(len(variables)):

        # Get all NetCDF files for this particular type.
        netcdf_files = get_all_files_in_path(outputs_path, 
                            file_name_substrings=outputs_substrings[index], file_extension='.nc')
        netcdf_files = get_netcdf_files_between_start_and_end_years(netcdf_files, start_year, end_year)
    
        # Use multiprocessing to extract data from all the files. Put the data from each file into DataFrame, store all such DataFrames in a list.
        arguments = list(zip(netcdf_files, [variables[index]]*len(netcdf_files), 
                    [calculation_types[index]]*len(netcdf_files)))
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            dataframes_for_each_nc_file = list(pool.starmap(extract_netcdf_file_into_dataframe, arguments))
        
        # Concatenate all DataFrames in the list together to form a single DataFrame for this NetCDF type. Sort by year and month.
        df = dataframes_for_each_nc_file[0]
        for df_index in range(1, len(dataframes_for_each_nc_file)):
            df = pd.concat([df, dataframes_for_each_nc_file[df_index]]).reset_index(drop=True)
        df.sort_values(['Year', 'Month'], inplace=True)
    
        # Add the DataFrame for this NetCDF type into the master list.
        dataframes.append(df)

    # Combine the DataFrames for all NetCDF types into a single large DataFrame by merging on the year and month columns.
    df = dataframes[0]
    for df_index in range(1, len(dataframes)):
        df = pd.merge(df, dataframes[df_index], on=['Year', 'Month'], how='inner')

    # Process the variables if specified to do so.
    if process_variables:
        df = process_dataframe(df)

    # Write the DataFrame to the specified output file.
    write_dataframe_to_fwf(extracted_outputs_file, df)


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file on the command line.
    if len(sys.argv) != 2:
        print('Usage: extract_time_series_data_from_netcdf_files.py `path/to/json/input/file\'')
        sys.exit()

    # Read and load the JSON file.
    input_file = sys.argv[1]
    with open(input_file) as f:
        inputs = json.load(f)

    # Process each dictionary (block) in the JSON file, one at a time, to produce an output file for each dictionary.
    for job_index in range(len(inputs)):
        start_time = time.time()
        extract_time_series_data_from_netcdf_files(**inputs[job_index])
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time for {inputs[job_index]['extracted_outputs_file']}: {elapsed_time:.2f} seconds")