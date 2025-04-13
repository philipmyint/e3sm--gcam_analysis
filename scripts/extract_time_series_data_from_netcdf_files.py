import json
import multiprocessing
import os
import pandas as pd
import sys
import time
import xarray as xr
from utility_functions import add_lists_elementwise, get_all_files_in_paths
from utility_dataframes import write_dataframe_to_fwf

def extract_year_and_month_from_name_of_netcdf_file(file):
    index = file.find('.nc')
    index -= 7
    year = int(file[index:index+4])
    month = int(file[index+5:index+7])
    return year, month

def get_netcdf_files_between_start_and_end_years(file_paths, start_year, end_year):
    new_file_paths = []
    for file in file_paths:
        year, _ = extract_year_and_month_from_name_of_netcdf_file(file)
        if year >= start_year and year <= end_year:
            new_file_paths.append(file)
    return new_file_paths

def put_column_means_of_netcdf_file_into_dataframe(file, variables=None, include_units_in_header=False):

    xrds = xr.open_dataset(file)
    if variables:
        xrds = xrds[variables]
    df = xrds.to_dataframe()

    if include_units_in_header:
        units = [xrds[var].attrs['units'] for var in xrds.data_vars]
        df.columns = add_lists_elementwise(xrds.data_vars, units, list2_are_units=True)

    df = df.mean().to_frame().T

    year, month = extract_year_and_month_from_name_of_netcdf_file(file)
    column_names_with_year_and_month_first = ['Year','Month']
    column_names_with_year_and_month_first.extend(df.columns)
    df['Year'] = year
    df['Month'] = month
    return df[column_names_with_year_and_month_first]

def extract_time_series_data_from_netcdf_files(output_paths, extracted_output_file, output_file_name_substrings=None, 
            start_year=2015, end_year=2100, variables=None, include_units_in_header=False):

    output_files = get_all_files_in_paths(output_paths, 
                    file_name_substrings=output_file_name_substrings, file_extension='.nc')
    output_files = get_netcdf_files_between_start_and_end_years(output_files, start_year, end_year)
    
    arguments = list(zip(output_files, [variables]*len(output_files), [include_units_in_header]*len(output_files)))
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        dataframes = list(pool.starmap(put_column_means_of_netcdf_file_into_dataframe, arguments))
    
    df_merged = dataframes[0]
    for df in dataframes[1:]:
        df_merged = pd.concat([df_merged,df]).reset_index(drop=True)
    df_merged.sort_values(['Year', 'Month'], inplace=True)

    write_dataframe_to_fwf(extracted_output_file, df_merged)

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Usage: extract_time_series_data_from_netcdf_files.py `path/to/json/input/file\'')
        sys.exit()

    input_file = sys.argv[1]
    with open(input_file) as f:
        inputs = json.load(f)

    for job_index in range(len(inputs)):
        start_time = time.time()
        extract_time_series_data_from_netcdf_files(**inputs[job_index])
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time for {inputs[job_index]['extracted_output_file']}: {elapsed_time:.2f} seconds")