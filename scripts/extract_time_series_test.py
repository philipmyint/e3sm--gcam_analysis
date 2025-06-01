import json
import multiprocessing
import numpy as np
import pandas as pd
import sys
import time
from utility_constants import *
from utility_dataframes import move_columns_next_to_each_other_in_dataframe, write_dataframe_to_fwf
from utility_functions import *
from utility_e3sm_netcdf import *

def process_dataframe(df):
    """ 
    Processes a Pandas DataFrame by adding new columns (e.g., total precipitation, mole fraction CO2) or changing the units (e.g., Pg instead of g). 

    Parameters:
        df: DataFrame to be processed.

    Returns:
        The processed DataFrame.
    """
    # The number of years will later be used in NumPy tiling operations to produce an array with length equal to the number of rows in the DataFrame.
    start_year = df['Year'].min()
    end_year = df['Year'].max()
    num_years = end_year - start_year + 1

    # Convert precipitation variables from m/s to mm/month, add a new column for the total precipitation, and update the column labels.
    substrings = ['PRECC', 'PRECL', 'PRECSC', 'PRECSL']
    if check_substrings_in_list(substrings, df.columns, all_or_any='all'):
        columns_to_modify = [label for label in df.columns for substring in substrings if substring in label]
        seconds_in_months_tiled = np.tile(NUM_SECONDS_IN_MONTHS, num_years).reshape(-1,1)
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
    
    # Convert fluxes from per second to per month.
    old_label = '/s)'
    columns_to_modify = [label for label in df.columns if old_label in label]
    seconds_in_months_tiled = np.tile(NUM_SECONDS_IN_MONTHS, num_years).reshape(-1,1)
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
    Extracts the specified variables from an E3SM-generated NetCDF file into a Pandas DataFrame and performs the indicated operation on the variables.
    Each NetCDF file contains simulation results for one month in a particular year from either EAM (atmosphere model) or ELM (land model). 
    One type of operation could be to perform an area-weighted mean over the latitude/longitude coordinates for the given month.

    Parameters:
        file: Complete path and name of the NetCDF file.
        variables: List of variables that we want to extract from the NetCDF file.
        calculation_type: String that indicates the operation type.

    Returns:
        DataFrame containing one column for each of the variables, plus year and month columns.
    """
    # Extract the area variables from the NetCDF file and transform the variables from an xarray Dataset into a DataFrame.
    areas, ds, landfrac, non_landfrac = find_gridcell_areas_in_netcdf_file(file)
    ds = ds[variables]
    df = ds.to_dataframe()

    # Add units to the DataFrame column header.
    units = [ds[variable].attrs['units'] for variable in ds.data_vars]
    df.columns = add_lists_elementwise(df.columns, units, list2_are_units=True)
    
    # 9 land units in order: vegetation, crop, ice, multiple ice, lake, wetland, urban tbd, urban hd, urban md. We want only vegetation (index = 0).
    if 'PCT_LANDUNIT' in variables:
        df = df.reset_index(level='ltype')
        df = df[df['ltype'] == 0].drop(columns=['ltype'])
        # Divide the percent by 100 to change it to a fraction, and update the column name to reflect that it is a vegetation fraction.
        df['PCT_LANDUNIT (%)'] /= 100
        df = df.rename(columns={'PCT_LANDUNIT (%)': 'FRAC_VEG'})

    # 17 plant-functional types (PFTs) in order: 1 bare, 8 tree, 3 shrub, 3 grass, 1 crop, 1 empty. These can be further subgrouped as follows:
    # Bare soil (index 0), forest (the 8 trees, indices 1--8); shrub (indices 9--11); grass (indices 12--14), crop (index 15). Ignore the empty PFT.
    if 'PCT_LANDUNIT' in variables and 'PCT_NAT_PFT' in variables:
        df = df.reset_index(level='natpft')
        pft_labels = ['BARE_AREA (km^2)', 'FOREST_AREA (km^2)', 'SHRUB_AREA (km^2)', 'GRASS_AREA (km^2)', 'CROP_AREA (km^2)']
        pft_min_max_indices = [(0, 0), (1, 8), (9, 11), (12, 14), (15, 15)]
        for index, pft_label in enumerate(pft_labels):
            # Select only the rows that pertain to this particular PFT subgroup and for each lat/lon coordinate, sum over all PFTS in the subgroup.
            pft_min_index, pft_max_index = pft_min_max_indices[index][0], pft_min_max_indices[index][1]
            df_this_pft = df[(df['natpft'] >= pft_min_index) & (df['natpft'] <= pft_max_index)]
            df_this_pft = df_this_pft.groupby(['lat', 'lon']).sum()
            # Divide all columns except for column containing the percent PFT by the total number of PFTs in the subgroup.
            columns_to_except_pctpft = [col for col in df.columns if col != 'PCT_NAT_PFT (%)']
            df_this_pft[columns_to_except_pctpft] /= (pft_max_index - pft_min_index + 1)
            # Divide the percent by 100 to change it to a fraction, and update the column so that it contains the area occupied by the PFT subgroup.
            df_this_pft['PCT_NAT_PFT (%)'] /= 100
            df_this_pft = df_this_pft.rename(columns={'PCT_NAT_PFT (%)': pft_label})
            df_this_pft['Area (km^2)'] = areas/km2_TO_m2
            df_this_pft[pft_label] = df_this_pft['Area (km^2)']*df_this_pft['FRAC_VEG']*df_this_pft[pft_label]
            # Drop these columns since we are now doing with processing.
            df_this_pft = df_this_pft.drop(columns=['natpft', 'Area (km^2)'])
            # Use the entire DataFrame for the first subgroup. For all other subgroups, add just the column for the PFT subgroup area to the DataFrame.
            if index == 0:
                df_with_columns_to_keep = df_this_pft.copy()
            else:
                df_with_columns_to_keep[pft_label] = df_this_pft[pft_label]
        df = df_with_columns_to_keep

    print(df.head(10))

    if calculation_type == 'area_weighted_mean_or_sum':
        # Calculate an area-weighted mean or sum over all latitude/longitude coordinates for each variable.

        # First, multiply all variables in the DataFrame by the grid cell areas at each latitude/longitude coordinate.
        df *= areas
        # For EAM variables that correspond specifically to land or non-land (ocean) quantities, multiply by the land or non-land fractions.
        columns_to_modify = [label for label in df.columns if '_LND' in label]
        df[columns_to_modify] *= landfrac
        columns_to_modify = [label for label in df.columns if '_OCN' in label]
        df[columns_to_modify] *= non_landfrac

        # Sum over all latitude/longitude coordinates to get an area-weighted sum for each variable.
        df = df.sum().to_frame().T

        # Variables that are not fluxes and stocks (so that they are not per-area quantities), should be global area-weighted means
        # rather than area-weighted sums, and therefore we need to divide these sums (which were computed a few lines above) by the total area.
        columns_to_modify = [label for label in df.columns if not check_substrings_in_string(['/m^2', '/m2'], label, all_or_any='any')]
        # LND and OCN refer to certain types of EAM output for which we need to multiply the total grid cell areas by the land or non-land fractions.
        columns_to_modify_LND = [label for label in columns_to_modify if '_LND' in label]
        total_area = np.sum(areas*landfrac)
        df[columns_to_modify_LND]/= total_area
        columns_to_modify_OCN = [label for label in columns_to_modify if '_OCN' in label]
        total_area = np.sum(areas*non_landfrac)
        df[columns_to_modify_OCN]/= total_area
        # Variables that are not LND or OCN refer to either ELM output or to EAM output where we need to multiply by the full area of each grid cell.
        columns_to_modify = [label for label in columns_to_modify if (label not in columns_to_modify_LND and label not in columns_to_modify_OCN)]
        total_area = np.sum(areas)
        df[columns_to_modify]/= total_area

        # Fluxes and stocks are global area-weighted sums and have been multiplied by areas, so we must update the labels to remove the '/m^2' part.
        for old_label in ['/m^2', '/m2']:
            condition = lambda x: old_label in x
            new_column_label_function = lambda x: x.replace(old_label, '')
            df.columns = modify_list_based_on_condition(df.columns, condition, new_column_label_function)    

    elif calculation_type == 'sum':
        # Perform a sum over all latitude/longitude coordinates for each variable.
        df = df.sum().to_frame().T
    elif calculation_type == 'mean':
        # Calculate a mean over all latitude/longitude coordinates for each variable.
        df = df.mean().to_frame().T

    # Add year and month columns.
    year, month = extract_year_and_month_from_name_of_netcdf_file(file)
    column_names_with_year_and_month_first = ['Year', 'Month']
    column_names_with_year_and_month_first.extend(df.columns)
    df['Year'] = year
    df['Month'] = month
    return df[column_names_with_year_and_month_first]

def extract_time_series_from_netcdf_files(simulation_path, output_file, netcdf_substrings, 
            variables, calculation_types, process_variables=True, start_year=2015, end_year=2100, write_to_csv=False):
    """ 
    Extracts time series data from E3SM-generated NetCDF files in a simulation directory into a Pandas DataFrame and writes it to an output file. 
    The NetCDF files can be of more than one type (e.g., one set generated from the ELM model and another set from the EAM model in E3SM).
    Each NetCDF file contains simulation results for one month in a particular year.

    Parameters:
        simulation_path: Complete path of the directory containing the NetCDF output files from running a simulation.
        output_file: Name of the output file where the contents of the DataFrame will be written.
        netcdf_substrings: List where each element is itself a list of substrings. Each list corresponds to a particular type of NetCDF file and
                           indicates the substrings that must be in the names of that NetCDF file type.
        variables: List where each element is itself a list of variables we want to extract for each type of NetCDF file. 
                   The aggregate of all variables contained in these lists will together form the columns of the DataFrame.
        calculation_types: List of strings that indicate what type of calculation we want to perform to produce the time series data 
                           for each type of NetCDF file. Current options include 'area_weighted_mean_or_sum', 'mean', and 'sum'.
        process_variables: Boolean indicating if further processing is to be done on the outputs in the DataFrame. This includes adding new variables
                           not in the NetCDF files (e.g., total precipitation, mole fraction CO2) or changing the units (e.g., Pg instead of g).
        start_year: First year in the extracted time series data.
        end_year: Last year in the extracted time series data.

    Returns:
        N/A.
    """
    # This list will store a DataFrame for each type of NetCDF file, and all elements of this list will later be merged into a single DataFrame.
    dataframes = []

    # Iterate over all NetCDF file types.
    for index in range(len(variables)):

        # Get all NetCDF files for this particular type that fall within the start and end years.
        netcdf_files = get_all_files_in_path(simulation_path, file_name_substrings=netcdf_substrings[index], file_extension='.nc')
        netcdf_files = get_netcdf_files_between_start_and_end_years(netcdf_files, start_year, end_year)
    
        # Use multiprocessing to extract data from all the files. Put the data from each file into DataFrame, store all such DataFrames in a list.
        arguments = list(zip(netcdf_files, [variables[index]]*len(netcdf_files), [calculation_types[index]]*len(netcdf_files)))
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            dataframes_for_each_nc_file = list(pool.starmap(extract_netcdf_file_into_dataframe, arguments))
        
        # Concatenate all DataFrames in the list together to form a single DataFrame for this NetCDF type. Sort by year and month.
        df = pd.concat(dataframes_for_each_nc_file)
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
    if write_to_csv or output_file.endswith('.csv'):
        df.to_csv(output_file, index=False)
    else:
        write_dataframe_to_fwf(output_file, df)


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file(s) on the command line.
    start_time_total = time.time()
    if len(sys.argv) < 2:
        print('Usage: plot_spatial_data.py `path/to/json/input/file(s)\'')
        sys.exit()

    # Read and load the JSON file(s) into a list of dictionaries.
    inputs = []
    for index in range(1, len(sys.argv)):
        input_file = sys.argv[index]
        with open(input_file) as f:
            inputs.extend(json.load(f))

    # Process each dictionary to produce a list of smaller dictionaries, each of which specifies data extraction options for a single output file.
    for index in range(len(inputs)):
        start_time = time.time()
        extract_time_series_from_netcdf_files(**inputs[index])
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time for {inputs[index]['output_file']}: {elapsed_time:.2f} seconds")
    
    # Print the total execution time needed to complete all data extraction operations.
    end_time = time.time()
    elapsed_time = end_time - start_time_total
    print(f"Elapsed time for extracting all time series data: {elapsed_time:.2f} seconds")