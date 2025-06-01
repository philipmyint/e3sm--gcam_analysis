import xarray as xr
from utility_constants import *
from utility_functions import create_numpy_array_from_ds

def extract_year_and_month_from_name_of_netcdf_file(file):
    """ 
    Finds the year and month from the name of a given E3SM-generated NetCDF file.

    Parameters:
        file: NetCDF file containing data for one month in a particular year.

    Returns:
        The year and month indicated in the name of the NetCDF file.
    """
    # Files are named *YYYY-MM.nc, so the year starts 7 characters before the .nc extension.
    index = file.find('.nc')
    index -= len('YYYY-MM')
    year = int(file[index:index+len('YYYY')])
    month = int(file[index+len('YYYY-'):index+len('YYYY-MM')])
    return year, month

def find_gridcell_areas_in_netcdf_file(file):
    """ 
    Obtains the grid cell areas of all latitude/longitude coordinates in an E3SM-generated (either EAM or ELM) NetCDF file.

    Parameters:
        file: NetCDF file.

    Returns:
        NumPy array containing the grid cell areas of all coordinates in units of m^2 and an xarray Dataset containing data from the file.
        Also returns the land and non-land (which is defined as ocean if in the case of EAM) fractions.
    """
    # Load the file into an xarray Dataset and put the grid cell areas into a NumPy array.
    ds = xr.open_dataset(file)
    areas = create_numpy_array_from_ds(ds, ['area'], [0])
    if 'elm.h0' in file:
        # If the NetCDF file is produced by the ELM model, multiply the areas by the land fraction, plus additional land and pft masks.
        variables = ['landfrac', 'landmask', 'pftmask']
        landfrac, landmask, pftmask = create_numpy_array_from_ds(ds, variables, [0]*len(variables))
        areas *= landfrac*landmask*pftmask*km2_TO_m2
        landfrac *= landmask*pftmask
        non_landfrac = 1 - landfrac
    elif 'eam.h0' in file:
        # If the NetCDF file is produced by the EAM model, multiply the areas by the Earth surface area so that the areas are in units of m^2.
        areas *= SURF_AREA
        # Find the land and ocean fractions.
        variables = ['LANDFRAC', 'OCNFRAC']
        landfrac, non_landfrac = create_numpy_array_from_ds(ds, variables, [0]*len(variables))
    return areas, ds, landfrac, non_landfrac

def get_netcdf_files_between_start_and_end_years(files, start_year, end_year):
    """ 
    Finds all E3SM-generated NetCDF files in a given list that fall between the start and end years.

    Parameters:
        files: List of NetCDF files.
        start_year: Start year.
        end_year: End year.

    Returns:
        List containing the NetCDF files that fall between the start and end years.
    """
    new_file_paths = []
    for file in files:
        year, _ = extract_year_and_month_from_name_of_netcdf_file(file)
        if year >= start_year and year <= end_year:
            new_file_paths.append(file)
    return new_file_paths