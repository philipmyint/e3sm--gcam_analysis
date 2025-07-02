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

def find_gridcell_areas_in_netcdf_file(file, region=None):
    """ 
    Obtains the grid cell areas of all latitude/longitude coordinates in an E3SM-generated (EAM or ELM or EHC) NetCDF file for the given region.

    Parameters:
        file: NetCDF file.
        region: String for the region of interest. If not specified or not recognized, then there will be no restrictions on the lat/lon coordinates. 

    Returns:
        NumPy array containing the grid cell areas of all coordinates in units of m^2 and an xarray Dataset containing data from the file.
        Also returns the land and non-land (which is defined as ocean if in the case of EAM) fractions.
    """
    # Load the file into an xarray Dataset and put the grid cell areas into a NumPy array.
    ds = xr.open_dataset(file)

    # If a region has been specified, get the bounds on the lat/lon coordinates and apply these bounds to restrict the lat/lon coordinates.
    if region:
        bounds = get_regional_bounds(region)
        lon_bounds = bounds[:2]
        lat_bounds = bounds[2:]
        if 'elm.h0' in file:
            ds = ds.where((ds.lon >= lon_bounds[0]), drop=True)
            ds = ds.where((ds.lon <= lon_bounds[1]), drop=True)
            ds = ds.where((ds.lat >= lat_bounds[0]), drop=True)
            ds = ds.where((ds.lat <= lat_bounds[1]), drop=True)
        elif 'eam.h0' in file:
            ds.load()
            ds = ds.where((ds.lon >= lon_bounds[0]), drop=True)
            ds = ds.where((ds.lon <= lon_bounds[1]), drop=True)
            ds = ds.where((ds.lat >= lat_bounds[0]), drop=True)
            ds = ds.where((ds.lat <= lat_bounds[1]), drop=True)
        elif 'surfdata_iESM_dyn' in file:
            ds = ds.where((ds.LONGXY >= lon_bounds[0]), drop=True)
            ds = ds.where((ds.LONGXY <= lon_bounds[1]), drop=True)
            ds = ds.where((ds.LATIXY >= lat_bounds[0]), drop=True)
            ds = ds.where((ds.LATIXY <= lat_bounds[1]), drop=True)

    if 'elm.h0' in file:
        # If the NetCDF file is produced by the ELM model, multiply the areas by the land fraction, plus additional land and pft masks.
        areas = create_numpy_array_from_ds(ds, ['area'], [0])
        variables = ['landfrac', 'landmask', 'pftmask']
        landfrac, landmask, pftmask = create_numpy_array_from_ds(ds, variables, [0]*len(variables))
        areas *= landfrac*landmask*pftmask*km2_TO_m2
        landfrac *= landmask*pftmask
        non_landfrac = 1 - landfrac
    elif 'eam.h0' in file:
        # If the NetCDF file is produced by the EAM model, multiply the areas by the Earth surface area so that the areas are in units of m^2.
        areas = create_numpy_array_from_ds(ds, ['area'], [0])
        areas *= SURF_AREA
        # Find the land and ocean fractions.
        variables = ['LANDFRAC', 'OCNFRAC']
        landfrac, non_landfrac = create_numpy_array_from_ds(ds, variables, [0]*len(variables))
    elif 'surfdata_iESM_dyn' in file:
        # Case where the NetCDF file is the land surface data file produced dynamically by the E3SM human component (EHC) model during run time.
        areas = create_numpy_array_from_ds(ds, ['AREA'], [0])
        variables = ['LANDFRAC_PFT', 'PFTDATA_MASK']
        landfrac, pftmask = create_numpy_array_from_ds(ds, variables, [0]*len(variables))
        areas *= landfrac*pftmask
        landfrac *= pftmask
        non_landfrac = 1 - landfrac
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

def get_regional_bounds(region):
    """
    Finds and returns a NumPy array that indicates the bounds on the longitude and latitude for given a geographical region.
    Adapted from a script by Daniel Ricciuto: /home/ac.eva.sinha/Sinha-etal-2025/workflow/e3sm_diags/diags_lineplots_hist_ZATM.py.

    Parameters:
        region: String for the geographical region.

    Returns:
        NumPy array with four numbers, [-min_lon, max_lon, min_lat, max_lat]. If the selected region is not recognized, returns [0, 360, -90, 90].
    """
    if region == 'noam':  # North America
        bounds = [-170.25, -45.25, 9.75, 79.75]
    elif region == 'bona':  # Boreal North America
        bounds = [-170.25, -60.25, 49.75, 79.75]
    elif region == 'tena':  # Temperate North America
        bounds = [-125.25, -66.25, 30.25, 49.75]
    elif region == 'conus':  # Continental US (CONUS)
        bounds = [-125.25, -66.25, 23.25, 54.75]
    elif region == 'columbia':   # Columbia River watershed
        bounds = [-126, -108, 40.0, 55.0]
    elif region == 'ceam':  # Central America
        bounds = [-115.25, -80.25, 9.75, 30.25]
    elif region == 'soam':    # South America
        bounds = [-80.25, -40.25, -59.75, 12.75]
    elif region == 'nhsa':  # Northern hemisphere South America
        bounds = [-80.25, -50.25, 0.25, 12.75]
    elif region == 'shsa':  # Southern hemisphere South America
        bounds = [-80.25, -40.25, -59.75, 0.25]
    elif region == 'amazon':  # Amazon
        bounds = [-85., -35., -25., 15.]
    elif region == 'euro':  # Europe
        bounds = [-10.25, 30.25, 35.25, 70.25]
    elif region == 'mide':  # Middle East
        bounds = [-10.25, 60.25, 20.24, 40.25]
    elif region == 'afrc':    # Africa, note can't span the Prime Meridian
     	bounds = [0.25, 45.25, -34.75, 20.25] 
    elif region == 'nhaf':  # Northern hemisphere Africa
        bounds = [-20.25, 45.25, 0.25, 20.25]
    elif region == 'shaf':  # Southern hemisphere Africa
        bounds = [10.25, 45.25, -34.75 ,0.25]
    elif region == 'asia':  # Asia
        bounds = [30.25, 179.75, -10.25, 70.25]
    elif region == 'boas':  # Boreal Asia
        bounds = [30.25, 179.75, 54.75, 70.25]
    elif region == 'ceas':  # Central Asia
        bounds = [30.25, 142.58, 30.25, 54.75]
    elif region == 'seas':  # Southeast Asia
        bounds = [65.25, 120.25, 5.25, 30.25]
    elif region == 'eqas':  # Equatorial Asia
        bounds = [99.75, 150.25, -10.25, 10.25]
    elif region == 'aust':  # Australia
        bounds = [112.00, 154.00, -41.25, -10.50]
    else:
        print(f'Did not recognize the selected region {region}! Setting bounds to global.\n')
        bounds = [-180, 180, -90, 90]   # Global
    bounds = np.array(bounds) 
    # Convert the longitudinal bounds to be between 0 and 360.
    bounds[:2] += 180
    return bounds