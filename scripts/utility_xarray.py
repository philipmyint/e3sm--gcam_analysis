from scipy import stats
import xarray as xr
import xarray.ufuncs as xrf

def calculate_mean_and_std_of_da_list(da_list, calculate_std=False):
    """
    Calculates the mean of a list of xarray DataArrays and optionally also the standard deviation of this list.

    Parameters:
        da_list: A list of DataArrays with the same dimensions.

    Returns:
        A DataArray containing the mean of the input DataArrays and optionally also the standard deviation of this list.
    """
    if not all(isinstance(da, xr.DataArray) for da in da_list):
        raise TypeError('All elements in the list must be xarray DataArrays.')

    if not all(da.dims == da_list[0].dims for da in da_list[1:]):
        raise ValueError('All DataArrays must have the same dimensions.')

    # Stack the DataArrays along a new dimension.
    stacked_da = xr.concat(da_list, dim='new_dim')

    # Calculate the mean along the new dimension and the standard deviation after reducing the dataArrays along all dimensions.
    mean_da = stacked_da.mean(dim='new_dim')
    if calculate_std:
        std = stacked_da.std().item()
        return mean_da, std
    else:
        return mean_da

def calculate_statistics_of_da(da, dim=None, skipna=True):
    """
    Calculates the min, mean, median, max, and the standard deviation of an xarray DataArray.

    Parameters:
        da: DataArray whose statistical properties we want to calculate.
        dim: Name of dimensions. If none, the reduction will be performed over all dimensions.
        skipna: If True, skip missing values (as marked by NaN). 

    Returns:
        min, mean, median, max, and standard deviation of the DataArray.
    """
    min = da.min(dim=dim, skipna=skipna).item()
    mean = da.mean(dim=dim, skipna=skipna).item()
    median = da.median(dim=dim, skipna=skipna).item()
    max = da.max(dim=dim, skipna=skipna).item()
    std = da.std(dim=dim, skipna=skipna).item()
    return min, mean, median, max, std
