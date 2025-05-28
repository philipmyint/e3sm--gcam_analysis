import uxarray as ux
import xarray as xr

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

def calculate_statistics_of_xarray(data, variable=None):
    """
    Calculates the min, mean, median, max, and the standard deviation of an xarray or uxarray object (DataArray, Dataset, or uxarray data structure).

    Parameters:
        data: Object of type xarray or uxarray whose statistical properties we want to calculate.
        variable: Variable of interest, used in the case of an xarray Dataset or uxarray data structure (UxDataArray or UxDataset).

    Returns:
        min, mean, median, max, and standard deviation of the xarray or uxarray object.
    """
    if isinstance(data, xr.DataArray):
        min = data.min().item()
        mean = data.mean().item()
        median = data.median().item()
        max = data.max().item()
        std = data.std().item()
    else:
        min = float(data[variable].min())
        mean = float(data[variable].mean())
        median = float(data[variable].median())
        max = float(data[variable].max())
        std = float(data[variable].std())
    return min, mean, median, max, std

def convert_xarray_to_uxarray(data, grid, variable=None, fillna=1):
    """
    Converts an xarray object (DataArray or Dataset) into an uxarray object (UxDataArray or UxDataset).

    Parameters:
        data: Object of type xarray.
        grid: Grid object for the unstructured grid contained within the uxarray object.
        variable: Variable of interest.
        fill_value: Value to fill NaNs.

    Returns:
        uxarray (UxDataArray or UxDataset) version of the given xarray object.
    """
    ds = xr.Dataset()
    if variable:
        ds[variable] = data
    else:
        ds = data
    if fillna:
        return ux.UxDataset.from_xarray(ds, grid).fillna(fillna)
    else:
        return ux.UxDataset.from_xarray(ds, grid)