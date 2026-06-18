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

def calculate_statistics_of_xarray(data, variable=None, weights=None):
    """
    Calculates the min, mean, median, max, and the standard deviation of an xarray or uxarray object (DataArray, Dataset, or uxarray data structure).

    Parameters:
        data: Object of type xarray or uxarray whose statistical properties we want to calculate.
        variable: Variable of interest, used in the case of an xarray Dataset or uxarray data structure (UxDataArray or UxDataset).
        weights: Optional xarray DataArray of weights to use when calculating the mean.

    Returns:
        min, mean, median, max, and standard deviation of the xarray or uxarray object.
    """
    if isinstance(data, xr.DataArray):
        min = data.min().item()
        if weights is not None:
            aligned_data, aligned_weights = xr.align(data, weights, join='inner')
            valid_weights = aligned_weights.where(aligned_data.notnull())
            weight_sum = valid_weights.sum()
            if weight_sum.item() == 0:
                mean = aligned_data.mean().item()
            else:
                mean = ((aligned_data * valid_weights).sum()/weight_sum).item()
        else:
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
    # Apply fillna on the xarray Dataset before constructing the UxDataset, to avoid the
    # uxarray/xarray >= 2025.7.1 incompatibility. Both ux.UxDataset.from_xarray and
    # ux.UxDataset.fillna internally reconstruct a UxDataset by passing a Dataset as the
    # first positional argument, which xarray >= 2025.7.1 no longer supports.
    # Applying fillna on the xr.Dataset first and passing dict(ds) to UxDataset avoids both issues.
    if fillna:
        ds = ds.fillna(fillna)
    # Rename the spatial dimension to match uxarray's face dimension name ('n_face') if needed.
    # When ux.open_dataset was used, this mapping was done automatically. With direct UxDataset
    # construction, we must do it manually — otherwise uxarray raises 'Data variable must be face
    # centered' at plot time because the data dimension (e.g. 'ncol') doesn't match 'n_face'.
    if hasattr(grid, 'n_face'):
        for dim in list(ds.dims):
            if dim not in ('time', 'year', 'lev', 'ilev') and ds.sizes[dim] == grid.n_face:
                ds = ds.rename({dim: 'n_face'})
                break
    return ux.UxDataset(dict(ds), uxgrid=grid)