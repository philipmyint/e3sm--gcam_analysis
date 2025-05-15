import json
import multiprocessing
import numpy as np
import pandas as pd
import sys
import time
from utility_constants import *
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_fwf
from utility_functions import *
from utility_e3sm_netcdf import *

def produce_synthetic_spatial_data(inputs):
    """ 
    Processes a Pandas DataFrame by adding new columns (e.g., total precipitation, mole fraction CO2) or changing the units (e.g., Pg instead of g). 

    Parameters:
        df: DataFrame to be processed.

    Returns:
        The processed DataFrame.
    """
    start_time = time.time()
    file = inputs[0]
    num_synthetic_sets_in_ensemble = inputs[1] - 1
    random_multipliers = np.random.uniform(low=0.99, high=1.02, size=num_synthetic_sets_in_ensemble)
    for index in range(len(random_multipliers)):
        ds = xr.open_dataset(file)
        variables = ds.data_vars
        for variable in variables:   
            ds[variable] *= random_multipliers[index]
        new_file = file.replace('.nc', f'_{index+1}.nc')
        ds.to_netcdf(new_file, mode='w')
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing the synthetic spatial data for {file}: {elapsed_time:.2f} seconds")

###---------------Begin execution---------------###
if __name__ == '__main__':

    start_time = time.time()

    files = ["./../2025_DiVittorio_et_al/control_spatial_data_elm.nc", "./../2025_DiVittorio_et_al/full_feedback_spatial_data_elm.nc",
        "./../2025_DiVittorio_et_al/carbon_scaling_spatial_data_elm.nc"]
    num_sets_in_ensemble = [5]*len(files)
    inputs = list(zip(files, num_sets_in_ensemble))

    # Process each dictionary to produce a list of smaller dictionaries, each of which specifies data extraction options for a single output file.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(produce_synthetic_spatial_data, inputs)
    
    # Print the total execution time needed to complete all data extraction operations.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all synthetic time series data: {elapsed_time:.2f} seconds")