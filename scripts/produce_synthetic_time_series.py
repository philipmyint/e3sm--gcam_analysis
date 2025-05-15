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

def produce_synthetic_time_series(inputs):
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
    df = read_file_into_dataframe(file)
    columns = [column for column in df.columns if column not in ['Year', 'Month']]
    base_multipliers = np.linspace(1.02, 1.05, num_synthetic_sets_in_ensemble)
    random_multipliers = np.random.uniform(low=-0.02, high=0.02, size=len(df))
    for index in range(len(base_multipliers)):
        df_new = df
        multipliers = base_multipliers[index] + random_multipliers
        df_new[columns] = df[columns].multiply(multipliers, axis='index')
        if file.endswith('.csv'):
            new_file = file.replace('.csv', f'_{index+1}.csv')
            df_new.to_csv(new_file, index=False)
        else:
            if '.' in file:
                new_file = file.replace('.', f'_{index+1}.')
            else:
                new_file = file + f'_{index+1}'
            write_dataframe_to_fwf(new_file, df_new)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing the synthetic time series data for {file}: {elapsed_time:.2f} seconds")

###---------------Begin execution---------------###
if __name__ == '__main__':

    start_time = time.time()

    files = ["./../2025_DiVittorio_et_al/control_time_series.csv", "./../2025_DiVittorio_et_al/full_feedback_time_series.csv",
        "./../2025_DiVittorio_et_al/ag_scaling_time_series.csv", "./../2025_DiVittorio_et_al/carbon_scaling_time_series.csv"]
    num_sets_in_ensemble = [5]*len(files)
    inputs = list(zip(files, num_sets_in_ensemble))

    # Process each dictionary to produce a list of smaller dictionaries, each of which specifies data extraction options for a single output file.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(produce_synthetic_time_series, inputs)
    
    # Print the total execution time needed to complete all data extraction operations.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all synthetic time series data: {elapsed_time:.2f} seconds")