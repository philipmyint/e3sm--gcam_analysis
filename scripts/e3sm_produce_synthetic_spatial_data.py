import multiprocessing
import numpy as np
import sys
import time
import xarray as xr
from utility_constants import MAX_PROCESSES
from utility_functions import validate_max_processes

def produce_synthetic_spatial_data(inputs):
    """ 
    Produces a synthetic ensemble set of spatial data using random numbers to introduce perturbations to the time series in a given file. 
    Each new synthetic spatial data map in the ensemble is output as a NetCDF file.

    Parameters:
        inputs: List with two items. The first item is the name of the file containing the base spatial data. 
                The second item is the total number of spatial data maps (including the base spatial data map) we want to include in the set.

    Returns:
        N/A.
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
        new_file = file.replace('.nc', f'_{index+2}.nc')
        ds.to_netcdf(new_file, mode='w')
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing the synthetic spatial data for {file}: {elapsed_time:.2f} seconds")


###---------------Begin execution---------------###
if __name__ == '__main__':

    # The ensemble will consist of a total of len(files)*num_files_in_each_set data sets.
    start_time = time.time()
    max_processes = MAX_PROCESSES
    args = sys.argv[1:]
    while args:
        arg = args.pop(0)
        if arg == '--max-processes':
            if not args:
                print("Error: --max-processes requires a positive integer")
                sys.exit(1)
            max_processes = validate_max_processes(args.pop(0))
        elif arg in ['-h', '--help']:
            print('Usage: python e3sm_produce_synthetic_spatial_data.py [--max-processes N]')
            sys.exit()
        else:
            print(f"Error: unknown option: '{arg}'")
            print('Usage: python e3sm_produce_synthetic_spatial_data.py [--max-processes N]')
            sys.exit(1)

    files = ["./../2025_DiVittorio_et_al_e3sm/control_spatial_data_elm.nc", 
             "./../2025_DiVittorio_et_al_e3sm/full_feedback_spatial_data_elm.nc", 
             "./../2025_DiVittorio_et_al_e3sm/carbon_scaling_spatial_data_elm.nc"]
    files.extend(["./../2025_DiVittorio_et_al_e3sm/control_spatial_data_eam.nc", 
                  "./../2025_DiVittorio_et_al_e3sm/full_feedback_spatial_data_eam.nc", 
                  "./../2025_DiVittorio_et_al_e3sm/carbon_scaling_spatial_data_eam.nc"])
    num_files_in_each_set = [5]*len(files)
    inputs = list(zip(files, num_files_in_each_set))

    num_processes = min(max_processes, max(1, len(inputs)))

    # Produce each set in parallel.
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(produce_synthetic_spatial_data, inputs)
    
    # Print the total execution time needed to produce all the sets.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all synthetic time series data: {elapsed_time:.2f} seconds")
