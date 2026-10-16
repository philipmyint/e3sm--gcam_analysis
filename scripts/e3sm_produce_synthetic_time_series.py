import multiprocessing
import numpy as np
import sys
import time
from utility_constants import MAX_PROCESSES
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_fwf
from utility_functions import validate_max_processes

def produce_synthetic_time_series(inputs):
    """ 
    Produces a synthetic ensemble set of time series using random numbers to introduce perturbations to the time series in a given file. 
    Each new synthetic time series in the ensemble is output as a .dat or .csv file.

    Parameters:
        inputs: List with two items. The first item is the name of the file containing the base time series. 
                The second item is the total number of time series (including the base time series) we want to include in the set.

    Returns:
        N/A.
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
            new_file = file.replace('.csv', f'_{index+2}.csv')
            df_new.to_csv(new_file, index=False)
        else:
            if '.dat' in file:
                new_file = file.replace('.dat', f'_{index+2}.dat')
            else:
                new_file = file + f'_{index+2}'
            write_dataframe_to_fwf(new_file, df_new)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing the synthetic time series data for {file}: {elapsed_time:.2f} seconds")


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
            print('Usage: python e3sm_produce_synthetic_time_series.py [--max-processes N]')
            sys.exit()
        else:
            print(f"Error: unknown option: '{arg}'")
            print('Usage: python e3sm_produce_synthetic_time_series.py [--max-processes N]')
            sys.exit(1)

    files = ["./../2025_DiVittorio_et_al_e3sm/control_time_series.dat", 
             "./../2025_DiVittorio_et_al_e3sm/full_feedback_time_series.dat", 
             "./../2025_DiVittorio_et_al_e3sm/ag_scaling_time_series.dat", 
             "./../2025_DiVittorio_et_al_e3sm/carbon_scaling_time_series.dat",
             "./../2025_DiVittorio_et_al_e3sm/control_time_series_amazon.dat"]
    files.extend(["./../2025_DiVittorio_et_al_e3sm/control_time_series_surfdata_iESM_dyn_20240730.dat",
        "./../2025_DiVittorio_et_al_e3sm/control_time_series_surfdata_iESM_dyn_20240730_amazon.dat", 
        "./../2025_DiVittorio_et_al_e3sm/full_feedback_time_series_surfdata_iESM_dyn_20240730.dat",
        "./../2025_DiVittorio_et_al_e3sm/control_time_series_surfdata_iESM_dyn_20240820.dat", 
        "./../2025_DiVittorio_et_al_e3sm/full_feedback_time_series_surfdata_iESM_dyn_20240820.dat"])
    num_files_in_each_set = [5]*len(files)
    inputs = list(zip(files, num_files_in_each_set))

    num_processes = min(max_processes, max(1, len(inputs)))

    # Produce each set in parallel.
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(produce_synthetic_time_series, inputs)
    
    # Print the total execution time needed to produce all the sets.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all synthetic time series data: {elapsed_time:.2f} seconds")
