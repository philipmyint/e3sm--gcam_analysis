import multiprocessing
import numpy as np
import pandas as pd
import time
from utility_dataframes import read_file_into_dataframe, write_dataframe_to_fwf

def produce_synthetic_time_series(inputs):
    """ 
    Produces a synthetic set of time series using random numbers to introduce perturbations to the time series in a given file. Each new
    synthetic time series is output as a .dat or .csv file.

    Parameters:
        inputs: List with five items. The first item is the name of the file containing the base time series. 
                The second item is a list of the scenarios in the file.
                The third item contains the label for the scenario column in the file.
                The fourth item is a list of the column labels representing the numerical data for each scenario.
                The fifth item is the total number of time series (including the base time series) we want to include in the set.

    Returns:
        N/A.
    """
    start_time = time.time()
    file = inputs[0]
    scenarios = inputs[1]
    scenario_label = inputs[2]
    columns_to_modify = inputs[3]
    num_synthetic_sets_in_ensemble = inputs[4] - 1
    df = read_file_into_dataframe(file)
    df_full = df.copy()
    base_multipliers = np.linspace(1.02, 1.05, num_synthetic_sets_in_ensemble)
    for scenario in scenarios:
        df_this_scenario = df[df[scenario_label] == scenario]
        for set_index in range(num_synthetic_sets_in_ensemble):
            df_new = df_this_scenario.copy()
            random_multipliers = np.random.uniform(low=-0.02, high=0.02, size=len(df_this_scenario))
            multipliers = base_multipliers[set_index] + random_multipliers
            df_new[columns_to_modify] = df_this_scenario[columns_to_modify].multiply(multipliers, axis='index')
            df_new[scenario_label] = f'{scenario}_{set_index+2}'
            df_full = pd.concat([df_full, df_new], axis=0)
    if file.endswith('.csv'):
        new_file = file.replace('.csv', f'_full.csv')
        df_full.to_csv(new_file, index=False)
    else:
        if '.dat' in file:
            new_file = file.replace('.dat', f'_full.dat')
        else:
            new_file = file + f'_full'
        write_dataframe_to_fwf(new_file, df_full)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing the synthetic data for {file}: {elapsed_time:.2f} seconds")


###---------------Begin execution---------------###
if __name__ == '__main__':

    # The ensemble will consist of a total of len(files)*num_variations_for_each_scenario data sets.
    start_time = time.time()

    # Uncomment these lines to make variations on the scenarios in the files listed below.
    '''
    files = ["./../2025_DiVittorio_et_al/gcam/ag_commodity_prices.csv", 
             "./../2025_DiVittorio_et_al/gcam/co2_emissions_regions.csv", 
             "./../2025_DiVittorio_et_al/gcam/co2_emissions_sectors.csv"]
    columns_to_modify = [['value'], ['value'], ['value']]
    scenario_labels = ['scenario', 'scenario', 'scenario']
    '''
    # Make variations on the scenarios for the vegetation and soil scalars.
    files = ["./../2025_DiVittorio_et_al/gcam/scalars/scalars_multiple_scenarios.csv"]
    columns_to_modify = [['Vegetation', 'Soil']]
    scenario_labels = ['Scenario']

    # Create inputs: list of scenarios for each file and scenario label for each file, columns with the numerical data (to later produce variations).
    all_scenarios = []
    for index, file in enumerate(files):
        df = read_file_into_dataframe(file)
        all_scenarios.append(df[scenario_labels[index]].unique())
    num_variations_for_each_scenario = [5]*len(files)
    inputs = list(zip(files, all_scenarios, scenario_labels, columns_to_modify, num_variations_for_each_scenario))

    # Produce data for each file in parallel.
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(produce_synthetic_time_series, inputs)
    
    # Print the total execution time needed to produce all the sets.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all synthetic data: {elapsed_time:.2f} seconds")