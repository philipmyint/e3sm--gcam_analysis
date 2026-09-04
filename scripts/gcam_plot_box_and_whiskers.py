import itertools
import json
import math
from matplotlib import pyplot as plt
import multiprocessing
import os
import pandas as pd
import seaborn as sns
import sys
import time
from utility_constants import *
from utility_dataframes import read_file_into_dataframe
from utility_functions import check_is_list_of_lists, transpose_scenarios_if_needed
from utility_gcam import gcam_landtype_groups, gcam_landtype_groups_original, produce_dataframe_for_landtype_group
from utility_plots import *

""" Dictionary of default input values for box plots (i.e., box-and-whisker plots). """
default_inputs = {
    'aggregation_function': 'sum',
    'basin_label': 'basin',
    'basins': None,
    'category_label': 'sector',
    'fill_boxes': True,
    'height': height_default,
    'hue': None,
    'key_columns': None,
    'landtype_groups': 'modified',
    'legend_label_size': legend_label_size_default,
    'legend_num_columns': 1,
    'legend_on': legend_on_default,
    'legend_place_outside': False, 
    'legend_x_offset': None,
    'linewidth': 1,
    'marker_size': 6,
    'multiplier': 1,
    'notify_scenarios_transposed': False,
    'plot_colors': plot_colors_default,
    'plot_directory': './',
    'plot_percent_difference': False,
    'plot_years': [2015, 2025, 2035, 2045],
    'plot_years_ncols': 2,
    'plot_type': 'ensemble',
    'produce_png': produce_png_default, 
    'region_label': 'region', 
    'regions': ['Global'],
    'scenario_label': 'scenario', 
    'scenario_sets': None,
    'use_latex': use_latex_default, 
    'value_label': 'value',
    'width': width_default,   
    'x_label': None,
    'x_label_size': axis_label_size_default,
    'x_scale': None,           
    'x_tick_labels': None,
    'x_tick_label_size': tick_label_size_default,   
    'y_label_size': axis_label_size_default,
    'y_limits': axis_limits_default,
    'y_scale': scale_default,
    'y_tick_label_size': tick_label_size_default,
    'year_label': 'year'            
}

def process_inputs(inputs):
    """ 
    Processes a dictionary of inputs (keys are options, values are choices for those options) for creating box (or box-and-whisker) plots.

    Parameters:
        inputs: Dictionary containing the user plotting choice inputs for different options. This dictionary may be incomplete or have invalid values.

    Returns:
        Dictionary that completely specifies all plotting options. 
        If the user did not make a choice for a particular option, the default choice for that plotting option will be selected.
    """
    # Check that the output file exists before trying to read it.
    if not os.path.exists(inputs['output_file']):
        print(f"Warning: output file not found, skipping: '{inputs['output_file']}'")
        return None
    if os.path.getsize(inputs['output_file']) == 0:
        print(f"Warning: output file is empty, skipping: '{inputs['output_file']}'")
        return None
    df = read_file_into_dataframe(inputs['output_file'])

    # If the category label (e.g., sector or landtype) has not been specified, use the default value.
    if 'category_label' not in inputs:
        category_label = default_inputs['category_label']
        inputs['category_label'] = category_label
    else:
        category_label = inputs['category_label']

    # If the list of categories (e.g., sectors or landtypes) have not been specified, populate the list with all categories except the excluded ones.
    if 'categories' not in inputs:
        if category_label in df.columns:
            inputs['categories'] = list(df[category_label].unique())
            if 'categories_to_exclude' in inputs:
                inputs['categories'] = [column for column in inputs['categories'] if column not in inputs['categories_to_exclude']]
        else:
            inputs['categories'] = 'All'
    categories = inputs['categories']
    # If the user entered a string indicating a single category, put that string in a list.
    if isinstance(categories, str):
        categories = [categories]
        inputs['categories'] = categories

    # Create the plot directory if it does not already exist.
    if 'plot_directory' not in inputs:
        inputs['plot_directory'] = default_inputs['plot_directory']
    if not os.path.exists(inputs['plot_directory']):
        os.makedirs(inputs['plot_directory'])

    # Use the name of the output file itself (without its path) to set defaults for the y-axis label and the name of the plot.
    index_of_last_backslash = inputs['output_file'].rfind('/')
    index_of_dot_csv = inputs['output_file'].find('.csv')
    if index_of_last_backslash == -1:
        output_file_name = inputs['output_file'][:index_of_dot_csv]
    else:
        output_file_name = inputs['output_file'][index_of_last_backslash+1:index_of_dot_csv]
    output_file_name = output_file_name.replace('_processed', '').replace('processed_', '')
    if 'y_label' not in inputs:
        inputs['y_label'] = output_file_name
    if 'plot_name' not in inputs:
        scenarios = inputs.get('scenarios', None)
        scenario_sets = inputs.get('scenario_sets', None)
        if scenarios is None:
            scenario_part = 'all_scenarios'
        elif check_is_list_of_lists(scenarios):
            num_scenario_sets = len(scenario_sets) if scenario_sets else len(scenarios[0])
            if num_scenario_sets == 1:
                scenario_part = str(scenario_sets[0] if scenario_sets else 'set0').replace(' ', '_')
            elif inputs.get('plot_percent_difference', False):
                scenario_part = 'pct_diff'
            else:
                scenario_part = 'comparison'
        elif len(scenarios) == 1:
            scenario_part = str(scenarios[0]).replace(' ', '_')
        elif inputs.get('plot_percent_difference', False):
            scenario_part = 'pct_diff'
        else:
            scenario_part = 'comparison'

        plot_years = inputs.get('plot_years', default_inputs['plot_years'])
        if plot_years is None:
            year_part = 'all_years'
        elif len(plot_years) == 1:
            year_part = str(plot_years[0])
        else:
            year_part = f'{min(plot_years)}-{max(plot_years)}'

        name = f'box_plot_{output_file_name}_{scenario_part}_{year_part}'
        inputs['plot_name'] = os.path.join(inputs['plot_directory'], name)
    elif 'plot_name' in inputs and '/' not in inputs['plot_name']:
        # If the user specified only a file name (no path) for the plot name, put the plot in the plot directory.
        inputs['plot_name'] = os.path.join(inputs['plot_directory'], inputs['plot_name'])

    # Add keys for plotting options that have not been specified in the inputs dictionary and use default values for them.
    for key in default_inputs.keys():
        if key not in inputs:
            inputs[key] = default_inputs[key]

    if inputs['hue'] == 'scenario' and inputs['plot_percent_difference']:
        print(f"Warning: skipping '{inputs['output_file']}' because hue='scenario' and "
              "plot_percent_difference=True cannot be used together.")
        return None

    # If the scenarios are specified as a list of lists (for ensemble plots), check if they need to be transposed.
    # Users can now specify scenarios in two formats:
    #   Format A (organized by ensemble member - original format):
    #       [["Control", "Full feedback"], ["Control_2", "Full feedback_2"], ...]
    #   Format B (organized by scenario set - new user-friendly format):
    #       [["Control", "Control_2", ...], ["Full feedback", "Full feedback_2", ...]]
    # The plotting functions expect Format A internally, so Format B will be automatically transposed.
    if check_is_list_of_lists(inputs['scenarios']):
        scenario_sets = inputs.get('scenario_sets', None)
        inputs['scenarios'], was_transposed = transpose_scenarios_if_needed(inputs['scenarios'], scenario_sets)
        if was_transposed and inputs['notify_scenarios_transposed']:
            print(f"Note: Scenarios were automatically transposed from 'organized by scenario set' format to 'organized by ensemble member' format.")

    # If the scenarios have not been specified, use all the scenarios in the Pandas DataFrame.
    if 'scenarios' not in inputs:
        scenario_label = inputs['scenario_label']
        inputs['scenarios'] = df[scenario_label].unique()
    
    # Check that the value_label column exists in the DataFrame.
    value_label = inputs['value_label']
    if value_label not in df.columns:
        raise KeyError(f"Error: value_label '{value_label}' not found in '{inputs['output_file']}'. "
                       f"Available columns: {list(df.columns)}")

    # If aggregation_function is area_weighted_mean, check that the 'area' column exists.
    if inputs['aggregation_function'] == 'area_weighted_mean' and 'area' not in df.columns:
        raise KeyError(f"Error: aggregation_function is 'area_weighted_mean' "
                       f"but no 'area' column found in '{inputs['output_file']}'. "
                       f"Consider using 'mean' or 'sum' instead, or add area data using gcam_add_areas_to_files.py.")

    # Check that hue and x_variable columns exist in the DataFrame if specified.
    hue = inputs.get('hue', None)
    if hue and hue not in df.columns:
        raise KeyError(f"Error: hue '{hue}' not found in '{inputs['output_file']}'. "
                       f"Available columns: {list(df.columns)}")
    x_variable = inputs.get('x_variable', inputs['category_label'])
    if x_variable not in df.columns:
        raise KeyError(f"Error: x_variable '{x_variable}' not found in '{inputs['output_file']}'. "
                       f"Available columns: {list(df.columns)}")

    # Warn if any specified scenario is not found in the DataFrame.
    scenario_label = inputs['scenario_label']
    if scenario_label in df.columns:
        available_scenarios = set(df[scenario_label].unique())
        specified_scenarios = inputs['scenarios']
        if check_is_list_of_lists(specified_scenarios):
            flat_scenarios = [s for member in specified_scenarios for s in member]
        else:
            flat_scenarios = list(specified_scenarios)
        missing_scenarios = [s for s in flat_scenarios if s not in available_scenarios]
        if missing_scenarios:
            print(f"Warning: the following scenarios were not found in '{inputs['output_file']}' "
                  f"and will produce empty plots: {missing_scenarios}")

    # Verify the plot directory is writable before doing any work.
    plot_name = inputs['plot_name']
    plot_dir = os.path.dirname(plot_name) or '.'
    if not os.access(plot_dir, os.W_OK):
        raise OSError(f"Error: plot directory is not writable: '{plot_dir}'")

    return inputs

def plot_box_and_whiskers(inputs):
    """ 
    Creates a box plot (box-and-whisker plot) and perform basic statistical analysis for a single output file. The data in the file are organized
    into scenarios or scenario sets, categories, and regions.
    Types of plots: 1) Plots in which each scenario is treated as an individual data set (not grouped into an ensemble) in the collection. 
    2) Ensemble plots in which the scenarios are grouped into sets (referred to as an ensemble) and results from of each set are plotted. 
    The scenarios must be specified as a list of lists (nested list) for ensemble plots.

    Parameters:
        input: Dictionary containing user inputs for different plotting options, where the keys are options and values are choices for those options.
               This dictionary is assumed to be complete (pre-processed).

    Returns:
        N/A.
    """
    # This function creates a box plot where the data (scenarios or scenario sets, categories, regions) all come from a single output file.
    start_time = time.time()
    output_file = inputs['output_file']
    
    # Extract all other plotting options.
    basin_label = inputs['basin_label']
    basins = inputs['basins']
    categories = inputs['categories']
    category_label = inputs['category_label']
    fill_boxes = inputs['fill_boxes']
    height = inputs['height']
    hue = inputs['hue']
    key_columns = inputs['key_columns']
    landtype_groups = inputs['landtype_groups']
    legend_label_size = inputs['legend_label_size']
    legend_num_columns = inputs['legend_num_columns']
    legend_on = inputs['legend_on']
    legend_place_outside = inputs['legend_place_outside']
    legend_x_offset = inputs['legend_x_offset']
    linewidth = inputs['linewidth']
    marker_size = inputs['marker_size']
    multiplier = inputs['multiplier']
    aggregation_function = inputs['aggregation_function']
    plot_colors = inputs['plot_colors']
    plot_name = inputs['plot_name']
    plot_percent_difference = inputs['plot_percent_difference']
    plot_years = inputs['plot_years']
    plot_type = inputs['plot_type']
    produce_png = inputs['produce_png']
    region_label = inputs['region_label']
    regions = inputs['regions']
    scenario_label = inputs['scenario_label']
    scenario_sets = inputs['scenario_sets']
    scenarios = inputs['scenarios']
    use_latex = inputs['use_latex']
    value_label = inputs['value_label']
    width = inputs['width']
    x_label = inputs['x_label']
    x_label_size = inputs['x_label_size']
    x_scale = inputs['x_scale']
    x_tick_labels = inputs['x_tick_labels'] or {}
    x_tick_label_size = inputs['x_tick_label_size']
    x_variable = inputs.get('x_variable', category_label)
    y_label = inputs['y_label']
    y_label_size = inputs['y_label_size']
    y_limits = inputs['y_limits']
    y_scale = inputs['y_scale']
    y_tick_label_size = inputs['y_tick_label_size']
    year_label = inputs['year_label']

    # Set the plotting options.
    plot_options = dict(width=width, height=height, name=plot_name, x_label=x_label, y_label=fr'{y_label}')
    plot_options.update(zip(['x_scale', 'y_scale', 'y_limits', 'use_latex'], [x_scale, y_scale, y_limits, use_latex]))
    plot_options.update(zip(['x_tick_label_size', 'y_tick_label_size', 'legend_on'], [x_tick_label_size, y_tick_label_size, legend_on]))
    plot_options.update(zip(['x_label_size', 'y_label_size', 'legend_label_size'], [x_label_size, y_label_size, legend_label_size]))
    plot_options.update(zip(['produce_png', 'legend_x_offset', 'legend_place_outside'], [produce_png, legend_x_offset, legend_place_outside]))
    plot_options.update(zip(['legend_num_columns'], [legend_num_columns]))

    # Use LaTeX fonts for figures and set font size of tick labels.
    setup_plot_params(plot_options)

    # Read the file, select requested years, and apply the user-specified multiplier.
    df = read_file_into_dataframe(output_file)
    if plot_years is not None:
        df = df[df[year_label].isin(plot_years)]
    df[value_label] *= multiplier

    # Create and concatenate Pandas DataFrames for the regions and basins; 'Global' involving creating a copy of the entire DataFrame.
    dataframes = []
    if 'Global' in regions:
        df_global = df.copy()
        df_global[region_label] = 'Global'
        dataframes.append(df_global) 
        regions.remove('Global')
    if regions:
        df_these_regions = df[df[region_label].isin(regions)]
        dataframes.append(df_these_regions) 
    if basins:
        df_these_basins = df[df[basin_label].isin(basins)]
        dataframes.append(df_these_basins) 
    df = pd.concat(dataframes)

    # Option 1: individual plots, in which each such box plot can include one or more individual (not grouped) data sets.
    if not check_is_list_of_lists(scenarios) or plot_type == 'individual':
        
        # If scenarios is a list of lists, but we do not actually want to group each inner list into an ensemble, 
        # unpack scenarios into a list of strings, with one string for each output file.
        if check_is_list_of_lists(scenarios) and plot_type == 'individual':
            scenarios = list(itertools.chain.from_iterable(scenarios))
            # Turn the corresponding plot_colors to list of lists with enough rows to match, then unpack them.
            plot_colors = list(itertools.chain.from_iterable(plot_colors))

        # Select the scenarios of interest.
        df = df[df[scenario_label].isin(scenarios)]
        scenario_list = scenarios
    
    # Option 2: ensemble plots, in which the outputs are subdivided into groups, where each group (ensemble) represents a set of scenarios.
    elif check_is_list_of_lists(scenarios) and plot_type == 'ensemble':
        # For each scenario set, collect all scenarios that belong to that set and put in a DataFrame, then concatenate the DataFrames for all sets.
        dataframes = []
        for index, scenario_set_name in enumerate(scenario_sets):
            scenarios_in_set = [scenarios[i][index] for i in range(len(scenarios))]
            df_this_set = df[df[scenario_label].isin(scenarios_in_set)].copy()
            # All scenarios in this set will be collected together and have the same label.
            df_this_set[scenario_label] = scenario_set_name
            dataframes.append(df_this_set)
        # Concatenate all DataFrames for all scenario sets by row. 
        # Averages over each scenario set will be automatically computed later when creating the box plot.
        df = pd.concat(dataframes)
        scenario_list = scenario_sets

    # The rest of the script involves creating DataFrames for different categories (e.g., landtypes or sectors) and producing the box plot.
    # This applies to both individual and ensemble plots.

    # Set the appropriate dictionary for the landtype group.
    if landtype_groups == 'modified':
        landtype_groups = gcam_landtype_groups
    elif landtype_groups == 'original':
        landtype_groups = gcam_landtype_groups_original

    # Create DataFrames for 1) all categories (create a copy of the entire DataFrame), 2) categories that correspond to
    # a group of landtypes (e.g., forest, crop, grass, shrub, pasture), and 3) a set of individual categories.
    dataframes = []
    if 'All' in categories:
        df_all = df.copy()
        df_all[category_label] = 'All'
        dataframes.append(df_all) 
        categories.remove('All')
    if any(item in landtype_groups for item in categories):
        for landtype_group in landtype_groups:
            if landtype_group in categories:
                if not df[category_label].isin(landtype_groups[landtype_group]).any():
                    continue
                df_this_landtype_group = produce_dataframe_for_landtype_group(df, landtype_group, category_label, 
                    value_label, landtype_groups, aggregation_function, key_columns)
                dataframes.append(df_this_landtype_group)
                categories.remove(landtype_group)
    if categories:
        df = df[df[category_label].isin(categories)]
        dataframes.append(df)
    df = pd.concat(dataframes).reset_index()

    units_by_category = {}
    if isinstance(category_label, str) and x_variable == category_label and 'units' in df.columns:
        for category, group in df.groupby(category_label):
            units = group['units'].dropna().unique()
            if len(units) == 1 and '/' in units[0]:
                units_by_category[str(category)] = units[0]
        numerators = {unit.split('/', 1)[0] for unit in units_by_category.values()}
        if units_by_category and len(numerators) == 1:
            numerator = numerators.pop()
            if numerator not in y_label:
                y_label += f' (% diff in {numerator})' if plot_percent_difference else f' ({numerator})'
                plot_options['y_label'] = fr'{y_label}'
    units = df['units'].dropna().unique() if 'units' in df.columns else []
    if not units_by_category and len(units) == 1 and units[0] not in y_label:
        unit_label = f'% diff in {units[0]}' if plot_percent_difference else units[0]
        y_label += f' ({unit_label})'
        plot_options['y_label'] = fr'{y_label}'

    years_to_plot = plot_years if plot_years is not None else sorted(df[year_label].unique())
    ncols = inputs['plot_years_ncols'] or len(years_to_plot)
    nrows = math.ceil(len(years_to_plot) / ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(width * ncols, height * nrows))
    axes_flat = list(axes.flat) if hasattr(axes, 'flat') else [axes]

    for panel_idx, year in enumerate(years_to_plot):
        ax = axes_flat[panel_idx]
        df_year = df[df[year_label] == year].copy()
        if df_year.empty:
            print(f"Warning: year {year} not found in data. Skipping panel.")
            ax.set_visible(False)
            continue
        if plot_percent_difference:
            reference = df_year[df_year[scenario_label] == scenario_list[0]][value_label].to_numpy()
            for scenario in scenario_list:
                scenario_filter = df_year[scenario_label] == scenario
                df_year.loc[scenario_filter, value_label] -= reference
                df_year.loc[scenario_filter, value_label] *= 100 / (reference + EPSILON)
            df_year = df_year[df_year[scenario_label] != scenario_list[0]]

        if df_year.empty or not df_year[value_label].notna().any():
            print(f"Warning: no plottable data for year {year}. Skipping panel.")
            ax.set_visible(False)
            continue

        num_x_variables = len(df_year[hue].unique()) if hue else len(df_year[x_variable].unique())
        palette = plot_colors[:num_x_variables] if not plot_percent_difference else plot_colors[1:num_x_variables + 1]
        box_properties = {'edgecolor': 'black'} if fill_boxes else {'color': 'black'}
        line_properties = {'color': 'black'}
        flier_properties = {'markeredgecolor': 'black', 'markerfacecolor': 'none'}
        if hue:
            sns.boxplot(df_year, x=x_variable, y=value_label, hue=hue, ax=ax, linewidth=linewidth,
                        palette=palette, fliersize=marker_size, fill=fill_boxes, boxprops=box_properties,
                        whiskerprops=line_properties, capprops=line_properties, medianprops=line_properties,
                        flierprops=flier_properties)
        else:
            sns.boxplot(df_year, x=x_variable, y=value_label, hue=x_variable, ax=ax, legend=False, linewidth=linewidth,
                        palette=palette, fliersize=marker_size, fill=fill_boxes, boxprops=box_properties,
                        whiskerprops=line_properties, capprops=line_properties, medianprops=line_properties,
                        flierprops=flier_properties)
        if units_by_category:
            ticks = ax.get_xticks()
            labels = [tick.get_text() for tick in ax.get_xticklabels()]
            ax.set_xticks(ticks, labels=[
                f"{x_tick_labels.get(label, label)} ({units_by_category[label].split('/', 1)[1]})"
                if label in units_by_category else x_tick_labels.get(label, label)
                for label in labels
            ])
        elif x_tick_labels:
            ticks = ax.get_xticks()
            labels = [tick.get_text() for tick in ax.get_xticklabels()]
            ax.set_xticks(ticks, labels=[x_tick_labels.get(label, label) for label in labels])
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
        ax.set_title(f'{chr(ord("a") + panel_idx)}) {year}')
        panel_options = plot_options.copy()
        panel_options['legend_on'] = legend_on if hue else False
        panel_options['y_label'] = ''
        set_figure_options(fig, ax, panel_options)

    for panel_idx in range(len(years_to_plot), nrows * ncols):
        axes_flat[panel_idx].set_visible(False)

    fig.set_size_inches(width * ncols, height * nrows)
    fig.supylabel(plot_options['y_label'], fontsize=y_label_size)
    fig.tight_layout(rect=(0.04, 0, 1, 1))
    plot_options['name'] = plot_name
    save_figure(plot_name, fig, plot_options)

    # Close the plot now that we are done with it. Record the elapsed time.
    plt.close(fig)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing plot {plot_name}: {elapsed_time:.2f} seconds")


###---------------Begin execution---------------###
if __name__ == '__main__':

    # Run this script together with the input JSON file(s) on the command line.
    start_time = time.time()
    if len(sys.argv) < 2:
        print('Usage: python gcam_plot_box_and_whiskers.py `path/to/json/input/file(s)\'')
        sys.exit()

    # Read and load the JSON file(s) into a list of dictionaries.
    inputs = []
    for i in range(1, len(sys.argv)):
        input_file = sys.argv[i]
        with open(input_file) as f:
            inputs.extend(json.load(f))

    # Process each dictionary so that each of them specifies a complete set of options (e.g., by adding default values) for a single plot.
    list_of_inputs = []
    for index in range(len(inputs)):
        result = process_inputs(inputs[index])
        if result is not None:
            list_of_inputs.append(result)

    # Create all of the box plots in parallel.
    with multiprocessing.Pool(processes=MAX_PROCESSES) as pool:
        pool.map(plot_box_and_whiskers, list_of_inputs)
    
    # Print the total execution time to produce all the plots.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all plots: {elapsed_time:.2f} seconds")