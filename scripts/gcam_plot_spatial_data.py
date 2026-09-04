import geopandas as gpd
import itertools
import json
import math
from matplotlib import pyplot as plt
import multiprocessing
import os
import pandas as pd
from scipy import stats
import sys
import time
from utility_constants import *
from utility_dataframes import perform_ttest, read_file_into_dataframe
from utility_functions import check_is_list_of_lists, print_p_values, sort_file, transpose_scenarios_if_needed
from utility_gcam import *
from utility_plots import *

""" Dictionary of default input values for spatial plots. """
default_inputs = {
    'aggregate_to_region': False,
    'basin_label': 'basin',
    'category_label': 'sector',
    'cbar_limits': None,
    'cbar_on': True,
    'cmap': 'viridis',
    'height': height_default,
    'key_columns': None,
    'landtype_groups': 'modified',
    'linewidth': 0.5,
    'multiplier': 1,
    'aggregation_function': 'sum',
    'notify_scenarios_transposed': False,
    'p_value_file': 'p_values.dat',
    'p_value_file_print_only_if_below_threshold': True,
    'p_value_threshold': 0.05,
    'plot_years': [2015, 2025, 2035, 2045],
    'plot_years_ncols': 2,
    'plot_categories': None,
    'plot_categories_ncols': 2,
    'plot_directory': './',
    'plot_type': 'mean',
    'produce_png': produce_png_default, 
    'region_label': 'region', 
    'scenario_label': 'scenario', 
    'scenario_sets': None,
    'shape_file_basin_label': 'glu_nm',
    'shape_file_region_label': 'reg_nm',
    'stippling_hatches': 'xxxx',
    'stippling_on': True,
    'units': None,
    'use_latex': use_latex_default, 
    'value_label': 'value',
    'width': width_default,   
    'x_tick_label_size': tick_label_size_default,   
    'y_tick_label_size': tick_label_size_default,
    'title_size': axis_label_size_default,
    'year_label': 'year'            
}

def _aggregate(df_rows, value_label, method):
    if method == 'sum':
        return df_rows[value_label].sum()
    if method == 'area_weighted_mean':
        total_area = df_rows['area'].sum()
        return (df_rows[value_label] * df_rows['area']).sum() / total_area if total_area != 0 else 0.0
    return df_rows[value_label].mean()


def _panel_units_str(df_panel, plot_type, default_units_str):
    """ Looks up the units for a single category panel, since different categories in the same file can have different units. """
    if plot_type not in ('absolute_difference', 'mean', 'percent_difference') or 'units' not in df_panel.columns:
        return default_units_str
    panel_units = df_panel['units'].dropna().unique()
    if len(panel_units) != 1:
        return default_units_str
    if plot_type == 'percent_difference':
        return f' (% diff in {panel_units[0]})'
    return f' ({panel_units[0]})'


def process_inputs(inputs):
    """ 
    Processes a dictionary of inputs (keys are options, values are choices for those options) for creating spatial plots.

    Parameters:
        inputs: Dictionary containing the user plotting choice inputs for different options. This dictionary may be incomplete or have invalid values.

    Returns:
        Dictionary that completely specifies all plotting options. 
        If the user did not select a plotting option fora particular category, the default choice for that plotting option will be selected.
    """
    # Check that the output file and shape file exist before trying to read them.
    if not os.path.exists(inputs['output_file']):
        print(f"Warning: output file not found, skipping: '{inputs['output_file']}'")
        return None
    if os.path.getsize(inputs['output_file']) == 0:
        print(f"Warning: output file is empty, skipping: '{inputs['output_file']}'")
        return None
    if not os.path.exists(inputs['shape_file']):
        raise FileNotFoundError(f"Error: shape file not found: '{inputs['shape_file']}'")
    df = read_file_into_dataframe(inputs['output_file'])

    # If the category label (e.g., sector or landtype) has not been specified, use the default value.
    if 'category_label' not in inputs:
        category_label = default_inputs['category_label']
        inputs['category_label'] = category_label
    else:
        category_label = inputs['category_label']

    # If the list of categories (e.g., sectors or landtypes) have not been specified, populate the list with all categories except the excluded ones.
    if 'categories' not in inputs:
        if isinstance(category_label, list) and all(c in df.columns for c in category_label):
            combos = df[category_label].drop_duplicates()
            inputs['categories'] = list(combos.itertuples(index=False, name=None))
            if 'categories_to_exclude' in inputs:
                exclude = {tuple(e) if isinstance(e, list) else e for e in inputs['categories_to_exclude']}
                inputs['categories'] = [c for c in inputs['categories'] if c not in exclude]
        elif not isinstance(category_label, list) and category_label in df.columns:
            if 'categories_to_exclude' in inputs:
                inputs['categories'] = [column for column in df[category_label].unique() if column not in inputs['categories_to_exclude']]
            else:
                # No categories were excluded, so this is equivalent to aggregating over 'All' categories; label it as such.
                inputs['categories'] = 'All'
        else:
            inputs['categories'] = 'All'
    categories = inputs['categories']
    # If categories is a dict keyed by column label, expand to the Cartesian product of the per-column value lists.
    if isinstance(categories, dict) and isinstance(category_label, list):
        per_column = [categories[col] for col in category_label]
        categories = list(itertools.product(*per_column))
        inputs['categories'] = categories
    elif isinstance(categories, str):
        categories = [categories]
        inputs['categories'] = categories
    elif isinstance(category_label, list):
        # Convert any list-type categories (from JSON) to tuples.
        categories = [tuple(c) if isinstance(c, list) else c for c in categories]
        inputs['categories'] = categories

    # If plot_categories is a dict keyed by column label, expand to the Cartesian product of the per-column value lists,
    # producing one panel per combination of values across the category_label entries.
    plot_categories = inputs.get('plot_categories')
    if isinstance(plot_categories, dict) and isinstance(category_label, list):
        per_column = [plot_categories[col] for col in category_label]
        plot_categories = list(itertools.product(*per_column))
        inputs['plot_categories'] = plot_categories

    # Panel categories must also be available in the prepared data. This is needed
    # when a panel requests an aggregate landtype group such as "crop".
    if isinstance(plot_categories, list):
        for category in plot_categories:
            normalized_category = tuple(category) if isinstance(category_label, list) and isinstance(category, list) else category
            if normalized_category not in categories:
                categories.append(normalized_category)
        inputs['categories'] = categories

    # If units have not been specified explicitly, fall back to the 'units' column in the data file, if present.
    # Filter to the selected categories first, since different categories in the same file can have different units.
    if isinstance(category_label, list):
        category_label_in_data = all(c in df.columns for c in category_label)
    else:
        category_label_in_data = category_label in df.columns
    if 'units' not in inputs and 'units' in df.columns:
        relevant_categories = plot_categories if plot_categories is not None else (categories if categories != 'All' else None)
        df_units_source = pd.DataFrame()
        if relevant_categories and category_label_in_data:
            if isinstance(category_label, list):
                mask = df[category_label].apply(tuple, axis=1).isin(set(relevant_categories))
            else:
                mask = df[category_label].isin(relevant_categories)
            df_units_source = df[mask]
        # Fall back to the whole file if none of the relevant categories matched (e.g., they are landtype group
        # names like 'crop' that are aggregated from raw values rather than literal column values).
        if df_units_source.empty:
            df_units_source = df
        data_units = df_units_source['units'].dropna().unique()
        if len(data_units) == 1:
            inputs['units'] = data_units[0]

    # Create the plot directory if it does not already exist. By default, put the name of the file containing p-value results in this directory.
    if 'plot_directory' not in inputs:
        inputs['plot_directory'] = default_inputs['plot_directory']
    if not os.path.exists(inputs['plot_directory']):
        os.makedirs(inputs['plot_directory'])
    if 'p_value_file' not in inputs:
        inputs['p_value_file'] = os.path.join(inputs['plot_directory'], default_inputs['p_value_file'])

    # Use the name of the output file itself (without its path) to set defaults for the title and the name of the plot.
    index_of_last_backslash = inputs['output_file'].rfind('/')
    index_of_dot_csv = inputs['output_file'].find('.csv')
    if index_of_last_backslash == -1:
        output_file_name = inputs['output_file'][:index_of_dot_csv]
    else:
        output_file_name = inputs['output_file'][index_of_last_backslash+1:index_of_dot_csv]
    if 'title' not in inputs:
        inputs['title'] = output_file_name
    if 'plot_name' not in inputs:
        # Build a descriptive default name from the key plot parameters; strip '_processed' from the file base.
        file_base = output_file_name.replace('_processed', '').replace('processed_', '')
        _plot_type_abbrevs = {'absolute_difference': 'abs_diff', 'percent_difference': 'pct_diff'}
        _raw_plot_type = inputs.get('plot_type', default_inputs['plot_type'])
        plot_type_str = _plot_type_abbrevs.get(_raw_plot_type, _raw_plot_type)

        raw_scenarios = inputs.get('scenarios', None)
        # Include scenario name only when there is no comparison (exactly one scenario).
        if raw_scenarios is None or check_is_list_of_lists(raw_scenarios) or len(list(raw_scenarios)) != 1:
            scen_part = None
        else:
            scen_part = str(list(raw_scenarios)[0]).replace(' ', '_')

        plot_years_val = inputs.get('plot_years', default_inputs['plot_years'])
        if plot_years_val is None:
            time_part = 'all_years'
        elif len(plot_years_val) == 1:
            time_part = str(plot_years_val[0])
        elif inputs.get('plot_categories') is not None:
            time_part = None
        else:
            time_part = f'{plot_years_val[0]}-{plot_years_val[-1]}'

        # For multi-year panel plots with a single selected category, include that category in the name.
        # 'All' should only be included when it represents a real category dimension being aggregated over.
        categories_val = inputs.get('categories')
        if isinstance(category_label, list):
            category_label_in_data = all(c in df.columns for c in category_label)
        else:
            category_label_in_data = category_label in df.columns
        if (inputs.get('plot_categories') is None and isinstance(categories_val, list) and len(categories_val) == 1
                and (categories_val[0] != 'All' or category_label_in_data)):
            category_part = str(categories_val[0]).replace(' ', '_')
        else:
            category_part = None

        name_parts = ['spatial', file_base]
        if category_part:
            name_parts.append(category_part)
        name_parts.append(plot_type_str)
        if scen_part:
            name_parts.append(scen_part)
        if time_part:
            name_parts.append(time_part)
        name = '_'.join(name_parts)
        inputs['plot_name'] = os.path.join(inputs['plot_directory'], name)
    elif 'plot_name' in inputs and '/' not in inputs['plot_name']:
        # If the user specified only a file name (no path) for the plot name, put the plot in the plot directory.
        inputs['plot_name'] = os.path.join(inputs['plot_directory'], inputs['plot_name'])

    # Add keys for plotting options that have not been specified in the inputs dictionary and use default values for them.
    for key in default_inputs.keys():
        if key not in inputs:
            inputs[key] = default_inputs[key]
    
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

    if inputs['aggregation_function'] == 'area_weighted_mean' and 'area' not in df.columns:
        raise KeyError(f"Error: aggregation_function is 'area_weighted_mean' "
                       f"but no 'area' column found in '{inputs['output_file']}'. "
                       f"Consider using 'mean' or 'sum' instead, or add area data using gcam_add_areas_to_files.py.")

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
                  f"and will produce empty/zero spatial plots: {missing_scenarios}")

    # Verify the plot directory is writable before doing any work.
    plot_name = inputs['plot_name']
    plot_dir = os.path.dirname(plot_name) or '.'
    if not os.access(plot_dir, os.W_OK):
        raise OSError(f"Error: plot directory is not writable: '{plot_dir}'")

    return inputs

def plot_spatial_data(inputs):
    """ 
    Creates a spatial plot and perform statistical analysis for a single output file. The data in the file are organized
    into individual scenarios or scenario sets (ensembles), categories, and regions.

    Parameters:
        input: Dictionary containing user inputs for different plotting options, where the keys are options and values are choices for those options.
               This dictionary is assumed to be complete (pre-processed).

    Returns:
        N/A.
    """
    # This function creates a spatial plot where the data (scenarios or scenario sets, categories, regions) all come from a single output file.
    start_time = time.time()
    output_file = inputs['output_file']
    
    # Extract all other plotting options.
    aggregate_to_region = inputs['aggregate_to_region']
    basin_label = inputs['basin_label']
    categories = inputs['categories']
    category_label = inputs['category_label']
    cbar_limits = inputs['cbar_limits']
    cbar_on = inputs['cbar_on']
    cmap_color = inputs['cmap']
    height = inputs['height']
    key_columns = inputs['key_columns']
    landtype_groups = inputs['landtype_groups']
    linewidth = inputs['linewidth']
    multiplier = inputs['multiplier']
    aggregation_function = inputs['aggregation_function']
    p_value_file = inputs['p_value_file']
    p_value_file_print_only_if_below_threshold = inputs['p_value_file_print_only_if_below_threshold']
    p_value_threshold = inputs['p_value_threshold']
    plot_name = inputs['plot_name']
    plot_type = inputs['plot_type']
    produce_png = inputs['produce_png']
    region_label = inputs['region_label']
    scenario_label = inputs['scenario_label']
    scenarios = inputs['scenarios']
    shape_file = inputs['shape_file']
    shape_file_basin_label = inputs['shape_file_basin_label']
    shape_file_region_label = inputs['shape_file_region_label']
    stippling_hatches = inputs['stippling_hatches']
    stippling_on = inputs['stippling_on']
    title = inputs['title']
    title_size = inputs['title_size']
    use_latex = inputs['use_latex']
    value_label = inputs['value_label']
    units = inputs.get('units', None)
    width = inputs['width']
    x_tick_label_size = inputs['x_tick_label_size']
    y_tick_label_size = inputs['y_tick_label_size']
    year_label = inputs['year_label']

    # Append type/units suffix to title; define units_str for panel titles.
    if plot_type == 'percent_difference':
        if title and '%' not in title and 'percent' not in title:
            title += r' (\% difference)'
        units_str = f' (% diff in {units})' if units else ' (%)'
    elif plot_type in ('absolute_difference', 'mean') and units:
        title += f' ({units})'
        units_str = f' ({units})'
    else:
        units_str = ''
    plot_options = dict(width=width, height=height, name=plot_name, produce_png=produce_png)
    plot_options.update(zip(['x_tick_label_size', 'y_tick_label_size', 'use_latex'], [x_tick_label_size, y_tick_label_size, use_latex]))

    # Use LaTeX fonts for figures and set font size of tick labels.
    setup_plot_params(plot_options)

    # Read the data file into a Pandas DataFrame and filter to the requested years if specified.
    df = read_file_into_dataframe(output_file)
    # 'All' should only be shown in panel titles when it represents a real category dimension being aggregated over.
    if isinstance(category_label, list):
        category_label_in_data = all(c in df.columns for c in category_label)
    else:
        category_label_in_data = category_label in df.columns
    selected_category = categories[0] if len(categories) == 1 and (categories[0] != 'All' or category_label_in_data) else None
    plot_years = inputs.get('plot_years', None)
    if plot_years is not None:
        df = df[df[year_label].isin(plot_years)]
    # Apply the multiplier to the value column (this could be used to change units, for example).
    df[value_label] *= multiplier

    # Set the appropriate dictionary for the landtype group.
    if landtype_groups == 'modified':
        landtype_groups = gcam_landtype_groups
    elif landtype_groups == 'original':
        landtype_groups = gcam_landtype_groups_original

    # Create separate DataFrames for the following cases: 1) all categories (create a copy of the entire DataFrame), 2) categories that correspond to
    # a group of landtypes (e.g., forest, crop, grass, shrub, pasture), and 3) a set of individual categories. Concatenate into a single DataFrame.
    dataframes = []
    if 'All' in categories:
        df_all = df.copy()
        df_all[category_label] = 'All'
        dataframes.append(df_all) 
        categories.remove('All')
    if any(item in landtype_groups for item in categories):
        col = category_label if isinstance(category_label, str) else category_label[0]
        for landtype_group in landtype_groups:
            if landtype_group in categories:
                # Only treat as a group if the sub-landtypes are actually present in the data.
                if not df[col].isin(landtype_groups[landtype_group]).any():
                    continue
                df_this_landtype_group = produce_dataframe_for_landtype_group(df, landtype_group, category_label,
                    value_label, landtype_groups, aggregation_function, key_columns)
                dataframes.append(df_this_landtype_group)
                categories.remove(landtype_group)
    if categories:
        if isinstance(category_label, list):
            df = df[df[category_label].apply(tuple, axis=1).isin(set(categories))]
        else:
            df = df[df[category_label].isin(categories)]
        dataframes.append(df)
    df = pd.concat(dataframes).reset_index()
    if aggregate_to_region and basin_label in df.columns:
        df = df.drop(columns=[basin_label])

    # Read the shape file into a GeoDataFrame from the GeoPandas library and get all regions in the file.
    gdf = gpd.read_file(shape_file)
    if basin_label not in df.columns:
        # The data has no sub-region (basin) granularity, so merge the basin polygons into region-level
        # polygons to avoid drawing meaningless sub-region boundary lines.
        gdf = gdf.dissolve(by=shape_file_region_label).reset_index()
    regions = gdf[shape_file_region_label].unique()

    # For both the individual plots and ensemble plots, create a common variable called scenario_list that will contain a 1D list of all scenarios.
    if not check_is_list_of_lists(scenarios):
        # Option 1: Individual plots, in which each such spatial plot can include one or more individual (not grouped) data sets.
        scenario_list = scenarios
        # Create column labels for these scenarios.
        scenario_columns = [f'scen={i}' for i in range(len(scenarios))]
    else:
        # Option 2: Ensemble plots, in which the ensemble is further subdivided into groups, where each group represents a set of scenarios.
        # The scenarios in this case are contained in a list of lists (i.e., a 2D or nested list). These will be unraveled into the 1D scenario_list.
        scenario_list = []
        scenario_columns = []
        num_scenario_sets = len(scenarios[0])
        num_scenarios_in_each_set = len(scenarios)
        # Group all scenarios that belong to the same set (they share the same column in the 2D list) and put them together in the 1D scenario_list.
        for index in range(num_scenario_sets):
            scenarios_in_set = [scenarios[i][index] for i in range(num_scenarios_in_each_set)]
            scenario_list.extend(scenarios_in_set)
            scenarios_indices_this_set = [f'scen={index}_{i}' for i in range(num_scenarios_in_each_set)]
            scenario_columns.extend(scenarios_indices_this_set)

    # If plot_years is set without plot_categories, produce one panel per year in a single multi-panel figure.
    plot_categories = inputs.get('plot_categories', None)
    if plot_years is not None and plot_categories is None:
        ncols = inputs.get('plot_years_ncols') or len(plot_years)
        nrows = math.ceil(len(plot_years) / ncols)
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(width * ncols, height * nrows))
        axes_flat = list(axes.flat) if hasattr(axes, 'flat') else [axes]
        num_scenario_sets_local = len(scenarios[0]) if check_is_list_of_lists(scenarios) else 0
        vmin, vmax = (cbar_limits[0], cbar_limits[1]) if cbar_limits else (None, None)

        for panel_idx, year in enumerate(plot_years):
            ax = axes_flat[panel_idx]
            df_year = df[df[year_label] == year]
            if df_year.empty:
                print(f"Warning: year {year} not found in data. Skipping panel.")
                ax.set_visible(False)
                continue
            gdf_panel = gdf.copy()
            for col in scenario_columns:
                gdf_panel[col] = float('nan')

            for index, scenario in enumerate(scenario_list):
                df_this_scenario = df_year[df_year[scenario_label] == scenario]
                for region in regions:
                    region_filter = gdf_panel[shape_file_region_label] == region
                    if basin_label in df_this_scenario.columns:
                        basins = gdf_panel[region_filter][shape_file_basin_label].unique()
                        for basin in basins:
                            basin_filter = gdf_panel[shape_file_basin_label] == basin
                            if basin not in gcam_basin_names_and_abbrevations:
                                gdf_panel.loc[region_filter & basin_filter, scenario_columns[index]] = 0
                            else:
                                basin_abbrv = gcam_basin_names_and_abbrevations[basin]
                                df_rows = df_this_scenario[
                                    (df_this_scenario[region_label] == region) &
                                    (df_this_scenario[basin_label] == basin_abbrv)
                                ]
                                gdf_panel.loc[region_filter & basin_filter, scenario_columns[index]] = _aggregate(
                                    df_rows, value_label, aggregation_function)
                    else:
                        df_rows = df_this_scenario[df_this_scenario[region_label] == region]
                        gdf_panel.loc[region_filter, scenario_columns[index]] = _aggregate(
                            df_rows, value_label, aggregation_function)

            gdf_panel.fillna(0, inplace=True)

            if len(scenario_list) == 1:
                gdf_panel['plot'] = gdf_panel[scenario_columns[0]]
            elif plot_type == 'mean':
                gdf_panel['plot'] = gdf_panel[scenario_columns].mean(axis=1)
            elif plot_type in ('absolute_difference', 'percent_difference'):
                if not check_is_list_of_lists(scenarios):
                    columns_control = [scenario_columns[0]]
                    columns_test = scenario_columns[1:]
                else:
                    columns_control = [x for x in scenario_columns if x.startswith('scen=0')]
                    columns_test = [x for x in scenario_columns if not x.startswith('scen=0')]
                    if num_scenario_sets_local == 2:
                        df_tmp = pd.DataFrame({c: gdf_panel[c] for c in columns_control + columns_test})
                        gdf_panel['p_value'] = df_tmp.apply(
                            perform_ttest, columns_set_1=columns_control, columns_set_2=columns_test, axis=1).fillna(1)
                control_data = gdf_panel[columns_control].mean(axis=1)
                test_data = gdf_panel[columns_test].mean(axis=1)
                if plot_type == 'percent_difference':
                    gdf_panel['plot'] = (test_data - control_data) / (control_data + EPSILON) * 100
                else:
                    gdf_panel['plot'] = test_data - control_data

            title_text = f'{selected_category} ({year}){units_str}' if selected_category is not None else f'{year}{units_str}'
            panel_title = f'{chr(ord("a") + panel_idx)}) {title_text}'
            ax.set_title(panel_title, fontdict={'fontsize': title_size})
            gdf_panel.plot('plot', ax=ax, legend=cbar_on, cmap=cmap_color, vmin=vmin, vmax=vmax,
                           legend_kwds={'shrink': .5}, edgecolor='k', linewidth=linewidth)
            if 'p_value' in gdf_panel.columns and stippling_on:
                gdf_panel[gdf_panel['p_value'] <= p_value_threshold].plot(
                    ax=ax, facecolor='none', color='none', hatch=stippling_hatches, linewidth=0)

        for panel_idx in range(len(plot_years), nrows * ncols):
            axes_flat[panel_idx].set_visible(False)

        plot_options['name'] = plot_name
        fig.set_size_inches(width * ncols, height * nrows)
        save_figure(plot_name, fig, plot_options)
        plt.close(fig)
        print(f"Elapsed time for producing plot {plot_name}: {time.time() - start_time:.2f} seconds")
        return

    # If plot_categories is set, produce one file per year in plot_years with one panel per category.
    if plot_categories is not None:
        years_to_loop = plot_years if plot_years is not None else sorted(df[year_label].unique())
        num_scenario_sets_local = len(scenarios[0]) if check_is_list_of_lists(scenarios) else 0
        for year in years_to_loop:
            df_year = df[df[year_label] == year]
            if df_year.empty:
                print(f"Warning: year {year} not found in data. Skipping.")
                continue
            ncols = inputs.get('plot_categories_ncols') or len(plot_categories)
            nrows = math.ceil(len(plot_categories) / ncols)
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(width * ncols, height * nrows))
            axes_flat = list(axes.flat) if hasattr(axes, 'flat') else [axes]
            vmin, vmax = (cbar_limits[0], cbar_limits[1]) if cbar_limits else (None, None)

            for panel_idx, category in enumerate(plot_categories):
                ax = axes_flat[panel_idx]
                if isinstance(category_label, list) and isinstance(category, (tuple, list)):
                    df_cat = df_year[(df_year[category_label] == list(category)).all(axis=1)]
                else:
                    df_cat = df_year[df_year[category_label] == category]
                if df_cat.empty:
                    print(f"Warning: category '{category}' not found for year {year}. Skipping panel.")
                    ax.set_visible(False)
                    continue
                gdf_panel = gdf.copy()
                for col in scenario_columns:
                    gdf_panel[col] = float('nan')

                for index, scenario in enumerate(scenario_list):
                    df_this_scenario = df_cat[df_cat[scenario_label] == scenario]
                    for region in regions:
                        region_filter = gdf_panel[shape_file_region_label] == region
                        if basin_label in df_this_scenario.columns:
                            basins = gdf_panel[region_filter][shape_file_basin_label].unique()
                            for basin in basins:
                                basin_filter = gdf_panel[shape_file_basin_label] == basin
                                if basin not in gcam_basin_names_and_abbrevations:
                                    gdf_panel.loc[region_filter & basin_filter, scenario_columns[index]] = 0
                                else:
                                    basin_abbrv = gcam_basin_names_and_abbrevations[basin]
                                    df_rows = df_this_scenario[
                                        (df_this_scenario[region_label] == region) &
                                        (df_this_scenario[basin_label] == basin_abbrv)
                                    ]
                                    gdf_panel.loc[region_filter & basin_filter, scenario_columns[index]] = _aggregate(
                                        df_rows, value_label, aggregation_function)
                        else:
                            df_rows = df_this_scenario[df_this_scenario[region_label] == region]
                            gdf_panel.loc[region_filter, scenario_columns[index]] = _aggregate(
                                df_rows, value_label, aggregation_function)

                gdf_panel.fillna(0, inplace=True)

                if len(scenario_list) == 1:
                    gdf_panel['plot'] = gdf_panel[scenario_columns[0]]
                elif plot_type == 'mean':
                    gdf_panel['plot'] = gdf_panel[scenario_columns].mean(axis=1)
                elif plot_type in ('absolute_difference', 'percent_difference'):
                    if not check_is_list_of_lists(scenarios):
                        columns_control = [scenario_columns[0]]
                        columns_test = scenario_columns[1:]
                    else:
                        columns_control = [x for x in scenario_columns if x.startswith('scen=0')]
                        columns_test = [x for x in scenario_columns if not x.startswith('scen=0')]
                        if num_scenario_sets_local == 2:
                            df_tmp = pd.DataFrame({c: gdf_panel[c] for c in columns_control + columns_test})
                            gdf_panel['p_value'] = df_tmp.apply(
                                perform_ttest, columns_set_1=columns_control, columns_set_2=columns_test, axis=1).fillna(1)
                    control_data = gdf_panel[columns_control].mean(axis=1)
                    test_data = gdf_panel[columns_test].mean(axis=1)
                    if plot_type == 'percent_difference':
                        gdf_panel['plot'] = (test_data - control_data) / (control_data + EPSILON) * 100
                    else:
                        gdf_panel['plot'] = test_data - control_data

                panel_title = f'{chr(ord("a") + panel_idx)}) {category} ({year}){_panel_units_str(df_cat, plot_type, units_str)}'
                ax.set_title(panel_title, fontdict={'fontsize': title_size})
                gdf_panel.plot('plot', ax=ax, legend=cbar_on, cmap=cmap_color, vmin=vmin, vmax=vmax,
                               legend_kwds={'shrink': .5}, edgecolor='k', linewidth=linewidth)
                if 'p_value' in gdf_panel.columns and stippling_on:
                    gdf_panel[gdf_panel['p_value'] <= p_value_threshold].plot(
                        ax=ax, facecolor='none', color='none', hatch=stippling_hatches, linewidth=0)

            for panel_idx in range(len(plot_categories), nrows * ncols):
                axes_flat[panel_idx].set_visible(False)

            year_plot_name = plot_name if plot_years is not None and len(plot_years) == 1 else f'{plot_name}_{year}'
            plot_options['name'] = year_plot_name
            fig.set_size_inches(width * ncols, height * nrows)
            save_figure(year_plot_name, fig, plot_options)
            plt.close(fig)
            print(f"Elapsed time for producing plot {year_plot_name}: {time.time() - start_time:.2f} seconds")
        return

    # For each scenario, region, and basin, calculate the mean or sum of all relevant categories over all years and put into the GeoDataFrame.
    for index, scenario in enumerate(scenario_list):
        df_this_scenario = df[df[scenario_label] == scenario]
        for region in regions:
            region_filter = gdf[shape_file_region_label] == region
            if basin_label in df_this_scenario.columns:
                # If the DataFrame containing the data of interest are further organized into basins, add basin-level information.
                basins = gdf[region_filter][shape_file_basin_label].unique()
                for basin in basins:
                    basin_filter = gdf[shape_file_basin_label] == basin
                    if basin not in gcam_basin_names_and_abbrevations:
                        # Set the value to 0 if there is no data for the current region-and-basin combination.
                        gdf.loc[region_filter & basin_filter, scenario_columns[index]] = 0
                    else:
                        # Basin names are fully written out in the GeoDataFrame (shape file), while they are abbreviated in the DataFrame (data file).
                        basin_abbrv = gcam_basin_names_and_abbrevations[basin]
                        df_rows = df_this_scenario[
                            (df_this_scenario[region_label] == region) &
                            (df_this_scenario[basin_label] == basin_abbrv)
                        ]
                        gdf.loc[region_filter & basin_filter, scenario_columns[index]] = _aggregate(
                            df_rows, value_label, aggregation_function)
            else:
                # Uniform value per region; all rows in region are aggregated together.
                df_rows = df_this_scenario[df_this_scenario[region_label] == region]
                gdf.loc[region_filter, scenario_columns[index]] = _aggregate(
                    df_rows, value_label, aggregation_function)
    # Fill any missing values in the GeoDataFrame with 0.
    gdf.fillna(0, inplace=True)

    if len(scenario_list) == 1:
        # If there is only one scenario in the list, simply plot that column.
        gdf['plot'] = gdf[scenario_columns[0]]
    elif plot_type == 'mean':
        # If there are multiple data columns, one option is to plot the mean across all columns.
        gdf['plot'] = gdf[scenario_columns].mean(axis=1)
    elif plot_type == 'absolute_difference' or plot_type == 'percent_difference':
        if not check_is_list_of_lists(scenarios):
            # If an individual plot, assume the first column is the reference for the absolute difference or percent difference calculation.
            columns_control = [scenario_columns[0]]
            columns_test = [x for x in scenario_columns[1:]] 
        else:
            # If an ensemble plot, group the scenarios into sets and calculate set means. Assume the first scenario set is the reference set.
            columns_control = [x for x in scenario_columns if x.startswith('scen=0')]
            columns_test = [x for x in scenario_columns if not x.startswith('scen=0')]
            if check_is_list_of_lists(scenarios) and num_scenario_sets == 2:
                # If there are two scenario sets, do a t-test at each individual region-and-basin combination and put results into the GeoDataFrame.
                df = pd.DataFrame()
                df[columns_control], df[columns_test] = gdf[columns_control], gdf[columns_test]
                gdf['p_value'] = df.apply(perform_ttest, columns_set_1=columns_control, columns_set_2=columns_test, axis=1).fillna(1)
        # Calculate either the absolute difference or percent difference between the means of the test and control data sets.
        control_data = gdf.loc[:, columns_control].mean(axis=1)
        test_data = gdf.loc[:, columns_test].mean(axis=1)
        if plot_type == 'percent_difference':
            gdf['plot'] = (test_data - control_data)/(control_data + EPSILON)*100
        else:
            # Absolute difference.
            gdf['plot'] = test_data - control_data
        # Do a t-test for the test and control data sets as a whole (global comparison). Print the output from this into a data file.
        ttest = stats.ttest_ind(control_data, test_data)
        print_p_values(ttest, columns_test, p_value_threshold, p_value_file, plot_name, p_value_file_print_only_if_below_threshold)

    # Create figure and axis objects for the plot, set the title and colorbar limits.
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.set_title(title, fontdict={'fontsize': title_size})
    if cbar_limits:
        vmin, vmax = cbar_limits[0], cbar_limits[1]
    else:
        vmin, vmax = None, None
        
    # Generate the plot and optionally add stippling to indicate statistically significant differences at individual regions and/or basins.
    gdf.plot('plot', ax=ax, legend=cbar_on, cmap=cmap_color, vmin=vmin, vmax=vmax, legend_kwds={"shrink": .5}, edgecolor='k', linewidth=linewidth)
    if 'p_value' in gdf.columns and stippling_on:
        mask = gdf['p_value'] <= p_value_threshold
        gdf[mask].plot(ax=ax, facecolor='none', color='none', hatch=stippling_hatches, linewidth=0)
    
    # Finalize the spatial plot now that all data sets have been processed.
    plot_options['name'] = plot_name
    fig.set_size_inches(width, height)
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
        print('Usage: python gcam_plot_spatial_data.py `path/to/json/input/file(s)\'')
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

    # Delete all the p-value files before we do any calculations to start a fresh run.
    for inputs in list_of_inputs:
        file = inputs['p_value_file']
        if os.path.exists(file): 
            os.remove(file)

    # Create all of the box plots in parallel.
    with multiprocessing.Pool(processes=MAX_PROCESSES) as pool:
        pool.map(plot_spatial_data, list_of_inputs)
    
    # Sort all the p-value files alphabetically.
    for inputs in list_of_inputs:
        file = inputs['p_value_file']
        if os.path.exists(file): 
            sort_file(file)
    
    # Print the total execution time to produce all the plots.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all plots: {elapsed_time:.2f} seconds")