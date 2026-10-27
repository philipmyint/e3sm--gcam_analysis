import itertools
import inspect
import json
import math
from matplotlib import pyplot as plt
import multiprocessing
import os
import pandas as pd
import re
import seaborn as sns
import sys
import time
from utility_constants import *
from utility_dataframes import read_file_into_dataframe
from utility_functions import check_is_list_of_lists, get_max_processes_from_inputs, transpose_scenarios_if_needed
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
    'multipanel': True,
    'multipanel_hspace': 0.35,
    'multipanel_wspace': 0.18,
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
    'title_size': axis_label_size_default,
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

def _safe_plot_name_part(value):
    """Return a filesystem-friendly name fragment for a year or category."""
    if isinstance(value, (list, tuple)):
        value = '_'.join(str(item) for item in value)
    value = str(value).strip().replace(' ', '_')
    return re.sub(r'[^A-Za-z0-9_.=-]+', '_', value).strip('_')


def _add_suffix_to_plot_name(plot_name, *suffix_parts):
    """Insert non-empty suffix parts before the plot file extension."""
    clean_parts = [_safe_plot_name_part(part) for part in suffix_parts if part is not None]
    clean_parts = [part for part in clean_parts if part]
    if not clean_parts:
        return plot_name
    directory, filename = os.path.split(plot_name)
    stem, extension = os.path.splitext(filename)
    updated_name = f"{stem}_{'_'.join(clean_parts)}{extension}"
    return os.path.join(directory, updated_name) if directory else updated_name


def _format_units_for_label(units, use_latex=False):
    """Return a readable units label, or None if the units value is empty/noisy."""
    if units is None:
        return None
    units_str = str(units).strip()
    if units_str.upper() in ('', 'NA', 'NONE', 'NAN', 'NULL') or units_str in ('(', ')', '()'):
        return None
    units_str = units_str.strip('()').strip()
    if units_str.upper() in ('', 'NA', 'NONE', 'NAN', 'NULL'):
        return None
    units_str = units_str.replace('thous km2', 'thousands km2')
    if use_latex:
        units_str = units_str.replace('^2', r'$^2$')
        units_str = units_str.replace('km2', r'km$^2$')
        units_str = units_str.replace('m2', r'm$^2$')
    return units_str


def _valid_units_fragment(fragment):
    """Return True if a piece of a unit label is worth displaying."""
    if fragment is None:
        return False
    fragment = str(fragment).strip().strip('()').strip()
    return bool(fragment) and fragment.upper() not in ('NA', 'NONE', 'NAN', 'NULL')


def _apply_multipanel_spacing(fig, wspace, hspace):
    """Keep year panels compact without clipping rotated labels or the shared y label."""
    fig.tight_layout(rect=(0.04, 0.02, 1, 0.98), pad=0.4, h_pad=0.6, w_pad=0.4)
    fig.subplots_adjust(wspace=wspace, hspace=hspace)


def _calculate_percent_difference(df_year, scenario_label, scenario_list, value_label):
    """
    Return rows for non-reference scenarios with values converted to percent difference.

    GCAM outputs are not always perfectly rectangular: a scenario may be missing
    rows that exist in the reference scenario, or vice versa. Pair rows by their
    shared identifier columns instead of relying on row order/array length.
    """
    reference_scenario = scenario_list[0]
    reference_df = df_year[df_year[scenario_label] == reference_scenario].copy()
    comparison_df = df_year[df_year[scenario_label] != reference_scenario].copy()

    if reference_df.empty:
        print(f"Warning: no reference data for scenario '{reference_scenario}'. Skipping percent-difference panel.")
        return df_year.iloc[0:0].copy()
    if comparison_df.empty:
        print("Warning: no comparison data found for percent-difference panel.")
        return comparison_df

    ignored_columns = {scenario_label, value_label, 'area', 'index', 'level_0'}
    merge_columns = [column for column in df_year.columns if column not in ignored_columns]
    if not merge_columns:
        if len(reference_df) != len(comparison_df):
            print("Warning: cannot pair reference and comparison rows for percent-difference panel. Skipping panel.")
            return df_year.iloc[0:0].copy()
        comparison_df[value_label] = (
            (comparison_df[value_label].to_numpy() - reference_df[value_label].to_numpy())
            / (reference_df[value_label].to_numpy() + EPSILON) * 100
        )
        return comparison_df

    reference_values = reference_df[merge_columns + [value_label]].rename(
        columns={value_label: '_reference_value'}
    )
    if reference_values.duplicated(merge_columns).any():
        reference_values = (
            reference_values
            .groupby(merge_columns, as_index=False, dropna=False)['_reference_value']
            .mean()
        )

    merged = comparison_df.merge(reference_values, on=merge_columns, how='inner')
    missing_rows = len(comparison_df) - len(merged)
    if missing_rows > 0:
        print(f"Warning: dropped {missing_rows} row(s) without matching reference data for percent-difference panel.")
    if merged.empty:
        print("Warning: no matched reference/comparison rows for percent-difference panel.")
        return comparison_df.iloc[0:0].copy()

    merged[value_label] = (
        (merged[value_label] - merged['_reference_value'])
        / (merged['_reference_value'] + EPSILON) * 100
    )
    return merged.drop(columns=['_reference_value'])


def _draw_boxplot(df_year, x_variable, value_label, hue, ax, linewidth, palette,
                  marker_size, fill_boxes, box_properties, line_properties,
                  flier_properties):
    """Draw a seaborn boxplot with kwargs supported by the installed seaborn."""
    boxplot_parameters = inspect.signature(sns.boxplot).parameters
    supports_fill = 'fill' in boxplot_parameters

    kwargs = dict(data=df_year, x=x_variable, y=value_label, ax=ax,
                  linewidth=linewidth, palette=palette, fliersize=marker_size,
                  boxprops=box_properties, whiskerprops=line_properties,
                  capprops=line_properties, medianprops=line_properties,
                  flierprops=flier_properties)
    if hue:
        kwargs['hue'] = hue
    if supports_fill:
        kwargs['fill'] = fill_boxes

    sns.boxplot(**kwargs)

    if not fill_boxes and not supports_fill:
        for patch in ax.artists:
            patch.set_facecolor('none')


def process_inputs(inputs):
    """ 
    Processes a dictionary of inputs (keys are options, values are choices for those options) for creating box (or box-and-whisker) plots.

    Parameters:
        inputs: Dictionary containing the user plotting choice inputs for different options. This dictionary may be incomplete or have invalid values.

    Returns:
        Dictionary that completely specifies all plotting options. 
        If the user did not make a choice for a particular option, the default choice for that plotting option will be selected.
    """
    # Check that the output data file exists before trying to read it.
    # Missing data files are skipped so one absent optional plot does not stop a long plotting run.
    if not os.path.exists(inputs['output_file']):
        print(f"Warning: output file not found, skipping: '{inputs['output_file']}'")
        return None
    if os.path.getsize(inputs['output_file']) == 0:
        print(f"Warning: output file is empty, skipping: '{inputs['output_file']}'")
        return None
    df = read_file_into_dataframe(inputs['output_file'])
    plot_years_was_specified = 'plot_years' in inputs

    # If the category label (e.g., sector or landtype) has not been specified, use the default value.
    if 'category_label' not in inputs:
        category_label = default_inputs['category_label']
        inputs['category_label'] = category_label
    else:
        category_label = inputs['category_label']

    # If categories were not specified, populate them from the data. If the
    # category column is absent, use "All" so the plot aggregates over the whole file.
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

    # Use the output CSV name as the default y-axis label and as the base for
    # generated plot names.
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
    inputs['_plot_years_was_specified'] = plot_years_was_specified

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

    # area_weighted_mean requires area metadata from gcam_process_extracted_data.py.
    if inputs['aggregation_function'] == 'area_weighted_mean' and 'area' not in df.columns:
        raise KeyError(f"Error: aggregation_function is 'area_weighted_mean' "
                       f"but no 'area' column found in '{inputs['output_file']}'. "
                       f"Consider using 'mean' or 'sum' instead, or add area data "
                       f"with gcam_process_extracted_data.py and a land_allocation_file entry.")

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

    # Verify the plot directory is writable before doing the plotting work.
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
    # This function creates one multi-panel box-and-whisker figure from a single
    # processed GCAM CSV. The data may be split by scenario, scenario set,
    # category, region, basin, hue, and year depending on the JSON options.
    start_time = time.time()
    output_file = inputs['output_file']
    
    # Extract all other plotting options.
    basin_label = inputs['basin_label']
    basins = list(inputs['basins']) if isinstance(inputs['basins'], list) else inputs['basins']
    categories = list(inputs['categories']) if isinstance(inputs['categories'], list) else inputs['categories']
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
    multipanel = inputs['multipanel']
    multipanel_hspace = inputs['multipanel_hspace']
    multipanel_wspace = inputs['multipanel_wspace']
    multiplier = inputs['multiplier']
    aggregation_function = inputs['aggregation_function']
    plot_colors = inputs['plot_colors']
    plot_name = inputs['plot_name']
    plot_percent_difference = inputs['plot_percent_difference']
    plot_years = inputs['plot_years']
    plot_years_was_specified = inputs.get('_plot_years_was_specified', True)
    plot_type = inputs['plot_type']
    produce_png = inputs['produce_png']
    region_label = inputs['region_label']
    regions = list(inputs['regions']) if isinstance(inputs['regions'], list) else inputs['regions']
    scenario_label = inputs['scenario_label']
    scenario_sets = inputs['scenario_sets']
    scenarios = inputs['scenarios']
    title_size = inputs['title_size']
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

    # Collect the plotting options that are shared across all year panels.
    plot_options = dict(width=width, height=height, name=plot_name, x_label=x_label, y_label=fr'{y_label}')
    plot_options.update(zip(['x_scale', 'y_scale', 'y_limits', 'use_latex'], [x_scale, y_scale, y_limits, use_latex]))
    plot_options.update(zip(['x_tick_label_size', 'y_tick_label_size', 'legend_on'], [x_tick_label_size, y_tick_label_size, legend_on]))
    plot_options.update(zip(['x_label_size', 'y_label_size', 'legend_label_size'], [x_label_size, y_label_size, legend_label_size]))
    plot_options.update(zip(['produce_png', 'legend_x_offset', 'legend_place_outside'], [produce_png, legend_x_offset, legend_place_outside]))
    plot_options.update(zip(['legend_num_columns'], [legend_num_columns]))

    # Use LaTeX fonts for figures and set font size of tick labels.
    setup_plot_params(plot_options)

    # Read the data file, filter to requested years, and apply the user-specified multiplier.
    df = read_file_into_dataframe(output_file)
    if not multipanel:
        years_to_loop = plot_years if plot_years is not None else sorted(df[year_label].unique())
        # If plot_years came from the script default rather than the JSON file,
        # quietly drop default years that are absent from this dataset. Explicit
        # user-requested years are kept so missing-year warnings still surface.
        if plot_years is not None and not plot_years_was_specified:
            available_years = set(df[year_label].unique())
            years_to_loop = [year for year in years_to_loop if year in available_years]
        include_year_in_name = plot_years is None or len(plot_years) != 1
        for year in years_to_loop:
            child_inputs = inputs.copy()
            child_inputs['multipanel'] = True
            child_inputs['plot_years'] = [year]
            if include_year_in_name:
                child_inputs['plot_name'] = _add_suffix_to_plot_name(plot_name, year)
            plot_box_and_whiskers(child_inputs)
        return
    if plot_years is not None:
        df = df[df[year_label].isin(plot_years)]
    df[value_label] *= multiplier

    # Create regional subsets for plotting. "Global" is represented by a copy
    # of the whole DataFrame with the region label overwritten, so it pools all regions.
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

    # Option 1: individual plots, where each scenario remains its own plotted data set.
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
    
    # Option 2: ensemble plots, where scenarios are grouped into named scenario sets.
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

    # From this point on, individual and ensemble plots share the same category
    # preparation and plotting logic.

    # Set the appropriate dictionary for the landtype group.
    if landtype_groups == 'modified':
        landtype_groups = gcam_landtype_groups
    elif landtype_groups == 'original':
        landtype_groups = gcam_landtype_groups_original

    # Create separate DataFrames for three category cases: all rows, named GCAM
    # landtype groups (e.g., forest/crop/grass), and exact category values.
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

    # Percent-difference plots are dimensionless, so label them directly as
    # percent differences instead of carrying through physical units such as
    # thousands km2. This also avoids dangling labels when source unit metadata
    # are empty or malformed.
    units_by_category = {}
    if plot_percent_difference:
        y_label_lower = y_label.lower()
        if '%' not in y_label and 'percent' not in y_label_lower and 'difference' not in y_label_lower:
            percent_label = r'(\% difference)' if use_latex else '(% difference)'
            y_label += f' {percent_label}'
            plot_options['y_label'] = fr'{y_label}'
    # If different plotted categories use related units such as kg/yr and kg/m2/yr,
    # keep the shared numerator on the y-axis and put each denominator on the x tick.
    elif isinstance(category_label, str) and x_variable == category_label and 'units' in df.columns:
        for category, group in df.groupby(category_label):
            units = [_format_units_for_label(unit, use_latex) for unit in group['units'].dropna().unique()]
            units = [unit for unit in units if unit]
            if len(units) == 1 and '/' in units[0]:
                numerator, denominator = units[0].split('/', 1)
                if _valid_units_fragment(numerator) and _valid_units_fragment(denominator):
                    units_by_category[str(category)] = units[0]
        numerators = {unit.split('/', 1)[0] for unit in units_by_category.values()}
        if units_by_category and len(numerators) == 1:
            numerator = numerators.pop()
            if _valid_units_fragment(numerator) and numerator not in y_label:
                y_label += f' (% diff in {numerator})' if plot_percent_difference else f' ({numerator})'
                plot_options['y_label'] = fr'{y_label}'
    if not plot_percent_difference:
        units = [_format_units_for_label(unit, use_latex) for unit in df['units'].dropna().unique()] if 'units' in df.columns else []
        units = [unit for unit in units if unit]
        if not units_by_category and len(units) == 1 and units[0] not in y_label:
            y_label += f' ({units[0]})'
            plot_options['y_label'] = fr'{y_label}'

    # Each requested year gets its own panel in the same output figure.
    years_to_plot = plot_years if plot_years is not None else sorted(df[year_label].unique())
    ncols = inputs['plot_years_ncols'] or len(years_to_plot)
    nrows = math.ceil(len(years_to_plot) / ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(width * ncols, height * nrows))
    axes_flat = list(axes.flat) if hasattr(axes, 'flat') else [axes]

    for panel_idx, year in enumerate(years_to_plot):
        ax = axes_flat[panel_idx]
        df_year = df[df[year_label] == year].copy()
        if df_year.empty:
            if plot_years_was_specified:
                print(f"Warning: year {year} not found in data. Skipping panel.")
            ax.set_visible(False)
            continue
        if plot_percent_difference:
            # Percent-difference box plots use the first scenario/scenario set as
            # the reference and then drop that reference from the plotted data.
            df_year = _calculate_percent_difference(df_year, scenario_label, scenario_list, value_label)

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
            _draw_boxplot(df_year, x_variable, value_label, hue, ax, linewidth,
                          palette, marker_size, fill_boxes, box_properties,
                          line_properties, flier_properties)
        else:
            _draw_boxplot(df_year, x_variable, value_label, None, ax, linewidth,
                          palette, marker_size, fill_boxes, box_properties,
                          line_properties, flier_properties)
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
        ax.set_title(str(year), fontsize=title_size)
        panel_options = plot_options.copy()
        panel_options['legend_on'] = legend_on if hue else False
        panel_options['y_label'] = ''
        set_figure_options(fig, ax, panel_options)

    for panel_idx in range(len(years_to_plot), nrows * ncols):
        axes_flat[panel_idx].set_visible(False)

    fig.set_size_inches(width * ncols, height * nrows)
    fig.supylabel(plot_options['y_label'], fontsize=y_label_size)
    _apply_multipanel_spacing(fig, multipanel_wspace, multipanel_hspace)
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

    # Use any max_processes value requested in the inputs, otherwise use the default.
    # Do not start more worker processes than there are independent plotting jobs.
    max_processes = get_max_processes_from_inputs(inputs, MAX_PROCESSES)
    num_processes = min(max_processes, max(1, len(list_of_inputs)))

    # Create all of the box plots in parallel.
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(plot_box_and_whiskers, list_of_inputs)
    
    # Print the total execution time to produce all the plots.
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all plots: {elapsed_time:.2f} seconds")
