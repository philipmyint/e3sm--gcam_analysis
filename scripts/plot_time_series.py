import json
from matplotlib import pyplot as plt
import pandas as pd
import sys
import time
from utility_plots import plot_colors, set_figure_options

def plot_time_series(plot_names, plot_variables, y_labels, output_files, plot_labels, multipliers, start_years, 
                     end_years, plot_widths, plot_heights, x_scales, y_scales, plot_Jan_July_options):

    for plot_index, plot_name in enumerate(plot_names):

        start_time = time.time()

        plot_variable = plot_variables[plot_index]
        y_label = y_labels[plot_index]
        output_files_set = output_files[plot_index]
        plot_labels_set = plot_labels[plot_index]
        multiplier = multipliers[plot_index]
        start_year = start_years[plot_index]
        end_year = end_years[plot_index]
        width = plot_widths[plot_index]
        height = plot_heights[plot_index]
        x_scale = x_scales[plot_index]
        y_scale = y_scales[plot_index]
        plot_Jan_July = plot_Jan_July_options[plot_index]

        fig, ax = plt.subplots(nrows=1, ncols=1)
        for file_index, file in enumerate(output_files_set):

            df = pd.read_fwf(file)
            df = df[(df['Year'] >= start_year) & (df['Year'] <= end_year)]

            for column in df.columns:
                if plot_variable in column:
                    plot_variable = column

            x = df.groupby('Year', as_index=False).mean()['Year']
            y_mean = df.groupby('Year', as_index=False).mean()[plot_variable]*multiplier
            ax.plot(x, y_mean, label=plot_labels_set[file_index] + ' (Annual mean)', color=plot_colors[file_index], linestyle='solid', linewidth=2)

            if plot_Jan_July:
                y_Jan = df[df['Month'] == 1][plot_variable]*multiplier
                y_July = df[df['Month'] == 7][plot_variable]*multiplier 
                ax.plot(x, y_Jan, label=plot_labels_set[file_index] + ' (January)', color=plot_colors[file_index], linestyle='dotted', linewidth=2)
                ax.plot(x, y_July, label=plot_labels_set[file_index] + ' (July)', color=plot_colors[file_index], linestyle='dashed', linewidth=2)
    
        plot_options = dict(width=width, height=height, name=plot_name, x_label='Year', y_label=fr'{y_label}')
        plot_options['x_scale'] = x_scale
        plot_options['y_scale'] = y_scale
        set_figure_options(fig,ax,plot_options)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time for {plot_name}: {elapsed_time:.2f} seconds")

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Usage: plot_time_series.py `path/to/json/input/file\'')
        sys.exit()

    start_time = time.time()
    input_file = sys.argv[1]
    with open(input_file) as f:
        inputs = json.load(f)

    for job_index in range(len(inputs)):

        keys = inputs[job_index].keys()
        values = inputs[job_index].values()
        num_plots = len(inputs[job_index]['plot_names'])
        values = [value*num_plots if len(value) == 1 else value for value in values]
        arguments = dict(zip(keys, values))
        plot_time_series(**arguments)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for producing all plots: {elapsed_time:.2f} seconds")