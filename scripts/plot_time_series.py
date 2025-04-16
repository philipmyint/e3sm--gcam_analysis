import json
from matplotlib import pyplot as plt
import pandas as pd
import sys
import time
from utility_plots import plot_colors, set_figure_options

def plot_time_series(output_files, plot_name, plot_labels, plot_variables, y_label, start_year=2015, 
                     end_year=2100, width=8.5, height=7, x_scale='linear', y_scale='linear', plot_Jan_July=True):

    fig, ax = plt.subplots(nrows=1, ncols=1)
    for index, file in enumerate(output_files):

        df = pd.read_fwf(file)
        df = df[(df['Year'] >= start_year) & (df['Year'] <= end_year)]

        x = df.groupby('Year', as_index=False).mean()['Year']
        y_mean = df.groupby('Year', as_index=False).mean()[plot_variables[0]]
        for var_index in range(1, len(plot_variables)):
            y_mean += df.groupby('Year', as_index=False).mean()[plot_variables[var_index]]
        ax.plot(x, y_mean, label=plot_labels[index]+ ' (Annual mean)', color=plot_colors[index], linestyle='solid', linewidth=2)

        if plot_Jan_July:
            y_Jan = df[df['Month'] == 1][plot_variables[0]]
            y_July = df[df['Month'] == 7][plot_variables[0]]
            for var_index in range(1, len(plot_variables)):
                y_Jan += df[df['Month'] == 1][plot_variables[var_index]]
                y_July += df[df['Month'] == 7][plot_variables[var_index]]
            ax.plot(x, y_Jan, label=plot_labels[index] + ' (January)', color=plot_colors[index], linestyle='dotted', linewidth=2)
            ax.plot(x, y_July, label=plot_labels[index] + ' (July)', color=plot_colors[index], linestyle='dashed', linewidth=2)
   
    plot_options = dict(width=width, height=height, name=plot_name, x_label='Year', y_label=fr'{y_label}')
    plot_options['x_scale'] = x_scale
    plot_options['y_scale'] = y_scale
    set_figure_options(fig,ax,plot_options)

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Usage: plot_time_series.py `path/to/json/input/file\'')
        sys.exit()

    input_file = sys.argv[1]
    with open(input_file) as f:
        inputs = json.load(f)

    for job_index in range(len(inputs)):
        start_time = time.time()
        plot_time_series(**inputs[job_index])
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time for {inputs[job_index]['plot_name']}: {elapsed_time:.2f} seconds")