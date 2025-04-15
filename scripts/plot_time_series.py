import json
from matplotlib import pyplot as plt
import multiprocessing
import numpy as np
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
        y_mean = pd.DataFrame(np.zeros((len(x), 1)))
        for variable in plot_variables:
            y_mean += df.groupby('Year', as_index=False).mean()[variable]
        ax.plot(x, y_mean, label=plot_labels[index]+ ' (Annual mean)', color=plot_colors[index], linestyle='solid', linewidth=2)

        if plot_Jan_July:
            y_Jan = pd.DataFrame(np.zeros((len(x), 1)))
            y_July = pd.DataFrame(np.zeros((len(x), 1)))
            for variable in plot_variables:
                y_Jan += df[df['Month'] == 1][variable]
                y_July += df[df['Month'] == 7][variable]
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

'''file_directory = './'
plot_labels = [r'Control', r'Full feedback']

fig, ax = plt.subplots(nrows=1, ncols=1)
files = ['control_eam_h0_time_series.dat', 'full_feedback_eam_h0_time_series.dat']
for index, file in enumerate(files):

    file = file_directory + file

    df = pd.read_fwf(file)
    x = df[df['Month'] == 1]['Year']
    y_mean = df.groupby('Year',as_index=False)['TREFHT (K)'].mean()['TREFHT (K)']
    y_Jan = df[df['Month'] == 1]['TREFHT (K)']
    y_July = df[df['Month'] == 7]['TREFHT (K)']

    ax.plot(x,y_mean, label=plot_labels[index]+ ' (Annual mean)', color=plot_colors[index], linestyle='solid', linewidth=2)
    ax.plot(x,y_Jan, label=plot_labels[index] + ' (January)', color=plot_colors[index], linestyle='dotted', linewidth=2)
    ax.plot(x,y_July, label=plot_labels[index] + ' (July)', color=plot_colors[index], linestyle='dashed', linewidth=2)

width = 8.5
height = 7
x_label = r'Year'
y_label = r'Reference height temperature (K)'
plot_directory = file_directory
file = file_directory + 'time_series_TREFHT'
plot_options = dict(width=width, height=height, name=file, x_label=x_label, y_label=y_label)
plot_options['x_scale'] = 'linear'
plot_options['y_scale'] = 'linear'
set_figure_options(fig,ax,plot_options)

fig, ax = plt.subplots(nrows=1, ncols=1)
files = ['control_eam_h0_time_series.dat', 'full_feedback_eam_h0_time_series.dat']
for index, file in enumerate(files):

    file = file_directory + file

    df = pd.read_fwf(file)
    x = df[df['Month'] == 1]['Year']
    y_mean = df.groupby('Year',as_index=False)['PRECC (m/s)'].mean()['PRECC (m/s)']
    y_Jan = df[df['Month'] == 1]['PRECC (m/s)']
    y_July = df[df['Month'] == 7]['PRECC (m/s)']

    ax.plot(x,y_mean, label=plot_labels[index]+ ' (Annual mean)', color=plot_colors[index], linestyle='solid', linewidth=2)
    ax.plot(x,y_Jan, label=plot_labels[index] + ' (January)', color=plot_colors[index], linestyle='dotted', linewidth=2)
    ax.plot(x,y_July, label=plot_labels[index] + ' (July)', color=plot_colors[index], linestyle='dashed', linewidth=2)

width = 8.5
height = 7
x_label = r'Year'
y_label = r'Convective precipitation rate (m/s)'
plot_directory = file_directory
file = file_directory + 'time_series_PRECC'
plot_options = dict(width=width, height=height, name=file, x_label=x_label, y_label=y_label)
plot_options['x_scale'] = 'linear'
plot_options['y_scale'] = 'linear'
set_figure_options(fig,ax,plot_options)


fig, ax = plt.subplots(nrows=1, ncols=1)
files = ['control_elm_h0_time_series.dat', 'full_feedback_elm_h0_time_series.dat']
for index, file in enumerate(files):

    file = file_directory + file

    df = pd.read_fwf(file)
    x = df[df['Month'] == 1]['Year']
    y_mean = df.groupby('Year',as_index=False)['NPP (gC/m^2/s)'].mean()['NPP (gC/m^2/s)']
    y_Jan = df[df['Month'] == 1]['NPP (gC/m^2/s)']
    y_July = df[df['Month'] == 7]['NPP (gC/m^2/s)']

    ax.plot(x,y_mean, label=plot_labels[index]+ ' (Annual mean)', color=plot_colors[index], linestyle='solid', linewidth=2)
    ax.plot(x,y_Jan, label=plot_labels[index] + ' (January)', color=plot_colors[index], linestyle='dotted', linewidth=2)
    ax.plot(x,y_July, label=plot_labels[index] + ' (July)', color=plot_colors[index], linestyle='dashed', linewidth=2)

width = 8.5
height = 7
x_label = r'Year'
y_label = r'Net primary production (g C/m$^2$/s)'
plot_directory = file_directory
file = file_directory + 'time_series_NPP'
plot_options = dict(width=width, height=height, name=file, x_label=x_label, y_label=y_label)
plot_options['x_scale'] = 'linear'
plot_options['y_scale'] = 'linear'
set_figure_options(fig,ax,plot_options)'''
