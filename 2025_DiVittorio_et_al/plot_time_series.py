import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sys

sys.path.append('/home/ac.myint1/scripts')
from utility_plots import plot_colors, set_figure_options
from utility_dataframes import clean_up_dataframe

file_directory = './'
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
    print(y_Jan)

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
set_figure_options(fig,ax,plot_options)