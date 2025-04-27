from matplotlib import pyplot as plt
import numpy as np
from matplotlib import cm, ticker

# Use LaTeX fonts for figures and set font size of tick labels.
plt.rc('text',usetex=True)
plt.rc('font',family='serif',weight='bold')
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20

""" Hex codes of Matplotlib Tableau color palette. """
""" Colors: blue, orange, green, red, purple, brown, pink, gray, olive, cyan. """
plot_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

""" Dictionary of default input values for time series plots. """
default_inputs_time_series = {'plot_directories': './',
                    'plot_names': 'time_series',
                    'calculation_types': 'mean',
                    'multipliers': 1,
                    'start_years': 2015,
                    'end_years': 2100,
                    'widths': 10.0,
                    'heights': 8.0,
                    'x_scales': 'linear',
                    'y_scales': 'linear',
                    'x_limits': 'default',
                    'y_limits': 'default',
                    'include_seasons': {'spring': False, 'summer': False, 'autumn': False, 'winter': False},
                    'seasons_to_plot_separately': {'spring': False, 'summer': False, 'autumn': False, 'winter': False},
                    'monthly_time_series_plot': False,
                    'monthly_time_series_start_year': 2071,
                    'monthly_time_series_end_year': 2090,
                    'monthly_time_series_x_limits': 'default',
                    'monthly_time_series_y_limits': 'default'
                    }

""" Tuples for different line styles in plots. """
linestyle_tuples = [
    ('solid',                 (0, ())),
    ('dashed',                (0, (5, 5))),
    ('dotted',                (0, (1, 5))),
    ('dashdot',               (0, (3, 5, 1, 5))),
    ('loosely dotted',        (0, (1, 10))),
    ('loosely dashed',        (0, (5, 10))),
    ('densely dashed',        (0, (5, 1))),
    ('loosely dashdotted',    (0, (3, 10, 1, 10))),
    ('densely dashdotted',    (0, (3, 1, 1, 1))),
    ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
    ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
    ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]

def create_contour_plot(x, y, z, options):
    """
    Creates a contour plot.

    Parameters:
        x: NumPy array or list containing the quantity to plot on the x axis.
        y: NumPy array or list containing the quantity to plot on the y axis.
        z: NumPy array or list containing the quantity to plot on the z axis.
        options: Dictionary that contains plotting options like number of contour levels or colorbar labels.

    Returns:
        N/A.
    """
    fig,ax = plt.subplots(nrows=1, ncols=1)
    levels = options.get('levels', 20)
    #levels = [0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300]
    plot = ax.contourf(x, y, z, levels=levels)
    cbar = fig.colorbar(plot)
    cbar.ax.tick_params(labelsize=20) 
    cbar_label = options.get('cbar_label', None)
    if cbar_label:
        cbar.set_label(fr'{cbar_label}', fontsize=24)
    cbar.set_ticks(cbar.locator.tick_values(z.min(), z.max()))
    set_figure_options(fig, ax, options)

def set_figure_options(fig, ax, options):
    """
    Sets options when creating a figure.

    Parameters:
        fig: Object for the figure of interest.
        ax: Axes object corresponding to the figure of interest. 
        options: Dictionary that contains plotting options like x and y labels, the name of the figure, etc.

    Returns:
        N/A.
    """
    _,labels = ax.get_legend_handles_labels()
    if labels:
        ax.legend(prop={'size': 14}, frameon=False, loc='best')
    ax.set_xlabel(options['x_label'], fontsize=24)
    ax.set_ylabel(options['y_label'], fontsize=24)
    x_scale = options.get('x_scale', None)
    if x_scale:
        ax.set_xscale(x_scale)
    y_scale = options.get('y_scale', None)
    if y_scale:
        ax.set_yscale(y_scale)
    x_limits = options.get('x_limits', 'default')
    if x_limits != 'default':
        ax.set_xlim(x_limits)
    y_limits = options.get('y_limits', 'default')
    if y_limits != 'default':
        ax.set_ylim(y_limits)
    width = options['width']
    height = options['height']
    fig.set_size_inches(width, height)
    name = options['name']
    fig.savefig(f'{name}.pdf', format='pdf')