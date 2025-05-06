from matplotlib import pyplot as plt

# Default values for different plotting options.
width_default = 10  # inches.
height_default = 8  # inches.
scale_default = 'linear'    # other options include 'linear' and 'log'.
axis_limits_default = None
num_contour_levels_default = 20
axis_label_size_default = 24
tick_label_size_default = 20
legend_label_size_default = 14
legend_on_default = True
linewidth_default = 2
produce_png_default = False
use_latex_default = False
""" Hex codes of Matplotlib Tableau color palette: blue, orange, green, red, purple, brown, pink, gray, olive, cyan. """
plot_colors_default = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
""" Tuples for different line styles in plots. """
linestyle_tuples_default = [
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

""" Dictionary of default input values for time series plots. """
default_inputs_time_series = {'plot_directory': './',
                    'calculation_type': 'mean',
                    'multiplier': 1,
                    'start_year': 2015,
                    'end_year': 2100,
                    'width': width_default,
                    'height': height_default,
                    'x_scale': scale_default,
                    'y_scale': scale_default,
                    'x_limits': axis_limits_default,
                    'y_limits': axis_limits_default,
                    'x_tick_label_size': tick_label_size_default,
                    'y_tick_label_size': tick_label_size_default,
                    'x_label_size': axis_label_size_default,
                    'y_label_size': axis_label_size_default,
                    'legend_label_size': legend_label_size_default,
                    'legend_on': legend_on_default,
                    'linewidth': linewidth_default,
                    'plot_colors': plot_colors_default,
                    'linestyle_tuples': linestyle_tuples_default,
                    'use_latex': use_latex_default, 
                    'produce_png': produce_png_default, 
                    'include_seasons': {'MAM': False, 'JJA': False, 'SON': False, 'DJF': False},
                    'seasons_to_plot_separately': {'MAM': False, 'JJA': False, 'SON': False, 'DJF': False},
                    'monthly_time_series_plot': False,
                    'monthly_time_series_start_year': 2071,
                    'monthly_time_series_end_year': 2090,
                    'monthly_time_series_x_limits': axis_limits_default,
                    'monthly_time_series_y_limits': axis_limits_default
}

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
    levels = options.get('levels', num_contour_levels_default)
    #levels = [0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300]
    plot = ax.contourf(x, y, z, levels=levels)
    cbar = fig.colorbar(plot)
    cbar.ax.tick_params(labelsize=tick_label_size_default) 
    cbar_label = options.get('cbar_label', None)
    if cbar_label:
        cbar.set_label(fr'{cbar_label}', fontsize=axis_label_size_default)
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
    if options.get('use_latex', use_latex_default):
        # Use LaTeX fonts for figures and set font size of tick labels.
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', weight='bold')
    legend_on = options.get('legend_on', legend_on_default)
    if legend_on:
        ax.legend(prop={'size': options.get('legend_label_size', legend_label_size_default)}, frameon=False, loc='best')
    else:
        ax.legend().set_visible(False)
    ax.set_xlabel(options['x_label'], fontsize=options.get('x_label_size', axis_label_size_default))
    plt.rcParams['xtick.labelsize'] = options.get('x_tick_label_size', tick_label_size_default)
    ax.set_ylabel(options['y_label'], fontsize=options.get('y_label_size', axis_label_size_default))
    plt.rcParams['ytick.labelsize'] = options.get('y_tick_label_size', tick_label_size_default)
    x_scale = options.get('x_scale', scale_default)
    ax.set_xscale(x_scale)
    y_scale = options.get('y_scale', scale_default)
    ax.set_yscale(y_scale)
    x_limits = options.get('x_limits', axis_limits_default)
    if not x_limits:
        ax.set_xlim(x_limits)
    y_limits = options.get('y_limits', axis_limits_default)
    if not y_limits:
        ax.set_ylim(y_limits)
    width = options.get('width', width_default)
    height = options.get('height', height_default)
    fig.set_size_inches(width, height)
    name = options['name']
    produce_png = options.get('produce_png', produce_png_default)
    if produce_png or name.endswith('.png'):
        if name.endswith('.png'):
            fig.savefig(f'{name}', format='png')
        else:
            fig.savefig(f'{name}.png', format='png')
    else:
        if name.endswith('.pdf'):
            fig.savefig(f'{name}', format='pdf')
        else:
            fig.savefig(f'{name}.pdf', format='pdf')