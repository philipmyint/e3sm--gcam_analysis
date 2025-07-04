from matplotlib import pyplot as plt

# Default values for different plotting options.
width_default = 10  # inches.
height_default = 8  # inches.
scale_default = 'linear'    # options include 'linear' and 'log'.
axis_limits_default = None
num_contour_levels_default = 20
axis_label_size_default = 24
tick_label_size_default = 20
legend_label_size_default = 14
legend_num_columns_default = 1
legend_on_default = True
legend_place_outside_default = False
linewidth_default = 2
produce_png_default = False
use_latex_default = False
bbox_inches_default = 'tight'

""" Hex codes of Matplotlib Tableau color palette: blue, orange, green, red, purple, brown, pink, gray, olive, cyan. """
plot_colors_default = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
""" XKCD colors (https://xkcd.com/color/rgb/): teal, magenta, sea green, bright pink, dark orange, goldenrod, forest, dirt, coral, baby blue, peach. """
plot_colors_default.extend(['#029386', '#c20078', '#53fca1', '#fe01b1', '#c65102', '#fac205', '#0b5509', '#8a6e45', '#fc5a50', '#a2cffe', '#ffb07c'])

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

""" Markers (https://matplotlib.org/stable/gallery/lines_bars_and_markers/marker_reference.html). """
markers_default = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X', '1', '2', '3', '4', '+', 'x', '|']

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
    fig, ax = plt.subplots(nrows=1, ncols=1)
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

def save_figure(name, fig, options):
    """
    Saves a figure as either a .pdf (default) or .png (if figure name ends with .png or if produce_png in the options dictionary is True).

    Parameters:
        name: Name of the figure.
        fig: Object for the figure of interest.
        options: Dictionary that contains plotting options, including the value of produce_png.

    Returns:
        N/A.
    """
    produce_png = options.get('produce_png', produce_png_default)
    bbox_inches = options.get('bbox_inches', bbox_inches_default)
    if produce_png or name.endswith('.png'):
        if name.endswith('.png'):
            fig.savefig(f'{name}', format='png', bbox_inches=bbox_inches)
        else:
            fig.savefig(f'{name}.png', format='png', bbox_inches=bbox_inches)
    else:
        if name.endswith('.pdf'):
            fig.savefig(f'{name}', format='pdf', bbox_inches=bbox_inches)
        else:
            fig.savefig(f'{name}.pdf', format='pdf', bbox_inches=bbox_inches)

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
        legend_num_columns = options.get('legend_num_columns', legend_num_columns_default)
        legend_place_outside = options.get('legend_place_outside', legend_place_outside_default)
        if not legend_place_outside:
            ax.legend(prop={'size': options.get('legend_label_size', legend_label_size_default)}, frameon=False, loc='best', ncol=legend_num_columns)
        else:
            legend_x_offset = options.get('legend_x_offset', None)
            if not legend_x_offset:
                if legend_num_columns == 1:
                    legend_x_offset = 1.3
                else:
                    legend_x_offset = 1 + 0.45*legend_num_columns
            plt.legend(prop={'size': options.get('legend_label_size', legend_label_size_default)}, 
                   frameon=False, loc='center right', bbox_to_anchor=(legend_x_offset, 0.5), ncol=legend_num_columns)
            
    ax.set_xlabel(options['x_label'], fontsize=options.get('x_label_size', axis_label_size_default))
    plt.rcParams['xtick.labelsize'] = options.get('x_tick_label_size', tick_label_size_default)
    ax.set_ylabel(options['y_label'], fontsize=options.get('y_label_size', axis_label_size_default))
    plt.rcParams['ytick.labelsize'] = options.get('y_tick_label_size', tick_label_size_default)

    x_scale = options.get('x_scale', scale_default)
    if x_scale:
        ax.set_xscale(x_scale)
    y_scale = options.get('y_scale', scale_default)
    if y_scale:
        ax.set_yscale(y_scale)
    x_limits = options.get('x_limits', axis_limits_default)
    if x_limits:
        ax.set_xlim(x_limits)
    y_limits = options.get('y_limits', axis_limits_default)
    if y_limits:
        ax.set_ylim(y_limits)

    width = options.get('width', width_default)
    height = options.get('height', height_default)
    fig.set_size_inches(width, height)
    name = options['name']
    save_figure(name, fig, options)