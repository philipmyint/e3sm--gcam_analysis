from matplotlib import pyplot as plt
import numpy as np
from matplotlib import cm, ticker

# Use LaTeX fonts for figures and set font size of tick labels.
plt.rc('text',usetex=True)
plt.rc('font',family='serif',weight='bold')
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20

""" List of colors expressed in terms of their hex codes. """
plot_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

def create_contour_plot(x, y, z, options):
    fig,ax = plt.subplots(nrows=1,ncols=1)
    levels = options.get('levels',20)
    #levels = [0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300]
    plot = ax.contourf(x,y,z,levels=levels)
    cbar = fig.colorbar(plot)
    cbar.ax.tick_params(labelsize=20) 
    cbar_label = options.get('cbar_label',None)
    if cbar_label:
        cbar.set_label(fr'{cbar_label}',fontsize=24)
    cbar.set_ticks(cbar.locator.tick_values(z.min(),z.max()))
    set_figure_options(fig,ax,options)

def set_figure_options(fig, ax, options):
    _,labels = ax.get_legend_handles_labels()
    if labels:
        ax.legend(prop={'size': 14},frameon=False,loc='best')
    ax.set_xlabel(options['x_label'],fontsize=24)
    ax.set_ylabel(options['y_label'],fontsize=24)
    x_scale = options.get('x_scale',None)
    if x_scale:
        ax.set_xscale(x_scale)
    y_scale = options.get('y_scale',None)
    if y_scale:
        ax.set_yscale(y_scale)
    width = options['width']
    height = options['height']
    fig.set_size_inches(width,height)
    name = options['name']
    fig.savefig(f'{name}.pdf',format='pdf')