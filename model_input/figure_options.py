"""
this file contains options to control Grompy output figures made using the make_model_fig.py script

"""

import matplotlib


use_gui = False

# figure file format, choose '.png', '.pdf' or '.jpg'. note the pdf file size can be quite large
figure_type = '.png'

# resolution
dpi = 200.0

debug = False

# matplotlib colormap to use, see here for all the options:
# https://matplotlib.org/examples/color/colormaps_reference.html
cmap = 'RdYlBu_r'

sea_lvl_color = matplotlib.cm.RdYlBu(0)

# plot the acual elements
# if False will make a scatter plot with a dot at centre of each element
show_elements = True

# show concentration only or a figure for each model parameter
conc_only = False

year = 365.25 * 24 * 60 * 60

# x and y resolution for showing contours (m)
dx = 0.0025
dy = 0.0025

# plot every x flow arrows
thin_arrows = 100

# density of streamlines
streamline_density = 1.0

# linewidth of streamlines. min and max are the streamline thickness for the lowest and highest flow velocity.
lw_min = 0.5
lw_max = 0.5

# limits of the y-axis. set to None for letting the script choose the limits
#ylim_fixed = [-140, 245]
ylim_fixed = None

# limits of the y-axis of the top panel. set to None for letting the script choose the limits
#tax_ylim_fixed = [-1.0, 6.0]
tax_ylim_fixed = None

# add a legend
add_legend = True
leg_fontsize = 'small'

# add a colorbar
add_colorbar = True

# add a figure title
add_title = False

# directories where make_model_figure.py searches for VTK files if no file is provided:
base_dir = 'model_output'
folder = None