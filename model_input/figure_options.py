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

# matplotlib colormap to use
cmap = 'RdYlBu_r'

sea_lvl_color = matplotlib.cm.RdYlBu(0)

# plot the acual elements
# if False will make a scatter plot with a dot at centre of each element
show_elements = True

# show concentration only
conc_only = False

#
#model_scenario = ''

year = 365.25 * 24 * 60 * 60

# resolution for showing contours
dx = 0.0025
dy = 0.0025

# plot every x flow arrows
thin_arrows = 100

streamline_density = 1.0
lw_min = 0.5
lw_max = 0.5

#ylim_fixed = [-140, 245]
ylim_fixed = None

#tax_ylim_fixed = [-1.0, 6.0]
tax_ylim_fixed = None

add_legend = True
add_colorbar = True
leg_fontsize = 'small'

add_title = False

#base_dir = '/home/elco/model_files/grompy_sgd_models'
base_dir = 'model_output'

#folder = 'model_files/grompy_sgd_models/base_case_highk'
folder = None