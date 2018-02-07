"""
Compare model results to salt wedge intrusion experiments by Goswami and Clement (2007) WRR 43.

Note the filename of the .csv file containing the model experiment data and the VTK files for the three model runs
are hardcoded in this python script. Use you favorite Python editor to change them.

This script generates a figure comparing the modeled to the measured position of the salt wedge, which is saved as
benchmark_data/model_vs_experimental_salt_wedge.pdf. It also reports the calculated freshwater flux to the screen.
"""

__author__ = 'Elco Luijendijk'

import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
import lib.read_vtu_file


def depth_sw_interface_Glover1959(x, Q, K, rho_f, rho_s):

    """
    Calculate depth of the fresh-salt water interface in a coastal aquifer, following Glover (1959).
    """

    gamma = (rho_s - rho_f) / rho_f

    y2 = 2 * Q / (gamma * K) * x + Q**2 / (gamma**2 * K**2)

    depth = y2**0.5

    return depth


# file containing the experimental results:
fn_exp = 'benchmark_data/table_a1_steady_state_salt_wedge_locations.csv'

vtk_files_init = ['benchmark_data/model_runs/runS0_specified_pressure_[0, 68.569]_final_output',
                  'benchmark_data/model_runs/runS1_specified_pressure_[0, 19.591]_final_output',
                  'benchmark_data/model_runs/runS2_specified_pressure_[0, 53.876]_final_output']

max_conc = 0.03624
target_cs = np.array([0.10, 0.50, 0.90]) * max_conc

# read experimental data
df_exp = pd.read_csv(fn_exp)

# read model results
vtk_files = [f + '_Elements.vtu' for f in vtk_files_init]
vtkf_files = [f + '_FaceElements.vtu' for f in vtk_files_init]

xs = []
ys = []
concs = []

freshwater_fluxes = []

for vtk_file, vtkf_file in zip(vtk_files, vtkf_files):
    print 'reading vtk file %s' % vtk_file

    xy, conn, pt_var_names, pt_var_arrays, cell_var_names, cell_var_arrays = \
        lib.read_vtu_file.read_vtu_file(vtk_file)

    #fnf = os.path.join(folder, vtkf_file)
    fnf = vtkf_file

    print 'reading vtk file %s' % fnf

    #xy_pts_face, conn_face, pt_var_names_face, pt_var_arrays_face, cell_var_names_face, cell_var_arrays_face = \
    #    lib.read_vtu_file.read_vtu_file(vtkf_file)

    c_ind = pt_var_names.index('concentration')
    conc = pt_var_arrays[c_ind]
    concs.append(conc)

    yi = np.unique(xy[:, 1])

    x_concs = np.zeros((target_cs.shape[0], yi.shape[0]))
    for i, yii in enumerate(yi):
        ind = np.where(xy[:, 1] == yii)
        x_interp = np.interp(target_cs, conc[ind][::-1], xy[:, 0][ind][::-1])
        x_concs[:, i] = x_interp

    xs.append(x_concs)
    ys.append(yi)

    # calculate discharge
    qx = cell_var_arrays[cell_var_names.index('qx')]
    qy = cell_var_arrays[cell_var_names.index('qy')]

    # calculate centre of elements
    elements = np.array([xy[conn_i] for conn_i in conn])
    element_x = np.mean(elements[:, :, 0], axis=1)
    element_y = np.mean(elements[:, :, 1], axis=1)

    # select elements on the left hand model boundary
    ind_left = (element_x == element_x.min()) #& (element_y <= 0.255)
    qx_left = qx[ind_left]
    qy_left = qy[ind_left]
    element_y_left = element_y[ind_left]

    # select elements on the right hand model boundary
    ind_right = element_x == element_x.max()
    qx_right = qx[ind_right]
    qy_right = qy[ind_right]
    element_y_right = element_y[ind_right]

    yi_sorted = np.sort(np.unique(xy[:, 1]))
    dy = yi_sorted[1] - yi_sorted[0]

    # convert flux to m2 s-1
    Qx_left = qx_left * dy
    Qx_right = qx_right * dy

    #
    total_flux_left = np.sum(Qx_left)
    total_flux_left_in = np.sum(Qx_left[Qx_left > 0])
    total_flux_left_out = np.sum(Qx_left[Qx_left < 0])

    total_flux_right = np.sum(Qx_right)
    total_flux_right_in = np.sum(Qx_right[Qx_right > 0])
    total_flux_right_out = np.sum(Qx_right[Qx_right < 0])

    freshwater_fluxes.append(total_flux_left_out)

# calculate width dispersion zone between fresh and salt water
# should be approx 1 wide acc. to experimental results:
print '-' * 20
print 'Width of modeled fresh-salt water mixing zone:'
width_mixing_zone = np.array([(xsi[0] - xsi[2]) for xsi in xs])
for i, w in enumerate(width_mixing_zone):
    print 'run %i, average width of mixing zone = %0.03f m' % (i + 1, np.mean(w[w > 0]))

# calculate analytical depth of transition fresh-salt water acc. to Glover (1959):
# note, fluxes reported in Goswami & Clement include third dimension with length of 2.7 cm
Qs = np.array([1.42, 0.59, 1.19]) / 2.7 * 1e-4
K = 1050.0 / (24.0 * 60.0 * 60.0)
dh = np.array([0.012, 0.007, 0.0105])
dx = 0.53
dhdx = dh / dx
b = 0.26
rho_f = 998.7
rho_s = 1026.0
xg = np.arange(0, 0.531, 0.001)
dg = [depth_sw_interface_Glover1959(xg, Q, K, rho_f, rho_s) for Q in Qs]

yg = 0.26 - np.array(dg)

###############
# report fluxes
###############
print '-' * 20
print 'Comparison modeled and observed fluxes, steady-state experiments 1, 2 and 3 (cm3/sec):'
Q_measured = Qs * 1e4 * 2.7
print 'measured: ', Q_measured
Q_modeled_transformed = -np.array(freshwater_fluxes) * 1e4 * 2.7
print 'modeled:  ', Q_modeled_transformed

##################################
# calculate misfit models and data
##################################
mean_errors = []
mean_abs_errors = []
rmses = []
for i in xrange(3):

    xc = 'x_ss%i' % (i + 1)
    yc = 'y_ss%i' % (i + 1)

    modeled_values = np.interp(df_exp[yc] / 100.0, ys[i], xs[i][1])
    df_exp['x_modeled_%i' % (i + 1)] = modeled_values

    model_error = modeled_values - df_exp[xc] / 100.0
    mean_errors.append(np.mean(model_error))
    mean_abs_errors.append(np.mean(np.abs(model_error)))
    rmses.append(np.sqrt(np.mean(model_error**2)))

print '-' * 20
print 'Misfit modeled and observed x coordinate salt wedge (m):'
for i in xrange(3):
    print 'experiment %i' % (i + 1)
    print '\tmean error = %0.4f' % mean_errors[i]
    print '\tmean abs error = %0.4f' % mean_abs_errors[i]
    print '\trmse = %0.4f' % rmses[i]
print '-' * 20

####################################
# make a figure of the model results
####################################
markers = ['s', 'd', '^']

fig = pl.figure()

ax = fig.add_subplot(1, 1, 1)
ax.grid(b=False)
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Elevation (m)')
ax.set_xlim(0.0, 0.53)
ax.set_ylim(0.0, 0.26)

leg_data_all = []

for i in xrange(3):

    xc = 'x_ss%i' % (i + 1)
    yc = 'y_ss%i' % (i + 1)

    leg_f = ax.fill_betweenx(ys[0], xs[i][0], xs[i][2], edgecolor='grey', facecolor='lightgrey', zorder=0)

    leg_data = ax.scatter(df_exp[xc] / 100.0, df_exp[yc] / 100.0, marker=markers[i],
                           facecolor='darkgray', edgecolor='black', s=60)
    leg_data_all.append(leg_data)

    leg_g, = ax.plot(xg, yg[i], ls=':', color='black')
    leg_m, = ax.plot(xs[i][1], ys[i], color='black')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

legs = leg_data_all + [leg_g, leg_m, leg_f]
dlabels = ['0.5 isochlor experiment %i' % i for i in range(1, 4)]
labels = dlabels + ['analytical solution, Glover (1959)', 'modeled 0.5 isochlor', 'modeled 0.1 - 0.9 isochlor']

ax.legend(legs, labels, frameon=False, fontsize='small')

fig.tight_layout()

fn_fig = 'benchmark_data/model_vs_experimental_salt_wedge.pdf'
print 'saving figure as %s' % fn_fig
fig.savefig(fn_fig)

fn_fig = 'benchmark_data/model_vs_experimental_salt_wedge.png'
print 'saving figure as %s' % fn_fig
fig.savefig(fn_fig)

print 'done'