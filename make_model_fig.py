import matplotlib
matplotlib.use('Agg')

"""
read vtk file and create model figure

usage:
    

"""

__author__ = 'Elco Luijendijk'

import sys
import os
import itertools
import random

import numpy as np
import matplotlib.pyplot as pl
import scipy.interpolate

import lib.read_vtu_file

from lib.grompy_lib import get_normal_flux_to_bnd

#import useful_functions

from model_input.figure_options import *

from matplotlib import ticker


def add_subplot_axes(ax, rect, axisbg='w'):
    """
    embed a new subplot in existing subplot

    based on:
    http://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
    """

    fig = pl.gcf()
    box = ax.get_position()

    width = box.width
    height = box.height

    inax_position = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)

    x = infig_position[0]
    y = infig_position[1]

    width *= rect[2]
    height *= rect[3]

    subax = fig.add_axes([x, y, width, height])#, axisbg=axisbg)

    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()

    x_labelsize *= rect[2] ** 0.5
    y_labelsize *= rect[3] ** 0.5

    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)

    return subax


def project_flux_to_raster(element_x, element_y, xy_pts, bnd_ind,
                           arrow_x_int=250.0, arrow_y_int=10.0):

    """
    interpolate fluid flux to a regular raster

    """

    # convert fluxes to regular grid to generate equally spaced flow arrows
    topxy = xy_pts[:, :][bnd_ind]

    # flow

    topx_grid = np.arange(topxy[:, 0].min(), topxy[:, 0].max() + arrow_x_int,
                          arrow_x_int)
    node_order = np.argsort(topxy[:, 0])
    topy_grid = np.interp(topx_grid,
                          topxy[:, 0][node_order],
                          topxy[:, 1][node_order])

    # find thickness
    xmin_ind = xy_pts[:, 0] == xy_pts[:, 0].min()
    ys_xmin = xy_pts[:, 1][xmin_ind]
    thickness = ys_xmin.max() - ys_xmin.min()

    # find surface for each element loc
    element_surface = np.interp(element_x,
                                topxy[:, 0][node_order],
                                topxy[:, 1][node_order])
    element_depth = element_y - element_surface

    # interpolate qx and qy to regular depth grid
    qxi = topx_grid

    qyi = np.arange(0, -thickness+arrow_y_int, -arrow_y_int)

    # interpolate pressure on regular grid
    xgq, ygq = np.meshgrid(qxi, qyi)
    xgqf, ygqf = xgq.flatten(), ygq.flatten()

    # interpolate u to grid
    xyq_pts = np.vstack((element_x, element_depth)).T
    qx_interp = scipy.interpolate.griddata(xyq_pts,
                                           qx,
                                           np.vstack((xgqf, ygqf)).T,
                                           method='linear')
    qy_interp = scipy.interpolate.griddata(xyq_pts,
                                           qy,
                                           np.vstack((xgqf, ygqf)).T,
                                           method='linear')

    # project back to surface
    ygqf_surface = np.interp(xgqf,
                             topxy[:, 0][node_order],
                             topxy[:, 1][node_order])
    ygqf_proj = ygqf + ygqf_surface

    return xgqf, ygqf_proj, qx_interp, qy_interp


def convert_to_grid(x, y, qx, qy, dx=1.0, dy=1.0):

    """
    interpolate a variable to a regular raster

    """

    # interpolate pressure on regular grid
    xi = np.arange(x.min(), x.max() + dx, dx)
    yi = np.arange(y.min(), y.max() + dy, dy)
    xg, yg = np.meshgrid(xi, yi)
    xgqf, ygqf = xg.flatten(), yg.flatten()

    # interpolate u to grid
    xyq_pts = np.vstack((x, y)).T
    qx_interp = scipy.interpolate.griddata(xyq_pts,
                                           qx,
                                           np.vstack((xgqf, ygqf)).T,
                                           method='linear')
    qy_interp = scipy.interpolate.griddata(xyq_pts,
                                           qy,
                                           np.vstack((xgqf, ygqf)).T,
                                           method='linear')

    qx_interp1 = np.resize(qx_interp, xg.shape)
    qy_interp1 = np.resize(qy_interp, xg.shape)

    return xg, yg, qx_interp1, qy_interp1


if 'vtu' in sys.argv[-1]:

    vtk_file = sys.argv[-1]

    try:
        assert vtk_file[-4:] == '.vtu'
    except AssertionError:
        raise NameError('file name does not end with .vtu, are you sure this is a grompy/escript VTK file?')

    if '_Elements.vtu' not in vtk_file:
        vtk_file = vtk_file.split('.vtu')[0] + '_Elements.vtu'

    fe_file = vtk_file.split('_Elements.vtu')[0] + '_FaceElements.vtu'
    vtk_files = [vtk_file]
    vtkf_files = [fe_file]

    folder = os.path.split(sys.argv[-1])[0]


elif folder is None:
    dirs = os.listdir(base_dir)
    dirs = [os.path.join(base_dir, directory) for directory in dirs]
    dirs.sort(key=os.path.getmtime)
    dirs = dirs[::-1]

    print('grompy-salt output directories from newest to oldest:')
    for i, directory in enumerate(dirs):
        print(i, directory)

    print('\nenter number or enter for the newest directory:')
    selection = input()

    if selection.isdigit() is True:
        folder = dirs[int(selection)]
    else:
        folder = dirs[0]
    folder = os.path.join(folder, 'vtk_files')

else:
    folder = os.path.join(folder, 'vtk_files')

print('enter f to make a figure for final results only or enter to include all timesteps')
selection = input()
if 'f' in selection:
    final_figs_only = True
else:
    final_figs_only = False

fns = os.listdir(folder)

vtk_files = [fn for fn in fns
             if fn[-4:] == '.vtu' and 'FaceElements' not in fn]
vtkf_files = [fn for fn in fns
             if fn[-4:] == '.vtu' and 'FaceElements' in fn]

if final_figs_only is True:
    vtk_files = [f for f in vtk_files if 'final' in f]
    vtkf_files = [f for f in vtkf_files if 'final' in f]
#if model_scenario is not '' and model_scenario is not None:
#    vtk_files = [f for f in vtk_files if model_scenario in f]
#    vtkf_files = [f for f in vtkf_files if model_scenario in f]

#files.sort(key=os.path.getmtime)
vtk_files = [os.path.join(folder, vtk_file) for vtk_file in vtk_files]
vtkf_files = [os.path.join(folder, vtkf_file) for vtkf_file in vtkf_files]

vtk_files.sort(key=os.path.getmtime)
vtkf_files.sort(key=os.path.getmtime)

vtk_files = vtk_files[::-1]
vtkf_files = vtkf_files[::-1]

# select file to show:
print('grompy-salt output files, from newest to oldest:')
for i, vtk_file in enumerate(vtk_files):
    print(i, os.path.split(vtk_file)[-1])

print('\nenter number of output file to show, two numbers separated by -, ' \
      'a for all files, or r for a random selection of 10 files or enter for the newest file:')
selection = input()

if selection.isdigit() is True:
    vtk_files = [vtk_files[int(selection)]]
    vtkf_files = [vtkf_files[int(selection)]]
elif '-' in selection:
    s2 = selection.split('-')
    start = int(s2[0])
    end = int(s2[1])
    vtk_files = vtk_files[start:end]
    vtkf_files = vtkf_files[start:end]
elif len(selection) > 0 and selection[0] == 'r':
    if len(selection) > 1:
        nrand = int(selection[1:])
    else:
        nrand = 10
    selection_list = [random.randint(0, len(vtk_files)) for i in range(nrand)]

    v = [vtk_files[i] for i in selection_list]
    vf = [vtkf_files[i] for i in selection_list]

    vtk_files = v
    vtkf_files = vf
elif selection == '':
    vtk_files = [vtk_files[0]]
    vtkf_files = [vtkf_files[0]]

print('make a figure of solute concentration (enter) only or all parameters (a)')

if 'a' in input():
    conc_only = False
else:
    conc_only = True

# find directory of files and create output dir:
folder_base = os.path.split(folder)[0]

fig_folder = os.path.join(folder_base, 'fig')
#fig_folder_alt = 'fig'
path = os.path.normpath(folder_base)
b = path.split(os.sep)
#fig_folder_alt = os.path.join('fig', b[-1])

if os.path.exists(fig_folder) is False:
    os.mkdir(fig_folder)
#if os.path.exists(fig_folder_alt) is False:
#    os.mkdir(fig_folder_alt)

for vtk_file, vtkf_file in zip(vtk_files[::-1], vtkf_files[::-1]):


    #fn = os.path.join(folder, vtk_file)
    fn = vtk_file

    print('reading vtk file %s' % fn)

    xy_pts, conn, pt_var_names, pt_var_arrays, cell_var_names, cell_var_arrays = \
        lib.read_vtu_file.read_vtu_file(fn)

    #fnf = os.path.join(folder, vtkf_file)
    fnf = vtkf_file

    print('reading vtk file %s' % fnf)

    xy_pts_face, conn_face, pt_var_names_face, pt_var_arrays_face, cell_var_names_face, cell_var_arrays_face = \
        lib.read_vtu_file.read_vtu_file(fnf)

    vtk_file = os.path.split(vtk_file)[-1]
    vtkf_file = os.path.split(vtkf_file)[-1]

    # construct label from filename
    f = vtk_file.split('.')[:-1]
    f = ''.join(f)
    f1 = f.split('_')
    ind = f1.index('Elements')
    f2 = f1[:ind]
    if f2[-1] == 'output':
        f2 = f2[:-1]

    f3 = f2[0].split('S')

    if add_title is True:
        plot_title = ' '.join(f3)
        plot_title += ', '
        plot_title += ' '.join(f2[1:-2])
        plot_title += '= %s' % f2[-2]
        if f2[-1] == 'final':
            plot_title += ', final timestep'
        else:
            f4 = f2[-1].split('s')
            plot_title += ', timestep %s' % f4[-1]

    # get fluxes
    qx_ind = cell_var_names.index('qx')
    qx = cell_var_arrays[qx_ind]
    qy_ind = cell_var_names.index('qy')
    qy = cell_var_arrays[qy_ind]

    # make pts 2d
    xy_pts = xy_pts[:, :2]
    xy_pts_face = xy_pts_face[:, :2]

    # make closed form of polygons
    conn_cl = np.zeros((conn.shape[0], conn.shape[1] + 1), dtype=int)
    conn_cl[:, :conn.shape[1]] = conn
    conn_cl[:, -1] = conn[:, 0]

    # set up figure
    #fig, ax, cax, tax = gwflow_output.init_figure()

    #fig = pl.figure(figsize=(8, 6))
    fig, axs = pl.subplots(2, 1, gridspec_kw={'height_ratios': [1, 3], 'hspace':0.0}, sharex=True)
    tax, ax = axs

    for axi in axs:
        axi.spines['top'].set_visible(False)
        axi.spines['right'].set_visible(False)
        axi.get_xaxis().tick_bottom()
        axi.get_yaxis().tick_left()

    if add_colorbar is True:
        rect = [0.75, 0.15, 0.3, 0.03]
        cax2 = add_subplot_axes(ax, rect)#, axisbg='w')

    # create polygons
    print('generate polygons for elements')
    polys = [xy_pts[conn_i] for conn_i in conn]
    patches = [matplotlib.patches.Polygon(poly)
               for poly in polys]

    print('generate polygons for face elements')
    polys_face = [xy_pts_face[conn_i] for conn_i in conn_face]
    #patches_face = [matplotlib.patches.Polygon(poly, True)
    #           for poly in polys_face]

    # calculate polygon centre
    polys_array = np.array(polys)
    polys_array_face = np.array(polys_face)

    #
    element_x = np.mean(polys_array[:, :, 0], axis=1)
    element_y = np.mean(polys_array[:, :, 1], axis=1)

    element_x_face = np.mean(polys_array_face[:, :, 0], axis=1)
    element_y_face = np.mean(polys_array_face[:, :, 1], axis=1)

    # find top bnd fluxes
    print('find top bnd fluxes')
    if 'nodal_flux_surface' in pt_var_names:
        f_ind = pt_var_names.index('nodal_flux')
    nfx = pt_var_arrays[f_ind][:, 0]
    nfy = pt_var_arrays[f_ind][:, 1]
    nodal_flux = pt_var_arrays[f_ind]

    #
    nodal_flux *= year

    nodata = -99999
    if 'surface' in pt_var_names:
        sf_ind = pt_var_names.index('surface')
        bnd_ind = pt_var_arrays[sf_ind] == 1
    else:
        sf_ind = pt_var_names.index('nodal_flux_surface')
        nfx = pt_var_arrays[sf_ind][:, 0]
        nfy = pt_var_arrays[sf_ind][:, 1]
        bnd_ind = (nfx > nodata) & (nfy > nodata)

    topxy = xy_pts[:, :][bnd_ind]
    x_order = np.argsort(topxy[:, 0])
    top_xy_sorted = topxy[x_order]
    nodal_flux_top_sorted = nodal_flux[bnd_ind][x_order]
    #nfy_top_sorted = nfy[bnd_ind][x_order]

    slope = np.diff(top_xy_sorted[:, 1]) / np.diff(top_xy_sorted[:, 0])
    slope_mean = slope.mean()

    # find h at top bnd
    h = pt_var_arrays[pt_var_names.index('h')]
    h_top = h[bnd_ind][x_order]

    print('rotate flux')
    nodal_flux_rotated = \
        get_normal_flux_to_bnd(nodal_flux_top_sorted,
                               top_xy_sorted)

    # correct for surpressed recharge at seepage bnd nodes
    est_rch_flux = nodal_flux_rotated[:, 1][-1]

    ind = (top_xy_sorted[:, 0] > 0) & (nodal_flux_rotated[:, 1] > (0.01 * est_rch_flux))
    nodal_flux_corrected = nodal_flux_rotated.copy()

    # convert flux to m/yr
    nfx *= year
    nfy *= year

    print('creating figs')

    show_watertable = False
    if show_watertable is True:
        print('calculating watertable position')

        p_ind = pt_var_names.index('pressure')
        pressure_array = pt_var_arrays[p_ind]

        xi = np.arange(xy_pts[:, 0].min(), xy_pts[:, 0].max() + dx, dx)
        yi = np.arange(xy_pts[:, 1].min(), xy_pts[:, 1].max() + dy, dy)
        xi_corr = np.arange(xy_pts[:, 0].min(), xy_pts[:, 0].max() + dx, dx)

        xy_pts_corr = xy_pts.copy()
        xy_pts_corr[:, 1] -= slope_mean * xy_pts_corr[:, 0]

        yi_corr = np.arange(xy_pts_corr[:, 1].min(), xy_pts_corr[:, 1].max() + dy, dy)

        xg, yg = np.meshgrid(xi, yi_corr)
        xgf, ygf = xg.flatten(), yg.flatten()

        # interpolate u to grid
        print('interpolating pressure to regular grid')
        pressure_array_regular = scipy.interpolate.griddata(xy_pts_corr,
                                                            pressure_array,
                                                            np.vstack((xgf, ygf)).T,
                                                            method='linear')

        # make a 2d grid again
        pressure_grid = np.resize(pressure_array_regular, xg.shape)
        ny, nx = xg.shape

        # find P=0 for each x
        print('finding point where P=0 for each column:')
        watertable_elevation = np.zeros(nx)
        watertable_x = xi.copy()
        for ix in range(nx):
            ind = np.isnan(pressure_grid[:, ix]) == False
            if np.min(pressure_grid[:, ix][ind]) > 0:
                # watertable exceeds surface, use approximate surface elevation
                # for watertable
                # TODO: find a way to get exact surface elevation...
                watertable_elevation[ix] = np.max(yi_corr[ind])
            else:
                # calculate position of watertable
                watertable_elevation[ix] = \
                    np.interp([0], pressure_grid[:, ix][ind][::-1],
                              yi_corr[ind][::-1])

        # only use if corrected for slope
        watertable_elevation = watertable_elevation + xi * slope_mean

    if conc_only is True:
        if 'concentration' in pt_var_names:
            c_ind = pt_var_names.index('concentration')
        elif 'conc' in pt_var_names:
            c_ind = pt_var_names.index('conc')

        pt_var_arrays_plot = [pt_var_arrays[c_ind]]
        pt_var_arrays_plot[0] = np.round(pt_var_arrays_plot[0], decimals=3)
        pt_var_names_plot = [pt_var_names[c_ind]]
    else:
        pt_var_arrays_plot = pt_var_arrays
        pt_var_names_plot = pt_var_names

    for nplot, pt_var_name, pt_var_array in zip(itertools.count(),
                                                pt_var_names_plot,
                                                pt_var_arrays_plot):

        print(pt_var_name)

        for vdim in range(len(pt_var_array.shape)):
            if len(pt_var_array.shape) > 1:
                # vector variable
                vtype = 'vector'
            else:
                vtype = 'scalar'

            if vtype is 'scalar':
                # calculate cell averages for model var
                cell_avg = np.array([pt_var_array[conn_i].mean()
                                     for conn_i in conn])
            else:
                cell_avg = np.array([pt_var_array[conn_i, vdim].mean()
                                     for conn_i in conn])

            if show_elements is True:
                # add element values:
                if nplot == 0:
                    elements = matplotlib.collections.PatchCollection(
                        patches,
                        linewidths=0.0,
                        antialiased=False)
                    ax.add_collection(elements)

                elements.set_array(cell_avg)
                elements.set_clim(pt_var_array.min(), pt_var_array.max())
                elements.set_edgecolor('None')
                elements.set_zorder(0.1)
                elements.set_antialiased(True)
                elements.set_cmap(cmap)

            else:
                if vtype == 'scalar':
                    if nplot == 0:
                        scatter_int = 1
                        leg_sc = ax.scatter(xy_pts[::scatter_int, 0],
                                            xy_pts[::scatter_int, 1],
                                            c=pt_var_array[::scatter_int],
                                            edgecolor='None', s=7)

                    else:
                        leg_sc.set_array(pt_var_array[::scatter_int])
                        if add_colorbar is True:
                            cb.set_clim(pt_var_array.min(),
                                        pt_var_array.max())

            xlim = [xy_pts[:, 0].min(), xy_pts[:, 0].max()]
            if ylim_fixed is None:
                ylim = [xy_pts[:, 1].min(), xy_pts[:, 1].max()]
            else:
                ylim = ylim_fixed

            if nplot == 0:
                ax.plot([xy_pts[:, 0].min(), 0], [0, 0],
                        color=sea_lvl_color, lw=0.25)

            if nplot == 0:
                va = (qx ** 2 + qy ** 2) ** 0.5
                scale = np.abs(va).max() * 10.0

                print('converting flux to regular grid')
                xg, yg, qxg, qyg = convert_to_grid(element_x, element_y,
                                                   qx * year,
                                                   qy * year, dx=dx, dy=dy)

                mask = np.zeros_like(qxg, dtype=bool)

                mask[np.isnan(qxg)] = True

                qxgm = np.ma.array(qxg, mask=mask)
                qygm = np.ma.array(qyg, mask=mask)
                vag = (qxgm ** 2 + qygm ** 2) ** 0.5
                lw = lw_min + lw_max * vag / vag.max()

                leg_splot = ax.streamplot(xg, yg, qxgm, qygm,
                                          density=streamline_density,
                                          color='k', linewidth=lw)

                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                ax.set_xlabel('Distance (m)')
                ax.set_ylabel('Elevation (m)')
                ax.grid()

            if nplot == 0:

                top_x = xy_pts[:, 0][bnd_ind]
                top_flux = nfy[bnd_ind]
                x_order = np.argsort(top_x)
                y0 = np.zeros_like(top_xy_sorted[:, 0])

                bnd_flux = \
                    tax.fill_between(top_xy_sorted[:, 0],
                                     y0,
                                     nodal_flux_corrected[:, 1],
                                     where=nodal_flux_corrected[:, 1] > 0,
                                     facecolor='lightblue',
                                     edgecolor='black',
                                     lw=0.5)
                bnd_flux2 = \
                    tax.fill_between(top_x[x_order],
                                     y0,
                                     nodal_flux_corrected[:, 1],
                                     where=nodal_flux_corrected[:, 1] < 0,
                                     facecolor='blue',
                                     edgecolor='black',
                                     lw=0.5)

                if debug is True:
                    if 'flux_surface_norm' in cell_var_names_face:
                        ind = cell_var_names_face.index('flux_surface_norm')
                        ind_sort = np.argsort(element_x_face)
                        bnd_flux3, = tax.plot(element_x_face[ind_sort],
                                              cell_var_arrays_face[ind][:, 1][ind_sort],
                                              color='green', lw=0.5, zorder=1001)

                    elif 'flux_surface_norm' in pt_var_names:
                        ind = pt_var_names.index('flux_surface_norm')
                        #ind_sort = np.argsort(element_x)
                        bnd_flux3 = tax.scatter(xy_pts[:, 0],
                                                pt_var_arrays[ind][:, 1],
                                                color='green', s=1, zorder=1001)

                # TODO show seepage bnd:
                if debug is True:
                    print('showing loc of seepage bnd')
                    ind = pt_var_names.index('active_seepage_bnd')
                    ind2 = pt_var_arrays[ind] > 0
                    xys = xy_pts[ind2]
                    zs = np.zeros((xys.shape[0]))
                    tax.scatter(xys[:, 0], zs, color='brown', s=2, zorder=1003)

                tax.set_xlim(xlim)
                tax.grid()
                tax.axhline(y=0, color='black', lw=0.5)
                tax.axvline(x=0, color='black', lw=0.5)
                tax.set_ylabel('Boundary flux (m/yr)')

                if add_title is True:
                    tax.set_title(plot_title, fontsize='small')

                if tax_ylim_fixed is not None:
                    tax.set_ylim(tax_ylim_fixed)
                else:
                    tax.set_ylim(nodal_flux_corrected[:, 1].min() * 1.25,
                                 nodal_flux_corrected[:, 1].max() * 1.25)

                if add_colorbar is True:
                    if show_elements is True:
                        cb = fig.colorbar(elements,
                                          cax=cax2,
                                          orientation='horizontal')
                    else:
                        cb = fig.colorbar(leg_sc,
                                          cax=cax2,
                                          orientation='horizontal')
                    cb.ax.tick_params(labelsize=leg_fontsize)

                    # (generate plot here)
                    tick_locator = ticker.MaxNLocator(nbins=3)
                    cb.locator = tick_locator
                    cb.update_ticks()

            if add_colorbar is True:
                if 'concentration' in pt_var_name:
                    cb.set_label('Solute concentration (kg/kg)'
                                 , fontsize='small')

                else:
                    cb.set_label(pt_var_name)

            if nplot == 0:
                leg_wt2, = ax.plot(top_xy_sorted[:, 0], h_top,
                                   color='gray', lw=1.5)

                if debug is True:
                    print('hydraulic head at top nodes: ', h_top)
                # add vertical line at coastline
                ax.axvline(x=0, color='black', lw=0.5, zorder=0)

            if nplot == 0 and add_legend is True:

                legs = [bnd_flux2, bnd_flux, leg_wt2]
                labels = ['recharge', 'discharge', 'watertable']
                tax.legend(legs, labels, frameon=False, loc='upper right', fontsize=leg_fontsize)

            # save regular figure:
            file_basename = ''.join(vtk_file.split('.')[:-1])
            if vtype is 'scalar':
                fig_fn = file_basename + '_' + pt_var_name + figure_type
            else:
                fig_fn = file_basename + '_' + pt_var_name \
                         + ['x', 'y', 'z'][vdim] + figure_type


            fn_out = os.path.join(fig_folder, fig_fn)
            #fn_out_alt = os.path.join(fig_folder_alt, fig_fn)
            print('saving ', fn_out)
            fig.savefig(fn_out, dpi=dpi)
            #print('saving ', fn_out_alt)
            #fig.savefig(fn_out_alt, dpi=dpi)

pl.clf()
print('done')