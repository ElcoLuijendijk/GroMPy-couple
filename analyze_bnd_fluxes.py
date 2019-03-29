"""
Read a series of Grompy model result VTK files from a single directory, analyze fluxes over upper boundary and save
the results to a .csv file.

The script will ask you to select a subdirectory where the VTK files are located. It will search for subdirectories
in the directory model_output
"""

__author__ = 'Elco Luijendijk'

import matplotlib
matplotlib.use('Agg')

import os
import re
import sys

import numpy as np
import matplotlib.pyplot as pl
import pandas as pd

import lib.read_vtu_file
from lib.grompy_lib import get_normal_flux_to_bnd, get_distance

year = 365.25 * 24 * 60 * 60

#default_dir = '/home/elco/model_files/grompy_sgd_models'
default_dir = 'model_output'

#folder = '/home/elco/model_files/grompy_sgd_models/base_case_highk/vtk_files'
#folder = None
#model_result_file = '/home/elco/model_files/grompy_sgd_models/combinations_final_merged/final_model_results_merged.csv'
#model_result_file = '/home/elco/model_files/grompy_sgd_models/combinations_final_topo2_jul17/final_model_results_merged.csv'
#model_result_file = '/home/elco/model_files/grompy_sgd_models/combinations_final_v2_jul17/final_model_results_merged.csv'

if '.csv' in sys.argv[-1]:
    model_result_file = sys.argv[-1]
else:
    model_result_file = None

if model_result_file is None:

    files = os.listdir(default_dir)
    files = [os.path.join(default_dir, directory) for directory in files]
    files.sort(key=os.path.getmtime)
    for i, f in enumerate(files):
        print i, f

    print 'select directory number'
    #folder = os.path.join(default_dir, files[int(raw_input())])
    folder = files[int(raw_input())]

    files = os.listdir(folder)
    files = [os.path.join(folder, f) for f in files if f[-4:] == '.csv']
    files.sort()
    for i, f in enumerate(files):
        print i, os.path.split(f)[-1]
    print 'select file number'
    model_result_file = files[int(raw_input())]


folder = os.path.split(model_result_file)[0]

df = pd.read_csv(model_result_file)

df = df.dropna(subset=['vtk_filename'])

vtk_files_raw = df['vtk_filename'].tolist()

vtk_files = [os.path.join('vtk_files', os.path.split(v)[-1][:-4] + '_Elements.vtu') for v in vtk_files_raw]

fig_folder = os.path.join(folder, 'fig')

if os.path.exists(fig_folder) is False:
    os.mkdir(fig_folder)

# add columns
newcols = ['file_number',
           'model_name',
           'vtk_flux_total',
           'vtk_flux_sea_total',
           'vtk_flux_sea_in',
           'vtk_flux_sea_out',
           'vtk_flux_land_total',
           'vtk_flux_land_in',
           'vtk_flux_land_out',
           'vtk_min_flux_sea',
           'vtk_max_flux_sea',
           'vtk_min_flux_land',
           'vtk_max_flux_land',
           'vtk_flux_sea_min_x',
           'vtk_flux_sea_max_x',
           'vtk_flux_out_sea_xmin',
           'vtk_flux_out_sea_xmax',
           'vtk_flux_out_land_xmin',
           'vtk_flux_out_land_xmax']

model_nos = df.index

for fileno, vtk_file in zip(model_nos, vtk_files):

    print 'reading vtk file %s' % vtk_file

    fnv = os.path.join(folder, vtk_file)
    xy_pts, conn, pt_var_names, pt_var_arrays, cell_var_names, cell_var_arrays = \
        lib.read_vtu_file.read_vtu_file(fnv)

    #fnf = os.path.join(folder, vtkf_file)
    #xy_pts_face, conn_face, pt_var_names_face, pt_var_arrays_face, cell_var_names_face, cell_var_arrays_face = \
    #    lib.read_vtu_file.read_vtu_file(fnf)

    use_old_method = False

    if use_old_method is True:
        # get fluxes
        qx_ind = cell_var_names.index('qx')
        qx = cell_var_arrays[qx_ind]
        qy_ind = cell_var_names.index('qy')
        qy = cell_var_arrays[qy_ind]

        # find top bnd fluxes
        if 'nodal_flux_surface' in pt_var_names:
            sf_ind = pt_var_names.index('nodal_flux_surface')
        bfx = pt_var_arrays[sf_ind][:, 0]
        bfy = pt_var_arrays[sf_ind][:, 1]

        nodata = -99999
        bnd_ind = (bfx > nodata) & (bfy > nodata)

        bfx = bfx[bnd_ind]
        bfy = bfy[bnd_ind]

        # convert flux to m/yr
        bfx *= year
        bfy *= year

        bxy = xy_pts[bnd_ind]
        #bxy = bxy[:2]
        bx = bxy[:, 0]
        by = bxy[:, 1]

        # sort by x-coordinate
        xsort = np.argsort(bx)

        bx = bx[xsort]
        by = by[xsort]
        bfx = bfx[xsort]
        bfy = bfy[xsort]

        dxi = np.diff(bx) / 2.0
        dx = np.zeros_like(bx)
        dx[:-1] += dxi
        dx[1:] += dxi

        nodal_flux_top_int = bfy * dx

        ind_sea = bx < 0
        ind_sea_out = ind_sea * bfy > 0
        ind_sea_in = ind_sea * bfy < 0

        ind_land = bx >= 0
        ind_land_out = ind_land * bfy > 0
        ind_land_in = ind_land * bfy < 0

    # new method from figure script:
    print 'find top bnd fluxes'
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

    print 'rotate flux'
    nodal_flux_top = \
        get_normal_flux_to_bnd(nodal_flux_top_sorted,
                               top_xy_sorted)

    dist = get_distance(top_xy_sorted)
    nodal_flux_top_int = nodal_flux_top[:, 1] * dist

    ind_sea = top_xy_sorted[:, 0] <= 0
    ind_land = top_xy_sorted[:, 0] > 0

    flux_in = nodal_flux_top[:, 1] < 0
    flux_out = nodal_flux_top[:, 1] > 0

    ind_sea_in = ind_sea & flux_in
    ind_sea_out = ind_sea & flux_out

    ind_land_in = ind_land & flux_in
    ind_land_out = ind_land & flux_out

    #####################################################################################################
    # new: calculate bnd fluxes based on fluxes in elements, and not fluxes projected to nodes by escript
    #####################################################################################################
    qx_ind = cell_var_names.index('qx')
    qx = cell_var_arrays[qx_ind]
    qy_ind = cell_var_names.index('qy')
    qy = cell_var_arrays[qy_ind]

    qx_yr = qx * year
    qy_yr = qy * year

    # find boundary nodes:
    nodata = -99999
    if 'surface' in pt_var_names:
        sf_ind = pt_var_names.index('surface')
        bnd_ind = pt_var_arrays[sf_ind] == 1
    else:
        sf_ind = pt_var_names.index('nodal_flux_surface')
        nfx = pt_var_arrays[sf_ind][:, 0]
        nfy = pt_var_arrays[sf_ind][:, 1]
        bnd_ind = (nfx > nodata) & (nfy > nodata)

    # find elements with two boundary nodes:
    bnd_ele = np.zeros_like(qx, dtype=bool)
    bnd_ele[:] = False
    conn_bnd = np.zeros_like(conn, dtype=bool)
    conn_bnd[:] = False
    for i, conni in enumerate(conn):
        bnd_ele_cnt = 0
        for j, conn_ii in enumerate(conni):
            if bnd_ind[conn_ii] == True:
                bnd_ele_cnt += 1
                conn_bnd[i, j] = True
                #print 'found bnd ele at element %i and node %i' % (i, conn_ii)
        if bnd_ele_cnt >= 2:
            bnd_ele[i] = True

    # get position of elements at top bnd
    top_xy_bnd_ele_list = [xy_pts[conni[conn_bnd_i], :2]
                           for conni, conn_bnd_i, bnd_ele_i in zip(conn, conn_bnd, bnd_ele)
                           if bnd_ele_i == True]
    top_xy_bnd_ele_nodes = np.array(top_xy_bnd_ele_list)
    top_xy_ele_unsorted = np.mean(top_xy_bnd_ele_nodes, axis=1)

    # sort nodes with increasing x
    xe_order = np.argsort(top_xy_ele_unsorted[:, 0])
    top_xy_ele = top_xy_ele_unsorted[xe_order]
    qx_top_sorted = qx_yr[bnd_ele][xe_order]
    qy_top_sorted = qy_yr[bnd_ele][xe_order]

    flux_top_sorted = np.array([qx_top_sorted, qy_top_sorted]).T

    print 'rotate flux'
    flux_top = \
        get_normal_flux_to_bnd(flux_top_sorted,
                               top_xy_ele)

    dist_ele = get_distance(top_xy_ele)
    flux_top_int = flux_top[:, 1] * dist_ele

    # get flux normal to boundary:
    ind_sea_ele = top_xy_ele[:, 0] <= 0
    ind_land_ele = top_xy_ele[:, 0] > 0

    flux_in_ele = flux_top[:, 1] < 0
    flux_out_ele = flux_top[:, 1] > 0

    ind_sea_in_ele = ind_sea_ele & flux_in_ele
    ind_sea_out_ele = ind_sea_ele & flux_out_ele

    ind_land_in_ele = ind_land_ele & flux_in_ele
    ind_land_out_ele = ind_land_ele & flux_out_ele

    ####################################
    # store results in pandas dataframe
    ####################################
    # store filename
    df.ix[fileno, 'model_name'] = vtk_file

    # calculate fluxes

    # total
    df.ix[fileno, 'vtk_flux_in'] = np.sum(flux_top_int[flux_in_ele])
    df.ix[fileno, 'vtk_flux_out'] = np.sum(flux_top_int[flux_out_ele])

    # sea
    df.ix[fileno, 'vtk_flux_sea_in'] = np.sum(flux_top_int[ind_sea_in_ele])
    df.ix[fileno, 'vtk_flux_sea_out'] = np.sum(flux_top_int[ind_sea_out_ele])
    df.ix[fileno, 'vtk_flux_sea_total'] = np.sum(flux_top_int[ind_sea_ele])
    df.ix[fileno, 'vtk_min_flux_sea'] = np.min(nodal_flux_top[:, 1][ind_sea])
    df.ix[fileno, 'vtk_max_flux_sea'] = np.max(nodal_flux_top[:, 1][ind_sea])

    # land
    df.ix[fileno, 'vtk_flux_land_in'] = np.sum(flux_top_int[ind_land_in_ele])
    df.ix[fileno, 'vtk_flux_land_out'] = np.sum(flux_top_int[ind_land_out_ele])
    df.ix[fileno, 'vtk_flux_land_total'] = np.sum(flux_top_int[ind_land_ele])
    df.ix[fileno, 'vtk_min_flux_land'] = np.min(flux_top[:, 1][ind_land_ele])
    df.ix[fileno, 'vtk_max_flux_land'] = np.max(flux_top[:, 1][ind_land_ele])
    df.ix[fileno, 'vtk_flux_land_right_hand_bnd'] = flux_top[:, 1][ind_land_ele][-1]

    df.ix[fileno, 'vtk_bnd_flux_left'] = flux_top[0, 1]
    df.ix[fileno, 'vtk_bnd_flux_right'] = flux_top[-1, 1]

    df.ix[fileno, 'vtk_submarine_discharge_area'] = np.sum(dist * ind_sea_out)
    df.ix[fileno, 'vtk_submarine_recharge_area'] = np.sum(dist * ind_sea_in)

    df.ix[fileno, 'vtk_land_discharge_area'] = np.sum(dist * ind_land_out)
    df.ix[fileno, 'vtk_land_recharge_area'] = np.sum(dist * ind_land_in)

    # calculate area where absolute flux exceeds 0.1 m/yr or 50% of the max flux
    flux_thresholds = [0.1, 0.5 * np.max(np.abs(flux_top[:, 1]))]

    for i, flux_threshold in enumerate(flux_thresholds):
        flux_in_threshold = flux_top[:, 1] < -flux_threshold
        flux_out_threshold = flux_top[:, 1] > flux_threshold

        ind_sea_in_threshold = ind_sea_ele & flux_in_threshold
        ind_sea_out_threshold = ind_sea_ele & flux_out_threshold

        ind_land_in_threshold = ind_land_ele & flux_in_threshold
        ind_land_out_threshold = ind_land_ele & flux_out_threshold

        df.ix[fileno, 'vtk_submarine_discharge_area_threshold_50perc_max'] = \
            np.sum(dist_ele * ind_sea_out_threshold)
        df.ix[fileno, 'vtk_land_discharge_area_threshold_50perc_max'] = \
            np.sum(dist_ele * ind_land_out_threshold)

    # calculate normalized cumulative flux
    cumulative_flux_top_land = np.cumsum(flux_top_int[ind_land_out_ele]) / np.sum(flux_top_int[ind_land_out_ele])
    cumulative_flux_top_sea = np.cumsum(flux_top_int[ind_sea_out_ele][::-1]) / np.sum(flux_top_int[ind_sea_out_ele])
    cumulative_flux_top_sea = cumulative_flux_top_sea[::-1]

    # calculate area where x % of the total terrestrial or marine discharge takes place:
    flux_thresholds2 = [0.25, 0.5, 0.75, 0.9, 0.95, 0.99]

    for flux_threshold2 in flux_thresholds2:
        try:
            ind_landc = np.where(cumulative_flux_top_land > flux_threshold2)[0][0]
            df.ix[fileno, 'vtk_land_discharge_area_%0.0fperc' % (flux_threshold2 * 100.0)] = \
                np.sum(dist_ele[ind_land_out_ele][:ind_landc])
        except:
            df.ix[fileno, 'vtk_land_discharge_area_%0.0fperc' % (flux_threshold2 * 100.0)] = 0.0

        try:
            ind_seac = np.where(cumulative_flux_top_sea > flux_threshold2)[0][-1]
            df.ix[fileno, 'vtk_sea_discharge_area_%0.0fperc' % (flux_threshold2 * 100.0)] = \
                np.sum(dist_ele[ind_sea_out_ele][ind_seac:])
        except:
            df.ix[fileno, 'vtk_sea_discharge_area_%0.0fperc' % (flux_threshold2 * 100.0)] = 0.0

    # flux from land to sea (= fresh submarine groundwater discharge)
    df.ix[fileno, 'vtk_flux_land_to_sea'] = df.ix[fileno, 'vtk_flux_sea_out'] \
            + df.ix[fileno, 'vtk_flux_sea_in']

    #
    df.ix[fileno, 'vtk_flux_total'] = np.sum(flux_top_int)

    # relative error in top bnd flux,
    df.ix[fileno, 'vtk_flux_error'] = np.abs(np.sum(flux_top_int)) / (np.abs(np.sum(flux_top_int[flux_top_int<0]))
                                                                       + np.abs(np.sum(flux_top_int[flux_top_int>0])))

    # add max extent seepage bnd:
    if 'active_seepage_bnd' in pt_var_names:
        si = pt_var_names.index('active_seepage_bnd')
        s = pt_var_arrays[si]
        if s.max() >= 1:
            xys = xy_pts[s == 1]
            df.ix[fileno, 'max_x_active_seepage_bnd'] = xys[:, 0].max()
        else:
            df.ix[fileno, 'max_x_active_seepage_bnd'] = np.nan

    print 'flux sea out (m2/yr) = ', df.ix[fileno, 'vtk_flux_sea_out']
    print 'flux sea in (m2/yr) = ', df.ix[fileno, 'vtk_flux_sea_in']
    print 'flux land to sea (m2/yr) = ', df.ix[fileno, 'vtk_flux_land_to_sea']
    print 'flux error (%) = ', df.ix[fileno, 'vtk_flux_error'] * 100.0

# join input and fluxes dataframe
df_final = df

fn_merged = model_result_file.split('.')[0] + '_with_bnd_stats.csv'

print 'saving input params and bnd fluxes to %s' % fn_merged
df_final.to_csv(fn_merged, index_label='model_run_no_v2')

print 'done'