from __future__ import print_function


"""

"""

import os
import numpy as np

import read_vtu_file
from grompy_lib import get_normal_flux_to_bnd


def get_vtk_files(base_dir, argv=[''], folder=None):

    if 'vtu' in argv[-1]:

        vtk_file = argv[-1]

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
        selection = raw_input()

        if selection.isdigit() is True:
            folder = dirs[int(selection)]
        else:
            folder = dirs[0]
        folder = os.path.join(folder, 'vtk_files')

    else:
        folder = os.path.join(folder, 'vtk_files')

    print('enter f to make a figure for final results only or enter to include all timesteps')
    selection = raw_input()
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
    print('output files, from newest to oldest:')
    for i, vtk_file in enumerate(vtk_files):
        print(i, os.path.split(vtk_file)[-1])

    print('\nenter number of output file to show, two numbers separated by -, ' 
          'a for all files, or r for a random selection of 10 files or enter for the newest file:')
    selection = raw_input()

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

    return vtk_files, vtkf_files


def extract_bnd_fluxes(fn, fnf):
    """
    extract location of top boundary and top bnd fluxes from a escript/grompy vtk file

    Parameters
    ----------
    fn : string
        filename of vtk element file
    fnf : string
        filename of vtk FaceElement file

    Returns
    -------
    bnd_xy : numpy array
        Nx3 array containing the x,y,z coordinates fo each node
    bnd_flux : numpy array
        array containing the fluxes over the top boundary for each node
    """

    print('reading vtk file %s' % fn)

    xy_pts, conn, pt_var_names, pt_var_arrays, cell_var_names, cell_var_arrays = \
        read_vtu_file.read_vtu_file(fn)

    #fnf = os.path.join(folder, vtkf_file)
    #fnf = vtkf_file

    print('reading vtk file %s' % fnf)

    xy_pts_face, conn_face, pt_var_names_face, pt_var_arrays_face, cell_var_names_face, cell_var_arrays_face = \
        read_vtu_file.read_vtu_file(fnf)

    #vtk_file = os.path.split(vtk_file)[-1]
    #vtkf_file = os.path.split(vtkf_file)[-1]

    nodata = -99999
    if 'surface' in pt_var_names:
        sf_ind = pt_var_names.index('surface')
        bnd_ind = pt_var_arrays[sf_ind] == 1
    else:
        sf_ind = pt_var_names.index('nodal_flux_surface')
        nfx = pt_var_arrays[sf_ind][:, 0]
        nfy = pt_var_arrays[sf_ind][:, 1]
        bnd_ind = (nfx > nodata) & (nfy > nodata)

    # find top bnd fluxes
    print('find top bnd fluxes')
    if 'nodal_flux_surface' in pt_var_names:
        f_ind = pt_var_names.index('nodal_flux')
    nfx = pt_var_arrays[f_ind][:, 0]
    nfy = pt_var_arrays[f_ind][:, 1]
    nodal_flux = pt_var_arrays[f_ind]

    #
    year = 365.24 * 24 * 3600
    nodal_flux *= year

    topxy = xy_pts[:, :][bnd_ind]
    x_order = np.argsort(topxy[:, 0])
    top_xy_sorted = topxy[x_order]
    nodal_flux_top_sorted = nodal_flux[bnd_ind][x_order]
    #nfy_top_sorted = nfy[bnd_ind][x_order]

    nodal_flux_rotated = get_normal_flux_to_bnd(nodal_flux_top_sorted, top_xy_sorted)

    #bnd_xys.append(top_xy_sorted)
    #bnd_fluxes.append(nodal_flux_rotated[:, 1])
    bnd_xy, bnd_flux = top_xy_sorted, nodal_flux_rotated[:, 1]

    return bnd_xy, bnd_flux