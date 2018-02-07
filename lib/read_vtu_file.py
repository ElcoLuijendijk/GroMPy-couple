"""
Read a a VTK (.vtu) file of an unstructured grid and return arrays of the variables
"""

__author__ = 'Elco Luijendijk'


import numpy as np


def find_name_after_index(inp_str, index_str):

    pt_ln2 = inp_str[inp_str.find(index_str)+len(index_str):]
    pt_ln3 = pt_ln2.split()[0]
    if pt_ln3.startswith('"'):
        inp_str = pt_ln3.split('"')[1]
    else:
        inp_str = pt_ln3

    return inp_str


def find_number_after_index(inp_str, index_str, dtype=int):

    pt_ln2 = inp_str[inp_str.find(index_str)+len(index_str):]
    pt_ln3 = pt_ln2.split()[0]
    if pt_ln3.startswith('"'):
        inp_nr = dtype(pt_ln3.split('"')[1])
    else:
        inp_nr = dtype(pt_ln3)

    return inp_nr


def find_line_index(inp_list, search_text):

    line_ind = [i for i, l in enumerate(inp_list) if search_text in l][0]

    return line_ind


def find_line(inp_list, search_text):

    inp_str = inp_list[find_line_index(inp_list, search_text)]

    return inp_str


def find_pointdata(di, npts):

    pt_var_names = []
    pt_var_arrays = []

    go = True

    ind_counter = 0

    if '<PointData>' not in di:
        return pt_var_names, pt_var_arrays

    pt_data_ind = find_line_index(di, '<PointData>') + 1
    pt_data_ind_end = find_line_index(di, '</PointData>')

    while go is True and pt_data_ind < pt_data_ind_end - 1:

        pt_data_ind = pt_data_ind + find_line_index(di[pt_data_ind:], '<DataArray')

        var_name = find_name_after_index(di[pt_data_ind], 'Name=')

        print 'found var name: ', var_name

        ncomponents = find_number_after_index(di[pt_data_ind], 'NumberOfComponents=')

        # check if variable is scalar or vector:
        if len(di[pt_data_ind+1].split()) > 1:
            # vector variable
            var_list = [dii.split() for dii in di[pt_data_ind+1:pt_data_ind+1+npts]]
            var_array = np.array(var_list, dtype=float)
        else:
            # scalar variable
            var_array = np.array(di[pt_data_ind+1:pt_data_ind+1+npts], dtype=float)

        pt_data_ind = pt_data_ind + find_line_index(di[pt_data_ind:], '</DataArray')

        #
        pt_var_names.append(var_name)
        pt_var_arrays.append(var_array)

    return pt_var_names, pt_var_arrays


def find_celldata(di, ncells):

    cell_var_names = []
    cell_var_arrays = []

    go = True

    ind_counter = 0

    cell_data_ind = find_line_index(di, '<CellData>') + 1
    cell_data_ind_end = find_line_index(di, '</CellData>')

    while go is True and cell_data_ind < cell_data_ind_end - 1:

        cell_data_ind = cell_data_ind + find_line_index(di[cell_data_ind:], '<DataArray')

        var_name = find_name_after_index(di[cell_data_ind], 'Name=')

        print 'found var name: ', var_name

        ncomponents = find_number_after_index(di[cell_data_ind], 'NumberOfComponents=')

        di_sel = di[cell_data_ind+1:cell_data_ind+1+ncells]
        dis = [dii.split() for dii in di_sel]

        var_array = np.array(dis, dtype=float)

        # check if 2nd dimension = len 1 and simplify array if yet
        if var_array.shape[1] == 1:
            var_array = var_array[:, 0]

        cell_data_ind = cell_data_ind + find_line_index(di[cell_data_ind:], '</DataArray')

        #
        cell_var_names.append(var_name)
        cell_var_arrays.append(var_array)

    return cell_var_names, cell_var_arrays


def read_vtu_file(fn):
    """
    """

    fin = open(fn, 'r')
    din = fin.readlines()
    fin.close()

    di = [dini.strip() for dini in din]

    # find number of points and cells

    #
    pt_ln = find_line(di, '<Piece')
    npts = find_number_after_index(pt_ln, 'NumberOfPoints=')
    ncells = find_number_after_index(pt_ln, 'NumberOfCells=')

    # find node locations
    xy_pt_index = find_line_index(di, '<Points>') + 2
    xyl = di[xy_pt_index:xy_pt_index+npts]
    xyl = [xyli.split() for xyli in xyl]
    xy_pts = np.array(xyl, dtype=float)

    # find connectivity data
    cell_index = find_line_index(di, '<Cells>') + 2
    cell_conn = di[cell_index:cell_index+ncells]
    cell_conn = [cell_conni.split() for cell_conni in cell_conn]
    conn = np.array(cell_conn, dtype=int)

    pt_var_names, pt_var_arrays = find_pointdata(di, npts)

    cell_var_names, cell_var_arrays = find_celldata(di, ncells)

    return xy_pts, conn, pt_var_names, pt_var_arrays, cell_var_names, cell_var_arrays
