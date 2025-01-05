

"""
supporting functions for Grompy

"""

import os
import math
import itertools
import pdb

import numpy as np

import lib.gwflow_lib as gwflow_lib

try:
    import esys.escript as es
    import esys.weipa
except ImportError:
    print('warning, could not find escript module')


def get_timestr(sec, limit_to_days=True):

    """
    Construct a string with the date and time, used for naming output files.
    """

    second = 1.0
    minute = 60.0
    hour = 60.0 * minute
    day = hour * 24
    year = day * 365.25

    years = math.floor(sec/year)
    sec_left = sec - years * year
    days = math.floor(sec_left / day)
    sec_left -= days * day
    hours = math.floor(sec_left / hour)
    sec_left -= hours * hour
    minutes = math.floor(sec_left / minute)
    seconds = sec_left - minutes * minute

    if limit_to_days is True:
        timestr = '%i yrs, %i days' % \
                  (years, days)
    else:
        timestr = '%i yrs, %i days, %i hrs, %i min, %i sec' % \
                  (years, days, hours, minutes, seconds)

    return timestr


def convert_to_array(u):

    """
    Return the x,y coordinates and the value of escript variable u as a numpy array.
    """

    coords = u.getFunctionSpace().getX()

    xy = np.array([coords[0].toListOfTuples(), coords[1].toListOfTuples()]).T

    u_array = np.array(u.toListOfTuples())

    #assert len(u_array.shape) == 1

    return xy, u_array


def calculate_angle(dx, dy):
    """
    Calculate the angle between two points.
    """

    degrees = np.arctan2(dy, dx) * 180 / np.pi

    return degrees


def get_variable_xy(variable):
    """
    Get x and y coordinate of an escript variable.
    """
    xy = variable.getFunctionSpace().getX()
    xya = np.array([xy[0].toListOfTuples(),
                    xy[1].toListOfTuples()]).T

    return xya


def get_sorted_flux_arrays(flux, bnd):
    """
    Return numpy arrays of flux and x,y coordinates of flux sorted by x-coordinate.
    """

    # calculate total flux over model boundary
    flux_xy = get_variable_xy(flux)
    flux_array = np.array([flux[0].toListOfTuples(),
                           flux[1].toListOfTuples()]).T

    # sort in order of increasing x
    #sort_order = np.lexsort((flux_xy[:, 0], flux_xy[:, 1]))
    #flux_array_sorted = flux_array[sort_order]
    #flux_xy_sorted = flux_xy[sort_order]

    bnd_array = np.array(bnd.toListOfTuples())
    bnd_xy = get_variable_xy(bnd)

    bnd_xy_filtered = bnd_xy[bnd_array > 0]

    # find if values are in bnd_array
    ind = np.array([(np.any(bnd_xy_filtered[:, 0] == xy[0]) &
                    np.any(bnd_xy_filtered[:, 1] == xy[1]))
                    for xy in flux_xy])

    flux_xy_bnd = flux_xy[ind]
    flux_array_bnd = flux_array[ind]

    sort_order = np.argsort(flux_xy_bnd[:, 0])

    return flux_array_bnd[sort_order], flux_xy_bnd[sort_order]


def get_dx_and_dy(xy):
    """
    Calculate dx and dy difference of a (Nx2) array with x and y coordinates.
    """

    if xy.shape[0] <= 1:
        return np.zeros(1)

    dxdy = np.diff(xy, axis=0)
    dxdy_ext = np.zeros((dxdy.shape[0] + 1, dxdy.shape[1]))
    dxdy_ext[1:-1] = (dxdy[1:] + dxdy[:-1]) / 2.0
    dxdy_ext[0] = dxdy[0]
    dxdy_ext[-1] = dxdy[-1]

    return dxdy_ext


def rotate_vector_np(v, angle):
    """
    Rotate a vector `v` by the given angle in radians
    """

    v_rot = v.copy()

    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)

    v_rot[:, 0] = v[:, 0] * cos_theta - v[:, 1] * sin_theta
    v_rot[:, 1] = v[:, 0] * sin_theta + v[:, 1] * cos_theta

    return v_rot


def get_normal_flux_to_bnd(flux_array, flux_bnd_xy):
    """
    Decompose flux vector into components normal and parallel to topography.
    """

    # find perpendicular vector for topo gradient
    # dxdy = np.diff(flux_bnd_xy, axis=0)

    # interpoalte to neighboring nodes again
    dxdy_ext = get_dx_and_dy(flux_bnd_xy)

    # rotate flux
    angle_topo = calculate_angle(dxdy_ext[:, 0], dxdy_ext[:, 1])

    flux_array_rot = rotate_vector_np(flux_array, -np.radians(angle_topo))

    return flux_array_rot


def get_distance(xy):
    """
    Calculate distance between series of points.
    """

    dist = np.zeros((xy.shape[0]))

    dxdy = np.diff(xy, axis=0)

    dist_diff = np.sqrt(dxdy[:, 0] ** 2 + dxdy[:, 1] ** 2)

    dist[1:] = 0.5 * dist_diff
    dist[:-1] += 0.5 * dist_diff

    return dist


def get_boundary_flux_simple(flux, bnd):
    """
    Calculate the flux normal to the model boundary.
    """

    # get flux arrays
    flux_array, flux_xy = get_sorted_flux_arrays(flux, bnd)

    if flux_array.shape[0] <= 1:
        return flux_xy, None

    # get flux normal and parallel to topography
    flux_normal = \
        get_normal_flux_to_bnd(flux_array, flux_xy)

    return flux_xy, flux_normal


def set_boundary_conditions(mesh, surface,  z_surface, Parameters):

    """
    Set boundary conditions.
    
    Parameters
    ----------
    mesh :
        escript mesh
    rho_f_0 : float
        fluid density at solute concentration = 0
    
    gamma : float
        coefficient of density increase
    g : float
        gravitational constant
    initial_condition_parameters : list
    
    Returns
    -------
    specified_pressure : escript variable
        specified pressure
    specified_pressure_bnd : escript variable
        location of specified pressure boundary
        (1=spec pressure, 0 is no spec pressure)
    specified_concentration : escript variable
        specified concentration 
    specified_concentration_bnd : escript variable
        lcoation of specified concentration boundary
    rch_bnd_loc : escript variable
        location of recharge boundary
    drain_bnd_loc : escript variable
        location of drain boundary
        
    """

    xy = mesh.getX()

    #rho_f_salt = Parameters.seawater_concentration * Parameters.gamma + \
    #    Parameters.rho_f_0
    #rho_f_salt = gwflow_lib.calculate_fluid_density(
    #    Parameters.seawater_concentration,
    #    Parameters.gamma,
    #    Parameters.rho_f_0)
    #rho_f_fresh = Parameters.freshwater_concentration * Parameters.gamma +
    # Parameters.rho_f_0

    # drain bnd location
    drain_bnd_loc = es.whereNegative(Parameters.drain_bnd_xmin - xy[0]) * \
        es.wherePositive(Parameters.drain_bnd_xmax - xy[0]) * \
        es.wherePositive(surface)

    # recharge bnd location
    try:
        # use tags to set bnd
        rch_bnd_loc = es.Scalar(0, es.FunctionOnBoundary(mesh))
        rch_bnd_loc.setTaggedValue("land_surface1", 1)
        rch_bnd_loc.setTaggedValue("land_surface2", 1)

        xy_rch = rch_bnd_loc.getFunctionSpace().getX()

        # make sure no recharge at sea or coastal node:
        rch_bnd_loc = rch_bnd_loc * es.wherePositive(xy_rch[0])

    except RuntimeError:
        print('could not set recharge bnd using tags, using location instead')
        rch_bnd_loc = (es.whereNegative(Parameters.recharge_mass_flux_xmin -
                                        xy[0])
                       * es.wherePositive(Parameters.recharge_mass_flux_xmax -
                                          xy[0])
                       * es.wherePositive(surface))

    rch_xy, rch_array = convert_to_array(rch_bnd_loc)
    print('number of active recharge nodes = %i' % rch_array.sum())

    #
    if Parameters.specified_concentration_surface is True:
        print('assigning specified concentration only to the surface up to x=0:')
        #specified_concentration_bnd = es.whereNegative(xy[0]) * surface
        #specified_concentration_bnd = es.whereNegative(Parameters.specified_concentration_xmin - xy[0]) * \
        #    es.wherePositive(Parameters.specified_concentration_xmax - xy[0]) * \
        #    es.wherePositive(surface)
        
        #xy_spec_conc = specified_concentration_bnd.getFunctionSpace().getX()

        # assign to user specified section:
        specified_concentration_bnd = surface * es.whereNegative(Parameters.specified_concentration_xmin - xy[0]) * \
            es.wherePositive(Parameters.specified_concentration_xmax - xy[0])

        
        specified_concentration = (es.whereNegative(xy[0]) * Parameters.specified_concentration
                                   + es.whereNonNegative(xy[0]) * Parameters.specified_concentration)
        
        spec_conc_xy, spec_conc_array = convert_to_array(specified_concentration_bnd)
        print(f"number of active specified concentration nodes = {spec_conc_array.sum()}")
        
    else:

        specified_concentration_bnd = surface * 0
        specified_concentration = surface * 0
        for xmin, xmax, ymin, ymax, spec_conc_segment in zip(Parameters.specified_concentration_xmin,
                                                             Parameters.specified_concentration_xmax,
                                                             Parameters.specified_concentration_ymin,
                                                             Parameters.specified_concentration_ymax,
                                                             Parameters.specified_concentration):

            print('assigning specified concentration of %0.4f at surface from x=%0.2f to %0.2f, y=%0.2f to %0.2f' %
                  (spec_conc_segment, xmin, xmax, ymin, ymax))

            specified_concentration_bnd_segment = (es.whereNonPositive(xmin - xy[0]) *
                                                   es.whereNonNegative(xmax - xy[0]) *
                                                   es.whereNonPositive(ymin - xy[1]) *
                                                   es.whereNonNegative(ymax - xy[1]))
            specified_concentration_bnd += specified_concentration_bnd_segment
            specified_concentration += (specified_concentration_bnd_segment * spec_conc_segment)

    # calculate product of concentration and density for bnd nodes
    specified_concentration_rho = \
        gwflow_lib.calculate_fluid_density(specified_concentration,
                                           Parameters.gamma,
                                           Parameters.rho_f_0)
    specified_concentration_rho_f = (specified_concentration *
                                     specified_concentration_rho)
    # set specified pressure boundary:
    print('setting specified pressure nodes at surface')
    if Parameters.specified_pressure_surface is True:
    
        specified_pressure_bnd = \
            es.whereNegative(Parameters.specified_pressure_xmin - xy[0]) * es.wherePositive(Parameters.specified_pressure_xmax - xy[0]) * es.wherePositive(surface)

        specified_pressure = (specified_pressure_bnd * Parameters.specified_pressure)
    else:
        # create scalars with value zero:
        specified_pressure_bnd = es.wherePositive(surface) * 0
        specified_pressure = es.wherePositive(surface) * 0

        for xmin, xmax, ymin, ymax, spec_pressure_segment in zip(Parameters.specified_pressure_xmin,
                                                                 Parameters.specified_pressure_xmax,
                                                                 Parameters.specified_pressure_ymin,
                                                                 Parameters.specified_pressure_ymax,
                                                                 Parameters.specified_pressure):

            print('assigning spec. pressure %0.3f to segment x=%0.3f-%0.3f, y=%0.3f-%0.3f'
                  % (spec_pressure_segment, xmin, xmax, ymin, ymax))

            spec_pressure_bound_segment = es.whereNonPositive(xmin - xy[0]) * es.whereNonNegative(xmax - xy[0]) * \
                                          es.whereNonPositive(ymin - xy[1]) * es.whereNonNegative(ymax - xy[1])

            specified_pressure_bnd += spec_pressure_bound_segment
            specified_pressure += spec_pressure_bound_segment * spec_pressure_segment

            # calculate hydrostatic pressure and add to pressure bnd:
            d = (ymax - xy[1])
            density_segment = es.sup(spec_pressure_bound_segment * specified_concentration_rho)
            dPh = spec_pressure_bound_segment * es.wherePositive(d) * d * density_segment * Parameters.g
            specified_pressure += spec_pressure_bound_segment * dPh

            print('added hydrostatic pressure to this segment: ', dPh)
            print('using density ', density_segment)
            print('warning: this assumes constant density in the fluid column at this location')
            print('hydrostatic pressure bnd condition with varying concentration/density is not implemented yet')

    #
    sp_xy, spec_pressure_array = convert_to_array(specified_pressure_bnd)
    print('number of specified pressure bnd nodes: %i' % spec_pressure_array.sum())

    su_xy, su_array = convert_to_array(surface)
    spc_xy, spec_c_array = convert_to_array(specified_concentration_bnd)
    print('number of surface nodes = %i' % su_array.sum())
    print('number of specified concentration bnd nodes: %i' % \
          spec_c_array.sum())
    print('specified concentration ', specified_concentration)
    ind = spec_c_array == 1
    max_x_sc = np.max(spc_xy[:, 0][ind])
    min_x_sc = np.min(spc_xy[:, 0][ind])
    print('active from x = %0.1f to %0.1f ' % (min_x_sc, max_x_sc))

    return (specified_pressure, specified_pressure_bnd,
            specified_concentration,
            specified_concentration_rho_f, specified_concentration_bnd,
            rch_bnd_loc, drain_bnd_loc)


def set_boundary_conditions_coastal_models(mesh, surface, sea_surface, z_surface, Parameters):

    """
    Set boundary conditions for coastal aquifer model.
    
    Parameters
    ----------
    mesh :
        escript mesh
    rho_f_0 : float
        fluid density at solute concentration = 0
    
    gamma : float
        coefficient of density increase
    g : float
        gravitational constant
    initial_condition_parameters : list
    
    Returns
    -------
    specified_pressure : escript variable
        specified pressure
    specified_pressure_bnd : escript variable
        location of specified pressure boundary
        (1=spec pressure, 0 is no spec pressure)
    specified_concentration : escript variable
        specified concentration 
    specified_concentration_bnd : escript variable
        lcoation of specified concentration boundary
    rch_bnd_loc : escript variable
        location of recharge boundary
    drain_bnd_loc : escript variable
        location of drain boundary
        
    """

    xy = mesh.getX()

    #rho_f_salt = Parameters.seawater_concentration * Parameters.gamma + \
    #    Parameters.rho_f_0
    rho_f_salt = gwflow_lib.calculate_fluid_density(
        Parameters.seawater_concentration,
        Parameters.gamma,
        Parameters.rho_f_0)
    #rho_f_fresh = Parameters.freshwater_concentration * Parameters.gamma +
    # Parameters.rho_f_0

    # drain bnd location
    drain_bnd_loc = es.whereNegative(Parameters.drain_bnd_xmin - xy[0]) * \
        es.wherePositive(Parameters.drain_bnd_xmax - xy[0]) * \
        es.wherePositive(surface)

    # recharge bnd location
    try:
        # use tags to set bnd
        rch_bnd_loc = es.Scalar(0, es.FunctionOnBoundary(mesh))
        rch_bnd_loc.setTaggedValue("land_surface1", 1)
        rch_bnd_loc.setTaggedValue("land_surface2", 1)

        xy_rch = rch_bnd_loc.getFunctionSpace().getX()

        # make sure no recharge at sea or coastal node:
        rch_bnd_loc = rch_bnd_loc * es.wherePositive(xy_rch[0])

    except RuntimeError:
        print('could not set recharge bnd using tags, using location instead')
        rch_bnd_loc = (es.whereNegative(Parameters.recharge_mass_flux_xmin -
                                        xy[0])
                       * es.wherePositive(Parameters.recharge_mass_flux_xmax -
                                          xy[0])
                       * es.wherePositive(surface))

    rch_xy, rch_array = convert_to_array(rch_bnd_loc)
    print('number of active recharge nodes = %i' % rch_array.sum())

    #
    if Parameters.specified_concentration_surface is True:
        print('assigning specified concentration only to the surface up to x=0:')
        if sea_surface is not None:
            specified_concentration_bnd = \
                es.whereNegative(xy[0]) * sea_surface + es.whereNonNegative(xy[0]) * surface
        else:
            #specified_concentration_bnd = es.whereNegative(xy[0]) * surface
            specified_concentration_bnd = es.whereNegative(Parameters.specified_concentration_xmin - xy[0]) * \
                es.wherePositive(Parameters.specified_concentration_xmax - xy[0]) * \
                es.wherePositive(surface)
            
        specified_concentration = (es.whereNegative(xy[0]) * Parameters.seawater_concentration
                                   + es.whereNonNegative(xy[0]) * Parameters.freshwater_concentration)

    else:

        specified_concentration_bnd = surface * 0
        specified_concentration = surface * 0
        for xmin, xmax, ymin, ymax, spec_conc_segment in zip(Parameters.specified_concentration_xmin,
                                                             Parameters.specified_concentration_xmax,
                                                             Parameters.specified_concentration_ymin,
                                                             Parameters.specified_concentration_ymax,
                                                             Parameters.specified_concentration):

            print('assigning specified concentration of %0.4f at surface from x=%0.2f to %0.2f, y=%0.2f to %0.2f' %
                  (spec_conc_segment, xmin, xmax, ymin, ymax))

            specified_concentration_bnd_segment = (es.whereNonPositive(xmin - xy[0]) *
                                                   es.whereNonNegative(xmax - xy[0]) *
                                                   es.whereNonPositive(ymin - xy[1]) *
                                                   es.whereNonNegative(ymax - xy[1]))
            specified_concentration_bnd += specified_concentration_bnd_segment
            specified_concentration += (specified_concentration_bnd_segment * spec_conc_segment)

    # calculate product of concentration and density for bnd nodes
    specified_concentration_rho = \
        gwflow_lib.calculate_fluid_density(specified_concentration,
                                           Parameters.gamma,
                                           Parameters.rho_f_0)
    specified_concentration_rho_f = (specified_concentration *
                                     specified_concentration_rho)
    # set specified pressure boundary:
    print('setting specified pressure nodes at surface')
    if Parameters.specified_pressure_surface is True:
        if sea_surface is not None:
            specified_pressure_bnd = \
                es.whereNegative(Parameters.specified_pressure_xmin - xy[0]) * \
                es.wherePositive(Parameters.specified_pressure_xmax - xy[0]) * \
                es.wherePositive(surface) +\
                es.whereNegative(Parameters.specified_pressure_xmin - xy[0]) * \
                es.wherePositive(Parameters.specified_pressure_xmax - xy[0]) * \
                es.wherePositive(sea_surface)
        else:
            specified_pressure_bnd = \
                es.whereNegative(Parameters.specified_pressure_xmin - xy[0]) * es.wherePositive(Parameters.specified_pressure_xmax - xy[0]) * es.wherePositive(surface)

        specified_pressure = (specified_pressure_bnd *
                      Parameters.specified_pressure)
    else:
        # create scalars with value zero:
        specified_pressure_bnd = es.wherePositive(surface) * 0
        specified_pressure = es.wherePositive(surface) * 0

        for xmin, xmax, ymin, ymax, spec_pressure_segment in zip(Parameters.specified_pressure_xmin,
                                                                 Parameters.specified_pressure_xmax,
                                                                 Parameters.specified_pressure_ymin,
                                                                 Parameters.specified_pressure_ymax,
                                                                 Parameters.specified_pressure):

            print('assigning spec. pressure %0.3f to segment x=%0.3f-%0.3f, y=%0.3f-%0.3f'
                  % (spec_pressure_segment, xmin, xmax, ymin, ymax))

            spec_pressure_bound_segment = es.whereNonPositive(xmin - xy[0]) * es.whereNonNegative(xmax - xy[0]) * \
                                          es.whereNonPositive(ymin - xy[1]) * es.whereNonNegative(ymax - xy[1])

            specified_pressure_bnd += spec_pressure_bound_segment
            specified_pressure += spec_pressure_bound_segment * spec_pressure_segment

            # calculate hydrostatic pressure and add to pressure bnd:
            d = (ymax - xy[1])
            density_segment = es.sup(spec_pressure_bound_segment * specified_concentration_rho)
            dPh = spec_pressure_bound_segment * es.wherePositive(d) * d * density_segment * Parameters.g
            specified_pressure += spec_pressure_bound_segment * dPh

            print('added hydrostatic pressure to this segment: ', dPh)
            print('using density ', density_segment)
            print('warning: this assumes constant density in the fluid column at this location')
            print('hydrostatic pressure bnd condition with varying concentration/density is not implemented yet')

    #
    sp_xy, spec_pressure_array = convert_to_array(specified_pressure_bnd)
    print('number of specified pressure bnd nodes: %i' % spec_pressure_array.sum())

    if Parameters.add_seawater_pressure is True:

        print('adding pressure of overlying seawater to specified pressure boundary')

        depth_under_sea = Parameters.sea_water_level - xy[1]
        sea_extent = es.whereNegative(z_surface - Parameters.sea_water_level)

        # calculate extra pressure by seawater
        specified_pressure_seawater = sea_extent * (
            depth_under_sea * Parameters.g * rho_f_salt)
        specified_pressure = specified_pressure + specified_pressure_seawater

        print('added pressure = ', specified_pressure_seawater)

        print('specified pressure sea bottom = ', specified_pressure * surface)

        if sea_surface is not None:
            print('specified pressure sealevel = ', specified_pressure * sea_surface)

    su_xy, su_array = convert_to_array(surface)
    spc_xy, spec_c_array = convert_to_array(specified_concentration_bnd)
    print('number of surface nodes = %i' % su_array.sum())
    print('number of specified concentration bnd nodes: %i' % \
          spec_c_array.sum())
    print('specified concentration ', specified_concentration)
    ind = spec_c_array == 1
    max_x_sc = np.max(spc_xy[:, 0][ind])
    min_x_sc = np.min(spc_xy[:, 0][ind])
    print('active from x = %0.1f to %0.1f ' % (min_x_sc, max_x_sc))

    return (specified_pressure, specified_pressure_bnd,
            specified_concentration,
            specified_concentration_rho_f, specified_concentration_bnd,
            rch_bnd_loc, drain_bnd_loc)


def depth_sw_interface_Glover1959(x, k, viscosity, hydr_gradient, thickness, rho_f, rho_s, gamma,
                                  g=9.81, Qmax=None,
                                  Q_correction_factor=0.5):

    """
    Calculate depth of the fresh-salt water interface in a coastal aquifer, following Glover (1959) JGR.
    """

    if hydr_gradient == 0.0:
        hydr_gradient = 1e-4

    K = k * rho_f * g / viscosity

    Q = K * thickness * hydr_gradient

    if Q_correction_factor is not None:
        print('correcting freshwater Q at shoreline by factor %0.2e to compensate for partitioning offshore ' \
              'and onshore discharge' % Q_correction_factor)
        Q = Q * Q_correction_factor

    if Qmax is not None and Qmax > 0:
        if Q > Qmax:
            print('calculated Darcy flux at shoreline exceeds maximum Q (ie recharge input)')
            print('calculated = %0.2e m2/s' % Q)
            print('max = %0.2e m2/s' % Qmax)
            print('using max Q in calculation of fresh-salt water interface following Glover (1959)')

            Q = Qmax

    gamma = (rho_s - rho_f) / rho_f

    if gamma == 0:
        gamma = (1025. - 998.7) / 998.7

    y2 = 2 * Q / (gamma * K) * x + Q**2 / (gamma**2 * K**2)

    # get rid of < 0 values:
    y2_fix = es.wherePositive(y2) * y2

    depth = y2_fix**0.5

    # transform to depth below sea lvl
    y = - depth

    # calculate intersection with sea bottom
    C = Q / (gamma * K)
    b = 2 * C
    c = C**2
    a = -hydr_gradient**2

    #intersect1 = (-b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)
    disc = b**2 - 4 * a * c
    #intersect_top = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
    intersect_top = (-b + disc**0.5) / (2 * a)

    # and intersection with bottom aquifer
    a = -hydr_gradient**2
    D = thickness
    b = - ( 2* D * hydr_gradient + 2 * C)
    c = D**2 - C**2
    disc = b ** 2 - 4 * a * c

    #intersect_bottom = (-b - np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    intersect_bottom = (-b - disc**0.5) / (2 * a)
    #intersect_bottom2 = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)

    return y, intersect_top, intersect_bottom


def set_initial_conditions(mesh, z_surface, Parameters):
    """
    Set initial conditions, default model setup
    
    Parameters
    ----------
    mesh :
        escript mesh
    Parameters : class
        class containing all model parameters
        
    Returns
    -------
    P0 : escript variable
        initial pressure (Pa)
    conc0 : escript variable
        initial solute concentration (kg/kg)
    rho_f_init : escript variable
        initial fluid density (kg/m^3)
    
    """

    # find top node
    xy = mesh.getX()

    if Parameters.analytical_solution_initial_h is True:
        print('calculating initial hydraulic head using analytical solution')
        # calculate initial hydraulic head using analytical solution
        R = Parameters.recharge_flux
        B = Parameters.thickness
        K = Parameters.k * Parameters.rho_f_0 * 9.81 / Parameters.viscosity
        L = Parameters.L

        # calculate hydraulic head using analytical solution for confined aq
        # with uniform recharge + Dupuit assumptions
        h = R / (K * B) * (L * xy[0] - 0.5 * xy[0]**2)

        h_adj = (es.wherePositive(xy[0]) * h
                 + es.wherePositive(-xy[0]) * z_surface)

        h_adj2 = (es.wherePositive(h_adj - z_surface) * z_surface +
                  es.whereNonPositive(h_adj - z_surface) * h_adj)

        z_surface = h_adj2

        print('analytical solution for initial h: ', h_adj)

    depth = z_surface - xy[1]
    
    # TODO: make this a bit less ugly:
    conc0 = es.wherePositive(xy[0]) * Parameters.initial_concentration
    conc0 += es.whereNonPositive(xy[0]) * Parameters.initial_concentration
    #conc0 = es.Scalar(0, es.Function(mesh))
    
    # initial fluid density
    rho_f_init = gwflow_lib.calculate_fluid_density(conc0,
                                                    Parameters.gamma,
                                                    Parameters.rho_f_0)

    # set initial hydrostatic pressure
    pressure0 = rho_f_init * Parameters.g * depth

    return pressure0, conc0, rho_f_init


def set_initial_conditions_coastal_models(mesh, z_surface, seawater, Parameters):
    """
    Set initial conditions for coastal aquifer models.
    
    Parameters
    ----------
    mesh :
        escript mesh
    Parameters : class
        class containing all model parameters
        
    Returns
    -------
    P0 : escript variable
        initial pressure (Pa)
    conc0 : escript variable
        initial solute concentration (kg/kg)
    rho_f_init : escript variable
        initial fluid density (kg/m^3)
    
    """

    # find top node
    xy = mesh.getX()

    if Parameters.analytical_solution_initial_h is True:
        print('calculating initial hydraulic head using analytical solution')
        # calculate initial hydraulic head using analytical solution
        R = Parameters.recharge_flux
        B = Parameters.thickness
        K = Parameters.k * Parameters.rho_f_0 * 9.81 / Parameters.viscosity
        L = Parameters.L

        # calculate hydraulic head using analytical solution for confined aq
        # with uniform recharge + Dupuit assumptions
        h = R / (K * B) * (L * xy[0] - 0.5 * xy[0]**2)

        h_adj = (es.wherePositive(xy[0]) * h
                 + es.wherePositive(-xy[0]) * z_surface)

        h_adj2 = (es.wherePositive(h_adj - z_surface) * z_surface +
                  es.whereNonPositive(h_adj - z_surface) * h_adj)

        z_surface = h_adj2

        print('analytical solution for initial h: ', h_adj)

    depth = z_surface - xy[1]
    sea_extent = es.whereNegative(z_surface - Parameters.sea_water_level)
    depth_under_sea = Parameters.sea_water_level - z_surface

    # initial salinity and fluid density

    # assume Ghyben-Herzberg initial conditions with P=0 at surface elevation
    if Parameters.ghyben_herzberg is True:

        #print 'initial concentration follow Ghyben-Herzberg relation'
        print('new initial cond: fresh-salt water interface now calculated using analytical solution Glover (1959) JGR')
        print('still called Ghyben-Herzberg in input file though....')

        old_style = False
        if old_style is True:
            conc0 = es.wherePositive(-xy[1] - z_surface * 40) * \
                Parameters.seawater_concentration
            conc0 += es.whereNonPositive(-xy[1] - z_surface * 40) * \
                Parameters.freshwater_concentration
        else:
            rho_f = Parameters.rho_f_0 * Parameters.freshwater_concentration * Parameters.gamma + Parameters.rho_f_0
            rho_s = Parameters.rho_f_0 * Parameters.seawater_concentration * Parameters.gamma + Parameters.rho_f_0
            Qmax = Parameters.recharge_flux * Parameters.L

            # calculate depth of sw interface following Glover (1959) JGR
            y_sw, int_sw_top, int_sw_bottom = depth_sw_interface_Glover1959(xy[0],
                                                                            Parameters.k,
                                                                            Parameters.viscosity,
                                                                            Parameters.topo_gradient,
                                                                            Parameters.thickness,
                                                                            rho_f, rho_s, Parameters.gamma,
                                                                            Qmax=Qmax)

            print('elevation fresh-saltwater interface: ', y_sw)

            conc0 = es.whereNegative(xy[1] - y_sw) * Parameters.seawater_concentration
            conc0 += es.whereNonNegative(xy[1] - y_sw) * Parameters.freshwater_concentration

            print('initial concentration: ', conc0)

    else:
        # TODO: make this a bit less ugly:
        conc0 = es.wherePositive(xy[0]) * Parameters.freshwater_concentration
        conc0 += es.whereNonPositive(xy[0]) * Parameters.seawater_concentration
        #conc0 = es.Scalar(0, es.Function(mesh))
        
    # make sure seawater salinity at sea:
    if seawater is not None:
        conc0 = es.wherePositive(seawater) * Parameters.seawater_concentration + es.whereZero(seawater) * conc0

    # initial fluid density
    rho_f_init = gwflow_lib.calculate_fluid_density(conc0,
                                                    Parameters.gamma,
                                                    Parameters.rho_f_0)

    rho_f_salt = gwflow_lib.calculate_fluid_density(Parameters.seawater_concentration,
                                                    Parameters.gamma,
                                                    Parameters.rho_f_0)

    # set initial hydrostatic pressure
    pressure0 = rho_f_init * Parameters.g * depth

    # add pressure of overlying seawater under the sea
    if Parameters.add_seawater_pressure is True:
        print('adding pressure of seawater to initial conditions')
        additional_pressure = sea_extent * (depth_under_sea * Parameters.g * rho_f_salt)
        print('additional pressure = ', additional_pressure)
        pressure0_adj = pressure0 + additional_pressure
    else:
        pressure0_adj = pressure0

    return pressure0_adj, conc0, rho_f_init


def add_cols_to_df(df, newcols):
    """
    Add a series of columns to a pandas dataframe.
    """

    for newcol in newcols:
        if newcol not in df.columns:
            df[newcol] = np.nan

    return df


def run_model_scenario(scenario_name,
                       model_output_dir,
                       model_file_adj,
                       ModelOptions,
                       Parameters,
                       mesh_function,
                       mesh_fn,
                       g=9.81):
    """
    Generate and run a single transient model scenario for coupled density-driven flow.
    
    """

    year = 365.25 * 24 * 60.0 * 60.0

    # setup model
    # set up mesh
    print('creating mesh')
    #mesh, surface, z_surface = mesh_function(Parameters, mesh_fn)
    mesh, surface, sea_surface, seawater, z_surface = mesh_function(Parameters,
                                                                    mesh_fn)

    print('set up PDEs')
    pressure_pde = es.linearPDEs.LinearPDE(mesh)
    solute_pde = es.linearPDEs.LinearPDE(mesh)

    pressure_pde.setSymmetryOn()
    solute_pde.setSymmetryOn()

    #pressure_pde.setReducedOrderOn()

    # set direct solver for pressure transport
    if ModelOptions.pressure_transport_solver is 'GMRES':
        print('using GMRES solver for pressure PDE')
        pressure_pde.getSolverOptions().setSolverMethod(es.SolverOptions.GMRES)
    elif ModelOptions.pressure_transport_solver is 'DIRECT':
        print('using direct solver for pressure PDE')
        pressure_pde.getSolverOptions().setSolverMethod(
            es.SolverOptions.DIRECT)
    elif ModelOptions.pressure_transport_solver is 'PCG':
        print('using PCG solver for pressure PDE')
        pressure_pde.getSolverOptions().setSolverMethod(
            es.SolverOptions.PCG)
    elif ModelOptions.pressure_transport_solver is 'TFQMR':
        print('using TFQMR solver for pressure PDE')
        pressure_pde.getSolverOptions().setSolverMethod(
            es.SolverOptions.TFQMR)

    # set direct solver for solute transport
    if ModelOptions.solute_transport_solver is 'GMRES':
        print('using GMRES solver for solute PDE')
        solute_pde.getSolverOptions().setSolverMethod(es.SolverOptions.GMRES)
    elif ModelOptions.solute_transport_solver is 'DIRECT':
        print('using direct solver for solute PDE')
        solute_pde.getSolverOptions().setSolverMethod(es.SolverOptions.DIRECT)
    elif ModelOptions.solute_transport_solver is 'PCG':
        print('using PCG solver for solute PDE')
        solute_pde.getSolverOptions().setSolverMethod(es.SolverOptions.PCG)
    elif ModelOptions.solute_transport_solver is 'PCG':
        print('using TFQMR solver for solute PDE')
        solute_pde.getSolverOptions().setSolverMethod(es.SolverOptions.TFQMR)

    # set gravity vector
    g_vector = es.Vector((0, -g), es.Function(mesh))

    print('set up boundary conditions')
    #(specified_pressure, specified_pressure_bnd,
    # specified_concentration, specified_concentration_rho_f,
    # specified_concentration_bnd, rch_bnd_loc, drain_bnd_loc) = \
    #    set_boundary_conditions(mesh, surface, z_surface, Parameters)
    if Parameters.mesh_type is 'coastal':
        (specified_pressure, specified_pressure_bnd,
        specified_concentration, specified_concentration_rho_f,
        specified_concentration_bnd, rch_bnd_loc, drain_bnd_loc) = \
            set_boundary_conditions_coastal_models(mesh, surface, sea_surface, z_surface,
                                    Parameters)
    else:
        (specified_pressure, specified_pressure_bnd,
        specified_concentration, specified_concentration_rho_f,
        specified_concentration_bnd, rch_bnd_loc, drain_bnd_loc) = \
            set_boundary_conditions(mesh, surface, z_surface, Parameters)

    print('set up initial conditions')
    if Parameters.mesh_type is 'coastal':
        pressure0, conc0, rho_f_init = set_initial_conditions_coastal_models(mesh,
                                                              z_surface,
                                                              seawater,
                                                              Parameters)
    else:
        pressure0, conc0, rho_f_init = set_initial_conditions(mesh,
                                                        z_surface,
                                                        Parameters)

    rho_f_init = gwflow_lib.calculate_fluid_density(conc0, Parameters.gamma,
                                                    Parameters.rho_f_0)

    # set initial conditions for PDE
    pressure = pressure0
    conc = conc0
    rho_f = rho_f_init
    phi = Parameters.porosity

    #
    dt = Parameters.dt0
    number_of_steps = 0
    runtime = 0
    output_runtime = 0
    output_step = 0

    # set bnd conditions
    print('setting boundary conditions')

    # calculate recharge mass flux
    recharge_mass_flux = (Parameters.recharge_flux *
                          Parameters.recharge_density)

    specified_flux = rch_bnd_loc * dt * recharge_mass_flux

    #if specified_pressure.getShape() == (1L,) or specified_pressure_bnd.getShape() == (1L,):
    #    msg = 'error, the specified pressure bnd parameter has the wrong shape.'
    #    msg += 'Most likely either the variable specified_pressure or specified_pressure_xmin in your input file '
    #    msg += 'is a list instead of a single number. When specified_pressure_surface = True these values should be ' \
    #           'floats instead of lists'
    #    raise ValueError(msg)

    pressure_pde.setValue(D=1, r=specified_pressure, q=specified_pressure_bnd, y=specified_flux)

    # create PDE to project element values to nodes:
    proj = es.linearPDEs.LinearPDESystem(mesh)
    proj.setSymmetryOn()

    # set up drain bnd
    active_seepage_bnd = None

    # construct permeability tensor
    k_tensor = es.Tensor(((0, 0), (0, 0)), es.Function(mesh))
    k_tensor[0, 0] = Parameters.k
    k_tensor[1, 1] = Parameters.k / Parameters.anisotropy

    # and vector, used to calculate flux. not sure why tensor doesnt work for
    # this, but it doesnt.
    k_vector = es.Vector((0, 0), es.Function(mesh))
    k_vector[0] = Parameters.k
    k_vector[1] = Parameters.k / Parameters.anisotropy

    # make sure low k in sea: no gw flow, just high diffusion of salt
    if seawater is not None:
        k_seawater = 1e-18
        k_tensor[0, 0] = es.whereZero(seawater) * k_tensor[0, 0] \
            + es.wherePositive(seawater) * k_seawater
        k_tensor[1, 1] = es.whereZero(seawater) * k_tensor[1, 1] \
            + es.wherePositive(seawater) * k_seawater
        k_vector[0] = es.whereZero(seawater) * k_vector[0] \
            + es.wherePositive(seawater) * k_seawater
        k_vector[1] = es.whereZero(seawater) * k_vector[1] \
            + es.wherePositive(seawater) * k_seawater

    # set up diffusivity
    if seawater is not None:
        diffusivity_seawater = Parameters.diffusivity_seawater
        diffusivity = es.whereZero(seawater) * Parameters.diffusivity \
            + es.wherePositive(seawater) * diffusivity_seawater
    else:
        diffusivity = Parameters.diffusivity

    dispersion_tensor = es.Tensor(((0, 0), (0, 0)), es.Function(mesh))

    # calculate transverse dispersivity
    t_disp = Parameters.l_disp * Parameters.disp_ratio

    # initial steady state run
    if ModelOptions.initial_steady_state_run is True:

        print('-' * 30)
        print('running initial steady state model')

        #pressure = gwflow_lib.solve_steady_state_pressure_eq(
        #    mesh, pressure_pde,
        #    rho_f, k_tensor, k_vector,
        #    Parameters.viscosity, g_vector,
        #    Parameters.Qf,
        #    rch_bnd_loc,
        #    recharge_mass_flux)

        # new steady-state function, with seepage bnd
        pressure, active_seepage_bnd = \
            gwflow_lib.solve_steady_state_pressure_eq_new(
                mesh, Parameters.topo_gradient,
                pressure_pde,
                rho_f, k_tensor, k_vector,
                Parameters.viscosity, g_vector,
                Parameters.Qf,
                rch_bnd_loc,
                recharge_mass_flux,
                specified_pressure_bnd,
                specified_pressure,
                drain_bnd_loc,
                rho_f_init,
                proj)


    # screen output for initial conditions
    q = pressure_pde.getFlux()
    ext_surf = es.integrate(surface,
                            where=es.FunctionOnBoundary(surface.getDomain()))
    ext_seepage = es.integrate(active_seepage_bnd,
                               where=es.FunctionOnBoundary(active_seepage_bnd.getDomain()))
    xy = pressure.getFunctionSpace().getX()
    h = pressure / (rho_f * g) + xy[1]

    # save mesh figure
    #pdb.set_trace()


    print('-' * 30)
    print('initial conditions:')
    print('pressure      ', pressure)
    print('h             ', h)
    print('concentration ', conc)
    print('flux          ', q * year)
    print('extent surface              ', ext_surf)
    print('extent active seepage bnd   ', ext_seepage)

    print('\n')

    ######################################
    print('-' * 30)
    print('starting transient iterations')

    reached_steady_state = False
    n_iterations = 0
    model_error = False

    pressure_differences_max = []
    pressure_differences_mean = []
    concentration_differences_max = []
    concentration_differences_mean = []
    dts = []
    runtimes = []
    max_CFL_number = 0

    go = True

    while go is True:

        #(runtime < Parameters.total_time
        #and reached_steady_state is False
        #and model_error is False):

        #sys.stdout.write('\rtimestep = %i / %i, iterations = %i    ' %
        #                 (number_of_steps, total_steps, n_iterations))
        #sys.stdout.flush()
        timestr = get_timestr(runtime)
        print('timestep = %i, t=%0.2e / %0.2e, %s, iterations=%i, max CFL=%0.2e'
              % (number_of_steps, runtime, Parameters.total_time, timestr, n_iterations, max_CFL_number))
        # store data at old timestep:
        pressure_t1 = pressure
        concentration_t1 = conc
        rho_f_t1 = rho_f

        # determine timestep length
        #dt = dt * Parameters.dt_inc
        #dts[number_of_steps]

        seepage_step = ((int(number_of_steps)
                        / Parameters.seepage_bnd_timestep_interval)
                            == (float(number_of_steps)
                                / Parameters.seepage_bnd_timestep_interval))

        # check if seepage bnd needs to be recalculated
        recalculate_seepage_bnd = True
        if (seepage_step == True or
                number_of_steps < Parameters.seepage_bnd_timestep_interval):
            recalculate_seepage_bnd = True
            #print('recalculating seepage bnd'
        if runtime > Parameters.seepage_bnd_max_time:
            recalculate_seepage_bnd = False

        # ignore convergence criterion for the first 20 timesteps
        # for some reason low-flux / low-k models have convergence problems initially
        if number_of_steps < 20:
            ignore_convergence_failure = True
        else:
            ignore_convergence_failure = False

        # iterate
        (pressure, conc, rho_f, viscosity, q, q_abs,
         active_concentration_bnd,
         active_seepage_bnd,
         active_rch_bnd,
         n_iterations,
         pressure_error, concentration_error,
         max_CFL_number,
         non_convergence,
         dt_real) = \
            gwflow_lib.iterate_coupled_flow_eqs(
                mesh, Parameters.topo_gradient,
                pressure_pde, solute_pde,
                Parameters.pressure_convergence_criterion,
                Parameters.concentration_convergence_criterion,
                Parameters.min_iterations,
                Parameters.max_iterations,
                dt, g_vector,
                pressure, conc, rho_f,
                phi,
                diffusivity,
                Parameters.l_disp,
                t_disp,
                Parameters.Qs,
                Parameters.specific_storage,
                k_tensor, k_vector,
                dispersion_tensor,
                Parameters.viscosity,
                Parameters.gamma,
                Parameters.alpha,
                Parameters.Qf,
                Parameters.rho_f_0,
                specified_pressure_bnd, specified_pressure,
                specified_concentration_bnd, specified_concentration,
                specified_concentration_rho_f,
                rch_bnd_loc,
                recharge_mass_flux,
                solute_transport=ModelOptions.solute_transport,
                steady_state=False,
                drain_loc=drain_bnd_loc,
                seepage_bnd=Parameters.seepage_bnd,
                recalculate_seepage_bnd=recalculate_seepage_bnd,
                proj=proj,
                active_seepage_bnd=active_seepage_bnd,
                concentration_bnd_inflow_only=
                Parameters.concentration_bnd_inflow_only,
                concentration_bnd_inflow_direction=
                Parameters.concentration_bnd_inflow_direction,
                coupled_iterations=ModelOptions.coupled_iterations,
                max_allowed_CFL_number=
                Parameters.max_allowed_CFL_number,
                force_CFL_timestep=Parameters.force_CFL_timestep,
                dt_max=Parameters.dt_max,
                iterate_seepage_in_one_timestep=
                Parameters.iterate_seepage_in_one_timestep,
                calculate_viscosity=Parameters.calculate_viscosity,
                verbose=ModelOptions.verbose,
                ignore_convergence_failure=ignore_convergence_failure)

        #except RuntimeError, msg:
        #    model_error = True
        #    print('solver failed, stopping this particular model scenario'
        #    print msg

        if non_convergence is True:
            model_error = True
            go = False
            print('non convergence, stopping this particular model scenario')

        # show iteration on screen
        #sys.stdout.write('%i / %i\r' % (number_of_steps + 1, total_steps))
        #sys.stdout.flush()

        #runtime += dts[number_of_steps]
        runtime += dt_real
        output_runtime += dt_real

        # calculate changes in pressure and concentration over 1 timestep
        pressure_difference = (pressure - pressure_t1) / (dt_real / year)
        concentration_difference = (conc - concentration_t1) / (dt_real / year)

        pressure_differences_max.append(es.Lsup(pressure_difference))
        pressure_differences_mean.append(es.integrate(pressure_difference))
        concentration_differences_max.append(es.Lsup(concentration_difference))
        concentration_differences_mean.append(es.integrate(concentration_difference))

        # store timestep size
        dts.append(dt_real)
        runtimes.append(runtime)

        # calculate if model has reached steady state
        if (es.Lsup(pressure_difference) <
                    Parameters.max_pressure_change_steady_state
                and es.Lsup(concentration_difference) <
                    Parameters.max_concentration_change_steady_state
                and Parameters.stop_when_steady_state is True
                and model_error is False
                and runtime >= Parameters.total_time):

            print('reached steady state at timestep %i' % number_of_steps)
            print('max abs. pressure change per year %0.2f' %
                  es.Lsup(pressure_difference))
            print('max abs. concentration change per year %0.2f' %
                  es.Lsup(concentration_difference))
            reached_steady_state = True
            go = False
        else:
            #print('no steady state yet'
            #print('max abs. pressure change per year %0.2f' % \
            #      es.Lsup(pressure_difference)
            pass
        # check if runtime exceeds max runtime
        if Parameters.stop_when_steady_state is False and runtime >= Parameters.total_time:
            print('exceeded runtime of ', Parameters.total_time)
            go = False

        if Parameters.stop_when_steady_state is True and (runtime >= Parameters.max_runtime
                                                          or number_of_steps >= Parameters.max_timesteps):
            print('exceeded maximum runtime without reaching steady state')
            print('max runtime = ', Parameters.max_runtime)
            go = False

        # give output to screen at regular interval or after last timestep
        if (output_runtime >= Parameters.output_interval
                or runtime == dt_real
                or go is False):

            print('')

            # reset output time counter
            if runtime != dt_real:
                output_runtime = 0

            timestr = get_timestr(runtime)

            print('model scenario: %s / %s' % (str(scenario_name), model_file_adj))
            print('\ntimestep %i, t = %i sec, (%s)' \
                  % (number_of_steps + 1, runtime, timestr))
            print('number of iterations: %i' % n_iterations)
            print('pressure iteration error = %0.3e (Pa)' % pressure_error)
            print('concentration iteration error = %0.3e (kg/kg)' \
                  % concentration_error)
            print('max CFL number = %0.2f' % max_CFL_number)
            print('timestep size = %0.2e sec or %0.2e days' \
                  % (dt_real, dt_real / (24 * 60 * 60.0)))

            print('\nmodel variables:\n')
            print('pressure (Pa)               ', pressure)
            print('concentration (kg/kg)       ', conc)
            print('fluid density (kg/m^3)      ', rho_f)
            print('qh (m/yr)                   ', q[0] * year)
            print('qv (m/yr)                   ', q[1] * year)
            print('pressure at surface         ', pressure * surface)
            print('concentration at surface    ', conc * surface)
            print('qv at surface (m/yr)        ', q[1] * surface * year)
            print('P change per year (Pa)    ', pressure_difference)
            print('C change per year (kg/kg) ', concentration_difference)

            # calculate hydraulic head
            xy = pressure.getFunctionSpace().getX()
            h = pressure / (rho_f * g) + xy[1]

            # interpolate flux to nodes
            proj.setValue(D=es.kronecker(mesh), Y=q)
            #proj.setValue(D=1.0, Y=q)
            nodal_flux = proj.getSolution()

            ext_surf = es.integrate(
                surface,
                where=es.FunctionOnBoundary(surface.getDomain()))
            ext_seepage = es.integrate(
                active_seepage_bnd,
                where=es.FunctionOnBoundary(active_seepage_bnd.getDomain()))
            ext_rch = es.integrate(
                active_rch_bnd,
                where=es.FunctionOnBoundary(active_rch_bnd.getDomain()))
            total_flux_over_surface = es.integrate(
                nodal_flux * surface * year,
                where=es.FunctionOnBoundary(nodal_flux.getDomain()))

            # normalize total flux
            # not sure if this works, this just seems to assume the top
            # boundary is horizontal....
            #flux_surface_norm = (nodal_flux * surface * year *
            #                     nodal_flux.getDomain().getNormal())
            #flux_surface_norm_old = (q * surface * year *
            #                         q.getDomain().getNormal())

            # TODO: figure out why the rotation doesnt work here, or at
            # least gives funny results...
            nodal_q_norm = gwflow_lib.rotate_vector_escript(nodal_flux,
                                                            Parameters.topo_gradient)
            #flux_surface_norm = nodal_flux * surface * year
            #flux_surface_norm = nodal_q_norm * surface * year
            flux_surface_norm = (nodal_q_norm * surface * year *
                                 nodal_q_norm.getDomain().getNormal())

            total_flux_over_surface_norm = es.integrate(
                flux_surface_norm,
                where=es.FunctionOnBoundary(flux_surface_norm.getDomain()))
            # flux at recharge bnd
            total_rch_flux = es.integrate(
                flux_surface_norm * active_rch_bnd,
                where=es.FunctionOnBoundary(flux_surface_norm.getDomain()))

            # flux at seepage bnd
            total_seepage_flux = es.integrate(
                flux_surface_norm * active_seepage_bnd,
                where=es.FunctionOnBoundary(flux_surface_norm.getDomain()))

            # flux at sea
            total_submarine_flux = es.integrate(
                flux_surface_norm * es.whereNegative(xy[0]),
                where=es.FunctionOnBoundary(flux_surface_norm.getDomain()))

            land_flux = flux_surface_norm * es.whereNonNegative(xy[0])
            land_flux_in = flux_surface_norm * es.whereNonNegative(xy[0]) \
                             * es.whereNegative(flux_surface_norm[1])
            land_flux_out = flux_surface_norm * es.whereNonNegative(xy[0]) \
                             * es.wherePositive(flux_surface_norm[1])

            total_land_flux_in = es.integrate(
                land_flux_in,
                where=es.FunctionOnBoundary(land_flux_in.getDomain()))
            total_land_flux_out = es.integrate(
                land_flux_out,
                where=es.FunctionOnBoundary(land_flux_out.getDomain()))

            submarine_flux = flux_surface_norm * es.whereNegative(xy[0])
            submarine_flux_in = flux_surface_norm * es.whereNegative(xy[0]) \
                                * es.whereNegative(flux_surface_norm[1])
            submarine_flux_out = flux_surface_norm * es.whereNegative(xy[0]) \
                                 * es.wherePositive(flux_surface_norm[1])

            total_submarine_flux_in = es.integrate(
                submarine_flux_in,
                where=es.FunctionOnBoundary(submarine_flux_in.getDomain()))
            total_submarine_flux_out = es.integrate(
                submarine_flux_out,
                where=es.FunctionOnBoundary(submarine_flux_out.getDomain()))

            # find transition point from recharge to seepage

            # min flux to be considered inflow or outflow zone. default limit
            # = 1 cm/yr
            flux_buffer = 1.0e-1 / year

            inflow_bool = es.whereNegative(flux_surface_norm[1])
            outflow_bool = es.wherePositive(flux_surface_norm[1])

            inflow_bool_buffer = es.whereNegative(flux_surface_norm[1] + flux_buffer)
            outflow_bool_buffer = es.wherePositive(flux_surface_norm[1] - flux_buffer)

            inflow_land = \
                inflow_bool * es.wherePositive(inflow_bool.getDomain().getX())
            outflow_land = \
                outflow_bool * es.wherePositive(outflow_bool.getDomain().getX())
            inflow_sea = \
                inflow_bool * es.whereNonPositive(inflow_bool.getDomain().getX())
            outflow_sea = \
                outflow_bool * es.whereNonPositive(outflow_bool.getDomain().getX())
            ext_outflow = \
                es.integrate(outflow_bool,
                    where=es.FunctionOnBoundary(outflow_bool.getDomain()))
            ext_outflow_land = \
                es.integrate(
                    outflow_land,
                    where=es.FunctionOnBoundary(outflow_bool.getDomain()))
            ext_outflow_sea = es.integrate(
                outflow_sea,
                where=es.FunctionOnBoundary(outflow_bool.getDomain()))
            ext_inflow = es.integrate(inflow_bool,
                      where=es.FunctionOnBoundary(inflow_bool.getDomain()))
            ext_inflow_land = es.integrate(
                inflow_land,
                where=es.FunctionOnBoundary(inflow_bool.getDomain()))
            ext_inflow_sea = es.integrate(
                inflow_sea,
                where=es.FunctionOnBoundary(inflow_bool.getDomain()))

            outflow_land_threshold = \
                outflow_bool_buffer * es.wherePositive(outflow_bool_buffer.getDomain().getX())
            outflow_sea_threshold = \
                outflow_bool_buffer * es.whereNonPositive(outflow_bool_buffer.getDomain().getX())
            ext_outflow_land_threshold = \
                es.integrate(
                    outflow_land_threshold,
                    where=es.FunctionOnBoundary(outflow_land_threshold.getDomain()))
            ext_outflow_sea_threshold = es.integrate(
                outflow_sea_threshold,
                where=es.FunctionOnBoundary(outflow_sea_threshold.getDomain()))

            min_land_flux = es.inf(land_flux[1])
            max_land_flux = es.sup(land_flux[1])

            max_seepage_flux = es.sup(active_seepage_bnd * flux_surface_norm[1])
            min_seepage_flux = es.inf(active_seepage_bnd * flux_surface_norm[1])

            min_submarine_flux = es.inf(submarine_flux[1])
            max_submarine_flux = es.sup(submarine_flux[1])

            print('extent surface               ', ext_surf)
            print('extent active seepage bnd    ', ext_seepage)
            print('extent active recharge bnd   ', ext_rch)
            print('total flux over surface      ', total_flux_over_surface)
            print('total flux normal to surface ', total_flux_over_surface_norm)
            print('total active recharche flux  ', total_rch_flux)
            print('total active seepage flux    ', total_seepage_flux)
            print('total land flux in           ', total_land_flux_in)
            print('total land flux out          ', total_land_flux_out)
            print('total submarine bnd flux     ', total_submarine_flux)
            print('total submarine bnd flux in  ', total_submarine_flux_in)
            print('total submarine bnd flux out ', total_submarine_flux_out)
            print('extent inflow                ', ext_inflow)
            print('extent outflow               ', ext_outflow)
            print('extent inflow land           ', ext_inflow_land)
            print('extent outflow land          ', ext_outflow_land)
            print('extent inflow sea            ', ext_inflow_sea)
            print('extent outflow sea           ', ext_outflow_sea)
            print('extent outflow land > threshold ', ext_outflow_land_threshold)
            print('extent outflow sea  > threshold ', ext_outflow_sea_threshold)
            print('min land bnd flux             ', min_land_flux)
            print('max land bnd flux             ', max_land_flux)
            print('min seepage bnd flux          ', min_seepage_flux)
            print('max seepage bnd flux          ', max_seepage_flux)
            print('min submarine bnd flux        ', min_submarine_flux)
            print('max submarine bnd flux        ', max_submarine_flux)

            spc_xy, spec_c_array = convert_to_array(specified_concentration_bnd)
            print('number of specified concentration bnd nodes: %i' % \
                spec_c_array.sum())
            print('specified concentration ', specified_concentration)

            #####################
            # boundary conditions
            #####################
            boundary_conditions = [specified_pressure_bnd,
                                   specified_pressure,
                                   specified_concentration_bnd,
                                   active_concentration_bnd,
                                   specified_concentration,
                                   specified_concentration_rho_f,
                                   rch_bnd_loc,
                                   active_seepage_bnd]

            boundary_fluxes = [flux_surface_norm,
                               land_flux_in,
                               land_flux_out,
                               submarine_flux,
                               submarine_flux_in,
                               submarine_flux_out]

            boundary_flux_stats = [total_flux_over_surface_norm,
                                   total_rch_flux,
                                   total_seepage_flux,
                                   total_land_flux_in,
                                   total_land_flux_out,
                                   total_submarine_flux,
                                   total_submarine_flux_in,
                                   total_submarine_flux_out,
                                   ext_rch,
                                   ext_seepage,
                                   ext_inflow,
                                   ext_outflow,
                                   ext_inflow_land,
                                   ext_outflow_land,
                                   ext_inflow_sea,
                                   ext_outflow_sea,
                                   ext_outflow_land_threshold,
                                   ext_outflow_sea_threshold,
                                   min_land_flux, max_land_flux,
                                   min_seepage_flux, max_seepage_flux,
                                   min_submarine_flux, max_submarine_flux]

            if ModelOptions.save_vtk_files_all_steps is True:
                # save model output to vtk file
                vtk_folder = os.path.join(model_output_dir, 'vtk_files')

                if not os.path.exists(vtk_folder):
                    os.makedirs(vtk_folder)

                fn_vtk = os.path.join(vtk_folder,
                                      '%s_ts%i.vtu' %
                                      (model_file_adj, output_step))

                # nodal flux with nodata value
                nodata = -99999
                flux_surface_plot = \
                    nodal_flux * surface + nodata * es.whereZero(surface)

                print('saving %s' % fn_vtk)

                if sea_surface is None:
                    sea_surface_save = surface
                else:
                    sea_surface_save = sea_surface

                esys.weipa.saveVTK(fn_vtk,
                                   pressure=pressure,
                                   concentration=conc,
                                   h=h,
                                   flux=q,
                                   qx=q[0],
                                   qy=q[1],
                                   kx=k_vector[0],
                                   ky=k_vector[1],
                                   nodal_flux=nodal_flux,
                                   nodal_flux_surface=flux_surface_plot,
                                   surface=surface,
                                   sea_surface=sea_surface_save,
                                   specified_pressure_bnd=specified_pressure_bnd,
                                   active_seepage_bnd=active_seepage_bnd,
                                   recharge_bnd=rch_bnd_loc,
                                   active_concentration_bnd=active_concentration_bnd,
                                   flux_surface_norm=flux_surface_norm)

            print('\ncontinuing iterations\n')
            # add increment to output step counter
            output_step += 1

        # adjust timestep
        # check if timestep was changed for CFL criterion
        if number_of_steps > 25:
            if dt == dt_real:
                dt = dt * Parameters.dt_inc
            elif Parameters.force_CFL_timestep is False:
                dt = dt_real * Parameters.dt_inc

        if dt > Parameters.dt_max:
            dt = Parameters.dt_max

        number_of_steps += 1

    return (mesh, surface, sea_surface,
            k_vector,
            pressure,
            conc, rho_f, viscosity, h, q, q_abs, nodal_flux,
            pressure_difference, concentration_difference,
            pressure_differences_max, concentration_differences_max,
            pressure_differences_mean, concentration_differences_mean,
            dts, runtimes, number_of_steps, output_step,
            boundary_conditions,
            boundary_fluxes,
            boundary_flux_stats,
            reached_steady_state)
