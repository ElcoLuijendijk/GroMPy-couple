from __future__ import print_function

"""
functions to solve groundwater flow, solute transport and heat
transport equations using the generic finite element model escript/Finley 

Elco Luijendijk
2012- 2017
McGill University & Goettingen University

<elco.luijendijk at geo.uni-goettingen.de>
"""

import math
import numpy as np
try:
    import esys.escript as es
except ImportError:
    print('warning, could not find escript module')


def calculate_q(k_vector, viscosity, pressure, density, g_vector):

    """
    Calculate flux vector (q).
    """

    q = -(k_vector / viscosity *
          (es.grad(pressure) - density * g_vector))

    return q


def solve_steady_state_pressure_eq_new(mesh, topo_gradient, pressure_pde,
                                       rho_f, k_tensor, k_vector,
                                       viscosity, g_vector,
                                       fluid_source,
                                       rch_bnd_loc,
                                       recharge_mass_flux,
                                       specified_pressure_bnd,
                                       specified_pressure,
                                       drain_bnd_loc,
                                       fluid_density,
                                       proj, debug=True):

    """
    Solve the steady-state fluid flow equation for fluid pressure using escript with optional seepage bnd condition.
    """

    year = 365.25 * 24 * 60 * 60.0

    # calculate boundary flux
    specified_flux = rch_bnd_loc * recharge_mass_flux

    a = rho_f * k_tensor / viscosity * es.kronecker(mesh)
    d = 0.0
    x = rho_f**2 * k_vector / viscosity * g_vector
    y = fluid_source

    #
    pressure_pde.setValue(A=a, D=d, X=x, Y=y, y=specified_flux,
                          q=specified_pressure_bnd, r=specified_pressure)

    # calculate pressure, without seepage bnd
    pressure = pressure_pde.getSolution()

    debug = True
    if debug is True:
        print('initial calculated steady-state pressure: ', pressure)

    if es.sup(drain_bnd_loc) == 0:
        print('no drain or seepage bnd')
        active_seepage_bnd = es.wherePositive(drain_bnd_loc)

        return pressure, active_seepage_bnd

    # check if h > surface elevation
    #pressure_threshold = 0.01 * 9.81 * 1000.0
    pressure_threshold = 0.0
    active_seepage_bnd = es.wherePositive(drain_bnd_loc * pressure - pressure_threshold)

    if debug is True:
        print('active seepage nodes: ', np.sum(np.array(active_seepage_bnd.toListOfTuples())))
        print('potential seepage nodes: ', np.sum(np.array(drain_bnd_loc.toListOfTuples())))

    n_seepage_change = 9999
    n_seepage_nodes = 99999
    n_iter = 0
    max_seepage_iter = 2000
    while n_seepage_change > 0 and n_iter < max_seepage_iter:

        # add seepage to specified pressure bnd and update bnd conditions
        specified_pressure_bnd_mod = \
        es.wherePositive(
            specified_pressure_bnd + active_seepage_bnd)

        #active_rch_bnd = rch_bnd_loc * es.whereZero(specified_pressure_bnd)
        active_rch_bnd = rch_bnd_loc
        specified_flux = active_rch_bnd * recharge_mass_flux
        pressure_pde.setValue(r=specified_pressure,
                              q=specified_pressure_bnd_mod,
                              y=specified_flux)

        # recalculate pressure
        pressure = pressure_pde.getSolution()

        if debug is True:
            print('new pressure: ', pressure)

        # calculate flux
        q = calculate_q(k_vector, viscosity, pressure,
                        fluid_density, g_vector)
        proj.setValue(D=es.kronecker(mesh), Y=q)
        try:
            nodal_q = proj.getSolution()
        except RuntimeError(msg):
            print('error, non-convergence')
            print(msg)
            non_convergence = True

        # calculate max vertical flux into the model domain at
        # drain bnd nodes

        # using outer norm to calculate correct bnd flux does not work because
        # flux cannot be interpolated to same functionspace as pressure variable
        #nodal_q_norm_bnd = nodal_q * nodal_q.getDomain().getNormal()
        #nodal_q_norm = es.interpolate(nodal_q_norm_bnd, active_seepage_bnd.getFunctionSpace())

        nodal_q_norm = rotate_vector_escript(nodal_q, topo_gradient)

        flux_seepage_bnd = active_seepage_bnd * nodal_q_norm[1]

        # remove inlfow nodes that are <= recharge
        flux_seepage_bnd_corr = flux_seepage_bnd + recharge_mass_flux / fluid_density
        # changed back to old method to speed up seepage bnd calc
        #flux_seepage_bnd_corr = flux_seepage_bnd

        # and remove only the worst x % of seepage nodes to avoid oscillations:
        seepage_threshold = es.inf(flux_seepage_bnd_corr) * 0.5
        #seepage_threshold = 0.0

        seepage_inflow_nodes = \
            es.whereNegative(flux_seepage_bnd_corr - seepage_threshold)

        removed_seepage_inflow_nodes = seepage_inflow_nodes

        if debug is True:
            print('number of seepage inflow nodes: ', \
                np.sum(np.array(seepage_inflow_nodes.toListOfTuples())))

            xmin_seep = es.inf(seepage_inflow_nodes * seepage_inflow_nodes.getDomain().getX()[0] + (1-seepage_inflow_nodes) * 999999)
            xmax_seep = es.sup(seepage_inflow_nodes * seepage_inflow_nodes.getDomain().getX()[0])

            print('from x= %0.2f m to x= %0.2f m' % (xmin_seep, xmax_seep))


        # add boundary nodes with P>0 to seepage bnd
        new_seepage_nodes = \
            es.wherePositive(drain_bnd_loc
                             * (1 - active_seepage_bnd)
                             * pressure)

        if debug is True:
            print('number of new seepage nodes: ', \
                np.sum(np.array(new_seepage_nodes.toListOfTuples())))

        # update the seepage bnd
        active_seepage_bnd = (active_seepage_bnd
                              + new_seepage_nodes
                              - removed_seepage_inflow_nodes)

        n_seepage_nodes_old = n_seepage_nodes
        n_seepage_nodes = np.sum(np.array(active_seepage_bnd.toListOfTuples()))

        n_seepage_change = np.abs(n_seepage_nodes_old - n_seepage_nodes)

        if debug is True:
            print('final active seepage nodes: ', np.sum(np.array(active_seepage_bnd.toListOfTuples())))
            print('potential seepage nodes: ', np.sum(np.array(drain_bnd_loc.toListOfTuples())))

        if n_seepage_nodes_old < n_seepage_nodes:
            print('lowest number of seepage nodes reached, stopping iterations')
            n_seepage_change = 0

        print('seepage iteration %i' % n_iter)
        print('seepage threshold ', seepage_threshold * year)
        print('change in seepage nodes from %0.0f to %0.0f' % (n_seepage_nodes_old,
                                                               n_seepage_nodes))

        n_iter += 1

    # update specified pressure bnd condition
    specified_pressure_bnd_mod = \
    es.wherePositive(
        specified_pressure_bnd + active_seepage_bnd)

    #active_rch_bnd = rch_bnd_loc * es.whereZero(specified_pressure_bnd)
    active_rch_bnd = rch_bnd_loc
    specified_flux = active_rch_bnd * recharge_mass_flux
    pressure_pde.setValue(r=specified_pressure,
                          q=specified_pressure_bnd_mod,
                          y=specified_flux)

    # recalculate pressure
    pressure = pressure_pde.getSolution()

    if debug is True:
        print('final pressure: ', pressure)

    return pressure, active_seepage_bnd


def update_pressure_pde(pressure_pde,
                        pressure_tmin1, phi, specific_storage,
                        k_tensor, k_vector, rho_f,
                        viscosity, dt,
                        rch_bnd_loc,
                        recharge_mass_flux,
                        fluid_source, g_vector,
                        gamma, concentration_change_rate, alpha,
                        temperature_change_rate):
    """
    Solve the transient groundwater flow equation with variable density.
    
    Parameters
    ----------
    pressure_pde
        escript partial differential equation
    mesh
        esys.escript mesh
    steady_state : bool
        solve transient or steady state eq.
    implicit_solution : bool
        implicit or explicit solution
    pressure_tmin1
        pressure at last timestep
    phi
        porosity
    specific_storage
        specific storage
    k
        permeability
    rho_f
        fluid density
    viscosity
        fluid viscosity
    dt
        timestep duration
    fluid_source
        fluid source term
    g_vector
        gravity vector
    gamma
        solute expansion coefficient (dimensionless?)
    concentration_change_rate
        change of concentration over time
    alpha
        thermal expansion coefficient (1/T)
    temperature_change_rate
        change of temperature over time
        
    Returns
    -------
    pressure : 
        Pressure
    """

    # calculate boundary flux
    specified_flux = rch_bnd_loc * dt * recharge_mass_flux

    # set coefficients
    a_coeff = dt * rho_f * k_tensor / viscosity
    d_coeff = specific_storage
    x_coeff = dt * rho_f**2 * k_vector / viscosity * g_vector
    y_coeff = \
        specific_storage * pressure_tmin1 - \
        dt * gamma * phi * concentration_change_rate - \
        dt * alpha * phi * temperature_change_rate + \
        dt * fluid_source

    pressure_pde.setValue(A=a_coeff, D=d_coeff, X=x_coeff, Y=y_coeff,
                          y=specified_flux)

    return pressure_pde


def update_solute_transport_pde(mesh, solute_pde,
                                concentration_old, v,
                                dt, solute_source,
                                dispersion_tensor,
                                diffusivity, l_disp, t_disp,
                                rho_f,
                                verbose=False):
    """
    Solve the solute transport equation.
    """

    # calculate hydrodynamic dispersivity tensor from molecular diffusivity
    # and longitudinal and transverse dispersivity

    # calculate absolute velocity
    v_abs = (v[0] ** 2 + v[1] ** 2) ** 0.5

    # get rid of 0 values of velocity
    v_abs_1 = es.whereZero(v_abs) * 1e-20 + es.whereNonZero(v_abs) * v_abs

    Dxx = l_disp * (v[0] ** 2) / v_abs_1 + diffusivity
    Dyy = t_disp * (v[1] ** 2) / v_abs_1 + diffusivity

    # and set dispersion tensor to 0 where 0 absolute velocity
    Dxx = (Dxx * es.whereNonZero(v_abs)
           + es.whereZero(v_abs) * diffusivity)
    Dyy = (Dyy * es.whereNonZero(v_abs)
           + es.whereZero(v_abs) * diffusivity)

    # horizontal and vertical values of dispersion (Dxx and Dyy)
    dispersion_tensor[0, 0] = Dxx
    dispersion_tensor[1, 1] = Dyy

    # off-diagonal terms in tensor (Dxy and Dyx)
    # set to 0 for now, model becomes numerically unstable
    # testing with fixed values did not result in significant cahnges in solute conc...
    # todo, figure out why this is unstable or check the potential error
    dispersion_tensor[0, 1] = (l_disp - t_disp) * (v[0] * v[1])/v_abs_1
    dispersion_tensor[1, 0] = (l_disp - t_disp) * (v[0] * v[1])/v_abs_1

    u_old = concentration_old

    a_coeff = dt * dispersion_tensor * es.kronecker(mesh)
    c_coeff = dt * v
    d_coeff = 1
    y_coeff = u_old + solute_source * dt

    if verbose is True:
        print('solute transport coefficients')
        print('A: ', a_coeff.getShape(), a_coeff)
        print('C: ', c_coeff.getShape(), c_coeff)
        print('D: ', d_coeff)
        print('Y: ', y_coeff.getShape(), y_coeff)

    solute_pde.setValue(A=a_coeff, C=c_coeff, D=d_coeff, Y=y_coeff)

    return solute_pde


def calculate_concentration(solute_mass, rho_f_0, gamma):

    """
    Calculate solute concentration using total solute mass and a concentration-density function.
    :param solute_mass:
    :param rho_f_0:
    :param gamma:
    :return:
    """

    # calculate new concentration, assuming no change in density:
    #concentration = u / rho_f
    a = rho_f_0 * gamma
    b = rho_f_0
    c = -solute_mass
    d = b**2 - 4 * a * c

    concentration = (-b + es.sqrt(d)) / (2 * a)

    return concentration


def calculate_fluid_density(concentration, gamma, rho_f_0):
    """
    Calculate fluid density.

    omitting density dependence on pressure and temperature for now

    Parameters
    ----------
    concentration
        solute concentration
    gamma
        coefficient of increase of density with increasing concentration
    rho_f_0
        fluid density at concentration zero
        
    Returns
    -------
    rho_f
        fluid density
        
    """

    # eq. 7 in seawat manual, lineralized version of eq.
    rho_f = rho_f_0 + rho_f_0 * gamma * concentration

    return rho_f


def calculate_angle(dx, dy):
    """
    Calculate the angle between two points.
    """

    degrees = np.arctan2(dy, dx) * 180 / np.pi

    return degrees


def rotate_vector_escript(v, gradient):
    """
    Rotate a vector `v` by angle calculated from the topographic gradient.
    """

    dx = 1.0
    dy = gradient * dx

    angle_degr = calculate_angle(dx, dy)

    angle = -np.radians(angle_degr)

    vr = v.copy()

    cos_theta = math.cos(angle)
    sin_theta = math.sin(angle)

    vr[0] = v[0]*cos_theta - v[1]*sin_theta
    vr[1] = v[0]*sin_theta + v[1]*cos_theta

    return vr


def iterate_coupled_flow_eqs(mesh, topo_gradient, pressure_pde, solute_pde,
                             pressure_convergence_criterion,
                             concentration_convergence_criterion,
                             min_iterations,
                             max_iterations,
                             dt, g_vector,
                             pressure, concentration, rho_f,
                             phi, diffusivity, l_disp, t_disp,
                             solute_source,
                             specific_storage,
                             k_tensor, k_vector,
                             dispersion_tensor,
                             viscosity,
                             gamma, alpha,
                             fluid_source,
                             rho_f_0,
                             specified_pressure_bnd,
                             specified_pressure,
                             specified_concentration_bnd,
                             specified_concentration,
                             specified_concentration_rho_f,
                             rch_bnd_loc,
                             recharge_mass_flux,
                             coupled_iterations=True,
                             solute_transport=True,
                             heat_transport=False,
                             steady_state=False,
                             proj=None,
                             drain_loc=None,
                             seepage_bnd=False,
                             recalculate_seepage_bnd=True,
                             active_seepage_bnd=None,
                             concentration_bnd_inflow_only=False,
                             concentration_bnd_inflow_direction='up',
                             max_allowed_CFL_number=None,
                             force_CFL_timestep=False,
                             dt_max=None,
                             calculate_viscosity=False,
                             verbose=False,
                             iterate_seepage_in_one_timestep=False,
                             max_seepage_iterations=50,
                             ignore_convergence_failure=False):

    """
    Iterative solve groundwater flow, solute transport and heat flow equations.
    
    solves either steady state or 1 timestep in implicit or explicit mode
    
    iterative coupling scheme of solute transport, pressure & flow eqs. and
    eqs of state follows Ackerer (2004), Geophysical Research Letters 31(12)
    
    Parameters
    ---------
    mesh : 
        escript mesh object
    pressure_pde : 
        groundwater flow PDE
    solute_pde
        solute transport PDE
    pressure_convergence_criterion : float
        convergence criterion groundwater flow eq. (Pa)
    concentration_convergence_criterion : float
        convergence criterion solute transport eq. (kg/kg)
    max_iterations : int
        max number of iterations
    dt : int
        timestep
    g_vector : 
        gravity vector (0,g)
    pressure : 
        pressure (Pa)
    concentration : 
        solute concentration (kg/kg)
    rho_f : 
        fluid density (kg / m3)
    phi :
        porosity
    D :
        solute diffusivity (...)
    l_disp :
        longitudinal dispersivity (...)
    t_disp :
        transverse dispersivity (...)
    solute_source :
        solute source (units...)
    specific_storage :
        specific storativity (...)
    k :
        permeability (m2)
    anisotropy :
        permeability anisotropy = horizontal/vertical permeability
        (dimensionless)
    viscosity : 
        viscosity (...)
    gamma :
        ?
    alpha :
        ?
    fluid_source :
        fluid source term (...)
    rho_f_0
        fluid density at solute concentration C=0 (kg/m3)
    specified_pressure_bnd
        location of specified pressure boundary
    specified_pressure
        specified pressure (Pa)
    specified_concentration_bnd
        location of specified concentration boundary
    specified_concentration
        specified concentration (kg/kg)
    rch_bnd_loc :

    recharge_mass_flux : float

    coupled_iterations : bool, optional
        couple groundwater and solute transport equations iteratively
        by adjusting density term
    solute_transport : bool, optional
        if True, simulate solute transport
    heat_transport : bool, optional
        if True, simulate heat transport
    steady_state : bool, optional
        True for steady state groundwater flow, False for transient
    verbose : bool, optional
        verbose text output
    drain_loc :
        location of drain boundary nodes
    debug : bool, optional
        debugging
    dt_max : float?
        =None                 ...
    proj : 
        escript PDE for projecting element data to nodes
    seepage_optimization_automated : boolean
        
    
    
    Returns 
    -------
    pressure_t2_i2 :
        pressure at next timestep (t2) and last iteration (i2)
    concentration_t2_i2 :
        solute concentration (kg/kg)
    rho_f_t2_i2 :
        fluid density
    iteration : int
        number of iterations
    dt_max :
        max timestep size
        
    """


    # calculate transverse dispersivity
    #t_disp = l_disp * disp_ratio

    year = 365.25 *24 * 60 * 60.

    if verbose is True:
        print('running iterative solver for pressure and concentration PDEs')
        if coupled_iterations is False:
            print('pressure and concentration are not coupled')

    #pressure_new = pressure
    pressure_old_ts = pressure
    concentration_old_ts = concentration
    fluid_density_new = fluid_density_old = rho_f
    #pressure_t1 = pressure.copy()
    #concentration_t1 = concentration.copy()

    # added 22 jun 2016, not sure if this is ok:
    active_rch_bnd = rch_bnd_loc

    if coupled_iterations is True and calculate_viscosity is True:
        viscosity_new = calculate_viscosity_simple(concentration)
    else:
        viscosity_new = viscosity

    active_specified_concentration_bnd = specified_concentration_bnd

    iteration = 0
    converged = False
    non_convergence = False
    ele_size = None
    q = None
    v = None

    while converged is False and non_convergence is False:

        if verbose is True:
            print('iteration ', iteration)
            if iteration > 0:
                print('pressure convergence ', es.Lsup(pressure_conv))

        if solute_transport is True:

            # get flux
            q = calculate_q(k_vector, viscosity_new, pressure_old_ts,
                            fluid_density_new, g_vector)

            v = q / phi

            # calculate new solute concentration
            concentration_old_iteration = concentration

            # finite element solute transport
            if concentration_bnd_inflow_only is True and iteration == 0:
                # only apply concentration bnd for inflow into model domain
                # assumes a horizontal model bnd
                # TODO: calculate flux normal to model boundary to account
                # for non-horizontal upper boundaries

                proj.setValue(D=es.kronecker(mesh), Y=q)
                try:
                    nodal_q = proj.getSolution()
                except RuntimeError(msg):
                        print('error, non-convergence')
                        print(msg)
                        non_convergence = True

                nodal_q_norm = rotate_vector_escript(nodal_q, topo_gradient)

                nodal_v = nodal_q / phi

                if concentration_bnd_inflow_direction == 'up':
                    inflow_bnd = (es.whereNegative(nodal_q_norm[1]) *
                                  specified_concentration_bnd)
                elif concentration_bnd_inflow_direction == 'down':
                    inflow_bnd = (es.wherePositive(nodal_q_norm[1]) *
                                  specified_concentration_bnd)
                elif concentration_bnd_inflow_direction == 'left':
                    inflow_bnd = (es.wherePositive(nodal_q[0]) *
                                  specified_concentration_bnd)
                elif concentration_bnd_inflow_direction == 'right':
                    inflow_bnd = (es.whereNegative(nodal_q[0]) *
                                  specified_concentration_bnd)

                if es.sup(inflow_bnd) > 0:
                    active_specified_concentration_bnd = inflow_bnd
                else:
                    min_x = es.inf(
                        specified_concentration_bnd *
                        specified_concentration_bnd.getDomain().getX()[0])

                    active_specified_concentration_bnd = \
                        (specified_concentration_bnd *
                         es.whereZero(
                             specified_concentration_bnd.getDomain().getX()[0]
                             - min_x))

                    if verbose is True:
                        print('warning, outflow for all specified ' \
                              'concentration boundary nodes')
                        #print 'using entire bnd instead'
                        #active_specified_concentration_bnd = \
                        #    specified_concentration_bnd
                        print('using first node as fixed conc bnd instead')
                        print('number of active conc bnd nodes:')
                        print(np.sum(np.array(
                            active_specified_concentration_bnd.
                            toListOfTuples())))

                if verbose is True:
                    import grompy_lib

                    xyi, ia = grompy_lib.convert_to_array(
                        active_specified_concentration_bnd)
                    xyc, ca = grompy_lib.convert_to_array(
                        specified_concentration_bnd)
                    print('inflow conc bnd nodes = %0.0f / %0.0f' \
                          % (ia.sum(), ca.sum()))
                    print('x = %0.3f - %0.3f' % (xyi[ia == 1, 0].min(),
                                                 xyi[ia == 1, 0].max()))
                    print('qv conc bnd: ', (nodal_q[1] *
                                            specified_concentration_bnd))

            #solute_pde.setValue(D=1,
            #                    r=specified_concentration_rho_f,
            #                    q=active_specified_concentration_bnd)
            solute_pde.setValue(D=1,
                                r=specified_concentration,
                                q=active_specified_concentration_bnd)

            solute_pde = update_solute_transport_pde(
                mesh, solute_pde,
                concentration_old_ts, v, dt, solute_source,
                dispersion_tensor,
                diffusivity, l_disp, t_disp, fluid_density_old)

            try:
                #solute_mass = solute_pde.getSolution()
                concentration = solute_pde.getSolution()
            except RuntimeError(error_msg):
                print('!! runtime error ', error_msg)
                print('solver options: ')
                print(solute_pde.getSolverOptions().getSummary())

                non_convergence = True

                #raise RuntimeError(error_msg)

            # calculate concentration, using new solute mass and eq of state
            #concentration_new = calculate_concentration(
            #    solute_mass, rho_f_0, gamma)

            #concentration_new = solve_solute_transport_v2(
            #        solute_pde, mesh,
            #        steady_state,
            #        concentration_t1, v, dt, solute_source,
            #        diffusivity, l_disp, t_disp, fluid_density_old,
            #        rho_f_0, gamma)

            concentration_change_rate = \
                (concentration - concentration_old_ts) / dt

        else:
            # no solute transport:
            concentration_change_rate = 0

        if heat_transport is True:
            # no temperature in models yet:
            temperature_change_rate = 0
        else:
            # no heat transport:
            temperature_change_rate = 0

        if coupled_iterations is True:
            if verbose is True:
                print('recalculating fluid density and viscosity')
            # recalculate fluid density            
            fluid_density_old = fluid_density_new
            fluid_density_new = \
                calculate_fluid_density(concentration, gamma, rho_f_0)

            if calculate_viscosity is True:
                viscosity_new = \
                    calculate_viscosity_simple(concentration)

        else:
            # leave fluid density unchanged
            concentration_change_rate = 0
            temperature_change_rate = 0

        # store old pressure
        pressure_old_iteration = pressure

        if drain_loc is None or es.sup(drain_loc) == 0:

            # calculate pressure, no drain or seepage bnd
            pressure_pde = \
                update_pressure_pde(pressure_pde,
                                    pressure_old_ts,
                                    phi, specific_storage,
                                    k_tensor, k_vector,
                                    fluid_density_new,
                                    viscosity_new, dt,
                                    rch_bnd_loc,
                                    recharge_mass_flux,
                                    fluid_source, g_vector,
                                    gamma, concentration_change_rate,
                                    alpha, temperature_change_rate)
            try:
                pressure = pressure_pde.getSolution()
            except RuntimeError(msg):
                print('error, non-convergence')
                print(msg)
                non_convergence = True
            #print 'no seepage bnd'
        else:
            # implement drain or seepage boundary
            if seepage_bnd is True:

                ## use seepage boundary:
                if active_seepage_bnd is None:
                    # calculate pressure without any drain boundary
                    pressure_pde.setValue(r=specified_pressure,
                                          q=specified_pressure_bnd)
                    active_rch_bnd = rch_bnd_loc
                else:
                    # incorporate active drain bnd of previous timestep
                    specified_pressure_bnd_mod = \
                        es.wherePositive(
                            specified_pressure_bnd + active_seepage_bnd)
                    pressure_pde.setValue(r=specified_pressure,
                                          q=specified_pressure_bnd_mod)

                    # do not change active rch bnd
                    active_rch_bnd = rch_bnd_loc

                    #active_rch_bnd = rch_bnd_loc * \
                    #                 es.whereZero(specified_pressure_bnd)
                    #specified_flux = rch_bnd_loc * dt * recharge_mass_flux

                # calculate pressure with existing seepage bnd
                pressure_pde = \
                    update_pressure_pde(pressure_pde,
                                        pressure_old_ts,
                                        phi, specific_storage,
                                        k_tensor, k_vector,
                                        fluid_density_new,
                                        viscosity_new, dt,
                                        active_rch_bnd, recharge_mass_flux,
                                        fluid_source, g_vector,
                                        gamma, concentration_change_rate,
                                        alpha, temperature_change_rate)
                try:
                    pressure = pressure_pde.getSolution()
                except RuntimeError:
                    print("error, pressure PDE solver failed")
                    converged = True
                    non_convergence = True
                    #if pressure_new not in locals():
                    #    pressure_new = pressure_t1

                # assign drain bnd nodes
                if active_seepage_bnd is None:
                    active_seepage_bnd = \
                        es.wherePositive(drain_loc * pressure)

                if iteration < max_seepage_iterations and recalculate_seepage_bnd is True:
                    # adjust seepage boundary, but only for first x iterations
                    # to avoid instability

                    if verbose is True:
                        seepage_xy = active_seepage_bnd.getDomain().getX()
                        seepage_nodes_xy = \
                            np.array(seepage_xy.toListOfTuples())
                        seepage_array = np.array(
                            active_seepage_bnd.toListOfTuples())
                        ind = seepage_array > 0
                        print('\tbefore adjustment:')
                        print('\tactive seepage bnd from x=%0.0f to %0.0f m' \
                              % (seepage_nodes_xy[ind, 0].min(),
                                 seepage_nodes_xy[ind, 0].max()))

                    # remove seepage nodes that have become source of water
                    q = calculate_q(k_vector, viscosity_new, pressure,
                                    fluid_density_new, g_vector)
                    proj.setValue(D=es.kronecker(mesh), Y=q)
                    try:
                        nodal_q = proj.getSolution()
                    except RuntimeError(msg):
                        print('error, non-convergence')
                        print(msg)
                        non_convergence = True

                    # calculate max vertical flux into the model domain at
                    # drain bnd nodes
                    # -> not possible, cannot mix face elements and normal elements
                    # later on to adjust seepage...
                    #nodal_q_norm = nodal_q * nodal_q.getDomain().getNormal()

                    #
                    nodal_q_norm = rotate_vector_escript(nodal_q, topo_gradient)

                    #flux_seepage_bnd = active_seepage_bnd * nodal_q[1]
                    flux_seepage_bnd = active_seepage_bnd * nodal_q_norm[1]

                    #flux_seepage_bnd_corr = flux_seepage_bnd +

                    seepage_change_buffer = 1e-3 / year

                    seepage_inflow_nodes = \
                        es.whereNegative(flux_seepage_bnd
                                         + recharge_mass_flux
                                         / fluid_density_new)

                    if verbose is True:

                        print('\tflux seepage bnd (m/yr): ', flux_seepage_bnd * year)
                        print('recharge ')
                        print('\tseepage inflow nodes: ', seepage_inflow_nodes)

                    #seepage_inflow_nodes = \
                    #    es.whereNegative(flux_seepage_bnd)

                    removed_seepage_inflow_nodes = seepage_inflow_nodes

                    # add boundary nodes with P>0 to seepage bnd
                    new_seepage_nodes = \
                        es.wherePositive(drain_loc
                                         * (1 - active_seepage_bnd)
                                         * pressure)

                    # update the seepage bnd
                    active_seepage_bnd = (active_seepage_bnd
                                          + new_seepage_nodes
                                          - removed_seepage_inflow_nodes)

                    if verbose is True:
                        seepage_xy = active_seepage_bnd.getDomain().getX()
                        seepage_nodes_xy = \
                            np.array(seepage_xy.toListOfTuples())
                        seepage_array = np.array(
                            active_seepage_bnd.toListOfTuples())
                        ind = np.array(seepage_array > 0)
                        print('\tafter adjustment:')
                        print('\tN=%i active seepage bnd x=%0.0f to %0.0f m' \
                              % (np.sum(ind.astype(int)),
                                 seepage_nodes_xy[ind, 0].min(),
                                 seepage_nodes_xy[ind, 0].max()))

                    if iterate_seepage_in_one_timestep is True:
                        # update the specified pressure boundary to include
                        # new seepage nodes
                        specified_pressure_bnd_mod = \
                            es.wherePositive(
                                specified_pressure_bnd + active_seepage_bnd)
                        #active_rch_bnd = rch_bnd_loc * es.whereZero(specified_pressure_bnd_mod)
                        # changed to have steady recharge bnd regardless of seepage bnd,
                        # 11 apr 2016, Elco
                        active_rch_bnd = rch_bnd_loc

                        # experiment, adjust recharge to have 0 rehcarge at seepage nodes
                        # not sure if this makes sense...
                        #specified_flux_adj = active_rch_bnd * dt * recharge_mass_flux
                        #
                        #pressure_pde.setValue(r=specified_pressure,
                        #                      q=specified_pressure_bnd_mod,
                        #                      y=specified_flux_adj)
                        pressure_pde.setValue(r=specified_pressure,
                                              q=specified_pressure_bnd_mod)

                        # recalculate pressure
                        #pressure_pde = \
                        #    update_pressure_pde(pressure_pde,
                        #                        pressure_t1,
                        #                        phi, specific_storage,
                        #                        k_tensor, k_vector,
                        #                        fluid_density_new,
                        #                        viscosity_new, dt,
                        #                        rch_bnd_loc, recharge_mass_flux,
                        #                        fluid_source, g_vector,
                        #                        gamma, concentration_change_rate,
                        #                        alpha, temperature_change_rate)

                        # recalculate pressure
                        try:
                            pressure = pressure_pde.getSolution()
                        except RuntimeError(msg):
                            print('error, non-convergence')
                            print(msg)
                            non_convergence = True

                        # remove inflow nodes again
                        #q = (k_vector / viscosity_new *
                        #     -(es.grad(pressure_new)
                        #       - fluid_density_new * g_vector)
                        #     / phi)
                        q = calculate_q(k_vector, viscosity_new, pressure,
                                        fluid_density_new, g_vector)
                        proj.setValue(D=es.kronecker(mesh), Y=q)
                        nodal_q = proj.getSolution()

                        # calculate max vertical flux into the model domain at
                        # drain bnd nodes
                        #nodal_q_norm = nodal_q * nodal_q.getDomain().getNormal()
                        nodal_q_norm = rotate_vector_escript(nodal_q, topo_gradient)
                        flux_seepage_bnd = active_seepage_bnd * nodal_q_norm[1]

                        #removed_seepage_inflow_nodes = \
                        #    es.whereNegative(flux_seepage_bnd -
                        #                     seepage_inflow_threshold)
                        #removed_seepage_inflow_nodes = \
                        #    es.whereNegative(flux_seepage_bnd
                        #                     + recharge_mass_flux
                        #                     / fluid_density_new)

                        removed_seepage_inflow_nodes = \
                            es.whereNegative(flux_seepage_bnd)

                        active_seepage_bnd = (active_seepage_bnd
                                              - removed_seepage_inflow_nodes)

                        if verbose is True:
                            seepage_xy = active_seepage_bnd.getDomain().getX()
                            seepage_nodes_xy = \
                                np.array(seepage_xy.toListOfTuples())
                            seepage_array = np.array(
                                active_seepage_bnd.toListOfTuples())
                            ind = seepage_array > 0
                            print('\tafter 2nd adjustment (removing inflow nodes):')
                            print('\tN=%i active seepage bnd from ' \
                                  'x=%0.0f to %0.0f m' \
                                  % (np.sum(ind.astype(int)),
                                     seepage_nodes_xy[ind, 0].min(),
                                     seepage_nodes_xy[ind, 0].max()))

                if iterate_seepage_in_one_timestep is True:
                    # assign updated seepage nodes:
                    specified_pressure_bnd_mod = \
                        es.wherePositive(
                            specified_pressure_bnd + active_seepage_bnd)

                    #active_rch_bnd = rch_bnd_loc * es.whereZero(specified_pressure_bnd_mod)
                    active_rch_bnd = rch_bnd_loc

                    # experiment, adjust recharge to have 0 rehcarge at seepage nodes
                    # not sure if this makes sense...
                    # ! probably the source of long timescale instability of seepage/rch bnd

                    #specified_flux_adj = active_rch_bnd * dt * recharge_mass_flux

                    #pressure_pde.setValue(r=specified_pressure,
                    #                      q=specified_pressure_bnd_mod,
                    #                      y=specified_flux_adj)
                    pressure_pde.setValue(r=specified_pressure,
                                          q=specified_pressure_bnd_mod)

                    # recalculate final pressure
                    #pressure_pde = \
                    #    update_pressure_pde(pressure_pde,
                    #                        pressure_t1,
                    #                        phi, specific_storage,
                    #                        k_tensor, k_vector,
                    #                        fluid_density_new,
                    #                        viscosity_new, dt,
                    #                        rch_bnd_loc, recharge_mass_flux,
                    #                        fluid_source, g_vector,
                    #                        gamma, concentration_change_rate,
                    #                        alpha, temperature_change_rate)

                    try:
                        pressure = pressure_pde.getSolution()
                    except RuntimeError(msg):
                        print('error, non-convergence')
                        print(msg)
                        non_convergence = True

        # calculate convergence criteria
        pressure_conv = pressure - pressure_old_iteration

        if solute_transport is True:
            conc_conv = concentration - concentration_old_iteration
        else:
            conc_conv = 0.0

        # check whether iterations have converged or not:
        if (es.Lsup(pressure_conv) < pressure_convergence_criterion) and \
                (es.Lsup(conc_conv) < concentration_convergence_criterion)\
                and iteration + 1 >= min_iterations:
            if iteration > 0 and verbose is True:
                print('iterations converged after %i iterations' % iteration)
            converged = True
        else:
            if verbose is True:
                print('iteration %i, max. pressure change %0.3e ' \
                      % (iteration, es.Lsup(pressure_conv)))
                print('              max. C change %0.3e ' \
                      % (es.Lsup(conc_conv)))

        if iteration + 1 >= max_iterations:
            print('warning, reached maximum number of %i iterations' \
                  % (iteration + 1))
            print('iteration %i, max. pressure change %0.3e Pa, ' \
                  'convergence at %0.2e Pa' \
                  % (iteration, es.Lsup(pressure_conv),
                     pressure_convergence_criterion))
            print('              max. C change %0.3e kg/kg, ' \
                  'convergence at %0.2e kg/kg' \
                  % (es.Lsup(conc_conv), concentration_convergence_criterion))
            converged = True
            non_convergence = True

        # check CFL number
        #max_CFL_number = calculate_CFL_number(q, dt)
        if ele_size is None:
            ele_size = q.getDomain().getSize()
        #print ele_size - q.getDomain().getSize()

        CFL_number = q * dt / ele_size
        max_CFL_number = es.Lsup(CFL_number)

        if max_CFL_number > 0.5 and verbose is True:
            print('warning, max CFL number = %0.2f, exceeds 0.5' \
                  % max_CFL_number)

        # recaclulcate timestep if max timestep exceeds CFL number
        if (max_allowed_CFL_number is not None
                and max_CFL_number > max_allowed_CFL_number
                and iteration == 0) \
                or (force_CFL_timestep is True and iteration <= 1):

            # make sure iteration is repeated
            converged = False

            #CFL_number / flux * flux.getDomain().getSize() = dtc /

            dtc = max_allowed_CFL_number / q * ele_size
            new_timestep = es.inf((dtc**2)**0.5)
            if dt_max is not None and new_timestep > dt_max:
                new_timestep = dt_max

            dt = new_timestep

            if verbose is True:
                print('max CFL number = ', max_CFL_number)
                print('changing timestep from %0.2e sec to %0.2e sec' \
                      % (dt, new_timestep))

        if coupled_iterations is False:
            converged = True

        iteration += 1

    return (pressure, concentration, fluid_density_new,
            viscosity_new, q, v,
            active_specified_concentration_bnd,
            active_seepage_bnd,
            active_rch_bnd,
            iteration,
            es.Lsup(pressure_conv),
            es.Lsup(conc_conv),
            max_CFL_number,
            non_convergence,
            dt)


def equations_of_state_batzle1992(P, T, C):

    """
    Calculate density using equation of state provided by Batzle and Wang (1992) Geophysics 57(11).

    Parameters
    ----------
    P : float or array
        pressure (Pa)
    T : float or array
        temperature (degrees C)
    C : float or array
        solute concentration (kg / kg)

    Returns
    -------
    rho_b : float or array
        water density (kg m^-3)

    """

    rho_w = 1 + 1e-6 * (-80*T - 3.3 * T**2 + 0.00175 * T**3 + 489 * P
                        - 2 * T * P + 0.016 * T**2 * P - 1.3e-5 * T**3 * P
                        - 0.333 * P**2 - 0.002 * T * P**2)

    rho_b = rho_w + C * (0.668 + 0.44 * C + 1e-6 * (300 * P - 2400 * P * C
                                                    + T * (80 + 3 * T
                                                           - 3300 * C
                                                           - 13 * P
                                                           + 47 * P * C)))

    return rho_b


def calculate_viscosity(C, T=20.0):

    """
    Calculate fluid viscosity following Batzle & Wang (1992).

    """

    viscosity = 1e-3 * (0.1 + 0.333 * C + (1.65 + 91.9*C**3)
                        * es.exp(-(0.42 * (C**0.8 - 0.17)**2
                                 + 0.045) * T**0.8))

    return viscosity


def calculate_viscosity_simple(C, T=20.0):

    """
    Calculate viscosity using a simplified function with fixed temperature at T=20 degr.

    linear fit to Batzle (1992) eq. over range of C=0.0 to 0.035 (ie, fresh to seawater range)
    """

    viscosity = 9.8080e-04 + 2.6515e-03 * C

    return viscosity


def calculate_CFL_number(flux, dt):

    """
    Calculate the CFL (Courant, Friedrich, Lewy) stability criterion.
    """

    CFL_number = flux * dt / flux.getDomain().getSize()

    max_CFL_number = es.Lsup(CFL_number)

    return max_CFL_number