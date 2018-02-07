"""
Model parameters for Grompy.
"""

import numpy as np


class ModelOptions():

    """
    Class containing various model options.
    """

    # flag to run single or multiple model scenarios
    run_multiple_scenarios = True

    # initial base run when running multiple model parameter sets:
    initial_base_run = False

    # run steady state or transient models
    steady_state = False

    # calculate initial pressure field using a steady-state solution
    initial_steady_state_run = True

    # turn solute transport on or off
    solute_transport = True

    # coupled fluid flow & solute transport
    coupled_iterations = True

    # save vtk file for final model results
    save_vtk_files = True

    # save vtk file for intermediate timesteps too
    save_vtk_files_all_steps = True

    # save pressure, concentration and flux variables to .csv files
    save_variables_to_csv = False

    # calculate fluxes over the model boundaries
    report_boundary_fluxes = True

    # calculate sgd statistics
    analyze_sgd = False

    # pressure and solute transport solvers
    # choices: 'DEFAULT', 'DIRECT', 'GMRES'
    # see escript manual for more info
    # DEFAULT = iterative solver, fast, but less stable
    # GMRES = stable but slow
    # DIRECT = very slow for large problems, but stable
    pressure_transport_solver = 'GMRES'
    solute_transport_solver = 'GMRES'

    # choose whether to run all combinations of input parameters or
    # vary parameters one by one to create a
    # sensitivity analysis
    # choose either 'combinations' or 'sensitivity'
    # or choose 'file' to read the parameter combinations from a .csv file
    model_scenario_list = 'sensitivity'

    # name of file containing parameter combinations for model scenarios
    # only used if model_scenario_list = 'file'
    model_scenario_file = ''

    # name for subfolder to store this set of model runs
    model_output_dir = 'model_output/salt_wedge_benchmark'

    scenario_name = 'sw'

    # option to overwrite model results or start numbering at last model found in folder
    # and append new model results to old ones
    overwrite_old_results = True

    # number of model runs at the same time
    # warning, each model run will use the number of processors set by the environment variable OMP_NUM_THREADS
    # which should be located in your .bashrc file
    max_proc = 4

    # option to run the set of model experiments in reverse order
    invert_model_run_order = False


class ModelParameters(dict):
    """
    Class that contains all default model variables.
    """

    g = 9.80665

    day = 24 * 60 * 60.0
    year = 365.25 * day

    ########
    # time:
    ########
    # initial timestep size
    dt0 = 0.25  # (sec)

    # increase in timestep length after each timestep
    dt_inc = 1.0

    # maximum timestep size
    dt_max = 10.0  # (sec)

    # total duration
    total_time = 80.0 * 60.0  # (sec)

    # time interval for which model stats are reported and figures are created
    output_interval = 5.0 * 60.0

    # maximum courant number allowed. Timestep sizes are reduced if
    # courant number is exceeded
    # set to None to turn this feature off
    max_allowed_CFL_number = 1.0
    #max_allowed_CFL_number = None

    # force timesteps to be recalculated to have CFL condition = max_allowed_CFL_number
    # turning this option on makes sure that the timestep size is always the max stable size
    force_CFL_timestep = False

    # stop model run automatically when steady-state is reached
    stop_when_steady_state = False

    # max. runtime and timesteps for model runs
    max_runtime = 5000.0 * year
    max_timesteps = 10000

    # threshold values for change in pressure and concentration for
    # assuming steady-state conditions and stopping model runs
    max_pressure_change_steady_state = 0.25  # (Pa / yr)
    max_concentration_change_steady_state = 1e-4  # (kg/kg / yr)

    ##################
    # mesh parameters
    ##################

    # mesh type, choices:
    # 'rectangle' : rectangular mesh
    # 'coastal' : coastal aquifer with 3 blocks, middle block at
    #  fresh-salt water interface with grid refinement
    #  interface calculated following Glover (1959)
    mesh_type = 'rectangle'

    # thickness of model domain
    thickness = 0.26  # (m)

    # grid cell size, not used when mesh type = 'rectangle'
    cellsize = 0.0025  # (m)

    # x and y cell size, only used when mesh type = 'rectangle'
    cellsize_x = 0.0025  # (m)
    cellsize_y = 0.0025  # (m)

    # length of model domain
    # = length of landward side if coastal mesh is used
    L = 0.53  # (m)

    #################################################################
    ## coastal mesh params, only used if mesh_type = 'coastal'
    #################################################################
    # topographic gradient
    topo_gradient = 0.00  # (m/m)
    # length of seaward side of model domain, only used for 'coastal' mesh
    L_sea = 0.0  # (m)
    # distance on either side of coastline with smaller cell size
    buffer_distance_sea = 500.0  # (m)
    buffer_distance_land = 250.0  # (m)
    # finer grid in middle block, only used for 3block mesh type:
    # factor of 0.1 means 10 times smaller grid cells
    grid_refinement_factor = 0.3

    # finer grid at sea to make sure recirc. sgd is ok
    grid_refinement_factor_sea = 1.0
    grid_refinement_factor_seawater = 1.0

    #####################
    # boundary conditions
    #####################

    # groundwater mass flux (density x flux) over top boundary
    recharge_flux = 0.0 / year  # (kg / (m^2 s)) (-> density * flux)
    # density of recharge fluid
    recharge_density = 998.7  # (kg / m^3)

    # specify min and max x-coordinates of zone that receives recharge
    # ignored if None is specified
    # use arbitrarily low and large xmin, xmax number to include entire model 
    # domain
    # warning: xmin and xmax not used for now due to problems with escript
    # boundary flux automatically assigned to entire land surface boundary.
    recharge_mass_flux_xmin = 1e6  # (m)
    recharge_mass_flux_xmax = 1e6  # (m)

    # spec pressure boundary value
    # additional P at upper right model bnd
    # SS-1: h= 0.267 -> dh = 0.007 m -> P = 68.569 Pa
    # SS-2: h= 0.262 -> dh = 0.002 m -> P = 19.591 Pa
    # SS-3: h= 0.2655 -> dh = 0.0055 m -> P = 53.876 Pa
    specified_pressure = [0, 68.569]  # (Pa)

    # location of specified pressure boundary
    # list of min locations
    specified_pressure_xmin = [0.0, 0.53]
    specified_pressure_xmax = [0.0, 0.53]

    # specified pressure bnd only applied to surface
    specified_pressure_surface = False

    # if specified_pressure_surface = False, specify y coordinates of
    # specified pressure boundary
    specified_pressure_ymin = [0.0, 0.0]
    specified_pressure_ymax = [0.255, 0.260]

    # salinity to use for calculating hydrostatic increase in pressure
    # at pressure boundary nodes
    specified_pressure_salinity = [0.03624, 0.0]

    # add pressure of overlying seawater in coastal aquifers:
    add_seawater_pressure = False

    # spec concentration boundary
    specified_concentration = [0.03624, 0.0]  # (kg / kg)
    # spec concentration boundary location
    specified_concentration_xmin = [0.0, 0.53]  # (m))
    specified_concentration_xmax = [0.0, 0.53]  # (m)

    # specified pressure bnd only applied to surface
    specified_concentration_surface = False

    # if specified_concentration_surface = False
    # set min and max y coordinate of specified concentration boundary
    specified_concentration_ymin = [0.0, 0.0]  # (m))
    specified_concentration_ymax = [0.255, 0.26]  # (m)

    # switch off concentration boundary for discharge nodes
    concentration_bnd_inflow_only = True

    # direction of inflow, choose 'left' 'right', 'up' or 'down'
    # 'left' means that any inflow is coming from the left
    # ie. concentration bnd will be switched off for nodes where qx > 0
    # TODO: modify the code so that grompy recognizes what the inflow direction is
    concentration_bnd_inflow_direction = 'left'

    # seepage boundary location
    drain_bnd_xmin = 1e6  # (m)
    drain_bnd_xmax = 1e6  # (m)

    # type of drain boundary
    # if True: convert drain nodes to specified head if hydraulic head
    # is above the surface
    # if False: regular flux-based drain boundary function
    seepage_bnd = False

    # option to run seepage bnd only once every x timesteps:
    seepage_bnd_timestep_interval = 1

    # full seepage iteration in 1 timestep or just remove inflow nodes for next timestep:
    iterate_seepage_in_one_timestep = False

    # option to stop recalculating seepage boundary after x time, to avoid
    # long computational time due to small oscillations
    seepage_bnd_max_time = 10000 * year

    # fluid source term
    Qf = 0.0  # (kg/(m^3 s)) (= fluid mass source)
    Qf_xmin = None
    Qf_xmax = None

    # solute source term
    Qs = 0.0
    Qs_xmin = None
    Qs_xmax = None

    ############################
    # initial conditions params
    ############################
    # add Ghyben-Herzberg initial salinity distribution
    # ignored if set to False
    ghyben_herzberg = False
    sea_water_level = 0  # (m)

    # seawater and freshwater concentrations
    seawater_concentration = 0.0  # (kg/kg)
    freshwater_concentration = 0.0  # (kg/kg)

    # derive initial pressure from analytical solution for 2D flow with
    # uniform recharge
    analytical_solution_initial_h = False

    # solve for pressure only using an initial steady-state model run
    # and use this as initial condition for subsequent runs
    initial_steady_state_run = True

    #################
    # other params
    ##################
    # permeability, estimate from hydraulic conductivity, using fresh water properties
    k = 1050 / (24.0 * 60.0 * 60.0) * 8.94e-4 / (998.7 * g)  # (m^2)
    #k = 1e-12

    # anisotropy (= horizontal permeability / vertical permeability)
    anisotropy = 1.0  # (dimensionless)

    #
    porosity = 0.385  # (V/V)

    # specific storativity
    specific_storage = 1.0e-4  # (...)

    # longitudinal dispersivity
    l_disp = 0.001  # (m)
    # transverse dispersivity ratio
    # used to calculate transverse dispersivity (transverse_disp = disp_ratio * l_disp)
    disp_ratio = 0.1  # (dimensionless)

    ## settings for iterative solver solute & pressure PDEs
    pressure_convergence_criterion = 1.0e-7  # [Pa]
    concentration_convergence_criterion = 1.0e-7  # [kg/kg]

    # maximum iterations for sequential iterative solving of solute transport and
    # fluid flow equations
    min_iterations = 4
    max_iterations = 200

    # fluid viscosity
    # note, grompy can also automatically calculate viscosity from concentration and temperature data following
    # Batzle and Wang (1992). This option is controlled by calculate_viscosity in the iterate_coupled_flow_eqs function
    # in the module gwflow_lib
    viscosity = 8.94e-4
    reference_porosity = 0.25  # (non dimensional)
    reference_density = 1000.0  # [kg m-3]

    ## eq. of state constants
    # thermal expansion coefficient
    alpha = 207e-6  # (1/K)
    #  density of water at solute concentration = 0 (and T=20 C)
    rho_f_0 = 998.872  # (kg/m^3)
    # solute expansion coefficient
    #gamma = 0.68412  # (dimensionless) -> linear fit to Batzle & Wang (1992)
    # calculate solute expansion coefficient for salt water intrusion
    # experiment by Goswami & Clement (2007)
    gamma = (1026.0 - rho_f_0) / (0.03624 * rho_f_0)  # (dimensionless)

    # solute diffusion coefficient
    diffusivity = 1.0e-9  # (m^2/s^1)

    # separate diffusivity for seawater,
    # used to simulate stable concentration of seawater by high
    # diffusivity and fixed concentration a sealevel
    #diffusivity_seawater = 1.0e-4


class ParameterRanges():

    """
    Specify the range of parameter values to include in model runs.

    variables can be either lists, tuples or numpy arrays

    Set value to None to ignore the parameter range and use the default value
    that is specified in the ModelParameters class

    add _s to any parameter name in the ``model_parameters.py`` file that
    you would like to change. For example topo_gradient_s = [0.01, 0.1] will
    change parameter topo_gradient in the ``model_parameters.py`` file and
    will first run a model scenario with a gradient of 0.01 and then one
    with a gradient of 0.1

    variables can be either lists, tuples or numpy arrays

    set value to None to ignore the parameter range and use the default value
    that is specified in the ModelParameters class

    """

    # specified pressure in model experiments Goswami and Clement (2007):
    # SS-1: h= 0.267 -> dh = 0.007 m -> P = 68.569 Pa
    # SS-2: h= 0.262 -> dh = 0.002 m -> P = 19.591 Pa
    # SS-3: h= 0.2655 -> dh = 0.0055 m -> P = 53.876 Pa
    specified_pressure_s = [[0, 68.569], [0, 19.591], [0, 53.876]]
