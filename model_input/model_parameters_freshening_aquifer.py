"""
Model parameters for Grompy.
"""

import numpy as np


class ModelOptions():

    """
    Class containing various model options.
    """

    verbose = False

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
    model_output_dir = 'model_output/freshening'

    scenario_name = 'ts_steady'

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
    dt0 = 5.0 * day  # (sec)
    # multiplier for increase in timestep length after each timestep
    dt_inc = 1.03
    # max timestep size
    dt_max = 100.0 * year  # (sec)
    # total duration (sec). this is the minimum runtime if the option stop_when_steady_state is True
    total_time = 9000.0 * year  # (sec)

    # number of timesteps for which model stats are reported and figures are
    # created
    output_interval = 500.0 * year

    # maximum Courant number allowed. Timestep sizes are reduced if
    # Courant number is exceeded
    # set to None to turn this feature off
    max_allowed_CFL_number = 1.0
    #max_allowed_CFL_number = None

    # force timesteps to be recalculated to have CFL condition = max_allowed_CFL_number
    # turning this option on makes sure that the timestep size is always the max stable size
    force_CFL_timestep = False

    # stop model run automatically when steady-state is reached
    stop_when_steady_state = True

    # max. runtime and timesteps for model runs
    max_runtime = 9000.0 * year
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
    # 'coastal' : coastal aquifer with 4 blocks.
    #             2 middle blocks with grid refinement at left and right of
    #             the fresh-salt water interface
    #             fresh-salt water bnd calculated using ghyben-herzberg,
    #             hydraulic head calculated using analytical solution gw flow
    #             in aquifer with constant recharge
    # "terrestrial" : simple aquifer with a stream on one side
    mesh_type = 'standard'

    # thickness of model domain
    thickness = 30.0  # (m)
    # grid cell size
    cellsize = 1.0  # (m)
    # x and y cell size, only used when mesh type = 'rectangle'
    cellsize_x = 1.0  # (m)
    cellsize_y = 1.0  # (m)

    # length of model domain
    L = 500  # (m)

    ## coastal mesh params:
    # topographic gradient proximal part of system
    topo_gradient = 0.025  # (m/m)

    # topographic gradient after x_topo_break
    topo_gradient_hinterland = 0.1

    topo_break = False

    x_topo_break = 250.0

    #####################
    # boundary functions
    #####################
    # groundwater mass flux (density x flux) over top boundary
    recharge_flux = 0.2 / year

    # density of recharge fluid
    recharge_density = 998.7  # (kg / m^3)

    # specify min and max x-coordinates of zone that receives recharge
    # ignored if None is specified
    # use arbitrarily low and large xmin, xmax number to include entire model 
    # domain
    # warning: xmin and xmax not used for now due to problems with escript
    # boundary flux automatically assigned to entire land surface boundary.
    recharge_mass_flux_xmin = 1e-6  # (m)
    recharge_mass_flux_xmax = 1e6  # (m)

    # spec pressure boundary value
    specified_pressure = 0  # (Pa)

    # location of specified pressure boundary
    # list of min locations
    specified_pressure_xmin = -0.01
    specified_pressure_xmax = 0.01

    # specified pressure bnd only applied to surface
    specified_pressure_surface = True

    # if specified_pressure_surface is False, specify y coordinates of
    # specified pressure boundary
    specified_pressure_ymins = [0.0]
    specified_pressure_ymaxs = [0.0]

    # salinity to use for calculating hydrostatic increase in pressure
    # at pressure boundary nodes
    specified_pressure_salinities = [0.0]

    # spec concentration boundary location
    specified_concentration_xmin = -0.01  # (m))
    specified_concentration_xmax = 1e6  # (m)

    # specified pressure bnd only applied to surface
    specified_concentration_surface = True

    # min and max y coordinate,
    # ignored if specified_concentration_surface = True
    specified_concentration_ymin = [-0.01]  # (m))
    specified_concentration_ymax = [1e6]  # (m)

    specified_concentration = 0.0

    # switch off concentration boundary for inflow nodes
    concentration_bnd_inflow_only = False

    # direction of inflow, choose 'left' 'right', 'up' or 'down'
    # 'left' means that any inflow is coming from the left
    # ie. concentration bnd will be switched off for nodes where
    # qx > 0
    concentration_bnd_inflow_direction = 'up'

    # seepage boundary location
    drain_bnd_xmin = -0.001  # (m)
    drain_bnd_xmax = 1e6 # (m)

    # type of drain boundary
    # if True: convert drain nodes to specified head if hydraulic head
    # is above the surface
    # if False: regular flux-based drain boundary function
    seepage_bnd = True

    # option to run seepage bnd only once every x timesteps:
    seepage_bnd_timestep_interval = 1

    # full seepage iteration in 1 timestep or just remove inflow nodes for next timestep:
    iterate_seepage_in_one_timestep = False

    # option to stop recalculating seepage boundary after x time, to avoid
    # long computational time due to small oscillations
    seepage_bnd_max_time = 10000 * year

    # fluid source term
    Qf = 0  # (kg/(m^3 s)) (= fluid mass source)
    Qf_xmin = None
    Qf_xmax = None

    # solute source term
    Qs = 0
    Qs_xmin = None
    Qs_xmax = None

    ############################
    # initial conditions params
    ############################

    # add Ghyben-Herzberg initial salinity distribution
    # ignored if set to False
    ghyben_herzberg = False
    #sea_water_level = 0  # (m)
    #seawater_concentration = 0.035  # (kg/kg)
    #freshwater_concentration = 0.0  # (kg/kg)

    # derive initial pressure from analytical solution for 2D flow with
    # uniform recharge. somewhat obsolete, first run is now a steady-state run to calculate initial pressure & flow
    analytical_solution_initial_h = False

    # solve for pressure only using an initial steady state run
    # and use this as initial condition for subsequent runs
    initial_steady_state_run = True

    # 
    initial_concentration = 0.035

    #################
    # other params
    #################
    # permeability
    # Tiller-Flotten: K = 1e-9
    k = 10**(-16.0)
    #k = 10**(-12.0)  # (m^2) = median k

    # anisotropy (= horizontal permeability / vertical permeability)
    anisotropy = 10.0

    #
    porosity = 0.4  # (V/V)

    # specific storage
    specific_storage = 1.0e-4  # (...)

    # longitudinal dispersivity
    l_disp = 50.0  # (m)
    # transverse dispersivity ratio
    # used to calculate transverse dispersivity (transverse_disp = disp_ratio * l_disp)
    disp_ratio = 0.1  # (dimensionless)

    # settings for iterative solver solute & pressure PDEs following Ackerer (2004) GRL 31(2):
    pressure_convergence_criterion = 1.0e-4  # [Pa]
    concentration_convergence_criterion = 1.0e-7  # [kg/kg]
    # maximum iterations for sequential iterative solving of solute and
    # pressure equations
    min_iterations = 3
    max_iterations = 200

    # fluid viscosity
    # note, grompy can also automatically calculate viscosity from concentration and temperature data following
    # Batzle and Wang (1992). This option is controlled by calculate_viscosity in the iterate_coupled_flow_eqs function
    # in the module gwflow_lib
    calculate_viscosity = False
    viscosity = 8.94e-4
    reference_porosity = 0.25  # (non dimensional)
    reference_density = 1000.0  # [kg m-3]

    # eq. of state constants
    # density of water at solute concentration = 0 (and T=20 C)
    rho_f_0 = 998.872  # [kg m^-3]
    # thermal expansion coefficient
    alpha = 207e-6  # (1/K)
    # solute expansion coefficient
    gamma = 0.68412  # (dimensionless)

    # solute diffusion coefficient
    diffusivity = 3.9e-11  # (m^2/s^1)

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
    """

    k_s = 10**np.arange(-16.0, -12.0, 1.0) 

    pass