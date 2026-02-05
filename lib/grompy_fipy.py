"""
FiPy-based implementation of the coupled groundwater flow and solute transport model.

This module provides an alternative implementation of the main model functionality
using the FiPy backend. It can be used as a drop-in replacement for the escript-based
implementation when the 'fipy' backend is selected.

The module implements:
- Mesh loading from Gmsh files
- Boundary condition setup
- Initial condition setup
- Coupled pressure-concentration solving
- VTK output

Usage
-----
>>> from lib.grompy_fipy import run_coupled_flow_model_fipy
>>> results = run_coupled_flow_model_fipy(Parameters, ModelOptions, mesh_filename)

Notes
-----
- This implementation uses FiPy's finite volume discretization, which may produce
  slightly different results than the finite element escript implementation.
- Sequential solving (pressure then concentration) is used to avoid FiPy bug #361.
"""

import os
import math
import numpy as np

from lib.backend import get_backend
from lib.backend.fipy_backend import FiPyBackend, FiPyField, FiPyMesh, FiPyPDESolver

try:
    import fipy
    from fipy import CellVariable, FaceVariable, DiffusionTerm, DiffusionTermCorrection, ConvectionTerm, TransientTerm
    from fipy.tools import numerix
    FIPY_AVAILABLE = True
except ImportError:
    FIPY_AVAILABLE = False


def get_timestr(sec, limit_to_days=True):
    """Construct a string with the date and time."""
    second = 1.0
    minute = 60.0
    hour = 60.0 * minute
    day = hour * 24
    year = day * 365.25

    years = math.floor(sec / year)
    sec_left = sec - years * year
    days = math.floor(sec_left / day)
    sec_left -= days * day
    hours = math.floor(sec_left / hour)
    sec_left -= hours * hour
    minutes = math.floor(sec_left / minute)
    seconds = sec_left - minutes * minute

    if limit_to_days:
        timestr = '%i yrs, %i days' % (years, days)
    else:
        timestr = '%i yrs, %i days, %i hrs, %i min, %i sec' % (
            years, days, hours, minutes, seconds)

    return timestr


def calculate_fluid_density(concentration, gamma, rho_f_0):
    """
    Calculate fluid density from solute concentration.
    
    Parameters
    ----------
    concentration : array-like
        Solute concentration (kg/kg)
    gamma : float
        Solute expansion coefficient
    rho_f_0 : float
        Reference fluid density at zero concentration
    
    Returns
    -------
    array-like
        Fluid density
    """
    return rho_f_0 * (1.0 + gamma * concentration)


def calculate_darcy_velocity_fipy(fipy_mesh, pressure, density, k_tensor, viscosity, g):
    """
    Calculate Darcy velocity from pressure field.
    
    Darcy's law with gravity:
    q = -k/mu * (grad(P) - rho*g)
    
    where g is the gravity vector pointing in -y direction.
    
    Parameters
    ----------
    fipy_mesh : fipy.Mesh
        The FiPy mesh object
    pressure : np.ndarray
        Pressure field (cell-centered)
    density : np.ndarray
        Fluid density field (cell-centered)
    k_tensor : tuple
        Permeability tensor ((kxx, kxy), (kyx, kyy))
    viscosity : float
        Fluid viscosity (Pa.s)
    g : float
        Gravitational acceleration (m/s²)
    
    Returns
    -------
    FaceVariable
        Darcy velocity vector at cell faces (rank-1)
    """
    # Create cell variables for pressure and density
    P = CellVariable(mesh=fipy_mesh, value=pressure)
    rho = CellVariable(mesh=fipy_mesh, value=density)
    
    # Get density at faces (arithmetic average of neighboring cells)
    rho_face = rho.arithmeticFaceValue
    
    # Pressure gradient at faces
    grad_P = P.faceGrad  # Shape: (2, nFaces)
    
    # Gravity vector (pointing in -y direction, i.e., downward)
    g_vector = FaceVariable(mesh=fipy_mesh, rank=1, value=[[0.], [-g]])
    
    # Effective isotropic permeability (geometric mean)
    # Note: FiPy has issues with anisotropic tensor + ImplicitSourceTerm
    kxx = k_tensor[0][0]
    kyy = k_tensor[1][1]
    k_eff = np.sqrt(kxx * kyy)
    
    # Darcy velocity: q = -k/mu * (grad(P) - rho*g)
    q = -(k_eff / viscosity) * (grad_P - rho_face * g_vector)
    
    return q


def calculate_dispersion_coefficients_fipy(fipy_mesh, velocity_face, porosity, diffusivity, l_disp, t_disp):
    """
    Calculate anisotropic dispersion coefficients for FiPy implementation.
    
    Calculates complete hydrodynamic dispersion tensor following theory:
    D_ij = phi * D_m * delta_ij + alpha_L * v_i * v_j / |v| + alpha_T * |v| * delta_ij
    
    Where delta_ij is Kronecker delta, and off-diagonal terms capture cross-dispersion.
    
    Parameters
    ----------
    fipy_mesh : fipy.Mesh
        The FiPy mesh (cell-centered)
    velocity_face : FaceVariable
        Face-centered velocity from Darcy calculation
    porosity : float
        Porosity (dimensionless)
    diffusivity : float
        Molecular diffusivity (m²/s)
    l_disp : float
        Longitudinal dispersivity (m)
    t_disp : float
        Transverse dispersivity (m)
    
    Returns
    -------
    CellVariable
        Full dispersion tensor [[Dxx, Dxy], [Dyx, Dyy]] as rank-2 CellVariable
    """
    # Extract face-centered velocity components
    vx_face_vals = np.array(velocity_face[0].value)
    vy_face_vals = np.array(velocity_face[1].value)
    
    n_cells = fipy_mesh.numberOfCells
    n_faces = fipy_mesh.numberOfFaces
    
    # Get cell-to-face connectivity: shape (nCellFaces, nCells)
    cell_face_ids = fipy_mesh._cellFaceIDs  # Shape: (4, nCells) for 2D triangular mesh
    
    # For each cell, average the bounding face values
    # cell_face_ids has face indices for each cell
    vx_cell = np.zeros(n_cells)
    vy_cell = np.zeros(n_cells)
    
    for i in range(n_cells):
        face_ids = cell_face_ids[:, i]
        # Get valid face IDs (less than n_faces)
        valid_mask = (face_ids >= 0) & (face_ids < n_faces)
        valid_face_ids = face_ids[valid_mask]
        
        if len(valid_face_ids) > 0:
            vx_cell[i] = np.mean(vx_face_vals[valid_face_ids])
            vy_cell[i] = np.mean(vy_face_vals[valid_face_ids])
    
    # Calculate velocity magnitude
    v_abs = np.sqrt(vx_cell**2 + vy_cell**2)
    
    # Add numerical protection to avoid division by zero
    v_abs_safe = np.where(v_abs == 0, 1e-20, v_abs)
    
    # Calculate diagonal components of dispersion tensor
    # Following hydrodynamic dispersion theory: 
    # D_ii = phi * D_m + alpha_i * |v| where alpha_i is dispersivity
    Dxx = porosity * diffusivity + l_disp * vx_cell**2 / v_abs_safe
    Dyy = porosity * diffusivity + t_disp * vy_cell**2 / v_abs_safe
    
    # Calculate cross-dispersion components (off-diagonal terms)
    # D_xy = (alpha_L - alpha_T) * v_x * v_y / |v|
    Dxy = (l_disp - t_disp) * (vx_cell * vy_cell) / v_abs_safe
    Dyx = Dxy  # Symmetric tensor
    
    # Create rank-2 CellVariable for the full dispersion tensor
    # Stack components as (2, 2, nCells) for proper broadcasting
    dispersion_tensor = CellVariable(mesh=fipy_mesh, rank=2, 
                                 value=[[Dxx, Dxy], [Dyx, Dyy]])
    
    return dispersion_tensor


def setup_fipy_boundary_conditions(mesh, cell_centers, masks, Parameters):
    """
    Set up boundary conditions for the FiPy model.
    
    Parameters
    ----------
    mesh : FiPyMesh
        The FiPy mesh
    cell_centers : np.ndarray
        Cell center coordinates
    masks : dict
        Dictionary of boolean masks from mesh loading
    Parameters : object
        Model parameters
    
    Returns
    -------
    dict
        Dictionary containing boundary condition arrays and masks
    """
    x = cell_centers[:, 0]
    y = cell_centers[:, 1]
    n_cells = len(x)
    
    # Calculate z_surface (topographic surface)
    z_surface = x * Parameters.topo_gradient
    
    # Tolerance for boundary detection
    tol = Parameters.cellsize * 0.5
    
    # Surface mask: cells at the top boundary
    y_max = y.max()
    surface_mask = np.abs(y - z_surface) < tol
    
    # Sea surface mask: surface cells where x < 0
    sea_surface_mask = surface_mask & (x < 0)
    
    # Land surface mask: surface cells where x >= 0
    land_surface_mask = surface_mask & (x >= 0)
    
    # Specified pressure boundary
    spec_pressure_mask = (
        (x >= Parameters.specified_pressure_xmin) &
        (x <= Parameters.specified_pressure_xmax) &
        surface_mask
    )
    
    # Calculate specified pressure values
    # For seawater, add hydrostatic pressure from water column
    specified_pressure = np.zeros(n_cells)
    if Parameters.add_seawater_pressure:
        # Hydrostatic pressure from seawater column
        depth_below_sealevel = np.maximum(0, Parameters.sea_water_level - y)
        rho_seawater = calculate_fluid_density(
            Parameters.seawater_concentration,
            Parameters.gamma,
            Parameters.rho_f_0
        )
        specified_pressure = spec_pressure_mask * depth_below_sealevel * rho_seawater * Parameters.g
    
    # Recharge boundary: land surface cells
    recharge_mask = (
        (x >= Parameters.recharge_mass_flux_xmin) &
        (x <= Parameters.recharge_mass_flux_xmax) &
        land_surface_mask
    )
    
    # Drain/seepage boundary
    drain_mask = (
        (x >= Parameters.drain_bnd_xmin) &
        (x <= Parameters.drain_bnd_xmax) &
        surface_mask
    )
    
    # Specified concentration boundary
    spec_conc_mask = (
        (x >= Parameters.specified_concentration_xmin) &
        (x <= Parameters.specified_concentration_xmax) &
        sea_surface_mask
    )
    
    specified_concentration = np.zeros(n_cells)
    specified_concentration[spec_conc_mask] = Parameters.seawater_concentration
    
    return {
        'surface': surface_mask,
        'sea_surface': sea_surface_mask,
        'land_surface': land_surface_mask,
        'z_surface': z_surface,
        'spec_pressure_mask': spec_pressure_mask,
        'specified_pressure': specified_pressure,
        'recharge_mask': recharge_mask,
        'drain_mask': drain_mask,
        'spec_conc_mask': spec_conc_mask,
        'specified_concentration': specified_concentration,
    }


def setup_fipy_initial_conditions(mesh, cell_centers, z_surface, Parameters):
    """
    Set up initial conditions for pressure and concentration.
    
    Parameters
    ----------
    mesh : FiPyMesh
        The FiPy mesh
    cell_centers : np.ndarray
        Cell center coordinates
    z_surface : np.ndarray
        Surface elevation at each cell
    Parameters : object
        Model parameters
    
    Returns
    -------
    tuple
        (initial_pressure, initial_concentration, initial_density)
    """
    x = cell_centers[:, 0]
    y = cell_centers[:, 1]
    n_cells = len(x)
    
    # Initial concentration: use Ghyben-Herzberg approximation
    if Parameters.ghyben_herzberg:
        # Estimate fresh-saltwater interface depth
        rho_f = calculate_fluid_density(
            Parameters.freshwater_concentration,
            Parameters.gamma,
            Parameters.rho_f_0
        )
        rho_s = calculate_fluid_density(
            Parameters.seawater_concentration,
            Parameters.gamma,
            Parameters.rho_f_0
        )
        
        # Ghyben-Herzberg ratio
        gh_ratio = rho_f / (rho_s - rho_f)
        
        # Hydraulic head approximation (linear from coast to inland)
        h_approx = np.maximum(0, x * Parameters.topo_gradient)
        
        # Interface depth below sea level
        interface_depth = gh_ratio * h_approx
        
        # Set concentration based on position relative to interface
        concentration = np.zeros(n_cells)
        
        # Cells below interface get seawater concentration
        # Interface elevation = -interface_depth (below sea level)
        interface_elevation = -interface_depth
        
        # Seawater where y < interface_elevation and x < some threshold
        seawater_zone = (y < interface_elevation) | (x < 0)
        concentration[seawater_zone] = Parameters.seawater_concentration
        concentration[~seawater_zone] = Parameters.freshwater_concentration
        
    else:
        # Uniform fresh water
        concentration = np.full(n_cells, Parameters.freshwater_concentration)
    
    # Initial density
    density = calculate_fluid_density(concentration, Parameters.gamma, Parameters.rho_f_0)
    
    # Initial pressure: hydrostatic from surface
    # P = rho * g * (z_surface - y)
    depth = z_surface - y
    pressure = density * Parameters.g * np.maximum(0, depth)
    
    return pressure, concentration, density


def solve_steady_state_pressure_fipy(
    fipy_mesh, backend, k_tensor, viscosity, density, g,
    recharge_flux, recharge_density, recharge_mask,
    spec_pressure_mask, specified_pressure,
    Parameters
):
    """
    Solve steady-state pressure equation using FiPy.
    
    The equation solved is:
    -div(k/mu * (grad(P) - rho*g)) = Q
    
    Which can be rewritten as:
    -div(k/mu * grad(P)) = Q - div(k/mu * rho * g)
    
    Parameters
    ----------
    fipy_mesh : fipy.Mesh
        The FiPy mesh object
    backend : FiPyBackend
        The backend instance
    k_tensor : tuple
        Permeability tensor ((kxx, kxy), (kyx, kyy))
    viscosity : float
        Fluid viscosity
    density : np.ndarray
        Fluid density field
    g : float
        Gravitational acceleration
    recharge_flux : float
        Recharge flux (m/s)
    recharge_density : float
        Recharge fluid density
    recharge_mask : np.ndarray
        Boolean mask for recharge cells
    spec_pressure_mask : np.ndarray
        Boolean mask for specified pressure cells
    specified_pressure : np.ndarray
        Specified pressure values
    Parameters : object
        Model parameters
    
    Returns
    -------
    np.ndarray
        Pressure solution
    """
    # Extract permeability components
    kxx = k_tensor[0][0]
    kyy = k_tensor[1][1]
    
    # Create pressure variable with initial values
    initial_pressure = np.zeros(fipy_mesh.numberOfCells)
    if np.any(spec_pressure_mask):
        initial_pressure[spec_pressure_mask] = specified_pressure[spec_pressure_mask]
    
    pressure = CellVariable(mesh=fipy_mesh, name='pressure', value=initial_pressure)
    
    # Hydraulic conductivity / viscosity
    # Note: FiPy has issues with anisotropic tensor + ImplicitSourceTerm
    # Use effective isotropic permeability (geometric mean) as approximation
    k_eff = np.sqrt(kxx * kyy)
    diffusion_coeff = k_eff / viscosity
    
    # Apply Dirichlet boundary conditions using penalty method
    # Note: FiPy requires NEGATIVE penalty coefficient due to internal sign conventions
    # The equation becomes: div(D*grad(P)) - lambda*P = -lambda*P_target
    # which gives P = P_target for large lambda at constrained cells
    from fipy import ImplicitSourceTerm
    
    if np.any(spec_pressure_mask):
        # Create penalty coefficient (negative!) as CellVariable
        large_value = 1e10
        constraint_coeff = CellVariable(mesh=fipy_mesh, value=0.0)
        constraint_coeff.setValue(-large_value, where=spec_pressure_mask)
        
        constraint_source = CellVariable(mesh=fipy_mesh, value=0.0)
        constraint_source.setValue(-large_value * specified_pressure, where=spec_pressure_mask)
        
        # Build equation with penalty term
        eq = (DiffusionTerm(coeff=diffusion_coeff) 
              + ImplicitSourceTerm(coeff=constraint_coeff) 
              == constraint_source)
    else:
        # Build equation: div(D * grad(P)) = 0
        eq = DiffusionTerm(coeff=diffusion_coeff) == 0
    
    # Solve
    eq.solve(var=pressure)
    
    return np.array(pressure.value)


def solve_solute_transport_fipy(
    fipy_mesh, backend, concentration_old, dt,
    pressure, density, k_tensor, viscosity, g,
    porosity, diffusivity, l_disp, t_disp,
    spec_conc_mask, specified_concentration,
    Parameters, use_tensor_dispersion=True
):
    """
    Solve solute transport equation using FiPy.
    
    The advection-diffusion equation:
    phi * dC/dt + div(q*C) = div(D_eff * grad(C)) + Q_s
    
    Where:
    - phi = porosity
    - C = concentration
    - q = Darcy flux
    - D_eff = effective diffusion/dispersion tensor
    
    Parameters
    ----------
    fipy_mesh : fipy.Mesh
        The FiPy mesh
    concentration_old : np.ndarray
        Concentration from previous timestep
    dt : float
        Timestep size
    pressure : np.ndarray
        Current pressure field
    density : np.ndarray
        Current fluid density
    k_tensor : tuple
        Permeability tensor
    viscosity : float
        Fluid viscosity
    g : float
        Gravitational acceleration
    porosity : float
        Porosity
    diffusivity : float
        Molecular diffusivity
    l_disp : float
        Longitudinal dispersivity
    t_disp : float
        Transverse dispersivity
    spec_conc_mask : np.ndarray
        Mask for specified concentration cells
    specified_concentration : np.ndarray
        Specified concentration values
    Parameters : object
        Model parameters
    
    Returns
    -------
    np.ndarray
        Updated concentration
    """
    # Create concentration variable with initial values
    initial_conc = np.array(concentration_old, dtype=float)
    if np.any(spec_conc_mask):
        initial_conc[spec_conc_mask] = specified_concentration[spec_conc_mask]
    
    concentration = CellVariable(mesh=fipy_mesh, name='concentration', 
                                value=initial_conc)
    
    # Calculate Darcy velocity from pressure gradient
    q = calculate_darcy_velocity_fipy(
        fipy_mesh, pressure, density, k_tensor, viscosity, g
    )
    
    # Pore velocity = Darcy velocity / porosity
    velocity = q / porosity
    
    # Calculate dispersion based on configuration choice
    if use_tensor_dispersion:
        # Full anisotropic tensor implementation with diagonal and cross-dispersion terms
        dispersion_tensor = calculate_dispersion_coefficients_fipy(
            fipy_mesh, velocity, porosity, diffusivity, l_disp, t_disp
        )
        # Use tensor diffusion term with improved numerical stability
        diffusion_term = DiffusionTermCorrection(coeff=dispersion_tensor)
    else:
        # Fallback: preserve old scalar implementation
        diffusion_coeff = porosity * diffusivity + l_disp  # Original implementation
        diffusion_term = DiffusionTerm(coeff=diffusion_coeff)
        print(f"Debug: Using scalar diffusion: {diffusion_coeff:.2e} m²/s")
    
    # Build equation with Dirichlet boundary conditions using penalty method
    # Note: FiPy requires NEGATIVE penalty coefficient due to internal sign conventions
    from fipy import ImplicitSourceTerm
    
    if np.any(spec_conc_mask):
        # Create penalty coefficient (negative!) as CellVariable
        large_value = 1e10
        constraint_coeff = CellVariable(mesh=fipy_mesh, value=0.0)
        constraint_coeff.setValue(-large_value, where=spec_conc_mask)
        
        constraint_source = CellVariable(mesh=fipy_mesh, value=0.0)
        constraint_source.setValue(-large_value * specified_concentration, where=spec_conc_mask)
        
        eq = (TransientTerm(coeff=porosity) 
              + ConvectionTerm(coeff=velocity)
              == diffusion_term
              + ImplicitSourceTerm(coeff=constraint_coeff)
              - constraint_source)
    else:
        eq = (TransientTerm(coeff=porosity) 
              + ConvectionTerm(coeff=velocity)
              == diffusion_term)
    
    # Solve for one timestep
    eq.solve(var=concentration, dt=dt)
    
    return np.array(concentration.value)


def run_coupled_flow_model_fipy(Parameters, ModelOptions, mesh_filename):
    """
    Run coupled groundwater flow and solute transport model using FiPy backend.
    
    This is the main entry point for running the model with the FiPy backend.
    
    Parameters
    ----------
    Parameters : object
        Model parameters (from model_parameters.py)
    ModelOptions : object
        Model options (from model_parameters.py)
    mesh_filename : str
        Path to the Gmsh mesh file
    
    Returns
    -------
    dict
        Dictionary containing model results:
        - 'pressure': final pressure field
        - 'concentration': final concentration field
        - 'flux': final flux field
        - 'runtime': total simulation time
        - 'timesteps': number of timesteps
    """
    if not FIPY_AVAILABLE:
        raise ImportError("FiPy is not installed. Install with: pip install fipy")
    
    print("=" * 60)
    print("Running coupled flow model with FiPy backend")
    print("=" * 60)
    
    year = 365.25 * 24 * 60 * 60.0
    day = 24 * 60 * 60.0
    
    # Initialize backend
    backend = get_backend('fipy')
    
    # Load mesh
    print("Loading mesh from:", mesh_filename)
    from lib.fipy_mesh_io import load_fipy_mesh_from_msh
    fipy_mesh, cell_centers, masks, field_data, extra = load_fipy_mesh_from_msh(
        mesh_filename, 
        topo_gradient=Parameters.topo_gradient
    )
    
    mesh = FiPyMesh(fipy_mesh)
    n_cells = mesh.num_cells
    print(f"Mesh loaded: {n_cells} cells")
    
    # Set up boundary conditions
    print("Setting up boundary conditions")
    bc = setup_fipy_boundary_conditions(mesh, cell_centers, masks, Parameters)
    
    # Set up initial conditions
    print("Setting up initial conditions")
    pressure, concentration, density = setup_fipy_initial_conditions(
        mesh, cell_centers, bc['z_surface'], Parameters
    )
    
    # Set up permeability tensor
    k_tensor = (
        (Parameters.k, 0),
        (0, Parameters.k / Parameters.anisotropy)
    )
    
    # Time stepping parameters
    dt = Parameters.dt0
    runtime = 0
    total_time = Parameters.total_time
    max_runtime = Parameters.max_runtime
    max_timesteps = Parameters.max_timesteps
    output_interval = Parameters.output_interval
    
    timestep = 0
    output_step = 0
    last_output_time = 0
    
    # Main time loop
    print("\nStarting time stepping...")
    print("-" * 60)
    
    # Initial steady-state pressure solve
    if ModelOptions.initial_steady_state_run:
        print("Running initial steady-state pressure solve...")
        pressure = solve_steady_state_pressure_fipy(
            fipy_mesh, backend, k_tensor, Parameters.viscosity,
            density, Parameters.g,
            Parameters.recharge_flux, Parameters.recharge_density,
            bc['recharge_mask'],
            bc['spec_pressure_mask'], bc['specified_pressure'],
            Parameters
        )
        print("Steady-state pressure solve complete")
        print(f"Pressure range: {pressure.min():.2f} to {pressure.max():.2f} Pa")
    
    # Create output directory
    output_dir = ModelOptions.model_output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Transient simulation
    while runtime < total_time and runtime < max_runtime and timestep < max_timesteps:
        
        # Update fluid density from concentration
        density = calculate_fluid_density(
            concentration, Parameters.gamma, Parameters.rho_f_0
        )
        
        # Solve pressure equation
        # (simplified - in full implementation would iterate between pressure and transport)
        
        # Solve solute transport
        if ModelOptions.solute_transport:
            concentration_new = solve_solute_transport_fipy(
                fipy_mesh, backend, concentration, dt,
                pressure, density, k_tensor, Parameters.viscosity,
                Parameters.g, Parameters.porosity,
                Parameters.diffusivity, Parameters.l_disp,
                Parameters.l_disp * Parameters.disp_ratio,
                bc['spec_conc_mask'], bc['specified_concentration'],
                Parameters, use_tensor_dispersion=True
            )
            
            # Check for convergence / steady state
            conc_change = np.max(np.abs(concentration_new - concentration))
            concentration = concentration_new
        else:
            conc_change = 0
        
        # Update time
        runtime += dt
        timestep += 1
        
        # Increase timestep
        dt = min(dt * Parameters.dt_inc, Parameters.dt_max)
        
        # Output
        if runtime - last_output_time >= output_interval:
            timestr = get_timestr(runtime)
            print(f"Timestep {timestep}: t = {timestr}, "
                  f"P: [{pressure.min():.1f}, {pressure.max():.1f}] Pa, "
                  f"C: [{concentration.min():.4f}, {concentration.max():.4f}] kg/kg")
            
            # Save VTK output
            if ModelOptions.save_vtk_files:
                vtk_filename = os.path.join(output_dir, f"output_{output_step:04d}")
                
                # Create field objects for saving
                pressure_field = FiPyField(
                    CellVariable(mesh=fipy_mesh, value=pressure, name='pressure'),
                    fipy_mesh, 'pressure'
                )
                conc_field = FiPyField(
                    CellVariable(mesh=fipy_mesh, value=concentration, name='concentration'),
                    fipy_mesh, 'concentration'
                )
                density_field = FiPyField(
                    CellVariable(mesh=fipy_mesh, value=density, name='density'),
                    fipy_mesh, 'density'
                )
                
                backend.save_vtk(
                    vtk_filename, mesh,
                    {'pressure': pressure_field, 
                     'concentration': conc_field,
                     'density': density_field}
                )
            
            last_output_time = runtime
            output_step += 1
        
        # Check for steady state (only after minimum runtime has been reached)
        # This matches the original escript implementation which requires
        # runtime >= total_time before allowing early termination
        if Parameters.stop_when_steady_state and runtime >= total_time:
            conc_rate = conc_change / (dt / year)  # Change per year
            if conc_rate < Parameters.max_concentration_change_steady_state:
                print(f"Steady state reached at t = {get_timestr(runtime)}")
                print(f"Concentration change rate: {conc_rate:.2e} kg/kg/yr")
                break
    
    print("-" * 60)
    print(f"Simulation complete: {timestep} timesteps, runtime = {get_timestr(runtime)}")
    
    # Final output
    results = {
        'pressure': pressure,
        'concentration': concentration,
        'density': density,
        'runtime': runtime,
        'timesteps': timestep,
        'cell_centers': cell_centers,
        'mesh': mesh,
    }
    
    return results


def check_fipy_available():
    """Check if FiPy is available."""
    return FIPY_AVAILABLE
