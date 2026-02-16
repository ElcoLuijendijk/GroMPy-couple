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
    from fipy import (
        UpwindConvectionTerm,
        PowerLawConvectionTerm,
        VanLeerConvectionTerm,
        ExponentialConvectionTerm
    )
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


def pressure_to_fresh_water_head(P, rho_f_0, g, y_coords):
    """
    Convert pressure field to equivalent fresh water hydraulic head.
    
    SEAWAT formulation (Langevin et al. 2007):
    
    h_f = P / (ρ_f * g) + y
    
    where:
    - h_f = equivalent fresh water hydraulic head (m)
    - P = pressure (Pa)
    - ρ_f = fresh water reference density (kg/m³), constant
    - g = gravitational acceleration (m/s²)
    - y = elevation (m)
    
    This transformation allows solving the flow equation directly for h_f,
    with all gravity effects eliminated from the differential equation.
    
    References:
    - Langevin et al. (2007) "SEAWAT Version 4: A Computer Program for 
      Simulation of Multi-Species Solute and Heat Transport"
      USGS Techniques and Methods 6-A7
    
    Parameters
    ----------
    P : np.ndarray
        Pressure field (Pa), shape (n_cells,)
    rho_f_0 : float
        Fresh water reference density (kg/m³)
    g : float
        Gravitational acceleration (m/s²)
    y_coords : np.ndarray
        Cell center y-coordinates (elevation, m), shape (n_cells,)
    
    Returns
    -------
    np.ndarray
        Fresh water hydraulic head (m), shape (n_cells,)
    """
    h_f = P / (rho_f_0 * g) + y_coords
    return h_f


def fresh_water_head_to_pressure(h_f, rho_f_0, g, y_coords):
    """
    Convert equivalent fresh water hydraulic head to pressure.
    
    Inverse of pressure_to_fresh_water_head():
    
    P = (h_f - y) * ρ_f * g
    
    This is used to convert the fresh water head solution back to pressure
    for output, analysis, and next iteration initialization.
    
    Parameters
    ----------
    h_f : np.ndarray
        Equivalent fresh water hydraulic head (m), shape (n_cells,)
    rho_f_0 : float
        Fresh water reference density (kg/m³)
    g : float
        Gravitational acceleration (m/s²)
    y_coords : np.ndarray
        Cell center y-coordinates (elevation, m), shape (n_cells,)
    
    Returns
    -------
    np.ndarray
        Pressure field (Pa), shape (n_cells,)
    """
    P = (h_f - y_coords) * rho_f_0 * g
    return P

 
def calculate_darcy_velocity_fipy(fipy_mesh, pressure, density, k_tensor, viscosity, g, cell_centers, rho_f_0):
    """
    Calculate Darcy velocity using SEAWAT fresh water head formulation.
    
    SEAWAT Approach (Langevin et al. 2007):
    
    Using equivalent fresh water head h_f = P/(ρ_f*g) + z, the Darcy flux becomes:
    
        q = -ρ * (k/μ) * ∇h_f
    
    where:
    - ρ = current density (variable with salt concentration)
    - k = permeability tensor
    - μ = fluid viscosity
    - ∇h_f = gradient of fresh water head
    
    Key insight: Gravity is completely eliminated from the flow equation.
    Buoyancy effects are implicitly handled through density variations in the 
    diffusion coefficient, not through an explicit gravity term in velocity.
    
    This is more stable numerically than the traditional formulation:
        q_traditional = -k/μ * (∇P - ρg)
    
    Parameters
    ----------
    fipy_mesh : fipy.Mesh
        The FiPy mesh object
    pressure : np.ndarray
        Pressure field (Pa), cell-centered
    density : np.ndarray
        Fluid density field (kg/m³), cell-centered
    k_tensor : tuple
        Permeability tensor ((kxx, kxy), (kyx, kyy))
    viscosity : float
        Fluid viscosity (Pa·s)
    g : float
        Gravitational acceleration (m/s²)
    cell_centers : np.ndarray
        Cell center coordinates, shape (n_cells, 2), used for elevation
    rho_f_0 : float
        Fresh water reference density (kg/m³)
    
    Returns
    -------
    FaceVariable
        Darcy velocity vector at cell faces (rank-1, m/s)
        
    References:
    -----------
    Langevin et al. (2007) SEAWAT Version 4 documentation.
    """
    # Get cell elevations
    y_coords = cell_centers[:, 1]
    
    # Convert pressure to fresh water head
    h_f = pressure_to_fresh_water_head(pressure, rho_f_0, g, y_coords)
    
    # Create FiPy variables
    h_f_var = CellVariable(mesh=fipy_mesh, value=h_f)
    rho = CellVariable(mesh=fipy_mesh, value=density)
    
    # Get density at faces (arithmetic average)
    rho_face = rho.arithmeticFaceValue
    
    # Gradient of fresh water head at faces
    grad_h_f = h_f_var.faceGrad  # Shape: (2, nFaces)
    
    # Effective isotropic permeability (geometric mean)
    kxx = k_tensor[0][0]
    kyy = k_tensor[1][1]
    k_eff = np.sqrt(kxx * kyy)
    
    # Darcy velocity with density-weighted coefficient:
    # q = -ρ * (k/μ) * ∇h_f
    # This includes all buoyancy effects implicitly through ρ
    q = -rho_face * (k_eff / viscosity) * grad_h_f
    
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



def calculate_boundary_fluxes_fipy(cell_centers, pressure, density, k_tensor, 
                                   viscosity, g, bc_dict, Parameters, year):
    """Calculate boundary fluxes matching escript implementation."""
    x = cell_centers[:, 0]
    y = cell_centers[:, 1]
    
    kxx = k_tensor[0][0]
    kyy = k_tensor[1][1]
    k_eff = np.sqrt(kxx * kyy)
    
    surface_mask = bc_dict['surface']
    sea_surface_mask = bc_dict['sea_surface']
    land_surface_mask = bc_dict['land_surface']
    
    n_cells = len(x)
    qx = np.zeros(n_cells)
    qy = np.zeros(n_cells)
    
    for i in range(n_cells):
        if surface_mask[i]:
            rho_local = density[i]
            qy[i] = -(k_eff / viscosity) * (0.0 - rho_local * g)
    
    flux_surface_norm = np.array([qx, qy])
    
    submarine_mask = sea_surface_mask & (x < 0)
    submarine_flux = flux_surface_norm[1] * submarine_mask.astype(float)
    submarine_flux_in = submarine_flux * (flux_surface_norm[1] < 0)
    submarine_flux_out = submarine_flux * (flux_surface_norm[1] > 0)
    
    total_submarine_flux = np.sum(submarine_flux)
    total_submarine_flux_in = np.sum(submarine_flux_in)
    total_submarine_flux_out = np.sum(submarine_flux_out)
    
    active_seepage_bnd = bc_dict.get('active_seepage_mask', np.zeros(n_cells))
    seepage_flux = flux_surface_norm[1] * active_seepage_bnd.astype(float)
    total_seepage_flux = np.sum(seepage_flux)
    
    land_mask = land_surface_mask & (x >= 0)
    land_flux = flux_surface_norm[1] * land_mask.astype(float)
    land_flux_in = land_flux * (flux_surface_norm[1] < 0)
    land_flux_out = land_flux * (flux_surface_norm[1] > 0)
    
    total_land_flux_in = np.sum(land_flux_in)
    total_land_flux_out = np.sum(land_flux_out)
    
    recharge_mask = bc_dict.get('recharge_mask', np.zeros(n_cells))
    recharge_flux = flux_surface_norm[1] * recharge_mask.astype(float)
    total_rch_flux = np.sum(recharge_flux)
    
    total_flux_over_surface_norm = np.array([0.0, total_land_flux_in + total_land_flux_out + total_submarine_flux])
    
    land_flux_vals = flux_surface_norm[1][land_mask]
    seepage_flux_vals = flux_surface_norm[1][active_seepage_bnd > 0]
    submarine_flux_vals = flux_surface_norm[1][submarine_mask]
    
    min_land_flux = np.min(land_flux_vals) if len(land_flux_vals) > 0 else 0.0
    max_land_flux = np.max(land_flux_vals) if len(land_flux_vals) > 0 else 0.0
    min_seepage_flux = np.min(seepage_flux_vals) if len(seepage_flux_vals) > 0 else 0.0
    max_seepage_flux = np.max(seepage_flux_vals) if len(seepage_flux_vals) > 0 else 0.0
    min_submarine_flux = np.min(submarine_flux_vals) if len(submarine_flux_vals) > 0 else 0.0
    max_submarine_flux = np.max(submarine_flux_vals) if len(submarine_flux_vals) > 0 else 0.0
    
    flux_buffer = 1.0e-1 / year
    
    outflow_mask = flux_surface_norm[1] > 0
    inflow_mask = flux_surface_norm[1] < 0
    
    inflow_land = inflow_mask & (x > 0)
    outflow_land = outflow_mask & (x > 0)
    inflow_sea = inflow_mask & (x <= 0)
    outflow_sea = outflow_mask & (x <= 0)
    
    ext_inflow_land = np.sum(inflow_land)
    ext_outflow_land = np.sum(outflow_land)
    ext_inflow_sea = np.sum(inflow_sea)
    ext_outflow_sea = np.sum(outflow_sea)
    
    outflow_land_threshold = outflow_mask & (x > 0) & (flux_surface_norm[1] > flux_buffer)
    outflow_sea_threshold = outflow_mask & (x <= 0) & (flux_surface_norm[1] > flux_buffer)
    
    ext_outflow_land_threshold = np.sum(outflow_land_threshold)
    ext_outflow_sea_threshold = np.sum(outflow_sea_threshold)
    
    boundary_fluxes = [flux_surface_norm, land_flux, land_flux_out, submarine_flux, submarine_flux_in, submarine_flux_out]
    
    boundary_flux_stats = [total_flux_over_surface_norm, total_rch_flux, total_seepage_flux,
                           total_land_flux_in, total_land_flux_out, total_submarine_flux,
                           total_submarine_flux_in, total_submarine_flux_out,
                           0, 0, 0, 0,
                           ext_inflow_land, ext_outflow_land, ext_inflow_sea, ext_outflow_sea,
                           ext_outflow_land_threshold, ext_outflow_sea_threshold,
                           min_land_flux, max_land_flux, min_seepage_flux, max_seepage_flux,
                           min_submarine_flux, max_submarine_flux]
    
    return boundary_fluxes, boundary_flux_stats

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
    
    # Helper function to get parameter with default
    def get_param(name, default=None):
        return getattr(Parameters, name, default)
    
    # Calculate z_surface (topographic surface)
    z_surface = x * get_param('topo_gradient', 0.0)
    
    # Tolerance for boundary detection
    cellsize = get_param('cellsize', None)
    if cellsize is None:
        # Calculate approximate cellsize from mesh
        x_range = x.max() - x.min()
        y_range = y.max() - y.min()
        estimated_area = (x_range * y_range) / len(x)
        cellsize = np.sqrt(estimated_area)
    tol = cellsize * 0.5
    
    # Surface mask: cells at the top boundary
    y_max = y.max()
    surface_mask = np.abs(y - z_surface) < tol
    
    # Sea surface mask: surface cells where x < 0
    sea_surface_mask = surface_mask & (x < 0)
    
    # Land surface mask: surface cells where x >= 0
    land_surface_mask = surface_mask & (x >= 0)
    
    # Specified pressure boundary
    spec_pressure_xmin = get_param('specified_pressure_xmin', [-10000])
    spec_pressure_xmax = get_param('specified_pressure_xmax', [0.01])
    
    # Handle both list and scalar values
    if isinstance(spec_pressure_xmin, (list, tuple)):
        spec_pressure_xmin = spec_pressure_xmin[0] if spec_pressure_xmin else -10000
    if isinstance(spec_pressure_xmax, (list, tuple)):
        spec_pressure_xmax = spec_pressure_xmax[0] if spec_pressure_xmax else 0.01
    
    spec_pressure_mask = (
        (x >= spec_pressure_xmin) &
        (x <= spec_pressure_xmax) &
        surface_mask
    )
    
    # Calculate specified pressure values
     # For seawater, add hydrostatic pressure from water column
    specified_pressure = np.zeros(n_cells)
    if getattr(Parameters, 'add_seawater_pressure', False):
        # Hydrostatic pressure from seawater column
        depth_below_sealevel = np.maximum(0, Parameters.sea_water_level - y)
        rho_seawater = calculate_fluid_density(
            Parameters.seawater_concentration,
            Parameters.gamma,
            Parameters.rho_f_0
        )
        specified_pressure = spec_pressure_mask * depth_below_sealevel * rho_seawater * Parameters.g
    
    # Recharge boundary: land surface cells
    recharge_xmin = get_param('recharge_mass_flux_xmin', 1e-6)
    recharge_xmax = get_param('recharge_mass_flux_xmax', 1e6)
    recharge_mask = (
        (x >= recharge_xmin) &
        (x <= recharge_xmax) &
        land_surface_mask
    )
    
    # Drain/seepage boundary
    drain_xmin = get_param('drain_bnd_xmin', -0.001)
    drain_xmax = get_param('drain_bnd_xmax', 1e6)
    drain_mask = (
        (x >= drain_xmin) &
        (x <= drain_xmax) &
        surface_mask
    )
    
    # Specified concentration boundary - loop through all regions
    spec_conc_xmin = get_param('specified_concentration_xmin', [-0.01])
    spec_conc_xmax = get_param('specified_concentration_xmax', [1e6])
    spec_conc_ymin = get_param('specified_concentration_ymin', [-0.01])
    spec_conc_ymax = get_param('specified_concentration_ymax', [1e6])
    spec_conc_values = get_param('specified_concentration', [0.0])
    specified_concentration_surface = get_param('specified_concentration_surface', False)
    
    # Ensure all are lists
    if not isinstance(spec_conc_xmin, (list, tuple)):
        spec_conc_xmin = [spec_conc_xmin]
    if not isinstance(spec_conc_xmax, (list, tuple)):
        spec_conc_xmax = [spec_conc_xmax]
    if not isinstance(spec_conc_ymin, (list, tuple)):
        spec_conc_ymin = [spec_conc_ymin]
    if not isinstance(spec_conc_ymax, (list, tuple)):
        spec_conc_ymax = [spec_conc_ymax]
    if not isinstance(spec_conc_values, (list, tuple)):
        spec_conc_values = [spec_conc_values]
    
    specified_concentration = np.zeros(n_cells)
    spec_conc_mask_regional = np.zeros(n_cells, dtype=bool)  # Track which cells get BC
    
    # Tolerance for boundary matching (allow cells within one cell size of boundary)
    unique_x = np.unique(x)
    unique_y = np.unique(y)
    if len(unique_x) > 1:
        dx = np.mean(np.diff(unique_x))
    else:
        dx = 0.001
    if len(unique_y) > 1:
        dy = np.mean(np.diff(unique_y))
    else:
        dy = 0.001
    tol = max(dx, dy)
    
    # Loop through all boundary regions
    for i_region in range(len(spec_conc_xmin)):
        xmin = spec_conc_xmin[i_region]
        xmax = spec_conc_xmax[i_region]
        ymin = spec_conc_ymin[i_region]
        ymax = spec_conc_ymax[i_region]
        conc_val = spec_conc_values[i_region] if i_region < len(spec_conc_values) else 0.0
        
        # Create mask for this region
        if specified_concentration_surface:
            # If applying to entire surface, use sea_surface_mask as before
            region_mask = (
                (x >= xmin) &
                (x <= xmax) &
                (y >= ymin) &
                (y <= ymax) &
                sea_surface_mask
            )
        else:
            # If applying to specific x,y regions, don't require sea_surface_mask
            # For vertical boundaries (xmin == xmax), use tolerance
            if abs(xmax - xmin) < tol:
                x_match = np.abs(x - xmin) <= tol
            else:
                x_match = (x >= xmin) & (x <= xmax)
            
            # For horizontal boundaries (ymin == ymax), use tolerance
            if abs(ymax - ymin) < tol:
                y_match = np.abs(y - ymin) <= tol
            else:
                y_match = (y >= ymin) & (y <= ymax)
            
            region_mask = x_match & y_match & surface_mask
        
        # Apply boundary condition for this region
        specified_concentration[region_mask] = conc_val
        spec_conc_mask_regional[region_mask] = True  # Mark cells that received BC
    
    # Create combined mask for return value
    # NOTE: We use the broader mask here for numerical stability with FiPy penalty method
    # The penalty method works better when applied to all zero-initialized cells
    spec_conc_mask = specified_concentration > -9999
    
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
        concentration[seawater_zone] = getattr(Parameters, 'seawater_concentration', 0.035)
        concentration[~seawater_zone] = getattr(Parameters, 'freshwater_concentration', 0.0)
        
    else:
        # Uniform fresh water
        freshwater_conc = getattr(Parameters, 'freshwater_concentration', 0.0)
        concentration = np.full(n_cells, freshwater_conc)
    
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


def solve_transient_pressure_fipy(
    fipy_mesh, backend, k_tensor, viscosity, density, g,
    recharge_flux, recharge_mask,
    spec_pressure_mask, specified_pressure,
    Parameters, dt, pressure_old, cell_centers, porosity=None, gamma=None,
    concentration_old=None, dC_dt_field=None
):
    """
    Solve transient flow equation using SEAWAT equivalent fresh water head formulation.
    
    SEAWAT Fresh Water Head Approach (Langevin et al. 2007):
    =========================================================
    
    Instead of solving for pressure P with explicit gravity terms in the velocity,
    we solve for equivalent fresh water hydraulic head h_f:
    
        h_f = P / (ρ_f * g) + z
    
    where ρ_f is the constant fresh water reference density.
    
    The flow equation in terms of h_f is:
    
        S_s * ∂h_f/∂t - div((ρ/ρ_f) * k/μ * ∇h_f) = Q_eff
    
    where:
    - S_s = specific storage coefficient
    - ρ = current density (varies with salt concentration)
    - ρ_f = reference fresh water density (constant)
    - k/μ = hydraulic conductivity / viscosity
    - Q_eff = recharge converted to head units
    
    Key advantages:
    1. Gravity completely eliminated from differential equation
    2. Density variations appear as coefficient (ρ/ρ_f) in diffusion term
    3. Numerically simpler and more stable than pressure formulation
    4. Industry standard used by SEAWAT, FEFLOW, and other codes
    
    Parameters
    ----------
    fipy_mesh : fipy.Mesh
        The FiPy mesh object
    backend : FiPyBackend
        The backend instance
    k_tensor : tuple
        Permeability tensor ((kxx, kxy), (kyx, kyy))
    viscosity : float
        Fluid viscosity (Pa·s)
    density : np.ndarray
        Current fluid density field (kg/m³), shape (n_cells,)
    g : float
        Gravitational acceleration (m/s²)
    recharge_flux : float
        Recharge flux (m/s)
    recharge_mask : np.ndarray
        Boolean mask for recharge cells
    spec_pressure_mask : np.ndarray
        Boolean mask for cells with specified pressure BCs
    specified_pressure : np.ndarray
        Specified pressure values (Pa) at BC cells
    Parameters : object
        Model parameters object
    dt : float
        Time step (s)
    pressure_old : np.ndarray
        Pressure field from previous timestep/iteration (Pa)
    cell_centers : np.ndarray
        Cell center coordinates, shape (n_cells, 2), used for elevation
    porosity : float, optional
        Porosity (not used in this formulation, kept for compatibility)
    gamma : float, optional
        Density expansion coefficient (not used in this formulation)
    concentration_old : np.ndarray, optional
        Concentration from previous iteration (not used, kept for compatibility)
    dC_dt_field : np.ndarray, optional
        Rate of concentration change (not used, kept for compatibility)
    
    Returns
    -------
    np.ndarray
        Pressure solution (Pa), converted back from fresh water head
        
    References:
    -----------
    Langevin, C. D., Thorne, D. T., Dausman, A. M., Sukop, M. C., & Guo, W. (2007).
    SEAWAT Version 4: A Computer Program for Simulation of Multi-Species Solute and 
    Heat Transport. U.S. Geological Survey Techniques and Methods 6-A7.
    """
    from fipy import TransientTerm, ImplicitSourceTerm
    
    # Extract permeability components
    kxx = k_tensor[0][0]
    kyy = k_tensor[1][1]
    k_eff = np.sqrt(kxx * kyy)
    
    # Get cell elevations (y-coordinates)
    y_coords = cell_centers[:, 1]
    
    # Convert pressure_old to fresh water head for initial guess
    h_f_old = pressure_to_fresh_water_head(pressure_old, Parameters.rho_f_0, g, y_coords)
    
    # Create fresh water head variable (this is what we solve for)
    h_f = CellVariable(mesh=fipy_mesh, name='fresh_water_head', value=h_f_old)
    
    # Storage coefficient from parameters
    S_s = Parameters.specific_storage
    
    # Effective diffusion coefficient with density-dependent permeability
    # D_eff = (ρ/ρ_f) * (k/μ)
    # This is where density variations are incorporated into the flow equation
    rho_rel = density / Parameters.rho_f_0  # Relative density ratio
    diffusion_coeff_array = rho_rel * k_eff / viscosity
    
    # Create as CellVariable for FiPy to handle properly
    diffusion_coeff = CellVariable(mesh=fipy_mesh, value=diffusion_coeff_array)
    
    # Recharge source term converted to head units
    # Q_eff = Q / (ρ_f * g)  [flux in head units, m/s]
    q_source = CellVariable(mesh=fipy_mesh, value=0.0)
    if np.any(recharge_mask):
        q_source.setValue(recharge_flux, where=recharge_mask)
    
    # Storage coefficient as CellVariable
    S_s_cell = CellVariable(mesh=fipy_mesh, value=S_s)
    
    # Convert specified pressure BCs to fresh water head
    specified_head = np.zeros_like(specified_pressure)
    if np.any(spec_pressure_mask):
        specified_head[spec_pressure_mask] = pressure_to_fresh_water_head(
            specified_pressure[spec_pressure_mask],
            Parameters.rho_f_0, g,
            y_coords[spec_pressure_mask]
        )
    
    # Apply Dirichlet boundary conditions using penalty method
    if np.any(spec_pressure_mask):
        large_value = 1e10
        constraint_coeff = CellVariable(mesh=fipy_mesh, value=0.0)
        constraint_coeff.setValue(-large_value, where=spec_pressure_mask)
        
        constraint_source = CellVariable(mesh=fipy_mesh, value=0.0)
        constraint_source.setValue(-large_value * specified_head, where=spec_pressure_mask)
        
        # Flow equation with penalty method for Dirichlet BC:
        # S_s * ∂h_f/∂t - div((ρ/ρ_f) * k/μ * ∇h_f) - λ*h_f = -λ*h_target + Q_eff
        eq = (TransientTerm(coeff=S_s_cell / dt) +
              DiffusionTerm(coeff=diffusion_coeff) +
              ImplicitSourceTerm(coeff=constraint_coeff)
              == q_source + constraint_source)
    else:
        # Flow equation without penalty:
        # S_s * ∂h_f/∂t - div((ρ/ρ_f) * k/μ * ∇h_f) = Q_eff
        eq = (TransientTerm(coeff=S_s_cell / dt) +
              DiffusionTerm(coeff=diffusion_coeff)
              == q_source)
    
    # Solve for fresh water head
    eq.solve(var=h_f, dt=dt)
    
    # Debug output: show relative density range
    rho_max = np.max(density)
    rho_min = np.min(density)
    if rho_max > rho_min:
        print(f"    Density range: {rho_min:.2f} to {rho_max:.2f} kg/m³ (rel: {rho_min/Parameters.rho_f_0:.6f} to {rho_max/Parameters.rho_f_0:.6f})")
    
    # Convert fresh water head solution back to pressure for output and next iteration
    h_f_solution = np.array(h_f.value)
    pressure_solution = fresh_water_head_to_pressure(h_f_solution, Parameters.rho_f_0, g, y_coords)
    
    return pressure_solution


def get_convection_term(velocity, scheme='exponential'):
    """
    Create a convection term using the specified FiPy advection scheme.
    
    Parameters
    ----------
    velocity : FaceVariable
        The velocity field (Darcy velocity / porosity)
    scheme : str, optional
        Convection scheme to use. Options:
        - 'central': ConvectionTerm (standard, can oscillate at high Pe)
        - 'upwind': UpwindConvectionTerm (stable, diffusive)
        - 'powerlaw': PowerLawConvectionTerm (good for Pe < 100)
        - 'vanleer': VanLeerConvectionTerm (TVD, excellent for sharp interfaces)
        - 'exponential': ExponentialConvectionTerm (default, exact for 1D, best for high Pe)
    
    Returns
    -------
    ConvectionTerm-like object
        The appropriate convection term for the advection-diffusion equation
    """
    scheme_lower = scheme.lower().strip()
    
    if scheme_lower == 'upwind':
        return UpwindConvectionTerm(coeff=velocity)
    elif scheme_lower == 'powerlaw':
        return PowerLawConvectionTerm(coeff=velocity)
    elif scheme_lower == 'vanleer':
        return VanLeerConvectionTerm(coeff=velocity)
    elif scheme_lower == 'exponential':
        return ExponentialConvectionTerm(coeff=velocity)
    elif scheme_lower == 'central':
        return ConvectionTerm(coeff=velocity)
    else:
        # Default to Exponential
        print(f"Warning: Unknown convection scheme '{scheme}', using 'exponential'")
        return ExponentialConvectionTerm(coeff=velocity)


def solve_solute_transport_fipy(
    fipy_mesh, backend, concentration_old, dt,
    pressure, density, k_tensor, viscosity, g,
    porosity, diffusivity, l_disp, t_disp,
    spec_conc_mask, specified_concentration,
    Parameters, cell_centers, use_tensor_dispersion=True, convection_scheme='exponential'
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
    cell_centers : np.ndarray
        Cell center coordinates, shape (n_cells, 2), used for elevation in fresh water head calculation
    use_tensor_dispersion : bool
        Whether to use full tensor dispersion
    convection_scheme : str
        Convection scheme to use
    
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
        fipy_mesh, pressure, density, k_tensor, viscosity, g, cell_centers, Parameters.rho_f_0
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
        
        convection_term = get_convection_term(velocity, scheme=convection_scheme)
        eq = (TransientTerm(coeff=porosity) 
              + convection_term
              == diffusion_term
              + ImplicitSourceTerm(coeff=constraint_coeff)
              - constraint_source)
    else:
        convection_term = get_convection_term(velocity, scheme=convection_scheme)
        eq = (TransientTerm(coeff=porosity) 
              + convection_term
              == diffusion_term)
    
    # Solve for one timestep
    eq.solve(var=concentration, dt=dt)
    
    return np.array(concentration.value)


def run_coupled_flow_model_fipy(Parameters, ModelOptions, mesh_filename, convection_scheme='exponential'):
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
    convection_scheme : str, optional
        Convection scheme for solute transport. Options:
        - 'exponential': ExponentialConvectionTerm (default, best for high Pe numbers)
        - 'upwind': UpwindConvectionTerm (stable, diffusive)
        - 'powerlaw': PowerLawConvectionTerm (good for Pe < 100)
        - 'vanleer': VanLeerConvectionTerm (TVD, excellent for sharp interfaces)
        - 'central': ConvectionTerm (standard, can oscillate at high Pe)
    
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
    # Initialize tracking arrays for metrics
    pressure_prev = pressure.copy()
    concentration_prev = concentration.copy()
    dts = []
    runtimes = []
    pressure_differences_max = []
    pressure_differences_mean = []
    concentration_differences_max = []
    concentration_differences_mean = []
    reached_steady_state = False
    
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
        
        # Enforce concentration boundary conditions after steady-state pressure solve
        spec_conc_xmin = bc['spec_conc_xmin'] if 'spec_conc_xmin' in bc else Parameters.specified_concentration_xmin
        spec_conc_xmax = bc['spec_conc_xmax'] if 'spec_conc_xmax' in bc else Parameters.specified_concentration_xmax
        spec_conc_ymin = Parameters.specified_concentration_ymin
        spec_conc_ymax = Parameters.specified_concentration_ymax
        spec_conc_values = Parameters.specified_concentration
        
        # Ensure all are lists
        if not isinstance(spec_conc_xmin, (list, tuple)):
            spec_conc_xmin = [spec_conc_xmin]
        if not isinstance(spec_conc_xmax, (list, tuple)):
            spec_conc_xmax = [spec_conc_xmax]
        if not isinstance(spec_conc_ymin, (list, tuple)):
            spec_conc_ymin = [spec_conc_ymin]
        if not isinstance(spec_conc_ymax, (list, tuple)):
            spec_conc_ymax = [spec_conc_ymax]
        if not isinstance(spec_conc_values, (list, tuple)):
            spec_conc_values = [spec_conc_values]
        
        # Get x and y cell centers
        x = cell_centers[:, 0]
        y = cell_centers[:, 1]
        
        # Tolerance for boundary matching (allow cells within one cell size of boundary)
        tol = max(np.mean(np.diff(np.sort(np.unique(x)))), 
                  np.mean(np.diff(np.sort(np.unique(y)))))
        
        # Apply concentration boundary conditions
        for i_region in range(len(spec_conc_xmin)):
            xmin = spec_conc_xmin[i_region]
            xmax = spec_conc_xmax[i_region]
            ymin = spec_conc_ymin[i_region]
            ymax = spec_conc_ymax[i_region]
            conc_val = spec_conc_values[i_region] if i_region < len(spec_conc_values) else 0.0
            
            # Create mask for this region with tolerance
            # For boundaries at exact locations (xmin == xmax), use tolerance to find nearby cells
            if abs(xmax - xmin) < tol:
                # This is a vertical line boundary - match cells within tolerance
                mask = (
                    (np.abs(x - xmin) <= tol) &
                    (y >= ymin) &
                    (y <= ymax)
                )
            else:
                # This is a region - match cells within bounds
                mask = (
                    (x >= xmin) &
                    (x <= xmax) &
                    (y >= ymin) &
                    (y <= ymax)
                )
            
            # Apply boundary condition
            concentration[mask] = conc_val
        
        print(f"Concentration range after BC enforcement: {concentration.min():.6f} to {concentration.max():.6f}")
    
    # Create output directory
    output_dir = ModelOptions.model_output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Transient simulation
    while runtime < total_time and runtime < max_runtime and timestep < max_timesteps:
        
        # Save state at start of timestep for convergence checking
        concentration_start_ts = concentration.copy()
        pressure_start_ts = pressure.copy()
        
        # Get convergence parameters
        max_coupled_iterations = Parameters.max_iterations if hasattr(Parameters, 'max_iterations') else 20
        pressure_tol = Parameters.pressure_convergence_criterion if hasattr(Parameters, 'pressure_convergence_criterion') else 1.0e-5
        concentration_tol = Parameters.concentration_convergence_criterion if hasattr(Parameters, 'concentration_convergence_criterion') else 1.0e-5
        
        # COUPLED ITERATION LOOP (NEW!)
        coupled_converged = False
        coupled_iter = 0
        print(f"    Starting coupled iterations (max={max_coupled_iterations})...")
        for coupled_iter in range(max_coupled_iterations):
            
            # [1] Update fluid density from current concentration
            density = calculate_fluid_density(
                concentration, Parameters.gamma, Parameters.rho_f_0
            )
            
            # [2] Calculate Darcy velocity from current pressure
            velocity_face = calculate_darcy_velocity_fipy(
                fipy_mesh, pressure, density, k_tensor, Parameters.viscosity, Parameters.g, cell_centers, Parameters.rho_f_0
            )
            
            # [3] Solve solute transport with updated velocity
            if ModelOptions.solute_transport:
                    concentration_new = solve_solute_transport_fipy(
                        fipy_mesh, backend, concentration, dt,
                        pressure, density, k_tensor, Parameters.viscosity,
                        Parameters.g, Parameters.porosity,
                        Parameters.diffusivity, Parameters.l_disp,
                        Parameters.l_disp * Parameters.disp_ratio,
                        bc['spec_conc_mask'], bc['specified_concentration'],
                        Parameters, cell_centers, use_tensor_dispersion=True, convection_scheme=convection_scheme
                    )
                    concentration = concentration_new
            
            # [3.5] Calculate concentration change rate (dC/dt) for coupling
            # This is used in the pressure equation to account for density-driven flow
            dC_dt_field = (concentration - concentration_start_ts) / dt
            
            # [4] Solve transient pressure with updated density and coupling
            # (This is where coupling happens - density changes affect pressure)
            pressure_new = solve_transient_pressure_fipy(
                fipy_mesh, backend, k_tensor, Parameters.viscosity,
                density, Parameters.g,
                Parameters.recharge_flux, bc['recharge_mask'],
                bc['spec_pressure_mask'], bc['specified_pressure'],
                Parameters, dt, pressure_start_ts, cell_centers,
                porosity=Parameters.porosity,
                gamma=Parameters.gamma,
                concentration_old=concentration_start_ts,
                dC_dt_field=dC_dt_field
            )
            
            # [5] Check convergence
            pressure_change = np.max(np.abs(pressure_new - pressure)) if len(pressure_new) > 0 else 0.0
            conc_change = np.max(np.abs(concentration_new - concentration)) if len(concentration_new) > 0 else 0.0
            
            pressure = pressure_new
            
            if (pressure_change < pressure_tol and conc_change < concentration_tol):
                coupled_converged = True
                if coupled_iter > 0:
                    print(f"  Coupled iteration {coupled_iter+1}: CONVERGED (P_Δ={pressure_change:.2e}, C_Δ={conc_change:.2e})")
                break
            else:
                if coupled_iter < 3 or (coupled_iter + 1) % 5 == 0:  # Print first 3 and every 5th
                    print(f"  Coupled iteration {coupled_iter+1}/{max_coupled_iterations}: P_Δ={pressure_change:.2e}, C_Δ={conc_change:.2e}")
        
        if not coupled_converged and coupled_iter == max_coupled_iterations - 1:
            print(f"  WARNING: Coupled iterations did not converge after {max_coupled_iterations} iterations")
        
        # Track differences
        pressure_diff = np.abs(pressure - pressure_prev)
        concentration_diff = np.abs(concentration - concentration_prev)
        
        pressure_differences_max.append(np.max(pressure_diff) if len(pressure_diff) > 0 else 0.0)
        pressure_differences_mean.append(np.mean(pressure_diff) if len(pressure_diff) > 0 else 0.0)
        concentration_differences_max.append(np.max(concentration_diff) if len(concentration_diff) > 0 else 0.0)
        concentration_differences_mean.append(np.mean(concentration_diff) if len(concentration_diff) > 0 else 0.0)
        
        pressure_prev = pressure.copy()
        concentration_prev = concentration.copy()
        
        # Update time
        runtime += dt
        timestep += 1
        
        # Increase timestep
        
        # Track timestep and runtime
        dts.append(dt)
        runtimes.append(runtime)
        
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
                reached_steady_state = True
                break
    
    print("-" * 60)
    print(f"Simulation complete: {timestep} timesteps, runtime = {get_timestr(runtime)}")
    
    # Final output
    # Compute final derived fields
    z_surface = bc['z_surface']
    h = (pressure / (density * Parameters.g)) + z_surface
    
    q_face = calculate_darcy_velocity_fipy(fipy_mesh, pressure, density, k_tensor, 
                                            Parameters.viscosity, Parameters.g, cell_centers, Parameters.rho_f_0)
    
    # Interpolate velocity from faces to cell centers
    # FiPy FaceVariable -> CellVariable conversion
    if hasattr(q_face, 'cellCenters'):
         # Get the cell-centered values from face values
         q_cell = np.array([q_face[0].cellCenters, q_face[1].cellCenters])
    else:
         # Fallback: use arithmetic average of neighboring faces
         q_vals = np.zeros((2, n_cells))
         for i in range(2):  # x and y components
             face_vals = q_face[i].value
             for cell_idx in range(n_cells):
                 # Get faces connected to this cell and average
                 cell_faces = fipy_mesh.cellFaceIDs[:, cell_idx]
                 q_vals[i, cell_idx] = np.mean(face_vals[cell_faces])
         q_cell = q_vals
    
    q = q_cell
    
    q_abs = np.sqrt(q[0]**2 + q[1]**2)
    nodal_flux = q[1] * bc['surface'].astype(float)
    
    Pdiff = np.array(pressure_differences_max[-1] if pressure_differences_max else 0.0)
    Cdiff = np.array(concentration_differences_max[-1] if concentration_differences_max else 0.0)
    
    dts = np.array(dts)
    runtimes = np.array(runtimes)
    pressure_differences_max = np.array(pressure_differences_max)
    concentration_differences_max = np.array(concentration_differences_max)
    pressure_differences_mean = np.array(pressure_differences_mean)
    concentration_differences_mean = np.array(concentration_differences_mean)
    
    boundary_fluxes, boundary_flux_stats = calculate_boundary_fluxes_fipy(
        cell_centers, pressure, density, k_tensor, 
        Parameters.viscosity, Parameters.g, bc, Parameters, year
    )
    
    boundary_conditions = [
        bc.get('spec_pressure_mask', np.zeros(n_cells)),
        bc.get('specified_pressure', np.zeros(n_cells)),
        bc.get('spec_conc_mask', np.zeros(n_cells)),
        bc.get('active_seepage_mask', np.zeros(n_cells)),
        bc.get('specified_concentration', np.zeros(n_cells)),
        bc.get('specified_concentration_rho_f', Parameters.rho_f_0),
        bc.get('recharge_mask', np.zeros(n_cells)),
        bc.get('active_seepage_mask', np.zeros(n_cells))
    ]
    
    results = (mesh, bc['surface'], bc.get('sea_surface', None), 
               k_tensor, pressure, concentration,
               Parameters.rho_f_0, Parameters.viscosity, h, q, q_abs, nodal_flux,
               Pdiff, Cdiff,
               pressure_differences_max, concentration_differences_max,
               pressure_differences_mean, concentration_differences_mean,
               dts, runtimes, timestep, output_step,
               boundary_conditions, boundary_fluxes, boundary_flux_stats,
               reached_steady_state)
    
    return results


def check_fipy_available():
    """Check if FiPy is available."""
    return FIPY_AVAILABLE
