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
    from fipy import LinearLUSolver
    FIPY_AVAILABLE = True
except ImportError:
    FIPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# MPI / parallel helpers
# ---------------------------------------------------------------------------
# Detect whether FiPy is running under MPI (i.e. launched via mpirun -np N).
# fipy.parallel.Nproc > 1 only when petsc4py is installed AND we are actually
# running inside an MPI job.  In serial / non-PETSc environments this block
# is a no-op and _FIPY_PARALLEL stays False.
try:
    from fipy import parallel as _fipy_parallel_comm
    _FIPY_PARALLEL = _fipy_parallel_comm.Nproc > 1
except (ImportError, AttributeError):
    _FIPY_PARALLEL = False


def _get_value(var):
    """Return the full value array of a FiPy variable as a numpy array.

    In serial runs (or when PETSc/MPI is not available) this is equivalent to
    ``np.array(var.value)``.  In MPI parallel runs ``var.globalValue`` performs
    an all-gather so that **every** rank receives the complete assembled array,
    avoiding the partial-partition data that ``var.value`` would return.
    """
    if _FIPY_PARALLEL:
        return np.array(var.globalValue)
    return np.array(var.value)


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

 
def calculate_darcy_flux_fipy(fipy_mesh, pressure, density, k_tensor, viscosity, g, cell_centers, rho_f_0):
    """
    Calculate face-centred Darcy velocity:

        q = −(k/μ) · (∇P − ρ · g · ẑ)

    where ρ is face-interpolated fluid density, P is total pressure, and ẑ is
    the upward unit vector (so the gravity term subtracts from the vertical
    pressure gradient, driving flow from high to low hydraulic head).

    Parameters
    ----------
    fipy_mesh : fipy.Mesh
    pressure : np.ndarray
        Total pressure field (Pa), shape (n_cells,)
    density : np.ndarray
        Fluid density (kg/m³), shape (n_cells,)
    k_tensor : tuple
        Permeability tensor ((kxx, kxy), (kyx, kyy))
    viscosity : float
        Fluid viscosity (Pa·s)
    g : float
        Gravitational acceleration (m/s²)
    cell_centers : np.ndarray
        Unused; kept for backward-compatible call signature.
    rho_f_0 : float
        Unused; kept for backward-compatible call signature.

    Returns
    -------
    FaceVariable
        Face-centred Darcy flux, shape (2, n_faces), m/s
    """
    kxx = k_tensor[0][0]
    kyy = k_tensor[1][1]
    k_eff = np.sqrt(kxx * kyy)
    mobility = k_eff / viscosity

    # Compute ∇P at faces via FiPy's built-in face-gradient interpolation.
    P_var = CellVariable(mesh=fipy_mesh, value=pressure)
    grad_P_face = P_var.faceGrad   # FaceVariable, shape (2, n_faces)

    # Face-interpolate density  ρ_face  (arithmetic average, numpy)
    _fci   = np.array(fipy_mesh.faceCellIDs)
    _own   = _fci[0]
    _nb    = _fci[1]
    _valid = _nb >= 0
    rho_face = density[_own].copy()
    rho_face[_valid] = 0.5 * (density[_own[_valid]] + density[_nb[_valid]])

    # q = −(k/μ) · (∇P − ρ_face · g · ẑ)
    qx_vals = -mobility * _get_value(grad_P_face[0])
    qy_vals = -mobility * (_get_value(grad_P_face[1]) - rho_face * g)

    return FaceVariable(mesh=fipy_mesh, value=np.array([qx_vals, qy_vals]))



def calculate_dispersion_coefficients_fipy(fipy_mesh, darcy_flux_face, porosity, diffusivity, l_disp, t_disp):
    """
    Calculate anisotropic dispersion coefficients for FiPy implementation.

    Uses the standard hydrodynamic dispersion tensor (Bear 1972):

        D_ij = (phi * D_m + alpha_T * |v|) * delta_ij
                + (alpha_L - alpha_T) * v_i * v_j / |v|

    where v = q/phi is the pore (seepage) velocity, |v| its magnitude, and
    q is the Darcy flux.  The factor phi is absorbed into the alpha terms so
    that the diffusion term in the PDE ``div(D · grad C)`` gives the correct
    macroscopic dispersion without an extra phi factor.

    Parameters
    ----------
    fipy_mesh : fipy.Mesh
        The FiPy mesh (cell-centered)
    darcy_velocity_face : FaceVariable
        Face-centered Darcy velocity q (NOT pore velocity).
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
    # Extract face-centered Darcy velocity components
    qx_face_vals = _get_value(darcy_flux_face[0])
    qy_face_vals = _get_value(darcy_flux_face[1])

    n_cells = fipy_mesh.numberOfCells
    n_faces = fipy_mesh.numberOfFaces

    # Get cell-to-face connectivity: shape (nCellFaces, nCells)
    cell_face_ids = fipy_mesh._cellFaceIDs

    # Vectorised face→cell average
    valid = (cell_face_ids >= 0) & (cell_face_ids < n_faces)   # (F, N) bool
    safe_ids = np.where(valid, cell_face_ids, 0)               # clamp invalids to 0
    qx_all = qx_face_vals[safe_ids]                            # (F, N)
    qy_all = qy_face_vals[safe_ids]                            # (F, N)
    valid_count = valid.sum(axis=0).astype(float)              # (N,)
    valid_count = np.where(valid_count == 0, 1.0, valid_count)
    qx_cell = (qx_all * valid).sum(axis=0) / valid_count       # cell-centred Darcy vel
    qy_cell = (qy_all * valid).sum(axis=0) / valid_count

    # Pore velocity v = q / phi
    vx_cell = qx_cell / porosity
    vy_cell = qy_cell / porosity

    # Magnitude of pore velocity
    v_abs = np.sqrt(vx_cell**2 + vy_cell**2)
    v_abs_safe = np.where(v_abs == 0, 1e-20, v_abs)

    # Correct hydrodynamic dispersion tensor (Bear 1972):
    #   D_xx = phi*Dm + alpha_T*|v| + (alpha_L - alpha_T)*vx^2 / |v|
    #   D_yy = phi*Dm + alpha_T*|v| + (alpha_L - alpha_T)*vy^2 / |v|
    #   D_xy = D_yx = (alpha_L - alpha_T)*vx*vy / |v|
    base = porosity * diffusivity + t_disp * v_abs
    Dxx = base + (l_disp - t_disp) * vx_cell**2 / v_abs_safe
    Dyy = base + (l_disp - t_disp) * vy_cell**2 / v_abs_safe
    Dxy = (l_disp - t_disp) * (vx_cell * vy_cell) / v_abs_safe
    Dyx = Dxy  # symmetric

    # Create rank-2 CellVariable for the full dispersion tensor
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
    # For flat domains (topo_gradient == 0), z_surface is the domain thickness
    # (the actual top boundary elevation), NOT zero.
    # Using z_surface = x * topo_gradient = 0 for a flat domain would incorrectly
    # identify the bottom boundary (y approx 0) as the surface.
    topo_gradient = get_param('topo_gradient', 0.0)
    thickness = get_param('thickness', y.max())
    if topo_gradient == 0.0:
        z_surface = np.full(n_cells, thickness)
    else:
        z_surface = x * topo_gradient

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
    surface_mask = np.abs(y - z_surface) < tol
    
    # Sea surface mask: surface cells where x < 0
    sea_surface_mask = surface_mask & (x < 0)
    
    # Land surface mask: surface cells where x >= 0
    land_surface_mask = surface_mask & (x >= 0)
    
    # Specified pressure boundary
    # Supports multiple regions (lists), vertical walls (xmin == xmax),
    # and the specified_pressure_surface = False path (ymin/ymax ranges).
    spec_pressure_xmin_raw = get_param('specified_pressure_xmin', [-10000])
    spec_pressure_xmax_raw = get_param('specified_pressure_xmax', [0.01])
    spec_pressure_ymin_raw = get_param('specified_pressure_ymin', [0.0])
    spec_pressure_ymax_raw = get_param('specified_pressure_ymax', [1e6])
    spec_pressure_values_raw = get_param('specified_pressure', [0.0])
    spec_pressure_surface_flag = get_param('specified_pressure_surface', True)
    spec_pressure_salinity_raw = get_param('specified_pressure_salinity', None)

    # Normalise to lists
    def _to_list(v):
        return list(v) if isinstance(v, (list, tuple)) else [v]

    spec_pressure_xmin_list = _to_list(spec_pressure_xmin_raw)
    spec_pressure_xmax_list = _to_list(spec_pressure_xmax_raw)
    spec_pressure_ymin_list = _to_list(spec_pressure_ymin_raw)
    spec_pressure_ymax_list = _to_list(spec_pressure_ymax_raw)
    spec_pressure_values_list = _to_list(spec_pressure_values_raw)
    if spec_pressure_salinity_raw is not None:
        spec_pressure_salinity_list = _to_list(spec_pressure_salinity_raw)
    else:
        spec_pressure_salinity_list = None

    specified_pressure = np.zeros(n_cells)
    spec_pressure_mask = np.zeros(n_cells, dtype=bool)

    n_pressure_regions = len(spec_pressure_xmin_list)
    for i_p in range(n_pressure_regions):
        pxmin = spec_pressure_xmin_list[i_p]
        pxmax = spec_pressure_xmax_list[i_p] if i_p < len(spec_pressure_xmax_list) else spec_pressure_xmin_list[i_p]
        pymin = spec_pressure_ymin_list[i_p] if i_p < len(spec_pressure_ymin_list) else 0.0
        pymax = spec_pressure_ymax_list[i_p] if i_p < len(spec_pressure_ymax_list) else 1e6
        pval = spec_pressure_values_list[i_p] if i_p < len(spec_pressure_values_list) else 0.0

        is_vertical = abs(pxmax - pxmin) < tol

        if spec_pressure_surface_flag and not is_vertical:
            # Apply only to surface cells in the x-range
            p_region_mask = (x >= pxmin) & (x <= pxmax) & surface_mask
        else:
            # Apply to full y-range (vertical walls, or surface flag off)
            if is_vertical:
                x_match = np.abs(x - pxmin) <= tol
            else:
                x_match = (x >= pxmin) & (x <= pxmax)
            y_match = (y >= pymin) & (y <= pymax)
            p_region_mask = x_match & y_match

        specified_pressure[p_region_mask] = pval
        spec_pressure_mask |= p_region_mask

        # Add hydrostatic pressure column below the top of the BC zone.
        # Mirrors the escript backend (grompy_lib.py lines 364-368):
        #   dPh = (ymax - y) * rho_segment * g
        # The specified_pressure value is the pressure at the TOP of the BC
        # zone (y == pymax); each cell lower in the column gets an additional
        # hydrostatic contribution.  The fluid density for this BC segment is
        # derived from its specified salinity (spec_pressure_salinity).
        if spec_pressure_salinity_list is not None and i_p < len(spec_pressure_salinity_list):
            sal_seg = spec_pressure_salinity_list[i_p]
            rho_seg = calculate_fluid_density(sal_seg, Parameters.gamma, Parameters.rho_f_0)
        else:
            # Fallback: use reference fresh water density
            rho_seg = Parameters.rho_f_0
        depth_below_top = np.maximum(0.0, pymax - y)
        specified_pressure[p_region_mask] += depth_below_top[p_region_mask] * rho_seg * Parameters.g

    # Add hydrostatic seawater pressure on top of specified values if requested
    if getattr(Parameters, 'add_seawater_pressure', False):
        depth_below_sealevel = np.maximum(0, Parameters.sea_water_level - y)
        rho_seawater = calculate_fluid_density(
            Parameters.seawater_concentration,
            Parameters.gamma,
            Parameters.rho_f_0
        )
        specified_pressure += spec_pressure_mask * depth_below_sealevel * rho_seawater * Parameters.g
    
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
    
    # Tolerance for boundary matching.
    # On unstructured (triangular) meshes, np.unique(x) has O(nCells) entries
    # with mean spacing ~L/nCells → O(10 µm), far smaller than a cell diameter.
    # Using the cellsize parameter (same as the pressure-BC section above) gives
    # a robust tolerance that catches the first layer of cells near each wall.
    # tol is already set above (line ~504) from cellsize; reuse it here.
    
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
            # If applying to specific x,y regions (vertical or horizontal boundaries)
            # For vertical boundaries (xmin == xmax), use tolerance to identify left/right walls
            if abs(xmax - xmin) < tol:
                x_match = np.abs(x - xmin) <= tol
            else:
                x_match = (x >= xmin) & (x <= xmax)
            
            # For horizontal boundaries (ymin == ymax), use tolerance to identify top/bottom
            if abs(ymax - ymin) < tol:
                y_match = np.abs(y - ymin) <= tol
            else:
                y_match = (y >= ymin) & (y <= ymax)
            
            # CRITICAL FIX: Don't require surface_mask for vertical boundaries (x fixed)
            # Vertical boundaries (like seawater inlet at x=0) need to apply to full depth
            # Only apply surface_mask for horizontal boundaries (y fixed)
            is_vertical_boundary = abs(xmax - xmin) < tol
            if is_vertical_boundary:
                # Vertical boundary: apply to full y-range without surface requirement
                region_mask = x_match & y_match
            else:
                # Horizontal boundary: apply to surface only
                region_mask = x_match & y_match & surface_mask
        
        # Apply boundary condition for this region
        specified_concentration[region_mask] = conc_val
        spec_conc_mask_regional[region_mask] = True  # Mark cells that received BC
        
        # Debug: Print how many cells received this BC
        n_cells_bc = np.sum(region_mask)
        if n_cells_bc > 0:
            print(f"  Region {i_region}: C={conc_val:.5f}, x=[{xmin:.4f}, {xmax:.4f}], "
                  f"y=[{ymin:.4f}, {ymax:.4f}], {n_cells_bc} cells")
    
    # Create combined mask for return value
    # Only mark cells that explicitly received a concentration BC (cells near the boundary walls).
    # Using a broader mask (e.g. specified_concentration > -9999, which is True for ALL cells)
    # would force C=0 everywhere via the penalty method and suppress salt-water intrusion.
    spec_conc_mask = spec_conc_mask_regional
    
    # Debug output: Boundary condition summary
    print("\n" + "="*70)
    print("BOUNDARY CONDITION SUMMARY")
    print("="*70)
    print(f"Surface cells: {np.sum(surface_mask)}")
    print(f"Sea surface cells: {np.sum(sea_surface_mask)}")
    print(f"Land surface cells: {np.sum(land_surface_mask)}")
    print(f"Pressure BC cells: {np.sum(spec_pressure_mask)}")
    print(f"Concentration BC cells: {np.sum(spec_conc_mask_regional)}")
    if np.sum(spec_conc_mask_regional) > 0:
        print(f"  Concentration BC range: {specified_concentration[spec_conc_mask_regional].min():.6f} to {specified_concentration[spec_conc_mask_regional].max():.6f}")
    print("="*70 + "\n")
    
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


def setup_fipy_initial_conditions(mesh, cell_centers, z_surface, Parameters, bc=None):
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
    bc : dict, optional
        Boundary condition dictionary containing 'specified_concentration' and 'spec_conc_mask'
    
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
    
    # OPTION B: Apply concentration boundary conditions to initial condition
    # This seeds salt water at inlet boundaries for faster intrusion
    if bc is not None and 'specified_concentration' in bc and 'spec_conc_mask' in bc:
        spec_conc_mask = bc['spec_conc_mask']
        specified_concentration = bc['specified_concentration']
        
        # Apply BC values to initial concentration at specified cells
        if np.any(spec_conc_mask):
            concentration[spec_conc_mask] = specified_concentration[spec_conc_mask]
            n_bc_cells = np.sum(spec_conc_mask)
            print(f"Applied initial concentration BC to {n_bc_cells} cells")
            print(f"  Concentration range: {concentration.min():.6f} to {concentration.max():.6f}")
    
    # Initial density
    density = calculate_fluid_density(concentration, Parameters.gamma, Parameters.rho_f_0)
    
    # Initial pressure: hydrostatic from surface
    # P = rho * g * (z_surface - y)
    depth = z_surface - y
    pressure = density * Parameters.g * np.maximum(0, depth)
    
    return pressure, concentration, density


def _build_face_pressure_bc(fipy_mesh, spec_pressure_mask, specified_pressure):
    """
    Build a (face_mask, face_pressure_values) pair for FiPy's .constrain() API.

    For each exterior face adjacent to a BC cell the total pressure value (Pa)
    is read directly from ``specified_pressure`` at the adjacent BC cell centre.
    The face is within ~dx/2 of the cell centre so the error is O(dx).

    Parameters
    ----------
    fipy_mesh : fipy.Mesh
    spec_pressure_mask : np.ndarray, bool, shape (n_cells,)
        True for cells that carry a specified-pressure Dirichlet BC.
    specified_pressure : np.ndarray, float, shape (n_cells,)
        Total pressure values (Pa) at BC cells.

    Returns
    -------
    bc_face_mask : np.ndarray, bool, shape (n_faces,)
    bc_face_vals : np.ndarray, float, shape (n_faces,)
    """
    face_cell_ids = np.array(fipy_mesh.faceCellIDs)   # (2, n_faces)
    n_faces = face_cell_ids.shape[1]

    # Identify exterior faces (works on both Grid2D and Gmsh2D meshes).
    try:
        ext_face_arr = np.array(fipy_mesh.exteriorFaces.value, dtype=bool)
    except Exception:
        ext_face_arr = np.zeros(n_faces, dtype=bool)

    owners    = face_cell_ids[0]   # shape (n_faces,)
    neighbors = face_cell_ids[1]   # shape (n_faces,); may be -1 or == owner on Gmsh BCs

    owner_is_bc = spec_pressure_mask[owners]

    # Clamp negative neighbour indices before indexing (masked out below).
    nb_safe = np.where(neighbors >= 0, neighbors, 0)
    nb_valid = neighbors >= 0
    nb_is_bc = nb_valid & spec_pressure_mask[nb_safe]

    # A face is a BC face iff exterior AND at least one adjacent cell is a BC cell.
    bc_face_mask = ext_face_arr & (owner_is_bc | nb_is_bc)

    # Set face pressure values directly from the adjacent BC cell's total pressure.
    bc_face_vals = np.zeros(n_faces)
    use_owner = bc_face_mask & owner_is_bc
    bc_face_vals[use_owner] = specified_pressure[owners[use_owner]]
    use_nb = bc_face_mask & (~owner_is_bc) & nb_is_bc
    bc_face_vals[use_nb] = specified_pressure[nb_safe[use_nb]]

    return bc_face_mask, bc_face_vals


def solve_steady_state_pressure_fipy(
    fipy_mesh, backend, k_tensor, viscosity, density, g,
    recharge_flux, recharge_density, recharge_mask,
    spec_pressure_mask, specified_pressure,
    Parameters, cell_centers=None
):
    """
    Solve the steady-state variable-density groundwater flow equation for
    total pressure P (Pa) using FiPy.

    Governing equation (Bear 1972 / SEAWAT):

        ∇·(ρ · k/μ · ∇P) = ∇·(ρ² · k/μ · g · ẑ)

    In FiPy operator form (moving the RHS source to the right):

        DiffusionTerm(ρ·k/μ) == ConvectionTerm((0, ρ²·k/μ·g)) + Q_recharge

    where the ConvectionTerm coefficient is a FaceVariable representing the
    gravity "velocity" (0, ρ²·k/μ·g) that drives the body-force flux.

    Dirichlet pressure BCs are applied via FiPy's .constrain() on the exterior
    faces adjacent to BC cells.  The specified_pressure array already contains
    physically correct total pressures in Pa (built with the full hydrostatic
    column in setup_fipy_boundary_conditions).

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
        Fluid density field (kg/m³), shape (n_cells,)
    g : float
        Gravitational acceleration (m/s²)
    recharge_flux : float
        Recharge source term (kg/(m³·s) or equivalent volumetric rate)
    recharge_density : float
        Density of recharge fluid (not used directly; kept for API compatibility)
    recharge_mask : np.ndarray
        Boolean mask for recharge cells
    spec_pressure_mask : np.ndarray
        Boolean mask for specified-pressure Dirichlet BC cells
    specified_pressure : np.ndarray
        Total pressure values (Pa) at BC cells
    Parameters : object
        Model parameters (must have .rho_f_0)
    cell_centers : np.ndarray, optional
        Unused; kept for backward-compatible call signature.

    Returns
    -------
    np.ndarray
        Total pressure solution P (Pa), shape (n_cells,)
    """
    # ------------------------------------------------------------------
    # Permeability / mobility
    # ------------------------------------------------------------------
    kxx = k_tensor[0][0]
    kyy = k_tensor[1][1]
    k_eff = np.sqrt(kxx * kyy)
    mobility = k_eff / viscosity          # scalar k/μ

    # ------------------------------------------------------------------
    # Face-interpolate density  ρ_face  (arithmetic average, numpy)
    # Avoids FiPy arithmeticFaceValue edge-cases on Gmsh meshes.
    # ------------------------------------------------------------------
    _fci   = np.array(fipy_mesh.faceCellIDs)   # (2, n_faces)
    _own   = _fci[0]
    _nb    = _fci[1]
    _valid = _nb >= 0
    rho_face = density[_own].copy()
    rho_face[_valid] = 0.5 * (density[_own[_valid]] + density[_nb[_valid]])

    # ------------------------------------------------------------------
    # DiffusionTerm coefficient:  ρ_face · k/μ  (FaceVariable)
    # ------------------------------------------------------------------
    diffusion_fv = FaceVariable(mesh=fipy_mesh, value=rho_face * mobility)

    # ------------------------------------------------------------------
    # Gravity source via ConvectionTerm:  ∇·(ρ²·k/μ · g · ẑ)
    # The ConvectionTerm in FiPy implements  ∇·(u · φ)  where u is the
    # coefficient FaceVariable.  Setting u = (0, ρ²·k/μ·g) with φ = 1
    # gives exactly the divergence of the gravity flux we need on the RHS.
    # ------------------------------------------------------------------
    n_faces = fipy_mesh.numberOfFaces
    gravity_vals = np.zeros((2, n_faces))
    gravity_vals[1] = rho_face**2 * mobility * g
    gravity_fv = FaceVariable(mesh=fipy_mesh, rank=1, value=gravity_vals)

    # ------------------------------------------------------------------
    # Recharge source term (cell-centred CellVariable)
    # ------------------------------------------------------------------
    recharge_array = np.zeros(fipy_mesh.numberOfCells)
    if np.any(recharge_mask):
        recharge_array[recharge_mask] = recharge_flux
    recharge_source = CellVariable(mesh=fipy_mesh, value=recharge_array)

    # ------------------------------------------------------------------
    # Pressure CellVariable + Dirichlet BCs via .constrain()
    # ------------------------------------------------------------------
    pressure = CellVariable(mesh=fipy_mesh, name='pressure', value=0.0)

    if np.any(spec_pressure_mask):
        bc_face_mask, bc_face_vals = _build_face_pressure_bc(
            fipy_mesh, spec_pressure_mask, specified_pressure
        )
        if bc_face_mask.any():
            bc_fv = FaceVariable(mesh=fipy_mesh, value=0.0)
            bc_fv.setValue(bc_face_vals, where=bc_face_mask)
            pressure.constrain(bc_fv, where=bc_face_mask)

    # ------------------------------------------------------------------
    # Assemble and solve:  ∇·(ρ·k/μ · ∇P) = ∇·(ρ²·k/μ · g · ẑ) + Q
    # ------------------------------------------------------------------
    eq = (DiffusionTerm(coeff=diffusion_fv)
          == ConvectionTerm(coeff=gravity_fv) + recharge_source)

    _lus = LinearLUSolver(tolerance=1e-12, iterations=1000)
    eq.solve(var=pressure, solver=_lus)

    return _get_value(pressure)


def solve_transient_pressure_fipy(
    fipy_mesh, backend, k_tensor, viscosity, density, g,
    recharge_flux, recharge_mask,
    spec_pressure_mask, specified_pressure,
    Parameters, dt, pressure_old, cell_centers, porosity=None, gamma=None,
    concentration_old=None, dC_dt_field=None
):
    """
    Solve the transient variable-density groundwater flow equation for total
    pressure P (Pa) using FiPy's implicit time discretisation.

    Governing equation (Bear 1972 / SEAWAT):

        ρ · S_s · ∂P/∂t = ∇·(ρ · k/μ · ∇P) − ∇·(ρ² · k/μ · g · ẑ)

    Rearranged into FiPy operator form:

        TransientTerm(ρ·S_s) == DiffusionTerm(ρ·k/μ) + ConvectionTerm((0, ρ²·k/μ·g)) + Q

    FiPy's TransientTerm applies a backward-Euler (implicit) time step
    internally when eq.solve(dt=dt) is called.

    Parameters
    ----------
    fipy_mesh : fipy.Mesh
    backend : FiPyBackend
    k_tensor : tuple
        Permeability tensor ((kxx, kxy), (kyx, kyy))
    viscosity : float
        Fluid viscosity (Pa·s)
    density : np.ndarray
        Current fluid density field (kg/m³), shape (n_cells,)
    g : float
        Gravitational acceleration (m/s²)
    recharge_flux : float
        Recharge source term
    recharge_mask : np.ndarray
        Boolean mask for recharge cells
    spec_pressure_mask : np.ndarray
        Boolean mask for specified-pressure Dirichlet BC cells
    specified_pressure : np.ndarray
        Total pressure values (Pa) at BC cells
    Parameters : object
        Model parameters (must have .specific_storage)
    dt : float
        Time step (s)
    pressure_old : np.ndarray
        Total pressure field from the previous time step (Pa)
    cell_centers : np.ndarray
        Unused; kept for backward-compatible call signature.
    porosity : float, optional
        Unused; kept for API compatibility.
    gamma : float, optional
        Unused; kept for API compatibility.
    concentration_old : np.ndarray, optional
        Unused; kept for API compatibility.
    dC_dt_field : np.ndarray, optional
        Unused; kept for API compatibility.

    Returns
    -------
    np.ndarray
        Total pressure solution P (Pa), shape (n_cells,)
    """
    # ------------------------------------------------------------------
    # Permeability / mobility
    # ------------------------------------------------------------------
    kxx = k_tensor[0][0]
    kyy = k_tensor[1][1]
    k_eff = np.sqrt(kxx * kyy)
    mobility = k_eff / viscosity          # scalar k/μ
    S_s = Parameters.specific_storage

    # ------------------------------------------------------------------
    # Face-interpolate density  ρ_face  (arithmetic average, numpy)
    # ------------------------------------------------------------------
    _fci   = np.array(fipy_mesh.faceCellIDs)   # (2, n_faces)
    _own   = _fci[0]
    _nb    = _fci[1]
    _valid = _nb >= 0
    rho_face = density[_own].copy()
    rho_face[_valid] = 0.5 * (density[_own[_valid]] + density[_nb[_valid]])

    # ------------------------------------------------------------------
    # DiffusionTerm coefficient:  ρ_face · k/μ  (FaceVariable)
    # ------------------------------------------------------------------
    diffusion_fv = FaceVariable(mesh=fipy_mesh, value=rho_face * mobility)

    # ------------------------------------------------------------------
    # Gravity source via ConvectionTerm:  ∇·(ρ²·k/μ · g · ẑ)
    # ------------------------------------------------------------------
    n_faces = fipy_mesh.numberOfFaces
    gravity_vals = np.zeros((2, n_faces))
    gravity_vals[1] = rho_face**2 * mobility * g
    gravity_fv = FaceVariable(mesh=fipy_mesh, rank=1, value=gravity_vals)

    # ------------------------------------------------------------------
    # TransientTerm coefficient:  ρ_cell · S_s  (CellVariable)
    # ------------------------------------------------------------------
    rho_cell_var = CellVariable(mesh=fipy_mesh, value=density * S_s)

    # ------------------------------------------------------------------
    # Recharge source term
    # ------------------------------------------------------------------
    recharge_array = np.zeros(fipy_mesh.numberOfCells)
    if np.any(recharge_mask):
        recharge_array[recharge_mask] = recharge_flux
    recharge_source = CellVariable(mesh=fipy_mesh, value=recharge_array)

    # ------------------------------------------------------------------
    # Pressure CellVariable initialised from previous timestep + Dirichlet BCs
    # ------------------------------------------------------------------
    pressure = CellVariable(mesh=fipy_mesh, name='pressure',
                            value=pressure_old.copy())

    if np.any(spec_pressure_mask):
        bc_face_mask, bc_face_vals = _build_face_pressure_bc(
            fipy_mesh, spec_pressure_mask, specified_pressure
        )
        if bc_face_mask.any():
            bc_fv = FaceVariable(mesh=fipy_mesh, value=0.0)
            bc_fv.setValue(bc_face_vals, where=bc_face_mask)
            pressure.constrain(bc_fv, where=bc_face_mask)

    # ------------------------------------------------------------------
    # Assemble and solve:
    #   ρ·S_s · ∂P/∂t = ∇·(ρ·k/μ · ∇P) + ∇·(ρ²·k/μ · g · ẑ) + Q
    # ------------------------------------------------------------------
    eq = (TransientTerm(coeff=rho_cell_var)
          == DiffusionTerm(coeff=diffusion_fv)
          + ConvectionTerm(coeff=gravity_fv)
          + recharge_source)

    _lus = LinearLUSolver(tolerance=1e-12, iterations=1000)
    eq.solve(var=pressure, dt=dt, solver=_lus)

    return _get_value(pressure)


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
    Parameters, cell_centers, use_tensor_dispersion=True, convection_scheme='exponential',
    darcy_velocity_face=None
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
    # Darcy flux q (FaceVariable). Use provided face values when supplied to
    # keep coupling consistent and avoid recomputation.
    if darcy_velocity_face is None:
        q = calculate_darcy_flux_fipy(
            fipy_mesh, pressure, density, k_tensor, viscosity, g, cell_centers, Parameters.rho_f_0
        )
    else:
        q = darcy_velocity_face
    
    # The transport PDE is:  phi * dC/dt + div(q * C) = div(D * grad(C))
    # Convection term uses Darcy flux q directly (NOT pore velocity q/phi).
    # The pore velocity (q/phi) is only used for mechanical dispersion calculation.
    velocity = q  # Darcy velocity for convection term

    # ------------------------------------------------------------------
    # concentration_bnd_inflow_only: restrict the concentration BC to
    # cells where flow is entering the domain (mirrors gwflow_lib.py:635-699)
    # ------------------------------------------------------------------
    active_conc_mask = np.array(spec_conc_mask, dtype=bool)

    conc_bnd_inflow_only = getattr(Parameters, 'concentration_bnd_inflow_only', False)
    conc_bnd_inflow_dir = getattr(Parameters, 'concentration_bnd_inflow_direction', 'left')

    if conc_bnd_inflow_only and np.any(active_conc_mask):
        # q is a FaceVariable (2, nFaces) — convert to cell-centred values for
        # the sign check.  We use the mesh face→cell connectivity to average
        # face values back to cell centres.  This is only used to determine
        # inflow vs. outflow direction, so a simple arithmetic average is fine.
        try:
            # FaceVariable global shape is (2, nFaces) across all MPI ranks.
            # fipy_mesh.cellFaceIDs has shape (nFacesPerCell, nCells).
            q_face_vals = _get_value(q)  # (2, nFaces)
            face_ids = np.array(fipy_mesh.cellFaceIDs)  # (nFacesPerCell, nCells)
            # Average the face values belonging to each cell
            qx_cell = np.mean(q_face_vals[0][face_ids], axis=0)
            qy_cell = np.mean(q_face_vals[1][face_ids], axis=0)
        except Exception:
            # Fallback: recompute a quick cell-centred Darcy velocity from
            # the fresh-water head gradient at cell centres.
            try:
                y_cc = cell_centers[:, 1]
                h_f_cc = pressure / (Parameters.rho_f_0 * g) + y_cc
                # Finite-difference gradient using cell-centre neighbours
                # (rough estimate; only the sign matters here)
                keff = np.sqrt(k_tensor[0][0] * k_tensor[1][1])
                dh_dx = np.gradient(h_f_cc.reshape(-1))  # 1-D approx
                qx_cell = -density * (keff / viscosity) * dh_dx
                qy_cell = np.zeros_like(qx_cell)
            except Exception:
                qx_cell = np.zeros(len(active_conc_mask))
                qy_cell = np.zeros(len(active_conc_mask))

        # Determine inflow sign convention (same as gwflow_lib.py)
        if conc_bnd_inflow_dir == 'left':
            inflow_flag = qx_cell > 0   # flow moving rightward = entering from left wall
        elif conc_bnd_inflow_dir == 'right':
            inflow_flag = qx_cell < 0   # flow moving leftward = entering from right wall
        elif conc_bnd_inflow_dir == 'up':
            inflow_flag = qy_cell < 0   # flow moving downward = entering from top wall
        elif conc_bnd_inflow_dir == 'down':
            inflow_flag = qy_cell > 0   # flow moving upward = entering from bottom wall
        else:
            inflow_flag = np.ones(len(active_conc_mask), dtype=bool)

        candidate_mask = active_conc_mask & inflow_flag

        if np.any(candidate_mask):
            active_conc_mask = candidate_mask
        else:
            # Fallback: all BC cells show outflow — pin only the cell at
            # minimum x (mirrors escript fallback in gwflow_lib.py:669-677)
            print('Warning: all concentration BC cells show outflow; '
                  'pinning minimum-x cell as fallback')
            bc_x = cell_centers[active_conc_mask, 0]
            min_x_val = bc_x.min()
            tol_x = np.abs(bc_x).mean() * 1e-6 + 1e-12
            fallback = active_conc_mask & (
                np.abs(cell_centers[:, 0] - min_x_val) <= tol_x
            )
            if np.any(fallback):
                active_conc_mask = fallback
            # else: keep active_conc_mask unchanged as last resort

    # Create concentration variable with initial values (no BC pinning yet)
    initial_conc = np.array(concentration_old, dtype=float)

    concentration = CellVariable(mesh=fipy_mesh, name='concentration',
                                 value=initial_conc)

    # ------------------------------------------------------------------
    # Dirichlet BCs via FiPy's native .constrain() API
    # (replaces the broken penalty method)
    # Find exterior faces adjacent to BC cells and pin them to the
    # specified concentration values.
    # ------------------------------------------------------------------
    if np.any(active_conc_mask):
        face_cell_ids_t = np.array(fipy_mesh.faceCellIDs)  # (2, nFaces)
        owners_t    = face_cell_ids_t[0]
        neighbors_t = face_cell_ids_t[1]
        try:
            ext_face_arr_t = np.array(fipy_mesh.exteriorFaces.value, dtype=bool)
        except Exception:
            ext_face_arr_t = np.zeros(face_cell_ids_t.shape[1], dtype=bool)

        nb_safe_t  = np.where(neighbors_t >= 0, neighbors_t, 0)
        nb_valid_t = neighbors_t >= 0

        owner_is_bc_t = active_conc_mask[owners_t]
        nb_is_bc_t    = nb_valid_t & active_conc_mask[nb_safe_t]
        bc_face_mask_t = ext_face_arr_t & (owner_is_bc_t | nb_is_bc_t)

        if bc_face_mask_t.any():
            bc_face_vals_t = np.zeros(face_cell_ids_t.shape[1])
            use_owner_t = bc_face_mask_t & owner_is_bc_t
            bc_face_vals_t[use_owner_t] = specified_concentration[owners_t[use_owner_t]]
            use_nb_t = bc_face_mask_t & (~owner_is_bc_t) & nb_is_bc_t
            bc_face_vals_t[use_nb_t] = specified_concentration[nb_safe_t[use_nb_t]]

            from fipy import FaceVariable as _FVC
            bc_cv = _FVC(mesh=fipy_mesh, value=0.0)
            bc_cv.setValue(bc_face_vals_t, where=bc_face_mask_t)
            concentration.constrain(bc_cv, where=bc_face_mask_t)
    
    # Calculate dispersion based on configuration choice
    if use_tensor_dispersion:
        # Full anisotropic tensor implementation with correct Bear (1972) formula.
        # Pass q (Darcy flux); the function internally divides by phi for pore velocity.
        dispersion_tensor = calculate_dispersion_coefficients_fipy(
            fipy_mesh, q, porosity, diffusivity, l_disp, t_disp
        )
        # Use tensor diffusion term
        diffusion_term = DiffusionTermCorrection(coeff=dispersion_tensor)
    else:
        # Scalar fallback: D_eff = phi*Dm + alpha_L * |v_mean|
        # Estimate mean pore velocity magnitude from face values
        qx_vals = _get_value(q[0])
        qy_vals = _get_value(q[1])
        v_mean = np.mean(np.sqrt(qx_vals**2 + qy_vals**2)) / porosity
        diffusion_coeff = porosity * diffusivity + l_disp * v_mean
        diffusion_term = DiffusionTerm(coeff=diffusion_coeff)

    # Build and solve the transport equation (no penalty terms needed)
    # velocity = q (Darcy), so the convection term correctly implements div(q*C)
    convection_term = get_convection_term(velocity, scheme=convection_scheme)
    eq = (TransientTerm(coeff=porosity)
          + convection_term
          == diffusion_term)
    
    # Solve for one timestep using LU solver (same reason as pressure solver:
    # penalty=1e10 is ill-conditioned with default iterative solver).
    _lus = LinearLUSolver(tolerance=1e-12, iterations=1000)
    eq.solve(var=concentration, dt=dt, solver=_lus)

    # Clip concentration to physical bounds [0, C_seawater] to suppress any
    # spurious numerical over/undershoot before the value is returned.
    C_max = getattr(Parameters, 'seawater_concentration', None)
    if C_max is None:
        # Derive from specified_concentration if seawater_concentration not set
        spec_conc = getattr(Parameters, 'specified_concentration', [0.03624])
        C_max = float(np.max(spec_conc))
    c_arr = _get_value(concentration)
    c_arr = np.clip(c_arr, 0.0, C_max)
    # Flush subnormal (denormal) values to zero to prevent gradual underflow
    # accumulation that shows up as non-zero C_min (e.g. 9.7e-320) in the log.
    c_arr[np.abs(c_arr) < 1e-300] = 0.0
    return c_arr


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
    # Pass bc=None so the initial concentration field is pure fresh water
    # everywhere (tank-full-of-freshwater start).  Salt is introduced solely
    # through the transient boundary conditions during the time loop.
    print("Setting up initial conditions")
    pressure, concentration, density = setup_fipy_initial_conditions(
        mesh, cell_centers, bc['z_surface'], Parameters, bc=None
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
            Parameters, cell_centers=cell_centers
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
        
        # Tolerance for boundary matching.
        # mean(diff(unique_x)) is O(L/nCells) ~11 µm on a 50k-cell unstructured mesh —
        # far smaller than the cell diameter (~2.5 mm).  Use cellsize instead.
        _cellsize_tol = getattr(Parameters, 'cellsize',
                                getattr(Parameters, 'cellsize_x', None))
        if _cellsize_tol is None:
            _cellsize_tol = 0.01  # fallback
        tol = _cellsize_tol * 0.6
        
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
    # Pre-compute mesh geometry needed for CFL checking
    _min_dx = getattr(Parameters, 'cellsize', getattr(Parameters, 'cellsize_x', None))
    if _min_dx is None:
        _min_dx = np.sqrt(np.array(fipy_mesh.cellVolumes).min())  # fallback
    _max_allowed_CFL = getattr(Parameters, 'max_allowed_CFL_number', 1.0)
    if _max_allowed_CFL is None:
        _max_allowed_CFL = 1.0
    _fn_y = np.array(fipy_mesh.faceNormals)[1]  # ny at each face
    # CFL check should use only INTERIOR faces — exterior boundary faces have
    # artificially large normal velocities (the Dirichlet pressure BC forces
    # a large flux through the face) that would over-constrain the timestep.
    try:
        _ext_face_mask = np.array(fipy_mesh.exteriorFaces.value, dtype=bool)
        _interior_face_mask = ~_ext_face_mask
    except Exception:
        _interior_face_mask = np.ones(fipy_mesh.numberOfFaces, dtype=bool)

    while runtime < total_time and runtime < max_runtime and timestep < max_timesteps:
        
        # ------------------------------------------------------------------
        # CFL-based timestep limiter
        # Before solving, estimate the maximum normal face velocity and clamp
        # dt so that CFL = |qn|_max * dt / dx <= max_allowed_CFL_number.
        # This prevents transport blow-up when large velocities develop.
        # ------------------------------------------------------------------
        density_for_cfl = calculate_fluid_density(
            concentration, Parameters.gamma, Parameters.rho_f_0
        )
        v_cfl = calculate_darcy_flux_fipy(
            fipy_mesh, pressure, density_for_cfl, k_tensor,
            Parameters.viscosity, Parameters.g, cell_centers, Parameters.rho_f_0
        )
        _vx_f = _get_value(v_cfl[0]); _vy_f = _get_value(v_cfl[1])
        _fn_x = np.array(fipy_mesh.faceNormals)[0]
        _qn_all = np.abs(_vx_f * _fn_x + _vy_f * _fn_y) / Parameters.porosity
        # Use only interior faces to avoid artificially constraining dt
        # due to large fluxes at Dirichlet pressure BC faces.
        _qn = _qn_all[_interior_face_mask]
        # Use the 99th percentile of interior face speeds rather than the
        # maximum.  A handful of faces adjacent to Dirichlet pressure BC cells
        # carry artificially large normal fluxes that would otherwise shrink
        # dt to ~0.01 s and require ~15 000 steps instead of ~170.
        _vmax = (np.percentile(_qn, 99) if len(_qn) > 0
                 else np.percentile(_qn_all, 99))
        if _vmax > 0:
            dt_cfl = _max_allowed_CFL * _min_dx / _vmax
            if dt > dt_cfl:
                dt = dt_cfl
        
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

        pressure_changes = []
        concentration_changes = []

        #print(f"    Starting coupled iterations (max={max_coupled_iterations})...")
        for coupled_iter in range(max_coupled_iterations):
            
            # [1] Update fluid density from current concentration
            density = calculate_fluid_density(
                concentration, Parameters.gamma, Parameters.rho_f_0
            )
            
            # [2] Calculate Darcy flux from current pressure
            darcy_flux_face = calculate_darcy_flux_fipy(
                fipy_mesh, pressure, density, k_tensor, Parameters.viscosity, Parameters.g, cell_centers, Parameters.rho_f_0
            )
            
            # [3] Solve solute transport with updated velocity
            if ModelOptions.solute_transport:
                    concentration_before_iter = concentration.copy()   # save for convergence check
                    concentration_new = solve_solute_transport_fipy(
                        fipy_mesh, backend, concentration, dt,
                        pressure, density, k_tensor, Parameters.viscosity,
                        Parameters.g, Parameters.porosity,
                        Parameters.diffusivity, Parameters.l_disp,
                        Parameters.l_disp * Parameters.disp_ratio,
                        bc['spec_conc_mask'], bc['specified_concentration'],
                        Parameters, cell_centers, use_tensor_dispersion=True, convection_scheme=convection_scheme,
                        darcy_velocity_face=darcy_flux_face
                    )
                    concentration = concentration_new
            
            # [3.5] Calculate concentration change rate (dC/dt) for coupling
            # This is used in the pressure equation to account for density-driven flow
            dC_dt_field = (concentration - concentration_start_ts) / dt
            
            # [4] Solve transient pressure with updated density and coupling
            # (This is where coupling happens - density changes affect pressure)
            # Pass the current best-estimate pressure (updated each Picard iter)
            # rather than the fixed pressure_start_ts so density-driven changes
            # can propagate within the coupled loop.
            pressure_new = solve_transient_pressure_fipy(
                fipy_mesh, backend, k_tensor, Parameters.viscosity,
                density, Parameters.g,
                Parameters.recharge_flux, bc['recharge_mask'],
                bc['spec_pressure_mask'], bc['specified_pressure'],
                Parameters, dt, pressure, cell_centers,
                porosity=Parameters.porosity,
                gamma=Parameters.gamma,
                concentration_old=concentration_start_ts,
                dC_dt_field=dC_dt_field
            )
            
            # [5] Check convergence
            pressure_change = np.max(np.abs(pressure_new - pressure)) if len(pressure_new) > 0 else 0.0
            # Use concentration_before_iter (saved before overwriting) so the diff is non-trivial
            conc_change = np.max(np.abs(concentration_new - concentration_before_iter)) if (ModelOptions.solute_transport and len(concentration_new) > 0) else 0.0
            
            pressure_changes.append(pressure_change)
            concentration_changes.append(conc_change)

            pressure = pressure_new
            
            if (pressure_change < pressure_tol and conc_change < concentration_tol):
                coupled_converged = True
                if coupled_iter > 0:
                    #print(f"  Coupled iteration {coupled_iter+1}: CONVERGED (P_Δ={pressure_change:.2e}, C_Δ={conc_change:.2e})")
                    pass
                break
            else:
                if coupled_iter < 3 or (coupled_iter + 1) % 5 == 0:  # Print first 3 and every 5th
                    #print(f"  Coupled iteration {coupled_iter+1}/{max_coupled_iterations}: P_Δ={pressure_change:.2e}, C_Δ={conc_change:.2e}")
                    pass
        
        if not coupled_converged and coupled_iter == max_coupled_iterations - 1:
            print(f"  WARNING: Coupled iterations did not converge after {max_coupled_iterations} iterations, P_Δ={pressure_change:.2e}, C_Δ={conc_change:.2e}")
            print(f" max P changes ", pressure_changes)
            print(f" max C change ", concentration_changes)
            
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
        if output_interval is None or (runtime - last_output_time >= output_interval):
            timestr = get_timestr(runtime)
            _qx = _get_value(darcy_flux_face[0])
            _qy = _get_value(darcy_flux_face[1])
            print(f"Step {timestep}: t = {timestr}, {runtime/total_time*100:.2f}% complete, "
                  f"P: [{pressure.min():.1e}, {pressure.mean():.1e}, {pressure.max():.1e}] Pa, "
                  f"C: [{concentration.min():.4e}, {concentration.mean():.4e}, {concentration.max():.4e}] kg/kg, "
                  f"qx: [{_qx.min():.3e}, {_qx.mean():.3e}, {_qx.max():.3e}] m/s, "
                  f"qy: [{_qy.min():.3e}, {_qy.mean():.3e}, {_qy.max():.3e}] m/s")
            
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
    
    q_face = calculate_darcy_flux_fipy(fipy_mesh, pressure, density, k_tensor, 
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
             face_vals = _get_value(q_face[i])
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
