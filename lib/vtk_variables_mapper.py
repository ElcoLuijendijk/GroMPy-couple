"""
VTK variables mapper for transforming FiPy model output to VTK format.

This module extracts and transforms simulation results from the FiPy backend's
25-element output tuple into a dictionary of 18 cell-centered variables suitable
for VTK visualization, plus metadata attributes.

The module handles:
- Extraction of data from the standard output tuple
- Precision conversion to float32 (for VTK compatibility)
- Creation of boundary condition masks
- Collection of model metadata for embedding in VTK files
- Data validation and error handling

Usage:
    from lib.vtk_variables_mapper import prepare_vtk_variables, prepare_vtk_metadata
    
    # After running FiPy model and getting 25-element tuple
    variables = prepare_vtk_variables(
        model_results=results_tuple,
        model_parameters=Parameters,
        mesh=fipy_mesh,
        cell_centers=cell_centers
    )
    
    metadata = prepare_vtk_metadata(
        model_parameters=Parameters,
        model_options=ModelOptions,
        boundary_flux_stats=results_tuple[24]
    )
"""

import numpy as np
import logging

logger = logging.getLogger(__name__)


def prepare_vtk_variables(
    model_results,
    model_parameters,
    mesh,
    cell_centers,
    boundary_info=None
):
    """
    Extract and prepare VTK variables from FiPy model output tuple.
    
    Maps the 25-element FiPy output tuple to 18 cell-centered variables
    plus metadata, matching the format of esys.weipa.saveVTK output.
    
    Parameters
    ----------
    model_results : tuple
        25-element tuple from run_coupled_flow_model_fipy():
        [0] mesh
        [1] surface
        [2] sea_surface
        [3] k_vector (permeability tensor)
        [4] P (pressure)
        [5] Conc (concentration)
        [6] rho_f (fluid density)
        [7] viscosity
        [8] h (hydraulic head)
        [9] q (velocity/flux vector [qx, qy])
        [10] q_abs (flux magnitude)
        [11] nodal_flux
        [12] Pdiff (pressure difference)
        [13] Cdiff (concentration difference)
        [14] pressure_differences_max
        [15] concentration_differences_max
        [16] pressure_differences_mean
        [17] concentration_differences_mean
        [18] dts
        [19] runtimes
        [20] nsteps
        [21] output_step
        [22] boundary_conditions (8-element list)
        [23] boundary_fluxes
        [24] boundary_flux_stats (dict)
        [25] reached_steady_state
    
    model_parameters : object
        Model parameters object with k, porosity, diffusivity, etc.
    
    mesh : FiPy mesh object
        Mesh used in simulation
    
    cell_centers : np.ndarray
        Cell center coordinates, shape (n_cells, 2)
    
    boundary_info : dict, optional
        Additional boundary information
    
    Returns
    -------
    dict
        Dictionary with 18 variables suitable for VTK:
        - pressure, concentration, hydraulic_head
        - velocity_x, velocity_y, flux_magnitude
        - permeability_xx, permeability_yy
        - nodal_flux
        - surface_mask, sea_surface_mask
        - specified_pressure_bnd, active_seepage_bnd, recharge_bnd
        - active_concentration_bnd
        
        All values converted to float32 for VTK compatibility.
    
    Examples
    --------
    >>> variables = prepare_vtk_variables(model_results, params, mesh, cc)
    >>> print(f"Extracted {len(variables)} variables")
    """
    
    variables = {}
    
    try:
        # Extract components from 25-element tuple
        n_cells = len(cell_centers)
        
        mesh_out = model_results[0]
        surface = model_results[1]
        sea_surface = model_results[2]
        k_tensor = model_results[3]
        P = model_results[4]
        Conc = model_results[5]
        h = model_results[8]
        q = model_results[9]
        q_abs = model_results[10]
        nodal_flux = model_results[11]
        boundary_conditions = model_results[22]
        
        # Core variables (8 variables)
        # ---------------------------------
        
        # Pressure (Pa)
        variables['pressure'] = _to_float32(P)
        
        # Concentration (ppt or g/kg)
        variables['concentration'] = _to_float32(Conc)
        
        # Hydraulic head (m)
        variables['hydraulic_head'] = _to_float32(h)
        
        # Velocity components (m/s)
        if isinstance(q, np.ndarray) and q.ndim == 2 and q.shape[0] == 2:
            variables['velocity_x'] = _to_float32(q[0])
            variables['velocity_y'] = _to_float32(q[1])
        else:
            logger.warning("Velocity data malformed, creating zero arrays")
            variables['velocity_x'] = np.zeros(n_cells, dtype=np.float32)
            variables['velocity_y'] = np.zeros(n_cells, dtype=np.float32)
        
        # Flux magnitude (m/s)
        variables['flux_magnitude'] = _to_float32(q_abs)
        
        # Permeability (m²) - diagonal only (kxx, kyy)
        if isinstance(k_tensor, (tuple, list)) and len(k_tensor) >= 2:
            if isinstance(k_tensor[0], (tuple, list)) and len(k_tensor[0]) >= 1:
                variables['permeability_xx'] = np.full(n_cells, float(k_tensor[0][0]), dtype=np.float32)
            else:
                variables['permeability_xx'] = np.full(n_cells, float(k_tensor[0]), dtype=np.float32)
            
            if isinstance(k_tensor[1], (tuple, list)) and len(k_tensor[1]) >= 2:
                variables['permeability_yy'] = np.full(n_cells, float(k_tensor[1][1]), dtype=np.float32)
            else:
                variables['permeability_yy'] = np.full(n_cells, float(k_tensor[1]), dtype=np.float32)
        else:
            logger.warning("Permeability tensor malformed")
            variables['permeability_xx'] = np.ones(n_cells, dtype=np.float32) * 1e-12
            variables['permeability_yy'] = np.ones(n_cells, dtype=np.float32) * 1e-12
        
        # Boundary & mask variables (5 variables)
        # -----------------------------------------
        
        # Nodal flux (m³/s) - surface flux
        variables['nodal_flux'] = _to_float32(nodal_flux)
        
        # Surface mask (0/1)
        variables['surface_mask'] = surface.astype(np.float32) if isinstance(surface, np.ndarray) else np.zeros(n_cells, dtype=np.float32)
        
        # Sea surface mask (0/1)
        if sea_surface is not None:
            variables['sea_surface_mask'] = sea_surface.astype(np.float32) if isinstance(sea_surface, np.ndarray) else np.zeros(n_cells, dtype=np.float32)
        else:
            variables['sea_surface_mask'] = np.zeros(n_cells, dtype=np.float32)
        
        # Boundary condition masks (derived from boundary_conditions list)
        # Format: [spec_pressure_mask, specified_pressure, spec_conc_mask,
        #          active_seepage_mask, specified_concentration, rho_f_0,
        #          recharge_mask, active_seepage_mask]
        
        # Specified pressure boundary condition mask
        if len(boundary_conditions) > 0:
            spec_p = boundary_conditions[0]
            variables['specified_pressure_bnd'] = spec_p.astype(np.float32) if isinstance(spec_p, np.ndarray) else np.zeros(n_cells, dtype=np.float32)
        else:
            variables['specified_pressure_bnd'] = np.zeros(n_cells, dtype=np.float32)
        
        # Active seepage boundary condition mask
        if len(boundary_conditions) > 3:
            active_seep = boundary_conditions[3]
            variables['active_seepage_bnd'] = active_seep.astype(np.float32) if isinstance(active_seep, np.ndarray) else np.zeros(n_cells, dtype=np.float32)
        else:
            variables['active_seepage_bnd'] = np.zeros(n_cells, dtype=np.float32)
        
        # Recharge boundary condition mask
        if len(boundary_conditions) > 6:
            recharge = boundary_conditions[6]
            variables['recharge_bnd'] = recharge.astype(np.float32) if isinstance(recharge, np.ndarray) else np.zeros(n_cells, dtype=np.float32)
        else:
            variables['recharge_bnd'] = np.zeros(n_cells, dtype=np.float32)
        
        # Active concentration boundary condition mask
        # Derived from spec_conc_mask if available
        if len(boundary_conditions) > 2:
            spec_c = boundary_conditions[2]
            variables['active_concentration_bnd'] = spec_c.astype(np.float32) if isinstance(spec_c, np.ndarray) else np.zeros(n_cells, dtype=np.float32)
        else:
            variables['active_concentration_bnd'] = np.zeros(n_cells, dtype=np.float32)
        
        logger.info(f"Successfully extracted {len(variables)} VTK variables")
        return variables
    
    except Exception as e:
        logger.error(f"Error preparing VTK variables: {e}")
        raise


def prepare_vtk_metadata(
    model_parameters,
    model_options=None,
    boundary_flux_stats=None
):
    """
    Create metadata dictionary for VTK file attributes.
    
    Metadata includes model parameters, boundary conditions, and flux
    statistics. These are stored as VTK field data attributes.
    
    Parameters
    ----------
    model_parameters : object
        Model parameters object with:
        - k or (kxx, kyy) : permeability
        - porosity : porosity
        - diffusivity : molecular diffusivity
        - l_disp : longitudinal dispersivity
        - t_disp : transverse dispersivity
        - viscosity : fluid viscosity
        - rho_f_0 : initial fluid density
        - g : gravitational acceleration
    
    model_options : object, optional
        Model options object with boundary condition settings:
        - left_pressure_head
        - right_pressure_head
        - initial_salinity
        - seepage_face, seepage_flow_in
        - recharge_flux
    
    boundary_flux_stats : dict, optional
        Dictionary with flux statistics from model run:
        - total_seepage_flux
        - total_submarine_flux_out
        - total_land_flux_in, total_land_flux_out
        - total_submarine_flux_in, total_submarine_flux_out
    
    Returns
    -------
    dict
        Metadata key-value pairs suitable for embedding in VTK file.
        Keys are strings, values are floats or strings.
    
    Examples
    --------
    >>> metadata = prepare_vtk_metadata(params, options, flux_stats)
    >>> print(f"Metadata: {list(metadata.keys())}")
    """
    
    metadata = {}
    
    try:
        # Model parameters
        # ----------------
        if hasattr(model_parameters, 'porosity'):
            metadata['porosity'] = float(model_parameters.porosity)
        
        if hasattr(model_parameters, 'k'):
            metadata['permeability_isotropic'] = float(model_parameters.k)
        
        if hasattr(model_parameters, 'diffusivity'):
            metadata['diffusivity'] = float(model_parameters.diffusivity)
        
        if hasattr(model_parameters, 'l_disp'):
            metadata['dispersivity_longitudinal'] = float(model_parameters.l_disp)
        
        if hasattr(model_parameters, 't_disp'):
            metadata['dispersivity_transverse'] = float(model_parameters.t_disp)
        
        if hasattr(model_parameters, 'viscosity'):
            metadata['viscosity'] = float(model_parameters.viscosity)
        
        if hasattr(model_parameters, 'rho_f_0'):
            metadata['fluid_density_initial'] = float(model_parameters.rho_f_0)
        
        if hasattr(model_parameters, 'g'):
            metadata['gravity'] = float(model_parameters.g)
        
        # Boundary conditions
        # -------------------
        if model_options is not None:
            if hasattr(model_options, 'left_pressure_head'):
                metadata['left_pressure_head'] = float(model_options.left_pressure_head)
            
            if hasattr(model_options, 'right_pressure_head'):
                metadata['right_pressure_head'] = float(model_options.right_pressure_head)
            
            if hasattr(model_options, 'initial_salinity'):
                metadata['initial_salinity'] = float(model_options.initial_salinity)
            
            if hasattr(model_options, 'seepage_face'):
                metadata['seepage_active'] = int(model_options.seepage_face)
            
            if hasattr(model_options, 'recharge_flux'):
                metadata['recharge_flux'] = float(model_options.recharge_flux)
        
        # Flux statistics (stored as metadata only, not cell variables)
        # -----------------------------------------------------------
        if boundary_flux_stats is not None and isinstance(boundary_flux_stats, dict):
            flux_keys = [
                'total_seepage_flux',
                'total_submarine_flux_out',
                'total_submarine_flux_in',
                'total_land_flux_in',
                'total_land_flux_out',
                'min_seepage_flux',
                'max_seepage_flux',
                'min_submarine_flux',
                'max_submarine_flux',
            ]
            
            for key in flux_keys:
                if key in boundary_flux_stats:
                    val = boundary_flux_stats[key]
                    if isinstance(val, (int, float)):
                        metadata[f'flux_{key}'] = float(val)
        
        # Simulation metadata
        # -------------------
        metadata['backend'] = 'fipy'
        metadata['vtk_format_version'] = '1.0'
        
        logger.info(f"Created metadata with {len(metadata)} attributes")
        return metadata
    
    except Exception as e:
        logger.warning(f"Error preparing metadata: {e}")
        return {}


def validate_variables(variables, n_cells):
    """
    Validate that all variables have correct shape and reasonable values.
    
    Parameters
    ----------
    variables : dict
        Dictionary of {name: np.ndarray(n_cells)}
    
    n_cells : int
        Expected number of cells
    
    Returns
    -------
    bool
        True if valid, False otherwise
    
    Logs warnings for issues found.
    """
    try:
        if not isinstance(variables, dict):
            logger.error("variables must be dict")
            return False
        
        if len(variables) == 0:
            logger.error("variables dict is empty")
            return False
        
        issues = []
        
        for name, arr in variables.items():
            if not isinstance(arr, np.ndarray):
                issues.append(f"{name}: not ndarray (got {type(arr)})")
                continue
            
            if arr.shape[0] != n_cells:
                issues.append(f"{name}: shape mismatch (got {arr.shape[0]}, expected {n_cells})")
                continue
            
            # Check for NaN in critical variables
            if name in ['pressure', 'concentration', 'hydraulic_head']:
                n_nan = np.isnan(arr).sum()
                if n_nan > 0:
                    issues.append(f"{name}: contains {n_nan} NaN values")
        
        if issues:
            for issue in issues:
                logger.warning(issue)
            return False
        
        logger.info(f"Validation passed for {len(variables)} variables")
        return True
    
    except Exception as e:
        logger.error(f"Validation error: {e}")
        return False


# Helper function
def _to_float32(arr):
    """
    Convert array to float32 with proper error handling.
    
    Parameters
    ----------
    arr : array-like
        Input array
    
    Returns
    -------
    np.ndarray
        Array converted to float32
    """
    try:
        if arr is None:
            return np.array([], dtype=np.float32)
        
        arr = np.asarray(arr, dtype=np.float64)
        
        # Handle infinities
        inf_mask = np.isinf(arr)
        if np.any(inf_mask):
            max_f32 = np.finfo(np.float32).max
            arr[inf_mask & (arr > 0)] = max_f32
            arr[inf_mask & (arr < 0)] = -max_f32
        
        return arr.astype(np.float32)
    
    except Exception as e:
        logger.warning(f"Could not convert to float32: {e}")
        return np.zeros_like(arr, dtype=np.float32)
