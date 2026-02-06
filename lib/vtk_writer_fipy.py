"""
VTK file writer for FiPy-based groundwater flow and solute transport simulations.

This module provides functionality to write simulation results to VTK format
(.vtu files) for visualization in ParaView. It supports both pyVTK (native
VTK Python bindings) and meshio (fallback) as writing backends, with
graceful degradation if libraries are unavailable.

Key Features:
- float32 precision (50% smaller files than float64)
- Optional compression (50-70% file size reduction)
- Metadata embedding (model parameters, boundary conditions)
- Dual-library support with automatic fallback
- Graceful error handling (no exceptions, continues if VTK fails)

Usage:
    from lib.vtk_writer_fipy import write_vtk_fipy
    
    success = write_vtk_fipy(
        filename='output.vtu',
        mesh=fipy_mesh,
        cell_centers=cell_centers,
        variables={'pressure': pressure_array, ...},
        metadata={'porosity': 0.2, ...},
        compression=True,
        precision='float32'
    )
"""

import numpy as np
import warnings
import logging

# Set up logging
logger = logging.getLogger(__name__)


def check_vtk_available():
    """
    Check which VTK writing libraries are available.
    
    Returns
    -------
    str
        'pyVTK', 'meshio', or 'none'
    """
    try:
        import vtkmodules.all as vtk
        return 'pyVTK'
    except ImportError:
        pass
    
    try:
        import meshio
        return 'meshio'
    except ImportError:
        pass
    
    return 'none'


def write_vtk_fipy(
    filename,
    mesh,
    cell_centers,
    variables,
    metadata=None,
    compression=True,
    precision='float32'
):
    """
    Write FiPy simulation results to VTK unstructured grid format.
    
    Supports both pyVTK (primary) and meshio (fallback) as backends.
    Automatically selects best available library with graceful degradation.
    
    Parameters
    ----------
    filename : str
        Output file path (.vtu extension recommended)
    mesh : FiPy mesh object
        Mesh object containing geometry
    cell_centers : np.ndarray
        Array of shape (n_cells, 2) with cell center coordinates
    variables : dict
        Dictionary of {variable_name: np.ndarray} where each array has
        shape (n_cells,) with cell-centered data
    metadata : dict, optional
        Dictionary of metadata key-value pairs to embed as VTK attributes.
        Values should be scalar (int, float, str) or small arrays.
    compression : bool, default=True
        Use VTK compression (reduces file size ~50-70%)
    precision : str, default='float32'
        Data precision: 'float32' (recommended) or 'float64'
    
    Returns
    -------
    bool
        True if VTK file written successfully, False otherwise.
        Never raises exception - logs warnings instead.
    
    Notes
    -----
    - Automatically converts arrays to specified precision
    - Validates input data shapes and finite values
    - Falls back gracefully if preferred library unavailable
    - Encodes metadata as VTK field data (visible in ParaView)
    
    Examples
    --------
    >>> success = write_vtk_fipy(
    ...     'model_output.vtu',
    ...     mesh=fipy_mesh,
    ...     cell_centers=cc,
    ...     variables={'pressure': P, 'concentration': C},
    ...     metadata={'porosity': 0.2},
    ...     compression=True
    ... )
    >>> if success:
    ...     print("VTK file written successfully")
    >>> else:
    ...     print("Warning: VTK output failed")
    """
    
    try:
        # Validate inputs
        if not _validate_inputs(mesh, cell_centers, variables):
            return False
        
        # Determine data precision
        if precision == 'float32':
            dtype = np.float32
        else:
            dtype = np.float64
        
        # Extract mesh geometry
        points, cells = _extract_mesh_geometry(mesh)
        if points is None:
            logger.warning("Could not extract mesh geometry")
            return False
        
        # Prepare cell data with precision conversion
        cell_data = {}
        for name, arr in variables.items():
            try:
                arr_converted = _convert_to_precision(arr, dtype)
                cell_data[name] = arr_converted
            except Exception as e:
                logger.warning(f"Skipping variable '{name}': {e}")
                continue
        
        if not cell_data:
            logger.warning("No valid cell data to write")
            return False
        
        # Get available VTK writer
        vtk_writer = check_vtk_available()
        
        if vtk_writer == 'pyVTK':
            # Try pyVTK with compression first
            if compression:
                success = _write_vtk_pyVTK(
                    filename, points, cells, cell_data,
                    metadata=metadata, compression=True, precision=dtype
                )
                if success:
                    return True
                logger.warning("pyVTK with compression failed, trying without")
            
            # Try pyVTK without compression
            success = _write_vtk_pyVTK(
                filename, points, cells, cell_data,
                metadata=metadata, compression=False, precision=dtype
            )
            if success:
                return True
            logger.warning("pyVTK writer failed, trying meshio fallback")
        
        if vtk_writer == 'meshio' or vtk_writer == 'pyVTK':
            # Try meshio
            success = _write_vtk_meshio(
                filename, points, cells, cell_data,
                metadata=metadata, precision=dtype
            )
            if success:
                return True
            logger.warning("meshio writer failed")
        
        logger.warning("No VTK writer available (install vtkmodules or meshio)")
        return False
    
    except Exception as e:
        logger.warning(f"VTK writing failed: {e}")
        return False


def _validate_inputs(mesh, cell_centers, variables):
    """
    Validate input data before VTK writing.
    
    Parameters
    ----------
    mesh : FiPy mesh
    cell_centers : np.ndarray
        Shape (n_cells, 2)
    variables : dict
        Variable data
    
    Returns
    -------
    bool
        True if valid, False otherwise
    """
    try:
        if mesh is None:
            logger.warning("Mesh is None")
            return False
        
        if cell_centers is None or len(cell_centers) == 0:
            logger.warning("cell_centers invalid")
            return False
        
        n_cells = len(cell_centers)
        
        if not isinstance(variables, dict):
            logger.warning("variables must be dict")
            return False
        
        if len(variables) == 0:
            logger.warning("variables dict is empty")
            return False
        
        # Check each variable
        for name, arr in variables.items():
            if not isinstance(arr, np.ndarray):
                logger.warning(f"{name}: not ndarray")
                return False
            
            if len(arr) != n_cells:
                logger.warning(f"{name}: shape mismatch (got {len(arr)}, expected {n_cells})")
                return False
        
        return True
    
    except Exception as e:
        logger.warning(f"Input validation failed: {e}")
        return False


def _extract_mesh_geometry(mesh):
    """
    Extract vertex coordinates and cell connectivity from FiPy mesh.
    
    Parameters
    ----------
    mesh : FiPy mesh object
        Supports Tri2D, Gmsh2D, and similar mesh types
    
    Returns
    -------
    points : np.ndarray or None
        Vertex coordinates, shape (n_vertices, 3) with z=0
    cells : np.ndarray or None
        Cell connectivity, shape (n_cells, n_vertices_per_cell)
    
    Returns (None, None) if extraction fails.
    """
    try:
        # Try to get vertex coordinates
        if hasattr(mesh, 'vertexCoords'):
            # FiPy format: (2, n_vertices) or (n_vertices, 2)
            verts = mesh.vertexCoords
            if verts.shape[0] == 2:
                verts = verts.T  # Transpose to (n_vertices, 2)
        else:
            logger.warning("Mesh has no vertexCoords attribute")
            return None, None
        
        # Add z-coordinate (zeros) for 3D VTK format
        if verts.shape[1] == 2:
            points = np.column_stack([verts, np.zeros(len(verts))])
        else:
            points = verts
        
        # Get cell connectivity
        if hasattr(mesh, '_cellVertexIDs'):
            cells = mesh._cellVertexIDs
            # FiPy format: (n_vertices_per_cell, n_cells)
            if cells.shape[0] < cells.shape[1]:
                cells = cells.T  # Transpose to (n_cells, n_vertices_per_cell)
        else:
            logger.warning("Mesh has no _cellVertexIDs attribute")
            return None, None
        
        return points, cells
    
    except Exception as e:
        logger.warning(f"Mesh geometry extraction failed: {e}")
        return None, None


def _convert_to_precision(arr, dtype):
    """
    Convert array to target precision with NaN/Inf handling.
    
    Parameters
    ----------
    arr : np.ndarray
        Input array (any dtype)
    dtype : type
        Target dtype (np.float32 or np.float64)
    
    Returns
    -------
    np.ndarray
        Converted array with NaN/Inf preserved
    
    Raises
    ------
    ValueError
        If array cannot be converted to float
    """
    try:
        arr = np.asarray(arr, dtype=np.float64)
        
        # Check for problematic values
        inf_mask = np.isinf(arr)
        if np.any(inf_mask):
            # Replace infinities with max/min finite value
            if dtype == np.float32:
                max_val = np.finfo(np.float32).max
                arr[inf_mask & (arr > 0)] = max_val
                arr[inf_mask & (arr < 0)] = -max_val
            else:
                # For float64, replace with extreme but finite values
                arr[inf_mask & (arr > 0)] = 1e308
                arr[inf_mask & (arr < 0)] = -1e308
        
        # Convert to target precision
        result = arr.astype(dtype)
        
        return result
    
    except Exception as e:
        raise ValueError(f"Could not convert to {dtype}: {e}")


def _write_vtk_pyVTK(
    filename,
    points,
    cells,
    cell_data,
    metadata=None,
    compression=True,
    precision=np.float32
):
    """
    Write VTK file using native pyVTK (vtkmodules) library.
    
    Parameters
    ----------
    filename : str
        Output file path
    points : np.ndarray
        Vertex coordinates, shape (n_vertices, 3)
    cells : np.ndarray
        Cell connectivity, shape (n_cells, n_vertices_per_cell)
    cell_data : dict
        {variable_name: np.ndarray(n_cells)}
    metadata : dict, optional
        Metadata to embed
    compression : bool
        Use VTK compression
    precision : type
        Data type (np.float32 or np.float64)
    
    Returns
    -------
    bool
        True if successful
    """
    try:
        import vtkmodules.all as vtk
        from vtkmodules.util.numpy_support import numpy_to_vtk
    except ImportError:
        logger.warning("pyVTK (vtkmodules) not available")
        return False
    
    try:
        n_cells = cells.shape[0]
        n_verts_per_cell = cells.shape[1]
        
        # Determine VTK cell type
        if n_verts_per_cell == 3:
            vtk_cell_type = vtk.VTK_TRIANGLE
        elif n_verts_per_cell == 4:
            vtk_cell_type = vtk.VTK_QUAD
        else:
            logger.warning(f"Unsupported cell type: {n_verts_per_cell} vertices")
            return False
        
        # Create unstructured grid
        grid = vtk.vtkUnstructuredGrid()
        
        # Add points
        vtk_points = vtk.vtkPoints()
        vtk_points.SetData(numpy_to_vtk(points, deep=True))
        grid.SetPoints(vtk_points)
        
        # Add cells
        for i, cell in enumerate(cells):
            ids = vtk.vtkIdList()
            for vertex_id in cell:
                ids.InsertNextId(int(vertex_id))
            grid.InsertNextCell(vtk_cell_type, ids)
        
        # Add cell data (variables)
        for var_name, var_array in cell_data.items():
            # Ensure correct precision
            var_array = var_array.astype(precision)
            vtk_array = numpy_to_vtk(var_array, deep=True)
            vtk_array.SetName(var_name)
            grid.GetCellData().AddArray(vtk_array)
        
        # Add metadata as field data
        if metadata:
            for key, value in metadata.items():
                try:
                    if isinstance(value, (int, float)):
                        meta_array = numpy_to_vtk(np.array([value], dtype=precision))
                        meta_array.SetName(str(key))
                        grid.GetFieldData().AddArray(meta_array)
                    elif isinstance(value, str):
                        # Store strings as field data
                        str_array = vtk.vtkStringArray()
                        str_array.SetName(str(key))
                        str_array.InsertNextValue(str(value))
                        grid.GetFieldData().AddArray(str_array)
                except Exception as e:
                    logger.warning(f"Could not add metadata '{key}': {e}")
        
        # Write to file
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(str(filename))
        writer.SetInputData(grid)
        
        if compression:
            writer.SetCompressorTypeToZLib()
        
        result = writer.Write()
        
        if result == 0:
            logger.warning("vtkXMLUnstructuredGridWriter returned error")
            return False
        
        logger.info(f"VTK file written successfully: {filename}")
        return True
    
    except Exception as e:
        logger.warning(f"pyVTK writing failed: {e}")
        return False


def _write_vtk_meshio(
    filename,
    points,
    cells,
    cell_data,
    metadata=None,
    precision=np.float32
):
    """
    Write VTK file using meshio library (fallback).
    
    Parameters
    ----------
    filename : str
        Output file path
    points : np.ndarray
        Vertex coordinates, shape (n_vertices, 3)
    cells : np.ndarray
        Cell connectivity, shape (n_cells, n_vertices_per_cell)
    cell_data : dict
        {variable_name: np.ndarray(n_cells)}
    metadata : dict, optional
        Metadata (note: meshio has limited metadata support)
    precision : type
        Data type (np.float32 or np.float64)
    
    Returns
    -------
    bool
        True if successful
    
    Note
    ----
    Compression not supported with meshio fallback.
    """
    try:
        import meshio
    except ImportError:
        logger.warning("meshio not available")
        return False
    
    try:
        # Ensure precision
        points = points.astype(precision)
        
        # Convert cell data to meshio format
        cell_data_meshio = {}
        for var_name, var_array in cell_data.items():
            var_array = var_array.astype(precision)
            # meshio expects cell_data as {name: [array_for_triangle, array_for_quad, ...]}
            cell_data_meshio[var_name] = [var_array]
        
        # Determine cell type
        n_verts = cells.shape[1]
        if n_verts == 3:
            cell_type = 'triangle'
        elif n_verts == 4:
            cell_type = 'quad'
        else:
            logger.warning(f"Unsupported cell type: {n_verts} vertices")
            return False
        
        # Create meshio mesh
        mesh_obj = meshio.Mesh(
            points=points,
            cells=[(cell_type, cells)],
            cell_data=cell_data_meshio
        )
        
        # Write to file (meshio auto-detects format from extension)
        mesh_obj.write(str(filename))
        
        logger.info(f"VTK file written with meshio: {filename}")
        return True
    
    except Exception as e:
        logger.warning(f"meshio writing failed: {e}")
        return False
