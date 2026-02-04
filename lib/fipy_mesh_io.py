"""
Helpers to load a Gmsh .msh into FiPy and construct simple cell masks

This module provides `load_fipy_mesh_from_msh`, which reads a Gmsh mesh
file (written by the existing `mesh_functions.py` workflow) and returns a
FiPy mesh plus cell-centered masks that approximate the escript masks used
in the rest of the codebase (``surface``, ``sea_surface``, ``seawater``)

Notes:
- Requires `meshio` and `fipy` installed in the environment.
- If FiPy's native Gmsh loader is available (`fipy.Gmsh2D`) it is used.
  Otherwise a triangular FiPy mesh is constructed from node + triangle
  connectivity using FiPy triangular mesh constructors.
- The returned masks are boolean numpy arrays aligned with FiPy cell order
  (cell-centred). They are not escript Field objects; downstream solver
  code needs to accept cell-centred masks or an adapter must be added.

"""

import os
import numpy as np


def load_fipy_mesh_from_msh(mesh_filename, topo_gradient=None, tol=1e-8):
    """
    Read a Gmsh .msh file and return a FiPy mesh and cell-centered masks.

    Parameters
    ----------
    mesh_filename : str
        Path to the Gmsh .msh file.
    topo_gradient : float or None
        If provided the function computes z_surface = x * topo_gradient and
        constructs `surface` and `seawater` masks similar to the escript
        implementation (based on cell-centre coordinates).
    tol : float
        Tolerance used for comparisons to identify surface cells.

    Returns
    -------
    fipy_mesh : fipy.mesh
        FiPy mesh instance (Gmsh2D or Tri2D depending on availability).
    cell_centers : ndarray (nCells, 2)
        Array with cell centre coordinates (x, y).
    masks : dict
        Dictionary with boolean arrays aligned to cells:
          - 'surface'
          - 'sea_surface'
          - 'seawater'
          - 'physical_groups' : dict(name -> boolean array)
    field_data : dict
        Mapping physical group name -> integer tag from the .msh (if present)
    extra : dict
        Additional information (e.g. 'z_surface_cells')

    Raises
    ------
    FileNotFoundError if the mesh file does not exist.
    RuntimeError if triangles are not present and a FiPy triangular mesh
    cannot be constructed.
    """
    if not os.path.exists(mesh_filename):
        raise FileNotFoundError(mesh_filename)

    # Try to import FiPy Gmsh loader
    fipy_mesh = None
    try:
        import fipy
    except Exception:
        fipy = None

    # Prefer FiPy Gmsh2D if present
    gmsh_loader = None
    if fipy is not None:
        try:
            # Some FiPy versions expose Gmsh2D in the top-level
            from fipy import Gmsh2D
            gmsh_loader = Gmsh2D
        except Exception:
            # try known alternative location
            try:
                from fipy.meshes.gmsh2D import Gmsh2D
                gmsh_loader = Gmsh2D
            except Exception:
                gmsh_loader = None

    # Read with meshio to obtain nodes, triangle cells and physical groups
    try:
        import meshio
    except Exception:
        raise RuntimeError("meshio is required to read .msh files")

    msh = meshio.read(mesh_filename)

    # field_data maps names to tags
    field_data = {}
    if hasattr(msh, 'field_data') and msh.field_data:
        for name, data in msh.field_data.items():
            try:
                field_data[name] = int(data[0])
            except Exception:
                field_data[name] = data

    # find triangle cell block
    tri_cells = None
    tri_block_type = None
    for block in msh.cells:
        if block.type in ('triangle', 'triangle3'):
            tri_cells = block.data
            tri_block_type = block.type
            break

    if tri_cells is None:
        # If no triangles found, try to triangulate quads or raise
        raise RuntimeError("No triangle cells found in mesh; FiPy triangular fallback requires triangles.")

    points = msh.points
    # ensure 2D points
    if points.shape[1] >= 2:
        xy = points[:, :2]
    else:
        raise RuntimeError("Mesh points do not contain x,y coordinates")

    # compute cell centers
    cell_centers = np.mean(xy[tri_cells], axis=1)

    # Attempt to construct FiPy mesh
    if gmsh_loader is not None:
        try:
            # Use Gmsh2D to load file directly into FiPy
            fipy_mesh = gmsh_loader(mesh_filename)
        except Exception:
            fipy_mesh = None

    if fipy_mesh is None:
        # Try to construct triangular mesh programmatically (Tri2D)
        try:
            # Import path can vary across FiPy versions
            try:
                from fipy.meshes import Tri2D
            except ImportError:
                try:
                    from fipy.meshes.triangularMesh import Tri2D
                except ImportError:
                    # older FiPy exposed Tri2D at a different path
                    from fipy.meshes.tri import Tri2D

            # Tri2D expects arrays of x, y and an element connectivity
            x = xy[:, 0]
            y = xy[:, 1]
            # tri_cells is of shape (nCells, 3)
            fipy_mesh = Tri2D(x, y, tri_cells)
        except Exception as e:
            raise RuntimeError(f"Could not construct a FiPy triangular mesh: {e}. Ensure FiPy supports Tri2D or Gmsh2D in your installation.")

    # Build simple masks similar to the escript outputs
    masks = {}
    x_centers = cell_centers[:, 0]
    y_centers = cell_centers[:, 1]

    if topo_gradient is not None:
        z_surface_cells = x_centers * topo_gradient
    else:
        z_surface_cells = np.zeros_like(x_centers)

    # surface: cells whose centre y ~= z_surface
    surface_mask = np.abs(y_centers - z_surface_cells) <= tol
    # sea surface: y == 0 and x < 0 (approx)
    sea_surface_mask = (np.abs(y_centers - 0.0) <= tol) & (x_centers < 0)
    # seawater: x < 0 and cell centre below surface
    seawater_mask = (x_centers < 0) & (y_centers <= z_surface_cells + tol)

    masks['surface'] = surface_mask
    masks['sea_surface'] = sea_surface_mask
    masks['seawater'] = seawater_mask

    # physical group masks: try to extract gmsh:physical from cell_data
    phys_groups = {}
    try:
        # meshio.cell_data_dict gives mapping e.g. {'gmsh:physical': {'triangle': array_of_ids}}
        cell_data_dict = msh.cell_data_dict
        if 'gmsh:physical' in cell_data_dict:
            phys = cell_data_dict['gmsh:physical']
            # phys may be a dict keyed by cell block type
            if isinstance(phys, dict) and tri_block_type in phys:
                phys_ids = np.asarray(phys[tri_block_type])
            else:
                # try to flatten cell_data list
                phys_ids = None
                for key, arr in phys.items():
                    if isinstance(arr, (list, tuple)):
                        # pick first
                        phys_ids = np.asarray(arr)
                        break
                if phys_ids is None:
                    phys_ids = None
        else:
            phys_ids = None
    except Exception:
        phys_ids = None

    if phys_ids is not None and field_data:
        for name, tag in field_data.items():
            phys_groups[name] = (phys_ids == tag)
    else:
        phys_groups = {}

    masks['physical_groups'] = phys_groups

    extra = {
        'z_surface_cells': z_surface_cells,
        'cell_centers': cell_centers
    }

    return fipy_mesh, cell_centers, masks, field_data, extra
