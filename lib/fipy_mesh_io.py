"""
Helpers to load a Gmsh mesh into FiPy and construct simple cell masks.

This module provides ``load_fipy_mesh_from_msh``, which reads a Gmsh mesh
file (written by the existing ``mesh_functions.py`` workflow) and returns a
FiPy mesh plus cell-centred masks that approximate the escript masks used
in the rest of the codebase (``surface``, ``sea_surface``, ``seawater``).

Notes
-----
- Requires ``meshio`` and ``fipy`` installed in the environment.
- When a sibling ``.geo`` file exists alongside the ``.msh`` (written by
  ``setup_rectangular_mesh_fipy``), it is passed directly to ``Gmsh2D``.
  ``Gmsh2D`` generates the mesh from the geometry description, which enables
  automatic MPI partitioning under ``mpirun -np N`` (Gmsh -part N).
- For coastal/standard meshes without a sibling ``.geo``, the ``.msh`` file
  is still passed to ``Gmsh2D`` (FiPy reads it directly).
- The returned masks are boolean numpy arrays aligned with FiPy cell order
  (cell-centred).  They are not escript Field objects.

"""

import os
import numpy as np


def load_fipy_mesh_from_msh(mesh_filename, topo_gradient=None, tol=1e-8):
    """
    Read a Gmsh mesh and return a FiPy mesh and cell-centred masks.

    Parameters
    ----------
    mesh_filename : str
        Path to the Gmsh .msh file.
    topo_gradient : float or None
        If provided the function computes z_surface = x * topo_gradient and
        constructs ``surface`` and ``seawater`` masks similar to the escript
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

        - ``'surface'``
        - ``'sea_surface'``
        - ``'seawater'``
        - ``'physical_groups'`` : dict(name -> boolean array)
    field_data : dict
        Mapping physical group name -> integer tag from the .msh (if present).
    extra : dict
        Additional information (e.g. ``'z_surface_cells'``).

    Raises
    ------
    FileNotFoundError
        If the mesh file does not exist.
    RuntimeError
        If triangles are not present and a FiPy triangular mesh cannot be
        constructed.
    """
    if not os.path.exists(mesh_filename):
        raise FileNotFoundError(mesh_filename)

    # ------------------------------------------------------------------
    # Resolve FiPy Gmsh2D loader
    # ------------------------------------------------------------------
    fipy_mesh = None
    try:
        import fipy  # noqa: F401 (side-effect: makes fipy CellVariable etc. available)
    except Exception:
        fipy = None

    gmsh_loader = None
    if fipy is not None:
        try:
            from fipy import Gmsh2D
            gmsh_loader = Gmsh2D
        except Exception:
            try:
                from fipy.meshes.gmsh2D import Gmsh2D
                gmsh_loader = Gmsh2D
            except Exception:
                gmsh_loader = None

    # ------------------------------------------------------------------
    # Always load from the pre-built .msh file via meshio (path B).
    #
    # A sibling .geo file may exist (written by setup_rectangular_mesh_fipy)
    # but we deliberately do NOT use Gmsh2D with an inline .geo string here.
    # Under MPI (mpirun -np N) that path causes Gmsh to re-mesh the domain
    # independently on every rank, partition the result, and write a format-2
    # file — producing gigabytes of redundant work, "Appending zeros" warnings,
    # and a broken global assembly.  Loading the pre-built .msh is faster,
    # deterministic, and works correctly in both serial and MPI runs.
    # ------------------------------------------------------------------
    if True:  # always use path B (.msh via meshio)
        # ---- path B: read .msh via meshio then build FiPy mesh ----
        try:
            import meshio
        except Exception:
            raise RuntimeError("meshio is required to read .msh files")

        msh = meshio.read(mesh_filename, file_format="gmsh")

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
            raise RuntimeError(
                "No triangle cells found in mesh; "
                "FiPy triangular fallback requires triangles."
            )

        points = msh.points
        if points.shape[1] >= 2:
            xy = points[:, :2]
        else:
            raise RuntimeError("Mesh points do not contain x,y coordinates")

        # Compute cell centres from meshio triangle connectivity
        cell_centers = np.mean(xy[tri_cells], axis=1)

        # Attempt to load FiPy mesh from .msh
        if gmsh_loader is not None:
            try:
                fipy_mesh = gmsh_loader(mesh_filename)
            except Exception:
                fipy_mesh = None

        if fipy_mesh is None:
            # Fallback: construct from nodes + connectivity
            try:
                try:
                    from fipy.meshes import Tri2D
                except ImportError:
                    try:
                        from fipy.meshes.triangularMesh import Tri2D
                    except ImportError:
                        from fipy.meshes.tri import Tri2D

                x = xy[:, 0]
                y = xy[:, 1]
                fipy_mesh = Tri2D(x, y, tri_cells)
            except Exception as e:
                raise RuntimeError(
                    f"Could not construct a FiPy triangular mesh: {e}. "
                    "Ensure FiPy supports Tri2D or Gmsh2D in your installation."
                )

        # physical group data (for .msh path only)
        phys_ids = None
        try:
            cell_data_dict = msh.cell_data_dict
            if 'gmsh:physical' in cell_data_dict:
                phys = cell_data_dict['gmsh:physical']
                if isinstance(phys, dict) and tri_block_type in phys:
                    phys_ids = np.asarray(phys[tri_block_type])
                else:
                    for arr in phys.values():
                        if isinstance(arr, (list, tuple, np.ndarray)):
                            phys_ids = np.asarray(arr)
                            break
        except Exception:
            phys_ids = None

    # ------------------------------------------------------------------
    # Build simple masks similar to the escript outputs
    # ------------------------------------------------------------------
    masks = {}
    x_centers = cell_centers[:, 0]
    y_centers = cell_centers[:, 1]

    if topo_gradient is not None:
        z_surface_cells = x_centers * topo_gradient
    else:
        z_surface_cells = np.zeros_like(x_centers)

    # For flat domains (topo_gradient == 0), z_surface = 0 but the domain top
    # is at y = thickness.  Use y_centers.max() as the effective surface level.
    if topo_gradient is None or topo_gradient == 0.0:
        z_surface_cells = np.full_like(x_centers, y_centers.max())

    # surface: cells whose centre y ~= z_surface
    y_range = y_centers.max() - y_centers.min()
    surface_tol = max(tol, y_range * 0.01)  # at least 1% of domain height
    surface_mask = np.abs(y_centers - z_surface_cells) <= surface_tol
    # sea surface: y == 0 and x < 0 (approx)
    sea_surface_mask = (np.abs(y_centers - 0.0) <= tol) & (x_centers < 0)
    # seawater: x < 0 and cell centre below surface
    seawater_mask = (x_centers < 0) & (y_centers <= z_surface_cells + tol)

    masks['surface'] = surface_mask
    masks['sea_surface'] = sea_surface_mask
    masks['seawater'] = seawater_mask

    # physical group masks (only available on .msh path)
    phys_groups = {}
    if phys_ids is not None and field_data:
        for name, tag in field_data.items():
            phys_groups[name] = (phys_ids == tag)
    masks['physical_groups'] = phys_groups

    extra = {
        'z_surface_cells': z_surface_cells,
        'cell_centers': cell_centers,
    }

    return fipy_mesh, cell_centers, masks, field_data, extra
