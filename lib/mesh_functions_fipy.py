"""
FiPy-native mesh creation functions for GroMPy.

This module provides mesh creation functions that use FiPy and the Gmsh Python API
directly, without requiring esys-escript. It mirrors the interface of mesh_functions.py
to allow drop-in replacement when using the FiPy backend.

Mesh Types:
- rectangular: Simple rectangular mesh using FiPy's Grid2D (no external dependencies)
- standard: Mesh with bi-linear surface topography, created via Gmsh Python API
- coastal: Complex coastal mesh with variable refinement, created via Gmsh Python API
- coastal_glover1959: Coastal mesh with salt-water interface from Glover (1959) analytical solution

Dependencies:
- fipy: For mesh representation and Grid2D
- gmsh: Python package for complex mesh generation (pip install gmsh)
- meshio: For reading Gmsh .msh files

Usage:
    from lib.mesh_functions_fipy import setup_rectangular_mesh_fipy
    mesh, surface, sea_surface, seawater, z_surface = setup_rectangular_mesh_fipy(Parameters, mesh_filename)
"""

import os
import math
import numpy as np

# Check for FiPy
FIPY_AVAILABLE = False
try:
    import fipy
    FIPY_AVAILABLE = True
except ImportError:
    pass

# Check for Gmsh Python API
GMSH_AVAILABLE = False
try:
    import gmsh
    GMSH_AVAILABLE = True
except ImportError:
    pass


def _check_fipy_available(func_name):
    """Check if FiPy is available."""
    if not FIPY_AVAILABLE:
        raise RuntimeError(
            f"{func_name} requires FiPy which is not installed.\n"
            f"Install with: pip install fipy"
        )


def _check_gmsh_available(func_name):
    """Check if Gmsh Python API is available."""
    if not GMSH_AVAILABLE:
        raise RuntimeError(
            f"{func_name} requires the Gmsh Python package for complex mesh generation.\n"
            f"Install with: pip install gmsh"
        )


def _load_mesh_from_msh(mesh_filename, topo_gradient=None, tol=1e-8):
    """
    Load a Gmsh .msh file into FiPy and create cell-centered masks.
    
    Uses FiPy's Gmsh2D loader directly for better compatibility.
    
    Parameters
    ----------
    mesh_filename : str
        Path to the .msh file (must be MSH format 2.2)
    topo_gradient : float, optional
        Topographic gradient for computing surface elevation
    tol : float
        Tolerance for mask comparisons
        
    Returns
    -------
    mesh : fipy mesh
        FiPy mesh object
    surface : ndarray
        Boolean mask for surface cells
    sea_surface : ndarray or None
        Boolean mask for sea surface cells
    seawater : ndarray or None
        Boolean mask for seawater cells
    z_surface : ndarray
        Surface elevation at each cell center
    """
    _check_fipy_available('_load_mesh_from_msh')
    
    # Load mesh using FiPy's Gmsh2D
    fipy_mesh = fipy.Gmsh2D(mesh_filename)
    
    # Get cell centers
    cell_centers = fipy_mesh.cellCenters
    x_centers = np.array(cell_centers[0])
    y_centers = np.array(cell_centers[1])
    n_cells = len(x_centers)
    
    # Compute z_surface based on topographic gradient
    if topo_gradient is not None:
        z_surface = x_centers * topo_gradient
    else:
        z_surface = np.zeros(n_cells)
    
    # Build simple masks
    # Surface: cells whose y-coordinate is close to z_surface
    surface = np.abs(y_centers - z_surface) < tol
    
    # Sea surface: y ~ 0 and x < 0
    sea_surface = (np.abs(y_centers) < tol) & (x_centers < 0)
    
    # Seawater: x < 0 and below surface
    seawater = (x_centers < 0) & (y_centers <= z_surface + tol)
    
    return fipy_mesh, surface, sea_surface, seawater, z_surface


def setup_rectangular_mesh_fipy(Parameters, mesh_filename):
    """
    Create a rectangular mesh for the FiPy backend.

    A Gmsh geometry (.geo) file is written alongside the .msh file so that
    ``load_fipy_mesh_from_msh`` can pass it directly to ``Gmsh2D``.  When
    running under ``mpirun -np N``, FiPy / Gmsh automatically partitions the
    mesh into N subdomains from the .geo description without any extra work.

    A triangulated .msh file is also written (via meshio) for caching and
    visualisation.  If the .msh already exists it is reused on subsequent
    runs (grompy.py reuses the most recent matching mesh file).

    Parameters
    ----------
    Parameters : object
        Model parameters object with attributes:
        - L: domain length (x direction)
        - thickness: domain thickness (y direction)
        - cellsize / cellsize_x / cellsize_y: target cell size (m)
    mesh_filename : str
        Path where the .msh file will be written.  The .geo file is written
        to the same directory with the same stem and a .geo extension.

    Returns
    -------
    mesh : fipy.Grid2D
        FiPy Grid2D mesh (used by grompy.py; the actual solver re-loads via
        load_fipy_mesh_from_msh which will prefer the .geo file).
    surface : ndarray
        Boolean mask for surface (top-boundary) cells.
    sea_surface : None
        Not applicable for a rectangular mesh.
    seawater : None
        Not applicable for a rectangular mesh.
    z_surface : ndarray
        Surface elevation array (constant = thickness).
    """
    _check_fipy_available('setup_rectangular_mesh_fipy')

    # Determine cell size (support both cellsize and cellsize_x/cellsize_y)
    cellsize = getattr(Parameters, 'cellsize',
                       getattr(Parameters, 'cellsize_x', None))
    if cellsize is None:
        raise AttributeError("Parameters must have 'cellsize' or 'cellsize_x'")
    cellsize_y = getattr(Parameters, 'cellsize_y', cellsize)

    nx = int(math.ceil(Parameters.L / cellsize))
    ny = int(math.ceil(Parameters.thickness / cellsize_y))
    dx = Parameters.L / nx
    dy = Parameters.thickness / ny

    # ------------------------------------------------------------------
    # Build Gmsh geometry string for the rectangle (triangles, no Recombine)
    # ------------------------------------------------------------------
    # characteristic length at each corner point; Gmsh will mesh with
    # approximately this edge length throughout the domain.
    cl = cellsize
    L = Parameters.L
    H = Parameters.thickness

    geo_str = (
        "// Rectangular domain for GroMPy FiPy backend\n"
        f"cl = {cl!r};\n"
        f"L  = {L!r};\n"
        f"H  = {H!r};\n"
        "Point(1) = {0, 0, 0, cl};\n"
        "Point(2) = {L, 0, 0, cl};\n"
        "Point(3) = {L, H, 0, cl};\n"
        "Point(4) = {0, H, 0, cl};\n"
        "Line(1) = {1, 2};\n"
        "Line(2) = {2, 3};\n"
        "Line(3) = {3, 4};\n"
        "Line(4) = {4, 1};\n"
        "Line Loop(1) = {1, 2, 3, 4};\n"
        "Plane Surface(1) = {1};\n"
        # No 'Recombine Surface' → Gmsh produces triangles (default)
    )

    # Ensure output directory exists
    mesh_dir = os.path.dirname(mesh_filename)
    if mesh_dir and not os.path.exists(mesh_dir):
        os.makedirs(mesh_dir)

    # Write .geo file (used by load_fipy_mesh_from_msh → Gmsh2D)
    geo_filename = os.path.splitext(mesh_filename)[0] + '.geo'
    try:
        with open(geo_filename, 'w') as fh:
            fh.write(geo_str)
        print(f"Gmsh geometry written to: {geo_filename}")
    except Exception as e:
        print(f"Warning: Failed to write .geo file: {e}")

    # ------------------------------------------------------------------
    # Also write a triangulated .msh for caching / visualisation (meshio)
    # ------------------------------------------------------------------
    # Create FiPy Grid2D mesh to obtain vertex/cell connectivity
    mesh = fipy.Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)

    try:
        import meshio

        vertex_coords = mesh.vertexCoords.T  # (n_vertices, 2)
        points = np.column_stack([vertex_coords,
                                  np.zeros(vertex_coords.shape[0])])

        quads = mesh._orderedCellVertexIDs.T  # (n_quads, 4)

        # Split each quad into two triangles
        tri_a = quads[:, [0, 1, 2]]
        tri_b = quads[:, [0, 2, 3]]
        triangles = np.vstack([tri_a, tri_b])

        mesh_data = meshio.Mesh(points, [("triangle", triangles)])
        # Use the gmsh22 writer directly with binary=False.
        # meshio.Mesh.write() does not reliably forward the binary=False kwarg
        # to the underlying format writer (it defaults to binary=True), which
        # produces a binary .msh that FiPy's Gmsh2D cannot parse.
        from meshio.gmsh._gmsh22 import write as _write_gmsh22
        _write_gmsh22(mesh_filename, mesh_data, binary=False)
        print(f"Rectangular mesh (.msh) saved to: {mesh_filename}")

    except ImportError:
        print("Warning: meshio not available; .msh file not written.")
    except Exception as e:
        print(f"Warning: Failed to save .msh file: {e}")

    # ------------------------------------------------------------------
    # Build masks from Grid2D cell centres
    # ------------------------------------------------------------------
    cell_centers = mesh.cellCenters
    x_centers = np.array(cell_centers[0])
    y_centers = np.array(cell_centers[1])

    # Surface mask: top row of cells (y centre ≈ thickness - dy/2)
    z_surface_value = Parameters.thickness
    surface = np.abs(y_centers - (z_surface_value - dy / 2)) < dy / 2

    z_surface = np.ones_like(x_centers) * z_surface_value

    sea_surface = None
    seawater = None

    return mesh, surface, sea_surface, seawater, z_surface


def setup_standard_mesh_fipy(Parameters, mesh_filename):
    """
    Create a mesh with bi-linear surface topography using Gmsh Python API.
    
    This replicates the geometry from setup_standard_mesh in mesh_functions.py
    but uses the Gmsh Python API instead of esys.pycad.
    
    Parameters
    ----------
    Parameters : object
        Model parameters object with attributes:
        - L: domain length
        - thickness: domain thickness
        - x_topo_break: x location of topographic break
        - topo_gradient: topographic gradient before break
        - topo_gradient_hinterland: topographic gradient after break (optional)
        - topo_break: boolean, whether there is a topographic break
        - cellsize: element size
    mesh_filename : str
        Path to save the generated .msh file
        
    Returns
    -------
    mesh : fipy mesh
        FiPy mesh object
    surface : ndarray
        Boolean mask for surface cells
    sea_surface : None
        Not applicable for standard mesh
    seawater : None  
        Not applicable for standard mesh
    z_surface : ndarray
        Surface elevation at each cell center
    """
    _check_fipy_available('setup_standard_mesh_fipy')
    _check_gmsh_available('setup_standard_mesh_fipy')
    
    # Handle topographic gradient
    if not getattr(Parameters, 'topo_break', False):
        topo_gradient_hinterland = Parameters.topo_gradient
    else:
        topo_gradient_hinterland = getattr(Parameters, 'topo_gradient_hinterland', Parameters.topo_gradient)
    
    # Calculate geometry points (same as mesh_functions.py)
    x_topo_break = getattr(Parameters, 'x_topo_break', Parameters.L / 2)
    
    z1 = x_topo_break * Parameters.topo_gradient
    z2 = z1 + (Parameters.L - x_topo_break) * topo_gradient_hinterland
    
    xs = np.array([0, x_topo_break, Parameters.L, Parameters.L, x_topo_break, 0])
    zs = np.array([0, z1, z2, z2 - Parameters.thickness, z1 - Parameters.thickness, -Parameters.thickness])
    
    # Create mesh using Gmsh Python API
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)  # Suppress output
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  # Use MSH 2.2 for FiPy compatibility
    gmsh.model.add("standard_mesh")
    
    try:
        # Create points
        points = []
        for i, (x, z) in enumerate(zip(xs, zs)):
            points.append(gmsh.model.geo.addPoint(x, z, 0, Parameters.cellsize))
        
        # Create lines
        lines = []
        for i in range(len(points)):
            p1 = points[i]
            p2 = points[(i + 1) % len(points)]
            lines.append(gmsh.model.geo.addLine(p1, p2))
        
        # Create curve loop and surface
        curve_loop = gmsh.model.geo.addCurveLoop(lines)
        surface = gmsh.model.geo.addPlaneSurface([curve_loop])
        
        # Add physical groups
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [surface], name="land")
        gmsh.model.addPhysicalGroup(1, [lines[0], lines[1]], name="land_surface")
        
        # Generate mesh
        gmsh.model.mesh.generate(2)
        
        # Ensure directory exists
        mesh_dir = os.path.dirname(mesh_filename)
        if mesh_dir and not os.path.exists(mesh_dir):
            os.makedirs(mesh_dir)
        
        # Save mesh
        gmsh.write(mesh_filename)
        
    finally:
        gmsh.finalize()
    
    # Load mesh into FiPy
    mesh, surface_mask, sea_surface, seawater, z_surface = _load_mesh_from_msh(
        mesh_filename, topo_gradient=None, tol=Parameters.cellsize / 2
    )
    
    # Recompute z_surface with bi-linear topography
    cell_centers = mesh.cellCenters
    x_centers = np.array(cell_centers[0])
    y_centers = np.array(cell_centers[1])
    
    # Compute z_surface based on bi-linear topography
    z_surface = np.where(
        x_centers <= x_topo_break,
        x_centers * Parameters.topo_gradient,
        z1 + (x_centers - x_topo_break) * topo_gradient_hinterland
    )
    
    # Surface mask: cells whose y-coordinate is close to z_surface
    tol = Parameters.cellsize / 2
    surface_mask = np.abs(y_centers - z_surface) < tol
    
    return mesh, surface_mask, None, None, z_surface


def setup_coastal_mesh_fipy(Parameters, mesh_filename, extend_domain=True, max_length=1e5):
    """
    Create a coastal mesh with variable refinement around the salt-water interface.
    
    Uses the Ghyben-Herzberg equation to estimate the salt-water toe location
    and creates a mesh with finer resolution around this interface.
    
    This replicates setup_coastal_mesh from mesh_functions.py using Gmsh Python API.
    
    Parameters
    ----------
    Parameters : object
        Model parameters with coastal aquifer properties
    mesh_filename : str
        Path to save the generated .msh file
    extend_domain : bool
        Whether to extend domain if salt-water toe exceeds initial length
    max_length : float
        Maximum allowed domain length
        
    Returns
    -------
    mesh : fipy mesh
        FiPy mesh object
    surface : ndarray
        Boolean mask for surface cells
    sea_surface : ndarray
        Boolean mask for sea surface cells
    seawater : ndarray
        Boolean mask for seawater cells
    z_surface : ndarray
        Surface elevation at each cell center
    """
    _check_fipy_available('setup_coastal_mesh_fipy')
    _check_gmsh_available('setup_coastal_mesh_fipy')
    
    # Calculate salt-water toe extent using Ghyben-Herzberg
    R = Parameters.recharge_flux
    B = Parameters.thickness
    K = Parameters.k * Parameters.rho_f_0 * 9.81 / Parameters.viscosity
    L = Parameters.L
    
    # Calculate hydraulic head using analytical solution
    xa = np.linspace(0, L, 101)
    h = R / (K * B) * (L * xa - 0.5 * xa**2)
    
    fine_mesh = False
    adjust_length = False
    
    if h[-1] < (B / 40):
        print('warning, calculated extent salt water toe exceeds model domain')
        print('calculated h at model bnd = %0.2f m' % h[-1])
        print('Ghyben-Herzberg depth of salt water interface = %0.2f m' % (h[-1] * 40))
        print('thickness = %0.2f m' % Parameters.thickness)
        
        if extend_domain is False:
            extent_salt_water = Parameters.L - Parameters.buffer_distance_land - 1.0
            print('entire top right triangle at landward side of model domain has fine discretization')
            fine_mesh = True
        else:
            adjust_length = True
            # Need to calculate extent for adjusted domain
            extent_salt_water = Parameters.L - Parameters.buffer_distance_land - 1.0
    else:
        # Salt water toe touches bottom of model domain
        a = 0.5 * R / (K * B)
        b = -(R * L) / (K * B)
        c = B / 40.0
        
        D = np.sqrt(b**2 - 4 * a * c)
        extent_salt_water = (-b - D) / (2 * a)
        print('calculated extent salt water toe = %0.2f m' % extent_salt_water)
    
    if adjust_length:
        L_land = extent_salt_water + Parameters.buffer_distance_land * 2
        if L_land > max_length:
            L_land = Parameters.L
            fine_mesh = True
        else:
            print('extending model domain size to %0.3e' % L_land)
    else:
        L_land = Parameters.L
    
    # Build mesh geometry (4 blocks like the escript version)
    L_sea = Parameters.L_sea
    buffer_sea = Parameters.buffer_distance_sea
    buffer_land = Parameters.buffer_distance_land
    topo_gradient = Parameters.topo_gradient
    thickness = Parameters.thickness
    
    # Define the 10 key points
    xs = np.array([
        -L_sea,                                    # 0: sea bottom-left
        extent_salt_water - buffer_sea,            # 1: salt wedge sea side bottom
        -buffer_sea,                               # 2: near coast bottom
        -L_sea,                                    # 3: sea top-left
        extent_salt_water,                         # 4: salt wedge bottom
        0,                                         # 5: coast bottom
        extent_salt_water + buffer_land,           # 6: salt wedge land side bottom
        buffer_land,                               # 7: near coast land side bottom
        L_land,                                    # 8: land bottom-right
        L_land,                                    # 9: land top-right
    ])
    
    zs = xs * topo_gradient
    zs[0:2] = zs[0:2] - thickness
    zs[4] = zs[4] - thickness
    zs[6] = zs[6] - thickness
    zs[8] = zs[8] - thickness
    
    # Create mesh using Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  # Use MSH 2.2 for FiPy compatibility
    gmsh.model.add("coastal_mesh")
    
    try:
        cellsize = Parameters.cellsize
        ref_factor = getattr(Parameters, 'grid_refinement_factor', 0.5)
        ref_factor_sea = getattr(Parameters, 'grid_refinement_factor_sea', 1.0)
        
        # Create points with mesh size
        points = []
        for i, (x, z) in enumerate(zip(xs, zs)):
            points.append(gmsh.model.geo.addPoint(x, z, 0, cellsize))
        
        # Create lines connecting the points for the 4 surfaces
        # Surface A (sea): 0-1-2-3-0
        line1 = gmsh.model.geo.addLine(points[0], points[1])
        line2 = gmsh.model.geo.addLine(points[1], points[2])
        line3 = gmsh.model.geo.addLine(points[2], points[3])
        line4 = gmsh.model.geo.addLine(points[3], points[0])
        
        # Surface B (salt wedge sea side): 1-4-5-2-1
        line5 = gmsh.model.geo.addLine(points[1], points[4])
        line6 = gmsh.model.geo.addLine(points[4], points[5])
        line7 = gmsh.model.geo.addLine(points[5], points[2])
        
        # Surface C (salt wedge land side): 4-6-7-5-4
        line8 = gmsh.model.geo.addLine(points[4], points[6])
        line9 = gmsh.model.geo.addLine(points[6], points[7])
        line10 = gmsh.model.geo.addLine(points[7], points[5])
        
        # Surface D (land): 6-8-9-7-6
        line11 = gmsh.model.geo.addLine(points[6], points[8])
        line12 = gmsh.model.geo.addLine(points[8], points[9])
        line13 = gmsh.model.geo.addLine(points[9], points[7])
        
        # Create curve loops and surfaces
        curve_a = gmsh.model.geo.addCurveLoop([line1, line2, line3, line4])
        curve_b = gmsh.model.geo.addCurveLoop([line5, line6, line7, -line2])
        curve_c = gmsh.model.geo.addCurveLoop([line8, line9, line10, -line6])
        curve_d = gmsh.model.geo.addCurveLoop([line11, line12, line13, -line9])
        
        surface_a = gmsh.model.geo.addPlaneSurface([curve_a])
        surface_b = gmsh.model.geo.addPlaneSurface([curve_b])
        surface_c = gmsh.model.geo.addPlaneSurface([curve_c])
        surface_d = gmsh.model.geo.addPlaneSurface([curve_d])
        
        gmsh.model.geo.synchronize()
        
        # Set mesh sizes for different regions
        # Get points in each surface for mesh refinement
        gmsh.model.mesh.field.add("MathEval", 1)
        gmsh.model.mesh.field.setString(1, "F", str(cellsize * ref_factor_sea))
        
        gmsh.model.mesh.field.add("MathEval", 2)
        gmsh.model.mesh.field.setString(2, "F", str(cellsize * ref_factor))
        
        # Use Box fields for regional refinement
        # Sea region (coarser)
        gmsh.model.mesh.field.add("Box", 3)
        gmsh.model.mesh.field.setNumber(3, "VIn", cellsize * ref_factor_sea)
        gmsh.model.mesh.field.setNumber(3, "VOut", cellsize)
        gmsh.model.mesh.field.setNumber(3, "XMin", -L_sea)
        gmsh.model.mesh.field.setNumber(3, "XMax", -buffer_sea)
        gmsh.model.mesh.field.setNumber(3, "YMin", -thickness * 2)
        gmsh.model.mesh.field.setNumber(3, "YMax", thickness)
        
        # Refined region around salt wedge
        gmsh.model.mesh.field.add("Box", 4)
        gmsh.model.mesh.field.setNumber(4, "VIn", cellsize * ref_factor)
        gmsh.model.mesh.field.setNumber(4, "VOut", cellsize)
        gmsh.model.mesh.field.setNumber(4, "XMin", -buffer_sea)
        gmsh.model.mesh.field.setNumber(4, "XMax", extent_salt_water + buffer_land)
        gmsh.model.mesh.field.setNumber(4, "YMin", -thickness * 2)
        gmsh.model.mesh.field.setNumber(4, "YMax", thickness)
        
        if fine_mesh:
            # Land region also refined
            gmsh.model.mesh.field.add("Box", 5)
            gmsh.model.mesh.field.setNumber(5, "VIn", cellsize * ref_factor)
            gmsh.model.mesh.field.setNumber(5, "VOut", cellsize)
            gmsh.model.mesh.field.setNumber(5, "XMin", extent_salt_water + buffer_land)
            gmsh.model.mesh.field.setNumber(5, "XMax", L_land)
            gmsh.model.mesh.field.setNumber(5, "YMin", -thickness * 2)
            gmsh.model.mesh.field.setNumber(5, "YMax", thickness)
            
            gmsh.model.mesh.field.add("Min", 6)
            gmsh.model.mesh.field.setNumbers(6, "FieldsList", [3, 4, 5])
            gmsh.model.mesh.field.setAsBackgroundMesh(6)
        else:
            gmsh.model.mesh.field.add("Min", 6)
            gmsh.model.mesh.field.setNumbers(6, "FieldsList", [3, 4])
            gmsh.model.mesh.field.setAsBackgroundMesh(6)
        
        # Add physical groups
        gmsh.model.addPhysicalGroup(2, [surface_a], name="sea")
        gmsh.model.addPhysicalGroup(2, [surface_b], name="salt_wedge_sea_side")
        gmsh.model.addPhysicalGroup(2, [surface_c], name="salt_wedge_land_side")
        gmsh.model.addPhysicalGroup(2, [surface_d], name="land")
        gmsh.model.addPhysicalGroup(1, [line3], name="sea_surface1")
        gmsh.model.addPhysicalGroup(1, [line7], name="sea_surface2")
        gmsh.model.addPhysicalGroup(1, [line10], name="land_surface1")
        gmsh.model.addPhysicalGroup(1, [line13], name="land_surface2")
        
        # Generate mesh
        gmsh.model.mesh.generate(2)
        
        # Ensure directory exists
        mesh_dir = os.path.dirname(mesh_filename)
        if mesh_dir and not os.path.exists(mesh_dir):
            os.makedirs(mesh_dir)
        
        gmsh.write(mesh_filename)
        
    finally:
        gmsh.finalize()
    
    # Load mesh into FiPy
    mesh, surface_mask, sea_surface, seawater, z_surface = _load_mesh_from_msh(
        mesh_filename, topo_gradient=topo_gradient, tol=cellsize / 2
    )
    
    # Recompute masks based on cell centers
    cell_centers = mesh.cellCenters
    x_centers = np.array(cell_centers[0])
    y_centers = np.array(cell_centers[1])
    
    # z_surface based on topographic gradient
    z_surface = x_centers * topo_gradient
    
    # Surface mask: cells at the surface
    tol = cellsize / 2
    surface_mask = np.abs(y_centers - z_surface) < tol
    
    # Sea surface: y ~ 0 and x < 0
    sea_surface = (np.abs(y_centers) < tol) & (x_centers < 0)
    
    # Seawater: x < 0 and below surface
    seawater = (x_centers < 0) & (y_centers <= z_surface + tol)
    
    return mesh, surface_mask, sea_surface, seawater, z_surface


def setup_coastal_mesh_glover1959_fipy(Parameters, mesh_filename):
    """
    Create a coastal mesh using the Glover (1959) analytical solution for the 
    fresh-salt water interface position.
    
    This replicates setup_coastal_mesh_glover1959 from mesh_functions.py using 
    Gmsh Python API.
    
    Parameters
    ----------
    Parameters : object
        Model parameters with coastal aquifer properties including:
        - ghyben_herzberg: bool, use analytical interface position
        - recharge_flux, thickness, k, viscosity, L, etc.
    mesh_filename : str
        Path to save the generated .msh file
        
    Returns
    -------
    mesh : fipy mesh
        FiPy mesh object
    surface : ndarray
        Boolean mask for surface cells
    sea_surface : None
        Not computed in original
    seawater : None
        Not computed in original
    z_surface : ndarray
        Surface elevation at each cell center
    """
    _check_fipy_available('setup_coastal_mesh_glover1959_fipy')
    _check_gmsh_available('setup_coastal_mesh_glover1959_fipy')
    
    # Calculate salt-water interface using Glover (1959) if requested
    if getattr(Parameters, 'ghyben_herzberg', False):
        R = Parameters.recharge_flux
        B = Parameters.thickness
        K = Parameters.k * Parameters.rho_f_0 * 9.81 / Parameters.viscosity
        L = Parameters.L
        
        xa = np.linspace(0, L, 1001)
        h = R / (K * B) * (L * xa - 0.5 * xa ** 2)
        
        # Import Glover 1959 function
        from .grompy_lib import depth_sw_interface_Glover1959
        
        rho_f = Parameters.rho_f_0 * Parameters.freshwater_concentration * Parameters.gamma + Parameters.rho_f_0
        rho_s = Parameters.rho_f_0 * Parameters.seawater_concentration * Parameters.gamma + Parameters.rho_f_0
        
        Qmax = Parameters.recharge_flux * L
        
        y_sw, int_sw_top, int_sw_bottom = depth_sw_interface_Glover1959(
            xa, Parameters.k, Parameters.viscosity,
            Parameters.topo_gradient, Parameters.thickness,
            rho_f, rho_s, Parameters.gamma, Qmax=Qmax
        )
        
        if Parameters.recharge_flux == 0.0:
            h = Parameters.topo_gradient * xa
        
        if int_sw_bottom > L:
            print('warning, calculated extent salt water toe exceeds model domain')
            print('calculated toe of fresh_salt water bnd = %0.2f m' % int_sw_bottom)
            extent_salt_water = Parameters.L - Parameters.buffer_distance_land - 1.0
            print('choosing maximum possible extent of %0.2f m for designing model grid' % extent_salt_water)
            fine_mesh = False
        else:
            fine_mesh = False
            extent_salt_water = int_sw_bottom
            print('calculated extent salt water toe = %0.2f m' % int_sw_bottom)
    else:
        print('assuming a vertical fresh-salt water interface')
        extent_salt_water = 0.0
        fine_mesh = False
    
    L_land = Parameters.L
    
    # Build mesh geometry (same structure as setup_coastal_mesh_fipy)
    L_sea = Parameters.L_sea
    buffer_sea = Parameters.buffer_distance_sea
    buffer_land = Parameters.buffer_distance_land
    topo_gradient = Parameters.topo_gradient
    thickness = Parameters.thickness
    
    xs = np.array([
        -L_sea,
        extent_salt_water - buffer_sea,
        -buffer_sea,
        -L_sea,
        extent_salt_water,
        0,
        extent_salt_water + buffer_land,
        buffer_land,
        L_land,
        L_land,
    ])
    
    zs = xs * topo_gradient
    zs[0:2] = zs[0:2] - thickness
    zs[4] = zs[4] - thickness
    zs[6] = zs[6] - thickness
    zs[8] = zs[8] - thickness
    
    # Create mesh using Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  # Use MSH 2.2 for FiPy compatibility
    gmsh.model.add("coastal_glover1959_mesh")
    
    try:
        cellsize = Parameters.cellsize
        ref_factor = getattr(Parameters, 'grid_refinement_factor', 0.5)
        ref_factor_sea = getattr(Parameters, 'grid_refinement_factor_sea', 1.0)
        
        # Create points
        points = []
        for i, (x, z) in enumerate(zip(xs, zs)):
            points.append(gmsh.model.geo.addPoint(x, z, 0, cellsize))
        
        # Create lines (same as setup_coastal_mesh_fipy)
        line1 = gmsh.model.geo.addLine(points[0], points[1])
        line2 = gmsh.model.geo.addLine(points[1], points[2])
        line3 = gmsh.model.geo.addLine(points[2], points[3])
        line4 = gmsh.model.geo.addLine(points[3], points[0])
        line5 = gmsh.model.geo.addLine(points[1], points[4])
        line6 = gmsh.model.geo.addLine(points[4], points[5])
        line7 = gmsh.model.geo.addLine(points[5], points[2])
        line8 = gmsh.model.geo.addLine(points[4], points[6])
        line9 = gmsh.model.geo.addLine(points[6], points[7])
        line10 = gmsh.model.geo.addLine(points[7], points[5])
        line11 = gmsh.model.geo.addLine(points[6], points[8])
        line12 = gmsh.model.geo.addLine(points[8], points[9])
        line13 = gmsh.model.geo.addLine(points[9], points[7])
        
        # Create surfaces
        curve_a = gmsh.model.geo.addCurveLoop([line1, line2, line3, line4])
        curve_b = gmsh.model.geo.addCurveLoop([line5, line6, line7, -line2])
        curve_c = gmsh.model.geo.addCurveLoop([line8, line9, line10, -line6])
        curve_d = gmsh.model.geo.addCurveLoop([line11, line12, line13, -line9])
        
        surface_a = gmsh.model.geo.addPlaneSurface([curve_a])
        surface_b = gmsh.model.geo.addPlaneSurface([curve_b])
        surface_c = gmsh.model.geo.addPlaneSurface([curve_c])
        surface_d = gmsh.model.geo.addPlaneSurface([curve_d])
        
        gmsh.model.geo.synchronize()
        
        # Set up mesh refinement fields (same as coastal mesh)
        gmsh.model.mesh.field.add("Box", 3)
        gmsh.model.mesh.field.setNumber(3, "VIn", cellsize * ref_factor_sea)
        gmsh.model.mesh.field.setNumber(3, "VOut", cellsize)
        gmsh.model.mesh.field.setNumber(3, "XMin", -L_sea)
        gmsh.model.mesh.field.setNumber(3, "XMax", -buffer_sea)
        gmsh.model.mesh.field.setNumber(3, "YMin", -thickness * 2)
        gmsh.model.mesh.field.setNumber(3, "YMax", thickness)
        
        gmsh.model.mesh.field.add("Box", 4)
        gmsh.model.mesh.field.setNumber(4, "VIn", cellsize * ref_factor)
        gmsh.model.mesh.field.setNumber(4, "VOut", cellsize)
        gmsh.model.mesh.field.setNumber(4, "XMin", -buffer_sea)
        gmsh.model.mesh.field.setNumber(4, "XMax", extent_salt_water + buffer_land)
        gmsh.model.mesh.field.setNumber(4, "YMin", -thickness * 2)
        gmsh.model.mesh.field.setNumber(4, "YMax", thickness)
        
        if fine_mesh:
            gmsh.model.mesh.field.add("Box", 5)
            gmsh.model.mesh.field.setNumber(5, "VIn", cellsize * ref_factor)
            gmsh.model.mesh.field.setNumber(5, "VOut", cellsize)
            gmsh.model.mesh.field.setNumber(5, "XMin", extent_salt_water + buffer_land)
            gmsh.model.mesh.field.setNumber(5, "XMax", L_land)
            gmsh.model.mesh.field.setNumber(5, "YMin", -thickness * 2)
            gmsh.model.mesh.field.setNumber(5, "YMax", thickness)
            
            gmsh.model.mesh.field.add("Min", 6)
            gmsh.model.mesh.field.setNumbers(6, "FieldsList", [3, 4, 5])
            gmsh.model.mesh.field.setAsBackgroundMesh(6)
        else:
            gmsh.model.mesh.field.add("Min", 6)
            gmsh.model.mesh.field.setNumbers(6, "FieldsList", [3, 4])
            gmsh.model.mesh.field.setAsBackgroundMesh(6)
        
        # Add physical groups
        gmsh.model.addPhysicalGroup(2, [surface_a], name="sea")
        gmsh.model.addPhysicalGroup(2, [surface_b], name="salt_wedge_sea_side")
        gmsh.model.addPhysicalGroup(2, [surface_c], name="salt_wedge_land_side")
        gmsh.model.addPhysicalGroup(2, [surface_d], name="land")
        gmsh.model.addPhysicalGroup(1, [line3], name="sea_surface1")
        gmsh.model.addPhysicalGroup(1, [line7], name="sea_surface2")
        gmsh.model.addPhysicalGroup(1, [line10], name="land_surface1")
        gmsh.model.addPhysicalGroup(1, [line13], name="land_surface2")
        
        # Generate mesh
        gmsh.model.mesh.generate(2)
        
        # Ensure directory exists
        mesh_dir = os.path.dirname(mesh_filename)
        if mesh_dir and not os.path.exists(mesh_dir):
            os.makedirs(mesh_dir)
        
        gmsh.write(mesh_filename)
        
    finally:
        gmsh.finalize()
    
    # Load mesh into FiPy
    mesh, surface_mask, sea_surface, seawater, z_surface = _load_mesh_from_msh(
        mesh_filename, topo_gradient=topo_gradient, tol=cellsize / 2
    )
    
    # Recompute masks
    cell_centers = mesh.cellCenters
    x_centers = np.array(cell_centers[0])
    y_centers = np.array(cell_centers[1])
    
    z_surface = x_centers * topo_gradient
    
    tol = cellsize / 2
    surface_mask = np.abs(y_centers - z_surface) < tol
    
    # sea_surface and seawater are None in original implementation
    sea_surface = None
    seawater = None
    
    return mesh, surface_mask, sea_surface, seawater, z_surface


# Aliases to match the original mesh_functions.py interface
# These allow the module to be used as a drop-in replacement
setup_rectangular_mesh = setup_rectangular_mesh_fipy
setup_standard_mesh = setup_standard_mesh_fipy
setup_coastal_mesh = setup_coastal_mesh_fipy
setup_coastal_mesh_glover1959 = setup_coastal_mesh_glover1959_fipy
