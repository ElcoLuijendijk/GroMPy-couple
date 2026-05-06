"""
functions to create escript meshes for grompy-salt
 
"""

import os
import sys
import math
import numpy as np
import esys.pycad as pc
import esys.pycad.gmsh as gmsh
import esys.finley as fl
import esys.escript as es



def setup_coastal_mesh_glover1959(Parameters,
                                  mesh_filename):
    """
    Create a mesh consisting of 3 blocks, with a smaller cell size in the middle block

    The middle block is centered around the fresh-salt water interface, which is calculated using an anlytical
    solution by Glover (1959) Journal of Geophys. Res.
    """

    # if Parameters.topo_gradient == 0:
    #    extent_salt_water = 0
    # else:
    #    extent_salt_water = (Parameters.thickness /
    #                         (41.0 * Parameters.topo_gradient))

    if Parameters.ghyben_herzberg is True:

        R = Parameters.recharge_flux
        B = Parameters.thickness
        K = Parameters.k * Parameters.rho_f_0 * 9.81 / Parameters.viscosity
        L = Parameters.L

        # calculate hydraulic head using analytical solution for confined aq
        # with uniform recharge + Dupuit assumptions
        xa = np.linspace(0, L, 1001)
        h = R / (K * B) * (L * xa - 0.5 * xa ** 2)


        # calculate depth salt water interface
        from .grompy_lib import depth_sw_interface_Glover1959

        rho_f = Parameters.rho_f_0 * Parameters.freshwater_concentration * Parameters.gamma + Parameters.rho_f_0
        rho_s = Parameters.rho_f_0 * Parameters.seawater_concentration * Parameters.gamma + Parameters.rho_f_0

        Qmax = Parameters.recharge_flux * L

        y_sw, int_sw_top, int_sw_bottom = depth_sw_interface_Glover1959(xa, Parameters.k, Parameters.viscosity,
                                                                        Parameters.topo_gradient, Parameters.thickness,
                                                                        rho_f, rho_s, Parameters.gamma,
                                                                        Qmax=Qmax)

        if Parameters.recharge_flux == 0.0:
            # assume with no rehcarge that the hydraulic head follows the land surface
            h = Parameters.topo_gradient * xa
            #y_sw, int_sw_top, int_sw_bottom

        if int_sw_bottom > L:
            print('warning, calculated extent salt water toe exceeds model domain')
            print('calculated toe of fresh_salt water bnd = %0.2f m' % int_sw_bottom)
            extent_salt_water = Parameters.L - Parameters.buffer_distance_land - 1.0
            print('choosing maximum possible extent of %0.2f m  for ' \
                  'designing model grid' % extent_salt_water)

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

    ###############################
    # use gmsh to construct domain
    ##############################
    xs = np.array([-Parameters.L_sea,
                   extent_salt_water - Parameters.buffer_distance_sea,
                   -Parameters.buffer_distance_sea,
                   -Parameters.L_sea,
                   extent_salt_water,
                   0,
                   extent_salt_water + Parameters.buffer_distance_land,
                   Parameters.buffer_distance_land,
                   L_land,
                   L_land])

    zs = xs * Parameters.topo_gradient

    zs[0:2] = zs[0:2] - Parameters.thickness
    zs[4] = zs[4] - Parameters.thickness
    zs[6] = zs[6] - Parameters.thickness
    zs[8] = zs[8] - Parameters.thickness

    # points = create_points(xs,zs)
    points = [pc.Point(x, z) for x, z in zip(xs, zs)]

    line1 = pc.Line(points[0], points[1])
    line2 = pc.Line(points[1], points[2])
    line3 = pc.Line(points[2], points[3])
    line4 = pc.Line(points[3], points[0])
    line5 = pc.Line(points[1], points[4])
    line6 = pc.Line(points[4], points[5])
    line7 = pc.Line(points[5], points[2])
    line8 = pc.Line(points[4], points[6])
    line9 = pc.Line(points[6], points[7])
    line10 = pc.Line(points[7], points[5])
    line11 = pc.Line(points[6], points[8])
    line12 = pc.Line(points[8], points[9])
    line13 = pc.Line(points[9], points[7])

    # finer grid cell size around fresh-salt water interface
    curve_a = pc.CurveLoop(line1, line2, line3, line4)
    curve_b = pc.CurveLoop(line5, line6, line7, -line2)
    curve_c = pc.CurveLoop(line8, line9, line10, -line6)
    curve_d = pc.CurveLoop(line11, line12, line13, -line9)

    surface_a = pc.PlaneSurface(curve_a)
    surface_b = pc.PlaneSurface(curve_b)
    surface_c = pc.PlaneSurface(curve_c)
    surface_d = pc.PlaneSurface(curve_d)

    surface_a.setLocalScale(factor=Parameters.grid_refinement_factor_sea)
    surface_b.setLocalScale(factor=Parameters.grid_refinement_factor)
    surface_c.setLocalScale(factor=Parameters.grid_refinement_factor)

    if fine_mesh is True:
        print('assigning refined grid to entire landward side of model domain')
        surface_d.setLocalScale(factor=Parameters.grid_refinement_factor)

    d = gmsh.Design(dim=2, element_size=Parameters.cellsize)

    ps1 = pc.PropertySet("sea_surface1", line3)
    ps2 = pc.PropertySet("sea_surface2", line7)
    ps3 = pc.PropertySet("land_surface1", line10)
    ps4 = pc.PropertySet("land_surface2", line13)

    d.addItems(pc.PropertySet('sea', surface_a),
               pc.PropertySet('salt_wedge_sea_side', surface_b),
               pc.PropertySet('salt_wedge_land_side', surface_c),
               pc.PropertySet('land', surface_d),
               ps1, ps2, ps3, ps4)

    d.setMeshFileName(mesh_filename)

    print('=' * 30)

    mesh = fl.MakeDomain(d, optimizeLabeling=True)

    # calculate surface
    xy = mesh.getX()
    z_surface = xy[0] * Parameters.topo_gradient
    surface = es.whereZero(xy[1] - z_surface)

    # sea surface
    # sea_surface = surface * es.whereNegative(xy[0])
    sea_surface = None
    seawater = None

    return mesh, surface, sea_surface, seawater, z_surface


def setup_coastal_mesh(Parameters,
                       mesh_filename,
                       extend_domain=True,
                       max_length=1e5):
    """
    Create a mesh consisting of 3 blocks, with a smaller cell size in the middle block

    The middle block is centered around the fresh-salt water interface, which is calculated using the
    Ghyben-Herzberg equation.
    """

    #if Parameters.topo_gradient == 0:
    #    extent_salt_water = 0
    #else:
    #    extent_salt_water = (Parameters.thickness /
    #                         (41.0 * Parameters.topo_gradient))

    R = Parameters.recharge_flux
    B = Parameters.thickness
    K = Parameters.k * Parameters.rho_f_0 * 9.81 / Parameters.viscosity
    L = Parameters.L

    # calculate hydraulic head using analytical solution for confined aq
    # with uniform recharge + Dupuit assumptions
    xa = np.linspace(0, L, 101)
    h = R / (K * B) * (L * xa - 0.5 * xa**2)

    if h[-1] < (B / 40):
        print('warning, calculated extent salt water toe exceeds model domain')
        print('calculated h at model bnd = %0.2f m' % h[-1])
        print('Ghyben-Herzberg depth of salt water interface = %0.2f m' % (h[-1] * 40))
        print('thickness = %0.2f m' % Parameters.thickness)

        if extend_domain is False:
            extent_salt_water = Parameters.L - Parameters.buffer_distance_land - 1.0
            #print 'choosing maximum possible extent of %0.2f m  for ' \
            #      'designing model grid' % extent_salt_water

            print('entire top right triangle at landward side of model domain has fine discretization')
            fine_mesh = True
            adjust_length = False
        else:
            adjust_length = True

    else:

        fine_mesh = False

        # salt water toe touches bottom of model domain
        a = 0.5 * R / (K * B)
        b = -(R * L) / (K * B)
        c = B / 40.0

        D = np.sqrt(b**2 - 4 * a * c)
        extent_salt_water = (-b - D) / (2 * a)

        hs1 = R / (K * B) * (L * extent_salt_water -
                             0.5 * extent_salt_water**2)

        print('calculated extent salt water toe = %0.2f m' % extent_salt_water)

        try:
            assert np.abs(hs1 - B / 40.0) < 1e-3
        except AssertionError:
            msg = 'error, something wrong with calculated extent ' \
                  'salt water toe'
            raise ValueError(msg)

    if adjust_length is True:
        L_land  = extent_salt_water + Parameters.buffer_distance_land * 2
        if L_land > max_length:
            L_land = Parameters.L
            fine_mesh = True
        else:
            print('extending model domain size to %0.3e' % L_land)
    else:
        L_land = Parameters.L

    ###############################
    # use gmsh to construct domain
    ##############################
    xs = np.array([-Parameters.L_sea,
                   extent_salt_water - Parameters.buffer_distance_sea,
                   -Parameters.buffer_distance_sea,
                   -Parameters.L_sea,
                   extent_salt_water,
                   0,
                   extent_salt_water + Parameters.buffer_distance_land,
                   Parameters.buffer_distance_land,
                   L_land,
                   L_land,
                   -Parameters.L_sea])

    zs = xs * Parameters.topo_gradient

    zs[0:2] = zs[0:2] - Parameters.thickness
    zs[4] = zs[4] - Parameters.thickness
    zs[6] = zs[6] - Parameters.thickness
    zs[8] = zs[8] - Parameters.thickness
    zs[10] = 0.0

    #points = create_points(xs,zs)
    points = [pc.Point(x, z) for x, z in zip(xs, zs)]

    line1 = pc.Line(points[0], points[1])
    line2 = pc.Line(points[1], points[2])
    line3 = pc.Line(points[2], points[3])
    line4 = pc.Line(points[3], points[0])
    line5 = pc.Line(points[1], points[4])
    line6 = pc.Line(points[4], points[5])
    line7 = pc.Line(points[5], points[2])
    line8 = pc.Line(points[4], points[6])
    line9 = pc.Line(points[6], points[7])
    line10 = pc.Line(points[7], points[5])
    line11 = pc.Line(points[6], points[8])
    line12 = pc.Line(points[8], points[9])
    line13 = pc.Line(points[9], points[7])

    # new lines for sea surface
    # seabottom, x=-buffer to x=0
    # -line3 (pt 3 to 2)
    # - line7 (pt 2 to pt 5)
    # coastline to x=0, z=0
    #line14 = pc.Line(points[5], points[10])
    # x=0, z=0 to x=0, z=sea bottom
    #line15 = pc.Line(points[10], points[3])
    # coastline to x=0, sea bottom

    # finer grid cell size around fresh-salt water interface
    curve_a = pc.CurveLoop(line1, line2, line3, line4)
    curve_b = pc.CurveLoop(line5, line6, line7, -line2)
    curve_c = pc.CurveLoop(line8, line9, line10, -line6)
    curve_d = pc.CurveLoop(line11, line12, line13, -line9)

    #curve_seawater = pc.CurveLoop(-line3, -line7, line14, line15)

    surface_a = pc.PlaneSurface(curve_a)
    surface_b = pc.PlaneSurface(curve_b)
    surface_c = pc.PlaneSurface(curve_c)
    surface_d = pc.PlaneSurface(curve_d)

    #surface_seawater = pc.PlaneSurface(curve_seawater)

    surface_a.setLocalScale(factor=Parameters.grid_refinement_factor_sea)
    surface_b.setLocalScale(factor=Parameters.grid_refinement_factor)
    surface_c.setLocalScale(factor=Parameters.grid_refinement_factor)

    #surface_seawater.setLocalScale(factor=Parameters.grid_refinement_factor_seawater)

    if fine_mesh is True:
        print('assigning refined grid from landward side of model domain')
        surface_d.setLocalScale(factor=Parameters.grid_refinement_factor)

    d = gmsh.Design(dim=2, element_size=Parameters.cellsize)

    ps1 = pc.PropertySet("sea_surface1", line3)
    ps2 = pc.PropertySet("sea_surface2", line7)
    ps3 = pc.PropertySet("land_surface1", line10)
    ps4 = pc.PropertySet("land_surface2", line13)
    #ps5 = pc.PropertySet("sea_surface", line14)

    #d.addItems(pc.PropertySet('sea', surface_a),
    #           pc.PropertySet('salt_wedge_sea_side', surface_b),
    #           pc.PropertySet('salt_wedge_land_side', surface_c),
    #           pc.PropertySet('land', surface_d),
    #           pc.PropertySet('seawater', surface_seawater),
    #           ps1, ps2, ps3, ps4, ps5)

    d.addItems(pc.PropertySet('sea', surface_a),
               pc.PropertySet('salt_wedge_sea_side', surface_b),
               pc.PropertySet('salt_wedge_land_side', surface_c),
               pc.PropertySet('land', surface_d),
               ps1, ps2, ps3, ps4)

    d.setMeshFileName(mesh_filename)

    mesh = fl.MakeDomain(d, optimizeLabeling=True)

    # calculate surface
    xy = mesh.getX()
    z_surface = xy[0] * Parameters.topo_gradient
    surface = es.whereZero(xy[1] - z_surface)

    # sea surface
    sea_surface = es.whereZero(xy[1]) * es.whereNegative(xy[0])

    seawater = es.whereNegative(xy[0]) * es.whereNegative(z_surface - xy[1])

    print(bla)

    return mesh, surface, sea_surface, seawater, z_surface


def setup_rectangular_mesh(Parameters,
                           mesh_filename):
    """
    Create a rectangular mesh.
    """
    nx = int(math.ceil(Parameters.L / Parameters.cellsize_x))
    ny = int(math.ceil(Parameters.thickness / Parameters.cellsize_y))
    mesh = fl.Rectangle(l0=Parameters.L, l1=Parameters.thickness,
                        n0=nx, n1=ny)

    # calculate surface
    xy = mesh.getX()
    z_surface = (xy[0] - xy[0] + 1) * Parameters.thickness
    surface = es.whereZero(xy[1] - z_surface)
    sea_surface = None
    seawater = None

    return mesh, surface, sea_surface, seawater, z_surface


def setup_standard_mesh_csv_topo(Parameters, mesh_filename):
    """
    Create a mesh with a land surface elevation profile read from a CSV file.

    The CSV file must have two columns named exactly 'distance' and 'elevation':
      - 'distance': horizontal distance from the left-hand boundary (m)
      - 'elevation': land surface elevation at that distance (m)

    The profile is linearly interpolated onto mesh node x-coordinates. The
    aquifer bottom is flat at z = -thickness (constant thickness domain).

    Required Parameters attributes:
      topo_csv_file  -- path to the CSV file (str)
      L              -- model domain length (m)
      thickness      -- aquifer thickness (m)
      cellsize       -- target element size (m)
    """
    import pandas as pd

    # ------------------------------------------------------------------
    # 1. Read and validate the CSV
    # ------------------------------------------------------------------
    topo_df = pd.read_csv(Parameters.topo_csv_file)

    required_cols = {'distance', 'elevation'}
    missing = required_cols - set(topo_df.columns)
    if missing:
        raise ValueError(
            "topo CSV file '%s' is missing required column(s): %s. "
            "The file must contain columns named 'distance' and 'elevation'."
            % (Parameters.topo_csv_file, ', '.join(sorted(missing)))
        )

    x_csv = topo_df['distance'].values.astype(float)
    z_csv = topo_df['elevation'].values.astype(float)

    # Sort by distance in case the CSV is not ordered
    sort_idx = np.argsort(x_csv)
    x_csv = x_csv[sort_idx]
    z_csv = z_csv[sort_idx]

    if x_csv[0] > 0.0:
        raise ValueError(
            "topo CSV file '%s' does not cover the full model domain: "
            "minimum distance is %g m but the domain starts at x=0."
            % (Parameters.topo_csv_file, x_csv[0])
        )
    if x_csv[-1] < Parameters.L:
        raise ValueError(
            "topo CSV file '%s' does not cover the full model domain: "
            "maximum distance is %g m but the domain ends at x=%g m (Parameters.L)."
            % (Parameters.topo_csv_file, x_csv[-1], Parameters.L)
        )

    # Restrict to [0, L] in case the CSV extends beyond the domain
    mask = (x_csv >= 0.0) & (x_csv <= Parameters.L)
    # Always include exact endpoints
    x_csv = x_csv[mask]
    z_csv = z_csv[mask]

    # Ensure endpoints at exactly 0 and L are present
    if x_csv[0] != 0.0:
        z0 = float(np.interp(0.0, x_csv, z_csv))
        x_csv = np.concatenate([[0.0], x_csv])
        z_csv = np.concatenate([[z0], z_csv])
    if x_csv[-1] != Parameters.L:
        zL = float(np.interp(Parameters.L, x_csv, z_csv))
        x_csv = np.concatenate([x_csv, [Parameters.L]])
        z_csv = np.concatenate([z_csv, [zL]])

    # ------------------------------------------------------------------
    # 2. Build the gmsh geometry
    # ------------------------------------------------------------------
    # Corner points (bottom-left, bottom-right)
    pt_bottom_left  = pc.Point(0.0,           -Parameters.thickness)
    pt_bottom_right = pc.Point(Parameters.L,  -Parameters.thickness)

    # Top-surface points from the CSV profile
    top_pts = [pc.Point(float(x), float(z)) for x, z in zip(x_csv, z_csv)]

    # Boundary lines
    # Bottom
    line_bottom = pc.Line(pt_bottom_left, pt_bottom_right)
    # Right wall
    line_right = pc.Line(pt_bottom_right, top_pts[-1])
    # Top surface segments (one per consecutive pair of CSV points)
    top_lines = [pc.Line(top_pts[i + 1], top_pts[i])
                 for i in range(len(top_pts) - 1)]
    # Left wall
    line_left = pc.Line(top_pts[0], pt_bottom_left)

    # Single curve loop: bottom -> right wall -> top (right to left) -> left wall
    curve = pc.CurveLoop(
        line_bottom,
        line_right,
        *top_lines,
        line_left
    )

    domain_surface = pc.PlaneSurface(curve)

    d = gmsh.Design(dim=2, element_size=Parameters.cellsize)

    # Tag all top-surface segments under one property set
    d.addItems(
        pc.PropertySet('land', domain_surface),
        pc.PropertySet('land_surface', *top_lines)
    )

    d.setMeshFileName(mesh_filename)

    print('=' * 30)

    mesh = fl.MakeDomain(d, optimizeLabeling=True)

    # ------------------------------------------------------------------
    # 3. Compute z_surface field on the mesh by interpolating the CSV profile
    # ------------------------------------------------------------------
    xy = mesh.getX()
    # es.Data objects support numpy via esys.escript; convert x-coords to numpy
    x_nodes = np.array(xy[0].toListOfTuples()).flatten()
    z_surface_np = np.interp(x_nodes, x_csv, z_csv)

    # Rebuild as an escript Data object on the same function space
    z_surface = es.Data(0.0, xy[0].getFunctionSpace())
    for i, zval in enumerate(z_surface_np):
        # use the interpolated numpy array to build an escript field
        pass
    # Simpler: use arithmetic on xy[0] — interpolate via a piecewise linear
    # construction using escript operations
    z_surface = es.Data(0.0, xy.getFunctionSpace())
    for i in range(len(x_csv) - 1):
        x0, x1 = float(x_csv[i]), float(x_csv[i + 1])
        z0, z1 = float(z_csv[i]), float(z_csv[i + 1])
        in_segment = es.whereNonNegative(xy[0] - x0) * es.whereNonPositive(xy[0] - x1)
        slope = (z1 - z0) / (x1 - x0)
        z_segment = z0 + slope * (xy[0] - x0)
        z_surface = z_surface + in_segment * z_segment

    surface = es.whereZero(xy[1] - z_surface)

    sea_surface = None
    seawater = None

    return mesh, surface, sea_surface, seawater, z_surface


def setup_standard_mesh(Parameters, mesh_filename):
    """
    Create a mesh with bi-linear surface topography

    """

    # if Parameters.topo_gradient == 0:
    #    extent_salt_water = 0
    # else:
    #    extent_salt_water = (Parameters.thickness /
    #                         (41.0 * Parameters.topo_gradient))


    if Parameters.topo_break == False:
        Parameters.topo_gradient_hinterland = Parameters.topo_gradient

    ###############################
    # use gmsh to construct domain
    ##############################
    xs = np.array([0, Parameters.x_topo_break, Parameters.L, Parameters.L, Parameters.x_topo_break, 0])

    z1 = Parameters.x_topo_break * Parameters.topo_gradient
    z2 = z1 + (Parameters.L - Parameters.x_topo_break) * Parameters.topo_gradient_hinterland
    zs = np.array([0, z1, z2, z2 - Parameters.thickness, z1 - Parameters.thickness, -Parameters.thickness ])

    # points = create_points(xs,zs)
    points = [pc.Point(x, z) for x, z in zip(xs, zs)]

    line1 = pc.Line(points[0], points[1])
    line2 = pc.Line(points[1], points[2])
    line3 = pc.Line(points[2], points[3])
    line4 = pc.Line(points[3], points[4])
    line5 = pc.Line(points[4], points[5])
    line6 = pc.Line(points[5], points[0])
    
    # finer grid cell size around fresh-salt water interface
    curve = pc.CurveLoop(line1, line2, line3, line4, line5, line6)
    
    surface = pc.PlaneSurface(curve)
    
    d = gmsh.Design(dim=2, element_size=Parameters.cellsize)

    ps1 = pc.PropertySet("land_surface1", line1)
    ps2 = pc.PropertySet("land_surface2", line2)

    d.addItems(pc.PropertySet('land', surface),
               ps1, ps2)

    d.setMeshFileName(mesh_filename)

    print('=' * 30)

    mesh = fl.MakeDomain(d, optimizeLabeling=True)

    # calculate surface
    xy = mesh.getX()
    after_break = es.wherePositive(xy[0] - Parameters.x_topo_break)
    x_after_break = xy[0] - Parameters.x_topo_break
    #z_extra_after_break = 
    z_surface1 = xy[0] * Parameters.topo_gradient
    z_surface2 = after_break * x_after_break * (Parameters.topo_gradient_hinterland - Parameters.topo_gradient)

    #z_surface[after_break] = z1 + (xy[0][after_break] - Parameters.x_topo_break) * Parameters.topo_gradient[1]
    z_surface = z_surface1 + z_surface2

    surface = es.whereZero(xy[1] - z_surface)

    # sea surface
    # sea_surface = surface * es.whereNegative(xy[0])
    sea_surface = None
    seawater = None

    return mesh, surface, sea_surface, seawater, z_surface