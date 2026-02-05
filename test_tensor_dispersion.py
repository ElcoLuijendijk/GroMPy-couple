#!/usr/bin/env python3
"""
Minimal test of FiPy tensor dispersion implementation
"""

import numpy as np
import sys
import os

# Add the lib directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'lib'))

try:
    import fipy
    from fipy import CellVariable, Grid2D, DiffusionTermCorrection, ConvectionTerm, TransientTerm
    from fipy.tools import numerix
    from lib.grompy_fipy import calculate_dispersion_coefficients_fipy
    print("FiPy and modules imported successfully")
except ImportError as e:
    print(f"Import error: {e}")
    sys.exit(1)

def test_tensor_dispersion():
    """Test that tensor dispersion creates concentration changes"""

    # Create a simple 2D grid mesh
    mesh = Grid2D(nx=20, ny=20, dx=5.0, dy=5.0)

    # Parameters
    porosity = 0.25
    diffusivity = 1e-9
    l_disp = 50.0  # longitudinal dispersivity
    t_disp = 5.0   # transverse dispersivity

    # Create initial concentration (high at left, low at right)
    initial_conc = np.zeros(mesh.numberOfCells)
    left_cells = mesh.cellCenters[0] < 50.0  # Left half
    initial_conc[left_cells] = 0.035  # seawater concentration
    initial_conc[~left_cells] = 0.0   # freshwater

    concentration = CellVariable(mesh=mesh, value=initial_conc, name='concentration')

    # Create constant velocity field (flowing rightward)
    vx_value = 1e-5  # m/s
    vy_value = 0.0   # no vertical flow

    vx = CellVariable(mesh=mesh, value=vx_value)
    vy = CellVariable(mesh=mesh, value=vy_value)
    velocity = [vx, vy]

    print(f"Initial concentration range: {concentration.min():.6f} to {concentration.max():.6f}")
    print(f"Velocity: vx={vx_value}, vy={vy_value}")

    # Calculate dispersion tensor
    dispersion_tensor = calculate_dispersion_coefficients_fipy(
        velocity, porosity, diffusivity, l_disp, t_disp
    )

    print("Dispersion tensor calculated")
    print(f"Dxx range: {dispersion_tensor.value[0,0].min():.2e} to {dispersion_tensor.value[0,0].max():.2e}")
    print(f"Dyy range: {dispersion_tensor.value[1,1].min():.2e} to {dispersion_tensor.value[1,1].max():.2e}")

    # Build the transport equation
    eq = (TransientTerm(coeff=porosity)
          + ConvectionTerm(coeff=velocity)
          == DiffusionTermCorrection(coeff=dispersion_tensor))

    # Take one timestep
    dt = 86400.0 * 30  # 30 days
    eq.solve(var=concentration, dt=dt)

    # Check concentration change
    conc_change = np.max(np.abs(np.array(concentration.value) - initial_conc))
    print(f"After {dt/86400:.0f} days: concentration change = {conc_change:.2e} kg/kg")

    if conc_change > 1e-10:
        print("SUCCESS: Tensor dispersion is working - concentration changed!")
        return True
    else:
        print("FAILURE: Concentration did not change significantly")
        return False

if __name__ == "__main__":
    success = test_tensor_dispersion()
    sys.exit(0 if success else 1)