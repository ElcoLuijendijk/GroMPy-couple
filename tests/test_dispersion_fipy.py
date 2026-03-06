"""
Tests for the FiPy tensor dispersion implementation.

Covers calculate_dispersion_coefficients_fipy from lib.grompy_fipy:
- Tensor structure (correct shape and symmetry)
- Longitudinal dominance along flow direction
- Non-trivial off-diagonal terms with angled flow
- Concentration change when the transport equation is solved
"""

import numpy as np
import pytest

fipy = pytest.importorskip('fipy')

from fipy import CellVariable, FaceVariable, Grid2D, DiffusionTermCorrection, ConvectionTerm, TransientTerm
from lib.grompy_fipy import calculate_dispersion_coefficients_fipy


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_mesh():
    """20×20 grid, 5 m cells."""
    return Grid2D(nx=20, ny=20, dx=5.0, dy=5.0)


@pytest.fixture
def uniform_horizontal_flow(simple_mesh):
    """Constant rightward Darcy flux (qx=1e-5 m/s, qy=0) as FaceVariables."""
    qx = FaceVariable(mesh=simple_mesh, value=1e-5)
    qy = FaceVariable(mesh=simple_mesh, value=0.0)
    return [qx, qy]


@pytest.fixture
def uniform_diagonal_flow(simple_mesh):
    """Constant 45° Darcy flux (qx=qy=1e-5/√2) as FaceVariables."""
    v = 1e-5 / np.sqrt(2)
    qx = FaceVariable(mesh=simple_mesh, value=v)
    qy = FaceVariable(mesh=simple_mesh, value=v)
    return [qx, qy]


@pytest.fixture
def uniform_horizontal_flow_vector(simple_mesh):
    """Rank-1 FaceVariable [qx, qy] needed by ConvectionTerm."""
    n_faces = simple_mesh.numberOfFaces
    q_vals = np.zeros((2, n_faces))
    q_vals[0, :] = 1e-5   # rightward
    return FaceVariable(mesh=simple_mesh, rank=1, value=q_vals)


@pytest.fixture
def dispersion_params():
    """Standard dispersion parameters."""
    return dict(porosity=0.25, diffusivity=1e-9, l_disp=50.0, t_disp=5.0)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestDispersionTensorStructure:
    """Verify the shape and basic properties of the returned tensor."""

    def test_tensor_shape(self, simple_mesh, uniform_horizontal_flow, dispersion_params):
        """Dispersion tensor must be 2×2 per cell."""
        tensor = calculate_dispersion_coefficients_fipy(
            simple_mesh, uniform_horizontal_flow, **dispersion_params
        )
        # FiPy rank-2 CellVariable: value shape is (2, 2, n_cells)
        assert tensor.value.shape == (2, 2, simple_mesh.numberOfCells), (
            f"Unexpected tensor shape: {tensor.value.shape}"
        )

    def test_tensor_symmetric(self, simple_mesh, uniform_horizontal_flow, dispersion_params):
        """Dispersion tensor must be symmetric (Dxy == Dyx)."""
        tensor = calculate_dispersion_coefficients_fipy(
            simple_mesh, uniform_horizontal_flow, **dispersion_params
        )
        Dxy = tensor.value[0, 1]
        Dyx = tensor.value[1, 0]
        np.testing.assert_allclose(
            Dxy, Dyx, rtol=1e-10,
            err_msg="Dispersion tensor is not symmetric: Dxy != Dyx"
        )

    def test_diagonal_entries_positive(self, simple_mesh, uniform_horizontal_flow, dispersion_params):
        """Diagonal entries Dxx and Dyy must be positive everywhere."""
        tensor = calculate_dispersion_coefficients_fipy(
            simple_mesh, uniform_horizontal_flow, **dispersion_params
        )
        assert np.all(tensor.value[0, 0] > 0), "Dxx has non-positive values"
        assert np.all(tensor.value[1, 1] > 0), "Dyy has non-positive values"


class TestDispersionTensorMagnitude:
    """Verify that dispersion magnitudes reflect the flow direction."""

    def test_longitudinal_dominates_horizontal_flow(
        self, simple_mesh, uniform_horizontal_flow, dispersion_params
    ):
        """With purely horizontal flow, Dxx (longitudinal) > Dyy (transverse)."""
        tensor = calculate_dispersion_coefficients_fipy(
            simple_mesh, uniform_horizontal_flow, **dispersion_params
        )
        Dxx_mean = tensor.value[0, 0].mean()
        Dyy_mean = tensor.value[1, 1].mean()
        assert Dxx_mean > Dyy_mean, (
            f"Expected Dxx ({Dxx_mean:.2e}) > Dyy ({Dyy_mean:.2e}) "
            "for horizontal flow"
        )

    def test_dxx_dyy_equal_diagonal_flow(
        self, simple_mesh, uniform_diagonal_flow, dispersion_params
    ):
        """With 45° flow, Dxx and Dyy should be equal by symmetry."""
        tensor = calculate_dispersion_coefficients_fipy(
            simple_mesh, uniform_diagonal_flow, **dispersion_params
        )
        Dxx_mean = tensor.value[0, 0].mean()
        Dyy_mean = tensor.value[1, 1].mean()
        np.testing.assert_allclose(
            Dxx_mean, Dyy_mean, rtol=1e-6,
            err_msg="Dxx and Dyy should be equal for 45° flow"
        )

    def test_off_diagonal_nonzero_diagonal_flow(
        self, simple_mesh, uniform_diagonal_flow, dispersion_params
    ):
        """With angled flow the off-diagonal term Dxy must be non-zero."""
        tensor = calculate_dispersion_coefficients_fipy(
            simple_mesh, uniform_diagonal_flow, **dispersion_params
        )
        Dxy_abs_mean = np.abs(tensor.value[0, 1]).mean()
        assert Dxy_abs_mean > 0, "Off-diagonal Dxy should be non-zero for diagonal flow"

    def test_off_diagonal_zero_axis_aligned_flow(
        self, simple_mesh, uniform_horizontal_flow, dispersion_params
    ):
        """With purely horizontal flow, off-diagonal terms must be zero."""
        tensor = calculate_dispersion_coefficients_fipy(
            simple_mesh, uniform_horizontal_flow, **dispersion_params
        )
        Dxy = tensor.value[0, 1]
        np.testing.assert_allclose(
            Dxy, 0.0, atol=1e-20,
            err_msg="Off-diagonal Dxy should be zero for axis-aligned flow"
        )


class TestDispersionTransport:
    """Verify that the dispersion tensor produces physical concentration changes."""

    def test_concentration_changes_after_timestep(
        self, simple_mesh, uniform_horizontal_flow, uniform_horizontal_flow_vector,
        dispersion_params
    ):
        """Solving one transport timestep must produce a measurable concentration change."""
        initial_conc = np.zeros(simple_mesh.numberOfCells)
        left_cells = np.array(simple_mesh.cellCenters[0]) < 50.0
        initial_conc[left_cells] = 0.035   # seawater
        initial_conc[~left_cells] = 0.0    # freshwater

        concentration = CellVariable(
            mesh=simple_mesh, value=initial_conc.copy(), name='concentration'
        )

        tensor = calculate_dispersion_coefficients_fipy(
            simple_mesh, uniform_horizontal_flow, **dispersion_params
        )

        # ConvectionTerm requires a rank-1 FaceVariable (2, n_faces)
        eq = (
            TransientTerm(coeff=dispersion_params['porosity'])
            + ConvectionTerm(coeff=uniform_horizontal_flow_vector)
            == DiffusionTermCorrection(coeff=tensor)
        )

        dt = 86400.0 * 30  # 30 days
        eq.solve(var=concentration, dt=dt)

        conc_change = np.max(np.abs(np.array(concentration.value) - initial_conc))
        assert conc_change > 1e-10, (
            f"Concentration change {conc_change:.2e} too small — "
            "dispersion tensor may not be applied correctly"
        )

    def test_concentration_stays_bounded(
        self, simple_mesh, uniform_horizontal_flow, uniform_horizontal_flow_vector,
        dispersion_params
    ):
        """Concentration must remain within [0, initial_max] after transport."""
        initial_max = 0.035
        initial_conc = np.zeros(simple_mesh.numberOfCells)
        left_cells = np.array(simple_mesh.cellCenters[0]) < 50.0
        initial_conc[left_cells] = initial_max

        concentration = CellVariable(
            mesh=simple_mesh, value=initial_conc.copy(), name='concentration'
        )

        tensor = calculate_dispersion_coefficients_fipy(
            simple_mesh, uniform_horizontal_flow, **dispersion_params
        )

        eq = (
            TransientTerm(coeff=dispersion_params['porosity'])
            + ConvectionTerm(coeff=uniform_horizontal_flow_vector)
            == DiffusionTermCorrection(coeff=tensor)
        )

        eq.solve(var=concentration, dt=86400.0 * 30)

        conc_vals = np.array(concentration.value)
        assert conc_vals.min() >= -1e-6, (
            f"Concentration went negative: {conc_vals.min():.2e}"
        )
        assert conc_vals.max() <= initial_max + 1e-6, (
            f"Concentration exceeded initial max: {conc_vals.max():.2e}"
        )
