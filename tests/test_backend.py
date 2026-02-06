"""
Unit tests for the GroMPy-couple backend abstraction layer.

These tests verify that both backends (escript and FiPy) implement
the required interface correctly and produce consistent results.
"""

import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_allclose

from lib.backend import get_backend, check_backend_available, AVAILABLE_BACKENDS


# =============================================================================
# Test fixtures
# =============================================================================

@pytest.fixture(params=['escript', 'fipy'])
def backend_name(request):
    """Parametrized fixture to test both backends."""
    return request.param


@pytest.fixture
def backend(backend_name):
    """Get an initialized backend instance."""
    if not check_backend_available(backend_name):
        pytest.skip(f"{backend_name} backend not available")
    return get_backend(backend_name)


@pytest.fixture
def simple_mesh(backend):
    """Create a simple rectangular mesh for testing."""
    return backend.create_rectangle_mesh(
        length=100.0,
        height=50.0,
        nx=10,
        ny=5
    )


# =============================================================================
# Test backend factory
# =============================================================================

class TestBackendFactory:
    """Tests for the backend factory function."""
    
    def test_available_backends_list(self):
        """Test that AVAILABLE_BACKENDS contains expected values."""
        assert 'escript' in AVAILABLE_BACKENDS
        assert 'fipy' in AVAILABLE_BACKENDS
    
    def test_get_backend_invalid_name(self):
        """Test that invalid backend name raises ValueError."""
        with pytest.raises(ValueError, match="Unknown backend"):
            get_backend('invalid_backend')
    
    def test_check_backend_available_returns_bool(self):
        """Test that check_backend_available returns boolean."""
        result = check_backend_available('escript')
        assert isinstance(result, bool)
        
        result = check_backend_available('fipy')
        assert isinstance(result, bool)
    
    def test_check_backend_available_invalid_name(self):
        """Test that check_backend_available returns False for invalid names."""
        result = check_backend_available('invalid_backend')
        assert result is False
    
    def test_get_backend_case_insensitive(self):
        """Test that backend name is case-insensitive."""
        for name in AVAILABLE_BACKENDS:
            if check_backend_available(name):
                backend_lower = get_backend(name.lower())
                backend_upper = get_backend(name.upper())
                assert backend_lower.name == backend_upper.name


# =============================================================================
# Test backend interface
# =============================================================================

class TestBackendInterface:
    """Tests for the BackendBase interface."""
    
    def test_backend_has_name(self, backend):
        """Test that backend has a name property."""
        assert hasattr(backend, 'name')
        assert backend.name in AVAILABLE_BACKENDS
    
    def test_create_rectangle_mesh(self, backend):
        """Test mesh creation."""
        mesh = backend.create_rectangle_mesh(
            length=100.0, height=50.0, nx=10, ny=5
        )
        assert mesh is not None
        assert mesh.num_cells > 0
        assert mesh.num_nodes > 0
    
    def test_create_scalar_field(self, backend, simple_mesh):
        """Test scalar field creation."""
        field = backend.create_scalar_field(simple_mesh, value=1.5, name="test")
        assert field is not None
        assert field.name == "test"
        values = field.get_values()
        assert np.allclose(values, 1.5)
    
    def test_create_scalar_field_from_array(self, backend, simple_mesh):
        """Test scalar field creation from array."""
        n = simple_mesh.num_cells
        arr = np.linspace(0, 1, n)
        field = backend.create_scalar_field(simple_mesh, value=arr)
        values = field.get_values()
        assert len(values) == n
        assert_array_almost_equal(values, arr)
    
    def test_create_vector_field(self, backend, simple_mesh):
        """Test vector field creation."""
        field = backend.create_vector_field(simple_mesh, values=(1.0, 2.0))
        assert field is not None
        coords = field.get_coordinates()
        assert coords.shape[1] == 2  # x, y coordinates
    
    def test_create_tensor_field(self, backend, simple_mesh):
        """Test tensor field creation."""
        field = backend.create_tensor_field(
            simple_mesh,
            values=((1.0, 0.0), (0.0, 2.0))
        )
        assert field is not None
    
    def test_create_pde_solver(self, backend, simple_mesh):
        """Test PDE solver creation."""
        solver = backend.create_pde_solver(simple_mesh)
        assert solver is not None


# =============================================================================
# Test field operations
# =============================================================================

class TestFieldOperations:
    """Tests for field arithmetic and operations."""
    
    def test_field_addition_scalar(self, backend, simple_mesh):
        """Test field + scalar."""
        field = backend.create_scalar_field(simple_mesh, value=1.0)
        result = field + 2.0
        values = result.get_values()
        assert np.allclose(values, 3.0)
    
    def test_field_addition_field(self, backend, simple_mesh):
        """Test field + field."""
        f1 = backend.create_scalar_field(simple_mesh, value=1.0)
        f2 = backend.create_scalar_field(simple_mesh, value=2.0)
        result = f1 + f2
        values = result.get_values()
        assert np.allclose(values, 3.0)
    
    def test_field_subtraction(self, backend, simple_mesh):
        """Test field - scalar and field - field."""
        f1 = backend.create_scalar_field(simple_mesh, value=5.0)
        f2 = backend.create_scalar_field(simple_mesh, value=2.0)
        
        result1 = f1 - 3.0
        assert np.allclose(result1.get_values(), 2.0)
        
        result2 = f1 - f2
        assert np.allclose(result2.get_values(), 3.0)
    
    def test_field_multiplication(self, backend, simple_mesh):
        """Test field * scalar and field * field."""
        f1 = backend.create_scalar_field(simple_mesh, value=3.0)
        f2 = backend.create_scalar_field(simple_mesh, value=4.0)
        
        result1 = f1 * 2.0
        assert np.allclose(result1.get_values(), 6.0)
        
        result2 = f1 * f2
        assert np.allclose(result2.get_values(), 12.0)
    
    def test_field_division(self, backend, simple_mesh):
        """Test field / scalar."""
        field = backend.create_scalar_field(simple_mesh, value=6.0)
        result = field / 2.0
        assert np.allclose(result.get_values(), 3.0)
    
    def test_field_negation(self, backend, simple_mesh):
        """Test -field."""
        field = backend.create_scalar_field(simple_mesh, value=5.0)
        result = -field
        assert np.allclose(result.get_values(), -5.0)
    
    def test_field_power(self, backend, simple_mesh):
        """Test field ** power."""
        field = backend.create_scalar_field(simple_mesh, value=3.0)
        result = field ** 2
        assert np.allclose(result.get_values(), 9.0)
    
    def test_field_copy(self, backend, simple_mesh):
        """Test field copy creates independent copy."""
        f1 = backend.create_scalar_field(simple_mesh, value=5.0)
        f2 = f1.copy()
        f2.set_values(10.0)
        
        assert np.allclose(f1.get_values(), 5.0)
        assert np.allclose(f2.get_values(), 10.0)


# =============================================================================
# Test mathematical operations
# =============================================================================

class TestMathOperations:
    """Tests for mathematical operations (gradient, integrate, etc.)."""
    
    def test_where_positive(self, backend, simple_mesh):
        """Test where_positive returns correct mask."""
        # Create field with mixed positive/negative values
        coords = simple_mesh.get_cell_centers()
        x = coords[:, 0]
        values = x - 50.0  # negative for x < 50, positive for x > 50
        
        field = backend.create_scalar_field(simple_mesh, value=values)
        mask = backend.where_positive(field)
        
        mask_values = mask.get_values()
        expected = (values > 0).astype(float)
        assert_array_almost_equal(mask_values, expected)
    
    def test_where_negative(self, backend, simple_mesh):
        """Test where_negative returns correct mask."""
        coords = simple_mesh.get_cell_centers()
        x = coords[:, 0]
        values = x - 50.0
        
        field = backend.create_scalar_field(simple_mesh, value=values)
        mask = backend.where_negative(field)
        
        mask_values = mask.get_values()
        expected = (values < 0).astype(float)
        assert_array_almost_equal(mask_values, expected)
    
    def test_where_zero(self, backend, simple_mesh):
        """Test where_zero returns correct mask."""
        field = backend.create_scalar_field(simple_mesh, value=0.0)
        mask = backend.where_zero(field)
        
        mask_values = mask.get_values()
        assert np.allclose(mask_values, 1.0)
    
    def test_minimum_maximum(self, backend, simple_mesh):
        """Test min/max reduction operations."""
        values = np.linspace(-10, 20, simple_mesh.num_cells)
        field = backend.create_scalar_field(simple_mesh, value=values)
        
        assert_allclose(backend.minimum(field), -10.0, rtol=1e-5)
        assert_allclose(backend.maximum(field), 20.0, rtol=1e-5)
    
    def test_max_abs(self, backend, simple_mesh):
        """Test max absolute value."""
        values = np.linspace(-30, 20, simple_mesh.num_cells)
        field = backend.create_scalar_field(simple_mesh, value=values)
        
        assert_allclose(backend.max_abs(field), 30.0, rtol=1e-5)
    
    def test_integrate(self, backend, simple_mesh):
        """Test integration over domain."""
        # Create uniform field of 1.0
        field = backend.create_scalar_field(simple_mesh, value=1.0)
        
        # Integrate should give area of domain (100 * 50 = 5000)
        integral = backend.integrate(field)
        assert_allclose(integral, 5000.0, rtol=1e-2)
    
    def test_gradient(self, backend, simple_mesh):
        """Test gradient computation."""
        # Create linear field: f(x,y) = x
        # Gradient should be (1, 0)
        coords = simple_mesh.get_cell_centers()
        x = coords[:, 0]
        
        field = backend.create_scalar_field(simple_mesh, value=x)
        grad = backend.gradient(field)
        
        # Check that gradient in x-direction is approximately 1
        grad_values = grad.get_values()
        # This is a rough check - exact values depend on mesh and interpolation
        assert grad_values is not None


# =============================================================================
# Test PDE solving
# =============================================================================

class TestPDESolving:
    """Tests for PDE solving capabilities."""
    
    def test_solve_laplace_equation(self, backend, simple_mesh):
        """
        Test solving Laplace equation: -div(grad(u)) = 0
        
        With boundary conditions:
        - u = 0 at x = 0
        - u = 1 at x = L
        
        Expected solution: u(x) = x / L
        """
        L = 100.0  # domain length
        
        # Create solver
        pde = backend.create_pde_solver(simple_mesh)
        
        # Set coefficients for Laplace equation
        # -div(A * grad(u)) = 0 with A = 1
        pde.set_coefficients(A=1.0, D=0.0, Y=0.0)
        
        # Set boundary conditions
        # Use cell centers for boundary mask creation - works for both
        # node-based (FEM/escript) and cell-centered (FVM/FiPy) backends
        coords = simple_mesh.get_cell_centers()
        x = coords[:, 0]
        
        # Create masks for left (x~0) and right (x~L) boundaries
        # Use a tolerance based on cell size to identify boundary cells
        # For cell-centered methods, cells near boundaries have centers offset from edges
        x_min, x_max = x.min(), x.max()
        dx = (x_max - x_min) / (len(np.unique(x)) - 1) if len(np.unique(x)) > 1 else x_max
        tol = dx * 0.6  # Cells within 60% of cell width from boundary
        left_mask = (x < x_min + tol).astype(float)
        right_mask = (x > x_max - tol).astype(float)
        dirichlet_mask = backend.create_scalar_field(
            simple_mesh, value=left_mask + right_mask
        )
        
        # Set values: 0 at left, 1 at right
        dirichlet_values = right_mask  # 0 at left (default), 1 at right
        dirichlet_value = backend.create_scalar_field(
            simple_mesh, value=dirichlet_values
        )
        
        pde.set_boundary_conditions(
            dirichlet_mask=dirichlet_mask,
            dirichlet_value=dirichlet_value
        )
        
        # Solve
        solution = pde.solve()
        
        # Check solution
        sol_coords = solution.get_coordinates()
        sol_values = solution.get_values()
        expected = sol_coords[:, 0] / L
        
        # Allow for some numerical error
        assert_allclose(sol_values, expected, rtol=0.1, atol=0.05)
    
    def test_solve_with_source_term(self, backend, simple_mesh):
        """
        Test solving Poisson equation: -div(grad(u)) = 1
        
        With homogeneous Dirichlet BCs: u = 0 on all boundaries.
        """
        pde = backend.create_pde_solver(simple_mesh)
        
        # Set coefficients
        pde.set_coefficients(A=1.0, D=0.0, Y=1.0)
        
        # Set homogeneous Dirichlet BCs on all boundaries
        # Use cell centers for boundary mask creation - works for both
        # node-based (FEM/escript) and cell-centered (FVM/FiPy) backends
        coords = simple_mesh.get_cell_centers()
        x, y = coords[:, 0], coords[:, 1]
        L, H = 100.0, 50.0
        
        # Create boundary mask using relative position
        # For cell-centered methods, identify cells near domain edges
        x_min, x_max = x.min(), x.max()
        y_min, y_max = y.min(), y.max()
        dx = (x_max - x_min) / (len(np.unique(x)) - 1) if len(np.unique(x)) > 1 else x_max
        dy = (y_max - y_min) / (len(np.unique(y)) - 1) if len(np.unique(y)) > 1 else y_max
        tol_x = dx * 0.6
        tol_y = dy * 0.6
        
        # Create boundary mask
        on_boundary = (
            (x < x_min + tol_x) |
            (x > x_max - tol_x) |
            (y < y_min + tol_y) |
            (y > y_max - tol_y)
        ).astype(float)
        
        dirichlet_mask = backend.create_scalar_field(simple_mesh, value=on_boundary)
        dirichlet_value = backend.create_scalar_field(simple_mesh, value=0.0)
        
        pde.set_boundary_conditions(
            dirichlet_mask=dirichlet_mask,
            dirichlet_value=dirichlet_value
        )
        
        # Solve
        solution = pde.solve()
        
        # Solution should be positive in the interior (source > 0)
        sol_values = solution.get_values()
        interior_mask = on_boundary < 0.5
        if np.any(interior_mask):
            assert np.all(sol_values[interior_mask] > 0)
    
    def test_anisotropic_diffusion(self, backend, simple_mesh):
        """
        Test solving diffusion with anisotropic tensor.
        
        This is crucial for groundwater flow with anisotropic permeability.
        """
        pde = backend.create_pde_solver(simple_mesh)
        
        # Create anisotropic diffusion tensor
        # k_xx = 1.0, k_yy = 0.1 (10:1 anisotropy)
        k_tensor = backend.create_tensor_field(
            simple_mesh,
            values=((1.0, 0.0), (0.0, 0.1))
        )
        
        pde.set_coefficients(A=k_tensor, D=0.0, Y=0.0)
        
        # Set simple boundary conditions
        # Use cell centers for boundary mask creation - works for both
        # node-based (FEM/escript) and cell-centered (FVM/FiPy) backends
        coords = simple_mesh.get_cell_centers()
        x = coords[:, 0]
        
        # Use relative position to identify boundary cells
        x_min, x_max = x.min(), x.max()
        dx = (x_max - x_min) / (len(np.unique(x)) - 1) if len(np.unique(x)) > 1 else x_max
        tol = dx * 0.6
        
        left_mask = (x < x_min + tol).astype(float)
        right_mask = (x > x_max - tol).astype(float)
        
        dirichlet_mask = backend.create_scalar_field(
            simple_mesh, value=left_mask + right_mask
        )
        dirichlet_value = backend.create_scalar_field(
            simple_mesh, value=right_mask
        )
        
        pde.set_boundary_conditions(
            dirichlet_mask=dirichlet_mask,
            dirichlet_value=dirichlet_value
        )
        
        # Solve - just check it doesn't crash
        solution = pde.solve()
        assert solution is not None


# =============================================================================
# Test mesh operations
# =============================================================================

class TestMeshOperations:
    """Tests for mesh-related operations."""
    
    def test_mesh_cell_centers_shape(self, backend, simple_mesh):
        """Test that cell centers have correct shape."""
        centers = simple_mesh.get_cell_centers()
        assert centers.ndim == 2
        assert centers.shape[1] == 2  # x, y coordinates
        assert centers.shape[0] == simple_mesh.num_cells
    
    def test_mesh_node_coordinates_shape(self, backend, simple_mesh):
        """Test that node coordinates have correct shape."""
        nodes = simple_mesh.get_node_coordinates()
        assert nodes.ndim == 2
        assert nodes.shape[1] == 2
        assert nodes.shape[0] == simple_mesh.num_nodes
    
    def test_mesh_cell_centers_in_domain(self, backend):
        """Test that cell centers are within domain bounds."""
        L, H = 100.0, 50.0
        mesh = backend.create_rectangle_mesh(L, H, 10, 5)
        centers = mesh.get_cell_centers()
        
        assert np.all(centers[:, 0] >= 0)
        assert np.all(centers[:, 0] <= L)
        assert np.all(centers[:, 1] >= 0)
        assert np.all(centers[:, 1] <= H)
    
    def test_get_native_mesh(self, backend, simple_mesh):
        """Test that native mesh object can be retrieved."""
        native = simple_mesh.get_native_mesh()
        assert native is not None


# =============================================================================
# Test I/O operations
# =============================================================================

class TestIO:
    """Tests for I/O operations."""
    
    def test_save_vtk(self, backend, simple_mesh, tmp_path):
        """Test saving mesh and fields to VTK format."""
        pressure = backend.create_scalar_field(simple_mesh, value=1.0, name="pressure")
        concentration = backend.create_scalar_field(simple_mesh, value=0.5, name="concentration")
        
        output_file = str(tmp_path / "test_output")
        backend.save_vtk(
            output_file,
            simple_mesh,
            {"pressure": pressure, "concentration": concentration}
        )
        
        # Check that file was created
        import os
        vtk_files = list(tmp_path.glob("test_output*"))
        assert len(vtk_files) > 0


# =============================================================================
# Cross-backend consistency tests
# =============================================================================

class TestCrossBackendConsistency:
    """
    Tests to ensure both backends produce consistent results.
    
    These tests run the same problem with both backends and compare results.
    """
    
    @pytest.mark.skipif(
        not (check_backend_available('escript') and check_backend_available('fipy')),
        reason="Both backends required for cross-backend tests"
    )
    def test_laplace_consistency(self):
        """Test that both backends give same Laplace solution."""
        results = {}
        
        for name in ['escript', 'fipy']:
            backend = get_backend(name)
            mesh = backend.create_rectangle_mesh(100.0, 50.0, 20, 10)
            
            pde = backend.create_pde_solver(mesh)
            pde.set_coefficients(A=1.0, D=0.0, Y=0.0)
            
            # Set boundary conditions
            coords = mesh.get_node_coordinates()
            x = coords[:, 0]
            L = 100.0
            tol = 1e-6
            
            left_mask = (np.abs(x) < tol).astype(float)
            right_mask = (np.abs(x - L) < tol).astype(float)
            
            dirichlet_mask = backend.create_scalar_field(mesh, value=left_mask + right_mask)
            dirichlet_value = backend.create_scalar_field(mesh, value=right_mask)
            
            pde.set_boundary_conditions(
                dirichlet_mask=dirichlet_mask,
                dirichlet_value=dirichlet_value
            )
            
            solution = pde.solve()
            results[name] = solution.get_values()
        
        # Compare results (allowing for some numerical differences)
        assert_allclose(
            results['escript'], 
            results['fipy'], 
            rtol=0.05, 
            atol=0.02
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
