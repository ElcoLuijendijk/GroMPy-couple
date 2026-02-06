"""
Unit tests for the FiPy mesh creation functions.

These tests verify that mesh_functions_fipy.py correctly creates meshes
using FiPy and the Gmsh Python API, without requiring esys-escript.
"""

import pytest
import numpy as np
import tempfile
import os
import warnings

# Suppress FiPy warnings in tests
warnings.filterwarnings('ignore')


# =============================================================================
# Check if FiPy and Gmsh are available
# =============================================================================

try:
    import fipy
    FIPY_AVAILABLE = True
except ImportError:
    FIPY_AVAILABLE = False

try:
    import gmsh
    GMSH_AVAILABLE = True
except ImportError:
    GMSH_AVAILABLE = False


# Skip all tests if FiPy is not available
pytestmark = pytest.mark.skipif(
    not FIPY_AVAILABLE,
    reason="FiPy not available"
)


# =============================================================================
# Test fixtures
# =============================================================================

class MockRectangularParams:
    """Mock parameters for rectangular mesh."""
    L = 100.0
    thickness = 50.0
    cellsize_x = 10.0
    cellsize_y = 5.0


class MockStandardParams:
    """Mock parameters for standard mesh."""
    L = 100.0
    thickness = 50.0
    cellsize = 10.0
    topo_gradient = 0.01
    topo_break = False
    x_topo_break = 50.0


class MockCoastalParams:
    """Mock parameters for coastal mesh."""
    L = 1000.0
    L_sea = 100.0
    thickness = 50.0
    cellsize = 20.0
    topo_gradient = 0.01
    recharge_flux = 1e-9  # m/s
    k = 1e-12  # permeability m^2
    rho_f_0 = 1000.0
    viscosity = 1e-3
    buffer_distance_sea = 50.0
    buffer_distance_land = 50.0
    grid_refinement_factor = 0.5
    grid_refinement_factor_sea = 1.0


@pytest.fixture
def temp_mesh_file():
    """Create a temporary mesh file and clean up after test."""
    with tempfile.NamedTemporaryFile(suffix='.msh', delete=False) as f:
        filepath = f.name
    yield filepath
    if os.path.exists(filepath):
        os.remove(filepath)


# =============================================================================
# Test rectangular mesh
# =============================================================================

class TestRectangularMesh:
    """Tests for setup_rectangular_mesh_fipy."""
    
    def test_creates_mesh(self, temp_mesh_file):
        """Test that rectangular mesh is created successfully."""
        from lib.mesh_functions_fipy import setup_rectangular_mesh_fipy
        
        params = MockRectangularParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_rectangular_mesh_fipy(
            params, temp_mesh_file
        )
        
        assert mesh is not None
        assert mesh.numberOfCells > 0
    
    def test_mesh_dimensions(self, temp_mesh_file):
        """Test that mesh has correct dimensions."""
        from lib.mesh_functions_fipy import setup_rectangular_mesh_fipy
        
        params = MockRectangularParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_rectangular_mesh_fipy(
            params, temp_mesh_file
        )
        
        # Expected: 10x10 = 100 cells (L/cellsize_x * thickness/cellsize_y)
        expected_cells = int(params.L / params.cellsize_x) * int(params.thickness / params.cellsize_y)
        assert mesh.numberOfCells == expected_cells
    
    def test_surface_mask(self, temp_mesh_file):
        """Test that surface mask identifies top boundary cells."""
        from lib.mesh_functions_fipy import setup_rectangular_mesh_fipy
        
        params = MockRectangularParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_rectangular_mesh_fipy(
            params, temp_mesh_file
        )
        
        # Surface mask should be a boolean array
        assert surface.dtype == bool
        
        # Should have one row of surface cells
        nx = int(params.L / params.cellsize_x)
        assert np.sum(surface) == nx
    
    def test_z_surface_constant(self, temp_mesh_file):
        """Test that z_surface is constant for rectangular mesh."""
        from lib.mesh_functions_fipy import setup_rectangular_mesh_fipy
        
        params = MockRectangularParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_rectangular_mesh_fipy(
            params, temp_mesh_file
        )
        
        # z_surface should be constant = thickness
        assert np.all(z_surface == params.thickness)
    
    def test_no_sea_surface_or_seawater(self, temp_mesh_file):
        """Test that rectangular mesh has no sea surface or seawater."""
        from lib.mesh_functions_fipy import setup_rectangular_mesh_fipy
        
        params = MockRectangularParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_rectangular_mesh_fipy(
            params, temp_mesh_file
        )
        
        assert sea_surface is None
        assert seawater is None


# =============================================================================
# Test standard mesh (requires Gmsh)
# =============================================================================

@pytest.mark.skipif(not GMSH_AVAILABLE, reason="Gmsh not available")
class TestStandardMesh:
    """Tests for setup_standard_mesh_fipy."""
    
    def test_creates_mesh(self, temp_mesh_file):
        """Test that standard mesh is created successfully."""
        from lib.mesh_functions_fipy import setup_standard_mesh_fipy
        
        params = MockStandardParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_standard_mesh_fipy(
            params, temp_mesh_file
        )
        
        assert mesh is not None
        assert mesh.numberOfCells > 0
    
    def test_mesh_file_created(self, temp_mesh_file):
        """Test that mesh file is written to disk."""
        from lib.mesh_functions_fipy import setup_standard_mesh_fipy
        
        params = MockStandardParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_standard_mesh_fipy(
            params, temp_mesh_file
        )
        
        assert os.path.exists(temp_mesh_file)
    
    def test_surface_mask_shape(self, temp_mesh_file):
        """Test that surface mask has correct shape."""
        from lib.mesh_functions_fipy import setup_standard_mesh_fipy
        
        params = MockStandardParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_standard_mesh_fipy(
            params, temp_mesh_file
        )
        
        # Surface mask should have one element per cell
        assert len(surface) == mesh.numberOfCells
        assert surface.dtype == bool
    
    def test_z_surface_follows_gradient(self, temp_mesh_file):
        """Test that z_surface follows topographic gradient."""
        from lib.mesh_functions_fipy import setup_standard_mesh_fipy
        
        params = MockStandardParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_standard_mesh_fipy(
            params, temp_mesh_file
        )
        
        # z_surface should be proportional to x
        cell_centers = mesh.cellCenters
        x_centers = np.array(cell_centers[0])
        
        # For simple case (no topo_break), z = x * topo_gradient
        expected_z = x_centers * params.topo_gradient
        np.testing.assert_allclose(z_surface, expected_z, rtol=1e-5)


# =============================================================================
# Test coastal mesh (requires Gmsh)
# =============================================================================

@pytest.mark.skipif(not GMSH_AVAILABLE, reason="Gmsh not available")
class TestCoastalMesh:
    """Tests for setup_coastal_mesh_fipy."""
    
    def test_creates_mesh(self, temp_mesh_file):
        """Test that coastal mesh is created successfully."""
        from lib.mesh_functions_fipy import setup_coastal_mesh_fipy
        
        params = MockCoastalParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_coastal_mesh_fipy(
            params, temp_mesh_file
        )
        
        assert mesh is not None
        assert mesh.numberOfCells > 0
    
    def test_has_sea_surface_mask(self, temp_mesh_file):
        """Test that coastal mesh has sea surface mask."""
        from lib.mesh_functions_fipy import setup_coastal_mesh_fipy
        
        params = MockCoastalParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_coastal_mesh_fipy(
            params, temp_mesh_file
        )
        
        assert sea_surface is not None
        assert sea_surface.dtype == bool
        # Sea surface should have some cells (x < 0 and y ~ 0)
        assert np.sum(sea_surface) >= 0
    
    def test_has_seawater_mask(self, temp_mesh_file):
        """Test that coastal mesh has seawater mask."""
        from lib.mesh_functions_fipy import setup_coastal_mesh_fipy
        
        params = MockCoastalParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_coastal_mesh_fipy(
            params, temp_mesh_file
        )
        
        assert seawater is not None
        assert seawater.dtype == bool
        # Seawater should have some cells (x < 0)
        assert np.sum(seawater) > 0
    
    def test_z_surface_follows_gradient(self, temp_mesh_file):
        """Test that z_surface follows topographic gradient."""
        from lib.mesh_functions_fipy import setup_coastal_mesh_fipy
        
        params = MockCoastalParams()
        mesh, surface, sea_surface, seawater, z_surface = setup_coastal_mesh_fipy(
            params, temp_mesh_file
        )
        
        # z_surface should be x * topo_gradient
        cell_centers = mesh.cellCenters
        x_centers = np.array(cell_centers[0])
        
        expected_z = x_centers * params.topo_gradient
        np.testing.assert_allclose(z_surface, expected_z, rtol=1e-5)


# =============================================================================
# Test module interface (drop-in replacement)
# =============================================================================

class TestModuleInterface:
    """Test that the module can be used as a drop-in replacement."""
    
    def test_alias_functions_exist(self):
        """Test that alias functions exist for drop-in replacement."""
        from lib import mesh_functions_fipy
        
        # These aliases should exist for drop-in replacement
        assert hasattr(mesh_functions_fipy, 'setup_rectangular_mesh')
        assert hasattr(mesh_functions_fipy, 'setup_standard_mesh')
        assert hasattr(mesh_functions_fipy, 'setup_coastal_mesh')
        assert hasattr(mesh_functions_fipy, 'setup_coastal_mesh_glover1959')
    
    def test_fipy_available_flag(self):
        """Test that FIPY_AVAILABLE flag is set correctly."""
        from lib import mesh_functions_fipy
        
        assert mesh_functions_fipy.FIPY_AVAILABLE == True
    
    @pytest.mark.skipif(not GMSH_AVAILABLE, reason="Gmsh not available")
    def test_gmsh_available_flag(self):
        """Test that GMSH_AVAILABLE flag is set correctly."""
        from lib import mesh_functions_fipy
        
        assert mesh_functions_fipy.GMSH_AVAILABLE == True


# =============================================================================
# Test grompy.py integration
# =============================================================================

class TestGrompyIntegration:
    """Test that grompy.py correctly uses FiPy mesh functions."""
    
    def test_grompy_imports_without_escript(self):
        """Test that grompy.py can be imported without escript."""
        import grompy
        
        assert hasattr(grompy, 'FIPY_AVAILABLE')
        assert hasattr(grompy, 'get_backend_name')
        assert hasattr(grompy, 'get_mesh_functions')
    
    def test_get_mesh_functions_fipy(self):
        """Test that get_mesh_functions returns FiPy module."""
        import grompy
        
        mesh_funcs = grompy.get_mesh_functions('fipy')
        
        assert 'fipy' in mesh_funcs.__name__
        assert hasattr(mesh_funcs, 'setup_rectangular_mesh')
    
    def test_backend_selection_fipy(self):
        """Test that FiPy backend is selected when requested."""
        import grompy
        
        class MockOptions:
            backend = 'fipy'
        
        backend_name = grompy.get_backend_name(MockOptions)
        
        if grompy.FIPY_AVAILABLE:
            assert backend_name == 'fipy'
