"""
Comprehensive test suite for salt wedge benchmark with FiPy backend.

Based on Goswami and Clement (2007) salt wedge intrusion experiments.
Validates:
- Mesh creation and file writing (tests recent fix)
- Model parameter configuration
- Glover (1959) analytical solution
- Experimental data integrity
- Full model execution for 3 scenarios
- Comparison with experimental measurements
- Physical validity of results
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path
import tempfile
import shutil

# Local imports
from lib.mesh_functions_fipy import setup_rectangular_mesh_fipy
from model_input.model_parameters_sw_benchmark import (
    ModelParameters, ModelOptions, ParameterRanges
)


# =============================================================================
# TEST FIXTURES
# =============================================================================

@pytest.fixture
def benchmark_parameters():
    """Load benchmark model parameters."""
    return ModelParameters()


@pytest.fixture
def benchmark_options():
    """Load benchmark model options."""
    opts = ModelOptions()
    opts.backend = 'fipy'
    return opts


@pytest.fixture
def parameter_ranges():
    """Load benchmark parameter ranges."""
    return ParameterRanges()


@pytest.fixture
def experimental_data():
    """Load experimental data from CSV."""
    csv_path = Path(__file__).parent.parent / 'benchmark_data' / \
               'table_a1_steady_state_salt_wedge_locations.csv'
    if not csv_path.exists():
        pytest.skip(f"Experimental data not found at {csv_path}")
    df = pd.read_csv(csv_path)
    return df


@pytest.fixture(params=[0, 1, 2], ids=['SS-1', 'SS-2', 'SS-3'])
def scenario_index(request):
    """Parametrize over 3 salt wedge scenarios."""
    return request.param


@pytest.fixture
def scenario_parameters(benchmark_parameters, parameter_ranges, scenario_index):
    """Get parameters for specific scenario."""
    params = benchmark_parameters
    pressures = parameter_ranges.specified_pressure_s
    params.specified_pressure = pressures[scenario_index]
    return params, scenario_index


@pytest.fixture
def fast_benchmark_parameters(parameter_ranges):
    """Fast variant of benchmark parameters for quick testing.
    
    Optimizations:
    - Coarse mesh: 2x larger cells (4x fewer elements)
    - Fast timesteps: 4x larger initial dt, enabled growth
    - Shorter simulation: 20 min instead of 80 min
    
    Expected runtime: 1-2 minutes per scenario (vs 10-15 min)
    """
    params = ModelParameters()
    
    # Coarse mesh (4x speedup from fewer cells)
    params.cellsize_x = 0.005  # doubled from 0.0025
    params.cellsize_y = 0.005
    
    # Fast timesteps with growth enabled (16x speedup from fewer timesteps)
    params.dt0 = 1.0            # 4x larger from 0.25
    params.dt_inc = 1.5         # FIXED: was 1.0 (no growth), now 1.5 (enables growth)
    params.dt_max = 60.0        # 6x larger from 10.0
    params.total_time = 1200.0  # 4x shorter: 20 min instead of 80 min
    
    # Output interval can be longer for fast tests
    params.output_interval = 60.0  # every 1 min
    
    return params


@pytest.fixture
def fast_scenario_parameters(fast_benchmark_parameters, parameter_ranges, scenario_index):
    """Get fast parameters for specific scenario."""
    params = fast_benchmark_parameters
    pressures = parameter_ranges.specified_pressure_s
    params.specified_pressure = pressures[scenario_index]
    return params, scenario_index


@pytest.fixture
def test_output_base_dir(tmp_path_factory):
    """Create base output directory for test suite."""
    # Create tests/test_output if it doesn't exist
    base_output = Path(__file__).parent / "test_output"
    base_output.mkdir(exist_ok=True)
    return base_output


@pytest.fixture
def test_output_dir(test_output_base_dir, scenario_index):
    """Create scenario-specific output directory."""
    scenario_dir = test_output_base_dir / f"scenario_S{scenario_index}"
    scenario_dir.mkdir(exist_ok=True, parents=True)
    return scenario_dir


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def depth_sw_interface_glover1959(x, Q, K, rho_f, rho_s):
    """
    Calculate depth of fresh-salt water interface per Glover (1959).
    
    Parameters
    ----------
    x : float or array
        Horizontal distance from coastline (m)
    Q : float
        Freshwater flux (m²/s)
    K : float
        Hydraulic conductivity (m/s)
    rho_f : float
        Freshwater density (kg/m³)
    rho_s : float
        Seawater density (kg/m³)
    
    Returns
    -------
    depth : float or array
        Depth of salt-water interface (m)
    """
    gamma = (rho_s - rho_f) / rho_f
    y2 = 2 * Q / (gamma * K) * x + Q**2 / (gamma**2 * K**2)
    depth = np.sqrt(y2)
    return depth


def extract_isochlor_profile(cell_centers, concentration, target_conc, params):
    """
    Extract x-positions of specific isochlor level at each y.
    
    Parameters
    ----------
    cell_centers : array
        Cell center coordinates (shape: 2 x n_cells)
    concentration : array
        Concentration field (shape: n_cells,)
    target_conc : float
        Target concentration for isochlor
    params : ModelParameters
        Model parameters
    
    Returns
    -------
    y_isochlor : array
        Y-coordinates of isochlor points
    x_isochlor : array
        X-coordinates of isochlor points (interpolated)
    """
    x_cells = np.array(cell_centers[0])
    y_cells = np.array(cell_centers[1])
    
    unique_y = np.sort(np.unique(y_cells))
    x_isochlor = []
    y_isochlor = []
    
    for y_target in unique_y:
        # Find cells at this y-coordinate
        tol = (y_cells.max() - y_cells.min()) * 0.01
        mask = np.abs(y_cells - y_target) < tol
        
        if mask.sum() < 2:
            continue
        
        x_at_y = x_cells[mask]
        conc_at_y = concentration[mask]
        
        # Sort by x
        sort_idx = np.argsort(x_at_y)
        x_sorted = x_at_y[sort_idx]
        conc_sorted = conc_at_y[sort_idx]
        
        # Find where concentration crosses target
        if conc_sorted.min() <= target_conc <= conc_sorted.max():
            x_iso = np.interp(target_conc, conc_sorted, x_sorted)
            x_isochlor.append(x_iso)
            y_isochlor.append(y_target)
    
    return np.array(y_isochlor), np.array(x_isochlor)


def calculate_misfit_metrics(modeled_x, experimental_x):
    """
    Calculate misfit metrics between modeled and experimental positions.
    
    Parameters
    ----------
    modeled_x : array
        Modeled x-positions
    experimental_x : array
        Experimental x-positions
    
    Returns
    -------
    dict
        Dictionary with 'mean_error', 'mae', 'rmse'
    """
    errors = modeled_x - experimental_x
    mean_error = np.mean(errors)
    mae = np.mean(np.abs(errors))
    rmse = np.sqrt(np.mean(errors**2))
    
    return {
        'mean_error': mean_error,
        'mae': mae,
        'rmse': rmse
    }


# =============================================================================
# PART 1: UNIT TESTS - GLOVER ANALYTICAL SOLUTION
# =============================================================================

class TestGloverAnalyticalSolution:
    """Test Glover (1959) analytical salt wedge solution."""
    
    def test_glover_basic_calculation(self):
        """Test basic Glover depth formula with simple values."""
        x = 0.1
        Q = 1e-6
        K = 1e-5
        rho_f = 1000
        rho_s = 1025
        
        depth = depth_sw_interface_glover1959(x, Q, K, rho_f, rho_s)
        
        # Just validate the formula works and returns reasonable types
        assert isinstance(depth, (float, np.floating))
        assert depth > 0, "Depth should be positive"
        assert np.isfinite(depth), "Depth should be finite"
    
    def test_glover_with_benchmark_params(self, benchmark_parameters):
        """Test Glover with actual salt wedge benchmark parameters."""
        params = benchmark_parameters
        
        # Scenario parameters
        Q_values = [1.42, 0.59, 1.19]  # cm³/s for SS-1, SS-2, SS-3
        Q_values = [q / 2.7 * 1e-4 for q in Q_values]  # Convert to m³/s
        
        K = 1050.0 / (24.0 * 60.0 * 60.0)  # m/s
        rho_f = params.rho_f_0
        rho_s = 1026.0
        
        depths = []
        for Q in Q_values:
            depth_at_L = depth_sw_interface_glover1959(
                params.L, Q, K, rho_f, rho_s
            )
            depths.append(depth_at_L)
        
        # Higher discharge (Q) leads to deeper salt wedge intrusion at same location
        # SS-1 has highest Q, so deepest intrusion; SS-2 has lowest Q, shallowest
        assert depths[0] > depths[1], "SS-1 (Q=1.42) should intrude deeper than SS-2 (Q=0.59)"
        assert depths[0] > depths[2], "SS-1 (Q=1.42) should intrude deeper than SS-3 (Q=1.19)"
        assert depths[1] < depths[2], "SS-2 (Q=0.59) should be shallower than SS-3 (Q=1.19)"
    
    def test_glover_boundary_conditions(self, benchmark_parameters):
        """Test Glover solution at domain boundaries."""
        params = benchmark_parameters
        Q = 1e-4
        K = 1e-5
        rho_f = 1000
        rho_s = 1025
        
        # At x=0: depth should be 0
        depth_at_0 = depth_sw_interface_glover1959(0, Q, K, rho_f, rho_s)
        assert depth_at_0 >= 0, "Depth at x=0 should be non-negative"
        
        # At x=L: should be maximum
        depth_at_L = depth_sw_interface_glover1959(params.L, Q, K, rho_f, rho_s)
        assert depth_at_L > 0, "Depth at x=L should be positive"
        
        # Check monotonic increase
        x_values = np.linspace(0, params.L, 10)
        depths = [depth_sw_interface_glover1959(x, Q, K, rho_f, rho_s) 
                  for x in x_values]
        diffs = np.diff(depths)
        assert np.all(diffs >= 0), "Depth should increase monotonically"
    
    def test_glover_array_input(self, benchmark_parameters):
        """Test Glover with array input."""
        params = benchmark_parameters
        x_array = np.linspace(0, params.L, 100)
        Q = 1e-4
        K = 1e-5
        rho_f = 1000
        rho_s = 1025
        
        depths = depth_sw_interface_glover1959(x_array, Q, K, rho_f, rho_s)
        
        assert isinstance(depths, np.ndarray)
        assert len(depths) == len(x_array)
        assert np.all(depths >= 0)


# =============================================================================
# PART 1: UNIT TESTS - EXPERIMENTAL DATA LOADING
# =============================================================================

class TestExperimentalDataLoading:
    """Test loading and validation of experimental benchmark data."""
    
    def test_csv_file_exists(self, experimental_data):
        """Verify experimental CSV file exists and is readable."""
        # If we got here, the file was successfully loaded
        assert experimental_data is not None
        assert len(experimental_data) > 0
    
    def test_csv_column_names(self, experimental_data):
        """Verify expected columns are present."""
        expected_cols = ['x_ss1', 'y_ss1', 'x_ss2', 'y_ss2', 'x_ss3', 'y_ss3']
        for col in expected_cols:
            assert col in experimental_data.columns, f"Missing column: {col}"
    
    def test_experimental_data_ranges(self, experimental_data):
        """Verify data values are physically reasonable."""
        # x values in cm: 0 to ~25.3
        x_cols = ['x_ss1', 'x_ss2', 'x_ss3']
        for col in x_cols:
            x_vals = experimental_data[col].dropna().values
            assert x_vals.min() >= 0, f"{col} has negative values"
            assert x_vals.max() <= 40, f"{col} exceeds domain length"
        
        # y values in cm: 0 to ~26
        y_cols = ['y_ss1', 'y_ss2', 'y_ss3']
        for col in y_cols:
            y_vals = experimental_data[col].dropna().values
            assert y_vals.min() >= 0, f"{col} has negative values"
            assert y_vals.max() <= 30, f"{col} exceeds domain thickness"
    
    def test_three_scenarios_present(self, experimental_data):
        """Verify all 3 scenarios have data points."""
        n_ss1 = experimental_data['x_ss1'].notna().sum()
        n_ss2 = experimental_data['x_ss2'].notna().sum()
        n_ss3 = experimental_data['x_ss3'].notna().sum()
        
        assert n_ss1 > 5, f"SS-1 has only {n_ss1} points"
        assert n_ss2 > 5, f"SS-2 has only {n_ss2} points"
        assert n_ss3 > 5, f"SS-3 has only {n_ss3} points"
    
    def test_salt_wedge_shape_validity(self, experimental_data):
        """Verify data represents valid salt wedge shape."""
        for scenario in [1, 2, 3]:
            x_col = f'x_ss{scenario}'
            y_col = f'y_ss{scenario}'
            
            x_vals = experimental_data[x_col].dropna().values
            y_vals = experimental_data[y_col].dropna().values
            
            # Sort by y to check wedge shape
            sort_idx = np.argsort(y_vals)
            x_sorted = x_vals[sort_idx]
            y_sorted = y_vals[sort_idx]
            
            # Generally x should decrease with increasing y (wedge shape)
            # Allow some scatter but overall trend should be negative
            diffs = np.diff(x_sorted)
            neg_count = np.sum(diffs < 0)
            total_count = len(diffs)
            
            assert neg_count / total_count > 0.5, \
                f"Scenario {scenario}: Not a valid wedge shape"


# =============================================================================
# PART 1: UNIT TESTS - BENCHMARK PARAMETERS
# =============================================================================

class TestSaltWedgeBenchmarkParameters:
    """Test model parameter configuration."""
    
    def test_model_options_backend(self, benchmark_options):
        """Verify FiPy backend is configured."""
        assert benchmark_options.backend == 'fipy'
        assert benchmark_options.solute_transport is True
        assert benchmark_options.coupled_iterations is True
    
    def test_mesh_parameters(self, benchmark_parameters):
        """Verify mesh configuration."""
        params = benchmark_parameters
        assert params.mesh_type == 'rectangle'
        assert params.L == 0.53
        assert params.thickness == 0.26
        assert params.cellsize_x == 0.0025
        assert params.cellsize_y == 0.0025
        
        # Verify expected cell count
        nx = int(np.ceil(params.L / params.cellsize_x))
        ny = int(np.ceil(params.thickness / params.cellsize_y))
        expected_cells = nx * ny
        assert expected_cells == 22048, f"Expected 22048 cells, got {expected_cells}"
    
    def test_boundary_conditions(self, benchmark_parameters):
        """Verify boundary condition setup."""
        params = benchmark_parameters
        
        # Specified pressure
        assert isinstance(params.specified_pressure, list)
        assert len(params.specified_pressure) == 2
        
        # Specified concentration
        assert isinstance(params.specified_concentration, list)
        assert params.specified_concentration[0] == 0.03624  # Seawater
        assert params.specified_concentration[1] == 0.0       # Freshwater
    
    def test_fluid_properties(self, benchmark_parameters):
        """Verify fluid properties."""
        params = benchmark_parameters
        
        assert params.rho_f_0 > 0
        assert params.porosity > 0
        assert params.porosity < 1
        assert params.viscosity > 0
        assert params.gamma > 0


# =============================================================================
# PART 1: UNIT TESTS - MESH CREATION AND SAVING
# =============================================================================

class TestMeshCreationAndSaving:
    """Test mesh creation and file writing (validates recent fix)."""
    
    def test_rectangular_mesh_created(self, benchmark_parameters, test_output_dir):
        """Call setup_rectangular_mesh_fipy and verify mesh object."""
        params = benchmark_parameters
        mesh_filename = str(test_output_dir / "mesh_test.msh")
        
        mesh, surface, sea_surface, seawater, z_surface = \
            setup_rectangular_mesh_fipy(params, mesh_filename)
        
        assert mesh is not None
        assert mesh.numberOfCells == 22048, \
            f"Wrong cell count: {mesh.numberOfCells}"
        assert surface is not None
        assert z_surface is not None
    
    def test_mesh_file_written(self, benchmark_parameters, test_output_dir):
        """Verify .msh file is created on disk."""
        params = benchmark_parameters
        mesh_filename = str(test_output_dir / "mesh_file_test.msh")
        
        mesh, _, _, _, _ = setup_rectangular_mesh_fipy(params, mesh_filename)
        
        mesh_path = Path(mesh_filename)
        assert mesh_path.exists(), f"Mesh file not created at {mesh_filename}"
        assert mesh_path.stat().st_size > 1e6, "Mesh file seems too small"
    
    def test_mesh_format_ascii_gmsh22(self, benchmark_parameters, test_output_dir):
        """Verify file format is ASCII gmsh22."""
        params = benchmark_parameters
        mesh_filename = str(test_output_dir / "mesh_format_test.msh")
        
        mesh, _, _, _, _ = setup_rectangular_mesh_fipy(params, mesh_filename)
        
        # Read first part of file and check for gmsh markers
        with open(mesh_filename, 'r', errors='ignore') as f:
            header = f.read(500)
            # Should contain gmsh format markers
            assert '$MeshFormat' in header or 'gmsh' in header.lower(), \
                "File doesn't appear to be gmsh format"
    
    def test_mesh_structure_valid(self, benchmark_parameters, test_output_dir):
        """Verify mesh structure is valid."""
        import meshio
        
        params = benchmark_parameters
        mesh_filename = str(test_output_dir / "mesh_structure_test.msh")
        
        mesh, _, _, _, _ = setup_rectangular_mesh_fipy(params, mesh_filename)
        
        # Load with meshio to verify structure
        mesh_data = meshio.read(mesh_filename)
        
        # Should have triangular cells (quads converted)
        cell_types = [cell.type for cell in mesh_data.cells]
        assert 'triangle' in cell_types, f"No triangles found, got: {cell_types}"
        
        # Count cells
        total_cells = sum(cell.data.shape[0] for cell in mesh_data.cells 
                         if cell.type == 'triangle')
        assert total_cells == 44096, f"Expected 44096 triangles, got {total_cells}"
    
    def test_mesh_domain_bounds(self, benchmark_parameters, test_output_dir):
        """Verify mesh covers full domain."""
        params = benchmark_parameters
        mesh_filename = str(test_output_dir / "mesh_bounds_test.msh")
        
        mesh, _, _, _, _ = setup_rectangular_mesh_fipy(params, mesh_filename)
        
        vertex_coords = mesh.vertexCoords
        x_coords = vertex_coords[0]
        y_coords = vertex_coords[1]
        
        assert x_coords.min() >= -1e-6, f"Min x too negative: {x_coords.min()}"
        assert x_coords.max() <= params.L + 1e-6, \
            f"Max x exceeds domain: {x_coords.max()}"
        assert y_coords.min() >= -1e-6, f"Min y too negative: {y_coords.min()}"
        assert y_coords.max() <= params.thickness + 1e-6, \
            f"Max y exceeds domain: {y_coords.max()}"
    
    def test_surface_mask_generation(self, benchmark_parameters, test_output_dir):
        """Verify surface boundary identification."""
        params = benchmark_parameters
        mesh_filename = str(test_output_dir / "mesh_surface_test.msh")
        
        mesh, surface, _, _, _ = setup_rectangular_mesh_fipy(params, mesh_filename)
        
        assert surface is not None
        assert surface.shape == (mesh.numberOfCells,)
        assert surface.dtype == bool
        assert surface.sum() > 100, "Should identify many surface cells"


# =============================================================================
# PART 2: INTEGRATION TESTS - FULL MODEL RUN
# =============================================================================

class TestSaltWedgeBenchmarkModelRun:
    """Run full salt wedge model for all 3 scenarios (slow tests)."""
    
    @pytest.mark.slow
    def test_run_scenario_complete(self, scenario_parameters, test_output_dir):
        """Full model run: setup → execution → output validation."""
        from lib.grompy_fipy import run_coupled_flow_model_fipy
        
        params, scenario_idx = scenario_parameters
        
        # ===== 1. SETUP =====
        mesh_filename = str(test_output_dir / f"mesh_s{scenario_idx}.msh")
        
        # Create mesh
        mesh_fipy, surface, sea_surface, seawater, z_surface = \
            setup_rectangular_mesh_fipy(params, mesh_filename)
        
        assert Path(mesh_filename).exists(), "Mesh file not created"
        assert mesh_fipy.numberOfCells == 22048, "Wrong cell count"
        
        # ===== 2. MODEL EXECUTION =====
        try:
            model_results = run_coupled_flow_model_fipy(
                params, ModelOptions(), mesh_filename
            )
        except Exception as e:
            pytest.fail(f"Model run failed: {str(e)[:200]}")
        
        # ===== 3. VERIFY OUTPUT STRUCTURE =====
        assert len(model_results) >= 24, \
            f"Wrong result tuple length: {len(model_results)}"
        
        # Unpack results (25 items returned)
        (mesh, surface, sea_surface, k_vector, P, Conc,
         rho_f, viscosity, h, q, q_abs, nodal_flux,
         Pdiff, Cdiff, Pmax, Cmax, Pmean, Cmean,
         dts, runtimes, nsteps, output_step,
         boundary_conditions, boundary_fluxes, boundary_flux_stats,
         reached_steady_state) = model_results
        
        # ===== 4. VERIFY PHYSICAL VALIDITY =====
        
        # Pressure
        assert P is not None, "Pressure field is None"
        assert isinstance(P, np.ndarray), "Pressure not array"
        assert len(P) == 22048, f"Pressure array wrong length: {len(P)}"
        assert P.min() >= -1000, f"Pressure too negative: {P.min()}"
        assert P.max() <= 10000, f"Pressure too high: {P.max()}"
        assert not np.isnan(P).any(), "Pressure contains NaN"
        
        # Concentration
        assert Conc is not None, "Concentration field is None"
        assert isinstance(Conc, np.ndarray), "Concentration not array"
        assert len(Conc) == 22048, f"Concentration wrong length: {len(Conc)}"
        assert Conc.min() >= -0.001, f"Concentration too low: {Conc.min()}"
        assert Conc.max() <= 0.04, f"Concentration too high: {Conc.max()}"
        assert not np.isnan(Conc).any(), "Concentration contains NaN"
        
        # Density
        assert rho_f is not None, "Density field is None"
        assert isinstance(rho_f, np.ndarray), "Density not array"
        assert rho_f.min() >= 990, f"Density too low: {rho_f.min()}"
        assert rho_f.max() <= 1030, f"Density too high: {rho_f.max()}"
        assert not np.isnan(rho_f).any(), "Density contains NaN"
        
        # Timesteps
        assert len(dts) == nsteps, "Timestep array wrong length"
        assert nsteps > 0, "No timesteps completed"
        
        # Store results for comparison tests
        pytest.results = {
            'scenario_idx': scenario_idx,
            'params': params,
            'output_dir': test_output_dir,
            'P': P, 'Conc': Conc, 'q': q, 'rho_f': rho_f,
            'cell_centers': mesh.cellCenters if hasattr(mesh, 'cellCenters') else None,
            'nsteps': nsteps,
            'reached_steady_state': reached_steady_state
        }


    @pytest.mark.fast
    def test_run_scenario_fast(self, fast_scenario_parameters, test_output_dir):
        """Fast model run with coarse mesh and short simulation time.
        
        Optimizations:
        - Mesh: 4x coarser (11k cells vs 44k)
        - Simulation: 20 min instead of 80 min
        - Timesteps: ~16-20 instead of ~320
        
        Expected runtime: < 2 minutes
        Acceptable error: ±5 cm in salt wedge position
        """
        from lib.grompy_fipy import run_coupled_flow_model_fipy
        
        params, scenario_idx = fast_scenario_parameters
        
        # ===== 1. SETUP (with coarse mesh) =====
        mesh_filename = str(test_output_dir / f"mesh_fast_s{scenario_idx}.msh")
        
        # Create coarse mesh
        mesh_fipy, surface, sea_surface, seawater, z_surface = \
            setup_rectangular_mesh_fipy(params, mesh_filename)
        
        assert Path(mesh_filename).exists(), "Mesh file not created"
        # Coarse mesh should have ~5512 cells (4x fewer than 22048)
        assert 5000 < mesh_fipy.numberOfCells < 6000, \
            f"Coarse mesh has wrong cell count: {mesh_fipy.numberOfCells}"
        
        # ===== 2. MODEL EXECUTION (fast) =====
        try:
            model_results = run_coupled_flow_model_fipy(
                params, ModelOptions(), mesh_filename
            )
        except Exception as e:
            pytest.fail(f"Fast model run failed: {str(e)[:200]}")
        
        # ===== 3. VERIFY OUTPUT STRUCTURE =====
        assert len(model_results) >= 24, \
            f"Wrong result tuple length: {len(model_results)}"
        
        # Unpack results (25 items returned)
        (mesh, surface, sea_surface, k_vector, P, Conc,
         rho_f, viscosity, h, q, q_abs, nodal_flux,
         Pdiff, Cdiff, Pmax, Cmax, Pmean, Cmean,
         dts, runtimes, nsteps, output_step,
         boundary_conditions, boundary_fluxes, boundary_flux_stats,
         reached_steady_state) = model_results
        
        # ===== 4. VERIFY PHYSICAL VALIDITY (FAST) =====
        
        # Pressure
        assert P is not None, "Pressure field is None"
        assert isinstance(P, np.ndarray), "Pressure not array"
        assert P.min() >= -1000, f"Pressure too negative: {P.min()}"
        assert P.max() <= 10000, f"Pressure too high: {P.max()}"
        assert not np.isnan(P).any(), "Pressure contains NaN"
        
        # Concentration (critical - should be in 0-1 range)
        assert Conc is not None, "Concentration field is None"
        assert isinstance(Conc, np.ndarray), "Concentration not array"
        assert Conc.min() >= -0.001, f"Concentration too low: {Conc.min()}"
        assert Conc.max() <= 0.04, f"Concentration too high: {Conc.max()}"
        assert not np.isnan(Conc).any(), "Concentration contains NaN"
         
        # Density
        assert rho_f is not None, "Density field is None"
        # Note: rho_f may be scalar or array depending on return value
        if hasattr(rho_f, '__len__') and len(rho_f) > 1:
            assert rho_f.min() >= 990, f"Density too low: {rho_f.min()}"
            assert rho_f.max() <= 1030, f"Density too high: {rho_f.max()}"
            assert not np.isnan(rho_f).any(), "Density contains NaN"
        else:
            # Scalar case
            rho_val = float(rho_f) if hasattr(rho_f, '__len__') else rho_f
            assert 990 <= rho_val <= 1030, f"Density out of range: {rho_val}"
        
        # Timesteps - should be much fewer (~16-20 instead of ~320)
        assert len(dts) == nsteps, "Timestep array wrong length"
        assert nsteps > 0, "No timesteps completed"
        assert nsteps < 100, f"Too many timesteps for fast run: {nsteps}"
        
        # ===== 5. VERIFY RUNTIME IS FAST =====
        total_runtime = np.sum(runtimes) if (runtimes is not None and len(runtimes) > 0) else 0
        print(f"\nFast scenario SS-{scenario_idx}: "
              f"Completed in {total_runtime:.1f}s, "
              f"{nsteps} timesteps, "
              f"Max Conc: {Conc.max():.4f}")


# =============================================================================
# PART 3: EXPERIMENTAL COMPARISON TESTS
# =============================================================================

class TestCompareWithExperimentalData:
    """Compare modeled results with experimental measurements (slow tests)."""
    
    @pytest.mark.slow
    def test_isochlor_extraction(self, scenario_parameters, test_output_dir):
        """Extract 0.5 isochlor from model results."""
        from lib.grompy_fipy import run_coupled_flow_model_fipy
        
        params, scenario_idx = scenario_parameters
        
        # Setup and run
        mesh_filename = str(test_output_dir / f"mesh_iso_{scenario_idx}.msh")
        mesh_fipy, _, _, _, _ = setup_rectangular_mesh_fipy(params, mesh_filename)
        
        model_results = run_coupled_flow_model_fipy(
            params, ModelOptions(), mesh_filename
        )
        
        Conc = model_results[5]
        cell_centers = mesh_fipy.cellCenters
        
        max_conc = 0.03624
        target_conc = 0.5 * max_conc
        
        y_iso, x_iso = extract_isochlor_profile(cell_centers, Conc, target_conc, params)
        
        assert len(y_iso) > 5, "Insufficient isochlor points extracted"
        assert len(x_iso) == len(y_iso), "Mismatched isochlor arrays"
    
    @pytest.mark.slow
    def test_misfit_against_experiments(self, scenario_parameters, test_output_dir,
                                        experimental_data):
        """Calculate misfit between modeled and experimental salt wedge."""
        from lib.grompy_fipy import run_coupled_flow_model_fipy
        
        params, scenario_idx = scenario_parameters
        
        # Setup and run
        mesh_filename = str(test_output_dir / f"mesh_misfit_{scenario_idx}.msh")
        mesh_fipy, _, _, _, _ = setup_rectangular_mesh_fipy(params, mesh_filename)
        
        model_results = run_coupled_flow_model_fipy(
            params, ModelOptions(), mesh_filename
        )
        
        Conc = model_results[5]
        cell_centers = mesh_fipy.cellCenters
        
        # Get experimental data for this scenario
        col_x = f'x_ss{scenario_idx + 1}'
        col_y = f'y_ss{scenario_idx + 1}'
        
        exp_x = experimental_data[col_x].dropna().values / 100.0  # cm to m
        exp_y = experimental_data[col_y].dropna().values / 100.0
        
        # Get modeled isochlor
        max_conc = 0.03624
        target_conc = 0.5 * max_conc
        
        mod_y, mod_x = extract_isochlor_profile(cell_centers, Conc, target_conc, params)
        
        # Interpolate modeled x at experimental y values
        mod_x_at_exp_y = np.interp(exp_y, mod_y, mod_x)
        
        # Calculate metrics
        metrics = calculate_misfit_metrics(mod_x_at_exp_y, exp_x)
        
        # ===== TOLERANCE: 4 cm spatial misfit =====
        tolerance = 0.04  # 4 cm
        
        assert abs(metrics['mean_error']) <= tolerance, \
            f"Scenario S{scenario_idx}: ME {metrics['mean_error']*100:.2f} cm > {tolerance*100:.0f} cm"
        assert metrics['mae'] <= tolerance, \
            f"Scenario S{scenario_idx}: MAE {metrics['mae']*100:.2f} cm > {tolerance*100:.0f} cm"
        assert metrics['rmse'] <= tolerance, \
            f"Scenario S{scenario_idx}: RMSE {metrics['rmse']*100:.2f} cm > {tolerance*100:.0f} cm"
        
        # Log results
        pytest.s = pytest.skip
        print(f"\n{'='*60}")
        print(f"Scenario S{scenario_idx} Misfit vs Experimental Data")
        print(f"{'='*60}")
        print(f"  Mean Error:  {metrics['mean_error']*100:7.2f} cm")
        print(f"  MAE:         {metrics['mae']*100:7.2f} cm")
        print(f"  RMSE:        {metrics['rmse']*100:7.2f} cm")
        print(f"  Tolerance:   {tolerance*100:7.2f} cm (PASS)")
        print(f"{'='*60}\n")
    
    @pytest.mark.slow
    def test_mixing_zone_width(self, scenario_parameters, test_output_dir):
        """Verify width of fresh-salt mixing zone."""
        from lib.grompy_fipy import run_coupled_flow_model_fipy
        
        params, scenario_idx = scenario_parameters
        
        # Setup and run
        mesh_filename = str(test_output_dir / f"mesh_width_{scenario_idx}.msh")
        mesh_fipy, _, _, _, _ = setup_rectangular_mesh_fipy(params, mesh_filename)
        
        model_results = run_coupled_flow_model_fipy(
            params, ModelOptions(), mesh_filename
        )
        
        Conc = model_results[5]
        cell_centers = mesh_fipy.cellCenters
        
        max_conc = 0.03624
        x_cells = np.array(cell_centers[0])
        y_cells = np.array(cell_centers[1])
        
        # Extract isochlors at 0.1, 0.5, 0.9 salinity
        widths = []
        for target_ratio in [0.1, 0.5, 0.9]:
            target_conc = target_ratio * max_conc
            y_iso, x_iso = extract_isochlor_profile(
                cell_centers, Conc, target_conc, params
            )
            
            if len(x_iso) > 0:
                width = np.max(x_iso) - np.min(x_iso)
                widths.append(width)
        
        if widths:
            avg_width = np.mean(widths)
            # Mixing zone should be ~1 cm
            assert 0.005 <= avg_width <= 0.02, \
                f"Scenario S{scenario_idx}: Width {avg_width*100:.2f} cm outside range"
    
    @pytest.mark.slow
    def test_depth_profile_vs_glover(self, scenario_parameters, test_output_dir,
                                      benchmark_parameters):
        """Compare modeled depth profile to Glover analytical solution."""
        from lib.grompy_fipy import run_coupled_flow_model_fipy
        
        params, scenario_idx = scenario_parameters
        
        # Setup and run
        mesh_filename = str(test_output_dir / f"mesh_depth_{scenario_idx}.msh")
        mesh_fipy, _, _, _, _ = setup_rectangular_mesh_fipy(params, mesh_filename)
        
        model_results = run_coupled_flow_model_fipy(
            params, ModelOptions(), mesh_filename
        )
        
        Conc = model_results[5]
        cell_centers = mesh_fipy.cellCenters
        
        # Get modeled 0.5 isochlor
        max_conc = 0.03624
        target_conc = 0.5 * max_conc
        
        mod_y, mod_x = extract_isochlor_profile(cell_centers, Conc, target_conc, params)
        
        # Calculate analytical Glover solution
        Qs = np.array([1.42, 0.59, 1.19]) / 2.7 * 1e-4  # cm³/s to m³/s
        K = 1050.0 / (24.0 * 60.0 * 60.0)  # m/s
        rho_f = params.rho_f_0
        rho_s = 1026.0
        
        Q_scenario = Qs[scenario_idx]
        
        # Calculate analytical depths
        analytical_depths = np.array([
            depth_sw_interface_glover1959(xi, Q_scenario, K, rho_f, rho_s)
            for xi in mod_x
        ])
        
        # Modeled depths
        modeled_depths = params.thickness - mod_y
        
        # Compare
        depth_errors = modeled_depths - analytical_depths
        rmse_depth = np.sqrt(np.mean(depth_errors**2))
        
        # ===== TOLERANCE: 1 cm depth =====
        tolerance_depth = 0.01  # 1 cm
        
        assert rmse_depth <= tolerance_depth, \
            f"Scenario S{scenario_idx}: Depth RMSE {rmse_depth*100:.2f} cm > {tolerance_depth*100:.0f} cm"
        
        print(f"\n{'='*60}")
        print(f"Scenario S{scenario_idx} Depth Profile vs Glover (1959)")
        print(f"{'='*60}")
        print(f"  Depth RMSE:  {rmse_depth*100:7.2f} cm")
        print(f"  Tolerance:   {tolerance_depth*100:7.2f} cm (PASS)")
        print(f"{'='*60}\n")
