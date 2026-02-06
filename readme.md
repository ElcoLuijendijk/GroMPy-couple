# GroMPy-couple

A 2D numerical model for coupled density-driven groundwater flow and solute transport in coastal aquifers.

## Overview

GroMPy-couple solves the coupled equations for:
- Groundwater flow (pressure equation)
- Solute transport (advection-diffusion equation)
- Variable-density effects

The model supports two numerical backends:
- **esys-escript** (Finite Element Method) - original implementation
- **FiPy** (Finite Volume Method) - new alternative backend

## Installation

### Dependencies

Core dependencies:
```
numpy
pandas
```

For FiPy backend:
```
fipy
gmsh
meshio
```

For escript backend:
```
esys-escript
```

### Installing FiPy Backend Dependencies

```bash
pip install fipy gmsh meshio numpy pandas
```

Note: The FiPy backend does not require esys-escript, making it easier to install and run on most systems.

## Usage

### Basic Usage

```bash
python grompy.py model_input/model_parameters.py
```

### Selecting the Backend

The backend is selected in your model parameters file:

```python
# In model_input/model_parameters.py

class ModelOptions:
    backend = 'fipy'  # Use FiPy backend (default for new installations)
    # backend = 'escript'  # Use esys-escript backend
```

### Running with FiPy Backend

1. Ensure FiPy and dependencies are installed:
   ```bash
   pip install fipy gmsh meshio
   ```

2. Set `backend = 'fipy'` in your model parameters file

3. Run the model:
   ```bash
   python grompy.py model_input/model_parameters.py
   ```

## Model Parameters

Key parameters are defined in `model_input/model_parameters.py`:

### Physical Parameters
- `k` - Permeability (m^2)
- `porosity` - Porosity (-)
- `diffusivity` - Molecular diffusivity (m^2/s)
- `l_disp` - Longitudinal dispersivity (m)
- `anisotropy` - Anisotropy ratio (kx/ky)

### Domain Parameters
- `L` - Domain length (m)
- `D` - Domain depth (m)
- `cellsize` - Mesh element size (m)

### Boundary Conditions
- `recharge_flux` - Recharge rate (m/s)
- `seawater_concentration` - Seawater salinity (kg/kg)

## VTK Output Visualization

Both escript and FiPy backends support VTK output for visualization in ParaView and other VTK-compatible tools.

### Requirements

**For escript backend:** Automatic via `esys.weipa` module

**For FiPy backend:** Install optional dependencies (recommended):
```bash
pip install vtkmodules    # Native VTK Python bindings (recommended)
pip install meshio        # Fallback VTK writer (fallback)
```

### Output Files

- **Location:** `model_output/vtk_files/`
- **Format:** `.vtu` (XML unstructured grid with compression)
- **Precision:** float32 (50% smaller files than float64)
- **Variables:** 18 cell-centered fields including:
  - Pressure, concentration, hydraulic head
  - Velocity components (x, y) and magnitude
  - Permeability (diagonal components)
  - Boundary condition masks and nodal flux
- **Metadata:** Embedded model parameters and boundary conditions (accessible in ParaView Object Inspector)

### Opening VTK Files

**ParaView (recommended, free and open-source):**
```bash
# Install ParaView if needed
# https://www.paraview.org/download/

paraview model_output/vtk_files/*.vtu
```

**Other options:**
- VisIt (https://visit.llnl.gov/)
- Gmsh (https://gmsh.info/)
- ICEM CFD, Salome, or other VTK-compatible viewers

### Variable Descriptions

| Variable | Units | Description |
|----------|-------|-------------|
| pressure | Pa | Fluid pressure |
| concentration | ppt | Solute concentration |
| hydraulic_head | m | Hydraulic head |
| velocity_x, velocity_y | m/s | Darcy velocity components |
| flux_magnitude | m/s | Magnitude of Darcy velocity |
| permeability_xx, permeability_yy | m² | Diagonal permeability tensor components |
| nodal_flux | m³/s | Surface nodal flux (non-zero on surface) |
| surface_mask, sea_surface_mask | 0/1 | Boolean masks for surface identification |
| specified_pressure_bnd | 0/1 | Specified pressure boundary mask |
| active_seepage_bnd | 0/1 | Active seepage boundary mask |
| recharge_bnd | 0/1 | Recharge boundary mask |
| active_concentration_bnd | 0/1 | Active concentration boundary mask |

### Metadata Attributes

The following metadata is embedded in each VTK file and visible in ParaView's Object Inspector:

**Model Parameters:**
- porosity, diffusivity, permeability (isotropic estimate)
- dispersivity_longitudinal, dispersivity_transverse
- viscosity, fluid_density_initial, gravity

**Boundary Conditions:**
- left_pressure_head, right_pressure_head
- initial_salinity, seepage_active
- recharge_flux

**Flux Statistics:**
- flux_total_seepage_flux, flux_total_submarine_flux_out
- flux_total_land_flux_in, flux_total_land_flux_out
- flux_total_submarine_flux_in

**Simulation Info:**
- backend ("fipy" or "escript")
- vtk_format_version, num_timesteps

### Troubleshooting VTK Output

**Problem:** "No VTK writer available"
- **Solution:** Install pyVTK: `pip install vtkmodules`
- **Fallback:** Also install meshio: `pip install meshio`
- VTK output will be skipped if no writer is available (model continues normally)

**Problem:** VTK file is too large
- **Solution:** Model uses float32 precision and compression by default
- For 50x100 cell mesh: typically 1-2 MB compressed vs 4-5 MB uncompressed
- Consider reducing mesh resolution for large models

**Problem:** File doesn't open in ParaView
- **Solution:** Ensure file was written successfully (check console output)
- Try with different ParaView version
- Verify vtkmodules installation: `python -c "import vtkmodules.all as vtk; print('OK')"`

## Backend Comparison

| Feature | esys-escript | FiPy |
|---------|--------------|------|
| Method | Finite Element | Finite Volume |
| Installation | Complex | Simple (pip) |
| Parallel | MPI support | Limited |
| Anisotropy | Full tensor | Isotropic approximation |

### FiPy Backend Notes

The FiPy backend uses:
- Penalty method for Dirichlet boundary conditions on internal cells
- Isotropic effective permeability (geometric mean) due to FiPy limitations with tensor diffusion
- Gmsh for mesh generation with MSH format 2.2 compatibility

## Project Structure

```
GroMPy-couple/
├── grompy.py              # Main entry point
├── lib/
│   ├── backend/           # Backend abstraction layer
│   │   ├── __init__.py    # Backend factory
│   │   ├── base.py        # Abstract base classes
│   │   ├── escript_backend.py
│   │   └── fipy_backend.py
│   ├── grompy_fipy.py     # FiPy solver implementation
│   ├── grompy_lib.py      # Shared utilities
│   ├── mesh_functions.py  # escript mesh creation
│   ├── mesh_functions_fipy.py  # FiPy mesh creation
│   └── fipy_mesh_io.py    # Mesh I/O for FiPy
├── model_input/
│   └── model_parameters.py
├── model_output/          # Output directory
└── tests/                 # Test suite
```

## Testing

Run the test suite:
```bash
python -m pytest tests/ -v
```

FiPy backend tests:
```bash
python -m pytest tests/test_backend.py tests/test_mesh_functions_fipy.py -v
```

## License

[Add license information here]

## Citation

[Add citation information here]

## Contributing

[Add contributing guidelines here]
