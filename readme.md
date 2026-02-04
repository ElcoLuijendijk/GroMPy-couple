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
