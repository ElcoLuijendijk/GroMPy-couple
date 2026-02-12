# GroMPy-couple: FiPy Advection Implementation Plan

**Date:** February 4, 2026  
**Status:** Ready for implementation

---

## Background

The FiPy backend for GroMPy-couple is functional but **missing the advection term** in the solute transport solver. This means solute is not transported with groundwater flow - only diffusion occurs.

### Current Problem

The FiPy transport equation (in `lib/grompy_fipy.py:439-445`) only has diffusion:

```python
eq = (TransientTerm(coeff=porosity) 
      == DiffusionTerm(coeff=D_eff)
      + ImplicitSourceTerm(coeff=constraint_coeff)
      - constraint_source)
```

The escript version (in `lib/gwflow_lib.py:329-331`) includes advection:

```python
a_coeff = dt * dispersion_tensor  # Diffusion
c_coeff = dt * v                   # ADVECTION (velocity) - MISSING IN FIPY
d_coeff = 1                        # Mass matrix
```

### Evidence of the Problem

- Concentration change rate is `0.00e+00 kg/kg/yr` even after 106 timesteps
- Initial Ghyben-Herzberg concentration distribution doesn't evolve
- Salt wedge remains static without advective transport

---

## Implementation Plan

### 1. Add FaceVariable to FiPy Imports

**File:** `lib/grompy_fipy.py`  
**Line:** 36

**Change from:**
```python
from fipy import CellVariable, DiffusionTerm, ConvectionTerm, TransientTerm
```

**Change to:**
```python
from fipy import CellVariable, FaceVariable, DiffusionTerm, ConvectionTerm, TransientTerm
```

---

### 2. Create Darcy Velocity Helper Function

**File:** `lib/grompy_fipy.py`  
**Location:** After line ~86 (after `calculate_fluid_density()` function)

**Add new function:**

```python
def calculate_darcy_velocity_fipy(fipy_mesh, pressure, density, k_tensor, viscosity, g):
    """
    Calculate Darcy velocity from pressure field.
    
    Darcy's law with gravity:
    q = -k/mu * (grad(P) - rho*g)
    
    where g is the gravity vector pointing in -y direction.
    
    Parameters
    ----------
    fipy_mesh : fipy.Mesh
        The FiPy mesh object
    pressure : np.ndarray
        Pressure field (cell-centered)
    density : np.ndarray
        Fluid density field (cell-centered)
    k_tensor : tuple
        Permeability tensor ((kxx, kxy), (kyx, kyy))
    viscosity : float
        Fluid viscosity (Pa.s)
    g : float
        Gravitational acceleration (m/s²)
    
    Returns
    -------
    FaceVariable
        Darcy velocity vector at cell faces (rank-1)
    """
    # Create cell variables for pressure and density
    P = CellVariable(mesh=fipy_mesh, value=pressure)
    rho = CellVariable(mesh=fipy_mesh, value=density)
    
    # Get density at faces (arithmetic average of neighboring cells)
    rho_face = rho.arithmeticFaceValue
    
    # Pressure gradient at faces
    grad_P = P.faceGrad  # Shape: (2, nFaces)
    
    # Gravity vector (pointing in -y direction, i.e., downward)
    g_vector = FaceVariable(mesh=fipy_mesh, rank=1, value=[[0.], [-g]])
    
    # Effective isotropic permeability (geometric mean)
    # Note: FiPy has issues with anisotropic tensor + ImplicitSourceTerm
    kxx = k_tensor[0][0]
    kyy = k_tensor[1][1]
    k_eff = np.sqrt(kxx * kyy)
    
    # Darcy velocity: q = -k/mu * (grad(P) - rho*g)
    q = -(k_eff / viscosity) * (grad_P - rho_face * g_vector)
    
    return q
```

---

### 3. Modify Transport Solver to Include Advection

**File:** `lib/grompy_fipy.py`  
**Function:** `solve_solute_transport_fipy()` (lines 354-450)

#### 3a. Add velocity calculation after line 417

After creating the concentration variable, add:

```python
# Calculate Darcy velocity from pressure gradient
q = calculate_darcy_velocity_fipy(
    fipy_mesh, pressure, density, k_tensor, viscosity, g
)

# Pore velocity = Darcy velocity / porosity
velocity = q / porosity
```

#### 3b. Modify equation to include ConvectionTerm

**Change from (lines 439-445):**
```python
eq = (TransientTerm(coeff=porosity) 
      == DiffusionTerm(coeff=D_eff)
      + ImplicitSourceTerm(coeff=constraint_coeff)
      - constraint_source)
```

**Change to:**
```python
eq = (TransientTerm(coeff=porosity) 
      + ConvectionTerm(coeff=velocity)
      == DiffusionTerm(coeff=D_eff)
      + ImplicitSourceTerm(coeff=constraint_coeff)
      - constraint_source)
```

Also update the `else` branch (lines 443-445) similarly:

**Change from:**
```python
eq = (TransientTerm(coeff=porosity) == 
      DiffusionTerm(coeff=D_eff))
```

**Change to:**
```python
eq = (TransientTerm(coeff=porosity) 
      + ConvectionTerm(coeff=velocity)
      == DiffusionTerm(coeff=D_eff))
```

---

## Technical Notes

### FiPy Sign Convention
- `ConvectionTerm` goes on the **left-hand side** (same side as `TransientTerm`)
- This represents: `∂C/∂t + ∇·(vC) = ∇·(D∇C)`
- FiPy's PowerLawConvectionTerm (default) is stable and handles upwinding automatically

### Velocity Calculation
- `P.faceGrad` gives face-centered gradient directly from cell-centered pressure
- `rho.arithmeticFaceValue` interpolates cell density to faces
- Gravity vector `[0, -g]` points downward (y-axis points up)

### Isotropic Permeability
- Using geometric mean `k_eff = sqrt(kxx * kyy)` as approximation
- Full anisotropic tensor causes issues with FiPy's `ImplicitSourceTerm`

### No CFL Check Needed
- FiPy's implicit ConvectionTerm is unconditionally stable
- No timestep restriction from advection

---

## Testing

### Unit Tests
```bash
cd /Users/elco/published_code/GroMPy-couple
python -m pytest tests/ -v
```

Expected: 54 tests pass, 31 skipped (escript-specific)

### Full Model Run
```bash
python grompy.py model_input/model_parameters.py
```

**Success criteria:**
- Concentration change rate is **non-zero** (e.g., `-1.5e-06 kg/kg/yr`)
- Salt wedge interface evolves from initial Ghyben-Herzberg position
- Model runs through full simulation time without numerical instabilities

---

## Files to Modify

| File | Changes |
|------|---------|
| `lib/grompy_fipy.py` | Add FaceVariable import, add velocity helper function, add ConvectionTerm |

---

## Future Enhancements (Not in This Plan)

1. **Velocity-dependent dispersion tensor** - Currently using simplified isotropic diffusion. Could implement full tensor like escript version.
2. **Seepage face boundary conditions** - Currently not implemented in FiPy backend.
3. **Anisotropic permeability** - Need to find workaround for FiPy's tensor handling.

---

## References

- FiPy documentation: https://www.ctcms.nist.gov/fipy/
- Original escript implementation: `lib/gwflow_lib.py`
- FiPy ConvectionTerm: Uses PowerLaw scheme by default (stable upwinding)
