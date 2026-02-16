# GroMPy-couple FiPy Backend - Comprehensive Code Analysis

**Date:** February 16, 2025  
**Analysis Focus:** Boundary Condition Enforcement and Flow-Transport Coupling  
**Status:** Critical Issues Identified

---

## Executive Summary

The FiPy backend implementation has **four major structural problems** causing zero concentration and constant pressure in simulation results:

### Root Causes (in priority order):

1. **NO PRESSURE ITERATION IN TRANSIENT LOOP** ⭐⭐⭐⭐⭐
   - Pressure is solved ONCE at t=0 (if `initial_steady_state_run=True`)
   - Never re-solved during entire time-stepping loop (lines 1018-1114)
   - Velocity field remains stale, preventing proper salt wedge development
   - Location: Line 1025-1026 (just comments, no actual solve call)

2. **MISSING INFLOW-ONLY CONCENTRATION BC LOGIC** ⭐⭐⭐⭐
   - Enforces concentration BCs at ALL boundary cells regardless of flow direction
   - Should only apply at inflow boundaries (where flux < 0)
   - At outflow, computed concentration should be free to vary
   - Location: Lines 813-827 (penalty method applied to all spec_conc_mask)
   - Compare to: gwflow_lib.py lines 635-704 (escript implementation with inflow logic)

3. **DECOUPLED DENSITY-PRESSURE FEEDBACK** ⭐⭐⭐⭐
   - Density is updated from concentration (line 1021-1023) ✓
   - But pressure equation NEVER sees this new density ✗
   - Without buoyancy forcing from density gradients, salt wedge can't intrude
   - Location: Line 674 (missing gravity source term)

4. **MISSING COUPLED ITERATION LOOP** ⭐⭐⭐
   - Escript version iterates pressure→transport→pressure until convergence
   - FiPy version does one-shot solve each timestep
   - Coupling iterations needed for accuracy and salt wedge formation
   - Location: Lines 1017-1114 (needs restructuring)

---

## Part 1: Concentration Boundary Condition Enforcement

### File Location: `/Users/elco/published_code/GroMPy-couple/lib/grompy_fipy.py`

### 1.1 Problem: Overly-Inclusive BC Mask (Line 498)

**Current Code:**
```python
# Line 498
spec_conc_mask = specified_concentration > -9999  # All non-zero entries
```

**Why It's Wrong:**
- Creates mask of ALL cells where `specified_concentration` was set to ANY value
- Initialized as all zeros (line 445): `specified_concentration = np.zeros(n_cells)`
- Then filled in loop (lines 461-495) for cells in boundary regions
- Result: `spec_conc_mask` = cells in the x,y regions of boundary conditions
- This includes many interior cells if the boundary region is large

**Expected Behavior (from escript):**
- Only cells on actual mesh boundary faces should be in the mask
- Interior cells should never be Dirichlet constrained

**Impact:**
Overly-inclusive mask combined with strong penalty method (next section) means concentration gets forced to BC values over large interior regions.

### 1.2 CRITICAL: Missing Inflow-Only Logic (Lines 813-827)

**Current Implementation (BROKEN):**
```python
# Lines 813-827
if np.any(spec_conc_mask):
    large_value = 1e10
    constraint_coeff = CellVariable(mesh=fipy_mesh, value=0.0)
    constraint_coeff.setValue(-large_value, where=spec_conc_mask)
    
    constraint_source = CellVariable(mesh=fipy_mesh, value=0.0)
    constraint_source.setValue(-large_value * specified_concentration, where=spec_conc_mask)
    
    convection_term = get_convection_term(velocity, scheme=convection_scheme)
    eq = (TransientTerm(coeff=porosity) 
          + convection_term
          == diffusion_term
          + ImplicitSourceTerm(coeff=constraint_coeff)
          - constraint_source)
```

**What's Missing:**

The escript implementation (gwflow_lib.py, lines 635-704) includes:

```python
if concentration_bnd_inflow_only is True and iteration == 0:
    # Project velocity to nodes to check flow direction
    proj.setValue(D=es.kronecker(mesh), Y=q)
    nodal_q = proj.getSolution()
    
    # Only apply BC where flow is ENTERING domain
    if concentration_bnd_inflow_direction == 'up':
        # whereNegative: q_y < 0 = upward flow into domain
        inflow_bnd = (es.whereNegative(nodal_q_norm[1]) * 
                      specified_concentration_bnd)
    elif concentration_bnd_inflow_direction == 'down':
        inflow_bnd = (es.wherePositive(nodal_q_norm[1]) *
                      specified_concentration_bnd)
    # ... etc for left/right ...
    
    # For pure outflow boundaries, use minimum x location as fallback
    if es.sup(inflow_bnd) > 0:
        active_specified_concentration_bnd = inflow_bnd
    else:
        # Use leftmost node to maintain reference concentration
        min_x = es.inf(specified_concentration_bnd * X[0])
        active_specified_concentration_bnd = (specified_concentration_bnd *
                                              es.whereZero(X[0] - min_x))
```

**Physical Meaning:**

For salt intrusion problems:
- **Inflow boundary (seawater side, upward flow):** Apply C = C_seawater BC
- **Outflow boundary (freshwater side, downward flow):** Let computed C vary freely
  - This prevents artificial seawater values at outflow
  - Allows transport equation to compute outflow concentrations

**Current FiPy Problem:**
- Applies concentration BC at ALL boundary cells uniformly
- Example: If sea surface is at y=0, and specified_concentration_surface=True:
  - Surface cells at x < 0 (seawater) get C_seawater (correct)
  - Surface cells at x > 0 (land) get C_seawater (WRONG - should be freshwater)
  - Interior cells near these boundaries also get over-constrained

**Result:**
Salt concentration becomes dominated by BC values, transport solution is over-written.

### 1.3 Equation Imbalance: Penalty Method Sign (Lines 814-827)

**Issue with Penalty Approach:**

The penalty method adds:
```
+ ImplicitSourceTerm(coeff=constraint_coeff)  # coeff = -1e10
- constraint_source                          # source = -1e10 * C_bc
```

This creates the equation:
```
phi*dC/dt + conv(C) = diff(C) - 1e10*C + 1e10*C_bc
```

Rearranged:
```
phi*dC/dt + conv(C) + 1e10*C = diff(C) + 1e10*C_bc
```

**Problems:**
1. The 1e10 coefficient is MUCH larger than diffusion/convection coefficients
   - Makes BC term dominate, over-constraining the solution
   - Numerical oscillations possible if system is ill-conditioned

2. Both penalty coefficient and source term are negative (double negative)
   - `constraint_coeff.setValue(-large_value, ...)` 
   - `constraint_source.setValue(-large_value * specified_concentration, ...)`
   - In implicit source term: `-lambda * C = -lambda * C_bc`
   - Solving: `C = C_bc` ✓ (correct, but by accident with double negative)

**Better Approach:**
Use FiPy's built-in `constrain()` method:
```python
concentration.constrain(specified_concentration, where=spec_conc_mask)
```

---

## Part 2: Pressure & Velocity Boundary Condition Enforcement

### File Location: `/Users/elco/published_code/GroMPy-couple/lib/grompy_fipy.py`

### 2.1 CRITICAL: No Pressure Iteration in Time Loop (Lines 1025-1026)

**Current Code:**
```python
# Lines 1017-1044 (time loop)
while runtime < total_time and runtime < max_runtime and timestep < max_timesteps:
    
    # Update fluid density from concentration
    density = calculate_fluid_density(
        concentration, Parameters.gamma, Parameters.rho_f_0
    )
    
    # Solve pressure equation
    # (simplified - in full implementation would iterate between pressure and transport)
    # ^ THIS IS THE PROBLEM! ^
    
    # Solve solute transport (uses STALE pressure from t=0)
    if ModelOptions.solute_transport:
        concentration_new = solve_solute_transport_fipy(...)
```

**What Should Happen:**

From escript (gwflow_lib.py, lines 600+):
```python
# For each timestep:
for iteration in range(number_of_iterations):
    # Step 1: Update pressure with current density
    pressure = solve_pressure_pde(...)
    
    # Step 2: Calculate velocity from new pressure
    v = q / phi
    
    # Step 3: Solve transport with new velocity
    concentration = solve_transport_pde(...)
    
    # Step 4: Update density from new concentration
    density = calculate_fluid_density(concentration, ...)
    
    # Step 5: Check convergence, break if converged
```

**Current FiPy Behavior:**
1. t=0: Solve pressure ONCE (if `initial_steady_state_run=True`)
   - Velocity = -(k/mu) * grad(P) at t=0
2. t=0 to t_end: Never solve pressure again
   - Velocity remains constant at t=0 value
   - Density changes but doesn't affect flow
   - Transport advects with stale velocity

**Why This Matters for Salt Intrusion:**

As saltwater (denser) enters and mixes:
1. Density increases at inflow regions
2. Density gradient should create pressure gradients (gravity-driven flow)
3. These pressure gradients drive more seawater intrusion
4. Salt wedge propagates inland

With constant pressure:
- No gravity-driven flow develops
- Seawater doesn't intrude unless initially prescribed
- No dynamic coupling between density and flow

### 2.2 Pressure Equation Missing Gravity Source Term (Line 674)

**Current Code:**
```python
# Line 674
# Build equation: div(D * grad(P)) = 0
eq = DiffusionTerm(coeff=diffusion_coeff) == 0
```

**What Should Be:**

The fundamental equation is:
```
-div(k/mu * (grad(P) - rho*g)) = Q
```

Expanding:
```
-div(k/mu * grad(P)) + div(k/mu * rho*g) = Q
```

Rearranging to FiPy form:
```
div(k/mu * grad(P)) = div(k/mu * rho*g) - Q
```

Or with P on LHS, gravity on RHS:
```
div(k/mu * grad(P)) = -k/mu * rho * g_y / (vertical unit)
```

**Why Missing Gravity is Critical:**

The gravity term creates the "buoyancy driving force":
- Without it: `q = -k/mu * grad(P)` only (just pressure gradient)
- With it: `q = -k/mu * (grad(P) - rho*g)` (pressure + gravity effects)

In a hydrostatic context:
- Denser fluid (saltwater) creates higher pressure at same depth
- This pressure difference drives seawater intrusion
- Gravity term accounts for this density effect on pressure

**Fix Required:**
```python
# Line 674, replace:
eq = DiffusionTerm(coeff=diffusion_coeff) == 0

# With:
# Gravity contribution: -k/mu * rho * g (where g points downward in -y)
gravity_coeff = -k_eff / viscosity * density * g  # cell-centered
eq = DiffusionTerm(coeff=diffusion_coeff) == CellVariable(
    mesh=fipy_mesh, value=gravity_coeff)
```

### 2.3 Pressure Penalty Method (Lines 659-677)

Similar issues to concentration penalty method:
```python
# Lines 659-677
if np.any(spec_pressure_mask):
    large_value = 1e10
    constraint_coeff = CellVariable(mesh=fipy_mesh, value=0.0)
    constraint_coeff.setValue(-large_value, where=spec_pressure_mask)
    
    constraint_source = CellVariable(mesh=fipy_mesh, value=0.0)
    constraint_source.setValue(-large_value * specified_pressure, where=spec_pressure_mask)
    
    eq = (DiffusionTerm(coeff=diffusion_coeff) 
          + ImplicitSourceTerm(coeff=constraint_coeff) 
          == constraint_source)
```

**Recommendations:**
1. Use `constrain()` instead of penalty method (cleaner)
2. Or reduce large_value to something more reasonable (1e6 or 1e8)
3. Or use FiPy's proper Dirichlet BC interface (different approach)

---

## Part 3: Flow-Transport Coupling

### File Location: `/Users/elco/published_code/GroMPy-couple/lib/grompy_fipy.py`

### 3.1 Missing Coupled Iteration Loop (Lines 1017-1044)

**Current Structure:**
```python
# Simple time-stepping without coupling iterations
while runtime < total_time:
    density = calculate_fluid_density(concentration, ...)  # Update
    # NO pressure solve!
    concentration = solve_transport(pressure_stale, ...)   # Single solve
    runtime += dt
```

**Required Structure:**
```python
# Coupled iteration approach (from escript)
while runtime < total_time:
    # Coupling iteration loop
    for iteration in range(ModelOptions.coupled_iterations):
        # Solve pressure with current density
        if iteration == 0 or coupled_iterations == True:
            pressure = solve_pressure(density, ...)
            v = q / porosity
        
        # Solve transport with current velocity
        concentration = solve_transport(v, ...)
        
        # Check convergence
        press_conv = max(abs(pressure - pressure_old))
        conc_conv = max(abs(concentration - concentration_old))
        if press_conv < tol and conc_conv < tol:
            break
    
    runtime += dt
```

**Why Coupling Iterations Matter:**

For salt intrusion:
1. Initial pressure drives flow
2. Flow carries salt, increasing density
3. Higher density creates new pressure gradients
4. New pressure drives more flow
5. Iterations ensure density changes are reflected in pressure solution

Without iterations:
- Density-driven circulation incomplete
- Salt wedge penetration underestimated
- Solution less accurate

### 3.2 Velocity Field Updates (Lines 787-793)

**Current Code:**
```python
# Lines 787-793
q = calculate_darcy_velocity_fipy(
    fipy_mesh, pressure, density, k_tensor, viscosity, g
)

velocity = q / porosity
```

**Status:** ✓ Formula is correct

**Problem:** This calculates from STALE pressure (not updated in loop)

**Fix:** Calculate after pressure solve in coupling loop

### 3.3 Darcy Velocity Calculation (Lines 96-147)

**Implementation:**
```python
# Line 144-145
# Darcy velocity: q = -k/mu * (grad(P) - rho*g)
q = -(k_eff / viscosity) * (grad_P - rho_face * g_vector)
```

**Status:** ✓ Formula is correct

**Issue:** `k_eff` uses geometric mean of tensor components
- Line 142: `k_eff = np.sqrt(kxx * kyy)`
- This approximates anisotropic permeability as isotropic
- Acceptable for groundwater applications but not exact

---

## Part 4: Initial Conditions

### File Location: `/Users/elco/published_code/GroMPy-couple/lib/grompy_fipy.py`

### 4.1 Initial Pressure (Lines 581-584)

**Implementation:**
```python
depth = z_surface - y
pressure = density * Parameters.g * np.maximum(0, depth)
```

**Status:** ✓ Correct hydrostatic initialization

### 4.2 Initial Concentration (Lines 538-576)

**Implementation:**
```python
if Parameters.ghyben_herzberg:
    # Ghyben-Herzberg approximation for salt wedge
    # ...
else:
    # Uniform freshwater
    concentration = np.full(n_cells, freshwater_conc)
```

**Status:** ✓ Correct initial conditions

### 4.3 Initial Steady-State Run (Lines 940-1010)

**Issue:**
- Pressure solve happens once at t=0
- Concentration BC enforcement happens at t=0
- But then nothing in transient loop

---

## Root Cause Analysis

### Why Zero Concentration?

**Mechanism:**
1. `specified_concentration` initialized to zeros (line 445)
2. Only boundary region cells get non-zero values (lines 461-495)
3. `spec_conc_mask` includes many cells due to region specification
4. Transport equation with large penalty term (1e10 * C = 1e10 * C_bc)
5. When C_bc = 0, penalty forces C = 0 everywhere in mask
6. Over-constrained equation: C → 0 (penalty dominates diffusion/convection)

**Why Not Advected:**
- Velocity field is stale from t=0
- Transport: `phi*dC/dt + div(q*C) = div(D*grad(C))`
- With small stale q, mostly diffusion acts
- Molecular diffusivity D_m << dispersivity-driven diffusion
- Tiny diffusive flux insufficient to overcome zero BC penalty

### Why Constant Pressure?

**Mechanism:**
1. Pressure solved only at t=0 (if initial_steady_state_run=True)
2. Time loop (1018+) never re-solves pressure
3. Density updates don't feed back to pressure equation
4. Pressure stays at initial hydrostatic values

**Why No Buoyancy:**
- Gravity term missing from pressure equation (line 674)
- Even if re-solved, without `- k/mu * rho * g` term:
  - Equation can't represent density-driven flow
  - Pressure gradients purely from boundary conditions
  - No internal circulation develops

---

## Summary of Issues

| Issue | Severity | Location | Root Cause |
|-------|----------|----------|-----------|
| No pressure iteration | ⭐⭐⭐⭐⭐ | Lines 1025-1026 | Comment instead of code |
| Missing inflow-only BC | ⭐⭐⭐⭐ | Lines 813-827 | Feature not implemented |
| Penalty method too strong | ⭐⭐⭐ | Lines 814-827 | 1e10 coefficient |
| Missing gravity term | ⭐⭐⭐ | Line 674 | Incorrect equation |
| No coupled iteration | ⭐⭐⭐ | Lines 1017-1044 | Single-pass design |
| Over-inclusive BC mask | ⭐⭐ | Line 498 | Region specification |
| Double-negative in penalty | ⭐⭐ | Lines 817-820 | Sign handling |

---

## What's Working Correctly

1. ✓ Mesh loading and parsing
2. ✓ Initial condition setup (hydrostatic, Ghyben-Herzberg)
3. ✓ Density calculation from concentration
4. ✓ Darcy velocity formula (when given updated pressure)
5. ✓ Dispersion tensor (anisotropic with cross-terms)
6. ✓ Convection term selection (exponential scheme available)
7. ✓ VTK output
8. ✓ Time stepping framework

---

## Recommended Fix Priority

### 1. Add Pressure Solve to Time Loop ⭐⭐⭐⭐⭐ (5 min)
```python
# Line 1025, replace comment with:
pressure = solve_steady_state_pressure_fipy(
    fipy_mesh, backend, k_tensor, Parameters.viscosity,
    density, Parameters.g,
    Parameters.recharge_flux, Parameters.recharge_density,
    bc['recharge_mask'],
    bc['spec_pressure_mask'], bc['specified_pressure'],
    Parameters
)
```

### 2. Add Coupled Iteration Loop ⭐⭐⭐⭐⭐ (30 min)
Restructure lines 1018-1045 to iterate pressure-concentration solving

### 3. Implement Inflow-Only BC Logic ⭐⭐⭐⭐ (45 min)
Add flow direction checking before applying concentration BC

### 4. Add Gravity Source Term ⭐⭐⭐ (10 min)
Modify line 674 to include `-k/mu * rho * g` source

### 5. Switch to Constrain Method ⭐⭐ (20 min)
Replace penalty method with `variable.constrain()`

---

## Testing Checklist

After fixes:
- [ ] Pressure changes during simulation (not constant at t > 0)
- [ ] Concentration changes during simulation (not constant)
- [ ] Salt concentration positive at inflow, zero at outflow
- [ ] Velocity field updates with time
- [ ] Density gradient creates pressure gradients (buoyancy works)
- [ ] Salt wedge penetrates inland (~40% domain for standard parameters)
- [ ] No artificial seawater at freshwater boundaries
- [ ] Coupled iterations improve solution in first few timesteps
- [ ] Total mass of salt conserved

---

## References to Original Implementations

**Escript Version (gwflow_lib.py):**
- Lines 635-704: Inflow-only BC logic
- Lines 600+: Coupled iteration structure  
- Lines 755-790: Density-driven coupling

**FiPy Backend Current:**
- Lines 96-147: Darcy velocity calculation (correct formula)
- Lines 150-226: Dispersion tensor (good implementation)
- Lines 589-679: Pressure solver (needs fixes)
- Lines 722-837: Transport solver (needs fixes)

