"""
FiPy backend implementation for GroMPy-couple.

This module provides the FiPy (Finite Volume Method) implementation of the
backend interface. It enables solving the coupled groundwater flow and
solute transport equations using FiPy's finite volume discretization.

Notes
-----
- FiPy uses a different discretization (FVM) than escript (FEM), which may
  lead to slightly different numerical results.
- Anisotropic diffusion tensors are supported via FiPy's DiffusionTerm.
- The sequential solving approach (pressure then concentration) avoids
  known FiPy bug #361 with coupled/vector equations.
"""

from typing import Tuple, Optional, Union, Any
import numpy as np

try:
    import fipy
    from fipy import (
        CellVariable, FaceVariable, 
        DiffusionTerm, ConvectionTerm, TransientTerm,
        ImplicitSourceTerm,
        Grid2D, Tri2D,
    )
    from fipy.solvers import LinearGMRESSolver, LinearPCGSolver
    from fipy.boundaryConditions import FixedValue, FixedFlux
    FIPY_AVAILABLE = True
except ImportError:
    FIPY_AVAILABLE = False

# Optional PETSc solvers — available when running under mpirun with petsc4py installed.
# Falls back silently to the default SciPy-based solvers when not available.
try:
    from fipy.solvers.petsc import LinearGMRESSolver as PETScGMRESSolver
    from fipy.solvers.petsc import LinearPCGSolver as PETScPCGSolver
    PETSC_AVAILABLE = True
except ImportError:
    PETSC_AVAILABLE = False

from lib.backend.base import BackendBase, MeshBase, FieldBase, PDESolverBase


class FiPyField(FieldBase):
    """
    FiPy implementation of FieldBase.
    
    Wraps FiPy CellVariable and FaceVariable objects.
    """
    
    def __init__(self, variable: Union['CellVariable', 'FaceVariable', np.ndarray],
                 mesh: 'fipy.meshes.mesh.Mesh' = None,
                 name: str = ""):
        """
        Initialize a FiPyField.
        
        Parameters
        ----------
        variable : CellVariable, FaceVariable, or np.ndarray
            FiPy variable or numpy array
        mesh : fipy.meshes.mesh.Mesh, optional
            FiPy mesh (required if variable is np.ndarray)
        name : str
            Name of the field
        """
        super().__init__(name)
        
        if isinstance(variable, np.ndarray):
            if mesh is None:
                raise ValueError("Mesh required when creating field from array")
            # Check if this is a 2x2 tensor (for anisotropic diffusion)
            if variable.ndim == 2 and variable.shape == (2, 2):
                # Store as raw numpy array for tensor operations
                self._variable = variable
                self._is_numpy = True
            else:
                # Create CellVariable from array
                self._variable = CellVariable(mesh=mesh, value=variable, name=name)
                self._is_numpy = False
        elif isinstance(variable, (int, float)):
            if mesh is None:
                raise ValueError("Mesh required when creating field from scalar")
            self._variable = CellVariable(mesh=mesh, value=variable, name=name)
            self._is_numpy = False
        else:
            self._variable = variable
            self._is_numpy = False
        
        self._mesh = mesh if mesh is not None else getattr(variable, 'mesh', None)
    
    @property
    def variable(self) -> Union['CellVariable', 'FaceVariable']:
        """Get the underlying FiPy variable."""
        return self._variable
    
    @property
    def mesh(self) -> 'fipy.meshes.mesh.Mesh':
        """Get the associated FiPy mesh."""
        return self._mesh
    
    def get_values(self) -> np.ndarray:
        """Get field values as numpy array."""
        if hasattr(self._variable, 'value'):
            val = np.array(self._variable.value)
        else:
            val = np.array(self._variable)
        return val.flatten() if val.ndim > 1 else val
    
    def set_values(self, values: Union[np.ndarray, float]) -> None:
        """Set field values."""
        if hasattr(self._variable, 'setValue'):
            self._variable.setValue(values)
        else:
            # For raw arrays, replace the data
            self._variable[:] = values
    
    def copy(self) -> 'FiPyField':
        """Create a copy of the field."""
        if hasattr(self._variable, 'copy'):
            new_var = self._variable.copy()
        else:
            new_var = CellVariable(mesh=self._mesh, value=self.get_values().copy())
        return FiPyField(new_var, self._mesh, self.name)
    
    def get_shape(self) -> Tuple[int, ...]:
        """Get the shape of the field data."""
        if hasattr(self._variable, 'shape'):
            return self._variable.shape
        return (len(self._variable),)
    
    def get_coordinates(self) -> np.ndarray:
        """Get x, y coordinates associated with field values."""
        if self._mesh is None:
            raise ValueError("No mesh associated with field")
        
        centers = self._mesh.cellCenters
        return np.array([centers[0], centers[1]]).T
    
    # Arithmetic operations
    def _ensure_compatible(self, other):
        """Ensure other operand is compatible."""
        if isinstance(other, FiPyField):
            return other._variable
        return other
    
    def __add__(self, other: Union['FiPyField', float]) -> 'FiPyField':
        other_val = self._ensure_compatible(other)
        result = self._variable + other_val
        return FiPyField(result, self._mesh, self.name)
    
    def __radd__(self, other: Union['FiPyField', float]) -> 'FiPyField':
        return self.__add__(other)
    
    def __sub__(self, other: Union['FiPyField', float]) -> 'FiPyField':
        other_val = self._ensure_compatible(other)
        result = self._variable - other_val
        return FiPyField(result, self._mesh, self.name)
    
    def __rsub__(self, other: Union['FiPyField', float]) -> 'FiPyField':
        other_val = self._ensure_compatible(other)
        result = other_val - self._variable
        return FiPyField(result, self._mesh, self.name)
    
    def __mul__(self, other: Union['FiPyField', float]) -> 'FiPyField':
        other_val = self._ensure_compatible(other)
        result = self._variable * other_val
        return FiPyField(result, self._mesh, self.name)
    
    def __rmul__(self, other: Union['FiPyField', float]) -> 'FiPyField':
        return self.__mul__(other)
    
    def __truediv__(self, other: Union['FiPyField', float]) -> 'FiPyField':
        other_val = self._ensure_compatible(other)
        result = self._variable / other_val
        return FiPyField(result, self._mesh, self.name)
    
    def __neg__(self) -> 'FiPyField':
        result = -self._variable
        return FiPyField(result, self._mesh, self.name)
    
    def __pow__(self, power: float) -> 'FiPyField':
        result = self._variable ** power
        return FiPyField(result, self._mesh, self.name)
    
    def __getitem__(self, index: int) -> 'FiPyField':
        """Get component of vector field."""
        if hasattr(self._variable, '__getitem__'):
            return FiPyField(self._variable[index], self._mesh, self.name)
        raise IndexError("Field does not support indexing")
    
    def __setitem__(self, index: int, value: Union['FiPyField', float]) -> None:
        """Set component of vector field."""
        if isinstance(value, FiPyField):
            self._variable[index] = value._variable
        else:
            self._variable[index] = value


class FiPyMesh(MeshBase):
    """
    FiPy implementation of MeshBase.
    
    Wraps FiPy mesh objects.
    """
    
    def __init__(self, mesh: 'fipy.meshes.mesh.Mesh'):
        """
        Initialize a FiPyMesh.
        
        Parameters
        ----------
        mesh : fipy.meshes.mesh.Mesh
            FiPy mesh object
        """
        self._mesh = mesh
    
    @property
    def fipy_mesh(self) -> 'fipy.meshes.mesh.Mesh':
        """Get the underlying FiPy mesh."""
        return self._mesh
    
    @property
    def num_cells(self) -> int:
        """Get number of cells."""
        return self._mesh.numberOfCells
    
    @property
    def num_nodes(self) -> int:
        """Get number of nodes (vertices)."""
        try:
            return self._mesh.vertexCoords.shape[1]
        except AttributeError:
            # Fallback for meshes without vertex coordinates
            return self._mesh.numberOfCells
    
    def get_cell_centers(self) -> np.ndarray:
        """Get coordinates of cell centers."""
        centers = self._mesh.cellCenters
        return np.array([centers[0], centers[1]]).T
    
    def get_node_coordinates(self) -> np.ndarray:
        """
        Get coordinates of mesh nodes.
        
        Note: FiPy doesn't directly expose vertex coordinates in the same way.
        We return face centers as an approximation for boundary operations.
        """
        # Get vertex coordinates if available
        try:
            verts = self._mesh.vertexCoords
            return np.array([verts[0], verts[1]]).T
        except:
            # Fallback to cell centers
            return self.get_cell_centers()
    
    def get_boundary_mask(self, boundary_name: str) -> np.ndarray:
        """Get boolean mask for boundary cells."""
        centers = self.get_cell_centers()
        x = centers[:, 0]
        y = centers[:, 1]
        
        tol = 1e-6
        x_min, x_max = x.min(), x.max()
        y_min, y_max = y.min(), y.max()
        
        # For FiPy, we identify cells near boundaries
        # This is an approximation since FiPy is cell-centered
        dx = (x_max - x_min) / np.sqrt(self.num_cells)  # Approximate cell size
        
        if boundary_name == 'left':
            return (x - x_min) < dx
        elif boundary_name == 'right':
            return (x_max - x) < dx
        elif boundary_name == 'bottom':
            return (y - y_min) < dx
        elif boundary_name == 'top':
            return (y_max - y) < dx
        else:
            return np.zeros(len(x), dtype=bool)
    
    def get_native_mesh(self) -> Any:
        """Get the native FiPy mesh object."""
        return self._mesh


class FiPyPDESolver(PDESolverBase):
    """
    FiPy implementation of PDESolverBase.
    
    Solves PDEs using FiPy's finite volume discretization.
    The general steady-state form being solved is:
    
        -div(A * grad(u)) + div(C * u) + D*u = Y
    
    This is assembled using FiPy's DiffusionTerm, ConvectionTerm, and
    ImplicitSourceTerm.
    """
    
    def __init__(self, mesh: FiPyMesh, symmetric: bool = True):
        """
        Initialize FiPyPDESolver.
        
        Parameters
        ----------
        mesh : FiPyMesh
            Mesh to solve on
        symmetric : bool
            Whether the PDE is symmetric (affects solver choice)
        """
        self._mesh = mesh
        self._fipy_mesh = mesh.fipy_mesh
        self._symmetric = symmetric
        
        # Create the solution variable
        self._solution = CellVariable(mesh=self._fipy_mesh, value=0.0, name="solution")
        
        # Store coefficients
        self._A = 1.0  # Diffusion coefficient
        self._C = None  # Advection coefficient (vector)
        self._D = 0.0  # Reaction/source coefficient
        self._Y = 0.0  # Source term
        
        # Boundary conditions
        self._dirichlet_mask = None
        self._dirichlet_value = None
        self._neumann_flux = None
        
        # Solver options
        self._solver_method = 'GMRES'
        self._tolerance = 1e-8
        self._max_iterations = 10000
    
    def set_coefficients(self,
                        A: Optional[Union[FiPyField, float, np.ndarray]] = None,
                        C: Optional[Union[FiPyField, np.ndarray]] = None,
                        D: Optional[Union[FiPyField, float]] = None,
                        X: Optional[Union[FiPyField, np.ndarray]] = None,
                        Y: Optional[Union[FiPyField, float]] = None) -> None:
        """Set PDE coefficients."""
        if A is not None:
            if isinstance(A, FiPyField):
                # Check if underlying variable is a tensor (2x2 array)
                raw = A._variable
                if isinstance(raw, np.ndarray) and raw.ndim == 2 and raw.shape == (2, 2):
                    # Store as 2x2 tensor for anisotropic diffusion
                    self._A = [[raw[0, 0], raw[0, 1]], [raw[1, 0], raw[1, 1]]]
                else:
                    self._A = A.get_values()
            elif isinstance(A, np.ndarray):
                # Check if it's a tensor (for anisotropic diffusion)
                if A.ndim == 2 and A.shape == (2, 2):
                    # Convert 2x2 tensor to FiPy format for anisotropic diffusion
                    self._A = [[A[0, 0], A[0, 1]], [A[1, 0], A[1, 1]]]
                else:
                    self._A = A
            else:
                self._A = A
        
        if C is not None:
            if isinstance(C, FiPyField):
                self._C = C.variable
            else:
                self._C = C
        
        if D is not None:
            if isinstance(D, FiPyField):
                self._D = D.get_values()
            else:
                self._D = D
        
        if Y is not None:
            if isinstance(Y, FiPyField):
                self._Y = Y.get_values()
            else:
                self._Y = Y
    
    def set_boundary_conditions(self,
                               dirichlet_mask: Optional[FiPyField] = None,
                               dirichlet_value: Optional[Union[FiPyField, float]] = None,
                               neumann_flux: Optional[Union[FiPyField, float]] = None) -> None:
        """Set boundary conditions."""
        if dirichlet_mask is not None:
            if isinstance(dirichlet_mask, FiPyField):
                self._dirichlet_mask = dirichlet_mask.get_values()
            else:
                self._dirichlet_mask = dirichlet_mask
        
        if dirichlet_value is not None:
            if isinstance(dirichlet_value, FiPyField):
                self._dirichlet_value = dirichlet_value.get_values()
            else:
                self._dirichlet_value = dirichlet_value
        
        if neumann_flux is not None:
            if isinstance(neumann_flux, FiPyField):
                self._neumann_flux = neumann_flux.get_values()
            else:
                self._neumann_flux = neumann_flux
    
    def solve(self) -> FiPyField:
        """Solve the PDE."""
        # Create a fresh solution variable to avoid constraint accumulation
        self._solution = CellVariable(mesh=self._fipy_mesh, value=0.0, name="solution")
        
        # Apply Dirichlet boundary conditions using face-based constraints
        # FiPy requires constraining at faces, not cells
        if self._dirichlet_mask is not None and self._dirichlet_value is not None:
            mask = np.array(self._dirichlet_mask) > 0.5
            if isinstance(self._dirichlet_value, (int, float)):
                bc_values = np.full(len(mask), self._dirichlet_value)
            else:
                bc_values = np.array(self._dirichlet_value)
            
            # For cells marked as boundary, we need to identify which faces
            # are on the external boundary and apply constraints there
            # Get cell centers and identify boundary faces
            cell_centers = self._mesh.get_cell_centers()
            x, y = cell_centers[:, 0], cell_centers[:, 1]
            x_min, x_max = x.min(), x.max()
            y_min, y_max = y.min(), y.max()
            
            # Compute approximate cell size from unique coordinate spacing
            unique_x = np.unique(x)
            unique_y = np.unique(y)
            dx = np.diff(unique_x).mean() if len(unique_x) > 1 else (x_max - x_min)
            dy = np.diff(unique_y).mean() if len(unique_y) > 1 else (y_max - y_min)
            
            # Use FiPy's built-in boundary face masks
            # Apply constraints based on which boundaries have marked cells
            # Only constrain a boundary if marked cells are predominantly at that boundary
            # (not corner cells that happen to be at multiple boundaries)
            
            # Check left boundary (x ~ x_min)
            # Only constrain if cells are specifically on left edge (not corners)
            left_cells = mask & (np.abs(x - x_min) < dx * 0.6)
            left_interior = left_cells & (y > y_min + dy * 0.6) & (y < y_max - dy * 0.6)
            if np.sum(left_cells) > 0 and np.sum(left_cells) >= 2:  # At least 2 cells marked
                # Use average of all left boundary cells (including corners)
                left_value = np.mean(bc_values[left_cells])
                self._solution.constrain(left_value, where=self._fipy_mesh.facesLeft)
            
            # Check right boundary (x ~ x_max)
            right_cells = mask & (np.abs(x - x_max) < dx * 0.6)
            if np.sum(right_cells) > 0 and np.sum(right_cells) >= 2:
                right_value = np.mean(bc_values[right_cells])
                self._solution.constrain(right_value, where=self._fipy_mesh.facesRight)
            
            # Check bottom boundary (y ~ y_min)
            # Only apply if there are interior bottom cells (not just corners)
            bottom_cells = mask & (np.abs(y - y_min) < dy * 0.6)
            bottom_interior = bottom_cells & (x > x_min + dx * 0.6) & (x < x_max - dx * 0.6)
            if np.sum(bottom_interior) > 0:  # Has interior bottom cells
                bottom_value = np.mean(bc_values[bottom_cells])
                self._solution.constrain(bottom_value, where=self._fipy_mesh.facesBottom)
            
            # Check top boundary (y ~ y_max)
            top_cells = mask & (np.abs(y - y_max) < dy * 0.6)
            top_interior = top_cells & (x > x_min + dx * 0.6) & (x < x_max - dx * 0.6)
            if np.sum(top_interior) > 0:  # Has interior top cells
                top_value = np.mean(bc_values[top_cells])
                self._solution.constrain(top_value, where=self._fipy_mesh.facesTop)
        
        # Build the equation
        # -div(A * grad(u)) + D*u = Y
        
        # Diffusion term - handle anisotropic tensors properly
        # FiPy's DiffusionTerm represents +div(coeff*grad(u)), NOT -div(coeff*grad(u))
        # For anisotropic diffusion, FiPy expects a FaceVariable with shape (dim, nFaces)
        if isinstance(self._A, list):
            a_tensor = np.array(self._A)
            if a_tensor.shape == (2, 2):
                # For diagonal anisotropy, use FaceVariable with rank=1
                # Each row contains the diagonal coefficient for that direction
                # For now, only support diagonal tensors (off-diagonal = 0)
                if np.abs(a_tensor[0, 1]) < 1e-14 and np.abs(a_tensor[1, 0]) < 1e-14:
                    # Diagonal anisotropy - use FaceVariable
                    k_diag = FaceVariable(mesh=self._fipy_mesh, 
                                         value=[[a_tensor[0, 0]], [a_tensor[1, 1]]], 
                                         rank=1)
                    diffusion = DiffusionTerm(coeff=k_diag)
                else:
                    # Full tensor - FiPy may not support this with constrain() BCs
                    # Fall back to scalar (isotropic) approximation with warning
                    import warnings
                    warnings.warn("FiPy backend: off-diagonal tensor terms not fully supported, "
                                  "using average of diagonal terms")
                    avg_k = 0.5 * (a_tensor[0, 0] + a_tensor[1, 1])
                    diffusion = DiffusionTerm(coeff=avg_k)
            else:
                diffusion = DiffusionTerm(coeff=self._A)
        elif isinstance(self._A, np.ndarray) and self._A.ndim == 2 and self._A.shape == (2, 2):
            # 2D array representing a tensor
            if np.abs(self._A[0, 1]) < 1e-14 and np.abs(self._A[1, 0]) < 1e-14:
                k_diag = FaceVariable(mesh=self._fipy_mesh, 
                                     value=[[self._A[0, 0]], [self._A[1, 1]]], 
                                     rank=1)
                diffusion = DiffusionTerm(coeff=k_diag)
            else:
                import warnings
                warnings.warn("FiPy backend: off-diagonal tensor terms not fully supported, "
                              "using average of diagonal terms")
                avg_k = 0.5 * (self._A[0, 0] + self._A[1, 1])
                diffusion = DiffusionTerm(coeff=avg_k)
        else:
            diffusion = DiffusionTerm(coeff=self._A)
        
        # Source terms
        eq = diffusion
        
        if self._D != 0:
            eq = eq + ImplicitSourceTerm(coeff=self._D)
        
        # Advection term (if present)
        if self._C is not None:
            eq = eq + ConvectionTerm(coeff=self._C)
        
        # RHS source term
        # Our model equation is: -div(A*grad(u)) + D*u = Y
        # FiPy's DiffusionTerm represents: +div(A*grad(u))
        # So the equation we're solving is: div(A*grad(u)) + D*u = -Y
        # Therefore we negate Y when adding to the equation
        if isinstance(self._Y, np.ndarray):
            source_var = CellVariable(mesh=self._fipy_mesh, value=self._Y)
            eq = eq + source_var  # +Y because FiPy uses opposite sign convention
        elif self._Y != 0:
            eq = eq + self._Y  # +Y because FiPy uses opposite sign convention
        
        # Choose solver — prefer PETSc (MPI-parallel) when available, fall back to
        # the default SciPy-based solvers for serial / non-PETSc environments.
        if PETSC_AVAILABLE:
            if self._solver_method.upper() == 'PCG' and self._symmetric:
                solver = PETScPCGSolver(tolerance=self._tolerance,
                                        iterations=self._max_iterations)
            else:
                solver = PETScGMRESSolver(tolerance=self._tolerance,
                                          iterations=self._max_iterations)
        else:
            if self._solver_method.upper() == 'PCG' and self._symmetric:
                solver = LinearPCGSolver(tolerance=self._tolerance,
                                        iterations=self._max_iterations)
            else:
                solver = LinearGMRESSolver(tolerance=self._tolerance,
                                           iterations=self._max_iterations)
        
        # Solve
        eq.solve(var=self._solution, solver=solver)
        
        return FiPyField(self._solution.copy(), self._fipy_mesh, "solution")
    
    def get_flux(self) -> FiPyField:
        """Get the flux from the solution."""
        # FiPy computes flux at faces
        # flux = -A * grad(u)
        grad_u = self._solution.grad
        
        if isinstance(self._A, list):
            # Anisotropic case
            A_tensor = np.array(self._A[0])
            flux_x = -(A_tensor[0, 0] * grad_u[0] + A_tensor[0, 1] * grad_u[1])
            flux_y = -(A_tensor[1, 0] * grad_u[0] + A_tensor[1, 1] * grad_u[1])
            # Return as cell-centered values
            flux = np.array([flux_x.cellVolumeAverage, flux_y.cellVolumeAverage])
        else:
            flux = -self._A * grad_u
        
        return FiPyField(flux, self._fipy_mesh, "flux")
    
    def set_solver_options(self,
                          method: str = 'GMRES',
                          tolerance: float = 1e-8,
                          max_iterations: int = 10000) -> None:
        """Set solver options."""
        self._solver_method = method
        self._tolerance = tolerance
        self._max_iterations = max_iterations


class FiPyBackend(BackendBase):
    """
    FiPy backend for GroMPy-couple.
    
    Provides the FiPy (Finite Volume Method) implementation of all
    backend operations for solving coupled groundwater flow and
    solute transport equations.
    """
    
    @property
    def name(self) -> str:
        """Return backend name."""
        return 'fipy'
    
    # =========================================================================
    # Mesh creation
    # =========================================================================
    
    def load_mesh(self, filename: str) -> FiPyMesh:
        """Load mesh from file."""
        from lib.fipy_mesh_io import load_fipy_mesh_from_msh
        
        fipy_mesh, cell_centers, masks, field_data, extra = \
            load_fipy_mesh_from_msh(filename)
        
        mesh = FiPyMesh(fipy_mesh)
        # Store additional info
        mesh._masks = masks
        mesh._field_data = field_data
        mesh._extra = extra
        
        return mesh
    
    def create_rectangle_mesh(self,
                             length: float,
                             height: float,
                             nx: int,
                             ny: int) -> FiPyMesh:
        """Create rectangular mesh."""
        dx = length / nx
        dy = height / ny
        mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)
        return FiPyMesh(mesh)
    
    # =========================================================================
    # Field creation
    # =========================================================================
    
    def create_scalar_field(self,
                           mesh: FiPyMesh,
                           value: Union[float, np.ndarray] = 0.0,
                           name: str = "",
                           on_boundary: bool = False) -> FiPyField:
        """Create scalar field."""
        fipy_mesh = mesh.fipy_mesh
        
        if isinstance(value, np.ndarray):
            var = CellVariable(mesh=fipy_mesh, value=value, name=name)
        else:
            var = CellVariable(mesh=fipy_mesh, value=value, name=name)
        
        return FiPyField(var, fipy_mesh, name)
    
    def create_vector_field(self,
                           mesh: FiPyMesh,
                           values: Union[Tuple[float, float], np.ndarray] = (0.0, 0.0),
                           name: str = "") -> FiPyField:
        """Create vector field."""
        fipy_mesh = mesh.fipy_mesh
        n_cells = mesh.num_cells
        
        if isinstance(values, tuple):
            data = np.zeros((2, n_cells))
            data[0, :] = values[0]
            data[1, :] = values[1]
        else:
            data = values
        
        # Create as a CellVariable with rank 1
        var = CellVariable(mesh=fipy_mesh, value=data, name=name)
        return FiPyField(var, fipy_mesh, name)
    
    def create_tensor_field(self,
                           mesh: FiPyMesh,
                           values: Union[np.ndarray, Tuple[Tuple[float, float],
                                                           Tuple[float, float]]] = None,
                           name: str = "") -> FiPyField:
        """Create 2x2 tensor field."""
        if values is None:
            values = ((0, 0), (0, 0))
        
        # Store tensor as nested array
        tensor = np.array(values, dtype=float)
        
        # For FiPy, tensors are typically stored as arrays and used
        # in the diffusion coefficient specification
        return FiPyField(tensor, mesh.fipy_mesh, name)
    
    # =========================================================================
    # PDE solver creation
    # =========================================================================
    
    def create_pde_solver(self, mesh: FiPyMesh, symmetric: bool = True) -> FiPyPDESolver:
        """Create PDE solver."""
        return FiPyPDESolver(mesh, symmetric)
    
    def create_projection_solver(self, mesh: FiPyMesh) -> FiPyPDESolver:
        """
        Create projection solver.
        
        Note: FiPy is cell-centered so projection is typically not needed.
        We return a standard solver that can be used if needed.
        """
        return FiPyPDESolver(mesh, symmetric=True)
    
    # =========================================================================
    # Mathematical operations
    # =========================================================================
    
    def gradient(self, field: FiPyField) -> FiPyField:
        """Compute gradient."""
        if hasattr(field.variable, 'grad'):
            grad = field.variable.grad
            return FiPyField(grad, field.mesh, "gradient")
        else:
            # For arrays, compute numerical gradient
            raise NotImplementedError("Gradient for raw arrays not implemented")
    
    def divergence(self, field: FiPyField) -> FiPyField:
        """Compute divergence."""
        if hasattr(field.variable, 'divergence'):
            div = field.variable.divergence
            return FiPyField(div, field.mesh, "divergence")
        else:
            raise NotImplementedError("Divergence for raw arrays not implemented")
    
    def integrate(self, field: FiPyField, over_boundary: bool = False) -> float:
        """Integrate field over domain."""
        if over_boundary:
            # FiPy doesn't have a direct boundary integral
            # This would need to sum over boundary faces
            raise NotImplementedError("Boundary integration not yet implemented")
        
        values = field.get_values()
        if field.mesh is not None:
            # Get cell volumes
            volumes = field.mesh.cellVolumes
            return float(np.sum(values * volumes))
        else:
            return float(np.sum(values))
    
    def interpolate(self, field: FiPyField, target_space: str) -> FiPyField:
        """
        Interpolate field to different function space.
        
        Note: FiPy is primarily cell-centered, so this is mostly a no-op.
        """
        # FiPy doesn't have the same function space concept as escript
        # Return a copy for now
        return field.copy()
    
    # =========================================================================
    # Conditional operations
    # =========================================================================
    
    def where_positive(self, field: FiPyField) -> FiPyField:
        """Return 1 where field > 0, 0 otherwise."""
        values = field.get_values()
        mask = (values > 0).astype(float)
        return FiPyField(mask, field.mesh, "mask")
    
    def where_negative(self, field: FiPyField) -> FiPyField:
        """Return 1 where field < 0, 0 otherwise."""
        values = field.get_values()
        mask = (values < 0).astype(float)
        return FiPyField(mask, field.mesh, "mask")
    
    def where_zero(self, field: FiPyField, tolerance: float = 1e-10) -> FiPyField:
        """Return 1 where field == 0 (within tolerance), 0 otherwise."""
        values = field.get_values()
        mask = (np.abs(values) < tolerance).astype(float)
        return FiPyField(mask, field.mesh, "mask")
    
    # =========================================================================
    # Reduction operations
    # =========================================================================
    
    def minimum(self, field: FiPyField) -> float:
        """Get minimum value."""
        return float(np.min(field.get_values()))
    
    def maximum(self, field: FiPyField) -> float:
        """Get maximum value."""
        return float(np.max(field.get_values()))
    
    def max_abs(self, field: FiPyField) -> float:
        """Get maximum absolute value."""
        return float(np.max(np.abs(field.get_values())))
    
    # =========================================================================
    # I/O operations
    # =========================================================================
    
    def save_vtk(self,
                filename: str,
                mesh: FiPyMesh,
                fields: dict) -> None:
        """Save to VTK format."""
        try:
            from fipy.viewers.vtkViewer import VTKCellViewer
            
            # Prepare variables for VTK export
            vars_to_save = []
            for name, field in fields.items():
                if isinstance(field, FiPyField):
                    var = field.variable
                    if hasattr(var, 'name'):
                        var.name = name
                    vars_to_save.append(var)
            
            if vars_to_save:
                viewer = VTKCellViewer(vars=vars_to_save)
                viewer.plot(filename=filename)
        except ImportError:
            # Fallback: use meshio to write VTK
            try:
                import meshio
                
                # Get mesh data
                fipy_mesh = mesh.fipy_mesh
                
                # Get vertices and cells
                try:
                    verts = fipy_mesh.vertexCoords.T
                    cells = fipy_mesh._cellVertexIDs.T
                    
                    # Determine cell type
                    if cells.shape[1] == 3:
                        cell_type = 'triangle'
                    elif cells.shape[1] == 4:
                        cell_type = 'quad'
                    else:
                        cell_type = 'polygon'
                    
                    # Create meshio mesh
                    points = np.column_stack([verts, np.zeros(len(verts))])
                    cells_meshio = [(cell_type, cells)]
                    
                    # Prepare cell data
                    cell_data = {}
                    for name, field in fields.items():
                        if isinstance(field, FiPyField):
                            values = field.get_values()
                            cell_data[name] = [values]
                    
                    mesh_obj = meshio.Mesh(points, cells_meshio, cell_data=cell_data)
                    mesh_obj.write(f"{filename}.vtk")
                except Exception as e:
                    print(f"Warning: Could not save VTK file: {e}")
            except ImportError:
                print("Warning: Neither FiPy VTK viewer nor meshio available for VTK output")
    
    # =========================================================================
    # Utility operations
    # =========================================================================
    
    def kronecker(self, mesh: FiPyMesh) -> FiPyField:
        """Return Kronecker delta (identity tensor)."""
        identity = np.array([[1.0, 0.0], [0.0, 1.0]])
        return FiPyField(identity, mesh.fipy_mesh, "kronecker")
    
    def get_domain_coordinates(self, mesh: FiPyMesh) -> FiPyField:
        """Get coordinate field."""
        centers = mesh.fipy_mesh.cellCenters
        coords = np.array([centers[0], centers[1]])
        return FiPyField(coords, mesh.fipy_mesh, "coordinates")
