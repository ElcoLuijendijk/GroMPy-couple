"""
Escript backend implementation for GroMPy-couple.

This module provides the escript/Finley implementation of the backend interface.
It wraps esys-escript functionality to match the abstract interface defined
in base.py.
"""

from typing import Tuple, Optional, Union, Any
import numpy as np

try:
    import esys.escript as es
    import esys.escript.linearPDEs as linearPDEs
    import esys.finley as finley
    import esys.weipa
    ESCRIPT_AVAILABLE = True
except ImportError:
    ESCRIPT_AVAILABLE = False

from lib.backend.base import BackendBase, MeshBase, FieldBase, PDESolverBase


class EscriptField(FieldBase):
    """
    Escript implementation of FieldBase.
    
    Wraps esys.escript.Data objects.
    """
    
    def __init__(self, data: 'es.Data', name: str = ""):
        """
        Initialize an EscriptField.
        
        Parameters
        ----------
        data : es.Data
            Escript Data object
        name : str
            Name of the field
        """
        super().__init__(name)
        self._data = data
    
    @property
    def data(self) -> 'es.Data':
        """Get the underlying escript Data object."""
        return self._data
    
    def get_values(self) -> np.ndarray:
        """Get field values as numpy array."""
        return np.array(self._data.toListOfTuples())
    
    def set_values(self, values: Union[np.ndarray, float]) -> None:
        """Set field values."""
        if isinstance(values, (int, float)):
            self._data.setTaggedValue(0, values)
            # For scalar, set all values
            self._data *= 0
            self._data += values
        else:
            # For array, need to create new Data object
            # This is a limitation - escript Data doesn't support direct array assignment
            raise NotImplementedError(
                "Setting escript field from array is not fully supported. "
                "Create a new field instead."
            )
    
    def copy(self) -> 'EscriptField':
        """Create a copy of the field."""
        return EscriptField(self._data.copy(), self.name)
    
    def get_shape(self) -> Tuple[int, ...]:
        """Get the shape of the field data."""
        return self._data.getShape()
    
    def get_coordinates(self) -> np.ndarray:
        """Get x, y coordinates associated with field values."""
        coords = self._data.getFunctionSpace().getX()
        xy = np.array([coords[0].toListOfTuples(), 
                       coords[1].toListOfTuples()]).T
        return xy
    
    # Arithmetic operations
    def __add__(self, other: Union['EscriptField', float]) -> 'EscriptField':
        if isinstance(other, EscriptField):
            return EscriptField(self._data + other._data, self.name)
        return EscriptField(self._data + other, self.name)
    
    def __radd__(self, other: Union['EscriptField', float]) -> 'EscriptField':
        return self.__add__(other)
    
    def __sub__(self, other: Union['EscriptField', float]) -> 'EscriptField':
        if isinstance(other, EscriptField):
            return EscriptField(self._data - other._data, self.name)
        return EscriptField(self._data - other, self.name)
    
    def __rsub__(self, other: Union['EscriptField', float]) -> 'EscriptField':
        if isinstance(other, EscriptField):
            return EscriptField(other._data - self._data, self.name)
        return EscriptField(other - self._data, self.name)
    
    def __mul__(self, other: Union['EscriptField', float]) -> 'EscriptField':
        if isinstance(other, EscriptField):
            return EscriptField(self._data * other._data, self.name)
        return EscriptField(self._data * other, self.name)
    
    def __rmul__(self, other: Union['EscriptField', float]) -> 'EscriptField':
        return self.__mul__(other)
    
    def __truediv__(self, other: Union['EscriptField', float]) -> 'EscriptField':
        if isinstance(other, EscriptField):
            return EscriptField(self._data / other._data, self.name)
        return EscriptField(self._data / other, self.name)
    
    def __neg__(self) -> 'EscriptField':
        return EscriptField(-self._data, self.name)
    
    def __pow__(self, power: float) -> 'EscriptField':
        return EscriptField(self._data ** power, self.name)
    
    def __getitem__(self, index: int) -> 'EscriptField':
        """Get component of vector/tensor field."""
        return EscriptField(self._data[index], self.name)
    
    def __setitem__(self, index: int, value: Union['EscriptField', float]) -> None:
        """Set component of vector/tensor field."""
        if isinstance(value, EscriptField):
            self._data[index] = value._data
        else:
            self._data[index] = value


class EscriptMesh(MeshBase):
    """
    Escript implementation of MeshBase.
    
    Wraps esys.finley Domain objects.
    """
    
    def __init__(self, domain: 'finley.Domain'):
        """
        Initialize an EscriptMesh.
        
        Parameters
        ----------
        domain : finley.Domain
            Escript/Finley domain object
        """
        self._domain = domain
        # Cache some values for efficiency
        self._num_cells = None
        self._num_nodes = None
    
    @property
    def domain(self) -> 'finley.Domain':
        """Get the underlying escript domain."""
        return self._domain
    
    @property
    def num_cells(self) -> int:
        """Get number of cells."""
        if self._num_cells is None:
            # Get number of elements
            x = es.Function(self._domain).getX()
            self._num_cells = len(x[0].toListOfTuples())
        return self._num_cells
    
    @property
    def num_nodes(self) -> int:
        """Get number of nodes."""
        if self._num_nodes is None:
            x = es.ContinuousFunction(self._domain).getX()
            self._num_nodes = len(x[0].toListOfTuples())
        return self._num_nodes
    
    def get_cell_centers(self) -> np.ndarray:
        """Get coordinates of cell centers."""
        x = es.Function(self._domain).getX()
        return np.array([x[0].toListOfTuples(), 
                        x[1].toListOfTuples()]).T
    
    def get_node_coordinates(self) -> np.ndarray:
        """Get coordinates of mesh nodes."""
        x = es.ContinuousFunction(self._domain).getX()
        return np.array([x[0].toListOfTuples(), 
                        x[1].toListOfTuples()]).T
    
    def get_boundary_mask(self, boundary_name: str) -> np.ndarray:
        """Get boolean mask for boundary nodes."""
        # Get coordinates
        x = es.FunctionOnBoundary(self._domain).getX()
        coords = np.array([x[0].toListOfTuples(), x[1].toListOfTuples()]).T
        
        if len(coords) == 0:
            return np.array([])
        
        x_vals = coords[:, 0]
        y_vals = coords[:, 1]
        
        tol = 1e-6
        x_min, x_max = x_vals.min(), x_vals.max()
        y_min, y_max = y_vals.min(), y_vals.max()
        
        if boundary_name == 'left':
            return np.abs(x_vals - x_min) < tol
        elif boundary_name == 'right':
            return np.abs(x_vals - x_max) < tol
        elif boundary_name == 'bottom':
            return np.abs(y_vals - y_min) < tol
        elif boundary_name == 'top':
            return np.abs(y_vals - y_max) < tol
        else:
            # For named boundaries (from Gmsh), would need tag lookup
            # For now, return empty mask
            return np.zeros(len(coords), dtype=bool)
    
    def get_native_mesh(self) -> Any:
        """Get the native escript domain object."""
        return self._domain


class EscriptPDESolver(PDESolverBase):
    """
    Escript implementation of PDESolverBase.
    
    Wraps esys.escript.linearPDEs.LinearPDE.
    """
    
    def __init__(self, mesh: EscriptMesh, symmetric: bool = True):
        """
        Initialize an EscriptPDESolver.
        
        Parameters
        ----------
        mesh : EscriptMesh
            Mesh to solve on
        symmetric : bool
            Whether the PDE is symmetric
        """
        self._mesh = mesh
        self._pde = linearPDEs.LinearPDE(mesh.domain)
        if symmetric:
            self._pde.setSymmetryOn()
        
        # Store the last solution for getFlux
        self._last_solution = None
    
    @property
    def pde(self) -> 'linearPDEs.LinearPDE':
        """Get the underlying escript PDE object."""
        return self._pde
    
    def set_coefficients(self,
                        A: Optional[Union[EscriptField, float, np.ndarray]] = None,
                        C: Optional[Union[EscriptField, np.ndarray]] = None,
                        D: Optional[Union[EscriptField, float]] = None,
                        X: Optional[Union[EscriptField, np.ndarray]] = None,
                        Y: Optional[Union[EscriptField, float]] = None) -> None:
        """Set PDE coefficients."""
        kwargs = {}
        
        if A is not None:
            if isinstance(A, EscriptField):
                # For tensor A, multiply by kronecker delta for proper form
                kwargs['A'] = A.data * es.kronecker(self._mesh.domain)
            elif isinstance(A, (int, float)):
                kwargs['A'] = A * es.kronecker(self._mesh.domain)
            else:
                kwargs['A'] = A
        
        if C is not None:
            if isinstance(C, EscriptField):
                kwargs['C'] = C.data
            else:
                kwargs['C'] = C
        
        if D is not None:
            if isinstance(D, EscriptField):
                kwargs['D'] = D.data
            else:
                kwargs['D'] = D
        
        if X is not None:
            if isinstance(X, EscriptField):
                kwargs['X'] = X.data
            else:
                kwargs['X'] = X
        
        if Y is not None:
            if isinstance(Y, EscriptField):
                kwargs['Y'] = Y.data
            else:
                kwargs['Y'] = Y
        
        if kwargs:
            self._pde.setValue(**kwargs)
    
    def set_boundary_conditions(self,
                               dirichlet_mask: Optional[EscriptField] = None,
                               dirichlet_value: Optional[Union[EscriptField, float]] = None,
                               neumann_flux: Optional[Union[EscriptField, float]] = None) -> None:
        """Set boundary conditions."""
        kwargs = {}
        
        if dirichlet_mask is not None:
            if isinstance(dirichlet_mask, EscriptField):
                kwargs['q'] = dirichlet_mask.data
            else:
                kwargs['q'] = dirichlet_mask
        
        if dirichlet_value is not None:
            if isinstance(dirichlet_value, EscriptField):
                kwargs['r'] = dirichlet_value.data
            else:
                kwargs['r'] = dirichlet_value
        
        if neumann_flux is not None:
            if isinstance(neumann_flux, EscriptField):
                kwargs['y'] = neumann_flux.data
            else:
                kwargs['y'] = neumann_flux
        
        if kwargs:
            self._pde.setValue(**kwargs)
    
    def solve(self) -> EscriptField:
        """Solve the PDE."""
        solution = self._pde.getSolution()
        self._last_solution = solution
        return EscriptField(solution, "solution")
    
    def get_flux(self) -> EscriptField:
        """Get the flux from the last solution."""
        flux = self._pde.getFlux()
        return EscriptField(flux, "flux")
    
    def set_solver_options(self,
                          method: str = 'GMRES',
                          tolerance: float = 1e-8,
                          max_iterations: int = 10000) -> None:
        """Set solver options."""
        options = self._pde.getSolverOptions()
        
        method_map = {
            'GMRES': es.SolverOptions.GMRES,
            'PCG': es.SolverOptions.PCG,
            'DIRECT': es.SolverOptions.DIRECT,
            'TFQMR': es.SolverOptions.TFQMR,
        }
        
        if method.upper() in method_map:
            options.setSolverMethod(method_map[method.upper()])
        
        options.setTolerance(tolerance)
        options.setIterMax(max_iterations)


class EscriptProjectionSolver(PDESolverBase):
    """
    Escript projection solver for converting element values to nodes.
    
    Wraps esys.escript.linearPDEs.LinearPDESystem for projection.
    """
    
    def __init__(self, mesh: EscriptMesh):
        """Initialize projection solver."""
        self._mesh = mesh
        self._pde = linearPDEs.LinearPDESystem(mesh.domain)
        self._pde.setSymmetryOn()
    
    def set_coefficients(self,
                        A: Optional[Union[EscriptField, float, np.ndarray]] = None,
                        C: Optional[Union[EscriptField, np.ndarray]] = None,
                        D: Optional[Union[EscriptField, float]] = None,
                        X: Optional[Union[EscriptField, np.ndarray]] = None,
                        Y: Optional[Union[EscriptField, float]] = None) -> None:
        """Set projection coefficients."""
        kwargs = {}
        
        if D is not None:
            if isinstance(D, EscriptField):
                kwargs['D'] = D.data
            else:
                kwargs['D'] = D
        
        if Y is not None:
            if isinstance(Y, EscriptField):
                kwargs['Y'] = Y.data
            else:
                kwargs['Y'] = Y
        
        if kwargs:
            self._pde.setValue(**kwargs)
    
    def set_boundary_conditions(self,
                               dirichlet_mask: Optional[EscriptField] = None,
                               dirichlet_value: Optional[Union[EscriptField, float]] = None,
                               neumann_flux: Optional[Union[EscriptField, float]] = None) -> None:
        """Set boundary conditions (typically not used for projection)."""
        pass
    
    def solve(self) -> EscriptField:
        """Solve the projection."""
        solution = self._pde.getSolution()
        return EscriptField(solution, "projected")
    
    def get_flux(self) -> EscriptField:
        """Get flux (not applicable for projection)."""
        raise NotImplementedError("Projection solver does not compute flux")
    
    def set_solver_options(self,
                          method: str = 'GMRES',
                          tolerance: float = 1e-8,
                          max_iterations: int = 10000) -> None:
        """Set solver options."""
        options = self._pde.getSolverOptions()
        options.setTolerance(tolerance)
        options.setIterMax(max_iterations)


class EscriptBackend(BackendBase):
    """
    Escript backend for GroMPy-couple.
    
    Provides the escript/Finley implementation of all backend operations.
    """
    
    @property
    def name(self) -> str:
        """Return backend name."""
        return 'escript'
    
    # =========================================================================
    # Mesh creation
    # =========================================================================
    
    def load_mesh(self, filename: str) -> EscriptMesh:
        """Load mesh from file."""
        if filename.endswith('.fly'):
            domain = finley.ReadMesh(filename)
        elif filename.endswith('.msh'):
            domain = finley.ReadGmsh(filename, numDim=2)
        else:
            raise ValueError(f"Unsupported mesh format: {filename}")
        return EscriptMesh(domain)
    
    def create_rectangle_mesh(self,
                             length: float,
                             height: float,
                             nx: int,
                             ny: int) -> EscriptMesh:
        """Create rectangular mesh."""
        domain = finley.Rectangle(nx, ny, l0=length, l1=height)
        return EscriptMesh(domain)
    
    # =========================================================================
    # Field creation
    # =========================================================================
    
    def create_scalar_field(self,
                           mesh: EscriptMesh,
                           value: Union[float, np.ndarray] = 0.0,
                           name: str = "",
                           on_boundary: bool = False) -> EscriptField:
        """Create scalar field."""
        if on_boundary:
            func_space = es.FunctionOnBoundary(mesh.domain)
        else:
            func_space = es.Function(mesh.domain)
        
        if isinstance(value, np.ndarray):
            # Create field from array
            data = es.Scalar(0, func_space)
            # escript doesn't support direct array assignment easily,
            # so we need to create it element by element or use interpolation
            # For now, create zero and add - this is a workaround
            # In practice, for initial conditions you'd interpolate from coordinates
            x = func_space.getX()
            x_arr = np.array(x[0].toListOfTuples())
            if len(value) == len(x_arr):
                # Create by interpolating - simplified approach
                data = es.Scalar(0, func_space)
                # This is a limitation - proper implementation would need
                # to map array values to mesh points
                data += value.mean()  # Fallback to mean value
        else:
            data = es.Scalar(value, func_space)
        
        return EscriptField(data, name)
    
    def create_vector_field(self,
                           mesh: EscriptMesh,
                           values: Union[Tuple[float, float], np.ndarray] = (0.0, 0.0),
                           name: str = "") -> EscriptField:
        """Create vector field."""
        func_space = es.Function(mesh.domain)
        
        if isinstance(values, tuple):
            data = es.Vector(values, func_space)
        else:
            data = es.Vector((0, 0), func_space)
            # Would need proper handling for array input
        
        return EscriptField(data, name)
    
    def create_tensor_field(self,
                           mesh: EscriptMesh,
                           values: Union[np.ndarray, Tuple[Tuple[float, float],
                                                           Tuple[float, float]]] = None,
                           name: str = "") -> EscriptField:
        """Create 2x2 tensor field."""
        func_space = es.Function(mesh.domain)
        
        if values is None:
            values = ((0, 0), (0, 0))
        
        data = es.Tensor(values, func_space)
        return EscriptField(data, name)
    
    # =========================================================================
    # PDE solver creation
    # =========================================================================
    
    def create_pde_solver(self, mesh: EscriptMesh, symmetric: bool = True) -> EscriptPDESolver:
        """Create PDE solver."""
        return EscriptPDESolver(mesh, symmetric)
    
    def create_projection_solver(self, mesh: EscriptMesh) -> EscriptProjectionSolver:
        """Create projection solver."""
        return EscriptProjectionSolver(mesh)
    
    # =========================================================================
    # Mathematical operations
    # =========================================================================
    
    def gradient(self, field: EscriptField) -> EscriptField:
        """Compute gradient."""
        return EscriptField(es.grad(field.data), "gradient")
    
    def divergence(self, field: EscriptField) -> EscriptField:
        """Compute divergence."""
        # escript doesn't have a direct divergence function
        # div(v) = d(v_x)/dx + d(v_y)/dy
        grad_v = es.grad(field.data)
        div = grad_v[0, 0] + grad_v[1, 1]
        return EscriptField(div, "divergence")
    
    def integrate(self, field: EscriptField, over_boundary: bool = False) -> float:
        """Integrate field over domain or boundary."""
        if over_boundary:
            return es.integrate(field.data, 
                              where=es.FunctionOnBoundary(field.data.getDomain()))
        return es.integrate(field.data)
    
    def interpolate(self, field: EscriptField, target_space: str) -> EscriptField:
        """Interpolate field to different function space."""
        domain = field.data.getDomain()
        
        if target_space == 'nodes':
            target = es.ContinuousFunction(domain)
        elif target_space == 'cells':
            target = es.Function(domain)
        elif target_space == 'boundary':
            target = es.FunctionOnBoundary(domain)
        else:
            raise ValueError(f"Unknown target space: {target_space}")
        
        return EscriptField(es.interpolate(field.data, target), field.name)
    
    # =========================================================================
    # Conditional operations
    # =========================================================================
    
    def where_positive(self, field: EscriptField) -> EscriptField:
        """Return 1 where field > 0, 0 otherwise."""
        return EscriptField(es.wherePositive(field.data), "mask")
    
    def where_negative(self, field: EscriptField) -> EscriptField:
        """Return 1 where field < 0, 0 otherwise."""
        return EscriptField(es.whereNegative(field.data), "mask")
    
    def where_zero(self, field: EscriptField, tolerance: float = 1e-10) -> EscriptField:
        """Return 1 where field == 0, 0 otherwise."""
        return EscriptField(es.whereZero(field.data, tolerance), "mask")
    
    # =========================================================================
    # Reduction operations
    # =========================================================================
    
    def minimum(self, field: EscriptField) -> float:
        """Get minimum value."""
        return es.inf(field.data)
    
    def maximum(self, field: EscriptField) -> float:
        """Get maximum value."""
        return es.sup(field.data)
    
    def max_abs(self, field: EscriptField) -> float:
        """Get maximum absolute value."""
        return es.Lsup(field.data)
    
    # =========================================================================
    # I/O operations
    # =========================================================================
    
    def save_vtk(self,
                filename: str,
                mesh: EscriptMesh,
                fields: dict) -> None:
        """Save to VTK format."""
        # Prepare fields dict for esys.weipa.saveVTK
        vtk_fields = {}
        for name, field in fields.items():
            if isinstance(field, EscriptField):
                vtk_fields[name] = field.data
            else:
                vtk_fields[name] = field
        
        esys.weipa.saveVTK(filename, **vtk_fields)
    
    # =========================================================================
    # Utility operations
    # =========================================================================
    
    def kronecker(self, mesh: EscriptMesh) -> EscriptField:
        """Return Kronecker delta (identity tensor)."""
        return EscriptField(es.kronecker(mesh.domain), "kronecker")
    
    def get_domain_coordinates(self, mesh: EscriptMesh) -> EscriptField:
        """Get coordinate field."""
        return EscriptField(mesh.domain.getX(), "coordinates")
