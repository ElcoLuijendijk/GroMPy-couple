"""
Abstract base classes for GroMPy-couple numerical backends.

This module defines the interface that all numerical backends must implement.
The interface is designed to abstract the key operations needed for solving
coupled groundwater flow and solute transport equations.

Classes
-------
BackendBase : Main backend interface
MeshBase : Mesh abstraction
FieldBase : Field/variable abstraction
PDESolverBase : PDE solver abstraction
"""

from abc import ABC, abstractmethod
from typing import Tuple, Optional, Union, Any
import numpy as np


class FieldBase(ABC):
    """
    Abstract base class for scalar and vector fields.
    
    A field represents a numerical variable defined on a mesh, such as
    pressure, concentration, or velocity. This class abstracts the
    differences between escript Data objects and FiPy CellVariable/FaceVariable.
    
    Attributes
    ----------
    name : str
        Name of the field for identification
    """
    
    def __init__(self, name: str = ""):
        self.name = name
    
    @abstractmethod
    def get_values(self) -> np.ndarray:
        """
        Get field values as a numpy array.
        
        Returns
        -------
        np.ndarray
            Field values
        """
        pass
    
    @abstractmethod
    def set_values(self, values: Union[np.ndarray, float]) -> None:
        """
        Set field values from a numpy array or scalar.
        
        Parameters
        ----------
        values : np.ndarray or float
            Values to set
        """
        pass
    
    @abstractmethod
    def copy(self) -> 'FieldBase':
        """
        Create a copy of the field.
        
        Returns
        -------
        FieldBase
            Copy of the field
        """
        pass
    
    @abstractmethod
    def get_shape(self) -> Tuple[int, ...]:
        """
        Get the shape of the field data.
        
        Returns
        -------
        tuple
            Shape of the field data
        """
        pass
    
    @abstractmethod
    def get_coordinates(self) -> np.ndarray:
        """
        Get the x, y coordinates associated with the field values.
        
        Returns
        -------
        np.ndarray
            Array of shape (N, 2) with x, y coordinates
        """
        pass
    
    # Arithmetic operations (all backends must support these)
    @abstractmethod
    def __add__(self, other: Union['FieldBase', float]) -> 'FieldBase':
        pass
    
    @abstractmethod
    def __radd__(self, other: Union['FieldBase', float]) -> 'FieldBase':
        pass
    
    @abstractmethod
    def __sub__(self, other: Union['FieldBase', float]) -> 'FieldBase':
        pass
    
    @abstractmethod
    def __rsub__(self, other: Union['FieldBase', float]) -> 'FieldBase':
        pass
    
    @abstractmethod
    def __mul__(self, other: Union['FieldBase', float]) -> 'FieldBase':
        pass
    
    @abstractmethod
    def __rmul__(self, other: Union['FieldBase', float]) -> 'FieldBase':
        pass
    
    @abstractmethod
    def __truediv__(self, other: Union['FieldBase', float]) -> 'FieldBase':
        pass
    
    @abstractmethod
    def __neg__(self) -> 'FieldBase':
        pass
    
    @abstractmethod
    def __pow__(self, power: float) -> 'FieldBase':
        pass


class MeshBase(ABC):
    """
    Abstract base class for mesh objects.
    
    A mesh represents the spatial discretization of the model domain.
    This class abstracts the differences between escript/Finley meshes
    and FiPy meshes.
    """
    
    @property
    @abstractmethod
    def num_cells(self) -> int:
        """
        Get the number of cells in the mesh.
        
        Returns
        -------
        int
            Number of cells
        """
        pass
    
    @property
    @abstractmethod
    def num_nodes(self) -> int:
        """
        Get the number of nodes in the mesh.
        
        Returns
        -------
        int
            Number of nodes
        """
        pass
    
    @abstractmethod
    def get_cell_centers(self) -> np.ndarray:
        """
        Get the coordinates of cell centers.
        
        Returns
        -------
        np.ndarray
            Array of shape (N, 2) with x, y coordinates
        """
        pass
    
    @abstractmethod
    def get_node_coordinates(self) -> np.ndarray:
        """
        Get the coordinates of mesh nodes.
        
        Returns
        -------
        np.ndarray
            Array of shape (N, 2) with x, y coordinates
        """
        pass
    
    @abstractmethod
    def get_boundary_mask(self, boundary_name: str) -> np.ndarray:
        """
        Get a boolean mask identifying boundary cells/nodes.
        
        Parameters
        ----------
        boundary_name : str
            Name of the boundary (e.g., 'top', 'bottom', 'left', 'right',
            'sea_surface', 'land_surface')
        
        Returns
        -------
        np.ndarray
            Boolean mask
        """
        pass
    
    @abstractmethod
    def get_native_mesh(self) -> Any:
        """
        Get the underlying native mesh object.
        
        Returns
        -------
        Any
            Native mesh object (escript Domain or FiPy Mesh)
        """
        pass


class PDESolverBase(ABC):
    """
    Abstract base class for PDE solvers.
    
    This class abstracts the linear PDE solver interface. Both pressure
    and solute transport equations are solved using this interface.
    
    The general form of the PDE being solved is:
    
        -div(A * grad(u)) + D*u = Y - div(X)
    
    With boundary conditions:
        - Dirichlet: u = r where q = 1
        - Neumann: n * A * grad(u) = y on boundaries
    
    For advection-diffusion (solute transport):
        -div(A * grad(u)) + div(C * u) + D*u = Y
    """
    
    @abstractmethod
    def set_coefficients(self,
                        A: Optional[Union[FieldBase, float, np.ndarray]] = None,
                        C: Optional[Union[FieldBase, np.ndarray]] = None,
                        D: Optional[Union[FieldBase, float]] = None,
                        X: Optional[Union[FieldBase, np.ndarray]] = None,
                        Y: Optional[Union[FieldBase, float]] = None) -> None:
        """
        Set the PDE coefficients.
        
        Parameters
        ----------
        A : FieldBase, float, or np.ndarray, optional
            Diffusion coefficient tensor (can be scalar for isotropic,
            2x2 array for anisotropic, or spatially varying field)
        C : FieldBase or np.ndarray, optional
            Advection coefficient vector
        D : FieldBase or float, optional
            Reaction/accumulation coefficient
        X : FieldBase or np.ndarray, optional
            Source vector
        Y : FieldBase or float, optional
            Source scalar
        """
        pass
    
    @abstractmethod
    def set_boundary_conditions(self,
                               dirichlet_mask: Optional[FieldBase] = None,
                               dirichlet_value: Optional[Union[FieldBase, float]] = None,
                               neumann_flux: Optional[Union[FieldBase, float]] = None) -> None:
        """
        Set boundary conditions for the PDE.
        
        Parameters
        ----------
        dirichlet_mask : FieldBase, optional
            Mask indicating where Dirichlet conditions apply (q in escript)
        dirichlet_value : FieldBase or float, optional
            Values for Dirichlet boundary conditions (r in escript)
        neumann_flux : FieldBase or float, optional
            Flux values for Neumann boundary conditions (y in escript)
        """
        pass
    
    @abstractmethod
    def solve(self) -> FieldBase:
        """
        Solve the PDE and return the solution.
        
        Returns
        -------
        FieldBase
            Solution field
        """
        pass
    
    @abstractmethod
    def get_flux(self) -> FieldBase:
        """
        Get the flux computed from the solution: -A * grad(u).
        
        Returns
        -------
        FieldBase
            Flux vector field
        """
        pass
    
    @abstractmethod
    def set_solver_options(self, 
                          method: str = 'GMRES',
                          tolerance: float = 1e-8,
                          max_iterations: int = 10000) -> None:
        """
        Set solver options.
        
        Parameters
        ----------
        method : str
            Solver method: 'GMRES', 'PCG', 'DIRECT', 'TFQMR'
        tolerance : float
            Convergence tolerance
        max_iterations : int
            Maximum number of iterations
        """
        pass


class BackendBase(ABC):
    """
    Abstract base class for numerical backends.
    
    This is the main interface that users interact with. It provides
    factory methods for creating meshes, fields, and PDE solvers.
    
    Attributes
    ----------
    name : str
        Name of the backend ('escript' or 'fipy')
    """
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Return the name of the backend."""
        pass
    
    # =========================================================================
    # Mesh creation
    # =========================================================================
    
    @abstractmethod
    def load_mesh(self, filename: str) -> MeshBase:
        """
        Load a mesh from a file.
        
        Parameters
        ----------
        filename : str
            Path to mesh file (Gmsh .msh or .fly format)
        
        Returns
        -------
        MeshBase
            Loaded mesh
        """
        pass
    
    @abstractmethod
    def create_rectangle_mesh(self,
                             length: float,
                             height: float,
                             nx: int,
                             ny: int) -> MeshBase:
        """
        Create a simple rectangular mesh.
        
        Parameters
        ----------
        length : float
            Length of the domain (x-direction)
        height : float
            Height of the domain (y-direction)
        nx : int
            Number of cells in x-direction
        ny : int
            Number of cells in y-direction
        
        Returns
        -------
        MeshBase
            Rectangular mesh
        """
        pass
    
    # =========================================================================
    # Field creation
    # =========================================================================
    
    @abstractmethod
    def create_scalar_field(self,
                           mesh: MeshBase,
                           value: Union[float, np.ndarray] = 0.0,
                           name: str = "",
                           on_boundary: bool = False) -> FieldBase:
        """
        Create a scalar field on the mesh.
        
        Parameters
        ----------
        mesh : MeshBase
            Mesh to create field on
        value : float or np.ndarray
            Initial value(s) for the field
        name : str
            Name for the field
        on_boundary : bool
            If True, create field on boundary (FunctionOnBoundary in escript)
        
        Returns
        -------
        FieldBase
            Scalar field
        """
        pass
    
    @abstractmethod
    def create_vector_field(self,
                           mesh: MeshBase,
                           values: Union[Tuple[float, float], np.ndarray] = (0.0, 0.0),
                           name: str = "") -> FieldBase:
        """
        Create a vector field on the mesh.
        
        Parameters
        ----------
        mesh : MeshBase
            Mesh to create field on
        values : tuple or np.ndarray
            Initial values for the vector components (vx, vy)
        name : str
            Name for the field
        
        Returns
        -------
        FieldBase
            Vector field
        """
        pass
    
    @abstractmethod
    def create_tensor_field(self,
                           mesh: MeshBase,
                           values: Union[np.ndarray, Tuple[Tuple[float, float], 
                                                           Tuple[float, float]]] = None,
                           name: str = "") -> FieldBase:
        """
        Create a 2x2 tensor field on the mesh.
        
        Parameters
        ----------
        mesh : MeshBase
            Mesh to create field on
        values : np.ndarray or nested tuple
            Initial values for the tensor [[a00, a01], [a10, a11]]
            If None, creates zero tensor
        name : str
            Name for the field
        
        Returns
        -------
        FieldBase
            Tensor field
        """
        pass
    
    # =========================================================================
    # PDE solver creation
    # =========================================================================
    
    @abstractmethod
    def create_pde_solver(self, mesh: MeshBase, symmetric: bool = True) -> PDESolverBase:
        """
        Create a PDE solver for a scalar equation.
        
        Parameters
        ----------
        mesh : MeshBase
            Mesh to solve on
        symmetric : bool
            Whether the PDE matrix is symmetric (enables faster solvers)
        
        Returns
        -------
        PDESolverBase
            PDE solver instance
        """
        pass
    
    @abstractmethod
    def create_projection_solver(self, mesh: MeshBase) -> PDESolverBase:
        """
        Create a solver for projecting element values to nodes.
        
        This is used to convert flux values (computed at element centers)
        to nodal values for output and visualization.
        
        Parameters
        ----------
        mesh : MeshBase
            Mesh to solve on
        
        Returns
        -------
        PDESolverBase
            Projection solver instance
        """
        pass
    
    # =========================================================================
    # Mathematical operations
    # =========================================================================
    
    @abstractmethod
    def gradient(self, field: FieldBase) -> FieldBase:
        """
        Compute the gradient of a scalar field.
        
        Parameters
        ----------
        field : FieldBase
            Scalar field
        
        Returns
        -------
        FieldBase
            Vector field containing gradient
        """
        pass
    
    @abstractmethod
    def divergence(self, field: FieldBase) -> FieldBase:
        """
        Compute the divergence of a vector field.
        
        Parameters
        ----------
        field : FieldBase
            Vector field
        
        Returns
        -------
        FieldBase
            Scalar field containing divergence
        """
        pass
    
    @abstractmethod
    def integrate(self, field: FieldBase, over_boundary: bool = False) -> float:
        """
        Integrate a field over the domain or boundary.
        
        Parameters
        ----------
        field : FieldBase
            Field to integrate
        over_boundary : bool
            If True, integrate over the boundary only
        
        Returns
        -------
        float
            Integral value
        """
        pass
    
    @abstractmethod
    def interpolate(self, field: FieldBase, target_space: str) -> FieldBase:
        """
        Interpolate a field to a different function space.
        
        Parameters
        ----------
        field : FieldBase
            Field to interpolate
        target_space : str
            Target space: 'nodes', 'cells', 'boundary'
        
        Returns
        -------
        FieldBase
            Interpolated field
        """
        pass
    
    # =========================================================================
    # Conditional operations (where functions)
    # =========================================================================
    
    @abstractmethod
    def where_positive(self, field: FieldBase) -> FieldBase:
        """
        Return 1 where field > 0, 0 otherwise.
        
        Equivalent to escript's wherePositive.
        
        Parameters
        ----------
        field : FieldBase
            Input field
        
        Returns
        -------
        FieldBase
            Binary mask field
        """
        pass
    
    @abstractmethod
    def where_negative(self, field: FieldBase) -> FieldBase:
        """
        Return 1 where field < 0, 0 otherwise.
        
        Equivalent to escript's whereNegative.
        
        Parameters
        ----------
        field : FieldBase
            Input field
        
        Returns
        -------
        FieldBase
            Binary mask field
        """
        pass
    
    @abstractmethod
    def where_zero(self, field: FieldBase, tolerance: float = 1e-10) -> FieldBase:
        """
        Return 1 where field == 0 (within tolerance), 0 otherwise.
        
        Equivalent to escript's whereZero.
        
        Parameters
        ----------
        field : FieldBase
            Input field
        tolerance : float
            Tolerance for zero comparison
        
        Returns
        -------
        FieldBase
            Binary mask field
        """
        pass
    
    # =========================================================================
    # Reduction operations
    # =========================================================================
    
    @abstractmethod
    def minimum(self, field: FieldBase) -> float:
        """
        Get the minimum value of a field.
        
        Equivalent to escript's inf.
        
        Parameters
        ----------
        field : FieldBase
            Input field
        
        Returns
        -------
        float
            Minimum value
        """
        pass
    
    @abstractmethod
    def maximum(self, field: FieldBase) -> float:
        """
        Get the maximum value of a field.
        
        Equivalent to escript's sup.
        
        Parameters
        ----------
        field : FieldBase
            Input field
        
        Returns
        -------
        float
            Maximum value
        """
        pass
    
    @abstractmethod
    def max_abs(self, field: FieldBase) -> float:
        """
        Get the maximum absolute value of a field.
        
        Equivalent to escript's Lsup.
        
        Parameters
        ----------
        field : FieldBase
            Input field
        
        Returns
        -------
        float
            Maximum absolute value
        """
        pass
    
    # =========================================================================
    # I/O operations
    # =========================================================================
    
    @abstractmethod
    def save_vtk(self,
                filename: str,
                mesh: MeshBase,
                fields: dict) -> None:
        """
        Save mesh and fields to VTK format for visualization.
        
        Parameters
        ----------
        filename : str
            Output filename (without extension)
        mesh : MeshBase
            Mesh to save
        fields : dict
            Dictionary of field names and FieldBase objects to save
        """
        pass
    
    # =========================================================================
    # Utility operations
    # =========================================================================
    
    @abstractmethod
    def kronecker(self, mesh: MeshBase) -> FieldBase:
        """
        Return the Kronecker delta (identity tensor) on the mesh.
        
        Parameters
        ----------
        mesh : MeshBase
            Mesh to create kronecker on
        
        Returns
        -------
        FieldBase
            Identity tensor field
        """
        pass
    
    @abstractmethod
    def get_domain_coordinates(self, mesh: MeshBase) -> FieldBase:
        """
        Get the coordinate field of the mesh domain.
        
        This is equivalent to mesh.getX() in escript.
        
        Parameters
        ----------
        mesh : MeshBase
            Mesh to get coordinates from
        
        Returns
        -------
        FieldBase
            Vector field containing x, y coordinates
        """
        pass
