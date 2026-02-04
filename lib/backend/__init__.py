"""
Backend abstraction layer for GroMPy-couple.

This module provides a unified interface for different numerical backends
(esys-escript and FiPy) to solve coupled groundwater flow and solute
transport equations.

Usage
-----
>>> from lib.backend import get_backend
>>> backend = get_backend('escript')  # or 'fipy'
>>> mesh = backend.create_mesh(...)
>>> pde = backend.create_pressure_pde(mesh)

Available Backends
------------------
- 'escript' : esys-escript (default) - Finite Element Method
- 'fipy' : FiPy - Finite Volume Method

"""

__all__ = ['get_backend', 'AVAILABLE_BACKENDS']

AVAILABLE_BACKENDS = ['escript', 'fipy']


def get_backend(name='escript'):
    """
    Factory function to get the appropriate numerical backend.
    
    Parameters
    ----------
    name : str, optional
        Backend name: 'escript' (default) or 'fipy'
    
    Returns
    -------
    BackendBase
        Instance of the requested backend
    
    Raises
    ------
    ValueError
        If the backend name is not recognized
    ImportError
        If the backend's dependencies are not installed
    
    Examples
    --------
    >>> backend = get_backend('escript')
    >>> backend = get_backend('fipy')
    """
    name = name.lower().strip()
    
    if name == 'escript':
        try:
            from lib.backend.escript_backend import EscriptBackend
            return EscriptBackend()
        except ImportError as e:
            raise ImportError(
                f"Cannot import escript backend. "
                f"Make sure esys-escript is installed. Original error: {e}"
            )
    
    elif name == 'fipy':
        try:
            from lib.backend.fipy_backend import FiPyBackend
            return FiPyBackend()
        except ImportError as e:
            raise ImportError(
                f"Cannot import FiPy backend. "
                f"Install with: pip install fipy meshio. Original error: {e}"
            )
    
    else:
        raise ValueError(
            f"Unknown backend: '{name}'. "
            f"Available backends: {AVAILABLE_BACKENDS}"
        )


def check_backend_available(name):
    """
    Check if a backend is available without instantiating it.
    
    Parameters
    ----------
    name : str
        Backend name to check
    
    Returns
    -------
    bool
        True if the backend can be imported, False otherwise
    """
    name = name.lower().strip()
    
    if name == 'escript':
        try:
            import esys.escript
            return True
        except ImportError:
            return False
    
    elif name == 'fipy':
        try:
            import fipy
            return True
        except ImportError:
            return False
    
    return False
