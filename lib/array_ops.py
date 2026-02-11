"""
Backend-agnostic array operations for escript Data objects and numpy arrays.

This module provides unified interface for array operations across different
numerical backends (escript, FiPy, numpy) to enable portable code that works
with any backend.

Functions:
    min_value: Get minimum value from array
    max_value: Get maximum value from array
    max_absolute_value: Get L-infinity norm (max absolute value)
    convert_to_numpy: Safe conversion of escript Data to numpy
    extract_flux_component: Extract scalar or vector component from flux values
"""

import numpy as np
import warnings


def convert_to_numpy(array):
    """
    Convert escript Data object to numpy array, or return numpy array as-is.
    
    Handles both escript Data objects and numpy arrays transparently.
    
    Parameters
    ----------
    array : numpy.ndarray or escript.Data
        Input array from either backend
    
    Returns
    -------
    numpy.ndarray
        Converted or original array as numpy array
    """
    # Check if it's an escript Data object
    if hasattr(array, 'asNumpy'):
        # escript Data object with asNumpy() method
        try:
            return array.asNumpy()
        except Exception as e:
            warnings.warn(f"Failed to convert escript Data with asNumpy(): {e}. Attempting numpy conversion.")
            try:
                return np.asarray(array)
            except Exception:
                raise ValueError(f"Could not convert array to numpy: {e}")
    
    # Already a numpy array or compatible object
    try:
        return np.asarray(array)
    except Exception as e:
        raise ValueError(f"Could not convert array to numpy: {e}")


def min_value(array):
    """
    Get minimum value from array (works with escript Data or numpy arrays).
    
    Silently falls back to numpy operations if escript operations unavailable,
    with a warning if conversion was necessary.
    
    Parameters
    ----------
    array : numpy.ndarray or escript.Data
        Input array
    
    Returns
    -------
    float
        Minimum value in array
    """
    # Try escript-specific method first (faster if available)
    if hasattr(array, 'inf'):
        try:
            result = array.inf()
            if hasattr(result, 'toFloat'):
                return float(result.toFloat())
            else:
                return float(result)
        except Exception as e:
            warnings.warn(f"escript.inf() failed, falling back to numpy.min(): {e}")
    
    # Fall back to numpy
    try:
        arr_np = convert_to_numpy(array)
        return float(np.min(arr_np))
    except Exception as e:
        raise ValueError(f"Could not compute minimum value: {e}")


def max_value(array):
    """
    Get maximum value from array (works with escript Data or numpy arrays).
    
    Silently falls back to numpy operations if escript operations unavailable,
    with a warning if conversion was necessary.
    
    Parameters
    ----------
    array : numpy.ndarray or escript.Data
        Input array
    
    Returns
    -------
    float
        Maximum value in array
    """
    # Try escript-specific method first (faster if available)
    if hasattr(array, 'sup'):
        try:
            result = array.sup()
            if hasattr(result, 'toFloat'):
                return float(result.toFloat())
            else:
                return float(result)
        except Exception as e:
            warnings.warn(f"escript.sup() failed, falling back to numpy.max(): {e}")
    
    # Fall back to numpy
    try:
        arr_np = convert_to_numpy(array)
        return float(np.max(arr_np))
    except Exception as e:
        raise ValueError(f"Could not compute maximum value: {e}")


def max_absolute_value(array):
    """
    Get L-infinity norm (maximum absolute value) from array.
    
    Works with escript Data or numpy arrays. This is equivalent to
    escript.Lsup() for computing the L-infinity norm.
    
    Silently falls back to numpy operations if escript operations unavailable,
    with a warning if conversion was necessary.
    
    Parameters
    ----------
    array : numpy.ndarray or escript.Data
        Input array
    
    Returns
    -------
    float
        Maximum absolute value in array (L-infinity norm)
    """
    # Try escript-specific method first (faster if available)
    if hasattr(array, 'Lsup'):
        try:
            result = array.Lsup()
            if hasattr(result, 'toFloat'):
                return float(result.toFloat())
            else:
                return float(result)
        except Exception as e:
            warnings.warn(f"escript.Lsup() failed, falling back to numpy.max(abs()): {e}")
    
    # Fall back to numpy
    try:
        arr_np = convert_to_numpy(array)
        return float(np.max(np.abs(arr_np)))
    except Exception as e:
        raise ValueError(f"Could not compute L-infinity norm: {e}")


def extract_flux_component(flux_value, component=1):
    """
    Safely extract flux component from escript tuple or numpy scalar/array.
    
    Handles both escript integration results (tuples with [0]=horizontal, [1]=vertical)
    and FiPy numpy results (scalars or arrays). This enables backend-agnostic code
    to work with flux data regardless of the numerical backend.
    
    Parameters
    ----------
    flux_value : scalar, tuple, array, or escript Data
        Flux value from either backend:
        - Escript: tuple (h_component, v_component)
        - FiPy scalar: single float value
        - FiPy array: numpy array with [0] or [1] indices
    component : int, default=1
        Component to extract (0=horizontal, 1=vertical)
        For scalar inputs, this parameter is ignored and the scalar is returned directly.
    
    Returns
    -------
    float
        Scalar value of the requested component
    
    Notes
    -----
    - component=1 extracts vertical flux (perpendicular to land surface or seabed)
    - component=0 extracts horizontal flux (rarely used, typically 0 for 1D integration)
    - For scalar inputs: returns float(flux_value) regardless of component argument
    - For indexable inputs (tuples, arrays): tries to access flux_value[component]
    - Debug warning logged when falling back from indexed to scalar handling
    
    Raises
    ------
    ValueError
        If component is not 0 or 1
    
    Examples
    --------
    >>> # Escript backend returns tuple (horizontal, vertical)
    >>> flux_escript = (0.0, 1.5)
    >>> extract_flux_component(flux_escript, component=1)
    1.5
    
    >>> # FiPy backend returns scalar (just the vertical component)
    >>> flux_fipy = 1.5
    >>> extract_flux_component(flux_fipy, component=1)
    1.5
    
    >>> # FiPy backend returns numpy array [horizontal, vertical]
    >>> flux_fipy_arr = np.array([0.0, 1.5])
    >>> extract_flux_component(flux_fipy_arr, component=1)
    1.5
    """
    if component not in (0, 1):
        raise ValueError(f"component must be 0 or 1, got {component}")
    
    # Try to treat as indexable first (tuples, arrays, etc.)
    try:
        # Check if it's indexable and has the requested component
        if hasattr(flux_value, '__getitem__') and hasattr(flux_value, '__len__'):
            length = len(flux_value)
            if length > component:
                result = flux_value[component]
                # Handle escript Data objects that might be returned from indexing
                if hasattr(result, 'toFloat'):
                    return float(result.toFloat())
                else:
                    return float(result)
    except (TypeError, IndexError):
        pass
    
    # Fall back to treating as scalar
    try:
        warnings.warn(
            f"Flux value is scalar or indexing failed for component {component}. "
            f"Returning value directly without component extraction.",
            stacklevel=2
        )
        return float(flux_value)
    except (TypeError, ValueError) as e:
        warnings.warn(
            f"Could not extract flux component {component}: {e}. Returning 0.0",
            stacklevel=2
        )
        return 0.0

