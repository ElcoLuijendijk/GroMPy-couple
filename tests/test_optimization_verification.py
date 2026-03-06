"""
Tests to verify benchmark ModelParameters are self-consistent and well-formed.

These tests are independent of specific values in model_parameters_sw_benchmark.py
so that the file can be adjusted freely without breaking the test suite.
"""

import pytest

from model_input.model_parameters_sw_benchmark import ModelParameters


class TestOptimizationParameters:
    """Verify that the benchmark parameters are self-consistent and well-formed.

    These tests do NOT assert specific values from model_parameters_sw_benchmark.py
    so that the file can be adjusted freely without breaking the test suite.
    """

    def test_dt_inc_positive(self):
        """Timestep growth factor must be >= 1 (no shrinking)."""
        params = ModelParameters()
        assert params.dt_inc >= 1.0, (
            f"dt_inc = {params.dt_inc} should be >= 1.0"
        )

    def test_pressure_convergence_criterion_positive(self):
        """Pressure convergence criterion must be a positive number."""
        params = ModelParameters()
        assert params.pressure_convergence_criterion > 0, (
            f"pressure_convergence_criterion = {params.pressure_convergence_criterion} "
            "must be positive"
        )

    def test_concentration_convergence_criterion_positive(self):
        """Concentration convergence criterion must be a positive number."""
        params = ModelParameters()
        assert params.concentration_convergence_criterion > 0, (
            f"concentration_convergence_criterion = "
            f"{params.concentration_convergence_criterion} "
            "must be positive"
        )

    def test_max_iterations_positive(self):
        """Max iterations must be a positive integer."""
        params = ModelParameters()
        assert params.max_iterations > 0, (
            f"max_iterations = {params.max_iterations} must be positive"
        )
