"""Unit tests for pyswmp.options module."""

import pytest
from pyswmp import TrackerOptions


class TestTrackerOptions:
    """Tests for TrackerOptions dataclass."""

    def test_default_options(self):
        """Test creating options with default values."""
        opts = TrackerOptions()

        # Check some default values
        assert opts.dt > 0
        assert opts.max_iterations > 0
        assert opts.coordinate_system in [1, 2]
        assert opts.max_arrivals >= 1

    def test_custom_options(self):
        """Test creating options with custom values."""
        opts = TrackerOptions(
            dt=0.1,
            max_iterations=5000,
            coordinate_system=2,
            max_arrivals=20
        )

        assert opts.dt == 0.1
        assert opts.max_iterations == 5000
        assert opts.coordinate_system == 2
        assert opts.max_arrivals == 20

    def test_coordinate_system_validation(self):
        """Test that coordinate_system is validated."""
        # Valid values should work
        opts1 = TrackerOptions(coordinate_system=1)
        opts2 = TrackerOptions(coordinate_system=2)

        assert opts1.coordinate_system == 1
        assert opts2.coordinate_system == 2

    def test_boolean_flags(self):
        """Test boolean configuration flags."""
        opts = TrackerOptions(
            extract_raypaths=True,
            extract_frechet=True,
            source_specific_receivers=False
        )

        assert opts.extract_raypaths is True
        assert opts.extract_frechet is True
        assert opts.source_specific_receivers is False

    def test_to_dict(self):
        """Test converting options to dictionary."""
        opts = TrackerOptions(
            dt=0.1,
            max_arrivals=15,
            coordinate_system=2
        )

        opts_dict = opts.to_dict()

        assert isinstance(opts_dict, dict)
        assert opts_dict['dt'] == 0.1
        assert opts_dict['max_arrivals'] == 15
        assert opts_dict['coordinate_system'] == 2

    def test_from_dict(self):
        """Test creating options from dictionary."""
        opts_dict = {
            'dt': 0.2,
            'max_arrivals': 25,
            'coordinate_system': 1
        }

        opts = TrackerOptions(**opts_dict)

        assert opts.dt == 0.2
        assert opts.max_arrivals == 25
        assert opts.coordinate_system == 1

    def test_ode_solver_options(self):
        """Test ODE solver configuration."""
        # Test different solver options (1=RK4, 2=RK5, 3=adaptive RK5)
        for solver in [1, 2, 3]:
            opts = TrackerOptions(ode_solver=solver)
            assert opts.ode_solver == solver

    def test_interpolator_options(self):
        """Test interpolator configuration."""
        # Test different interpolators (1=linear, 2=weighted average)
        for interp in [1, 2]:
            opts = TrackerOptions(interpolator=interp)
            assert opts.interpolator == interp

    def test_computation_mode(self):
        """Test computation mode (1=kinematic, 2=kinematic+spreading)."""
        opts1 = TrackerOptions(computation_mode=1)
        opts2 = TrackerOptions(computation_mode=2)

        assert opts1.computation_mode == 1
        assert opts2.computation_mode == 2
