"""Pytest configuration and fixtures for pyswmp tests."""

import numpy as np
import pytest


@pytest.fixture
def simple_velocity_model():
    """Create a simple constant velocity model for testing."""
    from pyswmp import create_constant_velocity_model

    return create_constant_velocity_model(
        nx=50, ny=40,
        x0=110.0, y0=-45.0,
        dx=0.5, dy=0.5,
        velocity=3.5
    )


@pytest.fixture
def gradient_velocity_model():
    """Create a gradient velocity model for testing."""
    from pyswmp import create_gradient_velocity_model

    return create_gradient_velocity_model(
        nx=50, ny=40,
        x0=110.0, y0=-45.0,
        dx=0.5, dy=0.5,
        v0=3.0,
        gradient_x=0.01
    )


@pytest.fixture
def simple_sources():
    """Create simple point sources for testing."""
    from pyswmp import Sources

    return Sources(
        positions=np.array([[135.0, -25.0], [140.0, -20.0]]),
        source_type=1
    )


@pytest.fixture
def simple_receivers():
    """Create simple receivers for testing."""
    from pyswmp import Receivers

    # Create a small grid of receivers
    x = np.linspace(120, 150, 10)
    y = np.linspace(-40, -20, 8)
    xx, yy = np.meshgrid(x, y)
    positions = np.column_stack([xx.ravel(), yy.ravel()])

    return Receivers(positions)


@pytest.fixture
def default_options():
    """Create default tracker options for testing."""
    from pyswmp import TrackerOptions

    return TrackerOptions(
        coordinate_system=2,  # Spherical
        max_arrivals=10
    )


@pytest.fixture
def sample_result():
    """Create a sample WaveFrontTrackerResult for testing."""
    from pyswmp import WaveFrontTrackerResult

    return WaveFrontTrackerResult(
        source_ids=np.array([1, 1, 2, 2]),
        receiver_ids=np.array([1, 2, 1, 2]),
        arrival_numbers=np.array([1, 1, 1, 1]),
        travel_times=np.array([10.5, 15.2, 12.3, 18.1]),
        azimuths=np.array([45.0, 90.0, 135.0, 180.0]),
        spreading_factors=np.array([1.0, 1.5, 1.2, 1.8])
    )
