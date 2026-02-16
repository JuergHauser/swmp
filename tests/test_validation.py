"""Integration tests for SWMP accuracy validation.

These tests validate SWMP's numerical accuracy against analytical solutions.
"""

import numpy as np
import pytest

import pyswmp


def great_circle_distance(lon1, lat1, lon2, lat2, radius=6371.0):
    """Compute great circle distance using haversine formula.

    Args:
        lon1, lat1: Source coordinates (degrees)
        lon2, lat2: Receiver coordinates (degrees)
        radius: Earth radius in km

    Returns:
        Distance in km
    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))

    return radius * c


@pytest.mark.slow
@pytest.mark.integration
def test_great_circle_accuracy_simple():
    """Test SWMP accuracy against great circle solution with simple setup.

    This is a simplified test with fewer receivers for faster execution.
    """
    # Create constant velocity model
    model = pyswmp.create_constant_velocity_model(
        nx=90, ny=45,  # Coarser grid for faster test
        x0=-180.0, y0=-90.0,
        dx=4.0, dy=4.0,  # 4 degree resolution
        velocity=4.0,
        cushion_nodes=3
    )

    # Single source
    source_lon, source_lat = 0.0, 0.0  # Equator at prime meridian
    sources = pyswmp.Sources(
        positions=np.array([[source_lon, source_lat]]),
        source_type=1
    )

    # Small grid of receivers
    lon_rec = np.linspace(-60, 60, 10)
    lat_rec = np.linspace(-30, 30, 8)
    lon_grid, lat_grid = np.meshgrid(lon_rec, lat_rec)
    receivers = pyswmp.Receivers(
        positions=np.column_stack([lon_grid.ravel(), lat_grid.ravel()])
    )

    # Compute analytical solution
    distances = great_circle_distance(
        source_lon, source_lat,
        receivers.positions[:, 0], receivers.positions[:, 1]
    )
    analytical_times = distances / 4.0

    # Run SWMP
    opts = pyswmp.TrackerOptions(
        dt=5.0,
        max_iterations=1000,
        coordinate_system=2,  # Spherical
        earth_radius=6371.0,
        computation_mode=1,
        max_arrivals=5,
        ode_solver=2
    )

    tracker = pyswmp.WaveFrontTracker(
        model=model,
        sources=sources,
        receivers=receivers,
        options=opts
    )

    result = tracker.forward()

    # Match results
    receiver_to_time = {}
    for rec_id, arr_num, tt in zip(result.receiver_ids, result.arrival_numbers, result.travel_times):
        if arr_num == 1:
            receiver_to_time[rec_id] = tt

    numerical_times = np.array([receiver_to_time.get(i+1, np.nan) for i in range(receivers.n_receivers)])

    # Compute errors
    valid_mask = ~np.isnan(numerical_times)
    errors = numerical_times[valid_mask] - analytical_times[valid_mask]
    relative_errors = (errors / analytical_times[valid_mask]) * 100

    # Assertions
    assert valid_mask.sum() > 0, "No arrivals computed"
    assert valid_mask.sum() / receivers.n_receivers > 0.8, "Too few arrivals (<80%)"

    # Check accuracy
    mean_abs_rel_error = np.abs(relative_errors).mean()
    max_abs_rel_error = np.abs(relative_errors).max()

    assert mean_abs_rel_error < 2.0, \
        f"Mean absolute relative error ({mean_abs_rel_error:.3f}%) exceeds 2%"
    assert max_abs_rel_error < 5.0, \
        f"Max absolute relative error ({max_abs_rel_error:.3f}%) exceeds 5%"

    # Check that most receivers are within 1% error
    within_1pct = (np.abs(relative_errors) < 1.0).sum()
    pct_within_1pct = 100 * within_1pct / valid_mask.sum()

    assert pct_within_1pct > 70, \
        f"Only {pct_within_1pct:.1f}% of receivers within 1% error (expected >70%)"


@pytest.mark.slow
@pytest.mark.integration
@pytest.mark.parametrize("velocity", [3.0, 4.0, 5.0])
def test_great_circle_different_velocities(velocity):
    """Test SWMP accuracy with different velocities.

    Args:
        velocity: Constant velocity to test (km/s)
    """
    # Small model for fast testing
    model = pyswmp.create_constant_velocity_model(
        nx=45, ny=23,
        x0=-90.0, y0=-45.0,
        dx=4.0, dy=4.0,
        velocity=velocity,
        cushion_nodes=3
    )

    # Source and receivers
    source_lon, source_lat = 0.0, 0.0
    sources = pyswmp.Sources(positions=np.array([[source_lon, source_lat]]), source_type=1)

    # Sparse receiver grid
    lon_rec = np.linspace(-60, 60, 6)
    lat_rec = np.linspace(-30, 30, 5)
    lon_grid, lat_grid = np.meshgrid(lon_rec, lat_rec)
    receivers = pyswmp.Receivers(
        positions=np.column_stack([lon_grid.ravel(), lat_grid.ravel()])
    )

    # Analytical solution
    distances = great_circle_distance(
        source_lon, source_lat,
        receivers.positions[:, 0], receivers.positions[:, 1]
    )
    analytical_times = distances / velocity

    # Run SWMP
    opts = pyswmp.TrackerOptions(
        dt=5.0,
        max_iterations=800,
        coordinate_system=2,
        earth_radius=6371.0,
        computation_mode=1,
        max_arrivals=3
    )

    tracker = pyswmp.WaveFrontTracker(
        model=model, sources=sources, receivers=receivers, options=opts
    )
    result = tracker.forward()

    # Match and compute errors
    receiver_to_time = {}
    for rec_id, arr_num, tt in zip(result.receiver_ids, result.arrival_numbers, result.travel_times):
        if arr_num == 1:
            receiver_to_time[rec_id] = tt

    numerical_times = np.array([receiver_to_time.get(i+1, np.nan) for i in range(receivers.n_receivers)])
    valid_mask = ~np.isnan(numerical_times)

    assert valid_mask.sum() > 0, f"No arrivals for velocity={velocity} km/s"

    errors = numerical_times[valid_mask] - analytical_times[valid_mask]
    relative_errors = (errors / analytical_times[valid_mask]) * 100

    # Accuracy should be velocity-independent (percentage error)
    mean_abs_rel_error = np.abs(relative_errors).mean()
    assert mean_abs_rel_error < 2.5, \
        f"Velocity {velocity} km/s: Mean error {mean_abs_rel_error:.3f}% exceeds 2.5%"


def test_haversine_formula():
    """Test the haversine formula implementation."""
    # Test known distances
    # Equator at 0° to equator at 90° longitude = 1/4 Earth circumference
    # Circumference = 2π × 6371 km ≈ 40030 km
    # Quarter circle = 10007.5 km
    d = great_circle_distance(0, 0, 90, 0, radius=6371.0)
    expected = np.pi / 2 * 6371.0
    assert np.abs(d - expected) < 1.0, "Haversine formula incorrect for equator"

    # Antipodal points (opposite sides of Earth)
    d = great_circle_distance(0, 0, 180, 0)
    expected = np.pi * 6371.0
    assert np.abs(d - expected) < 1.0, "Haversine formula incorrect for antipodal points"

    # Same point
    d = great_circle_distance(120, -30, 120, -30)
    assert d < 0.01, "Haversine formula should return ~0 for same point"


def test_constant_velocity_model_uniformity():
    """Test that constant velocity model is actually constant."""
    model = pyswmp.create_constant_velocity_model(
        nx=50, ny=40,
        x0=0.0, y0=0.0,
        dx=1.0, dy=1.0,
        velocity=3.5
    )

    # All velocities should be exactly 3.5
    assert np.all(model.velocities == 3.5), "Velocity model not constant"
    assert model.velocities.min() == 3.5
    assert model.velocities.max() == 3.5
    assert model.velocities.std() == 0.0


@pytest.mark.skip(reason="Requires long computation time - run manually for full validation")
def test_great_circle_accuracy_full():
    """Full validation test with dense receiver grid.

    This is the full validation from the examples/greatcircle/ notebook.
    Run manually when needed for comprehensive validation.
    """
    # Full resolution model
    model = pyswmp.create_constant_velocity_model(
        nx=180, ny=90,
        x0=-180.0, y0=-90.0,
        dx=2.0, dy=2.0,
        velocity=4.0,
        cushion_nodes=3
    )

    # Source
    source_lon, source_lat = 151.2, -33.9
    sources = pyswmp.Sources(positions=np.array([[source_lon, source_lat]]), source_type=1)

    # Dense receiver grid
    lon_rec = np.linspace(120, -120, 30)
    lat_rec = np.linspace(-60, 60, 20)
    lon_grid, lat_grid = np.meshgrid(lon_rec, lat_rec)
    receivers = pyswmp.Receivers(
        positions=np.column_stack([lon_grid.ravel(), lat_grid.ravel()])
    )

    # Analytical solution
    distances = great_circle_distance(
        source_lon, source_lat,
        receivers.positions[:, 0], receivers.positions[:, 1]
    )
    analytical_times = distances / 4.0

    # Run SWMP
    opts = pyswmp.TrackerOptions(
        dt=5.0,
        max_iterations=1000,
        coordinate_system=2,
        earth_radius=6371.0,
        computation_mode=1,
        max_arrivals=5,
        ode_solver=2
    )

    tracker = pyswmp.WaveFrontTracker(
        model=model, sources=sources, receivers=receivers, options=opts
    )
    result = tracker.forward()

    # Match and compute errors
    receiver_to_time = {}
    for rec_id, arr_num, tt in zip(result.receiver_ids, result.arrival_numbers, result.travel_times):
        if arr_num == 1:
            receiver_to_time[rec_id] = tt

    numerical_times = np.array([receiver_to_time.get(i+1, np.nan) for i in range(receivers.n_receivers)])
    valid_mask = ~np.isnan(numerical_times)

    errors = numerical_times[valid_mask] - analytical_times[valid_mask]
    relative_errors = (errors / analytical_times[valid_mask]) * 100

    # Strict accuracy requirements for full validation
    mean_abs_rel_error = np.abs(relative_errors).mean()
    within_1pct = (np.abs(relative_errors) < 1.0).sum() / valid_mask.sum() * 100

    assert mean_abs_rel_error < 1.0, \
        f"Full validation failed: mean error {mean_abs_rel_error:.3f}% exceeds 1%"
    assert within_1pct > 90, \
        f"Full validation failed: only {within_1pct:.1f}% within 1% (expected >90%)"

    print(f"\n✓ Full validation PASSED:")
    print(f"  Mean abs rel error: {mean_abs_rel_error:.4f}%")
    print(f"  Within 1% error: {within_1pct:.1f}%")
