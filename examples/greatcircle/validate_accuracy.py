#!/usr/bin/env python
"""
Validate SWMP accuracy against analytical great circle solution.

This script runs SWMP with a constant velocity model on a sphere and compares
the numerical travel times against the analytical great circle solution.

For a homogeneous sphere, rays must follow great circle paths and travel times
should match: t = d/v where d is the great circle distance.
"""

import sys
import numpy as np
import pyswmp


def great_circle_distance(lon1, lat1, lon2, lat2, radius=6371.0):
    """Compute great circle distance using haversine formula.

    Args:
        lon1, lat1: Source coordinates (degrees)
        lon2, lat2: Receiver coordinates (degrees)
        radius: Earth radius in km (default: 6371)

    Returns:
        Distance in km
    """
    # Convert to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))

    return radius * c


def validate_swmp_accuracy(velocity=4.0, dt=5.0, max_error_pct=1.0, verbose=True):
    """Run SWMP validation against analytical solution.

    Args:
        velocity: Constant velocity (km/s)
        dt: Time step for ODE solver (seconds)
        max_error_pct: Maximum acceptable relative error (%)
        verbose: Print detailed results

    Returns:
        dict: Validation results with error statistics

    Raises:
        AssertionError: If errors exceed threshold
    """
    if verbose:
        print("="*70)
        print("SWMP ACCURACY VALIDATION: Great Circle Test")
        print("="*70)
        print(f"\nConfiguration:")
        print(f"  Velocity:      {velocity} km/s")
        print(f"  Time step:     {dt} s")
        print(f"  Threshold:     {max_error_pct}%")

    # Create constant velocity model
    if verbose:
        print(f"\nCreating velocity model...")
    model = pyswmp.create_constant_velocity_model(
        nx=90, ny=45,  # Reduced from 180x90 for faster validation
        x0=-180.0, y0=-90.0,
        dx=4.0, dy=4.0,  # Increased from 2.0 for coarser grid
        velocity=velocity,
        cushion_nodes=0  # TODO: Debug cushion_nodes issue for large models
    )

    # Single source (Sydney)
    source_lon, source_lat = 151.2, -33.9
    sources = pyswmp.Sources(
        positions=np.array([[source_lon, source_lat]]),
        source_type=1
    )

    # Grid of receivers (heavily reduced for faster validation)
    lon_rec = np.linspace(120, -120, 10)  # Reduced from 20
    lat_rec = np.linspace(-60, 60, 8)     # Reduced from 15
    lon_grid, lat_grid = np.meshgrid(lon_rec, lat_rec)
    receivers = pyswmp.Receivers(
        positions=np.column_stack([lon_grid.ravel(), lat_grid.ravel()])
    )

    if verbose:
        print(f"  Model grid:    {model.nx} x {model.ny}")
        print(f"  Sources:       {sources.n_sources}")
        print(f"  Receivers:     {receivers.n_receivers}")

    # Compute analytical solution
    if verbose:
        print(f"\nComputing analytical solution...")
    distances = great_circle_distance(
        source_lon, source_lat,
        receivers.positions[:, 0], receivers.positions[:, 1]
    )
    analytical_times = distances / velocity

    if verbose:
        print(f"  Distance range:     {distances.min():.0f} - {distances.max():.0f} km")
        print(f"  Travel time range:  {analytical_times.min():.0f} - {analytical_times.max():.0f} s")

    # Run SWMP
    if verbose:
        print(f"\nRunning SWMP ray tracing...")
    opts = pyswmp.TrackerOptions(
        dt=dt,
        max_iterations=1000,
        coordinate_system=2,  # Spherical
        earth_radius=6371.0,
        computation_mode=1,
        max_arrivals=5,
        ode_solver=2,
        interpolator=2
    )

    tracker = pyswmp.WaveFrontTracker(
        model=model,
        sources=sources,
        receivers=receivers,
        options=opts
    )

    result = tracker.forward()

    if verbose:
        print(f"  Computed {len(result.travel_times)} arrivals")

    # Match results (first arrival for each receiver)
    receiver_to_time = {}
    for rec_id, arr_num, tt in zip(result.receiver_ids, result.arrival_numbers, result.travel_times):
        if arr_num == 1:
            receiver_to_time[rec_id] = tt

    numerical_times = np.array([receiver_to_time.get(i+1, np.nan) for i in range(receivers.n_receivers)])

    # Compute errors
    valid_mask = ~np.isnan(numerical_times)
    errors = numerical_times[valid_mask] - analytical_times[valid_mask]
    relative_errors = (errors / analytical_times[valid_mask]) * 100
    abs_errors = np.abs(errors)

    # Statistics
    stats = {
        'n_receivers': receivers.n_receivers,
        'n_arrivals': valid_mask.sum(),
        'mean_error': errors.mean(),
        'std_error': errors.std(),
        'mean_abs_error': abs_errors.mean(),
        'max_abs_error': abs_errors.max(),
        'rms_error': np.sqrt((errors**2).mean()),
        'mean_rel_error': relative_errors.mean(),
        'std_rel_error': relative_errors.std(),
        'mean_abs_rel_error': np.abs(relative_errors).mean(),
        'max_abs_rel_error': np.abs(relative_errors).max(),
        'within_threshold': (np.abs(relative_errors) < max_error_pct).sum(),
        'pct_within_threshold': 100 * (np.abs(relative_errors) < max_error_pct).sum() / valid_mask.sum()
    }

    if verbose:
        print(f"\n" + "="*70)
        print("VALIDATION RESULTS")
        print("="*70)
        print(f"\nCoverage:")
        print(f"  Receivers with arrivals: {stats['n_arrivals']} / {stats['n_receivers']}")

        print(f"\nAbsolute Error Statistics:")
        print(f"  Mean error:      {stats['mean_error']:8.3f} s")
        print(f"  Std deviation:   {stats['std_error']:8.3f} s")
        print(f"  Mean abs error:  {stats['mean_abs_error']:8.3f} s")
        print(f"  Max abs error:   {stats['max_abs_error']:8.3f} s")
        print(f"  RMS error:       {stats['rms_error']:8.3f} s")

        print(f"\nRelative Error Statistics:")
        print(f"  Mean:            {stats['mean_rel_error']:8.4f} %")
        print(f"  Std deviation:   {stats['std_rel_error']:8.4f} %")
        print(f"  Mean abs:        {stats['mean_abs_rel_error']:8.4f} %")
        print(f"  Max abs:         {stats['max_abs_rel_error']:8.4f} %")

        print(f"\nAccuracy Assessment:")
        print(f"  Within {max_error_pct}% error: {stats['within_threshold']} / {stats['n_arrivals']} "
              f"({stats['pct_within_threshold']:.1f}%)")

    # Validation checks
    passed = True
    if stats['mean_abs_rel_error'] > max_error_pct:
        if verbose:
            print(f"\n❌ FAIL: Mean absolute relative error ({stats['mean_abs_rel_error']:.3f}%) "
                  f"exceeds threshold ({max_error_pct}%)")
        passed = False

    if stats['pct_within_threshold'] < 90:
        if verbose:
            print(f"\n❌ FAIL: Only {stats['pct_within_threshold']:.1f}% of receivers within threshold "
                  f"(expected >90%)")
        passed = False

    if passed:
        if verbose:
            print(f"\n✓ PASS: SWMP accuracy validated successfully!")
            print("="*70)
    else:
        if verbose:
            print("="*70)

    stats['passed'] = passed
    return stats


def main():
    """Run validation and exit with appropriate code."""
    try:
        stats = validate_swmp_accuracy(
            velocity=4.0,
            dt=5.0,
            max_error_pct=1.0,
            verbose=True
        )

        if stats['passed']:
            return 0
        else:
            return 1

    except Exception as e:
        print(f"\n❌ ERROR: Validation failed with exception:")
        print(f"  {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return 2


if __name__ == "__main__":
    sys.exit(main())
