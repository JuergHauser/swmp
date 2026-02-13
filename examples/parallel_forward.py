#!/usr/bin/env python
"""Example: Parallel wavefront tracking with the file-free API.

Demonstrates two execution modes:
  1. Sequential (default)
  2. Parallel via user-provided ProcessPoolExecutor

Usage:
    python parallel_forward.py
"""

import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import pyswmp


def make_test_data():
    """Create a simple test dataset."""
    model = pyswmp.create_constant_velocity_model(
        nx=100, ny=80,
        x0=110.0, y0=-45.0,
        dx=0.5, dy=0.5,
        velocity=3.5,
    )

    sources = pyswmp.Sources(
        positions=np.array([
            [135.0, -25.0],
            [140.0, -20.0],
            [130.0, -30.0],
            [125.0, -35.0],
        ]),
        source_type=1,
    )

    receivers = pyswmp.Receivers(
        positions=np.array([
            [120.0, -30.0],
            [125.0, -25.0],
            [130.0, -35.0],
            [135.0, -40.0],
            [140.0, -25.0],
        ]),
    )

    return model, sources, receivers


def main():
    model, sources, receivers = make_test_data()
    opts = pyswmp.TrackerOptions(coordinate_system=2)
    tracker = pyswmp.WaveFrontTracker(model, opts)

    # --- Sequential (baseline) ---
    t0 = time.time()
    result_seq = tracker.forward(sources, receivers)
    t_seq = time.time() - t0
    print(f"Sequential:       {len(result_seq.travel_times)} arrivals in {t_seq:.3f}s")

    # --- Parallel via user-provided pool ---
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=2) as pool:
        result_pool = tracker.forward(sources, receivers, pool=pool)
    t_pool = time.time() - t0
    print(f"Pool (2 workers): {len(result_pool.travel_times)} arrivals in {t_pool:.3f}s")

    # --- Summary ---
    if t_seq > 0:
        print(f"\nSpeedup: {t_seq / t_pool:.2f}x")


if __name__ == "__main__":
    main()
