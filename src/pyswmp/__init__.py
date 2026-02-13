"""Python wrapper for SWMP: Surface Wave Multipathing ray tracing.

This package provides a modern Python interface to the SWMP Fortran library
with support for parallel execution and in-memory data handling.

Main classes:
    WaveFrontTracker: High-level solver with pool-based parallelism
    WaveFrontTrackerResult: Result data structure with aggregation support
    LibSWMP: Low-level ctypes wrapper (advanced users)

Legacy classes (from _pyswmp):
    VelocityModelGenerator: Generate synthetic velocity models
    ObservationGenerator: Generate synthetic observations
    Visualisation: Plotting utilities
"""

from .__version__ import __version__

# New parallel-capable implementation
from .solvers import WaveFrontTracker
from .results import WaveFrontTrackerResult
from ._libswmp import LibSWMP, SWMPError
from .options import TrackerOptions
from .data import (
    VelocityModel2D,
    Sources,
    Receivers,
    create_constant_velocity_model,
    create_gradient_velocity_model,
    write_velocity_model_file,
    write_sources_file,
    write_receivers_file,
    write_config_file,
)

# Legacy classes (for backward compatibility)
from ._pyswmp import (
    VelocityModel,
    TravelTimeData,
    VelocityModelGenerator,
    ObservationGenerator,
    Visualisation,
)

__all__ = [
    # Version
    "__version__",
    # Solver
    "WaveFrontTracker",
    "WaveFrontTrackerResult",
    "LibSWMP",
    "SWMPError",
    "TrackerOptions",
    # Data structures
    "VelocityModel2D",
    "Sources",
    "Receivers",
    "create_constant_velocity_model",
    "create_gradient_velocity_model",
    # Legacy classes (still useful)
    "VelocityModelGenerator",
    "ObservationGenerator",
    "Visualisation",
]
