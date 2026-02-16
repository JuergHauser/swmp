"""Data structures for SWMP models, sources, and receivers.

This module provides Pythonic data classes for configuring wavefront tracking
simulations programmatically without requiring input files.
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class VelocityModel2D:
    """2D velocity model on a regular grid.

    Attributes:
        velocities: 2D array of velocities (nx, ny)
        x0: X-coordinate of lower-left corner
        y0: Y-coordinate of lower-left corner
        dx: Grid spacing in X direction
        dy: Grid spacing in Y direction
        cushion_nodes: Number of cushion nodes for boundary (default: 3)

    Example:
        >>> # Create a simple constant velocity model
        >>> nx, ny = 100, 80
        >>> velocities = np.full((nx, ny), 3.5)  # 3.5 km/s everywhere
        >>> model = VelocityModel2D(
        ...     velocities=velocities,
        ...     x0=110.0, y0=-45.0,  # degrees
        ...     dx=0.5, dy=0.5
        ... )

        >>> # Or create from a function
        >>> x = np.linspace(110, 160, 100)
        >>> y = np.linspace(-45, -10, 80)
        >>> xx, yy = np.meshgrid(x, y, indexing='ij')
        >>> velocities = 3.0 + 0.01 * np.sqrt(xx**2 + yy**2)
        >>> model = VelocityModel2D(velocities, 110.0, -45.0, 0.5, 0.5)
    """

    velocities: np.ndarray
    x0: float
    y0: float
    dx: float
    dy: float
    cushion_nodes: int = 3

    def __post_init__(self):
        """Validate model parameters."""
        # Ensure velocities is numpy array
        if not isinstance(self.velocities, np.ndarray):
            self.velocities = np.array(self.velocities, dtype=np.float32)

        # Ensure correct dtype
        if self.velocities.dtype != np.float32:
            self.velocities = self.velocities.astype(np.float32)

        # Validate shape
        if self.velocities.ndim != 2:
            raise ValueError(f"velocities must be 2D array, got shape {self.velocities.shape}")

        # Validate positive spacings
        if self.dx <= 0 or self.dy <= 0:
            raise ValueError("Grid spacings dx and dy must be positive")

        if self.cushion_nodes < 0:
            raise ValueError("cushion_nodes must be non-negative")

    @property
    def nx(self) -> int:
        """Number of grid points in X direction."""
        return self.velocities.shape[0]

    @property
    def ny(self) -> int:
        """Number of grid points in Y direction."""
        return self.velocities.shape[1]

    @property
    def x1(self) -> float:
        """X-coordinate of upper-right corner."""
        return self.x0 + self.dx * self.nx

    @property
    def y1(self) -> float:
        """Y-coordinate of upper-right corner."""
        return self.y0 + self.dy * self.ny

    @property
    def extent(self):
        """Model extent as (x0, x1, y0, y1)."""
        return (self.x0, self.x1, self.y0, self.y1)

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"VelocityModel2D(\n"
            f"  shape=({self.nx}, {self.ny}),\n"
            f"  extent=({self.x0:.2f}, {self.x1:.2f}, {self.y0:.2f}, {self.y1:.2f}),\n"
            f"  spacing=({self.dx:.4f}, {self.dy:.4f}),\n"
            f"  velocity_range=({self.velocities.min():.3f}, {self.velocities.max():.3f})\n"
            f")"
        )


@dataclass
class Sources:
    """Source configuration for wavefront tracking.

    Attributes:
        positions: Array of source positions, shape (n_sources, 2) with columns [x, y]
        source_type: Source type (1=point source, 2=plane wave)
        angles: Optional array of angles for plane wave sources (degrees)

    Example:
        >>> # Single point source
        >>> sources = Sources(
        ...     positions=np.array([[135.0, -25.0]]),
        ...     source_type=1
        ... )

        >>> # Multiple point sources
        >>> positions = np.array([
        ...     [130.0, -30.0],
        ...     [135.0, -25.0],
        ...     [140.0, -20.0]
        ... ])
        >>> sources = Sources(positions, source_type=1)

        >>> # Plane wave sources with angles
        >>> sources = Sources(
        ...     positions=np.array([[135.0, -25.0]]),
        ...     source_type=2,
        ...     angles=np.array([45.0])  # degrees
        ... )
    """

    positions: np.ndarray
    source_type: int = 1
    angles: Optional[np.ndarray] = None

    def __post_init__(self):
        """Validate source parameters."""
        # Ensure positions is numpy array
        if not isinstance(self.positions, np.ndarray):
            self.positions = np.array(self.positions, dtype=np.float64)

        if self.positions.dtype != np.float64:
            self.positions = self.positions.astype(np.float64)

        # Validate shape
        if self.positions.ndim == 1:
            # Single source provided as 1D array [x, y]
            if len(self.positions) != 2:
                raise ValueError("Single source must be [x, y]")
            self.positions = self.positions.reshape(1, 2)
        elif self.positions.ndim == 2:
            if self.positions.shape[1] != 2:
                raise ValueError(f"positions must have shape (n, 2), got {self.positions.shape}")
        else:
            raise ValueError(f"positions must be 1D or 2D array, got shape {self.positions.shape}")

        # Validate source type
        if self.source_type not in [1, 2]:
            raise ValueError(f"source_type must be 1 (point) or 2 (plane wave), got {self.source_type}")

        # Validate angles if provided
        if self.angles is not None:
            if not isinstance(self.angles, np.ndarray):
                self.angles = np.array(self.angles, dtype=np.float64)

            if self.angles.dtype != np.float64:
                self.angles = self.angles.astype(np.float64)

            if len(self.angles) != self.n_sources:
                raise ValueError(
                    f"angles length ({len(self.angles)}) must match number of sources ({self.n_sources})"
                )
        else:
            # Create zero angles array
            self.angles = np.zeros(self.n_sources, dtype=np.float64)

    @property
    def n_sources(self) -> int:
        """Number of sources."""
        return len(self.positions)

    def __repr__(self) -> str:
        """String representation."""
        type_str = "point" if self.source_type == 1 else "plane wave"
        return (
            f"Sources(n={self.n_sources}, type={type_str}, "
            f"bounds=({self.positions[:, 0].min():.2f}-{self.positions[:, 0].max():.2f}, "
            f"{self.positions[:, 1].min():.2f}-{self.positions[:, 1].max():.2f}))"
        )


@dataclass
class Receivers:
    """Receiver configuration for wavefront tracking.

    Attributes:
        positions: Array of receiver positions, shape (n_receivers, 2) with columns [x, y]

    Example:
        >>> # Regular grid of receivers
        >>> x = np.linspace(110, 160, 50)
        >>> y = np.linspace(-45, -10, 40)
        >>> xx, yy = np.meshgrid(x, y)
        >>> positions = np.column_stack([xx.ravel(), yy.ravel()])
        >>> receivers = Receivers(positions)

        >>> # Or manually specified receivers
        >>> positions = np.array([
        ...     [120.0, -40.0],
        ...     [130.0, -30.0],
        ...     [140.0, -20.0]
        ... ])
        >>> receivers = Receivers(positions)
    """

    positions: np.ndarray

    def __post_init__(self):
        """Validate receiver parameters."""
        # Ensure positions is numpy array
        if not isinstance(self.positions, np.ndarray):
            self.positions = np.array(self.positions, dtype=np.float64)

        if self.positions.dtype != np.float64:
            self.positions = self.positions.astype(np.float64)

        # Validate shape
        if self.positions.ndim == 1:
            # Single receiver provided as 1D array [x, y]
            if len(self.positions) != 2:
                raise ValueError("Single receiver must be [x, y]")
            self.positions = self.positions.reshape(1, 2)
        elif self.positions.ndim == 2:
            if self.positions.shape[1] != 2:
                raise ValueError(f"positions must have shape (n, 2), got {self.positions.shape}")
        else:
            raise ValueError(f"positions must be 1D or 2D array, got shape {self.positions.shape}")

    @property
    def n_receivers(self) -> int:
        """Number of receivers."""
        return len(self.positions)

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"Receivers(n={self.n_receivers}, "
            f"bounds=({self.positions[:, 0].min():.2f}-{self.positions[:, 0].max():.2f}, "
            f"{self.positions[:, 1].min():.2f}-{self.positions[:, 1].max():.2f}))"
        )


def create_constant_velocity_model(
    nx: int,
    ny: int,
    x0: float,
    y0: float,
    dx: float,
    dy: float,
    velocity: float,
    cushion_nodes: int = 3
) -> VelocityModel2D:
    """Create a constant velocity model.

    Args:
        nx, ny: Grid dimensions
        x0, y0: Lower-left corner coordinates
        dx, dy: Grid spacings
        velocity: Constant velocity value
        cushion_nodes: Number of cushion nodes

    Returns:
        VelocityModel2D instance

    Example:
        >>> model = create_constant_velocity_model(
        ...     nx=100, ny=80,
        ...     x0=110.0, y0=-45.0,
        ...     dx=0.5, dy=0.5,
        ...     velocity=3.5
        ... )
    """
    velocities = np.full((nx, ny), velocity, dtype=np.float32)
    return VelocityModel2D(velocities, x0, y0, dx, dy, cushion_nodes)


def write_sources_file(sources: Sources, filename: str):
    """Write sources to binary file compatible with SWMP.

    Args:
        sources: Sources object
        filename: Output filename

    Raises:
        ValueError: If too many sources (>1 million)

    Example:
        >>> sources = Sources(positions=np.array([[135.0, -25.0]]))
        >>> write_sources_file(sources, 'sources.dat')
    """
    import struct
    from pathlib import Path

    # Security: Validate number of sources
    MAX_SOURCES = 1_000_000
    if sources.n_sources > MAX_SOURCES:
        raise ValueError(
            f"Too many sources: {sources.n_sources:,} (max {MAX_SOURCES:,})"
        )

    output_path = Path(filename).resolve()

    with open(output_path, 'wb') as f:
        # Write number of sources as integer (4 bytes)
        f.write(struct.pack('i', sources.n_sources))

        # Write source type as integer (4 bytes)
        f.write(struct.pack('i', sources.source_type))

        # Write positions as doubles (8 bytes each)
        for i in range(sources.n_sources):
            f.write(struct.pack('d', sources.positions[i, 0]))  # x/lon
            f.write(struct.pack('d', sources.positions[i, 1]))  # y/lat

        # Write angles as doubles (8 bytes each)
        for i in range(sources.n_sources):
            f.write(struct.pack('d', sources.angles[i]))


def write_receivers_file(receivers: Receivers, filename: str):
    """Write receivers to binary file compatible with SWMP.

    Args:
        receivers: Receivers object
        filename: Output filename

    Raises:
        ValueError: If too many receivers (>10 million)

    Example:
        >>> receivers = Receivers(positions=np.array([[120.0, -30.0]]))
        >>> write_receivers_file(receivers, 'receivers.dat')
    """
    import struct
    from pathlib import Path

    # Security: Validate number of receivers
    MAX_RECEIVERS = 10_000_000
    if receivers.n_receivers > MAX_RECEIVERS:
        raise ValueError(
            f"Too many receivers: {receivers.n_receivers:,} (max {MAX_RECEIVERS:,})"
        )

    output_path = Path(filename).resolve()

    with open(output_path, 'wb') as f:
        # Write number of receivers as integer (4 bytes)
        f.write(struct.pack('i', receivers.n_receivers))

        # Write positions as doubles (8 bytes each)
        for i in range(receivers.n_receivers):
            f.write(struct.pack('d', receivers.positions[i, 0]))  # x/lon
            f.write(struct.pack('d', receivers.positions[i, 1]))  # y/lat


def create_gradient_velocity_model(
    nx: int,
    ny: int,
    x0: float,
    y0: float,
    dx: float,
    dy: float,
    v0: float,
    gradient_x: float = 0.0,
    gradient_y: float = 0.0,
    cushion_nodes: int = 3
) -> VelocityModel2D:
    """Create a velocity model with linear gradient.

    Args:
        nx, ny: Grid dimensions
        x0, y0: Lower-left corner coordinates
        dx, dy: Grid spacings
        v0: Velocity at (x0, y0)
        gradient_x: Velocity gradient in X direction (dv/dx)
        gradient_y: Velocity gradient in Y direction (dv/dy)
        cushion_nodes: Number of cushion nodes

    Returns:
        VelocityModel2D instance

    Example:
        >>> # Velocity increasing to the east
        >>> model = create_gradient_velocity_model(
        ...     nx=100, ny=80,
        ...     x0=110.0, y0=-45.0,
        ...     dx=0.5, dy=0.5,
        ...     v0=3.0,
        ...     gradient_x=0.01  # +0.01 km/s per degree longitude
        ... )
    """
    x = np.arange(nx) * dx + x0
    y = np.arange(ny) * dy + y0
    xx, yy = np.meshgrid(x, y, indexing='ij')

    velocities = v0 + gradient_x * (xx - x0) + gradient_y * (yy - y0)
    return VelocityModel2D(velocities.astype(np.float32), x0, y0, dx, dy, cushion_nodes)


def write_velocity_model_file(model: VelocityModel2D, filename: str):
    """Write velocity model to binary file compatible with SWMP.

    Args:
        model: VelocityModel2D object
        filename: Output filename

    Raises:
        ValueError: If model is too large (>1GB)

    Example:
        >>> model = create_constant_velocity_model(100, 80, 110, -45, 0.5, 0.5, 3.5)
        >>> write_velocity_model_file(model, 'model.vel')
    """
    import struct
    from pathlib import Path

    # Security: Prevent path traversal and validate file size
    output_path = Path(filename).resolve()

    # Estimate file size: header + velocity data (float32)
    header_size = 4 * 8 + 3 * 4  # 4 doubles, 3 ints
    data_size = model.nx * model.ny * 4  # float32 = 4 bytes
    estimated_size = header_size + data_size

    # Check size limit (1GB = 1,073,741,824 bytes)
    MAX_FILE_SIZE = 1_073_741_824
    if estimated_size > MAX_FILE_SIZE:
        raise ValueError(
            f"Model too large: {estimated_size:,} bytes (max {MAX_FILE_SIZE:,} bytes). "
            f"Model has {model.nx} x {model.ny} grid points."
        )

    with open(output_path, 'wb') as f:
        # Write header: x0, y0, nx, ny, dx, dy, cn (all as doubles except nx, ny, cn as ints)
        f.write(struct.pack('d', model.x0))
        f.write(struct.pack('d', model.y0))
        f.write(struct.pack('i', model.nx))
        f.write(struct.pack('i', model.ny))
        f.write(struct.pack('d', model.dx))
        f.write(struct.pack('d', model.dy))
        f.write(struct.pack('i', model.cushion_nodes))

        # Write velocity values in Fortran order (column-major)
        vel_fortran = np.asfortranarray(model.velocities)
        for val in vel_fortran.ravel():
            f.write(struct.pack('f', val))  # float32


def write_config_file(
    filename: str,
    velocity_file: str,
    sources_file: str,
    receivers_file: str,
    output_dir: str = "output",
    options: Optional['TrackerOptions'] = None
):
    """Write a RAT configuration file.

    Args:
        filename: Output configuration filename
        velocity_file: Path to velocity model file
        sources_file: Path to sources file
        receivers_file: Path to receivers file
        output_dir: Directory for output files
        options: Optional TrackerOptions (uses defaults if None)

    Example:
        >>> write_config_file(
        ...     'config.in',
        ...     velocity_file='model.vel',
        ...     sources_file='sources.dat',
        ...     receivers_file='receivers.dat'
        ... )
    """
    from .options import TrackerOptions

    if options is None:
        options = TrackerOptions()

    # Create output directory if it doesn't exist
    from pathlib import Path
    Path(output_dir).mkdir(exist_ok=True)

    with open(filename, 'w') as f:
        f.write("!-------------------------------------------------------------------------------\n")
        f.write("! input files\n")
        f.write("!-------------------------------------------------------------------------------\n")
        f.write(f"{velocity_file}  ! velocity model\n")
        f.write(f"{sources_file}  ! sources\n")
        f.write(f"{receivers_file}  ! receivers\n")
        f.write("!-------------------------------------------------------------------------------\n")
        f.write("! ray tracing parameters\n")
        f.write("!-------------------------------------------------------------------------------\n")
        f.write(f"{options.dt}  ! ode solver time step size\n")
        f.write(f"{options.max_iterations}  ! maximum number of iterations\n")
        f.write(f"{options.n_bichar_nodes}  ! number of nodes on the bicharacteristic strip\n")
        f.write(f"{options.computation_mode}  ! kinematic (1) kinematic and spreading (2)\n")
        f.write(f"{options.ode_solver}  ! 4th ordRK (1) 5th ord RK (2) 5th ord adaptive RK (3)\n")
        f.write(f"{options.interpolator}  ! interpolator (1) linear) (2) weighted average\n")
        f.write(f"{options.max_arrivals}  ! maximum number of arrivals for a receiver\n")
        f.write(f"{options.coordinate_system}  ! (1) Cartesian cubic b splines (2) Spherical splines\n")
        f.write(f"{options.earth_radius}  ! radius of the earth\n")
        f.write(f"{1 if options.source_specific_receivers else 0}  ! source specific receiver configuration no/yes 0/1\n")

        # For now, we'll skip the receiver mode file if source_specific_receivers is False
        if options.source_specific_receivers and options.receiver_mode_file:
            f.write(f"{options.receiver_mode_file}  ! receiver configuration\n")
        else:
            f.write(f"recmode.dat  ! receiver configuration (not used)\n")

        f.write("!-------------------------------------------------------------------------------\n")
        f.write("! output\n")
        f.write("!-------------------------------------------------------------------------------\n")
        f.write(f"{options.wavefront_interval}  ! store every n-th wavefront\n")
        f.write(f"{1 if options.extract_raypaths else 0}  ! ray path extraction no/yes 0/1\n")
        f.write(f"{options.raypath_storage}  ! ray storage (0) no (1) one file (2) sorted by arrival (3) both\n")
        f.write(f"{1 if options.extract_frechet else 0}  ! frechet derivatives no/yes 0/1\n")
        f.write("!-------------------------------------------------------------------------------\n")
        f.write("! output files\n")
        f.write("!-------------------------------------------------------------------------------\n")
        f.write(f"{output_dir}/arrivals.dat  ! arrival time\n")
        f.write(f"{output_dir}/wafpos.dat  ! wavefront positions\n")
        f.write(f"{output_dir}/raypaths.dat  ! raypaths\n")
        f.write(f"{output_dir}/frechet.hdr  ! frechet header\n")
        f.write(f"{output_dir}/frechet.mat  ! frechet matrix\n")
        f.write(f"{output_dir}/frechet.rai  ! frechet ray index\n")
        f.write(f"{output_dir}/rat.sum  ! rat info file\n")
