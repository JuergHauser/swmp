"""Configuration options for SWMP wavefront tracking.

This module provides a dataclass-based configuration API similar to pyfm2d,
allowing users to configure all ray tracing parameters through Python instead
of configuration files.
"""

from dataclasses import dataclass
from typing import Optional
from pathlib import Path


@dataclass
class TrackerOptions:
    """Configuration options for wavefront tracking.

    This dataclass allows programmatic configuration of all parameters
    that are typically specified in RAT configuration files.

    Physical/Geometric Parameters:
        earth_radius: Radius of Earth in km (default: 6371.0)
        coordinate_system: Coordinate system type
            1 = Cartesian cubic B-splines (default)
            2 = Spherical splines

    Integration Parameters:
        dt: ODE solver time step size (default: 5.0)
        max_iterations: Maximum number of iterations (default: 500)
        n_bichar_nodes: Number of nodes on bicharacteristic strip (default: 75)
        ode_solver: ODE solver type (default: 1)
            1 = 4th order Runge-Kutta
            2 = 5th order Runge-Kutta
            3 = 5th order adaptive Runge-Kutta
        interpolator: Velocity interpolation method (default: 1)
            1 = Linear interpolation
            2 = Weighted average

    Computation Options:
        computation_mode: What to compute (default: 1)
            1 = Kinematic ray tracing only
            2 = Kinematic ray tracing + spreading factors
        max_arrivals: Maximum number of arrivals per receiver (default: 10)
        source_specific_receivers: Use source-specific receiver configuration (default: True)

    Output Options:
        extract_raypaths: Extract ray path data (default: False)
        extract_wavefronts: Extract wavefront data (default: False)
        extract_frechet: Compute Fréchet derivatives (default: False)
        raypath_storage: Ray storage mode (default: 2)
            0 = No storage
            1 = One file
            2 = Sorted by arrival
            3 = Both formats
        wavefront_interval: Store every n-th wavefront (default: 100)

    File Paths:
        velocity_file: Path to velocity model file (optional)
        sources_file: Path to sources file (optional)
        receivers_file: Path to receivers file (optional)
        receiver_mode_file: Path to receiver mode configuration (optional)

    Example:
        >>> # Create default options
        >>> opts = TrackerOptions()

        >>> # Customize parameters
        >>> opts = TrackerOptions(
        ...     dt=2.0,
        ...     max_iterations=1000,
        ...     extract_raypaths=True,
        ...     computation_mode=2  # Include spreading factors
        ... )

        >>> # Use with WaveFrontTracker
        >>> tracker = WaveFrontTracker(options=opts)
        >>> result = tracker.forward()
    """

    # Physical/geometric parameters
    earth_radius: float = 6371.0
    coordinate_system: int = 1

    # Integration parameters
    dt: float = 5.0
    max_iterations: int = 500
    n_bichar_nodes: int = 75
    ode_solver: int = 1
    interpolator: int = 1

    # Computation options
    computation_mode: int = 1
    max_arrivals: int = 10
    source_specific_receivers: bool = True

    # Output options
    extract_raypaths: bool = False
    extract_wavefronts: bool = False
    extract_frechet: bool = False
    raypath_storage: int = 2
    wavefront_interval: int = 100

    # File paths (optional - can be provided separately)
    velocity_file: Optional[str] = None
    sources_file: Optional[str] = None
    receivers_file: Optional[str] = None
    receiver_mode_file: Optional[str] = None

    def __post_init__(self):
        """Validate parameters after initialization."""
        # Validate coordinate system
        if self.coordinate_system not in [1, 2]:
            raise ValueError(f"coordinate_system must be 1 or 2, got {self.coordinate_system}")

        # Validate ODE solver
        if self.ode_solver not in [1, 2, 3]:
            raise ValueError(f"ode_solver must be 1, 2, or 3, got {self.ode_solver}")

        # Validate interpolator
        if self.interpolator not in [1, 2]:
            raise ValueError(f"interpolator must be 1 or 2, got {self.interpolator}")

        # Validate computation mode
        if self.computation_mode not in [1, 2]:
            raise ValueError(f"computation_mode must be 1 or 2, got {self.computation_mode}")

        # Validate raypath storage
        if self.raypath_storage not in [0, 1, 2, 3]:
            raise ValueError(f"raypath_storage must be 0, 1, 2, or 3, got {self.raypath_storage}")

        # Validate positive parameters
        if self.dt <= 0:
            raise ValueError(f"dt must be positive, got {self.dt}")
        if self.max_iterations <= 0:
            raise ValueError(f"max_iterations must be positive, got {self.max_iterations}")
        if self.n_bichar_nodes <= 0:
            raise ValueError(f"n_bichar_nodes must be positive, got {self.n_bichar_nodes}")
        if self.max_arrivals <= 0:
            raise ValueError(f"max_arrivals must be positive, got {self.max_arrivals}")
        if self.earth_radius <= 0:
            raise ValueError(f"earth_radius must be positive, got {self.earth_radius}")
        if self.wavefront_interval <= 0:
            raise ValueError(f"wavefront_interval must be positive, got {self.wavefront_interval}")

        # Convert file paths to Path objects if provided
        if self.velocity_file is not None:
            self.velocity_file = str(Path(self.velocity_file))
        if self.sources_file is not None:
            self.sources_file = str(Path(self.sources_file))
        if self.receivers_file is not None:
            self.receivers_file = str(Path(self.receivers_file))
        if self.receiver_mode_file is not None:
            self.receiver_mode_file = str(Path(self.receiver_mode_file))

    def to_dict(self) -> dict:
        """Convert options to dictionary for worker serialization.

        File path fields are excluded since they are not needed by
        in-memory workers and would cause errors during reconstruction.

        Returns:
            Dictionary with option values (excluding file paths)
        """
        return {
            'earth_radius': self.earth_radius,
            'coordinate_system': self.coordinate_system,
            'dt': self.dt,
            'max_iterations': self.max_iterations,
            'n_bichar_nodes': self.n_bichar_nodes,
            'ode_solver': self.ode_solver,
            'interpolator': self.interpolator,
            'computation_mode': self.computation_mode,
            'max_arrivals': self.max_arrivals,
            'source_specific_receivers': self.source_specific_receivers,
            'extract_raypaths': self.extract_raypaths,
            'extract_wavefronts': self.extract_wavefronts,
            'extract_frechet': self.extract_frechet,
            'raypath_storage': self.raypath_storage,
            'wavefront_interval': self.wavefront_interval,
        }

    @classmethod
    def from_config_file(cls, config_file: str) -> 'TrackerOptions':
        """Parse options from RAT configuration file.

        Args:
            config_file: Path to RAT configuration file

        Returns:
            TrackerOptions instance with parameters from file

        Example:
            >>> opts = TrackerOptions.from_config_file('config.in')
            >>> opts.dt = 2.0  # Override specific parameters
        """
        config_path = Path(config_file)
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_file}")

        # Read configuration file
        with open(config_path, 'r') as f:
            lines = [line.split('!')[0].strip() for line in f if line.strip() and not line.strip().startswith('!')]

        # Parse file paths (first 3 non-empty lines)
        velocity_file = lines[0] if len(lines) > 0 else None
        sources_file = lines[1] if len(lines) > 1 else None
        receivers_file = lines[2] if len(lines) > 2 else None

        # Parse numeric parameters (next section)
        dt = float(lines[3]) if len(lines) > 3 else 5.0
        max_iterations = int(lines[4]) if len(lines) > 4 else 500
        n_bichar_nodes = int(lines[5]) if len(lines) > 5 else 75
        computation_mode = int(lines[6]) if len(lines) > 6 else 1
        ode_solver = int(lines[7]) if len(lines) > 7 else 1
        interpolator = int(lines[8]) if len(lines) > 8 else 1
        max_arrivals = int(lines[9]) if len(lines) > 9 else 10
        coordinate_system = int(lines[10]) if len(lines) > 10 else 1
        earth_radius = float(lines[11]) if len(lines) > 11 else 6371.0
        source_specific_receivers = bool(int(lines[12])) if len(lines) > 12 else True
        receiver_mode_file = lines[13] if len(lines) > 13 else None

        # Parse output parameters
        wavefront_interval = int(lines[14]) if len(lines) > 14 else 100
        extract_raypaths = bool(int(lines[15])) if len(lines) > 15 else False
        raypath_storage = int(lines[16]) if len(lines) > 16 else 2
        extract_frechet = bool(int(lines[17])) if len(lines) > 17 else False

        # Resolve file paths relative to config file directory
        config_dir = config_path.parent
        if velocity_file:
            velocity_file = str(config_dir / velocity_file)
        if sources_file:
            sources_file = str(config_dir / sources_file)
        if receivers_file:
            receivers_file = str(config_dir / receivers_file)
        if receiver_mode_file:
            receiver_mode_file = str(config_dir / receiver_mode_file)

        return cls(
            earth_radius=earth_radius,
            coordinate_system=coordinate_system,
            dt=dt,
            max_iterations=max_iterations,
            n_bichar_nodes=n_bichar_nodes,
            ode_solver=ode_solver,
            interpolator=interpolator,
            computation_mode=computation_mode,
            max_arrivals=max_arrivals,
            source_specific_receivers=source_specific_receivers,
            extract_raypaths=extract_raypaths,
            extract_wavefronts=False,  # Not in config file
            extract_frechet=extract_frechet,
            raypath_storage=raypath_storage,
            wavefront_interval=wavefront_interval,
            velocity_file=velocity_file,
            sources_file=sources_file,
            receivers_file=receivers_file,
            receiver_mode_file=receiver_mode_file,
        )

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"TrackerOptions(\n"
            f"  dt={self.dt}, max_iterations={self.max_iterations},\n"
            f"  computation_mode={self.computation_mode}, ode_solver={self.ode_solver},\n"
            f"  extract_raypaths={self.extract_raypaths}, extract_frechet={self.extract_frechet}\n"
            f")"
        )
