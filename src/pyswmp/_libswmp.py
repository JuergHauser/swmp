"""Low-level ctypes wrapper for SWMP Fortran library.

This module provides the ctypes interface to the compiled Fortran library,
defining function signatures and basic error handling.
"""

import ctypes
import sys
from pathlib import Path
import numpy as np
import numpy.ctypeslib as npct

# Error codes from m_error.f90
SWMP_SUCCESS = 0
SWMP_ERROR_FILE = -1
SWMP_ERROR_MEMORY = -2
SWMP_ERROR_INVALID_PARAM = -3


class SWMPError(Exception):
    """Exception raised when SWMP library returns an error."""
    pass


def _find_library():
    """Locate the compiled SWMP library.

    Returns:
        Path to library file

    Raises:
        FileNotFoundError: If library cannot be located
    """
    # Look for library in package directory (installed case)
    pkg_dir = Path(__file__).parent

    # Platform-specific library names
    if sys.platform == 'darwin':
        lib_patterns = ['libswmp.dylib', 'libswmp.so']
    elif sys.platform == 'win32':
        lib_patterns = ['swmp.dll', 'libswmp.dll']
    else:  # Linux and others
        lib_patterns = ['libswmp.so']

    # Also check for Python extension module name
    from sysconfig import get_config_var
    ext_suffix = get_config_var('EXT_SUFFIX')
    if ext_suffix:
        lib_patterns.insert(0, f'libswmp{ext_suffix}')

    for pattern in lib_patterns:
        lib_path = pkg_dir / pattern
        if lib_path.exists():
            return str(lib_path)

    # Check build directory (development case)
    build_dirs = [
        pkg_dir.parent.parent / 'build',
        pkg_dir.parent.parent / 'build' / 'lib',
    ]

    for build_dir in build_dirs:
        if build_dir.exists():
            for pattern in lib_patterns:
                for lib_path in build_dir.rglob(pattern):
                    if lib_path.exists():
                        return str(lib_path)

    raise FileNotFoundError(
        f"Could not find SWMP library. Searched for: {lib_patterns}\n"
        f"Package dir: {pkg_dir}\n"
        "Please ensure the library is built and installed."
    )


class LibSWMP:
    """Wrapper for the SWMP Fortran library.

    This class loads the compiled library and defines ctypes signatures
    for all C-callable functions. Each instance maintains independent state.

    Example:
        >>> lib = LibSWMP()
        >>> lib.read_rat_conf(b"config.in")
        >>> lib.read_model(b"model.vel")
        >>> lib.forward()
    """

    def __init__(self, lib_path=None):
        """Initialize library wrapper.

        Args:
            lib_path: Optional explicit path to library file
        """
        if lib_path is None:
            lib_path = _find_library()

        self._lib = ctypes.CDLL(str(lib_path))
        self._setup_functions()

    def _setup_functions(self):
        """Define ctypes signatures for all library functions."""
        lib = self._lib

        # ================================================================
        # Error Handling
        # ================================================================
        lib.get_last_error.argtypes = [
            ctypes.POINTER(ctypes.c_int32),  # error code
            ctypes.c_char_p,                  # message buffer
            ctypes.c_int32,                   # buffer length
        ]
        lib.get_last_error.restype = ctypes.c_int32

        lib.clear_error.argtypes = []
        lib.clear_error.restype = ctypes.c_int32

        # ================================================================
        # Configuration and I/O
        # ================================================================
        # Note: rat_read_conf loads everything from config file
        # (model, sources, receivers, etc.)
        lib.rat_read_conf.argtypes = [ctypes.c_char_p, ctypes.c_int32]
        lib.rat_read_conf.restype = None

        # ================================================================
        # Forward Modeling
        # ================================================================
        lib.forward.argtypes = []
        lib.forward.restype = ctypes.c_int32

        lib.forward_single_source.argtypes = [ctypes.c_int32]
        lib.forward_single_source.restype = ctypes.c_int32

        # ================================================================
        # Control Functions
        # ================================================================
        lib.get_number_of_sources.argtypes = [ctypes.POINTER(ctypes.c_int32)]
        lib.get_number_of_sources.restype = ctypes.c_int32

        lib.set_enable_file_output.argtypes = [ctypes.c_int32]
        lib.set_enable_file_output.restype = ctypes.c_int32

        # Parameter setters
        lib.set_dt.argtypes = [ctypes.c_float]
        lib.set_dt.restype = None

        lib.set_maxit.argtypes = [ctypes.c_int32]
        lib.set_maxit.restype = None

        lib.set_nsenod.argtypes = [ctypes.c_int32]
        lib.set_nsenod.restype = None

        lib.set_mode.argtypes = [ctypes.c_int32]
        lib.set_mode.restype = None

        lib.set_solver.argtypes = [ctypes.c_int32]
        lib.set_solver.restype = None

        lib.set_interp.argtypes = [ctypes.c_int32]
        lib.set_interp.restype = None

        lib.set_mar.argtypes = [ctypes.c_int32]
        lib.set_mar.restype = None

        lib.set_velint.argtypes = [ctypes.c_int32]
        lib.set_velint.restype = None

        lib.set_rearth.argtypes = [ctypes.c_float]
        lib.set_rearth.restype = None

        # ================================================================
        # In-Memory Data Extraction - Arrivals
        # ================================================================
        lib.get_number_of_arrivals_in_memory.argtypes = [
            ctypes.POINTER(ctypes.c_int32)
        ]
        lib.get_number_of_arrivals_in_memory.restype = ctypes.c_int32

        lib.get_arrivals_from_memory.argtypes = [
            npct.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),  # receivers
            npct.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),  # arrival_nums
            npct.ndpointer(dtype=np.float32, flags='C_CONTIGUOUS'),  # times
            npct.ndpointer(dtype=np.float32, flags='C_CONTIGUOUS'),  # azimuths
            npct.ndpointer(dtype=np.float32, flags='C_CONTIGUOUS'),  # spreading
            ctypes.c_int32,  # array size
        ]
        lib.get_arrivals_from_memory.restype = ctypes.c_int32

        # ================================================================
        # In-Memory Data Extraction - Raypaths
        # ================================================================
        lib.get_raypath_count_in_memory.argtypes = [
            ctypes.POINTER(ctypes.c_int32),  # n_paths
            ctypes.POINTER(ctypes.c_int32),  # total_points
        ]
        lib.get_raypath_count_in_memory.restype = ctypes.c_int32

        lib.get_raypath_metadata_from_memory.argtypes = [
            npct.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),  # receivers
            npct.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),  # arrival_nums
            npct.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),  # npts
            ctypes.c_int32,  # n_paths
        ]
        lib.get_raypath_metadata_from_memory.restype = ctypes.c_int32

        lib.get_raypath_positions_from_memory.argtypes = [
            npct.ndpointer(dtype=np.float32, flags='C_CONTIGUOUS'),  # positions
            ctypes.c_int32,  # total_points
        ]
        lib.get_raypath_positions_from_memory.restype = ctypes.c_int32

        # ================================================================
        # In-Memory Data Extraction - Wavefronts
        # ================================================================
        # Note: Wavefront extraction functions not yet implemented in Fortran
        # Skip for now to allow basic functionality to work
        # TODO: Implement get_wavefront_count_in_memory, etc. in m_swmp.f90

        # ================================================================
        # In-Memory Data Setting - Model, Sources, Receivers (File-Free API)
        # ================================================================
        lib.set_velocity_model_from_memory.argtypes = [
            ctypes.c_float,  # x0
            ctypes.c_float,  # y0
            ctypes.c_int32,  # nx
            ctypes.c_int32,  # ny
            ctypes.c_float,  # dx
            ctypes.c_float,  # dy
            ctypes.c_int32,  # cn
            npct.ndpointer(dtype=np.float32, flags='C_CONTIGUOUS'),  # values
        ]
        lib.set_velocity_model_from_memory.restype = None

        lib.set_sources_from_memory.argtypes = [
            ctypes.c_int32,  # n_sources
            ctypes.c_int32,  # source_type
            npct.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),  # positions (Fortran order)
        ]
        lib.set_sources_from_memory.restype = None

        lib.set_receivers_from_memory.argtypes = [
            ctypes.c_int32,  # n_receivers
            ctypes.c_int32,  # max_arrivals
            npct.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),  # positions (Fortran order)
        ]
        lib.set_receivers_from_memory.restype = None

        # ================================================================
        # Ray/Frechet Control Setters
        # ================================================================
        lib.set_do_rays.argtypes = [ctypes.c_int32]
        lib.set_do_rays.restype = None

        lib.set_do_frechet.argtypes = [ctypes.c_int32]
        lib.set_do_frechet.restype = None

        # ================================================================
        # Jacobian (Frechet Matrix) - In-Memory Extraction
        # ================================================================
        lib.get_sparse_jacobian_size.argtypes = [
            ctypes.POINTER(ctypes.c_int32),  # nr (rows)
            ctypes.POINTER(ctypes.c_int32),  # nc (cols)
            ctypes.POINTER(ctypes.c_int32),  # nnz (non-zeros)
        ]
        lib.get_sparse_jacobian_size.restype = None

        lib.get_sparse_jacobian.argtypes = [
            npct.ndpointer(dtype=np.float32, flags='C_CONTIGUOUS'),  # jrow
            npct.ndpointer(dtype=np.float32, flags='C_CONTIGUOUS'),  # jcol
            npct.ndpointer(dtype=np.float32, flags='C_CONTIGUOUS'),  # jval
            ctypes.c_int32,  # n (nnz)
        ]
        lib.get_sparse_jacobian.restype = None

    def check_error(self, status):
        """Check return status and raise exception if error occurred.

        Args:
            status: Return code from library function

        Raises:
            SWMPError: If status indicates an error
        """
        if status != SWMP_SUCCESS:
            # Get detailed error message from library
            error_code = ctypes.c_int32()
            msg_buf = ctypes.create_string_buffer(256)
            self._lib.get_last_error(
                ctypes.byref(error_code),
                msg_buf,
                256
            )
            msg = msg_buf.value.decode('utf-8', errors='replace').strip()

            if not msg:
                msg = f"Unknown error (code {status})"

            raise SWMPError(f"SWMP error {error_code.value}: {msg}")

    # ================================================================
    # High-Level Wrappers
    # ================================================================

    def read_configuration(self, config_file):
        """Read configuration from file.

        This loads everything specified in the config file:
        model, sources, receivers, and all ray tracing parameters.

        Args:
            config_file: Path to configuration file
        """
        if isinstance(config_file, str):
            config_file = config_file.encode('utf-8')
        # rat_read_conf doesn't return status, but we can check for exceptions
        try:
            self._lib.rat_read_conf(config_file, len(config_file))
        except Exception as e:
            raise SWMPError(f"Failed to read configuration: {e}")

    def run_forward(self):
        """Run forward modeling for all sources (sequential).

        This is the backward-compatible sequential execution path.
        """
        status = self._lib.forward()
        self.check_error(status)

    def run_single_source(self, source_id):
        """Run forward modeling for a single source.

        Args:
            source_id: Source ID (1-based indexing)
        """
        status = self._lib.forward_single_source(ctypes.c_int32(source_id))
        self.check_error(status)

    def get_n_sources(self):
        """Get number of sources.

        Returns:
            Number of sources
        """
        n = ctypes.c_int32()
        status = self._lib.get_number_of_sources(ctypes.byref(n))
        self.check_error(status)
        return n.value

    def set_file_output(self, enable):
        """Enable or disable file output.

        Args:
            enable: True to enable file output, False for memory-only
        """
        status = self._lib.set_enable_file_output(ctypes.c_int32(1 if enable else 0))
        self.check_error(status)

    def get_arrivals(self, source_id):
        """Extract arrival data from memory for given source.

        Args:
            source_id: Source ID for which arrivals were computed

        Returns:
            Tuple of (source_ids, receiver_ids, arrival_nums, times, azimuths, spreading_factors)
        """
        # Get count
        n = ctypes.c_int32()
        status = self._lib.get_number_of_arrivals_in_memory(ctypes.byref(n))
        self.check_error(status)

        if n.value == 0:
            return (
                np.array([], dtype=np.int32),  # source_ids
                np.array([], dtype=np.int32),  # receiver_ids
                np.array([], dtype=np.int32),  # arrival_nums
                np.array([], dtype=np.float32),  # times
                np.array([], dtype=np.float32),  # azimuths
                np.array([], dtype=np.float32),  # spreading
            )

        # Allocate arrays
        receivers = np.zeros(n.value, dtype=np.int32)
        arrival_nums = np.zeros(n.value, dtype=np.int32)
        times = np.zeros(n.value, dtype=np.float32)
        azimuths = np.zeros(n.value, dtype=np.float32)
        spreading = np.zeros(n.value, dtype=np.float32)

        # Extract data
        status = self._lib.get_arrivals_from_memory(
            receivers, arrival_nums, times, azimuths, spreading, n.value
        )
        self.check_error(status)

        # Add source_id to all arrivals
        source_ids = np.full(n.value, source_id, dtype=np.int32)

        return source_ids, receivers, arrival_nums, times, azimuths, spreading

    def get_raypaths(self, source_id):
        """Extract raypath data from memory for given source.

        Args:
            source_id: Source ID for which raypaths were computed

        Returns:
            List of dicts with keys: source, receiver, arrival, path (ndarray)
        """
        # Get counts
        n_paths = ctypes.c_int32()
        total_pts = ctypes.c_int32()
        status = self._lib.get_raypath_count_in_memory(
            ctypes.byref(n_paths), ctypes.byref(total_pts)
        )
        self.check_error(status)

        if n_paths.value == 0:
            return []

        # Get metadata
        receivers = np.zeros(n_paths.value, dtype=np.int32)
        arrival_nums = np.zeros(n_paths.value, dtype=np.int32)
        npts = np.zeros(n_paths.value, dtype=np.int32)

        status = self._lib.get_raypath_metadata_from_memory(
            receivers, arrival_nums, npts, n_paths.value
        )
        self.check_error(status)

        # Get positions (flattened)
        positions = np.zeros(total_pts.value * 2, dtype=np.float32)
        status = self._lib.get_raypath_positions_from_memory(
            positions, total_pts.value
        )
        self.check_error(status)

        # Reconstruct raypaths
        raypaths = []
        pos_idx = 0
        for i in range(n_paths.value):
            n = npts[i]
            path = positions[pos_idx*2:(pos_idx+n)*2].reshape(n, 2)
            raypaths.append({
                'source': source_id,
                'receiver': int(receivers[i]),
                'arrival': int(arrival_nums[i]),
                'path': path,
            })
            pos_idx += n

        return raypaths

    def get_wavefronts(self, source_id):
        """Extract wavefront data from memory for given source.

        Args:
            source_id: Source ID for which wavefronts were computed

        Returns:
            List of wavefront arrays (each is ndarray of shape (n, 2))

        Note:
            Wavefront extraction not yet implemented in Fortran.
            This method always returns an empty list.
        """
        # TODO: Implement wavefront extraction functions in m_swmp.f90
        return []

    # ================================================================
    # Parameter Configuration Methods
    # ================================================================

    def set_dt(self, dt):
        """Set ODE solver time step.

        Args:
            dt: Time step size
        """
        self._lib.set_dt(ctypes.c_float(dt))

    def set_max_iterations(self, maxit):
        """Set maximum number of iterations.

        Args:
            maxit: Maximum iterations
        """
        self._lib.set_maxit(ctypes.c_int32(maxit))

    def set_n_bichar_nodes(self, nsenod):
        """Set number of nodes on bicharacteristic strip.

        Args:
            nsenod: Number of nodes
        """
        self._lib.set_nsenod(ctypes.c_int32(nsenod))

    def set_computation_mode(self, mode):
        """Set computation mode.

        Args:
            mode: Computation mode
                1 = Kinematic ray tracing only
                2 = Kinematic ray tracing + spreading factors
        """
        self._lib.set_mode(ctypes.c_int32(mode))

    def set_ode_solver(self, solver):
        """Set ODE solver type.

        Args:
            solver: Solver type
                1 = 4th order Runge-Kutta
                2 = 5th order Runge-Kutta
                3 = 5th order adaptive Runge-Kutta
        """
        self._lib.set_solver(ctypes.c_int32(solver))

    def set_interpolator(self, interp):
        """Set velocity interpolation method.

        Args:
            interp: Interpolation method
                1 = Linear interpolation
                2 = Weighted average
        """
        self._lib.set_interp(ctypes.c_int32(interp))

    def set_max_arrivals(self, mar):
        """Set maximum number of arrivals per receiver.

        Args:
            mar: Maximum number of arrivals
        """
        self._lib.set_mar(ctypes.c_int32(mar))

    def set_coordinate_system(self, velint):
        """Set coordinate system type.

        Args:
            velint: Coordinate system
                1 = Cartesian cubic B-splines
                2 = Spherical splines
        """
        self._lib.set_velint(ctypes.c_int32(velint))

    def set_earth_radius(self, rearth):
        """Set Earth radius for spherical coordinates.

        Args:
            rearth: Earth radius in km (default: 6371.0)
        """
        self._lib.set_rearth(ctypes.c_float(rearth))

    def set_do_rays(self, flag):
        """Enable or disable ray path computation in Fortran.

        Args:
            flag: True to compute raypaths, False to skip
        """
        self._lib.set_do_rays(ctypes.c_int32(1 if flag else 0))

    def set_do_frechet(self, flag):
        """Enable or disable Fréchet matrix computation in Fortran.

        Args:
            flag: True to compute Fréchet matrix, False to skip
        """
        self._lib.set_do_frechet(ctypes.c_int32(1 if flag else 0))

    def get_jacobian(self):
        """Extract sparse Jacobian (Fréchet matrix) from memory.

        Returns COO-format sparse data after forward_single_source() has been
        called with do_frechet enabled.

        Returns:
            Tuple of (rows, cols, vals, shape) where:
                rows: int32 array of row indices
                cols: int32 array of column indices
                vals: float64 array of values
                shape: (n_rows, n_cols) tuple
            Returns None if no Jacobian is available (nnz == 0).
        """
        nr = ctypes.c_int32()
        nc = ctypes.c_int32()
        nnz = ctypes.c_int32()
        self._lib.get_sparse_jacobian_size(
            ctypes.byref(nr), ctypes.byref(nc), ctypes.byref(nnz)
        )

        if nnz.value == 0:
            return None

        # Fortran returns all as float32 (COO format)
        jrow = np.zeros(nnz.value, dtype=np.float32)
        jcol = np.zeros(nnz.value, dtype=np.float32)
        jval = np.zeros(nnz.value, dtype=np.float32)

        self._lib.get_sparse_jacobian(jrow, jcol, jval, nnz.value)

        # Convert row/col to int, col from 1-based to 0-based
        rows = jrow.astype(np.int32) - 1  # 1-based to 0-based
        cols = jcol.astype(np.int32) - 1  # 1-based to 0-based

        return rows, cols, jval.astype(np.float64), (nr.value, nc.value)

    def apply_options(self, options):
        """Apply TrackerOptions to configure all parameters.

        Args:
            options: TrackerOptions instance with configuration

        Example:
            >>> lib = LibSWMP()
            >>> opts = TrackerOptions(dt=2.0, max_iterations=1000)
            >>> lib.apply_options(opts)
        """
        self.set_dt(options.dt)
        self.set_max_iterations(options.max_iterations)
        self.set_n_bichar_nodes(options.n_bichar_nodes)
        self.set_computation_mode(options.computation_mode)
        self.set_ode_solver(options.ode_solver)
        self.set_interpolator(options.interpolator)
        self.set_max_arrivals(options.max_arrivals)
        self.set_coordinate_system(options.coordinate_system)
        self.set_earth_radius(options.earth_radius)
        # Set ray/frechet computation flags
        self.set_do_rays(options.extract_raypaths or options.extract_frechet)
        self.set_do_frechet(options.extract_frechet)

    # ================================================================
    # Model, Sources, Receivers Configuration
    # ================================================================

    def set_velocity_model(self, model):
        """Set velocity model from VelocityModel2D object.

        This uses the new file-free API that initializes the complete
        velocity model structure in memory.

        Args:
            model: VelocityModel2D instance

        Example:
            >>> from pyswmp.data import VelocityModel2D
            >>> import numpy as np
            >>> velocities = np.full((100, 80), 3.5, dtype=np.float32)
            >>> model = VelocityModel2D(velocities, 110.0, -45.0, 0.5, 0.5)
            >>> lib.set_velocity_model(model)
        """
        from .data import VelocityModel2D

        if not isinstance(model, VelocityModel2D):
            raise TypeError("model must be a VelocityModel2D instance")

        # Ensure velocities are float32
        velocities = model.velocities.astype(np.float32)

        # Add cushion nodes if needed
        # Fortran expects array of size (nx+2*cn, ny+2*cn)
        if model.cushion_nodes > 0:
            # Pad with edge replication
            padded = np.pad(
                velocities,
                pad_width=model.cushion_nodes,
                mode='edge'
            )
        else:
            padded = velocities

        # Flatten velocity array in Fortran order (column-major)
        vel_flat = np.asfortranarray(padded).ravel(order='F')

        # Set complete model with metadata
        self._lib.set_velocity_model_from_memory(
            ctypes.c_float(model.x0),
            ctypes.c_float(model.y0),
            ctypes.c_int32(model.nx),
            ctypes.c_int32(model.ny),
            ctypes.c_float(model.dx),
            ctypes.c_float(model.dy),
            ctypes.c_int32(model.cushion_nodes),
            vel_flat
        )

    def set_sources(self, sources):
        """Set sources from Sources object (in-memory, no file).

        Args:
            sources: Sources instance

        Example:
            >>> from pyswmp.data import Sources
            >>> import numpy as np
            >>> positions = np.array([[135.0, -25.0], [140.0, -30.0]])
            >>> sources = Sources(positions, source_type=1)
            >>> lib.set_sources(sources)
        """
        from .data import Sources

        if not isinstance(sources, Sources):
            raise TypeError("sources must be a Sources instance")

        # Ensure positions are Fortran-contiguous (column-major)
        # Fortran expects dimension(n, 2) in column-major order
        positions = np.asfortranarray(sources.positions, dtype=np.float64)

        # Set sources in Fortran
        self._lib.set_sources_from_memory(
            ctypes.c_int32(sources.n_sources),
            ctypes.c_int32(sources.source_type),
            positions
        )

    def set_receivers(self, receivers, max_arrivals=10):
        """Set receivers from Receivers object (in-memory, no file).

        Args:
            receivers: Receivers instance
            max_arrivals: Maximum number of arrivals per receiver (default: 10)

        Example:
            >>> from pyswmp.data import Receivers
            >>> import numpy as np
            >>> positions = np.array([[120.0, -30.0], [130.0, -25.0]])
            >>> receivers = Receivers(positions)
            >>> lib.set_receivers(receivers, max_arrivals=10)
        """
        from .data import Receivers

        if not isinstance(receivers, Receivers):
            raise TypeError("receivers must be a Receivers instance")

        # Ensure positions are Fortran-contiguous (column-major)
        # Fortran expects dimension(n, 2) in column-major order
        positions = np.asfortranarray(receivers.positions, dtype=np.float64)

        # Set receivers in Fortran
        self._lib.set_receivers_from_memory(
            ctypes.c_int32(receivers.n_receivers),
            ctypes.c_int32(max_arrivals),
            positions
        )

    def get_velocity_model(self):
        """Get current velocity model as VelocityModel2D object.

        Returns:
            VelocityModel2D instance with current model data

        Example:
            >>> model = lib.get_velocity_model()
            >>> print(model.nx, model.ny)
            >>> print(model.velocities.shape)
        """
        from .data import VelocityModel2D

        # Get metadata
        x0 = ctypes.c_float()
        y0 = ctypes.c_float()
        nx = ctypes.c_int32()
        ny = ctypes.c_int32()
        dx = ctypes.c_float()
        dy = ctypes.c_float()
        cn = ctypes.c_int32()

        self._lib.get_model_meta_data(
            ctypes.byref(x0),
            ctypes.byref(y0),
            ctypes.byref(nx),
            ctypes.byref(ny),
            ctypes.byref(dx),
            ctypes.byref(dy),
            ctypes.byref(cn)
        )

        # Get velocity data
        total_size = (nx.value + 2 * cn.value) * (ny.value + 2 * cn.value)
        vel_flat = np.empty(total_size, dtype=np.float32)

        self._lib.get_model_vector(
            vel_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ctypes.byref(x0),
            ctypes.byref(y0),
            ctypes.byref(nx),
            ctypes.byref(ny),
            ctypes.byref(dx),
            ctypes.byref(dy),
            ctypes.byref(cn)
        )

        # Reshape to 2D (without cushion nodes for now)
        velocities = vel_flat[:nx.value * ny.value].reshape((nx.value, ny.value), order='F')

        return VelocityModel2D(
            velocities=velocities,
            x0=x0.value,
            y0=y0.value,
            dx=dx.value,
            dy=dy.value,
            cushion_nodes=cn.value
        )
