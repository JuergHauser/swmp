"""High-level solver classes for SWMP wavefront tracking.

This module provides the WaveFrontTracker class with support for
pool-based parallel execution following the pyfm2d pattern.
"""

from pathlib import Path
from functools import reduce
import operator
from typing import Optional
import warnings

import numpy as np

from ._libswmp import LibSWMP, SWMPError
from .results import WaveFrontTrackerResult


def _calc_single_source_worker(
    model_dict, sources_dict, receivers_dict, options_dict,
    source_id, enable_raypaths, enable_wavefronts, enable_frechet=False,
):
    """Worker function for parallel source processing.

    Each worker creates an independent LibSWMP instance to avoid shared memory
    issues with Fortran module-level state.

    Args:
        model_dict: Dict with model data (velocities, x0, y0, dx, dy, cushion_nodes)
        sources_dict: Dict with sources data (positions, source_type)
        receivers_dict: Dict with receivers data (positions)
        options_dict: Dict of TrackerOptions parameters
        source_id: Source ID to process (1-based)
        enable_raypaths: Whether to extract raypaths
        enable_wavefronts: Whether to extract wavefronts
        enable_frechet: Whether to extract Frechet/Jacobian matrix

    Returns:
        WaveFrontTrackerResult for this source
    """
    from .options import TrackerOptions
    from .data import VelocityModel2D, Sources, Receivers
    import numpy as np

    lib = LibSWMP()

    try:
        # Reconstruct data objects from dicts
        model = VelocityModel2D(
            velocities=np.asarray(model_dict['velocities'], dtype=np.float32),
            x0=model_dict['x0'],
            y0=model_dict['y0'],
            dx=model_dict['dx'],
            dy=model_dict['dy'],
            cushion_nodes=model_dict['cushion_nodes'],
        )
        lib.set_velocity_model(model)

        sources = Sources(
            positions=np.asarray(sources_dict['positions'], dtype=np.float64),
            source_type=sources_dict['source_type'],
        )
        lib.set_sources(sources)

        receivers = Receivers(
            positions=np.asarray(receivers_dict['positions'], dtype=np.float64),
        )
        max_arrivals = options_dict.get('max_arrivals', 10) if options_dict else 10
        lib.set_receivers(receivers, max_arrivals=max_arrivals)

        # Apply options
        if options_dict is not None:
            opts = TrackerOptions(**options_dict)
            lib.apply_options(opts)

        # Disable file output (memory-only mode)
        lib.set_file_output(False)

        # Ensure Fortran ray/frechet flags are set correctly
        if enable_raypaths or enable_frechet:
            lib.set_do_rays(True)
        if enable_frechet:
            lib.set_do_frechet(True)

        # Run forward model for this source
        lib.run_single_source(source_id)

        # Extract arrivals
        source_ids, receiver_ids, arrival_nums, times, azis, spfs = lib.get_arrivals(source_id)

        # Extract raypaths if requested
        raypaths = []
        if enable_raypaths:
            raypaths = lib.get_raypaths(source_id)

        # Extract wavefronts if requested
        wavefronts = {}
        if enable_wavefronts:
            wf_list = lib.get_wavefronts(source_id)
            if wf_list:
                wavefronts[source_id] = wf_list

        # Extract Jacobian if requested
        jacobian_rows = None
        jacobian_cols = None
        jacobian_vals = None
        jacobian_shape = None
        if enable_frechet:
            jac_data = lib.get_jacobian()
            if jac_data is not None:
                jacobian_rows, jacobian_cols, jacobian_vals, jacobian_shape = jac_data

        return WaveFrontTrackerResult(
            source_ids=source_ids,
            receiver_ids=receiver_ids,
            arrival_numbers=arrival_nums,
            travel_times=times,
            azimuths=azis,
            spreading_factors=spfs,
            raypaths=raypaths,
            wavefronts=wavefronts,
            jacobian_rows=jacobian_rows,
            jacobian_cols=jacobian_cols,
            jacobian_vals=jacobian_vals,
            jacobian_shape=jacobian_shape,
        )

    except Exception as e:
        raise RuntimeError(f"Error processing source {source_id}: {e}") from e


def _worker_wrapper(args):
    """Wrapper for schwimmbad/map-style pools that require a single-argument callable."""
    return _calc_single_source_worker(*args)


class WaveFrontTracker:
    """Wavefront tracker for computing seismic ray paths and travel times.

    This class provides a high-level interface to the SWMP Fortran library,
    with support for parallel execution over sources using pool executors.

    New API (recommended):
        >>> model = pyswmp.create_constant_velocity_model(100, 80, 110, -45, 0.5, 0.5, 3.5)
        >>> opts = pyswmp.TrackerOptions(coordinate_system=2, extract_frechet=True)
        >>> tracker = pyswmp.WaveFrontTracker(model, opts)
        >>> source = pyswmp.Sources(np.array([[135.0, -25.0]]), source_type=1)
        >>> receivers = pyswmp.Receivers(np.array([[120.0, -30.0]]))
        >>> result = tracker.forward(source, receivers)

    Legacy API (file-based, deprecated):
        >>> tracker = WaveFrontTracker(config_file='config.in')
        >>> result = tracker.forward()

    Legacy API (file-free, all at construction):
        >>> tracker = WaveFrontTracker(model=model, sources=sources, receivers=receivers)
        >>> result = tracker.forward()
    """

    def __init__(
        self,
        model_or_config=None,
        options: 'TrackerOptions' = None,
        *,
        config_file: str = None,
        model: 'VelocityModel2D' = None,
        sources: 'Sources' = None,
        receivers: 'Receivers' = None,
    ):
        """Initialize WaveFrontTracker.

        Three modes of initialization:

        1. New API (model + options):
           >>> tracker = WaveFrontTracker(model, options)
           >>> result = tracker.forward(source, receivers)

        2. Legacy file-based (deprecated):
           >>> tracker = WaveFrontTracker('config.in')

        3. Legacy file-free (all at construction):
           >>> tracker = WaveFrontTracker(model=model, sources=sources, receivers=receivers)

        Args:
            model_or_config: VelocityModel2D instance (new API) or config file path (legacy)
            options: TrackerOptions for configuration
            config_file: Path to RAT configuration file (legacy keyword)
            model: VelocityModel2D instance (legacy keyword)
            sources: Sources instance (legacy keyword, or use forward() argument)
            receivers: Receivers instance (legacy keyword, or use forward() argument)
        """
        from .options import TrackerOptions
        from .data import VelocityModel2D, Sources, Receivers

        # Resolve first positional argument
        if model_or_config is not None:
            if isinstance(model_or_config, (str, Path)):
                if config_file is not None:
                    raise ValueError("Cannot provide both positional config path and config_file keyword")
                config_file = str(Path(model_or_config).resolve())
            elif isinstance(model_or_config, VelocityModel2D):
                if model is not None:
                    raise ValueError("Cannot provide both positional model and model keyword")
                model = model_or_config
            else:
                raise TypeError(
                    f"First argument must be VelocityModel2D or config file path, "
                    f"got {type(model_or_config).__name__}"
                )

        self.config_file = str(Path(config_file).resolve()) if config_file else None
        self.options = options

        # Store options as dict for serialization to workers
        self._options_dict = options.to_dict() if options else None

        # Create library instance for metadata queries
        self._lib = LibSWMP()

        # Determine initialization mode
        if config_file:
            # Deprecated file-based initialization
            warnings.warn(
                "File-based initialization via config_file is deprecated. "
                "Use the file-free API with VelocityModel2D, Sources, and Receivers instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            self._lib.read_configuration(self.config_file)
            self._file_free = False
            self._model = None
            self._sources = None
            self._receivers = None
        elif model is not None:
            # File-free initialization: model required, sources/receivers optional
            self._lib.set_velocity_model(model)
            self._model = model
            self._file_free = True

            # Sources/receivers may be provided at construction (legacy) or via forward()
            if sources is not None and receivers is not None:
                self._lib.set_sources(sources)
                max_arrivals = options.max_arrivals if options else 10
                self._lib.set_receivers(receivers, max_arrivals=max_arrivals)
                self._sources = sources
                self._receivers = receivers
            else:
                self._sources = None
                self._receivers = None
        else:
            raise ValueError(
                "Must provide either:\n"
                "  1. model (new API: pass sources/receivers to forward()), or\n"
                "  2. config_file (legacy file-based mode), or\n"
                "  3. model + sources + receivers (legacy file-free mode)"
            )

        # Apply options if provided
        if self.options:
            self._lib.apply_options(self.options)

        # Cache number of sources (0 if sources not yet set)
        if self._sources is not None or config_file:
            self.n_sources = self._lib.get_n_sources()
        else:
            self.n_sources = 0

    def _require_file_free(self, operation):
        """Raise ValueError if not in file-free mode.

        Args:
            operation: Description of the operation requiring file-free mode
        """
        if not self._file_free:
            raise ValueError(
                f"{operation} requires the file-free API. "
                "Construct WaveFrontTracker with a VelocityModel2D instead of config_file."
            )

    def _set_sources_and_receivers(self, source, receivers):
        """Set sources and receivers on the library instance.

        Handles auto-wrapping of numpy arrays into Sources/Receivers objects.

        Args:
            source: Sources object or np.ndarray of positions
            receivers: Receivers object or np.ndarray of positions
        """
        from .data import Sources, Receivers

        # Auto-wrap numpy arrays
        if isinstance(source, np.ndarray):
            if source.ndim == 1:
                source = source.reshape(1, -1)
            source = Sources(positions=source, source_type=1)
        if isinstance(receivers, np.ndarray):
            if receivers.ndim == 1:
                receivers = receivers.reshape(1, -1)
            receivers = Receivers(positions=receivers)

        # Set on library
        self._lib.set_sources(source)
        max_arrivals = self.options.max_arrivals if self.options else 10
        self._lib.set_receivers(receivers, max_arrivals=max_arrivals)

        # Update stored references
        self._sources = source
        self._receivers = receivers
        self.n_sources = source.n_sources

    def _serialize_data(self):
        """Serialize model, sources, receivers for worker processes.

        Returns numpy arrays directly (they pickle efficiently).

        Returns:
            Tuple of (model_dict, sources_dict, receivers_dict)
        """
        model_dict = {
            'velocities': self._model.velocities,
            'x0': self._model.x0,
            'y0': self._model.y0,
            'dx': self._model.dx,
            'dy': self._model.dy,
            'cushion_nodes': self._model.cushion_nodes,
        }

        sources_dict = {
            'positions': self._sources.positions,
            'source_type': self._sources.source_type,
        }

        receivers_dict = {
            'positions': self._receivers.positions,
        }

        return model_dict, sources_dict, receivers_dict

    def _detect_pool_type(self, pool):
        """Detect pool executor type and validate compatibility.

        Args:
            pool: Pool executor instance

        Returns:
            Pool type string ('process', 'schwimmbad', or 'thread')

        Raises:
            ValueError: If pool type is incompatible (e.g., ThreadPoolExecutor)
        """
        pool_type = type(pool).__name__
        module = type(pool).__module__

        # Check for ThreadPoolExecutor (incompatible)
        if 'ThreadPool' in pool_type:
            raise ValueError(
                "ThreadPoolExecutor is not supported. Each worker needs an "
                "independent Fortran library instance, which requires separate "
                "processes. Use ProcessPoolExecutor or schwimmbad pools instead."
            )

        # Detect concurrent.futures ProcessPoolExecutor
        if 'concurrent.futures' in module and 'Process' in pool_type:
            return 'process'

        # Detect schwimmbad pools
        if 'schwimmbad' in module:
            return 'schwimmbad'

        # Unknown pool type - warn but try to proceed with map interface
        warnings.warn(
            f"Unknown pool type: {pool_type} from {module}. "
            "Attempting to use map() interface. "
            "If this fails, use ProcessPoolExecutor or schwimmbad pools."
        )
        return 'schwimmbad'

    def _resolve_extraction_flags(self, extract_raypaths, extract_wavefronts):
        """Resolve extraction flags from arguments vs options.

        Args:
            extract_raypaths: Explicit flag or None (use options)
            extract_wavefronts: Explicit flag or None (use options)

        Returns:
            Tuple of (extract_raypaths, extract_wavefronts, extract_frechet)
        """
        if extract_raypaths is None:
            extract_raypaths = self.options.extract_raypaths if self.options else False
        if extract_wavefronts is None:
            extract_wavefronts = self.options.extract_wavefronts if self.options else False
        extract_frechet = self.options.extract_frechet if self.options else False
        return extract_raypaths, extract_wavefronts, extract_frechet

    def forward(
        self,
        source=None,
        receivers=None,
        *,
        pool=None,
        extract_raypaths: bool = None,
        extract_wavefronts: bool = None,
    ) -> WaveFrontTrackerResult:
        """Compute wavefront tracking.

        New API: pass source and receivers as arguments.
        Legacy API: call with no positional arguments (uses stored sources/receivers).

        Args:
            source: Sources object or np.ndarray of shape (n, 2). Optional for legacy API.
            receivers: Receivers object or np.ndarray of shape (n, 2). Optional for legacy API.
            pool: Optional pool executor for parallel execution.
                Supports ProcessPoolExecutor and schwimmbad pools.
            extract_raypaths: Whether to extract raypath data. If None, uses options.
            extract_wavefronts: Whether to extract wavefront data. If None, uses options.

        Returns:
            WaveFrontTrackerResult with arrivals, raypaths, wavefronts, and Jacobian.
        """
        # Handle source/receiver arguments
        if source is not None and receivers is not None:
            self._set_sources_and_receivers(source, receivers)
        elif source is not None or receivers is not None:
            raise ValueError("Must provide both source and receivers, or neither.")
        elif self._sources is None and self.config_file is None:
            raise ValueError(
                "No sources/receivers configured. Either:\n"
                "  - Pass source and receivers to forward(), or\n"
                "  - Provide them at construction time"
            )

        # Resolve extraction flags
        do_raypaths, do_wavefronts, do_frechet = self._resolve_extraction_flags(
            extract_raypaths, extract_wavefronts
        )

        # Dispatch: pool or sequential
        if pool is not None:
            self._require_file_free("pool-based parallel execution")
            return self._forward_parallel(pool, do_raypaths, do_wavefronts, do_frechet)
        else:
            if self._file_free:
                return self._forward_sequential(do_raypaths, do_wavefronts, do_frechet)
            else:
                return self._forward_sequential_file_based(do_raypaths, do_wavefronts, do_frechet)

    def _forward_sequential(self, extract_raypaths, extract_wavefronts, extract_frechet):
        """Sequential execution using in-memory worker (file-free mode)."""
        model_dict, sources_dict, receivers_dict = self._serialize_data()

        results = []
        for source_id in range(1, self.n_sources + 1):
            result = _calc_single_source_worker(
                model_dict, sources_dict, receivers_dict, self._options_dict,
                source_id, extract_raypaths, extract_wavefronts, extract_frechet,
            )
            results.append(result)

        return reduce(operator.add, results)

    def _forward_sequential_file_based(self, extract_raypaths, extract_wavefronts, extract_frechet):
        """Sequential execution using the library instance directly (deprecated config_file mode)."""
        # Disable file output
        self._lib.set_file_output(False)

        # Ensure flags
        if extract_raypaths or extract_frechet:
            self._lib.set_do_rays(True)
        if extract_frechet:
            self._lib.set_do_frechet(True)

        results = []
        for source_id in range(1, self.n_sources + 1):
            self._lib.run_single_source(source_id)

            source_ids, receiver_ids, arrival_nums, times, azis, spfs = self._lib.get_arrivals(source_id)

            raypaths = []
            if extract_raypaths:
                raypaths = self._lib.get_raypaths(source_id)

            wavefronts = {}
            if extract_wavefronts:
                wf_list = self._lib.get_wavefronts(source_id)
                if wf_list:
                    wavefronts[source_id] = wf_list

            jacobian_rows = jacobian_cols = jacobian_vals = jacobian_shape = None
            if extract_frechet:
                jac_data = self._lib.get_jacobian()
                if jac_data is not None:
                    jacobian_rows, jacobian_cols, jacobian_vals, jacobian_shape = jac_data

            results.append(WaveFrontTrackerResult(
                source_ids=source_ids,
                receiver_ids=receiver_ids,
                arrival_numbers=arrival_nums,
                travel_times=times,
                azimuths=azis,
                spreading_factors=spfs,
                raypaths=raypaths,
                wavefronts=wavefronts,
                jacobian_rows=jacobian_rows,
                jacobian_cols=jacobian_cols,
                jacobian_vals=jacobian_vals,
                jacobian_shape=jacobian_shape,
            ))

        return reduce(operator.add, results)

    def _forward_parallel(self, pool, extract_raypaths, extract_wavefronts, extract_frechet):
        """Parallel execution using a user-provided pool executor."""
        model_dict, sources_dict, receivers_dict = self._serialize_data()

        worker_args = [
            (model_dict, sources_dict, receivers_dict, self._options_dict,
             sid, extract_raypaths, extract_wavefronts, extract_frechet)
            for sid in range(1, self.n_sources + 1)
        ]

        pool_type = self._detect_pool_type(pool)

        if pool_type == 'process':
            # ProcessPoolExecutor: use submit()
            futures = [pool.submit(_calc_single_source_worker, *a) for a in worker_args]
            results = [f.result() for f in futures]
        else:
            # schwimmbad / map-style pool
            results = list(pool.map(_worker_wrapper, worker_args))

        return reduce(operator.add, results)

    def forward_single_source(
        self,
        source_id: int,
        extract_raypaths: bool = None,
        extract_wavefronts: bool = None,
    ) -> WaveFrontTrackerResult:
        """Compute wavefront tracking for a single source.

        Args:
            source_id: Source ID to process (1-based indexing)
            extract_raypaths: Whether to extract raypath data. If None, uses options.
            extract_wavefronts: Whether to extract wavefront data. If None, uses options.

        Returns:
            WaveFrontTrackerResult for this source
        """
        if source_id < 1 or source_id > self.n_sources:
            raise ValueError(
                f"Source ID {source_id} out of range [1, {self.n_sources}]"
            )

        do_raypaths, do_wavefronts, do_frechet = self._resolve_extraction_flags(
            extract_raypaths, extract_wavefronts
        )

        if self._file_free:
            model_dict, sources_dict, receivers_dict = self._serialize_data()
            return _calc_single_source_worker(
                model_dict, sources_dict, receivers_dict, self._options_dict,
                source_id, do_raypaths, do_wavefronts, do_frechet,
            )
        else:
            # Deprecated file-based path: use library directly
            self._lib.set_file_output(False)
            if do_raypaths or do_frechet:
                self._lib.set_do_rays(True)
            if do_frechet:
                self._lib.set_do_frechet(True)

            self._lib.run_single_source(source_id)
            source_ids, receiver_ids, arrival_nums, times, azis, spfs = self._lib.get_arrivals(source_id)

            raypaths = []
            if do_raypaths:
                raypaths = self._lib.get_raypaths(source_id)

            wavefronts = {}
            if do_wavefronts:
                wf_list = self._lib.get_wavefronts(source_id)
                if wf_list:
                    wavefronts[source_id] = wf_list

            jacobian_rows = jacobian_cols = jacobian_vals = jacobian_shape = None
            if do_frechet:
                jac_data = self._lib.get_jacobian()
                if jac_data is not None:
                    jacobian_rows, jacobian_cols, jacobian_vals, jacobian_shape = jac_data

            return WaveFrontTrackerResult(
                source_ids=source_ids,
                receiver_ids=receiver_ids,
                arrival_numbers=arrival_nums,
                travel_times=times,
                azimuths=azis,
                spreading_factors=spfs,
                raypaths=raypaths,
                wavefronts=wavefronts,
                jacobian_rows=jacobian_rows,
                jacobian_cols=jacobian_cols,
                jacobian_vals=jacobian_vals,
                jacobian_shape=jacobian_shape,
            )

    def __repr__(self) -> str:
        """String representation."""
        if self.config_file:
            return f"WaveFrontTracker(config='{self.config_file}', n_sources={self.n_sources})"
        return f"WaveFrontTracker(model={self._model}, n_sources={self.n_sources})"
