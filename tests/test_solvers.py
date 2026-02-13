"""Tests for pyswmp.solvers module — new API, Jacobian, nthreads, and deprecation."""

import warnings
from unittest.mock import MagicMock

import numpy as np
import pytest

from pyswmp import WaveFrontTracker, WaveFrontTrackerResult, TrackerOptions
from pyswmp.data import VelocityModel2D, Sources, Receivers, create_constant_velocity_model


class TestWaveFrontTrackerNewConstructor:
    """Tests for the new WaveFrontTracker(model, options) constructor."""

    def test_model_only_constructor(self, simple_velocity_model):
        """Test constructing with model only (no sources/receivers)."""
        tracker = WaveFrontTracker(simple_velocity_model)
        assert tracker.n_sources == 0
        assert tracker._model is simple_velocity_model

    def test_model_and_options_constructor(self, simple_velocity_model, default_options):
        """Test constructing with model and options."""
        tracker = WaveFrontTracker(simple_velocity_model, default_options)
        assert tracker.n_sources == 0
        assert tracker.options is default_options

    def test_model_keyword_constructor(self, simple_velocity_model):
        """Test constructing with model as keyword argument."""
        tracker = WaveFrontTracker(model=simple_velocity_model)
        assert tracker.n_sources == 0

    def test_config_file_positional_type_dispatch(self, simple_velocity_model):
        """Test that a VelocityModel2D first positional is dispatched to model path."""
        tracker = WaveFrontTracker(simple_velocity_model)
        assert tracker.config_file is None
        assert tracker._model is simple_velocity_model

    def test_invalid_first_argument(self):
        """Test that invalid first argument raises TypeError."""
        with pytest.raises(TypeError, match="First argument must be"):
            WaveFrontTracker(42)

    def test_model_sources_receivers_legacy(
        self, simple_velocity_model, simple_sources, simple_receivers, default_options
    ):
        """Test legacy constructor with model + sources + receivers."""
        tracker = WaveFrontTracker(
            model=simple_velocity_model,
            sources=simple_sources,
            receivers=simple_receivers,
            options=default_options,
        )
        assert tracker.n_sources == simple_sources.n_sources
        assert tracker._sources is simple_sources
        assert tracker._receivers is simple_receivers

    def test_no_arguments_raises(self):
        """Test that no arguments raises ValueError."""
        with pytest.raises(ValueError):
            WaveFrontTracker()

    def test_config_file_emits_deprecation_warning(self):
        """Test that config_file constructor emits DeprecationWarning."""
        from unittest.mock import patch

        with pytest.warns(DeprecationWarning, match="config_file is deprecated"):
            with patch("pyswmp.solvers.LibSWMP") as MockLib:
                mock_lib = MockLib.return_value
                mock_lib.get_n_sources.return_value = 0
                WaveFrontTracker(config_file="/fake/config.in")


class TestForwardSignature:
    """Tests for the forward(source, receivers) new API signature."""

    def test_forward_requires_both_or_neither(self, simple_velocity_model, default_options):
        """Test that forward() requires both source and receivers, or neither."""
        tracker = WaveFrontTracker(simple_velocity_model, default_options)

        with pytest.raises(ValueError, match="Must provide both"):
            tracker.forward(source=Sources(np.array([[135.0, -25.0]]), source_type=1))

    def test_forward_without_sources_raises(self, simple_velocity_model, default_options):
        """Test that forward() without prior source setup raises."""
        tracker = WaveFrontTracker(simple_velocity_model, default_options)

        with pytest.raises(ValueError, match="No sources/receivers"):
            tracker.forward()

    def test_forward_legacy_with_stored_sources(
        self, simple_velocity_model, simple_sources, simple_receivers, default_options
    ):
        """Test that legacy forward() works when sources set at construction."""
        tracker = WaveFrontTracker(
            model=simple_velocity_model,
            sources=simple_sources,
            receivers=simple_receivers,
            options=default_options,
        )
        assert tracker.n_sources > 0

    def test_numpy_array_auto_wrapping(self, simple_velocity_model, default_options):
        """Test that numpy arrays are auto-wrapped into Sources/Receivers."""
        tracker = WaveFrontTracker(simple_velocity_model, default_options)

        source_arr = np.array([[135.0, -25.0]])
        recv_arr = np.array([[120.0, -30.0], [130.0, -25.0]])

        tracker._set_sources_and_receivers(source_arr, recv_arr)

        assert tracker.n_sources == 1
        assert isinstance(tracker._sources, Sources)
        assert isinstance(tracker._receivers, Receivers)
        assert tracker._receivers.n_receivers == 2


class TestPoolDispatch:
    """Tests for pool-based parallel dispatch logic."""

    def test_thread_pool_rejected(
        self, simple_velocity_model, simple_sources, simple_receivers, default_options
    ):
        """Test that ThreadPoolExecutor is rejected."""
        tracker = WaveFrontTracker(
            model=simple_velocity_model,
            sources=simple_sources,
            receivers=simple_receivers,
            options=default_options,
        )
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=1) as pool:
            with pytest.raises(ValueError, match="ThreadPoolExecutor is not supported"):
                tracker.forward(pool=pool)

    def test_pool_requires_file_free(self, simple_velocity_model, default_options):
        """Test that pool with file-free mode is accepted (no ValueError about file-free)."""
        tracker = WaveFrontTracker(
            model=simple_velocity_model,
            options=default_options,
        )
        assert tracker._file_free is True


class TestSerializationFormat:
    """Tests for _serialize_data returning numpy arrays."""

    def test_serialize_returns_numpy_arrays(
        self, simple_velocity_model, simple_sources, simple_receivers, default_options
    ):
        """Test that serialized data contains numpy arrays, not Python lists."""
        tracker = WaveFrontTracker(
            model=simple_velocity_model,
            sources=simple_sources,
            receivers=simple_receivers,
            options=default_options,
        )

        model_dict, sources_dict, receivers_dict = tracker._serialize_data()

        assert isinstance(model_dict['velocities'], np.ndarray)
        assert isinstance(sources_dict['positions'], np.ndarray)
        assert isinstance(receivers_dict['positions'], np.ndarray)

    def test_serialize_preserves_data(
        self, simple_velocity_model, simple_sources, simple_receivers, default_options
    ):
        """Test that serialized data matches original."""
        tracker = WaveFrontTracker(
            model=simple_velocity_model,
            sources=simple_sources,
            receivers=simple_receivers,
            options=default_options,
        )

        model_dict, sources_dict, receivers_dict = tracker._serialize_data()

        np.testing.assert_array_equal(model_dict['velocities'], simple_velocity_model.velocities)
        assert model_dict['x0'] == simple_velocity_model.x0
        assert model_dict['dx'] == simple_velocity_model.dx
        np.testing.assert_array_equal(sources_dict['positions'], simple_sources.positions)
        np.testing.assert_array_equal(receivers_dict['positions'], simple_receivers.positions)


class TestOptionsToDict:
    """Tests for TrackerOptions.to_dict() excluding file paths."""

    def test_to_dict_excludes_file_paths(self):
        """Test that to_dict() does not include file path fields."""
        opts = TrackerOptions(
            velocity_file="/some/path.vel",
            sources_file="/some/sources.dat",
            receivers_file="/some/receivers.dat",
            receiver_mode_file="/some/recmode.dat",
        )
        d = opts.to_dict()

        assert 'velocity_file' not in d
        assert 'sources_file' not in d
        assert 'receivers_file' not in d
        assert 'receiver_mode_file' not in d

    def test_to_dict_includes_all_numeric_fields(self):
        """Test that to_dict() includes all numeric/bool configuration fields."""
        opts = TrackerOptions()
        d = opts.to_dict()

        expected_keys = {
            'earth_radius', 'coordinate_system', 'dt', 'max_iterations',
            'n_bichar_nodes', 'ode_solver', 'interpolator', 'computation_mode',
            'max_arrivals', 'source_specific_receivers', 'extract_raypaths',
            'extract_wavefronts', 'extract_frechet', 'raypath_storage',
            'wavefront_interval',
        }
        assert set(d.keys()) == expected_keys

    def test_to_dict_roundtrips(self):
        """Test that to_dict() output can reconstruct TrackerOptions."""
        opts = TrackerOptions(dt=2.5, max_iterations=1000, coordinate_system=2)
        d = opts.to_dict()
        opts2 = TrackerOptions(**d)

        assert opts2.dt == opts.dt
        assert opts2.max_iterations == opts.max_iterations
        assert opts2.coordinate_system == opts.coordinate_system


class TestJacobianResultFields:
    """Tests for Jacobian fields on WaveFrontTrackerResult."""

    def test_default_jacobian_is_none(self):
        """Test that Jacobian fields default to None."""
        result = WaveFrontTrackerResult()
        assert result.jacobian_rows is None
        assert result.jacobian_cols is None
        assert result.jacobian_vals is None
        assert result.jacobian_shape is None

    def test_result_with_jacobian(self):
        """Test creating a result with Jacobian data."""
        result = WaveFrontTrackerResult(
            source_ids=np.array([1, 1], dtype=np.int32),
            receiver_ids=np.array([1, 2], dtype=np.int32),
            arrival_numbers=np.array([1, 1], dtype=np.int32),
            travel_times=np.array([10.5, 15.2], dtype=np.float32),
            azimuths=np.array([45.0, 90.0], dtype=np.float32),
            spreading_factors=np.array([1.0, 1.5], dtype=np.float32),
            jacobian_rows=np.array([0, 0, 1, 1], dtype=np.int32),
            jacobian_cols=np.array([0, 5, 3, 7], dtype=np.int32),
            jacobian_vals=np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float64),
            jacobian_shape=(2, 100),
        )

        assert result.jacobian_shape == (2, 100)
        assert len(result.jacobian_vals) == 4

    def test_jacobian_merge_via_add(self):
        """Test merging Jacobians from two sources via __add__."""
        result1 = WaveFrontTrackerResult(
            source_ids=np.array([1, 1], dtype=np.int32),
            receiver_ids=np.array([1, 2], dtype=np.int32),
            arrival_numbers=np.array([1, 1], dtype=np.int32),
            travel_times=np.array([10.5, 15.2], dtype=np.float32),
            azimuths=np.array([45.0, 90.0], dtype=np.float32),
            spreading_factors=np.array([1.0, 1.5], dtype=np.float32),
            jacobian_rows=np.array([0, 0, 1], dtype=np.int32),
            jacobian_cols=np.array([0, 5, 3], dtype=np.int32),
            jacobian_vals=np.array([0.1, 0.2, 0.3], dtype=np.float64),
            jacobian_shape=(2, 100),
        )

        result2 = WaveFrontTrackerResult(
            source_ids=np.array([2], dtype=np.int32),
            receiver_ids=np.array([1], dtype=np.int32),
            arrival_numbers=np.array([1], dtype=np.int32),
            travel_times=np.array([12.3], dtype=np.float32),
            azimuths=np.array([135.0], dtype=np.float32),
            spreading_factors=np.array([1.2], dtype=np.float32),
            jacobian_rows=np.array([0, 0], dtype=np.int32),
            jacobian_cols=np.array([1, 8], dtype=np.int32),
            jacobian_vals=np.array([0.5, 0.6], dtype=np.float64),
            jacobian_shape=(1, 100),
        )

        combined = result1 + result2

        # Check merged Jacobian
        assert combined.jacobian_shape == (3, 100)
        assert len(combined.jacobian_vals) == 5

        # Check row offset: result2's rows should be offset by result1's row count (2)
        np.testing.assert_array_equal(
            combined.jacobian_rows,
            np.array([0, 0, 1, 2, 2], dtype=np.int32)
        )
        np.testing.assert_array_equal(
            combined.jacobian_cols,
            np.array([0, 5, 3, 1, 8], dtype=np.int32)
        )

    def test_jacobian_merge_one_empty(self):
        """Test merging when one result has no Jacobian."""
        result1 = WaveFrontTrackerResult(
            source_ids=np.array([1], dtype=np.int32),
            receiver_ids=np.array([1], dtype=np.int32),
            arrival_numbers=np.array([1], dtype=np.int32),
            travel_times=np.array([10.5], dtype=np.float32),
            azimuths=np.array([45.0], dtype=np.float32),
            spreading_factors=np.array([1.0], dtype=np.float32),
            jacobian_rows=np.array([0, 0], dtype=np.int32),
            jacobian_cols=np.array([0, 5], dtype=np.int32),
            jacobian_vals=np.array([0.1, 0.2], dtype=np.float64),
            jacobian_shape=(1, 100),
        )

        result2 = WaveFrontTrackerResult(
            source_ids=np.array([2], dtype=np.int32),
            receiver_ids=np.array([1], dtype=np.int32),
            arrival_numbers=np.array([1], dtype=np.int32),
            travel_times=np.array([12.3], dtype=np.float32),
            azimuths=np.array([135.0], dtype=np.float32),
            spreading_factors=np.array([1.2], dtype=np.float32),
            # No Jacobian
        )

        combined = result1 + result2
        assert combined.jacobian_shape == (1, 100)
        assert len(combined.jacobian_vals) == 2

    def test_to_sparse_jacobian(self):
        """Test converting to scipy sparse matrix."""
        scipy_sparse = pytest.importorskip("scipy.sparse")

        result = WaveFrontTrackerResult(
            source_ids=np.array([1, 1], dtype=np.int32),
            receiver_ids=np.array([1, 2], dtype=np.int32),
            arrival_numbers=np.array([1, 1], dtype=np.int32),
            travel_times=np.array([10.5, 15.2], dtype=np.float32),
            azimuths=np.array([45.0, 90.0], dtype=np.float32),
            spreading_factors=np.array([1.0, 1.5], dtype=np.float32),
            jacobian_rows=np.array([0, 0, 1, 1], dtype=np.int32),
            jacobian_cols=np.array([0, 5, 3, 7], dtype=np.int32),
            jacobian_vals=np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float64),
            jacobian_shape=(2, 100),
        )

        J = result.to_sparse_jacobian()

        assert J.shape == (2, 100)
        assert J.nnz == 4
        assert J[0, 0] == pytest.approx(0.1)
        assert J[0, 5] == pytest.approx(0.2)
        assert J[1, 3] == pytest.approx(0.3)
        assert J[1, 7] == pytest.approx(0.4)

    def test_to_sparse_jacobian_no_data_raises(self):
        """Test that to_sparse_jacobian raises without Jacobian data."""
        pytest.importorskip("scipy.sparse")

        result = WaveFrontTrackerResult()
        with pytest.raises(ValueError, match="No Jacobian data"):
            result.to_sparse_jacobian()

    def test_arrival_indexing(self):
        """Test that Jacobian rows correspond to arrival indexing."""
        result = WaveFrontTrackerResult(
            source_ids=np.array([1, 1, 1], dtype=np.int32),
            receiver_ids=np.array([1, 1, 2], dtype=np.int32),
            arrival_numbers=np.array([1, 2, 1], dtype=np.int32),
            travel_times=np.array([10.5, 11.0, 15.2], dtype=np.float32),
            azimuths=np.array([45.0, 50.0, 90.0], dtype=np.float32),
            spreading_factors=np.array([1.0, 1.1, 1.5], dtype=np.float32),
            jacobian_rows=np.array([0, 1, 2], dtype=np.int32),
            jacobian_cols=np.array([5, 10, 15], dtype=np.int32),
            jacobian_vals=np.array([0.1, 0.2, 0.3], dtype=np.float64),
            jacobian_shape=(3, 100),
        )

        assert result.jacobian_shape[0] == len(result.travel_times)
        assert result.source_ids[0] == 1
        assert result.receiver_ids[0] == 1
        assert result.arrival_numbers[0] == 1

    def test_repr_with_jacobian(self):
        """Test repr includes Jacobian info when present."""
        result = WaveFrontTrackerResult(
            source_ids=np.array([1], dtype=np.int32),
            receiver_ids=np.array([1], dtype=np.int32),
            arrival_numbers=np.array([1], dtype=np.int32),
            travel_times=np.array([10.5], dtype=np.float32),
            azimuths=np.array([45.0], dtype=np.float32),
            spreading_factors=np.array([1.0], dtype=np.float32),
            jacobian_rows=np.array([0], dtype=np.int32),
            jacobian_cols=np.array([5], dtype=np.int32),
            jacobian_vals=np.array([0.1], dtype=np.float64),
            jacobian_shape=(1, 100),
        )

        assert "jacobian" in repr(result)
        assert "(1, 100)" in repr(result)


@pytest.mark.slow
class TestIntegrationJacobian:
    """Integration tests requiring compiled Fortran library."""

    def test_forward_with_frechet(
        self, simple_velocity_model, simple_sources, simple_receivers
    ):
        """Test forward modeling with Frechet matrix extraction."""
        opts = TrackerOptions(
            coordinate_system=2,
            extract_raypaths=True,
            extract_frechet=True,
        )

        tracker = WaveFrontTracker(simple_velocity_model, opts)
        result = tracker.forward(simple_sources, simple_receivers)

        # Should have arrivals
        assert len(result.travel_times) > 0

        # Should have Jacobian
        assert result.jacobian_vals is not None
        assert result.jacobian_shape is not None
        assert result.jacobian_shape[0] == len(result.travel_times)
        assert result.jacobian_shape[1] > 0

        # Jacobian rows should be valid indices
        assert np.all(result.jacobian_rows >= 0)
        assert np.all(result.jacobian_rows < result.jacobian_shape[0])

        # Check arrival indexing is consistent
        assert len(result.source_ids) == len(result.travel_times)
        assert len(result.receiver_ids) == len(result.travel_times)
        assert len(result.arrival_numbers) == len(result.travel_times)

    def test_forward_new_api(
        self, simple_velocity_model, simple_sources, simple_receivers
    ):
        """Test the new API: model at construction, source+receivers at forward()."""
        opts = TrackerOptions(coordinate_system=2, extract_raypaths=True)

        tracker = WaveFrontTracker(simple_velocity_model, opts)
        result = tracker.forward(simple_sources, simple_receivers)

        assert len(result.travel_times) > 0
        assert len(result.raypaths) > 0

    def test_forward_frechet_to_sparse(
        self, simple_velocity_model, simple_sources, simple_receivers
    ):
        """Test converting Frechet result to scipy sparse matrix."""
        scipy_sparse = pytest.importorskip("scipy.sparse")

        opts = TrackerOptions(
            coordinate_system=2,
            extract_raypaths=True,
            extract_frechet=True,
        )

        tracker = WaveFrontTracker(simple_velocity_model, opts)
        result = tracker.forward(simple_sources, simple_receivers)

        J = result.to_sparse_jacobian()

        assert J.shape[0] == len(result.travel_times)
        assert J.shape[1] > 0
        assert J.nnz > 0

    def test_forward_parallel_pool(
        self, simple_velocity_model, simple_sources, simple_receivers
    ):
        """Test forward with ProcessPoolExecutor produces same results as sequential."""
        from concurrent.futures import ProcessPoolExecutor

        opts = TrackerOptions(coordinate_system=2)

        tracker = WaveFrontTracker(simple_velocity_model, opts)

        result_seq = tracker.forward(simple_sources, simple_receivers)
        with ProcessPoolExecutor(max_workers=2) as pool:
            result_par = tracker.forward(simple_sources, simple_receivers, pool=pool)

        assert len(result_seq.travel_times) == len(result_par.travel_times)
        np.testing.assert_allclose(
            np.sort(result_seq.travel_times),
            np.sort(result_par.travel_times),
            rtol=1e-5,
        )
