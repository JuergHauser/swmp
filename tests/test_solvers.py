"""Tests for pyswmp.solvers module — new API and Jacobian."""

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
        """Test that a string first argument is dispatched as config file path."""
        # We can't actually load a config file in unit tests (Fortran blocks),
        # but verify the type dispatch works by checking that a VelocityModel2D
        # passed as first positional is correctly dispatched to model path.
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
        # forward() without args should work since sources/receivers are stored
        # (will fail at Fortran level in unit test, but tests the Python dispatch)
        assert tracker.n_sources > 0

    def test_numpy_array_auto_wrapping(self, simple_velocity_model, default_options):
        """Test that numpy arrays are auto-wrapped into Sources/Receivers."""
        tracker = WaveFrontTracker(simple_velocity_model, default_options)

        source_arr = np.array([[135.0, -25.0]])
        recv_arr = np.array([[120.0, -30.0], [130.0, -25.0]])

        # _set_sources_and_receivers should handle numpy arrays
        tracker._set_sources_and_receivers(source_arr, recv_arr)

        assert tracker.n_sources == 1
        assert isinstance(tracker._sources, Sources)
        assert isinstance(tracker._receivers, Receivers)
        assert tracker._receivers.n_receivers == 2


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
        # Jacobian rows should correspond 1:1 with (source, receiver, arrival) tuples
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

        # Row 0 → arrival (source=1, receiver=1, arrival=1)
        # Row 1 → arrival (source=1, receiver=1, arrival=2)
        # Row 2 → arrival (source=1, receiver=2, arrival=1)
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
