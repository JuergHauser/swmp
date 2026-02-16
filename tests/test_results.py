"""Unit tests for pyswmp.results module."""

import numpy as np
import pytest
from pyswmp import WaveFrontTrackerResult


class TestWaveFrontTrackerResult:
    """Tests for WaveFrontTrackerResult dataclass."""

    def test_empty_result(self):
        """Test creating an empty result."""
        result = WaveFrontTrackerResult()

        assert len(result.source_ids) == 0
        assert len(result.receiver_ids) == 0
        assert len(result.travel_times) == 0
        assert len(result.raypaths) == 0
        assert len(result.wavefronts) == 0

    def test_result_with_data(self):
        """Test creating a result with data."""
        result = WaveFrontTrackerResult(
            source_ids=np.array([1, 1, 2]),
            receiver_ids=np.array([1, 2, 1]),
            arrival_numbers=np.array([1, 1, 1]),
            travel_times=np.array([10.5, 15.2, 12.3]),
            azimuths=np.array([45.0, 90.0, 135.0]),
            spreading_factors=np.array([1.0, 1.5, 1.2])
        )

        assert len(result.source_ids) == 3
        assert len(result.receiver_ids) == 3
        assert len(result.travel_times) == 3
        assert result.travel_times[0] == 10.5

    def test_result_addition(self):
        """Test adding two results together."""
        result1 = WaveFrontTrackerResult(
            source_ids=np.array([1, 1]),
            receiver_ids=np.array([1, 2]),
            arrival_numbers=np.array([1, 1]),
            travel_times=np.array([10.5, 15.2]),
            azimuths=np.array([45.0, 90.0]),
            spreading_factors=np.array([1.0, 1.5])
        )

        result2 = WaveFrontTrackerResult(
            source_ids=np.array([2]),
            receiver_ids=np.array([1]),
            arrival_numbers=np.array([1]),
            travel_times=np.array([12.3]),
            azimuths=np.array([135.0]),
            spreading_factors=np.array([1.2])
        )

        combined = result1 + result2

        assert len(combined.source_ids) == 3
        assert len(combined.travel_times) == 3
        assert combined.travel_times[0] == 10.5
        assert combined.travel_times[2] == 12.3

    def test_result_addition_with_raypaths(self):
        """Test adding results with raypaths."""
        raypath1 = np.array([[0, 0], [1, 1], [2, 2]])
        raypath2 = np.array([[0, 0], [1.5, 1.5], [3, 3]])

        result1 = WaveFrontTrackerResult(
            source_ids=np.array([1]),
            receiver_ids=np.array([1]),
            arrival_numbers=np.array([1]),
            travel_times=np.array([10.5]),
            azimuths=np.array([45.0]),
            spreading_factors=np.array([1.0]),
            raypaths=[raypath1]
        )

        result2 = WaveFrontTrackerResult(
            source_ids=np.array([2]),
            receiver_ids=np.array([1]),
            arrival_numbers=np.array([1]),
            travel_times=np.array([12.3]),
            azimuths=np.array([135.0]),
            spreading_factors=np.array([1.2]),
            raypaths=[raypath2]
        )

        combined = result1 + result2

        assert len(combined.raypaths) == 2
        assert np.array_equal(combined.raypaths[0], raypath1)
        assert np.array_equal(combined.raypaths[1], raypath2)

    def test_to_dataframe(self):
        """Test converting result to pandas DataFrame."""
        pytest.importorskip("pandas")  # Skip test if pandas not installed

        result = WaveFrontTrackerResult(
            source_ids=np.array([1, 1, 2]),
            receiver_ids=np.array([1, 2, 1]),
            arrival_numbers=np.array([1, 1, 1]),
            travel_times=np.array([10.5, 15.2, 12.3]),
            azimuths=np.array([45.0, 90.0, 135.0]),
            spreading_factors=np.array([1.0, 1.5, 1.2])
        )

        df = result.to_dataframe()

        assert len(df) == 3
        assert 'source_id' in df.columns
        assert 'receiver_id' in df.columns
        assert 'travel_time' in df.columns
        assert df['travel_time'].iloc[0] == 10.5

    def test_to_dataframe_without_pandas(self):
        """Test that to_dataframe raises ImportError without pandas."""
        # Save pandas import and temporarily remove it
        import sys
        pandas_module = sys.modules.get('pandas')
        if pandas_module:
            sys.modules['pandas'] = None

        try:
            result = WaveFrontTrackerResult(
                source_ids=np.array([1]),
                receiver_ids=np.array([1]),
                arrival_numbers=np.array([1]),
                travel_times=np.array([10.5]),
                azimuths=np.array([45.0]),
                spreading_factors=np.array([1.0])
            )

            with pytest.raises(ImportError):
                result.to_dataframe()
        finally:
            # Restore pandas
            if pandas_module:
                sys.modules['pandas'] = pandas_module

    def test_repr(self):
        """Test string representation."""
        result = WaveFrontTrackerResult(
            source_ids=np.array([1, 1, 2]),
            receiver_ids=np.array([1, 2, 1]),
            arrival_numbers=np.array([1, 1, 1]),
            travel_times=np.array([10.5, 15.2, 12.3]),
            azimuths=np.array([45.0, 90.0, 135.0]),
            spreading_factors=np.array([1.0, 1.5, 1.2])
        )

        repr_str = repr(result)

        assert "WaveFrontTrackerResult" in repr_str
        assert "3 arrivals" in repr_str or "arrivals=3" in repr_str
