"""Result data structures for SWMP wavefront tracking.

This module provides result classes that support aggregation via __add__,
following the pyfm2d pattern for parallel execution.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple
import numpy as np
import numpy.typing as npt


@dataclass
class WaveFrontTrackerResult:
    """Results from wavefront tracking (single or multiple sources).

    This class stores arrival times, raypaths, and wavefronts from
    wavefront tracking computations. Data is stored in flattened arrays
    to naturally handle variable arrival counts per receiver.

    Attributes:
        source_ids: Source IDs for each arrival (1-based)
        receiver_ids: Receiver IDs for each arrival (1-based)
        arrival_numbers: Arrival numbers (1=first, 2=second, etc.)
        travel_times: Arrival times
        azimuths: Azimuths at receivers
        spreading_factors: Spreading factors (if mode=2)
        raypaths: List of raypath dicts with keys: source, receiver, arrival, path
        wavefronts: Dict mapping source_id to list of wavefront arrays

    Example:
        >>> result = WaveFrontTrackerResult(
        ...     source_ids=np.array([1, 1, 2]),
        ...     receiver_ids=np.array([1, 2, 1]),
        ...     arrival_numbers=np.array([1, 1, 1]),
        ...     travel_times=np.array([10.5, 12.3, 11.1]),
        ...     azimuths=np.array([45.0, 30.0, 50.0]),
        ...     spreading_factors=np.array([1.0, 0.9, 1.1])
        ... )
        >>> df = result.to_dataframe()
    """

    # Arrival data (flattened arrays)
    source_ids: npt.NDArray[np.int32] = field(
        default_factory=lambda: np.array([], dtype=np.int32)
    )
    receiver_ids: npt.NDArray[np.int32] = field(
        default_factory=lambda: np.array([], dtype=np.int32)
    )
    arrival_numbers: npt.NDArray[np.int32] = field(
        default_factory=lambda: np.array([], dtype=np.int32)
    )
    travel_times: npt.NDArray[np.float32] = field(
        default_factory=lambda: np.array([], dtype=np.float32)
    )
    azimuths: npt.NDArray[np.float32] = field(
        default_factory=lambda: np.array([], dtype=np.float32)
    )
    spreading_factors: npt.NDArray[np.float32] = field(
        default_factory=lambda: np.array([], dtype=np.float32)
    )

    # Raypaths (list of dicts)
    raypaths: List[Dict] = field(default_factory=list)
    # Each dict has keys: source, receiver, arrival, path (ndarray of shape (n, 2))

    # Wavefronts (dict by source)
    wavefronts: Dict[int, List[npt.NDArray]] = field(default_factory=dict)

    # Jacobian (Frechet matrix) in COO format
    # Row indices correspond 1:1 with arrivals (source_ids, receiver_ids, arrival_numbers)
    jacobian_rows: Optional[npt.NDArray[np.int32]] = None
    jacobian_cols: Optional[npt.NDArray[np.int32]] = None
    jacobian_vals: Optional[npt.NDArray[np.float64]] = None
    jacobian_shape: Optional[Tuple[int, int]] = None

    def __add__(self, other: 'WaveFrontTrackerResult') -> 'WaveFrontTrackerResult':
        """Merge two result objects (pyfm2d pattern).

        This enables reduction aggregation:
            combined = reduce(operator.add, results)

        Args:
            other: Another WaveFrontTrackerResult to merge

        Returns:
            New WaveFrontTrackerResult with concatenated data
        """
        if other is None:
            return self

        # Concatenate arrival arrays
        merged = WaveFrontTrackerResult(
            source_ids=np.concatenate([self.source_ids, other.source_ids]),
            receiver_ids=np.concatenate([self.receiver_ids, other.receiver_ids]),
            arrival_numbers=np.concatenate([self.arrival_numbers, other.arrival_numbers]),
            travel_times=np.concatenate([self.travel_times, other.travel_times]),
            azimuths=np.concatenate([self.azimuths, other.azimuths]),
            spreading_factors=np.concatenate([self.spreading_factors, other.spreading_factors]),
        )

        # Merge raypaths
        merged.raypaths = (self.raypaths or []) + (other.raypaths or [])

        # Merge wavefronts
        merged.wavefronts = {**(self.wavefronts or {}), **(other.wavefronts or {})}

        # Merge Jacobians (COO format — offset row indices for stacking)
        has_self_jac = self.jacobian_vals is not None and len(self.jacobian_vals) > 0
        has_other_jac = other.jacobian_vals is not None and len(other.jacobian_vals) > 0

        if has_self_jac and has_other_jac:
            # Offset other's row indices by self's row count
            self_nrows = self.jacobian_shape[0]
            ncols = max(self.jacobian_shape[1], other.jacobian_shape[1])
            total_rows = self_nrows + other.jacobian_shape[0]

            merged.jacobian_rows = np.concatenate([
                self.jacobian_rows,
                other.jacobian_rows + self_nrows,
            ])
            merged.jacobian_cols = np.concatenate([self.jacobian_cols, other.jacobian_cols])
            merged.jacobian_vals = np.concatenate([self.jacobian_vals, other.jacobian_vals])
            merged.jacobian_shape = (total_rows, ncols)
        elif has_self_jac:
            merged.jacobian_rows = self.jacobian_rows
            merged.jacobian_cols = self.jacobian_cols
            merged.jacobian_vals = self.jacobian_vals
            merged.jacobian_shape = self.jacobian_shape
        elif has_other_jac:
            merged.jacobian_rows = other.jacobian_rows
            merged.jacobian_cols = other.jacobian_cols
            merged.jacobian_vals = other.jacobian_vals
            merged.jacobian_shape = other.jacobian_shape

        return merged

    def to_dataframe(self):
        """Convert arrival data to pandas DataFrame.

        Returns:
            DataFrame with columns: source, receiver, arrival, time, azimuth, spreading

        Raises:
            ImportError: If pandas is not installed
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("pandas is required for to_dataframe(). Install with: pip install pandas")

        return pd.DataFrame({
            'source': self.source_ids,
            'receiver': self.receiver_ids,
            'arrival': self.arrival_numbers,
            'time': self.travel_times,
            'azimuth': self.azimuths,
            'spreading': self.spreading_factors,
        })

    def to_sparse_jacobian(self):
        """Convert Jacobian to scipy sparse CSR matrix.

        Returns:
            scipy.sparse.csr_matrix of shape (n_arrivals, n_model_params)

        Raises:
            ImportError: If scipy is not installed
            ValueError: If no Jacobian data is available
        """
        try:
            from scipy.sparse import coo_matrix
        except ImportError:
            raise ImportError(
                "scipy is required for to_sparse_jacobian(). "
                "Install with: pip install scipy"
            )

        if self.jacobian_vals is None or len(self.jacobian_vals) == 0:
            raise ValueError("No Jacobian data available. Run forward with extract_frechet=True.")

        return coo_matrix(
            (self.jacobian_vals, (self.jacobian_rows, self.jacobian_cols)),
            shape=self.jacobian_shape,
        ).tocsr()

    def __repr__(self) -> str:
        """String representation."""
        jac_str = ""
        if self.jacobian_vals is not None:
            jac_str = f", jacobian={self.jacobian_shape}"
        return (
            f"WaveFrontTrackerResult("
            f"arrivals={len(self.travel_times)}, "
            f"raypaths={len(self.raypaths)}, "
            f"sources={len(self.wavefronts)}"
            f"{jac_str})"
        )
