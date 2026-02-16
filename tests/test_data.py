"""Unit tests for pyswmp.data module."""

import numpy as np
import pytest

from pyswmp import VelocityModel2D, Sources, Receivers
from pyswmp import create_constant_velocity_model, create_gradient_velocity_model


class TestVelocityModel2D:
    """Tests for VelocityModel2D dataclass."""

    def test_valid_model_creation(self):
        """Test creating a valid velocity model."""
        velocities = np.ones((10, 20), dtype=np.float32)
        model = VelocityModel2D(
            velocities=velocities,
            x0=110.0,
            y0=-45.0,
            dx=0.5,
            dy=0.5,
            cushion_nodes=3
        )

        assert model.nx == 10
        assert model.ny == 20
        assert model.x0 == 110.0
        assert model.y0 == -45.0
        assert model.dx == 0.5
        assert model.dy == 0.5
        assert model.cushion_nodes == 3

    def test_automatic_dtype_conversion(self):
        """Test that velocities are automatically converted to float32."""
        velocities = np.ones((10, 20), dtype=np.float64)
        model = VelocityModel2D(velocities, 0, 0, 1, 1)

        assert model.velocities.dtype == np.float32

    def test_invalid_velocity_shape(self):
        """Test that 1D or 3D arrays raise ValueError."""
        with pytest.raises(ValueError, match="must be 2D array"):
            VelocityModel2D(np.ones(10), 0, 0, 1, 1)

        with pytest.raises(ValueError, match="must be 2D array"):
            VelocityModel2D(np.ones((10, 20, 30)), 0, 0, 1, 1)

    def test_invalid_grid_spacing(self):
        """Test that non-positive grid spacings raise ValueError."""
        velocities = np.ones((10, 20))

        with pytest.raises(ValueError, match="must be positive"):
            VelocityModel2D(velocities, 0, 0, -1, 1)

        with pytest.raises(ValueError, match="must be positive"):
            VelocityModel2D(velocities, 0, 0, 1, 0)

    def test_invalid_cushion_nodes(self):
        """Test that negative cushion_nodes raises ValueError."""
        velocities = np.ones((10, 20))

        with pytest.raises(ValueError, match="must be non-negative"):
            VelocityModel2D(velocities, 0, 0, 1, 1, cushion_nodes=-1)

    def test_computed_properties(self):
        """Test computed properties (x1, y1, extent)."""
        model = VelocityModel2D(
            velocities=np.ones((10, 20)),
            x0=110.0,
            y0=-45.0,
            dx=0.5,
            dy=0.4
        )

        assert model.x1 == 110.0 + 0.5 * 10
        assert model.y1 == -45.0 + 0.4 * 20
        assert model.extent == (110.0, 115.0, -45.0, -37.0)

    def test_repr(self):
        """Test string representation."""
        model = VelocityModel2D(np.ones((10, 20)) * 3.5, 0, 0, 1, 1)
        repr_str = repr(model)

        assert "VelocityModel2D" in repr_str
        assert "shape=(10, 20)" in repr_str
        assert "3.500" in repr_str


class TestSources:
    """Tests for Sources dataclass."""

    def test_single_point_source(self):
        """Test creating a single point source."""
        sources = Sources(
            positions=np.array([[135.0, -25.0]]),
            source_type=1
        )

        assert sources.n_sources == 1
        assert sources.source_type == 1
        assert sources.positions.shape == (1, 2)

    def test_multiple_point_sources(self):
        """Test creating multiple point sources."""
        positions = np.array([
            [130.0, -30.0],
            [135.0, -25.0],
            [140.0, -20.0]
        ])
        sources = Sources(positions, source_type=1)

        assert sources.n_sources == 3
        assert sources.positions.shape == (3, 2)

    def test_single_source_1d_input(self):
        """Test that 1D input is reshaped to (1, 2)."""
        sources = Sources(positions=np.array([135.0, -25.0]))

        assert sources.positions.shape == (1, 2)
        assert sources.n_sources == 1

    def test_plane_wave_sources_with_angles(self):
        """Test plane wave sources with angles."""
        sources = Sources(
            positions=np.array([[135.0, -25.0]]),
            source_type=2,
            angles=np.array([45.0])
        )

        assert sources.source_type == 2
        assert sources.angles[0] == 45.0

    def test_invalid_source_type(self):
        """Test that invalid source_type raises ValueError."""
        with pytest.raises(ValueError, match="must be 1.*or 2"):
            Sources(np.array([[135.0, -25.0]]), source_type=3)

    def test_mismatched_angles_length(self):
        """Test that mismatched angles length raises ValueError."""
        with pytest.raises(ValueError, match="angles length.*must match"):
            Sources(
                positions=np.array([[135.0, -25.0], [140.0, -20.0]]),
                source_type=2,
                angles=np.array([45.0])  # Only 1 angle for 2 sources
            )

    def test_automatic_zero_angles(self):
        """Test that angles are automatically created as zeros if not provided."""
        sources = Sources(positions=np.array([[135.0, -25.0]]), source_type=1)

        assert len(sources.angles) == 1
        assert sources.angles[0] == 0.0

    def test_invalid_positions_shape(self):
        """Test that invalid positions shape raises ValueError."""
        with pytest.raises(ValueError):
            Sources(positions=np.array([135.0]))  # Only 1 value, need 2

    def test_repr(self):
        """Test string representation."""
        sources = Sources(np.array([[135.0, -25.0]]))
        repr_str = repr(sources)

        assert "Sources" in repr_str
        assert "n=1" in repr_str
        assert "point" in repr_str


class TestReceivers:
    """Tests for Receivers dataclass."""

    def test_single_receiver(self):
        """Test creating a single receiver."""
        receivers = Receivers(positions=np.array([[120.0, -30.0]]))

        assert receivers.n_receivers == 1
        assert receivers.positions.shape == (1, 2)

    def test_multiple_receivers(self):
        """Test creating multiple receivers."""
        positions = np.array([
            [120.0, -40.0],
            [130.0, -30.0],
            [140.0, -20.0]
        ])
        receivers = Receivers(positions)

        assert receivers.n_receivers == 3
        assert receivers.positions.shape == (3, 2)

    def test_single_receiver_1d_input(self):
        """Test that 1D input is reshaped to (1, 2)."""
        receivers = Receivers(positions=np.array([120.0, -30.0]))

        assert receivers.positions.shape == (1, 2)
        assert receivers.n_receivers == 1

    def test_grid_receivers(self):
        """Test creating a grid of receivers."""
        x = np.linspace(110, 160, 10)
        y = np.linspace(-45, -10, 8)
        xx, yy = np.meshgrid(x, y)
        positions = np.column_stack([xx.ravel(), yy.ravel()])

        receivers = Receivers(positions)

        assert receivers.n_receivers == 80  # 10 * 8

    def test_invalid_positions_shape(self):
        """Test that invalid positions shape raises ValueError."""
        with pytest.raises(ValueError):
            Receivers(positions=np.array([120.0]))  # Only 1 value

    def test_repr(self):
        """Test string representation."""
        receivers = Receivers(np.array([[120.0, -30.0]]))
        repr_str = repr(receivers)

        assert "Receivers" in repr_str
        assert "n=1" in repr_str


class TestHelperFunctions:
    """Tests for helper functions."""

    def test_create_constant_velocity_model(self):
        """Test creating a constant velocity model."""
        model = create_constant_velocity_model(
            nx=100, ny=80,
            x0=110.0, y0=-45.0,
            dx=0.5, dy=0.5,
            velocity=3.5
        )

        assert model.nx == 100
        assert model.ny == 80
        assert np.all(model.velocities == 3.5)
        assert model.velocities.dtype == np.float32

    def test_create_gradient_velocity_model(self):
        """Test creating a gradient velocity model."""
        model = create_gradient_velocity_model(
            nx=10, ny=10,
            x0=0.0, y0=0.0,
            dx=1.0, dy=1.0,
            v0=3.0,
            gradient_x=0.1,
            gradient_y=0.0
        )

        # Check that velocity increases in X direction
        assert model.velocities[0, 0] == pytest.approx(3.0)
        assert model.velocities[-1, 0] == pytest.approx(3.9)  # 3.0 + 0.1 * 9

        # Check that velocity is constant in Y direction
        assert np.allclose(model.velocities[:, 0], model.velocities[:, -1])

    def test_create_gradient_velocity_model_y_gradient(self):
        """Test creating a gradient velocity model with Y gradient."""
        model = create_gradient_velocity_model(
            nx=10, ny=10,
            x0=0.0, y0=0.0,
            dx=1.0, dy=1.0,
            v0=3.0,
            gradient_x=0.0,
            gradient_y=0.2
        )

        # Check that velocity increases in Y direction
        assert model.velocities[0, 0] == pytest.approx(3.0)
        assert model.velocities[0, -1] == pytest.approx(4.8)  # 3.0 + 0.2 * 9


class TestFileIO:
    """Tests for file I/O functions (integration tests)."""

    def test_write_read_velocity_model_file(self, tmp_path):
        """Test writing and reading velocity model files."""
        from pyswmp import write_velocity_model_file

        # Create a test model
        model = create_constant_velocity_model(10, 10, 0, 0, 1, 1, 3.5)

        # Write to file
        output_file = tmp_path / "test_model.vel"
        write_velocity_model_file(model, str(output_file))

        # Check file exists and has content
        assert output_file.exists()
        assert output_file.stat().st_size > 0

    def test_write_sources_file(self, tmp_path):
        """Test writing sources file."""
        from pyswmp import write_sources_file

        sources = Sources(np.array([[135.0, -25.0], [140.0, -20.0]]))

        output_file = tmp_path / "test_sources.dat"
        write_sources_file(sources, str(output_file))

        assert output_file.exists()
        assert output_file.stat().st_size > 0

    def test_write_receivers_file(self, tmp_path):
        """Test writing receivers file."""
        from pyswmp import write_receivers_file

        receivers = Receivers(np.array([[120.0, -30.0], [130.0, -25.0]]))

        output_file = tmp_path / "test_receivers.dat"
        write_receivers_file(receivers, str(output_file))

        assert output_file.exists()
        assert output_file.stat().st_size > 0
