"""Plotting utilities for SWMP results.

Provides a Visualisation class that accepts in-memory VelocityModel2D and
WaveFrontTrackerResult objects, using cartopy for map projections.
"""

from collections import defaultdict

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

from .data import VelocityModel2D
from .results import WaveFrontTrackerResult

# tsitsul optimally distinct palette: http://tsitsul.in/blog/coloropt/
_TSITSUL_COLORS = [
    "#ebac23",
    "#b80058",
    "#008cf9",
    "#006e00",
    "#00bbad",
    "#d163e6",
    "#b24502",
    "#ff9287",
    "#5954d6",
    "#00c6f8",
    "#878500",
    "#00a76c",
]

# tsitsul tarnish 6 palette (muted): http://tsitsul.in/blog/coloropt/
_TARNISH_COLORS = [
    "#274d52",
    "#c7a2a6",
    "#818b70",
    "#604e3c",
    "#8c9fb7",
    "#796880",
]


class Visualisation:
    """Plot velocity models, raypaths, and wavefronts on cartopy maps.

    Accepts in-memory data objects from the file-free API.

    Args:
        model: 2D velocity model
        result: Wavefront tracker result (optional)
        interpolate: B-spline upsampling factor for display (default 4,
            set to 1 or None to show raw grid)

    Example:
        >>> vis = Visualisation(model, result)
        >>> fig = vis.plot_raypaths()
    """

    def __init__(self, model: VelocityModel2D, result: WaveFrontTrackerResult = None,
                 interpolate=4):
        self.model = model
        self.result = result
        # B-spline-interpolated model for display (matches what the ray tracer sees)
        if interpolate and interpolate > 1:
            self._display_model = model.interpolate(interpolate, interpolate)
        else:
            self._display_model = model

    def _make_axes(self, ax, coastlines):
        """Create or validate cartopy axes."""
        if ax is not None:
            return ax.figure, ax

        x0, x1, y0, y1 = self.model.extent
        xc = (x0 + x1) / 2.0
        projection = ccrs.Mercator(
            central_longitude=xc,
            min_latitude=y0,
            max_latitude=y1,
        )
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=projection)
        ax.set_extent([x0, x1, y0, y1], crs=ccrs.PlateCarree())
        if coastlines:
            ax.coastlines(resolution="10m", color="black")
        ax.gridlines(color="k", alpha=0.3, draw_labels=True)

        return fig, ax

    def _plot_model_background(self, ax, cmap, **kwargs):
        """Plot B-spline-interpolated velocity model as pcolormesh."""
        m = self._display_model
        x = np.linspace(m.x0, m.x0 + m.dx * (m.nx - 1), m.nx)
        y = np.linspace(m.y0, m.y0 + m.dy * (m.ny - 1), m.ny)
        yy, xx = np.meshgrid(y, x)

        plot_kwargs = dict(cmap=cmap, transform=ccrs.PlateCarree())
        plot_kwargs.update(kwargs)

        return ax.pcolormesh(xx, yy, m.velocities, **plot_kwargs)

    def plot_model(self, ax=None, cmap="Greys_r", coastlines=True, **kwargs):
        """Plot the velocity model.

        Args:
            ax: Cartopy GeoAxes (created automatically if None)
            cmap: Colormap name
            coastlines: Draw coastlines
            **kwargs: Passed to pcolormesh

        Returns:
            matplotlib.figure.Figure
        """
        fig, ax = self._make_axes(ax, coastlines)
        cm = self._plot_model_background(ax, cmap, **kwargs)
        fig.colorbar(cm, ax=ax, orientation="horizontal")
        return fig

    def plot_raypaths(self, ax=None, cmap="Greys_r", coastlines=True, **kwargs):
        """Plot raypaths overlaid on the velocity model.

        Raypaths are coloured by arrival number using the tarnish palette.
        First arrivals are plotted first (underneath), later arrivals on top.

        Args:
            ax: Cartopy GeoAxes (created automatically if None)
            cmap: Colormap name for velocity model
            coastlines: Draw coastlines
            **kwargs: Passed to pcolormesh

        Returns:
            matplotlib.figure.Figure
        """
        fig, ax = self._make_axes(ax, coastlines)
        cm = self._plot_model_background(ax, cmap, **kwargs)

        if self.result is not None and self.result.raypaths:
            transform = ccrs.PlateCarree()

            by_arrival = defaultdict(list)
            for rp in self.result.raypaths:
                by_arrival[rp["arrival"]].append(rp["path"])

            for arrival_num, paths in sorted(by_arrival.items()):
                color = _TARNISH_COLORS[(arrival_num - 1) % len(_TARNISH_COLORS)]
                for path in paths:
                    ax.plot(path[:, 0], path[:, 1], linewidth=1, color=color,
                            transform=transform)

        fig.colorbar(cm, ax=ax, orientation="horizontal")
        return fig

    def plot_wavefronts(self, ax=None, cmap="Greys_r", coastlines=True, **kwargs):
        """Plot wavefronts and raypaths overlaid on the velocity model.

        Wavefronts are drawn as black lines; raypaths are coloured by source.

        Args:
            ax: Cartopy GeoAxes (created automatically if None)
            cmap: Colormap name for velocity model
            coastlines: Draw coastlines
            **kwargs: Passed to pcolormesh

        Returns:
            matplotlib.figure.Figure
        """
        fig, ax = self._make_axes(ax, coastlines)
        cm = self._plot_model_background(ax, cmap, **kwargs)
        transform = ccrs.PlateCarree()

        # Wavefronts (black lines)
        if self.result is not None and self.result.wavefronts:
            for src_id, wf_list in self.result.wavefronts.items():
                for wf in wf_list:
                    ax.plot(wf[:, 0], wf[:, 1], linewidth=1, color="k",
                            transform=transform)

        # Raypaths (coloured by arrival number, first arrivals underneath)
        if self.result is not None and self.result.raypaths:
            by_arrival = defaultdict(list)
            for rp in self.result.raypaths:
                by_arrival[rp["arrival"]].append(rp["path"])

            for arrival_num, paths in sorted(by_arrival.items()):
                color = _TARNISH_COLORS[(arrival_num - 1) % len(_TARNISH_COLORS)]
                for path in paths:
                    ax.plot(path[:, 0], path[:, 1], linewidth=1, color=color,
                            transform=transform)

        fig.colorbar(cm, ax=ax, orientation="horizontal")
        return fig
