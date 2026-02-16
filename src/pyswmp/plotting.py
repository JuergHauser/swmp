"""Plotting utilities for SWMP results.

Provides a Visualisation class that accepts in-memory VelocityModel2D and
WaveFrontTrackerResult objects, with optional cartopy support.
"""

from collections import defaultdict

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


def _try_import_cartopy():
    try:
        import cartopy

        return cartopy
    except ImportError:
        return None


class Visualisation:
    """Plot velocity models, raypaths, and wavefronts.

    Accepts in-memory data objects from the file-free API.

    Args:
        model: 2D velocity model
        result: Wavefront tracker result (optional)

    Example:
        >>> vis = Visualisation(model, result)
        >>> fig = vis.plot_raypaths()
    """

    def __init__(self, model: VelocityModel2D, result: WaveFrontTrackerResult = None):
        self.model = model
        self.result = result

    def _make_axes(self, ax, coastlines):
        """Create or validate axes, using cartopy if available."""
        import matplotlib.pyplot as plt

        if ax is not None:
            return ax.figure, ax

        cartopy = _try_import_cartopy()
        if cartopy is not None:
            x0, x1, y0, y1 = self.model.extent
            xc = (x0 + x1) / 2.0
            fig = plt.figure()
            projection = cartopy.crs.Mercator(
                central_longitude=xc,
                min_latitude=y0,
                max_latitude=y1,
            )
            ax = fig.add_subplot(1, 1, 1, projection=projection)
            ax.set_extent([x0, x1, y0, y1], crs=cartopy.crs.PlateCarree())
            if coastlines:
                ax.coastlines(resolution="10m", color="black")
            ax.gridlines(color="k")
        else:
            fig, ax = plt.subplots()
            ax.set_xlim(self.model.x0, self.model.x1)
            ax.set_ylim(self.model.y0, self.model.y1)
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            ax.grid(True, alpha=0.3)

        return fig, ax

    def _plot_model_background(self, ax, cmap, **kwargs):
        """Plot velocity model as pcolormesh."""
        m = self.model
        x = np.linspace(m.x0, m.x1, m.nx)
        y = np.linspace(m.y0, m.y1, m.ny)
        yy, xx = np.meshgrid(y, x)

        cartopy = _try_import_cartopy()
        plot_kwargs = dict(cmap=cmap)
        plot_kwargs.update(kwargs)
        if cartopy is not None and hasattr(ax, "projection"):
            plot_kwargs["transform"] = cartopy.crs.PlateCarree()

        return ax.pcolormesh(xx, yy, m.velocities, **plot_kwargs)

    def _get_transform_kw(self, ax):
        """Return cartopy transform kwargs if applicable."""
        cartopy = _try_import_cartopy()
        if cartopy is not None and hasattr(ax, "projection"):
            return {"transform": cartopy.crs.PlateCarree()}
        return {}

    def plot_model(self, ax=None, cmap="Greys_r", coastlines=True, **kwargs):
        """Plot the velocity model.

        Args:
            ax: Matplotlib axes (created automatically if None)
            cmap: Colormap name
            coastlines: Draw coastlines (requires cartopy)
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

        Raypaths are coloured by source using the tsitsul palette.

        Args:
            ax: Matplotlib axes (created automatically if None)
            cmap: Colormap name for velocity model
            coastlines: Draw coastlines (requires cartopy)
            **kwargs: Passed to pcolormesh

        Returns:
            matplotlib.figure.Figure
        """
        fig, ax = self._make_axes(ax, coastlines)
        cm = self._plot_model_background(ax, cmap, **kwargs)

        if self.result is not None and self.result.raypaths:
            transform_kw = self._get_transform_kw(ax)

            by_source = defaultdict(list)
            for rp in self.result.raypaths:
                by_source[rp["source"]].append(rp["path"])

            for ic, (src_id, paths) in enumerate(sorted(by_source.items())):
                color = _TSITSUL_COLORS[ic % len(_TSITSUL_COLORS)]
                for path in paths:
                    ax.plot(path[:, 0], path[:, 1], linewidth=1, color=color, **transform_kw)

        fig.colorbar(cm, ax=ax, orientation="horizontal")
        return fig

    def plot_wavefronts(self, ax=None, cmap="Greys_r", coastlines=True, **kwargs):
        """Plot wavefronts and raypaths overlaid on the velocity model.

        Wavefronts are drawn as black lines; raypaths are coloured by source.

        Args:
            ax: Matplotlib axes (created automatically if None)
            cmap: Colormap name for velocity model
            coastlines: Draw coastlines (requires cartopy)
            **kwargs: Passed to pcolormesh

        Returns:
            matplotlib.figure.Figure
        """
        fig, ax = self._make_axes(ax, coastlines)
        cm = self._plot_model_background(ax, cmap, **kwargs)
        transform_kw = self._get_transform_kw(ax)

        # Wavefronts (black lines)
        if self.result is not None and self.result.wavefronts:
            for src_id, wf_list in self.result.wavefronts.items():
                for wf in wf_list:
                    ax.plot(wf[:, 0], wf[:, 1], linewidth=1, color="k", **transform_kw)

        # Raypaths (coloured by source)
        if self.result is not None and self.result.raypaths:
            by_source = defaultdict(list)
            for rp in self.result.raypaths:
                by_source[rp["source"]].append(rp["path"])

            for ic, (src_id, paths) in enumerate(sorted(by_source.items())):
                color = _TSITSUL_COLORS[ic % len(_TSITSUL_COLORS)]
                for path in paths:
                    ax.plot(path[:, 0], path[:, 1], linewidth=1, color=color, **transform_kw)

        fig.colorbar(cm, ax=ax, orientation="horizontal")
        return fig
