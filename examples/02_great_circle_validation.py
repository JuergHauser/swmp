import marimo

__generated_with = "0.19.7"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Great Circle Ray Tracing Validation

    On a constant-velocity spherical Earth, surface wave rays follow **great
    circle** paths. This provides an analytical reference solution:

    $$t = \frac{R \cdot \theta}{v}$$

    where $\theta$ is the angular distance (haversine formula), $R$ is the Earth
    radius, and $v$ is the constant velocity.

    This notebook compares pyswmp's numerical travel times against this
    analytical solution to quantify accuracy.
    """)
    return


@app.cell
def _():
    import numpy as np
    import matplotlib.pyplot as plt
    import pyswmp
    return np, plt, pyswmp


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Create a constant velocity model

    Global coverage at 2° resolution with a uniform 4 km/s velocity.
    """)
    return


@app.cell
def _(pyswmp):
    model = pyswmp.create_constant_velocity_model(
        nx=180, ny=90,
        x0=-180.0, y0=-90.0,
        dx=2.0, dy=2.0,
        velocity=4.0,
        cushion_nodes=3,
    )
    print(model)
    return (model,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Sources and receivers

    A single source near Sydney, with a grid of receivers across the Pacific.
    """)
    return


@app.cell
def _(np, pyswmp):
    source_lon, source_lat = 151.2, -33.9

    sources = pyswmp.Sources(
        positions=np.array([[source_lon, source_lat]]),
        source_type=1,
    )

    lon_rec = np.linspace(120, -120, 30)
    lat_rec = np.linspace(-60, 60, 20)
    lon_grid, lat_grid = np.meshgrid(lon_rec, lat_rec)

    receivers = pyswmp.Receivers(
        positions=np.column_stack([lon_grid.ravel(), lat_grid.ravel()])
    )
    print(f"{sources}\n{receivers}")
    return lat_grid, lat_rec, lon_grid, lon_rec, receivers, source_lat, source_lon, sources


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Analytical great circle travel times

    $$t = \frac{R \cdot \theta}{v}$$

    where $\theta$ is the angular distance from the haversine formula.
    """)
    return


@app.cell
def _(np, receivers, source_lat, source_lon):
    def great_circle_distance(lon1, lat1, lon2, lat2, radius=6371.0):
        """Great circle distance via haversine (km)."""
        lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
        return radius * 2 * np.arcsin(np.sqrt(a))

    distances = great_circle_distance(
        source_lon, source_lat,
        receivers.positions[:, 0], receivers.positions[:, 1],
    )
    analytical_times = distances / 4.0

    print(f"Distance range: {distances.min():.0f} – {distances.max():.0f} km")
    print(f"Travel time range: {analytical_times.min():.0f} – {analytical_times.max():.0f} s")
    return analytical_times, distances, great_circle_distance


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Run SWMP ray tracing
    """)
    return


@app.cell
def _(model, pyswmp, receivers, sources):
    opts = pyswmp.TrackerOptions(
        dt=5.0,
        max_iterations=200,
        coordinate_system=2,
        earth_radius=6371.0,
        computation_mode=1,
        max_arrivals=5,
        ode_solver=2,
        interpolator=2,
    )

    tracker = pyswmp.WaveFrontTracker(model, opts)
    result = tracker.forward(sources, receivers)

    print(f"Computed {len(result.travel_times)} arrivals")
    print(f"SWMP travel time range: {result.travel_times.min():.0f} – {result.travel_times.max():.0f} s")
    return (result,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Error analysis

    Compare first-arrival numerical travel times against the analytical solution.
    """)
    return


@app.cell
def _(analytical_times, np, receivers, result):
    # Map receiver_id -> first-arrival travel time
    receiver_to_time = {}
    for rec_id, arr_num, tt in zip(result.receiver_ids, result.arrival_numbers, result.travel_times):
        if arr_num == 1:
            receiver_to_time[rec_id] = tt

    numerical_times = np.array(
        [receiver_to_time.get(i + 1, np.nan) for i in range(receivers.n_receivers)]
    )

    valid = ~np.isnan(numerical_times)
    errors = numerical_times[valid] - analytical_times[valid]
    rel_errors = (errors / analytical_times[valid]) * 100

    print("=" * 55)
    print("SWMP vs Analytical Great Circle")
    print("=" * 55)
    print(f"Receivers with arrivals: {valid.sum()} / {receivers.n_receivers}")
    print(f"  Mean abs error:   {np.abs(errors).mean():.3f} s")
    print(f"  Max abs error:    {np.abs(errors).max():.3f} s")
    print(f"  RMS error:        {np.sqrt((errors ** 2).mean()):.3f} s")
    print(f"  Mean |rel error|: {np.abs(rel_errors).mean():.4f} %")
    print(f"  Max  |rel error|: {np.abs(rel_errors).max():.4f} %")

    within_1pct = (np.abs(rel_errors) < 1.0).sum()
    print(f"  Within 1% error:  {within_1pct}/{valid.sum()} ({100 * within_1pct / valid.sum():.0f}%)")
    print("=" * 55)
    return errors, numerical_times, rel_errors, valid


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Diagnostic plots

    Four views of the error:

    1. **1:1 comparison** — numerical vs analytical travel times with a ±1% band
    2. **Error vs distance** — coloured by relative error magnitude
    3. **Error distribution** — histogram of relative errors
    4. **Spatial error map** — contour plot of relative error across the receiver grid
    """)
    return


@app.cell
def _(
    analytical_times,
    distances,
    errors,
    lat_grid,
    lon_grid,
    np,
    numerical_times,
    plt,
    rel_errors,
    source_lat,
    source_lon,
    valid,
):
    _fig, _axes = plt.subplots(2, 2, figsize=(15, 12))

    # 1:1 comparison
    _ax = _axes[0, 0]
    _ax.scatter(analytical_times[valid], numerical_times[valid], c="blue", s=30, alpha=0.6,
                edgecolors="k", linewidth=0.5)
    _lo = min(analytical_times[valid].min(), numerical_times[valid].min())
    _hi = max(analytical_times[valid].max(), numerical_times[valid].max())
    _ax.plot([_lo, _hi], [_lo, _hi], "r--", lw=2, label="Perfect agreement")
    _ax.fill_between([_lo, _hi], [_lo * 0.99, _hi * 0.99], [_lo * 1.01, _hi * 1.01],
                     color="red", alpha=0.1, label="±1%")
    _ax.set_xlabel("Analytical (s)")
    _ax.set_ylabel("SWMP (s)")
    _ax.set_title("1:1 Comparison", fontweight="bold")
    _ax.legend()
    _ax.grid(True, alpha=0.3)
    _ax.set_aspect("equal")

    # Error vs distance
    _ax = _axes[0, 1]
    _sc = _ax.scatter(distances[valid] / 1000, errors, c=np.abs(rel_errors), s=30,
                      cmap="YlOrRd", alpha=0.7, edgecolors="k", linewidth=0.5)
    _ax.axhline(0, color="red", ls="--", lw=2)
    _ax.set_xlabel("Distance (x1000 km)")
    _ax.set_ylabel("Error (s)")
    _ax.set_title("Error vs Distance", fontweight="bold")
    _ax.grid(True, alpha=0.3)
    plt.colorbar(_sc, ax=_ax, label="|Relative Error| (%)")

    # Error distribution
    _ax = _axes[1, 0]
    _ax.hist(rel_errors, bins=30, color="steelblue", alpha=0.7, edgecolor="black")
    _ax.axvline(0, color="red", ls="--", lw=2)
    _ax.axvline(rel_errors.mean(), color="green", lw=2, label=f"Mean = {rel_errors.mean():.3f}%")
    _ax.set_xlabel("Relative Error (%)")
    _ax.set_ylabel("Count")
    _ax.set_title("Error Distribution", fontweight="bold")
    _ax.legend()
    _ax.grid(True, alpha=0.3, axis="y")

    # Spatial error map
    _ax = _axes[1, 1]
    _rel_grid = np.full(lon_grid.shape, np.nan)
    _valid_indices = np.where(valid)[0]
    for _k, _idx in enumerate(_valid_indices):
        _i, _j = np.unravel_index(_idx, lon_grid.shape)
        _rel_grid[_i, _j] = rel_errors[_k]
    _im = _ax.contourf(lon_grid, lat_grid, _rel_grid, levels=np.linspace(-1, 1, 21),
                       cmap="RdBu_r", extend="both")
    _ax.scatter(source_lon, source_lat, c="yellow", s=300, marker="*",
                edgecolors="black", linewidths=2, label="Source", zorder=5)
    _ax.set_xlabel("Longitude")
    _ax.set_ylabel("Latitude")
    _ax.set_title("Spatial Relative Error", fontweight="bold")
    _ax.legend()
    plt.colorbar(_im, ax=_ax, label="Relative Error (%)")
    _ax.grid(True, alpha=0.3)

    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Geometry and travel time map
    """)
    return


@app.cell
def _(analytical_times, distances, lat_grid, lat_rec, lon_grid, lon_rec,
      numerical_times, plt, source_lat, source_lon, valid):
    _fig, (_ax1, _ax2) = plt.subplots(1, 2, figsize=(14, 5))

    _ax1.scatter(lon_grid.ravel(), lat_grid.ravel(), c="blue", s=10, alpha=0.3, label="Receivers")
    _ax1.scatter(source_lon, source_lat, c="red", s=200, marker="*", label="Source", zorder=5)
    _ax1.set(xlabel="Longitude", ylabel="Latitude", title="Source / Receiver Layout",
             xlim=(-180, 180), ylim=(-90, 90))
    _ax1.legend()
    _ax1.grid(True, alpha=0.3)

    _ax2.scatter(distances[valid], analytical_times[valid], s=20, alpha=0.5, label="Analytical")
    _ax2.scatter(distances[valid], numerical_times[valid], s=10, alpha=0.5, marker="x", label="SWMP")
    _ax2.plot([distances.min(), distances.max()], [distances.min() / 4, distances.max() / 4],
              "k--", lw=1, label="v = 4 km/s")
    _ax2.set(xlabel="Distance (km)", ylabel="Travel time (s)", title="Travel Time vs Distance")
    _ax2.legend()
    _ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    _fig
    return


@app.cell
def _(analytical_times, lat_grid, lon_grid, lon_rec, lat_rec, plt, source_lat, source_lon):
    _times_grid = analytical_times.reshape(lon_grid.shape)

    _fig, _ax = plt.subplots(figsize=(12, 6))
    _c = _ax.contourf(lon_grid, lat_grid, _times_grid / 60, levels=20, cmap="viridis")
    _ax.contour(lon_grid, lat_grid, _times_grid / 60, levels=10, colors="white", alpha=0.3,
                linewidths=0.5)
    _ax.scatter(source_lon, source_lat, c="red", s=300, marker="*", edgecolors="white",
                linewidths=2, label="Source", zorder=5)
    plt.colorbar(_c, ax=_ax, label="Travel time (min)")
    _ax.set(xlabel="Longitude", ylabel="Latitude",
            title="Great Circle Travel Times from Sydney (v = 4 km/s)",
            xlim=(lon_rec.min(), lon_rec.max()), ylim=(lat_rec.min(), lat_rec.max()))
    _ax.legend()
    _ax.grid(True, alpha=0.3)
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Summary

    - **Model**: Constant velocity (4 km/s) spherical Earth at 2° resolution
    - **Validation**: Numerical travel times match the analytical great circle solution to < 1% relative error
    - **API**: Uses the file-free `WaveFrontTracker(model, opts)` + `forward(sources, receivers)` pattern

    Try varying `dt`, `ode_solver`, or grid resolution to observe convergence behaviour.
    """)
    return


if __name__ == "__main__":
    app.run()
