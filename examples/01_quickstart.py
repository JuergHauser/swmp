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
    # Getting started with pyswmp

    **pyswmp** is a Python-wrapped Fortran library for surface wave multipathing
    ray tracing on spherical and Cartesian coordinate systems
    (Hauser et al., 2008).

    This notebook walks through the basic workflow:

    1. Create a velocity model
    2. Define sources and receivers
    3. Configure ray tracing options
    4. Run forward modelling
    5. Inspect and visualise results

    All using the **file-free API** — no disk I/O required.
    """)
    return


@app.cell
def _():
    import numpy as np
    import pyswmp
    return np, pyswmp


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 1. Create a velocity model

    `create_constant_velocity_model()` builds a 2D grid with uniform velocity.
    The model covers a region of Australia at 0.5° spacing with a velocity of
    3.5 km/s (typical Rayleigh wave group velocity at ~20 s period).
    """)
    return


@app.cell
def _(pyswmp):
    model = pyswmp.create_constant_velocity_model(
        nx=100, ny=80,
        x0=110.0, y0=-45.0,
        dx=0.5, dy=0.5,
        velocity=3.5,
    )
    print(model)
    print(f"Domain: ({model.x0}, {model.y0}) to ({model.x1}, {model.y1})")
    print(f"Velocity: {model.velocities.mean():.2f} km/s")
    return (model,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 2. Define sources and receivers

    Sources and receivers are specified as `(longitude, latitude)` arrays.
    `source_type=1` indicates point sources.
    """)
    return


@app.cell
def _(np, pyswmp):
    # Single point source in central Australia
    sources = pyswmp.Sources(
        positions=np.array([[135.0, -25.0]]),
        source_type=1,
    )
    print(sources)

    # Grid of receivers across Australia
    x_rec = np.linspace(115, 155, 20)
    y_rec = np.linspace(-40, -15, 15)
    xx, yy = np.meshgrid(x_rec, y_rec)
    receivers = pyswmp.Receivers(
        positions=np.column_stack([xx.ravel(), yy.ravel()])
    )
    print(receivers)
    return receivers, sources


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3. Configure ray tracing options

    `TrackerOptions` controls the numerical solver, coordinate system, and
    output extraction. Key parameters:

    - `coordinate_system=2` — spherical coordinates (longitude/latitude)
    - `dt` — ODE time step size
    - `max_iterations` — maximum wavefront propagation steps
    - `extract_raypaths=True` — store individual ray paths for plotting
    """)
    return


@app.cell
def _(pyswmp):
    opts = pyswmp.TrackerOptions(
        coordinate_system=2,
        dt=5.0,
        max_iterations=500,
        extract_raypaths=True,
    )
    return (opts,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 4. Run forward modelling

    Create a `WaveFrontTracker` with the model and options, then call
    `forward(sources, receivers)`. The result contains travel times, azimuths,
    arrival numbers, and (optionally) ray paths.
    """)
    return


@app.cell
def _(model, opts, pyswmp, receivers, sources):
    tracker = pyswmp.WaveFrontTracker(model, opts)
    result = tracker.forward(sources, receivers)

    print(f"Arrivals:          {len(result.travel_times)}")
    print(f"Raypaths:          {len(result.raypaths)}")
    print(f"Travel time range: {result.travel_times.min():.1f} – {result.travel_times.max():.1f} s")
    return (result,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 5. Inspect results

    ### DataFrame view

    `to_dataframe()` converts results to a pandas DataFrame with columns:
    `source`, `receiver`, `arrival`, `time`, `azimuth`, `spreading`.
    """)
    return


@app.cell
def _(result):
    df = result.to_dataframe()
    print(f"Shape: {df.shape}")
    df.head(10)
    return (df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Ray path visualisation

    The `Visualisation` class plots ray paths overlaid on the velocity model.
    Each source is coloured differently.
    """)
    return


@app.cell
def _(model, pyswmp, result):
    _vis = pyswmp.Visualisation(model, result)
    _fig = _vis.plot_raypaths()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Travel time distribution

    For a constant velocity model and a single source, travel times
    should be proportional to source–receiver distance.
    """)
    return


@app.cell
def _(np, result):
    import matplotlib.pyplot as plt

    # First-arrival travel times per receiver
    _unique_recs = np.unique(result.receiver_ids)
    _first_times = []
    for _rid in _unique_recs:
        _mask = (result.receiver_ids == _rid) & (result.arrival_numbers == 1)
        if _mask.any():
            _first_times.append(result.travel_times[_mask][0])

    _fig, _ax = plt.subplots(figsize=(7, 4))
    _ax.hist(_first_times, bins=20, color="steelblue", edgecolor="black", alpha=0.7)
    _ax.set(xlabel="First-arrival travel time (s)", ylabel="Count",
            title="Travel time distribution")
    _ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Summary

    The complete file-free workflow in four lines:

    ```python
    model   = pyswmp.create_constant_velocity_model(...)
    opts    = pyswmp.TrackerOptions(coordinate_system=2, extract_raypaths=True)
    tracker = pyswmp.WaveFrontTracker(model, opts)
    result  = tracker.forward(sources, receivers)
    ```

    Next notebooks cover:
    - **02** — Validating accuracy against the analytical great circle solution
    - **03** — Real-world multipathing with the Australian velocity model
    - **04** — Frechet sensitivity kernels and the Jacobian matrix
    - **05** — Parallel execution for multi-source problems
    """)
    return


if __name__ == "__main__":
    app.run()
