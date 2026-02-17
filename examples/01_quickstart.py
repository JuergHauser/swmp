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

    1. Create a velocity model with anomalies
    2. Define sources and receivers
    3. Run forward modelling
    4. Inspect multipathing results

    All using the **file-free API** — no disk I/O required.
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
    ## 1. Create a velocity model

    Start from a constant background (3.5 km/s — typical Rayleigh wave group
    velocity at ~20 s period) and add strong Gaussian anomalies. A large slow
    anomaly causes wavefront triplications that produce **multiple arrivals**
    at individual receivers — the multipathing phenomenon that gives the library
    its name.
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

    # Large slow anomaly — strong enough to tripliccate the wavefront
    model = pyswmp.add_gaussian_anomaly(
        model, cx=130.0, cy=-30.0, sigma_x=5.0, sigma_y=5.0, amplitude=-1.5,
    )
    # Second slow anomaly in the south-east
    model = pyswmp.add_gaussian_anomaly(
        model, cx=140.0, cy=-38.0, sigma_x=4.0, sigma_y=3.0, amplitude=-1.0,
    )

    print(model)
    print(f"Velocity range: {model.velocities.min():.2f} – {model.velocities.max():.2f} km/s")
    return (model,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 2. Define sources and receivers

    Sources and receivers are specified as `(longitude, latitude)` arrays.
    `source_type=1` indicates point sources. A grid of receivers covers
    the domain so we can observe multipathing across the model.
    """)
    return


@app.cell
def _(np, pyswmp):
    sources = pyswmp.Sources(
        positions=np.array([[145.0, -18.0]]),
        source_type=1,
    )
    print(sources)

    # Grid of receivers across Australia
    _xx, _yy = np.meshgrid(np.linspace(115, 155, 20), np.linspace(-40, -15, 15))
    receivers = pyswmp.Receivers(
        positions=np.column_stack([_xx.ravel(), _yy.ravel()])
    )
    print(receivers)
    return receivers, sources


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3. Run forward modelling

    Create a `WaveFrontTracker` with the model and options, then call
    `forward(sources, receivers)`.
    """)
    return


@app.cell
def _(model, pyswmp, receivers, sources):
    opts = pyswmp.TrackerOptions(
        coordinate_system=2,
        dt=5.0,
        max_iterations=500,
        extract_raypaths=True,
    )

    tracker = pyswmp.WaveFrontTracker(model, opts)
    result = tracker.forward(sources, receivers)

    print(f"Arrivals:          {len(result.travel_times)}")
    print(f"Raypaths:          {len(result.raypaths)}")
    print(f"Travel time range: {result.travel_times.min():.1f} – {result.travel_times.max():.1f} s")
    return opts, result, tracker


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 4. Inspect results

    ### Ray paths overlaid on the velocity model

    The `Visualisation` class plots the B-spline-interpolated velocity field
    (matching what the ray tracer sees internally) with coastlines via cartopy.
    Rays bend around the slow anomaly — some receivers see multiple arrivals
    from different directions.
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
    ### Multipathing statistics

    How many arrivals does each receiver observe? In a homogeneous model
    every receiver gets exactly one. With strong velocity anomalies,
    wavefront triplications produce 2 or 3 arrivals at some receivers.
    """)
    return


@app.cell
def _(np, plt, result):
    _unique_recs, _counts = np.unique(result.receiver_ids, return_counts=True)
    _multi = (_counts > 1).sum()

    _fig, (_ax1, _ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Arrival count histogram
    _ax1.hist(_counts, bins=range(1, _counts.max() + 2), align="left",
              color="steelblue", edgecolor="black", alpha=0.7)
    _ax1.set(xlabel="Number of arrivals", ylabel="Number of receivers",
             title="Arrivals per receiver")
    _ax1.grid(True, alpha=0.3, axis="y")
    _ax1.annotate(f"{_multi}/{len(_unique_recs)} receivers\nwith multipathing",
                  xy=(0.95, 0.95), xycoords="axes fraction", ha="right", va="top",
                  fontsize=11, bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.8))

    # Travel time histogram (all arrivals)
    _ax2.hist(result.travel_times, bins=30, color="coral", edgecolor="black", alpha=0.7)
    _ax2.set(xlabel="Travel time (s)", ylabel="Count",
             title="Travel time distribution (all arrivals)")
    _ax2.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Summary

    The complete file-free workflow:

    ```python
    model   = pyswmp.create_constant_velocity_model(...)
    model   = pyswmp.add_gaussian_anomaly(model, cx=..., cy=..., ...)
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
