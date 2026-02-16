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
    # Surface wave multipathing across Australia

    When surface waves propagate through a **heterogeneous** velocity model,
    wavefronts can fold and develop triplications. This produces **multiple
    arrivals** at a single receiver — a phenomenon called multipathing.

    This notebook loads a real Australian velocity model and demonstrates
    multipathing ray tracing, inspired by Figure 17 from Hauser et al. (2008).

    ## References
    Hauser, J., Sambridge, M. and Rawlinson, N. (2008). Multiarrival wavefront
    tracking and its applications. *Geochem. Geophys. Geosyst.*, 9(11), Q11001.
    https://doi.org/10.1029/2008GC002069
    """)
    return


@app.cell
def _():
    from pathlib import Path

    import numpy as np
    import matplotlib.pyplot as plt
    import pyswmp
    return Path, np, plt, pyswmp


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Load the Australian velocity model

    The data files are stored in `data/` as binary files compatible with SWMP.
    We read them back into Python data objects for the file-free API.
    """)
    return


@app.cell
def _(Path, np, pyswmp):
    _data_dir = Path(__file__).parent / "data"

    # Read velocity model (Fortran text format)
    with open(_data_dir / "current.vel") as _f:
        _x0, _y0 = map(float, _f.readline().split())
        _nx, _ny = map(int, _f.readline().split())
        _dx, _dy = map(float, _f.readline().split())
        _cn = int(_f.readline().strip())
        # Grid includes cushion nodes: (nx + 2*cn) x (ny + 2*cn)
        _total_nx = _nx + 2 * _cn
        _total_ny = _ny + 2 * _cn
        _vals = [float(_f.readline().strip()) for _ in range(_total_nx * _total_ny)]
        _velocities = np.array(_vals, dtype=np.float32).reshape((_total_nx, _total_ny))

    model = pyswmp.VelocityModel2D(
        velocities=_velocities,
        x0=_x0, y0=_y0, dx=_dx, dy=_dy,
        cushion_nodes=_cn,
    )

    print(model)
    print(f"Velocity range: {_velocities.min():.3f} – {_velocities.max():.3f} km/s")
    return (model,)


@app.cell
def _(Path, np, pyswmp):
    _data_dir = Path(__file__).parent / "data"

    # Read sources (Fortran text format: n_sources, then lon lat per line)
    with open(_data_dir / "sources.dat") as _f:
        _n_src = int(_f.readline().strip())
        _src_pos = np.array(
            [list(map(float, _f.readline().split())) for _ in range(_n_src)], dtype=np.float64
        )

    sources = pyswmp.Sources(positions=_src_pos, source_type=1)

    # Read receivers (Fortran text format: n_receivers, then lon lat per line)
    with open(_data_dir / "receivers.dat") as _f:
        _n_rec = int(_f.readline().strip())
        _rec_pos = np.array(
            [list(map(float, _f.readline().split())) for _ in range(_n_rec)], dtype=np.float64
        )

    receivers = pyswmp.Receivers(positions=_rec_pos)

    print(sources)
    print(receivers)
    return receivers, sources


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Visualise the velocity model

    The velocity model shows lateral variations across Australia that will
    cause wavefront triplications and multipathing.
    """)
    return


@app.cell
def _(model, pyswmp):
    vis_model = pyswmp.Visualisation(model)
    fig_model = vis_model.plot_model()
    fig_model
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Run ray tracing with raypath extraction

    We enable `extract_raypaths=True` to visualise individual ray paths.
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

    print(f"Arrivals: {len(result.travel_times)}")
    print(f"Raypaths: {len(result.raypaths)}")
    print(f"Travel time range: {result.travel_times.min():.1f} – {result.travel_times.max():.1f} s")
    return (result,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Raypaths overlaid on velocity model

    Ray paths bend and deviate from great circles due to velocity
    heterogeneity, producing multiple arrivals at some receivers.
    """)
    return


@app.cell
def _(model, pyswmp, result):
    vis = pyswmp.Visualisation(model, result)
    fig_rays = vis.plot_raypaths()
    fig_rays
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Arrival statistics

    How many arrivals does each receiver observe? In a homogeneous model,
    every receiver gets exactly one arrival. With heterogeneity, some
    receivers see 2, 3, or more — that is multipathing.
    """)
    return


@app.cell
def _(np, plt, result):
    # Count arrivals per receiver
    _unique_recs, _arrival_counts = np.unique(result.receiver_ids, return_counts=True)

    _fig, (_ax1, _ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Arrival count histogram
    _ax1.hist(_arrival_counts, bins=range(1, _arrival_counts.max() + 2), align="left",
              color="steelblue", edgecolor="black", alpha=0.7)
    _ax1.set(xlabel="Number of arrivals", ylabel="Number of receivers",
             title="Arrivals per receiver")
    _ax1.grid(True, alpha=0.3, axis="y")

    _multi = (_arrival_counts > 1).sum()
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

    - **Model**: Real Australian surface wave velocity model (Hauser et al., 2008)
    - **Multipathing**: Velocity heterogeneity causes wavefront triplications,
      producing multiple arrivals at individual receivers
    - **API**: Binary data files loaded into Python objects, then passed to the
      file-free `WaveFrontTracker` API
    """)
    return


if __name__ == "__main__":
    app.run()
