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
    # Parallel Execution

    When tracing rays from multiple sources, each source can be processed
    independently. pyswmp supports parallel execution via a user-provided
    `concurrent.futures` pool.

    **Why processes, not threads?** The Fortran core uses **module-level state**
    (velocity model, receivers, wavefronts, etc.). Each worker must load its own
    independent `libswmp` instance, so `ProcessPoolExecutor` is required —
    `ThreadPoolExecutor` would share state and produce incorrect results.
    """)
    return


@app.cell
def _():
    import time
    from concurrent.futures import ProcessPoolExecutor

    import numpy as np
    import pyswmp
    return ProcessPoolExecutor, np, pyswmp, time


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Setup: multiple sources

    We create 8 point sources spread across Australia with a grid of receivers.
    More sources = more benefit from parallelism.
    """)
    return


@app.cell
def _(np, pyswmp):
    model = pyswmp.create_constant_velocity_model(
        nx=100, ny=80,
        x0=110.0, y0=-45.0,
        dx=0.5, dy=0.5,
        velocity=3.5,
    )

    sources = pyswmp.Sources(
        positions=np.array([
            [120.0, -30.0],
            [125.0, -25.0],
            [130.0, -35.0],
            [135.0, -20.0],
            [140.0, -30.0],
            [145.0, -25.0],
            [130.0, -40.0],
            [150.0, -35.0],
        ]),
        source_type=1,
    )

    x_rec = np.linspace(115, 155, 15)
    y_rec = np.linspace(-40, -15, 12)
    xx, yy = np.meshgrid(x_rec, y_rec)
    receivers = pyswmp.Receivers(
        positions=np.column_stack([xx.ravel(), yy.ravel()])
    )

    opts = pyswmp.TrackerOptions(coordinate_system=2, dt=5.0, max_iterations=500)
    tracker = pyswmp.WaveFrontTracker(model, opts)

    print(f"Sources:   {sources.n_sources}")
    print(f"Receivers: {receivers.n_receivers}")
    return model, opts, receivers, sources, tracker


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Sequential execution (baseline)
    """)
    return


@app.cell
def _(receivers, sources, time, tracker):
    _t0 = time.time()
    result_seq = tracker.forward(sources, receivers)
    t_seq = time.time() - _t0

    print(f"Sequential: {len(result_seq.travel_times)} arrivals in {t_seq:.3f}s")
    return result_seq, t_seq


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Parallel execution

    Pass a `ProcessPoolExecutor` to `forward(pool=...)`. Each source is
    dispatched to a separate worker process that creates its own Fortran
    library instance.
    """)
    return


@app.cell
def _(ProcessPoolExecutor, receivers, sources, time, tracker):
    _t0 = time.time()
    with ProcessPoolExecutor(max_workers=4) as _pool:
        result_par = tracker.forward(sources, receivers, pool=_pool)
    t_par = time.time() - _t0

    print(f"Parallel (4 workers): {len(result_par.travel_times)} arrivals in {t_par:.3f}s")
    return result_par, t_par


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Compare results
    """)
    return


@app.cell
def _(np, result_par, result_seq, t_par, t_seq):
    # Sort both result sets for comparison
    _seq_order = np.lexsort((result_seq.arrival_numbers, result_seq.receiver_ids,
                             result_seq.source_ids))
    _par_order = np.lexsort((result_par.arrival_numbers, result_par.receiver_ids,
                             result_par.source_ids))

    _times_match = np.allclose(
        result_seq.travel_times[_seq_order],
        result_par.travel_times[_par_order],
        atol=1e-6,
    )

    _speedup = t_seq / t_par if t_par > 0 else float("inf")

    print("=" * 50)
    print("Comparison")
    print("=" * 50)
    print(f"Sequential:    {t_seq:.3f}s  ({len(result_seq.travel_times)} arrivals)")
    print(f"Parallel (4):  {t_par:.3f}s  ({len(result_par.travel_times)} arrivals)")
    print(f"Speedup:       {_speedup:.2f}x")
    print(f"Results match: {_times_match}")
    print("=" * 50)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Summary

    - Pass a `ProcessPoolExecutor` (or schwimmbad pool for MPI) to
      `tracker.forward(sources, receivers, pool=pool)`
    - Each worker creates an independent Fortran library instance —
      `ThreadPoolExecutor` is **not supported**
    - Results are identical regardless of execution mode
    - Speedup depends on the number of sources and available cores

    For HPC environments with MPI, use a `schwimmbad.MPIPool` instead:

    ```python
    from schwimmbad import MPIPool
    with MPIPool() as pool:
        result = tracker.forward(sources, receivers, pool=pool)
    ```
    """)
    return


if __name__ == "__main__":
    app.run()
