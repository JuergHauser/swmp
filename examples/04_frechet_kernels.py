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
    # Frechet Derivatives and Sensitivity Kernels

    The **Jacobian** (or Frechet derivative matrix) $\mathbf{J}$ relates
    perturbations in the velocity model to changes in travel times:

    $$\delta \mathbf{t} = \mathbf{J} \, \delta \mathbf{v}$$

    Each row $J_{i,:}$ contains $\partial t_i / \partial v_j$ — the sensitivity
    of arrival $i$ to velocity parameter $j$. When reshaped onto the model
    grid, this gives the **sensitivity kernel** (or "fat ray") for that
    source–receiver pair.

    This notebook extracts the Jacobian using pyswmp's file-free API and
    visualises its structure.
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
    ## Setup: model, sources, receivers
    """)
    return


@app.cell
def _(np, pyswmp):
    model = pyswmp.create_constant_velocity_model(
        nx=50, ny=40,
        x0=110.0, y0=-45.0,
        dx=1.0, dy=1.0,
        velocity=3.5,
    )

    sources = pyswmp.Sources(
        positions=np.array([
            [135.0, -25.0],
            [125.0, -35.0],
        ]),
        source_type=1,
    )

    receivers = pyswmp.Receivers(
        positions=np.array([
            [120.0, -30.0],
            [130.0, -40.0],
            [145.0, -20.0],
            [150.0, -35.0],
        ]),
    )

    print(model)
    print(sources)
    print(receivers)
    return model, receivers, sources


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Forward model with Frechet extraction

    Enable `extract_frechet=True` and `extract_raypaths=True` in `TrackerOptions`.
    The Jacobian is returned as COO sparse arrays on `WaveFrontTrackerResult`.
    """)
    return


@app.cell
def _(model, pyswmp, receivers, sources):
    opts = pyswmp.TrackerOptions(
        coordinate_system=2,
        dt=5.0,
        max_iterations=500,
        extract_raypaths=True,
        extract_frechet=True,
    )

    tracker = pyswmp.WaveFrontTracker(model, opts)
    result = tracker.forward(sources, receivers)

    print(f"Arrivals:        {len(result.travel_times)}")
    print(f"Raypaths:        {len(result.raypaths)}")
    print(f"Jacobian shape:  {result.jacobian_shape}")
    print(f"Jacobian nnz:    {len(result.jacobian_vals)}")
    return (result,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Convert to scipy sparse matrix

    `result.to_sparse_jacobian()` returns a `scipy.sparse.coo_matrix` which can
    be converted to CSR/CSC for efficient arithmetic.
    """)
    return


@app.cell
def _(result):
    J = result.to_sparse_jacobian()

    print(f"Shape:    {J.shape}")
    print(f"NNZ:      {J.nnz}")
    print(f"Density:  {J.nnz / (J.shape[0] * J.shape[1]):.6f}")
    return (J,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Visualise Jacobian structure

    The sparsity pattern shows which model parameters each arrival is
    sensitive to — entries along the ray path from source to receiver.
    """)
    return


@app.cell
def _(J, plt):
    _fig, (_ax1, _ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Sparsity pattern
    _ax1.spy(J.tocsr(), markersize=1, aspect="auto")
    _ax1.set_xlabel("Model parameter index")
    _ax1.set_ylabel("Arrival index")
    _ax1.set_title("Jacobian sparsity pattern")

    # Row sums = total sensitivity per arrival
    _row_sums = abs(J.tocsr()).sum(axis=1).A1
    _ax2.barh(range(len(_row_sums)), _row_sums, color="steelblue", edgecolor="black")
    _ax2.set_xlabel("Sum |J| along row")
    _ax2.set_ylabel("Arrival index")
    _ax2.set_title("Total sensitivity per arrival")
    _ax2.invert_yaxis()

    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Sensitivity kernel maps

    Reshape a Jacobian row back onto the model grid to see the spatial
    sensitivity pattern — the "fat ray" or Frechet kernel for that arrival.
    """)
    return


@app.cell
def _(J, model, np, plt, result):
    _n_arrivals = J.shape[0]
    _n_to_show = min(_n_arrivals, 4)

    _fig, _axes = plt.subplots(1, _n_to_show, figsize=(4 * _n_to_show, 4), squeeze=False)

    # Jacobian columns include cushion nodes: (nx + 2*cn) x (ny + 2*cn)
    _cn = model.cushion_nodes
    _full_nx = model.nx + 2 * _cn
    _full_ny = model.ny + 2 * _cn

    for _k in range(_n_to_show):
        _ax = _axes[0, _k]
        _row = J.getrow(_k).toarray().ravel()

        # Reshape to full Fortran grid (including cushion), then trim
        _kernel_full = _row.reshape((_full_nx, _full_ny))
        _kernel = _kernel_full[_cn:_cn + model.nx, _cn:_cn + model.ny]

        _extent = [model.x0, model.x0 + model.dx * (model.nx - 1),
                   model.y0, model.y0 + model.dy * (model.ny - 1)]
        _im = _ax.imshow(
            _kernel.T, origin="lower", extent=_extent, aspect="auto",
            cmap="RdBu_r", vmin=-np.abs(_kernel).max(), vmax=np.abs(_kernel).max(),
        )
        plt.colorbar(_im, ax=_ax, shrink=0.8)

        _sid = result.source_ids[_k]
        _rid = result.receiver_ids[_k]
        _ax.set_title(f"src={_sid} -> rec={_rid}", fontsize=10)
        _ax.set_xlabel("Longitude")
        _ax.set_ylabel("Latitude")

    _fig.suptitle("Frechet sensitivity kernels", fontweight="bold", y=1.02)
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Raypath overlaid on sensitivity kernel

    The ray path sits at the centre of the sensitivity kernel — the Frechet
    derivative quantifies how much each grid cell contributes to the travel
    time along that ray.
    """)
    return


@app.cell
def _(J, model, np, plt, result):
    _fig, _ax = plt.subplots(figsize=(8, 6))

    # First arrival kernel — trim cushion nodes
    _cn = model.cushion_nodes
    _full_nx = model.nx + 2 * _cn
    _full_ny = model.ny + 2 * _cn
    _row = J.getrow(0).toarray().ravel()
    _kernel = _row.reshape((_full_nx, _full_ny))[_cn:_cn + model.nx, _cn:_cn + model.ny]
    _extent = [model.x0, model.x0 + model.dx * (model.nx - 1),
               model.y0, model.y0 + model.dy * (model.ny - 1)]
    _vmax = np.abs(_kernel).max()
    _im = _ax.imshow(_kernel.T, origin="lower", extent=_extent, aspect="auto",
                     cmap="RdBu_r", vmin=-_vmax, vmax=_vmax, alpha=0.8)
    plt.colorbar(_im, ax=_ax, label="dT/dV sensitivity")

    # Overlay ray path
    if result.raypaths:
        _path = result.raypaths[0]["path"]
        _ax.plot(_path[:, 0], _path[:, 1], "k-", lw=2, label="Ray path")
        _ax.plot(_path[0, 0], _path[0, 1], "r*", ms=15, zorder=5, label="Source")
        _ax.plot(_path[-1, 0], _path[-1, 1], "bv", ms=10, zorder=5, label="Receiver")

    _ax.set(xlabel="Longitude", ylabel="Latitude",
            title=f"Frechet kernel + ray path (arrival 0, t={result.travel_times[0]:.1f}s)")
    _ax.legend()
    _ax.grid(True, alpha=0.2)
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Summary

    - **`extract_frechet=True`** in `TrackerOptions` enables in-memory Jacobian computation
    - The Jacobian is stored as COO arrays on `WaveFrontTrackerResult` and converted
      to scipy sparse via `to_sparse_jacobian()`
    - Each row corresponds to one arrival (source-receiver-arrival triple)
    - Each column corresponds to one velocity model parameter (grid node)
    - Reshaping a row to the model grid shape gives the spatial sensitivity kernel
    """)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
