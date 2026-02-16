# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**pyswmp** is a Python-wrapped Fortran library for surface wave multipathing ray tracing on spherical and Cartesian coordinate systems (Hauser et al., 2008). Fortran core compiled via CMake, exposed to Python via ctypes, packaged with scikit-build-core.

## Build & Install

Requires: Python 3.9+, CMake 3.15+, gfortran (or Intel/NVIDIA Fortran compiler).

```bash
pip install -e ".[dev]"          # editable install with dev extras
pre-commit install               # set up pre-commit hooks
```

Manual Fortran build (already done in `build/`):
```bash
cd build && cmake .. && make
```

## Testing

```bash
pytest                                    # all tests (coverage enabled by default via pyproject.toml)
pytest tests/test_data.py                 # single file
pytest tests/test_data.py::TestVelocityModel2D::test_valid_model_creation  # single test
pytest -m slow                            # integration tests only (require compiled Fortran lib)
```

Coverage is configured in pyproject.toml: `--cov=pyswmp --cov-report=term-missing --cov-report=html`.

## Linting & Formatting

```bash
pre-commit run --all-files    # run all checks
black --check src/ tests/     # formatting (line length 100)
isort --check src/ tests/     # import sorting (profile=black)
ruff check src/ tests/        # linting
mypy src/                     # type checking (non-strict, non-blocking in CI)
```

## Architecture

### Fortran Core (`swmp/`)

Modules compiled into `libswmp` shared library with ISO C bindings:

- **m_types** — Data types (scalar_field_2d, sources, receivers, arrival, wavefront, node linked-list) and constants
- **m_error** — Error codes (SWMP_SUCCESS=0, etc.)
- **m_functions** — Math utilities (cubic B-spline interpolation, Gaussian noise)
- **m_solver** — ODE solvers (RK4, RK5, adaptive RK5) for bicharacteristic strip evolution
- **m_wavefronts** — Wavefront tracking via linked-list node insertion/removal
- **m_arrivals** — Arrival extraction from wavefronts at receiver locations
- **m_linalg** — Linear algebra (subspace inversion, SVD)
- **m_inverse** — Fréchet matrices, model updates for inverse problems
- **m_swmp** — Main C-callable API: contains `modgen`, `pred2obs`, `rat` modules with `forward()`, `forward_single_source()`, getters/setters, and file-free API entry points

**Critical detail:** The `rat` module uses module-level state variables (vmod, recs, sous, conf, wafc, jac). Each worker process must load its own independent `libswmp` instance — `ThreadPoolExecutor` is explicitly unsupported.

### Python Layer (`src/pyswmp/`)

Layered architecture from low-level to high-level:

1. **`_libswmp.py`** — ctypes wrapper defining all C function signatures. `LibSWMP` class. Library discovery searches package dir then `build/` dir.
2. **`_pyswmp.py`** — Legacy ctypes wrapper (backward compat). Has plotting via cartopy.
3. **`data.py`** — Dataclasses: `VelocityModel2D`, `Sources`, `Receivers` with validation in `__post_init__()`. File I/O helpers.
4. **`options.py`** — `TrackerOptions` dataclass with all ray tracing parameters.
5. **`results.py`** — `WaveFrontTrackerResult` with `__add__()` for merging per-source results via `functools.reduce()`.
6. **`solvers.py`** — `WaveFrontTracker`: supports file-based (config_file) or file-free (model+sources+receivers) init. `forward()` runs sequentially or parallel via `ProcessPoolExecutor`/schwimmbad.

### Key Design Patterns

- **File-free API** (v0.2.0): Complete in-memory workflow without disk I/O — `set_velocity_model_from_memory`, `set_sources_from_memory`, `set_receivers_from_memory`
- **Parallel execution**: Worker function `_calc_single_source_worker()` creates independent Fortran library instances per process
- **Dual initialization**: `WaveFrontTracker` accepts either config files or programmatic data objects

## Conventions

- Commit messages: Conventional Commits (`feat(scope): subject`)
- The `sundry/` directory contains archived legacy code and is excluded from all linting/formatting
- Wavefront extraction from memory is NOT yet implemented in Fortran (TODO in `_libswmp.py`)
