# Surface Wave Multipathing (SWMP)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

A Python-wrapped Fortran library for surface wave multipathing ray tracing, based on the methodology published in Hauser et al. (2008).

## Features

- 🚀 **High-performance Fortran core** with modern Python interface
- 🔄 **Parallel execution** via `ProcessPoolExecutor` and `schwimmbad` (MPI support)
- 💾 **File-free API** for in-memory workflows (no temporary files required)
- 🔧 **Integration ready** for CoFI/ESPRESSO inverse problem frameworks
- 📊 **Built-in visualization** with cartopy and matplotlib
- 🌍 **Multi-coordinate systems** supporting both Cartesian and spherical geometries
- ⚡ **Wavefront tracking** with multi-arrival detection
- 📐 **Frechet derivatives** for inverse problems

## Recent Improvements

This version includes several modernizations over the original implementation:

- **CMake-based build system** for cross-platform compilation (Linux, macOS, Windows)
- **CoFI/ESPRESSO integration** for inverse problem solving
- **Modern Python packaging** with PEP 517/518 compliance
- **File-free operation** allowing complete in-memory workflows
- **Migration to cartopy** from deprecated basemap
- **Parallel processing** with multiple backend support

## Quick Start

### Installation

#### From source (recommended for development):

```bash
# Requires: Python 3.9+, gfortran, CMake 3.15+
git clone https://github.com/yourusername/swmp.git
cd swmp
pip install -e .
```

#### Production installation:

```bash
pip install pyswmp
```

### Basic Usage

```python
import pyswmp
import numpy as np

# Create a simple velocity model
model = pyswmp.create_constant_velocity_model(
    nx=100, ny=80,
    x0=110.0, y0=-45.0,  # Lower-left corner (degrees)
    dx=0.5, dy=0.5,       # Grid spacing (degrees)
    velocity=3.5          # km/s
)

# Define source and receiver positions
sources = pyswmp.Sources(
    positions=np.array([[135.0, -25.0]]),  # [lon, lat]
    source_type=1  # Point source
)

receivers = pyswmp.Receivers(
    positions=np.array([[120.0, -30.0], [140.0, -20.0]])  # [lon, lat]
)

# Configure ray tracing options
options = pyswmp.TrackerOptions(
    coordinate_system=2,  # Spherical coordinates
    max_arrivals=10       # Maximum arrivals per receiver
)

# Run wavefront tracking
tracker = pyswmp.WaveFrontTracker(
    model=model,
    sources=sources,
    receivers=receivers,
    options=options
)

result = tracker.forward()

# Access results
print(f"Computed {len(result.travel_times)} arrivals")
print(f"Travel times: {result.travel_times}")
print(f"Azimuths: {result.azimuths}")

# Convert to pandas DataFrame
df = result.to_dataframe()
print(df.head())
```

### Parallel Execution

```python
from concurrent.futures import ProcessPoolExecutor

# Run with parallel processing
with ProcessPoolExecutor(max_workers=4) as pool:
    result = tracker.forward(pool=pool)
```

For MPI support:

```python
from schwimmbad import MPIPool

with MPIPool() as pool:
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    result = tracker.forward(pool=pool)
```

## Documentation

- **[API Guide](PROGRAMMATIC_API.md)** - Complete programmatic API reference
- **[File-Free API Status](FILE_FREE_API_STATUS.md)** - Implementation status and migration guide
- **[Examples](examples/)** - Working examples and Jupyter notebooks
  - `basic_usage.py` - Simple wavefront tracking
  - `file_free_example.py` - In-memory operation
  - `parallel_forward.py` - Parallel execution
  - `programmatic_example.py` - Complete workflow
- **[Scientific Paper](doc/Hauser_et_al._2008.pdf)** - Original methodology

## Examples

The `examples/` directory contains several demonstrations:

- **ausvel/** - Real-world surface wave propagation across Australia
- **greatcircle/** - Analytical validation with great circle paths

See [examples/README.md](examples/README.md) for detailed descriptions.

## Requirements

### Build Requirements
- Python 3.9 or higher
- CMake 3.15 or higher
- Fortran compiler (gfortran recommended)

### Runtime Dependencies
- numpy >= 1.21
- scipy >= 1.7
- matplotlib >= 3.4
- cartopy >= 0.20

### Optional Dependencies
- schwimmbad >= 0.3 (for parallel execution)
- pandas (for DataFrame conversion)
- pytest >= 7.0 (for development)

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Citation

If you use this software in your research, please cite:

```bibtex
@article{hauser2008multiarrival,
  title={Multiarrival wavefront tracking and its applications},
  author={Hauser, Juerg and Sambridge, Malcolm and Rawlinson, Nicholas},
  journal={Geochemistry, Geophysics, Geosystems},
  volume={9},
  number={11},
  year={2008},
  publisher={Wiley Online Library},
  doi={10.1029/2008GC002069}
}
```

## License

This project is licensed under the GNU General Public License v3.0 or later - see [LICENSE.txt](LICENSE.txt) for details.

## Authors

- **Juerg Hauser** - Original implementation and Python wrappers
  - CSIRO Mineral Resources
  - Email: juerg.hauser@csiro.au

## Acknowledgments

Original methodology developed by J. Hauser, M. Sambridge, and N. Rawlinson.

## References

Hauser, J., Sambridge, M. and Rawlinson, N. (2008). Multiarrival wavefront tracking and its applications. *Geochemistry, Geophysics, Geosystems*, 9(11), Q11001. https://doi.org/10.1029/2008GC002069
