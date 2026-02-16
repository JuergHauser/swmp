# Contributing to SWMP

Thank you for your interest in contributing to SWMP (Surface Wave Multipathing)! This document provides guidelines for contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Making Changes](#making-changes)
- [Testing](#testing)
- [Code Style](#code-style)
- [Submitting Changes](#submitting-changes)
- [Reporting Issues](#reporting-issues)

## Code of Conduct

This project follows standard open-source community guidelines. Please be respectful and constructive in all interactions.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/yourusername/swmp.git
   cd swmp
   ```
3. **Add the upstream repository**:
   ```bash
   git remote add upstream https://github.com/originalowner/swmp.git
   ```

## Development Setup

### Prerequisites

- Python 3.9 or higher
- CMake 3.15 or higher
- Fortran compiler (gfortran, ifort, or NVIDIA HPC SDK)
- Git

### Installation for Development

1. **Create a virtual environment**:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

2. **Install the package in editable mode with dev dependencies**:
   ```bash
   pip install -e ".[dev]"
   ```

3. **Install pre-commit hooks**:
   ```bash
   pre-commit install
   ```

   This will automatically run code quality checks before each commit.

### Verifying Installation

Run the test suite to ensure everything is working:

```bash
pytest
```

## Making Changes

### Workflow

1. **Create a new branch** for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b fix/issue-number-description
   ```

2. **Make your changes** following the code style guidelines below

3. **Test your changes** thoroughly:
   ```bash
   pytest
   ```

4. **Commit your changes** with clear, descriptive messages:
   ```bash
   git add .
   git commit -m "Add feature: brief description"
   ```

   The pre-commit hooks will automatically format your code.

### Commit Message Guidelines

Follow the conventional commits format:

```
<type>(<scope>): <subject>

<body>

<footer>
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting, etc.)
- `refactor`: Code refactoring
- `test`: Adding or updating tests
- `chore`: Maintenance tasks

**Examples:**
```
feat(solver): add support for adaptive time stepping

Implements adaptive RK5 solver with error control.

Closes #42
```

```
fix(data): validate file sizes before writing

Prevents memory issues with very large models.
```

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=pyswmp --cov-report=html

# Run specific test file
pytest tests/test_data.py

# Run specific test
pytest tests/test_data.py::TestVelocityModel2D::test_valid_model_creation
```

### Writing Tests

- Place tests in the `tests/` directory
- Name test files as `test_*.py`
- Name test functions as `test_*`
- Use pytest fixtures from `tests/conftest.py`
- Aim for high code coverage (>80%)

**Example test:**

```python
def test_velocity_model_validation():
    """Test that VelocityModel2D validates inputs correctly."""
    from pyswmp import VelocityModel2D
    import numpy as np

    # Valid model should work
    model = VelocityModel2D(
        velocities=np.ones((10, 10)),
        x0=0.0, y0=0.0, dx=1.0, dy=1.0
    )
    assert model.nx == 10

    # Invalid spacing should raise error
    with pytest.raises(ValueError):
        VelocityModel2D(
            velocities=np.ones((10, 10)),
            x0=0.0, y0=0.0, dx=-1.0, dy=1.0
        )
```

## Code Style

### Python

We follow PEP 8 with some modifications:

- **Line length**: 100 characters
- **Formatter**: Black
- **Import sorting**: isort
- **Linter**: Ruff
- **Type hints**: Encouraged but not required

The pre-commit hooks will automatically format your code. You can also run manually:

```bash
# Format code
black src/ tests/ examples/

# Sort imports
isort src/ tests/ examples/

# Lint
ruff check src/ tests/ examples/ --fix
```

### Python Style Guidelines

- Use descriptive variable names
- Add docstrings to all public functions/classes
- Include type hints where helpful
- Keep functions focused and small
- Use dataclasses for data structures

**Example:**

```python
def create_velocity_model(
    nx: int,
    ny: int,
    velocity: float
) -> VelocityModel2D:
    """Create a constant velocity model.

    Args:
        nx: Number of grid points in X direction
        ny: Number of grid points in Y direction
        velocity: Constant velocity value (km/s)

    Returns:
        VelocityModel2D instance

    Example:
        >>> model = create_velocity_model(100, 80, 3.5)
        >>> print(model.nx, model.ny)
        100 80
    """
    velocities = np.full((nx, ny), velocity, dtype=np.float32)
    return VelocityModel2D(velocities, 0, 0, 1, 1)
```

### Fortran

- Use Fortran 90+ features
- Keep subroutines focused
- Add comments for complex algorithms
- Use meaningful variable names
- Maintain ISO C bindings for Python interface

### CMake

- Use lowercase for commands
- Add comments for non-obvious logic
- Set compiler flags appropriately for each compiler

## Submitting Changes

### Before Submitting

1. **Ensure all tests pass**:
   ```bash
   pytest
   ```

2. **Check code quality**:
   ```bash
   pre-commit run --all-files
   ```

3. **Update documentation** if needed

4. **Add tests** for new functionality

### Creating a Pull Request

1. **Push your branch** to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

2. **Create a Pull Request** on GitHub:
   - Use a clear, descriptive title
   - Reference any related issues
   - Describe what changed and why
   - Include examples if applicable

3. **Pull Request Template**:
   ```markdown
   ## Description
   Brief description of changes

   ## Type of Change
   - [ ] Bug fix
   - [ ] New feature
   - [ ] Documentation update
   - [ ] Refactoring

   ## Testing
   - [ ] All existing tests pass
   - [ ] Added new tests for changes
   - [ ] Tested manually with examples

   ## Checklist
   - [ ] Code follows project style guidelines
   - [ ] Added docstrings for new functions
   - [ ] Updated relevant documentation
   - [ ] No new warnings from linters

   ## Related Issues
   Closes #42
   ```

4. **Address review feedback** promptly

## Reporting Issues

### Bug Reports

When reporting bugs, please include:

- **System information**: OS, Python version, compiler
- **Steps to reproduce**: Minimal example
- **Expected behavior**: What should happen
- **Actual behavior**: What actually happens
- **Error messages**: Full traceback if applicable

**Template:**

```markdown
**Environment:**
- OS: Ubuntu 22.04
- Python: 3.11
- Compiler: gfortran 12.2

**Description:**
Brief description of the issue

**Steps to Reproduce:**
1. Step one
2. Step two
3. ...

**Expected Behavior:**
What should happen

**Actual Behavior:**
What actually happens

**Error Message:**
```
Full error traceback here
```

**Additional Context:**
Any other relevant information
```

### Feature Requests

When requesting features:

- Describe the use case
- Explain why this would be useful
- Provide examples if possible
- Suggest implementation approach if you have ideas

## Questions?

If you have questions about contributing:

- Check existing issues and discussions
- Open a new issue with the "question" label
- Contact the maintainers: juerg.hauser@csiro.au

## License

By contributing, you agree that your contributions will be licensed under the same license as the project (GPL-3.0-or-later).

---

Thank you for contributing to SWMP!
