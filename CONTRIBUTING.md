# Contributing to PyForTool

Thank you for your interest in contributing to PyForTool!

This guide covers how to set up a development environment, coding standards, testing procedures, and how to submit changes.

## Table of Contents

- [Development Setup](#development-setup)
- [Code Standards](#code-standards)
- [Testing](#testing)
- [Submitting Changes](#submitting-changes)
- [Adding New Methods](#adding-new-methods)

## Development Setup

### Prerequisites

- Python 3.8+
- [pyfxtran](https://github.com/SebastienRietteMTO/pyfxtran)
- Git

### Clone and Install

```bash
git clone https://github.com/your-fork/pyfortool.git
cd pyfortool
pip install -e .
pip install flake8 pylint
```

### Run Tests

```bash
# From project root
./examples/tests.sh

# Or run linting separately
flake8 src/pyfortool/ bin/pyfortool_*

pylint -d R0912,C0209,R0915,R1702,C0302,R0913,R0914,W1202,R0904,R0902 \
    src/pyfortool/ bin/pyfortool_*
```

### Generate Documentation

```bash
# Doxygen API reference
cd doc/doxygen
doxygen doxygen_config

# View generated docs
open html/index.html
```

## Code Standards

### Style

- **Line length:** Maximum 100 characters
- **Naming:** Lower camelCase
- **Indentation:** 4 spaces

### Linting

PyForTool uses flake8 and pylint for code quality:

```bash
# No warnings from flake8
flake8 src/pyfortool/

# Pylint score above 9.8
pylint -d R0912,C0209,R0915,R1702,C0302,R0913,R0914,W1202,R0904,R0902 \
    src/pyfortool/
```

### Comment Style

Use NumPy-style docstrings:

```python
def myMethod(self, arg1, arg2):
    """
    Brief description of what the method does.

    Detailed description of the method's behavior, including
    any important implementation details.

    Parameters
    ----------
    arg1 : type
        Description of arg1.
    arg2 : type
        Description of arg2.

    Returns
    -------
    type
        Description of return value.

    Examples
    --------
    >>> pft = PYFT('input.F90')
    >>> pft.myMethod(x, y)
    result
    """
```

## Testing

### Test Files

Tests are in the `examples/` directory:

- `*_before.F90` - Input files
- `*_after.F90` - Expected output files
- `tests.sh` - Test runner script

### Running Tests

```bash
# Run all tests
./examples/tests.sh

# Run specific test
./examples/tests.sh test_name
```

### Adding Tests

1. Create `test_name_before.F90` with input code
2. Transform manually to create `test_name_after.F90`
3. Add command comment at top of before file:

```fortran
! pyfortool test_name --upperCase
PROGRAM test_name
...
```

4. Run tests to verify:

```bash
./examples/tests.sh test_name
```

## Submitting Changes

### Branch Naming

```
feature/description
bugfix/description
docs/description
refactor/description
```

### Commit Messages

- Use clear, concise descriptions
- Start with verb: "Add", "Fix", "Update", "Remove"
- Reference issues: "Fixes #123"

### Pull Request Process

1. Fork the repository
2. Create a feature branch
3. Make changes following code standards
4. Add/update tests
5. Run linting and tests
6. Submit pull request with description
7. Address review feedback

## Adding New Methods

### Step 1: Choose the Right Module

| Module | Purpose |
|--------|---------|
| `variables.py` | Variable declarations |
| `statements.py` | Statement manipulation |
| `cosmetics.py` | Code formatting |
| `applications.py` | Model-specific transforms |
| `cpp.py` | Preprocessor directives |
| `openacc.py` | GPU directives |

### Step 2: Apply Decorators

```python
class Variables:
    @debugDecor
    @noParallel           # If modifies XML
    @updateVarList        # If changes variables
    @updateTree('signal') # If affects dependencies
    def myNewMethod(self, arg):
        """
        Brief description.

        Parameters
        ----------
        arg : type
            Description.

        Examples
        --------
        >>> pft = PYFT('input.F90')
        >>> pft.myNewMethod(x)
        """
        # Implementation
        pass
```

### Step 3: Documentation

- Add docstring with NumPy style
- Include Parameters, Returns, Examples sections
- Add CLI option in `scripting.py` if appropriate

### Step 4: Add Tests

See [Adding Tests](#adding-tests) for instructions on adding test files.

## Project Structure

```
pyfortool/
├── src/pyfortool/       # Main package
│   ├── pyfortool.py     # PYFT class
│   ├── scope.py         # PYFTscope, ElementView
│   ├── variables.py     # VarList, Variables
│   ├── statements.py    # Statements
│   ├── cosmetics.py     # Cosmetics
│   ├── applications.py  # Applications
│   ├── cpp.py           # Preprocessor
│   ├── openacc.py       # OpenACC
│   ├── tree.py          # Tree
│   ├── expressions.py   # Expression helpers
│   └── util.py          # Utilities
├── bin/                 # CLI tools
├── examples/            # Test files
├── doc/                 # Documentation
│   ├── Documentation.md  # User guide
│   ├── doxygen/         # API reference
│   └── developer/        # Developer docs
└── CONTRIBUTING.md      # This file
```

## Getting Help

- [Documentation](doc/Documentation.html)
- [API Reference](doc/html/index.html)
- [Developer Guide](doc/developer/index.html)
- [Issue Tracker](https://github.com/SebastienRietteMTO/pyfortool/issues)

## License

By contributing, you agree that your contributions will be licensed under the project's license.
