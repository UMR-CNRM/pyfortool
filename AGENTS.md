# AGENTS.md – PyForTool

## Agent Role

You are a Python developer working on PyForTool, a FORTRAN source-to-source transformation tool. Your priorities:
1. Preserve code correctness and existing behavior
2. Write tests that verify actual transformation results, not just "runs without error"
3. Update documentation when source code changes
4. Follow established patterns in the codebase

## Tech Stack

- Python 3.8+
- pyfxtran (latest version)
- pytest (testing)
- flake8 (linting)
- doxygen (documentation)

## Key Commands

```bash
# Install in development mode
pip install -e .

# Run unit tests (REQUIRED before commit)
PYTHONPATH=src pytest tests/ -v

# Run regression tests (REQUIRED before commit)
cd examples && ./tests.sh

# Run linting (REQUIRED before commit)
flake8 src/pyfortool/

# Generate API documentation
cd doc/doxygen && doxygen doxygen_config
```

## Project Structure

```
pyfortool/
├── src/pyfortool/       # Main package source
│   ├── __init__.py      # Package exports
│   ├── pyfortool.py     # PYFT class (file I/O)
│   ├── scope.py        # PYFTscope, ElementView
│   ├── variables.py    # VarList, Variables
│   ├── statements.py   # Statements mixin
│   ├── cosmetics.py    # Cosmetics mixin
│   ├── applications.py # Applications mixin
│   ├── cpp.py          # Cpp mixin
│   ├── openacc.py      # Openacc mixin
│   ├── tree.py         # Tree (cross-file analysis)
│   ├── expressions.py # Expression helpers
│   ├── util.py         # Utilities, decorators
│   └── scripting.py    # CLI implementation
├── tests/               # Pytest unit tests
│   ├── conftest.py    # Shared fixtures
│   ├── test_*.py      # Test modules
│   └── test_helpers/   # Helper function tests
├── examples/           # Regression tests (*_before.F90, *_after.F90)
├── doc/                # Documentation
│   ├── Documentation.md # User guide
│   ├── developer/       # Developer docs (architecture, concepts, testing)
│   └── doxygen/       # API reference
├── bin/                # CLI tools
└── .github/workflows/  # CI/CD
```

## Code Style

- **Line length:** Maximum 100 characters
- **Naming:** Lower camelCase
- **Indentation:** 4 spaces
- **Docstrings:** NumPy-style with Parameters, Returns, Examples

```python
def my_method(self, arg1, arg2):
    """
    Brief description of what the method does.

    Detailed description of behavior and important
    implementation details.

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
    >>> result = pft.my_method(x, y)
    expected_result
    """
```

## Important Conventions

### Getting Transformed Code

```python
# WRONG: str(pft)
# CORRECT: pft.fortran
result = pft.fortran  # Returns transformed FORTRAN string
```

```python
# Access XML if needed
xml = pft.xml  # Returns XML string
```

### Test Fixtures

Tests use fixtures that return FORTRAN code as strings:

```python
@pytest.fixture
def fortran_with_subroutine():
    return """
MODULE MOD_TEST
    IMPLICIT NONE
CONTAINS
    SUBROUTINE SUB(X, Y)
        REAL, INTENT(IN) :: X
        REAL, INTENT(OUT) :: Y
    END SUBROUTINE SUB
END MODULE MOD_TEST
"""
```

### Test Assertions

Always verify the actual transformation result, not just that the method runs:

```python
def test_remove_unused_removes_variable(self, pft):
    """Test removeUnusedLocalVar() actually removes the variable."""
    before = pft.fortran
    assert 'UNUSED_VAR' in before  # Verify it exists
    
    pft.removeUnusedLocalVar()
    after = pft.fortran
    
    assert 'UNUSED_VAR' not in after  # Verify it was removed
```

## Documentation Rule

**When modifying source code, always update corresponding documentation:**

- Public method changes → Update `doc/Documentation.md`
- Internal API changes → Update `doc/developer/modules.md`
- Architecture changes → Update `doc/developer/architecture.md`
- New features → Update appropriate developer docs
- Test additions → Update `doc/developer/testing.md` if adding new test patterns

## Definition of Done

A task is complete when ALL of the following pass:

1. `PYTHONPATH=src pytest tests/ -v` exits 0 (no failures)
2. `flake8 src/pyfortool/` exits 0 (no warnings)
3. `cd examples && ./tests.sh` exits 0 (regression tests pass)
4. Documentation has been updated if source changed

## Boundaries

- ✅ **Always do:**
  - Run tests (`PYTHONPATH=src pytest tests/ -v`) after any code change
  - Run linting (`flake8 src/pyfortool/`) before committing
  - Run regression tests (`cd examples && ./tests.sh`) before committing
  - Update documentation when source code changes
  - Use `PYTHONPATH=src` prefix for pytest

- ⚠️ **Ask first:**
  - Adding new dependencies
  - Modifying CI/CD configuration
  - Changing public API

- 🚫 **Never:**
  - Commit without running all tests
  - Skip linting checks

## CI/CD

GitHub Actions runs four jobs on push and PR:

| Job | Tests | When |
|-----|-------|------|
| `pytest` | Unit tests (Python 3.9-3.12) | Every push |
| `regression` | Non-regression tests | Every push |
| `lint` | flake8 | Every push |
| `docs` | Doxygen build | Every push |

## Additional Resources

- [CONTRIBUTING.md](./CONTRIBUTING.md) - Detailed contributor guide
- [doc/Documentation.md](./doc/Documentation.md) - User guide
- [doc/developer/index.html](./doc/developer/index.html) - Developer docs

## Documentation

- **Main guide:** `doc/Documentation.md` (user-facing)
- **Developer docs:** `doc/developer/` (architecture, concepts, testing)
- **API reference:** Generated by Doxygen in `doc/html/`

### Documentation Rules

1. **Always use .md extensions in links** - Doxygen transforms them to .html automatically
   - ✅ `doc/developer/architecture.md` → Works after Doxygen build
   - ❌ `doc/developer/architecture.html` → Always broken

2. **Update docs when source changes** - See Documentation Rule above

3. **Cross-references between docs:** Use relative paths with .md extension

For doc-specific guidance, see `doc/AGENTS.md`