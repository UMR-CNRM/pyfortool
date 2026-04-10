# Testing Guide

This guide explains how to run and write tests for PyForTool.

## Test Suites

PyForTool has two test suites:

| Suite | Location | Framework | Purpose |
|-------|----------|----------|---------|
| Unit Tests | `tests/` | pytest | Test individual methods and classes |
| Regression Tests | `examples/` | Bash scripts | End-to-end transformation validation |

## Running Tests

### Unit Tests (pytest)

```bash
# Install pytest if not already installed
pip install pytest pytest-cov

# Run all unit tests
PYTHONPATH=src pytest tests/ -v

# Run specific test file
PYTHONPATH=src pytest tests/test_pyft.py -v

# Run specific test class
PYTHONPATH=src pytest tests/test_pyft.py::TestPYFTInit -v

# Run specific test
PYTHONPATH=src pytest tests/test_pyft.py::TestPYFTInit::test_simple_init -v

# Run with coverage
PYTHONPATH=src pytest tests/ --cov=pyfortool --cov-report=html

# Generate coverage report
PYTHONPATH=src pytest tests/ --cov=pyfortool --cov-report=term-missing
```

### Regression Tests

```bash
# Run all regression tests
./examples/tests.sh

# Run specific test
./examples/tests.sh test_name

# Update reference files (after manual verification)
./examples/tests.sh --update
```

**How it works:**
1. Copies all `*_before.F90` files to a temp directory
2. Runs pyfortool with the transformation command (from the `!#PYFT transfo:` comment)
3. Compares output with corresponding `*_after.F90` files

## Test Structure

```
tests/
├── conftest.py              # Shared pytest fixtures
│   ├── FORTRAN code fixtures (simple_fortran, module_with_subroutine, etc.)
│   ├── TempFortranFile context manager
│   └── PYFT fixtures (pft_simple, pft_module, etc.)
├── test_pyft.py             # PYFT class tests
│   ├── TestPYFTInit         # Initialization tests
│   ├── TestPYFTProperties  # Property access tests
│   ├── TestPYFTFileOperations # I/O tests
│   └── TestPYFTMultipleScopes # Scope handling tests
├── test_scope.py             # PYFTscope tests
│   ├── TestGetScopes        # getScopes() tests
│   ├── TestGetScopeNode     # getScopeNode() tests
│   ├── TestScopeProperties  # path, mainScope, parentScope
│   └── TestGetParent        # Parent navigation tests
├── test_varList.py           # VarList tests
│   ├── TestVarListFindVar   # findVar() tests
│   └── TestVarListRestrict  # restrict() tests
├── test_variables.py         # Variables mixin tests
│   ├── TestAttachArraySpecToEntity
│   ├── TestCheckImplicitNone
│   └── TestModifyAutomaticArrays
├── test_statements.py         # Statements mixin tests
│   ├── TestRemoveCall
│   ├── TestRemovePrints
│   └── TestSetFalseIfStmt
├── test_cosmetics.py          # Cosmetics mixin tests
│   ├── TestCosmeticsCase
│   ├── TestCosmeticsIndent
│   └── TestCosmeticsComments
├── test_cpp.py               # Cpp mixin tests
│   └── TestCppApplyIfdef
├── test_openacc.py            # Openacc mixin tests
│   ├── TestOpenaccRemoveACC
│   └── TestOpenaccAddACCData
├── test_applications.py       # Applications mixin tests
│   ├── TestApplicationsDrHook
│   └── TestApplicationsStack
└── test_helpers/
    ├── test_util.py          # Utility function tests
    │   ├── TestTag
    │   ├── TestN2Name
    │   ├── TestIsInt
    │   └── TestNonCode
    └── test_expressions.py   # Expression helper tests
        ├── TestCreateElem
        ├── TestCreateExprPart
        └── TestSimplifyExpr
```

## Writing Tests

### Test Fixtures

Shared test data is defined in `conftest.py`:

```python
@pytest.fixture
def module_with_subroutine():
    """Module with a subroutine."""
    return """
MODULE MOD_TEST
    IMPLICIT NONE
CONTAINS
    SUBROUTINE SUB(X, Y, Z)
        REAL, INTENT(IN) :: X(:)
        REAL, INTENT(OUT) :: Z(:)
        INTEGER :: I
        DO I = 1, SIZE(X)
            Z(I) = X(I) * 2.0
        END DO
    END SUBROUTINE SUB
END MODULE MOD_TEST
"""
```

### PYFT Fixtures

```python
@pytest.fixture
def pft_module(module_with_subroutine):
    """PYFT instance from module with subroutine."""
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(module_with_subroutine)
        return PYFT(fpath)
```

### Test Example

```python
class TestCosmeticsIndent:
    """Tests for indent method."""

    def test_indent_default(self, pft_module):
        """Test default indentation."""
        pft_module.indent()
        result = pft_module.fortran
        assert result is not None

    def test_indent_custom_params(self, pft_module):
        """Test indentation with custom parameters."""
        pft_module.indent(indentProgramunit=0, indentBranch=4)
        result = pft_module.fortran
        assert result is not None
```

## Best Practices

1. **Use descriptive test names**: `test_add_drhook_instruments_subroutine`
2. **One assertion per test**: Easier to diagnose failures
3. **Use fixtures for shared data**: Reduces duplication
4. **Test the result, not the implementation**: Test `pft.fortran` output
5. **Skip complex tests gracefully**: Use `pytest.skip()` when appropriate

```python
def test_complex_feature(self):
    """Test that requires specific setup."""
    try:
        # Try the transformation
        pft.complexMethod()
        result = pft.fortran
        assert 'expected' in result
    except (AttributeError, TypeError):
        pytest.skip("Complex feature requires full tree setup")
```

## CI/CD

Tests run automatically via GitHub Actions (`.github/workflows/test.yml`):

| Job | Tests | Python Version |
|-----|-------|----------------|
| `pytest` | Unit tests (pytest) | 3.9, 3.10, 3.11, 3.12 |
| `regression` | Non-regression tests | 3.11 |
| `lint` | flake8 | 3.11 |
| `docs` | Doxygen build | 3.11 |

### GitHub Actions Workflow

```yaml
jobs:
  pytest:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.11'
      - name: Install dependencies
        run: pip install pytest pytest-cov -e .
      - name: Run tests
        run: PYTHONPATH=src pytest tests/ -v

  regression:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.11'
      - name: Install dependencies
        run: pip install -e .
      - name: Run regression tests
        run: cd examples && ./tests.sh
```

## See Also

- [Architecture Guide](md__home_sriette_GIT_pyfortool_doc_developer_architecture.html) - How PyForTool is structured
- [Module Organization](md__home_sriette_GIT_pyfortool_doc_developer_modules.html) - What each module does
- [CONTRIBUTING.md](../../CONTRIBUTING.html#testing) - Contributing guidelines
