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
cd examples && ./tests.sh

# Run specific test (by filename prefix)
cd examples && ./tests.sh test_name

# Update reference files (after manual verification)
cd examples && ./tests.sh --update
```

**How it works:**
1. Copies all `*_before.F90` files to a temp directory
2. Runs pyfortool with the transformation command (from the `!#PYFT transfo:` comment)
3. Compares output with corresponding `*_after.F90` files

## Test Structure

```
tests/
├── conftest.py                    # Shared pytest fixtures
│   ├── fortran_simple             # Simple program
│   ├── fortran_with_subroutine    # Module with subroutine
│   ├── fortran_with_calls         # Module with CALL statements
│   ├── fortran_with_arrays        # Module with array syntax
│   ├── fortran_with_if            # Module with IF statements
│   ├── fortran_with_loops         # Module with DO loops
│   ├── TempFortranFile            # Context manager
│   └── pft_simple, pft_module, pft_calls, pft_arrays, pft_if, pft_loops  # Fixtures
│
├── test_pyft.py                   # PYFT class (17 tests)
│   ├── TestPYFTInit               # __init__, file loading
│   ├── TestPYFTProperties         # fortran, xml properties
│   ├── TestPYFTFileOperations     # write, writeXML, rename
│   ├── TestPYFTGetFileName        # getFileName method
│   ├── TestPYFTInheritance        # Conservative PYFT
│   └── TestPYFTMultipleScopes     # Multiple scopes
│
├── test_scope.py                  # PYFTscope class (22 tests)
│   ├── TestGetScopes              # getScopes() method
│   ├── TestGetScopeNode           # getScopeNode() method
│   ├── TestScopeProperties        # path, mainScope, parentScope
│   ├── TestGetParent              # Parent element navigation
│   ├── TestGetSiblings            # Sibling navigation
│   ├── TestNormalizeScope         # Scope path normalization
│   ├── TestGetScopePath           # getScopePath method
│   ├── TestGetParentScopeNode
│   ├── TestIsScopeNode
│   └── TestShowScopesList
│
├── test_varList.py                # VarList class (12 tests)
│   ├── TestVarListFindVar         # findVar() with array, exactScope
│   ├── TestVarListRestrict        # restrict() method
│   ├── TestVarListProperty        # __len__, __getitem__
│   └── TestVarListShowVarList
│
├── test_variables.py              # Variables mixin (27 tests)
│   ├── TestAttachArraySpecToEntity# DIMENSION to entity
│   ├── TestCheckImplicitNone      # IMPLICIT NONE check
│   ├── TestCheckIntent            # INTENT attribute check
│   ├── TestCheckOnly              # ONLY clause check
│   ├── TestCheckUnusedLocalVar    # Unused variable check
│   ├── TestShowUnusedVar          # Display unused vars
│   ├── TestAddExplicitArrayBounds
│   ├── TestAddArrayParentheses
│   ├── TestModifyAutomaticArrays
│   ├── TestRemoveUnusedLocalVar   # Actually removes variables
│   ├── TestRemoveVarIfUnused
│   ├── TestAddVar
│   ├── TestRemoveVar
│   ├── TestRenameVar
│   └── TestCheckKeyDimConsistency    # MERGE-based dims check
│
├── test_statements.py             # Statements mixin (14 tests)
│   ├── TestRemoveCall             # Remove CALL statements
│   ├── TestRemovePrints           # Remove PRINT statements
│   ├── TestRemoveArraySyntax      # Array to DO loops
│   ├── TestSetFalseIfStmt         # IF flags to .FALSE.
│   ├── TestChangeIfStatementsInIfConstructs
│   ├── TestEmpty
│   ├── TestIsNodeInProcedure
│   └── TestIsNodeInCall
│
├── test_cosmetics.py              # Cosmetics mixin (23 tests)
│   ├── TestCosmeticsCase          # upperCase, lowerCase
│   ├── TestCosmeticsIndent        # indent method
│   ├── TestCosmeticsComments      # removeComments
│   ├── TestCosmeticsEmptyLines
│   ├── TestCosmeticsContinuation
│   ├── TestCosmeticsSpaces
│   ├── TestCosmeticsIfStatements
│   ├── TestCosmeticsContains
│   └── TestCosmeticsFormatModuleUse  # formatModuleUse
│
├── test_cpp.py                    # Cpp mixin (12 tests)
│   ├── TestCppApplyIfdef          # #ifdef evaluation
│   ├── TestCppIfndef              # #ifndef handling
│   ├── TestCppNested              # Nested #ifdef
│   ├── TestCppElse                # #else handling
│   └── TestCppMultipleBlocks
│
├── test_openacc.py                # Openacc mixin (6 tests)
│   ├── TestOpenaccRemoveACC
│   ├── TestOpenaccRemoveBypass
│   ├── TestOpenaccAddACCData
│   ├── TestOpenaccAddRoutineSeq
│   ├── TestOpenaccAllocateHIP
│   └── TestOpenaccCrayBypass
│
├── test_applications.py           # Applications mixin (15 tests)
│   ├── TestApplicationsDrHook
│   ├── TestApplicationsDeleteCalls
│   ├── TestApplicationsBudgetDDH
│   ├── TestApplicationsMPPDB
│   ├── TestApplicationsPHYEX
│   ├── TestApplicationsConvertTypes
│   ├── TestApplicationsMesoNHGPU
│   ├── TestApplicationsStack
│   ├── TestApplicationsIJDim
│   ├── TestApplicationsModi
│   └── TestApplicationsRemoveExtraDO
│
└── test_helpers/
    ├── __init__.py
    ├── test_util.py               # Utility functions (26 tests)
    │   ├── TestTag
    │   ├── TestN2Name
    │   ├── TestAlltext
    │   ├── TestIsInt
    │   ├── TestIsFloat
    │   ├── TestNonCode
    │   ├── TestIsStmt
    │   ├── TestIsConstruct
    │   └── TestIsExecutable
    │
    └── test_expressions.py        # Expression helpers (22 tests)
        ├── TestCreateElem
        ├── TestCreateExpr
        ├── TestCreateExprPart
        └── TestSimplifyExpr
```

## Writing Tests

### Test Fixtures

Shared test data is defined in `conftest.py`. Each fixture returns FORTRAN code as a string:

```python
@pytest.fixture
def fortran_with_subroutine():
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

PYFT fixtures wrap the code and return a PYFT instance:

```python
@pytest.fixture
def pft_module(fortran_with_subroutine):
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(fortran_with_subroutine)
        return PYFT(fpath)
```

### Test Example

```python
class TestCosmeticsCase:
    """Tests for upperCase and lowerCase methods."""

    @pytest.fixture
    def fortran_lowercase(self):
        return """
module mod_test
    implicit none
contains
    subroutine sub_test(x, y)
        real, intent(in) :: x
        real, intent(out) :: y
    end subroutine sub_test
end module mod_test
"""

    @pytest.fixture
    def pft_lowercase(self, fortran_lowercase):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_lowercase)
            return PYFT(fpath)

    def test_uppercase_converts_keywords(self, pft_lowercase):
        """Test that upperCase() converts FORTRAN keywords."""
        pft_lowercase.upperCase()
        result = pft_lowercase.fortran
        assert 'MODULE' in result
        assert 'SUBROUTINE' in result
```

## Best Practices

1. **Use descriptive test names**: `test_remove_unused_removes_variable`
2. **Verify the result, not just execution**: Use `assert 'UNUSED_VAR' not in after`
3. **Use fixtures for shared data**: Reduces duplication
4. **Test actual transformation**: Check `pft.fortran` contains expected changes
5. **Skip complex tests gracefully**: Use `pytest.skip()` when appropriate

```python
def test_complex_feature(self):
    try:
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

## See Also

- Also in the Developer's Guide
  - [Architecture Guide](architecture.md) - Class hierarchy and data flow
  - [Core Concepts](concepts.md) - Detailed concept explanations
  - [Module Organization](modules.md) - What each module does
  - [CONTRIBUTING.md](../../CONTRIBUTING.md) - Contributing guidelines
- [User's Guide](../Documentation.md) - End-user documentation
- [API Reference](../html/index.html) - Auto-generated from docstrings
