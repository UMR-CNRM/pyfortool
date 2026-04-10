"""
Tests for VarList class in pyfortool.variables.
"""

import pytest
from pyfortool import PYFT


@pytest.fixture
def pft_with_vars():
    """PYFT instance with various variables."""
    code = """
MODULE MOD_VAR
    IMPLICIT NONE
    REAL :: MODULE_VAR
    INTEGER, PARAMETER :: CONST_VAL = 42
CONTAINS
    SUBROUTINE SUB(X, Y, Z, ARR)
        REAL, INTENT(IN) :: X
        INTEGER, INTENT(IN) :: Y
        REAL, INTENT(OUT) :: Z
        REAL, DIMENSION(10), INTENT(IN) :: ARR
        REAL :: LOCAL_REAL
        INTEGER :: LOOP_VAR
        Z = X + REAL(Y)
        LOCAL_REAL = X * 2.0
    END SUBROUTINE SUB
END MODULE MOD_VAR
"""
    import tempfile
    import os
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(code)
        yield PYFT(fpath)


@pytest.fixture
def pft_subroutine_vars():
    """PYFT instance with subroutine scope for testing varList."""
    code = """
MODULE MOD_VAR
    IMPLICIT NONE
CONTAINS
    SUBROUTINE TEST_SUB(X, Y, ARR)
        REAL, INTENT(IN) :: X
        INTEGER, INTENT(IN) :: Y
        REAL, DIMENSION(10), INTENT(IN) :: ARR
        REAL :: LOCAL_REAL
        REAL :: RESULT
        RESULT = X + REAL(Y)
    END SUBROUTINE TEST_SUB
END MODULE MOD_VAR
"""
    import tempfile
    import os
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(code)
        pft = PYFT(fpath)
        # Get the subroutine scope
        scopes = pft.getScopes()
        for scope in scopes:
            if 'sub:TEST_SUB' in scope.path:
                yield scope
                break


class TestVarListFindVar:
    """Tests for findVar() method."""

    def test_find_scalar_argument(self, pft_subroutine_vars):
        """Test finding scalar argument."""
        vl = pft_subroutine_vars.varList
        var = vl.findVar('X')
        assert var is not None
        assert var['n'] == 'X'
        assert var['arg'] is True
        assert var['as'] == []

    def test_find_integer_argument(self, pft_subroutine_vars):
        """Test finding integer argument."""
        vl = pft_subroutine_vars.varList
        var = vl.findVar('Y')
        assert var is not None
        assert var['n'] == 'Y'
        assert 'INTEGER' in var['t']

    def test_find_array_argument(self, pft_subroutine_vars):
        """Test finding array argument."""
        vl = pft_subroutine_vars.varList
        var = vl.findVar('ARR')
        assert var is not None
        assert var['n'] == 'ARR'
        assert var['as'] != []

    def test_find_local_variable(self, pft_subroutine_vars):
        """Test finding local variable."""
        vl = pft_subroutine_vars.varList
        var = vl.findVar('LOCAL_REAL')
        assert var is not None
        assert var['n'] == 'LOCAL_REAL'
        assert var['arg'] is False

    def test_find_nonexistent(self, pft_subroutine_vars):
        """Test finding nonexistent variable returns None."""
        vl = pft_subroutine_vars.varList
        var = vl.findVar('NONEXISTENT')
        assert var is None

    def test_find_with_array_true(self, pft_subroutine_vars):
        """Test finding only arrays."""
        vl = pft_subroutine_vars.varList
        var = vl.findVar('ARR', array=True)
        assert var is not None
        assert var['as'] != []

    def test_find_with_array_false(self, pft_subroutine_vars):
        """Test finding only scalars."""
        vl = pft_subroutine_vars.varList
        var = vl.findVar('X', array=False)
        assert var is not None
        assert var['as'] == []

    def test_find_with_exact_scope(self, pft_subroutine_vars):
        """Test finding with exact scope."""
        vl = pft_subroutine_vars.varList
        var = vl.findVar('X', exactScope=True)
        assert var is not None


class TestVarListRestrict:
    """Tests for restrict() method."""

    def test_restrict_to_scope(self, pft_with_vars):
        """Test restricting to a specific scope."""
        vl = pft_with_vars.varList
        sub_vl = vl.restrict('module:MOD_VAR/sub:SUB', excludeContains=True)
        assert sub_vl is not None

    def test_restrict_excludes_contains(self, pft_with_vars):
        """Test restrict with excludeContains."""
        vl = pft_with_vars.varList
        vl_with = vl.restrict('module:MOD_VAR', excludeContains=False)
        vl_without = vl.restrict('module:MOD_VAR', excludeContains=True)
        assert len(vl_with) >= len(vl_without)


class TestVarListProperty:
    """Tests for VarList property access."""

    def test_var_list_exists(self, pft_with_vars):
        """Test varList property exists."""
        assert hasattr(pft_with_vars, 'varList')

    def test_var_list_is_list(self, pft_with_vars):
        """Test varList returns accessible via indexing."""
        vl = pft_with_vars.varList
        assert hasattr(vl, '__getitem__')
        assert hasattr(vl, '__len__')

    def test_var_list_has_variables(self, pft_with_vars):
        """Test varList contains variables."""
        vl = pft_with_vars.varList
        assert len(vl) > 0


class TestVarListShowVarList:
    """Tests for showVarList() method."""

    def test_show_var_list_runs(self, pft_with_vars, capsys):
        """Test showVarList() runs without error."""
        vl = pft_with_vars.varList
        vl.showVarList()
        captured = capsys.readouterr()
        assert len(captured.out) > 0 or captured.err == ''
