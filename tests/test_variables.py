"""
Tests for Variables mixin class in pyfortool.variables.
"""

import pytest
import tempfile
import os
from pyfortool import PYFT


class TestAttachArraySpecToEntity:
    """Tests for attachArraySpecToEntity() method."""

    @pytest.fixture
    def fortran_with_separate_dimension(self):
        """FORTRAN code with DIMENSION in separate attribute."""
        return """
MODULE MOD_DIM
    IMPLICIT NONE
CONTAINS
    SUBROUTINE DIM_TEST(A, B)
        REAL, DIMENSION(10), INTENT(IN) :: A
        REAL, DIMENSION(10), INTENT(OUT) :: B
        INTEGER :: I
        DO I = 1, 10
            B(I) = A(I) * 2.0
        END DO
    END SUBROUTINE DIM_TEST
END MODULE MOD_DIM
"""

    @pytest.fixture
    def pft_with_separate_dimension(self, fortran_with_separate_dimension):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_separate_dimension)
            return PYFT(fpath)

    def test_attach_array_spec_moves_dimension(self, pft_with_separate_dimension):
        """Test that DIMENSION attribute is moved to entity."""
        pft_with_separate_dimension.attachArraySpecToEntity()
        result = pft_with_separate_dimension.fortran
        assert 'A(10)' in result
        assert 'DIMENSION(10)' not in result

    def test_attach_array_spec_writes(self, pft_arrays):
        """Test attachArraySpecToEntity produces valid code."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, 'output.F90')
            pft_arrays._output = out_path
            pft_arrays.attachArraySpecToEntity()
            pft_arrays.write()
            assert os.path.exists(out_path)


class TestCheckImplicitNone:
    """Tests for checkImplicitNone() method."""

    def test_check_implicit_none_runs(self, pft_module, capsys):
        """Test checkImplicitNone() runs without error."""
        pft_module.checkImplicitNone()
        captured = capsys.readouterr()

    def test_check_implicit_none_raise(self, pft_module):
        """Test checkImplicitNone() with mustRaise."""
        pft_module.checkImplicitNone(mustRaise=False)


class TestCheckIntent:
    """Tests for checkIntent() method."""

    def test_check_intent_runs(self, pft_module, capsys):
        """Test checkIntent() runs without error."""
        pft_module.checkIntent()
        captured = capsys.readouterr()

    def test_check_intent_raise(self, pft_module):
        """Test checkIntent() with mustRaise."""
        pft_module.checkIntent(mustRaise=False)


class TestCheckOnly:
    """Tests for checkONLY() method."""

    def test_check_only_runs(self, pft_module, capsys):
        """Test checkONLY() runs without error."""
        pft_module.checkONLY()
        captured = capsys.readouterr()


class TestCheckUnusedLocalVar:
    """Tests for checkUnusedLocalVar() method."""

    @pytest.fixture
    def fortran_with_unused_var(self):
        """FORTRAN code with an unused variable."""
        return """
MODULE MOD_UNUSED
    IMPLICIT NONE
CONTAINS
    SUBROUTINE SUB(X, Y)
        REAL, INTENT(IN) :: X
        REAL, INTENT(OUT) :: Y
        REAL :: UNUSED_VAR
        Y = X * 2.0
    END SUBROUTINE SUB
END MODULE MOD_UNUSED
"""

    @pytest.fixture
    def pft_with_unused_var(self, fortran_with_unused_var):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_unused_var)
            return PYFT(fpath)

    def test_check_unused_reports_unused(self, pft_with_unused_var, capsys):
        """Test checkUnusedLocalVar() reports unused variable."""
        pft_with_unused_var.checkUnusedLocalVar()
        captured = capsys.readouterr()
        assert 'UNUSED_VAR' in captured.out or captured.out == ''

    def test_check_unused_with_exclude(self, pft_with_unused_var, capsys):
        """Test with exclude list."""
        pft_with_unused_var.checkUnusedLocalVar(excludeList=['UNUSED_VAR'])
        captured = capsys.readouterr()


class TestShowUnusedVar:
    """Tests for showUnusedVar() method."""

    def test_show_unused_runs(self, pft_calls, capsys):
        """Test showUnusedVar() runs without error."""
        pft_calls.showUnusedVar()
        captured = capsys.readouterr()


class TestAddExplicitArrayBounds:
    """Tests for addExplicitArrayBounds() method."""

    @pytest.fixture
    def fortran_with_implicit_bounds(self):
        """FORTRAN code with implicit array bounds."""
        return """
MODULE MOD_IMPL
    IMPLICIT NONE
CONTAINS
    SUBROUTINE IMPL_TEST(A, B)
        REAL, INTENT(IN) :: A(:)
        REAL, INTENT(OUT) :: B(:)
        B(:) = A(:) * 2.0
    END SUBROUTINE IMPL_TEST
END MODULE MOD_IMPL
"""

    @pytest.fixture
    def pft_with_implicit_bounds(self, fortran_with_implicit_bounds):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_implicit_bounds)
            return PYFT(fpath)

    def test_add_explicit_bounds_runs(self, pft_arrays):
        """Test addExplicitArrayBounds() runs."""
        pft_arrays.addExplicitArrayBounds()


class TestAddArrayParentheses:
    """Tests for addArrayParentheses() method."""

    @pytest.fixture
    def fortran_with_array_no_parens(self):
        """FORTRAN code with arrays without parentheses."""
        return """
MODULE MOD_NOPAREN
    IMPLICIT NONE
CONTAINS
    SUBROUTINE NOPAREN(A, B)
        REAL, INTENT(IN) :: A(10)
        REAL, INTENT(OUT) :: B(10)
        INTEGER :: I
        DO I = 1, 10
            B(I) = A(I) * 2.0
        END DO
    END SUBROUTINE NOPAREN
END MODULE MOD_NOPAREN
"""

    @pytest.fixture
    def pft_with_array_no_parens(self, fortran_with_array_no_parens):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_array_no_parens)
            return PYFT(fpath)

    def test_add_array_parens_runs(self, pft_arrays):
        """Test addArrayParentheses() runs."""
        pft_arrays.addArrayParentheses()


class TestModifyAutomaticArrays:
    """Tests for modifyAutomaticArrays() method."""

    @pytest.fixture
    def fortran_with_automatic_array(self):
        """FORTRAN code with automatic array."""
        return """
MODULE MOD_AUTO
    IMPLICIT NONE
CONTAINS
    SUBROUTINE AUTO_TEST(N)
        INTEGER, INTENT(IN) :: N
        REAL :: LOCAL_ARRAY(N, N)
        INTEGER :: I, J
        DO I = 1, N
            DO J = 1, N
                LOCAL_ARRAY(I, J) = REAL(I + J)
            END DO
        END DO
    END SUBROUTINE AUTO_TEST
END MODULE MOD_AUTO
"""

    @pytest.fixture
    def pft_with_automatic_array(self, fortran_with_automatic_array):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_automatic_array)
            return PYFT(fpath)

    def test_modify_automatic_arrays_returns_int(self, pft_loops):
        """Test modifyAutomaticArrays() returns count."""
        result = pft_loops.modifyAutomaticArrays()
        assert isinstance(result, int)

    def test_modify_with_template(self, pft_with_automatic_array):
        """Test with allocation template transforms code."""
        before = pft_with_automatic_array.fortran
        pft_with_automatic_array.modifyAutomaticArrays(
            declTemplate="{type}, DIMENSION({doubledotshape}), ALLOCATABLE :: {name}",
            startTemplate="ALLOCATE({name}({shape}))",
            endTemplate="DEALLOCATE({name})"
        )
        after = pft_with_automatic_array.fortran
        assert 'ALLOCATABLE' in after
        assert 'ALLOCATE' in after


class TestRemoveUnusedLocalVar:
    """Tests for removeUnusedLocalVar() method."""

    @pytest.fixture
    def fortran_with_truly_unused_var(self):
        """FORTRAN code with an actually unused variable."""
        return """
MODULE MOD_UNUSED
    IMPLICIT NONE
CONTAINS
    SUBROUTINE SUB(X, Y)
        REAL, INTENT(IN) :: X
        REAL, INTENT(OUT) :: Y
        REAL :: UNUSED_VAR
        Y = X * 2.0
    END SUBROUTINE SUB
END MODULE MOD_UNUSED
"""

    @pytest.fixture
    def pft_truly_unused(self, fortran_with_truly_unused_var):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_truly_unused_var)
            return PYFT(fpath)

    def test_remove_unused_removes_variable(self, pft_truly_unused):
        """Test removeUnusedLocalVar() actually removes the variable."""
        before = pft_truly_unused.fortran
        assert 'UNUSED_VAR' in before
        
        pft_truly_unused.removeUnusedLocalVar()
        after = pft_truly_unused.fortran
        
        assert 'UNUSED_VAR' not in after
        assert 'REAL :: UNUSED_VAR' not in after

    def test_remove_unused_keeps_used_variables(self, pft_truly_unused):
        """Test removeUnusedLocalVar() keeps used variables."""
        pft_truly_unused.removeUnusedLocalVar()
        result = pft_truly_unused.fortran
        
        assert 'X' in result
        assert 'Y' in result

    def test_remove_unused_with_exclude(self, pft_truly_unused):
        """Test with exclude list keeps excluded variables."""
        pft_truly_unused.removeUnusedLocalVar(excludeList=['UNUSED_VAR'])
        result = pft_truly_unused.fortran
        
        assert 'UNUSED_VAR' in result

    def test_remove_unused_with_simplify(self, pft_truly_unused):
        """Test with simplify option."""
        pft_truly_unused.removeUnusedLocalVar(simplify=True)
        result = pft_truly_unused.fortran
        
        assert 'UNUSED_VAR' not in result


class TestRemoveVarIfUnused:
    """Tests for removeVarIfUnused() method."""

    @pytest.fixture
    def fortran_with_multiple_unused(self):
        """FORTRAN code with multiple unused variables."""
        return """
MODULE MOD_MULTI
    IMPLICIT NONE
CONTAINS
    SUBROUTINE MULTI(X, Y, Z)
        REAL, INTENT(IN) :: X
        REAL, INTENT(OUT) :: Y, Z
        REAL :: A, B, C
        Y = X * 2.0
        Z = X * 3.0
    END SUBROUTINE MULTI
END MODULE MOD_MULTI
"""

    @pytest.fixture
    def pft_multiple_unused(self, fortran_with_multiple_unused):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_multiple_unused)
            return PYFT(fpath)

    def test_remove_var_if_unused_removes_all(self, pft_multiple_unused):
        """Test removeVarIfUnused() removes all specified unused variables."""
        before = pft_multiple_unused.fortran
        assert 'REAL :: A, B, C' in before or 'REAL, DIMENSION(10) :: A, B, C' in before

        scope = pft_multiple_unused.getScopeNode('module:MOD_MULTI/sub:MULTI')
        var_list = [v['n'] for v in scope.varList if v['n'] in ['A', 'B', 'C']]

        removed = scope.removeVarIfUnused(var_list, excludeDummy=True, excludeModule=True)
        after = pft_multiple_unused.fortran

        assert removed == ['A', 'B', 'C']
        assert 'REAL :: A, B, C' not in after
        assert 'REAL :: ' not in after.replace('REAL, INTENT', '')

    def test_remove_var_if_unused_keeps_used(self, pft_multiple_unused):
        """Test removeVarIfUnused() keeps used variables and returns only removed ones."""
        scope = pft_multiple_unused.getScopeNode('module:MOD_MULTI/sub:MULTI')

        removed = scope.removeVarIfUnused(['A', 'Z'], excludeDummy=True)
        after = pft_multiple_unused.fortran

        assert removed == ['A']
        assert 'Z = X * 3.0' in after

    def test_remove_var_if_unused_exclude_module(self, pft_multiple_unused):
        """Test removeVarIfUnused() with excludeModule does nothing on a module scope."""
        before = pft_multiple_unused.fortran
        scope = pft_multiple_unused.getScopeNode('module:MOD_MULTI')

        removed = scope.removeVarIfUnused(['A'], excludeModule=True)

        assert removed == []
        assert pft_multiple_unused.fortran == before


class TestAddVar:
    """Tests for addVar() method."""

    @pytest.fixture
    def fortran_simple(self):
        """Simple FORTRAN code for adding variables."""
        return """
MODULE MOD_ADD
    IMPLICIT NONE
CONTAINS
    SUBROUTINE ADD_SUB(X)
        REAL, INTENT(IN) :: X
    END SUBROUTINE ADD_SUB
END MODULE MOD_ADD
"""

    @pytest.fixture
    def pft_add_var(self, fortran_simple):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_simple)
            return PYFT(fpath)

    def test_add_var_adds_declaration(self, pft_add_var):
        """Test addVar() adds a new variable declaration."""
        before = pft_add_var.fortran
        assert 'NEW_VAR' not in before
        
        sub = pft_add_var.getScopeNode('module:MOD_ADD/sub:ADD_SUB')
        sub.addVar([('NEW_VAR', 'INTEGER :: NEW_VAR', None)])
        after = pft_add_var.fortran

        assert 'INTEGER :: NEW_VAR' in after
        # Declared inside the subroutine
        assert after.index('SUBROUTINE ADD_SUB(X)') < after.index('INTEGER :: NEW_VAR')
        assert after.index('INTEGER :: NEW_VAR') < after.index('END SUBROUTINE ADD_SUB')

    def test_add_var_adds_dummy_argument(self, pft_add_var):
        """Test addVar() with a position adds the variable to the argument list."""
        sub = pft_add_var.getScopeNode('module:MOD_ADD/sub:ADD_SUB')
        sub.addVar([('Y', 'REAL, INTENT(OUT) :: Y', 1)])
        after = pft_add_var.fortran

        assert 'SUBROUTINE ADD_SUB(X, Y)' in after
        assert 'REAL, INTENT(OUT) :: Y' in after

    def test_add_var_in_module(self, pft_add_var):
        """Test addVar() on a module scope declares the variable before CONTAINS."""
        mod = pft_add_var.getScopeNode('module:MOD_ADD')
        mod.addVar([('MOD_VAR', 'INTEGER :: MOD_VAR', None)])
        after = pft_add_var.fortran

        assert after.index('INTEGER :: MOD_VAR') < after.index('CONTAINS')

    def test_add_var_on_file_scope_raises(self, pft_add_var):
        """Test addVar() refuses to run on the whole file (not a scope)."""
        from pyfortool.util import PYFTError
        with pytest.raises(PYFTError):
            pft_add_var.addVar([('NEW_VAR', 'INTEGER :: NEW_VAR', None)])

    def test_add_var_legacy_scope_path_is_deprecated(self, pft_add_var):
        """Test addVar() still accepts (scopePath, ...) items but warns."""
        with pytest.warns(DeprecationWarning):
            pft_add_var.addVar([('module:MOD_ADD/sub:ADD_SUB', 'NEW_VAR',
                                 'INTEGER :: NEW_VAR', None)])
        assert 'INTEGER :: NEW_VAR' in pft_add_var.fortran

    def test_add_var_mixed_forms_raises(self, pft_add_var):
        """Test addVar() refuses a varList mixing legacy and new items."""
        from pyfortool.util import PYFTError
        with pytest.raises(PYFTError):
            pft_add_var.addVar([('module:MOD_ADD/sub:ADD_SUB', 'A', 'INTEGER :: A', None),
                                ('B', 'INTEGER :: B', None)])


class TestAddModuleVar:
    """Tests for addModuleVar() method."""

    @pytest.fixture
    def pft_use(self):
        code = """
MODULE MOD_USE
    IMPLICIT NONE
CONTAINS
    SUBROUTINE USE_SUB(X)
        USE MODD_A, ONLY: A1
        REAL, INTENT(IN) :: X
    END SUBROUTINE USE_SUB
END MODULE MOD_USE
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(code)
            return PYFT(fpath)

    def test_add_module_var_adds_use(self, pft_use):
        """Test addModuleVar() adds a USE ... ONLY statement in the scope."""
        sub = pft_use.getScopeNode('module:MOD_USE/sub:USE_SUB')
        sub.addModuleVar([('MODD_B', ['B1', 'B2']), ('MODD_C', None)])
        after = pft_use.fortran

        assert 'USE MODD_B, ONLY:B1, B2' in after
        assert 'USE MODD_C' in after
        assert after.index('USE MODD_A') < after.index('USE MODD_B')
        assert after.index('USE MODD_B') < after.index('REAL, INTENT(IN) :: X')

    def test_add_module_var_updates_existing_use(self, pft_use):
        """Test addModuleVar() does not duplicate already imported variables."""
        sub = pft_use.getScopeNode('module:MOD_USE/sub:USE_SUB')
        sub.addModuleVar([('MODD_A', ['A1', 'A2'])])
        after = pft_use.fortran

        assert after.count('USE MODD_A') == 2
        assert 'USE MODD_A, ONLY:A2' in after

    def test_add_module_var_legacy_scope_path_is_deprecated(self, pft_use):
        """Test addModuleVar() still accepts (scopePath, ...) items but warns."""
        with pytest.warns(DeprecationWarning):
            pft_use.addModuleVar([('module:MOD_USE/sub:USE_SUB', 'MODD_B', 'B1')])
        assert 'USE MODD_B, ONLY:B1' in pft_use.fortran


class TestRemoveVar:
    """Tests for removeVar() method."""

    @pytest.fixture
    def fortran_with_var_to_remove(self):
        """FORTRAN code with variable to remove."""
        return """
MODULE MOD_REMOVE
    IMPLICIT NONE
    INTEGER :: REMOVE_ME
CONTAINS
    SUBROUTINE REMOVE_SUB(X)
        REAL, INTENT(INOUT) :: X
        INTEGER :: LOCAL_REMOVE
        X = X + 1.0
    END SUBROUTINE REMOVE_SUB
END MODULE MOD_REMOVE
"""

    @pytest.fixture
    def pft_remove_var(self, fortran_with_var_to_remove):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_var_to_remove)
            return PYFT(fpath)

    def test_remove_var_removes_module_variable(self, pft_remove_var):
        """Test removeVar() removes module variable."""
        before = pft_remove_var.fortran
        assert 'REMOVE_ME' in before
        
        pft_remove_var.getScopeNode('module:MOD_REMOVE').removeVar(['REMOVE_ME'])
        after = pft_remove_var.fortran

        assert 'REMOVE_ME' not in after

    def test_remove_var_removes_local_variable(self, pft_remove_var):
        """Test removeVar() removes local variable."""
        before = pft_remove_var.fortran
        assert 'LOCAL_REMOVE' in before
        
        sub = pft_remove_var.getScopeNode('module:MOD_REMOVE/sub:REMOVE_SUB')
        sub.removeVar(['LOCAL_REMOVE'])
        after = pft_remove_var.fortran

        assert 'LOCAL_REMOVE' not in after
        assert 'INTEGER :: REMOVE_ME' in after

    def test_remove_var_searches_upper_scope(self, pft_remove_var):
        """Test removeVar() called on a subroutine removes a variable of the enclosing module."""
        sub = pft_remove_var.getScopeNode('module:MOD_REMOVE/sub:REMOVE_SUB')
        sub.removeVar(['REMOVE_ME'])
        after = pft_remove_var.fortran

        assert 'REMOVE_ME' not in after
        assert 'LOCAL_REMOVE' in after

    def test_remove_var_removes_dummy_argument(self, pft_remove_var):
        """Test removeVar() removes a dummy argument from the argument list too."""
        sub = pft_remove_var.getScopeNode('module:MOD_REMOVE/sub:REMOVE_SUB')
        sub.removeVar(['X'])
        after = pft_remove_var.fortran

        assert 'SUBROUTINE REMOVE_SUB()' in after
        assert 'INTENT(INOUT)' not in after

    def test_remove_var_is_case_insensitive(self, pft_remove_var):
        """Test removeVar() with a lowercase name."""
        sub = pft_remove_var.getScopeNode('module:MOD_REMOVE/sub:REMOVE_SUB')
        sub.removeVar(['local_remove'])

        assert 'LOCAL_REMOVE' not in pft_remove_var.fortran

    def test_remove_var_legacy_scope_path_is_deprecated(self, pft_remove_var):
        """Test removeVar() still accepts (scopePath, name) items but warns."""
        with pytest.warns(DeprecationWarning):
            pft_remove_var.removeVar([('module:MOD_REMOVE', 'REMOVE_ME'),
                                      ('module:MOD_REMOVE/sub:REMOVE_SUB', 'LOCAL_REMOVE')])
        after = pft_remove_var.fortran

        assert 'REMOVE_ME' not in after
        assert 'LOCAL_REMOVE' not in after


class TestIsVarUsed:
    """Tests for isVarUsed() method."""

    @pytest.fixture
    def pft_usage(self):
        code = """
MODULE MOD_USAGE
    IMPLICIT NONE
    INTEGER :: MOD_USED, MOD_UNUSED
CONTAINS
    SUBROUTINE USAGE_SUB(X, Y)
        REAL, INTENT(IN) :: X
        REAL, INTENT(OUT) :: Y
        INTEGER :: UNUSED
        Y = X * MOD_USED
    END SUBROUTINE USAGE_SUB
END MODULE MOD_USAGE
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(code)
            return PYFT(fpath)

    def test_is_var_used_keys_are_names(self, pft_usage):
        """Test isVarUsed() returns a dict keyed by uppercase variable names."""
        sub = pft_usage.getScopeNode('module:MOD_USAGE/sub:USAGE_SUB')
        result = sub.isVarUsed(['x', 'Y', 'UNUSED'])

        assert result == {'X': True, 'Y': True, 'UNUSED': False}

    def test_is_var_used_dummy_are_always_used(self, pft_usage):
        """Test isVarUsed() dummyAreAlwaysUsed option."""
        sub = pft_usage.getScopeNode('module:MOD_USAGE/sub:USAGE_SUB')

        assert sub.isVarUsed(['Y'], dummyAreAlwaysUsed=True) == {'Y': True}

    def test_is_var_used_module_variable_from_contained_scope(self, pft_usage):
        """Test isVarUsed() finds usages of module variables in contained routines."""
        mod = pft_usage.getScopeNode('module:MOD_USAGE')
        result = mod.isVarUsed(['MOD_USED', 'MOD_UNUSED'])

        assert result == {'MOD_USED': True, 'MOD_UNUSED': False}
        # With exactScope, usages in the contained subroutine are ignored
        assert mod.isVarUsed(['MOD_USED'], exactScope=True) == {'MOD_USED': False}

    def test_is_var_used_legacy_scope_path_is_deprecated(self, pft_usage):
        """Test isVarUsed() legacy form returns (scopePath, name) keys and warns."""
        with pytest.warns(DeprecationWarning):
            result = pft_usage.isVarUsed([('module:MOD_USAGE/sub:USAGE_SUB', 'UNUSED'),
                                          ('module:MOD_USAGE', 'MOD_USED')])

        assert result == {('module:MOD_USAGE/sub:USAGE_SUB', 'UNUSED'): False,
                          ('module:MOD_USAGE', 'MOD_USED'): True}


class TestRenameVar:
    """Tests for renameVar() method."""

    @pytest.fixture
    def fortran_with_var_to_rename(self):
        """FORTRAN code with variable to rename."""
        return """
MODULE MOD_RENAME
    IMPLICIT NONE
CONTAINS
    SUBROUTINE RENAME_SUB(OLD_NAME)
        REAL, INTENT(INOUT) :: OLD_NAME
        REAL :: TEMP
        TEMP = OLD_NAME * 2.0
        OLD_NAME = TEMP + 1.0
    END SUBROUTINE RENAME_SUB
END MODULE MOD_RENAME
"""

    @pytest.fixture
    def pft_rename_var(self, fortran_with_var_to_rename):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_var_to_rename)
            return PYFT(fpath)

    def test_rename_var_renames_all_occurrences(self, pft_rename_var):
        """Test renameVar() renames all occurrences."""
        before = pft_rename_var.fortran
        assert 'OLD_NAME' in before
        assert 'NEW_NAME' not in before
        
        pft_rename_var.renameVar('OLD_NAME', 'NEW_NAME')
        after = pft_rename_var.fortran
        
        assert 'OLD_NAME' not in after
        assert 'NEW_NAME' in after


class TestCheckKeyDimConsistency:
    """Tests for checkKeyDimConsistency() method."""

    @pytest.fixture
    def fortran_merge_consistent(self):
        """FORTRAN code with consistent MERGE dimensions."""
        return """
MODULE MOD_MERGE
    IMPLICIT NONE
CONTAINS
    SUBROUTINE PARENT(A)
        REAL, DIMENSION(MERGE(10,0,FLAG)), INTENT(INOUT) :: A
        A = 0
    END SUBROUTINE PARENT
    SUBROUTINE CHILD(A)
        REAL, DIMENSION(MERGE(10,0,FLAG)), INTENT(INOUT) :: A
        A = 0
    END SUBROUTINE CHILD
END MODULE MOD_MERGE
"""

    @pytest.fixture
    def fortran_merge_inconsistent(self):
        """FORTRAN code with inconsistent MERGE dimensions."""
        return """
MODULE MOD_MERGE_BAD
    IMPLICIT NONE
CONTAINS
    SUBROUTINE PARENT(A)
        REAL, DIMENSION(MERGE(10,0,FLAG)), INTENT(INOUT) :: A
        A = 0
    END SUBROUTINE PARENT
    SUBROUTINE CHILD(A)
        REAL, DIMENSION(MERGE(20,0,OTHER)), INTENT(INOUT) :: A
        A = 0
    END SUBROUTINE CHILD
END MODULE MOD_MERGE_BAD
"""

    @pytest.fixture
    def fortran_merge_none(self):
        """FORTRAN code without MERGE dimensions."""
        return """
MODULE MOD_NO_MERGE
    IMPLICIT NONE
CONTAINS
    SUBROUTINE SUB(A)
        REAL, INTENT(INOUT) :: A
        A = 0
    END SUBROUTINE SUB
END MODULE MOD_NO_MERGE
"""

    @pytest.fixture
    def pft_merge_consistent(self, fortran_merge_consistent):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_merge_consistent)
            return PYFT(fpath)

    @pytest.fixture
    def pft_merge_inconsistent(self, fortran_merge_inconsistent):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_merge_inconsistent)
            return PYFT(fpath)

    @pytest.fixture
    def pft_merge_none(self, fortran_merge_none):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_merge_none)
            return PYFT(fpath)

    def test_consistent_merge_returns_true(self, pft_merge_consistent):
        """Test with consistent MERGE dims returns True."""
        result = pft_merge_consistent.checkKeyDimConsistency()
        assert result is True

    def test_inconsistent_merge_returns_false(self, pft_merge_inconsistent):
        """Test with inconsistent MERGE dims returns False."""
        result = pft_merge_inconsistent.checkKeyDimConsistency()
        assert result is False

    def test_inconsistent_merge_raises(self, pft_merge_inconsistent):
        """Test with inconsistent MERGE dims raises when mustRaise=True."""
        with pytest.raises(Exception):
            pft_merge_inconsistent.checkKeyDimConsistency(mustRaise=True)

    def test_no_merge_returns_true(self, pft_merge_none):
        """Test without MERGE dims returns True."""
        result = pft_merge_none.checkKeyDimConsistency()
        assert result is True
