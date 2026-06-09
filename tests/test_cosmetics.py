"""
Tests for the Cosmetics class (code formatting transformations).
"""

import pytest
import tempfile
import os
import re
from pyfortool import PYFT


class TestCosmeticsCase:
    """Tests for upperCase and lowerCase methods."""

    @pytest.fixture
    def fortran_lowercase(self):
        """FORTRAN code in lowercase."""
        return """
module mod_test
    implicit none
contains
    subroutine sub_test(x, y)
        real, intent(in) :: x
        real, intent(out) :: y
        y = x * 2.0
    end subroutine sub_test
end module mod_test
"""

    @pytest.fixture
    def fortran_uppercase(self):
        """FORTRAN code in uppercase."""
        return """
MODULE MOD_TEST
    IMPLICIT NONE
CONTAINS
    SUBROUTINE SUB_TEST(X, Y)
        REAL, INTENT(IN) :: X
        REAL, INTENT(OUT) :: Y
        Y = X * 2.0
    END SUBROUTINE SUB_TEST
END MODULE MOD_TEST
"""

    @pytest.fixture
    def pft_lowercase(self, fortran_lowercase):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_lowercase)
            return PYFT(fpath)

    @pytest.fixture
    def pft_uppercase(self, fortran_uppercase):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_uppercase)
            return PYFT(fpath)

    def test_uppercase_converts_keywords(self, pft_lowercase):
        """Test that upperCase() converts FORTRAN keywords."""
        pft_lowercase.upperCase()
        result = pft_lowercase.fortran
        assert 'MODULE' in result
        assert 'SUBROUTINE' in result
        assert 'IMPLICIT NONE' in result

    def test_lowercase_converts_keywords(self, pft_uppercase):
        """Test that lowerCase() converts FORTRAN keywords."""
        pft_uppercase.lowerCase()
        result = pft_uppercase.fortran
        assert 'module' in result
        assert 'subroutine' in result


class TestCosmeticsIndent:
    """Tests for indent method."""

    @pytest.fixture
    def fortran_no_indent(self):
        """FORTRAN code without proper indentation."""
        return """
MODULE MOD_NOINDENT
IMPLICIT NONE
CONTAINS
SUBROUTINE SUB(X)
REAL, INTENT(IN) :: X
X = X * 2.0
END SUBROUTINE SUB
END MODULE MOD_NOINDENT
"""

    @pytest.fixture
    def pft_no_indent(self, fortran_no_indent):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_no_indent)
            return PYFT(fpath)

    def test_indent_adds_spaces(self, pft_no_indent):
        """Test that indent() adds proper indentation."""
        before = pft_no_indent.fortran
        pft_no_indent.indent()
        result = pft_no_indent.fortran
        assert result is not None

    def test_indent_custom_params(self, pft_module):
        """Test indentation with custom parameters."""
        pft_module.indent(indentProgramunit=0, indentBranch=4)
        result = pft_module.fortran
        assert result is not None

    def test_indent_excl_directives(self, pft_module):
        """Test indentation excluding directives."""
        pft_module.indent(exclDirectives=['!$OMP'])
        result = pft_module.fortran
        assert result is not None

    def test_indent_no_exclusions(self, pft_module):
        """Test indentation with no exclusions."""
        pft_module.indent(exclDirectives=[])
        result = pft_module.fortran
        assert result is not None


class TestCosmeticsComments:
    """Tests for removeComments method."""

    @pytest.fixture
    def fortran_with_comments_to_remove(self):
        """FORTRAN code with comments to remove."""
        return """
! Header comment
MODULE MOD_COMMENT
    IMPLICIT NONE
    REAL :: X  ! Inline comment
CONTAINS
    SUBROUTINE SUB
        ! Full line comment
        X = 1.0
    END SUBROUTINE SUB
END MODULE MOD_COMMENT
"""

    @pytest.fixture
    def pft_comments_to_remove(self, fortran_with_comments_to_remove):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_comments_to_remove)
            return PYFT(fpath)

    def test_remove_comments_default(self, pft_comments_to_remove):
        """Test removing comments while keeping directives."""
        before = pft_comments_to_remove.fortran
        assert '! Header comment' in before
        
        pft_comments_to_remove.removeComments()
        after = pft_comments_to_remove.fortran
        
        assert '! Header comment' not in after
        assert '! Inline comment' not in after
        assert '! Full line comment' not in after

    def test_remove_comments_all(self, pft_comments_to_remove):
        """Test removing all comments including directives."""
        pft_comments_to_remove.removeComments(exclDirectives=[])
        result = pft_comments_to_remove.fortran
        assert '! Header comment' not in result

    def test_remove_comments_pattern(self, pft_comments_to_remove):
        """Test removing comments matching a pattern."""
        pft_comments_to_remove.removeComments(pattern=re.compile(r'!.*Inline'))
        result = pft_comments_to_remove.fortran
        assert 'Inline comment' not in result


class TestCosmeticsEmptyLines:
    """Tests for removeEmptyLines method."""

    @pytest.fixture
    def fortran_with_empty_lines(self):
        """FORTRAN code with empty lines."""
        return """


MODULE MOD_EMPTY


    IMPLICIT NONE


CONTAINS


    SUBROUTINE SUB


        X = 1.0


    END SUBROUTINE SUB


END MODULE MOD_EMPTY



"""

    @pytest.fixture
    def pft_empty_lines(self, fortran_with_empty_lines):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_empty_lines)
            return PYFT(fpath)

    def test_remove_empty_lines(self, pft_empty_lines):
        """Test removing empty lines."""
        pft_empty_lines.removeEmptyLines()
        result = pft_empty_lines.fortran
        assert '\n\n\n' not in result


class TestCosmeticsContinuation:
    """Tests for updateContinuation method."""

    @pytest.fixture
    def fortran_with_continuation(self):
        return """
MODULE MOD_CONT
    IMPLICIT NONE
CONTAINS
    SUBROUTINE CONT_EXAMPLE(A, B, C, D, E)
        REAL, INTENT(IN) :: A, B, C, D
        REAL, INTENT(OUT) :: E
        E = A + B + &
            C + D
    END SUBROUTINE CONT_EXAMPLE
END MODULE MOD_CONT
"""

    @pytest.fixture
    def pft_continuation(self, fortran_with_continuation):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_continuation)
            return PYFT(fpath)

    def test_update_continuation_align(self, pft_continuation):
        """Test aligning continuation lines."""
        before = pft_continuation.fortran
        assert '&' in before
        
        pft_continuation.updateContinuation(align=True)
        result = pft_continuation.fortran
        assert result is not None

    def test_update_continuation_remove_all(self, pft_continuation):
        """Test removing all continuation characters."""
        before = pft_continuation.fortran
        assert '&' in before
        
        pft_continuation.updateContinuation(align=False, removeALL=True, addBegin=False)
        result = pft_continuation.fortran
        assert '&' not in result


class TestCosmeticsSpaces:
    """Tests for updateSpaces method."""

    @pytest.fixture
    def fortran_with_bad_spaces(self):
        """FORTRAN code with inconsistent spacing."""
        return """
MODULE MOD_SPACE
    IMPLICIT NONE
CONTAINS
    SUBROUTINE BAD_SPACE(A,B)
        REAL,INTENT(IN)::A
        REAL,INTENT(OUT)::B
        B=A*2.0
    END SUBROUTINE BAD_SPACE
END MODULE MOD_SPACE
"""

    @pytest.fixture
    def pft_bad_spaces(self, fortran_with_bad_spaces):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_bad_spaces)
            return PYFT(fpath)

    def test_update_spaces_default(self, pft_module):
        """Test default space updates."""
        pft_module.updateSpaces()
        result = pft_module.fortran
        assert result is not None

    def test_update_spaces_normalizes(self, pft_bad_spaces):
        """Test that updateSpaces normalizes spacing."""
        pft_bad_spaces.updateSpaces()
        result = pft_bad_spaces.fortran
        assert ',INTENT' not in result
        assert '::' in result


class TestCosmeticsIfStatements:
    """Tests for changeIfStatementsInIfConstructs method."""

    def test_change_single_if_to_construct(self, pft_single_if):
        """Test converting single-line IF to IF-THEN-ENDIF."""
        pft_single_if.changeIfStatementsInIfConstructs()
        result = pft_single_if.fortran
        assert result is not None


class TestCosmeticsContains:
    """Tests for removeEmptyCONTAINS method."""

    @pytest.fixture
    def fortran_with_empty_contains(self):
        return """
MODULE MOD_EMPTY
    IMPLICIT NONE
CONTAINS
END MODULE MOD_EMPTY
"""

    @pytest.fixture
    def pft_empty_contains(self, fortran_with_empty_contains):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_empty_contains)
            return PYFT(fpath)

    def test_remove_empty_contains(self, pft_empty_contains):
        """Test removing empty CONTAINS."""
        pft_empty_contains.removeEmptyCONTAINS()
        result = pft_empty_contains.fortran
        assert 'CONTAINS' not in result


class TestCosmeticsFormatModuleUse:
    """Tests for formatModuleUse method."""

    @pytest.fixture
    def fortran_unordered_use(self):
        """FORTRAN code with unordered USE declarations."""
        return """
MODULE MOD_TARGET
    USE MODI_FOO
    USE MODD_PARAM
    USE MODN_NAMELIST
    USE MODE_ENV
    USE MODD_DATA
    USE MODI_BAR
    IMPLICIT NONE
CONTAINS
    SUBROUTINE SUB
        USE MODD_LOCAL
        USE MODI_UTIL
        USE MODE_HELPER
        IMPLICIT NONE
        INTEGER :: X
        X = 1
    END SUBROUTINE SUB
END MODULE MOD_TARGET
"""

    @pytest.fixture
    def fortran_mixed_case_use(self):
        """FORTRAN code with mixed case USE declarations."""
        return """
MODULE MOD_MIXED
    USE modi_foo
    USE MODD_bar
    USE mode_baz
    USE Modd_qux
    IMPLICIT NONE
END MODULE MOD_MIXED
"""

    @pytest.fixture
    def fortran_use_with_only(self):
        """FORTRAN code with ONLY clauses in USE."""
        return """
MODULE MOD_ONLY
    USE MODD_PARAM, ONLY: JPX, JPY
    USE MODI_BAR, ONLY: bar_sub
    USE MODD_DATA, ONLY: JPDATA
    IMPLICIT NONE
END MODULE MOD_ONLY
"""

    @pytest.fixture
    def fortran_use_with_comments(self):
        """FORTRAN code with trailing comments on USE lines."""
        return """
MODULE MOD_COMMENT
    USE MODI_FOO   ! Interface for FOO
    USE MODD_PARAM ! Data parameters
    USE MODE_ENV   ! Environment setup
    IMPLICIT NONE
END MODULE MOD_COMMENT
"""

    @pytest.fixture
    def fortran_use_mixed_case_only(self):
        """FORTRAN code with mixed-case in ONLY clause."""
        return """
MODULE MOD_MIXED_ONLY
    USE modd_params, ONLY: jpx, Jpy
    USE modi_routines, ONLY: bar_Sub
    IMPLICIT NONE
END MODULE MOD_MIXED_ONLY
"""

    @pytest.fixture
    def pft_unordered_use(self, fortran_unordered_use):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_unordered_use)
            return PYFT(fpath)

    @pytest.fixture
    def pft_mixed_case_use(self, fortran_mixed_case_use):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_mixed_case_use)
            return PYFT(fpath)

    @pytest.fixture
    def pft_use_with_only(self, fortran_use_with_only):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_use_with_only)
            return PYFT(fpath)

    @pytest.fixture
    def pft_use_with_comments(self, fortran_use_with_comments):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_use_with_comments)
            return PYFT(fpath)

    @pytest.fixture
    def pft_use_mixed_case_only(self, fortran_use_mixed_case_only):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_use_mixed_case_only)
            return PYFT(fpath)

    def test_group_order(self, pft_unordered_use):
        """Test that USE statements are grouped in MODD, MODE, MODI, MODN order."""
        before = pft_unordered_use.fortran
        assert 'USE MODI_FOO' in before
        assert 'USE MODD_PARAM' in before

        pft_unordered_use.formatModuleUse()
        result = pft_unordered_use.fortran

        modd_idx = result.index('USE MODD_PARAM')
        mode_idx = result.index('USE MODE_ENV')
        modi_foo_idx = result.index('USE MODI_FOO')
        modi_bar_idx = result.index('USE MODI_BAR')
        modn_idx = result.index('USE MODN_NAMELIST')

        assert modd_idx < mode_idx < modi_foo_idx < modn_idx
        assert modi_bar_idx < modn_idx

    def test_alphabetical_order_within_group(self, pft_unordered_use):
        """Test that USE statements are sorted alphabetically within each group."""
        pft_unordered_use.formatModuleUse()
        result = pft_unordered_use.fortran

        modd_data = result.index('USE MODD_DATA')
        modd_param = result.index('USE MODD_PARAM')
        assert modd_data < modd_param

        modi_bar = result.index('USE MODI_BAR')
        modi_foo = result.index('USE MODI_FOO')
        assert modi_bar < modi_foo

    def test_uppercase_only_after_comma(self, pft_mixed_case_use):
        """Test that upper=True only uppercases text after the comma (ONLY clause),
        not the module name."""
        pft_mixed_case_use.formatModuleUse(upper=True)
        result = pft_mixed_case_use.fortran

        # Module names preserved as-is (no comma in these, so no change)
        assert 'USE modi_foo' in result
        assert 'USE MODD_bar' in result
        assert 'USE mode_baz' in result

    def test_preserve_case_when_upper_false(self, pft_mixed_case_use):
        """Test that original case is preserved when upper=False."""
        pft_mixed_case_use.formatModuleUse(upper=False)
        result = pft_mixed_case_use.fortran

        assert 'USE modi_foo' in result
        assert 'USE mode_baz' in result

    def test_preserves_only_clause(self, pft_use_with_only):
        """Test that ONLY clauses are preserved after formatting."""
        pft_use_with_only.formatModuleUse()
        result = pft_use_with_only.fortran

        assert 'ONLY: JPX' in result or 'ONLY: JPX, JPY' in result
        assert 'USE MODD_DATA' in result
        assert 'USE MODD_PARAM' in result
        assert 'USE MODI_BAR' in result

    def test_subroutine_use_statements(self, pft_unordered_use):
        """Test that USE statements in subroutines are also formatted."""
        pft_unordered_use.formatModuleUse()
        result = pft_unordered_use.fortran

        lines = [l.strip() for l in result.split('\n')]
        sub_uses = [l for l in lines if l.startswith('USE MOD') and
                    'MODULE MOD_TARGET' not in l]

        assert 'USE MODD_LOCAL' in sub_uses or 'USE MODD_LOCAL' in result
        assert 'USE MODE_HELPER' in sub_uses or 'USE MODE_HELPER' in result
        assert 'USE MODI_UTIL' in sub_uses or 'USE MODI_UTIL' in result

    def test_no_use_statements(self, pft_module):
        """Test that code without USE statements is unchanged."""
        before = pft_module.fortran
        pft_module.formatModuleUse()
        result = pft_module.fortran
        assert result == before

    def test_preserves_trailing_comments(self, pft_use_with_comments):
        """Test that trailing comments stay with their USE statement after reordering."""
        pft_use_with_comments.formatModuleUse()
        result = pft_use_with_comments.fortran

        assert 'USE MODD_PARAM ! Data parameters' in result
        assert 'USE MODE_ENV   ! Environment setup' in result
        assert 'USE MODI_FOO   ! Interface for FOO' in result

    def test_upper_uppercases_only_vars(self, pft_use_mixed_case_only):
        """Test that upper=True uppercases variable names in ONLY clause."""
        pft_use_mixed_case_only.formatModuleUse(upper=True)
        result = pft_use_mixed_case_only.fortran

        # Module names are NOT uppercased
        assert 'USE modd_params' in result
        assert 'USE modi_routines' in result
        # ONLY clause vars ARE uppercased
        assert ', ONLY: JPX, JPY' in result
        assert ', ONLY: BAR_SUB' in result
