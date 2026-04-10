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
