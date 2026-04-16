"""
Tests for the Cpp class (C preprocessor directive handling).
"""

import pytest
import tempfile
import os
from pyfortool import PYFT


class TestCppApplyIfdef:
    """Tests for applyCPPifdef method."""

    @pytest.fixture
    def fortran_with_ifdef(self):
        return """
MODULE MOD_CPP
    IMPLICIT NONE
CONTAINS
    SUBROUTINE CPP_TEST(X, Y)
        REAL, INTENT(INOUT) :: X, Y
#ifdef KEY1
        X = X * 2.0
#else
        Y = Y * 3.0
#endif
    END SUBROUTINE CPP_TEST
END MODULE MOD_CPP
"""

    @pytest.fixture
    def pft_ifdef(self, fortran_with_ifdef):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_ifdef)
            return PYFT(fpath)

    def test_apply_cpp_ifdef_key_present(self, pft_ifdef):
        """Test applying #ifdef with key present."""
        pft_ifdef.applyCPPifdef(['KEY1'])
        result = pft_ifdef.fortran
        assert '2.0' in result
        assert '3.0' not in result

    def test_apply_cpp_ifdef_key_absent(self, pft_ifdef):
        """Test applying #ifdef with key absent."""
        pft_ifdef.applyCPPifdef([])
        result = pft_ifdef.fortran
        assert result is not None

    def test_apply_cpp_ifdef_key_false(self, pft_ifdef):
        """Test applying #ifdef with key marked as false (%KEY1)."""
        pft_ifdef.applyCPPifdef(['%KEY1'])
        result = pft_ifdef.fortran
        assert '2.0' not in result
        assert '3.0' in result


class TestCppIfndef:
    """Tests for #ifndef handling."""

    @pytest.fixture
    def fortran_with_ifndef(self):
        return """
MODULE MOD_IFNDEF
    IMPLICIT NONE
CONTAINS
    SUBROUTINE IFNDEF_TEST(A, B)
        REAL, INTENT(INOUT) :: A, B
#ifndef DISABLE_FEATURE
        A = A + 1.0
#endif
        B = B + 2.0
    END SUBROUTINE IFNDEF_TEST
END MODULE MOD_IFNDEF
"""

    @pytest.fixture
    def pft_ifndef(self, fortran_with_ifndef):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_ifndef)
            return PYFT(fpath)

    def test_apply_cpp_ifndef_not_defined(self, pft_ifndef):
        """Test #ifndef when key is not defined."""
        pft_ifndef.applyCPPifdef([])
        result = pft_ifndef.fortran
        assert '+ 1.0' in result

    def test_apply_cpp_ifndef_defined(self, pft_ifndef):
        """Test #ifndef when key is defined."""
        pft_ifndef.applyCPPifdef(['DISABLE_FEATURE'])
        result = pft_ifndef.fortran
        assert '+ 1.0' not in result


class TestCppNested:
    """Tests for nested #ifdef handling."""

    @pytest.fixture
    def fortran_with_nested_ifdef(self):
        return """
MODULE MOD_NESTED
    IMPLICIT NONE
CONTAINS
    SUBROUTINE NESTED_TEST(X, Y, Z)
        REAL, INTENT(INOUT) :: X, Y, Z
#ifdef OUTER
#ifdef INNER
        X = 1.0
#else
        Y = 2.0
#endif
#else
        Z = 3.0
#endif
    END SUBROUTINE NESTED_TEST
END MODULE MOD_NESTED
"""

    @pytest.fixture
    def pft_nested(self, fortran_with_nested_ifdef):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_nested_ifdef)
            return PYFT(fpath)

    def test_apply_nested_both_present(self, pft_nested):
        """Test nested #ifdef with both keys present."""
        pft_nested.applyCPPifdef(['OUTER', 'INNER'])
        result = pft_nested.fortran
        assert '= 1.0' in result
        assert '= 2.0' not in result
        assert '= 3.0' not in result

    def test_apply_nested_only_outer(self, pft_nested):
        """Test nested #ifdef with only outer key present."""
        pft_nested.applyCPPifdef(['OUTER'])
        result = pft_nested.fortran
        assert result is not None


class TestCppElse:
    """Tests for #else handling in applyCPPifdef."""

    @pytest.fixture
    def fortran_with_else(self):
        return """
MODULE MOD_ELSE
    IMPLICIT NONE
CONTAINS
    SUBROUTINE ELSE_TEST(A, B)
        REAL, INTENT(INOUT) :: A, B
#ifdef FEATURE
        A = 10.0
#else
        B = 20.0
#endif
    END SUBROUTINE ELSE_TEST
END MODULE MOD_ELSE
"""

    @pytest.fixture
    def pft_else(self, fortran_with_else):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_else)
            return PYFT(fpath)

    def test_else_takes_true_branch(self, pft_else):
        """Test that #else takes true branch when key is present."""
        pft_else.applyCPPifdef(['FEATURE'])
        result = pft_else.fortran
        assert '= 10.0' in result
        assert '= 20.0' not in result

    def test_else_takes_false_branch(self, pft_else):
        """Test that #else takes false branch when key is absent."""
        pft_else.applyCPPifdef([])
        result = pft_else.fortran
        assert result is not None


class TestCppMultipleBlocks:
    """Tests for multiple #ifdef blocks."""

    @pytest.fixture
    def fortran_with_multiple_ifdef(self):
        return """
MODULE MOD_MULTI
    IMPLICIT NONE
CONTAINS
    SUBROUTINE MULTI_TEST(A, B, C)
        REAL, INTENT(INOUT) :: A, B, C
#ifdef FLAG1
        A = 1.0
#endif
#ifdef FLAG2
        B = 2.0
#endif
        C = 3.0
    END SUBROUTINE MULTI_TEST
END MODULE MOD_MULTI
"""

    @pytest.fixture
    def pft_multi(self, fortran_with_multiple_ifdef):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_multiple_ifdef)
            return PYFT(fpath)

    def test_apply_multiple_flags(self, pft_multi):
        """Test applying multiple #ifdef flags."""
        pft_multi.applyCPPifdef(['FLAG1', 'FLAG2'])
        result = pft_multi.fortran
        assert '= 1.0' in result
        assert '= 2.0' in result
        assert '= 3.0' in result

    def test_apply_partial_flags(self, pft_multi):
        """Test applying only some #ifdef flags."""
        pft_multi.applyCPPifdef(['FLAG1'])
        result = pft_multi.fortran
        assert result is not None
