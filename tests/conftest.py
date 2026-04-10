"""
Shared pytest fixtures for PyForTool tests.
"""

import pytest
import tempfile
import os
from pathlib import Path


@pytest.fixture
def simple_fortran():
    """Simple FORTRAN code for basic tests."""
    return """
PROGRAM TEST
    INTEGER :: X, Y
    X = 1
    Y = X + 2
    PRINT *, Y
END PROGRAM TEST
"""


@pytest.fixture
def module_with_subroutine():
    """Module with a subroutine."""
    return """
MODULE MOD_TEST
    IMPLICIT NONE
CONTAINS
    SUBROUTINE SUB(X, Y, Z)
        REAL, INTENT(IN) :: X(:)
        REAL, INTENT(IN) :: Y(:)
        REAL, INTENT(OUT) :: Z(:)
        INTEGER :: I
        DO I = 1, SIZE(X)
            Z(I) = X(I) + Y(I)
        END DO
    END SUBROUTINE SUB
END MODULE MOD_TEST
"""


@pytest.fixture
def module_with_functions():
    """Module with subroutines and functions."""
    return """
MODULE MOD_FUNC
    IMPLICIT NONE
CONTAINS
    SUBROUTINE SUB1(A, B)
        REAL, INTENT(IN) :: A
        REAL, INTENT(OUT) :: B
        B = A * 2.0
    END SUBROUTINE SUB1

    FUNCTION FUNC1(X) RESULT(Y)
        REAL, INTENT(IN) :: X
        REAL :: Y
        Y = X + 1.0
    END FUNCTION FUNC1

    SUBROUTINE SUB2(C, D)
        REAL, INTENT(INOUT) :: C
        REAL, INTENT(OUT) :: D
        D = C ** 2
    END SUBROUTINE SUB2
END MODULE MOD_FUNC
"""


@pytest.fixture
def fortran_with_arrays():
    """FORTRAN code with array syntax."""
    return """
MODULE MOD_ARRAY
    IMPLICIT NONE
CONTAINS
    SUBROUTINE ARRAY_OPS(A, B, C)
        REAL, INTENT(IN) :: A(:), B(:)
        REAL, INTENT(OUT) :: C(:)
        C(:) = A(:) + B(:)
        C(:) = C(:) * 2.0
    END SUBROUTINE ARRAY_OPS
END MODULE MOD_ARRAY
"""


@pytest.fixture
def fortran_with_conditionals():
    """FORTRAN code with IF statements."""
    return """
MODULE MOD_IF
    IMPLICIT NONE
CONTAINS
    SUBROUTINE CHECK_VAL(X, FLAG)
        REAL, INTENT(IN) :: X
        LOGICAL, INTENT(OUT) :: FLAG
        IF (X > 0.0) THEN
            FLAG = .TRUE.
        ELSE
            FLAG = .FALSE.
        END IF
    END SUBROUTINE CHECK_VAL

    SUBROUTINE MULTI_IF(Y, Z)
        REAL, INTENT(IN) :: Y
        REAL, INTENT(OUT) :: Z
        IF (Y > 10.0) THEN
            Z = 100.0
        ELSE IF (Y > 5.0) THEN
            Z = 50.0
        ELSE
            Z = 10.0
        END IF
    END SUBROUTINE MULTI_IF
END MODULE MOD_IF
"""


@pytest.fixture
def fortran_with_calls():
    """FORTRAN code with CALL statements."""
    return """
MODULE MOD_CALLS
    IMPLICIT NONE
CONTAINS
    SUBROUTINE HELPER(A, B)
        REAL, INTENT(INOUT) :: A
        REAL, INTENT(OUT) :: B
        A = A + 1.0
        B = A * 2.0
    END SUBROUTINE HELPER

    SUBROUTINE MAIN_SUB(X, Y)
        REAL, INTENT(INOUT) :: X
        REAL, INTENT(OUT) :: Y
        CALL HELPER(X, Y)
    END SUBROUTINE MAIN_SUB
END MODULE MOD_CALLS
"""


@pytest.fixture
def fortran_with_comments():
    """FORTRAN code with comments."""
    return """
! This is a header comment
MODULE MOD_COMMENT
    IMPLICIT NONE
CONTAINS
    SUBROUTINE WITH_COMMENTS(A, B)
        REAL, INTENT(IN) :: A  ! Input value
        REAL, INTENT(OUT) :: B
        ! Perform calculation
        B = A * 2.0  ! Multiply by 2
    END SUBROUTINE WITH_COMMENTS
END MODULE MOD_COMMENT
"""


@pytest.fixture
def fortran_with_types():
    """FORTRAN code with derived types."""
    return """
MODULE MOD_TYPE
    IMPLICIT NONE
    TYPE :: MY_TYPE
        REAL :: X
        REAL :: Y
        INTEGER :: I
    END TYPE MY_TYPE
CONTAINS
    SUBROUTINE TYPE_OPS(OBJ, VAL)
        TYPE(MY_TYPE), INTENT(INOUT) :: OBJ
        REAL, INTENT(IN) :: VAL
        OBJ%X = VAL
        OBJ%Y = VAL * 2.0
        OBJ%I = INT(VAL)
    END SUBROUTINE TYPE_OPS
END MODULE MOD_TYPE
"""


@pytest.fixture
def fortran_with_use():
    """FORTRAN code with USE statements."""
    return """
MODULE MOD_USE
    USE MOD_BASE, ONLY: BASE_VAL
    IMPLICIT NONE
CONTAINS
    SUBROUTINE USE_EXAMPLE(A)
        REAL, INTENT(IN) :: A
        REAL :: B
        B = A + BASE_VAL
    END SUBROUTINE USE_EXAMPLE
END MODULE MOD_USE
"""


@pytest.fixture
def fortran_with_loops():
    """FORTRAN code with various loops."""
    return """
MODULE MOD_LOOPS
    IMPLICIT NONE
CONTAINS
    SUBROUTINE DO_LOOP(A, B, N)
        INTEGER, INTENT(IN) :: N
        REAL, INTENT(IN) :: A(N)
        REAL, INTENT(OUT) :: B(N)
        INTEGER :: I
        DO I = 1, N
            B(I) = A(I) * 2.0
        END DO
    END SUBROUTINE DO_LOOP

    SUBROUTINE NESTED_LOOP(A, B, M, N)
        INTEGER, INTENT(IN) :: M, N
        REAL, INTENT(IN) :: A(M, N)
        REAL, INTENT(OUT) :: B(M, N)
        INTEGER :: I, J
        DO I = 1, M
            DO J = 1, N
                B(I, J) = A(I, J) + 1.0
            END DO
        END DO
    END SUBROUTINE NESTED_LOOP
END MODULE MOD_LOOPS
"""


@pytest.fixture
def fortran_with_declarations():
    """FORTRAN code with various variable declarations."""
    return """
MODULE MOD_DECL
    IMPLICIT NONE
    REAL, PARAMETER :: PI = 3.14159
    INTEGER, DIMENSION(10) :: ARRAY1
    REAL, DIMENSION(5, 5) :: ARRAY2
CONTAINS
    SUBROUTINE DECL_EXAMPLE()
        REAL :: LOCAL_VAL
        INTEGER :: I
        DO I = 1, 10
            ARRAY1(I) = I
        END DO
    END SUBROUTINE DECL_EXAMPLE
END MODULE MOD_DECL
"""


@pytest.fixture
def fortran_with_where():
    """FORTRAN code with WHERE constructs."""
    return """
MODULE MOD_WHERE
    IMPLICIT NONE
CONTAINS
    SUBROUTINE WHERE_EXAMPLE(A, B, C)
        REAL, INTENT(IN) :: A(:), B(:)
        REAL, INTENT(OUT) :: C(:)
        WHERE (A(:) > 0.0)
            C(:) = A(:) + B(:)
        ELSEWHERE
            C(:) = 0.0
        END WHERE
    END SUBROUTINE WHERE_EXAMPLE
END MODULE MOD_WHERE
"""


@pytest.fixture
def fortran_with_if_single():
    """FORTRAN code with single-line IF statements."""
    return """
MODULE MOD_SINGLE_IF
    IMPLICIT NONE
CONTAINS
    SUBROUTINE SINGLE_IF(X)
        REAL, INTENT(INOUT) :: X
        IF (X > 0.0) X = X * 2.0
    END SUBROUTINE SINGLE_IF
END MODULE MOD_SINGLE_IF
"""


class TempFortranFile:
    """Context manager for creating temporary FORTRAN files."""

    def __init__(self, content, suffix='.F90'):
        self.content = content
        self.suffix = suffix
        self.path = None
        self._tmp_dir = None

    def __enter__(self):
        self._tmp_dir = tempfile.TemporaryDirectory()
        self.path = os.path.join(self._tmp_dir.name, f'test{self.suffix}')
        with open(self.path, 'w') as f:
            f.write(self.content)
        return self.path

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._tmp_dir:
            self._tmp_dir.cleanup()


@pytest.fixture
def tmp_fortran_file(simple_fortran):
    """Create a temporary FORTRAN file with simple content."""
    return TempFortranFile(simple_fortran)


@pytest.fixture
def tmp_module_file(module_with_subroutine):
    """Create a temporary FORTRAN file with module and subroutine."""
    return TempFortranFile(module_with_subroutine)


@pytest.fixture
def tmp_array_file(fortran_with_arrays):
    """Create a temporary FORTRAN file with array syntax."""
    return TempFortranFile(fortran_with_arrays)


@pytest.fixture
def tmp_if_file(fortran_with_conditionals):
    """Create a temporary FORTRAN file with IF statements."""
    return TempFortranFile(fortran_with_conditionals)


@pytest.fixture
def tmp_calls_file(fortran_with_calls):
    """Create a temporary FORTRAN file with CALL statements."""
    return TempFortranFile(fortran_with_calls)


@pytest.fixture
def tmp_comment_file(fortran_with_comments):
    """Create a temporary FORTRAN file with comments."""
    return TempFortranFile(fortran_with_comments)


@pytest.fixture
def tmp_type_file(fortran_with_types):
    """Create a temporary FORTRAN file with types."""
    return TempFortranFile(fortran_with_types)


@pytest.fixture
def tmp_loop_file(fortran_with_loops):
    """Create a temporary FORTRAN file with loops."""
    return TempFortranFile(fortran_with_loops)


@pytest.fixture
def tmp_where_file(fortran_with_where):
    """Create a temporary FORTRAN file with WHERE."""
    return TempFortranFile(fortran_with_where)


@pytest.fixture
def tmp_single_if_file(fortran_with_if_single):
    """Create a temporary FORTRAN file with single-line IF."""
    return TempFortranFile(fortran_with_if_single)


@pytest.fixture
def tmp_decl_file(fortran_with_declarations):
    """Create a temporary FORTRAN file with declarations."""
    return TempFortranFile(fortran_with_declarations)


def pft_from_fixture(fixture_name, tmp_path_factory, request):
    """Helper to create PYFT instance from fixture content."""
    from pyfortool import PYFT

    content = request.getfixturevalue(fixture_name.replace('tmp_', '').replace('_file', ''))
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(content)
        return PYFT(fpath)


@pytest.fixture
def pft_simple(simple_fortran):
    """PYFT instance from simple FORTRAN code."""
    from pyfortool import PYFT
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(simple_fortran)
        return PYFT(fpath)


@pytest.fixture
def pft_module(module_with_subroutine):
    """PYFT instance from module with subroutine."""
    from pyfortool import PYFT
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(module_with_subroutine)
        return PYFT(fpath)


@pytest.fixture
def pft_arrays(fortran_with_arrays):
    """PYFT instance from code with arrays."""
    from pyfortool import PYFT
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(fortran_with_arrays)
        return PYFT(fpath)


@pytest.fixture
def pft_calls(fortran_with_calls):
    """PYFT instance from code with CALL statements."""
    from pyfortool import PYFT
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(fortran_with_calls)
        return PYFT(fpath)


@pytest.fixture
def pft_comments(fortran_with_comments):
    """PYFT instance from code with comments."""
    from pyfortool import PYFT
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(fortran_with_comments)
        return PYFT(fpath)


@pytest.fixture
def pft_loops(fortran_with_loops):
    """PYFT instance from code with loops."""
    from pyfortool import PYFT
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(fortran_with_loops)
        return PYFT(fpath)


@pytest.fixture
def pft_if(fortran_with_conditionals):
    """PYFT instance from code with IF statements."""
    from pyfortool import PYFT
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(fortran_with_conditionals)
        return PYFT(fpath)


@pytest.fixture
def pft_single_if(fortran_with_if_single):
    """PYFT instance from code with single-line IF."""
    from pyfortool import PYFT
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(fortran_with_if_single)
        return PYFT(fpath)


@pytest.fixture
def pft_where(fortran_with_where):
    """PYFT instance from code with WHERE."""
    from pyfortool import PYFT
    with tempfile.TemporaryDirectory() as tmpdir:
        fpath = os.path.join(tmpdir, 'test.F90')
        with open(fpath, 'w') as f:
            f.write(fortran_with_where)
        return PYFT(fpath)
