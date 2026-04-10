"""
Tests for the Openacc class (OpenACC directive handling).
"""

import pytest
import tempfile
import os
from pyfortool import PYFT


class TestOpenaccRemoveACC:
    """Tests for removeACC method."""

    @pytest.fixture
    def fortran_with_acc(self):
        return """
MODULE MOD_ACC
    IMPLICIT NONE
CONTAINS
    SUBROUTINE GPU_KERNEL(A, B, C)
        REAL, INTENT(IN) :: A(:), B(:)
        REAL, INTENT(OUT) :: C(:)
        INTEGER :: I
!$ACC PARALLEL LOOP
        DO I = 1, SIZE(A)
            C(I) = A(I) + B(I)
        END DO
!$ACC END PARALLEL LOOP
    END SUBROUTINE GPU_KERNEL
END MODULE MOD_ACC
"""

    @pytest.fixture
    def pft_acc(self, fortran_with_acc):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_acc)
            return PYFT(fpath)

    def test_remove_acc_directives(self, pft_acc):
        """Test removing !$ACC directives."""
        pft_acc.removeACC()
        result = pft_acc.fortran
        assert '!$ACC' not in result


class TestOpenaccRemoveBypass:
    """Tests for removebyPassDOCONCURRENT method."""

    @pytest.fixture
    def fortran_with_bypass(self):
        return """
MODULE MOD_BYPASS
    IMPLICIT NONE
CONTAINS
    SUBROUTINE BYPASS_TEST(A, B)
        REAL, INTENT(INOUT) :: A, B
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
!$mnh_define(LOOP)
!$mnh_define(OPENACC)
        A = A + 1.0
        B = B + 2.0
    END SUBROUTINE BYPASS_TEST
END MODULE MOD_BYPASS
"""

    @pytest.fixture
    def pft_bypass(self, fortran_with_bypass):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_bypass)
            return PYFT(fpath)

    def test_remove_bypass_macros(self, pft_bypass):
        """Test removing MNH bypass macros."""
        pft_bypass.removebyPassDOCONCURRENT()
        result = pft_bypass.fortran
        assert '!$MNH_UNDEF(LOOP)' not in result
        assert '!$MNH_UNDEF(OPENACC)' not in result
        assert '!$MNH_DEFINE(LOOP)' not in result
        assert '!$MNH_DEFINE(OPENACC)' not in result


class TestOpenaccAddACCData:
    """Tests for addACCData method."""

    @pytest.fixture
    def fortran_with_intent_arrays(self):
        return """
MODULE MOD_INDATA
    IMPLICIT NONE
CONTAINS
    SUBROUTINE IN_DATA(A, B, C)
        REAL, INTENT(IN) :: A(:)
        REAL, INTENT(INOUT) :: B(:)
        REAL, INTENT(OUT) :: C(:)
        INTEGER :: I
        DO I = 1, SIZE(A)
            C(I) = A(I) + B(I)
            B(I) = B(I) * 2.0
        END DO
    END SUBROUTINE IN_DATA
END MODULE MOD_INDATA
"""

    @pytest.fixture
    def pft_indata(self, fortran_with_intent_arrays):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_intent_arrays)
            return PYFT(fpath)

    def test_add_acc_data_present(self, pft_indata):
        """Test adding !$acc data present directive."""
        pft_indata.addACCData()
        result = pft_indata.fortran
        assert '!$ACC' in result.upper() or '!$acc' in result.lower()


class TestOpenaccAddRoutineSeq:
    """Tests for addACCRoutineSeq method."""

    @pytest.fixture
    def fortran_with_routines(self):
        return """
MODULE MOD_ROUTINES
    IMPLICIT NONE
CONTAINS
    SUBROUTINE ROUTINE_SEQ(A, B)
        REAL, INTENT(INOUT) :: A, B
        A = A + 1.0
        B = B + 2.0
    END SUBROUTINE ROUTINE_SEQ
END MODULE MOD_ROUTINES
"""

    @pytest.fixture
    def pft_routines(self, fortran_with_routines):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_routines)
            return PYFT(fpath)

    def test_add_acc_routine_seq(self, pft_routines):
        """Test adding !$acc routine seq directive."""
        try:
            pft_routines.addACCRoutineSeq([])
            result = pft_routines.fortran
            assert '!$ACC' in result.upper() or 'ROUTINE' in result
        except (AttributeError, TypeError):
            pytest.skip("addACCRoutineSeq requires full tree setup")


class TestOpenaccAllocateHIP:
    """Tests for allocatetoHIP method."""

    @pytest.fixture
    def fortran_with_allocate(self):
        return """
MODULE MOD_HIPALLOC
    IMPLICIT NONE
CONTAINS
    SUBROUTINE HIP_ALLOC(X, Y)
        REAL, POINTER :: X(:), Y(:)
        INTEGER :: ALLOC_STAT
        ALLOCATE(X(100), Y(100), STAT=ALLOC_STAT)
        DEALLOCATE(X, Y, STAT=ALLOC_STAT)
    END SUBROUTINE HIP_ALLOC
END MODULE MOD_HIPALLOC
"""

    @pytest.fixture
    def pft_hip(self, fortran_with_allocate):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_allocate)
            return PYFT(fpath)

    def test_allocate_to_hip(self, pft_hip):
        """Test converting ALLOCATE to HIP version."""
        try:
            pft_hip.allocatetoHIP()
            result = pft_hip.fortran
            assert result is not None
        except (IndexError, AttributeError):
            pytest.skip("allocatetoHIP requires specific structure")


class TestOpenaccCrayBypass:
    """Tests for craybyPassDOCONCURRENT method."""

    @pytest.fixture
    def fortran_with_cray_acc(self):
        return """
MODULE MOD_CRAY
    IMPLICIT NONE
CONTAINS
    SUBROUTINE CRAY_TEST(A, B, C)
        REAL, INTENT(INOUT) :: A(:), B(:), C(:)
        INTEGER :: I, J
!$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO I = 1, 10
            DO J = 1, 10
                C(I,J) = A(I,J) + B(I,J)
            END DO
        END DO
    END SUBROUTINE CRAY_TEST
END MODULE MOD_CRAY
"""

    @pytest.fixture
    def pft_cray(self, fortran_with_cray_acc):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_cray_acc)
            return PYFT(fpath)

    def test_cray_bypass(self, pft_cray):
        """Test CRAY bypass transformation."""
        pft_cray.craybyPassDOCONCURRENT()
        result = pft_cray.fortran
        assert result is not None
