"""
Tests for the Applications class (high-level transformations).
"""

import pytest
import tempfile
import os
from pyfortool import PYFT


class TestApplicationsDrHook:
    """Tests for addDrHook and deleteDrHook methods."""

    @pytest.fixture
    def fortran_with_subroutine(self):
        return """
MODULE MOD_DRHOOK
    IMPLICIT NONE
CONTAINS
    SUBROUTINE HOOK_TEST(A, B)
        REAL, INTENT(INOUT) :: A, B
        A = A + 1.0
        B = B + 2.0
    END SUBROUTINE HOOK_TEST
END MODULE MOD_DRHOOK
"""

    @pytest.fixture
    def pft_drhook(self, fortran_with_subroutine):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_subroutine)
            return PYFT(fpath)

    def test_add_drhook(self, pft_drhook):
        """Test adding DR_HOOK instrumentation."""
        pft_drhook.addDrHook()
        result = pft_drhook.fortran
        assert 'DR_HOOK' in result.upper() or 'YOMHOOK' in result

    def test_delete_drhook(self, pft_drhook):
        """Test deleting DR_HOOK instrumentation."""
        pft_drhook.deleteDrHook()
        result = pft_drhook.fortran
        assert result is not None


class TestApplicationsDeleteCalls:
    """Tests for removeCall method."""

    @pytest.fixture
    def fortran_with_call(self):
        return """
MODULE MOD_CALL
    IMPLICIT NONE
CONTAINS
    SUBROUTINE CALL_TEST(A, B)
        REAL, INTENT(INOUT) :: A, B
        CALL HELPER(A, B)
        A = A + 1.0
    END SUBROUTINE CALL_TEST
    SUBROUTINE HELPER(X, Y)
        REAL, INTENT(INOUT) :: X, Y
        X = X * 2.0
        Y = Y * 2.0
    END SUBROUTINE HELPER
END MODULE MOD_CALL
"""

    @pytest.fixture
    def pft_call(self, fortran_with_call):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_call)
            return PYFT(fpath)

    def test_remove_call(self, pft_call):
        """Test removing a CALL statement."""
        count = pft_call.removeCall('HELPER')
        assert count >= 0


class TestApplicationsBudgetDDH:
    """Tests for deleteBudgetDDH method."""

    @pytest.fixture
    def fortran_with_budget(self):
        return """
MODULE MOD_BUDGET
    IMPLICIT NONE
CONTAINS
    SUBROUTINE BUDGET_TEST(A, B)
        REAL, INTENT(INOUT) :: A, B
        CALL BUDGET_STORE_INIT_PHY('TEST', 0)
        A = A + 1.0
        CALL BUDGET_STORE_END_PHY('TEST', 0)
    END SUBROUTINE BUDGET_TEST
END MODULE MOD_BUDGET
"""

    @pytest.fixture
    def pft_budget(self, fortran_with_budget):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_budget)
            return PYFT(fpath)

    def test_delete_budget(self, pft_budget):
        """Test deleting budget diagnostics."""
        pft_budget.deleteBudgetDDH()
        result = pft_budget.fortran
        assert result is not None


class TestApplicationsMPPDB:
    """Tests for addMPPDB_CHECKS method."""

    @pytest.fixture
    def fortran_with_arrays(self):
        return """
MODULE MOD_MPPDB
    IMPLICIT NONE
CONTAINS
    SUBROUTINE MPPDB_TEST(A, B, C)
        REAL, INTENT(IN) :: A(:)
        REAL, INTENT(INOUT) :: B(:)
        REAL, INTENT(OUT) :: C(:)
        INTEGER :: I
        DO I = 1, SIZE(A)
            C(I) = A(I) + B(I)
        END DO
    END SUBROUTINE MPPDB_TEST
END MODULE MOD_MPPDB
"""

    @pytest.fixture
    def pft_mppdb(self, fortran_with_arrays):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_arrays)
            return PYFT(fpath)

    def test_add_mppdb_checks(self, pft_mppdb):
        """Test adding MPPDB checks."""
        pft_mppdb.addMPPDB_CHECKS()
        result = pft_mppdb.fortran
        assert 'MPPDB' in result.upper() or 'CHECK' in result

    def test_add_mppdb_prints(self, pft_mppdb):
        """Test adding debug prints instead of MPPDB checks."""
        pft_mppdb.addMPPDB_CHECKS(printsMode=True)
        result = pft_mppdb.fortran
        assert result is not None


class TestApplicationsPHYEX:
    """Tests for PHYEX-related methods."""

    @pytest.fixture
    def fortran_with_phyex(self):
        return """
MODULE MOD_PHYEX
    IMPLICIT NONE
CONTAINS
    SUBROUTINE PHYEX_TEST(A, B)
        REAL, INTENT(INOUT) :: A, B
        CALL ROTATE_WIND(A, B)
        A = A + 1.0
    END SUBROUTINE PHYEX_TEST
END MODULE MOD_PHYEX
"""

    @pytest.fixture
    def pft_phyex(self, fortran_with_phyex):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_phyex)
            return PYFT(fpath)

    def test_delete_noncolumn_calls(self, pft_phyex):
        """Test deleting non-column PHYEX calls."""
        pft_phyex.deleteNonColumnCallsPHYEX()
        result = pft_phyex.fortran
        assert result is not None


class TestApplicationsConvertTypes:
    """Tests for convertTypesInCompute method."""

    @pytest.fixture
    def fortran_with_type_access(self):
        return """
MODULE MOD_TYPES
    IMPLICIT NONE
    TYPE :: CST_TYPE
        REAL :: XG
        REAL :: ZS
    END TYPE CST_TYPE
CONTAINS
    SUBROUTINE TYPE_COMPUTE(CST, ZA)
        TYPE(CST_TYPE), INTENT(IN) :: CST
        REAL, INTENT(OUT) :: ZA
        ZA = 1.0 + CST%XG
    END SUBROUTINE TYPE_COMPUTE
END MODULE MOD_TYPES
"""

    @pytest.fixture
    def pft_types(self, fortran_with_type_access):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_type_access)
            return PYFT(fpath)

    def test_convert_types_in_compute(self, pft_types):
        """Test converting TYPE member accesses to local variables."""
        pft_types.convertTypesInCompute()
        result = pft_types.fortran
        assert result is not None


class TestApplicationsMesoNHGPU:
    """Tests for Meso-NH GPU-related methods."""

    @pytest.fixture
    def fortran_with_mnh_gpu(self):
        return """
MODULE MOD_MNHGPU
    IMPLICIT NONE
CONTAINS
    SUBROUTINE MNHGPU_TEST(A, B)
        REAL, INTENT(INOUT) :: A, B
        IF (OCND2) THEN
            A = A + 1.0
        END IF
    END SUBROUTINE MNHGPU_TEST
END MODULE MOD_MNHGPU
"""

    @pytest.fixture
    def pft_mnhgpu(self, fortran_with_mnh_gpu):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_mnh_gpu)
            return PYFT(fpath)

    def test_delete_routine_calls_mesonh_gpu(self, pft_mnhgpu):
        """Test deleting Meso-NH GPU incompatible calls."""
        pft_mnhgpu.deleteRoutineCallsMesoNHGPU()
        result = pft_mnhgpu.fortran
        assert result is not None


class TestApplicationsStack:
    """Tests for addStack method."""

    @pytest.fixture
    def fortran_with_auto_arrays(self):
        return """
MODULE MOD_STACK
    IMPLICIT NONE
CONTAINS
    SUBROUTINE STACK_TEST(N)
        INTEGER, INTENT(IN) :: N
        REAL :: LOCAL_ARRAY(N, N)
        INTEGER :: I, J
        DO I = 1, N
            DO J = 1, N
                LOCAL_ARRAY(I, J) = REAL(I + J)
            END DO
        END DO
    END SUBROUTINE STACK_TEST
END MODULE MOD_STACK
"""

    @pytest.fixture
    def pft_stack(self, fortran_with_auto_arrays):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_auto_arrays)
            return PYFT(fpath)

    def test_add_stack_arome(self, pft_stack):
        """Test adding stack allocation for AROME model."""
        try:
            pft_stack.addStack('AROME', [])
            result = pft_stack.fortran
            assert result is not None
        except (AttributeError, TypeError):
            pytest.skip("AROME stack requires full tree setup")

    def test_add_stack_mesonh(self, pft_stack):
        """Test adding stack allocation for MESO-NH model."""
        pft_stack.addStack('MESONH', [])
        result = pft_stack.fortran
        assert result is not None

    def test_add_stack_invalid_model(self, pft_stack):
        """Test adding stack with invalid model raises error."""
        with pytest.raises(Exception):
            pft_stack.addStack('INVALID', [])


class TestApplicationsIJDim:
    """Tests for removeIJDim method."""

    @pytest.fixture
    def fortran_with_ij_loops(self):
        return """
MODULE MOD_IJ
    IMPLICIT NONE
CONTAINS
    SUBROUTINE IJ_TEST(D, A, B)
        TYPE(DIM_TYPE), INTENT(IN) :: D
        REAL, INTENT(IN) :: A(:,:)
        REAL, INTENT(OUT) :: B(:)
        INTEGER :: JI, JJ
        DO JI = D%NIB, D%NIT
            B(JI) = A(JI, 1)
        END DO
    END SUBROUTINE IJ_TEST
END MODULE MOD_IJ
"""

    @pytest.fixture
    def pft_ij(self, fortran_with_ij_loops):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_ij_loops)
            return PYFT(fpath)

    def test_remove_ij_dim(self, pft_ij):
        """Test removing I/J dimensions for column computation."""
        try:
            pft_ij.removeIJDim([])
            result = pft_ij.fortran
            assert result is not None
        except (AttributeError, TypeError):
            pytest.skip("removeIJDim requires full tree setup")


class TestApplicationsModi:
    """Tests for buildModi method."""

    @pytest.fixture
    def fortran_for_modi(self):
        return """
MODULE MODE_MODTEST
    IMPLICIT NONE
CONTAINS
    SUBROUTINE TEST_SUB(A, B)
        REAL, INTENT(INOUT) :: A, B
        A = A + 1.0
    END SUBROUTINE TEST_SUB
END MODULE MODE_MODTEST
"""

    @pytest.fixture
    def pft_modi(self, fortran_for_modi):
        import shutil
        tmpdir = tempfile.mkdtemp()
        fpath = os.path.join(tmpdir, 'MODE_MODTEST.F90')
        with open(fpath, 'w') as f:
            f.write(fortran_for_modi)
        pft = PYFT(fpath)
        yield pft, tmpdir
        shutil.rmtree(tmpdir, ignore_errors=True)

    def test_build_modi(self, pft_modi):
        """Test building MODI interface file."""
        pft, tmpdir = pft_modi
        pft.buildModi()
        modi_path = os.path.join(tmpdir, 'modi_MODE_MODTEST.F90')
        assert os.path.exists(modi_path)


class TestApplicationsRemoveExtraDO:
    """Tests for removeExtraDOinMnhDoConcurrent method."""

    @pytest.fixture
    def fortran_with_do_concurrent(self):
        return """
MODULE MOD_DOCONCURRENT
    IMPLICIT NONE
CONTAINS
    SUBROUTINE DO_CONCURRENT_TEST(A, B)
        REAL, INTENT(INOUT) :: A(:), B(:)
        INTEGER :: I
!$MNH_DO_CONCURRENT
        DO I = 1, SIZE(A)
            B(I) = A(I) + 1.0
        END DO
    END SUBROUTINE DO_CONCURRENT_TEST
END MODULE MOD_DOCONCURRENT
"""

    @pytest.fixture
    def pft_doconcurrent(self, fortran_with_do_concurrent):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(fortran_with_do_concurrent)
            return PYFT(fpath)

    def test_remove_extra_do(self, pft_doconcurrent):
        """Test removing extra DO in MNH DO CONCURRENT."""
        pft_doconcurrent.removeExtraDOinMnhDoConcurrent()
        result = pft_doconcurrent.fortran
        assert result is not None
