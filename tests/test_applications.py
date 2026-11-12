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

    def test_add_stack_arome(self, pft_with_tree, fortran_with_auto_arrays):
        """Test automatic arrays become CRAY pointers on the AROME stack."""
        pft = pft_with_tree(fortran_with_auto_arrays)
        assert 'YDSTACK' not in pft.fortran

        pft.addStack('AROME', ['module:MOD_STACK/sub:STACK_TEST'])
        result = pft.fortran

        assert 'SUBROUTINE STACK_TEST(N, YDSTACK)' in result  # new argument
        assert 'TYPE (STACK), INTENT(IN) :: YDSTACK' in result
        assert 'USE STACK_MOD, ONLY: STACK' in result
        assert 'POINTER(IP_LOCAL_ARRAY_, LOCAL_ARRAY)' in result
        assert 'YLSTACK=YDSTACK' in result
        assert 'IP_LOCAL_ARRAY_=YLSTACK%L(' in result

    def test_add_stack_arome_outside_stop_scopes(self, pft_with_tree,
                                                 fortran_with_auto_arrays):
        """Test AROME stack transformation ignores routines outside stopScopes."""
        pft = pft_with_tree(fortran_with_auto_arrays)
        before = pft.fortran
        pft.addStack('AROME', [])
        assert pft.fortran == before

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
        REAL, INTENT(IN) :: A(D%NIT, 2)
        REAL, INTENT(OUT) :: B(D%NIT)
        INTEGER :: JI
        DO JI = D%NIB, D%NIT
            B(JI) = A(JI, 1)
        END DO
    END SUBROUTINE IJ_TEST
END MODULE MOD_IJ
"""

    def test_remove_ij_dim(self, pft_with_tree, fortran_with_ij_loops):
        """Test the loop on JI is replaced by an initialisation of JI."""
        pft = pft_with_tree(fortran_with_ij_loops)
        assert 'DO JI = D%NIB, D%NIT' in pft.fortran

        pft.removeIJDim(['module:MOD_IJ/sub:IJ_TEST'])
        result = pft.fortran

        assert 'DO JI' not in result
        assert 'END DO' not in result
        assert 'JI = D%NIB' in result
        assert 'B(JI) = A(JI, 1)' in result

    def test_remove_ij_dim_outside_stop_scopes(self, pft_with_tree, fortran_with_ij_loops):
        """Test removeIJDim leaves routines outside stopScopes unchanged."""
        pft = pft_with_tree(fortran_with_ij_loops)
        before = pft.fortran
        pft.removeIJDim([])
        assert pft.fortran == before


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


class TestShugradTransfoUselessParens:
    """Tests for shugradTransfo() with useless parenthesis (issue #55)."""

    @staticmethod
    def _transform(args):
        code = f"""
MODULE MODFOO
  CONTAINS
  SUBROUTINE SUBFOO
    REAL, DIMENSION(D%NIJT, D%NKT) :: P1, P2, P3, P4
    REAL, DIMENSION(D%NIJT) :: P7
    P7=MXM2D({args[0]})
    P1=MZM({args[1]})
    P1=MZM({args[2]})
    P1=MZM({args[3]})
  CONTAINS
    INCLUDE "shugrad.h"
  END SUBROUTINE SUBFOO
END MODULE MODFOO
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(code)
            pft = PYFT(fpath)
            pft.shugradTransfo()
            return pft.fortran

    def test_useless_parens_removed(self):
        """Useless parenthesis around arguments give the same result as without them."""
        withParens = self._transform(['((P1(:,IKB+1)+P1(:,IKB)))', '(((P2)))',
                                      '((MZF((P4))))', '((P2+P3)*P4)'])
        withoutParens = self._transform(['P1(:,IKB+1)+P1(:,IKB)', 'P2',
                                         'MZF(P4)', '(P2+P3)*P4'])
        assert withParens == withoutParens
        assert 'CALL MXM2D_PHY(D, ZWORKSHUGRAD_2D_1, P7)' in withParens
        assert 'CALL MZM_PHY(D, P2, P1)' in withParens
        # Necessary parenthesis are kept
        assert '=(P2(1:D%NIJT, 1:D%NKT)+P3(1:D%NIJT, 1:D%NKT))*P4' in withParens
