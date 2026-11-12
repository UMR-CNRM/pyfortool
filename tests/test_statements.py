"""
Tests for Statements mixin class in pyfortool.statements.
"""

import pytest
import tempfile
import os
from pyfortool import PYFT


class TestRemoveCall:
    """Tests for removeCall() method."""

    def test_remove_call_returns_int(self, pft_calls):
        """Test removeCall() returns count."""
        result = pft_calls.removeCall('HELPER')
        assert isinstance(result, int)
        assert result >= 0

    def test_remove_nonexistent_call(self, pft_calls):
        """Test removing nonexistent call."""
        result = pft_calls.removeCall('NONEXISTENT')
        assert result == 0

    def test_remove_call_with_simplify(self, pft_calls):
        """Test removeCall() with simplify."""
        result = pft_calls.removeCall('HELPER', simplify=True)
        assert isinstance(result, int)

    def test_remove_call_modifies_code(self, pft_calls):
        """Test removeCall() actually modifies code."""
        original = pft_calls.fortran
        pft_calls.removeCall('HELPER')
        modified = pft_calls.fortran
        assert 'HELPER' not in modified or 'CALL HELPER' not in modified


class TestRemovePrints:
    """Tests for removePrints() method."""

    def test_remove_prints_runs(self, pft_simple):
        """Test removePrints() runs without error."""
        pft_simple.removePrints()

    def test_remove_prints_no_prints(self, pft_calls):
        """Test removePrints() on code without prints."""
        pft_calls.removePrints()


class TestRemoveArraySyntax:
    """Tests for removeArraySyntax() method."""

    def test_remove_array_syntax_runs(self, pft_arrays):
        """Test removeArraySyntax() runs."""
        pft_arrays.removeArraySyntax()
        # Should complete without error

    def test_remove_array_syntax_concurrent(self, pft_arrays):
        """Test with concurrent=True."""
        pft_arrays.removeArraySyntax(concurrent=True)

    def test_remove_array_syntax_everywhere(self, pft_arrays):
        """Test with everywhere parameter."""
        pft_arrays.removeArraySyntax(everywhere=True)

    def test_remove_array_syntax_no_everywhere(self, pft_arrays):
        """Test with everywhere=False."""
        pft_arrays.removeArraySyntax(everywhere=False)


class TestSetFalseIfStmt:
    """Tests for setFalseIfStmt() method."""

    def test_set_false_if_stmt_runs(self, pft_if):
        """Test setFalseIfStmt() runs."""
        pft_if.setFalseIfStmt('LFLAG')

    def test_set_false_if_stmt_list(self, pft_if):
        """Test with list of flags."""
        pft_if.setFalseIfStmt(['FLAG1', 'FLAG2'])

    def test_set_false_if_stmt_with_simplify(self, pft_if):
        """Test with simplify."""
        pft_if.setFalseIfStmt('FLAG', simplify=True)


class TestChangeIfStatementsInIfConstructs:
    """Tests for changeIfStatementsInIfConstructs() method."""

    def test_change_single_if_runs(self, pft_single_if):
        """Test on single-line IF."""
        pft_single_if.changeIfStatementsInIfConstructs()

    def test_change_all_if_runs(self, pft_if):
        """Test on all IF statements."""
        pft_if.changeIfStatementsInIfConstructs()


class TestEmpty:
    """Tests for empty() method."""

    def test_empty_runs(self, pft_calls):
        """Test empty() runs."""
        pft_calls.empty()


class TestIsNodeInProcedure:
    """Tests for isNodeInProcedure() method."""

    def test_is_node_in_procedure_runs(self, pft_calls):
        """Test isNodeInProcedure() runs."""
        scope = pft_calls.getScopeNode('module:MOD_CALLS/sub:HELPER')
        nodes = scope.findall('.//{*}named-E')
        if nodes:
            result = scope.isNodeInProcedure(nodes[0], ['ALLOCATED'])
            assert isinstance(result, bool)


class TestIsNodeInCall:
    """Tests for isNodeInCall() method."""

    def test_is_node_in_call_runs(self, pft_calls):
        """Test isNodeInCall() runs."""
        scope = pft_calls.getScopeNode('module:MOD_CALLS/sub:MAIN_SUB')
        nodes = scope.findall('.//{*}named-E')
        if nodes:
            result = scope.isNodeInCall(nodes[0])
            assert isinstance(result, bool)


class TestCheckContinuationInArraySection:
    """Tests for checkContinuationInArraySection() method."""

    @staticmethod
    def _pft(code):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(code)
            return PYFT(fpath)

    def test_continuation_in_section_detected(self, caplog):
        """A continuation line inside an array section is reported."""
        pft = self._pft("""
SUBROUTINE FOO(A, B)
  REAL, DIMENSION(:,:,:) :: A, B
  A(:,:,:)=B(:, &
  &:,:)
END SUBROUTINE FOO
""")
        assert pft.checkContinuationInArraySection() is False
        assert len(caplog.records) == 1
        assert 'B(:,' in caplog.records[0].getMessage()

    def test_continuation_in_section_raise(self):
        """mustRaise=True raises a PYFTError."""
        from pyfortool.util import PYFTError
        pft = self._pft("""
SUBROUTINE FOO(A, B)
  REAL, DIMENSION(:,:) :: A, B
  A(:,:)=B(1:&
  &2, :)
END SUBROUTINE FOO
""")
        with pytest.raises(PYFTError):
            pft.checkContinuationInArraySection(mustRaise=True)

    def test_nested_section_reported_once(self, caplog):
        """Nested array sections are reported only once (outermost)."""
        pft = self._pft("""
SUBROUTINE FOO(A, B, C)
  REAL, DIMENSION(:,:) :: A, B, C
  A(:,:)=B(:, C(:, &
  &1))
END SUBROUTINE FOO
""")
        assert pft.checkContinuationInArraySection() is False
        assert len(caplog.records) == 1

    def test_continuation_in_element_detected(self, caplog):
        """A continuation line inside an array element selection is reported."""
        pft = self._pft("""
SUBROUTINE FOO(A)
  REAL, DIMENSION(:,:,:) :: A
  A(1, &
   &2, 3)=1.
END SUBROUTINE FOO
""")
        assert pft.checkContinuationInArraySection() is False
        assert len(caplog.records) == 1
        assert 'A(1,' in caplog.records[0].getMessage()

    def test_continuation_in_function_ignored(self, caplog):
        """A function call (name not declared as an array) is not reported."""
        pft = self._pft("""
SUBROUTINE FOO(A)
  REAL, DIMENSION(:,:,:) :: A
  A(1,1,1)=F(1, &
            &2)
END SUBROUTINE FOO
""")
        assert pft.checkContinuationInArraySection(mustRaise=True) is True
        assert len(caplog.records) == 0

    def test_array_declared_in_host_scope(self, caplog):
        """Arrays declared in a parent scope are found, function names are not."""
        pft = self._pft("""
MODULE MODTEST
  REAL, DIMENSION(10, 10) :: XA
CONTAINS
  SUBROUTINE SUB1(P)
    REAL :: P
    P=XA(1, &
        &2)
  END SUBROUTINE SUB1
  SUBROUTINE SUB2(P)
    REAL :: P
    REAL, DIMENSION(10, 10) :: G
    P=G(1, &
       &2)
  END SUBROUTINE SUB2
  SUBROUTINE SUB3(P)
    REAL :: P
    P=G(1, &
       &2)
  END SUBROUTINE SUB3
END MODULE MODTEST
""")
        assert pft.checkContinuationInArraySection() is False
        messages = [rec.getMessage() for rec in caplog.records]
        assert len(messages) == 2
        assert any('XA(1,' in msg for msg in messages)
        assert any('G(1,' in msg for msg in messages)

    def test_structure_component(self, caplog):
        """Component selections are reported only when they contain a ':'."""
        pft = self._pft("""
SUBROUTINE FOO(X, P, Q)
  USE MODD_TYPE, ONLY: TYP
  TYPE(TYP) :: X
  REAL :: P
  REAL, DIMENSION(:) :: Q
  P=X%Y(1, &
       &2)
  Q(:)=X%Z(:, &
          &2)
END SUBROUTINE FOO
""")
        assert pft.checkContinuationInArraySection() is False
        assert len(caplog.records) == 1
        assert 'X%Z(:,' in caplog.records[0].getMessage()

    def test_valid_code(self, caplog):
        """Continuations between operands or in CALL arguments are accepted."""
        pft = self._pft("""
SUBROUTINE FOO(A, B)
  REAL, DIMENSION(:,:,:) :: A, B
  A(:,:,:)=B(:,:,:) + &
          &B(:,:,:)
  CALL F(A, &
        &B)
END SUBROUTINE FOO
""")
        assert pft.checkContinuationInArraySection(mustRaise=True) is True
        assert len(caplog.records) == 0
