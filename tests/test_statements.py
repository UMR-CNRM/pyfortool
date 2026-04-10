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
