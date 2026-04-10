"""
Tests for PYFTscope class in pyfortool.scope.
"""

import pytest
from pyfortool import PYFT


class TestGetScopes:
    """Tests for getScopes() method."""

    def test_get_all_scopes(self, pft_module):
        """Test getting all scopes."""
        scopes = pft_module.getScopes()
        assert len(scopes) > 0

    def test_get_scopes_with_level_1(self, pft_module):
        """Test getting scopes with level=1 (direct children only)."""
        scopes = pft_module.getScopes(level=1)
        assert len(scopes) >= 1

    def test_get_scopes_with_level_minus_1(self, pft_module):
        """Test getting scopes with level=-1 (all descendants)."""
        scopes = pft_module.getScopes(level=-1)
        assert len(scopes) >= 1

    def test_exclude_contains(self, pft_module):
        """Test excludeContains parameter."""
        scopes_with_contains = pft_module.getScopes(excludeContains=False)
        scopes_without_contains = pft_module.getScopes(excludeContains=True)
        assert len(scopes_with_contains) >= len(scopes_without_contains)

    def test_include_itself(self, pft_module):
        """Test includeItself parameter."""
        scopes_with_self = pft_module.getScopes(includeItself=True)
        scopes_without_self = pft_module.getScopes(includeItself=False)
        assert len(scopes_with_self) >= len(scopes_without_self)

    def test_exclude_kinds(self, pft_module):
        """Test excludeKinds parameter."""
        scopes_all = pft_module.getScopes()
        scopes_no_sub = pft_module.getScopes(excludeKinds=['sub'])
        assert len(scopes_all) > len(scopes_no_sub)


class TestGetScopeNode:
    """Tests for getScopeNode() method."""

    def test_get_module_scope(self, pft_module):
        """Test getting module scope."""
        scope = pft_module.getScopeNode('module:MOD_TEST')
        assert scope is not None
        assert 'MOD_TEST' in scope.path

    def test_get_subroutine_scope(self, pft_module):
        """Test getting subroutine scope."""
        scope = pft_module.getScopeNode('module:MOD_TEST/sub:SUB')
        assert scope is not None
        assert 'SUB' in scope.path

    def test_get_nonexistent_scope(self, pft_module):
        """Test getting nonexistent scope raises error."""
        from pyfortool.util import PYFTError
        with pytest.raises(PYFTError):
            pft_module.getScopeNode('module:NONEXISTENT')


class TestScopeProperties:
    """Tests for scope properties."""

    def test_path_property(self, pft_module):
        """Test path property."""
        scope = pft_module.getScopeNode('module:MOD_TEST')
        assert 'module:MOD_TEST' in scope.path

    def test_main_scope_property(self, pft_module):
        """Test mainScope property."""
        scope = pft_module.getScopeNode('module:MOD_TEST/sub:SUB')
        assert scope.mainScope is not None

    def test_parent_scope_property(self, pft_module):
        """Test parentScope property."""
        scope = pft_module.getScopeNode('module:MOD_TEST/sub:SUB')
        parent = scope.parentScope
        assert parent is not None
        assert hasattr(parent, 'path')


class TestGetParent:
    """Tests for getParent() method."""

    def test_get_direct_parent(self, pft_module):
        """Test getting direct parent."""
        scope = pft_module.getScopeNode('module:MOD_TEST/sub:SUB')
        stmts = scope.findall('.//{*}a-stmt')
        if stmts:
            parent = scope.getParent(stmts[0])
            assert parent is not None

    def test_get_grandparent(self, pft_module):
        """Test getting grandparent with level=2."""
        scope = pft_module.getScopeNode('module:MOD_TEST/sub:SUB')
        stmts = scope.findall('.//{*}a-stmt')
        if stmts:
            grandparent = scope.getParent(stmts[0], level=2)
            assert grandparent is not None


class TestGetSiblings:
    """Tests for getSiblings() method."""

    def test_get_all_siblings(self, pft_module):
        """Test getting all siblings."""
        scope = pft_module.getScopeNode('module:MOD_TEST/sub:SUB')
        stmts = scope.findall('.//{*}a-stmt')
        if len(stmts) >= 2:
            siblings = scope.getSiblings(stmts[1])
            assert isinstance(siblings, list)

    def test_get_siblings_before(self, pft_module):
        """Test getting siblings before only."""
        scope = pft_module.getScopeNode('module:MOD_TEST/sub:SUB')
        stmts = scope.findall('.//{*}a-stmt')
        if len(stmts) >= 2:
            siblings = scope.getSiblings(stmts[1], after=False)
            assert isinstance(siblings, list)


class TestNormalizeScope:
    """Tests for normalizeScope() static method."""

    def test_normalize_mixed_case(self, pft_module):
        """Test normalizing mixed case scope path."""
        from pyfortool.scope import PYFTscope
        result = PYFTscope.normalizeScope('module:Test/sub:Sub')
        assert 'TEST' in result
        assert 'SUB' in result

    def test_normalize_already_normal(self, pft_module):
        """Test normalizing already normal scope path."""
        from pyfortool.scope import PYFTscope
        result = PYFTscope.normalizeScope('module:TEST/sub:SUB')
        assert result == 'module:TEST/sub:SUB'


class TestGetScopePath:
    """Tests for getScopePath() method."""

    def test_get_scope_path_for_element(self, pft_module):
        """Test getting scope path for an element."""
        scope = pft_module.getScopeNode('module:MOD_TEST/sub:SUB')
        stmts = scope.findall('.//{*}a-stmt')
        if stmts:
            path = scope.getScopePath(stmts[0])
            assert 'module:' in path or 'sub:' in path


class TestGetParentScopeNode:
    """Tests for getParentScopeNode() method."""

    def test_get_parent_scope_node(self, pft_module):
        """Test getting parent scope node."""
        scope = pft_module.getScopeNode('module:MOD_TEST/sub:SUB')
        stmts = scope.findall('.//{*}a-stmt')
        if stmts:
            parent_scope = scope.getParentScopeNode(stmts[0])
            assert parent_scope is not None


class TestIsScopeNode:
    """Tests for isScopeNode() method."""

    def test_subroutine_is_scope(self, pft_module):
        """Test subroutine element is recognized as scope."""
        scope = pft_module.getScopeNode('module:MOD_TEST/sub:SUB')
        assert scope.isScopeNode(scope[0]) or scope.isScopeNode(scope)

    def test_statement_not_scope(self, pft_module):
        """Test statement element is not recognized as scope."""
        scope = pft_module.getScopeNode('module:MOD_TEST/sub:SUB')
        stmts = scope.findall('.//{*}a-stmt')
        if stmts:
            assert not scope.isScopeNode(stmts[0])


class TestShowScopesList:
    """Tests for showScopesList() method."""

    def test_show_scopes_list_runs(self, pft_module, capsys):
        """Test showScopesList() runs without error."""
        pft_module.showScopesList()
        captured = capsys.readouterr()
        assert 'MOD_TEST' in captured.out or 'scope' in captured.out.lower()
