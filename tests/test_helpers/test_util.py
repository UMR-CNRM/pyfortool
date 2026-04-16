"""
Tests for utility functions in pyfortool.util.
"""

import pytest
import xml.etree.ElementTree as ET
from pyfortool.util import (
    tag, n2name, alltext, isint, isfloat,
    nonCode, isExecutable, isStmt, isConstruct
)


class TestTag:
    """Tests for tag() function."""

    def test_simple_tag(self):
        """Test extracting tag without namespace."""
        elem = ET.Element('{http://fxtran.net/#syntax}subroutine-stmt')
        assert tag(elem) == 'subroutine-stmt'

    def test_tag_with_namespace(self):
        """Test extracting tag with different namespace."""
        elem = ET.Element('{custom}my-tag')
        assert tag(elem) == 'my-tag'

    def test_tag_without_namespace(self):
        """Test tag with namespace prefix."""
        elem = ET.Element('{http://fxtran.net/#syntax}simple')
        assert tag(elem) == 'simple'


class TestN2Name:
    """Tests for n2name() function."""

    def test_single_n(self):
        """Test extracting name from single n element."""
        from pyfortool.expressions import createElem
        nodeN = createElem('N')
        nodeN.append(createElem('n', text='X'))
        assert n2name(nodeN) == 'X'

    def test_multiple_n(self):
        """Test extracting name from multiple n elements."""
        from pyfortool.expressions import createElem
        nodeN = createElem('N')
        nodeN.append(createElem('n', text='AR'))
        nodeN.append(createElem('n', text='RAY'))
        assert n2name(nodeN) == 'ARRAY'

    def test_empty_n(self):
        """Test with empty n elements."""
        from pyfortool.expressions import createElem
        nodeN = createElem('N')
        nodeN.append(createElem('n', text=''))
        assert n2name(nodeN) == ''


class TestAllText:
    """Tests for alltext() function."""

    def test_simple_text(self):
        """Test extracting all text from simple element."""
        elem = ET.Element('test')
        elem.text = 'hello'
        assert alltext(elem) == 'hello'

    def test_with_tail(self):
        """Test extracting text including tail."""
        parent = ET.Element('parent')
        child = ET.SubElement(parent, 'child')
        child.text = 'hello'
        child.tail = ' world'
        assert alltext(parent) == 'hello world'

    def test_nested_elements(self):
        """Test extracting text from nested elements."""
        parent = ET.Element('parent')
        child1 = ET.SubElement(parent, 'child1')
        child1.text = 'a'
        child2 = ET.SubElement(parent, 'child2')
        child2.text = 'b'
        assert alltext(parent) == 'ab'

    def test_deeply_nested(self):
        """Test extracting text from deeply nested elements."""
        root = ET.Element('root')
        level1 = ET.SubElement(root, 'l1')
        level1.text = 'text'
        level2 = ET.SubElement(level1, 'l2')
        level2.text = 'more'
        assert 'textmore' in alltext(root)


class TestIsInt:
    """Tests for isint() function."""

    def test_positive_int(self):
        """Test positive integers."""
        assert isint('42') is True
        assert isint('0') is True
        assert isint('123456') is True

    def test_negative_int(self):
        """Test negative integers."""
        assert isint('-42') is True
        assert isint('-123') is True

    def test_float_not_int(self):
        """Test that floats are not integers."""
        assert isint('3.14') is False
        assert isint('0.0') is False
        assert isint('-1.5') is False

    def test_string_not_int(self):
        """Test that non-numeric strings are not integers."""
        assert isint('hello') is False
        assert isint('42abc') is False
        assert isint('') is False

    def test_scientific_notation(self):
        """Test scientific notation is not integer."""
        assert isint('1e10') is False


class TestIsFloat:
    """Tests for isfloat() function."""

    def test_positive_float(self):
        """Test positive floats."""
        assert isfloat('3.14') is True
        assert isfloat('0.0') is True
        assert isfloat('123.456') is True

    def test_negative_float(self):
        """Test negative floats."""
        assert isfloat('-3.14') is True
        assert isfloat('-0.001') is True

    def test_int_is_float(self):
        """Test that integers are also valid floats."""
        assert isfloat('42') is True
        assert isfloat('0') is True

    def test_scientific_notation_float(self):
        """Test scientific notation is valid float."""
        assert isfloat('1e10') is True
        assert isfloat('1.5e-3') is True

    def test_string_not_float(self):
        """Test that non-numeric strings are not floats."""
        assert isfloat('hello') is False
        assert isfloat('3.14.15') is False
        assert isfloat('') is False


class TestNonCode:
    """Tests for nonCode() function."""

    def test_comment(self):
        """Test that comment elements are non-code."""
        from pyfortool.expressions import createElem
        elem = createElem('C')
        assert nonCode(elem) is True

    def test_cnt(self):
        """Test that cnt elements are non-code."""
        from pyfortool.expressions import createElem
        elem = createElem('cnt')
        assert nonCode(elem) is True

    def test_cpp(self):
        """Test that cpp elements are non-code."""
        from pyfortool.expressions import createElem
        elem = createElem('cpp')
        assert nonCode(elem) is True

    def test_code_element(self):
        """Test that code elements are not non-code."""
        from pyfortool.expressions import createElem
        elem = createElem('a-stmt')
        assert nonCode(elem) is False

    def test_subroutine_stmt(self):
        """Test subroutine-stmt is not non-code."""
        from pyfortool.expressions import createElem
        elem = createElem('subroutine-stmt')
        assert nonCode(elem) is False


class TestIsStmt:
    """Tests for isStmt() function."""

    def test_stmt_with_suffix(self):
        """Test elements ending with -stmt."""
        from pyfortool.expressions import createElem
        elem = createElem('a-stmt')
        assert isStmt(elem) is True

    def test_call_stmt(self):
        """Test call-stmt."""
        from pyfortool.expressions import createElem
        elem = createElem('call-stmt')
        assert isStmt(elem) is True

    def test_construct(self):
        """Test that constructs are not statements."""
        from pyfortool.expressions import createElem
        elem = createElem('if-construct')
        assert isStmt(elem) is False

    def test_no_suffix(self):
        """Test elements without -stmt suffix."""
        from pyfortool.expressions import createElem
        elem = createElem('named-E')
        assert isStmt(elem) is False


class TestIsConstruct:
    """Tests for isConstruct() function."""

    def test_construct_with_suffix(self):
        """Test elements ending with -construct."""
        from pyfortool.expressions import createElem
        elem = createElem('if-construct')
        assert isConstruct(elem) is True

    def test_do_construct(self):
        """Test do-construct."""
        from pyfortool.expressions import createElem
        elem = createElem('do-construct')
        assert isConstruct(elem) is True

    def test_stmt_not_construct(self):
        """Test that statements are not constructs."""
        from pyfortool.expressions import createElem
        elem = createElem('a-stmt')
        assert isConstruct(elem) is False

    def test_no_suffix(self):
        """Test elements without -construct suffix."""
        from pyfortool.expressions import createElem
        elem = createElem('named-E')
        assert isConstruct(elem) is False


class TestIsExecutable:
    """Tests for isExecutable() function."""

    def test_call_stmt(self):
        """Test call statement is executable."""
        from pyfortool.expressions import createElem
        elem = createElem('call-stmt')
        assert isExecutable(elem) is True

    def test_a_stmt(self):
        """Test assignment statement is executable."""
        from pyfortool.expressions import createElem
        elem = createElem('a-stmt')
        assert isExecutable(elem) is True

    def test_if_construct(self):
        """Test if-construct is executable."""
        from pyfortool.expressions import createElem
        elem = createElem('if-construct')
        assert isExecutable(elem) is True

    def test_decl_stmt_not_executable(self):
        """Test declaration statement is not executable."""
        from pyfortool.expressions import createElem
        elem = createElem('T-decl-stmt')
        assert isExecutable(elem) is False

    def test_subroutine_stmt_not_executable(self):
        """Test subroutine statement is not executable."""
        from pyfortool.expressions import createElem
        elem = createElem('subroutine-stmt')
        assert isExecutable(elem) is False

    def test_use_stmt_not_executable(self):
        """Test use statement is not executable."""
        from pyfortool.expressions import createElem
        elem = createElem('use-stmt')
        assert isExecutable(elem) is False
