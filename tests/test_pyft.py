"""
Tests for the PYFT class in pyfortool.pyfortool.
"""

import pytest
import os
import tempfile
from pyfortool import PYFT


class TestPYFTInit:
    """Tests for PYFT initialization."""

    def test_simple_init(self, simple_fortran):
        """Test basic PYFT initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(simple_fortran)
            pft = PYFT(fpath)
            assert pft is not None
            pft.close()

    def test_init_with_output(self, simple_fortran):
        """Test PYFT initialization with output file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, 'input.F90')
            out_path = os.path.join(tmpdir, 'output.F90')
            with open(in_path, 'w') as f:
                f.write(simple_fortran)
            pft = PYFT(in_path, output=out_path)
            assert pft is not None
            pft.close()

    def test_context_manager(self, simple_fortran):
        """Test PYFT as context manager."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(simple_fortran)
            with PYFT(fpath) as pft:
                assert pft is not None


class TestPYFTProperties:
    """Tests for PYFT properties."""

    def test_xml_property(self, pft_simple):
        """Test xml property returns string."""
        xml = pft_simple.xml
        assert isinstance(xml, str)
        assert '<f:subroutine-stmt>' in xml or '<f:program-stmt>' in xml

    def test_fortran_property(self, pft_simple):
        """Test fortran property returns string."""
        fortran = pft_simple.fortran
        assert isinstance(fortran, str)
        assert 'PROGRAM' in fortran or 'SUBROUTINE' in fortran


class TestPYFTFileOperations:
    """Tests for PYFT file operations."""

    def test_write(self, pft_simple):
        """Test write() method."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, 'output.F90')
            pft_simple._output = out_path
            pft_simple.write()
            assert os.path.exists(out_path)
            with open(out_path) as f:
                content = f.read()
            assert 'PROGRAM' in content or 'SUBROUTINE' in content

    def test_writexml(self, pft_simple):
        """Test writeXML() method."""
        with tempfile.TemporaryDirectory() as tmpdir:
            xml_path = os.path.join(tmpdir, 'output.xml')
            pft_simple.writeXML(xml_path)
            assert os.path.exists(xml_path)
            with open(xml_path) as f:
                content = f.read()
            assert 'subroutine-stmt' in content or 'program-stmt' in content

    def test_rename_upper(self, simple_fortran):
        """Test renameUpper() method."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.f90')
            with open(fpath, 'w') as f:
                f.write(simple_fortran)
            pft = PYFT(fpath)
            pft.renameUpper()
            assert pft._filename.endswith('.F90')
            pft.close()

    def test_rename_lower(self, simple_fortran):
        """Test renameLower() method."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = os.path.join(tmpdir, 'test.F90')
            with open(fpath, 'w') as f:
                f.write(simple_fortran)
            pft = PYFT(fpath)
            pft.renameLower()
            assert pft._filename.endswith('.f90')
            pft.close()


class TestPYFTGetFileName:
    """Tests for getFileName() method."""

    def test_get_filename(self, pft_simple):
        """Test getFileName() returns path."""
        fname = pft_simple.getFileName()
        assert isinstance(fname, str)
        assert 'test' in fname

    def test_filename_normalized(self, pft_simple):
        """Test getFileName() returns normalized path."""
        fname = pft_simple.getFileName()
        assert os.path.isabs(fname) or fname.startswith('/')


class TestPYFTInheritance:
    """Tests for PYFT inheritance from PYFTscope."""

    def test_has_varList(self, pft_module):
        """Test PYFT has varList property."""
        assert hasattr(pft_module, 'varList')

    def test_has_getScopes(self, pft_module):
        """Test PYFT has getScopes() method."""
        assert hasattr(pft_module, 'getScopes')
        scopes = pft_module.getScopes()
        assert isinstance(scopes, list)

    def test_has_uppercase(self, pft_simple):
        """Test PYFT has upperCase() method."""
        assert hasattr(pft_simple, 'upperCase')

    def test_has_lower_case(self, pft_simple):
        """Test PYFT has lowerCase() method."""
        assert hasattr(pft_simple, 'lowerCase')


class TestPYFTMultipleScopes:
    """Tests for PYFT with multiple scopes."""

    def test_get_multiple_scopes(self, pft_module):
        """Test getting multiple scopes."""
        scopes = pft_module.getScopes()
        assert len(scopes) >= 2

    def test_scope_paths(self, pft_module):
        """Test scope paths are correct."""
        scopes = pft_module.getScopes()
        paths = [s.path for s in scopes]
        assert any('module:' in p for p in paths)
        assert any('sub:' in p for p in paths)

    def test_get_scope_node(self, pft_module):
        """Test getting specific scope by path."""
        scope = pft_module.getScopeNode('module:MOD_TEST')
        assert scope is not None
        assert 'MOD_TEST' in scope.path

    def test_get_nested_scope(self, pft_module):
        """Test getting nested scope."""
        scope = pft_module.getScopeNode('module:MOD_TEST/sub:SUB')
        assert scope is not None
        assert 'SUB' in scope.path
