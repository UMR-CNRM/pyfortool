"""
Tests for the pyfortool_parallel command (scripting.mainParallel).

The examples/tests.sh regression suite covers the pyfortool command; these
tests cover the parallel command, which applies transformations on all the
files of a tree, in several processes sharing the same Tree instance.
"""

import json
import pickle

import pytest

from pyfortool import PYFT
from pyfortool.scripting import mainParallel, TreeManager, poolInit
from pyfortool.util import PYFTError


MOD_A = """MODULE MOD_A
IMPLICIT NONE
CONTAINS
SUBROUTINE SUB_A(N, X)
INTEGER, INTENT(IN) :: N
REAL, INTENT(INOUT) :: X(N)
CALL SUB_B(N, X)
END SUBROUTINE SUB_A
END MODULE MOD_A
"""

# No IMPLICIT NONE in SUB_B and SUB_C: used by the error test
SUB_B = """SUBROUTINE SUB_B(N, Y)
INTEGER, INTENT(IN) :: N
REAL, INTENT(INOUT) :: Y(N)
REAL :: ZWORK(N)
ZWORK(:) = 1.
Y(:) = Y(:) + ZWORK(:)
END SUBROUTINE SUB_B
"""

SUB_C = """SUBROUTINE SUB_C(N, Y)
INTEGER, INTENT(IN) :: N
REAL, INTENT(INOUT) :: Y(N)
REAL :: ZWORK(N)
ZWORK(:) = 1.
Y(:) = ZWORK(:)
END SUBROUTINE SUB_C
"""


@pytest.fixture
def src_tree(tmp_path, monkeypatch):
    """
    Source tree spread over sub-directories; cwd is set to its parent.

    Call graph: SUB_A (src/mod_a.F90) calls SUB_B (src/sub/sub_b.F90).
    SUB_C (src/sub/sub_c.F90) is not called.
    """
    (tmp_path / 'src' / 'sub').mkdir(parents=True)
    files = {'a': tmp_path / 'src' / 'mod_a.F90',
             'b': tmp_path / 'src' / 'sub' / 'sub_b.F90',
             'c': tmp_path / 'src' / 'sub' / 'sub_c.F90'}
    files['a'].write_text(MOD_A)
    files['b'].write_text(SUB_B)
    files['c'].write_text(SUB_C)
    monkeypatch.chdir(tmp_path)
    return files


def run_parallel(*options):
    """Run pyfortool_parallel on the 'src' tree with the given options."""
    mainParallel(['pyfortool_parallel', '--tree', 'src', '--descTree', 'desc.json',
                  '--nbPar', '2', '--logLevel', 'error'] + list(options))


class TestMainParallel:
    """Tests for the pyfortool_parallel command."""

    def test_transforms_all_files(self, src_tree):
        """Test a transformation is applied, in place, on every file of the tree."""
        for path in src_tree.values():
            path.write_text(path.read_text().lower())

        run_parallel('--upperCase')

        assert src_tree['a'].read_text() == MOD_A
        assert src_tree['b'].read_text() == SUB_B
        assert src_tree['c'].read_text() == SUB_C

    def test_writes_desc_tree(self, src_tree):
        """Test the description of the tree is written in the --descTree file."""
        run_parallel('--upperCase')

        with open('desc.json', encoding='utf-8') as file:
            desc = json.load(file)
        assert desc['scopes'] == {'src/mod_a.F90': ['module:MOD_A', 'module:MOD_A/sub:SUB_A'],
                                  'src/sub/sub_b.F90': ['sub:SUB_B'],
                                  'src/sub/sub_c.F90': ['sub:SUB_C']}

    def test_dry_run(self, src_tree):
        """Test --dryRun leaves the files untouched."""
        for path in src_tree.values():
            path.write_text(path.read_text().lower())
        before = {key: path.read_text() for key, path in src_tree.items()}

        run_parallel('--upperCase', '--dryRun')

        assert {key: path.read_text() for key, path in src_tree.items()} == before

    def test_cross_file_transformation(self, src_tree):
        """
        Test a transformation modifying other files through the shared tree.

        addStack adds the YDSTACK argument to SUB_B, and also to its caller
        SUB_A, which is in another file (processed by another process).
        """
        run_parallel('--addStack', 'AROME', '--stopScopes', 'module:MOD_A/sub:SUB_A')

        resultA = src_tree['a'].read_text()
        resultB = src_tree['b'].read_text()
        assert 'SUBROUTINE SUB_B(N, Y, YDSTACK)' in resultB
        assert 'POINTER(IP_ZWORK_, ZWORK)' in resultB
        # Caller, in another file, updated
        assert 'SUBROUTINE SUB_A(N, X, YDSTACK)' in resultA
        assert 'CALL SUB_B(N, X, YDSTACK=YDSTACK)' in resultA
        # Not called from the stop scope
        assert src_tree['c'].read_text() == SUB_C

    def test_opts_by_env(self, src_tree, monkeypatch):
        """Test --optsByEnv applies file-specific options."""
        # The regular expression is matched (re.match) against the file path
        monkeypatch.setenv('PYFT_OPTS', r'.*/sub_b\.F90:=:--lowerCase')

        run_parallel('--optsByEnv', 'PYFT_OPTS')

        assert src_tree['a'].read_text() == MOD_A
        assert src_tree['b'].read_text() == SUB_B.lower()
        assert src_tree['c'].read_text() == SUB_C

    def test_error_reported(self, src_tree):
        """Test an error in some files is reported with the number of files in error."""
        with pytest.raises(PYFTError, match='Errors have been reported in 2 file'):
            run_parallel('--checkIMPLICIT', 'Err')

    def test_parallel_unset_after_error(self, src_tree):
        """Test the parallel configuration of PYFT is removed even after an error."""
        with pytest.raises(PYFTError):
            run_parallel('--checkIMPLICIT', 'Err')

        assert PYFT.SHARED_TREE is None
        assert PYFT.PARALLEL_FILE_LOCKS is None
        assert PYFT.NO_PARALLEL_LOCK is None
        # A PYFT instance can still be created (it used the dead shared tree before)
        assert 'SUBROUTINE SUB_C' in PYFT(str(src_tree['c'])).fortran


class TestParallelSetup:
    """
    Tests for the objects used to set up parallel processing.

    The 'forkserver' and 'spawn' start methods send the process target through a
    pipe, so these objects must be picklable. Pickle serialises classes and
    functions by qualified name, which means they cannot be defined inside
    mainParallel(). Only the 'fork' start method, which inherits the parent
    memory, tolerates that.
    """

    def test_manager_is_picklable(self):
        """TreeManager must be reachable by qualified name."""
        assert pickle.loads(pickle.dumps(TreeManager)) is TreeManager

    def test_pool_initializer_is_picklable(self):
        """poolInit must be reachable by qualified name."""
        assert pickle.loads(pickle.dumps(poolInit)) is poolInit

    def test_tree_is_registered_on_manager(self):
        """The Tree type must be registered on TreeManager at import time."""
        assert 'Tree' in TreeManager._registry  # pylint: disable=protected-access
