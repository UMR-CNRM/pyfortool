"""
Tests for the pyfortool.scripting module.
"""

import pickle

from pyfortool.scripting import MyManager, poolInit


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
        """MyManager must be reachable by qualified name."""
        assert pickle.loads(pickle.dumps(MyManager)) is MyManager

    def test_pool_initializer_is_picklable(self):
        """poolInit must be reachable by qualified name."""
        assert pickle.loads(pickle.dumps(poolInit)) is poolInit

    def test_tree_is_registered_on_manager(self):
        """The Tree type must be registered on MyManager at import time."""
        assert 'Tree' in MyManager._registry  # pylint: disable=protected-access
