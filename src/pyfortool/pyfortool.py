#!/usr/bin/env python3

"""
This module contains the main PYFT class.

PYFT (Python FORTRAN Tool) is the primary class for reading, manipulating,
and writing FORTRAN source files. It provides file-level operations and
wraps the scope-level functionality from PYFTscope.

Examples
--------
>>> from pyfortool import PYFT

# Read and modify a file
>>> pft = PYFT('input.F90')
>>> pft.removePrints()
>>> pft.write()

# Use context manager
>>> with PYFT('input.F90') as pft:
...     pft.removeComments()
...     pft.write()
"""

import os
import sys
from multiprocessing import Lock, RLock

from pyfortool.scope import PYFTscope
from pyfortool.tree import Tree, updateTree
from pyfortool.util import (debugDecor, tostring, tofortran, fortran2xml,
                            setVerbosity, printInfos, PYFTError)


@debugDecor
def conservativePYFT(filename, parserOptions, wrapH,
                     tree=None, verbosity=None, clsPYFT=None):
    """
    Create a PYFT object with conservative parsing options.

    Similar to PYFT constructor but forces the '-no-include' option
    to prevent automatic inclusion of header files.

    Parameters
    ----------
    filename : str
        Path to the FORTRAN source file.
    parserOptions : list or None
        Parser options passed to PYFT. If None, uses default options.
    wrapH : bool
        Whether to wrap .h file content. See PYFT class.
    tree : Tree, optional
        Tree instance for cross-file analysis.
    verbosity : str or int, optional
        Logging verbosity level.
    clsPYFT : type, optional
        PYFT subclass to use instead of default PYFT.

    Returns
    -------
    PYFT
        A PYFT instance configured for conservative tree manipulation.

    See Also
    --------
    PYFT : The main PYFT class.
    """
    options = PYFT.DEFAULT_FXTRAN_OPTIONS if parserOptions is None else parserOptions
    options = options.copy()
    if len(set(options).intersection(('-no-include', '-noinclude'))) == 0:
        # We add the option to not include 'include files' when analysing the tree
        options.append('-no-include')
    if clsPYFT is None:
        clsPYFT = PYFT
    pft = clsPYFT(filename, parserOptions=options, wrapH=wrapH,
                  tree=tree, verbosity=verbosity)
    return pft


def generateEmptyPYFT(filename, fortran=None, **kwargs):
    """
    Generate a new PYFT object from a file path, creating the file if needed.

    Parameters
    ----------
    filename : str
        Path for the new file.
    fortran : str, optional
        FORTRAN source code to write to the file.
        If None, creates a stub subroutine "SUBROUTINE FOO\nEND".
    **kwargs
        Additional arguments passed to PYFT constructor.

    Returns
    -------
    PYFT
        PYFT instance for the newly created file.

    Examples
    --------
    >>> pft = generateEmptyPYFT('new.F90', 'MODULE TEST\nEND MODULE')
    >>> pft.addVar([('module:TEST', 'X', 'INTEGER :: X', None)])
    >>> pft.write()
    """
    with open(filename, 'w', encoding='utf-8') as fo:
        fo.write('SUBROUTINE FOO\nEND' if fortran is None else fortran)
    pft = PYFT(filename, **kwargs)
    if fortran is None:
        pft.remove(pft.find('.//{*}program-unit'))
    return pft


class PYFT(PYFTscope):
    """
    Main class for FORTRAN file manipulation.

    PYFT extends PYFTscope with file-level operations (read/write) and
    integrates with the fxtran parser for FORTRAN source code analysis.

    Class Attributes
    ----------------
    DEFAULT_FXTRAN_OPTIONS : list
        Default options for fxtran parser: ['-construct-tag', '-no-include',
        '-no-cpp', '-line-length', '9999']
    MANDATORY_FXTRAN_OPTIONS : list
        Required options: ['-construct-tag']

    Examples
    --------
    Basic usage:
    >>> pft = PYFT('myfile.F90')
    >>> pft.upperCase()
    >>> pft.write()

    Using context manager:
    >>> with PYFT('myfile.F90', output='output.F90') as pft:
    ...     pft.removeComments()
    ...     pft.indent()

    Parallel processing:
    >>> tree = Tree(['/path/to/src'])
    >>> PYFT.setParallel(tree)
    """
    DEFAULT_FXTRAN_OPTIONS = ['-construct-tag', '-no-include', '-no-cpp', '-line-length', '9999']
    MANDATORY_FXTRAN_OPTIONS = ['-construct-tag']
    SHARED_TREE = None  # Can be updated by setParallel
    NO_PARALLEL_LOCK = None  # Can be updated by setParallel
    PARALLEL_FILE_LOCKS = None  # Can be updated by setParallel

    @updateTree('signal')
    def __init__(self, filename, output=None, parserOptions=None, verbosity=None,
                 wrapH=False, tree=None, enableCache=False):
        """
        Initialize a PYFT instance from a FORTRAN source file.

        Parameters
        ----------
        filename : str
            Path to the input FORTRAN file.
        output : str, optional
            Path for output file. If None, overwrites the input file.
        parserOptions : list, optional
            Options for the fxtran parser. See DEFAULT_FXTRAN_OPTIONS.
        verbosity : str or int, optional
            Logging level (e.g., 'DEBUG', 'INFO', 'WARNING').
        wrapH : bool, optional
            If True, wrap .h file content in a MODULE to enable parsing
            as free-form FORTRAN.
        tree : Tree, optional
            Tree instance for cross-file analysis. If None, creates a new Tree.
        enableCache : bool, optional
            If True, cache parent nodes for faster traversal.

        Raises
        ------
        PYFTError
            If Python version < 3.8 or file does not exist.

        Examples
        --------
        >>> pft = PYFT('myfile.F90')
        >>> pft = PYFT('myfile.F90', output='newfile.F90')
        >>> pft = PYFT('code.h', wrapH=True)
        """
        self.__class__.lockFile(filename)
        if not sys.version_info >= (3, 8):
            # At least version 3.7 for ordered dictionary
            # At least verison 3.8 for namsepace wildcard (use of '{*}' in find or findall)
            raise PYFTError("PyForTool needs at least version 3.8 of python")
        self._filename = filename
        self._originalName = filename
        assert os.path.exists(filename), f'Input filename ({filename})must exist'
        self._output = output
        tree = Tree() if tree is None else tree
        if self.SHARED_TREE is not None:
            assert tree is not None, 'tree must be None if setParallel has been called'
            tree.copyFromOtherTree(self.SHARED_TREE)
        if parserOptions is None:
            self._parserOptions = self.DEFAULT_FXTRAN_OPTIONS.copy()
        else:
            self._parserOptions = parserOptions.copy()
        for tDir in tree.getDirs():
            self._parserOptions.extend(['-I', tDir])
        for option in self.MANDATORY_FXTRAN_OPTIONS:
            if option not in self._parserOptions:
                self._parserOptions.append(option)
        includesRemoved, xml = fortran2xml(self._filename, self._parserOptions, wrapH)
        super().__init__(xml, enableCache=enableCache, tree=tree)
        if includesRemoved:
            self.tree.signal(self)
        if verbosity is not None:
            setVerbosity(verbosity)

    @classmethod
    def setParallel(cls, tree, clsLock=None, clsRLock=None):
        """
        Configure PYFT for parallel processing.

        Must be called before creating PYFT instances for parallel execution.
        Sets up shared tree and file locking mechanisms.

        Parameters
        ----------
        tree : Tree
            Tree object shared among all processes.
        clsLock : type, optional
            Lock class to use. Defaults to multiprocessing.Lock.
        clsRLock : type, optional
            Recursive lock class. Defaults to multiprocessing.RLock.

        Examples
        --------
        >>> tree = Tree(['/path/to/src'], descTreeFile='tree.json')
        >>> PYFT.setParallel(tree)
        >>> # Now create PYFT instances in parallel processes
        """
        if clsLock is None:
            clsLock = Lock
        if clsRLock is None:
            clsRLock = RLock
        cls.NO_PARALLEL_LOCK = clsRLock()
        cls.SHARED_TREE = tree
        cls.PARALLEL_FILE_LOCKS = {os.path.normpath(file): clsLock()
                                   for file in tree.getFiles()}

    @classmethod
    def unsetParallel(cls):
        """
        Remove parallel processing configuration
        """
        cls.NO_PARALLEL_LOCK = None
        cls.SHARED_TREE = None
        cls.PARALLEL_FILE_LOCKS = None

    @classmethod
    def lockFile(cls, filename):
        """
        Acquire file lock for parallel processing.

        Parameters
        ----------
        filename : str
            Path to the file to lock.
        """
        filename = os.path.normpath(filename)
        # pylint: disable-next=unsupported-membership-test
        if cls.PARALLEL_FILE_LOCKS is not None and filename in cls.PARALLEL_FILE_LOCKS:
            # pylint: disable-next=unsubscriptable-object
            cls.PARALLEL_FILE_LOCKS[filename].acquire()

    @classmethod
    def unlockFile(cls, filename, silence=False):
        """
        Release file lock for parallel processing.

        Parameters
        ----------
        filename : str
            Path to the file to unlock.
        silence : bool, optional
            If True, suppress ValueError when file is not locked.
        """
        filename = os.path.normpath(filename)
        # pylint: disable-next=unsupported-membership-test
        if cls.PARALLEL_FILE_LOCKS is not None and filename in cls.PARALLEL_FILE_LOCKS:
            try:
                # pylint: disable-next=unsubscriptable-object
                cls.PARALLEL_FILE_LOCKS[filename].release()
            except ValueError:
                if not silence:
                    raise

    def __enter__(self):
        """
        Enter context manager.

        Returns
        -------
        PYFT
            Self reference for use in with statement.
        """
        return self

    def __exit__(self, excType, excVal, excTb):
        """
        Exit context manager and close file.
        """
        self.close()

    def close(self):
        """
        Close the FORTRAN file and release resources.

        Prints debug statistics and releases file locks.
        Automatically called when exiting context manager.
        """
        printInfos()
        self.__class__.unlockFile(self.getFileName())

    @property
    def xml(self):
        """
        Get the XML representation of the parsed code.

        Returns
        -------
        str
            XML string representation of the FORTRAN source.
        """
        return tostring(self)

    @property
    def fortran(self):
        """
        Get the FORTRAN source code representation.

        Returns
        -------
        str
            FORTRAN source code string.
        """
        return tofortran(self)

    def renameUpper(self):
        """
        Set output file extension to uppercase.

        Examples
        --------
        >>> pft = PYFT('file.F90')
        >>> pft.renameUpper()  # Output will be file.F90
        """
        self._rename(str.upper)

    def renameLower(self):
        """
        Set output file extension to lowercase.

        Examples
        --------
        >>> pft = PYFT('file.F90')
        >>> pft.renameLower()  # Output will be file.f90
        """
        self._rename(str.lower)

    def _rename(self, mod):
        """
        Apply a transformation function to the file extension.

        Parameters
        ----------
        mod : callable
            Function to apply to file extension (e.g., str.upper, str.lower).
        """
        def _transExt(path, mod):
            filename, ext = os.path.splitext(path)
            return filename + mod(ext)
        if self._output is None:
            self._filename = _transExt(self._filename, mod)
        else:
            self._output = _transExt(self._output, mod)

    def write(self):
        """
        Write the transformed FORTRAN source to file.

        Writes the current state of the code tree as FORTRAN source
        to the output file (or overwrites input if no output specified).
        """
        with open(self._filename if self._output is None else self._output, 'w',
                  encoding='utf-8') as fo:
            fo.write(self.fortran)
            fo.flush()  # ensuring all the existing buffers are written
            os.fsync(fo.fileno())  # forcefully writing to the file

        if self._output is None and self._filename != self._originalName:
            # We must perform an in-place update of the file, but the output file
            # name has been updated. Then, we must remove the original file.
            os.unlink(self._originalName)

    def writeXML(self, filename):
        """
        Write the internal XML representation to a file.

        Parameters
        ----------
        filename : str
            Path for the output XML file.

        Examples
        --------
        >>> pft = PYFT('input.F90')
        >>> pft.writeXML('output.xml')
        """
        with open(filename, 'w', encoding='utf-8') as fo:
            fo.write(self.xml)
