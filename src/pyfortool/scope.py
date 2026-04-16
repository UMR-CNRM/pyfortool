#!/usr/bin/env python3

"""
Scope-level operations for FORTRAN code.

Provides PYFTscope class for navigating and manipulating FORTRAN scopes
(modules, subroutines, functions, types) with integrated support for
variables, statements, cosmetics, applications, C++ directives, and OpenACC.

Key Features
------------
- Scope path navigation (e.g., 'module:MOD/sub:SUB')
- XML tree traversal with CONTAINS section filtering
- Parent/sibling element lookup with caching
- Integration of Variables, Statements, Cosmetics, Applications, Cpp, Openacc

Classes
-------
PYFTscope : Core scope wrapper extending ElementView
ElementView : XML tree view with optional CONTAINS exclusion

Examples
--------
>>> pft = PYFT('input.F90')
>>> scopes = pft.getScopes()  # Get all scopes
>>> sub = pft.getScopeNode('module:MOD/sub:SUB')  # Get specific scope
>>> for scope in pft.getScopes(excludeKinds=['type']):
...     print(scope.path)
"""

import copy
import os

from pyfortool.variables import Variables, updateVarList
from pyfortool.cosmetics import Cosmetics
from pyfortool.applications import Applications
from pyfortool.statements import Statements
from pyfortool.cpp import Cpp
from pyfortool.openacc import Openacc
from pyfortool.util import PYFTError, debugDecor, n2name, tag
from pyfortool.tree import Tree, updateTree
from pyfortool.expressions import createElem, createExpr


class ElementView():
    """
    View of an ElementTree exposing a subset of subelements.

    Provides filtering capabilities for XML trees, particularly for
    excluding CONTAINS sections from scope traversal.
    """

    def __init__(self, xml, excludeContains=False):
        """
        Initialize ElementView.

        Parameters
        ----------
        xml : xml element
            Root XML element for this view.
        excludeContains : bool, optional
            If True, ignore elements after CONTAINS statement.
        """
        super().__init__()
        self._excludeContains = excludeContains
        self._xml = xml

    @property
    def _virtual(self):
        """
        :param xml: xml corresponding to a scope
        :return: a node (possibly the xml node) containing only the relevant subelements
        """
        if self._excludeContains:
            contains = self._xml.find('./{*}contains-stmt')
            if contains is None:
                return self._xml
            indexContains = list(self._xml).index(contains)
            childNode = createElem('virtual')
            childNode.extend(self._xml[:indexContains] + [self._xml[-1]])
            return childNode
        return self._xml

    # PROPERTIES

    @property
    def tag(self):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.tag
        """
        return self._xml.tag

    @property
    def tail(self):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.tail
        """
        return self._xml.tail

    @property
    def text(self):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.text
        """
        return self._xml.text

    # READ-ONLY METHODS, they can always use the virtual approach

    def findtext(self, *args, **kwargs):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.findtext
        """
        return self._virtual.findtext(*args, **kwargs)

    def iterfind(self, *args, **kwargs):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.iterfind
        """
        return self._virtual.iterfind(*args, **kwargs)

    def itertext(self, *args, **kwargs):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.itertext
        """
        return self._virtual.itertext(*args, **kwargs)

    def __getitem__(self, *args, **kwargs):
        return self._virtual.__getitem__(*args, **kwargs)

    def __len__(self, *args, **kwargs):
        return self._virtual.__len__(*args, **kwargs)

    def __iter__(self):
        return list(self._virtual).__iter__()

    def find(self, *args, **kwargs):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.find
        """
        return self._virtual.find(*args, **kwargs)

    def findall(self, *args, **kwargs):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.findall
        """
        return self._virtual.findall(*args, **kwargs)

    def iter(self, *args, **kwargs):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.iter
        """
        return self._virtual.iter(*args, **kwargs)

    def items(self, *args, **kwargs):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.items
        """
        return self._virtual.items(*args, **kwargs)

    # WRITE METHODS

    @updateVarList
    def clear(self):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.clear
        """
        for item in self:
            self.remove(item)
        self.text = None
        self.tail = None

    def append(self, *args, **kwargs):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.append
        """
        # Append after the 'END SUBROUTINE' statement
        return self._xml.append(*args, **kwargs)

    def extend(self, *args, **kwargs):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.extend
        """
        # Extend after the 'END SUBROUTINE' statement
        return self._xml.extend(*args, **kwargs)

    def _getIndex(self, index):
        """
        :param index: index in the virtual node
        :return: index in the _xml node
        """
        if not self._excludeContains:
            return index
        contains = self._xml.find('./{*}contains-stmt')
        if contains is None:
            return index
        indexContains = list(self._xml).index(contains)
        # Checks
        if index > indexContains or index < -indexContains - 1:
            raise IndexError('list index out of range')
        # Converts negative index into positive index
        if index == -1:
            # END SUBROUTINE
            return index
        if index < -1:
            index = indexContains + index + 1

        return len(self._xml) if index == indexContains else index

    def __setitem__(self, index, item):
        return self._xml.__setitem__(self._getIndex(index), item)

    @updateVarList
    def __delitem__(self, index):
        return self._xml.__delitem__(self._getIndex(index))

    def insert(self, index, item):
        """
        https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element.insert
        """
        return self._xml.insert(0 if index == 0 else (self._getIndex(index - 1) + 1), item)

    @updateVarList
    def remove(self, node):
        """
        Remove node from the xml
        """
        if isinstance(node, ElementView):
            node = node._xml  # pylint: disable=protected-access
        self.getParent(node).remove(node)


class PYFTscope(ElementView, Variables, Cosmetics, Applications, Statements, Cpp, Openacc):
    """
    Wrap an XML node representing a FORTRAN scope.

    PYFTscope provides methods to navigate, query, and modify FORTRAN
    source code at the scope level (modules, subroutines, functions, types).

    Scope Path Format
    -----------------
    Scope paths are '/' separated strings identifying the location in the
    source tree. Examples:
    - 'module:MODULE' - a module
    - 'module:MOD/sub:SUB' - subroutine SUB in module MODULE
    - 'module:MOD/type:TYPE' - type TYPE in module MODULE
    - 'module:MOD/sub:SUB/func:FUNC' - function FUNC in subroutine SUB

    Examples
    --------
    >>> pft = PYFT('myfile.F90')
    >>> scopes = pft.getScopes()  # Get all scopes
    >>> sub = pft.getScopeNode('module:MOD/sub:SUB')  # Get specific scope
    >>> for scope in pft.getScopes(excludeKinds=['type']):
    ...     print(scope.path)
    """
    SCOPE_STMT = {'module': 'module-stmt',
                  'func': 'function-stmt',
                  'sub': 'subroutine-stmt',
                  'type': 'T-stmt',
                  'prog': 'program-stmt',
                  'interface': 'interface-stmt',
                  'submodule': 'submodule-stmt'}
    SCOPE_CONSTRUCT = {'module': 'program-unit',
                       'func': 'program-unit',
                       'sub': 'program-unit',
                       'type': 'T-construct',
                       'prog': 'program-unit',
                       'interface': 'interface-construct',
                       'submodule': 'program-unit'}

    def __init__(self, xml, scopePath='/', parentScope=None,
                 enableCache=False, tree=None, excludeContains=False):
        """
        Initialize a PYFTscope instance.

        Parameters
        ----------
        xml : xml element
            XML element representing the scope.
        scopePath : str, optional
            Path string identifying this scope (e.g., 'module:MOD/sub:SUB').
        parentScope : PYFTscope, optional
            Parent scope instance.
        enableCache : bool, optional
            If True, cache parent nodes for faster traversal.
        tree : Tree, optional
            Tree instance for cross-file analysis.
        excludeContains : bool, optional
            If True, ignore CONTAINS sections in scope traversal.
        """
        super().__init__(xml=xml, excludeContains=excludeContains)
        self._mainScope = self if parentScope is None else parentScope._mainScope
        self._path = scopePath
        self._parentScope = parentScope
        self.tree = Tree() if tree is None else tree
        self._cacheParent = {}

        if enableCache and parentScope is None:
            # parent cache associated to the main scope
            for node in self.iter():
                for subNode in node:
                    self._cacheParent[id(subNode)] = node

    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for key, val in self.__dict__.items():
            setattr(result, key, copy.deepcopy(val, memo))
        return result

    def __getattr__(self, attr):
        """
        Get some attributes defined in the PYFT class
        """
        if attr in ('SHARED_TREE', 'NO_PARALLEL_LOCK', 'PARALLEL_FILE_LOCKS '):
            return getattr(self._parentScope, attr)
        raise AttributeError(f"{attr} doesn't exist")

    @property
    def path(self):
        """
        Get the scope path.

        Returns
        -------
        str
            Scope path string (e.g., 'module:MOD/sub:SUB').
        """
        return self._path

    @property
    def mainScope(self):
        """
        Get the main (root) scope.

        Returns
        -------
        PYFTscope
            The top-level scope in the file.
        """
        return self._mainScope

    @property
    def parentScope(self):
        """
        Get the parent scope.

        Returns
        -------
        PYFTscope or None
            Parent scope, or None if this is the root scope.
        """
        return self._parentScope

    # No @debugDecor for this low-level method
    def getParent(self, item, level=1):
        """
        Get the parent element of an XML node.

        Parameters
        ----------
        item : xml element
            Element whose parent to find.
        level : int, optional
            Number of levels to traverse (1 = direct parent, 2 = grandparent, etc.).

        Returns
        -------
        xml element or None
            Parent element at the specified level.
        """
        # pylint: disable=protected-access
        def check(node):
            # We check if the registered parent is still the right one
            # node must be in its parent, and the parent chain must go to the root node
            return node in self.mainScope._cacheParent.get(id(node), []) and \
                   (self.mainScope._cacheParent[id(node)] == self.mainScope._xml or
                    check(self.mainScope._cacheParent[id(node)]))

        assert level >= 1
        parent = None
        if check(item):
            parent = self.mainScope._cacheParent[id(item)]
        else:
            for node in self.mainScope.iter():
                if item in list(node):
                    parent = node
                    if self.mainScope._cacheParent:
                        # _cacheParent not empty means that we want to use the caching system
                        self.mainScope._cacheParent[id(item)] = node
                    break
        return parent if level == 1 else self.getParent(parent, level - 1)

    def getSiblings(self, item, before=True, after=True):
        """
        Get sibling elements.

        Parameters
        ----------
        item : xml element
            Element whose siblings to find.
        before : bool, optional
            If True, include siblings before the item. Default is True.
        after : bool, optional
            If True, include siblings after the item. Default is True.

        Returns
        -------
        list
            List of sibling elements.

        Examples
        --------
        >>> siblings = scope.getSiblings(node)  # All siblings
        >>> before = scope.getSiblings(node, after=False)  # Only before
        """

        siblings = self.getParent(item).findall('./{*}*')
        if not after:
            siblings = siblings[:siblings.index(item)]
        if not before:
            siblings = siblings[siblings.index(item) + 1:]
        return [s for s in siblings if s != item]

    # No @debugDecor for this low-level method
    @staticmethod
    def _getNodeName(node):
        """
        Internal methode to compute the name of a scope
        :param node: program-unit node
        :return: name
        """
        if tag(node) == 'T-stmt':
            # To not capture the extension name in "TYPE, EXTENDS(FOO) :: FOO2"
            nodeN = node.find('./{*}T-N/{*}N')
        elif tag(node) == 'submodule-stmt':
            nodeN = node.find('./{*}submodule-module-N')
        else:
            nodeN = node.find('.//{*}N')
            if nodeN is not None and nodeN.find('.//{*}N') is not None:
                # As of 7 Jul 2023, this is the case for interface statements
                nodeN = nodeN.find('.//{*}N')
        if nodeN is not None:
            name = n2name(nodeN).upper()
        else:
            name = '--UNKNOWN--'
        return name

    # No @debugDecor for this low-level method
    def _getNodePath(self, node):
        """
        Internal methode to compute a path part from a node
        :param node: program-unit node
        :return: path part (e.g. module:MODU)
        """
        stmt = tag(node[0])
        name = self._getNodeName(node[0])
        return {v: k for (k, v) in self.SCOPE_STMT.items()}[stmt] + ':' + name

    # No @debugDecor for this low-level method
    @staticmethod
    def normalizeScope(scopePath):
        """
        Normalize a scope path to standard format.

        Converts scope path to lowercase prefix and uppercase names.

        Parameters
        ----------
        scopePath : str
            Scope path to normalize.

        Returns
        -------
        str
            Normalized scope path.

        Examples
        --------
        >>> PYFTscope.normalizeScope('module:Test/sub:Sub')
        'module:TEST/sub:SUB'
        """
        return '/'.join([(k.lower() + ':' + w.upper())
                         for (k, w) in [component.split(':')
                         for component in scopePath.split('/')]])

    @debugDecor
    def showScopesList(self, includeItself=False):
        """
        Display all scopes found in the source code.

        Parameters
        ----------
        includeItself : bool, optional
            If True, include the current scope in the output.
            Default is False.

        Examples
        --------
        >>> pft = PYFT('myfile.F90')
        >>> pft.showScopesList()
        These scopes have been found in the source code:
          - /module:MOD
          - /module:MOD/sub:SUB
          - /module:MOD/func:FUNC
        """
        print("These scopes have been found in the source code:")
        print("\n".join(['  - ' + scope.path
                         for scope in self.getScopes(includeItself=includeItself)]))

    @debugDecor
    def getScopes(self, level=-1, excludeContains=True, excludeKinds=None, includeItself=True):
        """
        Get child scopes from the current scope.

        Parameters
        ----------
        level : int, optional
            Depth of scope traversal:
            - -1 (default): All child scopes recursively.
            - 1: Direct children only.
            - 2: Children and grandchildren, etc.
        excludeContains : bool, optional
            If True (default), exclude scopes from CONTAINS sections
            (nested subroutines/functions).
        excludeKinds : list of str, optional
            Scope kinds to exclude. Options: 'module', 'sub', 'func',
            'type', 'prog', 'interface', 'submodule'.
        includeItself : bool, optional
            If True (default), include the current scope in results.

        Returns
        -------
        list of PYFTscope
            List of scope instances matching the criteria.

        Examples
        --------
        >>> pft = PYFT('myfile.F90')
        >>> all_scopes = pft.getScopes()  # All scopes
        >>> subs = pft.getScopes(excludeKinds=['module', 'func', 'type'])
        >>> direct = pft.getScopes(level=1)
        """
        assert level == -1 or level > 0, 'level must be -1 or a positive int'

        def _getRecur(node, level, basePath=''):
            # If node is the entire xml
            if tag(node) == 'object':
                usenode = node.find('./{*}file')
            else:
                usenode = node
            results = []
            for child in [child for child in usenode
                          if tag(child) in self.SCOPE_CONSTRUCT.values()]:
                nodePath = self._getNodePath(child)
                if excludeKinds is None or nodePath.split(':')[0] not in excludeKinds:
                    scopePath = nodePath if basePath in ('', '/') \
                                else basePath + '/' + nodePath
                    results.append(PYFTscope(child,
                                             scopePath=scopePath, parentScope=self,
                                             enableCache=False, tree=self.tree,
                                             excludeContains=excludeContains))
                    if level != 1:
                        results.extend(_getRecur(child, level - 1, scopePath))
            return results

        if includeItself and tag(self) in self.SCOPE_CONSTRUCT.values():
            itself = [self]
        else:
            itself = []

        return _getRecur(self, level, self.path) + itself

    @debugDecor
    def getScopeNode(self, scopePath, excludeContains=True, includeItself=True):
        """
        Get a specific scope by path.

        Parameters
        ----------
        scopePath : str
            Scope path to search for (e.g., 'module:MOD/sub:SUB').
        excludeContains : bool, optional
            See getScopes. Default is True.
        includeItself : bool, optional
            See getScopes. Default is True.

        Returns
        -------
        PYFTscope
            Scope instance matching the path.

        Raises
        ------
        PYFTError
            If scope not found or found multiple times.

        Examples
        --------
        >>> pft = PYFT('myfile.F90')
        >>> sub = pft.getScopeNode('module:MOD/sub:SUB')
        >>> func = pft.getScopeNode('/module:MOD/func:FUNC')
        """
        scope = [scope for scope in self.getScopes(excludeContains=excludeContains,
                                                   includeItself=includeItself)
                 if scope.path == scopePath]
        if len(scope) == 0:
            raise PYFTError(f'{scopePath} not found')
        if len(scope) > 1:
            raise PYFTError(f'{scopePath} found several times')
        return scope[0]

    @debugDecor
    def isScopeNode(self, node):
        """
        Check if a node is a scope-level construct.

        Parameters
        ----------
        node : xml element
            Node to check.

        Returns
        -------
        bool
            True if node is a scope construct (program-unit, interface-construct,
            or T-construct).

        Examples
        --------
        >>> node = pft.find('.//{*}subroutine-stmt/..')
        >>> pft.isScopeNode(node)
        True
        """
        return tag(node) in self.SCOPE_CONSTRUCT.values()

    @debugDecor
    def getParentScopeNode(self, item, mustRaise=True):
        """
        Get the scope containing an element.

        Parameters
        ----------
        item : xml element
            Element whose containing scope to find.
        mustRaise : bool, optional
            If True (default), raise PYFTError if scope not found.

        Returns
        -------
        xml element or None
            The scope node containing the item.

        Examples
        --------
        >>> call = pft.find('.//{*}call-stmt')
        >>> scope = pft.getParentScopeNode(call)
        """
        result = self.getParent(item)
        while result is not None and not self.isScopeNode(result):
            result = self.getParent(result)
        if result is None and mustRaise:
            raise PYFTError("The scope parent has not been found.")
        return result

    @debugDecor
    def getScopePath(self, item, includeItself=True):
        """
        Get the scope path for an element.

        Parameters
        ----------
        item : xml element
            Element whose scope path to determine.
        includeItself : bool, optional
            If True (default) and item is a scope node, include it.

        Returns
        -------
        str
            Full scope path of the element's containing scope.

        Examples
        --------
        >>> node = pft.find('.//{*}a-stmt')
        >>> pft.getScopePath(node)
        '/module:MOD/sub:SUB'
        """
        if includeItself and self.isScopeNode(item):
            result = [self._getNodePath(item)]
        else:
            result = []
        item = self.getParentScopeNode(item, mustRaise=False)
        while item is not None:
            result = [self._getNodePath(item)] + result
            item = self.getParentScopeNode(item, mustRaise=False)
        return '/'.join(result)

    def getFileName(self):
        """
        Get the source filename.

        Returns
        -------
        str
            Normalized path to the source file.

        Examples
        --------
        >>> pft = PYFT('/path/to/file.F90')
        >>> pft.getFileName()
        '/path/to/file.F90'
        """
        return os.path.normpath(self.mainScope.find('.//{*}file').attrib['name'])

    @updateTree()
    @updateVarList
    def empty(self, addStmt=None, simplify=False):
        """
        Remove all statements from scopes.

        Removes all executable statements from scopes while preserving:
        - Dummy argument declarations
        - USE statements
        - Scope declarations and endings

        Parameters
        ----------
        addStmt : str or list, optional
            Statement(s) to insert into the body of emptied scopes.
        simplify : bool, optional
            If True, also remove unused local variables.

        Examples
        --------
        >>> pft = PYFT('myfile.F90')
        >>> pft.empty()  # Remove all statements

        Empty and add a comment:
        >>> pft.empty(addStmt='! TODO: Implement')
        """
        scopes = []  # list of scopes to empty
        for scope in self.getScopes(level=1, excludeContains=False):
            if scope.path.split('/')[-1].split(':')[0] == 'module':
                scopes.extend(scope.getScopes(level=1, excludeContains=False, includeItself=False))
            else:
                scopes.append(scope)
        tagExcluded = (list(self.SCOPE_STMT.values()) +
                       ['end-' + decl for decl in self.SCOPE_STMT.values()] +
                       ['T-decl-stmt', 'use-stmt', 'C'])
        for scope in scopes:
            for node in list(scope):
                if tag(node) not in tagExcluded:
                    scope.remove(node)
            scope.removeUnusedLocalVar(simplify=simplify)
            if addStmt is not None:
                if isinstance(addStmt, str):
                    addStmt = createExpr(addStmt)
                elif not isinstance(addStmt, list):
                    addStmt = [addStmt]
                for stmt in addStmt:
                    scope.insertStatement(stmt, False)
        if simplify:
            # Apllied on self (and not scope) to remove lines before the first scope
            # in case self represents an entire file
            self.removeComments()
            self.removeEmptyLines()
