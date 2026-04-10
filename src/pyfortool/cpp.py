"""
This module implements the Cpp class containing the methods for dealing with cpp directives
"""

from pyfortool.util import debugDecor, alltext, PYFTError, tag, noParallel
from pyfortool.tree import updateTree
from pyfortool.variables import updateVarList


class Cpp:
    """
    C preprocessor directive handling methods.

    Provides utilities for evaluating and manipulating C preprocessor
    (#ifdef, #ifndef, #else, #endif) directives in FORTRAN code.
    """

    @debugDecor
    @noParallel
    @updateVarList
    @updateTree('signal')
    def applyCPPifdef(self, keys):
        """
        Evaluate and reduce #ifdef / #ifndef blocks based on specified keys.

        Parameters
        ----------
        keys : list of str
            List of preprocessor keys to evaluate.
            Keys preceded by '%' are treated as undefined (False).
            Example: ['KEY1', '%KEY2'] means KEY1=True, KEY2=False.

        Transformation Examples
        ----------------------
        If 'K' is in keys (K=True):
            #ifdef K
              A
            #else
              B
            #endif
        becomes simply: A

        If '%K' is in keys (K=False):
            #ifdef K
              A
            #else
              B
            #endif
        becomes simply: B

        Notes
        -----
        - Only handles #ifdef and #ifndef (not "#if defined ...")
        - Nested conditionals are supported
        - Code not matching any conditional is kept unchanged
        """

        # We make the hypothesis that #ifdef, #else and #endif have the same parent
        # Get all nodes containing #ifdef or #ifndef
        parents = set(self.getParent(cppNode) for cppNode in self.findall('.//{*}cpp')
                      if (cppNode.text.startswith('#ifdef ') or
                          cppNode.text.startswith('#ifndef ')))

        # Iteration over nodes contained in each parent
        toRemove = []
        for par in parents:
            # we deal with nested #ifdef #ifndef and #if
            # We need to track #if cpp directives to discard #else and #endif related to these #if
            # Each time we enter an #ifdef, #ifndef or #if, we add a value to the keep list
            # True or False to keep or discard it, None not to touch it
            keep = [True]
            for node in par:
                if tag(node) == 'cpp':
                    if node.text.startswith('#ifdef '):
                        k = alltext(node).split(' ')[1].strip()
                        if k in keys:
                            toRemove.append((node, par))
                            keep.append(True)
                        elif '%' + k in keys:
                            toRemove.append((node, par))
                            keep.append(False)
                        else:
                            keep.append(None)
                            if False in keep:
                                toRemove.append((node, par))
                    elif node.text.startswith('#ifndef '):
                        k = alltext(node).split(' ')[1].strip()
                        if k in keys:
                            toRemove.append((node, par))
                            keep.append(False)
                        elif '%' + k in keys:
                            toRemove.append((node, par))
                            keep.append(True)
                        else:
                            keep.append(None)
                            if False in keep:
                                toRemove.append((node, par))
                    elif node.text.startswith('#if '):
                        if False in keep:
                            toRemove.append((node, par))
                        # We are in a #if,following #else / #endif  is associated to this #if
                        keep.append(None)
                    elif node.text.startswith('#else'):
                        if keep[-1] is not None:
                            toRemove.append((node, par))
                            keep[-1] = not keep[-1]
                        elif False in keep:
                            toRemove.append((node, par))
                    elif node.text.startswith('#endif'):
                        if keep[-1] is not None or False in keep:
                            toRemove.append((node, par))
                        keep.pop()
                    elif node.text.startswith('#elifdef') or node.text.startswith('#elifndef'):
                        raise NotImplementedError("#elifdef and #elifndef not (yet?) implemented")
                    else:
                        if False in keep:
                            toRemove.append((node, par))
                else:
                    if False in keep:
                        toRemove.append((node, par))
            if len(keep) != 1:
                # We check the hypothesis done at the beginning
                raise PYFTError("#else or #endif hasn't the same parent as #ifdef " +
                                "or #ifndef in {f}".format(f=self.getFileName()))
        # Suppress node in reverse order to attach tail to previous node
        if len(toRemove) != 0:
            self.tree.signal(self)  # Tree may need to be updated
        for node, par in toRemove[::-1]:
            index = list(par).index(node)
            if index != 0:
                if node.tail is not None:
                    if par[index - 1].tail is None:
                        par[index - 1].tail = ""
                    # We only keep '\n' and spaces at the end (indentation)
                    par[index - 1].tail += (node.tail.count('\n') * '\n' +
                                            (len(node.tail) - len(node.tail.rstrip(' '))) * ' ')
            par.remove(node)
