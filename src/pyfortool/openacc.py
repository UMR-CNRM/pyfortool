"""
This module implements the Openacc class containing the methods relative to openacc
"""

import re

from pyfortool.util import debugDecor, n2name, alltext, tag
from pyfortool.expressions import createElem, createExpr


class Openacc():
    """
    Methods relative to openacc
    """
    @debugDecor
    def removeACC(self):
        """
        Remove openACC directives
        """
        self.removeComments(exclDirectives=[],
                            pattern=re.compile(r'^\!\$ACC', re.IGNORECASE))

    @debugDecor
    def removebyPassDOCONCURRENT(self):
        """
        Remove macro !$mnh_(un)def(OPENACC) and !$mnh_(un)def(LOOP) directives
        for other compiler than Cray
        """

        self.removeComments(exclDirectives=[],
                            pattern=re.compile(r'^\!\$mnh_undef\(LOOP\)'))
        self.removeComments(exclDirectives=[],
                            pattern=re.compile(r'^\!\$mnh_undef\(OPENACC\)'))
        self.removeComments(exclDirectives=[],
                            pattern=re.compile(r'^\!\$mnh_define\(LOOP\)'))
        self.removeComments(exclDirectives=[],
                            pattern=re.compile(r'^\!\$mnh_define\(OPENACC\)'))

    @debugDecor
    def craybyPassDOCONCURRENT(self):
        """
        By pass a bug of the CRAY compiler in which the vectorisation is not done with
        BR_ fonctions use or locally.
        On all expanded compute kernels with !$acc loop independent collapse(X) placed:
            - if BR_ fonction is used : !$acc loop independent collapse(X) is removed and
              the nested DO loops are factorised into DO CONCURRENT
            - if a mnh_undef(OPENACC) macro is in place, !$acc loop collapse independant(X)
              is removed
            - if a mnh_undef(LOOP) macro is in place the nested DO loops are factorised into
              DO CONCURRENT
        """
        def checkPresenceofBR(node):
            """Return True if a BR_ (math BIT-REPRODUCTIBILITY) function is present in the node"""
            mathBRList = ['ALOG', 'LOG', 'EXP', 'COS', 'SIN', 'ASIN', 'ATAN', 'ATAN2',
                          'P2', 'P3', 'P4']
            namedE = node.findall('.//{*}named-E/{*}N/{*}n')
            for el in namedE:
                if alltext(el) in ['BR_' + e for e in mathBRList]:
                    return True
            return False

        def getStatementsInDoConstruct(node, savelist):
            for sNode in node:
                if 'do-construct' in tag(sNode):
                    getStatementsInDoConstruct(sNode, savelist)
                elif 'do-stmt' not in tag(sNode):
                    savelist.append(sNode)

        useNestedLoops = True  # Cray compiler needs nested loops by default
        useAccLoopIndependent = True  # It needs nested !$acc loop independent collapse(X)
        toremove = []  # list of nodes to remove
        comments = self.findall('.//{*}C')

        for comment in comments:
            if comment.text.startswith('!$mnh_undef(LOOP)'):
                useNestedLoops = False  # Use DO CONCURRENT
            if comment.text.startswith('!$mnh_undef(OPENACC)'):
                useAccLoopIndependent = False

            if comment.text.startswith('!$mnh_define(LOOP)'):
                useNestedLoops = True  # Use DO CONCURRENT
            if comment.text.startswith('!$mnh_define(OPENACC)'):
                useAccLoopIndependent = True

            if comment.text.startswith('!$acc loop independent collapse('):
                # Get the statements content in the DO-construct
                par = self.getParent(comment)
                ind = list(par).index(comment)
                nestedLoop = par[ind + 1]
                statements = []
                getStatementsInDoConstruct(nestedLoop, statements)

                # Check presence of BR_ within the statements
                isBRPresent = False
                for stmt in statements:
                    if checkPresenceofBR(stmt):
                        isBRPresent = True
                        break

                # Remove !$acc loop independent collapse
                # if BR_ is present or if !$mnh_undef(OPENACC)
                if not useAccLoopIndependent or isBRPresent:
                    toremove.append((self, comment))

                # Use DO CONCURRENT instead of nested-loop if BR_ is present or if !$mnh_undef(LOOP)
                if not useNestedLoops or isBRPresent:
                    # Determine the table of indices
                    doStmt = nestedLoop.findall('.//{*}do-stmt')
                    table = {}
                    for do in reversed(doStmt):
                        table[do.find('.//{*}do-V/{*}named-E/{*}N/{*}n').text] = \
                            [alltext(do.find('.//{*}lower-bound')),
                             alltext(do.find('.//{*}upper-bound'))]

                    # Create the do-construct
                    inner, outer, _ = self.createDoConstruct(table,
                                                             indent=len(nestedLoop.tail),
                                                             concurrent=True)

                    # Insert the statements in the new do-construct
                    for stmt in statements:
                        inner.insert(-1, stmt)

                    # Insert the new do-construct and delete all the old do-construct
                    par.insert(ind, outer)
                    toremove.append((par, nestedLoop))

        # Suppression of nodes
        for parent, elem in toremove:
            parent.remove(elem)

    @debugDecor
    def allocatetoHIP(self):
        """
        Convert (DE)ALLOCATE to (DE)ALLOCATE_HIP on variables only sent to the GPU via
        !$acc enter data copyin
        or
        !$acc enter data create
        This is necessary for using the managed memory with GPU AMD MI250X (on Adastra)
        """
        scopes = self.getScopes()
        for scope in scopes:
            varsToChange = []
            comments = scope.findall('.//{*}C')
            pointers = scope.findall('.//{*}pointer-a-stmt')

            for coms in comments:
                if ('!$acc enter data' in coms.text or '!$acc exit data' in coms.text) \
                        and coms.text.count('!') == 1:
                    if coms.text.count('(') == 1:
                        # !$acc enter data copyin( XRRS, XRRS_CLD ) ==> [' XRRS, XRRS_CLD ']
                        varsToChange.extend(coms.text.split(')')[0].split('(')[1:][0].split(','))
                    else:
                        # $acc exit data delete(xb_mg(level,m)%st) ==> xb_mg(level,m)%st
                        variableName = re.search(r'\(.*\)', coms.text).group(0)[1:-1]
                        varsToChange.insert(0, variableName)

            if len(pointers) > 0:
                for i, var in enumerate(varsToChange):
                    for pointer in pointers:
                        if alltext(pointer.find('.//{*}E-1/{*}named-E')) == var:
                            varsToChange[i] = alltext(pointer.find('.//{*}E-2/{*}named-E'))

            for i, var in enumerate(varsToChange):
                varsToChange[i] = var.replace(' ', '')

            if len(varsToChange) > 0:
                scope.addModuleVar([(scope.path, 'MODE_MNH_HIPFORT', None)])

            allocateStmts = scope.findall('.//{*}allocate-stmt')
            allocateStmts.extend(scope.findall('.//{*}deallocate-stmt'))

            for stmt in allocateStmts:
                varsChecking = ''
                allocateArg = stmt.find('.//{*}arg-spec')
                arrayR = allocateArg.find('.//{*}array-R')
                if arrayR is None:
                    # Particular case of multiple parensR such as Tjacobi(level)%r(nz)
                    parensR = allocateArg.findall('.//{*}parens-R')
                    if len(parensR) == 2:
                        # It's a component case
                        varsChecking = alltext(allocateArg).split('%')[0]+'%' + \
                            alltext(allocateArg).split('%')[1].split('(')[0]
                    elif len(parensR) == 1:
                        # Tjacobi(level)%Sr case where Sr is a type itself
                        if allocateArg.find('.//{*}component-R'):
                            varsChecking = alltext(allocateArg)
                        else:
                            varsChecking = alltext(allocateArg).split('(')[0]
                    elif len(parensR) == 0:
                        varsChecking = alltext(allocateArg).split('(')[0]
                else:
                    varsChecking = alltext(allocateArg).split('(')[0]

                if varsChecking in varsToChange:
                    removeLastComma = True
                    stmt.text = "CALL MNH_HIP" + stmt.text  # For allocate/deallocate statements
                    if tag(stmt) == 'allocate-stmt':
                        # Replace upper:lower into upper,lower
                        lowerBounds = allocateArg.findall('.//{*}lower-bound')
                        if lowerBounds is not None:
                            for bounds in lowerBounds:
                                bounds.tail = ','
                        arrayR = allocateArg.find('.//{*}array-R')
                        if arrayR is None:
                            # Particular case of multiple parensR such as Tjacobi(level)%r(nz)
                            parensR = allocateArg.findall('.//{*}parens-R')
                            if len(parensR) > 1:
                                parensR[-1].text = ','
                            else:
                                if allocateArg.find('.//{*}component-R'):
                                    removeLastComma = False
                                else:
                                    parensR[0].text = ','
                        else:
                            arrayR.text = ','
                        if removeLastComma:
                            allocateArg.tail = ''

    @debugDecor
    def addACCData(self):
        """
        1) Add after declaration:
        !$acc data present ( list of intent arrays)
        2) Add at the end of the routine
        !$acc end data
        """
        scopes = self.getScopes()
        if scopes[0].path.split('/')[-1].split(':')[1][:4] == 'MODD':
            return
        for scope in scopes:
            # Do not add !$acc data directives to :
            # - MODULE or FUNCTION object,
            # - interface subroutine from a MODI
            # but only to SUBROUTINES
            if 'sub:' in scope.path and 'func' not in scope.path and 'interface' not in scope.path:
                # Look for all intent arrays only
                arraysIntent = []
                for var in scope.varList:
                    # intent arrays, not of type TYPE (only REAL, INTEGER, CHARACTER)
                    if var['arg'] and var['as'] and 'TYPE' not in var['t'] and \
                       var['scopePath'] == scope.path:
                        arraysIntent.append(var['n'])
                # Check if there is any intent variables
                if len(arraysIntent) == 0:
                    break

                # 1) First !$acc data present()
                listVar = "!$acc data present ( "
                count = 0
                for var in arraysIntent:
                    if count > 6:
                        listVar = listVar + '\n!$acc &              '
                        count = 0
                    listVar = listVar + var + ", &"
                    count += 1
                listVarEnd = listVar[:-3]  # remove last comma and &
                accAddMultipleLines = createExpr(listVarEnd + ')')
                idx = scope.insertStatement(scope.indent(accAddMultipleLines[0]), first=True)

                # 2) multi-lines !$acc &
                for iLine, line in enumerate(accAddMultipleLines[1:]):
                    scope.insert(idx + 1 + iLine, line)

                # 3) !$acc end data
                comment = createElem('C', text='!$acc end data', tail='\n')
                scope.insertStatement(scope.indent(comment), first=False)

    @debugDecor
    def addACCRoutineSeq(self, stopScopes):
        """
        Adds the '!$acc routine (<name>) seq' directive
        :param stopScopes: scope paths where we stop to add the directive
        """
        for scope in self.getScopes():
            if self.tree.isUnderStopScopes(scope.path, stopScopes,
                                           includeInterfaces=True,
                                           includeStopScopes=True):
                name = n2name(scope[0].find('.//{*}N')).upper()
                acc = createElem('C', text=f'!$acc routine ({name}) seq',
                                 tail=scope[0].tail)
                scope[0].tail = '\n'
                scope.insert(1, acc)
