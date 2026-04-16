"""
Expression manipulation functions.

These functions are independent of PYFT and PYFTscope objects and provide
low-level utilities for creating and manipulating FORTRAN expression XML nodes.
"""

import re
from functools import lru_cache
import copy
import xml.etree.ElementTree as ET

from pyfortool.util import debugDecor, isint, isfloat, fortran2xml, PYFTError
from pyfortool import NAMESPACE


def createElem(tagName, text=None, tail=None, childs=None):
    """
    Create an XML element with the given tag and attributes.

    Parameters
    ----------
    tagName : str
        XML tag name (without namespace).
    text : str, optional
        Text content for the element.
    tail : str, optional
        Tail text (text after element).
    childs : Element or list, optional
        Child element(s) to append.

    Returns
    -------
    Element
        Created XML element.

    Examples
    --------
    >>> elem = createElem('named-E')
    >>> elem = createElem('literal-E', text='42')
    >>> elem = createElem('n', text='X', tail='\\n')
    """
    node = ET.Element(f'{{{NAMESPACE}}}{tagName}')
    if text is not None:
        node.text = text
    if tail is not None:
        node.tail = tail
    if childs is not None:
        if isinstance(childs, list):
            node.extend(childs)
        else:
            node.append(childs)
    return node


@lru_cache
def _cachedCreateExprPart(value):
    """
    :param value: expression part to put in a *-E node

    If value is:
      - a FORTRAN string (python sting containing a ' or a "), returns
        <f:string-E><f:S>...
      - a FORTRAN value (python string convertible in real or int, or .FALSE./.TRUE.), returns
        <f:literal-E><f:l>...
      - a FORTRAN variable name (pyhon string with only alphanumerical characters and _), returns
        <named-E/><N><n>...
      - a FORTRAN operation (other python string), returns the right part of
        the X affectation statement of the code:
        "SUBROUTINE T; X=" + value + "; END". The xml is obtained by calling fxtran.
    """

    # Allowed characters in a FORTRAN variable name
    allowed = "abcdefghijklmnopqrstuvwxyz"
    allowed += allowed.upper() + '0123456789_'

    if isint(value) or isfloat(value) or value.upper() in ('.TRUE.', '.FALSE.'):
        node = createElem('literal-E')
        node.append(createElem('l', text=str(value)))
    elif "'" in value or '"' in value:
        node = createElem('string-E')
        node.append(createElem('S', text=value))
    elif all(c in allowed for c in value):
        nodeN = createElem('N')
        nodeN.append(createElem('n', text=value))
        node = createElem('named-E')
        node.append(nodeN)
    elif re.match(r'[a-zA-Z_][a-zA-Z0-9_]*%[a-zA-Z_][a-zA-Z0-9_]*$', value):
        # A%B
        nodeN = createElem('N')
        nodeN.append(createElem('n', text=value.split('%')[0]))
        ct = createElem('ct', text=value.split('%')[1])
        componentR = createElem('component-R', text='%')
        componentR.append(ct)
        nodeRLT = createElem('R-LT')
        nodeRLT.append(componentR)
        node = createElem('named-E')
        node.append(nodeN)
        node.append(nodeRLT)
    else:
        _, xml = fortran2xml(f"SUBROUTINE T; X={value}; END")
        node = xml.find('.//{*}E-2')[0]
    return node


@debugDecor
def createExprPart(value):
    """
    Create an XML node from a FORTRAN expression part.

    Parameters
    ----------
    value : str
        Expression part value to convert.

    Returns
    -------
    Element
        XML element representing the expression:
        - Integer/float: <literal-E><l>value</l></literal-E>
        - String: <string-E><S>value</S></string-E>
        - Variable: <named-E><N><n>value</n></N></named-E>
        - Structure member:
          <named-E><N><n>A</n></N><R-LT><component-R>%B</component-R></R-LT></named-E>
        - Expression: parsed via fxtran

    Examples
    --------
    >>> createExprPart('42')  # Literal
    >>> createExprPart('X')   # Variable
    >>> createExprPart('A%B') # Structure member
    """
    return copy.deepcopy(_cachedCreateExprPart(value))


@lru_cache
def _cachedCreateExpr(value):
    """
    Internal cached function for createExpr.

    Parameters
    ----------
    value : str
        FORTRAN statement(s) to convert.

    Returns
    -------
    list
        List of XML nodes from the statement.
    """
    return fortran2xml(f"SUBROUTINE T\n{value}\nEND")[1].find('.//{*}program-unit')[1:-1]


@debugDecor
def createExpr(value):
    """
    Convert FORTRAN statements to XML nodes.

    Parameters
    ----------
    value : str
        One or more FORTRAN statements to convert.

    Returns
    -------
    list
        List of XML nodes representing the statements.

    Examples
    --------
    >>> nodes = createExpr('X = 42')
    >>> nodes = createExpr('CALL SUB(X, Y)')
    >>> nodes = createExpr('IF (A > B) THEN\\n  X = 1\\nEND IF')
    """
    return copy.deepcopy(_cachedCreateExpr(value))


@debugDecor
def simplifyExpr(expr, add=None, sub=None):
    """
    Simplify a numeric expression by combining constants.

    Parameters
    ----------
    expr : str
        Expression to simplify (e.g., '1+I+2+JI-I').
    add : str, optional
        Expression to add to the result.
    sub : str, optional
        Expression to subtract from the result.

    Returns
    -------
    str
        Simplified expression string.

    Examples
    --------
    >>> simplifyExpr('1+1+I+JI-I')
    '2+JI'
    >>> simplifyExpr('X+1', add='Y')
    'X+Y+1'

    Notes
    -----
    - Only handles addition and subtraction.
    - Does not simplify expressions within parentheses.
    """
    # We could have used external module, such as sympy, but this routine
    # (as long as it's sufficient) avoids introducing dependencies.
    if re.search(r'\([^()]*[+-][^()]*\)', expr):
        raise NotImplementedError("Expression cannot (yet) contain + or - sign inside " +
                                  f"parenthesis: {expr}")

    def split(expr):
        """
        :param s: expression
        :return: a list of (sign, abs(value))
        """
        # splt is ['1', '+', '1', '+', 'I', '+', 'JI', '-', 'I']
        splt = re.split('([+-])', expr.replace(' ', '').upper())
        if splt[0] == '':
            # '-1' returns [
            splt = splt[1:]
        if len(splt) % 2 == 1:
            # expr doesn't start with a sign
            splt = ['+'] + splt  # ['+', '1', '+', '1', '+', 'I', '+', 'JI', '-', 'I']
        # group sign and operand [('+', '1'), ('+', '1'), ('+', 'I'), ('+', 'JI'), ('-', 'I')]
        splt = [(splt[2 * i], splt[2 * i + 1]) for i in range(len(splt) // 2)]
        return splt

    splt = split(expr)
    if add is not None:
        splt += split(add)
    if sub is not None:
        splt += [('-' if sign == '+' else '+', elem) for (sign, elem) in split(sub)]
    # Suppress elements with opposite signs
    for sign, elem in splt.copy():
        if ('+', elem) in splt and ('-', elem) in splt:
            splt.remove(('+', elem))
            splt.remove(('-', elem))
    # Pre-compute integer additions/substractions
    found = -1
    for i, (sign, elem) in enumerate(splt.copy()):
        if isint(elem):
            if found == -1:
                found = i
            else:
                result = str((1 if splt[found][0] == '+' else -1) * int(splt[found][1]) +
                             (1 if sign == '+' else -1) * int(elem))
                splt[found] = split(str(result))[0]
                splt.pop(i)
    # Order (no matter what ordering is done but we need to order to allow comparisons)
    splt.sort(key=''.join)
    # Empty e.g. '1-1'
    if len(splt) == 0:
        splt = [('+', '0')]
    # Concatenate
    result = ' '.join(s[0] + ' ' + s[1] for s in splt)
    if result.startswith('+'):
        result = result[1:]
    return result.lstrip(' ')


@debugDecor
def createArrayBounds(lowerBoundstr, upperBoundstr, context):
    """
    Return a lower-bound and upper-bound node
    :param lowerBoundstr: string for the fortran lower bound of an array
    :param upperBoundstr: string for the fortran upper bound of an array
    :param context: 'DO' for DO loops
                    'DOCONCURRENT' for DO CONCURRENT loops
                    'ARRAY' for arrays
    """
    lowerBound = createElem('lower-bound')
    lowerBound.insert(0, createExprPart(lowerBoundstr))
    upperBound = createElem('upper-bound')
    upperBound.insert(0, createExprPart(upperBoundstr))
    if context == 'DO':
        lowerBound.tail = ', '
    elif context in ('DOCONCURRENT', 'ARRAY'):
        lowerBound.tail = ':'
    else:
        raise PYFTError(f'Context unknown in createArrayBounds: {context}')
    return lowerBound, upperBound
