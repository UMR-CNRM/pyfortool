# Core Concepts

This document explains the key concepts in PyForTool.

## Scope Paths

Scope paths identify locations in FORTRAN source code. They're used to specify where transformations should be applied.

### Format

```
module:MODULE/sub:SUB/func:FUNC/type:TYPE
```

Components are separated by `/` and each component has the form `kind:NAME`.

### Component Types

| Prefix | Meaning | Example |
|--------|---------|---------|
| `module:` | FORTRAN module | `module:MODULE_NAME` |
| `sub:` | Subroutine | `sub:SUB_NAME` |
| `func:` | Function | `func:FUNC_NAME` |
| `type:` | Derived type definition | `type:TYPE_NAME` |
| `prog:` | Program | `prog:MAIN` |
| `interface:` | Interface block | `interface:IF_NAME` |
| `submodule:` | Submodule | `submodule:PARENT_CHILD` |

### Examples

```python
pft = PYFT('input.F90')

# Get module
mod = pft.getScopeNode('module:MOD')

# Get subroutine in module
sub = pft.getScopeNode('module:MOD/sub:SUB')

# Get function in subroutine
func = pft.getScopeNode('module:MOD/sub:SUB/func:INNER')

# Get type definition in module
tp = pft.getScopeNode('module:MOD/type:MYTYPE')

# Get all subroutines (any module)
for scope in pft.getScopes(excludeKinds=['module', 'func', 'type']):
    print(scope.path)
```

### Path Normalization

```python
# Paths are normalized to lowercase prefix, uppercase name
PYFTscope.normalizeScope('module:Test/sub:Sub')
# Returns: 'module:TEST/sub:SUB'
```

### Path Matching

```python
# Get scopes under a path
pft.getScopes(level=1)  # Direct children only
pft.getScopes(level=-1)  # All descendants (default)

# Check if path is under another
tree.isUnderStopScopes('module:MOD/sub:SUB', 
                      stopScopes=['module:MOD'])
# Returns: True
```

## VarList

VarList stores information about variables declared in a scope.

### Variable Descriptor Dictionary

```python
var = {
    'n': 'X',              # Variable name (uppercase)
    't': 'REAL',           # Type specification
    'as': [],              # Array specifications
    'i': 'IN',             # INTENT attribute
    'arg': True,           # Is dummy argument
    'argorder': 0,         # Position in argument list
    'use': 'MODULE',       # Module if imported via USE
    'opt': False,          # OPTIONAL attribute
    'allocatable': False,  # ALLOCATABLE attribute
    'pointer': False,      # POINTER attribute
    'parameter': False,    # PARAMETER attribute
    'result': False,       # Function result variable
    'init': None,          # Initial value
    'scopePath': 'module:MOD/sub:SUB'  # Declaration scope
}
```

### Array Specifications

```python
# 1D array 1:10
var['as'] = [(None, '10')]

# 2D array 1:10, 1:20
var['as'] = [(None, '10'), (None, '20')]

# Explicit bounds
var['as'] = [('1', '10'), ('LOWER', 'UPPER')]
```

### Accessing Variables

```python
vl = pft.varList

# Find any variable
var = vl.findVar('X')

# Find only arrays
array = vl.findVar('Y', array=True)

# Find only scalars
scalar = vl.findVar('Z', array=False)

# Find in specific scope
var = vl.findVar('X', exactScope=True)

# Restrict to a scope
sub_vl = vl.restrict('module:MOD/sub:SUB', excludeContains=True)
```

## XML Structure

PyForTool uses fxtran to convert FORTRAN to XML.

### Namespace

```python
NAMESPACE = 'http://fxtran.net/#syntax'
```

### Common Element Tags

| Tag | FORTRAN Element |
|-----|-----------------|
| `subroutine-stmt` | SUBROUTINE statement |
| `end-subroutine-stmt` | END SUBROUTINE |
| `function-stmt` | FUNCTION statement |
| `T-decl-stmt` | Type declaration |
| `T-stmt` | TYPE definition |
| `call-stmt` | CALL statement |
| `a-stmt` | Assignment |
| `use-stmt` | USE statement |
| `C` | Comment |
| `if-construct` | IF block |
| `do-construct` | DO loop |

### Querying XML

```python
# Find single element
node = pft.find('.//{*}subroutine-stmt')

# Find all matching elements
nodes = pft.findall('.//{*}call-stmt')

# XPath with namespace
nodes = pft.findall('.//{*}named-E/{*}N/{*}n')
```

### Element Properties

```python
elem.tag      # Tag name (without namespace)
elem.text     # Text content
elem.tail     # Text after element
elem.attrib   # Attributes dictionary
elem.find()   # Find child
elem.findall() # Find all children
elem.iter()   # Iterate all descendants
```

## Tree (Cross-File Dependencies)

Tree tracks compilation and execution dependencies.

### Creating a Tree

```python
from pyfortool.tree import Tree

# Scan directories
tree = Tree(['/path/to/src'], descTreeFile='tree.json')

# Load existing tree
tree = Tree(tree='tree.json')
```

### Data Structures

```python
# Files in tree
tree.getFiles()

# Compilation dependencies (USE statements)
tree._useList[file][scopePath] = [(module, onlyList), ...]

# Execution dependencies (CALL statements)
tree._callList[file][scopePath] = [callName, ...]

# Functions
tree._funcList[file][scopePath] = [funcName, ...]

# Includes
tree._includeList[file][scopePath] = [include, ...]
```

### Querying Dependencies

```python
# What does file.F90 need?
deps = tree.needsFile('file.F90')
# Returns: {'needs': {file: [...scopes...]}, 'needsFromDir': [...]}

# Who needs file.F90?
users = tree.usedByFile('file.F90')
# Returns: {'usedBy': {file: [...scopes...]}, 'usedByFromDir': [...]}

# Who calls SUB?
callers = tree.whoCall('SUB')

# What scopes call SCOPE?
callees = tree.scopeCallees('module:MOD/sub:SUB')
```

### Visualization

```python
# Generate DOT file for dependency graph
tree.plotCompilTreeFromFile('main.F90', 'deps.dot')

# With limits
tree.plotExecTreeFromFile('main.F90', 'calls.png', 
                          maxUpper=2, maxLower=2)
```

## Decorators

### @debugDecor

Located in `util.py`. Provides tracing and profiling.

```python
@debugDecor
def myMethod(self, arg):
    """Method documentation."""
    pass

# With logging level INFO:
# - Shows call count and execution time
# INFO: Function myMethod called 42 times, total time: 0.12s

# With logging level DEBUG:
# - Shows arguments and return values
# DEBUG: myMethod(arg='value') -> {'result': 42}
```

**Performance Note:** Low overhead except when called many times in tight loops. Consider removing for performance-critical code.

### @updateVarList

Located in `variables.py`. Invalidates variable cache.

```python
@updateVarList
def modifyVariables(self):
    # Changes to declarations
    self.addVar([...])
```

**Important:** Must be applied to any method that:
- Adds or removes variables
- Modifies variable declarations
- Changes argument lists

### @updateTree

Located in `tree.py`. Updates dependency tracking.

```python
@updateTree('file')    # Current file only
@updateTree('scan')    # Scan for new/removed files
@updateTree('signal')  # Process signaled files
```

**Usage:** Apply to methods that affect cross-file dependencies.

### @noParallel

Located in `util.py`. Serializes access for XML-modifying methods.

```python
@noParallel
@updateTree('signal')
def modifyXML(self):
    # XML modifications
```

**Important:** Decorator order matters:

```python
# CORRECT order:
@noParallel          # First
@updateTree('signal') # Second

# WRONG - will cause issues:
@updateTree('signal')
@noParallel
```

## Expression Helpers

Located in `expressions.py`.

### createElem

Create XML elements.

```python
from pyfortool.expressions import createElem

elem = createElem('named-E')
elem = createElem('literal-E', text='42')
elem = createElem('n', text='X', tail='\n')
```

### createExpr

Parse FORTRAN statements to XML.

```python
from pyfortool.expressions import createExpr

nodes = createExpr('X = 42')
nodes = createExpr('CALL SUB(A, B)')
nodes = createExpr('IF (A > B) THEN\n  X = 1\nEND IF')
```

### createExprPart

Create single expression nodes.

```python
from pyfortool.expressions import createExprPart

# Variables
createExprPart('X')
# -> <named-E><N><n>X</n></N></named-E>

# Literals
createExprPart('42')
# -> <literal-E><l>42</l></literal-E>

# Structure members
createExprPart('A%B')
# -> <named-E><N><n>A</n></N><R-LT><component-R>%B</component-R></R-LT></named-E>
```

### simplifyExpr

Combine numeric constants.

```python
from pyfortool.expressions import simplifyExpr

simplifyExpr('1+1+I+JI-I')  # -> '2+JI'
simplifyExpr('X+1', add='Y')  # -> 'X+Y+1'
```

## Error Handling

### PYFTError

```python
from pyfortool.util import PYFTError

try:
    scope = pft.getScopeNode('nonexistent:SCOPE')
except PYFTError as e:
    print(f"Error: {e}")
```

### Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `scopePath not found` | Invalid scope path | Check path format |
| `scopePath found several times` | Duplicate scope | Use more specific path |
| `file not found` | Invalid filename | Check file exists |

## See Also

- [Architecture Guide](md__home_sriette_GIT_pyfortool_doc_developer_architecture.html) - Class hierarchy and data flow
- [Module Organization](md__home_sriette_GIT_pyfortool_doc_developer_modules.html) - What each module contains
- [API Reference](../html/index.html) - Method documentation
