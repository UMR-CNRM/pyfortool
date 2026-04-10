# Module Organization

This guide explains what each module in PyForTool contains and when to use it.

## Module Overview

```
pyfortool/
├── __init__.py          # Package initialization, exports PYFT
├── pyfortool.py         # PYFT class (file I/O)
├── scope.py             # PYFTscope, ElementView classes
├── variables.py        # VarList, Variables classes
├── statements.py        # Statements class
├── cosmetics.py         # Cosmetics class
├── applications.py     # Applications class
├── cpp.py              # Cpp class
├── openacc.py          # Openacc class
├── tree.py             # Tree class
├── util.py             # Utility functions
├── expressions.py      # Expression helpers
└── scripting.py        # CLI tools implementation
```

## pyfortool.py - File Operations

**Main class:** `PYFT`

Entry point for file-level operations.

### Key Methods

| Method | Description |
|--------|-------------|
| `__init__(filename, output)` | Open FORTRAN file |
| `write()` | Write transformed code |
| `writeXML(filename)` | Write XML representation |
| `renameUpper()` / `renameLower()` | Change file extension case |
| `setParallel(tree)` | Enable parallel processing |
| `conservativePYFT()` | Create conservative parser instance |
| `generateEmptyPYFT()` | Create new file from scratch |

### Usage

```python
from pyfortool import PYFT

# Basic usage
pft = PYFT('input.F90')
pft.upperCase()
pft.write()

# With output file
pft = PYFT('input.F90', output='output.F90')
pft.write()

# Context manager
with PYFT('input.F90') as pft:
    pft.removeComments()
    pft.write()
```

## scope.py - Scope Operations

**Main classes:** `PYFTscope`, `ElementView`

Provides scope-level navigation and XML tree access.

### PYFTscope

| Method | Description |
|--------|-------------|
| `getScopes()` | Get all scopes |
| `getScopeNode(path)` | Get scope by path |
| `getScopePath(item)` | Get path for element |
| `getParent(item)` | Get parent element |
| `getFileName()` | Get source filename |
| `empty()` | Remove all statements |

### ElementView (XML Access)

| Method | Description |
|--------|-------------|
| `find()` | Find element (xpath) |
| `findall()` | Find all elements |
| `iter()` | Iterate elements |
| `iterfind()` | Iterate matching elements |
| `insert()` | Insert element |
| `remove()` | Remove element |
| `insertStatement()` | Insert FORTRAN statement |

### Usage

```python
# Navigate scopes
scopes = pft.getScopes()
sub = pft.getScopeNode('module:MOD/sub:SUB')

# Query XML
node = pft.find('.//{*}call-stmt')
nodes = pft.findall('.//{*}named-E')

# Get parent
parent = sub.getParent(node)

# Insert statement
from pyfortool.expressions import createExpr
stmt = createExpr('X = 42')[0]
sub.insertStatement(stmt, first=True)
```

## variables.py - Variable Management

**Main classes:** `VarList`, `Variables`

### VarList

Manages variable declarations for a scope.

| Method | Description |
|--------|-------------|
| `findVar(name)` | Find variable |
| `restrict(path, excludeContains)` | Filter by scope |
| `showVarList()` | Display all variables |

### Variables (Mixin)

| Method | Description |
|--------|-------------|
| `addVar(varList)` | Add variable declaration |
| `removeVar(varList)` | Remove variable |
| `addModuleVar(moduleVarList)` | Add USE statement |
| `attachArraySpecToEntity()` | Move DIMENSION to entities |
| `checkImplicitNone()` | Check for IMPLICIT NONE |
| `checkIntent()` | Check INTENT attributes |
| `checkUnusedLocalVar()` | Find unused variables |
| `removeUnusedLocalVar()` | Remove unused variables |
| `addExplicitArrayBounds()` | Expand A(:) to A(lbound:ubound) |
| `addArrayParentheses()` | Add A to A(:) |
| `modifyAutomaticArrays()` | Transform automatic arrays |

### Usage

```python
# Query variables
vl = pft.varList
var = vl.findVar('X')

# Add variable
pft.addVar([('module:MOD/sub:SUB', 'NEW_VAR', 
             'INTEGER :: NEW_VAR', None)])

# Remove variable
pft.removeVar([('module:MOD/sub:SUB', 'OLD_VAR')])

# Add USE statement
pft.addModuleVar([('module:MOD/sub:SUB', 'MODD_XX', 'VAR')])

# Transform automatic arrays
pft.modifyAutomaticArrays(
    declTemplate="{type}, DIMENSION({doubledotshape}), ALLOCATABLE :: {name}",
    startTemplate="ALLOCATE({name}({shape}))",
    endTemplate="DEALLOCATE({name})"
)
```

## statements.py - Statement Manipulation

**Main class:** `Statements`

### Methods

| Method | Description |
|--------|-------------|
| `removeCall(name)` | Remove CALL statements |
| `removePrints()` | Remove PRINT statements |
| `removeArraySyntax()` | Convert A(:) = B(:) to DO loops |
| `inlineContainedSubroutines()` | Inline helper subroutines |
| `setFalseIfStmt(flags)` | Set flags to .FALSE. |

### Usage

```python
# Remove calls
pft.removeCall('FOO')  # All CALL FOO
pft.removeCall('FOO', simplify=True)  # + remove unused vars

# Array syntax expansion
pft.removeArraySyntax()  # Standard DO loops
pft.removeArraySyntax(concurrent=True)  # DO CONCURRENT
pft.removeArraySyntax(everywhere=False)  # Only marked sections

# Inline subroutines
pft.inlineContainedSubroutines()
pft.inlineContainedSubroutines(simplify=True)

# Disable conditional code
pft.setFalseIfStmt('LDEBUG')  # IF (LDEBUG) becomes .FALSE.
pft.setFalseIfStmt(['LDEBUG', 'LVERBOSE'], simplify=True)
```

## cosmetics.py - Code Formatting

**Main class:** `Cosmetics`

### Methods

| Method | Description |
|--------|-------------|
| `upperCase()` | UPPER CASE keywords |
| `lowerCase()` | lower case keywords |
| `indent()` | Fix indentation |
| `removeComments()` | Remove comments |
| `removeEmptyLines()` | Remove blank lines |
| `changeIfStatementsInIfConstructs()` | IF x → IF x THEN |
| `updateSpaces()` | Normalize whitespace |
| `updateContinuation()` | Handle line continuations |
| `prettify()` | Full formatting |
| `minify()` | Compact formatting |

### Usage

```python
# Case
pft.upperCase()
pft.lowerCase()

# Indentation
pft.indent()  # Default (2 spaces)
pft.indent(indentProgramunit=0, indentBranch=4)

# Comments
pft.removeComments()  # Keep directives
pft.removeComments(exclDirectives=[])  # Remove all
pft.removeComments(pattern=re.compile(r'!.*TODO'))

# Prettify (all formatting)
pft.prettify()
```

## applications.py - High-Level Transforms

**Main class:** `Applications`

Model-specific and application-level transformations.

### Methods

| Method | Description |
|--------|-------------|
| `addDrHook()` / `deleteDrHook()` | Profiling instrumentation |
| `deleteBudgetDDH()` | Remove budget diagnostics |
| `addStack(model, stopScopes)` | GPU stack allocation |
| `convertTypesInCompute()` | Inline structure members |
| `deleteNonColumnCallsPHYEX()` | Remove multi-column dependencies |
| `removeIJDim()` | Remove I,J dimensions |
| `expandAllArrays()` | Expand all array syntax |
| `inlineContainedSubroutinesPHYEX()` | PHYEX-style inlining |
| `buildModi()` | Generate *_modi.f90 files |
| `splitModuleRoutineFile()` | Split file by scopes |

### Usage

```python
# Profiling
pft.addDrHook()    # Add timing
pft.deleteDrHook() # Remove timing

# GPU preparation
pft.addStack('AROME', ['module:MOD/sub:DRIVER'])
pft.addStack('MESONH', ['module:MOD/sub:DRIVER'])

# Optimization
pft.convertTypesInCompute()  # STR%VAR → LOCAL_VAR

# Split file
pft.splitModuleRoutineFile()  # Creates one file per module/subroutine
```

## cpp.py - Preprocessor

**Main class:** `Cpp`

### Methods

| Method | Description |
|--------|-------------|
| `applyCPPifdef(keys)` | Evaluate #ifdef blocks |

### Usage

```python
# Evaluate specific keys
pft.applyCPPifdef(['KEY1', '%KEY2'])
# - KEY1 is defined (True)
# - KEY2 is undefined (False, via % prefix)
```

## openacc.py - GPU Directives

**Main class:** `Openacc`

### Methods

| Method | Description |
|--------|-------------|
| `addACCData()` | Add !$acc data directives |
| `addACCRoutineSeq(stopScopes)` | Add !$acc routine seq |
| `removeACC()` | Remove !$acc directives |
| `craybyPassDOCONCURRENT()` | Handle CRAY compiler issues |
| `allocatetoHIP()` | Convert to AMD HIP |
| `removebyPassDOCONCURRENT()` | Remove MNH bypass macros |

### Usage

```python
# Add directives
pft.addACCData()  # !$acc data present() for INTENT arrays

# Remove directives
pft.removeACC()

# AMD GPU preparation
pft.allocatetoHIP()  # ALLOCATE → CALL MNH_HIPALLOCATE
```

## tree.py - Cross-File Analysis

**Main class:** `Tree`

### Methods

| Method | Description |
|--------|-------------|
| `getFiles()` | List all source files |
| `needsFile()` | Get compilation dependencies |
| `usedByFile()` | Get files depending on file |
| `whoCall()` | Find callers of routine |
| `scopeCallees()` | Find routines called by scope |
| `plotCompilTreeFromFile()` | Generate compilation graph |
| `plotExecTreeFromFile()` | Generate call graph |

### Usage

```python
from pyfortool.tree import Tree

tree = Tree(['/path/to/src'], descTreeFile='tree.json')

# Query dependencies
deps = tree.needsFile('file.F90')
users = tree.usedByFile('file.F90')

# Visualization
tree.plotCompilTreeFromFile('main.F90', 'deps.dot')
tree.plotExecTreeFromFile('main.F90', 'calls.dot')
```

## util.py - Utilities

Helper functions and decorators.

### Functions

| Function | Description |
|----------|-------------|
| `fortran2xml(source)` | Parse FORTRAN to XML |
| `tostring(xml)` | XML to string |
| `tofortran(xml)` | XML to FORTRAN |
| `tag(elem)` | Get tag without namespace |
| `n2name(nodeN)` | Extract name from N element |
| `alltext(elem)` | Get all text content |
| `isint()` / `isfloat()` | Type checking |
| `nonCode()` | Check if comment/whitespace |
| `isExecutable()` | Check if executable |
| `setVerbosity(level)` | Set logging level |

### Decorators

| Decorator | Purpose |
|-----------|---------|
| `@debugDecor` | Tracing and profiling |
| `@noParallel` | Prevent parallel execution |

## expressions.py - Expression Helpers

Low-level XML construction.

### Functions

| Function | Description |
|----------|-------------|
| `createElem()` | Create XML element |
| `createExpr()` | Parse statement to XML |
| `createExprPart()` | Create expression node |
| `simplifyExpr()` | Simplify numeric expressions |
| `createArrayBounds()` | Create array bounds |

### Usage

```python
from pyfortool.expressions import createExpr, createElem

# Create statement
nodes = createExpr('X = 42')

# Insert into tree
pft.insertStatement(nodes[0], first=True)
```

## scripting.py - CLI Implementation

Contains functions for the command-line tools.

### Functions

| Function | Description |
|---------|-------------|
| `main()` | Single-file transformation |
| `mainParallel()` | Parallel transformation |
| `getArgs()` | Parse command-line arguments |
| `applyTransfo()` | Apply transformations |
| `updateParser*()` | Configure argument parser |

**Not typically used directly** - use via CLI tools:
```bash
pyfortool input.F90 --upperCase --indent
pyfortool_parallel --tree /path/to/src ...
```

## Choosing the Right Module

| Task | Module |
|------|--------|
| Open file, read, write | `pyfortool.py` (PYFT) |
| Navigate scopes | `scope.py` |
| Modify variables | `variables.py` |
| Remove/add statements | `statements.py` |
| Format code | `cosmetics.py` |
| Model-specific transforms | `applications.py` |
| Preprocessor handling | `cpp.py` |
| GPU directives | `openacc.py` |
| Cross-file analysis | `tree.py` |
| Build XML nodes | `expressions.py` |
| Debug/tracing | `util.py` |

## See Also

- [Architecture Guide](md__home_sriette_GIT_pyfortool_doc_developer_architecture.html) - How modules fit together
- [Core Concepts](md__home_sriette_GIT_pyfortool_doc_developer_concepts.html) - Detailed concept explanations
- [API Reference](../html/index.html) - Method documentation
