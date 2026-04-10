# PyForTool API Reference

This section contains the auto-generated API reference from source code docstrings.

For general documentation, see:

- [User's Guide](md__home_sriette_GIT_pyfortool_doc_Documentation.html) - End-user documentation
- [Architecture Guide](md__home_sriette_GIT_pyfortool_doc_developer_architecture.html) - Class hierarchy and patterns
- [Core Concepts](md__home_sriette_GIT_pyfortool_doc_developer_concepts.html) - Scope paths, VarList, decorators
- [Module Organization](md__home_sriette_GIT_pyfortool_doc_developer_modules.html) - What each module contains

## Main Classes

| Class | Description |
|-------|-------------|
| [PYFT](@ref pyfortool.pyfortool.PYFT) | Main class for file operations |
| [PYFTscope](@ref pyfortool.scope.PYFTscope) | Scope-level operations |
| [VarList](@ref pyfortool.variables.VarList) | Variable management |
| [Tree](@ref pyfortool.tree.Tree) | Cross-file dependency tracking |

## Module Reference

| Module | Description |
|--------|-------------|
| [pyfortool.variables](@ref pyfortool.variables) | Variable declaration management |
| [pyfortool.statements](@ref pyfortool.statements) | Statement manipulation |
| [pyfortool.cosmetics](@ref pyfortool.cosmetics) | Code formatting |
| [pyfortool.applications](@ref pyfortool.applications) | High-level transformations |
| [pyfortool.cpp](@ref pyfortool.cpp) | Preprocessor directives |
| [pyfortool.openacc](@ref pyfortool.openacc) | OpenACC directives |
| [pyfortool.tree](@ref pyfortool.tree) | Dependency analysis |
| [pyfortool.expressions](@ref pyfortool.expressions) | Expression helpers |
| [pyfortool.util](@ref pyfortool.util) | Utilities and decorators |

## Quick Examples

### Basic File Operations

```python
from pyfortool import PYFT

# Read, transform, write
pft = PYFT('input.F90')
pft.upperCase()
pft.write()
```

### Scope Navigation

```python
# Get specific scope
sub = pft.getScopeNode('module:MOD/sub:SUB')

# Find variables
var = sub.varList.findVar('X')
```

### Code Transformations

```python
# Remove calls and simplify
pft.removeCall('FOO', simplify=True)

# Expand array syntax
pft.removeArraySyntax(concurrent=True)

# Add profiling
pft.addDrHook()
```
