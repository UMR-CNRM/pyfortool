"""
PyForTool - Python FORTRAN source-to-source transformation tool.

A Python package built on top of fxtran for manipulating FORTRAN source files.
Provides capabilities for:

- Parsing FORTRAN source code into XML
- Analyzing and querying code structure (scopes, variables, statements)
- Applying transformations (code optimization, GPU preparation, instrumentation)
- Writing transformed code back to FORTRAN

Basic Usage
----------
>>> from pyfortool import PYFT
>>> pft = PYFT('input.F90')
>>> pft.upperCase()
>>> pft.write()

Main Classes
-----------
PYFT : Main class for file-level operations
PYFTscope : Scope-level operations (modules, subroutines, functions)

See Also
--------
fxtran : External FORTRAN parser used by PyForTool
"""

__version__ = "0.2.15"

# NAMESPACE should be defined first
NAMESPACE = 'http://fxtran.net/#syntax'

# Import the public part of the package
from . import util  # pylint: disable=wrong-import-position # noqa: F401 E402
from . import expressions  # pylint: disable=wrong-import-position # noqa: F401 E402
# pylint: disable-next=wrong-import-position
from .pyfortool import PYFT, conservativePYFT  # noqa: F401 E402
from . import tree  # pylint: disable=wrong-import-position # noqa: F401 E402
