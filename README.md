# Python-fortran-tool

This set of tools is primarily designed to be applied on the [PHYEX](https://github.com/UMR-CNRM/PHYEX)
repository. But they are normally generic enough to be used in other contexts.

The tools can be used with the command line (the pyfortool.py and pyfortool\_parallel.py utilities)
or as a python package. The documentation is in the doc directory and some examples are in the
examples directory. 

Prerequisites:
  - python > 3.8 (but only tested with version 3.10)

Other prerequisites are linked to the pyfxtran installation:
  - a C compiler (tested with cc 11.4.0)
  - make and git to compile fxtran (on the first use only)

There are two installation methods, depending on your use.
On some systems, you may need to create a python virtual environment to use these methods.

Method 1 (suitable for developpers and end users):
  - open a terminal on a system satisfying the prerequisites and enter the following commands
  - if you don't have a github ssh key or don't know what it is:
    > git clone https://github.com/UMR-CNRM/pyfortool.git
  - if you have a github ssh key:
    > git clone git@github.com:UMR-CNRM/pyfortool.git
  - > cd pyfortool; pip install -e .

Method 2 (suitable for users who do not develop PyForTool):
  - open a terminal on a system satisfying the prerequisites and enter the following commands
  - pip install pyfortool
