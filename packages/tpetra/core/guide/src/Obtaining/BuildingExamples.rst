
.. _building_examples:

Building and Running Tpetra Examples
####################################

All of the examples in this guide can be built and run.  The source code for all
examples is contained in this documentation's ``Source/Examples/SourceCode`` directory.  To
build the examples, Tpetra must first be built as outlined in
:ref:`building_tpetra`.  It is recommended that an "installed" version of Tpetra be
built by adding the following option to the Trilinos configure script:

.. code-block:: sh

   -D CMAKE_INSTALL_PREFIX:PATH="/path/to/TrilinosBuild"

Making Examples
===============

The ``Make.py`` script in ``Source/Examples/SourceCode`` will build each example in the
correct environment.  ``Make.py`` takes a single argument - the name of an
example.  The ``TRILINOS_INSTALL`` environment variable must exist in the
environment and point to the directory to which Trilinos was installed.  For
example, to build the power method example, do:

.. code-block:: sh

   TRILINOS_INSTALL=/path/to/TrilinosInstall ./Make.py power_method_1.cpp

When completed, the ``power_method_1.exe`` executable will have been built.

Running Examples
================

Execute examples at the command line:

.. code-block:: sh

   mpirun -np <N> <example_name>.exe

Or, for non-MPI examples:

.. code-block:: sh

   ./<example_name>.exe
