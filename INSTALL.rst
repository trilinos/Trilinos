==================================================================
Quick configure, build and install hints for Trilinos
==================================================================

:Contact: trilinos-framework@software.sandia.gov

This document is intended to be a very concise set of examples of how to
configure, build, install, and build downstream code against Trilinos. The
intended audience is those who need a quick refresher on Trilinos CMake build
system or those wanting a quick install without worrying about all the
features and options that are available.  For a more in-depth document on what
features and options are available and how to use them, see the document:

  https://trilinos.org/docs/files/TrilinosBuildReference.html

.. sectnum::

.. contents::


Requirements
============

* CMake 3.23.0 or newer
* C and C++ compiler
* Optionally a Fortran compiler
* Optionally an installation of MPI


Instructions
============

Following are a few examples of simple configurations for Trilinos. Anything in
<> should be replaced with the appropriate path or value and excluding the <>.

It is recommended that you put your configure options in a script (e..g
``do-configure``) so you can repeat the configure if necessary.

Note: all examples assume a unix like command line and the CMake Makefile
Generator.


Simple MPI instructions (enables most packages)
-----------------------------------------------

::

  cmake \
  -DTPL_ENABLE_MPI=ON \
  -DMPI_BASE_DIR=<path to mpi install> \
  -DTrilinos_ENABLE_ALL_PACKAGES=ON \
  -DCMAKE_INSTALL_PREFIX=<path to install Trilinos into> \
  <path to Trilinos source>
  
  make -j<n> install

NOTE: Enabling all packages will trigger the enable of several third-party
libraries (TPLs).  If CMake can't find these TPLs on your system, then you
will need to either point to them or disable them as specified in the CMake
output.


Simple non-MPI instructions (enables most packages)
---------------------------------------------------

::

  cmake \
  -DCMAKE_C_COMPILER=<path to C compiler> \
  -DCMAKE_CXX_COMPILER=<path to C++ compiler> \
  -DCMAKE_Fortran_COMPILER=<path to Fortran compiler> \
  -DTrilinos_ENABLE_ALL_PACKAGES=ON \
  -DCMAKE_INSTALL_PREFIX=<path to install Trilinos into> \
  <path to Trilinos source>
  
  make -j<n> install


Intermediate MPI instructions (enables a few packages)
------------------------------------------------------

::

  cmake \
  -DTPL_ENABLE_MPI=ON \
  -DMPI_BASE_DIR=<path to mpi install> \
  -DTrilinos_ENABLE_Epetra=ON \
  -DTrilinos_ENABLE_AztecOO=ON \
  -DTrilinos_ENABLE_Ifpack=ON \
  -DCMAKE_INSTALL_PREFIX=<path to install Trilinos into> \
  <path to Trilinos source>
  
  make -j<n> install


Intermediate MPI instructions (enables a few packages, explicit compilers)
-------------------------------------------------------------------------

::

  cmake \
  -DTPL_ENABLE_MPI=ON \
  -DCMAKE_C_COMPILER=<path to MPI C compiler wrapper> \
  -DCMAKE_CXX_COMPILER=<path to MPI C++ compiler wrapper> \
  -DCMAKE_Fortran_COMPILER=<path to MPI Fortran compiler wrapper> \
  -DTrilinos_ENABLE_Epetra=ON \
  -DTrilinos_ENABLE_AztecOO=ON \
  -DTrilinos_ENABLE_Ifpack=ON \
  -DCMAKE_INSTALL_PREFIX=<path to install Trilinos into> \
  <path to Trilinos source>
  
  make -j<n> install


Useful Options
==============

To generate Ninja build files (Ninja 1.10+) instead of Makefiles use::

  -GNinja

To use shared libraries (much smaller executables and faster linking) use::

  -DBUILD_SHARED_LIBS=ON

To enable support for the ``float`` scalar type use::

  -DTrilinos_ENABLE_FLOAT=ON

To enable support for ``std::complex<T>`` scalar types use::

  -DTrilinos_ENABLE_COMPLEX=ON

To disable Fortran use the following::

  -DTrilinos_ENABLE_Fortran=OFF

To enable a package::

  -DTrilinos_ENABLE_<package name>=ON

To get the list of packages that can be enabled, run::

  cmake <path to Trilinos source> 2>&1 \
    | grep "Final set of non-enabled SE packages"

To enable tests::

  -DTrilinos_ENABLE_TESTS=ON


Building against installed Trilinos
===================================

For information on how to build against an installation of Trilinos, see
`demos/simpleBuildAgainstTrilinos`_

.. _demos/simpleBuildAgainstTrilinos: demos/simpleBuildAgainstTrilinos/README.md


Support for Different Compilers/MPIs
====================================

Trilinos tests with all Compiler/MPI combinations listed at:

  https://github.com/trilinos/Trilinos/wiki/Pull-Request-Testing-Interface

Compilers/MPIs that are not part of our automated process will not receive support.
We are happy to accept Pull Requests enhancing support for other compilers/MPIs as
needed by our customers, with the understanding that said configurations cannot be
guaranteed to work.
