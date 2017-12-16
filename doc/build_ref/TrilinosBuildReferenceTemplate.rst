============================================================
Trilinos Configure, Build, Test, and Install Reference Guide
============================================================

:Author: Roscoe A. Bartlett (bartlettra@orn.gov)
:Contact: trilinos-framework@software.sandia.gov

:Abstract: This document contains reference information on how to configure, build, test, and install Trilinos using the TriBITS CMake build system.  The primary audience are users of Trilinos that need to configure and build the software.  The secondary audience are actual developers of Trilinos.

.. sectnum::

.. contents::

Introduction
============

Trilinos contains a large number of packages that can be enabled and there is a fairly complex dependency tree of required and optional package enables.  The following sections contain fairly generic information on how to configure, build, test, and install Trilinos that addresses a wide range of issues.

This is not the first document that a user should read when trying to set up to install Trilinos.  For that, see the INSTALL.* file in the base Trilinos source directory.  There is a lot of information and activities mentioned in this reference that most users (and even some Trilinos developers) will never need to know about.

Also, this particular reference has no information at all on what is actually in Trilinos.  For that, go to:

  http://trilinos.org

to get started.

Trilinos-specific options
=========================

Below, configure options specific to Trilinos are given.  The later sections
give more generic options that are the same for all TriBITS projects.


Enabling float and complex Scalar types
----------------------------------------

Many of the packages in Trilinos are implemented using C++ templates and
therefore support a variety of data-types.  In addition to the default scalar
type ``double``, many of the data-structures and solvers are tested with the
types ``float``, ``std::complex<float>``, and ``std::complex<double>``.  In
addition, these packages support the explicit template instantiation
(i.e. ``-DTrilinos_EXPLICIT_TEMPLATE_INSTANTIATION=ON``) of these types as
well.  However, support and explicit instantiations for these types are off by
default since most users don't need these extra types and enabling them
greatly increases the compilation time for Trilinos libraries and tests and
can consume a great deal more disk space.  But support for these types in the
various Trilinos packages can be enabled using the following options:

  ``-DTrilinos_ENABLE_FLOAT=ON``

    Enables suppport and explicit instantiations for the ``float`` scalar
    data-type in all supported Trilinos packages.

  ``-DTrilinos_ENABLE_COMPLEX=ON``

    Enables suppport and explicit instantiations for the ``std::complex<T>``
    scalar data-type in all supported Trilinos packages.

  ``-DTrilinos_ENABLE_COMPLEX_FLOAT=ON``

    Enables suppport and explicit instantiations for the
    ``std::complex<float>`` scalar data-type in all supported Trilinos
    packages.  This is set to ``ON`` by default when
    ``-DTrilinos_ENABLE_FLOAT=ON`` and ``-DTrilinos_ENABLE_COMPLEX=ON`` are
    set.

  ``-DTrilinos_ENABLE_COMPLEX_DOUBLE=ON``

    Enables suppport and explicit instantiations for the
    ``std::complex<double>`` scalar data-type in all supported Trilinos
    packages.  This is set to ``ON`` by default when
    ``-DTrilinos_ENABLE_COMPLEX=ON`` is set.


Enabling/disabling thread safety
--------------------------------

By default, many Trilinos classes are not thread-safe.  However, some of these
classes can be made thread safe by configuring with::

  -D Trilinos_ENABLE_THREAD_SAFE:BOOL=ON
  
This will set the default value ``Teuchos_ENABLE_THREAD_SAFE=ON`` which makes
the Teuchos Memory Management classes (Teuchos::RCP, Teuchos::Ptr,
Teuchos::Array, Teuchos::ArrayView, and Teuchos::ArrayRCP) thread-safe.  See
documentation for other Trilinos packages for what parts become thread safe
when setting this option.


Enabling/disabling time monitors
--------------------------------

I order to enable instrumentation of select code to generate timing statistics, set::

 -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON

This will enable Teuchos time monitors by default in all Trilinos packages
that support them.  To print the timers at the end of the program, call
``Teuchos::TimeMonitor::summarize()``.

Select different TriBITS implementation
----------------------------------------

In order do co-development of TriBTS and Trilinos (see
http://http://trac.trilinos.org/wiki/TriBITSTrilinosDev), set, for example::

   -D Trilinos_TRIBITS_DIR:STRING=TriBITS/tribits

(NOTE: You have to use the data-type ``STRING`` with ``Trilinos_TRIBITS_DIR``
or CMake will automatically assume it is relative to the build dir!)


Configuring with Kokkos and advanced back-ends
----------------------------------------------

Kokkos (https://github.com/kokkos/kokkos) is a C++ implementation of a
cross-platform shared-memory parallel programming model. Many Trilinos packages,
and other stand-alone applications, use it to implement parallel algorithms.

If the Kokkos package is enabled (e.g. ``-DTrilinos_ENABLE_Kokkos=ON``), then
the following CMake cache variables can be used to get the included Kokkos
configuration system to select compiler and other build related flags for the
target machine.  These build-related flags are selected to create correct and
perforamnt code and for C++ software that uses Kokkos.

============================    ======================================
Functionality                   CMake Cache Varaible
============================    ======================================
Specify architecture            ``KOKKOS_HOST_ARCH``
Debug builds                    ``KOKKOS_DEBUG``
Device options:
* Enable Cuda                   ``TPL_ENABLE_CUDA``
* Enable OpenMP                 ``Trilinos_ENABLE_OpenMP``
* Enable Pthread                ``TPL_ENABLE_PThread``
* Specify Serial                ``TPL_ENABLE_MPI=FALSE``
Advanced options:
* Enable compiler warnings      ``KOKKOS_ENABLE_COMPILER_WARNINGS``
* Aggressive Vectorization      ``KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION``
* Profiling                     ``KOKKOS_ENABLE_PROFILING``
* Enable profile load print     ``KOKKOS_ENABLE_PROFILE_LOAD_PRINT``
* Enable dualview modify chk    ``KOKKOS_ENABLE_DUALVIEW_MODIFY_CHECK``
Kokkos TPLs:                 
* Use hwloc library             ``TPL_ENABLE_HWLOC``
* Use memkind library           ``KOKKOS_ENABLE_MEMKIND``
* Use librt                     ``KOKKOS_ENABLE_LIBRT``
CUDA Options:                
* Enable CUDA LDG               ``KOKKOS_ENABLE_CUDA_LDG_INTRINSIC`` (global mem load)
* Enable CUDA UVM               ``KOKKOS_ENABLE_CUDA_UVM`` (unified virtual mem)
* Enable CUDA RDC               ``KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE``
* Enable CUDA LAMBDA            ``KOKKOS_ENABLE_CUDA_LAMBDA``
============================    ======================================

If the cache var ``KOKKOS_HOST_ARCH`` is not set (or is set to ``None``) then
the Kokkos settings are not used and the default Trilinos CMake configuration
is used as described below.

If ``KOKKOS_HOST_ARCH != None`` is set, then the correct compiler flags for
C++11 and OpenMP are selected by the Kokkos system and the values of the cache
vars ``Trilinos_CXX11_FLAGS`` and ``OpenMP_CXX_FLAGS`` set by the user will be
ignored.

To see more documentation for each of these options, run a configure with
``-DTrilinos_ENABLE_Kokkos=ON`` and then look in the CMakeCache.txt file (as
raw text or using the CMake QT GUI or ``ccmake``).
