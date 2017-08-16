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

  -D <Project>_ENABLE_THREAD_SAFE:BOOL=ON
  
This will set the default value ``Teuchos_ENABLE_THREAD_SAFE=ON`` which makes
the Teuchos Memory Management classes (Teuchos::RCP, Teuchos::Ptr,
Teuchos::Array, Teuchos::ArrayView, and Teuchos::ArrayRCP) thread-safe.  See
documentation for other Trilinos packages for what parts become thread safe
when setting this option.


Enabling/disabling time monitors
--------------------------------

I order to enable instrumentation of select code to generate timing statistics, set::

 -D <Project>_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON

This will enable Teuchos time monitors by default in all Trilinos packages
that support them.  To print the timers at the end of the program, call
``Teuchos::TimeMonitor::summarize()``.

In order do co-development of TriBTS and Trilinos (see http://http://trac.trilinos.org/wiki/TriBITSTrilinosDev), set::

   -D <Project>_TRIBITS_DIR:STRING=TriBITS \
   -D <Project>_TRIBITS_PACKAGE_USE_TRIBITS_DIR=TRUE

(NOTE: You have to use the data-type ``STRING`` with ``Trilinos_TRIBITS_DIR``
or CMake will automatically assume it is relative to the build dir!)
