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

    Enables support and explicit instantiations for the ``float`` scalar
    data-type in all supported Trilinos packages.

  ``-DTrilinos_ENABLE_COMPLEX=ON``

    Enables support and explicit instantiations for the ``std::complex<T>``
    scalar data-type in all supported Trilinos packages.

  ``-DTrilinos_ENABLE_COMPLEX_FLOAT=ON``

    Enables support and explicit instantiations for the
    ``std::complex<float>`` scalar data-type in all supported Trilinos
    packages.  This is set to ``ON`` by default when
    ``-DTrilinos_ENABLE_FLOAT=ON`` and ``-DTrilinos_ENABLE_COMPLEX=ON`` are
    set.

  ``-DTrilinos_ENABLE_COMPLEX_DOUBLE=ON``

    Enables support and explicit instantiations for the
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

When co-developing TriBTS and Trilinos (after cloning the TriBITS repo
https://github.com/TriBITSPub/TriBITS under the local Trilinos git repo)
configure Trilinos to use that TriBITS implementation using, for example::

   -D Trilinos_TRIBITS_DIR:STRING=TriBITS/tribits

(NOTE: You have to use the data-type ``STRING`` with ``Trilinos_TRIBITS_DIR``
or CMake will automatically assume it is relative to the build dir!)


Configuring with Kokkos and advanced back-ends
----------------------------------------------

Kokkos (https://github.com/kokkos/kokkos) is a C++ implementation of a
cross-platform shared-memory parallel programming model. Many Trilinos packages,
and other stand-alone applications, use it to implement parallel algorithms.

The Kokkos package is enabled with ``-DTrilinos_ENABLE_Kokkos=ON``), then
the native configuration option of the Kokkos package are available such as
``-DKokkos_ENABLE_OPENMP=ON`` or ``-DKokkos_ENABLE_CUDA=ON``. We recommend that
you visit https://kokkos.org/kokkos-core-wiki/get-started/building-from-source.html
and https://kokkos.org/kokkos-core-wiki/get-started/configuration-guide.html
for further information on how to select the various configuration options
for your build with Kokkos.

The following Trilinos variables will pass their value to the equivalent
Kokkos variable

============================    ======================================   ============================
Functionality                   CMake Cache Variable                     Kokkos Variable
============================    ======================================   ============================
Device options:
* Enable Cuda                   ``TPL_ENABLE_CUDA``                      ``Kokkos_ENABLE_CUDA``
* Enable OpenMP                 ``Trilinos_ENABLE_OpenMP``               ``Kokkos_ENABLE_OPENMP``
* Enable Pthread                ``TPL_ENABLE_PThread``                   ``Kokkos_ENABLE_THREADS``
============================    ======================================   ============================

Note: Trilinos always turns on the Kokkos Serial backend even if it was disabled
explicitly using ``-D Kokkos_ENABLE_SERIAL=OFF``.

To see more documentation for each of these options, run a configure with
``-DTrilinos_ENABLE_Kokkos=ON`` and then look in the ``CMakeCache.txt`` file
(as raw text or using the CMake QT GUI or ``ccmake``).


Setting the C++ language standard for Trilinos
----------------------------------------------

Trilinos currently supports building with the C++17 language standard as
supported by a wide range of C++ compilers.  In addition, the library targets
imported from the installed ``<Package>Config.cmake`` files (also pulled in
through ``TrilinosConfig.cmake``) will automatically require downstream CMake
projects turn on C++17 or later standard support in the compiler options
(using the CMake ``INTERFACE_COMPILE_FEATURES`` properties of the Trilinos
library targets).  Building Trilinos with C++14 or lower C++ language
standards is not supported.

However, to try building Trilinos with a higher C++ language standard (with a
supporting compiler), set the CMake cache variable ``CMAKE_CXX_STANDARD`` to
an appropriate value.  For example, to try building Trilinos with C++20 turned
on, configure with::

  -D CMAKE_CXX_STANDARD:STRING=20

As mentioned above, that will also result in all downstream C++ software built
with CMake to be built with C++20 compiler options turned on as well.

However, Trilinos is currently only rigorously tested with C++17 compiler
options so trying to build and use with a higher language standard may not
give satisfactory results.


Addressing problems with large builds of Trilinos
-------------------------------------------------

Trilinos is a large collection of complex software.  Depending on what gets
enabled when building Trlinos, one can experience build and installation
problems due to this large size.

When running into problems like these, the first thing that should be tried is
to **upgrade to and use the newest supported version of CMake!** In some
cases, newer versions of CMake may automatically fix problems with building
and installing Trilinos.  Typically, Trilinos is kept current with new CMake
releases as they come out.

Otherwise, some problems that can arise when and solutions to those problems
are mentioned below.

**Command-line too long errors:**

When turning on some options and enabling some set of package's one may
encounter command-lines that are too long for the OS shell or the tool being
called.  For example, on some systems, enabling CUDA and COMPLEX variable
types (e.g. ``-D TPL_ENABLE_CUDA=ON -D Trilinos_ENABLE_COMPLEX=ON``) can
result in "File 127" errors when trying to create libraries due to large
numbers of ``*.o`` object files getting passed to create some libraries.

Also, on some systems, the list of include directories may become so long that
one gets "Command-line too long" errors during compilation.

These and other cases can be addressed by explicitly enabling built-in CMake
support for ``*.rsp`` resource files as described in the section `Enabling the
usage of resource files to reduce length of build lines`_.

**Large Object file errors:**

Depending on settings and which packages are enabled, some of the ``*.o``
files can become very large, so large that it overwhelms the system tools to
create libraries.  One known case is older versions of the ``ar`` tool used to
create static libraries (i.e. ``-D BUILD_SHARED_LIBS=OFF``) on some systems.
Versions of ``ar`` that come with the BinUtils package **before** version 2.27
may generate "File Truncated" failures when trying to create static libraries
involving these large object files.

The solution to that problem is to use a newer version of BinUtils 2.27+ for
which ``ar`` can handle these large object files to create static libraries.
Just put that newer version of ``ar`` in the default path and CMake will use
it or configure with::

  -D CMAKE_AR=<path-to-updated-binutils>/bin/ar

**Long make logic times:**

On some systems with slower disk operations (e.g. NFS mounted disks), the time
that the ``make`` program with the ``Unix Makefiles`` generator to do
dependency analysis can be excessively long (e.g. cases of more than 2 minutes
to do dependency analysis have been reported to determine if a single target
needs to be rebuilt).  The solution is to switch from the default ``Unix
Makefiles`` generator to the ``Ninja`` generator (see `Enabling support for
Ninja`_).


Enabling and viewing build statistics
-------------------------------------

The Trilinos project has portable built-in support for generating and
reporting build statistics such high-watermark for RAM, wall clock time, file
size, and many other statistics used to build each and every object file,
library, and executable target in the project (and report that information to
CDash).  To enable support for these build statistics, configure with::

  -D Trilinos_ENABLE_BUILD_STATS=ON \

This will do the following:

* Generate wrappers ``build_stats_<op>_wrapper.sh`` for C, C++, and Fortran
  (and for static builds also ``ar``, ``randlib`` and ``ld``) in the build
  tree that will compute statics as a byproduct of every invocation of these
  commands.  (The wrappers create a file ``<output-file>.timing`` for every
  generated object, library and executable ``<output-file>`` file.)

* Define a build target called ``generate-build-stats`` that when run will
  gather up all of the generated build statistics into a single CSV file
  ``build_stats.csv`` in the base build directory.  (This target also runs at
  the end of the ``ALL`` target so a raw ``make`` will automatically create an
  up-to-date ``build_stats.csv`` file.)

* By default, enable the package ``TrilinosBuildStats`` (and when
  ``-DTrilinos_ENABLE_TESTS=ON`` or ``-DTrilinosBuildStats_ENABLE_TESTS=ON``
  are also set) will define the test ``TrilinosBuildStats_Results`` to
  summarize and report the build statistics.  When run, this test calls the
  tools ``gather_build_stats.py`` and ``summarize_build_stats.py`` to gather
  and report summary build stats to STDOUT and will also upload the file
  ``build_stats.csv`` to CDash as using the CTest property ``ATTACHED_FILES``
  when submitting test results to CDash.

The default for the cache variable ``Trilinos_ENABLE_BUILD_STATS`` is
determined as follows:

* If the variable ``Trilinos_ENABLE_BUILD_STATS`` is set in the environment
  (e.g. with ``export Trilinos_ENABLE_BUILD_STATS=ON``), then it will be used
  as the default value.

* Else if the CMake variable ``Trilinos_ENABLE_BUILD_STATS_DEFAULT`` is set in
  a ``*.cmake`` file included using
  ``-DTrilinos_CONFIGURE_OPTIONS_FILE=<config_file>.cmake``, then it will be
  used as the default value.

* Else, the default value is set to ``OFF``.

Otherwise, if ``Trilinos_ENABLE_BUILD_STATS`` is explicitly set in the cache
with ``-DTrilinos_ENABLE_BUILD_STATS=ON|OFF``, then that value will be used.

When the test ``TrilinosBuildStats_Results`` is run, it produces summary
statistics to STDOUT like shown below::

  Full Project: sum(max_resident_size_size_mb) = ??? (??? entries)
  Full Project: max(max_resident_size_size_mb) = ??? (<file-name>)
  Full Project: max(elapsed_real_time_sec) = ??? (<file-name>)
  Full Project: sum(elapsed_real_time_sec) = ??? (??? entries)
  Full Project: sum(file_size_mb) = ??? (??? entries)
  Full Project: max(file_size_mb) = ??? (<file-name>)

  <package1>: sum(max_resident_size_mb) = ??? (??? entries)
  <package1>: max(max_resident_size_mb) = ??? (<file-name>)
  <package1>: max(elapsed_real_time_sec) = ??? (<file-name>)
  <package1>: sum(elapsed_real_time_sec) = ??? (??? entries)
  <package1>: sum(file_size_mb) = ??? (??? entries)
  <package1>: max(file_size_mb) = ??? (<file-name>)

  ...

  <packagen>: sum(max_resident_size_mb) = ??? (??? entries)
  <packagen>: max(max_resident_size_mb) = ??? (<file-name>)
  <packagen>: max(elapsed_real_time_sec) = ??? (<file-name>)
  <packagen>: sum(elapsed_real_time_sec) = ??? (??? entries)
  <packagen>: sum(file_size_mb) = ??? (??? entries)
  <packagen>: max(file_size_mb) = ??? (<file-name>)

where:

* ``max_resident_size_size_mb`` is the high watermark for RAM usage to build a
  given target measured in MB.
* ``elapsed_real_time_sec`` is the wall clock time used to build a given
  target measured in seconds.
* ``file_size_mb`` is the file size of a given build target (i.e. object file,
  library, or executable) measured in MB.
* ``Full Project`` are the stats for all of the enabled Trilinos packages.
* ``<packagei>`` are the build stats for the ``<packagei>`` subdirectory under
  the base directories ``commonTools`` and ``packages``.  (These map to Trilinos
  packages is most cases.)

This output format makes it easy to query and view these statistics directly
on CDash using the "Test Output" filter on the ``cdash/queryTests.php`` page.
(This allows viewing and comparing these statistics across many different
compilers, platforms, and build configurations and even across the same builds
over days, weeks, and months.)

The generated ``build_stats.csv`` file contains many other types of useful
build stats as well but the above three are some of the more significant build
statistics.

To avoid situations where a full rebuild does not occur (e.g. any build target
fails) and an old obsolete ``build_stats.csv`` file is hanging around, one can
cause that file to get deleted on every (re)configure by setting::

  -D Trilinos_REMOVE_BUILD_STATS_ON_CONFIGURE=ON

This will remove the file ``build_stats.csv`` very early in the configure
process and therefore will usually remove the file even of later configure
operations fail.

Finally, to make rebuilds more robust and to restrict build stats to only new
targets getting (re)built after an initial configure, then configure with::

  -D Trilinos_REMOVE_BUILD_STATS_TIMING_FILES_ON_FRESH_CONFIGURE=ON

This will remove **all** of the ``*.timing`` files under the base build
directory during a fresh configure (i.e. where the ``CMakeCache.txt`` file
does not exist).  But this will not remove ``*.timing`` files on reconfigures
(i..e where a ``CMakeCache.txt`` file is preserved).  Timing stats for targets
that are already built and don't need to be rebuilt after the last fresh
configure will not get reported.  (But this can be useful for CI builds where
one only wants to see build stats for the files updated in the last PR
iteration.

NOTES:

* The underlying compilers must already be specified in the cache variables
  ``CMAKE_C_COMPILER``, ``CMAKE_CXX_COMPILER``, and ``CMAKE_Fortran_COMPILER``
  and not left up to CMake to determine.  The best way to do that is, for
  example ``-DCMAKE_C_COMPILER=$(which mpicc)`` on the ``cmake`` command-line.

* The tool ``gather_build_stats.py`` is very robust and will discard data from
  any invalid or corrupted ``*.timing`` files and can deal with ``*.timing``
  files with different sets and ordering of the data fields from different
  versions of the build stats wrapper tool.  (Therefore, one can keep
  rebuilding in an existing build directory with old ``*.timing`` files
  hanging around and never have to worry about being able to create an updated
  ``build_stats.csv`` file.)

* The installed ``TrilinosConfig.cmake`` and ``<Package>Config.cmake`` files
  list the original underlying C, C++, and Fortran compilers, **not** the
  build stats compiler wrappers.

* The ``generate-build-stats`` target has dependencies on every object,
  library, and executable build target in the project so it will always only
  run after all of those targets are up to date.

* After uploading the test results to CDash, the file ``build_stats.csv`` can
  be downloaded off CDash from the ``TrilinosBuildStats_Results`` test results
  details page.  (The file is downloaded as a compressed
  ``build_stats.csv.tgz`` file which will then need to be uncompressed using
  ``tar -xzvf build_stats.csv.tgz`` before viewing.)

* Any ``build_stats.csv`` file can be viewed and queried by uploading it to
  the site ``https://jjellio.github.io/build_stats/index.html``.
