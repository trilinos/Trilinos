TriBITS/tribits/ Directory Contents
+++++++++++++++++++++++++++++++++++

This directory contains the implementation for the varous parts of TriBITS
that are used by TriBITS projects to implement TriBITS functionality.  It also
contains basic documentation in the subdirectory doc/ that is very close to
the TriBITS implementation.  Files and directories from here are what get
installed on the system or will be snapshotted into
``<projectDir>/cmake/tribits/``.  Each TriBITS Project decides what parts of
TriBITS it wants to install or shapshot using the script
``tribits/snapshot_tribits.py`` (which takes arguments for what dirs to
snapshot). This directory contains no tests at all. All of the tests for
TriBITS are in the ``test/`` directory in the parent TriBITS repository.

The breakdown of the contents of this directory are described below:

**TriBITS.cmake**: The one file that needs to be included in order to use
TriBITS in a CMake project. This one file insulates clients from future
TriBITS refactorings of TriBITS.

**Version.cmake**: Version of TriBITS.  This gets included by TriBITS.cmake

.. _TriBITS Core:

**core/**: Core TriBITS package-based architecture for CMake projects. This
only depends on raw CMake and contains just the minimal support for building,
testing, installing, and deployment.  Only depends on CMake and nothing else.

**python_utils/**: Some basic Python utilities that are not specific to
TriBITS but are used in TriBITS CI and testing support software.  There are
some very useful python scripts here like ``gitdist`` and ``snapshot-dir.py``.

**ci_support/**: Support code for pre-push continuous integration testing.
This contains the ``checkin-test.py`` script and its supporting Python
modules.

**ctest_driver/**: Support for package-by-package testing driven by CTest
submitting to CDash (to CDash project ``<Project>``).  This contains the file
``TribitsCTestDriveCore.cmake`` and some supporting modules.

**dashboard_driver/**: TriBITS Dashboard Driver system which uses CTest to
drive individual project builds in parallel and submits results to a separate
``<Project>Driver`` CDash project.  WARNING: This was written by a contractor
and is least well tested, more confusing, and least desirable part of TriBITS.
If you have a better way to manage multiple builds (e.g. Jenkins) then use
that instead.

**common_tpls/**: TPLs that are very common and are used by several different
TriBITS projects but are not built into the TriBITS system itslef. Having some
of these common TPLs in a central location enhances uniformity, reuse, and
makes it easier to pull TriBITS packages out of a repo and build them
independently.

**doc/**: Basic TriBITS documentation built using docutils. The generated
documents are stored in the git repo but instead are built on command using
the ``doc/build_docs.sh`` script.

**examples/**: Example TriBITS projects and TPLs. These can be copied out and
used as examples

**devtools_install/**: Basic install scripts for tools like CMake, GCC,
OpenMPI, MPICH, Git, etc. By default, these all download tarballs from the
``github.com/TriBITSPub/`` site and the repos are named
``devtools-<toolname>-<version>-base``. This makes it easy to set up a new dev
environment for projects that uses TriBITS (or don't use TriBITS for that
matter).

**win_interface/**: Some non-Windows C header files ported to Windows to make
porting to Windows easier.

The script ``snapshot_tribits.py`` install the different pieces for of this
``tribits/`` directory into a project's ``<projectDir>/cmake/tribits/``
subdirectory. It supports the argument ``--components`` with values ``core``,
``python_utils``, ``ci_support``, ``ctest_driver``, ``dashboard_driver``,
``common_tpls``, ``doc``, ``examples``, ``win_interface``, and
``devtools_install``. These snapshot components have the dependencies:

* ``core`` => (external CMake)
* ``python_utils`` => (external Python 2.4)
* ``win_interface`` => (external C compiler)
* ``TriBITS.cmake`` => ``core``
* ``ci_support`` => ``core``, ``python_utils``
* ``ctest_driver`` => ``core``, ``ci_support``
* ``dashboard_driver`` => ``ctest_driver``
* ``common_tpls`` => ``core``
* ``examples`` => (external tribits installation)
* ``doc`` => ``core``, ``ci_support``, ``examples``
* ``devtools_install`` => ``python_utils``
