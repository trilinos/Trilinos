Welcome to TriBITS's documentation!
===================================

:Author: Roscoe A. Bartlett (rabartl@sandia.gov)
:Date: |date|

.. |date| date::

.. toctree::
   :maxdepth: 3
   :hidden:
   :caption: Project, Build, Test, and Install:

   build_ref

`Project, Build, Test, and Install: <./build_ref.html>`_
--------------------------------------------------------

:Abstract: Provides a general project-independent reference on how to configure, build, test, and install a project that uses the TriBITS CMake build system.  The primary audience of this particular build of this document are TriBITS project developers themselves.  A project-specific version of this document should be created and accessed by users of a particular TriBITS-based project.

.. toctree::
   :maxdepth: 3
   :hidden:
   :caption: Users Guide:

   users_guide

`Users Guide: <./users_guide.html>`_
------------------------------------

:Abstract: This document describes the usage of TriBITS to build, test, and deploy complex software.  The primary audience are those individuals who develop on a software project which uses TriBITS as the framework for the project's CMake build system.  The overall structure of a TriBITS project is described including all of the various project- and package-specific files that TriBITS requires or can use and what order these files are processed.  It also contains detailed reference information on all of the various TriBITS macros and functions directly used in a TriBITS project's CMakeLists.txt files.  Many other topics of interest to a TriBITS project developer and architect are discussed as well.

.. toctree::
   :maxdepth: 3
   :hidden:
   :caption: Maintainers Guide:

   maintainers_guide

`Maintainers Guide: <./maintainers_guide.html>`_
------------------------------------------------

:Abstract: This document describes the internal implementation and the maintenance of the TriBITS project itself.  The primary audience are those individuals who will make changes and contributions to the TriBITS project or just want to understand its implementation details.  Included is all of the same information as the TriBITS Users Guide but also include file names and line numbers for all of the documented TriBITS macros and functions.
