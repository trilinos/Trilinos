=================================
TriBITS Users Guide and Reference
=================================

:Author: Roscoe A. Bartlett (rabartl@sandia.gov)
:Date: |date|
:Version: .. include:: ../TribitsGitVersion.txt

.. |date| date::

:Abstract: This document describes the usage of TriBITS to build, test, and deploy complex software.  The primary audience are those individuals who develop on a software project which uses TriBITS as the framework for the project's CMake build system.  The overall structure of a TriBITS project is described including all of the various project- and package-specific files that TriBITS requires or can use and what order these files are processed.  It also contains detailed reference information on all of the various TriBITS macros and functions directly used in a TriBITS project's CMakeLists.txt files.  Many other topics of interest to a TriBITS project developer and architect are discussed as well.

.. sectnum::
   :depth: 2

.. Above, the depth of the TOC is set to just 2 because I don't want the
.. TriBITS function/macro names to have section numbers appearing before them.
.. Also, some of them are long and I don't want them to go off the page of the
.. PDF document.

.. Sections in this document use the underlines:
..
.. Level-1 ==================
.. Level-2 ------------------
.. Level-3 ++++++++++++++++++
.. Level-4 ..................

.. contents::


Introduction
=============

This document describes the usage of the TriBITS (Tribal Build, Integration,
Test System) to develop software projects.  An initial overview of TriBITS is
provided in the `TriBITS Overview`_ document which contains the big picture
and provides a high-level road map to what TriBITS provides.  This particular
document, however, describes the details on how to use the TriBITS system to
create a CMake build system for a set of compiled software packages.  Also
described are an extended set of tools and processes to create a complete
software development, testing, and deployment environment consistent with
modern agile software development best practices.

TriBITS is a fairly extensive framework that is built on top of the
open-source CMake/CTest/CPack/CDash system (which in itself is a very
extensive system of software and tools).  The most important thing to remember
is that a software project that uses TriBITS is really just a CMake project.
TriBITS makes no attempt to hide that fact either from the TriBITS project
developers or from the users that need to configure and build the software.
Therefore, to make effective usage of TriBITS, one must learn the basics of
CMake (both as a developer and as a user).  In particular, CMake is a
Turning-complete programming language with local variables, global variables, macros, functions, targets, commands, and other
features.  One needs to understand how to define and use variables, macros,
and functions in CMake.  One needs to know how to debug CMakeLists.txt files
and CMake code in general (i.e. using ``message()`` print statements).  One
needs to understand how CMake defines and uses targets for various qualities
like libraries, executables, etc.  Without this basic understanding of CMake,
one will have trouble resolving problems when they occur.

The remainder of this documented is structured as follows.  First, there is
some additional `Background`_ material provided.  Then, a detailed
specification of `TriBITS Project Structure`_ is given which lists and defines
all of the files that a TriBITS project contains and how they are processed.
This is followed up by short descriptions of `Example TriBITS Projects`_ that
are provided with the TriBITS source tree that are used throughout this
document.  The topic of `Package Dependencies and Enable/Disable Logic`_ is
then discussed which is the backbone of the TriBITS system.  An overview of
the foundations for `TriBITS Automated Testing`_ is then given.  The topic of
TriBITS `Multi-Repository Support`_ is examined next.  `Development
Workflows`_ using TriBITS is then explored.  This is followed by a set of
detailed `Howtos`_.  Later some `Miscellaneous Topics`_ are presented that
don't fit well into other sections.  Then the main bulk of the detailed
reference material for TriBITS is given in the section `TriBITS Detailed
Reference Documentation`_.  Finally, several bits of information are provided
in the `Appendix`_.


.. include:: ../TribitsGuidesBody.rst


.. include:: TribitsCoreDetailedReference.rst


.. include:: ../TribitsGuidesReferences.rst


.. include:: ../TribitsFAQ.rst


Appendix
========

.. include:: ../TribitsCMakeLanguageOverviewAndGotchas.rst

.. include:: ../TribitsHistory.rst

.. include:: ../TribitsPackageNotCMakePackage.rst

.. include:: ../TribitsDesignConsiderations.rst

.. include:: ../TribitsToolsDocumentation.rst


.. ***
.. *** Extra References
.. ***

.. _tribits_read_all_project_deps_files_create_deps_graph(): TribitsMaintainersGuide.html#tribits-read-all-project-deps-files-create-deps-graph

.. _${PACKAGE_NAME}_LIB_DEFINED_DEPENDENCIES: TribitsMaintainersGuide.html#package-name-lib-enabled-dependencies

.. _${PACKAGE_NAME}_TEST_DEFINED_DEPENDENCIES: TribitsMaintainersGuide.html#package-name-test-enabled-dependencies
