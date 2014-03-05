========================
TriBITS Developers Guide
========================

:Author: Roscoe A. Bartlett (bartlettra@ornl.gov)

.. sectnum::

.. contents::


Introduction
=============

This document describes the usage of the TriBITS (Tribal Build, Integration,
Test System) to develop software projects.  An initial overview of TriBITS is
provided in the `TriBITS Overview <../overview/TribitsOverview.pdf>`_ document
which contains the big picture and provides a high-level road map to what
TriBITS provides.  This particular document, however, describes the details on
how to use the TriBITS system to create a CMake build system for a set of
compiled software packages.

TriBITS is a fairly extensive framework that is build on
CMake/CTest/CPack/CDash which in of itself is a very extensive system of
software and tools.  The most important thing to remember is that software
project that use TriBITS are really just CMake projects.  TriBITS makes no
attent to hide that either from the TriBITS project developers or from the
users that need to configure and build the software.  Therefore, to make
effective usage of TriBITS, one must learn the basics of CMake.  In
particular, CMake is a Turning complete programming lanauge with local and
global variables (with strange scoping rules), macros, functions, targets,
commands, and other features.  One needs to understand how to define and use
variables, macros, functions in CMake.  One needs to know how to debug
CMakeLists.txt files and CMake code in general (i.e. using ``MESSAGE()` print
statements).  One needs to understand how CMake defines and uses targets for
various qualities like libraries, executables, etc.  Without this basic
understanding of CMake, one will have trouble resolving problems when they
might occur.


Structure of a TriBITS Project
==============================

???

TriBITS Global Project Settings
===============================

TriBITS defines a number of global project-level settings that can be set by
the user and can have their default determined by each individual TriBITS
project.  If a given TriBITS project does not define its own default, a
reasonble default is set by the TriBITS system automatically.

ToDo: Document what parameters influence the entire TriBITS project, what
parameters can have project-specific defaults, etc.


TriBITS Macros and Functions
============================

.. include:: TribitsDetailedMacroFunctionDoc.rst


References
==========

.. [SCALE] http://scale.ornl.gov/


Appendix
========


History of TriBITS
------------------

TriBITS started development in November 2007 as a set of helper macros to
provide a CMake build system for a small subset of packages in Trilinos.  The
initial goal was to just to support a native Windows build (using Visual C++)
to compile and install these few Trilinos packages on Windows for usage by
another project (the Sandia Titan project which included VTK).  At that time,
Trilinos was using a highly customized autotools build system.  Initially,
this CMake system was just a set of macros to streamline creating executables
and tests.  Some of the conventions started in that early effort (e.g. naming
conventions of variables and macros where functions use upper case like old
FORTRAN and variables are mixed case) were continued in later efforts and are
reflected in the current.  Then, stating in early 2008, a more detailed
evaluation was performed to see if Trilinos should stitch over to CMake as the
default (and soon only) supported build and test system (see "Why CMake?" in
`TriBITS Overview <../overview/TribitsOverview.pdf>`_).  This lead to the
initial implementation of a scale-able package-based architecture
(PackageArch) for the Trilinos CMake project in late 2008.  This Trilinos
CMake PackageArch system evolved over the next few years with development in
the system slowing into 2010.  This Trilinos CMake build system was then
adopted as the build infrastructure for the CASL VERA effort in 2011 where
CASL VERA packages were treated as add-on Trilinos packages (see Section ???).
Over the next year, there was significant development of the system to support
larger multi-repo projects in support of CASL VERA.  That lead to the decision
to formally generalize the Trilinos CMake PackageArch build system outside of
Trilinos and the name TriBITS was formally adopted in November 2011.  Work to
refactor the Trilinos CMake system into a general reusable stand-alone
CMake-based build system started in October 2011 and an initial implementation
was complete in December 2011 when it was used for the CASL VERA build system.
In early 2012, the ORNL CASL-related projects Denovo and SCALE ([SCALE]_)
adopted TriBITS as their native development build systems.  Shortly after
TriBITS was adopted the native build system for the the CASL-related
University of Michigan code MPACT.  In addition to being used in CASL, all of
these codes also had a significant life outside of CASL.  Because they used
the same TriBITS build system, it proved relatively easy to keep these various
codes integrated together in the CASL VERA code meta-build.  At the same time,
TriBITS well served the independent development teams and non-CASL projects
independent from CASL VERA.  Since the initial extraction of TriBITS from
Trilinos, the TriBITS system was further extended and refined, driven by CASL
VERA development and expansion.  Independently, an early version of TriBITS
from 2012 was adopted by the LiveV
project\footnote{https://github.com/lifev/cmake} which was forked and extended
independently.
