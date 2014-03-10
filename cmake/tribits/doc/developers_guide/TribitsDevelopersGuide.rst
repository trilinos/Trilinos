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
CMakeLists.txt files and CMake code in general (i.e. using ``MESSAGE()`` print
statements).  One needs to understand how CMake defines and uses targets for
various qualities like libraries, executables, etc.  Without this basic
understanding of CMake, one will have trouble resolving problems when they
might occur.


TriBITS Develoepr and User Roles
================================

ToDo: Discuss the three different primary roles for related to TriBITS (core
TriBITS system development, TriBITS project architect, TriBITS project
developer, TriBITS project user).


Brief CMake Language Tutorial
==============================

TriBITS removes a lot of the boiler plate code needed to write a CMake
project.  As a result, many people can come into a project that uses TriBITS
and quickly start to contribute by adding new source files, adding new
libraries, adding new tests, and even adding new TriBITS packages and TPLs;
all without really having learned anything about CMake.  One just needs to
copy-and-paste existing example CMake code and files as basically "monkey see,
monkey do".  As long as nothing out of the ordinary happens, many people can
get along just fine in this mode for a time.

However, we have observed that most mistakes that people make when using
TriBITS, and most of the problems they have when using the sytem, are due to a
basic lack of knowlege of the CMake language.  One can find basic tutorials
and references on the CMake language in various locations online for free.
One can also purchase the `offical CMake reference book`_.  Therefore, this
document will not even attempt to provide a first reference to CMake (which is
a large topic in itself).  However, what we try to provide below is a short
overivew of the CMake langauge and a description of its unique features in
order to help avoid some of these common mistakes and provide greater
understanding of how TriBITS works.

.. _Offical CMake reference book: http://www.cmake.org/cmake/help/book.html

The CMake language that is used to write CMake projects with TriBITS (and that
core TriBITS itself is implemented in) is a fairly simply programming languge
with fairly simple rules (for the most part).  However, compared to other
programming lanuages, there are a few peculiar aspects to the CMake language
like strange varible scoping rules, arguments to macros and function, that can
make working with it difficult if you don't understand these.  Also, CMake has
some interesting gotchas.  In order to effectively use TriBITS (or just raw
CMake) to construct and maintain a project's CMake files, one must know the
basic rules of CMake.  

The first thing to understand about the CMake language is that everthing is
just a string (or an array of strings) and functions that operate on strings.
An array argument is just a single with elements separated by semi-colons
"<str0>;<str1>;...".

Varibles are set using a built-in CMke function that just takes string
arguments like::

  SET(SOME_VARIABLE "some_value")

In CMake, the above is idential, in every way, to::

  SET(SOME_VARIABLE some_value)

or::

  SET("SOME_VARIABLE;"some_value")

The function ``SET()`` simply interprets the first argument to as the name of
a varible to set in the local scope.  Many other built-in and user-defined
CMake functions work the same way.

CMake offers a rich assortment of built-in functions for doing all sorts of
things.  As part of these functions are the built-in ``MACRO()`` and the
``FUNCTION()`` functions which allow you to create user-defined macros and
function.  All of these built-in and user-defined macros and functions work
exactly the same way; they take in an array of string arguments.  Some
functions take in positional arguments but most actually take a combination of
positional and keyword arguments.

Varible names are translated into their stored values using
``${SOME_VARIABLE}``.  The value that is extracted depends on if the varible
is set in the local or global (cache) scope.  The local scopes for CMake start
in the base project directory in its ``CMakeLists.txt`` file.  Any varibles
that are created by macros in that base local scope are seen across an entire
project but are *not* persistent across ``cmake`` configure invocations.

The handling of variables is one area where CMake is radically different from
most other languages.  First, a varible that is not defined simply returns
nothing.  What is surprising to most peoople about this is that it does not
even return an empty string.  For example, the following set statement::

   SET(SOME_VAR a ${SOME_UNDEFINED_VAR} c)

produces ``SOME_VAR='a;c'`` and *not* ``'a;;c'``!  The same thing occurs when
an empty varible is dereferenced such as with::

   SET(EMPTY_VAR "")
   SET(SOME_VAR a ${EMPTY_VAR} c)

which produces ``SOME_VAR='a;c'`` and *not* ``'a;;c'``.  In order to always
produce an element in the array even if the varible is empty, one must quote
the argument as with::

   SET(EMPTY_VAR "")
   SET(SOME_VAR a "${EMPTY_VAR}" c)

which produces ``SOME_VAR='a;;c'``, or three elements as one might assue.

This is a common error the people make when they call functions (built-in or
TriBITS-defined) involving varibles that might be undefined or set to the
empty string.  For example, for the macro::

   MACRO(SOME_MACRO A_ARG B_ARG C_ARG)
      ...
   ENDMACRO()

If someone trys to call it with::

  SOME_MACRO(a ${SOME_OHTER_VAR} c)

and if ``SOME_OHTER_VAR=""`` or if it is undefined, then CMake will error out
with the error message saying that the macro ``SOME_MACRO()`` takes 3
arguments but only 2 were provided.  If a varible might be empty but that is
still a valid argument to a function (or element in a general array variable,
then it must be quoted as::

  SOME_MACRO(a "${SOME_OHTER_VAR}" c)

Related to this problem is that if you mispell the name of a variable in a
CMake if statement like::

   IF (SOME_VARBLE)
     ...
   ENDIF()

then it will always be false.  To avoid this problem, use the utility function
`ASSERT_DEFINED()`_ as::

   ASSERT_DEFINED(SOME_VARBLE)
   IF (SOME_VARBLE)
     ...
   ENDIF()

In this case, the mispelled variable would be caught.

A quick note of some strange CMake langauge behavior is case sensitivity:

* Calls of built-in and user-defined functions is *case insensitive*!  That is
  ``set(...)``, ``SET(...)``, ``Set()``, and all other combinations of upper
  and lower case characters for 'S', 'E', 'T' all call the bulit-in `SET()``
  function.  The convention in TriBITS is to use all caps for functions and
  macros.  The convention in CMake literature from Kitware seems to use
  lower-case for functions and macros.

* The names of CMake varables (local or cache/global) are *case sensitive*!
  That is, ``SOME_VAR`` and ``some_var`` are *different* variables.  Built-in
  CMake varibles tend use all caps with underscores
  (e.g. ``CMAKE_CURRENT_SOURCE_DIR``) but other built-in CMake varibles tend
  to use mixed case wtih underscores (e.g. ``CMAKE_Fortran_FLAGS``).  TriBITS
  tends to use a similar naming convention where most varibles have mostly
  upper-case letters except for proper nouns like the project, package or TPL
  name (e.g. ``TribitsProj_TRIBITS_DIR``, ``TriBITS_SOURCE_DIR``,
  ``Boost_INCLUDE_DIRS``).

I don't now of any other language that uses different case senstivity rules
for varibles verses functions.  However, because we must parse macro and
function arguments when writing user-defined macros and functions, it is a
good thing that CMake varibles are not case insensitive.  Case insenstivity
would make it much harder and more expensive to parse argument lists (see 

The other mistakes that people make is not understanding how CMake scopes
variables and other entities.  CMake defaults a global scope (i.e. "cache"
varibles) and several nested local scopes that are created by
``ADD_SUBDIRECTORY()`` and entering FUNCTIONS.  See `DUAL_SCOPE_SET()`_ for a
short discussion of these scoping rules.







???


Structure of a TriBITS Project
==============================

???


Processing of TriBITS Files
===========================

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

The following subsections give detailed documentation for the CMake macros and
functions that make up the core TriBITS system.  These are what are used by
TriBITS project developers in their ``CMakeLists.txt`` and other files.  These
are listed in approximately the order they will be encounted in a project or
packages ``CMakeLists.txt`` and other files.

.. include:: TribitsMacroFunctionDoc.rst


General Utility Macros and Functions
=====================================

The following subsections give detailed documentation for some CMake macros
and functions which are *not* a core part of the TriBITS system but are
included in the TriBITS system that are used inside of the TriBITS system and
are provided as a convenience to TriBITS project developers.  One will see
many of these functions and macros used throughout the implementation of
TriBITS and even in the ``CMakeLists.txt`` files for projects that use
TriBITS.

These macros and functions are *not* prefixed with ``TRIBITS_``.  There is
really not a large risk to defining and using these non-namespaces utility
functions and macros.  It turns out that CMake allows you to redefine any
macro or function, even built-in ones, inside of your project so even if CMake
did add new comamnds that clashed with these names, there would be no
conflicit.  When overridding a built-in command ``some_bultin_command()``, you
can always access the original built-in command as ``_some_bultin_command()``.


.. include:: UtilsMacroFunctionDoc.rst


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
