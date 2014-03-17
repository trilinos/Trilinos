=======================================
TriBITS Developers Guide and Reference
=======================================

:Author: Roscoe A. Bartlett (bartlettra@ornl.gov)

:Abstract: This document describes the usage of TriBITS to build, test, and deploy complex software.  The primary audience are those individuals who develop on a software project which uses TriBITS.  The overall structure of a TriBITS project is descrided including all of the various project- and package-specific files that TriBITS requires or can use and how and what order these files are processed.  It also contains detailed reference information on all of the various TriBITS macros and functions.  Many other topics of interest to a TriBITS project developer and archetect are also discussed.

.. sectnum::
   :depth: 1

.. Above, the depth of the TOC is set to just 1 because I don't want the
.. TriBITS function/macro names to have section numbers appearing before them.
.. Also, some of them are long and I don't want them to go off the page of the
.. PDF document.

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

The remainder of this documented is structured as follows.  First, there is a
discussion of the various `TriBITS Developer and User Roles`_.  Then a brief
discussion of `CMake Language Overivew and Gotchas`_ is provided.

ToDo: Finish outline of the document.

The final sections are detailed documentation of `TriBITS Global Project
Settings`_, `TriBITS Macros and Functions`_, and `General Utility Macros and
Functions`_.


TriBITS Developer and User Roles
================================

There are approximately five different types roles with respect to TriBITS.
These different roles require different levels of expertise and knowlege of
CMake and knowledge of the TriBITS system.  The primary roles are 1) *TriBITS
Project User*, 2) *TriBITS Project Developer*, 3) *TriBITS Project Architect*,
4) *TriBITS System Developer*, and 5) *TriBITS System Architect*.  Each of
these roles builds on the necessary knolwege of the lower-level roles.
 
The first role is that of a **TriBITS Project User** who only needs to be able
to configure, build, and test a project that uses TriBITS as its build system.
A person acting in this role needs to know little about CMake other than
basics about how to run the ``cmake`` and ``ctest`` exectuables, how to set
CMake cache variables, and the basics of building software and running tests
with ``ctest``.  The proper reference for a TriBITS Project User is the
`Project-Specific Build Quick Reference`_.  Also, the `TriBITS
Overview <../overview/TribitsOverview.pdf>`_ document may be of some help
also.  A TriBITS project user does not need to know anything about the CMake
langauge itself or any of the TriBITS macros or functions described in
`TriBITS Macros and Functions`_ or really anything else described in this
current document.

A **TriBITS Project Developer** is someone who contributes to a software
project that uses TriBITS.  They will add source files, libraries and
exectuables, add test executables and define tests run with ``ctest``.  They
have to configure and build the project code in order to be able to develop
and run tests and therefore this role includes all of the necessary knowledge
and functions of a TriBITS Project User.  A casual TriBITS Project Developer
typically does not need to know a lot about CMake and really only need to know
a subset of the `TriBITS Macros and Functions`_ defined in this document.  A
slightly more sophsiticated TriBITS Project Developer will also add new
packages, add new package dependencies, and define new TPLs.  This current
TriBITS Developers Guide should supply everything such a developer needs to
know and more.  Only a smaller part of this document needs to be understood
and accessed by people assuming this role.

The next level of roles is a **TriBITS Project Architect**.  This is someone
(perhaps only one person on a project development team) that knows the usage
and functioning of TriBITS in great detail.  They understand how to set up a
TriBITS project from scrach, how to set up automated testing using the TriBITS
system, and know how to use TriBITS to implement the overall software
development process.  A person in this role is also likely to be the one who
makes the initial technical decision for their project to adopt TriBITS as is
native build and test system.  This document (along with detailed
CMake/CTest/CDash documentation provided by Kitware and the larger community)
should provide most of what a person in this role needs to know.  A person
assuming this role is the primary audience for this document.

The last two roles **TriBITS System Developer** and **TriBITS System
Architect** are for those individuals that actually extend and modify the
TriBITS system itself.  A TriBITS System Developer needs to know how to add
new functionlity while maintaining backward compatibility, how to add new unit
tests to the TriBITS system, and perform other related tasks.  Such a
developer needs to be very knowledgeable of the basic functioning of CMake and
know how TriBITS is implemented in the CMake language.  A TriBITS System
Architect is someone who must be consusted on almost all non-trivial changes
or additions to the TriBITS system.  A TriBITS System Architect in addition
needs to know the entire TriBITS system, the design philosophy that provides
the foundation for TriBIITS and be an expert in CMake, CTest, and CDash.
Everything that needs to be known by a TriBITS System Developer and a TriBITS
System Architect is not contained in this document.  Instead, the primary
documentation will be in the TriBITS CMake source code and various unit tests
itself.  At the time of this writing, there is currently there is only one
TriBITS System Architect (who also happens to be the primary author of this
document).

An explicit goal of this document is to make new TriBITS Project System
Archetects (i.e. those would make the decision to adopt TriBITS), and new
TriBITS System Developers to help extend and maintain the system.  As TriBITS
matures and its development stabilizes, the need for a TriBITS System
Architect will be diminished.

So depending on the particular role that a reader falls into, this documnet
may or may not be necessary but instead the TriBITS Overview or the
<Project>BuildQuickRef documents may be more appropriate.


CMake Language Overivew and Gotchas
====================================

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
One can also purchase the `offical CMake reference book`_.  Also, documenation
for any built-in CMake command is available locally by running::

   $ cmake --help-command <CMAKE_COMMAND>

Because tutorials and detailed documentation for the CMake language already
exists, this document will not even attempt to provide a first reference to
CMake (which is a large topic in itself).  However, what we try to provide
below is a short overivew of the more quarky or supprising aspects of the
CMake language that an programmer experienced in another language might get
tripped up by or surprised by.  Some of the more unique features of the
language are described in order to help avoid some of these common mistakes
and provide greater understanding of how TriBITS works.

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

The first thing to understand about the CMake language is that everthing line
of CMake code is just a command taking a string (or an array of strings) and
functions that operate on strings.  An array argument is just a single with
elements separated by semi-colons ``"<str0>;<str1>;..."``.  CMake is a bit odd
in how it deals with these arrays (which just represented as a string with
elements separated with semi-colons ``';'``).  For example, all of the
following are equivalent and pass in a CMake array with 3 elements [``A``],
[``B``], and [``C``]::

  SOME_FUNC(A B C)
  SOME_FUNC("A" "B" "C")
  SOME_FUNC("A;B;C")

However, the above is *not* the same as::

  SOME_FUNC("A B C")

which just passes in a single element with value [``A B C``].  Raw quotes in
CMake basically escapes the interpetation of space characters as array element
boundaries.  Quotes around arguments with no spaces does nothing (as seen
above).  In order to get a quote char [``"``] into string, you must escape it
as::

  SOME_FUNC(\"A\")

which passes an array with the single argument [``\"A\"``].

Varibles are set using a built-in CMke function that just takes string
arguments like::

  SET(SOME_VARIABLE "some_value")

In CMake, the above is idential, in every way, to::

  SET(SOME_VARIABLE some_value)
  SET("SOME_VARIABLE;"some_value")
  SET("SOME_VARIABLE;some_value")

The function ``SET()`` simply interprets the first argument to as the name of
a varible to set in the local scope.  Many other built-in and user-defined
CMake functions work the same way.  That is some of the string argumnets are
interpreted as the names of variables.

However, CMake appears to parse arguments differently for built-in CMake
control structure functions like ``FOREACH()`` and ``IF()`` and does not just
interpret them as a string array.  For example::

  FOREACH (SOME_VAR "a;b;c")
    MESSAGE("SOME_VAR='${SOME_VAR}'")
  ENDFOREACH()

prints ```SOME_VAR='a;b;c'`` instead of printing ``SOME_VAR='a'`` followed by
``SOME_VAR='b'``, etc., as you would otherwise expect.  Therefore, this simple
rule for the handling of function arguments as string arrays does not hold for
CMake logic control commands.  Just follow the CMake documentation for these
control structures..

CMake offers a rich assortment of built-in functions for doing all sorts of
things.  As part of these functions are the built-in ``MACRO()`` and the
``FUNCTION()`` functions which allow you to create user-defined macros and
functions (which is what TriBITS is built on).  All of these built-in and
user-defined macros and functions work exactly the same way; they take in an
array of string arguments.  Some functions take in positional arguments but
most actually take a combination of positional and keyword arguments (see
`PARSE_ARGUMENTS()`_).

Varible names are translated into their stored values using
``${SOME_VARIABLE}``.  The value that is extracted depends on if the varible
is set in the local or global (cache) scope.  The local scopes for CMake start
in the base project directory in its base ``CMakeLists.txt`` file.  Any
varibles that are created by macros in that base local scope are seen across
an entire project but are *not* persistent across ``cmake`` configure
invocations.

The handling of variables is one area where CMake is radically different from
most other languages.  First, a varible that is not defined simply returns
nothing.  What is surprising to most peoople about this is that it does not
even return an empty string!  For example, the following set statement::

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

which produces ``SOME_VAR='a;;c'``, or three elements as one might assume.

This is a common error that people make when they call CMake functions
(built-in or TriBITS-defined) involving varibles that might be undefined or
empty.  For example, for the macro::

   MACRO(SOME_MACRO  A_ARG  B_ARG  C_ARG)
      ...
   ENDMACRO()

if someone trys to call it with::

  SOME_MACRO(a ${SOME_OHTER_VAR} c)

and if ``SOME_OHTER_VAR=""`` or if it is undefined, then CMake will error out
with the error message saying that the macro ``SOME_MACRO()`` takes 3
arguments but only 2 were provided.  If a varible might be empty but that is
still a valid argument to the command, then it must be quoted as::

  SOME_MACRO(a "${SOME_OHTER_VAR}" c)

Related to this problem is that if you mispell the name of a variable in a
CMake ``IF()`` statement like::

   IF (SOME_VARBLE)
     ...
   ENDIF()

then it will always be false and the code inside the if statement will never
be executed!  To avoid this problem, use the utility function
`ASSERT_DEFINED()`_ as::

   ASSERT_DEFINED(SOME_VARBLE)
   IF (SOME_VARBLE)
     ...
   ENDIF()

In this case, the mispelled variable would be caught.

While on the subject of ``IF()`` statements, CMake has a strange convention.
When you say::

  IF (SOME_VAR)
    DO_SOMETHING()
  ENDIF()

then ``SOME_VAR` is interpreted as a variable and will be considered true and
``DO_SOMETHING()`` will be called if ``${SOME_VAR}`` does *not* evaluate to
``0``, ``OFF``, ``NO``, ``FALSE``, ``N``, ``IGNORE``, ``""``, or ends in the
suffix ``-NOTFOUND``.  How about that for a true/false rule!  To be safe, use
``ON/OFF`` and ``TRUE/FASLE`` pairs for setting variables.  Look up native
CMake documentation on ``IF()``.



CMake langauge behavior with respect to case sensitivity is also strange:

* Calls of built-in and user-defined macros and functions is *case
  insensitive*!  That is ``set(...)``, ``SET(...)``, ``Set()``, and all other
  combinations of upper and lower case characters for 'S', 'E', 'T' all call
  the bulit-in `SET()`` function.  The convention in TriBITS is to use all
  caps for functions and macros (was adopted by following the conventions used
  in the early versions of TriBITS, see `History of TriBITS`_).  The
  convention in CMake literature from Kitware seems to use lower-case letters
  for functions and macros.

* The names of CMake varables (local or cache/global) are *case sensitive*!
  That is, ``SOME_VAR`` and ``some_var`` are *different* variables.  Built-in
  CMake varibles tend use all caps with underscores
  (e.g. ``CMAKE_CURRENT_SOURCE_DIR``) but other built-in CMake varibles tend
  to use mixed case wtih underscores (e.g. ``CMAKE_Fortran_FLAGS``).  TriBITS
  tends to use a similar naming convention where most varibles have mostly
  upper-case letters except for proper nouns like the project, package or TPL
  name (e.g. ``TribitsProj_TRIBITS_DIR``, ``TriBITS_SOURCE_DIR``,
  ``Boost_INCLUDE_DIRS``).

I don't know of any other programming language that uses different case
senstivity rules for varibles verses functions.  However, because we must
parse macro and function arguments when writing user-defined macros and
functions, it is a good thing that CMake varibles are case sensitive.  Case
insenstivity would make it much harder and more expensive to parse argument
lists that take keyword-based arguments (see `PARSE_ARGUMENTS()`_).

Other mistakes that people make result from not understanding how CMake scopes
variables and other entities.  CMake defaults a global scope (i.e. "cache"
varibles) and several nested local scopes that are created by
``ADD_SUBDIRECTORY()`` and entering FUNCTIONS.  See `DUAL_SCOPE_SET()`_ for a
short discussion of these scoping rules.  It is not just varibles that can
have local and global scoping rules.  Other entities, like defines set with
the built-in command ``ADD_DEFINITIONS()`` only apply to the local scope and
child scopes.  That means that if you call ``ADD_DEFINITIONS()`` to set a
define that affects the meaning of a header-file in C or C++, for example,
that definition will *not* carry over to a peer subdirectory and those
definitions will not be set (see warning in `TRIBITS_ADD_LIBRARY()`_).

Now that some CMake basics and common gotchas have been reviewed, we now get
into the meat of TriBITS starting with the overall structure of a TriBITS
project.


Structure of a TriBITS Project
==============================

ToDo: Fill in!


Processing of TriBITS Files
===========================

ToDo: Fill in!

Automated testing
=================

ToDo: Define test group PT, ST, and EX

ToDo: Define test category BASIC, CONTINUOUS, NIGHTLY, WEEKLY, and PERFORMANCE

ToDo: Define repo category Continuous, Nightly, and Experimental which also map to CDash tracks

ToDo: Discuss the propery usage of these test categories and why NIGHTLY
testing should be the default.

ToDo: Fill in!

Multi-Repository Support
========================

ToDo: Discuss 'egdist', ExtraRepositoriesList.cmake, and the rep clone script.


Basic Development Workflow
==========================

ToDo: Fill in!

Multi-Repository Development Workflow
=====================================

ToDo: Fill in!

Project-Specific Build Quick Reference
======================================

TriBITS provides a mechanisms to quickly create a project-specific build quick
reference document in restructured text (RST) format and with HTML and
LaTeX/PDF outputs.  These document are generally created in the base project
source tree and given then name ``<Project>BuildQuickRef.[rst,html,pdf]``.
This document consists of two parts.  One part is a generic template document
``TribitsBuildQuickRefBody.rst`` that uses the place-holder ``<Project>`` that
is substituted for the for the real project name (read from the project's
``ProjectName.cmake`` file by default).  In order to produce this document, a
project must have the template file::

  <projectBaseDir>/cmake/<Project>BuildQuickRefTemplate.rst

defined which provides the outer RST doucment (with title, authors, abstract,
introduction, other introductory sections).  From this, the script::

  tribits/doc/build_quick_ref/create-project-build-quickref.py

is used to generate the read-only files::

  <projectBaseDir>/
    <Project>BuildQuickRef.rst
    <Project>BuildQuickRef.html
    <Project>BuildQuickRef.pdf

To see a simple example of this, see::

  tribits/doc/examples/TribitsExampleProject/cmake/create-build-quickref.sh

A project-indepenent version of this file is provided in the
`TribitsBuildQuickRef.[rsts,html,pdf]
<../build_quick_ref/TribitsBuildQuickRef.html>`_


Creating Source Distributions
=============================

ToDo: Fill in!


Multi-Repository Almost Continuous Integration
==============================================

ToDo: Fill in!


Regulated Backward Compatibility and Deprecated Code
====================================================

ToDo: Fill in!


TriBITS Global Project Settings
===============================

TriBITS defines a number of global project-level settings that can be set by
the user and can have their default determined by each individual TriBITS
project.  If a given TriBITS project does not define its own default, a
reasonble default is set by the TriBITS system automatically.  These options
are defined and are set, for the most part, in the internal TriBITS function
``TRIBITS_DEFINE_GLOBAL_OPTIONS_AND_DEFINE_EXTRA_REPOS()`` in the TriBITS
CMake code file ``TribitsGlobalMacros.cmake`` which gets called inside of the
``TRIBITS_PROJECT()`` macro.  That function and that file are the definitive
source the options that a TriBITS project takes and what the default values
are but we strive to document them here as well.  Many of these global options
(i.e. cache variables) such as ``${PROJECT_NAME}_<SOME_OPTION>`` allow the
project to define a default by setting a local varible
``${PROJECT_NAME}_<SOME_OPTION>_DEFAULT`` as::

  SET(${PROJECT_NAME}_<SOME_OPTION>_DEFAULT <someDefault>)

either in its top-level ``CMakeLists.txt`` file or in its
``ProjectName.cmake`` file.  If ``${PROJECT_NAME}_<SOME_OPTION>_DEFAULT`` is
not set by the project, then TriBITS provides a reasonable default value.  The
TriBITS code for this looks like::

  IF ("${${PROJECT_NAME}_<SOME_OPTION>_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_<SOME_OPTION>_DEFAULT <someDefault>)
  ENDIF()

  ADVANCED_SET( ${PROJECT_NAME}_<SOME_OPTION>
    ${PROJECT_NAME}_<SOME_OPTION>_DEFAULT}
    CACHE BOOL "[documentation]."
    )

where ``<SOME_OPTION>`` is the option name like ``TEST_CATEGORIES`` and
``<someDefault>`` is the default set by TriBITS if the project does not define
a default.  In this way, if the project sets the variable
``${PROJECT_NAME}_<SOME_OPTION>_DEFAULT`` before this code exeutates, then
``${${PROJECT_NAME}_<SOME_OPTION>_DEFAULT}`` will be used as the default for
the cache varible ``${PROJECT_NAME}_<SOME_OPTION>`` which, of course, can be
overridden by the user when calling ``cmake`` in a number of ways.

Most of these global options that can be overridden externally by setting the
cache variable ``${PROJECT_NAME}_<SOME_OPTION>`` should be documented in the
`Project-Specific Build Quick Reference`_ document.  A generic version of this
document is found in `TribitsBuildQuickRef.[rsts,html,pdf]
<../build_quick_ref/TribitsBuildQuickRef.html>`_.  Some of the more unusual
options that might only be of interest to developers mentioned below may not
be documented in ``<Project>BuildQuickRef.[rst,html,pdf]``.

The global project-level TriBITS options for which defaults can be provided by
a given TriBITS project are:

* `${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES`_
* `${PROJECT_NAME}_ENABLE_Fortran`_
* `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_
* `${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES`_
* `${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES`_
* `${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES`_
* `${PROJECT_NAME}_ELEVATE_ST_TO_PT`_
* `${PROJECT_NAME}_ENABLE_CPACK_PACKAGING`_
* `${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION`_
* `${PROJECT_NAME}_CPACK_SOURCE_GENERATOR`_
* `${PROJECT_NAME}_TEST_CATEGORIES`_
* `MPI_EXEC_MAX_NUMPROCS`_

These options are described below.

.. _${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES:

**${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES**

  If ``${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES`` is ``ON`` (the
  TriBITS default value), then any explicitly enabled packages that have
  disabled upstream required packages or TPLs will be disabled.  If ``OFF``,
  then an configure error will occur (for more details see
  `TribitsBuildQuickRef.* <../build_quick_ref/TribitsBuildQuickRef.html>`_).
  A project define a different default value by setting::
  
    SET(${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT FALSE)
  
  .. _${PROJECT_NAME}_ENABLE_Fortran:
  
**${PROJECT_NAME}_ENABLE_Fortran**
  
  If ``${PROJECT_NAME}_ENABLE_Fortran`` is ``ON``, then Fortran support for the
  project will be enabled and the Fortran compiler(s) must be found.  By
  default, TriBITS sets this to ``ON`` for non-Windows systems (i.e. ``WIN32``
  is not set by CMake) but is ``OFF`` for a Windows system.  A project always
  requires Fortran, for example, it can set the default::
  
    SET(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT TRUE)
  
  If a project does not have any native Fortran code a good default would be::
  
    SET(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT OFF)
  
  NOTE: It is usually not a good idea to always force off Fortran, or any
  compiler, because extra repositories and packages might be added by someone
  that might require the compiler and we don't want to unnecessarily limit the
  generality of a given TriBITS build.  Setting the default for all platforms
  should be sufficient.
  
.. _${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS:

**${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS**

  If ``${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`` is set to ``ON``, then
  any defined libraries or header files that are listed in calls to
  `TRIBITS_ADD_LIBRARY()`_ will be installed (unless options are passed into
  `TRIBITS_ADD_LIBRARY()`_ that disable installs).  If set to ``OFF``, then
  headers and librareis will be installed by default and only ``INSTALLABLE``
  executables added with `TRIBITS_ADD_EXECUTABLE()`_ will be installed.
  However, as described in `TribitsBuildQuickRef.*
  <../build_quick_ref/TribitsBuildQuickRef.html>`_, shared libraries will
  still be always be installed if enabled since they are needed by the
  installed executables.  The TriBITS default is to set this to ``ON``.
  
  For a TriBITS project that primarily is delivering libraries (e.g. Trilinos),
  then it makes sense to leave the TriBITS default or explicitly set::
  
    SET(${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT ON)
  
  For a TriBITS project that is primarily delivering executablers (e.g. VERA),
  then it makes sense to set the default to::
  
    SET(${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT OFF)
  
.. _${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES:

**${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES**
  
  If ``${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES`` is ``ON``, then
  ``Makefile.export.<PACKAGE_NAME>`` will get created at configure time in the
  build tree and installed into the install tree.  See `TribitsBuildQuickRef.*
  <../build_quick_ref/TribitsBuildQuickRef.html>`_ for details.  The TriBITS
  default is ``ON`` but a project can decide to turn this off by default by
  setting::
  
    SET(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES_DEFAULT OFF)
  
  A project might want to disable the generation of export makefiles by default
  if its main purpose is to provide executables.  There is no reason to provide
  an export makefile if libraies and headers are not actaully installed (see
  `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_)

.. _${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES:

**${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES**

  If ``${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES`` is set to ``ON``,
  then ``<PACKAGE_NAME>Config.cmake`` files are created at configure time in
  the build tree and installed into the install tree.  These files are used by
  external CMkae projects to pull in the list of compilers, compiler options,
  include directories and libraries.  The TriBITS default is ``ON``.  A
  project can change the default by setting, for example::

    SET(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES_DEFAULT OFF)

  A project would want to turn off the creation and installation of
  ``<PACKAGE_NAME>Config.cmake`` files if it was only installing and providing
  executables. See `TribitsBuildQuickRef.*
  <../build_quick_ref/TribitsBuildQuickRef.html>`_ for details.

.. _${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES:

**${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES**

  If ``${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES`` is ``ON``, then the
  data-structures needed to generate ``Makefile.export.<PACKAGE_NAME>`` and
  ``<PACKAGE_NAEM>Config.cmake`` are created.  These data structures are also
  needed in order to generate export makefiles on demand using the function
  `TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES()`_.  The default in
  TriBITS is to turn this ``ON`` automatically by default if
  ``${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES`` or
  ``${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES`` are ``ON``.  Else, by
  default, TriBITS sets this to ``OFF``.  The only reason for the project to
  override the default is to set it to ``ON`` as with:

    SET(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT ON)

  is so that the necessary data-structures are generated in order to use the
  function `TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES()`_.

.. _${PROJECT_NAME}_ELEVATE_ST_TO_PT:

**${PROJECT_NAME}_ELEVATE_ST_TO_PT**

  If ``${PROJECT_NAME}_ELEVATE_ST_TO_PT`` is set to ``ON``, then all ``ST`` SE
  packages will be elevated to ``PT`` packages.  The TriBITS default is
  obviously ``OFF``.  The default can be changed by setting::

    SET(${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT ON)

  There are projects, especially meta-projects, where the distiction between
  ``PT`` and ``ST`` code is not helpful or the assignment of ``PT`` and ``ST``
  packages in a repository is not appropriate.  An example project like this
  CASL VERA.  Changing the default to ``ON`` allows any packages to be
  considered in pre-push testing.

.. _${PROJECT_NAME}_ENABLE_CPACK_PACKAGING:

**${PROJECT_NAME}_ENABLE_CPACK_PACKAGING**

  If ``${PROJECT_NAME}_ENABLE_CPACK_PACKAGING`` is ``ON``, then CPack support
  is enabled and some TriBITS code is avoided that is needed to set up
  data-structures that are used by the built-in CMake target
  ``package_source``.  The TriBITS default is ``OFF`` with the idea that the
  average developer or user will not be wanting to create source distributions
  with CPack.  However, this default can be changed by setting::

    SET(${PROJECT_NAME}_ENABLE_CPACK_PACKAGING ON)

.. _${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION:

**${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION**

  If ``${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION`` is
  ``TRUE``, then the directories for subpackages that are not enabled are left
  out of the source tarball.  This reduces the size of the tarball as much as
  possible but does require that the TriBITS packages and subpackages be
  properly set up to allow disabled subpackages from being excluded.  The
  TriBITS default is ``TRUE`` but this can be changed by setting::

    SET(${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION_DEFAULT FALSE)

.. _${PROJECT_NAME}_CPACK_SOURCE_GENERATOR:

**${PROJECT_NAME}_CPACK_SOURCE_GENERATOR**

  The variable ``${PROJECT_NAME}_CPACK_SOURCE_GENERATOR`` determines the CPack
  source generation types that are created when the ``package_source`` target
  is run.  The TriBITS default is set to ``TGZ``.  However, this default can
  be overridded by setting, for example::

    SET(${PROJECT_NAME}_CPACK_SOURCE_GENERATOR_DEFAULT "TGZ;TBZ2")

  This variable should generally be set in the file::

     <projectDir>/cmake/CallbackDefineProjectPackaging.cmake

  instead of in the base-level ``CMakeLists.txt`` file so that it goes along
  with rest of the project-specific CPack packaging options.

.. _${PROJECT_NAME}_TEST_CATEGORIES:

**${PROJECT_NAME}_TEST_CATEGORIES**

  The cache varaible ``${PROJECT_NAME}_TEST_CATEGORIES`` detemines what tests
  defined using `TRIBITS_ADD_TEST()`_ and `TRIBITS_ADD_ADVANCED_TEST()`_ will
  be added for ``ctest`` to run (see `Automated testing`) for
  discussion of test categories).  The TriBITS default is ``NIGHTLY`` for a
  standard local build.  The ``checkin-test.py`` script sets this to
  ``BASIC``.  A TriBITS project can override the default for a basic configure
  using, for example::

    SET(${PROJECT_NAME}_TEST_CATEGORIES BASIC)

  The justification for having the default test category be ``NIGHTLY``
  instead of ``BASIC`` is that when someone is enabling a package to develop
  on it or install it, we want them by default to be seeing the full version
  of the test suite (shy of the ``WEEKLY`` tests which can be very expensive)
  for the packages they are explictly enabling.  Typically they will not be
  enabling forward (downstream) dependent packages so the cost of running the
  test suite should not be too prohibitive.  This all depends on how good of a
  job the development teams do in making their test suites run fast and
  keeping the cost of running the tests down.  See the section `Automated
  testing`_ for a more detailed discussion.

.. _MPI_EXEC_MAX_NUMPROCS:

**MPI_EXEC_MAX_NUMPROCS**

  The varaible ``MPI_EXEC_MAX_NUMPROCS`` gives the maximum number of processes
  for an MPI test that will be allowed as defined by `TRIBITS_ADD_TEST()`_ and
  `TRIBITS_ADD_ADVANCED_TEST()`_.  The TriBITS default is set to be ``4`` (for
  no good reason really but it needs to stay that way for backward
  compatibility).  This default can be changed by setting::

    SET(MPI_EXEC_MAX_NUMPROCS_DEFAULT <newDefaultMax>)

  While this default can be changed for the project as a whole on all
  platforms, it is likely better to change this default on a
  machine-by-machine basis to correspond to the loat that can be accomidated
  by a given machine (or class of machines).  For example if a given machine
  has 64 cores, a reasonble number for ``MPI_EXEC_MAX_NUMPROCS_DEFAULT`` is
  64.


TriBITS Macros and Functions
============================

The following subsections give detailed documentation for the CMake macros and
functions that make up the core TriBITS system.  These are what are used by
TriBITS project developers in their ``CMakeLists.txt`` and other files.  These
are listed in approximately the order they will be encounted in a project or
packages ``CMakeLists.txt`` and other files.  All of these functions and
macros should be aviable when processing the project's and package's variables
files if used properly.  Therefore, no explicit ``INCLUDE()`` statements
should be needed other than the initial include of the
``TribitsProject.cmake`` file in the top-level CMakeLists.txt file so the
command `TRIBITS_PROJECT()`_ can be executed.

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
conflicit.  When overridding a built-in command ``some_builtin_command()``,
you can always access the original built-in command as
``_some_builtin_command()``.


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
CASL VERA packages were treated as add-on Trilinos packages (see Section
`Multi-Repository Support`_).  Over the next year, there was significant
development of the system to support larger multi-repo projects in support of
CASL VERA.  That lead to the decision to formally generalize the Trilinos
CMake PackageArch build system outside of Trilinos and the name TriBITS was
formally adopted in November 2011.  Work to refactor the Trilinos CMake system
into a general reusable stand-alone CMake-based build system started in
October 2011 and an initial implementation was complete in December 2011 when
it was used for the CASL VERA build system.  In early 2012, the ORNL
CASL-related projects Denovo and SCALE ([SCALE]_) adopted TriBITS as their
native development build systems.  Shortly after TriBITS was adopted the
native build system for the the CASL-related University of Michigan code
MPACT.  In addition to being used in CASL, all of these codes also had a
significant life outside of CASL.  Because they used the same TriBITS build
system, it proved relatively easy to keep these various codes integrated
together in the CASL VERA code meta-build.  At the same time, TriBITS well
served the independent development teams and non-CASL projects independent
from CASL VERA.  Since the initial extraction of TriBITS from Trilinos, the
TriBITS system was further extended and refined, driven by CASL VERA
development and expansion.  Independently, an early version of TriBITS from
2012 was adopted by the LiveV project\footnote{https://github.com/lifev/cmake}
which was forked and extended independently.
