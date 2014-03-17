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

TriBITS Project Structure
=========================

TriBITS is a framework, implemented in CMake to create CMake projects.  As a
framework, it defines the the overall structure of a CMake build system for a
project and processes project, repository, and package specific files in a
specified order.  All of this processing takes place in the
`TRIBITS_PROJECT()`_ macro.

TriBITS Structural Units
------------------------

A CMake project that uses TriBITS as its build and test system is composed of
a single *TriBITS Project*, one or more *TriBITS Repositories* and one or more
*TriBITS Packages*.  In addition, a TriBITS Package can be broken up into
*TriBITS Subpackages*.  Together, the collection of TriBITS Packages and
TriBITS Subpackages are called *TriBITS Software Engineering Packages*, or
*TriBITS SE Packages* for short.

First, to establish the basic nomenclature, the key structural TriBITS units
are:

* `TriBITS Package`_: A collection of related software that typically includes
  one or more source files built into one or more libraries and has assoicated
  tests to help define and protect the functionality provided by the software.
  A package also typically defines a unit of documentation and testing (see
  `Automated testing`_).  A TriBITS package may or may not be broken down into
  multiple subpackages. Examples of TriBITS packages in
  ``TribitsExampleProject`` include ``SimpleCXX``, ``MixedLanguage`` and
  ``PackageWithSubpackages``.

* `TriBITS Subpackage`_: A part of a parent package that also typically has
  source files built into libraries and tests but is documented and tested
  along with the other subpackages the parent package.  The primary purpose
  for supporting subpackages to provide finer-grained of control software
  dependencies.  In ``TribitsExampleProject``, ``PackageWithSubpackages`` is
  an example of a package with subpackages ``SubpackaeA``, ``SubpackaeB``, and
  ``SubpackageC``.  The full subpackage name has the parent package name
  prefixed the the subpackage name
  (e.g. ``PackageWithSubpackagesSubpackageA``).  The parent package is always
  implicitly dependent on its subpackages.

* **TriBITS SE Package**: The combined set of TriBITS packages and
  subpackages. SE packages are the basis for setting dependencies in the
  TriBITS system.  For example, the SE Packages provided by the example
  ``PackageWithSubpackages`` is (in order of increasing dependencies)
  ``PackageWithSubpackagesSubpackaeA``, ``PackageWithSubpackagesSubpackaeB``,
  ``PackageWithSubpackagesSubpackageC``, and ``PackageWithSubpackages``.

* `TriBITS TPL`_: The specification for a particular external dependency that
  is required or can be used in one or more TriBITS SE Packages.  A TPL (a
  Third Party Library) typically provides a list of libraries or a list
  include directories for header files but can also be manifisted in order
  ways as well.  Examples of basic TPLs include ``BLAS``, ``LAPACK``, and
  ``Boost``.

* `TriBITS Repository`_: A collection of one or more TriBITS packages
  specified in a ``PackagesList.cmake`` file.

* `TriBITS Project`_: A collection of TriBITS Repositories and Packages that
  defines a CMake ``PROJECT`` and can be configured, built, and tested.

The following subsections define the major structural units of a TriBITS
project in more detail.  Each structural unit is described along with the
files and directories assoicated with each.  In addition, a key set of TriBITS
CMake variables for each are defined as well.

In the next major section following this one, some `Example TriBITS Projects`_
are described.  For those who just want to jump in and learn best by example,
these exmaple projects are a good way to start.  These example projects will
be referenced in the more detailed descriptions given in this document.

The CMake varaibles defined by TriBITS described in the structural using below
below fall into one of two types:

* *Local Fixed-Name Variables* are used temporarily in the processing of a
  TriBITS unit.  These include variables such as ``PROJECT_NAME``,
  ``REPOSITORY_NAME``, ``PACKAGE_NAME``, and ``PARENT_PACKAGE_NAME`` These are
  distinguished by having a fixed/constant name.  They are typically part of
  TriBITS reflection system, allowing subordinate units to determine the
  encapsulating unit in which they are participating.  For, example, a TriBITS
  subpackage can determine its name, its parent package's name and
  directories, its parent repository name and directories, and the enclosing
  project's name and directories.
* *Namespaced Variables* are used to refer to properties of a named TriBITS
  unit.  These include variables such as ``${REPOSTORY_NAME}_SOURCE_DIR``
  (e.g. ``TribitsExProj_SOURCE_DIR``) and ``${PACKAGE_NAME}_BINARY_DIR``
  (e.g. ``SimpleCXX_BINARY_DIR``).  They are available after processing the
  unit, for use by downstream or subordinate units.  They are part of the
  TriBITS dependency system, allowing downstream units to access properties of
  their known upstream dependencies.

More information about these various files is described in section `Processing
of TriBITS Files: Ordering and Details`_.


TriBITS Project
+++++++++++++++

* This defines a complete CMake project build by calling
  ``PROJECT(${PROJECT_NAME} ...)``.
* Consists of one or more TriBITS Repositories (see `TriBITS Repository`_).
* Defines ``PROJECT_NAME`` CMake variable (defined in ``ProjectName.cmake``)
* Defines a set of native Repositories (see below) that define packages and
  TPLs.  This set of native Repositories defines the official packages list
  and dependencies data-structure.  The list of native Repositories can be
  empty!
* Allows for extra Repositories to be added on before and after the set of
  native Repositories (given in ``ExtraRepositoriesList.cmake`` or by CMake
  variables)
* Defines a default CDash server and default project name on the server (the
  project name on the CDash server must be the same as ``${PROJECT_NAME})``.

**TriBITS Project Files**

The core files making up a TriBITS Project are::

  ${${PROJECT_SOURCE_DIR}}/
     CMakeLists.txt
     CTestConfig.cmake
     ProjectName.cmake    # Defines PACAKGE_NAME
     Version.cmake        # [Optional] Dev mode, Project version, VC branch
     cmake/
       NativeRepositoriesList.cmake    # [Optional] Used for some meta-projects
       ExtraRepositoriesList.cmake     # [Optional] Lists repos and VC URLs 
       ProjectDependenicesSetup.cmake  # [Optional] Project deps overrides
       CallbackDefineProjectPackaging.cmake  # [Optional] CPack settings
       tribits/    # [Optional] Or provide ${PROJECT_NAME}_TRIBITS_DIR
       ctest/
         CTestCustom.ctest.in  # [Optional] Custom ctest settings

ToDo: Document each file in detail here!

**TriBITS Project Core Variables**

The following local variables are defined in the top-level Project
``CMakeLists.txt`` file scople and are therefore acceessible by all files
prcessed by TriBITS:


  ``PROJECT_NAME``

    The name of the TriBTS Project.  This exists to support, among other
    things, the ability for subordinate units (Repositories and Packages) to
    determine the Project in which is participating.  This is typically read
    from a ``SET()`` statement in the project's ``ProjectName.cmake`` file.

  ``PROJECT_SOURCE_DIR``

    The absolute path to the base Project source directory.  This is set
    automatically by TriBITS given the directory passed into ``cmake`` at
    configure time at the beginning of the `TRIBITS_PROJECT()`_ macro.

  ``PROJECT_BINARY_DIR``

    The absolute path to the base Project binary/build directory.  This is set
    automatically by TriBITS and is the directory where ``cmake`` is run from
    and is set at the beginning of the `TRIBITS_PROJECT()`_ macro

  ``${PROJECT_NAME}_SOURCE_DIR``

    Set to ``${PROJECT_SOURCE_DIR}`` automatically by TriBITS at the beginning
    of the `TRIBITS_PROJECT()`_ macro.

  ``${PROJECT_NAME}_BINARY_DIR``

    Set to ``${PROJECT_BINARY_DIR}`` automatically by TriBITS at the beginning
    of the `TRIBITS_PROJECT()`_ macro.

  ``${PACKAGE_NAME}_ENABLE_TESTS``

    CMake cache varaibles that if set to ``ON``, then tests for all explicitly
    enabled packages will be turned on.

  ``${PACKAGE_NAME}_ENABLE_EXAMPLES``

    CMake cache varaibles that if set to ``ON``, then examples for all
    explicitly enabled packages will be turned on.

TriBITS Repository
++++++++++++++++++

A TriBITS Repository is the basic unit of ready-made composition between
different collections of software that use the TriBITS CMake build and system.

In short, a TriBITS Repository:

* Is a collection of related TriBITS Packages and TPLs
* Defines an association between one or more packages
* Defines a CMake variable specific to the collection referred to the
  in the variable ``REPOSITORY_NAME``.
* Defines the base source and binary directories for the Repository
  ``${REPOSITORY_NAME}_SORUCE_DIR`` and ``${REPOSITORY_NAME}_BINARY_DIR``.
* [Optional] Defines a list of add-on packages and TPLs (in ``PackagesList.cmake``
  and ``TPLsList.cmake``)
* [Optional] Defines a common set of initializations and other hooks for a
  collection of projects.

**TriBITS Repository Core Files**

The core files making up a TriBITS Repository are::

  ${${REPOSITORY_NAME}_SOURCE_DIR}/
    PackagesList.cmake
    TPLsList.cmake
    Copyright.txt  # [Optional] Only needed if creating version header file
    Version.cmake  # [Optional] Info inserted into ${REPO_NAME}_version.h
    cmake/
       RepositoryDependenciesSetup.cmake # [Optional]
       CallbackSetupExtraOptions.cmake # [Optional] Called after tribits options
       CallbackDefineRepositoryPackaging.cmake # [Optioinal] CPack packaging

ToDo: Document each file in detail here!

**TriBITS Repository Core Varaibles**

The following local variables are defined automatically by TriBITS before
processing a given TriBITS repositories files (e.g. ``PackagesList.cmake``,
``TPLsList.cmake``, etc.):

  ``REPOSITORY_NAME``

    The name of the current TriBITS repository.

  ``REPOSITORY_DIR``

    Path of the current Repository relative to the Project directory.  This is
    typically just the repository name but can be an arbitrary directory if
    specified through a ``ExtraRepositoriesList.cmake`` file.

The following base project-scope local variables are available once the list
of TriBITS repositories are defined:

  ``${REPOSITORY_NAME}_SOURCE_DIR``

    The absolute path to the base of a given Repository source directory.
    CMake code, for example in a packages's ``CMakeLists.txt`` file, typically
    refers to this by the raw name like ``RepoX_SOURCE_DIR``.  This makes such
    CMake code independent of where the various TriBITS repos are in relation
    to each other or the Project.

  ``${REPOSITORY_NAME}_BINARY_DIR``

    The absolute path to the base of a given Repository binary directory.
    CMake code, for example in packages, refer to this by the raw name like
    ``RepoX_SOURCE_DIR``.  This makes such CMake code independent of where the
    various TriBITS repos are in relation to each other or the Project.

Before processing the various repositories files, the following internal
project-scope local CMake variables are defined by TriBITS:


  ``${PROJECT_NAME}_NATIVE_REPOSITORIES``

     The list of Native Repositories for a given TriBITS project
     (i.e. Repositories that are always present when configuring the Project).
     If the file ``${PROJECT_SOURCE_DIR}/NativeRepositoriesList.cmake``
     exists, then the list of native repositories will be read from that file.
     If the file ``NativeRepositoriesList.cmake`` does not exist, then the
     project is assumed to also be a repository and the list of native
     repositories is just the local project directory
     ``${PROJECT_SOURCE_DIR}/.``.  In this case, the
     ``${PROJECT_SOURCE_DIR}/`` must contain at a minumum a
     ``PackagesList.cmake``, and a ``TPLsList.cmake`` file.

  ``${PROJECT_NAME}_EXTRA_REPOSITORIES``

     The list of Extra Repositories that the project is being configured with.
     The packages in these repositories are *not* listed in the main project
     dependency files but are listed in the dependenicy files in other
     contexts.  This list of repositories either comes from the project's
     ``ExtraRepositoriesList.cmake`` file or comes from the CMake variables
     ${PROJECT_NAME}_EXTRA_REPOSITORIES.  See `Enabling extra repositories
     with add-on packages
     <../build_quick_ref/TribitsBuildQuickRef.html#enabling-extra-repositories-with-add-on-packages>`_
     for details.

  ``${PROJECT_NAME}_ALL_REPOSITORIES``

    Concatenation of all the repos listed in
    ``${PROJECT_NAME}_NATIVE_REPOSITORIES`` and
    ``${PROJECT_NAME}_EXTRA_REPOSITORIES`` in the order that the project is
    being configured with.

TriBITS Package
+++++++++++++++

A TriBITS Package:

* Typically defines a set of libraries and/or headers and/or executables
  and/or tests with cmake build targets for building these and publishes the
  list of include directories and libraries that are created (along with CMake
  dependencies).
* Is declared in its parent repositories ``PackagesList.cmake`` file.
* Defines dependencies on upstream TPLs and/or other SE packages by just
  naming the dependencies in the file ``cmake/Dependencies.cmake`` using the
  macro `TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_.
* Can optionally have subpackages listed in the argument
  ``SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS`` to
  `TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_.
* Publishes a list of header file include paths and/or libraries for other
  downstream dependnet packages to use.

.. _TriBITS Package Core Files:

**TriBITS Package Core Files**

Teh core files that make up TriBITS Package (or Subpackage) are::

  ${${PACKAGE_NAME}_SOURCE_DIR}/
    CMakeLists.txt  # Only processed if the package is enabled
    cmake/
      Dependencies.cmake  # Always processed if the package is listed in the
                           # enclosing Repository

ToDo: Document each file in detail here!

**TriBITS Package Core Variables**

These are top-level (non-cache) CMake variables that are available once an SE
Package (Package or Subpackage) is processed.

  ``PACKAGE_NAME``

    The name of the current TriBITS SE package.  This is set automatically by
    TriBITS before the packages's ``CMakeLists.txt`` file is processed.

  ``PACKAGE_SOURCE_DIR``

    The absolute path to the base package's base source directory.  This is
    set automatically by TriBITS in the macro `TRIBITS_PACKAGE()`_.

  ``PACKAGE_BINARY_DIR``

    The absolute path to the base package's base binary/build directory.  This
    is set automatically by TriBITS in the macro `TRIBITS_PACKAGE()`_.

  ``${PACKAGE_NAME}_SOURCE_DIR``

    The absolute path to the base of a given package's source directory.
    CMake code, for example in other packages, refer to this by the raw name
    like ``PackageX_SOURCE_DIR``.  This makes such CMake code independent of
    where the package is in relation to other packages.  This variable is
    defined for all processed packages, independent of whether they are
    enabled.

  ``${PACKAGE_NAME}_BINARY_DIR``

    The absolute path to the base of a given package's binary directory.
    CMake code, for example in other packages, refer to this by the raw name
    like ``PackageX_BINARY_DIR``.  This makes such CMake code independent of
    where the package is in relation to other packages.  This variable is only
    defined if the package is enabled.

  ``PARENT_PACKAGE_SOURCE_DIR``

    The absolute path to the parent package's source directory.  This this
    only defined for a subpackage.

  ``PARENT_PACKAGE_BINARY_DIR``

    The absolute path to the parent package's binary directory.  This this
    only defined for a subpackage.

  ``${PACKAGE_NAME}_PARENT_REPOSITORY``

    The name of the package's parent repository.  This can be used by a
    package to access information about its parent repository.  For example,
    the variable ``${${PACKAGE_NAME}_PARENT_REPOSITORY}_SOURCE_DIR`` can be
    dereferenced.

  ``${PACKAGE_NAME}_ENABLE_TESTS``

    Set to ``ON`` if the package's tests are to be enabled.  This will enable
    a package's tests and all of its subpackage's tests.

  ``${PACKAGE_NAME}_ENABLE_EXAMPLES``

    Set to ``ON`` if the package's examples are to be enabled.  This will
    enable a package's examples and all of its subpackage's examples.

Currently, a Package can refer to its containing Repository and refer to its
source and binary directories.  This is so that it can refer to
repository-level resources (e.g. like the ``Trilinos_version.h`` file for
Trilinos packages).  This may be undesirable because it will make it very hard
to pull a package out of one Repository and place it in another repository for
a different use.  However, a package can indirectly refer to its own
repository without concern for what it is call by reading the variable
``${PACKAGE_NAME}_PARENT_REPOSITORY``.

TriBITS Subpackage
++++++++++++++++++

A TriBITS Subpackage:

* Typically defines a set of libraries and/or headers and/or executables
  and/or tests with cmake build targets for building these and publishes the
  list of include directories and libraries that are created (along with CMake
  dependencies).
* Is declared in its parent packages's ``Dependencies.cmake`` file in a call
  to `TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_ using the argument
  ``SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS``.
* Defines dependencies on upstream TPLs and/or other SE packages by just
  naming the dependencies in the file ``cmake/Dependencies.cmake`` using the
  macro `TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_.
* Publishes a list of header file include paths and/or libraries for other
  downstream dependnet packages to use.
* Can **NOT** have its own subpakcages defined (only top-level packages can
  have subpackages).

A subpakages set of core files is identical to the `TriBITS Package Core
Files`_ and the variables defined are almost identicial.  The only difference
is that a subpakcage may need to refer to its parent package where a top-level
package has not parent package.

TriBITS TPL
+++++++++++

A TriBITS TPL:

* Defines a set of pre-built libraries and/or header files and/or executables
  and/or some other resoruces used by one or more TriBITS Packages.
* Is listed as an explicit optional and required dependency in one or more
  package's ``Dependencies.cmake`` file using the macro
  `TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_.
* Publishes the include directories and/or libraries and/or executables
  provided by the TPL.

Using a TriBITS TPL is to be preferred over using a raw CMake
``FIND_PACKAGE(${TPL_NAME})`` because the TriBITS system guarantees that only
a single unique version of TPL will be used by multiple packages and by
declaring a TPL using TriBITS, automatical enable/disable logic will be
applied as described in `Package Enable/Disable Logic`_.

For each TPL referenced in a ``TPLsList.cmake`` file using the macro
`TRIBITS_DEFINE_REPOSITORY_TPLS()`_, there should exist a file, typically
called ``FindTPL${TPL_NAME}.cmake``, that once processed, produces the
variables ``${TPL_NAME}_LIBRARIES`` and ``${TPL_NAME}_INCLUDE_DIRS``.  A file
``FindTPL${TPL_NAME}.cmake`` is typically defined using the function
`TRIBITS_TPL_DECLARE_LIBRARIES()`_.  A simple example of such a file is the
standard ``FindTPLPETSC.cmake`` module which is:

.. include:: ../../tpls/FindTPLPETSC.cmake
   :literal:

However, the TriBITS system does not require the usage of of the function
`TRIBITS_TPL_DECLARE_LIBRARIES()`_ and does not even care about the TPL module
name ``FindTPL${TPL_NAME}.cmake``.  All that is required is that some CMake
file fragment exist that once included, will define the varaibles
``${TPL_NAME}_LIBRARIES`` and ``${TPL_NAME}_INCLUDE_DIRS``.  However, to be
user friendly, such a CMake file should respond to the same variables as
accepted by the standard `TRIBITS_TPL_DECLARE_LIBRARIES()`_ funtion.

The only core variables related to an enabled TPL ``${TPL_NAME}_LIBRARIES``
and ``${TPL_NAME}_INCLUDE_DIRS`` as defined in
`TRIBITS_TPL_DECLARE_LIBRARIES()`_ need to be defined.  For more details, see
`TRIBITS_DEFINE_REPOSITORY_TPLS()`_.

Processing of TriBITS Files: Ordering and Details
--------------------------------------------------

ToDo: Fill in!

Coexisting Projects, Repositories, and Packages
-----------------------------------------------

Certain simplifications area allowed when defining TriBITS projects,
repositories and packages.  The known allowed simplifications are described
below.

**TriBITS Repository == TriBITS Project**: It is allowed for a TriBITS Project
and a TriBITS Repository to be the same directory and have the same names with
no ``NativeRepositoriesList.cmake`` file.  In this case, the repository name
and the project name and the same name.  This is quite common.  This is the
case, for example, with the the Trilinos and the TribitsExampleProject
projects and repositories.  In this case, the Project's and the Repository's
``Version.cmake`` and ``Copyright.txt`` files are one and the same.

**TriBITS Package == TriBITS Repository**: It is also allowed for a TriBITS
Repository to have only one package and to have that package be the base
repository directory and have the same names.  The TriBITS Pepository and the
single TriBITS Package would typically have the same name in this case (but
that is actually not required but is confusing if one does not make it so).
For example, in the TriBITS test project ``MockTrilinos``, the repostory and
package ``extraRepoOnePackage`` are one in the same.  In this case, the file
``extraRepoOnePackage/PackagesList.cmake`` looks like:

.. include:: ../../package_arch/UnitTests/MockTrilinos/extraRepoOnePackage/PackagesList.cmake
   :literal:

This is used in the real TriBITS repository `DataTransferKit`_.

.. _DataTransferKit: https://github.com/CNERG/DataTransferKit

However, to maximize flexibility, it is recommended that a TriBITS package and
TriBITS repository not share the same directory.

**TriBITS Package, TriBITS Repository, TriBITS Package sharing the same source
directory**: In the extreme, it is posisble to collapase a single TriBITS
dackage, repository, and project into the same base source directory.
However, in this case, the TriBITS Project name and the TriBITS Package name
can not be the same and some modifications and tricks are needed to allow this
to work.  One example of this is the TriBITS project and TriBITS package
themselves, which are both rooted in the base ``tribits`` directory.  There
are a few restructions and modifications needed to get this to work:

* The Project and Package names can not be the same: In the case of TriBTS,
  the project name is ``TriBITSProj`` (as defined in
  ``tribits/ProjectName.cmake``) and the package name is ``TriBITS`` (as
  defined in ``tribits/PackagesList.cmake``).
* The base ``CMakeLists.txt`` file must me modified to allow it to be
  processed both as the base project file and as the package file: In the case
  of ``tribits/CMakeLists.txt``, a big if statement based in the local
  variable ``TriBITS_PROCESSED_BASED_PROJECT`` is used.
* The TriBITS package's binary directory must be modified.  Because of
  directory level targets like ``${PROJECT_NAME}_libs`` and
  ``${PACKAGE_NAME}_libs``, a subdirectory for package's the binary directory
  must be created.  This is simply done by overriding the binary directory
  name ``${PACKAGE_NAME}_SPECIFIED_BINARY_DIR``.  In the case of TriBITS, this
  is set to ``tribits`` in the ``tribits/PackagesList.cmake`` file.

Other than those modifications, a TriBITS project, repository, and package can
all be rooted in the same source directory.  However, to do so is a little
messing and is not recommended.  It was only done this way with the base
TriBITS directory in order to maintain backward compatibility.

However, one possible use case for collapsing a project, repository, and
package into a single base source directory would be to support the
stand-alone build of a TriBITS package as its own entity that uses and
installation of the TriBITS CMake modules.  If a given upstream TriBITS
package has no required upstream TriBITS package dependencies and minimal TPL
dependencies (or only uses "standard" TriBITS TPLs already defined in the
``tribits/tpls/`` directory), then creating a stand-alone project build of a
loan TriBITS package requires fairly little extra overhead.

Standard TriBITS TPLs
---------------------

ToDo: Fill in!


Example TriBITS Projects
=========================

In this section, a few different example TriBITS projects and packages are
previewed.  All of these examples exist in the TriBITS source tree directory
``tribits`` itself so they are available to all users of TriBITS.  These
examples also provide a means to test the TriBITS system itself.

The first example covered in detail is `TribitsExampleProject`_.  This first
example covers all the basics for setting up a simple multi-package TriBITS
project.  The second example outlined is `MockTrilinos` which mostly exists to
test the TirBITS system itself.  The last example mentioned is the the
``TriBITS`` package itself which allows the TriBITS system to be tested and
installed from any TriBITS project that lists it, including the
``TriBITSProj`` project itself (see `Coexisting Projects, Repositories, and
Packages`_).

The ``tribits/doc/examples/`` directory contains some other example TriBITS
projects and repositories as well that are refered to in this and other
documents.

TribitsExampleProject
----------------------

Here, the example project TribitsExampleProject that is provided with the
TriBITS directory tree is outlined.  From this simple example, one can quickly
see how the basic structural elements a TriBITS example project are put
together.

ToDo: Fill in!


MockTrilinos
-------------

ToDo: Fill in!


The TriBITS Package Itself
--------------------------

ToDo: Fill in!

Package Enable/Disable Logic
============================

ToDo: Fill in!

The following TriBITS repository-related variables alter what packages in a
given TriBITS repository get enabled implicitly or not:

  ``${REPOSITORY_NAME}_NO_IMPLICIT_PACKAGE_ENABLE``

    If set to ``ON``, then the packages in Repository ``${REPOSITORY_NAME}``
    will not be implicitly enabled in any of the package adjustment logic.

  ``${REPOSITORY_NAME}_NO_IMPLICIT_PACKAGE_ENABLE_EXCEPT``

    List of packages in the Repository ``${REPOSITORY_NAME}`` that will be
    allowed to be implicitly enabled.  Only checked if
    ``${REPOSITORY_NAME}_NO_IMPLICIT_PACKAGE_ENABLE`` is true.

The above variables typically are defined in the outer TriBITS Project's
``PackageName.cmake`` file in order to adjust how its listed repositories are
handled.

ToDo: Fill in!

Automated testing
=================

Much of the value provided by the TriBITS system is related to the support of
testing of a complex project.  Many different types of testing is required in
a complex project and development effort.  In addition a large project with
lots of repositories and packages provides a number of testing and developmetn
challanges but also provides a number of opportunities to do testing in an
efficient way; expecially pre-push and post-push continuous integration (CI)
testing.  In addition, a number of post-push automated nightly test cases must
be managed.  TriBITS takes full advantage of the features of raw CMake, CTest,
and CDash in support of testing and where gaps exist, TriBITS provides tools
and customizations.

The following subsections describe several aspects to the TriBITS support for
testing.  ToDo: outline the following subsections.


Testing categories for Repositories, Packages, and Tests
--------------------------------------------------------

ToDo: Define repo category Continuous, Nightly, and Experimental which also map to CDash tracks

ToDo: Define SE package test group PT, ST, and EX

ToDo: Define test category BASIC, CONTINUOUS, NIGHTLY, WEEKLY, and PERFORMANCE

ToDo: Discuss the propery usage of these test categories and why NIGHTLY
testing should be the default.

ToDo: Fill in!

Pre-push Testing using checkin-test.py
--------------------------------------

Describe the checkin-test.py script

TriBITS Package-by-Package CTest/Dash Driver
--------------------------------------------

ToDo: Fill in!

TriBITS CDash Customizations
---------------------------------

ToDo: Fill in!


CDash regression email addresses
++++++++++++++++++++++++++++++++++

Every TriBITS Package has a regression email address associated with it that
gets uploaded to a CDash project on a CDash server that is used to determine
what email address to use when a package has configure, build, or test
failures.  Because of the complex organizational nature of different projects
and different integration models, a single static email address for a given
package in every project build is not practical.

The TriBITS system allows for a Package's regression email to be specified in
the following order.:

1) REGRESSION_EMAIL_LIST (defined in Dependneices.cmake): Package-specific
email address overrides specified in the Dependencies.cmake file (the local
variable REGRESSION_EMAIL_LIST).  This package-specific email address will be
overridden if ${REPOSITORY_NAME}_REPOSITORY_OVERRIDE_PACKAGE_EMAIL_LIST using
the single Project or Repository master regression emails lists described
above.

2) ${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESSS_BASE: A base email address
specified at the Repository level creating package-specific email addresses
(e.g. <lower-case-package-name>-regression@some.repo.gov, where
${PROJECT_NAME}_REPOSITORY_EMAIL_URL_ADDRESSS_BASE=some.repo.gov).  This
variable is set in the Repositories cmake/DependenciesSetup.cmake file to
provide a default for the repository but can be overridden (i.e. by the
project).

3) ${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESSS: A Single email address
for all packages specified at the Repository level
(e.g. my-repo-regression@some.repo.gov).  This variable is set in the
Repositories cmake/DependenciesSetup.cmake file to provide a default for the
repository but can be overridden (i.e. by the project).

4) ${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESSS_BASE: A base email address
specified at the Project level creating package-specific email addresses
(e.g. <lower-case-package-name>-regression@some.project.gov, where
${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESSS_BASE=some.project.gov).  If not
already set, this variable will be set to
${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESSS_BASE for the first repostory
processed that has this set.

5) ${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESSS: A Single email address for
all packages specified at the Project level
(e.g. my-project-regression@some.project.gov).  If not already set, this
variable will be set to ${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESSS
for the first repostory processed that has this set.

If any of the email lists or URL string variables listed above are set to
"OFF" or "FALSE" (or some other value that CMake interprests as false) then
the varaibles are treated as empty and not set.

If a TriBITS project does not use CDash, then no email address needed to be
assigned to packages at all (which will be the case if none of the above
variables are set).

As a general rule, repository-level settings override project-level settings
and package-level settings override both.  Also, a project can redefine a
reposiotry's regression email list settings.

All of the email dependency managment logic must be accessable by just running
the macro:

    TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()

The above email address configuration variables are read from the Repository
and Project files RepositoryDependenciesSetup.cmake and
ProjectDependenciesSetup.cmake respectively.  The
RepositoryDependenciesSetup.cmake are read first in the specified order
followed up reading the ProjectDependenciesSetup.cmake file.  In this way, the
project can override any of the repository settings.

Here is the precedence order for how regression email addresses are selected
for a given package:

1) Package-specific email list is selected if defined (unless an override is
in place).

2) Repository-level option is selected over a project-level option.

3) Default email form with repository or project address base is selected over
single repository or project email address.

4) If none of the above are selected, then no email address is assigned.

What the above setup does is it results in the TriBITS system
creating a file called CDashSubprojectDependencies.xml that gets sent to
CDash. CDash then takes this file and creates, or updates, a set of CDash
users and sets up a mapping of Labels (which are used for TriBITS packages) to
CDash user emails addresses. CDash is automatically set up to process this XML
file and create and updates CDash users. It is not, however, set up to remove
labels from existing users.  Therefore, if you change a TriBITS package's
CDash regression email list (using one of the methods described above), then
you need to manually remove the associated labels from the old email address.
CDash will not remove them for you.

Therefore, to change the mapping of CDash regression email addresses to
TriBITS packages, you must perform the actions:

1) Change the TriBITS CMake files as described above that result in the
desired CDashSubprojectDependeinces.xml file. You can debug this by running
the checkin-test.py script and seeing what gets written in the
${PROJECT_NAME}Dependencies.xml file in the CHECKIN directory.

2) Log onto the CDash server using an administrator account and then remove
the auto-generated account for the CDash user email address for which labels
are being removed (i.e. no longer associated with a TriBITS package).  This is
needed since CDash seems to be unable to remove labels from an existing CDash
user (however this might be fixed in a current version of CDash).

3) The next time a CDash submit is performed by the
TribitsCTestDriverCore.cmake script, the CDash user associated with the mail
list with labels being removed will get automatically recreated with the right
list of labels (according to the current CDashSubprojectDependencies.xml
file).  Also, any new CDash users for new email addresses will be created.

Hopefully that should be enough clues to manage the mapping of CDash
regression email lists to TriBITS packages.


Multi-Repository Support
========================

ToDo: Discuss 'egdist', ExtraRepositoriesList.cmake, and the rep clone script.


Multi-Repository Almost Continuous Integration
----------------------------------------------

ToDo: Fill in!

Development Workflows
======================

Basic Development Workflow
--------------------------

ToDo: Fill in!

Multi-Repository Development Workflow
-------------------------------------

ToDo: Fill in!

Additional Topics
=================

In this section, a number of miscellaneous topics and TriBITS features are
discussed.  These features and topics are not considered the primary fetures
of TriBITS but can be very useful in many situations.

Project-Specific Build Quick Reference
--------------------------------------

If a project that uses TriBITS is going to have a significnat user base that
will configure, build, and test the project, then having some documentation
that explains how to do this would be useful.  For this purpose, TriBITS
provides a mechanism to quickly create a project-specific build quick
reference document in restructured text (RST) format and with HTML and
LaTeX/PDF outputs.  These documents are generally created in the base project
source tree and given then name ``<Project>BuildQuickRef.[rst,html,pdf]``.
This document consists of two parts.  One part is a generic template document::

  tribits/doc/build_quick_ref/TribitsBuildQuickRefBody.rst

provided in the TriBITS source tree that uses the place-holder ``<Project>``
for the for the real project name.  The second part is a project-specific
template file::

  <projectBaseDir>/cmake/<Project>BuildQuickRefTemplate.rst

which provides the outer RST doucment (with title, authors, abstract,
introduction, other introductory sections).  From these two files, the
script::

  tribits/doc/build_quick_ref/create-project-build-quickref.py

is used to replace ``<Project>`` in the ``TribitsBuildQuickRefBody.rst`` file
with the real project name (read from the project's ``ProjectName.cmake`` file
by default) and then generates the read-only files::

  <projectBaseDir>/
    <Project>BuildQuickRef.rst
    <Project>BuildQuickRef.html
    <Project>BuildQuickRef.pdf

To see a simple example of this, see::

  tribits/doc/examples/TribitsExampleProject/cmake/create-build-quickref.sh

A project-indepenent version of this file is provided in the
`TribitsBuildQuickRef.[rsts,html,pdf]
<../build_quick_ref/TribitsBuildQuickRef.html>`_ which is referred to many
times in this developers guide.


Creating Source Distributions
-----------------------------

ToDo: Fill in!


Regulated Backward Compatibility and Deprecated Code
----------------------------------------------------

ToDo: Fill in!

TriBITS Dashboard Driver
------------------------

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
