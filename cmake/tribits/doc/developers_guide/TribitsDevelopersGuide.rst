=======================================
TriBITS Developers Guide and Reference
=======================================

:Author: Roscoe A. Bartlett (bartlettra@ornl.gov)

:Abstract: This document describes the usage of TriBITS to build, test, and deploy complex software.  The primary audience are those individuals who develop on a software project which uses TriBITS.  The overall structure of a TriBITS project is descrided including all of the various project- and package-specific files that TriBITS requires or can use and how and what order these files are processed.  It also contains detailed reference information on all of the various TriBITS macros and functions directly used in TriBITS project CMkae files.  Many other topics of interest to a TriBITS project developer and archetect are discussed as well.

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


.. Common references to TribitsBuildQuickRef document

.. _Selecting the list of packages to enable: ../build_quick_ref/TribitsBuildQuickRef.html#selecting-the-list-of-packages-to-enable

.. _<Project>_EXTRAREPOS_FILE: ../build_quick_ref/TribitsBuildQuickRef.html#project-extrarepos-file

.. Common references to the TribitsOverview document

.. _TriBITS Overview: ../overview/TribitsOverview.pdf


.. Common references to the TribitsLifecycleModel document

.. _TriBITS Lifecycle Model: ../lifecycle_model/TribitsLifecycleModel.pdf

Introduction
=============

This document describes the usage of the TriBITS (Tribal Build, Integration,
Test System) to develop software projects.  An initial overview of TriBITS is
provided in the `TriBITS Overview`_ document which contains the big picture
and provides a high-level road map to what TriBITS provides.  This particular
document, however, describes the details on how to use the TriBITS system to
create a CMake build system for a set of compiled software packages.

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

The remainder of this documented is structured as follows.  First, there is
some additional `Background`_ material provided.  Then, a detailed
specification of `TriBITS Project Structure`_ is given.  This is followed up
by short descriptions of `Example TriBITS Projects`_ that are provided with
the TriBITS source tree that are used throughout this document.  The topic of
`Package Dependencies and Enable/Disable Logic`_ is then discussed. An
overview of the foundations for `TriBITS Automated Testing`_ is then given.
The topic of TriBITS `Multi-Repository Support`_ is examined next.
`Development Workflows`_ using TriBITS is then explored.  This is followed by
a set of detailed `Howtos`_.  Later some `Additional Topics`_ are presented
that don't fit well into other sections.  Then the main bulk of the detailed
reference material for TriBITS is given in the section `TriBITS Detailed
Reference Documentation`_.  Finally, several bits of information is provided
in the `Appendix`_.

Background
==========

Before diving into the details about TriBITS in the following sections, first
some background is in order.  First, a discussion of `TriBITS Developer and
User Roles`_ for TriBITS is provided to help the reader identify their own
roles and to help guide the reader to the appropriate documentation.  Then,
section `CMake Language Overview and Gotchas`_ tries to orient readers with
little to no CMake knowledge or experience on where to start and provide some
warnings about non-obvious CMake behavior that often trip up new users of
TriBITS.


TriBITS Developer and User Roles
--------------------------------

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
`Project-Specific Build Quick Reference`_.  Also, the `TriBITS Overview`_
document may be of some help also.  A TriBITS project user does not need to
know anything about the CMake langauge itself or any of the TriBITS macros or
functions described in `TriBITS Macros and Functions`_ or really anything else
described in this current document.

A **TriBITS Project Developer** is someone who contributes to a software
project that uses TriBITS.  They will add source files, libraries and
exectuables, test executables and define tests run with ``ctest``.  They have
to configure and build the project code in order to be able to develop and run
tests and therefore this role includes all of the necessary knowledge and
functions of a TriBITS Project User.  A casual TriBITS Project Developer
typically does not need to know a lot about CMake and really only need to know
a subset of the `TriBITS Macros and Functions`_ defined in this document in
addition to the genertic `TriBITS Build Quick Reference
<../build_quick_ref/TribitsBuildQuickRef.html>`_ document.  A slightly more
sophsiticated TriBITS Project Developer will also add new packages, add new
package dependencies, and define new TPLs.  This current TriBITS Developers
Guide and Reference document should supply everything such a developer needs
to know and more.  Only a smaller part of this document needs to be understood
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


CMake Language Overview and Gotchas
-----------------------------------

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
a single *TriBITS Project*, one or more `TriBITS Repositories`_ and one or
more `TriBITS Packages`_.  In addition, a TriBITS Package can be broken up
into `TriBITS Subpackages`_.  Together, the collection of TriBITS Packages and
TriBITS Subpackages are called *TriBITS Software Engineering Packages*, or
`TriBITS SE Packages`_ for short.

First, to establish the basic nomenclature, the key structural TriBITS units
are:

.. _TriBITS Packages:

* `TriBITS Package`_: A collection of related software that typically includes
  one or more source files built into one or more libraries and has assoicated
  tests to help define and protect the functionality provided by the software.
  A package also typically defines a unit of documentation and testing (see
  `TriBITS Automated Testing`_).  A TriBITS package may or may not be broken down into
  multiple subpackages. Examples of TriBITS packages in
  ``TribitsExampleProject`` include ``SimpleCXX``, ``MixedLanguage`` and
  ``PackageWithSubpackages``.  (Don't confuse a TriBTS "Package" with a raw
  CMake "Package" (see `History of TriBITS`_).  A raw CMake "Package" actually
  maps to a `TriBITS TPL`_.)

.. _TriBITS Subpackages:

* `TriBITS Subpackage`_: A part of a parent `TriBITS Package`_ that also
  typically has source files built into libraries and tests but is documented
  and tested along with the other subpackages the parent package.  The primary
  purpose for supporting subpackages to provide finer-grained of control
  software dependencies.  In ``TribitsExampleProject``,
  ``PackageWithSubpackages`` is an example of a package with subpackages
  ``SubpackageA``, ``SubpackaeB``, and ``SubpackageC``.  The full subpackage
  name has the parent package name prefixed the the subpackage name
  (e.g. ``PackageWithSubpackagesSubpackageA``).  The parent package is always
  implicitly dependent on its subpackages.

.. _TriBITS SE Package:

.. _TriBITS SE Packages:

* **TriBITS SE Package**: The combined set of `TriBITS Packages`_ and `TriBITS
  Subpackages`_ that constitute the basice *Softare Engineering* packages (see
  ???) of a TriBITS proejct: SE packages are the basis for setting
  dependencies in the TriBITS system.  For example, the SE Packages provided
  by the top-level example ``PackageWithSubpackages`` is (in order of
  increasing dependencies) ``PackageWithSubpackagesSubpackageA``,
  ``PackageWithSubpackagesSubpackaeB``, ``PackageWithSubpackagesSubpackageC``,
  and ``PackageWithSubpackages`` (see `TribitsExampleProject`_).

.. _TriBITS TPLs:

* `TriBITS TPL`_: The specification for a particular external dependency that
  is required or can be used in one or more `TriBITS SE Packages`_.  A TPL (a
  Third Party Library) typically provides a list of libraries or a list
  include directories for header files but can also be manifisted in order
  ways as well.  Examples of basic TPLs include ``BLAS``, ``LAPACK``, and
  ``Boost``.

.. _TriBITS Repositories:

* `TriBITS Repository`_: A collection of one or more `TriBITS Packages`_
  specified in a `<repoDir>/PackagesList.cmake`_ file.

.. _TriBITS Projects:

* `TriBITS Project`_: A collection of `TriBITS Repositories`_ and `TriBITS
  Packages`_ that defines a CMake ``PROJECT`` defining soiftware which can be
  directly configured with ``cmake`` and then built, tested, and installed.

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

A TriBITS Project:

* Defines a complete CMake project and calls ``PROJECT(${PROJECT_NAME} ...)``.
* Consists of one or more TriBITS Repositories (see `TriBITS Repository`_).
* Defines the ``PROJECT_NAME`` CMake variable (defined in
  `<projectDir>/ProjectName.cmake`_)
* Defines a set of native Repositories (see below) that define packages and
  TPLs.
* Allows for extra Repositories to be added on before or after the set of
  native Repositories (specified in
  `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ or by CMake variables)
* Defines a default CDash server and default project name on the server (the
  project name on the CDash server must be the same as ``${PROJECT_NAME}``).

For more details on the definition of a TriBITS Project, see:

* `TriBITS Project Core Files`_
* `TriBITS Project Core Variables`_

TriBITS Project Core Files
..........................

The core files making up a TriBITS Project (where ``<projectDir> =
${PROJECT_SOURCE_DIR}``) are::

  <projectDir>/
     ProjectName.cmake    # Defines PACAKGE_NAME
     CMakeLists.txt       # Base project CMakeLists.txt file
     CTestConfig.cmake    # [Optional] Needed for CDash submits
     Version.cmake        # [Optional] Dev mode, Project version, VC branch
     cmake/
       NativeRepositoriesList.cmake    # [Optional] Used for some meta-projects
       ExtraRepositoriesList.cmake     # [Optional] Lists repos and VC URLs 
       ProjectDependenciesSetup.cmake  # [Optional] Project deps overrides
       CallbackDefineProjectPackaging.cmake  # [Optional] CPack settings
       tribits/    # [Optional] Or provide ${PROJECT_NAME}_TRIBITS_DIR
       ctest/
         CTestCustom.cmake.in  # [Optional] Custom ctest settings

These TriBITS Project files are documented in more detail below:

* `<projectDir>/ProjectName.cmake`_
* `<projectDir>/CMakeLists.txt`_
* `<projectDir>/CTestConfig.cmake`_
* `<projectDir>/Version.cmake`_
* `<projectDir>/cmake/NativeRepositoriesList.cmake`_
* `<projectDir>/cmake/ExtraRepositoriesList.cmake`_
* `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_
* `<projectDir>/cmake/CallbackDefineProjectPackaging.cmake`_
* `<projectDir>/cmake/tribits/`_
* `<projectDir>/cmake/ctest/CTestCustom.cmake.in`_

.. _<projectDir>/ProjectName.cmake:

**<projectDir>/ProjectName.cmake**: [Required] At a minimum provides a
``SET()`` statement to set the local variable ``PROJECT_NAME``.  This file is
the first file is read by a number of tools in order to get the TriBITS
project's name.  This file is read first in every context that involves
processing the TriBITS project's files, including processes and tools that
just need to build a package and TPL dependnecy tree (see `Package
Dependencies and Enable/Disable Logic`_).  Being this is the first file read
in for a TriBITS project and that it is read in first at the top level scope
in every context, this is a good file to put other universal static project
options in that affect dependency handling.  Note that this is a project, not
a repository file so no general repository-specific settings should go in this
file.  A simple example of this file is
`TribitsExampleProject`_/``PackageName.cmake``:

.. include:: ../examples/TribitsExampleProject/ProjectName.cmake
   :literal:


.. _<projectDir>/CMakeLists.txt:

**<projectDir>/CMakeLists.txt**: [Required] The top-level CMake project file.
This is the first file that the ``cmake`` exectuable processes that starts
everything off and is the base level scope for local CMake variables. Due to a
number of CMake limitations and quarks, a project's top-level
``CMakeLists.txt`` file is not a clean as one might otherwise hope would be
but it is not too bad.  A simple, but representative, example is
`TribitsExampleProject`_/``CMakeLists.txt``:

.. include:: ../examples/TribitsExampleProject/CMakeLists.txt
   :literal:

A couple of CMake and TriBITS quarks that that above example
``CMakeLists.txt`` addresses are worth some discussion.  First, to avoid
duplication, the project's ``ProjectName.cmake`` file is read in with an
``INCLUDE()`` that defines the local variable ``PROJECT_NAME``.  Right after
this initial include, the built-in CMake command ``PROJECT(${PROJECT_NAME}
NONE)`` is run.  This command must be explicitly called with ``NONE`` so as to
avoid default CMake behavior for defining compilers.  The definition of
compilers comes later as part of the TriBITS system inside of the
`TRIBITS_PROJECT()`_ command (see `Full Processing of TriBITS Project
Files`_).

As noted in the above example file, the only project defaults that should be
set in this top-level ``CMakeLists.txt`` file are those that do not impact the
list of package enables/disables.  The latter type of defaults should set in
other files (see below).

In this example project, a CMake cache variable
``${PROJECT_NAME}_TRIBITS_DIR`` must be set by the user to define where the
base ``tribits`` source directory is located.  With this variable set
(i.e. passed into ``cmake`` command-line use
``-DTribitsExProj_TRIBITS_DIR=<someDir>``), one just includes a single file to
pull in the TriBITS system::

  INCLUDE("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")

With the ``TriBITS.cmake`` file included, the configuration of the project
using TriBITS occurs with a single call to `TRIBITS_PROJECT()`_.

Some projects, like Trilinos, actually snapshot the ``tribits`` directory into
their source tree and therefore don't need to have this variable set.  In
Trilinos, the include line is just::

  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits/TriBITS.cmake)

The minimum CMake version must also be delcared in the top-level
``CMakeLists.txt`` file as shown.  Explicitly setting the minimum CMake
version avoids strange errors that can occur when someone tries to build the
project using a version of CMake that is too old.  If the given project
requires a version of CMake newer than what is required by TriBITS itelf (as
defined in the variable ``TRIBITS_CMAKE_MINIMUM_REQUIRED`` which was set when
the ``TriBITS.cmake`` file was included), then that version can be passed
instead of using ``${TRIBITS_CMAKE_MINIMUM_REQUIRED}`` (the current minimum
version of CMake required by TriBITS is given at in `TribitsBuildQuickRef
<../build_quick_ref/TribitsBuildQuickRef.html#getting-set-up-to-use-cmake>`_)
.  For example, the ``VERA/CMakeLists.txt`` file lists as its first line::

  SET(VERA_TRIBITS_CMAKE_MINIMUM_REQUIRED 2.8.5)
  CMAKE_MINIMUM_REQUIRED(VERSION ${VERA_TRIBITS_CMAKE_MINIMUM_REQUIRED})

.. _<projectDir>/CTestConfig.cmake:

**<projectDir>/CTestConfig.cmake**: [Optional] Specifies the CDash site and
project to submit results to when doing an automated build driven by CTest (it
is read by the scripting code in ``TribitsCTestDriverCore.cmake``).  This file
only needs to be present when the CTest driver code in
``TribitsCTestDriverCore.cmake`` is used.  This file is also required to use
the TriBITS-generated ``dashboard`` target (see `Dashboard Submissions
<../build_quick_ref/TribitsBuildQuickRef.html#dashboard-submissions>`_).  An
example of this file is `TribitsExampleProject`_/``CTestConfig.cmake``:

.. include:: ../examples/TribitsExampleProject/CTestConfig.cmake
   :literal:

All of the varibles set in this file are directly understood by raw ``ctest``
and will not be explained here further (see documentation for the standard
CMake module ``CTest``).  The usage of the function
`SET_DEFAULT_AND_FROM_ENV()`_ allows the varaibles to be overridded both as
CMake cache variables and in the environment.  The latter is needed when
running using ``ctest`` as the driver.  Given that all of these variables are
nicely namespaced, overriding them on the shell environment is not as
dangerous as might otherwise be the case but this is what had to be done to
get around limitations for older versions of CMake/CTest.

.. _<projectDir>/Version.cmake:

**<projectDir>/Version.cmake**: If defined, gives the project's version and
determines development/release mode (see `Project and Repositiory Versioning
and Release Mode`_).  This file is read in (using ``INCLUDE()``) in the
project's base-level ``<projectDir>/CMakeLists.txt`` file scope so local
variables set in this file are seen by the entire CMake project.  For example,
`TribitsExampleProject`_/``Version.cmake``, this looks like:

.. include:: ../examples/TribitsExampleProject/Version.cmake
   :literal:

Note that the prefix ``${REPOSITORY_NAME}_`` is used instead of hard-coding
the project name.  This is so that the same ``Version.txt`` file can be used
as the the `<repoDir>/Version.cmake`_ file and have the repository name be
flexible.  TriBITS sets ``REPOSITORY_NAME = ${PROJECT_NAME}`` when it reads in
this file at the project-level scope.

It is strongly recommended that every TriBITS project contain a
``Version.cmake`` file, even if a release has never occured.  Otherwise, the
project needs to define the variable
``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT`` at the global project
scope (perhaps in ``<projectDir>/ProjectName.cmake``) to get right development
mode of behavior (see `${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE`_).

.. _<projectDir>/cmake/NativeRepositoriesList.cmake:

**<projectDir>/cmake/NativeRepositoriesList.cmake**: [Optional] If present,
this file gives the list of native repositories to this TriBITS project.  The
file must contain a ``SET()`` statement defining the variable
``${PROJECT_NAME}_NATIVE_REPOSITORES`` which is just a flat list of repository
names that must also be directory names under ``<projectDir>/``.  For example,
if this file contains::

  SET(${PROJECT_NAME}_NATIVE_REPOSITORES Repo0 Repo1)

then the directories ``<projectDir>/Repo0`` and ``<projectDir>/Repo1`` must
exist and must be valid TriBITS repositories (see `TriBITS Repository`).  

There are no examples for the usage of this file in any of the TriBITS
examples or test projects.  However, support for this file is maintained for
backward compatibility since there are some TriBITS projects that use it.  It
is recommended instead to define multiple repositories using the
`<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file as it allows for more
flexibility in how extra repositories are specified and how they are accessed.

If this file ``NativeRepositoriesList.cmake`` does not exist, then TriBITS
sets ```${PROJECT_NAME}_NATIVE_REPOSITORES`` eqaual to ".", or the base
project directory (i.e. ``<projectDir>/.``).  In this case, the file
``<projectDir>/PackagesList.cmake`` and ``<projectDir>/TPLsList.cmake`` must
exist.  However, if the project has no native packages or TPLs, then these
files be set up to list no packages or TPLs.  This is the case for
meta-projects like VERA that have only extra repositories specified in the
file `<projectDir>/cmake/ExtraRepositoriesList.cmake`_.

.. _<projectDir>/cmake/ExtraRepositoriesList.cmake:

**<projectDir>/cmake/ExtraRepositoriesList.cmake**: [Optional] If present,
this file defines a list of extra repositories that are added on to the
project's native repositories.  The list of repositories is defined using the
macro `TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES()`_.  For example, the extra
repos file:

.. include:: ../../python/UnitTests/ExtraReposList.cmake
   :literal:

shows the speification of both TriBITS Repositories and non-TriBITS VC
Repositories.  In the above file, the repositories ``ExtraRepo1``,
``ExtraRepo3``, and ``ExtraRepo4`` are VC repositories that are cloned into
directories under ``<projectDir>`` of the same names from the URLs
``someurl.com:/ExtraRepo1``, ``someurl3.com:/ExtraRepo3``, and
``someurl4.com:/ExtraRepo4``, respectively.  However, the repository
``ExtraRepo2`` is **not** a `TriBITS Repository`_ because it is marked as
``NOPACKAGES``.  In this case, it gets cloned as the directory::

  <projectDir>/packages/SomePackage/Blah

However, the code in the tools ``checkin-test.py`` and
``TribitsCTetsDriverCore.cmake`` will consider this repository and directory
and any changes to this repository will be listed as changes to
``somePackage``.

.. _<projectDir>/cmake/ProjectDependenciesSetup.cmake:

**<projectDir>/cmake/ProjectDependenciesSetup.cmake**: [Optional] If present,
this file is included a single time as part of the generation of the project's
dependency data-structure (see `Reduced Package Dependency Processing`_).  It
gets included at the top project level scope after all of the
`<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ files have been included
but before all of the package `<packageDir>/cmake/Dependencies.cmake`_ files
are included.  Any local variables set in this file have project-wide scope.
The primary purpose for this file is to set variables that will impact the
processing of project's package's ``Dependencies.cmake`` files.

The typical usage of this file is to set the default CDash email base address
that will be the default for all of the defined packages (see `CDash
regression email addresses`_).  For example, to set the default email address
for all of the packages, one would set in this file::

   SET_DEFAULT(${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESSS
       projectx-regressions@somemailserver.org)

The repository email address varaibles
``${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESSS_BASE`` and
``${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESSS`` possibly set in the
just processed `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ files can
also be overridded here.  The CASL VERA meta-project uses this file to
override several of the repository-specific email addresses for its
constituent CDash email addresses.

In general, variables that affect how package dependnecies are defined or
affect package and TPL enable/disable logic should be defined in this file.

.. _<projectDir>/cmake/CallbackDefineProjectPackaging.cmake:

**<projectDir>/cmake/CallbackDefineProjectPackaging.cmake**: [Optional] If
exists, defines the CPack settings for the project (see `Offical CPack
Documentation <http://www.cmake.org/cmake/help/documentation.html>`_ and
`Online CPack Wiki <http://www.cmake.org/Wiki/CMake:Packaging_With_CPack>`_).
This file must define a macro called ``TRIBITS_PROJECT_DEFINE_PACKAGING()``
which is then invoked by TriBITS.  The file:

  `TribitsExampleProject`_/``cmake/CallbackDefineProjectPackaging.cmake``

provides a good example which is:

.. include:: ../examples/TribitsExampleProject/cmake/CallbackDefineProjectPackaging.cmake
   :literal:

The CPack variables that should be defined at the project-level should be
described in the `Offical CPack Documentation`_.

Settings that are general for all distributions (like non-package repository
files to exclude from the tarball) should be set at the in the file
`<repoDir>/cmake/CallbackDefineRepositoryPackaging.cmake`_.

.. _<projectDir>/cmake/tribits/:

**<projectDir>/cmake/tribits/**: [Optional] The typical location of the
``tribits`` soruce tree for projects that choose to snapshot or checkout
TriBITS into their source tree.  Trilinos, for example, currently snapshots
the TriBITS source tree into this directory.

.. _<projectDir>/cmake/ctest/CTestCustom.cmake.in:

**<projectDir>/cmake/ctest/CTestCustom.cmake.in**: [Optional] If this file
exists, it is processed using a ``CONFIGURE_FILE()`` command to write the file
``CTestCustom.cmake`` in the project base build tree.  This file is picked up
automatically by ``ctest`` (see `CTest documentation
<http://www.cmake.org/Wiki/CMake/Testing_With_CTest>`_).  This file is
typically used to change the maximum size of test output.  For example, the
`TribitsExampleProject`_/cmake/ctest/CTestCustom.cmake.in`` looks like:

.. include:: ../examples/TribitsExampleProject/cmake/ctest/CTestCustom.cmake.in
   :literal:

which sets the output size for each test submitted to CDash be unlimited
(which is not really recommended).  These variables used by Trilinos at one
time were::

  SET(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE 50000)
  SET(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE 5000000)

which sets the max output for passed and failed tests to 50000k and 5000000k,
respectively.

For documentation of the options one can change for CTest, see `This Online
Wiki Page <http://www.cmake.org/Wiki/CMake/Testing_With_CTest>`_.


TriBITS Project Core Variables
..............................

The following local variables are defined in the top-level Project
``CMakeLists.txt`` file scople and are therefore acceessible by all files
prcessed by TriBITS:

  ``PROJECT_NAME``

    The name of the TriBTS Project.  This exists to support, among other
    things, the ability for subordinate units (Repositories and Packages) to
    determine the Project in which is participating.  This is typically read
    from a ``SET()`` statement in the project's
    `<projectDir>/ProjectName.cmake`_ file.

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

The following internal project-scope local CMake variables are defined by
TriBITS for the project's repositoreis.:

  ``${PROJECT_NAME}_NATIVE_REPOSITORIES``

     The list of Native Repositories for a given TriBITS project
     (i.e. Repositories that are always present when configuring the Project
     and are managed in the same VC repo typically).  If the file
     ``${PROJECT_SOURCE_DIR}/NativeRepositoriesList.cmake`` exists, then the
     list of native repositories will be read from that file.  If the file
     ``NativeRepositoriesList.cmake`` does not exist, then the project is
     assumed to also be a repository and the list of native repositories is
     just the local project directory ``${PROJECT_SOURCE_DIR}/.``.  In this
     case, the ``${PROJECT_SOURCE_DIR}/`` must contain at a minumum a
     ``PackagesList.cmake``, and a ``TPLsList.cmake`` file (see `TriBITS
     Repository`_).

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


TriBITS Repository
++++++++++++++++++

A TriBITS Repository is the basic unit of ready-made composition between
different collections of software that use the TriBITS CMake build and system.

In short, a TriBITS Repository:

* Is a named collection of related TriBITS Packages and TPLs (defined in
  `<repoDir>/PackagesList.cmake`_ and `<repoDir>/TPLsList.cmake`_
  respectively)
* Defines the base source and binary directories for the Repository
  ``${REPOSITORY_NAME}_SORUCE_DIR`` and ``${REPOSITORY_NAME}_BINARY_DIR``.
* [Optional] Defines a common set of initializations and other hooks for a all
  the packages in the repository.

For more details on the definition of a TriBITS Repository, see:

* `TriBITS Repository Core Files`_
* `TriBITS Repository Core Variables`_

TriBITS Repository Core Files
.............................

The core files making up a TriBITS Repository (where ``<reposDir> =
${${REPOSITORY_NAME}_SOURCE_DIR}``) are::

  <repoDir>/
    PackagesList.cmake
    TPLsList.cmake
    Copyright.txt  # [Optional] Only needed if creating version header file
    Version.cmake  # [Optional] Info inserted into ${REPO_NAME}_version.h
    cmake/
       RepositoryDependenciesSetup.cmake # [Optional]
       CallbackSetupExtraOptions.cmake # [Optional] Called after tribits options
       CallbackDefineRepositoryPackaging.cmake # [Optioinal] CPack packaging

These TriBITS Repository files are documented in more detail below:

* `<repoDir>/PackagesList.cmake`_
* `<repoDir>/TPLsList.cmake`_
* `<repoDir>/Copyright.txt`_
* `<repoDir>/Version.cmake`_
* `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_
* `<repoDir>/cmake/CallbackSetupExtraOptions.cmake`_
* `<repoDir>/cmake/CallbackDefineRepositoryPackaging.cmake`_

.. _<repoDir>/PackagesList.cmake:

**<repoDir>/PackagesList.cmake**: [Required] Provides the list of top-level
packages defined by the repository.  This file typically just calls the macro
`TRIBITS_DEFINE_REPOSITORY_PACKAGES()`_ to define the list of packages along
with their directories and other properties.  For example, the file
`TribitsExampleProject`_/``PackagesList.cmake`` looks like:

.. include:: ../examples/TribitsExampleProject/PackagesList.cmake
   :literal:

Other comamnds that are appropriate to use in this file include
`TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS()`_.  Also, if the binary directory for
any package ``<packageName>`` needs to be changed from the default, then the
variable ``<packageName>_SPECIFIED_BINARY_DIR`` can be set.  One can see an
example of this in the file ``tribits/PackageList.cmake`` which shows

.. include:: ../../PackagesList.cmake
   :literal:

(see `TriBITS Package, TriBITS Repository, TriBITS Package sharing the same
source directory`_).

It is perfectly legal for a TriBITS repjository to define no packages at all
with::

  TRIBITS_DEFINE_REPOSITORY_PACKAGES()

and this would be the case for a TriBITS meta-project that has no native
packages, only extra repositories.

.. ToDo: Point to example meta-project.

.. _<repoDir>/TPLsList.cmake:

**<repoDir>/TPLsList.cmake**: [Required] Provides the list of TPLs that are
listed as TPLs in the repository's SE packages
`<packageDir>/cmake/Dependencies.cmake`_ files (see `TriBITS TPL`).  This file
typically just calls the macro `TRIBITS_DEFINE_REPOSITORY_TPLS()`_ to define
the TPLs along with their find modules and other properties.  See an example
from `TribitsExampleProject`_/``TPLsList.cmake`` which shows:

.. include:: ../examples/TribitsExampleProject/TPLsList.cmake
   :literal:

It is perfectly fine to specify no TPLs for a repository with::

  TRIBITS_DEFINE_REPOSITORY_TPLS()

but the macro ``TRIBITS_DEFINE_REPOSITORY_TPLS()`` has to be called, even if
there are no TPLs.  See `TRIBITS_DEFINE_REPOSITORY_TPLS()`_ for further
details.

.. _<repoDir>/Copyright.txt:

**<repoDir>/Copyright.txt**: [Optional] Gives the default copyright and
license declaration for all of the software in the TriBITS repository
``<repoDir>``.  This file is read into a string and then used to configure the
repository's version file (see `Project and Repositiory Versioning and Release
Mode`_).

.. _<repoDir>/Version.cmake:

**<repoDir>/Version.cmake**: Contains version information for the repository
(and the project also if this is also the base project).  For example,
`TribitsExampleProject`_/``Version.cmake``, this looks like:

.. include:: ../examples/TribitsExampleProject/Version.cmake
   :literal:

Note that the prefix ``${REPOSITORY_NAME}_`` is used instead of hard-coding
the repository's name to allow flexibility in what a meta-project names a
given TriBITS repository.

The local variables in these set statements are processed in the base project
directory's local scope and are therefore seen by the entire CMake project.
When this file is read in repository mode, the variable
``${REPOSITORY_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT`` is ignored.

.. _<repoDir>/cmake/RepositoryDependenciesSetup.cmake:

**<repoDir>/cmake/RepositoryDependenciesSetup.cmake**: [Optional] If present,
this file is included a single time as part of the generation of the project
dependency data-structure (see `Reduced Package Dependency Processing`_).  It
gets included at the top project level scope in order with the other
repositories listed in ``${PROJECT_NAME}_ALL_REPOSITORIES``.  Any local
variables set in this file have project-wide scope.  The primary purpose for
this file is to set variables that will impact the processing of project's
package's ``Dependencies.cmake`` files and take care of other enable/disable
issues that are not clearly handled by the TriBITS system automatically.

The typical usage of this file is to set the default CDash email base address
that will be the default for all of the defined packages (see `CDash
regression email addresses`_).  For example, to set the default email address
for all of the packages in this repository, one would set in this file::

   SET_DEFAULT(${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESSS
       repox-regressions@somemailserver.org)

.. _<repoDir>/cmake/CallbackSetupExtraOptions.cmake:

**<repoDir>/cmake/CallbackSetupExtraOptions.cmake**: [Optional] If defined,
this file is processed (included) for each repo in order right after the basic
TriBITS options are defined in the macro
``TRIBITS_DEFINE_GLOBAL_OPTIONS_AND_DEFINE_EXTRA_REPOS()``.  This file must
define the macro ``TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS()`` which is then
called by the TriBITS system.  This file is only processed when doing a basic
configuration of the project and **not** when it is just building up the
dependency data-structures as part of other tools (like ``checkin-test.py``,
``TribitsCTestDriverCore.cmake``, etc.).  Any local variables set in this file
or macro have project-wide scope.

A few additional variables are defined by the time this file is procesed and
can be ued in the logic in these files.  The variables that should already be
defined, in addition to all of the basic user TriBITS cache variables, include
``CMAKE_HOST_SYSTEM_NAME``, ``${PROJECT_NAME}_HOSTNAME``, and
``PYTHON_EXECUTABLE``.  The types of logic to put in this file includes:

* Setting additional user cache variable options that are used by multiple
  packages into the TriBITS Repository.  For example, Trilinos defines a
  ``Trilinos_DATA_DIR`` user cache varaible that several Trilinos packages
  use.
* Disabling packages in the TriBITS Repository when conditions will not allow
  them to be enabled.  For example, Trilinos disables the package ForTrilinos
  when Fortran is disabled and disables PyTrilinos when Python support is
  disabled.

An example of this file is:

  `TribitsExampleProject`_/``/cmake/CallbackSetupExtraOptions.cmake``

which currently looks like:

.. include:: ../examples/TribitsExampleProject/cmake/CallbackSetupExtraOptions.cmake
   :literal:

.. _<repoDir>/cmake/CallbackDefineRepositoryPackaging.cmake:

**<repoDir>/cmake/CallbackDefineRepositoryPackaging.cmake**: [Optional] If
this file exists, then this file defines extra CPack-related options that are
specific to the TriBITS Repository.  This file must define the macro
``TRIBITS_REPOSITORY_DEFINE_PACKAGING()`` which is called by TriBITS.  This
file is procdessed as the top project-level scope so any local variables set
have project-wide effect.  This file is processed before the project's
`<projectDir>/cmake/CallbackDefineProjectPackaging.cmake`_ so any options
defined in the repositories file are overridden by the project.  This file
typically this just involves setting extra excludes to remove files from the
tarball.  The file:

  `TribitsExampleProject`_/``cmake/CallbackDefineRepositoryPackaging.cmake``

provides a good example which is:

.. include:: ../examples/TribitsExampleProject/cmake/CallbackDefineRepositoryPackaging.cmake
   :literal:


TriBITS Repository Core Variables
.................................

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

TriBITS Package
+++++++++++++++

A TriBITS Package:

* Typically defines a set of libraries and/or header files and/or executables
  and/or tests with CMake build targets for building these and exports the
  list of include directories, libraries, and targets that are created (along
  with CMake dependencies).
* Is declared in its parent repository's `<repoDir>/PackagesList.cmake`_ file.
* Defines dependencies on upstream TPLs and/or other SE packages by just
  naming the dependencies in the file
  `<packageDir>/cmake/Dependencies.cmake`_..
* Can optionally have subpackages listed in the argument
  ``SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS`` to
  `TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_.
* Is the fundamental unit of software partitioning and aggregation and must
  have a unique packae name that is globally unique, not only within its
  defined repository, but cross all possible repositories that might be
  coupled together some day using TriBITS.
* Is the unit of testing as drived by ``TribitsCTestDriverCore.cmake`` and
  displayed on CDash.

**WARNING:** As noted above, one must be very careful to pick package names
that will be globally unique

For more details on the definition of a TriBITS Package (or subpackage), see:

* `TriBITS Package Core Files`_
* `TriBITS Package Core Variables`_

TriBITS Package Core Files
..........................

The core files that make up TriBITS Package (where ``<packageDir> =
${${PACKAGE_NAME}_SOURCE_DIR}``) are::

  <packageDir>/
    CMakeLists.txt  # Only processed if the package is enabled
    cmake/
      Dependencies.cmake  # Always processed if the package is listed in the
                           # enclosing Repository

**NOTE:** Before a TriBITS Package's files are described in more detail, it is
important to note that all of the package's files that define its behavior and
tests should strictly be contained under the package's base source directory
``<packageDir>/`` if at all possible.  While this is not a requirement for the
basic TriBITS build system, the approach for automatically detecting when a
package has changed by looking at what files have changed (which is used in
``checkin-test.py`` and ``TribitsCTestDriverCore.cmake``) requires that the
package's files be listed under ``<packageDir>/`` (see `Pre-push Testing using
checkin-test.py`_).  Without this, these development tetsing tools will not be
able to effectively determine what needs to be rebuilt and retested which can
clean to pushing broken software.

These TriBITS Package files are documented in more detail below:

* `<packageDir>/cmake/Dependencies.cmake`_
* `<packageDir>/CMakeLists.txt`_

.. _<packageDir>/cmake/Dependencies.cmake:

**<packageDir>/cmake/Dependencies.cmake**: [Required] Defines the dependencies
of a given SE package using the macro
`TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_.  This file is processed at the
top-level project scope (using an ``INCLUDE()``) so any local variables set
will be seen by the entire project.  This file is always processed, including
when just building the project's `<Project>PackageDependencies.xml`_ file.

An example of a ``Dependencies.cmake`` file for a package with optional and
required dependencies is for the mock ``Panzer`` package in `MockTrilinos`_:

.. include:: ../../package_arch/UnitTests/MockTrilinos/packages/panzer/cmake/Dependencies.cmake
   :literal:

An example of a package with subpackages is ``PackageWithSubpackages`` in
`TribitsExampleProject`_ with the ``Dependencies.cmake`` file:

.. _package_with_subpackages/cmake/Dependencies.cmake:

.. include:: ../examples/TribitsExampleProject/packages/package_with_subpackages/cmake/Dependencies.cmake
   :literal:

The last case defines three subpackages which creates three new SE packages
with names ``PackageWithSubpackagesSubpackageA``,
``PackageWithSubpackagesSubpackageB``, and
``PackageWithSubpackagesSubpackageC``.

.. _<packageDir>/CMakeLists.txt:

**<packageDir>/CMakeLists.txt**: [Required] The package's top-level
``CMakeLists.txt`` file that defines the libraries, include directories, and
contains the tests for the package.

The basic structure of this file for a **package without subpackages** is
shown in:

  `TribitsExampleProject`_/``packages/simple_cxx/CMakeLists.txt``

which is:

.. include:: ../examples/TribitsExampleProject/packages/simple_cxx/CMakeLists.txt
   :literal:

The first command at the top of the file in `TRIBITS_PACKAGE()`_ where the
package name is passed in in addition to a few other options.  TriBITS
obviously already knows the package name.  The purpose for repeating it is as
documentation for the developer's sake.  Then a set of platform confgiure-time
tests is typically performed (if there are any).  In this example, the
existance of the C++ ``__int64`` data-type is checked using the module
`CheckFor__int64.cmake`` (which is in the ``cmake/`` directory of this
package.  (CMake has great support for configure-time tests, see
`Configure-time System Tests`_.)  This is followed by package-specific
options.  In this case, the standard TriBITS options for debug checking and
deprecated warnings are added using the standard macros
`TRIBITS_ADD_DEBUG_OPTION()`_ and
`TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()`_.  After all of this up front
stuff (which will be present in any moderately complex CMake-configured
project) the source and the test subdirectories are added that actually define
the library and the tests.  In this case, the standard
`TRIBITS_ADD_TEST_DIRECTORIES()`_ macro is used which only conditionally adds
the tests for the package.

The final command in the package's ``CMakeLists.txt`` file must always be
`TRIBITS_PACKAGE_POSTPROCESS()`_.  This is needed in order to perform some
necessary post-processing by TriBITS.

It is also possible to for the package's top-level ``CMakeLists.txt`` to be
the only such file in a package.  Such an example can be see in the example
project `TribitsHelloWorld`_.

When a TriBITS package is broken up into subpackages (see `TriBITS
Subpackage`_), its ``CMakeLists.txt`` file looks a little different.  The
basic structure of this file for a **package with subpackages** is shown in:

  `TribitsExampleProject`_/``packages/package_with_subpackages/CMakeLists.txt``

which contains:

.. include:: ../examples/TribitsExampleProject/packages/package_with_subpackages/CMakeLists.txt
   :literal:

What is different about ``CMakeLists.txt`` files for packages without
subpackages is that the `TRIBITS_PACKAGE()`_ command is broken up into two
parts `TRIBITS_PACKAGE_DECL()`_ and `TRIBITS_PACKAGE_DEF()`_.  In between
these two comamnds, the parent package can define the common package options
and then calls the commnand `TRIBITS_PROCESS_SUBPACKAGES()`_ which fully
processes the packages.  If the parent package has libraries and/or
tests/example of its own, it can define those after calling
`TRIBITS_PACKAGE_DEF()`_, just like with a regular package.  However, it is
rare for a package broken up into subpackages to have its own libraries and/or
tests and examples.  As always, the final command called inside of the
``CMakeLists.txt`` is `TRIBITS_PACKAGE_POSTPROCESS()`_.

NOTE: The package's ``CMakeLists.txt`` file only gets processed if the package
is actually enabled with ``${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=ON``.  This
is an important design feature of TriBITS is that the contents of non-enabled
package's can't damnage the configure, build, and test of the enabled packages
based on errors in non-enabled packages.  This is critical to allow
experimental ``EX`` test packages and lower-maturity packages to exist in the
same soruce repositories safely.


TriBITS Package Core Variables
..............................

.. _TriBITS Package Local Varaibles:

The following locally scoped **TriBITS Package Local Varaibles** are defined
when the files for a given TriBITS Package (or any SE packaeg for that matter)
are being processed:

  ``PACKAGE_NAME``

    The name of the current TriBITS SE package.  This is set automatically by
    TriBITS before the packages's ``CMakeLists.txt`` file is processed.
    **WARNING:** The TriBITS packae name must be globally unique cross all
    possible repositories that might be coupled together at some point using
    TriBITS!

  ``PACKAGE_SOURCE_DIR``

    The absolute path to the base package's base source directory.  This is
    set automatically by TriBITS in the macro `TRIBITS_PACKAGE()`_.

  ``PACKAGE_BINARY_DIR``

    The absolute path to the base package's base binary/build directory.  This
    is set automatically by TriBITS in the macro `TRIBITS_PACKAGE()`_.

  ``PACKAGE_NAME_UC``

    This is set to the upper-case versionof ``${PACKAGE_NAME}``.  This is set
    automatically by TriBITS in the macro `TRIBITS_PACKAGE()`_.

.. _TriBITS Package Top-Level Local Variables:

Once all of the TriBITS SE package's ``Dependencies.cmake`` files have been
processed, the following **TriBITS Package Top-Level Local Variables** are
defined:

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

  ``${PACKAGE_NAME}_PARENT_REPOSITORY``

    The name of the package's parent repository.  This can be used by a
    package to access information about its parent repository.  For example,
    the variable ``${${PACKAGE_NAME}_PARENT_REPOSITORY}_SOURCE_DIR`` can be
    dereferenced.

.. _TriBITS Package Cache Variables:

In addition, the following user-settable **TriBITS Package Cache Variables**
are defined before an SE Package's ``CMakeLists.txt`` file is processed:

  ``${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}``

    Set to ``ON`` if the package is to be enabled.

  ``${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE_NAME}``

    Set to ``ON`` if support for the optional upstream dependent package
    ``${OPTIONAL_DEP_PACKAGE_NAME}`` is enabled in package
    ``${PACKAGE_NAME}``.  Here ``${OPTIONAL_DEP_PACKAGE_NAME}`` corresponds to
    each optional upstream SE package listed in the ``LIB_OPTIONAL_PACKAGES``
    and ``TEST_OPTIONAL_PACKAGES`` arguments to the
    `TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_ macro.

  ``${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEPENDENT_TPL_NAME}``

    Set to ``ON`` if support for the optional upstream dependent TPL
    ``${OPTIONAL_DEPENDENT_TPL_NAME}`` is enabled in package
    ``${PACKAGE_NAME}``.  Here ``${OPTIONAL_DEPENDENT_TPL_NAME}`` corresponds
    each to the optional upstream TPL listed in the ``LIB_OPTIONAL_TPLS`` and
    ``TEST_OPTIONAL_TPLS`` arguments to the
    `TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_ macro.

  ``${PACKAGE_NAME}_ENABLE_TESTS``

    Set to ``ON`` if the package's tests are to be enabled.  This will enable
    a package's tests and all of its subpackage's tests.

  ``${PACKAGE_NAME}_ENABLE_EXAMPLES``

    Set to ``ON`` if the package's examples are to be enabled.  This will
    enable a package's examples and all of its subpackage's examples.

The above varaibles can be set by the user or may be set automatically as part
of the `Package Dependencies and Enable/Disable Logic`_.

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

* Typically defines a set of libraries and/or header files and/or executables
  and/or tests with CMake build targets for building these and exports the
  list of include directories, libraries, and targets that are created (along
  with CMake dependencies).
* Is declared in its parent packages's
  `<packageDir>/cmake/Dependencies.cmake`_ file in a call to
  `TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_ using the argument
  ``SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS``.
* Defines dependencies on upstream TPLs and/or other SE packages by just
  naming the dependencies in the file ``cmake/Dependencies.cmake`` using the
  macro `TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_.
* Can **NOT** have its own subpackages defined (only top-level packages can
  have subpackages).
* Is enabled or disabled along with all other subpackages in the parent
  package automatically if it's parent package is enabled or disabled with
  ``${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME}`` set to ``ON`` or ``OFF``
  respectively.
* Has tests turned on automatically if
  ``${PARENT_PACKAGE_NAME}_ENABLE_TESTS==ON``.

The contents of a TriBITS Subpackage are almost idential to those of a TriBITS
Package.  The differences are described below.

For more details on the definition of a TriBITS Package (or subpackage), see:

* `TriBITS Subpackage Core Files`_
* `TriBITS Subpackage Core Variables`_

TriBITS Subpackage Core Files
..............................

The set of core files for a subpackage are identical to the `TriBITS Package
Core Files`_.  The core files that make up a TriBITS Subpackage (where
``<packageDir> = ${${PARENT_PACKAGE_NAME}_SOURCE_DIR}`` and ``<spkgDir>`` is
the subpackage directory listed in the
``SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS`` to
`TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_) are::

  <packageDir>/<spkgDir>/
    CMakeLists.txt  # Only processed if this subpackage is enabled
    cmake/
      Dependencies.cmake  # Always processed if the parent package
                            # is listed in the enclosing Repository

These TriBITS Subpackage files are documented in more detail below:

* `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_
* `<packageDir>/<spkgDir>/CMakeLists.txt`_
* `How is TriBITS Subpackage is different from a TriBITS Package?`_

.. _<packageDir>/<spkgDir>/cmake/Dependencies.cmake:

**<packageDir>/<spkgDir>/cmake/Dependencies.cmake**: The contents of this file
for subpackages is idential as for top-level packages.  It just contains a
call to the macro `TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_ to define this SE
package's upstream TPL and SE package dependencies.  A simple example is for
``SubpackageB`` declared in
`package_with_subpackages/cmake/Dependencies.cmake`_ shown shown in:

  `TribitsExampleProject`_/``packages/package_with_subpackages/B/cmake/Dependneices.cmake``

which is:

.. include:: ../examples/TribitsExampleProject/packages/package_with_subpackages/B/cmake/Dependencies.cmake
   :literal:

What this shows is that subpackages must list their dependenices on each other
(if such dependencies exist) using the full SE package name
``${PARENT_PACKAGE_NAME}${SUBPACKGE_NAME}`` or in this case
``PackageWithSubpackages`` + ``SubpackageA`` =
``PackageWithSubpackagesSubpackageA``.

Note that the parent SE package depends on its subpackages, not the other way
around.  For example, the ``PackageWithSubpackages`` parent SE package depends
its SE subpackages ``PackageWithSubpackagesSubpackageA``,
``PackageWithSubpackagesSubpackageC``, and
``PackageWithSubpackagesSubpackageC``.  As such all (direct) dependencies for
a subpackage must be listed in its own ``Dependencies.cmake`` file.  For
example, the ``PackageWithSubpackages`` subpackage ``SubpackageA`` depends on
the ``SimpleCxx`` package and is declared as such as shown in:

  `TribitsExampleProject`_/``packages/package_with_subpackages/A/cmake/Dependneices.cmake``

which is:

.. include:: ../examples/TribitsExampleProject/packages/package_with_subpackages/A/cmake/Dependencies.cmake
   :literal:

What this means is that any TPL or library depencenices listed in the parent
package's `<packageDir>/cmake/Dependencies.cmake`_ file are **NOT**
dependencies of its subpackages.  For example, if
`package_with_subpackages/cmake/Dependencies.cmake`_ where changed to be::

  TRIBITS_DEFINE_PACKAGE_DEPENDENCIES(
    LIB_REQUIRED_TPLS Boost
    SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
      SubpackageA   A   PT  REQUIRED
      ...
    )

then the ``Boost`` TPL would **NOT** be a dependency of the SE package
``PackageWithSubpackagesSubpackageA`` but instead would be listed as a
dependency of the parent SE package ``PackageWithSubpackages``.  (And in this
case, this TPL dependency is pretty worthless since the SE package
``PackageWithSubpackages`` does not even define any libraries or tests of its
own.)

.. _<packageDir>/<spkgDir>/CMakeLists.txt:

**<packageDir>/<spkgDir>/CMakeLists.txt**: [Required] The subpackage's
top-level ``CMakeLists.txt`` file that defines the libraries, include
directoies, and contains the tests for the subpackage.  The contents of a
subpackage's top-level ``CMakeLists.txt`` is almost indentical to a top-level
package's `<packageDir>/CMakeLists.txt`_ file.  The primary difference is that
the commmands `TRIBITS_PACKAGE()`_ and `TRIBITS_PACKAGE_POSTPROCESS()`_ and
replaced with `TRIBITS_SUBPACKAGE()`_ and `TRIBITS_SUBPACKAGE_POSTPROCESS()`_
as shown in the file:


  `TribitsExampleProject`_/``packages/package_with_subpackages/A/CMakeLists.txt``

which contains:

.. include:: ../examples/TribitsExampleProject/packages/package_with_subpackages/A/CMakeLists.txt
   :literal:

Unlike `TRIBITS_PACKAGE()`_, `TRIBITS_SUBPACKAGE()`_ does not take any extra
arguments.  Those extra settings are assumed to be defined by the top-level
parent package.  Like top-level packages, subpackages are free to define
user-setable options and configure-time tests but typically don't.  The idea
is that subpackages should be lighter weight than top-level packages.  OTher
than using `TRIBITS_SUBPACKAGE()`_ and `TRIBITS_SUBPACKAGE_POSTPROCESS()`_, a
subpackage can be layed out just like any other package and can call on any
other commands to add libraries, add executables, add test, etc.


TriBITS Subpackage Core Variables
.................................

The core varaibles assoicated with a subpackage are identicial to the `TriBITS
Package Core Variables`_.  The only difference is that a subpackage may need
to refer to its parent package where a top-level package does not have a
parent package.  The extra varaibles that are defined when processing a
subpackages files are:

  ``PARENT_PACKAGE_NAME``

    The name of the parent package.

  ``PARENT_PACKAGE_SOURCE_DIR``

    The absolute path to the parent package's source directory.  This this
    only defined for a subpackage.

  ``PARENT_PACKAGE_BINARY_DIR``

    The absolute path to the parent package's binary directory.  This this
    only defined for a subpackage.

How is TriBITS Subpackage is different from a TriBITS Package?
..............................................................

A common question this is natural to ask is how a TriBITS Subpackage is
different from a TriBITS Package?  They contain the same basic files (i.e. a
``cmake/Dependencies.cmake``, a top-level ``CMakeList.txt`` file, source
files, test files, etc.).  They both are included in the list of TriBITS SE
Packages and therefore can both be enabled/disabled by the user or in
automatica dependnecy logic.  The primary difference is that a subpackage is
meant to involve less overhead in defining and is to be used to partition the
parent package's software into chunks according to software engineering
packaging principles.  Also, the dependency logic treats a parent package's
subpackages as part of itself so when the parent package is explicitly enabled
or disabled, it is identical to explicitly enabling or disabling all of its
subpackages.  Other differences and issues between packages as subpackages are
discussed throughout this guide.

.. ToDo: List out the basic SE packaging design principles from "Agile
   Software Design".

.. ToDo: Finish this section!

TriBITS TPL
+++++++++++

A TriBITS TPL:

* Defines a set of pre-built libraries and/or header files and/or executables
  and/or some other resoruces used by one or more TriBITS Packages and
  publishes the list of include directories and/or libraries and/or
  executables provided by the TPL to the TriBITS project.
* Is declared in a `<repoDir>/TPLsList.cmake` file.  
* Is listed as an explicit optional or required dependency in one or more
  TriBITS SE package's `<packageDir>/cmake/Dependencies.cmake`_ files.

Using a TriBITS TPL is to be preferred over using a raw CMake
``FIND_PACKAGE(<someCMakePackage>)`` because the TriBITS system guarantees
that only a single unique version of TPL will be used by multiple packages and
by declaring a TPL using TriBITS, automatical enable/disable logic will be
applied as described in `Package Dependencies and Enable/Disable Logic`_.

For each TPL referenced in a ``TPLsList.cmake`` file using the macro
`TRIBITS_DEFINE_REPOSITORY_TPLS()`_, there should exist a file, typically
called ``FindTPL${TPL_NAME}.cmake``, that once processed, produces the
variables ``${TPL_NAME}_LIBRARIES`` and ``${TPL_NAME}_INCLUDE_DIRS``.  Most
``FindTPL${TPL_NAME}.cmake`` files just use the the function
`TRIBITS_TPL_DECLARE_LIBRARIES()`_ the define the TriBITS TPL.  A simple
example of such a file is the standard ``FindTPLPETSC.cmake`` module which is:

.. include:: ../../tpls/FindTPLPETSC.cmake
   :literal:

Some concrete ``FindTPL${TPL_NAME}.cmake`` files actually do use
``FIND_PACKAGE()`` and a standard CMake package find modulue to fill in the
guts of finding at TPL.

Note that the TriBITS system does not require the usage of of the function
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

One of the most important things to know about TriBITS is what files it
processes, in what order, and in what context.  This is critical to being able
to understand what impact (if any) setting a variable or otherwise changing
the CMake runtime state will have on configuring a CMake project which uses
TriBITS.  While the different files that make up a `TriBITS Project`_,
`TriBITS Repository`_, `TriBITS Package`_, and `TriBITS TPL`_ were defined in
the section `TriBITS Project Structure`_, that material did not fully describe
the context and in what order these files are processed by the TriBITS
framework.

The TriBITS system processes the project's files in one of two use cases.  The
first use case is in the basic configuration of the project with a standard
``cmake`` command invocation in order to set up the build files in the binary
directory (see `Full TriBITS Project Configuration`_).  The second use case is
in reading the project's dependency-related files in order to build a package
dependency datastructure (e.g. the `<Project>PackageDependencies.xml`_ file,
see `Reduced Package Dependency Processing`_)).  The second use case of
reading the project's dependency files is largely a subset of the first.

Another factor that is important to understand is the scoping in which the
various files are processed (with ``INCLUDE()`` or ``ADD_SUBDIRECTORY()``).
This scoping has a large impact on the configuration of the project and what
effect the processing of files and setting varaibles have on the project as a
whole.  Some of the strange scoping rules for CMake are discussed in `CMake
Language Overview and Gotchas`_ and should be understood before trying to
debug issues with processesing.  Many of the basic files are processed
(included) in the base project `<projectDir>/CMakeLists.txt`_ scope and
therefore any local variables set in these files are accessible to the entire
CMake project (after the file is processed, of course).  Other files get
processed inside of functions which have their own local scope and therefore
only impact the rest of the project in more purposeful ways.

Full TriBITS Project Configuration
++++++++++++++++++++++++++++++++++

The first use case to describe is the full processing of all of the TriBITS
project's files starting with the base `<projectDir>/CMakeLists.txt`_ file.
This begins with the invocation of the command::

  $ cmake [options] <projectDir>

Below, is a short pseudo-code algorithm for the TriBITS framework processing
and callbacks that begin in the `<projectDir>/CMakeLists.txt`_ and proceed
through the call to `TRIBITS_PROJECT()`_.

.. _Full Processing of TriBITS Project Files:

**Full Processing of TriBITS Project Files:**

| 1.  Read `<projectDir>/ProjectName.cmake`_ (sets ``PROJECT_NAME``)
| 2.  Call ``PROJECT(${PROJECT_NAME} NONE)`` (sets ``${PROJECT_NAME}_SOURCE_DIR``
|     and ``${PROJECT_NAME}_BINARY_DIR``)
| 3.  Execute `TRIBITS_PROJECT()`_:
|   1)  Set ``PROJECT_SOURCE_DIR`` and ``PROJECT_BINARY_DIR``
|   2)  For each ``<optFilei>`` in ``${${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE}``:
|       * ``INCLUDE(<optFilei>)``
|   3)  Set variables ``CMAKE_HOST_SYSTEM_NAME`` and ``${PROJECT_NAME}_HOSTNAME``
|       (both of these can be overridden in the cache by the user)
|   4)  Find Python (sets ``PYTHON_EXECUTABLE``)
|   5)  ``INCLUDE(`` `<projectDir>/Version.cmake`_ ``)``
|   6)  Define primary TriBITS options and read in the list of extra repositories
|       (calls ``TRIBITS_DEFINE_GLOBAL_OPTIONS_AND_DEFINE_EXTRA_REPOS()``)
|       * ``INCLUDE(`` `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ ``)``      
|   7)  For each ``<repoDir>`` in all defined TriBITS repositories:
|       * ``INCLUDE(`` `<repoDir>/cmake/CallbackSetupExtraOptions.cmake`_ ``)``
|       * Call ``TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS()``
|   9)  Call ``TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()``:
|     a)  For each ``<repoDir>`` in all defined TriBITS repositories:
|         * ``INCLUDE(`` `<repoDir>/PackagesList.cmake`_ ``)``
|         * ``INCLUDE(`` `<repoDir>/TPLsList.cmake`_ ``)``
|     b)  For each ``<repoDir>`` in all defined TriBITS repositories:
|         * ``INCLUDE(`` `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ ``)``
|     c)  ``INCLUDE(`` `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_ ``)``
|     d)  For each ``<packageDir>`` in all defined top-level packages:
|         * ``INCLUDE(`` `<packageDir>/cmake/Dependencies.cmake`_ ``)``
|           - Sets all package-specific options (see `TriBITS Package Cache Variables`_)
|         * For each ``<spkgDir>`` in all subpackages for this package:
|           * ``INCLUDE(`` `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_ ``)``
|             - Sets all subpackage-specific options
|   10) Adjust SE packae and TPLs enable/disable
|       (see `Package Dependencies and Enable/Disable Logic`_)
|   11) `Probe and set up the environment`_ (finds MPI, compilers, etc.)
|       (see `TriBITS Environment Probing and Setup`_)
|   12) For each enabled TPL, ``INCLUDE(FindTPL<tplName>.cmake`` (see `TriBITS TPL`_)
|   13) For each ``<repoDir>`` in all defined TriBITS repositories:
|       * Read `<repoDir>/Copyright.txt`_
|       * ``INCLUDE(`` `<repoDir>/Version.cmake`_ ``)``
|       (see `Project and Repositiory Versioning and Release Mode`_)
|   14) For each ``<packageDir>`` in all enabled top-level packages
|       * ``ADD_SUBDIRECTORY(`` `<packageDir>/CMakeLists.txt`_ ``)``
|       * For each ``<spkgDir>`` in all enabled subpackages for this package:
|         * ``ADD_SUBDIRECTORY(`` `<packageDir>/<spkgDir>/CMakeLists.txt`_ ``)``
|   15) For each ``<repoDir>`` in all defined TriBITS repositories:
|       * ``INCLUDE(`` `<repoDir>/cmake/CallbackDefineRepositoryPackaging.cmake`_ ``)``
|       * Call ``TRIBITS_REPOSITORY_DEFINE_PACKAGING()``
|   17) ``INCLUDE(`` `<projectDir>/cmake/CallbackDefineProjectPackaging.cmake`_ ``)``
|       * Call ``TRIBITS_PROJECT_DEFINE_PACKAGING()``

The TriBITS Framework obviously does a lot more that what is described above
but the basic trace of major operations and ordering and the processing of
project, repository, package, and subpackage files should be clear.  All of
this information should also be clear when enabling `File Processing
Tracing`_. 

Reduced Package Dependency Processing 
++++++++++++++++++++++++++++++++++++++

In addition to the full processing that occurs as part of the `Full TriBITS
Project Configuration`_, there are also TriBITS tools that only process as
subset of project's file.  This reduced processing is performed in order to
build up the project's package dependenices data-structure (see `TriBITS
Environment Probing and Setup`_) and to write the file
`<Project>PackageDependencies.xml`_.  For example, the tool ``checkn-test.py``
and the script ``TribitsCTestDriverCore.cmake`` both drive this type of
processing.  In particular, the CMake -P script
``TribitsDumpDepsXmlScript.cmake`` reads all of the project's
dependency-related files and dumps out the `<Project>PackageDependencies.xml`_
file a defined set of native and extra repositories defined for the project.
This reduced processing is given below.

.. _Reduced Dependency Processing of TriBITS Project Files:

**Reduced Dependency Processing of TriBITS Project:**

| 1.  Read `<projectDir>/ProjectName.cmake`_ (sets ``PROJECT_NAME``)
| 2. ``INCLUDE(`` `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ ``)``      
| 3.  Call ``TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()``:
|   a)  For each ``<repoDir>`` in all defined TriBITS repositories:
|       * ``INCLUDE(`` `<repoDir>/PackagesList.cmake`_ ``)``
|       * ``INCLUDE(`` `<repoDir>/TPLsList.cmake`_ ``)``
|   b)  For each ``<repoDir>`` in all defined TriBITS repositories:
|       * ``INCLUDE(`` `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ ``)``
|   c)  ``INCLUDE(`` `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_ ``)``
|   d)  For each ``<packageDir>`` in all defined top-level packages:
|       * ``INCLUDE(`` `<packageDir>/cmake/Dependencies.cmake`_ ``)``
|         - Sets all package-specific options (see `TriBITS Package Cache Variables`_)
|       * For each ``<spkgDir>`` in all subpackages for this package:
|         * ``INCLUDE(`` `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_ ``)``
|           - Sets all subpackage-specific options

When comparing the above reduced dependency processing to the `Full Processing
of TriBITS Project Files`_ it is important to note that that several files are
**not** processed in these cases.  The files that are not processed include
`<projectDir>/Version.cmake`_, `<repoDir>/Version.cmake`_ and
`<repoDir>/cmake/CallbackSetupExtraOptions.cmake`_.  Therefore, you can't put
anything in these files that would impact the definition of TriBITS
repositories, packages, TPLs, etc.  Anything that would affect the
dependencies data-structure that gets written out as
`<Project>PackageDependencies.xml`_ must be contained in the files that are
processed shown above.

Debugging issues with `Reduced Dependency Processing of TriBITS Project Files`
is more difficult because one can not easily turn on `File Processing Tracing`
like they can when doing the full CMake configure.  However, options may be
added to the various tools to show this file processing and help debug
problems.

File Processing Tracing
+++++++++++++++++++++++

In order to aid in debugging problems with configuration, TriBITS defines the
CMake cache option ``${PROJECT_NAME}_TRACE_FILE_PROCESSING``.  When enabled,
TriBITS will print out when any of the project-related, repository-related, or
package-related file is being processed by TriBITS.  When
``${PROJECT_NAME}_TRACE_FILE_PROCESSING=ON``, lines starting with ``"-- File
Trace:"`` are printed in the ``cmake`` STDOUT for files that TriBITS
automatically processes where there may be any confusion about what files are
processed and when. 

For example, for `TribitsExampleProject`_, the configure file trace for the
configure command::

  $ cmake \
    -DTribitsExProj_TRIBITS_DIR=<tribitsDir> \
    -DTribitsExProj_ENABLE_MPI=ON \
    -DTribitsExProj_ENABLE_ALL_PACKAGES=ON \
    -DTribitsExProj_ENABLE_TESTS=ON \
    -DTribitsExProj_TRACE_FILE_PROCESSING=ON \
    -DTribitsExProj_ENABLE_CPACK_PACKAGING=ON \
    -DTribitsExProj_DUMP_CPACK_SOURCE_IGNORE_FILES=ON \
    <tribitsDir>/doc/examples/TribitsExampleProject \
    | grep "^-- File Trace:"

looks something like::

  -- File Trace: PROJECT    INCLUDE    [...]/Version.cmake
  -- File Trace: REPOSITORY INCLUDE    [...]/cmake/CallbackSetupExtraOptions.cmake
  -- File Trace: REPOSITORY INCLUDE    [...]/PackagesList.cmake
  -- File Trace: REPOSITORY INCLUDE    [...]/TPLsList.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/simple_cxx/cmake/Dependencies.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/mixed_language/cmake/Dependencies.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/package_with_subpackages/cmake/Dependencies.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/package_with_subpackages/A/cmake/Dependencies.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/package_with_subpackages/B/cmake/Dependencies.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/package_with_subpackages/C/cmake/Dependencies.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/wrap_external/cmake/Dependencies.cmake
  -- File Trace: PROJECT    CONFIGURE  [...]/cmake/ctest/CTestCustom.cmake.in
  -- File Trace: REPOSITORY READ       [...]/Copyright.txt
  -- File Trace: REPOSITORY INCLUDE    [...]/Version.cmake
      "${TPL_MPI_FILE_TRACE}
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/simple_cxx/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/simple_cxx/test/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/mixed_language/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/mixed_language/test/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/package_with_subpackages/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/package_with_subpackages/A/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/package_with_subpackages/A/tests/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/package_with_subpackages/B/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/package_with_subpackages/B/tests/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/package_with_subpackages/C/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/package_with_subpackages/C/tests/CMakeLists.txt
  -- File Trace: REPOSITORY INCLUDE    [...]/cmake/CallbackDefineRepositoryPackaging.cmake
  -- File Trace: PROJECT    INCLUDE    [...]/cmake/CallbackDefineProjectPackaging.cmake

However, every file that TriBITS processes is not printed in this file trace
if it should be obvious that the file is being processed.  For example, the
package's configured header file created using `TRIBITS_CONFIGURE_FILE()`_
does not result in a file trace print statment because this is an
unconditional command that is explicitly called in one of the the package's
``CMakeLists.txt`` files so it should be clear that this file is being
processed.


Coexisting Projects, Repositories, and Packages
-----------------------------------------------

Certain simplifications are allowed when defining TriBITS projects,
repositories and packages.  The known allowed simplifications are described
below.

**TriBITS Repository == TriBITS Project**: It is allowed for a TriBITS Project
and a TriBITS Repository to be the same source directory and in fact this is
the default for every TriBITS project (unless the the
`<projectDir>/cmake/NativeRepositoriesList.cmake`_ is defined).  In this case,
the repository name and the project name are the same as well.  This is quite
common and is in fact the default that every TriBITS Project is also a TriBITS
repository (and therefore must contain `<repoDir>/PackagesList.cmake`_ and
`<repoDir>/TPLsList.cmake`_ files).  This is the case, for example, with the
the Trilinos and the `TribitsExampleProject`_ projects and repositories.  In
this case, the Project's and the Repository's ``Version.cmake`` and
``Copyright.txt`` files are one and the same, as they should be (see `Project
and Repositiory Versioning and Release Mode`_).

**TriBITS Package == TriBITS Repository**: It is also allowed for a TriBITS
Repository to have only one package and to have that package be the base
repository directory.  The TriBITS Repository and the single TriBITS Package
would typically have the same name in this case (but that is actually not
required but it is confusing if they are not the same).  For example, in the
TriBITS test project ``MockTrilinos``, the repostory and package
``extraRepoOnePackage`` are one in the same.  In this case, the file
``extraRepoOnePackage/PackagesList.cmake`` looks like:

.. include:: ../../package_arch/UnitTests/MockTrilinos/extraRepoOnePackage/PackagesList.cmake
   :literal:

This is used in the real TriBITS repository `DataTransferKit`_.

.. _DataTransferKit: https://github.com/CNERG/DataTransferKit

However, to maximize flexibility, it is recommended that a TriBITS package and
TriBITS repository not share the same directory.

.. _TriBITS Package, TriBITS Repository, TriBITS Package sharing the same source directory:

**TriBITS Package, TriBITS Repository, TriBITS Package sharing the same source
directory**: In the extreme, it is posisble to collapase a single TriBITS
package, repository, and project into the same base source directory.
However, in this case, the TriBITS Project name and the TriBITS Package name
cannot be the same and some modifications and tricks are needed to allow this
to work.  One example of this is the TriBITS project and `The TriBITS Test
Package`_ themselves, which are both rooted in the base ``tribits`` source
directory.  There are a few restructions and modifications needed to get this
to work:

* The Project and Package names cannot be the same: In the case of the TriBITS
  project, the project name is ``TriBITSProj`` (as defined in
  ``tribits/ProjectName.cmake``) and the package name is ``TriBITS`` (as
  defined in ``tribits/PackagesList.cmake``).
* The base ``CMakeLists.txt`` file must be modified to allow it to be
  processed both as the base project ``CMakeLists.txt`` file and as the
  package's base ``CMakeLists.txt`` file: In the case of
  ``tribits/CMakeLists.txt``, a big if statement is used.
* An extra subdirectory must be created for TriBITS package's binary
  directory: Because of directory-level targets like ``${PROJECT_NAME}_libs``
  and ``${PACKAGE_NAME}_libs``, a subdirectory for package's the binary
  directory must be created.  This is simply done by overriding the binary
  directory name ``${PACKAGE_NAME}_SPECIFIED_BINARY_DIR``.  In the case of
  TriBITS, this is set to ``tribits`` in the ``tribits/PackagesList.cmake``
  file.

Other than those modifications, a TriBITS project, repository, and package can
all be rooted in the same source directory.  However, as one can see above, to
do so is a little messy and is not recommended.  It was only done this way
with the base TriBITS directory in order to maintain backward compatibility
for the usage of TriBITS in existing TriBITS projects.

However, one possible use case for collapsing a project, repository, and
package into a single base source directory would be to support the
stand-alone build of a TriBITS package as its own entity that uses an
installation of the TriBITS.  If a given TriBITS package has no required
upstream TriBITS package dependencies and minimal TPL dependencies (or only
uses `Standard TriBITS TPLs`_ already defined in the ``tribits/tpls/``
directory), then creating a stand-alone project build of a loan TriBITS
package requires fairly little extra overhead or duplication.  However, as
mentioned above, one cannot use the same name for the package and the project.

.. NOTE: We could modify the TriBITS system to allow having the project and
.. package names be the same if we disable one of the sets of XXX_libs and
.. other targets created by TriBITS.  We could also provide a standard wrapper
.. CMakeLists.txt file to make it easy for packages to support stand-alone
.. builds.  However, that is an effort for later.


Standard TriBITS TPLs
---------------------

TriBITS contains find modules for a few standard TPLs that are either integral
to the TriBITS system or are likely to be used across many independent TriBITS
repositories.  The goal of maintaining a few of these in the later case under
TriBITS is to enforce conformity in case these independent repositories are
combined into a single metra-project.

The standard TriBITS TPLs are contained under the directory::

  tribits/tpls/

The current list of standard TriBITS TPLs is:

.. include:: TribitsStandardTPLsList.txt
   :literal:

The TPLs ``MPI`` and ``CUDA`` are standard because they are special in that
they define compilers and other special tools that are used in
`TRIBITS_ADD_LIBRARY()`_, `TRIBITS_ADD_EXECUTABLE()`_, `TRIBITS_ADD_TEST()`_
and other commands.

These standard TPLs are used in a `<repoDir>/TPLsList.cmake`_ file as::

  TRIBITS_DEFINE_REPOSITORY_TPLS(
    MPI   "${${PROJECT_NAME}_TRIBITS_DIR}/tpls/"  PT
    CUDA  "${${PROJECT_NAME}_TRIBITS_DIR}/tpls/"  ST
    ...
    )

Other than the special TPLs ``MPI`` and ``CUDA``, other TPLs that are
candidates to put into TriBITS are those that are likley to be used by
different stand-alone TriBITS repositories that need to be combined into a
single TriBITS meta-project.  By using a standard TPL definition, it is
guaranteed that the TPL used will be consistent with all of the repositories.

Note that just because packages in two repositories reference the same TPL
does not necessarily mean that it needs to be a standard TriBITS TPL.  For
example, if the TPL ``BLAS`` is defined in an upstream repository
(e.g. Trilinos), then a package in a downstream repository can list a
dependency on the TPL ``BLAS`` without having to define its own ``BLAS`` TPL
in its repository's `<repoDir>/TPLsList.cmake`_ file.  For more details on
TPLs, see `TriBITS TPL`_.


.. Where to set variables?
.. -----------------------

.. A TriBITS project, repository, package have a number of files in which
   varaibles can be set and commands can be called in order to affect how the
   project is defined and configured.  With so many files that can be included by
   TriBITS, it can be difficult to know where the right place is to set a given
   set of variables.  The primary considerations for where to set a variable
   depend on:
 
.. ToDo: Describe considerations on where to set variables ...


Example TriBITS Projects
=========================

In this section, a few different example TriBITS projects and packages are
previewed.  All of these examples exist in the TriBITS source directory
``tribits`` itself so they are available to all users of TriBITS.  These
examples also provide a means to test the TriBITS system itself (see `The
TriBITS Test Package`_).

The first example covered is the bare bones `TribitsHelloWorld`_ example
project.  The second example covered in detail is `TribitsExampleProject`_.
This example covers all the basics for setting up a simple multi-package
TriBITS project.  The third example outlined is `MockTrilinos` which mostly
exists to test the TriBITS system itself but use contains some nice examples
of a few different TriBITS features and behaviors.  The last example mentioned
is `The TriBITS Test Package`_ itself which allows the TriBITS system to be
tested and installed from any TriBITS project that lists it, including the
``TriBITSProj`` project itself (see `Coexisting Projects, Repositories, and
Packages`_).

The directory ``tribits/doc/examples/`` contains some other example TriBITS
projects and repositories as well that are refered to in this and other
documents.


TribitsHelloWorld
-----------------

This is the simplest possible TriBITS project that you can imagine and is
contained under the directory::

  tribits/doc/examples/TribitsHelloWorld/

It contains only a single TriBITS package and no frills at all (does not
support MPI or Fortran).  However, it does show how minimal a `TriBITS
Project`_ (which is also a `TriBITS Repository`_) and a `TriBITS Package`_ can
be and still show the value of TriBITS over raw CMake.  The simple
`HelloWorld` package is used to compare with the raw CMakeList.txt file in the
``RawHeloWorld`` example project in the `TriBITS Overview`_ document.

The directory structure for this examples shows what is necessary for a
minimal TriBITS project:

.. include:: TribitsHelloWorldDirAndFiles.txt
   :literal:

This has all of the required `TriBITS Project Core Files`_, `TriBITS
Repository Core Files`_, and `TriBITS Package Core Files`_.  It just build a
simle library, a simple exectuable, a test exectuable, and the tests them as
shown by the file ``TribitsHelloWorld/hello_world/CMakeLists.txt`` which is:

.. include:: ../examples/TribitsHelloWorld/hello_world/CMakeLists.txt
   :literal:

The build and test of this simple project is tested in the `TriBITS Package`_
testing file::

  tribits/doc/examples/UnitTests/CMakeLists.txt

Note that this little example is a fully functional `TriBITS Repository`_ and
can be embedded in to a larger TriBITS metra-project and be seamlessly built
along with any other such TriBITS-based software.

.. ToDo: Put in reference to the example meta-project.


TribitsExampleProject
----------------------

TribitsExampleProject in an example `TriBITS Project`_ and `TriBITS
Repository`_ contained in the TriBITS source tree under::

  tribits/doc/examples/TribitsExampleProject/

When this used as the base TriBITS project, this is the directory coresponds
to ``<projectDir>`` and ``<repoDir>`` referenced in `TriBITS Project Core
Files`_ and `TriBITS Repository Core Files`_, respectively.

Several files from this project were used as examples in the section `TriBITS
Project Structure`_.  Here, a fuller description is given of this project and
how TriBITS works using it.  From this simple example project, one can quickly
see how the basic structural elements a TriBITS project, repository, and
package (and subpackage) are pulled together.

The name of this project ``PROJECT_NAME`` given in its
``TribitsExampleProject/ProjectName.cmake`` file:

.. include:: ../examples/TribitsExampleProject/ProjectName.cmake
   :literal:

The variable ``PROJECT_NAME=TribitsExProj`` is used to prefix with
``"${PROJECT_NAME}_"`` all of the projects global TriBITS variables like
``TribitsExProj_ENABLE_TESTS``, ``TribitsExProj_ENABLE_ALL_PACKAGES``, etc.
Note, as shown in this example, the project name and the base project
directory name do **not** need to match.

.. _TribitsExampleProject Files and Directories:

The directory structure and key files for this example project is shown in
this partial list of **TribitsExampleProject Files and Directories**::

  TribitsExampleProject/
    CMakeLists.txt
    Copyright.txt
    PackagesList.cmake
    ProjectName.cmake
    TPLsList.cmake
    Version.cmake
    ...
    cmake/
      CallbackDefineProjectPackaging.cmake
      CallbackDefineRepositoryPackaging.cmake
      CallbackSetupExtraOptions.cmake
    packages/
      simple_cxx/
        CMakeLists.txt
        cmake/
          CheckFor__int64.cmake
          Dependencies.cmake
          SimpleCxx_config.h.in
        src/
          CMakeLists.txt
          SimpleCxx_HelloWorld.cpp
          SimpleCxx_HelloWorld.hpp
        test/
          CMakeLists.txt
          SimpleCxx_HelloWorld_Tests.cpp
      mixed_language/ ...
      package_with_subpackages/
        CMakeLists.txt
        cmake/
          Dependencies.cmake
        A/
          CMakeLists.txt
          cmake/
            Dependencies.cmake
          ...
        B/ ...
        C/ ...
      wrap_external/ ...


Above, the subdirectoires under ``packages/`` are sorted according to the
order listed in the ``TribitsExampleProject/PackagesList.cmake`` file:

.. include:: ../examples/TribitsExampleProject/PackagesList.cmake
   :literal:

From this file, we get the list of top-level packages ``SimpleCxx``,
``MixedLanguage``, ``PackageWithSubpackages``, and ``WrapExternal`` (and their
base package directories and testing group, see
`<repoDir>/PackagesList.cmake`_).

The full listing of package files in `TribitsExampleProject Files and
Directories`_ is only shown for the ``SimpleCxx`` package directory
``packages/simple_cxx/``.  This gives ``<packageDir> =
<repoDir>/packages/simple_cxx`` for the package ``PACKAGE_NAME = SimpleCxx``
referenced in `TriBITS Package Core Files`_.  As explained there, the files
`<packageDir>/cmake/Dependencies.cmake`_ and `<packageDir>/CMakeLists.txt`_
must exist for every package directory listed in
`<repoDir>/PackagesList.cmake`_ and we see these files under in the directory
``packages/simple_cxx/``.  The package ``SimpleCxx`` does not have any
upstream SE package dependencies.

Now consider the example top-level package ``PackageWithSubpackages`` which,
as the name suggests, is broken down into subpackages.  The
``PackageWithSubpackages`` dependencies file::

  TribitsExampleProject/packages/package_with_subpackages/cmake/Dependencies.cmake

with contents:

.. include:: ../examples/TribitsExampleProject/packages/package_with_subpackages/cmake/Dependencies.cmake
   :literal:

references the three subpackage with sub-directories ``<spkgDir>`` = ``A``,
``B``, and ``C`` under the parent package directory
``packages/package_with_packages/`` which are shown in `TribitsExampleProject
Files and Directories`_.  This gives another set of three SE packages
``PackageWithSubpackagesSubpackageA``, ``PackageWithSubpackagesSubpackaeB``,
and ``PackageWithSubpackagesSubpackageC``.  Combining ``<packageDir> =
packages/package_with_packages`` and ``<spkgDir>`` for each subpackage gives
the subpackage directories::

  TribitsExampleProject/packages/package_with_subpackages/A/
  TribitsExampleProject/packages/package_with_subpackages/B/
  TribitsExampleProject/packages/package_with_subpackages/C/

Together with the top-level parent SE package ``PackageWithSubpackages``
itself, this top-level package provides four SE packages giving the final list
of SE packages provided by this TriBITS repo as::

  SimpleCxx MixedLanguage PackageWithSubpackagesSubpackageA \
    PackageWithSubpackagesSubpackaeB PackageWithSubpackagesSubpackaeC \
    PackageWithSubpackages WrapExternal 7

The above list of SE packages is shown formatted this way since this is the
format that the SE packages are printed by TriBITS in the ``cmake`` STDOUT on
the line starting with ``"Final set of non-enabled SE packages:"`` when no
packages are enabled (see `Selecting the list of packages to enable`_).
TriBITS defines enable/disable cache variables for each of the defined SE
packages like ``TribitsExProj_ENBLE_SimpleCxx``,
``TribitsExProj_ENBLE_PackageWithSubpackagesSubpackageA``, and defines all the
variables listed in `TriBITS Package Cache Variables`_ that are settable by
the users or by the dependnecy logic described in section `Package
Dependencies and Enable/Disable Logic`_.

Hopefully this simple project shows how what is listed in files:

* `<repoDir>/PackagesList.cmake`_,
* `<packageDir>/cmake/Dependencies.cmake`_, and
* `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_

is used to specify the packages and SE packages in a TriBITS project and
repository.  More details about the contents of the ``Dependencies.cmake``
files is described in the section `Package Dependencies and Enable/Disable
Logic`_.

When starting a new TriBITS project, repository, or package, one should
consider basing these on the examples in this project.  In fact, the skeletons
for any of the

* `TriBITS Project Core Files`_,
* `TriBITS Repository Core Files`_,
* `TriBITS Package Core Files`_, or
* `TriBITS Subpackage Core Files`_

should be copied from this example project as they represent best practice
when using TriBITS for the typical use cases.


MockTrilinos
-------------

The TriBITS project ``MockTrilinos`` is contained under the directory::

  tribits/package_arch/UnitTests/MockTrilinos/

This TriBITS project is not a full TriBITS project (i.e. it does not build
anything).  Instead, it is used to test the TriBITS system using tests defined
in the `The TriBITS Test Package`_.  The ``MockTrilinos`` project is actually
given the name ``PROJECT_NAME = Trilinos`` and contains a subset of packages
with slightly modified dependencies from a snapshot of the real Trilinos
project from May 2009.  The list of packages in::

  tribits/package_arch/UnitTests/MockTrilinos/PackagesList.cmake

is:

.. include:: ../../package_arch/UnitTests/MockTrilinos/PackagesList.cmake
   :literal:

All of the package directories listed above have ``cmake/Dependenices.cmake``
files but generally do not have ``CMakeLists.txt`` files since most of the
testing of ``MockTrilinos`` just involves dependency handling.

``MockTrilinos`` also contains a number of extra TriBITS repositories used in
various tests.  These extra repositories offer examples of different types of
TriBITS repositories like:

* ``extraRepoOnePackage``: Contains just the single package
  ``extraRepoOnePackage`` which is defined in the base repository directory.

* ``extraRepoOnePackageThreeSubpackages``: Contains just the single package
  ``extraRepoOnePackageThreeSubpackages`` which is defined in the base
  repository directory but is broken up into subpackages.

* ``extraRepoTwoPackages``: Contains just two packages but provides an exmaple
  of defining multiple repositories with possible missing required and
  optional upstream packages (see `Multi-Repository Support`_).

* ``extraTrilinosRepo``: Just a typical extra repo with add-on packages and
  new TPLs defined that depends on a few ``MockTrilinos`` packages.

New test extra repostories are added when new types of tests are needed that
would require new package and TPL dependency structures since existing
dependency tests based on ``MockTrilinos`` are expensive to change by their
very nature.

The reason that the ``MockTrilinos`` test project is mentioned in this
developers guide is because it contains a greater variety of packages,
subpackages, and TPLs with a greater variety of different types of
dependencies.  This variety is is needed to fully test the TriBITS system but
this project and the tests also serve as examples and extra documentation for
the behavior of the TriBITS system.  Several of the examples referenced in
this document come from ``MockTrilinos``.

Most of the dependency tests involving ``MockTrilinos`` are specified in::

  tribits/package_arch/UnitTests/DependencyUnitTests/CMakeLists.txt

A great deal about the current behavior of TriBITS `Package Dependencies and
Enable/Disable Logic`_ can be learned from inspecting the tests defined in
this ``CMakeLists.txt`` file.  There are also some faster-running unit tests
involving ``MockTrilinos`` defined in the file::

  tribits/package_arch/UnitTests/TribitsAdjustPackageEnables_UnitTests.cmake


ReducedMockTrilinos
-------------------

The TriBITS project ``ReducedMockTrilinos`` is contained under the directory::

  tribits/package_arch/UnitTests/ReducedMockTrilinos/

It is a scaled-down version of the `MockTrilinos`_ test project with just a
handful of packages and some modified dependencies.  Its primary purpose is to
be used for examples in the section `Package Dependencies and Enable/Disable
Logic`_ and to test a few features of the TriBITS system not tested in other
tests.

The list of packages in::

  tribits/package_arch/UnitTests/ReducedMockTrilinos/PackagesList.cmake

is:

.. include:: ../../package_arch/UnitTests/ReducedMockTrilinos/PackagesList.cmake
   :literal:

All of the listed packages are standard TriBITS packages except for the mock
``Thyra`` package which is broken down into subpackages.  More details of this
example project are described in `Package Dependencies and Enable/Disable
Logic`_.


The TriBITS Test Package
------------------------

The last TriBITS example mentioned here is the TriBITS test package named
(appropriately) ``TriBITS`` itself.  The directory for the ``TriBITS`` package
is the base TriBITS source directory ``tribits``.  This allows any TriBITS
project to add testing for the TriBITS system by just listing this package and
its directory in its repository's `<repoDir>/PackagesList.cmake`_ file.  For
example, the Trilinos repository which currently snaphsots the TriBITS source
tree lists the ``TriBITS`` package with::

  TRIBITS_DEFINE_REPOSITORY_PACKAGES(
    TriBITS   cmake/tribits  PT   # Only tests, no libraries/capabilities!
    ...
    )

No downstream packages list a dependency on ``TriBITS`` in their
`<packageDir>/cmake/Dependencies.cmake`_ files.  Listing the ``TriBITS``
package in only done in the ``PackagesList.cmake`` file for testing TriBITS.

Other TriBITS projects/repositories that don't snapshot TriBITS but also want
to test TriBITS (perhaps just to mine the running tests for examples) can do
so by including the ``TriBITS`` test package in their ``PackagesList.cmake``
file using::

  TRIBITS_DEFINE_REPOSITORY_PACKAGES(
    TriBITS   ${${PROJECT_NAME}_TRIBITS_DIR}   PT
    ...
    )

Once the ``TriBITS`` test package is added to the list of project/repository
packages, it can be enabled just like any other package by adding the
following to the ``cmake`` command-line options::

  -D <Project>_ENABLE_TriBITS=ON \
  -D <Project>_ENABLE_TESTS=ON

One can then inspect the added tests prefixed by ``"TriBITS_"`` to see what
tests are defined and how they are run.  There is a wealth of information
about the TriBITS system embedded in these tests and where documentation and
these tests disagreed, believe the tests!


Package Dependencies and Enable/Disable Logic
=============================================

Arguably, the more important feature/aspect of the TriBITS system is the
partitioning of a large software project into pakages and managing the
dependencies between these packages to support building, testing, and depoying
different pieces as needed.  This is especially useful in incremental CI
testing of large projects.  However, this is also a critical component in
creating and maintaining Self Sustaining Software (see the `TriBITS Lifecycle
Model`_).  The fundamental mechanism for breaking up a large software into
manageable pieces is to partition the software into different `TriBITS
Packages`_ and then define the dependencies between these packages (which are
defined inside of the `<packageDir>/cmake/Dependencies.cmake`_ files for each
package).

Example ReducedMockTrilinos Project Dependency Structure
--------------------------------------------------------

To demonstrate the TriBITS package and TPL dependency handling system, the
small simple `ReducedMockTrilinos`_ project is used.  The list of packages for
this project is defined in ``ReducedMockTrilinos/PackagesList.cmake`` (see
`<repoDir>/PackagesList.cmake`_) which contents:

.. include:: ../../package_arch/UnitTests/ReducedMockTrilinos/PackagesList.cmake
   :literal:

This gives the full list of top-level packages::

  Teuchos RTOp Epetra Triutils EpetraExt Thyra

All of the listed packages are standard TriBITS packages except for the mock
``Thyra`` package which is broken down into subpackages as shown in
``thyra/cmake/Dependnecies.cmake`` (see
`<packageDir>/cmake/Dependencies.cmake`_) which is:

.. include:: ../../package_arch/UnitTests/ReducedMockTrilinos/packages/thyra/cmake/Dependencies.cmake
   :literal:

Adding in the subpackages defined in the top-level ``Thyra`` package, the full
set of `TriBITS SE Packages`_ for this project is::

  Teuchos RTOp Epetra Triutils EpetraExt ThyraCoreLibs ThyraGoodStuff \
    ThyraCrazyStuff ThyraEpetra ThyraEpetraExt Thyra

The list of `TriBITS TPLs`_ for this example project given in
``ReducedMockTrilinos/TPLsList.cmake`` (see `<repoDir>/TPLsList.cmake`_) is:

.. include:: ../../package_arch/UnitTests/ReducedMockTrilinos/TPLsList.cmake
   :literal:

Take note of the testing group (i.e. ``PT``, ``ST``, or ``EX``) assigned to
each SE package and TPL plays a signficiant role in how the TriBITS dependency
system handles enables and disable of the SE packages and TPLs.

The dependency structure of this simple TriBITS project is shown below in
`ReducedMockTrilinos Dependencies`_.

.. _ReducedMockTrilinos Dependencies:

**ReducedMockTrilinos Dependencies:**

.. include:: ../../package_arch/UnitTests/DependencyUnitTests/ReducedModelTrilinos_ExpectedDependencies.txt
   :literal:

The above dependency structure printout is produced by configuring with
``${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES=ON`` (which also results in more
dependency information than what is shown above, e.g. like computed forward
package dependencies).  Note that the top-level SE package ``Thyra`` is shown
to depend on its subpackages (not the other way around).  (Many people are
confused about this the nature of the dependencies between packages and
subpackages.)

.. ToDo: Show diagram with this dependency structure.

A number of user-setable cache variables determine what SE packages (and TPLs)
and what tests and examples get enabled.  These cache variables are described
in `Selecting the list of packages to enable`_ and are described below.  Also,
the assigned `SE Package Test Group`_ (i.e. ``PT``, ``ST``, and ``EX``) also
affects what packages get enabled or disabled.

Any of these SE packages can be enabled or disabled with
``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=(ON|OFF)`` (the default enable is
typically ``""``, see `unset enable/disable state by default`_ below).  For
``ReducedMockTrilinos``, this gives the enable/disable cache variables (with
the initial default values)::

  Trilinos_ENABLE_Teuchos=""
  Trilinos_ENABLE_RTOp""
  Trilinos_ENABLE_Epetra=""      
  Trilinos_ENABLE_Triutils=""
  Trilinos_ENABLE_EpetraExt=""
  Trilinos_ENABLE_ThyraCore=""
  Trilinos_ENABLE_ThyraGoodStuff=""
  Trilinos_ENABLE_ThyraCrazyStuff="OFF"  # Because it is 'EX'
  Trilinos_ENABLE_ThyraEpetra=""
  Trilinos_ENABLE_ThyraEpetraExt=""

Every TriBITS SE package is assumed to have tests and/or examples so TriBITS
defines the following cache varaibles as well (with the initial default
values)::

  Teuchos_ENABLE_TESTS=""
  RTOp_ENABLE_TESTS=""
  Epetra_ENABLE_TESTS=""
  Triutils_ENABLE_TESTS=""
  EpetraExt_ENABLE_TESTS=""
  ThyraCoreLibs_ENABLE_TESTS=""
  ThyraGoodStuff_ENABLE_TESTS=""
  ThyraEpetra_ENABLE_TESTS=""
  ThyraEpetraExt_ENABLE_TESTS=""
  Thyra_ENABLE_TESTS=""

NOTE: TriBITS only sets the variables ``<TRIBITS_PACKAGE>_ENABLE_TESTS`` into
the cache if the SE package ``<TRIBITS_PACKAGE>`` becomes enabled at some
point.  This cuts down the clutter in the CMake cache for large projects with
lots of packages where the user only enbles a subset of the pakages.

NOTE: TriBITS also defines the cache varaibles
``<TRIBITS_PACKAGE>_ENABLE_EXAMPLES`` for each enabled TriBITS package which
is handled the same way as the TEST variables.

Also, every defined TPL is given its own ``TPL_ENABLE_<TRIBITS_TPL>``
enable/disable cache variable.  For the TPLs in ``ReducedMockTrilinos``, this
gives the enable/disable cache variables (with default values)::

  TPL_ENABLE_MPI=""
  TPL_ENABLE_BLAS=""
  TPL_ENABLE_LAPACK=""
  TPL_ENABLE_Boost=""
  TPL_ENABLE_UMFPACK=""
  TPL_ENABLE_AMD=""
  TPL_ENABLE_PETSC=""

In addition, for every optional SE package and TPL dependency, TriBITS defines
a cache variable ``<TRIBITS_PACAKGE>_ENABLE_<OPTIONAL_DEP>``.  For the
optional dependenices shown in `ReducedMockTrilinos Dependencies`_, that gives
the additional cache variables (with default values)::

  Teuchos_ENABLE_Boost=""
  Teuchos_ENABLE_MPI=""
  Teuchos_ENABLE_Boost=""
  Epetra_ENABLE_MPI=""
  EpetraExt_ENABLE_Triutils=""
  EpetraExt_ENABLE_UMFPACK=""
  EpetraExt_ENABLE_AMD=""
  EpetraExt_ENABLE_PETSC=""
  Thyra_ENABLE_ThyraGoodStuff=""
  Thyra_ENABLE_ThyraCrazytuff=""
  Thyra_ENABLE_ThyraEpetra=""
  Thyra_ENABLE_ThyraEpetraExt=""

The above optional package-specific cache variables allow one to control
whether or not support for upstream dependency X is turned on in package Y
independent of whether or not X and Y are themselves both enabled.  For
example, if the packages ``Triutils`` and ``EpetraExt`` are both enabled, one
can explicitly disable support for the optional dependency ``Triutils`` in
``EpetraExt`` by setting ``EpetraExt_ENABLE_Triutils=OFF``.  One may want to
do this for several reasons but the bottom line is that this gives the user
more detailed control over package dependenices.  See the `TriBITS Dependency
Handling Rules/Behaviors`_ and `Explicit disable of an optional package
dependency`_ for more discussion and examples.

Before getting into specific examples for some `Standard Enable/Disable Use
Cases`_, some of the `TriBITS Dependency Handling Rules/Behaviors`_ are
defined below.

TriBITS Dependency Handling Rules/Behaviors
-------------------------------------------

Below, some of the rules and behaviors of the TriBITS dependency management
system are described.  Examples refer to the `Example ReducedMockTrilinos
Project Dependency Structure`_.  More detailed examples of these behaviors are
given in the section `Standard Enable/Disable Use Cases`_.

.. _unset enable/disable state by default:

1) An SE package ``<TRIBITS_PACKAGE>`` with testing group ``PT`` or ``ST`` is
   given an **unset enable/disable state by default**
   (i.e. ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=""``).  For example, the
   ``PT`` package ``Teuchos`` is not enabled or disabled by default and is
   given the initial value ``Trilinos_ENABLE_Teuchos=""``.  This allows
   ``PT``, and ``ST`` packages to be enabled or disabled using other logic
   defined by TriBITS which is described below.

.. _EX SE package disabled by default:

2) An SE package ``<TRIBITS_PACKAGE>`` with testing group ``EX`` is **disabled
   by default** (i.e. ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=OFF``).  This
   results in all required downstream SE packages to be disabled by default.
   However, the user can explicitly set
   ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=ON`` for an ``EX`` package and
   it will be enabled (unless one of its required dependencies are not enabled
   for some reason).

3) Any SE package ``<TRIBITS_PACKAGE>`` can be explictly **enabled** by the
   user by setting the cache variable
   ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=ON``
   (e.g. ``Trilinos_ENABLE_EpetraExt=ON``).  When a package is enabled in this
   way, the TriBITS system will try to enable all of the required upstream SE
   packages and TPLs defined by the package (specified in its
   ``Dependencies.cmake`` file).  If an enabled SE package can't be enabled
   and has to be disabled, either a warning is printed or processing will stop
   with an error (depending on the value of
   ``${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES``, see below).  In
   addition, if ``${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES=ON``, then
   TriBITS will try to enable all of the specified optional SE packages as
   well.

4) Any SE package ``<TRIBITS_PACKAGE>`` can be explictly **disabled** by the
   user by setting the cache variable
   ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=OFF``
   (e.g. ``Trilinos_ENABLE_Teuchos=OFF``).  When an SE package is explicitly
   disabled, it will result in the disable of all downstream SE packages that
   have required dependency on it.  It will also disable optional support for
   the disabled packages in downstream packages that list it as an optional
   dependency.  For an example, see `Explicit disable of a package`_.

5) A TriBITS TPL ``<TRIBITS_TPL>`` with testing group ``PT`` or ``ST`` is
   given an **unset enable/disable state by default**
   (i.e. ``TPL_ENABLE_<TRIBITS_TPL>=""``).  For example, the ``PT`` TPL
   ``BLAS`` is not enabled or disabled by default
   (i.e. ``TPL_ENABLE_BLAS=""``) .  This allows ``PT``, and ``ST`` TPLs to be
   enabled or disabled using other logic.

6) A TriBITS TPL ``<TRIBITS_TPL>`` with testing group ``EX`` (as is for `PT``
   and ``EX`` TPLs), is given an **unset enable/disable state by default**
   (i.e. ``TPL_ENABLE_<TRIBITS_TPL>""``).  This is different behavior than for
   ``EX`` SE pakages described above which provides an initial hard disable.
   However, since TriBITS will never automatically enable an optional TPL (see
   below) and since only downstream ``EX`` SE packages are allowed to have a
   required dependenices on an ``EX`` TPL (see ???), there is no need to set
   the default enable for an ``EX`` TPL to ``OFF``.

7) All TPLs listed as required TPL dependencies for the final set of enabled
   SE packages are **set to enabled** (i.e. ``TPL_ENABLE_<TRIBITS_TPL>=ON``)
   for the final set of enabled packages (unless the listed TPLs are already
   explicit disabled).  For example, if the ``Epetra`` package is enabled,
   then that will trigger the enable of its required TPLs ``BLAS`` and
   ``LAPACK``.

8) TPLs with testing group ``PT`` and ``ST`` will only be enabled if they are
   explicitly enabled by the user.  For example, just because the package
   ``Teuchos`` is enabled, the optional TPLs ``Boost`` and ``MPI`` will
   **not** be enabled by default.  To enable the optional TPL ``Boost``, for
   example, and enable support for ``Boost`` in the ``Teuchos`` package, the
   user must explicitly set ``TPL_ENABLE_Boost=ON``.

9) Any TPLs that are explicitly disabled
   (i.e. ``TPL_ENABLE_<TRIBITS_TPL>=OFF``) will result in the disable of all
   downstream dependent SE packages that have a required dependency on the
   TPL.  For example, if the user sets ``TPL_ENABLE_LAPACK=OFF``, then this
   will result in the disable of SE packages ``Teuchos`` and ``Epetra``, and
   all of the required SE packages downstream from them (see ???).  Also, the
   explicitly disabled TPL will result in disable of optional support in all
   downstream SE packages.  For example, if the user sets
   ``TPL_ENABLE_MPI=OFF``, then TriBITS will automatically set
   ``Teuchos_ENABLE_MPI=OFF`` and ``Epetra_ENABLE_MPI=OFF`` (see ???).

10) Disables trump enables where there is a conflicit and TriBITS will never
    override a disable in order to satisfy some dependency.  For example, if
    the user sets ``Trilinos_ENABLE_Teuchos=OFF`` and
    ``Trilinos_ENABLE_RTOp=ON``, then TriBITS will **not** override the
    disable of ``Teuchos`` in order to satisfy the required dependency of
    ``RTOp``.  In cases such as this, the behavior of the TriBITS dependency
    adjustment system will depend on the setting of the top-level user cache
    variable ``${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES``:

    .. _${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON:
 
    * If ``${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON``: TriBITS
      will disable explicit enable and continue on.  i.e., TriBITS will
      override ``Trilinos_ENABLE_RTOp=ON`` and set
      ``Trilinos_ENABLE_RTOp=OFF`` and print a verbose warning to STDOUT (see
      ???).

    .. _${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=OFF:

    * If ``${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=OFF``: TriBITS
      will generate a detailed error message and abort processing.  i.e.,
      TriBITS will report that ``RTOp`` is enabled but the required SE package
      ``Teuchos`` is disabled and therefore ``RTOp`` can't be enabled and
      processing must stop (see ???).

11) An explicit enable/disable of a top-level parent package with subpackages
    with ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=(ON|OFF)`` is equivanent
    to the explicit enable/disable of all of the parent package's subpakages.
    For example, explicitly setting ``Trilinos_ENABLE_Thyra=ON`` is equivalent
    to explicitly setting::

      Trilinos_ENABLE_ThyraCoreLibs=ON
      Trilinos_ENABLE_ThyraGoodStuff=ON   # Only if enabling ST code!
      Trilinos_ENABLE_ThyraEpetra=ON
      Trilinos_ENABLE_ThyraEpetraExt=ON   # Only if enabling ST code!
   
    (Note that ``Trilinos_ENABLE_ThyraCrazyStuff`` is **not** set to ``ON``
    because it is aready set to ``OFF`` by default, see `EX SE package
    disabled by default`_.)  Likewise, explicitly setting
    ``Trilinos_ENABLE_Thyra=OFF`` is equivalent to explicitly setting all of
    the ``Thyra`` subpackages to ``OFF`` at the outset.  See ??? for an
    example.

12) TriBITS will only enable an optional ``ST`` SE package when
    ``${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES=ON`` if
    ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=ON`` is also set.  If an optional
    ``ST`` upstream dependent SE package is not enabled due to
    ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=OFF``, then a one-line warning is
    printed to STDOUT.  The TriBITS default is
    ``${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES=ON``.  This helps to avoid
    problems when users try to set a permuation of enables/disables which is
    not regularly tested.

13) The user cache-variable ``${PROJECT_NAME}_ENABLE_ALL_PACKAGES`` will
    result in the enable of all ``PT`` SE packages when
    ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=OFF`` and all ``PT`` and ``ST`` SE
    packages when ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=ON``.  For an
    example, see ???.

14) Setting ``${PROJECT_NAME}_ENABLE_TESTS=ON`` will **only enable tests for
    explicitly enabled SE packages**.  For example, configuring with
    ``Trilinos_ENABLE_RTOp=ON`` and ``Trilinos_ENABLE_TESTS=ON`` will only
    result in the enable of tests for ``RTOp``, not ``Teuchos`` (even through
    TriBITS will enable ``Teuchos`` because it is a required dependency of
    ``RTOp``).  See ?? for an example.

15) TriBITS prints out (to ``cmake`` STDOUT) the initial set of
    enables/disables on input, prints a line whenever it sets (or overrides)
    an enable or disable, and prints out the final set of enables/disables.
    Therefore, the user just needs to grep the ``cmake`` STDOUT to find out
    why any particular SE package or TPL is enabled or disabled in the end.
    In addtion, will print out when tests/examples for a given SE package gets
    enabled and when support for optional SE packages and TPLs is enabled or
    not.  A detailed description of this output is given in all of the below
    examples but in particular see `Explicit enable of a package (and its
    tests)`_.

16) TriBITS setting (or overrides) of enable/disable cache variables are done
    by setting local non-cache variables at the top project-level scope.  This
    is done so they don't get set in the cache and so that the same dependency
    enable/disable logic is redone, from scratch, with each re-configure which
    results in the same enable/disable logic output as for the initial
    configure.  This is to avoid confusion by the user about why some SE
    packages and TPLs are enabled and some are not on subsequent reconfigures.

TriBITS prints out a lot of information about the enable/disable logic as it
applies the above rules/behaviors.  For a large TriBITS project with lots of
packages, this can produce a lot of output to STDOUT.  One just needs to
understand what TriBITS is printing out and where to look in the output for
different information.  The example `Standard Enable/Disable Use Cases`_ given
below show what this output looks like for the various enable/disable
scenarios and tries to explains in more detail the reasons for why the given
behavior is implemented the way that it is.  Given this output, the rule
definitions given above, and the detailed example `Standard Enable/Disable Use
Cases`_ described below, one should always be able to figure out exactly why
the final set of enables/disables is the way it is, even in the largest and
most complex of TriBITS projects.  (NOTE: The same can *not* be said for many
other large software configuration and depolyment systems where basic
decisions about what to enable and disable are hidden from the user and can be
very difficult to debug).


Standard Enable/Disable Use Cases
---------------------------------

Below, a few of the standard enable/disable use cases for a TriBITS project
are given using the `Example ReducedMockTrilinos Project Dependency
Structure`_.

* `Default configure with no packages enabled on input`_
* `Explicit enable of a package (and its tests)`_
* `Explicit disable of a package`_
* `Explicit enable of an optional package dependency`_
* `Explicit disable of an optional package dependency`_
* `Enable all packages`_

.. _Default configure with no packages enabled on input:

**Default configure with no packages enabled on input**

ToDo: Discuss parts of the printout.

ToDo: Fill in!

.. _Explicit enable of a package (and its tests):

**Explicit enable of a package (and its tests)**

ToDo: Fill in!

.. _Explicit disable of a package:

**Explicit disable of a package**

ToDo: Enable RTOp, Disable Teuchos

ToDo: Show what happens when
${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON

ToDo: Show what happens when
${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=OFF

ToDo: Fill in!

.. _Explicit enable of an optional package dependency:

**Explicit enable of an optional package dependency**

ToDo: Fill in!

.. _Explicit disable of an optional package dependency:

**Explicit disable of an optional package dependency**

ToDo: Fill in!

.. _Enable all packages:

**Enable all packages**

Enabling all packages the with::

   $ cmake -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON \
      -DTrilinos_DUMP_PACKAGE_DEPENDENCIES:BOOL=ON \
      <projectDir>

produces the dependency-related output:

.. include:: ../../package_arch/UnitTests/DependencyUnitTests/ReducedModelTrilinos_EnableAllPackages.txt
   :literal:

Just this one configure case provides a lot information and an oppurtunity to
address a number of issues.

1) The first thing to note is the set of initially enabled and disabled
packages is printed.  By default, ``PT`` and ``ST`` packages are not
explicitly enabled or disabled by default and therefore they don't appear in
the list of initially enabled or disabled packages.  This allows ``PT``, and
``ST`` packages to be enabled or disabled using other logic.  However, ``EX``
SE packages are explicited disabled by default by TriBITS.  This explicit
disabling of the ``EX`` SE packages is shown in::

  Explicitly disabled SE packages on input (by user or by default):  ThyraCrazyStuff 1

2) The option ``Trilinos_ENABLE_ALL_PACKAGES=ON`` enables only ``PT`` SE
packages.  If ``Trilinos_ENABLE_SECONDARY_TESTED_CODE=ON`` is set, then the
``ST`` packages would be enabled also.

3) Note the warnings about not enabling the SE packages ``ThyraGoodStuff`` and
``ThyraEpetraExt`` which are optional dependencies of ``Thyra`` due to
``Trilinos_ENABLE_SECONDARY_TESTED_CODE=OFF``.  These warnings are printed to
avoid confusion about why these SE packages are not enabled even through
``Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON``.

4) TriBITS will only enable required TPLs for enabled SE packages, as shown
for the TPLs ```BLAS`` and ``LAPACK`` and *not* for the optional TPLs listed.
The reason for this is that while all optional packages are present (because
their source exists in the project source tree) not all TPLs can expected to
be on a given system.  Optional TPLs are only enabled if they are explicitly
enabled.

5) The final set of enabled packages and subpackages is shown on the
configuration of the packages.









As shown above, only the "Primary Tested" ``PT`` packages are enabled by
default.

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

The basic ideas of breaking up a large set of software into pieces, defining
dependencies between the pieces, and then applying algorithms to manipulate
the dependency data-structures is nothing new.  If fact, nearly every binary
package deployment system provided in various Linux OS distributions have the
concept of packages and dependnecies and will automatically install all of the
necessary upstream dependencies when a downstream dependency install is
requested.

<Project>PackageDependencies.xml
--------------------------------

ToDo: Fill in!


TriBITS Automated Testing
=========================

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

ToDo: Define repo category Continuous, Nightly, and Experimental which also
map to CDash tracks.

.. _SE Package Test Group:

ToDo: Define SE package test group PT, ST, and EX.

ToDo: Define test category BASIC, CONTINUOUS, NIGHTLY, WEEKLY, and PERFORMANCE.

ToDo: Discuss the propery usage of these test categories and why NIGHTLY
testing should be the default.

ToDo: Fill in!


Pre-push Testing using checkin-test.py
--------------------------------------

ToDo: Describe the checkin-test.py script

ToDo: Describe the system for mapping changed files to changed packages.


TriBITS Package-by-Package CTest/Dash Driver
--------------------------------------------

ToDo: Fill in!

ToDo: Document CTEST_TEST_TIMEOUT and DART_TESTING_TIMEOUT and how these
interact.


TriBITS CDash Customizations
----------------------------

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
the following order of precedence:

1) **${REPOSITORY_NAME}_REPOSITORY_OVERRIDE_PACKAGE_EMAIL_LIST** (typically
defined in `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_): Defines a
single email address for all packages for the repository
``${REPOSITORY_NAME}`` and overrides all other package email regression
specification varaibles.  This is typically used by a meta-project to
redefined the regression email addresses for the packages in an externally
developed repository.

2) **REGRESSION_EMAIL_LIST** (defined in
`<packageDir>/cmake/Dependencies.cmake`_): Package-specific email address
specified in the packages's ``Dependencies.cmake`` file using
`TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()`_.

3) **${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESSS_BASE** (set in
`<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_): A base email address
specified at the Repository level creating package-specific email addresses
(e.g. ``<lower-case-package-name>-regression@some.repo.gov``, where
``${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESSS_BASE=some.repo.gov``).
This variable is used, for example, by the Trilinos project to provide
automatic regression email addresses for packages.

4) **${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESSS** (set in
`<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_): A single email address
for all packages specified at the Repository level
(e.g. ``my-repo-regression@some.repo.gov``).  This varaible is used for
smaller repositories with smaller development groups who just want all
regression emails for the repository's packages going to a single email
address.  This reduces the overhead of managing a bunch of individual package
email addresses but at the expense of spamming too many people with CDash
failure emails.

5) **${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESSS_BASE** (set in
`<projectDir>/cmake/ProjectDependenciesSetup.cmake`_): A base email address
specified at the Project level creating package-specific email addresses
(e.g. ``<lower-case-package-name>-regression@some.project.gov``, where
``${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESSS_BASE=some.project.gov``).  If not
already set, this variable will be set to
``${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESSS_BASE`` for the first
repostory processed that has this set. This behavior is used, for example by
the Trilinos project to automatically assign email addresses for add-on
packages and was added to maintain backward compatibility.

6) **${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESSS** (set in
`<projectDir>/cmake/ProjectDependenciesSetup.cmake`_): A single default email
address for all packages specified at the Project level
(e.g. ``my-project-regression@some.project.gov``).  If not already set, this
variable will be set to
``${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESSS`` for the first
repostory processed that has this set.  Every meta-project should set this
varaible so that it will be the defualt email address for any new package
added.

WARNING: If any of the email lists or URL string variables listed above are
set to ``"OFF"`` or ``"FALSE"`` (or some other value that CMake interprests as
false, see ``CMake Language Overivew and Gotchas``) then the varaibles are
treated as empty and not set.

If a TriBITS project does not use CDash, then no email address needed to be
assigned to packages at all (which will be the case if none of the above
variables are set).

As a general rule, repository-level settings override project-level settings
and package-level settings override both.  Also, a project can redefine a
reposiotry's regression email list settings by resetting the varibles in the
project's `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_ file.

All of the email dependency managment logic must be accessable by just running
the macro::

    TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()

The above email address configuration variables are read from the Repository
and Project files `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ and
`<projectDir>/cmake/ProjectDependenciesSetup.cmake`_, respectively.  The
``RepositoryDependenciesSetup.cmake`` files are read first in the specified
repository order followed up by reading the ``ProjectDependenciesSetup.cmake``
file.  In this way, the project can override any of the repository settings.

Here is a short review of the precedence order for how regression email
addresses are selected for a given package:

1) Package-specific email list is selected if defined (unless an override is
in place).

2) Repository-level option is selected over a project-level option.

3) Default email form with repository or project address base is selected over
single repository or project email address.

4) If none of the above are selected, then no email address is assigned.

What the above setup does is it results in the TriBITS system (in the
``TribitsCTestDriverCore.cmake`` file called under ``ctest``) creating a file
called ``CDashSubprojectDependencies.xml`` that gets sent to the CDash
server. CDash then takes this file and creates, or updates, a set of CDash
users and sets up a mapping of Labels (which are used for TriBITS package
names) to CDash user emails addresses. CDash is automatically set up to
process this XML file and create and updates CDash users. It is not, however,
set up to remove labels from existing users.  Therefore, if you change a
TriBITS package's CDash regression email list (using one of the methods
described above), then you need to manually remove the associated labels from
the old email address.  CDash will not remove them for you.

Therefore, to change the mapping of CDash regression email addresses to
TriBITS packages, you must perform the actions:

1) Change the TriBITS CMake files as described above that will result in the
desired email addresses in the ``CDashSubprojectDependeinces.xml`` file. You
can debug this by running the ``checkin-test.py`` script and seeing what gets
written in the generated `<Project>PackageDependencies.xml`_ file in the
``CHECKIN`` directory.

2) Log onto the CDash server using an administrator account and then remove
the auto-generated account for the CDash user email address for which labels
are being removed (i.e. no longer associated with a TriBITS package).  This is
needed since CDash seems to be unable to remove labels from an existing CDash
user (however this might be fixed in a current version of CDash).

3) The next time a CDash submit is performed by the
``TribitsCTestDriverCore.cmake`` script, the CDash user associated with the
mail list with labels being removed will get automatically recreated with the
right list of labels (according to the current
``CDashSubprojectDependencies.xml`` file).  Also, any new CDash users for new
email addresses will be created.

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

In this section, the typical development workflows for a TriBITS project are
described.  First, the `Basic Development Workflow`_ for a sinlge-repository
TriBITS project is described.  This is followed up with a slightly more
complex `Multi-Repository Development Workflow`_.


Basic Development Workflow
--------------------------

ToDo: Fill in!


Multi-Repository Development Workflow
-------------------------------------

ToDo: Discuss 'egdist' and the rep clone script.


Howtos
======

ToDo: Fill in!


How to Add a new TriBITS Package
--------------------------------

ToDo: Fill in!


How to Add a new TriBITS Package with Subpackages
-------------------------------------------------

ToDo: Fill in!


How to Add a new TriBITS TPL
----------------------------

ToDo: Fill in!


Additional Topics
=================

In this section, a number of miscellaneous topics and TriBITS features are
discussed.  These features and topics are either not considered primary
fetures of TriBITS (but can be very useful in many situations) or don't neatly
fit into one of the other sections.


TriBITS System Project Dependencies
-----------------------------------

The basic TriBITS system itself which is used to configure, built, test,
create tarballs, and install software that uses the TriBITS system has no
dependencies other than a basic installation of CMake (which typically
includes the exectuables ``cmake``, ``ctest``, and ``cpack``).  Great effort
has been expended to implement all of the core functionality of TriBITS just
using raw CMake.  That means that anyone who needs to configure, build, and
install the software just needs a compable CMake implementation.  TriBITS is
purposfully maintained to require an older version of CMake.  At the time of
this writing, the mimimum required version of CMake needed to use TriBITS is
CMake 2.8.1 (relasesed in March 2010, see `CMake Release Wiki
<http://www.cmake.org/Wiki/CMake_Released_Versions>`_).  CMake is becoming
iniquitous enough that many clients will already have a current-enough version
of CMake installed by default on their systems and will therefore not need to
download or install any extra software when building and installing a project
that uses TriBITS (assuming the necessary compilers etc. required by the
project are also installed).  If a current-enough version of CMake is not
installed on a given system, it is easy to download the source code and all it
needs is a basic C++ compiler to build and install.

However, note that a specific TriBITS project is free to use any newer CMake
features it wants and therefore these projects will require newer versions of
CMake than what is required by TriBITS (see discussion of
``CMAKE_MINIMUM_REQUIRED()`` in `<projectDir>/CMakeLists.txt`_).  But also
note that specific TriBITS projects and packages will also require additional
tools like compilers, Python, Perl, or many other such dependencies.  It is
just that TriBITS itself does not require any of these.  The goal of TriBITS
is not to amke the portability of software that uses it any worse than it
already is but instead to make it easier in most cases (that after all is the
whole goal of CMake).

While the core TriBITS functionality is just written using raw CMake, the more
sophisticated development tools needed to implement the full TriBITS
development environment requires Python 2.4 (or higher, but not Python 3.x).
Python is needed for tools like ``checkin-test.py`` and ``egdist``.  In
addition, these python tools are used in ``TribitsCTestDriverCore.cmake`` to
drive automated testing and submittals to CDash.  In addition, git is the
chosen version control tool for TriBITS and all of the VC related
functionality requires git support.  But none of this is required for doing
the most basic building, testing, or installation of a TriBITS project.


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


Project and Repositiory Versioning and Release Mode
----------------------------------------------------

TriBITS has built-in support for project and repository versioning and release
mode control.  When the project contains the file
`<projectDir>/Version.cmake`_, it is used to define the project's offical
version.  The idea is that when it is time to branch for a release, the only
file that needs to be changed is the file is `<projectDir>/Version.cmake`_

Each TriBITS repository can also contain a `<repoDir>/Version.cmake`_ file
that sets varaibles which TriBITS packages in that repository can use to
derive development and release version information.  If the TriBITS repository
also contains a `<repoDir>/Copyright.txt`_ file, then the information in
``<repoDir>/Version.cmake`` and ``<repoDir>/Copyright.txt`` are used to
configure a repository version header file::

  ${${REPOSITORY_NAME}_BINARY_DIR}/${REPOSITORY_NAME}_version.h

The configured header file ``${REPOSITORY_NAME}_version.h`` gives the
repository version number in several formats, which allows C/C++ code (or any
software that uses the C preprocessor) to write conditional code like::

  #if Trilinos_MAJOR_MINOR_VERSION > 100200
    /* Contains feature X */
    ...
  #else
    /* Does not contain feature X */
    ...
  #endif

Of course when the TriBITS project and the TriBITS repository are the same
directory, the ``<projectDir>/Version.cmake`` and ``<repoDir>/Version.cmake``
files are the same file, which works just fine.


TriBITS Environment Probing and Setup
-------------------------------------

Part of the TriBITS Framework is to probe the environment, set up the
compilers, and get ready to compile code.  This was mentioned in the step
"Probe and set up the environment" in `Full Processing of TriBITS Project
Files`_.  This is exectued by the TriBITS macro ``TRIBITS_SETUP_ENV()``.  Some
of to things this macro does are:

.. _Probe and set up the environment:

**Probe and set up the environment:**

* Set ``CMAKE_BUILD_TYPE``
* Set up for MPI (MPI compilers, etc.)
* Set up C, C++, and Fortran compiler
* Find Perl (sets ``PERL_EXECUTABLE``)
* Determine mixed langauge C/Fortran linking
* Set up C++11, OpenMP, and Windows issues
* Find Doxygen
* Perform some other configure-time tests (see output)

At the completion of this part of the processing, the TriBITS CMake project is
ready to compile code.


Configure-time System Tests
---------------------------

CMake has very nice support for defining configure-time checks of the system
to help inconfiguring the project.  One can check for whether a header file
exists or not, if the compiler supports a given data-type or language feature
or perform almost any other type of check that one can imagine that can be
done using the configured compilers, libraries, system tools, etc.

ToDo: Fill in!


Creating Source Distributions
-----------------------------

ToDo: Fill in!


Regulated Backward Compatibility and Deprecated Code
----------------------------------------------------

ToDo: Fill in!


Wrapping Exterally Configured/Built Software
--------------------------------------------

ToDo: Fill in!


TriBITS Dashboard Driver
------------------------

ToDo: Fill in!


References
==========

.. _SCALE:

**SCALE** http://scale.ornl.gov/


TriBITS Detailed Reference Documentation
========================================

The following subsections contain detailed reference documentation for the
various TriBITS variables and functions and macros that are used by TriBITS.


TriBITS Global Project Settings
-------------------------------

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

  If `${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON`_ (the TriBITS
  default value), then any explicitly enabled packages that have disabled
  upstream required packages or TPLs will be disabled.  If
  `${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=OFF`_, then an
  configure error will occur.  For more details also see
  `TribitsBuildQuickRef.* <../build_quick_ref/TribitsBuildQuickRef.html>`_).
  A project can define a different default value by setting::
  
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
  be added for ``ctest`` to run (see `TriBITS Automated Testing`) for
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
  keeping the cost of running the tests down.  See the section `TriBITS
  Automated Testing`_ for a more detailed discussion.

.. _${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE:

**${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE**

  The variable ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE`` switches the
  TriBITS project from development mode to release mode.  The default for this
  variable ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT`` should be set
  in the project's `<projectDir>/Version.cmake`_ file and switched from ``ON``
  to ``OFF`` when creating a release (see `Project and Repositiory Versioning
  and Release Mode`_).  When ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE = ON``,
  several other variables are given defaults appropriate for development mode.
  For example, ``${PROJECT_NAME}_ASSERT_MISSING_PACKAGES`` is set to ``ON`` by
  default in development mode but is set to ``OFF`` by default in release
  mode.  In addition, strong compiler warnings are enabled by default in
  development mode but are disabled by default in release mode.  Thsi variable
  also affects the behavior of `TRIBITS_SET_ST_FOR_DEV_MODE()`_.

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
----------------------------

The following subsections give detailed documentation for the CMake macros and
functions that make up the core TriBITS system.  These are what are used by
TriBITS project developers in their ``CMakeLists.txt`` and other files.  All
of these functions and macros should be available when processing the
project's and package's variables files if used properly.  Therefore, no
explicit ``INCLUDE()`` statements should be needed other than the initial
include of the ``TribitsProject.cmake`` file in the top-level
`<projectDir>/CMakeLists.txt`_ file so the command `TRIBITS_PROJECT()`_ can be
executed.

.. include:: TribitsMacroFunctionDoc.rst


General Utility Macros and Functions
------------------------------------

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
`TriBITS Overview`_).  This lead to the initial implementation of a scale-able
package-based architecture (PackageArch) for the Trilinos CMake project in
late 2008.  This Trilinos CMake PackageArch system evolved over the next few
years with development in the system slowing into 2010.  This Trilinos CMake
build system was then adopted as the build infrastructure for the CASL VERA
effort in 2011 where CASL VERA packages were treated as add-on Trilinos
packages (see Section `Multi-Repository Support`_).  Over the next year, there
was significant development of the system to support larger multi-repo
projects in support of CASL VERA.  That lead to the decision to formally
generalize the Trilinos CMake PackageArch build system outside of Trilinos and
the name TriBITS was formally adopted in November 2011.  Work to refactor the
Trilinos CMake system into a general reusable stand-alone CMake-based build
system started in October 2011 and an initial implementation was complete in
December 2011 when it was used for the CASL VERA build system.  In early 2012,
the ORNL CASL-related projects Denovo and SCALE (see [`SCALE`_]) adopted
TriBITS as their native development build systems.  Shortly after TriBITS was
adopted the native build system for the the CASL-related University of
Michigan code MPACT.  In addition to being used in CASL, all of these codes
also had a significant life outside of CASL.  Because they used the same
TriBITS build system, it proved relatively easy to keep these various codes
integrated together in the CASL VERA code meta-build.  At the same time,
TriBITS well served the independent development teams and non-CASL projects
independent from CASL VERA.  Since the initial extraction of TriBITS from
Trilinos, the TriBITS system was further extended and refined, driven by CASL
VERA development and expansion.  Independently, an early version of TriBITS
from 2012 was adopted by the LiveV
project\footnote{https://github.com/lifev/cmake} which was forked and extended
independently.

Note that a TriBITS "Package" is not the same thing as a "Package" in raw
CMake terminology.  In raw CMake, a "Package" is some externally provided bit
of software or other utiltiy for which the current CMake project has an
optional or required dependency.  Therefore, a raw CMake "Package" actually
maps to a `TriBITS TPL`_.  A raw CMake "Package" (e.g. Boost, CUDA, etc.)  can
be found using a standard CMake find module ``Find<rawPackageName>.cmake``
using the built-in command ``FIND_PACKAGE(<rawPackageName>)``.  It is
unfortunate that the TriBITS and the raw CMake defintions of the term
"Package" are not the same.  However, the term "Package" was coined by the
Trilinos project long ago before CMake was adopted as the Trilinos build
system and Trilinos' definition of "Package" (going back to 1998) pre-dates
the development of CMake and therefore dictated the terminology of TriBITS


Design Considerations for TriBITS
---------------------------------

ToDo: Discuss design requirements.

ToDo: Discuss why it is a good idea to explicitly list packages instead of
just searching for them.


checkin-test.py --help
----------------------

Below is a snapshot of the output from ``checkin-test.py --help``.  This
documentation contains a lot of information about the recommended development
workflow (mostly related to pushing commits) and outlines a number of
different use cases for using the script.

.. include:: checkin-test-help.txt
   :literal:


egdist --help
-------------

Below is a snapshot of the output from ``egdist --help``.  For more details on
the usage of ``egdist``, see `Multi-Repository Support`_ and `Multi-Repository
Development Workflow`_.

.. include:: egdist-help.txt
   :literal:

