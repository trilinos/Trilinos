=======================================
TriBITS Developers Guide and Reference
=======================================

:Author: Roscoe A. Bartlett (bartlettra@ornl.gov)
:Date: |date|
:Version: .. include:: TribitsGitVersion.txt

.. |date| date::

:Abstract: This document describes the usage of TriBITS to build, test, and deploy complex software.  The primary audience are those individuals who develop on a software project which uses TriBITS.  The overall structure of a TriBITS project is described including all of the various project- and package-specific files that TriBITS requires or can use and how and what order these files are processed.  It also contains detailed reference information on all of the various TriBITS macros and functions directly used in TriBITS project CMake files.  Many other topics of interest to a TriBITS project developer and architect are discussed as well.

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
and CMake code in general (i.e. using ``MESSAGE()`` print statements).  One
needs to understand how CMake defines and uses targets for various qualities
like libraries, executables, etc.  Without this basic understanding of CMake,
one will have trouble resolving problems when they occur.

.. jfrye: Below paragraph is redundant with table of contents above?

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
detailed `Howtos`_.  Later some `Additional Topics`_ are presented that don't
fit well into other sections.  Then the main bulk of the detailed reference
material for TriBITS is given in the section `TriBITS Detailed Reference
Documentation`_.  Finally, several bits of information are provided in the
`Appendix`_.

Background
==========

In order to easily find the most appropriate documentation, see the 'TriBITS 
Developer and User Roles`_ guide.  This guide describes the different roles 
that users of TriBITS may play and offers links to relevant sections of the
documentation.  Additionally, the reader may wish to review the `CMake Language 
Overview and Gotchas`_ section which is meant for users that are new to both 
CMake and TriBITS. This section gives a brief overview of getting started with
CMake and provides some warnings about non-obvious CMake behavior that often 
trips up new users of TriBITS.


TriBITS Developer and User Roles
--------------------------------

There are approximately five different types roles related to TriBITS.  These
different roles require different levels of expertise and knowledge of CMake
and knowledge of the TriBITS system.  The primary roles are 1) `TriBITS
Project User`_, 2) `TriBITS Project Developer`_, 3) `TriBITS Project
Architect`_, 4) `TriBITS System Developer`_, and 5) `TriBITS System
Architect`_.  Each of these roles builds on the necessary knowledge of the
lower-level roles.
 
.. _TriBITS Project User:
.. _TriBITS Project Users:

The first role is that of a **TriBITS Project User** who only needs to be able
to configure, build, and test a project that uses TriBITS as its build system.
A person acting in this role needs to know little about CMake other than
basics about how to run the ``cmake`` and ``ctest`` executables, how to set
CMake cache variables, and the basics of building software by typing ``make``
and running tests with ``ctest``.  The proper reference for a TriBITS Project
User is the `Project-Specific Build Reference`_.  The `TriBITS
Overview`_ document may also be of some help.  A TriBITS project user may also
need to consult `Package Dependencies and Enable/Disable Logic`_.

.. _TriBITS Project Developer:
.. _TriBITS Project Developers:

A **TriBITS Project Developer** is someone who contributes to a software
project that uses TriBITS.  They will add source files, libraries and
executables , test executables and define tests run with ``ctest``.  They have
to configure and build the project code in order to be able to develop and run
tests and therefore this role includes all of the necessary knowledge and
functions of a TriBITS Project User.  A casual TriBITS Project Developer
typically does not need to know a lot about CMake and really only needs to
know a subset of the `TriBITS Macros and Functions`_ defined in this document
in addition to the genetic `TriBITS Build Reference
<TribitsBuildReference.html>`_ document.  A slightly more
sophisticated TriBITS Project Developer will also add new packages, add new
package dependencies, and define new TPLs.  This current TriBITS Developers
Guide and Reference document should supply everything such a developer needs
to know and more.  Only a smaller part of this document needs to be understood
and accessed by people assuming this role.

.. jfrye: Add links to parts of document above?

.. _TriBITS Project Architect:
.. _TriBITS Project Architects:

The next level of roles is a **TriBITS Project Architect**.  This is someone
(perhaps only one person on a project development team) that knows the usage
and functioning of TriBITS in great detail.  They understand how to set up a
TriBITS project from scratch, how to set up automated testing using the TriBITS
system, and know how to use TriBITS to implement the overall software
development process.  A person in this role is also likely to be the one who
makes the initial technical decision for their project to adopt TriBITS for
its native build and test system.  This document (along with detailed
CMake/CTest/CDash documentation provided by Kitware and the larger community)
should provide most of what a person in this role needs to know.  A person
assuming this role is the primary audience for a lot of the more advanced
material in this document.

.. _TriBITS System Developer:
.. _TriBITS System Developers:
.. _TriBITS System Architect:
.. _TriBITS System Architects:

The last two roles **TriBITS System Developer** and **TriBITS System
Architect** are for those individuals that actually extend and modify the
TriBITS system itself.  A TriBITS System Developer needs to know how to add
new TriBITS functionality while maintaining backward compatibility, know how to
add new unit tests for the TriBITS system (see `The TriBITS Test Package`_),
and perform other related tasks.  Such a developer needs to be very
knowledgeable of the basic functioning of CMake and know how TriBITS is
implemented in the CMake language.  A TriBITS System Architect is someone who
must be consulted on almost all non-trivial changes or additions to the
TriBITS system.  A *TriBITS System Architect* in addition needs to know the
entire TriBITS system, the design philosophy that provides the foundation for
TriBITS and be an expert in CMake, CTest, and CDash.  Everything that needs to
be known by a TriBITS System Developer and a TriBITS System Architect is not
contained in this document.  Instead, the primary documentation will be in the
TriBITS CMake source code and various unit tests itself defined in `The
TriBITS Test Package`_.  At the time of this writing, there is currently there
is only one TriBITS System Architect (who also happens to be the primary
author of this document).

An explicit goal of this document is to make new `TriBITS Project Architects`_
(i.e. those who would make the decision to adopt TriBITS for their projects),
and new `TriBITS System Developers`_ to help extend and maintain TriBITS.

Depending on the particular role that a reader falls into, this document
may not be necessary and the `TriBITS Overview`_ or the
`<Project>BuildReference`_ documents may be more appropriate.  Hopefully the
above roles and discussion help the reader select the right document to
start with.


CMake Language Overview and Gotchas
-----------------------------------

TriBITS removes a lot of the boiler plate code needed to write a CMake
project.  As a result, many people can come into a project that uses TriBITS
and quickly start to contribute by adding new source files, adding new
libraries, adding new tests, and even adding new TriBITS packages and TPLs;
all without really having learned anything about CMake.  Often one can use existing
example CMake code as a guide and be successful using basic functionality. As long 
as nothing out of the ordinary happens, many people can get along just fine in this 
mode for a time.

However, we have observed that most mistakes and problems that people run into when 
using TriBITS are due to lack of basic knowledge of the CMake language.  One can find 
basic tutorials and references on the CMake language in various locations online for free.
One can also purchase the `official CMake reference book`_.  Also, documentation
for any built-in CMake command is available locally by running::

   $ cmake --help-command <CMAKE_COMMAND>

Because tutorials and detailed documentation for the CMake language already
exists, this document does not attempt to provide a first reference to CMake
(which is a large topic in itself).  However, what we try to provide below is
a short overview of the more quirky or surprising aspects of the CMake
language that a programmer experienced in another language might get tripped
up or surprised by.  Some of the more unique features of the language are
described in order to help avoid some of these common mistakes and provide
greater understanding of how TriBITS works.

.. _Official CMake reference book: http://www.cmake.org/cmake/help/book.html

The CMake language is used to write CMake projects with TriBITS. In fact the
core TriBITS functionality itself is implemented in the CMake language (see 
`TriBITS System Project Dependencies`_). CMake is a fairly simple programming 
language with relatively simple rules (for the most part).  However, compared 
to other programming languages, there are a few peculiar aspects to the CMake 
language that can make working with it difficult if you don't understand these 
rules.  For example there are unexpected variable scoping rules and how arguments 
are passed to macros and functions can be tricky. Also, CMake has some interesting 
gotchas.  In order to effectively use TriBITS (or just raw CMake) to construct 
and maintain a project's CMake files, one must know the basic rules of CMake 
and be aware of these gotchas.

The first thing to understand about the CMake language is that nearly every
line of CMake code is just a command taking a string (or an array of strings)
and functions that operate on strings.  An array argument is just a single
string literal with elements separated by semi-colons ``"<str0>;<str1>;..."``.
CMake is a bit odd in how it deals with these arrays, which are just
represented as a string with elements separated with semi-colons ``';'``.  For
example, all of the following are equivalent and pass in a CMake array with 3
elements [``A``], [``B``], and [``C``]::

  SOME_FUNC(A B C)
  SOME_FUNC("A" "B" "C")
  SOME_FUNC("A;B;C")

However, the above is *not* the same as::

  SOME_FUNC("A B C")

which just passes in a single element with value [``A B C``].  Raw quotes in
CMake basically escape the interpretation of space characters as array element
boundaries.  Quotes around arguments with no spaces does nothing (as seen
above, except for the interpretation as variable names in an ``IF()``
statement).  In order to get a quote char [``"``] into string, you must escape
it as::

  SOME_FUNC(\"A\")

which passes an array with the single argument [``\"A\"``].

Variables are set using the built-in CMake ``SET()`` command that just takes
string arguments like::

  SET(SOME_VARIABLE "some_value")

In CMake, the above is identical, in every way, to::

  SET(SOME_VARIABLE some_value)
  SET("SOME_VARIABLE";"some_value")
  SET("SOME_VARIABLE;some_value")

The function ``SET()`` simply interprets the first argument to as the name of
a variable to set in the local scope.  Many other built-in and user-defined
CMake functions work the same way.  That is, some of the string arguments are
interpreted as the names of variables.  There is no special language feature
that interprets them as variables (except in an ``IF()`` statement).

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
control structures (i.e. see ``cmake --help-command if`` and ``cmake
--help-command foreach``).

CMake offers a rich assortment of built-in commands for doing all sorts of
things.  Two of these are the built-in ``MACRO()`` and the ``FUNCTION()``
commands which allow you to create user-defined macros and functions. TriBITS
is actually built on CMake functions and macros.  All of the built-in and 
user-defined macros, and some functions take an  array of string arguments.  
Some functions take in positional arguments. In fact,  most functions take a 
combination of positional and keyword arguments.

Variable names are translated into their stored values using
``${SOME_VARIABLE}``.  The value that is extracted depends on if the variable
is set in the local or global (cache) scope.  The local scopes for CMake start
in the base project directory in its base ``CMakeLists.txt`` file.  Any
variables that are created by macros in that base local scope are seen across
an entire project but are *not* persistent across multiple successive
``cmake`` configure invocations where the cache file ``CMakeCache.txt`` is not
deleted in between.

The handling of variables is one area where CMake is radically different from
most other languages.  First, a variable that is not defined simply returns
nothing.  What is surprising to most people about this is that it does not
even return an empty string that would register as an array element!  For
example, the following set statement::

   SET(SOME_VAR a ${SOME_UNDEFINED_VAR} c)

(where ``SOME_UNDEFINED_VAR`` is an undefined variable) produces
``SOME_VAR='a;c'`` and *not* ``'a;;c'``!  The same thing occurs when an empty
variable is de-references such as with::

   SET(EMPTY_VAR "")
   SET(SOME_VAR a ${EMPTY_VAR} c)

which produces ``SOME_VAR='a;c'`` and *not* ``'a;;c'``.  In order to always
produce an element in the array even if the variable is empty, one must quote
the argument as with::

   SET(EMPTY_VAR "")
   SET(SOME_VAR a "${EMPTY_VAR}" c)

which produces ``SOME_VAR='a;;c'``, or three elements as one might assume.

This is a common error that people make when they call CMake functions
(built-in or TriBITS-defined) involving variables that might be undefined or
empty.  For example, for the macro::

   MACRO(SOME_MACRO  A_ARG  B_ARG  C_ARG)
      ...
   ENDMACRO()

if someone tries to call it with (misspelled variable?)::

  SOME_MACRO(a ${SOME_OHTER_VAR} c)

and if ``SOME_OHTER_VAR=""`` or if it is undefined, then CMake will error out
with the error message saying that the macro ``SOME_MACRO()`` takes 3
arguments but only 2 were provided.  If a variable might be empty but that is
still a valid argument to the command, then it must be quoted as::

  SOME_MACRO(a "${SOME_OHTER_VAR}" c)

Related to this problem is that if you misspell the name of a variable in a
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

In this case, the misspelled variable would be caught.

While on the subject of ``IF()`` statements, CMake has a strange convention.
When you say::

  IF (SOME_VAR)
    DO_SOMETHING()
  ENDIF()

then ``SOME_VAR`` is interpreted as a variable and will be considered true and
``DO_SOMETHING()`` will be called if ``${SOME_VAR}`` does *not* evaluate to
``0``, ``OFF``, ``NO``, ``FALSE``, ``N``, ``IGNORE``, ``""``, or ends in the
suffix ``-NOTFOUND``.  How about that for a true/false rule!  To be safe, use
``ON/OFF`` and ``TRUE/FALSE`` pairs for setting variables.  Look up native
CMake documentation on ``IF()`` for all the interesting details and all the
magical things it can do.

**WARNING:** If you mistype ``"ON"`` as ``"NO"``, it evaluates to
``FALSE``/``OFF``!  (That is a fun defect to track down!)

CMake language behavior with respect to case sensitivity is also strange:

* Calls of built-in and user-defined macros and functions is *case
  insensitive*!  That is ``set(...)``, ``SET(...)``, ``Set()``, and all other
  combinations of upper and lower case characters for 'S', 'E', 'T' all call
  the built-in ``SET()`` function.  The convention in TriBITS is to use all
  caps for functions and macros (which was adopted by following the
  conventions used in the early versions of TriBITS, see the `History of
  TriBITS`_).  The convention in CMake literature from Kitware seems to use
  lower-case letters for functions and macros.

* However, the names of CMake (local or cache/global) variables are *case
  sensitive*!  That is, ``SOME_VAR`` and ``some_var`` are *different*
  variables.  Built-in CMake variables tend use all caps with underscores
  (e.g. ``CMAKE_CURRENT_SOURCE_DIR``) but other built-in CMake variables tend
  to use mixed case with underscores (e.g. ``CMAKE_Fortran_FLAGS``).  TriBITS
  tends to use a similar naming convention where variables have mostly
  upper-case letters except for parts that are proper nouns like the project,
  package or TPL name (e.g. ``TribitsProj_TRIBITS_DIR``,
  ``TriBITS_SOURCE_DIR``, ``Boost_INCLUDE_DIRS``).

I don't know of any other programming language that uses different case
sensitivity rules for variables and functions.  However, because we must parse
macro and function arguments when writing user-defined macros and functions,
it is a good thing that CMake variables are case sensitive.  Case insensitivity
would make it much harder and more expensive to parse argument lists that take
keyword-based arguments.

Other mistakes that people make result from not understanding how CMake scopes
variables and other entities.  CMake defines a global scope (i.e. "cache"
variables) and several nested local scopes that are created by
``ADD_SUBDIRECTORY()`` and entering FUNCTIONS.  See `DUAL_SCOPE_SET()`_ for a
short discussion of these scoping rules.  And it is not just variables that
can have local and global scoping rules.  Other entities, like defines set
with the built-in command ``ADD_DEFINITIONS()`` only apply to the local scope
and child scopes.  That means that if you call ``ADD_DEFINITIONS()`` to set a
define that affects the meaning of a header-file in C or C++, for example,
that definition will *not* carry over to a peer subdirectory and those
definitions will not be set (see warning in `Miscellaneous Notes
(TRIBITS_ADD_LIBRARY())`_).

Now that some CMake basics and common gotchas have been reviewed, we now get
into the meat of TriBITS starting with the overall structure of a TriBITS
project in the next section.

Software Engineering Packaging Principles
-----------------------------------------

The term "software engineering" package is adopted in TriBITS nomenclature in
order to highlight and the role of defining and managing "Packages" according
to standard software engineering packaging principles.  In his book [`Agile
Software Development, 2003`_], Robert Martin defines several object-oriented
(OO) software engineering principles related to packaging software which are
listed below:

* *Package Cohesion OO Principles*:

  1) *REP (Release-Reuse Equivalency Principle)*: The granule of reuse is the
     granule of release.

  2) *CCP (Common Closure Principle)*: The classes in a package should be
     closed together against the same kinds of changes.  A change that affects
     a closed package affects all the classes in that package and no other
     packages.

  3) *CRP (Common Reuse Principle)*: The classes in a package are used
     together.  If you reuse one of the classes in a package, you reuse them
     all.

* *Package Coupling OO Principles*:

  4) *ADP (Acyclic Dependencies Principle)*: Allow no cycles in the package
     dependency graph.

  5) *SDP (Stable Dependencies Principle)*: Depend in the direction of
     stability.

  6) SAP (Stable Abstractions Principle): A package should be as abstract as
     it is stable.

Any of these six OO packaging principles (and other issues) may be considered
when deciding how to partition software into different `TriBITS SE Packages`_.

NOTE: The purpose of this TriBITS Developers Guide document is not teach
basic software engineering so these various principles will not be expanded on
further.  However, interested readers are strongly encouraged to read
[`Agile Software Development, 2003`_] as one of the better software
engineering books out there (see
http://web.ornl.gov/~8vt/readingList.html#Most_Recommended_SE_Books).

TriBITS Project Structure
=========================

.. jfrye: What does the "F" in framework mean?

TriBITS is a Framework, implemented in CMake, to create CMake projects.  As a
`Software Framework`_, TriBITS defines the overall structure of a CMake build
system for a project and it processes the various project-, repository-, and
package-specific files in a specified order.  Almost all of this processing
takes place in the `TRIBITS_PROJECT()`_ macro (or macros and functions it
calls).  The following subsections define the essence of the TriBITS framework
in some detail.  Later sections cover specific topics and the various sections
link to each other.  Within this section, the subsection `TriBITS Structural
Units`_ defines the basic units `TriBITS Project`_, `TriBITS Repository`_,
`TriBITS Package`_, `TriBITS TPL`_ and other related structural units.  The
subsection `Processing of TriBITS Files: Ordering and Details`_ defines
exactly what files TriBITS processes and in what order.  It also shows how to
get TriBITS to show exactly what files it is processing to help in debugging
issues.  The subsection `Coexisting Projects, Repositories, and Packages`_
gives some of the rules and constrains for how the different structure units
can co-exist in the same directories.  The last two subsections in this
section cover `Standard TriBITS TPLs`_ and `Common TriBITS TPLs`_.

TriBITS Structural Units
------------------------

A CMake project that uses TriBITS as its build and test system is composed of
a single `TriBITS Project`_, one or more `TriBITS Repositories`_ and one or
more `TriBITS Packages`_.  In addition, a TriBITS Package can be broken up
into `TriBITS Subpackages`_.  Together, the collection of TriBITS Packages and
TriBITS Subpackages are called *TriBITS Software Engineering Packages*, or
`TriBITS SE Packages`_ for short.

First, to better establish the basic nomenclature, the key structural TriBITS
units are:

.. _TriBITS Packages:

* `TriBITS Package`_: A collection of related software that typically includes
  one or more source files built into one or more libraries and has associated
  tests to help define and protect the functionality provided by the software.
  A package also typically defines a single integrated unit of documentation
  and testing (see `TriBITS Automated Testing`_).  A TriBITS package may or
  may not be broken down into multiple subpackages. Examples of TriBITS
  packages in `TribitsExampleProject`_ include ``SimpleCXX``,
  ``MixedLang`` and ``WithSubpackages``.  (Don't confuse a TriBITS
  "Package" with a raw CMake "Package". A raw CMake "Package" actually maps to
  a `TriBITS TPL`_; see `Why a TriBITS Package is not a CMake Package`_.)

.. _TriBITS Subpackages:

* `TriBITS Subpackage`_: A part of a parent `TriBITS Package`_ that also
  typically has source files built into libraries and tests but is documented
  and tested along with the other subpackages in the parent package.  The primary
  purpose for supporting subpackages is to provide finer-grained control of
  software dependencies.  In `TribitsExampleProject`_, ``WithSubpackages`` is
  an example of a package with subpackages ``'A'``, ``'B'``, and ``'C'``.  The
  full subpackage name is prefixed by the parent package name (e.g. the full
  name for subpackage ``'A'`` is ``WithSubpackagesA``).  The parent package is
  always implicitly dependent on its subpackages (e.g. parent package
  ``WithSubpackages`` depends on its subpackages ``WithSubpackagesA``,
  ``WithSubpackagesB``, and ``WithSubpackagesC``).

.. _TriBITS SE Package:
.. _TriBITS SE Packages:

* **TriBITS SE Package**: The combined set of `TriBITS Packages`_ and `TriBITS
  Subpackages`_ that constitute the basic *Software Engineering* packages of a
  TriBITS project (see `Software Engineering Packaging Principles`_): SE
  packages are the basis for setting dependencies in the TriBITS system.  For
  example, the SE Packages provided by the top-level example package
  ``WithSubpackages`` (in order of increasing dependencies) are
  ``WithSubpackagesA``, ``WithSubpackagesB``, ``WithSubpackagesC``, and
  ``WithSubpackages`` (see `TribitsExampleProject`_).

.. _TriBITS TPLs:

* `TriBITS TPL`_: The specification for a particular external dependency that
  is required or can be used in one or more `TriBITS SE Packages`_.  A TPL (a
  Third Party Library) typically provides a list of libraries or a list
  include directories for header files but can also be manifested in order
  ways as well.  Examples of basic TPLs include ``BLAS``, ``LAPACK``, and
  ``Boost``.

.. _TriBITS Repositories:

* `TriBITS Repository`_: A collection of one or more `TriBITS Packages`_
  specified in a `<repoDir>/PackagesList.cmake`_ file and zero or more TPL
  declarations specified in a `<repoDir>/TPLsList.cmake`_ file.
  As discussed below, Repositories can include Native Repositories and Version
  Control (VC) Repositories.

.. _TriBITS Projects:

* `TriBITS Project`_: A collection of `TriBITS Repositories`_ and `TriBITS
  Packages`_ that defines a complete CMake ``PROJECT`` defining software which can be
  directly configured with ``cmake`` and then be built, tested, installed,
  etc.

.. _TriBITS Meta-Project:
.. _TriBITS Meta-Projects:

* **TriBITS Meta-Project**: A `TriBITS Project`_ that contains no native
  `TriBITS packages`_ or `TriBITS TPLs`_ but is composed out packages from
  other `TriBITS Repositories`_.

In this document, dependencies are described as either being *upstream* or
*downstream/forward* defined as:

.. _upstream:
.. _upstream dependency:
.. _upstream dependencies:

* If unit "B" requires (or can use, or comes before) unit "A", then "A" is an
  **upstream dependency** of "B".

.. _downstream:
.. _downstream dependency:
.. _downstream dependencies:

* If unit "B" requires (or can use, or comes before) unit "A", then "B" is a
  **downstream dependency** or a **forward dependency** of "A".

The following subsections define the major structural units of a TriBITS
project in more detail.  Each structural unit is described along with the
files and directories associated with each.  In addition, a key set of TriBITS
CMake variables for each are defined as well.

In the next major section following this one, some `Example TriBITS Projects`_
are described.  For those who just want to jump in and learn best by example,
these example projects are a good way to start.  These example projects will
be referenced in the more detailed descriptions given in this document.

The last issue to touch on before getting into the detailed descriptions of
the different TriBITS structural units is the issue of how CMake variables are
defined and used by TriBITS. The CMake variables described in the TriBITS
structural units below fall into one of two major types:

* *Local Fixed-Name Variables* are used temporarily in the processing of a
  TriBITS unit.  These include variables such as ``PROJECT_NAME``,
  ``REPOSITORY_NAME``, ``PACKAGE_NAME``, and ``PARENT_PACKAGE_NAME``.  These
  are distinguished by having a non-namespaced fixed/constant name.  They are
  typically part of TriBITS reflection system, allowing subordinate units to
  determine the encapsulating unit in which they are participating.  For,
  example, a TriBITS subpackage can determine its name, its parent package's
  name and directories, its parent repository name and directories, and the
  enclosing project's name and directories without having to refer to any
  specific names.

* *Globally Scoped Namespaced Variables* are used to refer to properties of a
  named TriBITS unit that are seen globally by the entire TriBITS project.
  These include variables such as ``${REPOSITORY_NAME}_SOURCE_DIR``
  (e.g. ``TribitsExProj_SOURCE_DIR``) and ``${PACKAGE_NAME}_BINARY_DIR``
  (e.g. ``SimpleCXX_BINARY_DIR``).  They are available after processing the
  unit, for use by `downstream`_ or subordinate units.  They are part of the
  TriBITS dependency system, allowing downstream units to access properties of
  their known `upstream dependencies`_.

More information about these various files is described in section `Processing
of TriBITS Files: Ordering and Details`_.


TriBITS Project
+++++++++++++++

A TriBITS Project:

* Defines the ``PROJECT_NAME`` CMake variable (defined in
  `<projectDir>/ProjectName.cmake`_)
* Defines a complete CMake project which calls ``PROJECT(${PROJECT_NAME}
  ...)`` and can be directly configured, built, tested, installed, etc.
* Consists of one or more `TriBITS Repositories`_ (and may itself be a
  `TriBITS Repository`_) which can include native and extra repositories.
* Allows for extra Repositories to be added on before or after the set of
  native Repositories (specified in
  `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ or by CMake cache
  variables)
* Defines a default CDash server and default project name on the server (the
  project name on the CDash server must be the same as ``${PROJECT_NAME}``).
* Defines pre-push testing standard builds (using `checkin-test.py`_)

For more details on the definition of a TriBITS Project, see:

* `TriBITS Project Core Files`_
* `TriBITS Project Core Variables`_

TriBITS Project Core Files
..........................

The core files making up a TriBITS Project (where ``<projectDir> =
${PROJECT_SOURCE_DIR}``) are::

  <projectDir>/
     ProjectName.cmake    # Defines PACKAGE_NAME
     CMakeLists.txt       # Base project CMakeLists.txt file
     CTestConfig.cmake    # [Optional] Needed for CDash submits
     Version.cmake        # [Optional] Dev mode, Project version, VC branch
     project-checkin-test-config.py # [Optional] checkin-test.py config
     cmake/
       NativeRepositoriesList.cmake    # [Optional] Rarely used
       ExtraRepositoriesList.cmake     # [Optional] Lists repos and VC URLs 
       ProjectCiFileChangeLogic.py     # [Optional] CI global change/test logic
       ProjectCompilerPostConfig.cmake # [Optional] Override/tweak build flags
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
* `<projectDir>/project-checkin-test-config.py`_
* `<projectDir>/cmake/NativeRepositoriesList.cmake`_
* `<projectDir>/cmake/ExtraRepositoriesList.cmake`_
* `<projectDir>/cmake/ProjectCiFileChangeLogic.py`_
* `<projectDir>/cmake/ProjectCompilerPostConfig.cmake`_
* `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_
* `<projectDir>/cmake/CallbackDefineProjectPackaging.cmake`_
* `<projectDir>/cmake/tribits/`_
* `<projectDir>/cmake/ctest/CTestCustom.cmake.in`_

.. _<projectDir>/ProjectName.cmake:

**<projectDir>/ProjectName.cmake**: [Required] At a minimum provides a
``SET()`` statement to set the local variable ``PROJECT_NAME``.  This file is
the first file that is read by a number of tools in order to get the TriBITS
project's name.  This file is read first in every context that involves
processing the TriBITS project's files, including processes and tools that
just need to build the package and TPL dependency tree (see `Reduced Package
Dependency Processing`_).  Being this is the first file read in for a TriBITS
project and that it is read in first at the top level scope in every context,
this is a good file to put in other universal static project options.  Note
that this is a project, not a repository file so no general
repository-specific settings should go in this file.  A simple example of this
file is `TribitsExampleProject`_/``ProjectName.cmake``:

.. include:: ../../examples/TribitsExampleProject/ProjectName.cmake
   :literal:

A meta-project's ``ProjectName.cmake`` file might have a number of other
variables set that define how the various different TriBITS repos are cobbled
together into a single TriBITS meta-project (usually because the different
TriBITS repos and packages are a bit messy and have other issues).  For
example, the CASL VERA TriBITS meta-project at one point had a very extensive
collection of set statements in this file.

.. _<projectDir>/CMakeLists.txt:

**<projectDir>/CMakeLists.txt**: [Required] The top-level CMake project file.
This is the first file that the ``cmake`` executable processes that starts
everything off and is the base-level scope for local (non-cache) CMake
variables.  Due to a few CMake limitations and quirks, a project's top-level
``CMakeLists.txt`` file is not quit as clean as one might otherwise hope would
be but it is not too bad.  A simple, but representative, example is
`TribitsExampleProject`_/``CMakeLists.txt``:

.. include:: ../../examples/TribitsExampleProject/CMakeLists.txt
   :literal:

A couple of CMake and TriBITS quirks that the above example
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

In this example project, a CMake cache variable `${PROJECT_NAME}_TRIBITS_DIR`_
must be set by the user to define where the base ``tribits`` source directory
is located.  With this variable set (i.e. passed into ``cmake`` command-line
use ``-DTribitsExProj_TRIBITS_DIR=<someDir>``), one just includes a single
file to pull in the TriBITS system::

  INCLUDE("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")

With the ``TriBITS.cmake`` file included, the configuration of the project
using TriBITS occurs with a single call to `TRIBITS_PROJECT()`_.

Some projects, like Trilinos, actually snapshot the `TriBITS/tribits/`_
directory into their source tree `<projectDir>/cmake/tribits/`_ and therefore
don't need to have this variable set.  In Trilinos, the include line is just::

  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits/TriBITS.cmake)

The minimum CMake version must also be declared in the top-level
``CMakeLists.txt`` file as shown.  Explicitly setting the minimum CMake
version avoids strange errors that can occur when someone tries to build the
project using a version of CMake that is too old.  The project should set the
minimum CMake version based on the CMake features used in that project's own
CMake files.  The minimum CMake version required by TriBITS is defined in in
the variable ``TRIBITS_CMAKE_MINIMUM_REQUIRED`` (the current minimum version
of CMake required by TriBITS is given at in `Getting set up to use CMake`_) .
For example, the ``VERA/CMakeLists.txt`` file lists as its first line::

  SET(VERA_TRIBITS_CMAKE_MINIMUM_REQUIRED 3.10.0)
  CMAKE_MINIMUM_REQUIRED(VERSION ${VERA_TRIBITS_CMAKE_MINIMUM_REQUIRED}
    FATAL_ERROR)

.. _<projectDir>/CTestConfig.cmake:

**<projectDir>/CTestConfig.cmake**: [Optional] Specifies the CDash site and
project to submit results to when doing an automated build driven by the CTest
driver function `TRIBITS_CTEST_DRIVER()`_ (see `TriBITS CTest/CDash Driver`_).
This file is also required to use the TriBITS-generated ``dashboard`` target
(see `Dashboard Submissions`_).  An example of this file is
`TribitsExampleProject`_/``CTestConfig.cmake``:

.. include:: ../../examples/TribitsExampleProject/CTestConfig.cmake
   :literal:

Most of the variables set in this file are directly understood by raw
``ctest`` and those variables not be explained here further (see documentation
for the standard CMake module ``CTest``).  The usage of the function
`SET_DEFAULT_AND_FROM_ENV()`_ allows the variables to be overridden both as
CMake cache variables and in the environment.  The latter is needed when
running using ``ctest`` as the driver (since older versions of ``ctest`` did
not support ``-D<var-name>:<type>=<value>`` command-line arguments like for
``cmake``).  Given that all of these variables are nicely namespaced,
overriding them in the shell environment is not as dangerous as might
otherwise be the case but this is what had to be done to get around
limitations for older versions of CMake/CTest.

NOTE: One can also set::

  SET_DEFAULT_AND_FROM_ENV(TRIBITS_2ND_CTEST_DROP_SITE ...)
  SET_DEFAULT_AND_FROM_ENV(TRIBITS_2ND_CTEST_DROP_LOCATION ...)

in this file in order to submit to a second CDash site/location.  For details,
see `Dashboard Submissions`_.  This is useful when considering a CDash upgrade
and/or implementing new CDash features or tweaks.

.. _<projectDir>/Version.cmake:

**<projectDir>/Version.cmake**: If defined, gives the project's version and
determines development/release mode (see `Project and Repository Versioning
and Release Mode`_).  This file is read in (using ``INCLUDE()``) in the
project's base-level ``<projectDir>/CMakeLists.txt`` file scope so local
variables set in this file are seen by the entire CMake project.  For example,
`TribitsExampleProject`_/``Version.cmake``, looks like:

.. include:: ../../examples/TribitsExampleProject/Version.cmake
   :literal:

When this file exists in the base project, these will be used to create
standard SOVERSION symlinks to shared libs.  For example, on Linux, in
addition to the real shared lib ``lib<libname>.so``, the standard SOVERSION
symlinks are created like::

  lib<libname>.so.01
  lib<libname>.so.1.1

When this file exists at the repository level, the prefix
``${REPOSITORY_NAME}_`` is used instead of hard-coding the project name.  This
is so that the same ``Version.txt`` file can be used as the
`<repoDir>/Version.cmake`_ file and have the repository name be flexible.
TriBITS sets ``REPOSITORY_NAME = ${PROJECT_NAME}`` when it reads in this file
at the project-level scope.

It is strongly recommended that every TriBITS project contain a
``Version.cmake`` file, even if a release has never occurred.  Otherwise, the
project needs to define the variable
`${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT`_ at the global project scope
(perhaps in ``<projectDir>/ProjectName.cmake``) to get right development mode
behavior.

.. _<projectDir>/project-checkin-test-config.py:

**<projectDir>/project-checkin-test-config.py**: [Optional] Used to define the
``--default-builds`` and other project-level configuration options for the
project's usage of the `checkin-test.py`_ tool.  Machine or package-specific
options should **not** be placed in this file.  An example of this file for
`TribitsExampleProject`_/``project-checkin-test-config.py`` is shown below:

.. include:: ../../examples/TribitsExampleProject/project-checkin-test-config.py
   :literal:

The contents of the file ``project-checkin-test-config.py`` show above are
pretty self explanatory.  This file defines a single python dictionary
data-structure called ``configuration`` which gives some default arguments in
``defaults``, and then cmake options that define the projects
``--default-builds``.  For more details, see the section `Pre-push Testing
using checkin-test.py`_.

.. _<projectDir>/cmake/NativeRepositoriesList.cmake:

**<projectDir>/cmake/NativeRepositoriesList.cmake**: [Deprecated] If present,
this file gives the list of native repositories for the TriBITS project.  The
file must contain a ``SET()`` statement defining the variable
``${PROJECT_NAME}_NATIVE_REPOSITORIES`` which is just a flat list of repository
names that must also be directory names under ``<projectDir>/``.  For example,
if this file contains::

  SET(${PROJECT_NAME}_NATIVE_REPOSITORIES Repo0 Repo1)

then the directories ``<projectDir>/Repo0/`` and ``<projectDir>/Repo1/`` must
exist and must be valid TriBITS repositories (see `TriBITS Repository`_).

There are no examples for the usage of this file in any of the TriBITS
examples or test projects.  However, support for this file is maintained for
backward compatibility since there may be some TriBITS projects that still use
it.  It is recommended instead to define multiple repositories using the
`<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file as it allows for more
flexibility in how extra repositories are specified and how they are accessed.
The latter file allows the various tools to perform version control (VC)
activities with these repos while "native repositories" do not.

If this file ``NativeRepositoriesList.cmake`` does not exist, then TriBITS
sets ``${PROJECT_NAME}_NATIVE_REPOSITORIES`` equal to ".", or the base project
directory (i.e. ``<projectDir>/.``).  In this case, the file
``<projectDir>/PackagesList.cmake`` and ``<projectDir>/TPLsList.cmake`` must
exist.  However, if the project has no native packages or TPLs, then these
files can be set up to list no packages or TPLs.  This is the case for
meta-projects like VERA that have only extra repositories specified in the
file `<projectDir>/cmake/ExtraRepositoriesList.cmake`_.

.. _<projectDir>/cmake/ExtraRepositoriesList.cmake:

**<projectDir>/cmake/ExtraRepositoriesList.cmake**: [Optional] If present,
this file defines a list of extra repositories that are added on to the
project's native repositories.  The list of repositories is defined using the
macro `TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES()`_.  For example, the extra
repos file:

.. include:: ExtraReposList.cmake
   :literal:

shows the specification of both TriBITS Repositories and non-TriBITS VC
Repositories.  In the above file, the repositories ``ExtraRepo1``,
``ExtraRepo3``, and ``ExtraRepo4`` are both TriBITS and VC repositories that
are cloned into directories under ``<projectDir>`` of the same names from the
URLs ``someurl.com:/ExtraRepo1``, ``someurl3.com:/ExtraRepo3``, and
``someurl4.com:/ExtraRepo4``, respectively.  However, the repository
``ExtraRepo2`` is **not** a `TriBITS Repository`_ because it is marked as
``NOPACKAGES``.  In this case, it gets cloned as the directory::

  <projectDir>/packages/SomePackage/Blah

However, the code in the tools `checkin-test.py`_ and
`TRIBITS_CTEST_DRIVER()`_ will consider non-TriBITS VC repos like
``ExtraRepo2`` and any changes to this repository will be listed as changes to
``somePackage`` (see `Pre-push Testing using checkin-test.py`_).

NOTE: This file can be overridden by setting the cache variable
`<Project>_EXTRAREPOS_FILE`_.

.. ToDO: I am not really sure what repos get cloned the first time based on
.. the list of extra repos.  From looking at TribitsCTestDriverCore.cmake, it
.. looks like only the selected repos will be cloned.  I need to add some unit
.. tests that really show what the real behavior is and then document that
.. behavior here.


.. _<projectDir>/cmake/ProjectCiFileChangeLogic.py:

**<projectDir>/cmake/ProjectCiFileChangeLogic.py**: [Optional] If present,
then this Python module is imported and the Python class defined there
ProjectCiFileChangeLogic there is used to determine which files need to
trigger a global rebuild of the project enabling all packages.

An example of this given in the file
``TribitsExampleProject/cmake/ProjectCiFileChangeLogic.py``:

.. include:: ../../examples/TribitsExampleProject/cmake/ProjectCiFileChangeLogic.py
   :literal:

This logic is used in all code that is used in CI testing including
`checkin-test.py`_, `TRIBITS_CTEST_DRIVER()`_ and
`get-tribits-packages-from-files-list.py`_.  If this file does not exist, then
TriBITS has some default logic which may or may not be sufficient for the
needs of a given project.


.. _<projectDir>/cmake/ProjectCompilerPostConfig.cmake:

**<projectDir>/cmake/ProjectCompilerPostConfig.cmake**: [Optional] If present,
then this file is read using ``INCLUDE()`` at the top-level CMakeLists.txt
file scope right after the compilers for the languages ``<LANG>`` = ``C``,
``CXX``, and ``Fortran`` are determined and checked using
``ENABLE_LANGUAGE(<LANG>)`` but before any other checks are performed.  This
file can contain logic for the project to adjust the flags set in
``CMAKE_<LANG>_FLAGS`` and changes to other aspects of the build flags
(including link flags, etc.).

One example of the usage of this file is the Trilinos project where this file
is (or was) used to apply specialized logic implemented in the Kokkos build
system to select compiler options and to determine how C++11 and OpenMP flags
are set.  This file in Trilinos looked like::

  IF (${Trilinos_ENABLE_Kokkos})

    ...

    include(${Kokkos_GEN_DIR}/kokkos_generated_settings.cmake)
  
    IF (NOT KOKKOS_ARCH STREQUAL "None")
    
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${KOKKOS_CXX_FLAGS}")
    
      MESSAGE("-- " "Skip adding flags for C++11 because Kokkos flags does that ...")
      SET(${PROJECT_NAME}_CXX11_FLAGS " ")
    
      MESSAGE("-- " "Skip adding flags for OpenMP because Kokkos flags does that ...")
      SET(OpenMP_CXX_FLAGS_OVERRIDE " ")
    
    ENDIF()
  
  ENDIF()

The exact context where this file is processed (if it exists) is described in
`Full Processing of TriBITS Project Files`_ and `TriBITS Environment Probing
and Setup`_.


.. _<projectDir>/cmake/ProjectDependenciesSetup.cmake:

**<projectDir>/cmake/ProjectDependenciesSetup.cmake**: [Optional] If present,
this file is included a single time as part of the generation of the project's
dependency data-structure (see `Reduced Package Dependency Processing`_).  It
gets included at the top project level scope after all of the
`<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ files have been included
but before all of the package `<packageDir>/cmake/Dependencies.cmake`_ files
are included.  Any local variables set in this file have project-wide scope.
The primary purpose for this file is to set variables that will impact the
processing of project's package ``Dependencies.cmake`` files.

The typical usage of this file is to set the default CDash email address for
all packages or override the email addresses for all of a repository's package
CDash regression email addresses (see `CDash regression email addresses`_).
For example, to set the default email address for all of the packages, one
would set in this file::

   SET_DEFAULT(${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESS
       projectx-regressions@somemailserver.org)

The repository email address variables
`${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESS_BASE`_ and
`${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESS`_ possibly set in the just
processed `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ files can also
be overridden in this file.  The CASL VERA meta-project uses this file to
override several of the repository-specific email addresses for its
constituent repositories.

In general, variables that affect how package dependencies are defined or
affect package and TPL enable/disable logic for *only this particular project*
should be defined in this file.

.. _<projectDir>/cmake/CallbackDefineProjectPackaging.cmake:

**<projectDir>/cmake/CallbackDefineProjectPackaging.cmake**: [Optional] If
exists, defines the CPack settings for the project (see `Official CPack
Documentation <http://www.cmake.org/cmake/help/documentation.html>`_ and
`Online CPack Wiki <http://www.cmake.org/Wiki/CMake:Packaging_With_CPack>`_).
This file must define a macro called ``TRIBITS_PROJECT_DEFINE_PACKAGING()``
which is then invoked by TriBITS.  The file:

  `TribitsExampleProject`_/``cmake/CallbackDefineProjectPackaging.cmake``

provides a good example which is:

.. include:: ../../examples/TribitsExampleProject/cmake/CallbackDefineProjectPackaging.cmake
   :literal:

The CPack variables show above that should be defined at the project-level are
described in the `Official CPack Documentation`_.

Settings that are general for all distributions (like non-package repository
files to exclude from the tarball) should be set at the in the file
`<repoDir>/cmake/CallbackDefineRepositoryPackaging.cmake`_.  See `Creating
Source Distributions` for more details.

.. _<projectDir>/cmake/tribits/:

**<projectDir>/cmake/tribits/**: [Optional] This is the typical location of
the `TriBITS/tribits/`_ source tree for projects that choose to snapshot
TriBITS into their source tree.  In fact, TriBITS assumes this is the default
location for the TriBITS source tree if ``${PROJECT_NAME}_TRIBITS_DIR`` is not
otherwise specified.  Trilinos, for example, currently snapshots the TriBITS
source tree into this directory.  See `TriBITS directory snapshotting`_ for
more details.

.. _<projectDir>/cmake/ctest/CTestCustom.cmake.in:

**<projectDir>/cmake/ctest/CTestCustom.cmake.in**: [Optional] If this file
exists, it is processed using a ``CONFIGURE_FILE()`` command to write the file
``CTestCustom.cmake`` in the project base build directory
```${PROJECT_BINARY_DIR}/``.  This file is picked up automatically by
``ctest`` (see `CTest documentation`_).  This file is typically used to change
the maximum size of test output.  For example, the
`TribitsExampleProject`_/``cmake/ctest/CTestCustom.cmake.in`` looks like:

.. include:: ../../examples/TribitsExampleProject/cmake/ctest/CTestCustom.cmake.in
   :literal:

which sets the output size for each test submitted to CDash be unlimited
(which is not really recommended).  These variables used by Trilinos at one
time were::

  SET(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE 50000)
  SET(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE 5000000)

which sets the max output for passed and failed tests to 50000k and 5000000k,
respectively.

For documentation of the options one can change for CTest, see online `CTest
documentation`_.


TriBITS Project Core Variables
..............................

The following `local variables` are defined in the top-level Project
``CMakeLists.txt`` file scope and are therefore accessible by all files
processed by TriBITS:

  .. _PROJECT_NAME:

  ``PROJECT_NAME``

    The name of the TriBITS Project.  This exists to support, among other
    things, the ability for subordinate units (Repositories and Packages) to
    determine the Project in which is participating.  This is typically read
    from a ``SET()`` statement in the project's
    `<projectDir>/ProjectName.cmake`_ file.

  .. _PROJECT_SOURCE_DIR:

  ``PROJECT_SOURCE_DIR``

    The absolute path to the base Project source directory.  This is set
    automatically by TriBITS given the directory passed into ``cmake`` at
    configure time at the beginning of the `TRIBITS_PROJECT()`_ macro.

  .. _PROJECT_BINARY_DIR:

  ``PROJECT_BINARY_DIR``

    The absolute path to the base Project binary/build directory.  This is set
    automatically by TriBITS and is the directory where ``cmake`` is run from
    and is set at the beginning of the `TRIBITS_PROJECT()`_ macro.

  .. _${PROJECT_NAME}_SOURCE_DIR:

  ``${PROJECT_NAME}_SOURCE_DIR``

    Set to the same directory as ``${PROJECT_SOURCE_DIR}`` automatically by
    the built-in ``PROJECT()`` command called in the top-level
    `<projectDir>/CMakeLists.txt`_ file..

  .. _${PROJECT_NAME}_BINARY_DIR:

  ``${PROJECT_NAME}_BINARY_DIR``

    Set to the same directory as ``${PROJECT_BINARY_DIR}`` automatically by
    the built-in ``PROJECT()`` command called in the top-level
    `<projectDir>/CMakeLists.txt`_ file..

The following `cache variables` are defined for every TriBITS project:

  .. _${PROJECT_NAME}_TRIBITS_DIR:

  ``${PROJECT_NAME}_TRIBITS_DIR``

    CMake cache variable that gives the path to the TriBITS implementation
    directory.  When set to a relative path (set as type ``STRING``, see
    below), this is taken relative to ``${CMAKE_CURRENT_SOURCE_DIR}/`` (the
    project base source dir).  When an absolute path is given, it is used
    without modification.  If this variable is not set in the
    `<projectDir>/CMakeLists.txt`_ file, then it will be automatically set as
    a ``PATH`` cache variable by the include of ``TriBITS.cmake`` by the
    statement ::

      SET( ${PROJECT_NAME}_TRIBITS_DIR
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits" CACHE PATH "...")

    Therefore, projects that snapshot TriBITS into
    `<projectDir>/cmake/tribits/`_ don't need to explicitly set
    ``${PROJECT_NAME}_TRIBITS_DIR``.  In addition, one can also point to a
    different TriBITS implementation just by setting the absolute path::

      -D <Project>_TRIBITS_DIR=<some-abs-dir>

    or to a relative path using, for example::

      -D <Project>_TRIBITS_DIR:STRING=TriBITS/tribits

    Note that when the ``TriBITS`` git repo itself is cloned by a TriBITS
    project, then ``${PROJECT_NAME}_TRIBITS_DIR`` should be set to the
    directory ``TriBITS/tribits`` (see `TriBITS/tribits/`_) as shown above.

  .. _${PROJECT_NAME}_ENABLE_TESTS:

  ``${PROJECT_NAME}_ENABLE_TESTS``

    CMake cache variable that if set to ``ON``, then tests for all explicitly
    enabled packages will be turned on.  This has a default value of ``OFF``.
    This is used in logic to enable individual SE package tests (see
    `<Project>_ENABLE_TESTS only enables explicitly enabled SE package
    tests`_).

  ``${PACKAGE_NAME}_ENABLE_EXAMPLES``

    CMake cache variable that if set to ``ON``, then examples for all
    explicitly enabled packages will be turned on.  This has a default value
    of ``OFF``.

The following `internal project-scope local` (non-cache) CMake variables are
defined by TriBITS giving the project's TriBITS repositories.:

  .. _${PROJECT_NAME}_NATIVE_REPOSITORIES:

  ``${PROJECT_NAME}_NATIVE_REPOSITORIES``

     The list of Native Repositories for a given TriBITS project
     (i.e. Repositories that are always present when configuring the Project
     and are managed in the same VC repo typically).  This variable is set in
     the file `<projectDir>/cmake/NativeRepositoriesList.cmake`_ if it exists.
     If the file ``NativeRepositoriesList.cmake`` does not exist, then the
     project is assumed to also be a repository and the list of native
     repositories is just the local project directory
     ``${PROJECT_SOURCE_DIR}/.``.  In this case, the
     ``${PROJECT_SOURCE_DIR}/`` must contain at a minimum a
     ``PackagesList.cmake`` file, and a ``TPLsList.cmake`` file (see `TriBITS
     Repository`_).

  .. _${PROJECT_NAME}_EXTRA_REPOSITORIES:

  ``${PROJECT_NAME}_EXTRA_REPOSITORIES``

     The list of Extra Repositories that the project is being configured with.
     This list of repositories either comes from processing the project's
     `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file or comes from the
     CMake cache variable ``${PROJECT_NAME}_EXTRA_REPOSITORIES``.  See
     `Enabling extra repositories with add-on packages`_ for details.

  .. _${PROJECT_NAME}_ALL_REPOSITORIES:

  ``${PROJECT_NAME}_ALL_REPOSITORIES``

    Concatenation of all the repos listed in
    `${PROJECT_NAME}_NATIVE_REPOSITORIES`_ and
    `${PROJECT_NAME}_EXTRA_REPOSITORIES`_ in the order they are processed.


TriBITS Repository
++++++++++++++++++

A TriBITS Repository is the basic unit of ready-made composition between
different collections of software that use the TriBITS CMake build and system.

In short, a TriBITS Repository:

* Is a named collection of related TriBITS Packages and TPLs (defined in
  `<repoDir>/PackagesList.cmake`_ and `<repoDir>/TPLsList.cmake`_
  respectively)
* Defines the base source and binary directories for the Repository
  ``${REPOSITORY_NAME}_SOURCE_DIR`` and ``${REPOSITORY_NAME}_BINARY_DIR``.
* Defines a common set of initializations and other hooks for all the
  packages in the repository.
* Typically maps to a VC (i.e. git) repository and therefore represents a unit
  of integration, versioning and reuse.  (But core TriBITS has no dependency
  on any VC tool.)

For more details on the definition of a TriBITS Repository, see:

* `TriBITS Repository Core Files`_
* `TriBITS Repository Core Variables`_

TriBITS Repository Core Files
.............................

The core files making up a TriBITS Repository (where ``<repoDir> =
${${REPOSITORY_NAME}_SOURCE_DIR}``) are::

  <repoDir>/
    PackagesList.cmake
    TPLsList.cmake
    Copyright.txt  # [Optional] Only needed if creating version header file
    Version.cmake  # [Optional] Info inserted into ${REPO_NAME}_version.h
    cmake/
       RepositoryDependenciesSetup.cmake # [Optional] CDash email addresses?
       CallbackSetupExtraOptions.cmake # [Optional] Called after main options
       CallbackDefineRepositoryPackaging.cmake # [Optional] CPack packaging

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
`TRIBITS_REPOSITORY_DEFINE_PACKAGES()`_ to define the list of packages along
with their directories and other properties.  For example, the file
`TribitsExampleProject`_/``PackagesList.cmake`` looks like:

.. include:: ../../examples/TribitsExampleProject/PackagesList.cmake
   :literal:

Other commands that are appropriate to use in this file include
`TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS()`_ and
`TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES()`_.  Also, if the binary directory
for any package ``<packageName>`` needs to be changed from the default, then
the variable ``<packageName>_SPECIFIED_BINARY_DIR`` can be set.  (see `TriBITS
Package == TriBITS Repository == TriBITS Project`_).

It is perfectly legal for a TriBITS repository to define no packages at all
with::

  TRIBITS_REPOSITORY_DEFINE_PACKAGES()

and this would be the case for a TriBITS meta-project that has no native
packages, only extra repositories.

.. ToDo: Point to example meta-project.

.. _<repoDir>/TPLsList.cmake:

**<repoDir>/TPLsList.cmake**: [Required] Provides the list of TPLs that are
referenced as TPL dependencies in the repository's SE package's
`<packageDir>/cmake/Dependencies.cmake`_ files (see `TriBITS TPL`_).  This
file typically just calls the macro `TRIBITS_REPOSITORY_DEFINE_TPLS()`_ to
define the TPLs along with their find modules and other properties.  An
example is `ReducedMockTrilinos`_/``TPLsList.cmake`` which shows:

.. include:: ../../examples/ReducedMockTrilinos/TPLsList.cmake
   :literal:

See `TriBITS TPL`_ for details on what gets defined for each TriBITS TPL once
this file is processed.

It is perfectly fine to specify no TPLs at all for a repository with::

  TRIBITS_REPOSITORY_DEFINE_TPLS()

but the macro ``TRIBITS_REPOSITORY_DEFINE_TPLS()`` has to be called, even if
there are no TPLs.  See `TRIBITS_REPOSITORY_DEFINE_TPLS()`_ for further
details and constraints.

.. _<repoDir>/Copyright.txt:

**<repoDir>/Copyright.txt**: [Optional] Gives the default copyright and
license declaration for all of the software in the TriBITS repository
directory ``<repoDir>/``.  This file is read into a string and then used to
configure the repository's version header file (see `Project and Repository
Versioning and Release Mode`_).  Even if a repository version header file is
not produced, it is a good idea for every TriBITS repository to define this
file, just for legal purposes.  For a good open-source license, one should
consider copying the ``TriBITS/Copyright.txt`` file which is a simple 3-clause
BSD-like license like:

.. include:: ../../Copyright.txt
   :literal:

.. _<repoDir>/Version.cmake:

**<repoDir>/Version.cmake**: [Optional] Contains version information for the
repository (and the project also if this is also the base project).  For
example, `TribitsExampleProject`_/``Version.cmake``, this looks like:

.. include:: ../../examples/TribitsExampleProject/Version.cmake
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
gets included in the order listed in `${PROJECT_NAME}_ALL_REPOSITORIES`_.  Any
local variables set in this file have project-wide scope.  The primary purpose
for this file is to set variables that will impact the processing of project's
package's ``Dependencies.cmake`` files and take care of other enable/disable
issues that are not otherwise cleanly handled by the TriBITS system
automatically.

The typical usage of this file is to set the default CDash email address for
all of the defined packages (see `CDash regression email addresses`_).  For
example, to set the default email address for all of the packages in this
repository, one would set in this file::

   SET_DEFAULT(${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESS
       repox-regressions@somemailserver.org)

Note that the prefix ``${REPOSITORY_NAME}_`` is used instead of hard-coding
the repo name to allow greater flexibility in how meta-projects refer to this
TriBITS repo.

.. _<repoDir>/cmake/CallbackSetupExtraOptions.cmake:

**<repoDir>/cmake/CallbackSetupExtraOptions.cmake**: [Optional] If defined,
this file is processed (included) for each repo in order right after the basic
TriBITS options are defined in the macro
``TRIBITS_DEFINE_GLOBAL_OPTIONS_AND_DEFINE_EXTRA_REPOS()``.  This file must
define the macro ``TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS()`` which is then
called by the TriBITS system.  This file is only processed when doing a basic
configuration of the project and **not** when it is just building up the
dependency data-structures (i.e. it is **not** processed in the `Reduced
Package Dependency Processing`_).  Any local variables set have project-wide
scope.

A few additional variables are defined by the time this file is processed and
can be used in the logic in these files.  Some of the variables that should
already be defined (in addition to all of the basic user TriBITS cache
variables set in ``TRIBITS_DEFINE_GLOBAL_OPTIONS_AND_DEFINE_EXTRA_REPOS()``)
include ``CMAKE_HOST_SYSTEM_NAME``, ``${PROJECT_NAME}_HOSTNAME``, and
``PYTHON_EXECUTABLE`` (see `Python Support`_).  The types of commands and
logic to put in this file include:

* Setting additional user cache variable options that are used by multiple
  packages in the TriBITS Repository.  For example, Trilinos defines a
  ``Trilinos_DATA_DIR`` user cache variable that several Trilinos packages use
  to get extra test data to define extra tests.
* Disabling packages in the TriBITS Repository when conditions will not allow
  them to be enabled.  For example, Trilinos disables the package
  ``ForTrilinos`` when ``Fortran`` is disabled and disables the package
  ``PyTrilinos`` when Python support is disabled.

An example of this file is:

  `TribitsExampleProject`_/``/cmake/CallbackSetupExtraOptions.cmake``

which currently looks like:

.. include:: ../../examples/TribitsExampleProject/cmake/CallbackSetupExtraOptions.cmake
   :literal:

.. _<repoDir>/cmake/CallbackDefineRepositoryPackaging.cmake:

**<repoDir>/cmake/CallbackDefineRepositoryPackaging.cmake**: [Optional] If
this file exists, then it defines extra CPack-related options that are
specific to this TriBITS Repository.  This file must define the macro
``TRIBITS_REPOSITORY_DEFINE_PACKAGING()`` which is called by TriBITS.  This
file is processed as the top project-level scope so any local variables set
have project-wide effect.  This file is processed after the project's
`<projectDir>/cmake/CallbackDefineProjectPackaging.cmake`_ file so any project
CPACK variables are defined for the repository-level options and commands are
created.  This file typically just sets extra excludes to remove files from
the tarball.  The file:

  `TribitsExampleProject`_/``cmake/CallbackDefineRepositoryPackaging.cmake``

provides a good example which is:

.. include:: ../../examples/TribitsExampleProject/cmake/CallbackDefineRepositoryPackaging.cmake
   :literal:

As shown in the above example, it is important to prefix the excluded files
and directories with the repository base directory
``${${REPOSITORY_NAME}_SOURCE_DIR}/`` since these are interpreted by CPack as
regular-expressions.


TriBITS Repository Core Variables
.................................

The following temporary local variables are defined automatically by TriBITS
before processing a given TriBITS repository's files
(e.g. ``PackagesList.cmake``, ``TPLsList.cmake``, etc.):

  .. _REPOSITORY_NAME:

  ``REPOSITORY_NAME``

    The name of the current TriBITS repository.  This name will be the
    repository name listed in
    `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file or if this
    repository directory is the project base directory, ``REPOSITORY_NAME``
    will be set to ``${PROJECT_NAME}``.

  ``REPOSITORY_DIR``

    Path of the current Repository *relative* to the Project's base source
    directory `${PROJECT_NAME}_SOURCE_DIR`_..  This is typically just the
    repository name but can be an arbitrary directory if specified through the
    `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file.

The following project-scope (non-cache) local variables are set once the list
of TriBITS repositories is processed and before any of the repository's files
are processed:

  ``${REPOSITORY_NAME}_SOURCE_DIR``

    The absolute path to the base of a given TriBITS Repository's source
    directory.  CMake code, for example in a packages' ``CMakeLists.txt``
    file, typically refers to this by the raw name like ``RepoX_SOURCE_DIR``.
    This makes such CMake code independent of where the various TriBITS repos
    are in relation to each other or the TriBITS Project (but does hard-code
    the repository name which is not ideal).

  ``${REPOSITORY_NAME}_BINARY_DIR``

    The absolute path to the base of a given TriBITS Repository's binary
    directory.  CMake code, for example in packages, refer to this by the raw
    name like ``RepoX_SOURCE_DIR``.  This makes such CMake code independent of
    where the various TriBITS repos are in relation to each other or the
    Project.

The following project-level local variables can be defined by the project or
the user to help define the what packages from the repository
``${REPOSITORY_NAME}`` contribute to the primary meta-project packages (PMPP):

  .. _${REPOSITORY_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES:

  ``${REPOSITORY_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES``

    If set to ``TRUE``, then the package's in the TriBITS repository are not
    considered to be part of the primary meta-project packages.  This affects
    what packages get enabled by default when enabling all packages with
    setting ``${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON`` and what tests and
    examples get enabled by default when setting
    ``${PROJECT_NAME}_ENABLE_TESTS=ON``.  See `TriBITS Dependency Handling
    Behaviors`_ for more details.

  .. _${REPOSITORY_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES_EXCEPT:

  ``${REPOSITORY_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES_EXCEPT``

     When the above variable is set to ``TRUE``, this variable is read by
     TriBITS to find the list of TriBITS packages selected packages in the
     repository ``${REPOSITORY_NAME}`` which are considered to be part of the
     set of the project's primary meta-project package when the above variable
     is set to ``ON``.  NOTE: It is not necessary to list all of the
     subpackages in a given parent package.  Only the parent package need be
     listed and it will be equivalent to listing all of the subpackages.  See
     `TriBITS Dependency Handling Behaviors`_ for more details.

The above primary meta-project variables should be set in the meta-project's
`<projectDir>/ProjectName.cmake`_ file so that they will be set in all
situations.


TriBITS Package
+++++++++++++++

A TriBITS Package:

* Is the fundamental TriBITS structural unit of software partitioning and
  aggregation.
* Must have a unique package name (``PACKAGE_NAME``) that is globally unique
  (see `Globally unique TriBITS package names`_).
* Typically defines a set of libraries and/or header files and/or executables
  and/or tests with CMake build targets for building these for which TriBITS
  exports the list of include directories, libraries, and targets that are
  created (along with CMake dependencies).
* Is declared in its parent repository's `<repoDir>/PackagesList.cmake`_ file.
* Declares dependencies on `upstream`_ TPLs and/or other SE packages by just
  naming the dependencies in the file
  `<packageDir>/cmake/Dependencies.cmake`_.
* Can optionally have subpackages listed in the argument
  `SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS`_ to the macro call
  `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_.
* Is the unit of testing as driven by `TRIBITS_CTEST_DRIVER()`_ and displayed
  on CDash.

.. _Globally unique TriBITS package names:

**WARNING:** As noted above, one must be very careful to pick **globally
unique TriBITS package names**.  This name must be unique not only within its
defined TriBITS repository but also across all SE packages in all TriBITS
repositories that ever might be cobbled together into a single TriBITS (meta)
project!  Choosing a good package name is the single most important decision
when it comes to defining a TriBITS package.  One must be careful **not** to
pick names like "Debug" or "StandardUtils" that have a high chance of clashing
with poorly named TriBITS packages from other TriBITS repositories.

For more details on the definition of a TriBITS Package (or subpackage), see:

* `TriBITS Package Core Files`_
* `TriBITS Package Core Variables`_

TriBITS Package Core Files
..........................

The core files that make up a *TriBITS Package* (where ``<packageDir> =
${${PACKAGE_NAME}_SOURCE_DIR}``) are::

  <packageDir>/
    CMakeLists.txt  # Only processed if the package is enabled
    cmake/
      Dependencies.cmake  # Always processed if its repo is processed
      <packageName>_config.h.in  # [Optional], name is not fixed

There are a few simple rules for the location and the contents of the
``<packageDir>/`` directory:

* The directory ``<packageDir>/`` must not be a subdirectory of the package
  directory of any other SE package (e.g. not ``pkga/pkgb``).
* All of the source files, test files, etc. for the package should be included
  under ``<packageDir>/``.

The above rules are not needed for basic building and testing but are needed
for extended features like automatically detecting when a package has changed
by looking at what files have changed (see `Pre-push Testing using
checkin-test.py`_) and for creating source tarballs correctly (see `Creating
Source Distributions`_).  Therefore, it would be wise to abide by the above
rules when defining packages.

The following TriBITS Package files are documented in more detail below:

* `<packageDir>/cmake/Dependencies.cmake`_
* `<packageDir>/cmake/<packageName>_config.h.in`_
* `<packageDir>/CMakeLists.txt`_

.. _<packageDir>/cmake/Dependencies.cmake:

**<packageDir>/cmake/Dependencies.cmake**: [Required] Defines the dependencies
for a given TriBITS (SE) package using the macro
`TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_.  This file is processed at the
top-level project scope (using an ``INCLUDE()``) so any local variables set
will be seen by the entire project.  This file is always processed, including
when just building the project's dependency data-structure (see `Reduced
Package Dependency Processing`_).

An example of a ``Dependencies.cmake`` file for a package with optional and
required dependencies is for the mock ``Panzer`` package in `MockTrilinos`_:

.. include:: ../../examples/MockTrilinos/packages/panzer/cmake/Dependencies.cmake
   :literal:

.. _with_subpackages/cmake/Dependencies.cmake:

An example of a package with subpackages is ``WithSubpackages`` which has the
dependencies file:

  `TribitsExampleProject`_/``packages/with_subpackages/cmake/Dependencies.cmake``

which is:

.. include:: ../../examples/TribitsExampleProject/packages/with_subpackages/cmake/Dependencies.cmake
   :literal:

``WithSubpackages`` defines three subpackages which creates three new SE
packages with names ``WithSubpackagesA``, ``WithSubpackagesB``, and
``WithSubpackagesC``.

if a TriBITS Package or Subpackage has no dependencies, it still has to call
``TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`` but it is called with no arguments
such as with:

  `TribitsHelloWorld`_/``hello_world/cmake/Dependencies.cmake:``

which contains:

.. include:: ../../examples/TribitsHelloWorld/hello_world/cmake/Dependencies.cmake
   :literal:

Other TriBITS macros/functions that can be called in this file include
`TRIBITS_TPL_TENTATIVELY_ENABLE()`_ and
`TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES()`_.

.. _<packageDir>/cmake/<packageName>_config.h.in:

**<packageDir>/cmake/<packageName>_config.h.in**: [Optional] The package's
configured header file.  This file will contain placeholders for variables
that will be substitute at configure time with `TRIBITS_CONFIGURE_FILE()`_.
This includes usage of ``#cmakedefine <varName>`` and other standard CMake
file configuration features.

**NOTE:** The file name ``<packageName>_config.h.in`` is not at all fixed and
the package can call this file anything it wants.  Also, a package can
configure multiple header files in different directories for different
purposes using `TRIBITS_CONFIGURE_FILE()`_ or even calls to the raw CMake
function `CONFIGURE_FILE()`_.

.. _<packageDir>/CMakeLists.txt:

**<packageDir>/CMakeLists.txt**: [Required] The package's top-level
``CMakeLists.txt`` file that defines the libraries, include directories, and
contains the tests for the package.

The basic structure of this file for a **package without subpackages** is
shown in:

  `TribitsExampleProject`_/``packages/simple_cxx/CMakeLists.txt``

which is:

.. include:: ../../examples/TribitsExampleProject/packages/simple_cxx/CMakeLists.txt
   :literal:

The first command at the top of the file is a call to `TRIBITS_PACKAGE()`_
which takes the package name (``SimpleCxx`` in this case) in addition to a few
other options.  While TriBITS obviously already knows the package name (since
it read it from the `<repoDir>/PackagesList.cmake`_ file), the purpose for
repeating it in this call is as documentation for the developer's sake (and
this name is checked against the expected package name).  Then a set of
configure-time tests are typically performed (if the package needs any of
these).  In this example, the existence of the C++ ``__int64`` data-type is
checked using the module ``CheckFor__int64.cmake`` (which is in the ``cmake/``
directory of this package.  (CMake has great support `Configure-time System
Tests`_.)  This is followed by package-specific options.  In this case, the
standard TriBITS options for debug checking and deprecated warnings are added
using the standard macros `TRIBITS_ADD_DEBUG_OPTION()`_ and
`TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()`_.  After all of this up-front
stuff is complete (which will be present in any moderately complex
CMake-configured project) the source and the test sub-directories are added
that actually define the library and the tests.  In this case, the standard
`TRIBITS_ADD_TEST_DIRECTORIES()`_ macro is used which only conditionally adds
the tests for the package.

The final command in the package's base ``CMakeLists.txt`` file must always be
`TRIBITS_PACKAGE_POSTPROCESS()`_.  This is needed in order to perform some
necessary post-processing by TriBITS.

It is also possible for the package's top-level ``CMakeLists.txt`` to be the
only ``CMakeLists.txt`` file for a package.  Such an example can be seen in the
example project `TribitsHelloWorld`_ in the ``HelloWorld`` package.

When a TriBITS package is broken up into subpackages (see `TriBITS
Subpackage`_), its ``CMakeLists.txt`` file looks a little different from a
package with no subpackages as shown above.  The basic structure of this file
for a **package with subpackages** is shown in:

  `TribitsExampleProject`_/``packages/with_subpackages/CMakeLists.txt``

which contains:

.. include:: ../../examples/TribitsExampleProject/packages/with_subpackages/CMakeLists.txt
   :literal:

What is different about ``CMakeLists.txt`` files for packages without
subpackages is that the `TRIBITS_PACKAGE()`_ command is broken up into two
parts `TRIBITS_PACKAGE_DECL()`_ and `TRIBITS_PACKAGE_DEF()`_.  In between
these two commands, the parent package can define the common package options
and then calls the command `TRIBITS_PROCESS_SUBPACKAGES()`_ which fully
processes the packages.  If the parent package has libraries and/or
tests/example of its own, it can define those after calling
`TRIBITS_PACKAGE_DEF()`_, just like with a regular package.  However, it is
rare for a package broken up into subpackages to have its own libraries and/or
tests and examples.  As always, the final command called inside of a package's
base ``CMakeLists.txt`` file is `TRIBITS_PACKAGE_POSTPROCESS()`_.

NOTE: The package's base ``CMakeLists.txt`` file only gets processed if the
package is actually enabled
(i.e. ``${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=ON``).  This is an important
design feature of TriBITS in that the contents of non-enabled packages can't
damage the configure, build, and test of the enabled packages based on errors
in non-enabled packages.  This is critical to allow experimental `EX`_
test-group packages and lower-maturity packages to exist in the same source
repositories safely with higher-maturity and more important packages.


TriBITS Package Core Variables
..............................

A packages' core variables are broken down into the following categories:

* `TriBITS Package Local Variables`_
* `TriBITS Package Top-Level Local Variables`_
* `TriBITS Package Cache Variables`_
* `TriBITS Package Optional Dependency Macro Variables`_

.. _TriBITS Package Local Variables:

The following locally scoped *TriBITS Package Local Variables* are defined
when the files for a given TriBITS Package (or any SE package for that matter)
are being processed:

  ``PACKAGE_NAME``

    The name of the current TriBITS SE package.  This is set automatically by
    TriBITS before the packages' ``CMakeLists.txt`` file is processed.
    **WARNING:** This name must be globally unique across the entire project
    (see `Globally unique TriBITS package names`_).

  ``PACKAGE_SOURCE_DIR``

    The absolute path to the package's base source directory.    This is
    set automatically by TriBITS in the macro `TRIBITS_PACKAGE()`_.

  ``PACKAGE_BINARY_DIR``

    The absolute path to the package's base binary/build directory.  This is
    set automatically by TriBITS in the macro `TRIBITS_PACKAGE()`_.

  ``PACKAGE_NAME_UC``

    This is set to the upper-case version of ``${PACKAGE_NAME}``.  This is set
    automatically by TriBITS in the macro `TRIBITS_PACKAGE()`_.

.. _TriBITS Package Top-Level Local Variables:

Once all of the TriBITS SE package's ``Dependencies.cmake`` files have been
processed, the following *TriBITS Package Top-Level Local Variables* are
defined:

  ``${PACKAGE_NAME}_SOURCE_DIR``

    The absolute path to the package's base source directory.  CMake code, for
    example in other packages, refer to this by the raw name like
    ``PackageX_SOURCE_DIR``.  This makes such CMake code independent of where
    the package is in relation to other packages.  NOTE: This variable is
    defined for all declared packages that exist, independent of whether they
    are enabled or not.  This variable is set as soon as it is known if the
    given package exists or not.

  ``${PACKAGE_NAME}_BINARY_DIR``

    The absolute path to the package's base binary directory.  CMake code, for
    example in other packages, refer to this by the raw name like
    ``PackageX_BINARY_DIR``.  This makes such CMake code independent of where
    the package is in relation to other packages.  NOTE: This variable is
    **only** defined if the package is actually enabled!

  ``${PACKAGE_NAME}_PARENT_REPOSITORY``

    The name of the package's parent repository.  This can be used by a
    package to access information about its parent repository.  For example,
    the variable ``${${PACKAGE_NAME}_PARENT_REPOSITORY}_SOURCE_DIR`` can be
    dereferenced and read of needed (but it is not recommended that packages
    be aware of their parent repository in general)..

  .. _${PACKAGE_NAME}_TESTGROUP:

  ``${PACKAGE_NAME}_TESTGROUP``

    Defines the `SE Package Test Group`_ for the package.  This determines in
    what contexts the package is enabled or not for testing-related purposes
    (see `Nested Layers of TriBITS Project Testing`_)

.. _TriBITS Package Cache Variables:

In addition, the following user-settable *TriBITS Package Cache Variables* are
defined before a (SE) Package's ``CMakeLists.txt`` file is processed:

  .. _${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}:

  ``${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}``

    Set to ``ON`` if the package is enabled and is to be processed.

  .. _${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE_NAME}:

  ``${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE_NAME}``

    Set to ``ON`` if support for the optional `upstream`_ dependent package
    ``${OPTIONAL_DEP_PACKAGE_NAME}`` is enabled in package
    ``${PACKAGE_NAME}``.  Here ``${OPTIONAL_DEP_PACKAGE_NAME}`` corresponds to
    each optional upstream SE package listed in the ``LIB_OPTIONAL_PACKAGES``
    and ``TEST_OPTIONAL_PACKAGES`` arguments to the
    `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_ macro.

    **NOTE:** It is important that the CMake code in the package's
    ``CMakeLists.txt`` files key off of this variable and **not** the
    project-level variable
    ``${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE_NAME}`` because the
    package-level variable
    ``${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE_NAME}`` can be explicitly
    turned off by the user even through the packages ``${PACKAGE_NAME}`` and
    ``${OPTIONAL_DEP_PACKAGE_NAME}`` are both enabled at the project level!
    See `Support for optional SE package/TPL can be explicitly disabled`_.

  .. _${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL_NAME}:

  ``${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL_NAME}``

    Set to ``ON`` if support for the optional `upstream`_ dependent TPL
    ``${OPTIONAL_DEP_TPL_NAME}`` is enabled in package ``${PACKAGE_NAME}``.
    Here ``${OPTIONAL_DEP_TPL_NAME}`` corresponds to each optional upstream
    TPL listed in the ``LIB_OPTIONAL_TPLS`` and ``TEST_OPTIONAL_TPLS``
    arguments to the `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_ macro.
    **NOTE:** It is important that the CMake code in the package
    ``${PACKAGE_NAME}`` key off of this variable and **not** the global
    ``TPL_ENABLE_${OPTIONAL_DEP_TPL_NAME}`` variable because
    ``${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL_NAME}`` can be explicitly
    turned off by the user even through the package ``${PACKAGE_NAME}`` and
    the TPL ``${OPTIONAL_DEP_TPL_NAME}`` are both enabled (see `Support for
    optional SE package/TPL can be explicitly disabled`_).

  .. _${PACKAGE_NAME}_ENABLE_TESTS:

  ``${PACKAGE_NAME}_ENABLE_TESTS``

    Set to ``ON`` if the package's tests are to be enabled.  This will enable
    a package's tests and all of its subpackage's tests.

  .. _${PACKAGE_NAME}_ENABLE_EXAMPLES:

  ``${PACKAGE_NAME}_ENABLE_EXAMPLES``

    Set to ``ON`` if the package's examples are to be enabled.  This will
    enable a package's examples and all of its subpackage's examples.

The above global cache variables can be explicitly set by the user or may be
set automatically as part of the `Package Dependencies and Enable/Disable
Logic`_.

.. _TriBITS Package Optional Dependency Macro Variables:

The following local *TriBITS Package Optional Dependency Macro Variables* are
defined in the top-level project scope before a (SE) Package's
``CMakeLists.txt`` file is processed:

  .. _HAVE_<PACKAGE_NAME_UC>_<OPTIONAL_DEP_PACKAGE_NAME_UC>:

  ``HAVE_<PACKAGE_NAME_UC>_<OPTIONAL_DEP_PACKAGE_NAME_UC>``

    Set to ``ON`` if support for optional upstream package
    ``${OPTIONAL_DEP_PACKAGE}`` is enabled in downstream package
    ``${PACKAGE_NAME}``
    (i.e. `${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE_NAME}`_ = ``ON``) and
    is set to ``FALSE`` otherwise.  Here, ``<PACKAGE_NAME_UC>`` and
    ``<OPTIONAL_DEP_PACKAGE_NAME_UC>`` are the upper-case names for the packages
    ``${PACKAGE_NAME}`` and ``${OPTIONAL_DEP_PACKAGE_NAME}``, respectively.
    For example, if optional support for upstream package ``Triutils`` is
    enabled in downstream package ``EpetraExt`` in `ReducedMockTrilinos`_,
    then ``EpetraExt_ENABLE_TriUtils=ON`` and ``HAVE_EPETRAEXT_TRIUTILS=ON``.
    This variable is meant to be used in::

      #cmakedefine HAVE_<PACKAGE_NAME_UC>_<OPTIONAL_DEP_PACKAGE_NAME_UC>

    in configured header files
    (e.g. `<packageDir>/cmake/<packageName>_config.h.in`_).  For example, for
    the ``EpetraExt`` and ``Triutils`` example, this would be::

      #cmakedefine HAVE_EPETRAEXT_TRIUTILS

  .. _HAVE_<PACKAGE_NAME_UC>_<OPTIONAL_DEP_TPL_NAME_UC>:

  ``HAVE_<PACKAGE_NAME_UC>_<OPTIONAL_DEP_TPL_NAME_UC>``

    Set to ``ON`` if support for the optional TPL ``${OPTIONAL_DEP_TPL_NAME}``
    is enabled in downstream package ``${PACKAGE_NAME}``
    (i.e. `${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL_NAME}`_ = ``ON``) and is
    set to ``FALSE`` otherwise.  Here, ``<PACKAGE_NAME_UC>`` and
    ``<OPTIONAL_DEP_TPL_NAME_UC>`` are the upper-case names for the packages
    ``${PACKAGE_NAME}`` and ``${OPTIONAL_DEP_TPL_NAME}``, respectively.  For
    example, if optional support for the TPL ``Boost`` is enabled in the
    package ``Teuchos`` in `ReducedMockTrilinos`_, then
    ``Teuchos_ENABLE_Boost=ON`` and ``HAVE_TUECHOS_BOOST=ON``.  This variable
    is meant to be used in::

      #cmakedefine HAVE_<PACKAGE_NAME_UC>_<OPTIONAL_DEP_TPL_NAME_UC>

    in configured header files
    (e.g. `<packageDir>/cmake/<packageName>_config.h.in`_).  For example, for
    the ``Teuchos`` and ``Boost`` example, this would be::

      #cmakedefine HAVE_TEUCHOS_BOOST

Currently, a Package can refer to its containing Repository and refer to its
source and binary directories.  This is so that it can refer to
repository-level resources (e.g. like the ``Trilinos_version.h`` file for
Trilinos packages).  However, this may be undesirable because it will make it
hard to pull a package out of one TriBITS repository and place it in another
repository for a different use.  However, a package can indirectly refer to
its own repository without loss of generality by reading the variable
``${PACKAGE_NAME}_PARENT_REPOSITORY``.  The problem is referring to other
TriBITS repositories explicitly.


TriBITS Subpackage
++++++++++++++++++

A TriBITS Subpackage:

* Is a compartmentalization of a parent `TriBITS Package`_ according to
  `Software Engineering Packaging Principles`_.
* Typically defines a set of libraries and/or header files and/or executables
  and/or tests with CMake build targets for building these for which TriBITS
  exports the list of include directories, libraries, and targets that are
  created (along with CMake dependencies).
* Is declared in its parent package's `<packageDir>/cmake/Dependencies.cmake`_
  file in a call to `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_ using the
  argument `SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS`_.
* Defines dependencies on `upstream`_ TPLs and/or other SE packages by just
  naming the dependencies in the file
  `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_ using the macro
  `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_.
* Can **NOT** have its own subpackages defined (only top-level packages can
  have subpackages).
* Is enabled or disabled along with all other subpackages in the parent
  package automatically if it's parent package is enabled or disabled with
  ``${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME}`` set to ``ON`` or ``OFF``
  respectively (see `Enable/disable of parent package is enable/disable for
  subpackages`_ and `PARENT_PACKAGE_NAME`_).
* Has tests turned on automatically if
  ``${PARENT_PACKAGE_NAME}_ENABLE_TESTS=ON``.

The contents of a TriBITS Subpackage are almost identical to those of a
TriBITS Package.  The differences are described below and in `How is a TriBITS 
Subpackage different from a TriBITS Package?`_.

For more details on the definition of a TriBITS Package (or subpackage), see:

* `TriBITS Subpackage Core Files`_
* `TriBITS Subpackage Core Variables`_


TriBITS Subpackage Core Files
..............................

The set of core files for a subpackage are identical to the `TriBITS Package
Core Files`_ and are::

  <packageDir>/<spkgDir>/
    CMakeLists.txt  # Only processed if this subpackage is enabled
    cmake/
      Dependencies.cmake  # Always processed if the parent package
                            # is listed in the enclosing Repository

(where ``<packageDir> = ${${PARENT_PACKAGE_NAME}_SOURCE_DIR}`` and
``<spkgDir>`` is the subpackage directory listed in the
`SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS`_ to
`TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_).

There are a few simple rules for the location and the contents of the
``<spkgDir>/`` directory:

* The relative directory ``<spkgDir>/`` should be a strict subdirectory of
  ``<packageDir>/`` (e.g. not ``../../somewhereelse``).
* The directory ``<spkgDir>/`` must not be a subdirectory of the package
  directory of any other subpackage (e.g. not ``spkga/spkgb``).
* All of the source files, test files, etc. for the subpackage should be
  included under ``<spkgDir>/``.

The above rules are not needed for basic building and testing but are needed
for extended features like automatically detecting when a package has changed
by looking at what files have changed (see `Pre-push Testing using
checkin-test.py`_) and for creating source tarballs correctly (see `Creating
Source Distributions`_).  Therefore, it would be wise to abide by the above
rules when defining subpackages.

These TriBITS Subpackage files are documented in more detail below:

* `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_
* `<packageDir>/<spkgDir>/CMakeLists.txt`_

.. _<packageDir>/<spkgDir>/cmake/Dependencies.cmake:

**<packageDir>/<spkgDir>/cmake/Dependencies.cmake**: [Required] The contents
of this file for subpackages are identical as for top-level packages.  It just
contains a call to the macro `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_ to
define this SE package's `upstream`_ TPL and SE package dependencies.  A
simple example is for the example subpackage ``WithSubpackagesB`` (declared in
`with_subpackages/cmake/Dependencies.cmake`_) with the file:

  `TribitsExampleProject`_/``packages/with_subpackages/b/cmake/Dependencies.cmake``

which is:

.. include:: ../../examples/TribitsExampleProject/packages/with_subpackages/b/cmake/Dependencies.cmake
   :literal:

What this shows is that subpackages must list their dependencies on each other
(if such dependencies exist) using the full SE package name
``${PARENT_PACKAGE_NAME}${SUBPACKAGE_NAME}`` or in this case::

   'WithSubpackagesA' = 'WithSubpackages' + 'A'

Note that the parent SE package depends on its subpackages, not the other way
around.  For example, the ``WithSubpackages`` parent SE package
automatically depends its SE subpackages
``WithSubpackagesA``, ``WithSubpackagesC``,
and ``WithSubpackagesC``.  As such all (direct) dependencies
for a subpackage must be listed in its own ``Dependencies.cmake`` file.  For
example, the ``WithSubpackages`` subpackage ``A`` depends on
the ``SimpleCxx`` package and is declared as such as shown in:

  `TribitsExampleProject`_/``packages/with_subpackages/a/cmake/Dependencies.cmake``

which is:

.. include:: ../../examples/TribitsExampleProject/packages/with_subpackages/a/cmake/Dependencies.cmake
   :literal:

What this means is that any TPL or library dependencies listed in the parent
package's `<packageDir>/cmake/Dependencies.cmake`_ file are **NOT**
dependencies of its subpackages.  For example, if
`with_subpackages/cmake/Dependencies.cmake`_ where changed to be::

  TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
    LIB_REQUIRED_TPLS Boost
    SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
      A   A   PT  REQUIRED
      ...
    )

then the ``Boost`` TPL would **NOT** be a dependency of the SE package
``WithSubpackagesA`` but instead would be listed as a
dependency of the parent SE package ``WithSubpackages``.  (And in this
case, this TPL dependency is pretty worthless since the SE package
``WithSubpackages`` does not even define any libraries or tests of its
own.)

.. _<packageDir>/<spkgDir>/CMakeLists.txt:

**<packageDir>/<spkgDir>/CMakeLists.txt**: [Required] The subpackage's
top-level ``CMakeLists.txt`` file that defines the libraries, include
directories, and contains the tests for the subpackage.  The contents of a
subpackage's top-level ``CMakeLists.txt`` file are almost identical to a
top-level package's `<packageDir>/CMakeLists.txt`_ file.  The primary
difference is that the commands `TRIBITS_PACKAGE()`_ and
`TRIBITS_PACKAGE_POSTPROCESS()`_ and replaced with `TRIBITS_SUBPACKAGE()`_ and
`TRIBITS_SUBPACKAGE_POSTPROCESS()`_ as shown in the file:


  `TribitsExampleProject`_/``packages/with_subpackages/a/CMakeLists.txt``

which contains:

.. include:: ../../examples/TribitsExampleProject/packages/with_subpackages/a/CMakeLists.txt
   :literal:

Unlike `TRIBITS_PACKAGE()`_, `TRIBITS_SUBPACKAGE()`_ does not take any extra
arguments.  Those extra settings are assumed to be defined by the top-level
parent package.  Like top-level packages, subpackages are free to define
user-settable options and configure-time tests but typically don't.  The idea
is that subpackages should be lighter weight than top-level packages.  Other
than using `TRIBITS_SUBPACKAGE()`_ and `TRIBITS_SUBPACKAGE_POSTPROCESS()`_, a
subpackage can be laid out just like any other package and can call on any
other commands to add libraries, add executables, add test, etc.


TriBITS Subpackage Core Variables
.................................

The core variables associated with a subpackage are identical to the `TriBITS
Package Core Variables`_.  The only difference is that a subpackage may need
to refer to its parent package where a top-level package does not have a
parent package.  The extra variables that are defined when processing a
subpackage's files are:

  .. _PARENT_PACKAGE_NAME:

  ``PARENT_PACKAGE_NAME``

    The name of the parent package.

  ``PARENT_PACKAGE_SOURCE_DIR``

    The absolute path to the parent package's base source directory.

  ``PARENT_PACKAGE_BINARY_DIR``

    The absolute path to the parent package's base binary directory.


How is a TriBITS Subpackage different from a TriBITS Package?
.............................................................

A common question this is natural to ask is how a TriBITS Subpackage is
different from a TriBITS Package?  They contain the same basic files (i.e. a
``cmake/Dependencies.cmake`` file, a top-level ``CMakeLists.txt`` file, source
files, test files, etc.).  They both are included in the list of `TriBITS SE
Packages`_ and therefore can both be enabled/disabled by the user or in
automatic dependency logic (see `Package Dependencies and Enable/Disable
Logic`_).  The primary difference is that a subpackage is meant to involve
less overhead and is to be used to partition the parent package's software
into chunks according to `Software Engineering Packaging Principles`_.  Also,
the dependency logic treats a parent package's subpackages as part of itself
so when the parent package is explicitly enabled or disabled, it is identical
to explicitly enabling or disabling all of its subpackages (see
`Enable/disable of parent package is enable/disable for subpackages`_).  Also,
subpackages are tested along with their peer subpackages with the parent
package as part of `TriBITS CTest/CDash Driver`_.  This effectively means that if a
build failure is detected in any subpackage, then that will effectively
disable the parent package and all of its other subpackages in downstream
testing.  This is a type of "all for one and all for one" when it comes to the
relationship between the subpackages within a single parent package.  These
are some of the issues to consider when breaking up software into packages and
subpackages that will be mentioned in other sections as well.

TriBITS TPL
+++++++++++

A *TriBITS TPL*:

* Defines a set of pre-built libraries and/or header files and/or executables
  and/or some other resources used by one or more TriBITS Packages for which
  TriBITS publishes the list of include directories and/or libraries and/or
  executables provided by the TPL to the TriBITS project.
* Has a globally unique name (see `Globally unique TriBITS TPL names`_) that
  is declared in a `<repoDir>/TPLsList.cmake`_ file.
* Is listed as an explicit optional or required dependency in one or more
  TriBITS SE package's `<packageDir>/cmake/Dependencies.cmake`_ files.

.. _Globally unique TriBITS TPL names:

**WARNING:** One must be very careful to pick **Globally unique TriBITS TPL
names** across all TPLs in all TriBITS repositories that ever might be cobbled
together into a single TriBITS (meta) project!  However, choosing TPL names is
usually much easier and less risky than choosing `Globally unique TriBITS
package names`_ because the native TPLs themselves tend to be uniquely named.
For example, the TPL names ``BLAS`` and ``LAPACK`` are well defined in the
applied math and computational science community and are not likely to clash.

Using a TriBITS TPL is to be preferred over using raw CMake
``FIND_PACKAGE(<someCMakePackage>)`` because the TriBITS system guarantees
that only a single unique version of TPL of the same version will be used by
multiple packages.  Also, by defining a TriBITS TPL, automatic enable/disable
logic will be applied as described in `Package Dependencies and Enable/Disable
Logic`_.  For example, if a TPL is explicitly disabled, all of the downstream
packages that depend on that TPL will be automatically disabled as well (see
`TPL disable triggers auto-disables of downstream dependencies`_).

For each TPL referenced in a `<repoDir>/TPLsList.cmake`_ file using the macro
`TRIBITS_REPOSITORY_DEFINE_TPLS()`_, there must exist a file, typically called
``FindTPL${TPL_NAME}.cmake``, that once processed, produces the variables
``TPL_${TPL_NAME}_LIBRARIES`` and ``TPL_${TPL_NAME}_INCLUDE_DIRS``.  Most
``FindTPL${TPL_NAME}.cmake`` files just use the function
`TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()`_ the define the TriBITS TPL.
A simple example of such a file is the common TriBITS ``FindTPLPETSC.cmake``
module which is currently:

.. include:: ../../common_tpls/FindTPLPETSC.cmake
   :literal:

Some concrete ``FindTPL${TPL_NAME}.cmake`` files actually do use
``FIND_PACKAGE()`` and a standard CMake package find module to fill in the
guts of finding at TPL which is perfectly fine.  In this case, the purpose for
the wrapping ``FindTPL${TPL_NAME}.cmake`` is to standardize the output
variables ``TPL_${TPL_NAME}_INCLUDE_DIRS`` and ``TPL_${TPL_NAME}_LIBRARIES``.
For more details on properly using ``FIND_PACKAGE()`` to define a
``FindTPL${TPL_NAME}.cmake`` file, see `How to use FIND_PACKAGE() for a
TriBITS TPL`_.

Once the `<repoDir>/TPLsList.cmake`_ files are all processed, then each
defined TPL ``TPL_NAME`` is assigned the following global non-cache variables:

  .. _${TPL_NAME}_FINDMOD:

  ``${TPL_NAME}_FINDMOD``

    Gives the location for the TPL's find module.  This is set using the
    ``FINDMOD`` field in the call to `TRIBITS_REPOSITORY_DEFINE_TPLS()`_.  The
    final value of the variable is defined by the *last*
    `<repoDir>/TPLsList.cmake`_ file that is processed that declares the TPL
    ``TPL_NAME``.  For example, if ``Repo1/TPLsList.cmake`` and
    ``Repo2/TPLsList.cmake`` both list the TPL ``SomeTpl``, then if ``Repo2``
    is processed after ``Repo1``, then ``SomeTpl_FINDMOD`` is determined by
    ``Repo2/TPLsList.cmake`` and the find module listed in
    ``Repo1/TPLsList.cmake`` is ignored.

  .. _${TPL_NAME}_TESTGROUP:

  ``${TPL_NAME}_TESTGROUP``

    Gives the TPL's `SE Package Test Group`_. This is set using the
    ``CLASSIFICATION`` field in the call to
    `TRIBITS_REPOSITORY_DEFINE_TPLS()`_.  If multiple repos define a given
    TPL, then the *first* `<repoDir>/TPLsList.cmake`_ file that is processed
    that declares the TPL ``TPL_NAME`` specifies the test group.  For example,
    if ``Repo1/TPLsList.cmake`` and ``Repo2/TPLsList.cmake`` both list the TPL
    ``SomeTpl``, then if ``Repo2`` is processed after ``Repo1``, then
    ``SomeTpl_TESTGROUP`` is determined by ``Repo1/TPLsList.cmake`` and the
    test group in ``Repo2/TPLsList.cmake`` is ignored.  However, if
    ``${TPL_NAME}_TESTGROUP`` is already set before the
    `<repoDir>/TPLsList.cmake`_ files are processed, then that test group will
    be used.  Therefore, the project can override the test group for a given
    TPL if desired.

The specification given in `Enabling support for an optional Third-Party
Library (TPL)`_ and `TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()`_ describes
how the a ``FindTPL${TPL_NAME}.cmake`` module should behave and allow users to
override and specialize how a TPL's include directories and libraries are
determined.  However, note that the TriBITS system does not require the usage
of the function ``TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()`` and does not
even care about the TPL module name ``FindTPL${TPL_NAME}.cmake``.  All that is
required is that some CMake file fragment exist that once included, will
define the variables ``${TPL_NAME}_LIBRARIES`` and
``${TPL_NAME}_INCLUDE_DIRS``.  However, to be user friendly, such a CMake file
should respond to the same variables as accepted by the standard
``TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()`` function.

The core variables related to an enabled TPL are ``${TPL_NAME}_LIBRARIES``,
``${TPL_NAME}_INCLUDE_DIRS``, and ``${TPL_NAME}_TESTGROUP`` as defined in
`TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()`_ need to be defined.  For more
details, see `TRIBITS_REPOSITORY_DEFINE_TPLS()`_.

Processing of TriBITS Files: Ordering and Details
--------------------------------------------------

One of the most important things to know about TriBITS is what files it
processes, in what order, and in what context.  This is critical to being able
to understand what impact (if any) setting a variable or otherwise changing
the CMake run-time state will have on configuring a CMake project which uses
TriBITS.  While the different files that make up a `TriBITS Project`_,
`TriBITS Repository`_, `TriBITS Package`_, `TriBITS Subpackage`_, and `TriBITS
TPL`_ were defined in the section `TriBITS Structural Units`_, that material
did not fully describe the context and in what order these files are processed
by the TriBITS framework.

The TriBITS system processes the project's files in one of two general use
cases.  The first use case is in the basic configuration of the project with a
standard ``cmake`` command invocation in order to set up the build files in
the binary directory (see `Full TriBITS Project Configuration`_).  The second
use case is in reading the project's dependency-related files in order to
build the package dependency data-structure (e.g. the
`<Project>PackageDependencies.xml`_ file, see `Reduced Package Dependency
Processing`_).  The second use case of reading the project's dependency files
is largely a subset of the first.

Another factor that is important to understand is the scoping in which the
various files are processed (with ``INCLUDE()`` or ``ADD_SUBDIRECTORY()``).
This scoping has a large impact on the configuration of the project and what
effect the processing of files and setting variables have on the project as a
whole.  Some of the strange scoping rules for CMake are discussed in `CMake
Language Overview and Gotchas`_ and should be understood before trying to
debug issues with processing.  Many of the basic files are processed
(included) in the base project `<projectDir>/CMakeLists.txt`_ scope and
therefore any local variables set in these files are accessible to the entire
CMake project (after the file is processed, of course).  Other files get
processed inside of functions which have their own local scope and therefore
only impact the rest of the project in more purposeful ways.  And of course
all of the package `<packageDir>/CMakeLists.txt`_ files that are processed
using ``ADD_SUBDIRECTORY()`` create a new local scope for that given package.


Full TriBITS Project Configuration
++++++++++++++++++++++++++++++++++

The first use case to describe is the full processing of all of the TriBITS
project's files starting with the base `<projectDir>/CMakeLists.txt`_ file.
This begins with the invocation of the following CMake command to generate the
project's build files::

  $ cmake [options] <projectDir>

Below, is a short pseudo-code algorithm for the TriBITS framework processing
and callbacks that begins in the `<projectDir>/CMakeLists.txt`_ file and
proceeds through the call to `TRIBITS_PROJECT()`_.

.. _Full Processing of TriBITS Project Files:

**Full Processing of TriBITS Project Files:**

| 1.  Read `<projectDir>/ProjectName.cmake`_ (sets `PROJECT_NAME`_)
| 2.  Call ``PROJECT(${PROJECT_NAME} NONE)`` (sets `${PROJECT_NAME}_SOURCE_DIR`_
|     and `${PROJECT_NAME}_BINARY_DIR`_)
| 3.  Call `TRIBITS_PROJECT()`_:
|   1)  Set `PROJECT_SOURCE_DIR`_ and `PROJECT_BINARY_DIR`_
|   2)  For each ``<optFile>`` in ${`${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE`_}
|         then in ${`${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_APPEND`_}
        :
|       * ``INCLUDE(<optFile>)``
|   3)  Set variables ``CMAKE_HOST_SYSTEM_NAME`` and ``${PROJECT_NAME}_HOSTNAME``
|       (both of these can be overridden in the cache by the user)
|   4)  Find Python (sets ``PYTHON_EXECUTABLE``, see `Python Support`_)
|   5)  ``INCLUDE(`` `<projectDir>/Version.cmake`_ ``)``
|   6)  Define primary TriBITS options and read in the list of extra repositories
|       (calls ``TRIBITS_DEFINE_GLOBAL_OPTIONS_AND_DEFINE_EXTRA_REPOS()``)
|       * ``INCLUDE(`` `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ ``)``      
|   7)  For each ``<repoDir>`` in all defined TriBITS repositories:
|       * ``INCLUDE(`` `<repoDir>/cmake/CallbackSetupExtraOptions.cmake`_ ``)``
|       * Call macro ``TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS()``
|   9)  Call ``TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()``:
|     a)  For each ``<repoDir>`` in all defined TriBITS repositories:
|         * ``INCLUDE(`` `<repoDir>/PackagesList.cmake`_ ``)`` and process list
|         * ``INCLUDE(`` `<repoDir>/TPLsList.cmake`_ ``)`` and process list
|     b)  For each ``<repoDir>`` in all defined TriBITS repositories:
|         * ``INCLUDE(`` `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ ``)``
|     c)  ``INCLUDE(`` `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_ ``)``
|     d)  For each ``<packageDir>`` in all defined top-level packages:
|         * ``INCLUDE(`` `<packageDir>/cmake/Dependencies.cmake`_ ``)``
|           - Sets all package-specific options (see `TriBITS Package Cache Variables`_)
|         * For each ``<spkgDir>`` in all subpackages for this package:
|           * ``INCLUDE(`` `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_ ``)``
|             - Sets all subpackage-specific options
|   10) Adjust SE package and TPLs enables and disables
|       (see `Package Dependencies and Enable/Disable Logic`_)
|   11) `Probe and set up the environment`_ (finds MPI, compilers, etc.)
|       (see `TriBITS Environment Probing and Setup`_)
|       * ``INCLUDE(`` `<projectDir>/cmake/ProjectCompilerPostConfig.cmake`_ ``)``
|   12) For each ``<tplName>`` in the set of enabled TPLs:
|       * ``INCLUDE(${<tplName>_FINDMOD})`` (see `TriBITS TPL`_)
|   13) For each ``<repoDir>`` in all defined TriBITS repositories:
|       * Read `<repoDir>/Copyright.txt`_
|       * ``INCLUDE(`` `<repoDir>/Version.cmake`_ ``)``
|       (see `Project and Repository Versioning and Release Mode`_)
|   14) For each ``<packageDir>`` in all enabled top-level packages
|       * ``ADD_SUBDIRECTORY(`` `<packageDir>/CMakeLists.txt`_ ``)``
|       * For each ``<spkgDir>`` in all enabled subpackages for this package:
|         * ``ADD_SUBDIRECTORY(`` `<packageDir>/<spkgDir>/CMakeLists.txt`_ ``)``
|   16) ``INCLUDE(`` `<projectDir>/cmake/CallbackDefineProjectPackaging.cmake`_ ``)``
|       * Call ``TRIBITS_PROJECT_DEFINE_PACKAGING()``
|   16) For each ``<repoDir>`` in all defined TriBITS repositories:
|       * ``INCLUDE(`` `<repoDir>/cmake/CallbackDefineRepositoryPackaging.cmake`_ ``)``
|       * Call ``TRIBITS_REPOSITORY_DEFINE_PACKAGING()``

The TriBITS Framework obviously does a lot more than what is described above
but the basic trace of major operations and ordering and the processing of
project, repository, package, and subpackage files should be clear.  All of
this information should also be clear when enabling `File Processing
Tracing`_ and watching the output from the ``cmake`` configure STDOUT.

Reduced Package Dependency Processing 
++++++++++++++++++++++++++++++++++++++

In addition to the full processing that occurs as part of the `Full TriBITS
Project Configuration`_, there are also TriBITS tools that only process as
subset of project's files.  This reduced processing is performed in order to
build up the project's package dependencies data-structure and to write the
file `<Project>PackageDependencies.xml`_.  For example, the tool
`checkin-test.py`_ and the function `TRIBITS_CTEST_DRIVER()`_ both drive this
type of processing.  In particular, the CMake -P script
`TribitsDumpDepsXmlScript.cmake`_ reads all of the project's
dependency-related files and dumps out the `<Project>PackageDependencies.xml`_
file (see `TriBITS Project Dependencies XML file and tools`_).  This reduced
processing (e.g. as executed in ``cmake -P TribitsDumpDepsXmlScript.cmake``)
is described below.

.. _Reduced Dependency Processing of TriBITS Project Files:

**Reduced Dependency Processing of TriBITS Project:**

| 1.  Read `<projectDir>/ProjectName.cmake`_ (sets `PROJECT_NAME`_)
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
|   e) Write the file `<Project>PackageDependencies.xml`_
|      (specified using `${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE`_)

When comparing the above reduced dependency processing to the `Full Processing
of TriBITS Project Files`_ it is important to note that that several files are
**not** processed in the reduced algorithm shown above.  The files that are
**not** processed include `<projectDir>/Version.cmake`_,
`<repoDir>/Version.cmake`_ and
`<repoDir>/cmake/CallbackSetupExtraOptions.cmake`_ (in addition to not
processing any of the ``CMakeLists.txt`` files obviously).  Therefore, one
cannot put anything in these non-processed files that would impact the
definition of TriBITS repositories, packages, TPLs, etc.  Anything that would
affect the dependencies data-structure that gets written out as
`<Project>PackageDependencies.xml`_ must be contained in the files that are
processed shown in the reduced processing above.

Debugging issues with `Reduced Dependency Processing of TriBITS Project
Files`_ is more difficult because one cannot easily turn on `File Processing
Tracing`_ like one can when doing the full CMake configure.  However, options
may be added to the various tools to show this file processing and help debug
problems.

File Processing Tracing
+++++++++++++++++++++++

In order to aid in debugging problems with `Full TriBITS Project
Configuration`_ and `Reduced Package Dependency Processing`_, TriBITS defines
the CMake cache option `${PROJECT_NAME}_TRACE_FILE_PROCESSING`_.  When
enabled, TriBITS will print out when any of the project-related,
repository-related, or package-related file is being processed by TriBITS.
When ``${PROJECT_NAME}_TRACE_FILE_PROCESSING=ON``, lines starting with ``"--
File Trace:"`` are printed to ``cmake`` stdout for files that TriBITS
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
    <tribitsDir>/doc/TribitsExampleProject \
    | grep "^-- File Trace:"

looks something like::

  -- File Trace: PROJECT    INCLUDE    [...]/Version.cmake
  -- File Trace: REPOSITORY INCLUDE    [...]/cmake/CallbackSetupExtraOptions.cmake
  -- File Trace: REPOSITORY INCLUDE    [...]/PackagesList.cmake
  -- File Trace: REPOSITORY INCLUDE    [...]/TPLsList.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/simple_cxx/cmake/Dependencies.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/mixed_lang/cmake/Dependencies.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/with_subpackages/cmake/Dependencies.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/with_subpackages/a/cmake/Dependencies.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/with_subpackages/b/cmake/Dependencies.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/with_subpackages/c/cmake/Dependencies.cmake
  -- File Trace: PACKAGE    INCLUDE    [...]/packages/wrap_external/cmake/Dependencies.cmake
  -- File Trace: PROJECT    CONFIGURE  [...]/cmake/ctest/CTestCustom.cmake.in
  -- File Trace: REPOSITORY READ       [...]/Copyright.txt
  -- File Trace: REPOSITORY INCLUDE    [...]/Version.cmake
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/simple_cxx/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/simple_cxx/test/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/mixed_lang/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/mixed_lang/test/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/with_subpackages/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/with_subpackages/a/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/with_subpackages/a/tests/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/with_subpackages/b/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/with_subpackages/b/tests/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/with_subpackages/c/CMakeLists.txt
  -- File Trace: PACKAGE    ADD_SUBDIR [...]/packages/with_subpackages/c/tests/CMakeLists.txt
  -- File Trace: PROJECT    INCLUDE    [...]/cmake/CallbackDefineProjectPackaging.cmake
  -- File Trace: REPOSITORY INCLUDE    [...]/cmake/CallbackDefineRepositoryPackaging.cmake

However, every file that TriBITS processes is not printed in this file trace
if it should be obvious that the file is being processed.  For example, the
package's configured header file created using `TRIBITS_CONFIGURE_FILE()`_
does not result in a file trace print statement because this is an
unconditional command that is explicitly called in one of the package's
``CMakeLists.txt`` files so it should be clear that this file is being
processed and exactly when it is processed.


Coexisting Projects, Repositories, and Packages
-----------------------------------------------

Certain simplifications are allowed when defining TriBITS projects,
repositories and packages.  The known allowed simplifications are described
below.

**TriBITS Repository Dir == TriBITS Project Dir**: It is allowed for a TriBITS
Project and a TriBITS Repository to be the same source directory and in fact
this is the default for every TriBITS project (unless the
`<projectDir>/cmake/NativeRepositoriesList.cmake`_ is defined).  In this case,
the repository name `REPOSITORY_NAME`_ and the project name `PROJECT_NAME`_
are the same as well.  This is quite common and is in fact the default that
every TriBITS Project is also a TriBITS repository (and therefore must contain
`<repoDir>/PackagesList.cmake`_ and `<repoDir>/TPLsList.cmake`_ files).  This
is the case, for example, with the `Trilinos`_ and the
`TribitsExampleProject`_ projects and repositories.  In this case, the
Project's and the Repository's ``Version.cmake`` and ``Copyright.txt`` files
are also one and the same, as they should be (see `Project and Repository
Versioning and Release Mode`_).

**TriBITS Package Dir == TriBITS Repository Dir**: It is also allowed for a
TriBITS Repository to have only one package and to have that package be the
base repository directory.  The TriBITS Repository and the single TriBITS
Package would typically have the same name in this case (but that is actually
not required but it is confusing if they are not the same).  For example, in
the TriBITS test project `MockTrilinos`_, the repository and package
``extraRepoOnePackage`` are the same directory.  In this case, the file
``extraRepoOnePackage/PackagesList.cmake`` looks like:

.. include:: ../../examples/MockTrilinos/extraRepoOnePackage/PackagesList.cmake
   :literal:

(Note the dot ``'.'`` for the package directory.)

This is also how the real TriBITS repository and package `DataTransferKit`_ is
set up (at least that is the way it was when this document was first written).

.. _DataTransferKit: https://github.com/CNERG/DataTransferKit

However, to maximize flexibility, it is recommended that a TriBITS package and
its TriBITS repository **not** share the same directory or the same name.
This allows a TriBITS repository to define more packages later.

.. _TriBITS Package == TriBITS Repository == TriBITS Project:

**TriBITS Package Dir == TriBITS Repository Dir == TriBITS Project Dir**: In
the extreme, it is possible to collapse a single TriBITS package, repository,
and project into the same base source directory.  They can also share the same
name for the package, repository and package.  One example of this is the
``TriBITS`` project and `The TriBITS Test Package`_ themselves, which are both
rooted in the base ``TriBITS/`` source directory of the stand-alone TriBITS
repository.  There are a few restrictions and modifications needed to get this
to work:

* The base ``CMakeLists.txt`` file must be modified to allow it to be
  processed both as the base project ``CMakeLists.txt`` file and as the
  package's base ``CMakeLists.txt`` file.  In the case of
  ``TriBITS/CMakeLists.txt``, a big if statement is used.

.. ToDo: Create a new macro called TRIBITS_PROJECT_AND_PACKAGE() that will
.. automatically take care of the details of a TriBITS package also being a
.. TriBITS project.

Other than that simple modification to the top-level ``CMakeLists.txt`` file,
a TriBITS project, repository, and package can all be rooted in the same
source directory.

The primary use case for collapsing a project, repository, and package into a
single base source directory would be to support the stand-alone build of a
TriBITS package as its own entity that uses an independent installation of the
TriBITS (or a minimal snapshot of TriBITS).  If a given TriBITS package has no
required `upstream`_ TriBITS package dependencies and minimal TPL dependencies
(or only uses `Standard TriBITS TPLs`_ or `Common TriBITS TPLs`_ already
defined in the ``tribits/core/std_tpls/`` or ``tribits/common_tpls/``
directories), then creating a stand-alone project build of a single TriBITS
package requires fairly little extra overhead or duplication.


Standard and Common TPLs
------------------------

While a TriBITS Repository can define their own TPLs and their own TPL find
modules (see `TriBITS TPL`_), the TriBITS source tree contains the find
modules for a few different standard TPLs and common TPLs.  `Standard TriBITS
TPLs`_ are integral to the TriBITS system itself while `Common TriBITS TPLs`_
are TPL that are used in several different TriBITS Repositories and are
contained in TriBITS for convenience and uniformity.


Standard TriBITS TPLs
+++++++++++++++++++++

TriBITS contains find modules for a few standard TPLs integral to the TriBITS
system.  The standard TriBITS TPLs are contained under the directory::

  tribits/core/std_tpls/

The current list of standard TriBITS TPL find modules is:

.. include:: TribitsStandardTPLsList.txt
   :literal:

The TPLs ``MPI`` and ``CUDA`` are standard because they are special in that
they define compilers and other special tools that are used in
`TRIBITS_ADD_LIBRARY()`_, `TRIBITS_ADD_EXECUTABLE()`_, `TRIBITS_ADD_TEST()`_
and other commands.

These standard TPLs are used in a `<repoDir>/TPLsList.cmake`_ file as::

  TRIBITS_REPOSITORY_DEFINE_TPLS(
    MPI   "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"  PT
    CUDA  "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"  ST
    ...
    )


Common TriBITS TPLs
+++++++++++++++++++

TriBITS also contains find modules for several TPLs that are used across many
independent TriBITS repositories.  The goal of maintaining these under TriBITS
is to enforce conformity in case these independent repositories are combined
into a single meta-project.

The common TriBITS TPLs are contained under the directory::

  tribits/common_tpls/

The current list of common TriBITS TPL find modules is:

.. include:: TribitsCommonTPLsList.txt
   :literal:

Common TPLs are used in a `<repoDir>/TPLsList.cmake`_ file as::

  TRIBITS_REPOSITORY_DEFINE_TPLS(
    BLAS   "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"  PT
    LAPACK  "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"  PT
    ...
    )

By using a standard TPL definition, it is guaranteed that the TPL used will be
consistent with all of the TriBITS packages that depend on these TPLs in case
they are combined into a single project.

Note that just because packages in two different TriBITS repositories
reference the same TPL does not necessarily mean that it needs to be moved
into the TriBITS source tree under ``tribits/common_tpls/``.  For example, if
the TPL ``QT`` is defined in an `upstream`_ repository (e.g. Trilinos), then a
package in a `downstream`_ repository can list a dependency on the TPL ``QT``
without having to define its own ``QT`` TPL in its repository's
`<repoDir>/TPLsList.cmake`_ file.  For more details, see `TriBITS TPL`_.

.. Where to set variables?
.. -----------------------

.. A TriBITS project, repository, package have a number of files in which
   variables can be set and commands can be called in order to affect how the
   project is defined and configured.  With so many files that can be included by
   TriBITS, it can be difficult to know where the right place is to set a given
   set of variables.  The primary considerations for where to set a variable
   depend on:
 
.. ToDo: Describe considerations on where to set variables ...


Example TriBITS Projects
=========================

In this section, a few different example TriBITS projects and packages are
previewed.  Most of these examples exist in the TriBITS source directory
``tribits`` itself so they are available to all users of TriBITS.  These
examples also provide a means to test the TriBITS system itself (see `The
TriBITS Test Package`_).

The first example covered is the bare bones `TribitsHelloWorld`_ example
project.  The second example covered in detail is `TribitsExampleProject`_.
This example covers all the basics for setting up a simple multi-package
TriBITS project.  The third example outlined is `MockTrilinos`_ which mostly
exists to test the TriBITS system itself but also contains some nice examples
of a few different TriBITS features and behaviors.  The forth example is the
`ReducedMockTrilinos`_ project which is used to demonstrate TriBITS behavior
in this document.  Also mentioned is the `Trilinos`_ project itself which can
be a useful example of the usage of TriBITS (see disclaimers in the section
`Trilinos`_).  The last example mentioned is `The TriBITS Test Package`_
itself which allows the TriBITS system to be tested and installed from any
TriBITS project that lists it, including the ``TriBITS`` project itself
(see `Coexisting Projects, Repositories, and Packages`_).

The directory ``tribits/examples/`` contains some other example TriBITS
projects and repositories as well that are referred to in this and other
documents.


TribitsHelloWorld
-----------------

``TribitsHelloWorld`` is about the simplest possible TriBITS project that you
can imagine and is contained under the directory::

  tribits/examples/TribitsHelloWorld/

This example project contains only a single TriBITS package and no frills at
all (does not support MPI or Fortran).  However, it does show how minimal a
`TriBITS Project`_ (which is also a `TriBITS Repository`_) and a `TriBITS
Package`_ can be and still demonstrates some of the value of TriBITS over raw
CMake.  The simple ``HelloWorld`` package is used to compare with the raw
CMakeLists.txt file in the ``RawHeloWorld`` example project in the `TriBITS
Overview`_ document.

The directory structure for this example his shown below:

.. include:: TribitsHelloWorldDirAndFiles.txt
   :literal:

This has all of the required `TriBITS Project Core Files`_, `TriBITS
Repository Core Files`_, and `TriBITS Package Core Files`_.  It just builds a
simple library, a simple executable, a test executable, and the tests them as
shown by the file ``TribitsHelloWorld/hello_world/CMakeLists.txt`` which is:

.. include:: ../../examples/TribitsHelloWorld/hello_world/CMakeLists.txt
   :literal:

The build and test of this simple project is tested in the `The TriBITS Test
Package`_ file::

  TriBITS/test/core/ExamplesUnitTests/CMakeLists.txt

Note that this little example is a fully functional `TriBITS Repository`_ and
can be embedded in to a larger TriBITS meta-project and be seamlessly built
along with any other such TriBITS-based software.

.. ToDo: Put in reference to the example meta-project.


TribitsExampleProject
----------------------

``TribitsExampleProject`` in an example `TriBITS Project`_ and `TriBITS
Repository`_ contained in the TriBITS source tree under::

  tribits/examples/TribitsExampleProject/

When this is used as the base TriBITS project, this is the directory
corresponds to ``<projectDir>`` and ``<repoDir>`` referenced in `TriBITS
Project Core Files`_ and `TriBITS Repository Core Files`_, respectively.

Several files from this project are used as examples in the section `TriBITS
Project Structure`_.  Here, a fuller description is given of this project and
a demonstration of how TriBITS works.  From this simple example project, one
can quickly see how the basic structural elements of a TriBITS project,
repository, and package (and subpackage) are pulled together.

This simple project shows how what is listed in files:

* `<repoDir>/PackagesList.cmake`_,
* `<packageDir>/cmake/Dependencies.cmake`_, and
* `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_

are used to specify the packages and SE packages in a TriBITS project and
repository.  More details about the contents of the ``Dependencies.cmake``
files is described in the section `Package Dependencies and Enable/Disable
Logic`_.

The name of this project ``PROJECT_NAME`` is given in its
``TribitsExampleProject/ProjectName.cmake`` file:

.. include:: ../../examples/TribitsExampleProject/ProjectName.cmake
   :literal:

The variable ``PROJECT_NAME=TribitsExProj`` is used to prefix (using
``"${PROJECT_NAME}_"``) all of the project's global TriBITS variables like
``TribitsExProj_ENABLE_TESTS``, ``TribitsExProj_ENABLE_ALL_PACKAGES``, etc.

.. _TribitsExampleProject Files and Directories:

The directory structure and key files for this example project are shown in
the partial list of `TribitsExampleProject Files and Directories`_ below::

  TribitsExampleProject/
    CMakeLists.txt
    Copyright.txt
    PackagesList.cmake
    ProjectName.cmake
    project-checkin-test-config.py
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
      mixed_lang/ ...
      with_subpackages/
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

Above, the sub-directories under ``packages/`` are sorted according to the
order listed in the ``TribitsExampleProject/PackagesList.cmake`` file:

.. include:: ../../examples/TribitsExampleProject/PackagesList.cmake
   :literal:

From this file, we get the list of top-level packages ``SimpleCxx``,
``MixedLang``, ``WithSubpackages``, and ``WrapExternal`` (and their base
package directories and testing group, see `<repoDir>/PackagesList.cmake`_).
(NOTE: By default the package ``InsertedPkg`` is not defined because its
directory is missing, see `How to insert a package into an upstream repo`_.)

A full listing of package files in `TribitsExampleProject Files and
Directories`_ is only shown for the ``SimpleCxx`` package directory
``packages/simple_cxx/``.  For this package, ``<packageDir> =
<repoDir>/packages/simple_cxx`` and ``PACKAGE_NAME = SimpleCxx``.  As
explained in `TriBITS Package Core Files`_, the files
`<packageDir>/cmake/Dependencies.cmake`_ and `<packageDir>/CMakeLists.txt`_
must exist for every package directory listed in
`<repoDir>/PackagesList.cmake`_ and we see these files under in the directory
``packages/simple_cxx/``.  The package ``SimpleCxx`` does not have any
`upstream`_ SE package dependencies.

Now consider the example top-level package ``WithSubpackages`` which,
as the name suggests, is broken down into subpackages.  The
``WithSubpackages`` dependencies file::

  TribitsExampleProject/packages/with_subpackages/cmake/Dependencies.cmake

with contents:

.. include:: ../../examples/TribitsExampleProject/packages/with_subpackages/cmake/Dependencies.cmake
   :literal:

references the three subpackages with sub-directories ``<spkgDir>`` = ``A``,
``B``, and ``C`` under the parent package directory
``packages/package_with_packages/`` which are shown in `TribitsExampleProject
Files and Directories`_.  This gives another set of three SE packages
``WithSubpackagesA``, ``WithSubpackagesB``,
and ``WithSubpackagesC``.  Combining ``<packageDir> =
packages/package_with_packages`` and ``<spkgDir>`` for each subpackage gives
the subpackage directories::

  TribitsExampleProject/packages/with_subpackages/a/
  TribitsExampleProject/packages/with_subpackages/b/
  TribitsExampleProject/packages/with_subpackages/c/

Together with the top-level parent SE package ``WithSubpackages``
itself, this top-level package provides four `TriBITS SE Packages`_ giving the
final list of SE packages provided by this TriBITS repo as::

  SimpleCxx MixedLang WithSubpackagesA WithSubpackagesB WithSubpackagesC \
    WithSubpackages WrapExternal 7

The above list of SE packages is printed (with the number of packages printed
at the end) by TriBITS to the ``cmake`` stdout on the line starting with
``"Final set of non-enabled SE packages:"`` when no packages are enabled (see
`Selecting the list of packages to enable`_).  (Note that TriBITS does not put
in line-brakes with continuation characters ``"\"`` as shown above.)  TriBITS
defines enable/disable cache variables for each of these defined SE packages
like ``TribitsExProj_ENABLE_SimpleCxx`` and
``TribitsExProj_ENABLE_WithSubpackagesA``, and defines all
the variables listed in `TriBITS Package Cache Variables`_ that are settable
by the users or by the dependency logic described in section `Package
Dependencies and Enable/Disable Logic`_.

When starting a new TriBITS project, repository, or package, one should
consider basing these on the examples in ``TribitsExampleProject``.  In fact,
the skeletons for any of the

* `TriBITS Project Core Files`_,
* `TriBITS Repository Core Files`_,
* `TriBITS Package Core Files`_, or
* `TriBITS Subpackage Core Files`_

should be copied from this example project as they represent best practice
when using TriBITS for the typical use cases.


MockTrilinos
-------------

The TriBITS project ``MockTrilinos`` is contained under the directory::

  tribits/examples/MockTrilinos/

This TriBITS project is not a full TriBITS project (i.e. it does not build
anything).  Instead, it is primarily used to test the TriBITS system using
tests defined in the `The TriBITS Test Package`_.  The ``MockTrilinos``
project is actually given the name ``PROJECT_NAME = Trilinos`` and contains a
subset of Trilinos packages with slightly modified dependencies from a real
version of the Trilinos project from May 2009.  The list of packages in::

  tribits/examples/MockTrilinos/PackagesList.cmake

is:

.. include:: ../../examples/MockTrilinos/PackagesList.cmake
   :literal:

All of the package directories listed above have
`<packageDir>/cmake/Dependencies.cmake`_ files but generally do not have
`<packageDir>/CMakeLists.txt`_ files since most of usage of ``MockTrilinos``
just involves the testing of the algorithms and behaviors described in the
section `Package Dependencies and Enable/Disable Logic`_.

``MockTrilinos`` also contains a number of extra TriBITS repositories used in
various tests.  These extra repositories offer examples of different types of
TriBITS repositories like:

* ``extraRepoOnePackage``: Contains just the single package
  ``extraRepoOnePackage`` which is defined in the base repository directory.

* ``extraRepoOnePackageThreeSubpackages``: Contains just the single package
  ``extraRepoOnePackageThreeSubpackages`` which is defined in the base
  repository directory but is broken up into subpackages.

* ``extraRepoTwoPackages``: Contains just two packages but provides an example
  of defining multiple repositories with possible missing required and
  optional `upstream`_ packages (see `Multi-Repository Support`_).

* ``extraTrilinosRepo``: Just a typical extra repo with add-on packages and
  TPLs that depends on a few upstream ``MockTrilinos`` packages.

New test extra repositories are added when new types of tests are needed that
would require new package and TPL dependency structures since existing
dependency tests based on ``MockTrilinos`` are expensive to change by their
very nature.

The primary reason that the ``MockTrilinos`` test project is mentioned in this
developers guide is because it contains a variety of packages, subpackages,
and TPLs with a variety of different types of dependencies.  This variety is
is needed to more fully test the TriBITS system but this project and the tests
also serve as examples and extra documentation for the behavior of the TriBITS
system.  Several of the dependency-related examples referenced in this
document come from ``MockTrilinos``.

Most of the dependency tests involving ``MockTrilinos`` are specified in::

  TriBITS/test/core/DependencyUnitTests/CMakeLists.txt

A great deal about the current behavior of TriBITS `Package Dependencies and
Enable/Disable Logic`_ can be learned from inspecting these tests.  There are
also some faster-running unit tests involving ``MockTrilinos`` defined in the
file::

  TriBITS/test/core/TribitsAdjustPackageEnables_UnitTests.cmake


ReducedMockTrilinos
-------------------

The TriBITS project ``ReducedMockTrilinos`` is contained under the directory::

  tribits/examples/ReducedMockTrilinos/

It is a scaled-down version of the `MockTrilinos`_ test project with just a
handful of packages and some modified dependencies.  Its primary purpose for
this example project is to be used for examples in the section `Package
Dependencies and Enable/Disable Logic`_ and to test a few features of the
TriBITS system not tested in other tests.

The list of packages in::

  tribits/examples/ReducedMockTrilinos/PackagesList.cmake

is:

.. include:: ../../examples/ReducedMockTrilinos/PackagesList.cmake
   :literal:

All of the listed packages are standard TriBITS packages except for the mock
``Thyra`` package which is broken down into subpackages.  More details of this
example project are described in `Package Dependencies and Enable/Disable
Logic`_.

Trilinos
--------

The real Trilinos project and repository itself is an advanced example for the
usage of TriBITS.  Almost every single-repository use case for TriBITS is
demonstrated somewhere in Trilinos.  While some of the usage of TriBITS in
Trilinos may not be not exemplary (e.g., because it represents old usage, or
was written by CMake/TriBITS beginners) it does represent real working usage.
Given that Trilinos is a widely available software repository, anyone should be
able to access a newer version of Trilinos and mine it for CMake and TriBITS
examples.


The TriBITS Test Package
------------------------

The last TriBITS example mentioned here is the TriBITS test package named
(appropriately) ``TriBITS`` itself defined in the ``TriBITS`` repository.  The
directory for the ``TriBITS`` test package is the base TriBITS source
directory ``tribits``.  This allows any TriBITS project to add testing for the
TriBITS system by just listing the TriBITS repository in its
`<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file.  Trilinos lists the
TriBITS repository in its ``ExtraRepositoriesList.cmake`` file as:

  TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES(
    TriBITS  ""  GIT  https://github.com/TriBITSPub/TriBITS  ""  Continuous
    ...
    )

No `downstream`_ TriBITS packages list a dependency on ``TriBITS`` in their
`<packageDir>/cmake/Dependencies.cmake`_ files.  Defining the ``TriBITS``
TriBITS package in only done for running the TriBITS tests.

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
partitioning of a large software project into packages and managing the
dependencies between these packages to support building, testing, and
deploying different pieces as needed (see discussion of `Software Engineering
Packaging Principles`_).  This is especially useful in incremental CI testing
of large projects.  However, maintaining such dependencies is also a critical
component in creating and maintaining Self-Sustaining Software (see the
`TriBITS Lifecycle Model`_).  The fundamental mechanism for breaking up a
large software into manageable pieces is to partition the software into
different `TriBITS Packages`_ and then define the dependencies between these
packages (which are defined inside of the
`<packageDir>/cmake/Dependencies.cmake`_ files for each package).

Note that the basic idea of breaking up a large set of software into pieces,
defining dependencies between the pieces, and then applying algorithms to
manipulate the dependency data-structures is nothing new.  If fact, nearly
every binary package deployment system provided in various Linux OS
distributions have the concept of packages and dependencies and will
automatically install all of the necessary `upstream dependencies`_ when a
`downstream dependency`_ install is requested.  The main difference (and the
added complexity) with TriBITS is that it can handle both required and
optional dependencies since it can build from source.  A binary package
installation system, however, typically can't support optional dependencies
because only pre-built binary libraries and tools are available to install.

This section is organized and broken-down as follows.  First, the subsection
`Example ReducedMockTrilinos Project Dependency Structure`_ presents the
`ReducedMockTrilinos`_ example project, describes its dependency structure,
and uses it to begin to describe how TriBITS sets up and manages package and
TPL dependencies.  The packages in this `ReducedMockTrilinos`_ example project
are used in the following subsections so one will be constantly referred back
to this subsection.  The following subsection `TriBITS Dependency Handling
Behaviors`_ defines and describes the nitty-gritty details of the TriBITS
package and TPL dependency structure and the algorithms that manipulate the
various package, TPL, and test enables and disables.  Specific examples for
the TriBITS dependency handling algorithms are given in the subsection
`Example Enable/Disable Use Cases`_.  Finally, the subsection
`<Project>PackageDependencies.xml`_ describes the standard XML output
data-structure that gets created by TriBITS that defines a project's package
and TPL dependencies.


Example ReducedMockTrilinos Project Dependency Structure
--------------------------------------------------------

To demonstrate the TriBITS package and TPL dependency handling system, the
small simple `ReducedMockTrilinos`_ project is used.  The list of packages for
this project is defined in the file ``ReducedMockTrilinos/PackagesList.cmake``
(see `<repoDir>/PackagesList.cmake`_) which contents:

.. include:: ../../examples/ReducedMockTrilinos/PackagesList.cmake
   :literal:

All of the listed packages are standard TriBITS packages except for the mock
``Thyra`` package which is broken down into subpackages as shown in
``packages/thyra/cmake/Dependnecies.cmake`` (see
`<packageDir>/cmake/Dependencies.cmake`_) which is:

.. include:: ../../examples/ReducedMockTrilinos/packages/thyra/cmake/Dependencies.cmake
   :literal:

This gives the full list of top-level TriBITS packages::

  Teuchos RTOp Epetra Triutils EpetraExt Thyra

Adding in the subpackages defined in the top-level ``Thyra`` package, the full
set of `TriBITS SE Packages`_ for this project is::

  Teuchos RTOp Epetra Triutils EpetraExt ThyraCoreLibs ThyraGoodStuff \
    ThyraCrazyStuff ThyraEpetra ThyraEpetraExt Thyra

Note that one can see this full list of top-level packages and SE packages in
the ``cmake`` configure output lines starting with::

  Final set of non-enabled packages:
  Final set of non-enabled SE packages:

respectively, when configuring with no package enables as shown in the example
`Default configure with no packages enabled on input`_.

The list of `TriBITS TPLs`_ for this example project given in the file
``ReducedMockTrilinos/TPLsList.cmake`` (see `<repoDir>/TPLsList.cmake`_) which
is:

.. include:: ../../examples/ReducedMockTrilinos/TPLsList.cmake
   :literal:

Take note of the `SE Package Test Group`_ (i.e. `PT`_, `ST`_, or `EX`_)
assigned to each SE package and TPL as it plays a significant role in how the
TriBITS dependency system handles enables and disables.

The dependency structure of this simple TriBITS project is shown below in
`ReducedMockTrilinos Dependencies`_.

.. _ReducedMockTrilinos Dependencies:

**ReducedMockTrilinos Dependencies:**

.. include:: ReducedMockTrilinosOutput/ExpectedDependencies.txt
   :literal:

The above dependency structure printout is produced by configuring with
``${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES=ON`` (which also results in more
dependency information than what is shown above, e.g. like computed forward
package dependencies).  Note that the top-level SE package ``Thyra`` is shown
to depend on its subpackages (not the other way around).  (Many people are
confused about the nature of the dependencies between packages and
subpackages.  See `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_ for more
discussion.)

.. ToDo: Show diagram with this dependency structure.

A number of user-settable CMake cache variables determine what SE packages
(and TPLs) and what tests and examples get enabled.  These cache variables are
described in `Selecting the list of packages to enable`_ and are described
below.  Also, the assigned `SE Package Test Group`_ (i.e. `PT`_, `ST`_, and
`EX`_) also affects what packages get enabled or disabled.

Any of these SE packages can be enabled or disabled with
``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=(ON|OFF)`` (the default enable is
typically empty ``""``, see `PT/ST SE packages given default unset
enable/disable state`_).  For ``ReducedMockTrilinos``, this gives the
enable/disable cache variables (with the initial default values)::

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
defines the following cache variables as well (with the initial default
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
lots of packages where the user only enables a subset of the packages.

NOTE: TriBITS also defines the cache variables
``<TRIBITS_PACKAGE>_ENABLE_EXAMPLES`` for each enabled TriBITS package which
is handled the same way as the ``<TRIBITS_PACKAGE>_ENABLE_TEST`` variables.

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
a cache variable ``<TRIBITS_PACKAGE>_ENABLE_<OPTIONAL_DEP>``.  For the
optional dependencies shown in `ReducedMockTrilinos Dependencies`_, that gives
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
  Thyra_ENABLE_ThyraCrazyStuff=""
  Thyra_ENABLE_ThyraEpetra=""
  Thyra_ENABLE_ThyraEpetraExt=""

The above optional package-specific cache variables allow one to control
whether or not support for `upstream`_ dependency X is turned on in package Y
independent of whether or not X and Y are themselves both enabled.  For
example, if the packages ``Triutils`` and ``EpetraExt`` are both enabled, one
can explicitly disable support for the optional dependency ``Triutils`` in
``EpetraExt`` by setting ``EpetraExt_ENABLE_Triutils=OFF``.  One may want to
do this for several reasons but the bottom line is that this gives the user
more detailed control over package dependencies.  See the `TriBITS Dependency
Handling Behaviors`_ and `Explicit disable of an optional package
dependency`_ for more discussion and examples.

Before getting into specific `Example Enable/Disable Use Cases`_, some of the
`TriBITS Dependency Handling Behaviors`_ are first defined below.

TriBITS Dependency Handling Behaviors
-------------------------------------

Below, some of the rules and behaviors of the TriBITS dependency management
system are described.  Examples refer to the `Example ReducedMockTrilinos
Project Dependency Structure`_.  More detailed examples of these behaviors are
given in the section `Example Enable/Disable Use Cases`_.

In brief, the rules/behaviors of the TriBITS package and TPL dependency
management system are:

0)  `No circular dependencies of any kind are allowed`_
1)  `PT/ST SE packages given default unset enable/disable state`_
2)  `EX SE packages disabled by default`_
3)  `SE package enable triggers auto-enables of upstream dependencies`_
4)  `SE package disable triggers auto-disables of downstream dependencies`_
5)  `PT/ST TPLs given default unset enable/disable state`_
6)  `EX TPLs given default unset enable/disable state`_
7)  `Required TPLs are auto-enabled for enabled SE packages`_
8)  `Optional TPLs only enabled explicitly by the user`_
9)  `TPL disable triggers auto-disables of downstream dependencies`_
10) `Disables trump enables where there is a conflict`_
11) `Enable/disable of parent package is enable/disable for subpackages`_
12) `Subpackage enable does not auto-enable the parent package`_
13) `Support for optional SE package/TPL is enabled by default`_
14) `Support for optional SE package/TPL can be explicitly disabled`_
15) `Explicit enable of optional SE package/TPL support auto-enables SE package/TPL`_
16) `ST SE packages only auto-enabled if ST code is enabled`_
17) `<Project>_ENABLE_ALL_FORWARD_DEP_PACKAGES enables downstream packages/tests`_
18) `<Project>_ENABLE_ALL_PACKAGES enables all PT (cond. ST) SE packages`_
19) `<Project>_ENABLE_TESTS only enables explicitly enabled SE package tests`_
20) `If no SE packages are enabled, nothing will get built`_
21) `TriBITS prints all enables and disables to stdout`_
22) `TriBITS auto-enables/disables done using non-cache local variables`_

In more detail, these rules/behaviors are:

.. _No circular dependencies of any kind are allowed:

0) **No circular dependencies of any kind are allowed**: The zeroth rule of
   the TriBITS dependency management system is that no TriBITS SE package (or
   its tests) can declare a dependency on a `downstream`_ SE package, period!
   This is one of the most basic `Software Engineering Packaging Principles`_
   stated as the *ADP (Acyclic Dependencies Principle)*, i.e. "Allow no cycles
   in the package dependency graph".  By default, the TriBITS system will
   check for circular dependencies and will fail the configure if any are
   found. For a more detailed discussion of circular package dependencies, see
   `Design Considerations for TriBITS`_.

.. _PT/ST SE packages given default unset enable/disable state:

1) **PT/ST SE packages given default unset enable/disable state**: An SE
   package ``<TRIBITS_PACKAGE>`` with testing group ``PT`` or ``ST`` is given
   an **unset enable/disable state by default**
   (i.e. ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=""``).  For example, the
   ``PT`` package ``Teuchos`` is not enabled or disabled by default and is
   given the initial value ``Trilinos_ENABLE_Teuchos=""``.  For an example,
   see `Default configure with no packages enabled on input`_.  This allows
   ``PT``, and ``ST`` packages to be enabled or disabled using other logic
   defined by TriBITS which is described below.  TriBITS defines persistent
   cache variables for these with the default value of empty "".  Therefore,
   if the user or other CMake code does not hard enable or disable one of
   these variables, then on future configures it will be defined but have the
   value of empty "".

.. _EX SE packages disabled by default:

2) **EX SE packages disabled by default**: An SE package ``<TRIBITS_PACKAGE>``
   with testing group ``EX`` is **disabled by default**
   (i.e. ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=OFF``).  For an example,
   see `Default configure with no packages enabled on input`_.  This results
   in all required `downstream`_ SE packages to be disabled by default.
   However, the user can explicitly set
   ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=ON`` for an ``EX`` package and
   it will be enabled (unless one of its required dependencies are not enabled
   for some reason).  In this case, the cache variable is given the cache
   value of ``OFF``.

.. _SE package enable triggers auto-enables of upstream dependencies:

3) **SE package enable triggers auto-enables of upstream dependencies**: Any
   SE package ``<TRIBITS_PACKAGE>`` can be explicitly **enabled** by the user
   by setting the cache variable
   ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=ON``
   (e.g. ``Trilinos_ENABLE_EpetraExt=ON``).  When an SE package is enabled in
   this way, the TriBITS system will try to enable all of the required
   upstream SE packages and TPLs defined by the package (specified in its
   ``Dependencies.cmake`` file).  If an enabled SE package can't be enabled
   and has to be disabled, either a warning is printed or processing will stop
   with an error (depending on the value of
   ``${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES``, see `Disables
   trump enables where there is a conflict`_).  In addition, if
   ``${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES=ON``, then TriBITS will try
   to enable all of the specified optional `PT`_ SE packages as well (and also
   optional upstream `ST`_ SE packages if
   ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=ON``).  For an example, see
   `Explicit enable of a package, its tests, an optional TPL, with ST
   enabled`_.

.. _SE package disable triggers auto-disables of downstream dependencies:

4) **SE package disable triggers auto-disables of downstream dependencies**:
   Any SE package ``<TRIBITS_PACKAGE>`` can be explicitly **disabled** by the
   user by setting the cache variable
   ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=OFF``
   (e.g. ``Trilinos_ENABLE_Teuchos=OFF``).  When an SE package is explicitly
   disabled, it will result in the disable of all `downstream`_ SE packages
   that have required dependency on it.  It will also disable optional support
   for the disabled packages in downstream packages that list it as an
   optional dependency.  For an example, see `Explicit disable of a package`_.

.. _PT/ST TPLs given default unset enable/disable state:

5) **PT/ST TPLs given default unset enable/disable state**: A TriBITS TPL
   ``<TRIBITS_TPL>`` with testing group ``PT`` or ``ST`` is given an **unset
   enable/disable state by default** (i.e. ``TPL_ENABLE_<TRIBITS_TPL>=""``).
   For example, the ``PT`` TPL ``BLAS`` is not enabled or disabled by default
   (i.e. ``TPL_ENABLE_BLAS=""``).  For an example, see `Default configure with
   no packages enabled on input`_.  This allows ``PT``, and ``ST`` TPLs to be
   enabled or disabled using other logic.

.. _EX TPLs given default unset enable/disable state:

6) **EX TPLs given default unset enable/disable state**: A TriBITS TPL
   ``<TRIBITS_TPL>`` with testing group ``EX``, is given an **unset
   enable/disable state by default** (i.e. ``TPL_ENABLE_<TRIBITS_TPL>=""``,
   same as for ``PT`` and ``EX`` TPLs).  For an example, see `Default
   configure with no packages enabled on input`_.  This is different behavior
   than for ``EX`` SE packages described above which provides an initial hard
   disable.  However, since TriBITS will never automatically enable an
   optional TPL (see `Optional TPLs only enabled explicitly by the user`_) and
   since only `downstream`_ ``EX`` SE packages are allowed to have a required
   dependencies on an ``EX`` TPL, there is no need to set the default enable
   for an ``EX`` TPL to ``OFF``.  If an ``EX`` SE package has a required
   dependency on an ``EX`` TPL, just enabling the ``EX`` SE package should
   automatically enable the ``EX`` TPL as described in `Required TPLs are
   auto-enabled for enabled SE packages`_.

.. _Required TPLs are auto-enabled for enabled SE packages:

7) **Required TPLs are auto-enabled for enabled SE packages**: All TPLs listed
   as required TPL dependencies for the final set of enabled SE packages are
   **set to enabled** (i.e. ``TPL_ENABLE_<TRIBITS_TPL>=ON``), unless the
   listed TPLs are already explicit disabled (in which case the SE package
   would be disabled or an error would occur, see `Disables trump enables
   where there is a conflict`_).  For example, if the ``Teuchos`` package is
   enabled, then that will trigger the enable of its required TPLs ``BLAS``
   and ``LAPACK``.  For an example, see `Explicit enable of a package and its
   tests`_.

.. _Optional TPLs only enabled explicitly by the user:

8) **Optional TPLs only enabled explicitly by the user**: Optional TPLs,
   regardless of their testing group ``PT``, ``ST`` or ``EX``, will only be
   enabled if they are explicitly enabled by the user.  For example, just
   because the package ``Teuchos`` is enabled, the optional TPLs ``Boost`` and
   ``MPI`` will **not** be enabled by default, even if
   ``${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES=ON``.  To enable the
   optional TPL ``Boost``, for example, and enable support for ``Boost`` in
   the ``Teuchos`` package, the user must explicitly set
   ``TPL_ENABLE_Boost=ON``.  For an example, see `Explicit enable of a
   package, its tests, an optional TPL, with ST enabled`_.

.. _TPL disable triggers auto-disables of downstream dependencies:

9) **TPL disable triggers auto-disables of downstream dependencies**: Any TPLs
   that are explicitly disabled (i.e. ``TPL_ENABLE_<TRIBITS_TPL>=OFF``) will
   result in the disable of all `downstream`_ dependent SE packages that have
   a required dependency on the TPL.  For example, if the user sets
   ``TPL_ENABLE_LAPACK=OFF``, then this will result in the disable of SE
   packages ``Teuchos`` and ``Epetra``, and all of the required SE packages
   downstream from them.  Also, the explicitly disabled TPL will result in the
   disable of optional support in all downstream SE packages.  For example, if
   the user sets ``TPL_ENABLE_MPI=OFF``, then TriBITS will automatically set
   ``Teuchos_ENABLE_MPI=OFF`` and ``Epetra_ENABLE_MPI=OFF``.  For examples,
   see `Explicit disable of an optional TPL`_ and `Explicit disable of a
   required TPL`_.

.. _Disables trump enables where there is a conflict:

10) **Disables trump enables where there is a conflict** and TriBITS will
    never override a disable in order to satisfy some dependency.  For
    example, if the user sets ``Trilinos_ENABLE_Teuchos=OFF`` and
    ``Trilinos_ENABLE_RTOp=ON``, then TriBITS will **not** override the
    disable of ``Teuchos`` in order to satisfy the required dependency of
    ``RTOp``.  In cases such as this, the behavior of the TriBITS dependency
    adjustment system will depend on the setting of the top-level user cache
    variable `${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES`_:

    .. _${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON:
 
    * If ``${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON``, then
      TriBITS will disable the explicit enable and continue on.  In the above
      example, TriBITS will override ``Trilinos_ENABLE_RTOp=ON`` and set
      ``Trilinos_ENABLE_RTOp=OFF`` and print a verbose warning to the
      ``cmake`` stdout.

    .. _${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=OFF:

    * If ``${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=OFF``, then
      TriBITS will generate a detailed error message printed to ``cmake``
      stdout and then abort processing.  In the above example, TriBITS will
      report that ``RTOp`` is enabled but the required SE package ``Teuchos``
      is disabled and therefore ``RTOp`` can't be enabled and processing must
      stop.

    For an example of both behaviors, see `Conflicting explicit enable and
    disable`_.

.. _Enable/disable of parent package is enable/disable for subpackages:

11) **Enable/disable of parent package is enable/disable for subpackages**: An
    explicit enable/disable of a top-level parent package with subpackages
    with ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=(ON|OFF)`` is equivalent
    to the explicit enable/disable of all of the parent package's subpackages.
    For example, explicitly setting ``Trilinos_ENABLE_Thyra=ON`` is equivalent
    to explicitly setting::

      Trilinos_ENABLE_ThyraCoreLibs=ON
      Trilinos_ENABLE_ThyraGoodStuff=ON   # Only if enabling ST code!
      Trilinos_ENABLE_ThyraEpetra=ON
      Trilinos_ENABLE_ThyraEpetraExt=ON   # Only if enabling ST code!
   
    (Note that ``Trilinos_ENABLE_ThyraCrazyStuff`` is **not** set to ``ON``
    because it is already set to ``OFF`` by default, see `EX SE packages
    disabled by default`_.)  Likewise, explicitly setting
    ``Trilinos_ENABLE_Thyra=OFF`` is equivalent to explicitly setting all of
    the ``Thyra`` subpackages to ``OFF`` at the outset.  For a ``PT`` example,
    see `Explicit enable of a package and its tests`_.  For a ``ST`` example,
    see `Explicit enable of a package, its tests, an optional TPL, with ST
    enabled`_.

.. _Subpackage enable does not auto-enable the parent package:

12) **Subpackage enable does not auto-enable the parent package**: Enabling an
    SE package that is a subpackage does **not** automatically enable the
    parent package (except for at the very end, mostly just for show).  For
    example, enabling the SE package ``ThyraEpetra`` does not result in enable
    of the parent ``Thyra`` package, (except when
    ``${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON`` for course).  For
    an example, see `Explicit enable of a subpackage`_. This means that if a
    `downstream`_ package declares a dependency on the SE package
    ``ThyraEpetra``, but not the parent package ``Thyra``, then the ``Thyra``
    package (and its other subpackages and their dependencies) will not get
    auto-enabled.  This is a key aspect of the TriBITS SE package management
    system.

.. _Support for optional SE package/TPL is enabled by default:

13) **Support for optional SE package/TPL is enabled by default**: For an SE
    package ``<TRIBITS_PACKAGE>`` with an optional dependency on an
    `upstream`_ SE package or TPL ``<TRIBITS_DEP_PACKAGE_OR_TPL>``, TriBITS
    will automatically set the intra-enable variable
    ``<TRIBITS_PACKAGE>_ENABLE_<TRIBITS_DEP_PACKAGE_OR_TPL>=ON`` if
    ``<TRIBITS_PACKAGE>`` and ``<TRIBITS_DEP_PACKAGE_OR_TPL>`` are **both
    enabled**.  For example, if the packages ``Triutils`` and ``EpetraExt``
    are both enabled, then TriBITS will automatically set
    ``EpetraExt_ENABLE_Triutils=ON`` by default and enables support for
    ``Triutils`` in ``EpetraExt``.  Likewise, if ``Teuchos`` and the optional
    TPL ``Boost`` are both enabled, then TriBITS will automatically set
    ``Teuchos_ENABLE_Boost=ON`` by default.  This is obviously the logical
    behavior.  See the full context for these examples in `Explicit enable of
    a package, its tests, an optional TPL, with ST enabled`_.

.. _Support for optional SE package/TPL can be explicitly disabled:

14) **Support for optional SE package/TPL can be explicitly disabled:** Even
    though TriBITS will automatically set
    ``<TRIBITS_PACKAGE>_ENABLE_<TRIBITS_DEP_PACKAGE_OR_TPL>=ON`` by default if
    ``<TRIBITS_PACKAGE>`` and ``<TRIBITS_DEP_PACKAGE_OR_TPL>`` are both
    enabled at the project level (as described above), the user can explicitly
    set ``<TRIBITS_PACKAGE>_ENABLE_<TRIBITS_DEP_PACKAGE_OR_TPL>=OFF`` which
    will turn off optional package-level support for the upstream SE package
    or TPL ``<TRIBITS_DEP_PACKAGE_OR_TPL>`` in the downstream SE package
    ``<TRIBITS_PACKAGE>``.  For example, the user can enable ``EpetraExt`` and
    ``Triutils`` at the project level, but set
    ``EpetraExt_ENABLE_Triutils=OFF`` which will turn off package-level
    support for ``Triutils`` in the ``EpetraExt`` package.  Likewise, the user
    can enable ``Teuchos``, ``Epetra`` and the optional TPL ``Boost``, but set
    ``Epetra_ENABLE_Boost=ON``.  This would provide support for ``Boost`` in
    ``Teuchos`` but not in ``Epetra``.  For examples, see `Explicit disable of
    an optional package dependency`_ and `Explicit disable of an optional TPL
    dependency`_.

.. _Explicit enable of optional SE package/TPL support auto-enables SE package/TPL:

15) **Explicit enable of optional SE package/TPL support auto-enables SE
    package/TPL**: If the user explicitly enables the TriBITS SE package
    ``<TRIBITS_PACKAGE>`` and explicitly sets
    ``<TRIBITS_PACKAGE>_ENABLE_<TRIBITS_DEP_PACKAGE_OR_TPL>=ON`` on input,
    then that will automatically enable the SE package or TPL
    ``<TRIBITS_DEP_PACKAGE_OR_TPL>`` (and all of its upstream dependencies
    accordingly).  For example, if the user sets
    ``Trilinos_ENABLE_EpetraExt=ON`` and ``EpetraExt_ENABLE_Triutils=ON``,
    then that will result in the auto-enable of
    ``Trilinos_ENABLE_Triutils=ON`` regardless of the value of
    ``${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES`` or
    ``${PROJECT_NAME}_SECONDARY_TESTED_CODE``. This true even if the optional
    SE package or TPL is ``EX``.  For example, setting
    ``Thyra_ENABLE_ThyraCrazyStuff=ON`` will result in the enabling of the
    ``EX`` package ``ThyraCrazyStuff``.  However, always remember that
    `Disables trump enables where there is a conflict`_.  For example, see
    `Explicit enable of an optional package dependency`_.

.. _ST SE packages only auto-enabled if ST code is enabled:

16) **ST SE packages only auto-enabled if ST code is enabled**: TriBITS will
    only enable an optional ``ST`` SE package when
    ``${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES=ON`` if
    ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=ON`` is also set.  Otherwise, when
    an optional ``ST`` `upstream`_ dependent SE package is not enabled due to
    ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=OFF``, then a one-line warning is
    printed to stdout.  For example, if the ``EpetraExt`` package is enabled
    but ``ST`` code is not enabled, then the optional SE package ``Triutils``
    will not be enabled (for an example, see `Explicit enable of a package and
    its tests`_).  However, when ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=ON``
    and ``EpetraExt`` is enabled, then ``Triutils`` will be enabled too (for an
    example, see `Explicit enable of a package, its tests, an optional TPL,
    with ST enabled`_).  The TriBITS default is
    ``${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES=ON``.  This helps to avoid
    problems when users try to set a permutation of enables/disables which is
    not regularly tested.

.. _<Project>_ENABLE_ALL_FORWARD_DEP_PACKAGES enables downstream packages/tests:

17) **<Project>_ENABLE_ALL_FORWARD_DEP_PACKAGES enables downstream packages/tests**:
    Setting the user cache-variable
    ``${PROJECT_NAME}_ENABLE_ALL_FORWARD_PACKAGES=ON`` will result in the
    `downstream`_ ``PT`` SE packages and tests to be enabled (and all ``PT``
    and ``ST`` SE packages and tests when
    ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=ON``) for all explicitly enabled
    SE packages.  For example, configuring with ``Trilinos_ENABLE_Epetra=ON``,
    ``Trilinos_ENABLE_TESTS=ON``, and
    ``Trilinos_ENABLE_ALL_FORWARD_PACKAGES=ON`` will result the package
    enables (and test and example enables) for the SE packages ``Triutils``,
    ``EpetraExt``, ``ThyraCoreLibs``, ``ThyraEpetra`` and ``Thyra``.  For an
    example, see `Explicit enable of a package and downstream packages and
    tests`_.

.. _${PROJECT_NAME}_ENABLE_ALL_PACKAGES:

.. _<Project>_ENABLE_ALL_PACKAGES enables all PT (cond. ST) SE packages:

18) **<Project>_ENABLE_ALL_PACKAGES enables all PT (cond. ST) SE packages**:
    Setting the user cache-variable ``${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON``
    will result in the enable of all ``PT`` SE packages when
    ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=OFF`` and all ``PT`` and ``ST`` SE
    packages when ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=ON``.  For an
    example, see `Enable all packages`_.  When the project is a meta-project,
    this will only enable the project's primary meta-project packages (PMPP).
    That is, a package will only be enabled due to
    ``${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON`` when its parent repository does
    not have `${REPOSITORY_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES`_ set to
    ``TRUE``.  However, if::

       ${REPOSITORY_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES=TRUE

    then the SE package may be enabled if it (or its parent package) is listed
    in `${REPOSITORY_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES_EXCEPT`_.

.. _<Project>_ENABLE_TESTS only enables explicitly enabled SE package tests:

19) **<Project>_ENABLE_TESTS only enables explicitly enabled SE package
    tests**: Setting ``${PROJECT_NAME}_ENABLE_TESTS=ON`` will **only enable
    tests for explicitly enabled SE packages** on input.  For example,
    configuring with ``Trilinos_ENABLE_RTOp=ON`` and
    ``Trilinos_ENABLE_TESTS=ON`` will only result in the enable of tests for
    ``RTOp``, not ``Teuchos`` (even through TriBITS will enable ``Teuchos``
    because it is a required dependency of ``RTOp``).  See an example, see
    `Explicit enable of a package and its tests`_.  When the project is a
    meta-project, this will only enable tests for the project's primary
    meta-project packages (PMPP).  The uses the same logic as for
    ``${PROJECT_NAME}_ENABLE_ALL_PACKAGES``, see
    `<Project>_ENABLE_ALL_PACKAGES enables all PT (cond. ST) SE packages`_.

.. _If no SE packages are enabled, nothing will get built:

20) **If no SE packages are enabled, nothing will get built**: Most TriBITS
    projects are set up such that if the user does not explicitly enable at
    least one SE package in some way, then nothing will be enabled or built.
    In this case, when ``${PROJECT_NAME}_ALLOW_NO_PACKAGES=TRUE`` a warning
    will be printed and configuration will complete.  However, if
    ``${PROJECT_NAME}_ALLOW_NO_PACKAGES=FALSE``, then the configure will die
    with an error message.  For example, the ``checkin-test.py`` tool sets
    ``${PROJECT_NAME}_ALLOW_NO_PACKAGES=OFF`` to make sure that something gets
    enabled and tested in order to accept the results of the test and allow a
    push.  For an example, see `Default configure with no packages enabled on
    input`_.

.. _TriBITS prints all enables and disables to stdout:

21) **TriBITS prints all enables and disables to stdout**: TriBITS prints out
    (to ``cmake`` stdout) the initial set of enables/disables on input, prints
    a line whenever it sets (or overrides) an enable or disable, and prints
    out the final set of enables/disables.  Therefore, the user just needs to
    grep the ``cmake`` stdout to find out why any particular SE package or TPL
    is enabled or disabled in the end.  In addition, will print out when
    tests/examples for a given SE package gets enabled and when support for
    optional SE packages and TPLs is enabled or not.  Examples of this output
    is given in all of the below examples but a detailed description of this
    output is given in `Explicit enable of a package and its tests`_.

.. _TriBITS auto-enables/disables done using non-cache local variables:

22) **TriBITS auto-enables/disables done using non-cache local variables**:
    TriBITS setting (or overrides) of enable/disable cache variables are done
    by setting local non-cache variables at the top project-level scope.  This
    is done so they don't get set in the cache and so that the same dependency
    enable/disable logic is redone, from scratch, with each re-configure.
    This results in the same enable/disable logic output as for the initial
    configure.  This is to avoid confusion by the user about why some SE
    packages and TPLs are enabled and some are not on subsequent reconfigures.

TriBITS prints out a lot of information about the enable/disable logic as it
applies the above rules/behaviors.  For a large TriBITS project with lots of
packages, this can produce a lot of output to stdout.  One just needs to
understand what TriBITS is printing out and where to look in the output for
different information.  The examples in the section `Example Enable/Disable
Use Cases`_ show what this output looks like for the various enable/disable
scenarios and tries to explain in more detail the reasons for why the given
behavior is implemented the way that it is.  Given this output, the rule
definitions given above, and the detailed `Example Enable/Disable Use Cases`_,
one should always be able to figure out exactly why the final set of
enables/disables is the way it is, even in the largest and most complex of
TriBITS projects.  (NOTE: The same can *not* be said for many other large
software configuration and deployment systems where basic decisions about what
to enable and disable are hidden from the user and can be very difficult to
track down and debug.)

The above behaviors address the majority of the functionality of the TriBITS
dependency management system.  However, when dealing with TriBITS projects
with multiple repositories, some other behaviors are supported through the
definition of a few more variables.  The following TriBITS repository-related
variables alter what packages in a given TriBITS repository get enabled
implicitly or not by TriBITS:

  ``${REPOSITORY_NAME}_NO_IMPLICIT_PACKAGE_ENABLE``

    If set to ``ON``, then the packages in Repository ``${REPOSITORY_NAME}``
    will not be implicitly enabled in any of the package adjustment logic.

  ``${REPOSITORY_NAME}_NO_IMPLICIT_PACKAGE_ENABLE_EXCEPT``

    List of packages in the Repository ``${REPOSITORY_NAME}`` that will be
    allowed to be implicitly enabled.  Only checked if
    ``${REPOSITORY_NAME}_NO_IMPLICIT_PACKAGE_ENABLE`` is true.

The above variables typically are defined in the outer TriBITS Project's CTest
driver scripts or even in top-level project files in order to adjust how
packages in its listed repositories are handled.  What these variable do is to
allow a large project to turn off the auto-enable of optional SE packages in a
given TriBITS repository to provide more detailed control of what gets used
from a given TriBITS repository.  This, for example, is used in the CASL VERA
project to manage some of its extra repositories and packages to further reduce
the number of packages that get auto-enabled.

.. ToDo: We should likely change this so that only tests and example don't get
.. enabled by default, not the packages itself.  This is mainly used to avoid
.. certain packages from getting enabled in the CI server when updates are
.. received.  This needs to be explained in the section on the CI server.


Example Enable/Disable Use Cases
---------------------------------

Below, a few of the standard enable/disable use cases for a TriBITS project
are given using the `Example ReducedMockTrilinos Project Dependency
Structure`_ that demonstrate the `TriBITS Dependency Handling Behaviors`_.

The use cases covered are:

* `Default configure with no packages enabled on input`_
* `Explicit enable of a package and its tests`_
* `Explicit enable of a package, its tests, an optional TPL, with ST enabled`_
* `Explicit disable of a package`_
* `Conflicting explicit enable and disable`_
* `Explicit enable of an optional TPL`_
* `Explicit disable of an optional TPL`_
* `Explicit disable of a required TPL`_
* `Explicit enable of a subpackage`_
* `Explicit enable of an optional package dependency`_
* `Explicit disable of an optional package dependency`_
* `Explicit enable of an optional TPL dependency`_
* `Explicit disable of an optional TPL dependency`_
* `Explicit enable of a package and downstream packages and tests`_
* `Enable all packages`_

All of these use cases and more can be easily run from the command-line by
first setting:

  $ export REDUCED_MOCK_TRILINOS=<base-dir>/tribits/examples/ReducedMockTrilinos

and then copy and pasting the ``cmake`` commands shown below.  Just make sure
to run these in a temp directory because this actually configures a CMake
project in the local directory.  Just make sure and run::

  $ rm -r CMake*

before each run to clear the CMake cache.

These use cases are now described in detail below.

.. _Default configure with no packages enabled on input:

**Default configure with no packages enabled on input**

The first use-case to consider is the configure of a TriBITS project without
enabling any packages.  For the ``ReducedMockTrilinos`` project, this is done
with::

   $ cmake ${REDUCED_MOCK_TRILINOS}

which produces the relevant dependency-related output:

.. include:: ReducedMockTrilinosOutput/NoEnables.txt
   :literal:

The above example demonstrates the following behaviors of the TriBITS
dependency handling system:

* `PT/ST SE packages given default unset enable/disable state`_
* `EX SE packages disabled by default`_ (i.e., the ``EX`` SE package
  ``ThyraCrazyStuff`` is set to ``OFF`` by default at the very beginning).
* `PT/ST TPLs given default unset enable/disable state`_
* `EX TPLs given default unset enable/disable state`_
* `If no SE packages are enabled, nothing will get built`_

.. _Explicit enable of a package and its tests:

**Explicit enable of a package and its tests**

One of the most typical use cases is for the user to explicitly enable one or
more top-level TriBITS package and enable its tests.  This configuration would
be used to drive local development on a specific set of packages (i.e. tests
do not need to be enabled for packages not being changed).

Consider the configure of the ``ReducedMockTrilinos`` project enabling the
top-level ``Thyra`` package and its tests with::

   $ cmake -DTrilinos_ENABLE_Thyra:BOOL=ON \
      -DTrilinos_ENABLE_TESTS:BOOL=ON \
      ${REDUCED_MOCK_TRILINOS}

which produces the relevant dependency-related output:

.. include:: ReducedMockTrilinosOutput/EnableThyra_EnableTests.txt
   :literal:

This is a configuration that a developer would use to develop on the ``Thyra``
package and its subpackages for example.  There is no need to be enabling the
tests and examples for upstream packages unless those packages are going to be
changed as well.

This case demonstrates a number of TriBITS dependency handling behaviors that
are worth some discussion.

First, note that enabling the parent package ``Thyra`` with
``Trilinos_ENABLE_Thyra=ON`` right away results in the auto-enable of its
``PT`` subpackages ``ThyraCoreLibs`` and ``ThyraEpetra`` which demonstrates
the behavior `Enable/disable of parent package is enable/disable for
subpackages`_.  Note that the ``ST`` subpackages ``ThyraGoodStuff`` and
``ThyraEpetraExt`` where *not* enabled because
``${PROJECT_NAME}_SECONDARY_TESTED_CODE=OFF`` (which is off by default) which
demonstrates the behavior `ST SE packages only auto-enabled if ST code is
enabled`_.

Second, note the auto-enable of required upstream SE packages ``Epetra``,
``RTOp`` and ``Teuchos`` shown in lines like::

  -- Setting Trilinos_ENABLE_Teuchos=ON because ThyraCoreLibs has a required dependence on Teuchos

Lastly, note that the final set of enabled packages, SE packages,
tests/examples and TPLs can be clearly seen when processing the TPLs and
top-level packages in the lines::

  Getting information for all enabled TPLs ...
  
  -- Processing enabled TPL: BLAS
  -- Processing enabled TPL: LAPACK
  
  Configuring individual enabled Trilinos packages ...
  
  Processing enabled package: Teuchos (Libs)
  Processing enabled package: RTOp (Libs)
  Processing enabled package: Epetra (Libs)
  Processing enabled package: Thyra (CoreLibs, Epetra, Tests, Examples)

Note that subpackage enables are listed with their parent packages along with
if the tests and/or examples are enabled.  Top-level level packages that don't
have subpackages just show if the ``Libs`` or also the ``Tests`` and
``Examples`` have been enabled as well.

.. _Explicit enable of a package, its tests, an optional TPL, with ST enabled:

**Explicit enable of a package, its tests, an optional TPL, with ST enabled**

An extended use case shown here is for the explicit enable of a package and its
tests along with the enable of an optional TPL and with ``ST`` code enabled.
This is a configuration that would be used to support the local development of
a TriBITS package that involves modifying ``ST`` software.

Consider the configure of the ``ReducedMockTrilinos`` project with::

   $ cmake -DTPL_ENABLE_Boost:BOOL=ON \
      -DTrilinos_ENABLE_Thyra:BOOL=ON \
      -DTrilinos_ENABLE_TESTS:BOOL=ON \
      -DTrilinos_ENABLE_SECONDARY_TESTED_CODE:BOOL=ON \
      -DTrilinos_ENABLE_TESTS:BOOL=ON \
      ${REDUCED_MOCK_TRILINOS}

which produces the relevant dependency-related output:

.. include:: ReducedMockTrilinosOutput/EnableThyra_EnableTests_EnableBoost_ST.txt
   :literal:

A few more behaviors of the TriBITS system that this particular configuration
use-case shows are described below.

First, note the enable of the ``ST`` ``Thyra`` subpackages in lines like::

  -- Setting subpackage enable Trilinos_ENABLE_ThyraGoodStuff=ON because parent package \
   Trilinos_ENABLE_Thyra=ON

Second, note the auto-enable of support for optional SE packages in lines
like::

  -- Setting EpetraExt_ENABLE_Triutils=ON since Trilinos_ENABLE_EpetraExt=ON \
   AND Trilinos_ENABLE_Triutils=ON

Third, note the auto-enable of support for the optional TPL ``Boost`` in the
line::

  -- Setting Teuchos_ENABLE_Boost=ON since TPL_ENABLE_Boost=ON

.. _Explicit disable of a package:

**Explicit disable of a package**

Another common use case is to enable a package but to disable an optional
upstream package.  This type of configuration would be used as part of a
"black list" approach to enabling only a subset of packages and optional
support.  The "black list" approach is to enable a package with
``${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES=ON`` (the TriBITS default) but
then to turn off a specific set of packages that you don't want.  This is
contrasted with a "white list" approach where you would configure with
``${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES=OFF`` and then have to manually
enable all of the optional packages you want.  Experience with projects like
Trilinos show that the "black list" approach is generally to be preferred for
a few reasons.

Consider the configure of the ``ReducedMockTrilinos`` project enabling ``Thyra`` but disabling ``Epetra`` with::

   $ cmake  -DTrilinos_ENABLE_Thyra:BOOL=ON \
      -DTrilinos_ENABLE_Epetra:BOOL=OFF \
      -DTrilinos_ENABLE_TESTS:BOOL=ON \
      ${REDUCED_MOCK_TRILINOS}

which produces the relevant dependency-related output:

.. include:: ReducedMockTrilinosOutput/EnableThyra_DisableEpetra_EnableTests.txt
   :literal:

Note how the disable of ``Epetra`` wipes out all of the required and optional
SE packages and intra-package dependencies that depend on ``Epetra``.  What is
left is only the ``ThyraCoreLibs`` and its upstream dependencies that don't
depend on ``Epetra`` (which is only ``RTOp`` and ``Teuchos``).

.. _Conflicting explicit enable and disable:

**Conflicting explicit enable and disable**

One use case that occasionally comes up is when a set of inconsistent enables
and disables are set.  While this seems illogical that anyone would ever do
this, when it comes to larger more complex projects with lots of packages and
lots of dependencies, this can happen very easily.  In some cases, someone is
enabling a set of packages they want and is trying to weed out as many of the
(what they think) are optional dependencies that they don't need and
accidentally disables a package that is an indirect required dependency of one
of the packages they want (in which case the configure should likely fail and
provide a good error message).  The other use case where conflicting
enables/disables can occur is in CTest drivers using `TRIBITS_CTEST_DRIVER()`_
where an upstream package has failed and is explicitly disabled (in which case
it should gracefully disable downstream dependent packages).  TriBITS can
either be set up to have the disable override the explicit enable or stop the
configure in error depending on the value of the cache variable
`${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES`_ (see `Disables trump
enables where there is a conflict`_).

For example, consider what happens with the ``ReducedMockTrilinos`` project if
someone tries to enable the ``RTOp`` package and disable the ``Teuchos``
package.  This is not consistent because ``RTOp`` has a required dependency on
``Teuchos``.  The default behavior of TriBITS is this case is shown in the
below configure::

   $ cmake -DTrilinos_ENABLE_Epetra:BOOL=ON \
      -DTrilinos_ENABLE_RTOp:BOOL=ON \
      -DTrilinos_ENABLE_Teuchos:BOOL=OFF \
      ${REDUCED_MOCK_TRILINOS}

which produces the relevant dependency-related output:

.. include:: ReducedMockTrilinosOutput/EnableEpetra_EnableRTOp_DisableTeuchos.txt
   :literal:

As shown above, the TriBITS default (which is
``${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=OFF``) results in a
configure-time error with a good error message.

However, if one sets
``${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON`` and configures
with::

   $ cmake -DTrilinos_ENABLE_Epetra:BOOL=ON \
      -DTrilinos_ENABLE_RTOp:BOOL=ON \
      -DTrilinos_ENABLE_Teuchos:BOOL=OFF \
      -DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES:BOOL=ON \
      ${REDUCED_MOCK_TRILINOS}

then the disable trumps the enable and results in a successful configure as
shown in the following relevant dependency-related output:

.. include:: ReducedMockTrilinosOutput/EnableEpetra_EnableRTOp_DisableTeuchos_DisableEnabledFwdDepPackages.txt
   :literal:

As shown above, what you end up with is just the enabled package ``Epetra``
which does not have a required dependency on the disabled package ``Teuchos``.
Developers of large complex TriBITS projects would be wise to set the default
for `${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES`_ to ``ON``,
especially in automated builds and testing.

.. _Explicit enable of an optional TPL:

**Explicit enable of an optional TPL**:

ToDo: Set Trilinos_ENABLE_Thyra=ON and TPL_ENABLE_MPI=ON

.. _Explicit disable of an optional TPL:

**Explicit disable of an optional TPL**:

ToDo: Set Trilinos_ENABLE_Thyra=ON and TPL_ENABLE_MPI=OFF

.. _Explicit disable of a required TPL:

**Explicit disable of a required TPL**

ToDo: Set Trilinos_ENABLE_Epetra=ON and Trilinos_ENABLE_BLAS=OFF

.. _Explicit enable of a subpackage:

**Explicit enable of a subpackage**

ToDo: Enable ThyraEpetra and show how it enables other SE packages and at the
end, enables the Thyra package (just for show).

.. _Explicit enable of an optional package dependency:

**Explicit enable of an optional package dependency**

ToDo: Set Trilinos_ENABLE_EpetraExt=ON and EpetraExt_ENABLE_Triutils=ON and
shows how it enables Trilinos_ENABLE_Triutils=ON even through ST code is not
enabled.

.. _Explicit disable of an optional package dependency:

**Explicit disable of an optional package dependency**

ToDo: Set Trilinos_ENABLE_EpetraExt=ON, Trilinos_ENABLE_Triutils=ON, and
EpetraExt_ENABLE_Triutils=OFF.  Discuss how EpetraExt's and ThyraEpetraExt's
CMakeLists.txt files might turn off some features if they detects that
EpetraExt/Triutils support is turned off.

.. _Explicit enable of an optional TPL dependency:

**Explicit enable of an optional TPL dependency**

ToDo: The current ReducedMockTrilinos is not set up to give a good example of
this.  We should add an optional Boost dependency to say, Epetra.  Then we
could show the enable of Teuchos and Epetra and Epetra_ENABLE_Boost=ON.  That
would enable Boost and enable support for Boost in Epetra but would not
provide support for Boost in Teuchos.

.. _Explicit disable of an optional TPL dependency:

**Explicit disable of an optional TPL dependency**

ToDo: The current ReducedMockTrilinos is not set up to give a good example of
this.  We should add an optional Boost dependency to say, Epetra.  Then we
could show the enable of Teuchos and Epetra and TPL_ENABLE_Boost=ON but set
Epetra_ENABLE_Boost=OFF.  That would provide support for Boost in Teuchos but
not in Epetra.

.. _Explicit enable of a package and downstream packages and tests:

**Explicit enable of a package and downstream packages and tests**

ToDo: Set Trilinos_ENABLE_RTOp=ON,
Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON, and Trilinos_ENABLE_TESTS=ON and
show what packages and tests/examples get enabled.  This is the use case for
the checkin-test.py tool for PT enabled code.

.. _Enable all packages:

**Enable all packages**

The last use case to consider is enabling all defined packages.  This
configuration would be used for either doing a full test of all of the
packages defined or to create a distribution of the project.

Enabling all ``PT`` packages the with::

   $ cmake -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON \
      -DTrilinos_DUMP_PACKAGE_DEPENDENCIES:BOOL=ON \
      ${REDUCED_MOCK_TRILINOS}

produces the relevant dependency-related output:

.. include:: ReducedMockTrilinosOutput/EnableAllPackages.txt
   :literal:

As shown above, only the ``PT`` SE packages get enabled.  To also enable the
``ST`` packages as well, one additionally set
``${PROJECT_NAME}_SECONDARY_TESTED_CODE=ON`` at configure time.

.. _<Project>PackageDependencies.xml:


TriBITS Project Dependencies XML file and tools
-----------------------------------------------

The TriBITS CMake configure system can write out the project's package
dependencies into a file ``<Project>Dependencies.xml`` (or any name one wants
to give it).  This file is used by a number of the TriBITS SE-related tools.
The structure of this XML file, showing one of the more interesting mock
packages from the `MockTrilinos`_ project is shown below::

  <PackageDependencies project="Trilinos">
    ...
    <Package name="Amesos" dir="packages/amesos" type="PT">
      <LIB_REQUIRED_DEP_PACKAGES value="Teuchos,Epetra"/>
      <LIB_OPTIONAL_DEP_PACKAGES value="EpetraExt"/>
      <TEST_REQUIRED_DEP_PACKAGES/>
      <TEST_OPTIONAL_DEP_PACKAGES value="Triutils,Galeri"/>
      <LIB_REQUIRED_DEP_TPLS/>
      <LIB_OPTIONAL_DEP_TPLS value="SuperLUDist,ParMETIS,UMFPACK,SuperLU,MUMPS"/>
      <TEST_REQUIRED_DEP_TPLS/>
      <TEST_OPTIONAL_DEP_TPLS/>
      <EmailAddresses>
        <Regression address="amesos-regression@repo.site.gov"/>
      </EmailAddresses>
      <ParentPackage value=""/>
    </Package>
    ...
  </PackageDependencies>

This XML file contains the names, directories, `Test Test Category`_
(i.e. ``type``), CDash email address, and all of the SE package and TPL
dependencies for every SE package in the TriBITS project (including add-on
repositories if specified).  There are several python tools under
``tribits/ci_support/`` that read in this file and use the created
data-structure for various tasks.  This file and these tools are used by
`checkin-test.py`_ and `TRIBITS_CTEST_DRIVER()`_.  But these tools can also be
used to construct other workflows and tools.

.. _TribitsDumpDepsXmlScript.cmake:

A TriBITS project configure can create this file as a byproduct of
configuration by setting the configure option (see `Outputting package
dependency information`_), or the CMake -P script
**TribitsDumpDepsXmlScript.cmake** can be used to create this file on the fly
without having to configure a TriBITS project.  To create this file outside of
configuration, one can run::

  cmake \
    [-D PROJECT_SOURCE_DIR=<projectSourceDir>] \
    [-D <Project>_PRE_REPOSITORIES=<prepo0>,<prepo1>,...] \
    [-D <Project>_EXTRA_REPOSITORIES=<erepo0>,<erepo1>,...] \
    -D <Project>_DEPS_XML_OUTPUT_FILE=<projectDepsFileOut> \
    -P <tribitsDir>/ci_support/TribitsDumpDepsXmlScript.cmake

If TriBITS is snashotted into the project in the standard location
``<projectDir>/cmake/tribits`` or the entire TriBITS repo is cloned under
``<projectDir>/TriBITS`` (so that the ``tribits`` dir is
``<projectDir>/TriBITS/tribits``) then one can leave off
``-DPROJECT_SOURCE_DIR=<projectSourceDir>`` and (if not wanting to include
extra repos) just run::

  cmake \
    -D <Project>_DEPS_XML_OUTPUT_FILE=<projectDepsFileOut> \
    -P <projectSourceDir>/cmake/tribits/ci_support/TribitsDumpDepsXmlScript.cmake

Once the XML file ``<projectDepsFileOut>`` is created, it can be used in
various types of analysis and used with different tools and commands.

The tool `get-tribits-packages-from-files-list.py`_ can be used to determine
the list of TriBITS SE packages that need to be tested given a list of changed
files (e.g. as returned from ``git diff --name-only <from>..<to> >
changed-files.txt``).  This is used in the `checkin-test.py`_ tool and the
`TRIBITS_CTEST_DRIVER()`_ function to determine what TriBITS packages need to
be tested based on what files have been changed.

The tool `get-tribits-packages-from-last-tests-failed.py`_ can be used to
extract the list of TriBITS SE packages that correspond to the failings tests
listed in the CTest-generated
``<build-dir>/Testing/Temporary/LastTestsFailed*.log`` file.  This tool is
used in the `TRIBITS_CTEST_DRIVER()`_ function in CI-testing mode to determine
what packages must be re-tested if they failed in the last CI iteration.

The tool `filter-packages-list.py`_ takes in a list of TriBITS SE package
names and then filters the list according the `Test Test Category`_ of the
packages.  This is used in testing workflows that only test a subset of
packages according to the Test Test Category at different stages in the
workflow.  For example, the `checkin-test.py`_ tool and the
`TRIBITS_CTEST_DRIVER()`_ function use this filtering to only test Primary
Tested (PT) or Secondary Tested (ST) packages for a given set of changed files
in a continuous integration workflow (see `Nested Layers of TriBITS Project
Testing`_).


TriBITS Automated Testing
=========================

Much of the value provided by the TriBITS system is support for testing of
complex projects.  Many different types of testing are required in a
complex project and development effort.  A large project with lots of
repositories and packages provides a number of testing and development
challenges but also provides a number of opportunities to do testing in an
efficient way; especially pre-push and post-push continuous integration (CI)
testing.  In addition, a number of post-push automated nightly test cases must
be managed.  TriBITS takes full advantage of the features of raw CMake, CTest,
and CDash in support of testing and where gaps exist, TriBITS provides tools
and customizations.

The following subsections describe several aspects to the TriBITS support for
testing.  The subsection `Test Classifications for Repositories, Packages, and
Tests`_ defines the different types of test-related classifications that are
defined by TriBITS.  These different test classifications are then used to
define a number of different standard `Nested Layers of TriBITS Project
Testing`_ which include different types of CI testing as well as nightly and
other tests.  One of the most important types of CI testing, pre-push testing,
is then described in more detail in the subsection `Pre-push Testing using
checkin-test.py`_.  The subsection `TriBITS CTest/CDash Driver`_ describes the usage
of the advanced `TRIBITS_CTEST_DRIVER()`_ function to do incremental project
testing of a projects using advanced ``ctest -S`` scripts.  The final
subsection `TriBITS CDash Customizations`_ describes how projects can use a
CDash server to more effectively display test results and provide
notifications for failures that are compartmentalized for a large project.


Test Classifications for Repositories, Packages, and Tests
----------------------------------------------------------

TriBITS defines a few different testing-related classifications for a TriBITS
project.  These different classifications are used to select subsets of the
project's repositories, packages (and code within these packages), and tests
to be included in a given project build and test definition.  These different
classification are:

* `Repository Test Classification`_
* `SE Package Test Group`_
* `Test Test Category`_

These different test-related classifications are used to defined several
different `Nested Layers of TriBITS Project Testing`_.  First, the `Repository
Test Classification`_ determines what repositories are even processed in order
for their SE packages to even consider being enabled.  Second, if a repository
is selected, then the `SE Package Test Group`_ determines what SE packages
(and optional code in those packages) are even enabled such that their
`<packageDir>/CMakeLists.txt`_ files are even processed (i.e. according to
`TriBITS Dependency Handling Behaviors`_).  Lastly, if an SE package gets
enabled, then the `Test Test Category`_ determines what test executables and
test cases get defined using the functions `TRIBITS_ADD_EXECUTABLE()`_,
`TRIBITS_ADD_TEST()`_ and `TRIBITS_ADD_ADVANCED_TEST()`_.

More detailed descriptions of `Repository Test Classifications`_, `SE Package
Test Groups`_, and `Test Test Categories`_ are given in the following
subsections.

.. _Repository Test Classification:
.. _Repository Test Classifications:

**Repository Test Classification**

The first type of test-related classification is for extra repositories
defined in the file `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ (pulled
in trough the `${PROJECT_NAME}_EXTRAREPOS_FILE`_ cache variable) using the
``REPO_CLASSIFICATION`` field in the macro call
`TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES()`_.  These classifications map to
the standard CTest dashboard types ``Continuous``, ``Nightly``, and
``Experimental`` (see `CTest documentation`_ and `TriBITS CTest/CDash Driver`_
for details).

.. _Repository Test Continuous:

* Repositories marked **Continuous** match the standard CTest dashboard type
  ``Continuous``.  These repositories are pulled in when
  ``${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Continuous``, or
  ``Nightly``.  Repositories marked as ``Continuous`` are cloned, updated, and
  processed by default in all project automated testing described in `Nested
  Layers of TriBITS Project Testing`_.  NOTE: One should **not** confuse this
  with the `Test Test Category CONTINUOUS`_.

.. _Repository Test Nightly:

* Repositories marked **Nightly** match the standard CTest dashboard type
  ``Nightly``.  These repositories are pulled in by default when
  ``${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Nightly``.  Repositories
  marked as ``Nightly`` are not processed by default as part of either
  `Pre-Push CI Testing`_ or `Post-Push CI Testing`_.  One would mark a
  repository as ``Nightly`` for few reasons.  First, an extra repo may be
  marked as ``Nightly`` if it may not be available to all developers to clone
  and only the nightly testing processes and machines may have access.  Also,
  an extra repo may also be marked as ``Nightly`` if it does not contain any
  packages that the project wants to pay the cost to include in even
  `Post-Push CI Testing`_.  NOTE: One should **not** confuse this with the
  `Test Test Category NIGHTLY`_.

.. _Repository Test Experimental:

* Repositories marked **Experimental** match the standard CTest dashboard type
  ``Experimental``.  Repositories marked as ``Experimental`` are not processed
  by default as part of any automated testing described in the subsection
  `Nested Layers of TriBITS Project Testing`_ (except for perhaps some
  experimental builds).  The main reason that an extra repo may be marked as
  ``Experimental`` is that it may only contain ``EX`` SE packages and
  therefore none of these packages would be enabled by default anyway.  Also,
  a repo may be marked as ``Experimental`` if it is developed in a very sloppy
  way such that one cannot even assume that the repository's
  `<repoDir>/PackagesList.cmake`_, `<repoDir>/TPlsList.cmake`_, and
  `<packageDir>/cmake/Dependencies.cmake`_ files are not without errors.
  Since these files are always processed if the repository is included, then
  any errors in these files will cause the entire configure to fail!

.. ToDo: Determine, add unit tests for, and document exactly how
.. 'Experimental' repos are handled
.. w.r.t. ${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE.  Currently it
.. looks like there are not processed if
.. ${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Experimental.  I am not
.. sure why we bother with Experimental repos other than to just ignore them.
.. However, it does seem they get cloned no matter what their classification
.. in the TribitsCTestDriverCore.cmake script.

.. _SE Package Test Group:
.. _SE Package Test Groups:

**SE Package Test Group**

Once a set of TriBITS repositories are selected in accordance with their
`Repository Test Classification`_, that determines the set of SE packages and
TPLs defined for the TriBITS project.  Given the set of defined SE packages
and TPLs, the set of SE packages that get enabled is determined by the **SE
Package Test Group** which is defined and described here.

Every `TriBITS SE Package`_ and `TriBITS TPL`_ is assigned a test group. These
test groups are for *Primary Tested* (`PT`_) code, *Secondary Tested* (`ST`_)
code, and *Experimental* (`EX`_) code.  The test group defines *what* SE
packages and TPL get selected (or are excluded from being selected) to include
in a given build for testing-related purposes.  SE packages may also
conditionally build in additional code based on the testing group.  The
detailed rules for when an SE package or TPL are selected or excluded from the
build based on the test group is given in `TriBITS Dependency Handling
Behaviors`_.  We only summarize those rules here.

More detailed descriptions of the test groups are given below.

.. _PT:

* **Primary Tested (PT)** Code is of the highest priority to keep working for
  the current development effort.  SE packages and TPLs may be selected to be
  ``PT`` for a number of reasons.  First, if the capability provided by the
  code is mature and if a regression would cause major harm to a customer, the
  code should likely be marked as ``PT``.  Also, if the build and correct
  functioning of the code is needed by other development team members to
  support their day-to-day development activities, then the code should be
  marked as ``PT`` as well.  A TPL, on the other hand, is marked as ``PT`` if
  it is required by a ``PT`` SE package.  Every project developer is expected
  to have every ``PT`` TPL installed on every machine where they do
  development on and from which they push to the global repo (see
  `checkin-test.py`_ tool).  ``PT`` SE packages and TPLs are the foundation
  for `Pre-Push CI Testing`_.

.. _ST:

* **Secondary Tested (ST)** Code is still very important code for the project
  and represents important capability to maintain but is excluded from the
  ``PT`` set of code for one of a few different reasons.  First, code may be
  marked as ``ST`` if is not critical to drive most day-to-day development
  activities.  If ``ST`` code breaks, it usually will not cause immediate and
  major harm to most developers.  Also, code may be marked as ``ST`` if it has
  required dependencies on ``ST`` TPLs which are either hard to install or may
  not be available on all platforms where developers do their development and
  from where they push changes to the global repo.  In addition, code may be
  marked as ``ST`` if the project is just too big and developers can't be
  expected to build and test all of this code with every push (so a decision
  is made to only make some code as ``PT`` so that pushes don't take too
  long).  ``ST`` code can be included in the TriBITS auto-enable algorithms by
  setting the variable ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=ON`` (see
  `TriBITS Dependency Handling Behaviors`_).  Otherwise, ``ST`` code **is
  not** enabled by auto-enable algorithms.  Typically, ``ST`` code is excluded
  from the default builds in `Pre-Push CI Testing`_ but ``ST`` code **is**
  typically tested in `Post-Push CI Testing`_, `Nightly Testing`_ as well as
  in other testing cases.

.. _EX:

* **Experimental (EX)** Code is usually too unstable, buggy, or non-portable
  to be maintained as part of the project's automated testing processes.  Or
  the code may not be important enough to the project to bother paying the
  cost to test it in the project's automated testing processes.  The ability
  to mark some code as ``EX`` allows the developers of that code to include
  their code in the project's VC repos along with the rest of the project's
  code and be able to take advantage of all of the development tools provided
  by the project (using TriBITS) but at the same time not having to "sign up"
  for all of the responsibilities of maintaining working software that every
  developer has to take a part in helping to keep working.

The test group for each type of entity is assigned in the following places:

* The top-level `TriBITS Package`_'s test group is assigned using the
  ``CLASSIFICATION`` field in the macro call
  `TRIBITS_REPOSITORY_DEFINE_PACKAGES()`_ in its parent repository's
  `<repoDir>/PackagesList.cmake`_ file.

* A `TriBITS Subpackage`_'s test group is assigned using the
  ``CLASSIFICATIONS`` field of the `SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS`_
  argument in the macro call `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_ in its
  parent package's `<packageDir>/cmake/Dependencies.cmake`_ file.

* A `TriBITS TPL`_'s test group is assigned using the ``CLASSIFICATION`` field
  in the macro call `TRIBITS_REPOSITORY_DEFINE_TPLS()`_ in its parent
  repository's `<repoDir>/TPLsList.cmake`_ file.

After these files are processed, the variable `${PACKAGE_NAME}_TESTGROUP`_
gives the test group for each defined SE Package while the variable
``${TPL_NAME}_TESTGROUP`` gives the test group for each defined TPL.

Note that the test group classification ``PT``/``ST``/``EX`` is *not* to be
confused with the *maturity level* of the SE package as discussed in the
`TriBITS Lifecycle Model`_.  The test group classification in no way implies
the maturity of the given TriBITS SE Package or piece of code.  Instead, the
test group is just used to sub-select packages (and pieces of code within
those packages) that are the most important to sustain for the various current
development group's activities.  While more-mature code would typically never
be classified as ``EX`` (Experimental), there are cases were immature packages
may be classified as ``ST`` or even ``PT``.  For example, a very important
research project may be driving the development of a very new algorithm with
the low *maturity level* of *Research Stable* (``RS``) or even *Exploratory*
(``EP``) because keeping that code working may be critical to keeping the
research project on track.

In addition to just selecting ``PT`` and ``ST`` SE packages as a whole, a
TriBITS ``PT`` SE package can also contain conditional code and test
directories that get enabled when ``${PROJECT_NAME}_SECONDARY_TESTED_CODE=ON``
and therefore represents more ``ST`` code.  The package's
`<packageDir>/CMakeLists.txt`_ files can contain simple if statements and can
use the `TRIBITS_SET_ST_FOR_DEV_MODE()`_ function to automatically select
extra code to enable when ``ST`` is enabled or when the project is in release
mode.

.. ToDo: Provide examples in TribitsExampleProject of ST and EX SE packages
.. and of conditional ST and EX code in an PT SE package.  This is important
.. to round out the examples.

.. _Test Test Category:
.. _Test Test Categories:

**Test Test Category**

Once a package is even defined (due to its parent repository's selection
consistent with its `Repository Test Classification`_) and is the package is
enabled (consistent with its `SE Package Test Group`_) then the set of
individual test executables and test cases that are included or not in that
package depends on the ``CATEGORIES`` argument in the functions
`TRIBITS_ADD_EXECUTABLE()`_, `TRIBITS_ADD_TEST()`_ and
`TRIBITS_ADD_ADVANCED_TEST()`_, and the `${PROJECT_NAME}_TEST_CATEGORIES`_
variable.  This **Test Test Category** defines the last "knob" that the
development team has in controlling what tests get run in a particular test
scenario as described in the section `Nested Layers of TriBITS Project
Testing`_.

The currently allowed values for the *Test Test Category* are ``BASIC``,
``CONTINUOUS``, ``NIGHTLY``, ``HEAVY``, and ``PERFORMANCE``.  Tests are
enabled based on their assigned test test category matching the categories set
in the CMake cache variable `${PROJECT_NAME}_TEST_CATEGORIES`_.  The test test
categories ``BASIC``, ``CONTINUOUS``, ``NIGHTLY``, and ``HEAVY`` are subsets
of each other.  That is, a ``BASIC`` test is automatically included in the set
of ``CONTINUOUS``, ``NIGHTLY``, and ``HEAVY`` tests (as set using
``${PROJECT_NAME}_TEST_CATEGORIES``).

The different test test categories are described below in more detail:

.. _Test Test Category BASIC:

* Tests marked **BASIC** represent key functionality that is needed by nearly
  every developer that works on the project and so must be protected at all
  times and are therefore included in `Pre-Push CI Testing`_.  Tests marked as
  ``BASIC`` are enabled for the values of ``${PROJECT_NAME}_TEST_CATEGORIES``
  of ``BASIC``, ``CONTINUOUS``, ``NIGHT``, and ``HEAVY``.  The category
  ``BASIC`` is the default test test category given to all test executables
  and tests that don't specify the ``CATEGORIES`` argument.

.. _Test Test Category CONTINUOUS:

* Tests marked **CONTINUOUS** also represent importantly functionality but are
  typically not run in `Pre-Push CI Testing`_ testing but instead are run in
  `Post-Push CI Testing`_, `Nightly Testing`_, and other types of testing.
  Tests marked as ``CONTINUOUS`` are enabled for the values of
  ``${PROJECT_NAME}_TEST_CATEGORIES`` equal to ``CONTINUOUS``, ``NIGHT``, and
  ``HEAVY``.  A test may be marked ``CONTINUOUS`` and not ``BASIC`` for a few
  different reasons.  For example, the code needed to run the test may take
  too long to build or the test itself may take too long to run in order to
  afford including it in `Pre-Push CI Testing`_.

.. _Test Test Category NIGHTLY:

* Tests marked **NIGHTLY** usually take even longer to build and/or run than
  ``CONTINUOUS`` tests and therefore are too expensive to include in
  `Post-Push CI Testing`_.  Tests may also be marked as ``NIGHTLY`` even if
  they might run relatively fast if there is a desire to not cause the CI
  server to fail if these tests fail.  In this case, the decision is to take
  the testing and maintenance of these tests and the capabilities they
  represent "offline" so that they don't influence the daily development cycle
  for the project but instead are addressed in a "secondary feedback loop".
  Tests marked as ``NIGHTLY`` are enabled for the values of
  ``${PROJECT_NAME}_TEST_CATEGORIES`` equal to ``NIGHT``, and ``HEAVY``.

.. _Test Test Category HEAVY:

* Tests marked **HEAVY** are usually reserved for very expensive tests that
  are too expensive to run nightly.  ``HEAVY`` tests require more testing
  resources and therefore may only be run on a fully optimized build and/or
  run less frequently.  Tests marked as ``HEAVY`` are enabled only for the
  value of ``${PROJECT_NAME}_TEST_CATEGORIES`` equal to ``HEAVY``.

.. _Test Test Category PERFORMANCE:

* Tests marked **PERFORMANCE** are a special category of tests that are
  specially designed to measure the serial (non-parallel) run-time performance
  of parts of the software (see `Performance Testing`_).  Tests marked as
  ``PERFORMANCE`` are enabled only for the value of
  ``${PROJECT_NAME}_TEST_CATEGORIES`` equal to ``PERFORMANCE``.

Every TriBITS project has a default setting for
``${PROJECT_NAME}_TEST_CATEGORIES`` that is set for a basic ``cmake``
configure of the project (see `${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT`_ for
more details).  In addition, the different testing processes described in the
section `Nested Layers of TriBITS Project Testing`_ set this to different
values.


Nested Layers of TriBITS Project Testing
----------------------------------------

Now that the different types of `Test Classifications for Repositories,
Packages, and Tests`_ have been defined, this section describes how these
different test-related classifications are used to select repositories,
packages (and code) and tests to run in the standard project testing
processes.  More than any other section in this document, this section will
describe and assume a certain class of software development processes (namely
agile processes) where testing and *continuous integration* (CI) are critical
components.  However, detailed descriptions of these processes are deferred to
the later sections `Pre-push Testing using checkin-test.py`_ and `TriBITS
CTest/CDash Driver`_.

The standard TriBITS-supported project testing processes are:

* `Pre-Push CI Testing`_
* `Post-Push CI Testing`_
* `Nightly Testing`_
* `Heavy Testing`_
* `Performance Testing`_

.. ToDo: Discuss why we would want to create standardized test cases?  The
.. answer is that it standardized how testing is done across multiple TriBITS
.. projects so that larger meta-projects and run tests in a consistent way.

.. ToDo: Show a diagram for these different testing processes and show how
.. they nest inside of each other for the most part.

These standard testing processes are outlined in more detail below and show
how the different test-related categories are used to define each of these.

.. _Pre-Push CI Testing:

**Pre-Push CI Testing**

The first level of testing is *Pre-Push CI Testing* that is performed before
changes to the project are pushed to the master branch(es) in the global
repository(s).  With TriBITS, this type of testing and the following push is
typically done using the `checkin-test.py`_ tool.  This category of testing is
described in much more detail in `Pre-push Testing using checkin-test.py`_.
All of the "default builds" used with the ``checkin-test.py`` tool select
repositories, SE packages and code, and individual tests using the following
test-related classifications:

=========================  ==================  ====================================
   Classification Type        Classification           (See Reference)
=========================  ==================  ====================================
Repository Test Classif.   ``Continuous``      (`Repository Test Continuous`_)
SE Package Test Group      ``PT``              (`PT`_)
Test Test Category         ``BASIC``           (`Test Test Category BASIC`_)
=========================  ==================  ====================================

Typically a TriBITS project will define a "standard development environment"
which is comprised of a standard compiler (e.g. GCC 4.6.1), TPL versions
(e.g. OpenMPI 1.4.2, Boost 4.9, etc.), and other tools (e.g. cmake 3.10.0, git
1.8.2, etc.).  This standard development environment is expected to be used to
test changes to the project's code before any push.  By using a standard
development environment, if the code builds and all the tests pass for the
"default" pre-push builds for one developer, then that maximizes the
probability that the code will also build and all tests will pass for every
other developer using the same development environment.  This is critical to
keep the development team maximally productive.  Portability is also important
for most projects but portability testing is best done in a secondary feedback
look using `Nightly Testing`_ builds.  TriBITS has some support for helping to
set up a standard software development environment as described in section
`TriBITS Development Toolset`_.

The basic assumption of all CI processes (including the one described here) is
that if anyone pulls the project's development sources at any time, then all
of the code will build and all of the tests will pass for the "default" build
cases.  For a TriBITS project, this means that the project's
``--default-builds`` (see above) will all pass for every ``PT`` package.  All
of these software development processes make this basic assumption and agile
software development methods fall apart if this is not true.

.. _Post-Push CI Testing:

**Post-Push CI Testing**

After changes are pushed to the master branch(es) in the global repository(s),
*Post-Push CI Testing* is performed where a CI server detects the changes and
immediately fires off a CI build using CTest to test the changes and the
results are posted to a CDash server (in the "Continuous" section on the
project's dashboard page).  This process is driven by CTest driver code that
calls `TRIBITS_CTEST_DRIVER()`_ as described in the section `TriBITS
CTest/CDash Driver`_.  Various types of specific CI builds can be constructed
and run (see `CTest/CDash CI Server`_) but these post-push CI builds typically
select repositories, SE packages and code, and individual tests using the
following test-related classifications:

=========================  ==================  ====================================
   Classification Type        Classification           (See Reference)
=========================  ==================  ====================================
Repository Test Classif.   ``Continuous``      (`Repository Test Continuous`_)
SE Package Test Group      ``PT`` & ``ST``     (`PT`_ and `ST`_)
Test Test Category         ``CONTINUOUS``      (`Test Test Category CONTINUOUS`_)
=========================  ==================  ====================================

Post-push CI testing would assume to use the same standard development
environment as used for `Pre-Push CI Testing`_.  Also, the project may also
choose to run additional automated post-push CI builds that exactly match the
pre-push CI default builds to help check on the health of these builds
continuously and not just rely on the development team to always perform the
pre-push CI builds correctly before pushing.

.. _Nightly Testing:

**Nightly Testing**

In addition to pre-push and post-push CI testing, a typical TriBITS project
will set up multiple *Nightly Testing* builds (or once-a-day builds, they
don't need to only run at night).  These builds are also driven by CTest
driver scripts described in the section `TriBITS CTest/CDash Driver`_ and post
results to the project's CDash server (in the "Nightly" section on the
project's dashboard page).  Nightly builds don't run in a continuous loop but
instead are run once a day (e.g. driven by a cron job) and there tends to be
many different nightly build cases that test the project using different
compilers (e.g.  GCC, Intel, Microsoft, etc., and different versions of each),
different TPL versions (e.g. different OpenMPI versions, different MPICH
versions, etc.), different platforms (e.g. Linux, Windows, etc.), and varying
many other options and settings on these different platforms.  What all
nightly builds have in common is that they tend to select repositories, SE
packages and code, and individual tests using the following test-related
classifications:

=========================  ==================  ====================================
   Classification Type        Classification           (See Reference)
=========================  ==================  ====================================
Repository Test Classif.   ``Nightly``         (`Repository Test Nightly`_)
SE Package Test Group      ``PT`` & ``ST``     (`PT`_ and `ST`_)
Test Test Category         ``NIGHTLY``         (`Test Test Category NIGHTLY`_)
=========================  ==================  ====================================

The nightly builds comprise the basic "heart beat" for the project.

.. _Heavy Testing:

**Heavy Testing**

*Heavy Testing* builds are just an extension to the `Nightly Testing`_ builds
that add on more expensive tests marked using the `Test Test Category HEAVY`_.
For projects that define heavy tests and heavy builds, individual test cases
may be allowed to take 24 hours or longer to run so they can't even be run
every day in nightly testing.  What standard heavy builds have in common is
that they tend to select repositories, SE packages and code, and individual
tests using the following test-related classifications:

=========================  ==================  ====================================
   Classification Type        Classification           (See Reference)
=========================  ==================  ====================================
Repository Test Classif.   ``Nightly``         (`Repository Test Nightly`_)
SE Package Test Group      ``PT`` & ``ST``     (`PT`_ and `ST`_)
Test Test Category         ``HEAVY``           (`Test Test Category HEAVY`_)
=========================  ==================  ====================================

Project developer teams should strive to limit the number of test cases that
are marked as ``HEAVY`` since these tests will typically *not* get run in very
may builds or may not be run every day and developers will tend to never
enable them when doing more extensive testing using ``--st-extra-builds`` with
the `checkin-test.py`_ tool in extended pre-push testing.

.. _Performance Testing:

**Performance Testing** 

*Performance Testing* builds are a special class of builds that have tests
that are specifically designed to test the run-time performance of a particular
piece of code or algorithm.  These tests tend to be sensitive to loads on the
machine and therefore typically need to be run on an unloaded machine for
reliable results.  Details on how to write good performance tests with hard
pass/fail time limits is beyond the scope of this document.  All TriBITS does
is to define the special `Test Test Category PERFORMANCE`_ to allow TriBITS
packages to declare these tests in a consistent way so that they can be run
along with performance tests defined in other TriBITS packages.  From a
TriBITS standpoint, all performance testing builds would tend to select
repositories, SE packages and code, and individual tests using the following
test-related classifications:

=========================  ==================  ====================================
   Classification Type        Classification           (See Reference)
=========================  ==================  ====================================
Repository Test Classif.   ``Nightly``         (`Repository Test Nightly`_)
SE Package Test Group      ``PT`` & ``ST``     (`PT`_ and `ST`_)
Test Test Category         ``PERFORMANCE``     (`Test Test Category PERFORMANCE`_)
=========================  ==================  ====================================

.. ToDo: I need to set up automated testing for TriBITS to use as the example
.. for all of these types of testing.  There is no better example that one
.. that actually works.  It would also be nice to have a snapshot repo of
.. TribitsExampleProject that also had this testing enabled for it but I am
.. not sure that really makes sense.


.. _checkin-test.py:

Pre-push Testing using checkin-test.py
--------------------------------------

CMake provides the integrated tool CTest (executable ``ctest``) which is used
to define and run different tests.  However, a lot more needs to be done to
effectively test changes for a large project before pushing to the master
branch(es) in the main repository(s).  Things get especially complicated and
tricky when multiple version-control (VC) repositories are involved.  The
TriBITS system provides the tool ``checkin-test.py`` for automating the
process of:

1) Determining if the local VC repo(s) are ready to integrate with the remote
   master branch(es),

2) Pulling and integrating the most current changes from the remote VC repo(s),

3) Figuring out what TriBITS packages need to be enabled and testing (by
   examining VC file diffs),

4) Configuring only the necessary TriBITS packages and tests (and their
   downstream dependencies by default) and building, running tests, and
   reporting results (via email), and

5) Only if all specified builds and tests pass, amending the last commit
   message with the test results, then pushing local commits to the remove VC
   repo(s) and sending out summary emails.

There are several advantages to using a project's ``checkin-test.py`` tool for
pushing changes to the main development branch which include:

a) provides a consistent definition for "okay to push" for all developers

b) protects other developers from pulling badly broken code

c) reduces the number of packages that need to be tested by automatically
   determining the set based on changed files and package dependencies

d) avoids developer mistakes in performing repetitive tasks and forgetting
   important steps in the process

e) marks a set of working commits that are safe to search with ``git bisect``
   to find problems (see `Using Git Bisect with checkin-test.py workflows`_)

When using the ``checkin-test.py`` tool, every TriBITS project defines one or
more "default builds" (specified through the ``--default-builds`` argument)
for pre-push CI testing that form the criteria for if it is okay to push code
changes or not.  The "default builds" select repositories, SE packages and
code, and individual tests as described in `Pre-Push CI Testing`_.  A TriBITS
project defines its default pre-push builds using the file
`<projectDir>/project-checkin-test-config.py`_.  For an example, the file
`TribitsExampleProject`_/``project-checkin-test-config.py`` is shown below:

.. include:: ../../examples/TribitsExampleProject/project-checkin-test-config.py
   :literal:

This gives ``--default-builds=MPI_DEBUG,SERIAL_RELEASE``.  As shown, typically
two default builds are defined so that various options can be toggled between
the two builds.  Typical options to toggle include enabling/disabling MPI and
enabling/disabling run-time debug mode checking (i.e. toggle
``${PROJECT_NAME}_ENABLE_DEBUG``).  Typically, other important options will
also be toggled between these two builds.

Note that both of the default builds shown above, including the ``MPI_DEBUG``
build, actually set optimized compiler flags with
``-DCMAKE_BUILD_TYPE:STRING=RELEASE``.  What makes the ``MPI_DEBUG`` build a
"debug" build is turning on optional run-time debug-mode checking, not
disabling optimized code.  This is important so that the defined tests run
fast.  For most projects, the default pre-push builds should **not** be used
to debug-enabled code which is suitable to run through a debugger
(e.g. ``gdb``).  Instead, these "debug" builds are designed to test changes to
the project's code efficiently before pushing changes.  Typically, a
development team should not have to test the chosen compiler's ability to
generate non-optimized debug code and suffer slower test times before pushing.

Note that turning on ``-DTribitsExProj_ENABLE_CHECKED_STL=ON`` as shown above
can only be used when the TPLs have no C++ code using the C++ STL or if that
particular build points to C++ TPLs also compiled with checked STL enabled.
The `TribitsExampleProject`_ default builds do not depend on any C++ TPLs that
might use the C++ STL so enabling this option adds additional positive
debug-mode checking for C++ code.

The ``checkin-test.py`` tool is a fairly sophisticated piece of software that
is well tested and very robust.  The level of testing of this tool is likely
greater than any of the software that it will be used to test (unless the
project is a real-time flight control system or nuclear reactor control system
or something).  This is needed so as to provide confidence in the developers
that the tool will only push their changes if everything checks out as it
should.  There are a lot of details and boundary cases that one has to
consider and a number of use cases that need to be supported by such a tool.
For more detailed documentation, see `checkin-test.py --help`_.

Note that the ``checkin-test.py`` tool can also be used to implement
"poor-man's" post-push testing processes as described in `Post-Push CI and
Nightly Testing using checkin-test.py`_.  However, most software projects will
want to go with the more elaborate and more feature-full CTest/CDash system
described in `TriBITS CTest/CDash Driver`_.

.. ToDo: Describe the standard workflow for using the checkin-test.py tool.

.. ToDo: Describe why the --default-builds must only include PT code and not
.. ST due to changing the behavior of the PT software.  As an example, discuss
.. Boost usage in Teuchos and how that changes behavior of a few things.
.. Another example is BinUtils.

.. ToDo: Describe the system for mapping changed files to changed packages.
.. Discuss how important the directory structure is.

.. ToDo: Describe the local test, remote push process by pulling from the
.. Trilinos developers webpage.

.. ToDo: Describe the role that checkin-test.py plays in multi-repo ACI


TriBITS CTest/CDash Driver
--------------------------

The TriBITS system uses a sophisticated and highly customized CTest -S driver
script to test TriBITS projects and submit results to a CDash server.  The
primary code for driving this is contained in the CTest function
`TRIBITS_CTEST_DRIVER()`_ contained in the file
``TribitsCTestDriverCore.cmake``.  This script loops through all of the
specified TriBITS packages for a given TriBITS project and does a configure,
built, and test and then submits results to the specified CDash server
incrementally.  If the configure or library build of any upstream TriBITS
package fails, then that TriBITS package is disabled in all downstream TriBITS
package builds so as not to propagate already known failures.  Each TriBITS
top-level package is assigned its own CDash regression email address (see
`CDash regression email addresses`_) and each package configure/build/test is
given its own row for the package build in the CDash server.  A CTest script
using `TRIBITS_CTEST_DRIVER()`_ is run in one of three different modes.
First, it can run standard once-a-day, from-scratch builds as described in
`CTest/CDash Nightly Testing`_.  Second, it can run as a CI server as
described in `CTest/CDash CI Server`_.  Third, it can run in experimental mode
testing a local repository using the TriBITS-defined `make dashboard`_ target.

.. ToDo: Discuss the behavior of Continuous, Nightly, and Experimental
.. CTest/CDash tracks and what that means in terms of sending out CDash error
.. emails and allowing for multiple builds of the same name in a single test
.. day.

.. ToDo: Document CTEST_TEST_TIMEOUT and DART_TESTING_TIMEOUT and how these
.. interact.

.. ToDo: Import and format the contents of tribits/ctest/README related to
.. this topic.


CTest/CDash Nightly Testing
+++++++++++++++++++++++++++

When a TriBITS CTest script using `TRIBITS_CTEST_DRIVER()`_ is run in
"Nightly" testing mode, it builds the project from scratch package-by-package
and submits results to the TriBITS project's CDash project on the designated
CDash server.

.. ToDo: Give example of the package-by-package build.

CTest/CDash CI Server
+++++++++++++++++++++

When a TriBITS ctest driver script is used in continuous integration (CI)
mode, it starts every day with a clean from-scratch build and then performs
incremental rebuilds as new commits are pulled from the master branch in the
main repository(s).  In this mode, a continuous loop is performed after the
initial baseline build constantly pulling commits from the master git
repository(s).  If any package changes are detected (looking at git file
diffs), then the tests and examples for those packages and all downstream
packages are enabled and run using a reconfigure/rebuild.  Since all of the
upstream package libraries are already built, this rebuild and retest can take
place in a fraction of the time of a complete from-scratch build and test of
the project.

.. ToDo Give some examples of the behavior of this process for the
.. ReducedMockTrilins project.  We could show example output from the script
.. that shows how it behaves.

.. ToDo: Discuss the usage of ${REPOSITORY_NAME}_NO_IMPLICIT_PACKAGE_ENABLE
.. and ${REPOSITORY_NAME}_NO_IMPLICIT_PACKAGE_ENABLE_EXCEPT to allow the CI
.. server to directly processing and running tests for upstream packages which
.. are not primary packages in the TriBITS meta-project.

.. ToDo: Mention that failing packages are always explicitly processed even if
.. the current CI iteration does not pull any changes.

.. ToDo: Give some examples from ReducedMockTrilinos of how file changes
.. result in package enables and rebuilds.


TriBITS CDash Customizations
----------------------------

CDash is not currently designed to accommodate multi-package, multi-repository
VC projects in the way they are supported by TriBITS.  However, CDash provides
some ability to customize a CDash project and submits to address missing
features.  Each TriBITS package is given a CTest/CDash "Label" with the name
of the TriBITS package.  CDash will then aggregate the different package
configure/build/test runs for each package into aggregated "builds".  The
commits pulled for each of the extra VC repos listed in the
`<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file are shown in an
uploaded CDash "Notes" file for each TriBITS package configure/build/test
submit.  This uploaded "Notes" file also contains a cleaned up version of the
CMakeCache.txt file as well as a copy of the ``ctest -s`` script that ran the
case.

CDash offers numerous features such as the ability to construct a number of
different types of queries and is extremely helpful in using past test data.

.. ToDo: Describe more on what CDash provides and how TriBITS uses it.  This
.. could be a *long* section with a lot of screen shots.


CDash regression email addresses
++++++++++++++++++++++++++++++++

Every TriBITS Package has a regression email address associated with it that
gets uploaded to a CDash project on a CDash server that is used to determine
what email address to use when a package has configure, build, or test
failures.  Because of the complex organizational nature of different projects
and different integration models, a single static email address for a given
package in every project build is not practical.

The TriBITS system allows for a package's regression email address to be
specified in the following order of precedence:

.. _${REPOSITORY_NAME}_REPOSITORY_OVERRIDE_PACKAGE_EMAIL_LIST:

1) **${REPOSITORY_NAME}_REPOSITORY_OVERRIDE_PACKAGE_EMAIL_LIST** (typically
defined in `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_): Defines a
single email address for all packages for the repository
``${REPOSITORY_NAME}`` and overrides all other package email regression
specification variables.  This is typically used by a meta-project to redefine
the regression email addresses for all of the packages in an externally
developed repository.

2) **REGRESSION_EMAIL_LIST** (defined in
`<packageDir>/cmake/Dependencies.cmake`_): Package-specific email address
specified in the package's ``Dependencies.cmake`` file using
`TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_.

.. _${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESS_BASE:

3) **${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESS_BASE** (set in
`<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_): A base email address
specified at the Repository level creating package-specific email addresses
(e.g. ``<lower-case-package-name>-regression@some.repo.gov``, where
``${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESS_BASE=some.repo.gov``).
This variable is used, for example, by the Trilinos project to provide
automatic regression email addresses for packages.

.. _${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESS:

4) **${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESS** (set in
`<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_): A single email address
for all packages specified at the Repository level
(e.g. ``my-repo-regression@some.repo.gov``).  This variable is used for
smaller repositories with smaller development groups who just want all
regression emails for the repository's packages going to a single email
address.  This reduces the overhead of managing a bunch of individual package
email addresses but at the expense of spamming too many people with CDash
failure emails.

.. _${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESS_BASE:

5) **${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESS_BASE** (set in
`<projectDir>/cmake/ProjectDependenciesSetup.cmake`_): A base email address
specified at the Project level creating package-specific email addresses
(e.g. ``<lower-case-package-name>-regression@some.project.gov``, where
``${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESS_BASE=some.project.gov``).  If not
already set, this variable will be set to
``${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESS_BASE`` for the first
repository processed that has this set. This behavior is used, for example by
the Trilinos project to automatically assign email addresses for add-on
packages and was added to maintain backward compatibility.

.. _${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESS:

6) **${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESS** (set in
`<projectDir>/cmake/ProjectDependenciesSetup.cmake`_): A single default email
address for all packages specified at the Project level
(e.g. ``my-project-regression@some.project.gov``).  If not already set, this
variable will be set to
``${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESS`` for the first
repository processed that has this set.  Every meta-project should set this
variable so that it will be the default email address for any new package
added.

WARNING: If any of the email lists or URL string variables listed above are
set to ``"OFF"`` or ``"FALSE"`` (or some other value that CMake interprets as
false, see `CMake Language Overview and Gotchas`_) then the variables are
treated as empty and not set.

If a TriBITS project does not use CDash, then no email address needs to be
assigned to packages at all (and therefore none of the above variables need be
set).

As a general rule, repository-level settings override project-level settings
and package-level settings override both.  Also, a project can redefine a
repository's regression email list settings by resetting the variables in the
project's `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_ file.

All of the email dependency management logic must be accessible by just running
the macro::

    TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()

The above email address configuration variables are read from the Repository
and Project files `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ and
`<projectDir>/cmake/ProjectDependenciesSetup.cmake`_, respectively.  The
``RepositoryDependenciesSetup.cmake`` files are read first in the specified
repository order followed up by reading the ``ProjectDependenciesSetup.cmake``
file.  In this way, the project can override any of the repository settings.

In review, the precedence order for how regression email addresses are
selected for a given package is:

1) Package-specific email list is selected if defined (unless an override is
   in place).

2) Repository-level option is selected over a project-level option.

3) Default email form with repository or project address base is selected over
   single repository or project email address.

4) If none of the above are selected, then no email address is assigned for a
   given package.

What the above setup does is it results in the TriBITS system (in the
`TRIBITS_CTEST_DRIVER()`_ function called in a ``ctest -S`` script) creating a
file called ``CDashSubprojectDependencies.xml`` (which contains the list of
TriBITS packages, which CDash calls "subprojects", and email address for each
package to send regression emails to) and that file gets sent to the CDash
server. CDash then takes this file and creates, or updates, a set of CDash
users (same name and password as the email list address) and sets up a mapping
of Labels (which are used for TriBITS package names) to CDash user emails
addresses. CDash is automatically set up to process this XML file and create
and update CDash users.  There are several consequences of this implementation
of which project maintainers need to be aware.

First, **one should not list the email address for a CDash user account
already on CDash**.  This is because labels will be added for the TriBITS
packages that this email address is associated with and CDash emails will not
be sent out for any other TriBITS packages, no matter what setting that CDash
user account has.  Therefore, one should only list email addresses as CDash
regression email lists that are not already CDash user accounts and wish to be
maintained separately.  For example, if there is an email list that one wants
to have CDash emails sent to but is already a CDash user account, then one can
create another email list (e.g. using Mailman) which can then be registered
with the TriBITS packages in the ``CDashSubprojectDependencies.xml`` file and
then that new email list forward email to the target email list.

Second, the CDash implementation currently is not set up to remove labels from
existing users when an email address is disassociated with a TriBITS package
in the ``CDashSubprojectDependencies.xml`` file.  Therefore, **if one changes
a TriBITS package's CDash regression email address then one needs to manually
remove the associated labels from the old email address**.  CDash will not
remove them automatically.  Otherwise, email will continue to be sent to that
email address for that package.

Therefore, **to change the mapping of CDash regression email addresses to
TriBITS packages, one must perform the following actions**:

1) Change the TriBITS CMake files as described above that will result in the
   desired email addresses in the ``CDashSubprojectDependencies.xml``
   file. One can debug this by generating the file
   `<Project>PackageDependencies.xml`_ by using the cmake -P script
   `TribitsDumpDepsXmlScript.cmake`_.

2) Log onto the CDash server using an administrator account and then remove
   the auto-generated account for the CDash user email address for which
   labels are being removed (i.e. no longer associated with a TriBITS
   package).  This is needed since CDash seems to be unable to remove labels
   from an existing CDash user (however this might be fixed in a current
   version of CDash).

3) The next time a CDash submit is performed by a CTest driver script calling
   `TRIBITS_CTEST_DRIVER()`_, the CDash user associated with the mail list
   with labels being removed will get automatically recreated with the right
   list of labels (according to the current
   ``CDashSubprojectDependencies.xml`` file).  Also, any new CDash users for
   new email addresses will be created.

Hopefully that should be enough information to manage the mapping of CDash
regression email lists to TriBITS packages for single and multi-repository
TriBITS projects.


Multi-Repository Support
========================

TriBITS has built-in support for projects involving multiple `TriBITS
Repositories`_ which contain multiple `TriBITS Packages`_ (see `How to set up
multi-repository support`_).  The basic configuration, build, and test of such
projects requires only raw CMake/CTest, just like any other CMake project (see
`TriBITS System Project Dependencies`_).  Every TriBITS project automatically
supports tacking on add-on TriBITS packages and TPLs through the
`${PROJECT_NAME}_EXTRA_REPOSITORIES`_ cmake cache variable as described in
`Enabling extra repositories with add-on packages`_.  In addition, a TriBITS
project can be set up to pull in other TriBITS Repositories using the
`<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file.  A special form of
this type of project is a `TriBITS Meta-Project`_ that contains no native
packages or TPLs of its own.  The ability to create meta-projects out of
individual TriBITS repositories allows TriBITS to be used to provide
coordinated builds (or meta-builds) of large aggregations of software.

To help set up a full-featured development environment (i.e. not just the
basic configure, build, test, and install) for TriBITS projects with multiple
repositories, TriBITS provides some extra development tools implemented using
Python which are provided in the "extended" parts of TriBITS (see
`TriBITS/tribits/ Directory Contents`_).  The primary tools supporting
multi-repository projects are the Python tools `clone_extra_repos.py`_,
`gitdist`_, and `checkin-test.py`_.

To demonstrate, consider the TriBITS meta-project with the following
``ExtraRepositoriesList.cmake`` file::

  TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES(
    ExtraRepo1  ""  GIT  git@someurl.com:ExtraRepo1  ""  Continuous
    ExtraRepo2  "ExtraRepo1/ExtraRepos2"  GIT  git@someurl.com:ExtraRepo2
                                                   NOPACKAGES  Continuous
    ExtraRepo3  ""  GIT  git@someurl.com:ExtraRepo3  ""  Nightly
    )

Once cloned, the directories would be laid out as::

  MetaProject/
    .git/
    .gitignore
    ExtraRepo1/
      .git/
      ExrraRepo2/
        .git/
    ExtraRepo3/
      .git/

.. _clone_extra_repos.py:

The tool **clone_extra_repos.py** is used to clone the extra repositories for
a multi-repositories TriBITS proejct.  It reads the repository URLs and
destination directories from the file
`<projectDir>/cmake/ExtraRepositoriesList.cmake`_ and does the clones.  For
example, to clone all the repos for the ``MetaProject`` project, one would use
the commands::

  $ git clone git@someurl.com:MetaProject
  $ cd MetaProject/
  $ ./cmake/tribits/ci_support/clone_extra_repos.py 

which produces the output like::

  ...  
  
  ***
  *** Clone the selected extra repos:
  ***
  
  Cloning repo ExtraRepo1 ...
  
  Running: git clone git@someurl.com:ExtraRepo1 ExtraRepo1
  
  Cloning repo ExtraRepo2 ...
  
  Running: git clone git@someurl.com:ExtraRepo2 ExtraRepo1/ExtraRepo2
  
  Cloning repo ExtraRepo3 ...
  
  Running: git clone git@someurl.com:ExtraRepo3 ExtraRepo3

See `clone_extra_repos.py --help`_ for more details.

.. _gitdist:

Once cloned, one needs to work with the multiple repositories to perform basic
VC operations.  For this, TriBITS provides the tool **gitdist** which is a
simple stand-alone Python script that distributes a git command across a set
of git repos.  This tool is not specific to TriBITS but it is very useful for
dealing with TriBITS projects with multiple repositories.  It only requires a
local base git repo and a set of zero or more git repos cloned under it.

To use ``gitdist`` with this aggregate meta-project, one would first set up
the file ``MetaProject/.gitdist`` (or a version-controlled
``MetaProject.gitdist.default`` file) which would contain the lines::

  ExtraRepo1
  ExtraRepo1/ExtraRepo2
  ExtraRepo3

and one would set up the tracked ignore file ``MetaProject/.gitignore`` which
contains the lines::

  /ExtraRepo1/
  /ExtraRepo1/ExtraRepo2/
  /ExtraRepo3/

To use ``gitdist``, one would put ``gitdist`` into their path and also set up
the command-line shell aliases ``gitdist-status`` and ``gitdist-mod`` (see
`gitdist --dist-help=aliases`_).

Some of the aggregate commands that one would typically run under the base
``MetaProject/`` directory are::

  # See status of all repos at once
  gitdist-status

  # Pull updates to all
  gitdist pull

  # Push local commits to tracking branches
  gitdist push

The tool ``gitdist`` is provided under TriBITS directory::

  cmake/tribits/python_utils/gitidst   

and can be installed by the `install_devtools.py`_ tool (see `TriBITS
Development Toolset`_).  See `gitdist documentation`_ for more details.

For projects with a standard set of extra repositories defined in the
`<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file, the
``checkin-test.py`` tool only requires passing in the option
``--extra-repos-file=project`` and ``--extra-repos-type=Continuous`` (or
``Nightly``, see `Repository Test Classification`_) and it will automatically
perform all of the various actions for all of the selected repositories.  See
`checkin-test.py`_ and `checkin-test.py --help`_ for more details.

To keep track of compatible versions of the git repos, TriBITS provides
support for a ``<Project>RepoVersion.txt`` file.  Any TriBITS project can
generate this file automatically by setting the option
`${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE`_.  For the above example
``MetaProject``, this file looks like::

  *** Base Git Repo: MetaProject
  e102e27 [Mon Sep 23 11:34:59 2013 -0400] <author0@someurl.com>
  First summary message
  *** Git Repo: ExtraRepo1
  b894b9c [Fri Aug 30 09:55:07 2013 -0400] <author1@someurl.com>
  Second summary message
  *** Git Repo: ExtraRepo1/ExtraRepo2
  97cf1ac [Thu Dec 1 23:34:06 2012 -0500] <author2someurl.com>
  Third summary message
  *** Git Repo: ExtraRepo3
  cd4a3af [Mon Mar 9 19:39:06 2013 -0400] <author3someurl.com>
  Fourth summary message

This file gets created in the build directory, gets echoed in the configure
output, gets installed into the install directory, and get added to the source
distributions tarball.  It also gets pushed up to CDash for all automated
builds.  The tool `gitdist`_ can then use this file to checkout and tag
compatible versions, difference two versions of the meta-project, etc. (see
`gitdist documentation`_ for more details on git operations).

The TriBITS approach to managing multiple VC repos described above works well
for around 20 or 30 VC repos but is likely not a good solution for many more
git repos.  For larger numbers of VC repos, one should consider nested
integration creating snapshot git repos (e.g. using the tool
`snapshot-dir.py`_) that aggregate several related repositories into a single
git repo.  Another approach might be to use git submodules.  (However, note
that the TriBITS tools and processes described here are **not** currently set
up to support aggregate VC repos that use git submodules.)  The design
decision with TriBITS was to explicitly handle the different git VC repos by
listing them in the `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file and
then using the simple, easy to understand, tools `clone_extra_repos.py`_ and
`gitdist`_.  There are advantages and disadvantages to explicitly handling the
different git repos as is currently employed by the TriBITS software
development tools verses using git submodules.  It is possible that TriBITS
will add support for aggregate git repos using git submodules in the future
but only if there are important projects that choose to use them.  The
discussion of these various approaches and strategies to dealing with
aggregate repos is beyond the scope of this document.


Development Workflows
======================

In this section, the typical development workflows for a TriBITS project are
described.  First, the `Basic Development Workflow`_ for a single-repository
TriBITS project is described.  This is followed up with a slightly more
complex `Multi-Repository Development Workflow`_.


Basic Development Workflow
--------------------------

The basic development workflow of a TriBITS project is not much different than
with any other CMake project that uses CTest to define and run tests.  One
pulls updates from the master VC repo then configures with ``cmake``, and
iteratively builds, runs tests, adds files, changes files, does a final test,
then pushes updates.  The major difference is that a well constructed
development process will use the `checkin-test.py`_ tool to test and push all
changes that affect the build or the tests.  The basic steps in configuring,
building, running tests, etc., are given in the project's
`<Project>BuildReference`_. file (see `Project-Specific Build Reference`_).

Multi-Repository Development Workflow
-------------------------------------

The development workflow for a project with multiple VC repos is very similar
to a project with just a single VC repo if the project provides a standard
`<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file.  The major difference
is in making changes, creating commits, etc.  The `gitdist`_ tool makes these
steps easier and has been shown to work fairly well for up to 20 extra VC
repos (as used in the CASL VERA project).  The `checkin-test.py`_ tool
automatically handles all of the details of pulling, diffing, pushing etc. to
all the VC repos.



.. ToDo: Discuss usage of 'gitdist' and the repo clone script.


Howtos
======

This section provides short, succinct lists of the steps to accomplish a few
common tasks.  Extra details are referenced.


How to add a new TriBITS Package
--------------------------------

To add a new TriBITS package, it is recommended to take the template from one
of the `TribitsExampleProject`_ packages that most closely fits the needs of
the new package and modify it for the new package.  For example, the files for
the ``SimpleCxx`` package can be copied one at a time and modified for the new
package.

To add a new TriBITS package (with no subpackages), do the following:

1) Chose a name ``<packageName>`` for the new package and which TriBITS
   repository (``<repoDir>``) to put the package into.  **WARNING!** The
   chosen name ``<packageName>`` must be unique across all TriBITS
   repositories (see `Globally unique TriBITS package names`_).

2) Create the directory ``<repoDir>/<packageDir>/`` for the new package and
   put in skeleton files for `<packageDir>/cmake/Dependencies.cmake`_ and
   `<packageDir>/CMakeLists.txt`_.  Set the desired upstream TPL and SE
   package dependencies in the new ``Dependencies.cmake`` file but initially
   comment out everything in the ``CMakeLists.txt`` file except for the
   `TRIBITS_PACKAGE()`_ and `TRIBITS_PACKAGE_POSTPROCESS()`_ commands.

3) Add a row for the new package to the `<repoDir>/PackagesList.cmake`_ file
   after all of its upstream dependent packages.  If a mistake is made and it
   is listed before one of its upstream dependent packages, the TriBITS CMake
   code will catch this and issue an error.

4) Configure the TriBITS project enabling the new empty package
   ``<packageName>``.  This will enable the listed dependencies.

5) Incrementally fill in the package's ``CMakeLists.txt`` files defining
   libraries, executables, tests and examples.  The project should be built
   and tests run as new pieces are added.

Once the new package is defined, downstream SE packages can define
dependencies on this new package.

.. ToDo: Expand on the above bullets a lot!


How to add a new TriBITS Package with Subpackages
-------------------------------------------------

Adding a new package with subpackages is similar to adding a new regular
package described in `How to add a new TriBITS Package`_.  Again, it is
recommended that one copies an example package from `TribitsExampleProject`_.
For example, one could copy files and directories from the example package
``WithSubpackages``.

To add a new TriBITS package with packages, do the following:

1) Chose a name ``<packageName>`` for the new package and which TriBITS
   repository (``<repoDir>``) to put the package into.  **WARNING!** The
   chosen name ``<packageName>`` must be unique across all TriBITS
   repositories (see `Globally unique TriBITS package names`_).

2) Create the directory ``<repoDir>/<packageDir>/`` for the new package and
   put in skeleton files for `<packageDir>/cmake/Dependencies.cmake`_ and
   `<packageDir>/CMakeLists.txt`_.  Initially don't define any subpackages and
   comment out everything in the ``CMakeLists.txt`` file except for the
   `TRIBITS_PACKAGE()`_ and `TRIBITS_PACKAGE_POSTPROCESS()`_ commands.

3) Add a new row for the new package to the `<repoDir>/PackagesList.cmake`_
   file after all of the upstream dependencies of its to-be-defined
   subpackages.

4) Configure the TriBITS project enabling the new empty package
   ``<packageName>``.

5) Incrementally add the subpackages as described in `How to add a new TriBITS
   Subpackage`_, filling out the various ``CMakeLists.txt`` files defining
   libraries, executables, tests and examples.

Once the new SE packages are defined, downstream SE packages can define
dependencies on these.

How to add a new TriBITS Subpackage
-----------------------------------

Given an existing top-level TriBITS package that is already broken down into
subpackages (see `How to add a new TriBITS Package with Subpackages`_), adding
a new subpackage does not require changing any project-level or
repository-level files.  One only needs to add the declaration for the new
subpackages in its parent's `<packageDir>/cmake/Dependencies.cmake`_ file then
fill out the pieces of the new subpackage defined in the section `TriBITS
Subpackage Core Files`_.  It is recommended to copy files from one of the
`TribitsExampleProject`_ subpackages in the ``WithSubpackages``
package.

To add a new TriBITS subpackage to a top-level package that already has
subpackages, do the following:

1) Chose a name ``<spkgName>`` for the new subpackage which only has to be
   different than the other subpackages in the parent package.  This name gets
   appended to the parent package's name ``<packageName>`` to form the SE
   package name ``<packageName><spkgName>``.

2) Create the directory ``<packageDir><spkgDir>/`` for the new package and put
   in skeleton files for `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_
   and `<packageDir>/<spkgDir>/CMakeLists.txt`_.  Set the desired upstream TPL
   and SE package dependencies in the new ``Dependencies.cmake`` file but
   initially comment out everything in the ``CMakeLists.txt`` file except for
   the `TRIBITS_SUBPACKAGE()`_ and `TRIBITS_SUBPACKAGE_POSTPROCESS()`_
   commands.

3) Add a row for the new subpackage to the argument
   `SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS`_ in the macro call
   `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_ in the parent package's
   `<packageDir>/cmake/Dependencies.cmake`_ file after all of its upstream
   dependent subpackages.  If a mistake is made and it is listed before one of
   its upstream dependent subpackages, the TriBITS CMake code will catch this
   and issue an error.

4) Configure the TriBITS project enabling the new empty SE package
   ``<packageName><spkgName>``.  This will enable the listed dependencies.

5) Incrementally fill in the subpackage's ``CMakeLists.txt`` files defining
   libraries, executables, tests and examples.  The project should be built
   and tests run as new pieces are added.


How to add a new TriBITS TPL
----------------------------

In order for an SE package to define a dependency on a new TPL (i.e. one that
has not already been declared in the current repo's or an upstream repo's
`<repoDir>/TPLsList.cmake`_ file), one must add and modify a few
repository-level files.

To add a new TriBITS TPL, do the following:

1) Chose a name ``<tplName>`` for the new TPL (must be globally unique across
   all TriBITS repos, see `Globally unique TriBITS TPL names`_) and which
   TriBITS repository (``<repoDir>``) to define the new TPL in.  The right
   repo is usually the one where the package exists that needs the new TPL
   dependency.

2) Create the TPL find module, e.g. ``<repoDir>/tpls/FindTPL<tplName>.cmake``
   (see `TriBITS TPL`_ and `TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()`_
   for details).  List the default required header files and/or libraries that
   must be provided by the TPL.  NOTE: This find module can also (optionally)
   use ``FIND_PACKAGE(<tplName> ...)`` for the default find operation.  For
   details, see `How to use FIND_PACKAGE() for a TriBITS TPL`_.

3) Add a row for the new TPL to the `<repoDir>/TPLsList.cmake`_ file after any
   TPLs that this new TPL may depend on.

4) Configure the TriBITS project enabling the new TPL with
   ``TPL_ENABLE_<tplName>=ON`` and see that the new find module finds the TPL
   correctly at configure time.

5) List the new TPL in the package(s) that need this dependency in the
   package's `<packageDir>/cmake/Dependencies.cmake`_ file (or
   `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_ for a subpackage).

6) If the TPL is an optional TPL for the package, then::

    #cmakedefine HAVE_<PACKAGE_NAME_UC>_<TPL_NAME_UC>

   should be added to the package's ``<packageName>_config.h.in`` file (see
   `TRIBITS_CONFIGURE_FILE()`_) so that the code knows if the TPL is defined
   or not.

7) Use the TPL in the packages that define the dependency on the new TPL,
   configure, test, etc.

.. ToDo: Provide more details on various things.

.. ToDo: Add an example of an optional TPL to some TribitsExampleProject
.. package and refer to it here.


How to use FIND_PACKAGE() for a TriBITS TPL
-------------------------------------------

When defining a ``FindTPL<tplName>.cmake`` file, it is possible (and
encouraged) to utilize ``FIND_PACKAGE(<tplName> ...)`` to provide the default
find operation.  In order for the resulting ``FindTPL<tplName>.cmake`` to
behave consistently between all the various TriBITS TPLs (and allow the
standard TriBITS TPL find overrides) one must use the TriBITS function
`TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE()`_ in combination with the function
`TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()`_.  The basic form of the
resulting TriBITS TPL module file ``FindTPL<tplName>.cmake`` looks like::

  # First, set up the variables for the (backward-compatible) TriBITS way of
  # finding <tplName>.  These are used in case FIND_PACKAGE(<tplName> ...) is
  # not called or does not find <tplName>.  Also, these variables need to be
  # non-null in order to trigger the right behavior in the function
  # TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES().
  SET(REQUIRED_HEADERS <header0> <header1> ...)
  SET(REQUIRED_LIBS_NAMES <libname0> <libname1> ...)
  
  # Second, search for <tplName> components (if allowed) using the standard
  # FIND_PACKAGE(<tplName> ...).
  TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE(<tplName>  <tplName>_ALLOW_PREFIND)
  IF (<tplName>_ALLOW_PREFIND)
    MESSAGE("-- Using FIND_PACKAGE(<tplName> ...) ...") 
    FIND_PACKAGE(<tplName>)
    IF (<tplName>_FOUND)
      # Tell TriBITS that we found <tplName> and there no need to look any further!
      SET(TPL_<tplName>_INCLUDE_DIRS ${<tplName>_INCLUDE_DIRS} CACHE PATH "...")
      SET(TPL_<tplName>_LIBRARIES ${<tplName>_LIBRARIES} CACHE FILEPATH "...")
      SET(TPL_<tplName>_LIBRARY_DIRS ${<tplName>_LIBRARY_DIRS} CACHE PATH "...")
    ENDIF()
  ENDIF()
  
  # Third, call TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()
  TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( <tplName>
    REQUIRED_HEADERS ${REQUIRED_HEADERS}
    REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
    )
  # NOTE: If FIND_PACKAGE(<tplName> ...) was called and successfully found
  # <tplName>, then TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES() will use the
  # already-set variables and just print them out.  This is the final "hook"
  # into the TriBITS TPL system.

With this approach, the ``FindTPL<tplName>.cmake`` module preserves all of the
user behavior described in `Enabling support for an optional Third-Party
Library (TPL)`_ for overriding what TPL components to look for, where to look
and finally to override what is actually used.  That is, if the user sets the
cache variables ``TPL_<tplName>_INCLUDE_DIRS``, ``TPL_<tplName>_LIBRARIES``,
or ``TPL_<tplName>_LIBRARY_DIRS``, then they should be used without question
(which is why the ``SET( ... CACHE ...)`` calls in the above example do not
use ``FORCE``).

If one wants to skip and ignore the standard TriBITS TPL override variables
``<tplName>_INCLUDE_DIRS``, ``<tplName>_LIBRARY_NAMES``, or
``<tplName>_LIBRARY_DIRS``, then one can set::

  SET(<tplName>_FORCE_PRE_FIND_PACKAGE TRUE CACHE BOOL
    "Always first call FIND_PACKAGE(<tplName> ...) unless explicit override")

at the top of the ``FindTPL<tplName>.cmake`` and it will ignore these
variables.  This avoids name classes with the variables
``<tplName>_INCLUDE_DIRS`` and ``<tplName>_LIBRARY_DIRS`` which are often used
in the concrete CMake ``Find<tplName>.cmake`` module files themselves.

For a slightly more complex (but real-life) example, see ``FindTPLHDF5.cmake``
which is:

.. include:: ../../common_tpls/FindTPLHDF5.cmake
   :literal:

Note that some specialized ``Find<tplName>.cmake`` modules do more than just
return a list of include directories and libraries.  Some, like
``FindQt4.cmake`` also return other variables that are used in downstream
packages. therefore, in these cases, ``FIND_PACKAGE(Qt4 ...)`` must be called
on every configure.  In specialized cases such as this, one must write a more
specialized ``FindTPL<tplName>.cmake`` file and can't use the
`TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE()`_ function like shown above.  Such
find modules cannot completely adhere to the standard behavior described in
`Enabling support for an optional Third-Party Library (TPL)`_.


How to add a new TriBITS Repository
-----------------------------------

To add a new TriBITS and/ git VC repository to a TriBITS project that already
contains other extra repositories, do the following:

1) Add a row for the new TriBITS and/or git repo to the file
   `<projectDir>/cmake/ExtraRepositoriesList.cmake`_.  Commit this file.

2) Add an ignore for the extra repo name to the base project's git repo's
   ``.gitignore`` file (see `Multi-Repository Support`_).  Commit this file.

3) Set up the new package dependencies for the new package in the new TriBITS
   repo or make other adjustments.

4) Consider the potential for missing upstream repos and packages by using
   `TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES()`_.

See `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ for more details and
links.

.. ToDo: Provide more detail!


How to add a new TriBITS SE Package dependency
----------------------------------------------

It is often the case where one will want to add a new dependency for an
existing `downstream`_ package to an existing `upstream`_ `TriBITS SE
Package`_.  This can either be a required dependency or an optional
dependency.  Here, we will refer to the downstream SE package as
``<packageName>`` with base directory ``<packageDir>`` and will refer to the
upstream package as ``<upstreamPackageName>``.

The process for adding a new dependency to an existing upstream SE package is
as follows:

1) **Add the name of the upstream package to the downstream package's
   Dependencies.cmake file:** Add ``<upstreamPackagename>`` to the call of
   `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_ in the downstream package's
   `<packageDir>/cmake/Dependencies.cmake`_ file.  If this is to be a required
   library dependency, then ``<upstreamPackagename>`` is added to the
   ``LIB_REQUIRED_PACKAGES`` argument.  Alternatively, if this is to be an
   optional library dependency, then ``<upstreamPackagename>`` is added to the
   ``LIB_OPTIONAL_PACKAGES`` argument.  (For example, see the file
   ``packages/EpetaExt/cmake/Dependencies.cmake`` file in the
   `ReducedMockTrilinos`_ project.)  If only the test and/or example sources,
   and not the package's core library sources, will have the required or
   optional dependency, then ``<upstreamPackagename>`` is added to the
   arguments ``TEST_REQUIRED_PACKAGES`` or ``TEST_OPTIONAL_PACKAGES``,
   respectively.

2) **For an optional dependency, add `HAVE_` preprocessor macro to the
   package's configured header file:** If this is an optional dependency,
   typically a C/C++ processor macro will be added to the package's configured
   `<packageDir>/cmake/<packageName>_config.h.in`_ file using the line::

     #cmakedefine HAVE_<PACKAGE_NAME_UC>_<OPTIONAL_DEP_PACKAGE_NAME_UC>

   (see `HAVE_<PACKAGE_NAME_UC>_<OPTIONAL_DEP_PACKAGE_NAME_UC>`_.)

   .. _Warning, do not add optional defines for tests/examples to configured header files:

   **Warning, do not add optional defines for tests/examples to configured
   header files:** If this is a test-only and/or example-only dependency then
   please do **not** add a ``#cmakedefine`` to the package's core
   `<packageDir>/cmake/<packageName>_config.h.in`_ file.  Instead, add the
   ``#cmakedefine`` line to a configured header that is only included by
   sources for the tests/examples or just add a define on the compiler command
   line (see the ``DEFINES`` argument to `TRIBITS_ADD_LIBRARY()`_ and
   `TRIBITS_ADD_EXECUTABLE()`_, but see the warning about problems with
   ``ADD_DEFINITIONS()`` in `Miscellaneous Notes (TRIBITS_ADD_LIBRARY())`_).
   We don't want the package's header files to change or libraries to have to
   be rebuilt if tests/examples get enabled or disabled.  Otherwise, the
   `TriBITS CTest/CDash Driver`_ process will result in unnecessary rebuilds
   of software over and over again.

3) **Use the features of the upstream package in the source files of the
   downstream package sources and/or tests/examples:** Usage of the features
   of the upstream package ``<upstreamPackageName>`` in the downstream package
   ``<packageName>`` will typically involve adding ``#include
   <upstreamPackageName>_<fileName>`` in the package's C/C++ source (or
   test/example) files (or the equivalent in Fortran).  If it is an optional
   dependency, then these includes will typically be protected using
   preprocessor ifdefs, for example, as::

     #include "<packageName>_config.h"

     #if HAVE_<PACKAGE_NAME_UC>_<OPTIONAL_DEP_PACKAGE_NAME_UC>
     #  include "<upstreamPackageName>_<fileName>"
     #endif 

4) **For an optional dependency, use CMake IF() statements based on
   ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE_NAME}:** When a package
   ``PACKAGE_NAME`` has an optional dependency on an upstream package
   ``OPTIONAL_DEP_PACKAGE_NAME`` and needs to put in optional logic in a
   CMakeLists.txt file, then the IF() statements should use the variable
   `${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE_NAME}`_ and **not** the
   variable ``${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE_NAME}``.  For
   example, to optionally enable a test that depends on the enable of the
   optional upstream dep package, one would use::

     IF (${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE_NAME})
       TRIBITS_ADD_TEST( ... )
     ENDIF()

  .. ToDo: Find an example to point to in TribitsExampleProject.
  
NOTE: TriBITS will automatically add the include directories for the upstream
package to the compile lines for the downstream package source builds and will
add the libraries for the upstream package to the link lines to the downstream
package library and executable links.  See documentation in the functions
`TRIBITS_ADD_LIBRARY()`_ and `TRIBITS_ADD_EXECUTABLE()`_, and the ``DEPLIBS``
argument to these functions, for more details.


How to add a new TriBITS TPL dependency
---------------------------------------

It is often the case where one will want to add a new dependency for an
existing `downstream`_ package to an existing `upstream`_ `TriBITS TPL`_.
This can either be a required dependency or an optional dependency.  Here, we
will refer to the downstream SE package as ``<packageName>`` with base
directory ``<packageDir>`` and will refer to the upstream TPL as
``<tplName>``.

The process for adding a new dependency to an existing upstream TPL is as
follows:

1) **Add the name of the upstream TPL to the downstream package's
   Dependencies.cmake file:** Add ``<tplName>`` to the call of
   `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_ in the downstream package's
   `<packageDir>/cmake/Dependencies.cmake`_ file.  If this is to be a required
   library dependency, then ``<tplName>`` is added to the
   ``LIB_REQUIRED_TPLs`` argument.  Alternatively, if this is to be an
   optional library dependency, then ``<tplName>`` is added to the
   ``LIB_OPTIONAL_TPL`` argument.  (For example, see the file
   ``packages/Teuchos/cmake/Dependencies.cmake`` file in the
   `ReducedMockTrilinos`_ project.)  If only the test and/or example sources,
   and not the package's core library sources, will have the required or
   optional dependency, then ``<tplName>`` is added to the arguments
   ``TEST_REQUIRED_TPLs`` or ``TEST_OPTIONAL_TPLS``, respectively.

2) **For an optional dependency, add `HAVE_` preprocessor macro to the
   package's configured header file:** If this is an optional dependency,
   typically a C/C++ processor macro will be added to the package's configured
   `<packageDir>/cmake/<packageName>_config.h.in`_ file using the line::

     #cmakedefine HAVE_<PACKAGE_NAME_UC>_<OPTIONAL_DEP_TPL_NAME_UC>

   (see `HAVE_<PACKAGE_NAME_UC>_<OPTIONAL_DEP_TPL_NAME_UC>`_.)

   **WARNING:** If this is a test-only and/or example-only dependency then
   please do **not** add a ``#cmakedefine`` to the package's core
   `<packageDir>/cmake/<packageName>_config.h.in`_ file.  See `Warning, do not
   add optional defines for tests/examples to configured header files`_.

3) **Use the features of the upstream TPL in the source files of the
   downstream package sources and/or tests/examples:** Usage of the features
   of the upstream package ``<tplName>`` in the downstream package
   ``<packageName>`` will typically involve adding ``#include
   <tplName>_<fileName>`` in the package's C/C++ source (or test/example)
   files (or the equivalent in Fortran).  If it is an optional dependency,
   then these includes will typically be protected using preprocessor ifdefs,
   for example, as::

     #include "<packageName>_config.h"

     #if HAVE_<PACKAGE_NAME_UC>_<OPTIONAL_DEP_TPL_NAME_UC>
     #  include "<upstreamPackageName>_<fileName>"
     #endif 

4) **For an optional dependency, use CMake IF() statements based on
   ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL_NAME}:** When a package
   ``PACKAGE_NAME`` has an optional dependency on TPL
   ``OPTIONAL_DEP_TPL_NAME`` and needs to put in optional logic in a
   CMakeLists.txt file, then the IF() statements should use the variable
   `${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL_NAME}`_ and **not** the variable
   ``${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_TPL_NAME}``.  For example, to
   optionally enable a test that depends on the enable of the optional TPL,
   one would use::

     IF (${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL_NAME})
       TRIBITS_ADD_TEST( ... )
     ENDIF()

  .. ToDo: Find an example to point to in TribitsExampleProject.
  
NOTE: TriBITS will automatically add the include directories for the upstream
TPL to the compile lines for the downstream package source builds and will add
the libraries for the upstream TPL to the link lines to the downstream package
library and executable links.  See documentation in the functions
`TRIBITS_ADD_LIBRARY()`_ and `TRIBITS_ADD_EXECUTABLE()`_, and the ``DEPLIBS``
argument to these functions, for more details.

How to tentatively enable a TPL
-------------------------------

A TriBITS package can request the tentative enable of any of its optional TPLs
(see `How to add a new TriBITS TPL dependency`_).  This is done by calling
`TRIBITS_TPL_TENTATIVELY_ENABLE()`_ in the SE package's
`<packageDir>/cmake/Dependencies.cmake`_ file.  For example::

  TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
    ...
    LIB_OPTIONAL_TPLS  SomeTpl
    ...
    )

  TRIBITS_TPL_TENTATIVELY_ENABLE(SomeTpl)

This will result is an attempt to find the components for the TPL ``SomeTpl``.
But if that attempt fails, then the TPL will be disabled and
``${PACKAGE_NAME}_ENABLE_SomeTpl`` will be set to ``OFF``.


How to insert a package into an upstream repo
---------------------------------------------

Sometimes it is desired to insert a package from a downstream VC repo into an
upstream `TriBITS Repository`_ in order for one or more packages in the
upstream repo to define a dependency on that package.  The way this is
supported in TriBITS is to just list the inserted package into the
``PackagesList.cmake`` file of the upstream TriBITS repo after the packages it
depends on and before the packages that will use it then call the
`TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES()`_ function to allow the package to
be missing.  This is demonstrated in `TribitsExampleProject`_ with the package
``InsertedPkg`` which is **not** included in the default
``TribitsExampleProject`` source tree.  The
`TribitsExampleProject`_/``PackagesList.cmake`` file looks like:

.. include:: ../../examples/TribitsExampleProject/PackagesList.cmake
   :literal:

In this example, ``InsertedPkg`` has a required dependency on ``SimpleCxx``
and the SE package ``WithSubpackagesB`` has an optional dependency on
``InsertedPkg``.  Therefore, the inserted package ``InsertedPkg`` has upstream
and downstream dependencies on packages in the ``TribitsExampleProject`` repo.

The function ``TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES()`` tells TriBITS to
treat ``InsertedPkg`` the same as any other package if the directory
``TribitsExampleProject/InsertedPkg`` exists or to completely ignore the
package ``InsertedPkg`` otherwise.  In addition, TriBITS will automatically
disable of all downstream package dependencies for the missing package (and
print a note about the disables).  NOTE: By default TriBITS will silently
ignore missing inserted packages and disable optional support for the missing
package.  To see what packages are missing and being ignored, configure with::

  -D <Project>_WARN_ABOUT_MISSING_EXTERNAL_PACKAGES=TRUE

The way one would set up ``TribitsExampleProject`` to enable ``InsertedPkg``,
if these were in separate VC (e.g. git) repos for example, would be to do::

  $ git clone <some-url-base>/TribitsExampleProject
  $ cd TribitsExampleProject
  $ git clone <some-other-url-base>/ExteranlPkg
  $ echo /InsertedPkg/ >> .git/info/excludes

Then, when you configure ``TribitsExampleProject``, the package
``InsertedPkg`` would automatically appear and could then be enabled or
disabled like any other TriBITS package.  The TriBITS test
``Tribits_TribitsExampleProject_InsertedPkg`` demonstrates this.

Assuming that one would put the (new) external package in a separate VC repo,
one would perform the following steps:

1) Pick a name for the new inserted external package ``<insertedPackageName>``
   (see `Globally unique TriBITS package names`_).  (NOTE: The external
   package may already exist in which case the name would already be
   selected).

2) Create a new VC repo containing the new package ``<insertedPackageName>``
   or add it to an existing extra VC repo. (Or, add the new package to an
   existing downstream repo but don't add it to the ``PackagesList.cmake``
   file to downstream repo.  That would define the package twice!)

3) Clone the downstream VC repo under the upstream TriBITS repository.  (This
   is currently needed since TriBITS only allows package dirs to be contained
   under the repository directory.)

4) Insert the package into the `<repoDir>/PackagesList.cmake`_ file as
   described in `How to add a new TriBITS Package`_ except one must also call
   ``TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(<insertedPackageName>)`` as
   described above.

5) Flesh out the new package and use it in downstream SE packages just like it
   was any other package.  But note that any downstream SE package that has a
   required dependency on ``<insertedPackageName>`` will always be hard
   disabled when the source for ``<insertedPackageName>`` is missing.

6) When configuring and building to get the package working, add
   ``-D<insertedPackageName>_ALLOW_MISSING_EXTERNAL_PACKAGE=FALSE`` so that
   TriBITS will catch mistakes in specifying the package directory.
   Otherwise, to see notes about ignoring missing inserted/external packages,
   set the variable ``-D<Project>_WARN_ABOUT_MISSING_EXTERNAL_PACKAGES=TRUE``
   and TriBITS will print warnings about missing external packages.

How to check for and tweak TriBITS "ENABLE" cache variables
-----------------------------------------------------------

TriBITS defines a number of special ``<XXX>_ENABLE_<YYY>`` variables for
enabling/disabling various entities that allow for a default "undefined" empty
``""`` enable status.  Examples of these special variables include:

* ``${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>`` ((SE) packages)
* ``TPL_ENABLE_<tplName>`` (TPLs)
* ``<TRIBITS_PACKAGE>_ENABLE_<TRIBITS_DEP_PACKAGE_OR_TPL>`` (Optional support
  for a (SE) package or TPL in a downstream package)
* ``<TRIBITS_PACKAGE>_ENABLE_TESTS`` (Package tests)
* ``<TRIBITS_PACKAGE>_ENABLE_EXAMPLES`` (Package examples)
* ``${PROJECT_NAME}_ENABLE_TESTS`` (Tests for explicitly enabled packages)
* ``${PROJECT_NAME}_ENABLE_EXAMPLES`` (Examples for explicitly enabled
  packages)

(see `TriBITS Dependency Handling Behaviors`_).

To check for and tweak these special "ENABLE" variables, perform the
following:

1) To check to see if an ``ENABLE`` variable has been enabled or disabled
   (either explicitly or through auto enable/disable logic), use::

     IF ("${<XXX>_ENABLE_<YYY>}" STREQUAL "")
       # Variable has not been set to 'ON' or 'OFF' yet
       ...
     ENDIF()

   This will work correctly independent of if the cache variable has been
   default defined or not.

2) To tweak the enable/disable of one or more of these variables after user
   input but before auto-enable/disable logic:

  a) To tweak the enables/disables for a TriBITS Repository (i.e. affecting
     all TriBITS projects) add enable/disable code to the file
     `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_.

  b) To tweak the enables/disables for a specific TriBITS Project
     (i.e. affecting only that TriBITS project) add enable/disable code to the
     file `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_.

  For example, one might default disable a package if it has not been
  explicitly enabled (or disabled) in one of these files using logic like::

    IF (NOT ${PROJECT_NAME}_ENABLE_Fortran)
      IF ("${${PROJECT_NAME}_ENABLE_<SomeFortranPackage>}" STREQUAL "")
        MESSAGE("-- " "NOTE: Setting ${PROJECT_NAME}_ENABLE_<SomeFortranPackage>=OFF because"
          "${PROJECT_NAME}_ENABLE_Fortran = '${${PROJECT_NAME}_ENABLE_Fortran}'")
        SET(${PROJECT_NAME}_ENABLE_<SomeFortranPackage> OFF)
      ENDIF()
    ENDIF()

In order to understand the above steps for properly querying and tweaking
these ``ENABLE`` variables, one must understand how TriBITS CMake code defines
and interprets variables of this type.

First, note that most of these particular ``ENABLE`` variables are not
``BOOL`` cache variables but are actually ``STRING`` variables with the
possible values of ``ON``, ``OFF`` and empty ``""`` (see the macro
`SET_CACHE_ON_OFF_EMPTY()`_).  Therefore, just because the value of a
``<XXX>_ENABLE_<YYY>`` variable is defined (e.g. ``IF (DEFINED
<XXX>_ENABLE_<YYY>) ... ENDIF()``) does not mean that it has been set to
``ON`` or ``OFF`` yet (or any non-empty values that evaluates to true or false
in CMake).  To see if an ``ENABLE`` variable is one of these variables, look
in the CMakeCache.txt file for the type.  If the type is ``STRING``, then it
is most likely this type of variable with a default value of empty ``""``.
However, if the cache type is ``BOOL`` then it is likely a standard bool
variable that is not allowed to have a value of empty ``""``.

Second, note that the value of empty ``""`` evaluates to ``FALSE`` in CMake
``IF()`` statements.  Therefore, if one just wants to know if one of these
variables evaluates to true, then just use ``IF (<XXX>_ENABLE_<YYY>)
... ENDIF()``.

Third, note that TriBITS will not define these cache variables until TriBITS
processes the ``Dependencies.cmake`` files on the first configure (see `Full
TriBITS Project Configuration`_).  On future reconfigures, these variables are
all defined (but most will have a default value of empty ``""`` stored in the
cache).

The reason the files ``RepositoryDependenciesSetup.cmake`` and
``ProjectDependenciesSetup.cmake`` are the best places to put in these tweaks
is because, as shown in `Full Processing of TriBITS Project Files`_, they get
processed after all of the user input has been read (in CMake cache variables
set with ``-D<variable>=<value>`` and read in from
`${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE`_ files) but before TriBITS adjusts
the SE package and TPLs enables and disables (see `Package Dependencies and
Enable/Disable Logic`_).  Also, these files get processed in `Reduced Package
Dependency Processing`_ as well so they get processed in all contexts where
enable/disable logic is applied.


How to set up multi-repository support
--------------------------------------

The following steps describe how to set up support for TriBITS project
managing multiple version control and TriBITS repositories by default (see
`Multi-Repository Support`_).

1) Add file `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ and list out
extra repos

  For example, this file would contain something like::

    TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES(
      ExtraRepo1  ""  GIT  git@someurl.com:ExtraRepo1  ""  Continuous
      ExtraRepo2  "ExtraRepo1/ExtraRepos2"  GIT  git@someurl.com:ExtraRepo2
                                                     NOPACKAGES  Continuous
      ExtraRepo3  ""  GIT  git@someurl.com:ExtraRepo3  ""  Nightly
      )

  NOTE: If one will not be using the `checkin-test.py`_ tool, or
  `clone_extra_repos.py`_ tool, or the `TriBITS CTest/CDash Driver`_ system,
  then one can leave the **REPO_VCTYPE** and **REPO_URL** fields empty (see
  `TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES()`_ for details).  (TriBITS Core
  does not have any dependencies on any specific VC tool.  These fields are
  listed here to avoid duplicating the list of repos in another file when
  using these additional TriBITS tools.)

2) Set default values for `${PROJECT_NAME}_EXTRAREPOS_FILE`_ and
``${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE`` (and possibly
``${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES``) in the file
`<projectDir>/ProjectName.cmake`_

  For example, add::

    SET(${PROJECT_NAME}_EXTRAREPOS_FILE  cmake/ExtraRepositoriesList.cmake
      CACHE  FILEPATH  "Set in ProjectName.cmake")
    SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE  Continuous
      CACHE  STRING  "Set in ProjectName.cmake")

  to the `<projectDir>/ProjectName.cmake`_ file.  Otherwise, no extra repos
  will be defined or processed by default when configuring the project.

  And if the project can operate without all of its extra repos, the project
  can set the following default in this file as well with::

    SET(${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES  TRUE
      CACHE  STRING  "Set in ProjectName.cmake")

  Otherwise, all of the extra repos need to be present or the project
  configure will fail.

3) If using git as the VC tool, then set the variable
``${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT`` in the
`<projectDir>/ProjectName.cmake`_ file

  For example::

    SET(${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT  TRUE)

4) If wanting a clone tool with git repos, set up a link to the
`clone_extra_repos.py`_ tool in the base ``<projectDir>/`` directory

  Create a symlink to the script `clone_extra_repos.py`_ in the base project
  repo, for example with::

    cd <projecDir>/
    ln -s cmake/tribits/ci_support/clone_extra_repos.py .
    git add clone_extra_repos.py
    git commit


How to submit testing results to a CDash site
---------------------------------------------

The following steps describe how to submit results to a CDash site using the
`TriBITS CTest/CDash Driver`_ support.

1) Create a CDash project ``<ProjectName>`` on the targeted CDash site.

  To do this, one must have an account on the CDash site and the permissions
  to create a new CDash project.  The name of the project on CDash should
  generally match the name of the TriBITS project ``PROJECT_NAME`` but it does
  not have to.  In fact, one can generally submit to any CDash project with
  any name so creating a new CDash project is actually optional.

  NOTE: For open-source projects, Kitware provides the free CDash site
  my.cdash.org that allows a limited number of submits and data per day.  But
  it should be enough to test out submitting to CDash.

2) Create ``CTestConfig.cmake`` and ``CTestCustom.cmake.in`` files for the
   project.

  * The file `<projectDir>/CTestConfig.cmake`_ can be copied and pasted from
    TribitsExampleProject/CTestConfig.cmake.  To customize for your project,
    you generally just need to update the variables ``CTEST_DROP_SITE``,
    ``CTEST_PROJECT_NAME``, and ``CTEST_DROP_LOCATION``.

  * The file `<projectDir>/cmake/ctest/CTestCustom.cmake.in`_ can be copied
    and pasted from ``TribitsExampleProject/cmake/ctest/CTestCustom.cmake.in``
    and them modified as desired.

3) Test experimental submits with ``make dashboard``.

  Once the CDash project and the `<projectDir>/CTestConfig.cmake`_ and
  `<projectDir>/cmake/ctest/CTestCustom.cmake.in`_ files are created, one
  perform an experimental submission by just configuring the project as normal
  (except configuring additionally with ``-DCTEST_BUILD_FLAGS=-j8`` and
  ``-DCTEST_PARALLEL_LEVEL=8`` to use parallelism in the build and testing in
  the ``ctest -S`` script) and then running the build target::

    make dashboard

  That will configure, build, test, and submit results to the CDash site to
  the ``Experimental`` track of the target CDash project (see `Dashboard
  Submissions`_).

  To work out problems locally without spamming the CDash site, one can run
  with::

    env CTEST_DO_SUBMIT=OFF make dashboard

  To submit to a different CDash site and project, change the cache vars
  ``CTEST_DROP_SITE``, ``CTEST_PROJECT_NAME``, and ``CTEST_DROP_LOCATION`` at
  configure time. For example, to submit TribitsExampleProject results to a
  different CDash site, configure with::

    cmake \
    -DCTEST_DROP_SITE=testing.sandia.gov/cdash \
    -DCTEST_PROJECT_NAME=TribitsExProj \
    -DCTEST_DROP_LOCATION="/submit.php?project=TribitsExProj" \
    [other cmake options] \
    <baseDir>/TribitsExamplProject

  and then run ``make dashboard``.
     
4) Add custom CTest -S driver scripts.

  For driving different builds and tests, one needs to set up one or more
  CTest -S driver scripts.  There are various ways to do this but a simple
  approach that avoids duplication is to first create a file like
  ``TribitsExampleProject/cmake/ctest/TribitsExProjCTestDriver.cmake``:

  .. include:: ../../examples/TribitsExampleProject/cmake/ctest/TribitsExProjCTestDriver.cmake
     :literal:

  and then create a set of CTest -S driver scripts that uses that file.  One
  example is the file
  ``TribitsExampleProject/cmake/ctest/general_gcc/ctest_serial_debug.cmake``:

  .. include:: ../../examples/TribitsExampleProject/cmake/ctest/general_gcc/ctest_serial_debug.cmake
     :literal:

5) Test CTest -S driver scripts

  Once a CTest -S driver script (like the ``ctest_serial_debug.cmake`` example
  shown above) is created, one can test it locally and then test a submit to
  CDash.  To test the exact state of the repository locally, one can create a
  temporary base directory, symbolically link in the local project source
  directory, and then run the CTest -S script by setting
  ``CTEST_DO_SUBMIT=OFF``.  For example, the TribitsExampleProject CTest -S
  script can be run and tested locally by doing::

    $ mkdir MOCK_TRIBITSEXPROJ_SERIAL_DEBUG
    $ cd MOCK_TRIBITSEXPROJ_SERIAL_DEBUG/
    $ ln -s $TRIBITS_DIR/examples/TribitsExampleProject .
    $ env  CTEST_DASHBOARD_ROOT=$PWD \
         CTEST_DROP_SITE=testing.sandia.gov/cdash \
         CTEST_PROJECT_NAME=TribitsExProj \
         CTEST_DROP_LOCATION="/submit.php?project=TribitsExProj" \
         CTEST_TEST_TYPE=Experimental \
         CTEST_DO_SUBMIT=OFF \
         CTEST_DO_UPDATES=OFF \
         CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=TRUE \
         ctest -V -S \
          $TRIBITS_DIR/examples/TribitsExampleProject/cmake/ctest/general_gcc/ctest_serial_debug.cmake \
         &> console.out

  where ``TRIBITS_DIR`` is an env var that points to the location of the
  ``TriBITS/tribits`` directory on the local machine (and the location of the
  CDash site and project is changed, since the free my.cdash.org site can only
  accept as small number of builds each day).

  Once that CTest -S driver script is working correctly without submitting to
  CDash, the above ``ctest -S`` command can be run with ``CTEST_DO_SUBMIT=ON``
  which submit the CDash and then print the location on CDash for the
  submitted configure, build, and test results.

6) Set up automated runs of CTest -S driver scripts

  The custom CTest -S driver scripts created above can be run and used to
  submit to CDash in a variety of ways:

  * Cron jobs can be set up to run them at the same time every day.
 
  * Jenkins jobs can be set up to run them based on various criteria.

  * Travis CI can run them to respond to GitHub pushes.

  * Use the legacy `TriBITS Dashboard Driver`_ system (not recommended).

  The setup of Jenkins, Travis CI and other more sophisticated automated
  testing systems will not be described here.  What will be briefly outlined
  below is the setup using cron jobs on a Linux machine.  That is sufficient
  for most smaller projects and provides tremendous value.

  To set up an automated build using a cron job, one will typically create a
  shell driver script that sets the env and then calls the ``ctest -S
  <script>`` command.  Then one just adds a call to that shell driver script
  using ``crontab -e``.  That is about all there is to it.

  .. ToDo: Add some more details and examples on how to do this, using
  .. TribitsExampleProject implementations.  In particular, you need two
  .. clones of your base git repo.  One to provide the CTest -S script, and
  .. the other cloned and updated by the CTest driver script.

Additional Topics
=================

In this section, a number of miscellaneous topics and TriBITS features are
discussed.  These features and topics are either not considered primary
features of TriBITS (but can be very useful in many situations) or don't neatly
fit into one of the other sections.

TriBITS Repository Contents
---------------------------

The TriBITS git repository is organized as a `TriBITS Project`_, `TriBITS
Repository`_, and `TriBITS Package`_ all in the same base directory.  The base
contents are described in the file::

   TriBITS/README.DIRECTORY_CONTENTS.rst

The part of TriBITS that is installed or snapshotted in contained in the
subdirectory `TriBITS/tribits/`_ and is described in the following section.

.. include:: TriBITS.README.DIRECTORY_CONTENTS.rst

.. _TriBITS/tribits/:

.. include:: ../../README.DIRECTORY_CONTENTS.rst


TriBITS System Project Dependencies
-----------------------------------

The core TriBITS system itself (see ``tribts/core/`` in `TriBITS/tribits/`_)
which is used to configure, built, test, create tarballs, and install software
has no dependencies other than a basic installation of CMake (which typically
includes the executables ``cmake``, ``ctest``, and ``cpack``).  Great effort
has been expended to implement all of this core functionality of TriBITS just
using raw CMake.  That means that anyone who needs to configure, build, and
install software that uses TriBITS just needs a compatible CMake
implementation.  CMake is becoming iniquitous enough that many machines will
already have a current-enough version of CMake installed by default on their
systems and therefore no one will need to download or install any extra
software when building and installing a project that uses TriBITS (assuming
the necessary compilers etc. required by the project are also installed).  If
a current-enough version of CMake is not installed on a given system, it is
easy to download the source code and all it needs is a basic C++ compiler to
build and install.

However, note that a specific TriBITS project is free to use any newer CMake
features it wants and therefore these projects will require newer versions of
CMake than what is required by TriBITS (see discussion of
``CMAKE_MINIMUM_REQUIRED()`` in `<projectDir>/CMakeLists.txt`_).  But also
note that specific TriBITS projects and packages will also require additional
tools like compilers, Python (see `Python Support`_), Perl, or many other such
dependencies.  It is just that TriBITS itself does not require any of these in
order to perform the basic configure, build, test, and install of software.
The goal of TriBITS is not to make the portability of software that uses it
any worse than it already is but instead to make it easier in most cases (that
after all is the whole goal of CMake).

While the TriBITS Core functionality to configure, build, test, and install
software is written using only raw CMake, the more sophisticated development
tools needed to implement the full TriBITS development environment require
Python 2.4 (or higher, but not Python 3.x) (see `Python Support`_).  Python is
needed for tools like `checkin-test.py`_ and `gitdist`_.  In addition, these
python tools are used in `TRIBITS_CTEST_DRIVER()`_ to drive automated testing
and submits to CDash.  Also note that ``git`` is the chosen version control
tool for the TriBITS software development tools and all the VC-related
functionality in TriBITS.  But none of this is required for doing the most
basic building, testing, or installation of a project using TriBITS Core.


Python Support
--------------

TriBITS Core does not require anything other than raw CMake.  However, Python
Utils, TriBITS CI Support, and other extended TriBITS components require
Python.  These extra TriBITS tools only require Python 2.4+.  By default, when
a TriBITS project starts to configure using CMake, it will try to find Python
2.4+ on the system (see `Full Processing of TriBITS Project Files`_).  If
Python is found, it will set the global cache variable ``PYTHON_EXECUTABLE``.
If it is not found, then it will print a warning and ``PYTHON_EXECUTABLE``
will be empty.  With this default behavior, if Python is found, then the
TriBITS project can use it.  Otherwise, it can do without it.

While the default behavior for finding Python described above is useful for
many TriBITS project (such as Trilinos), some TriBITS projects need different
behavior such as:

1. The TriBITS project may not ever use Python so there is no need to look for
   it at all.  In this case, the TriBITS project would set
   `${PROJECT_NAME}_USES_PYTHON`_ to ``FALSE``.

2. Some TriBITS projects require Python and should not even configure if it
   can't be found.  In this case, the TriBITS project would set
   `${PROJECT_NAME}_REQUIRES_PYTHON`_ to ``TRUE``.

3. Some TriBITS projects may require a version of Python more recent than 2.4.
   In this case, the TriBITS project would set `PythonInterp_FIND_VERSION`_ to
   some value higher than ``2.4``.  For example, may newer systems have Python
   2.6.6 or higher versions installed by default and projects developed on
   such a system typically requires this version or higher.


Project-Specific Build Reference
--------------------------------

If a project that uses TriBITS is going to have a significant user base that
will configure, build, and test the project, then having some documentation
that explains how to do this would be useful.  For this purpose, TriBITS
provides a mechanism to quickly create a project-specific build reference
document in restructured text (RST) format and with HTML and LaTeX/PDF
outputs.  This document are generally created in the base project source tree
and given then name ``<Project>BuildReference.[rst,html,pdf]``.  This document
consists of two parts.  One part is a generic template document::

  tribits/doc/TribitsBuildReferenceBody.rst

provided in the TriBITS source tree that uses the place-holder ``<Project>``
for the for the real project name.  The second part is a project-specific
template file::

  <projectDir>/cmake/<Project>BuildReferenceTemplate.rst

which provides the outer RST document (with title, authors, abstract,
introduction, other introductory sections).  From these two files, the
script::

  tribits/doc/build_ref/create-project-build-quickref.py

is used to replace ``<Project>`` in the ``TribitsBuildReferenceBody.rst`` file
with the real project name (read from the project's ``ProjectName.cmake`` file
by default) and then generates the read-only files::

  <projectDir>/
    <Project>BuildReference.rst
    <Project>BuildReference.html
    <Project>BuildReference.pdf

For a simple example of this, see::

  tribits/doc/build_ref/create-build-ref.sh

A project-independent version of this file is provided in the
`TribitsBuildReference`_.[rst,html,pdf] which is referred to many times in this
developers guide.


Project and Repository Versioning and Release Mode
----------------------------------------------------

TriBITS has built-in support for project and repository versioning and
release-mode control.  When the project contains the file
`<projectDir>/Version.cmake`_, it is used to define the project's official
version.  The idea is that when it is time to branch for a release, the only
file that needs to be changed is the file `<projectDir>/Version.cmake`_.

Each TriBITS repository can also contain a `<repoDir>/Version.cmake`_ file
that sets version-related variables which TriBITS packages in that repository
can use to derive development and release version information.  If the TriBITS
repository also contains a `<repoDir>/Copyright.txt`_ file, then the
information in ``<repoDir>/Version.cmake`` and ``<repoDir>/Copyright.txt`` are
used to configure a repository version header file::

  ${${REPOSITORY_NAME}_BINARY_DIR}/${REPOSITORY_NAME}_version.h

The configured header file ``${REPOSITORY_NAME}_version.h`` defines C
pre-processor macros that give the repository version number in several
formats, which allows C/C++ code (or any software that uses the C
preprocessor) to write conditional code like::

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
compilers, and get ready to compile code.  This was mentioned in `Full
Processing of TriBITS Project Files`_.  This is executed by the TriBITS macro
``TRIBITS_SETUP_ENV()``.  Some of the things this macro does are:

.. _Probe and set up the environment:

**Probe and set up the environment:**

* Set ``CMAKE_BUILD_TYPE`` default (if not already set)
* Set up for MPI (MPI compilers, etc.)
* Set up C, C++, and Fortran compilers using ``ENABLE_LANGUAGE(<LANG>)``
* ``INCLUDE(`` `<projectDir>/cmake/ProjectCompilerPostConfig.cmake`_ ``)``
* Find Perl (sets ``PERL_EXECUTABLE``)
* Determine mixed language C/Fortran linking
* Set up C++11 (`${PROJECT_NAME}_ENABLE_CXX11`_)
* Set up OpenMP (with ``FIND_PACKAGE(OpenMP)``)
* Set up optional Windows support
* Find Doxygen (sets ``DOXYGEN_EXECUTABLE``)
* Perform some other configure-time tests (see ``cmake`` configure output)

At the completion of this part of the processing, the TriBITS CMake project is
ready to compile code.  All of the major variables set as part of this process
are printed to the ``cmake`` stdout when the project is configured.


RPATH Handling
--------------

As explained in `Setting install RPATH`_, TriBITS changes the CMake defaults
to write in the RPATH for shared libraries and executables so that they run
right out of the install directory without needing to set paths in the
environment (e.g. ``LD_LIBRARY_PATH``).  However, these defaults can be
changed by changing setting different project defaults for the variables
`${PROJECT_NAME}_SET_INSTALL_RPATH`_ and `CMAKE_INSTALL_RPATH_USE_LINK_PATH`_.
But most projects should likely keep these defaults in place since they make
it so that doing builds and installations on a single machine work correctly
by default out of the box.  For other installation/distribution use cases, the
user is told how to manipulate CMake variables for those cases in `Setting
install RPATH`_.


Configure-time System Tests
---------------------------

CMake has good support for defining configure-time checks of the system to
help in configuring the project.  One can check for whether a header file
exists, if the compiler supports a given data-type or language feature, or
perform almost any other type of check that one can imagine that can be done
using the configured compilers, libraries, system tools, etc.  An example was
given in `TribitsExampleProject`_.  Just follow that example, look at some of
the built-in CMake configure-time test modules, and consult provided on-line
CMake documentation in order to learn how to create a configure-time test for
almost anything.

.. ToDo: Provide more detail how to integrate a configure-time test into a
.. TriBITS project and TriBITS package.


Creating Source Distributions
-----------------------------

The TriBITS system uses CMake's built-in CPack support to create source
distributions in a variety of zipped and tarred formats.  (Note that the term
"source tarball" or just "tarball" may be used below but should be interpreted
as "source distribution".)  TriBITS will automatically add support for CPack
when the variable `${PROJECT_NAME}_ENABLE_CPACK_PACKAGING`_ is set to ``ON``.
The commands for creating a source distribution are described in `Creating a
tarball of the source tree`_ using the built-in ``package_source`` build
target. The value added by TriBITS is that TriBITS will automatically exclude
the source for any defined packages that are not enabled and TriBITS provides
a framework for systematically excluding files and directories from individual
repositories and packages.  In addition, the source for non-enabled
subpackages can also be excluded depending on the value of
`${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION`_.  All of
this allows one to create distributions which only includes subsets of a
larger project (even a single package in some cases).

Unlike other build systems (like autotools), CMake will put **EVERYTHING**
into the source distribution (e.g. tarball) that is sitting in the source tree
by default.  Therefore, setting up for a source distribution usually means
deciding what extra files and directories should be excluded.  Beyond the
directories for non-enabled packages, further files can be selected to be
excluded on a package-by-package based and at the repository level (see
below).

Individual packages can list additional files/directories under the package's
source tree to be excluded from the project's source distribution using a call
to `TRIBITS_EXCLUDE_FILES()`_ in their `<packageDir>/CMakeLists.txt`_ file.
Note that if the package is not enabled, these excludes will never get added!
That is no problem if these excludes only apply to the given package since
TriBITS will add an exclude for the entire package but is a problem if a
package lists excludes for files outside of the package's source tree.

Additional files and entire directories can also be excluded at the repository
level by listing them in the
`<repoDir>/cmake/CallbackDefineRepositoryPackaging.cmake`_ file which
typically just appends the built-in CMake variable
``CPACK_SOURCE_IGNORE_FILES``.  However, if a repository is not processed,
this file is never processed and therefore no files will be excluded from a
repository sitting in the source tree that is not processed (see below).

There are a number of project-level settings that need to be defined and these
are specified in the file
`<projectDir>/cmake/CallbackDefineProjectPackaging.cmake`_.

The `TribitsExampleProject`_ is set up for creating source distributions and
this is demonstrated in one of the tests defined in::

  TriBITS/test/core/ExamplesUnitTests/CMakeLists.txt

There are a few points of caution to note about creating source distributions.

**NOTE:** It is worth stressing again that **EVERY** file that is in the
source tree will be included in the source distribution (tarball) unless there
is an exclude regex matching it appended to the variable
``CPACK_SOURCE_IGNORE_FILES``.  TriBITS can only add excludes for defined
non-enabled packaged.  Every other file listed in the source tree will be
included in the tarball.

**NOTE:** The entries in ``CPACK_SOURCE_IGNORE_FILES`` are interpreted as
**REGULAR EXPRESSIONS** not globs so if you add ``"someFile.*"`` as an
exclude, it will exclude every file in the entire source tree that has
``"someFile"`` in the name!  This is because, in regex terminology, the
trailing ``".*"`` means "match any character zero or more times" and
``"someFile"`` can match anywhere in the file name path.  Also, note that if
you add in an exclude like ``"*.pyc"`` (i.e. trying to exclude all of the
generated Python bite code files) that it will exclude every file that has
``"pyc"`` in the name and **not** just those with the file extension
``"pyc"``.  For example, the exclude ``".pyc"`` would exclude the files
``"puppyc"``, ``"lpycso"``, etc.  If you want to exclude all files with
extension ``"pyc"``, you have to add the exclude regex ``".*[.]pyc$"``!  One's
lack of understanding of this fact will cost someone hours of lost time
debugging what happens when random files are missing when one tries to
configure what is left.  Somethings, what is left will actually configure and
might almost build!

**NOTE:** As warned in `TriBITS Package Core Files`_ and `TriBITS Subpackage
Core Files`_, SE Packages must have directories that are strictly independent
of the directories of other SE packages.  If they don't, then the source
directory for an enabled package will get excluded from the source
distribution if its directory is under the directory of a package that is not
enabled.  For example, if ``PackageA`` is enabled but its package directory
``packageb/packagea/`` is under the package directory ``packageb/`` for the
disabled package ``PackageB``, then every file and directory under
``packageb/`` will be excluded from the source distribution (tarball),
including everything under ``packageb/packagea/``!  It would be too expensive
to put in an automated check for cases like this so package developers should
just take care not to nest the directories of packages inside of each other to
avoid problems like this.

**NOTE:** Extra repositories that are sitting in the source tree but not
processed by TriBITS for some reason (e.g. due to explicitly listing in the
variable `${PROJECT_NAME}_EXTRA_REPOSITORIES`_ only a subset of the repositories
listed in `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ ) will get added
to the source distribution in full by default even though there are no enabled
packages from these repos.  These non-processed repo dirs are like any other
random directory sitting in the source tree, they will get copied over into
the source distribution!.

**NOTE:** When debugging tarball creation problems, always configure with the
variable ``<Project>_DUMP_CPACK_SOURCE_IGNORE_FILES=ON``.  If you don't see a
regex listed for the file or directory you expect to be excluded, then that
file/directory it will be included in the source distribution!


Using Git Bisect with checkin-test.py workflows
-----------------------------------------------

There are cases where a customer will do an update of an upstream project from
a git repo and then find out that some feature or behavior is broken with
respect to their usage.  This can happen even if the upstream project's own
test suite passes all of its tests.  Depending on the situation, there may be
hundreds to thousands of commits between the last known "good" working version
of the code and pulled "bad" version.  To help customers find the first commit
that contains the changes which are causing the breakage, git supplies ``git
bisect``.  This set of commands does a binary search of the commits in the
range ``<good-sha>..<bad_sha>`` and finds the first commit that is "bad" (or a
range of commits which contains the first "bad" commit if commits are skipped,
as described below).

But the ``git bisect`` commands require that all of the commits in the range
``<good-sha>..<bad_sha>`` be complete commits that provide full working
versions of the upstream project.  However, many beginning git developers and
even many experienced developers don't always create individual git commits
that build and pass all of the upstream project's tests and therefore can
create false "bad" commits during the binary search.  This can happen when
developers create intermediate "check-point" commits during the development
process but did not squash the intermediate commits together to create
cohesive working commits.  This can also happen when experienced developers
have a set of what they believe are all working commits but do not actually
test all of the commits to know that they pass *all* of the upstream project's
tests before pushing these commits to the main development branch.  This lack
of detailed testing of each and every individual commit can give rise to false
"bad" commits which will result in ``git bisect`` reporting the wrong first
"bad" commit.

Projects that use the `checkin-test.py`_ tool to push sets of commits to the
main development branch have an advantage in the usage of ``git bisect``.
This is because the default mode of the ``checkin-test.py`` script is to amend
the top commit message with a summary of what was tested and therefore marks a
set of commits that are known to have more complete testing.  For example, the
``checkin-test.py`` tool amends the top commit (after the final pull and
rebase by default) as shown in the following Trilinos commit::

  commit 71ce56bd2d268922fda7b8eca74fad0ffbd7d807
  Author: Roscoe A. Bartlett <bartlettra@ornl.gov>
  Date:   Thu Feb 19 12:04:11 2015 -0500
  
      Define HAVE_TEUCHOSCORE_CXX11 in TeuchosCore_config.h
      
      This makes TeuchosCore a good example for how Trilinos (or any TriBITS)
      subpackages should put in an optional dependency on C++11.
      
      Build/Test Cases Summary
      Enabled Packages: TeuchosCore
      Disabled Packages: [...]
      0) MPI_DEBUG => passed: passed=44,notpassed=0 (2.61 min)
      1) SERIAL_RELEASE => passed: passed=43,notpassed=0 (1.08 min)

Therefore, these special known-tested commits can be flagged by grepping the
``git log -1 HEAD`` output for the string ``"Build/Test Cases Summary"``.  By
bisecting on these commits, one has a lower chance of encountering false "bad"
commits and has a higher chance of finding a smaller range of commits where
the first true "bad" commit might be found.  To aid in performing ``git
bisect`` and only checking ``checkin-test.py``-tested commits, the tool
`is_checkin_tested_commit.py`_ is provided.

To demonstrate how the ``is_checkin_tested_commit.py`` tool can be used with
``git bisect``, suppose that someone writes a customized script
``build_and_test_customer_code.sh`` that will build the upstream project and
the downstream customer's code and then run a set of tests to see if the "bad"
behavior seen by the customer code is the current ``HEAD`` version of the
upstream project.  Assume this script is copied into the upstream project's
local git repo using::

  $ cd <upstream-repo>/
  $ cp ~/build_and_test_customer_code.sh .
  $ cat /build_and_test_customer_code.sh >> .git/info/exclude

Now, one could use ``build_and_test_customer_code.sh`` directly with::

  $ git bisect run ./build_and_test_customer_code.sh

but that would result in testing *all* the commits which may have a high
chance of producing false "bad" commits as described above and fail to
correctly bracket the location of the true first "bad" commit.

So instead, one can write a filtered testing script
``safe_build_and_test_customer_code.sh`` which calls
``is_checkin_tested_commit.py`` and ``build_and_test_customer_code.sh`` as
follows::

  #!/bin/bash
  #
  # Script: safe_build_and_test_customer_code.sh

  $TRIBITS_DIR/ci_support/is_checkin_tested_commit.py
  IS_CHECKIN_TESTED_COMMIT_RTN=$?
  if [ "$IS_CHECKIN_TESTED_COMMIT_RTN" != "0" ] ; then
    exit 125 # Skip the commit because HEAD is not known to be tested!
  fi

  ./build_and_test_customer_code.sh  # Rtn 0 "good", or [1, 124] if "bad"

The above test script ``safe_build_and_test_customer_code.sh`` will skip the
testing of commits that are not marked by the ``checkin-test.py`` tool.

To demonstrate how to use the ``is_checkin_tested_commit.py`` script with
``git bisect``, an example from Trilinos is used below.  (Note that the
current Trilinos public repository may have been filtered so the commit SHA1s
shown below may not match what is in the current Trilinos repository.  But one
can use the commit summary message, author, and author date to find the
updated SHA1s and then to update this example for the current repository.)

Consider a scenario where a customer application updates Trilinos from the
commit::

  d44c17d "Merge branch 'master' of software.sandia.gov:/space/git/Trilinos"
  Author: Roscoe A. Bartlett <xxx@ornl.gov>
  Date:   Tue May 26 12:43:25 2015 -0400

to the commit::

  605b91b "Merge branch 'master' of software.sandia.gov:/git/Trilinos"
  Author: Vitus Leung <xxx@sandia.gov>
  Date:   Tue Sep 29 20:18:54 2015 -0600

and it is found that some critical feature broke or is not behaving acceptably
for the customer code (but all of the tests for Trilinos pass just fine).
This range of commits ``d44c17d..605b91b`` gives 2257 commits to search as
shown by::

  $ cd Trilinos/
  $ git log-oneline d44c17d..605b91b | wc -l
  2257

However, as described above, it is likely that doing ``git bisect`` on that
full set of 2257 commits may result in hitting false "bad" commits and
therefore result in a false bracketing of the first "bad" commit.  This is
where the usage of the ``checkin-test.py`` tool helps which is used by many
(but not currently all) Trilinos developers to push changes to the Trilinos
'master' branch in the current single-branch workflow.  The commits marked
with the ``checkin-test.py`` tool are known (with some caveats mentioned
below) to be working commits and for this the range of commits
``d44c17d..605b91b`` yields 166 commits as shown by::

  $ git log-oneline --grep="Build/Test Cases Summary" d44c17d..605b91b | wc -l
  166

That is an average of 2257/166 = 13.6 commits between commits pushed with the
``checkin-test.py`` tool.  So bisecting on just the commits marked by
``checkin-test.py`` should bound the "bad" commit in a set of 13.6 commits on
average.  Bisecting on this set of 166 commits should likely give no false
“bad” commits, and therefore result in the correct bracketing of the first
"bad" commit.

Using the ``safe_build_and_test_customer_code.sh`` shown above, one would
search for the first bad commit over this range using::

  $ git bisect start 605b91b d44c17d
  $ env DUMMY_TEST_COMMIT_BAD_SHA=83f05e8 \
      time git bisect run ./safe_build_and_test_customer_code.sh

and this would return the range of commits that contains the first "bad"
commit (listed at the end of ``git bisect log`` output, see example below).

To provide a concrete example, suppose the commit that first introduced the
problem in the range of commits ``d44c17d..605b91b`` was::

  83f05e8 "MueLu: stop semi-coarsening if no z-layers are left."
  Author: Tobias Wiesner <tawiesn@sandia.gov>
  Date:   Wed Jul 1 14:54:20 2015 -0600

And instead of using the script ``safe_build_and_test_customer_code.sh``, we
use a dummy driver script ``dummy_test_commit.sh`` to simulate this which is
provided in the set of TriBITS documentation::

  $TRIBITS_DIR/doc/developers_guide/scripts/dummy_test_commit.sh

as:

.. include:: scripts/dummy_test_commit.sh
   :literal:

This driver script allows one to simulate the usage of ``git bisect`` to
understand how it works without having to actually build and test code.  It is
a useful training and experimentation tool.

Using ``git bisect`` (with git version 2.1.0) over the range of commits
``d44c17d..605b91b`` searching for the first "bad" commit is done by running
the commands::

  $ git bisect start 605b91b d44c17d
  $ env DUMMY_TEST_COMMIT_BAD_SHA=83f05e8 \
      time git bisect run \
      $TRIBITS_DIR/doc/developers_guide/scripts/dummy_test_commit.sh \
      &> ../git_bisect_run.log
  $ git2 bisect log &> ../git_bisect_log.log
  $ cat ../git_bisect_log.log | grep "possible first bad commit" | \
      sed "s/possible first bad commit://g"  | sed "s/[a-z0-9]\{30\}\]/]/g"
  $ git bisect reset
  
This set of commands yield the output::

  Bisecting: 1128 revisions left to test after this (roughly 10 steps)
  [9634d462dba77704b598e89ba69ba3ffa5a71471] Revert "Trilinos: remove _def.hpp [...]"
  
  real	1m22.961s
  user	0m57.157s
  sys	3m40.376s
  
  #  [165067ce53] MueLu: SemiCoarsenPFactory. Use strided maps to properly transfer [...]
  #  [ada21a95a9] MueLu: refurbish LineDetectionFactory
  #  [83f05e8970] MueLu: stop semi-coarsening if no z-layers are left.

  Previous HEAD position was 83f05e8... MueLu: stop semi-coarsening if no z-layers are left.
  Switched to branch 'master'

This output shows the dummy bad commit 83f05e8 in a set of just 3 commits,
bounded in the set of commits ``8b79832..165067c``::

  165067c "MueLu: SemiCoarsenPFactory. Use strided maps to properly [...]."
  Author: Tobias Wiesner <tawiesn@sandia.gov>
  Date:   Thu Jul 2 12:11:24 2015 -0600
  
  8b79832 "Ifpack2: RBILUK: adding additional ETI types"
  Author: Jonathan Hu <jhu@sandia.gov>
  Date:   Thu Jul 2 14:17:40 2015 -0700

The set of commits that were actually tested by ``git bisect run <script>`` is
shown by::

  $ cat ../git_bisect_log.log | grep "\(good:\|bad:\)" | sed "s/[a-z0-9]\{30\}\]/]/g"
  # bad: [605b91b012] Merge branch 'master' of software.sandia.gov:/git/Trilinos
  # good: [d44c17d5d2] Merge branch 'master' of software.sandia.gov:/space/git/Trilinos
  # good: [7e13a95774] Ifpack2: If the user does not provide the bandwith of the banded [...]
  # bad: [7335d8bc92] MueLu: fix documentation
  # bad: [9997ecf0ba] Belos::LSQRSolMgr: Fixed bug in setParameters.
  # bad: [b6e0453224] MueLu: add a nightly test for the combination of semicoarsening [...]
  # bad: [165067ce53] MueLu: SemiCoarsenPFactory. Use strided maps to properly [...]
  # good: [3b5453962e] Ifpack2: Nuking the old ETI system
  # good: [8b79832f1d] Ifpack2: RBILUK: adding additional ETI types

This is only 9 commits out of the possible set of 166 ``checkin-test.py``
marked commits which is out of the total set of 2257 possible commits.  With
just 9 build/test cycles, it bounded the first "bad" commit in a set of 3
commits in this case.  And it does not matter how sloppy or broken the
intermediate commits are in Trilinos.  All that matters is the usage of the
``checkin-test.py`` tool (another motivation for the usage of the
``checkin-test.py`` tool, see `Pre-push Testing using checkin-test.py`_ for
others as well).
 
Note that above, we grep the output from ``git bisect log`` for the set of
possible "bad" commits instead of just looking at the output from the ``git
bisect run <script>`` command (which also lists the set of possible "bad"
commits).  This is because the direct output from the ``git bisect run
<script>`` command (shown in the log file ``git_bisect_run.log``) shows the
set of possible bad commits at the end of the output but they are unsorted and
give no other git commit information::

  There are only 'skip'ped commits left to test.
  The first bad commit could be any of:
  83f05e89706590c4b384dd191f51ef4ab00ce9bb
  ada21a95a991cd238581e5a6a96800d209a57924
  165067ce538af2cd0bd403e2664171726ec86f3f
  We cannot bisect more!
  bisect run cannot continue any more

The problem with unsorted commits is that it is difficult to use an unsorted
set to do further bisection.  However, the output of the set of commits from
``git bisect log`` is sorted and also shows the commit summary message and
therefore is much more useful.  (But note that older versions of git don’t
show this set of commits at the end of ``git bisect log`` so make sure and use
an updated version of git, preferably >= 2.1.0.)

Now that one has the range of possible "bad" commits (just 3 in this example)
doing a further manual bisection or even manual inspection of these commits
may be enough to find the change that is causing the problem for the
downstream customer application.

Without the usage of the ``checkin-test.py`` tool, one would not have an
automated way to ensure that ``git bisect`` avoids false "bad" commits.  This
allows for less experienced developers to create commits and push to the main
development branch but still ensure effective usage of ``git bisect``.  (This
is another example where automated tools in TriBITS help to overcome lacking
developer experience and discipline.)


Multi-Repository Almost Continuous Integration
----------------------------------------------

The `checkin-test.py`_ tool can be used to the implement staged integration of
the various repositories in a multi-repo TriBITS project (see
`Multi-Repository Support`_) .  This is referred to here as Almost Continuous
Integration (ACI).  The basic concept of Almost Continuous Integration (ACI)
is defined and described in the paper [`Integration Strategies for CSE,
2009`_].

This topic is broken down into the following subsections:

* `ACI Multi-Git/TriBITS Repo Integration Example`_
* `ACI Local Sync Git Repo Setup`_
* `ACI Integration Build Directory Setup`_
* `ACI Sync Driver Script`_
* `ACI Cron Job Setup`_
* `Addressing ACI Failures and Summary`_


ACI Introduction
++++++++++++++++

The TriBITS system allows for setting up composite meta-builds of large
collections of software pulled in from many different git/TriBITS code
repositories as described in the section `Multi-Repository Support`_.  The
`checkin-test.py`_ tool is a key tool to enable the testing of a set of
packages in different git/TriBITS repos before pushing to remote tracking
branches for the set of git repos; all in one robust command invocation.

While the ``checkin-test.py`` tool was originally designed and its default
behavior is to test a set of local commits created by a developer before
pushing changes to one or more (public) git repos, it can also be used to set
up an Almost Continuous Integration (ACI) process to keep these various
git/TriBITS repos in sync thereby integrating the work of various disconnected
development teams and projects.  To use the ``checkin-test.py`` tool for ACI
requires some setup and changing what the tool does a little by passing in
additional options that a regular developer typically never uses.

The following subsections describe how to use the `checkin-test.py`_ tool to
implement an ACI process for a given set of git/TriBITS repositories and also
provides a little background and context behind ACI.


ACI Multi-Git/TriBITS Repo Integration Example
++++++++++++++++++++++++++++++++++++++++++++++

In order to set up the context for the ACI process, consider the following
simple TriBITS project with two extra repositories::

  BaseProj/
      ExtraRepo1
      ExtraRepo2

Here, ``BaseProj`` is the base TriBITS project/repository and ``ExtraRepo1``
and ``ExtraRepo2`` are extra repositories that supply additional TriBITS
packages that are appended to the TriBITS packages defined in ``BaseProj``
(see `<projectDir>/cmake/ExtraRepositoriesList.cmake`_).  Also, assume that
``BaseProj``, ``ExtraRepo1``, and ``ExtraRepo2`` are developed by three
different development teams that all have different funding sources and
different priorities so they tend not to work closely together or consider the
other efforts too much when developing their software.  However, in this
example, there is great value in combining all of this software into a single
integrated TriBITS meta-project.  This combined meta-build is driven by a 4th
integration/development team.  In this case, the core developers for each of
these three different git/TriBITS repos do not test compatibility with the
other git/TriBITS repos when pushing commits to their own git/TriBITS repos.
This gives three different git repos on three different machines:

* ``BaseProj`` main repo:  Pushed to by the core BaseProj team::

    url1.gov:/git/BaseProj

* ``ExtraRepo1`` main repo:  Pushed to by the core ExtraRepo1 team::

    url2.gov:/git/ExtraRepo1

* ``ExtraRepo2`` main repo:  Pushed to by the core ExtraRepo2 team::

    url3.gov:/git/ExtraRepo2

Because of the independent development processes of these three teams, unless
these development teams maintain 100% backward compatibility w.r.t. the
interfaces and behavior of the combined software, one cannot at any time pull
the code from these three different git repos and expect to be able to
successfully build all of the code and have all of the tests pass.  Therefore,
how does the 4th integration team expect to be able to build, test, and
possibly extend the combined software?  In this case, the integration team
would set up their own clones of all three git/TriBITS repos on their own
machine such as:

**Integration project mirrored git repos**::

  url4.gov:/git/BaseProj
  url4.gov:/git/ExtraRepo1
  url4.gov:/git/ExtraRepo2

Once an initial collaboration effort between the integration team and the
three other development teams is able to get a version of all three
git/TriBITS repos to work correctly in the combined meta-project, these
versions (assume the ``master`` branches) would be pushed to the git repos on
the git integration server ``url4.gov``.  The state where the TriBITS packages
in the three different git/TriBITS repos in the ``master`` branch on
``url4.gov`` all work together correctly constitutes the initial condition for
the ACI process described below.  From that initial condition, the ACI
processes ensures that updates the ``master`` branches for the git/TriBITS
repos on ``url4.gov`` do not break any builds or tests of the integrated
software.

In order to describe how to set up an ACI process using the
``checkin-test.py`` tool, the following subsections will focus on the update
of the git/TriBITS ``ExtraRepo1`` repo keeping the other two git/TriBITS repos
``BaseProj`` and ``ExtraRepo2`` constant as the ACI use case.


ACI Local Sync Git Repo Setup
+++++++++++++++++++++++++++++

In order to set up an ACI process for the multi-git/TriBITS repo example
outlined above, first local repos are created by cloning the repos on the
integration server url4.gov as follows (all of which become 'origin')::

  $ cd $SYNC_BASE_DIR
  $ git clone url4.gov:/git/BaseProj
  $ cd BaseProj
  $ git clone url4.gov:/git/ExtraRepo1
  $ git clone url4.gov:/git/ExtraRepo2

where, ``SYNC_BASE_DIR=~/sync_base_dir`` for example, which must already be
created.

Next, one defines a remote to pull changes for the ``ExtraRepo1`` from the
main development repo:

  $ cd $SYNC_BASE_DIR/BaseProj/ExtraRepo1
  $ git remote add public url2.gov:/git/ExtraRepo1

Here, one should pick a name for the remote repo for ``ExtraRepo1`` that is
most descriptive for that particular situation.  In this case, the name
``public`` is chosen to signify the main public development repo.

This gives the remotes::

  $ cd $SYNC_BASE_DIR/BaseProj
  $ gitdist remote -v | grep -v push | grep -v "^$"
  *** Base Git Repo: BaseProj
  origin	        url4.gov:/git/BaseProj (fetch)
  *** Git Repo: ExtraRepo1
  origin	        url4.gov:/git/ExtraRepo1 (fetch)
  public		url2.gov:/git/ExtraRepo1 (fetch)
  *** Git Repo: ExtraRepo2
  origin	        url4.gov:/git/ExtraRepo2 (fetch)

The remote ``public`` is used by the ``checkin-test.py`` wrapper script (see
`ACI Sync Driver Script`_ below) to pull and merge in additional changes that
will be tested and pushed to the 'origin' repos on ``url4.gov``.  In this
case, the ``ExtraRepo1`` remote ``public`` will result in updates being pulled
from the main development repo on ``url2.gov``, thereby facilitating the
update of ``ExtraRepo1`` in the integrated meta-project.


ACI Integration Build Directory Setup
+++++++++++++++++++++++++++++++++++++

After the git repos are cloned and the remotes are set up as described above,
a build base directory is set up as::

  $ cd $SYNC_BASE_DIR
  $ mkdir BUILDS
  $ mkdir BUILDS/CHECKIN

An ACI wrapper script for ``checkin-test.py`` is created to drive the syncing
process.  It is assumed that this script would be called only once a day and
not continuously in a loop (but that is possible as well but is not documented
here).

NOTE: Other build directory structures are possible, it all depends how one
writes the ``checkin-test.py`` wrapper scripts but the above directory
structure is fairly standard in the usage of the ``checkin-test.py`` script.


ACI Sync Driver Script
++++++++++++++++++++++

The sync driver script for this example should be called something like
``sync_ExtraRepo1.sh``, placed under version control, and would look something
like::

  #!/bin/bash -e
  
  # Set up the environment (i.e. PATH; needed for cron jobs)
  ...

  SYNC_BASE_DIR=~/sync_base_dir
  CHECKIN_TEST_WRAPPER=$SYNC_BASE_DIR/BaseProj/sampleScripts/checkin-test-foo.sh

  cd $SYNC_BASE_DIR/BUILDS/CHECKIN

  $CHECKIN_TEST_WRAPPER \
    --extra-pull-from=ExtraRepo1:public:master \
    --abort-gracefully-if-no-changes-to-push \
    --enable-extra-packages=Package1A \
    --send-build-case-email=only-on-failure \
    --send-email-to=base-proj-integrators@url4.gov \
    --send-email-to-on-push=base-proj-integrators@url4.gov \
    --no-append-test-results --no-rebase \
    --do-all --push \
    -j16 \
    --wipe-clean \
    "$@"

NOTE, in the above example ``sync_ExtraRepo1.sh`` script, the variable
``CHECKIN_TEST_WRAPPER`` is set to a wrapper script::

   BaseProj/sampleScripts/checkin-test-foo.sh

which would be set up to call the project's `checkin-test.py`_ tool with
configure options for the specific machine.  The location and the nature of
the wrapper script will vary from project to project and machine to machine.
In some simple cases, ``CHECKIN_TEST_WRAPPER`` might just be set to be the raw
machine-independent ``checkin-test.py`` tool for the project.

A description of each option passed into this invocation of the
`checkin-test.py`_ tool is given below (see `checkin-test.py --help`_ for more
details):

  ``--extra-pull-from=ExtraRepo1:public:master``
  
    This option instructs the ``checkin-test.py`` tool to pull and merge in
    commits that define the integration.  One could do the pull(s) manually of
    doing so has the disadvantage that if they fail for some reason, they will
    not be seen by the ``checkin-test.py`` tool and no notification email
    would go out.
  
  ``--abort-gracefully-if-no-changes-to-push``
  
    The option ``--abort-gracefully-if-no-changes-to-push`` makes the
    ``checkin-test.py`` tool gracefully terminate without sending out any
    emails if after all the pulls, there are no local changes to push to the
    'origin' repos.  This can happen, for example, if no commits were pushed
    to the main development git repo for ``ExtraRepo1`` at
    ``url2.gov:/git/ExtraRepo1`` since the last time this sync process was
    run.  This avoids getting confusing and annoying emails like ``"PUSH
    FAILED"``.  The reason this option is not generally needed for local
    developer usage of the ``checkin-test.py`` tool is that in general a
    developer will not run the ``checkin-test.py`` tool with ``--push`` unless
    they have made local changes; it just does not make any sense at all to do
    that and if they do by accident, they should get an error email.  However,
    for an automated ACI sync process, there is no easy way to know a-priori
    if changes need to be synced so the script supports this option to deal
    with that case gracefully.

  ``--enable-extra-packages=Package1A``

    This option should be set if one wants to ensure that all commits get
    synced, even when these changes don't impact the build or the tests of the
    project.  If not setting ``--enable-extra-packages=<some-package>`` , then
    the ``checkin-test.py`` tool will only decide on its own what packages to
    test just based on what packages have changed files in the ``ExtraRepo1``
    repo and if no modified files map to a package, then no packages will be
    auto-enabled and therefore no packages will be enabled at all.  For
    example, if a top-level README file in the base ``ExtraRepo1`` repo gets
    modified that does not sit under a package directory, then the automatic
    logic in the checkin-test.py tool will not trigger a package enable. In
    that case, no configure, build, testing, or push will take place (must run
    at least some tests in order to assume it is safe to push) and therefore
    the sync will not occur.  Therefore, if one wants to ensure that every
    commit gets safely synced over on every invocation, then the safest way to
    that is to always enable at least one or more packages by specify
    ``--enable-extra-packages=<pkg0>,<pkg1>``.  **WARNING:** it is not
    advisable to manually set ``--enable-packages=<package-list>`` since it
    turns off the auto-enable logic for changed files.  This is because if
    there are changes to other packages, then these packages will not get
    enabled and not get tested, which could break the global build and tests.
    Also, this is fragile if new packages are added to ``ExtraRepo1`` later
    that are not listed in ``--enable-packages=<pkg0>,<pkg1>,...`` as they
    will not be included in the testing.  Also, if someone makes local commits
    in other local git repos before running the sync script again, then these
    packages will not get enabled and tested.  Therefore, in general, don't
    set ``--enable-packages=<pkg0>,<pkg1>,...`` in a sync script, only set
    ``--enable-extra-packages=<pkg0>,<pkg1>,...`` to be robust and safe.

  ``--send-build-case-email=only-on-failure``

    This makes the checkin-test.py tool skip sending email about a build case
    (e.g. ``MPI_DEBUG``) unless it fails.  That way, if everything passes,
    then only a final ``DID PUSH`` email will go out.  But if a build case
    does fail (i.e. configure, build, or tests fail), then an "early warning"
    email will still go out.  However, if one wants to never get the early
    build-case emails, one can turn this off by setting
    ``--send-build-case-email=never``.
 
  ``--send-email-to=base-proj-integrators@url4.gov``
 
    The results of the builds will be sent this email address.  If you only
    want an email sent when a push actually happens, you can set
    ``--send-email-to=''`` and rely on ``--send-email-to-on-push``.
  
  ``--send-email-to-on-push=base-proj-integrators@url4.gov``
  
    A confirmation and summary email will be sent to this address if the push
    happens.  This can be a different email address than set by the
    ``--send-email-to`` option.  It is highly recommended that a mail list be
    used for this email address since this will be the only more permanent
    logging of the ACI process.
  
  ``--no-append-test-results --no-rebase``
  
    These options are needed to stop the ``checkin-test.py`` tool from
    modifying the commits being tested and pushed from one public git repo to
    another.  The option ``--no-append-test-results`` is needed to instruct
    the ``checkin-test.py`` tool to **NOT** amend the last commit with the
    test results.  The option ``--no-rebase`` is needed to avoid rebasing the
    new commits pulled.  While the default behavior of the ``checkin-test.py``
    tool is to amend the last commit message and rebase the local commits
    (which is considered best practice when testing local commits), this is a
    very bad thing to do when a ACI sync server is only testing and moving
    commits between public repos.  Amending the last commit would change the
    SHA1 of the commit (just as a rebase would) and would fork the history and
    mess up a number of workflows that otherwise should work smoothly.  Since
    an email logging what was tested will go out if a push happens due to the
    ``--send-email-to-on-push`` argument, there is no value in appending the
    test results to the last commit pulled and merged (which will generally
    not be a merge commit but just a fast-forward).  There are cases, however,
    where appending the test results in an ACI process might be acceptable but
    they are not discussed here.
  
  ``--do-all --push -j16``
  
    These are standard options that always get used when invoking the
    ``checkin-test.py`` tool and need no further explanation.
  
  ``--wipe-clean``
  
    This option is added if you want to make the sync server more robust to
    changes that might require a clean configure from script.  If you care
    more about using less computer resources and testing that rebuilds work
    smoothly, remove this option.

The sync script can be created and tested locally to ensure that it works
correctly first, before setting it as a cron job as described next.  Also, the
sync script should be version controlled in one of the project's git repos.
This ensures that changes to script pushed to the repos will get invoked
automatically when they are pulled (but only on the second invocation of the
script).  If a change to script is critical in order to do the pull, then one
must manually pull the updated commit to the local sync repo.

Note, if using this in a continuous sync server that runs many times in a day
in a loop, you also want to set the option
``--abort-gracefully-if-no-changes-pulled`` in addition to the option
``--abort-gracefully-if-no-changes-to-push``.  That is because if the updated
repos are in a broken state such that there are always local changes at every
CI iteration (because they have not been pushed to origin), you don't want to
do a new CI build unless something has changed that would otherwise perhaps
make the error go away.  That allows the CI server to sit ready to try out any
change that gets pulled that might allow the integrated build to work and then
push the updates.


ACI Cron Job Setup
++++++++++++++++++

Once the sync script ``sync_ExtraRepo1.sh`` has been locally tested, then it
should be committed to a version control git repo and then run automatically
as a cron job.  For example, the cron script shown below would fire off the
daily ACI process at 8pm local time every night::

  # ----------------- minute (0 - 59)
  # |  -------------- hour (0 - 23)
  # |  |  ----------- day of month (1 - 31)
  # |  |  |  -------- month (1 - 12)
  # |  |  |  |  ----- day of week (0 - 7) (Sunday=0 or 7)
  # |  |  |  |  |
  # *  *  *  *  *  command to be executed
   00 20  *  *  *  ~/sync_base_dir/sync_ExtraRepo1.sh &> ~/sync_base_dir/sync_ExtraRepo1.out
  
In the above crontab file (set with ``'crontab -e'`` or ``'crontab
my-crontab-file.txt'``), the script::

  ~/sync_base_dir/sync_ExtraRepo1.sh
 
is assumed to be a soft symbolic link to some version controlled copy of the
ACI sync script.  For example, it might make sense for this script to be
version controlled in the ``BaseProj`` repo and therefore the symbolic link
would be created with something like::

  $ cd ~/sync_base_dir/
  $ ln -s BaseProj/sampleScripts/sync_ExtraRepo1.sh .

Such a setup would ensure that sync scripts would always be up-to-date due to
the git pulls part of the ACI process.


Addressing ACI Failures and Summary
+++++++++++++++++++++++++++++++++++

After the above cron job starts running (setup described above), the
``checkin-test.py`` tool will send out emails to the email addresses passed
into the underlying ``checkin-test.py`` tool.  If the emails report an update,
configure, build, or test failure, then someone will need to log onto the
machine where the ACI sync server is running and investigate what went wrong,
just like they would if they were running the ``checkin-test.py`` tool for
testing locally modified changes before pushing.

In the above example, only a single git/TriBITS repo is integrated in this ACI
sync scripts.  For a complete system, other ACI sync scripts would be written
to sync the two other git/TriBITS repos in order to maintain some
independence.  Or, a single ACI sync script that tries to update all three
git/TriBITS repos at could would be written and used.  The pattern of
integrations chosen will depend many different factors and these patterns can
change over time according to current needs and circumstances.

In summary, the `checkin-test.py`_ tool can be used to set up robust and
effective Almost Continuous Integration (ACI) sync servers that can be used to
integrate large collections of software in logical configurations at any
logical frequency.  Such an approach, together with the practice of `Regulated
Backward Compatibility and Deprecated Code`_, can allow very large collections
of software to be kept integrated in short time intervals using ACI.

Post-Push CI and Nightly Testing using checkin-test.py
------------------------------------------------------

While the post-push CI and Nightly testing processes using ``ctest -S``
scripts using `TRIBITS_CTEST_DRIVER()`_ (see `TriBITS CTest/CDash Driver`_)
which posts results to a CDash server (see `TriBITS CDash Customizations`_) is
a very attractive testing system with many advantages, setting up a CDash
server can be a bit difficult and a CDash server can require non-trivial
storage and CPU resources (due to the MySQL DB of test results) and requires
some amount of extra maintenance.  As an intermediate approach, one can
consider just using the project's `checkin-test.py`_ tool to implement basic
post-push CI and/or Nightly testing servers using simple cron jobs and some
other helper scripts.  The ``checkin-test.py`` tool will robustly pull new
commits, configure the project, build, run tests, and send out emails with
results and pass/fail.  A bunch of builds can be run at once using multiple
builds specified in the ``--default-builds``, ``--st-extra-builds``, and
``--extra-build`` arguments, or different invocations of the
``checkin-test.py`` tool can be run in parallel for better machine
utilization.

What one gives up with this approach over the full-blow CTest/CDash
implementation is:

1) Test results over multiple builds/days are not archived in a database for
   later retrieval and query.  Only the latest build's detailed results are
   accessible.

2) Detailed test results (pass or fail) cannot be easily looked up on a
   website.  Only detailed results from the last build can be found and one
   has to go to the machine where the build was run from to see the results.

The process for setting up a basic nightly build using ``checkin-test.py`` as
a cron job is a subset of the steps needed to set up a driver and cron job for
the more complex ACI process described in `ACI Sync Driver Script`_ and `ACI
Cron Job Setup`_.  Setting up a CI server using ``checkin-test.py`` is a
little more involved than setting up a nightly build (and less desirable
because it blows away old CI iteration results pretty fast) but can be done
using existing tools.

.. ToDo: Show how to set up a cron job using checkin-test.py.

For smaller, private, less-critical projects with few developers, setting up
CI and Nightly testing using ``checkin-test.py`` may be quite adequate.  In
fact, post-push testing processes implemented with ``checkin-test.py`` are
much more solid and feature-full than have been employed in many software
projects that we have seen over the years that were larger, more public, had
many developers, and were quite important to many users and development
teams.


TriBITS Dashboard Driver
------------------------

TriBITS also includes a system based on CMake/CTest/CDash to drive the builds
of a TriBITS project and post results to another CDash project.  This system
is contained under the directory::

  tribits/ctest/tdd/

If the TriBITS project name is ``<projectName>``, the TDD driver CDash project
is typically called ``<projectName>Driver``.  Using CTest to drive the Nightly
builds of a TriBITS project makes sense because CTest can run different builds
in parallel, can time-out builds that are taking too long, and will report
results to a dashboard and submit notification emails when things fail.
However, this is the most confusing and immature part of the TriBITS system.
The `TriBITS CTest/CDash Driver`_ system using the `TRIBITS_CTEST_DRIVER()`_
can be used without this TriBITS Dashboard Driver (TDD) system.

However, this TriBITS subsystem is not well tested with automated tests, is
difficult to extend and test manually, and has other problems.  Therefore, it
is not recommended that projects adopt the usage of this subsystem.  A simple
set of cron jobs or a tool like Jenkins is likely a better option (if done
carefully).

.. ToDo: Import and format the contents of tribits/ctest/README related to
.. this topic into this section.


Regulated Backward Compatibility and Deprecated Code
----------------------------------------------------

The motivation and ideas behind `Regulated Backward Compatibility` and
deprecated code are described in the `TriBITS Lifecycle Model`_ document.
Here, the details of the implementation in TriBITS are given and how
transitions between non-backward compatible major versions is accomplished.

This section describes the process for deprecating code within a major
backward-compatible version sequence and finally removing deprecated code when
transitioning to a new major version X to X+1 for the semantic versioning
numbering scheme X.Y.Z.  This information is given as a process with different
phases.

The processing for managing deprecated code is as follows and more details are
given below:

* `Setting up support for deprecated code handling`_

* `Deprecating code but maintaining backward compatibility`_

  * `Deprecating C/C++ classes, structs, functions, typdefs, etc.`_
  * `Deprecating preprocessor macros`_
  * `Deprecating an entire header and/or source file`_

* `Hiding deprecated code to certify and facilitate later removal`_

  * `Hiding C/C++ entities`_
  * `Hiding entire deprecated header and source files`_

* `Physically removing deprecated code`_

  * `Removing entire deprecated header and source files`_
  * `Removing deprecated code from remaining files`_


Setting up support for deprecated code handling
+++++++++++++++++++++++++++++++++++++++++++++++

Setting up support for managing deprecated code in a TriBITS package requires
just two simple changes to the TriBITS-related files in a package.  First, the
top-level package `<packageDir>/CMakeLists.txt`_ file needs to have a call
to::

  TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

Second, the package's configure header::

  <packageDir>/cmake/<PACKAGE_UCNAME>_config.h.in

file needs to have::

  @<PACKAGE_UCNAME>_DEPRECATED_DECLARATIONS@

where ``<PACKAGE_UCNAME>`` is the upper-case name of the TriBITS package.

That is all that is needed to provide support for TriBITS deprecated code.

When a TriBITS project is configured with::

   -D<Project>_SHOW_DEPRECATED_WARNINGS=ON

by default, all packages will show deprecated warnings.  These deprecated
warnings can be turned on and off on a package-by-package basis by setting::

   -D<PackageName>_SHOW_DEPRECATED_WARNINGS=[ON|OFF]

This gives developers and users a little extra control over what deprecated
warnings are shown.

In addition, deprecated code can be hidden from the build to help certify that
downstream code is clean by configuring with::

   -D<Project>_HIDE_DEPRECATED_CODE=ON

This will remove the deprecated code from the build and install (see details
below) so that other software in the TriBITS project can be shown to build
clean without deprecated code and so that outside client code can be shown to
be clean of deprecated code.

As with deprecated warnings, showing or hiding deprecated code can be
controlled on a package-by-package basis by setting::

   -D<PackageName>_HIDE_DEPRECATED_CODE=[ON|OFF]

In this case, hiding deprecated code on a package-by-package basis may not
work because deprecated code in a downstream package might rely on deprecated
code in an upstream package (which might have its deprecated code hidden).


Deprecating code but maintaining backward compatibility
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

One of the most important aspects of the `TriBITS Lifecycle Model`_ for
later-stage Production Growth and Production Maintenance code is to provide
backward compatibility between a continuous stream of versions of the software
within a major version number X in the version numbering scheme X.Y.Z.  In all
cases, if a piece of client code builds and runs correctly against version
X0.Y0.Z0, it should also build, without modification, against versions
X0.Y1,Z1 for all Y1 >= Y0 and all Z1 and up to (but not including) X0+1.0.0.
There are many types of constructs that one will want to deprecate and
therefore later remove.  When deprecating code, one wants to give users
compile-time warnings of the usage of deprecated features so that they know
what they need to remove.  One also wants to allow them to certify that their
code is free of deprecated warnings by hiding deprecated code.  Below, the
different types of entities that one wants to deprecate and how to support
hiding code (which also facilitates removing it later) are described.


Deprecating C/C++ classes, structs, functions, typdefs, etc.
............................................................

To deprecate standard C/C++ constructs, one can just use the standard TriBITS
compile-time macro ``<PACKAGE_UCNAME>_DEPRECATED`` which is properly ifdefed
by the TriBITS system to add a GCC/Intel deprecated attribute or not.  For
example, one would deprecate standard C/C++ constructs for the package
``SomePackage`` with::

  // Deprecate a class (or struct)
  class SOMEPACKAGE_DEPRECATED SomeClass { ... };

  // Deprecate a function
  SOMEPACKAGE_DEPRECATED int someFunc(...);

  // Deprecate a typedef
  SOMEPACKAGE_DEPRECATED typedef someTypeDef int;

The GCC (version 3.1 and newer) and Intel C and C++ compilers both support
adding extra attributes including the ``__deprecated__`` attribute.  When this
attribute is applied to a given C/C++ entity, it produces a compiler warning
that can be searched for in the compiler output and elevated to an error (when
``-Werror`` is also passed to the compiler).

In addition to the basic deprecated warning, one can also add an optional
deprecated warning message using the macro
``<PACKAGE_UCNAME>_DEPRECATED_MSG()``.  For example, if a new function is
replacing an old function, one might use::

  // Deprecated old unsafe function taking raw pointers
  SOMEPACKAGE_DEPRECATED_MSG(
    "Please use the safe someFunc(const Ptr<const std::string>&) instead!")
  void someFunc(const std::string *str);

  // New version take does not take raw pointers (and is therefore safer)
  void someFunc(const Teuchos::Ptr<const std::string> &str);

Then, if user code calls the version ``someFunc(const std::string*)`` they
will get the string::

  "Please use the safe someFunc(const Ptr<const std::string>&) instead!"

printed as well by the compiler.  Note that the custom message will only be
printed for GCC versions 4.5 and newer.  If an older version of GCC is used,
then the message string is ignored.


Deprecating preprocessor macros
...............................

A C/C++ preprocessor macro is not an entity seen by the C/C++ compiler and
therefore cannot directly take advantage of a feature such as the
``__deprecated__`` attribute of the GCC/Intel compilers.  However, in some
cases, for function-like macros can such as::

  // The definition of the macro
  #define SOME_OLD_FUNC_MACRO(ARG1, ARG2, ...) ...
  ...
  // A use of the macro
  SOME_OLD_FUNC_MACRO(a, b, ...)

there is a strategy where one can define the macro to call a dummy deprecated
function such as with::

  SOMEPACKAGE_DEPRECATED
  inline void SOME_OLD_FUNC_MACRO_is_deprecated()
  {}

  // The definition of the macro
  #define SOME_OLD_FUNC_MACRO(ARG1, ARG2, ...) \
    { \
      SOME_OLD_FUNC_MACRO_is_deprecated(); \
      ... \
    }

  ...

  // A use of the macro
  SOME_OLD_FUNC_MACRO(a, b, ...)

In the above example, client code calling ``SOME_OLD_FUNC_MACRO()`` will now
result in a deprecated compiler warning which should make it clear what is
being deprecated.

Note that this approach will not work in cases were a macro is only providing
a single value such as a constant (but one should not use macros for providing
constants anyway).


Deprecating an entire header and/or source file
...............................................

There are times when one wants to deprecate an entire set of files and all of
the contents in those files.  In addition to deprecating the contents of the
files one will want to deprecate the entire file as well.  There are a few
steps to this.  First, one wants to put a warning in the file such as::

  #ifdef __GNUC__
  #  warning "This header <THIS_HEADER> is deprecated!  Use <NEW_HEADER> instead!"
  #endif

The above ``ifdef`` on ``__GNUC__`` is needed in order avoid the non-standard
``#warning`` preprocessor directive on non-compliant compilers (but should
work on all later version GCC and Intel compilers).


Hiding deprecated code to certify and facilitate later removal
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In addition to adding deprecation warnings at preprocessing or compile-time,
it is also highly desirable to allow the deprecated code to be removed from
the build to help certify that client code indeed no longer needs the
deprecated code.  The following subsections describe how to hide deprecated
code from existing files and how to hide deprecated files entirely.


Hiding C/C++ entities
.....................

In the case when various C/C++ entities will be removed from an existing file,
but the file will remain, then the deprecated code can be ifdefed out, for the
package ``SomePackage`` for example, using::

  #ifndef SOMEPACKAGE_HIDE_DEPRECATED_CODE
    // Deprecate a class (or struct)
    class SOMEPACKAGE_DEPRECATED SomeClass { ... };
    ...
  #endif /* SOMEPACKAGE_HIDE_DEPRECATED_CODE */

In this way, when the CMake variable ``SomePackae_HIDE_DEPRECATED_CODE=ON``,
then the deprecated code will be completely removed resulting in compile
errors for any downstream code still using them.


Hiding entire deprecated header and source files
................................................

In order to hide entire deprecated header and source files when the CMake
variable ``<PackageName>_HIDE_DEPRECATED_CODE=ON`` is set, one needs to move
the headers and sources to another directory and provide for conditional
inclusion in the TriBITS build of the library.  For example, suppose one wants
to remove support for the deprecated files ``SomeOldStuff.hpp`` and
``SomeOldStuff.cpp``.  In this case, one would move the files onto a new
``deprecated/`` sub-directory and then write the ``CMakeLists.txt`` file
like::

  SET(HEADERS "")
  SET(SOURCES "")

  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
  APPEND_GLOB(HEADERS ${DIR}/*.hpp)
  APPEND_GLOB(SOURCES ${DIR}/*.cpp)
  
  IF (NOT ${PACKAGE_NAME}_HIDE_DEPRECATED_CODE)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/deprecated)
    APPEND_SET(HEADERS
      SomeOldStuff.hpp
      )
    APPEND_SET(SOURCES
      SomeOldStuff.cpp
      )
  ENDIF()

  ...

  TRIBITS_ADD_LIBRARY(
    <LIBRARY_NAME>
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )

In this way, when ``${PACKAGE_NAME}_HIDE_DEPRECATED_CODE=TRUE``, then the
directory for the deprecated headers will not be in the include path (so
downstream clients will not even be able to see them) and they will not be
installed so external clients will not be able to see them either.  However,
when ``${PACKAGE_NAME}_HIDE_DEPRECATED_CODE=FALSE``, then these files will be
included in the build and include path and downstream clients can use them.

Once these files need to be permanently removed, one just then needs to remove
them from the version control repository (i.e. ``git rm <files_to_remove>``)
and then remove them from the above ``CMakeLists.txt`` code.


Physically removing deprecated code
+++++++++++++++++++++++++++++++++++

The final step in the code deprecation cycle is to actually remove the
deprecated code.  This is necessary to clean house, remove clutter and finally
get the payoff in the reduction of technical debt that occurs when removing
what is no longer needed or wanted.

It is recommended to remove deprecated files first, then remove deprecated
file fragments from remaining files.  Also, it is recommended to create git
commits after each step.

It is also recommended that some time before deprecated code is actually
removed, that a TriBITS repo change the default of
``<Project>_HIDE_DEPRECATED_CODE`` from ``OFF`` to ``ON`` so that downstream
clients will see the effects of hiding the deprecated code before the code is
actually removed.  In fact, this default should be changed several days to a
week or more before the code is actually removed.  This way, downstream code
developers will get a "shock" about removal of the deprecated code but can
manually configure with ``-D<Project>_HIDE_DEPRECATED_CODE=OFF`` to keep
building in the short-term until they can remove their usage of deprecated
code.


Removing entire deprecated header and source files
..................................................

To remove entire deprecated header and source files one just needs to first
remove them from the version control repository and local directories
(e.g. ``git rm deprecated/*``) and then remove any traces of them from the
``CMakeLists.txt`` file.  For the example in `Hiding entire deprecated header
and source files`_, one would just remove the files ``SomeOldStuff.hpp`` and
``SomeOldStuff.cpp`` from the ``CMakeLists.txt`` file leaving::
  
  IF (NOT ${PACKAGE_NAME}_HIDE_DEPRECATED_CODE)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/deprecated)
    APPEND_SET(HEADERS
      )
    APPEND_SET(SOURCES
      )
  ENDIF()

Since more files may be deprecated later, it may be a good idea to leave the
machinery for conditionally including deprecated files by leaving the above
empty CMake code or just commenting it out.

To find which ``CMakeLists.txt`` files need to be modified, do a search like::

  $ find . -name CMakeLists.txt -exec grep -nH HIDE_DEPRECATED_CODE {} \;

After removing the files, create a local commit of the removed files and the
updated ``CMakeLists.txt`` files before removing deprecated fragments from the
source files.  In other words, do::

  $ emacs -nw CMakeLists.txt  # Remove the references to the deprecated files
  $ git rm SomeOldStuff.hpp SomeOldStuff.cpp
  $ git commit -m "Removing deprecated files"


Removing deprecated code from remaining files
.............................................

The deprecated ifdefed blocks described in `Hiding C/C++ entities`_ can be
removed manually but it is generally preferred to use a tool.  One simple tool
that can do this is called ``unifdef``, that can be downloaded and it is
documented at:

    http://dotat.at/prog/unifdef/

Just download, build, and install the program ``unifdef`` (see
``unifdef/INSTALL`` in untarred source) and then run it as described below.
In the example below, assume the program is installed in the user's home
directory under::

  ~/install/bin/unifdef

For a given TriBITS package, the program is then run as::

  $ find . -name "*pp" -exec ~/install/bin/unifdef \
    -DSomePackage_HIDE_DEPRECATED_CODE {} -o {} \;

After the above command runs, look at the diffs to make sure the ifdef
deprecated code blocks were removed correctly.  For example, run::

  $ git diff -- .

If the diffs look correct, commit the changes::

  $ git commit -m "Removing deprecated code blocks" -- .

Then test everything and push using the `checkin-test.py`_ tool.

After that, all deprecated code is removed and the next period of incremental
change and deprecation begins.


Installation and Backward Compatibility Testing
-----------------------------------------------

TriBITS has some built-in support for installation testing and backward
compatibility testing.  The way it works is that one can install the headers,
libraries, and executables for a TriBITS project and then configure the tests
and examples in the TriBITS project against the installed
headers/libraries/executables.  In this mode, the TriBITS project's libraries
and executables are **not** build and the header file locations to local
source are not included.

When the same version of the project sources are used to build the
tests/examples against the installed headers/libraries/executables, then this
constitutes *installation testing*.  When an older version of the project is
used to build and run tests and examples against the
headers/libraries/executables for a version of the project, then this
constitutes *backward compatibility testing* which also includes installation
testing of course.

ToDo: Describe how the installation testing and backward compatibility testing
process works with some examples.

ToDo: Discuss how this fits into the TriBITS lifecycle model.


Wrapping Externally Configured/Built Software
---------------------------------------------

It is possible to take an external piece of software that uses any arbitrary
build system and wrap it as a TriBITS package and have it integrate in with
the package dependency infrastructure.  The `TribitsExampleProject`_ package
``WrapExternal`` shows how this can be done.  Support for this in TriBITS is
slowly evolving but some key TriBITS features that have been added to support
the arbitrary case include:

* `TRIBITS_DETERMINE_IF_CURRENT_PACKAGE_NEEDS_REBUILT()`_: Uses brute-force
  searches for recently modified files in upstream SE packages to determine if
  the external piece of software needs to be rebuilt.

* `TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES()`_: Write an export
  makefile or a ``XXXConfig.cmake`` file for usage by the wrapped externally
  configured and built software.

.. ToDo: Show how this is done in detail once TriBITS has nice support for
.. external software and the WrapExternal example looks cleaner.

While it is possible to wrap nearly any externally configured and built piece
of software as a TriBITS package, in most cases, it is usually better to just
create a TriBITS build system for the software.  For projects that use a raw
CMake build system, a TriBITS build system can be created side-by-side with
the existing raw CMake build using a number of approaches. The common approach
that is not too invasive is to create a ``CMakeLists.tribits.txt`` file along
side every native ``CMakeLists.txt`` file in the external software project and
have the native ``CMakeLists.txt`` file defined like::

  IF (DOING_A_TRIBITS_BUILD)
    INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.tribits.txt")
    RETURN()
  ENDIF()

  # Rest of native CMakeLists.txt file ...

Experience from the CASL VERA project has shown that, overall, there is less
hassle, less work, and better portability when creating a native TriBITS
build, even if it is a secondary build system for a given piece of software.

The decision whether to just wrap the build system for the existing software
or to create a (secondary) TriBITS build system for it depends on a number of
factors.

Note that while it is generally recommended to create a TriBITS build for an
existing piece of software, it is generally not recommended with switch over
all of the tests to use CTest (unless the existing software is going to ditch
their current test driver and reporting system).  Instead, the external
software's native test system can just be called by the wrapped TriBITS
package in one or more CTest tests.


TriBITS directory snapshotting
------------------------------

Some TriBITS projects choose to snapshot the `TriBITS/tribits/`_ directory
source tree into their project's source tree, typically under
`<projectDir>/cmake/tribits/`_.  The independent ``TriBITS/tribts/`` source
tree contains the tool ``snapshot_tribits.py`` (calls `snapshot-dir.py`_) that
allows one to update the snapshot of the TriBITS source tree as simply as::

  $ cd <projectDir>/cmake/tribits/
  $ <some-base-dir>/TriBITS/tribits/snapshot_tribits.py

This will create a git commit in the local ``<projectDir>/`` git repo that
looks like::

    Automatic snapshot commit from tribits at f8c1682
    
    Origin repo remote tracking branch: 'casl-dev-collab/tribits_reorg_26'
    Origin repo remote repo URL: 'casl-dev-collab = git@casl-dev:collaboration/TriBITS'
    
    At commit:
    
    f8c1682 Assert TriBITS min CMake version in TriBITS itself
    Author: Roscoe A. Bartlett <bartlettra@ornl.gov>
    Date: Fri Dec 5 05:40:49 2014 -0500

This, of course, assumes that ``<projectDir>/`` is a local git repo (or is in
local git repo).  If that is not the case, then one cannot use the script
``snapshot_tribits.py`` or must use it with the ``--skip-commit`` option.

See `snapshot-dir.py --help`_ for more details.  Note the guidance on using a
different branch for the snapshot sync followed by a merge.  This allows for
one to maintain local changes to TriBITS and use git to manage the merges.
However, this will increase the changes of merge conflicts so one should
consider just directly snapshotting into the master branch to avoid merge
conflicts.


TriBITS Development Toolset
---------------------------

Most TriBITS projects need git, a compiler (e.g. GCC), MPI, and a number of
other standard TPLs and other tools in order to develop on and test the
project code.  To this end, TriBITS contains some helper scripts for
downloading, configuring, building, and installing packages like git, cmake,
GCC, MPICH, and others needed to set up a development environment for a
typical computational science software project.  These tools are used to set
up development environments on new machines for projects like Trilinos and
CASL VERA.  Scripts with names like ``install-gcc.py`` are defined which pull
sources from public git repos then configure, build, and install into
specified installation directories.

.. _install_devtools.py:

The script **install_devtools.py** is provided in the directory::

  tribits/devtools_install/

To use this script, one just needs to create some scratch directory like::

  $ mkdir scratch
  $ cd scratch/

then install the tools using, for example::

  $ install_devtools.py --install-dir=~/install/tribits_devtools \
    --parallel=16 --do-all

Then to access installed development environment, one just needs to source the
script::

  ~/install/tribits_devtools/load_dev_env.sh

and then the installed versions of GCC, MPICH, CMake, and `gitdist`_ are
placed in one's path.

See `install_devtools.py --help`_ for more details.


References
==========

.. _Agile Software Development, 2003:

Martin, Robert. *Agile Software Development (Principles, Patterns, and
Practices)*. Prentice Hall. 2003. 

.. _Integration Strategies for CSE, 2009:

Bartlett, Roscoe. *Integration Strategies for Computational Science &
Engineering Software*.  2009-0655, Second International Workshop on Software
Engineering for Computational Science and Engineering, 2009.  http://web.ornl.gov/~8vt/CSE_SoftwareIntegration_Strategies.pdf.

.. _SCALE, 2011:

*SCALE: A Comprehensive Modeling and Simulation Suite for Nuclear Safety Analysis and Design*, ORNL/TM-2005/39, Version 6.1, Oak Ridge National Laboratory, Oak Ridge, Tennessee, June 2011. Available from Radiation Safety Information Computational Center at Oak Ridge National Laboratory as CCC-785.  http://scale.ornl.gov/

.. _LiveV:

*LiveV*  https://github.com/lifev/cmake

TriBITS Detailed Reference Documentation
========================================

The following subsections contain detailed reference documentation for the
various TriBITS variables, functions, and macros that are used by TriBITS
projects that `TriBITS Project Developers`_ need to know about.  Variables,
functions and macros that are used only internally in TriBITS are generally
not documented here (see the TriBITS ``*.cmake`` source files).


TriBITS Global Project Settings
-------------------------------

TriBITS defines a number of global project-level settings that can be set by
the user and can have their default determined by each individual TriBITS
project.  If a given TriBITS project does not define its own default, a
reasonable default is set by the TriBITS system automatically.  These options
are defined and are set, for the most part, in the internal TriBITS function
``TRIBITS_DEFINE_GLOBAL_OPTIONS_AND_DEFINE_EXTRA_REPOS()`` in the TriBITS
CMake code file ``TribitsGlobalMacros.cmake`` which gets called inside of the
`TRIBITS_PROJECT()`_ macro.  That function and that file are the definitive
source the options that a TriBITS project takes and what the default values
are but we strive to document them here as well.  Many of these global options
(i.e. cache variables) such as ``${PROJECT_NAME}_<SOME_OPTION>`` allow the
project to define a default by setting a local variable
``${PROJECT_NAME}_<SOME_OPTION>_DEFAULT`` as::

  SET(${PROJECT_NAME}_<SOME_OPTION>_DEFAULT <someDefault>)

either in its top-level ``CMakeLists.txt`` file or in its
``ProjectName.cmake`` file (depends on what variable it is as to where it
should be set).  If ``${PROJECT_NAME}_<SOME_OPTION>_DEFAULT`` is not set by
the project, then TriBITS provides a reasonable default value.  The TriBITS
code that uses these defaults for this looks like::

  IF ("${${PROJECT_NAME}_<SOME_OPTION>_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_<SOME_OPTION>_DEFAULT <someDefault>)
  ENDIF()

  ADVANCED_SET( ${PROJECT_NAME}_<SOME_OPTION>
    ${PROJECT_NAME}_<SOME_OPTION>_DEFAULT}
    CACHE BOOL "[documentation]."
    )

where ``<SOME_OPTION>`` is an option name, for example like
``TEST_CATEGORIES``, and ``<someDefault>`` is the default set by TriBITS if
the project does not define a default.  In this way, if the project sets the
variable ``${PROJECT_NAME}_<SOME_OPTION>_DEFAULT`` before this code executes,
then ``${${PROJECT_NAME}_<SOME_OPTION>_DEFAULT}`` will be used as the default
for the cache variable ``${PROJECT_NAME}_<SOME_OPTION>`` which, of course, can
be overridden by the user when calling ``cmake`` in a number of ways.

Most of these global options that can be overridden externally by setting the
cache variable ``${PROJECT_NAME}_<SOME_OPTION>`` should be documented in the
`Project-Specific Build Reference`_ document.  A generic version of this
document is found in `TribitsBuildReference`_.  Some of the more unusual
options that might only be of interest to developers mentioned below may not
be documented in `TribitsBuildReference`_.

The global project-level TriBITS options for which defaults can be provided by
a given TriBITS project are:

* `${PROJECT_NAME}_C_Standard`_
* `${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS`_
* `${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_APPEND`_
* `${PROJECT_NAME}_CPACK_SOURCE_GENERATOR`_
* `${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE`_
* `${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES`_
* `${PROJECT_NAME}_ELEVATE_ST_TO_PT`_
* `${PROJECT_NAME}_ENABLE_CPACK_PACKAGING`_
* `${PROJECT_NAME}_ENABLE_CXX11`_
* `${PROJECT_NAME}_ENABLE_CXX`_
* `${PROJECT_NAME}_ENABLE_C`_
* `${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE`_
* `${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES`_
* `${PROJECT_NAME}_ENABLE_Fortran`_
* `${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES`_
* `${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE`_
* `${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION`_
* `${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES`_
* `${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE`_
* `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_
* `${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS`_
* `${PROJECT_NAME}_REQUIRES_PYTHON`_
* `${PROJECT_NAME}_SET_INSTALL_RPATH`_
* `${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME`_
* `${PROJECT_NAME}_TEST_CATEGORIES`_
* `${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS`_
* `${PROJECT_NAME}_TRACE_ADD_TEST`_
* `${PROJECT_NAME}_USE_GNUINSTALLDIRS`_
* `${PROJECT_NAME}_USES_PYTHON`_
* `DART_TESTING_TIMEOUT`_
* `CMAKE_INSTALL_RPATH_USE_LINK_PATH`_
* `MPI_EXEC_MAX_NUMPROCS`_
* `PythonInterp_FIND_VERSION`_

These options are described below.

.. _${PROJECT_NAME}_C_Standard:

**${PROJECT_NAME}_C_Standard**

  The variable ``${PROJECT_NAME}_C_Standard`` is used define the C standard
  pass to the compiler in ``--std=<cstd>`` for GCC builds of the project.
  TriBITS sets the default as ``c99`` but the project can set a new default in
  the project's base `<projectDir>/CMakeLists.txt`_ file with, for example::

    SET(${PROJECT_NAME}_C_Standard_DEFAULT c11)

.. _${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS:

**${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS**

  The variable ``${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS`` determines how
  unparsed and otherwise ignored arguments are handled in TriBITS functions
  that are called by the client TriBITS projects.  These are arguments that
  are left over from parsing input options to functions and macros that take
  both positional arguments and keyword arguments/options handled with the
  ``CMAKE_PARSE_ARGUMENTS()`` function.  For example, for the a TriBITS
  function declared like::

    TRIBITS_COPY_FILES_TO_BINARY_DIR(
      <targetName>
      [SOURCE_FILES <file1> <file2> ...]
      [SOURCE_DIR <sourceDir>]
      ...
      )

  the arguments ``SOURCE_FILES <file1> <file2> ...`` and those that follow are
  parsed by the ``CMAKE_PARSE_ARGUMENTS()`` function while the argument
  ``<targetName>`` is a positional argument.  The problem is that any
  arguments passed between the first ``<targetName>`` argument and the
  specified keyword arguments like ``SOURCE_FILES`` and ``SOURCE_DIR`` are
  returned as unparsed arguments and are basically ignored (which is what
  happened in earlier versions of TriBITS).  For example, calling the function
  as::

    TRIBITS_COPY_FILES_TO_BINARY_DIR( FooTestCopyFiles
      ThisArgumentIsNotParsedAndIsIgnored
      SOURCE_FILES file1.cpp file2.cpp ...
      ...
      )

  would result in the unparsed argument
  ``ThisArgumentIsNotParsedAndIsIgnored``.

  The value of ``${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS`` determines how
  that ignored argument is handled.  If the value is ``WARNING``, then it will
  just result in a ``MESSAGE(WARNING ...)`` command that states the warning
  but configure is allowed to be completed. This would be the right value to
  allow an old TriBITS project to keep configuring until the warnings can be
  cleaned up.  If the value is ``SEND_ERROR``, then ``MESSAGE(SEND_ERROR
  ...)`` is called.  This will result in the configure failing but will allow
  configure to continue until the end (or a ``FATAL_ERROR`` is raised).  This
  would be the right value when trying to upgrade a TriBITS project where you
  wanted to see all of the warnings when upgrading TriBITS (so you could fix
  them all in one shot).  Finally, the value of ``FATAL_ERROR`` will result in
  ``MESSAGE(FATAL_ERROR ...)`` being called which will halt configure right
  away.  This is the best value when developing on a TriBITS project that is
  already clean but you want to catch new developer-inserted errors right
  away.

  The default value for ``${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS`` is
  ``WARNING``, so that it will be backward compatible for TriBITS projects
  that might have previously undetected unparased and therefore ignored
  argument .  However, a project can change the default by setting, for
  example::

    SET(${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS_DEFAULT FATAL_ERROR)

  in the `<projectDir>/ProjectName.cmake`_ file.

  The user of a TriBITS project should not be able to trigger this unparsed
  arguments condition so this variable is not documented in the `TriBITS Build
  Reference`_.  But it is still a CMake cache var that is documented in the
  CMakeCache.txt file and can be set by the user or developer if desired.

.. _${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_APPEND:

**${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_APPEND**

  The variable ``${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_APPEND`` is used to
  define the absolute path to a file (or a list of files) that should be
  included after the files listed in
  ``${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE``.  This variable can be used by
  the TriBITS project to define, for example, a standard set of development
  environments in the base `<projectDir>/CMakeLists.txt`_ file with::

    SET(${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_APPEND_DEFAULT
      "${CMAKE_CURRENT_LIST_DIR}/cmake/StdDevEnvs.cmake")

  **before** the `TRIBITS_PROJECT()`_ command.  By including this file(s)
  after the file(s) listed in ``${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE``, the
  user can override the variables set in this appended file(s).  But it is
  important that these variables best set after the users options have been
  set but before the Package and TPL dependency analysis is done (because this
  might enable or disable some TPLs).

.. _${PROJECT_NAME}_CPACK_SOURCE_GENERATOR:

**${PROJECT_NAME}_CPACK_SOURCE_GENERATOR**

  The variable ``${PROJECT_NAME}_CPACK_SOURCE_GENERATOR`` determines the CPack
  source generation types that are created when the ``package_source`` target
  is run.  The TriBITS default is set to ``TGZ``.  However, this default can
  be overridden by setting, for example::

    SET(${PROJECT_NAME}_CPACK_SOURCE_GENERATOR_DEFAULT "TGZ;TBZ2")

  This variable should generally be set in the file::

     <projectDir>/cmake/CallbackDefineProjectPackaging.cmake

  instead of in the base-level ``CMakeLists.txt`` file so that it goes along
  with rest of the project-specific CPack packaging options.

.. _${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE:

**${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE**

  The variable ``${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE`` determines if the
  CTest driver scripts using `TRIBITS_CTEST_DRIVER()`_ configure, build, test
  and submit results to CDash all-at-once for all of the packages being tested
  or if instead is done package-by-package.  Currently, the default is set to
  ``FALSE`` for the package-by-package mode (for historical reasons) but the
  default can be set to ``TRUE`` by setting:

    SET(${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE_DEFAULT "TRUE")

  in the project's `<projectDir>/ProjectName.cmake`_ file.  (This default must
  be changed in the ``<projectDir>/ProjectName.cmake`` file and **NOT** the
  `<projectDir>/CMakeLists.txt`_ file because the latter is not directly
  processed in CTest -S driver scripts using ``TRIBITS_CTEST_DRIVER()``.)

  In general, a project should change the default to ``TRUE`` when the minimum
  CMake version being used with the project is CMake 3.10+ and when using a
  newer CDash installation that can accomidate the results coming from ctest
  -S and display them package-by-package very nicely.  Otherwise, most
  projects are better off with package-by-package mode since it results in
  nicer display on CDash.

.. _${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES:

**${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES**

  If `${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON`_ (the TriBITS
  default value), then any explicitly enabled packages that have disabled
  `upstream`_ required packages or TPLs will be disabled.  If
  `${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=OFF`_, then an
  configure error will occur.  For more details also see
  `TribitsBuildReference`_ and `Disables trump enables where there is a
  conflict`_.  A project can define a different default value by setting::
  
    SET(${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT FALSE)

.. _${PROJECT_NAME}_ELEVATE_ST_TO_PT:

**${PROJECT_NAME}_ELEVATE_ST_TO_PT**

  If ``${PROJECT_NAME}_ELEVATE_ST_TO_PT`` is set to ``ON``, then all ``ST`` SE
  packages will be elevated to ``PT`` packages.  The TriBITS default is
  obviously ``OFF``.  The default can be changed by setting::

    SET(${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT ON)

  There are projects, especially meta-projects, where the distinction between
  ``PT`` and ``ST`` code is not helpful or the assignment of ``PT`` and ``ST``
  packages in a repository is not appropriate with respect to the outer
  meta-project.  An example project like this CASL VERA.  Changing the default
  to ``ON`` allows any and packages to be considered in pre-push testing.

.. _${PROJECT_NAME}_ENABLE_CPACK_PACKAGING:

**${PROJECT_NAME}_ENABLE_CPACK_PACKAGING**

  If ``${PROJECT_NAME}_ENABLE_CPACK_PACKAGING`` is ``ON``, then CPack support
  is enabled and some TriBITS code is run that is needed to set up
  data-structures that are used by the built-in CMake target
  ``package_source``.  The TriBITS default is ``OFF`` with the idea that the
  average developer or user will not be wanting to create source distributions
  with CPack.  However, this default can be changed by setting::

    SET(${PROJECT_NAME}_ENABLE_CPACK_PACKAGING_DEFAULT ON)

.. _${PROJECT_NAME}_ENABLE_CXX11:
  
**${PROJECT_NAME}_ENABLE_CXX11**
  
  If ``${PROJECT_NAME}_ENABLE_CXX11`` is ``ON``, then C++ compiler options
  that turn on C++11 support will be searched for.  By default, TriBITS sets
  this to ``OFF`` for all systems.  However, if project requires C++11 support
  by default, then the project should set the default:
  
    SET(${PROJECT_NAME}_ENABLE_CXX11_DEFAULT TRUE)

.. _${PROJECT_NAME}_ENABLE_CXX:
  
**${PROJECT_NAME}_ENABLE_CXX**
  
  If ``${PROJECT_NAME}_ENABLE_CXX`` is ``ON``, then C++ language support for
  the project will be enabled and the C++ compiler must be found.  By default,
  TriBITS sets this to ``ON`` for all systems.  A project never requires C++
  can set this to off by default by setting:
  
    SET(${PROJECT_NAME}_ENABLE_CXX_DEFAULT FALSE)

.. _${PROJECT_NAME}_ENABLE_C:
  
**${PROJECT_NAME}_ENABLE_C**
  
  If ``${PROJECT_NAME}_ENABLE_C`` is ``ON``, then C language support for the
  project will be enabled and the C compiler must be found.  By default,
  TriBITS sets this to ``ON`` for all systems.  A project never requires C can
  set this to off by default by setting:
  
    SET(${PROJECT_NAME}_ENABLE_C_DEFAULT FALSE)
  
  If a project does not have any native C code a good default would be::
  
    SET(${PROJECT_NAME}_ENABLE_C_DEFAULT FALSE)
  
  NOTE: It is usually not a good idea to always force off C, or any compiler,
  because extra repositories and packages might be added by someone that might
  require the compiler and we don't want to unnecessarily limit the generality
  of a given TriBITS build.  Setting the default for all platforms should be
  sufficient.

.. _${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE:
.. _${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT:

**${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE**

  The variable ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE`` switches the
  TriBITS project from development mode to release mode.  The default for this
  variable ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT`` should be set
  in the project's `<projectDir>/Version.cmake`_ file and switched from ``ON``
  to ``OFF`` when creating a release (see `Project and Repository Versioning
  and Release Mode`_).  When ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE`` is
  ``ON``, several other variables are given defaults appropriate for
  development mode.  For example, ``${PROJECT_NAME}_ASSERT_MISSING_PACKAGES``
  is set to ``ON`` by default in development mode but is set to ``OFF`` by
  default in release mode.  In addition, strong compiler warnings are enabled
  by default in development mode but are disabled by default in release mode.
  This variable also affects the behavior of `TRIBITS_SET_ST_FOR_DEV_MODE()`_.
  
.. _${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES:

**${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES**
  
  If ``${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES`` is ``ON``, then
  ``Makefile.export.<PackageName>`` will get created at configure time in the
  build tree and installed into the install tree.  See `TribitsBuildReference`_
  for details.  The TriBITS default is ``ON`` but a project can decide to turn
  this off by default by setting::
  
    SET(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES_DEFAULT OFF)
  
  A project might want to disable the generation of export makefiles by default
  if its main purpose is to provide executables.  There is no reason to provide
  an export makefile if libraries and headers are not actually installed (see
  `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_)
 
.. _${PROJECT_NAME}_ENABLE_Fortran:
  
**${PROJECT_NAME}_ENABLE_Fortran**
  
  If ``${PROJECT_NAME}_ENABLE_Fortran`` is ``ON``, then Fortran support for
  the project will be enabled and the Fortran compiler(s) must be found.  By
  default, TriBITS sets this to ``ON`` .
  
  If a project does not have any native Fortran code a good default would be::
  
    SET(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT OFF)

  This default can be set in `<projectDir>/ProjectName.cmake`_ or `<projectDir>/CMakeLists.txt`_.

  Given that a native Fortran compiler is not supported by default on Windows
  and on most Mac OSX systems, projects that have optional Fortran code may
  decide to set the default depending on the platform by setting, for example::

    IF ( (WIN32 AND NOT CYGWIN) OR (CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin") )
      MESSAGE(STATUS "Warning: Setting ${PROJECT_NAME}_ENABLE_Fortran=OFF by default"
       " because this is Windows (not cygwin) and we assume to not have Fortran!")
      SET(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT OFF)
    ELSE()
      SET(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT ON)
    ENDIF()
  
  NOTE: It is usually not a good idea to always force off Fortran, or any
  compiler, because extra repositories and packages might be added by someone
  that might require the compiler and we don't want to unnecessarily limit the
  generality of a given TriBITS build.  Setting the default for all platforms
  should be sufficient.

.. _${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES:

**${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES**

  If ``${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES`` is set to ``ON``,
  then ``<PackageName>Config.cmake`` files are created at configure time in
  the build tree and installed into the install tree.  These files are used by
  external CMake projects to pull in the list of compilers, compiler options,
  include directories and libraries.  The TriBITS default is ``ON`` but a
  project can change the default by setting, for example::

    SET(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES_DEFAULT OFF)

  A project would want to turn off the creation and installation of
  ``<PackageName>Config.cmake`` files if it was only installing and providing
  executables (see `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_).

.. _${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE:

**${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE**

  If ``${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE`` is ``ON``, then packages
  and subpackages marked as ``ST`` in the `<repoDir>/PackagesList.cmake`_ file
  will be implicitly enabled along with the ``PT`` packages.  Additional code
  and tests may also be enabled using this option.  The TriBITS default is
  ``OFF`` but this can be changed by setting::

    SET(${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE_DEFAULT ON)

  in the `<projectDir>/ProjectName.cmake`_ file.

.. _${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION:

**${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION**

  If ``${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION`` is
  ``TRUE``, then the directories for subpackages that are not enabled are left
  out of the source tarball.  This reduces the size of the tarball as much as
  possible but does require that the TriBITS packages and subpackages be
  properly set up to allow disabled subpackages from being excluded.  The
  TriBITS default is ``TRUE`` but this can be changed by setting::

    SET(${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION_DEFAULT FALSE)

.. _${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE:

**${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE**

  If ``${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE`` is ``ON``, then the file
  ``<Project>RepoVersion.txt`` will get generated as a byproduct of
  configuring with CMake.  See `Multi-Repository Support`_ and
  `<Project>_GENERATE_REPO_VERSION_FILE`_.  The default is ``OFF`` but the
  project can change that by setting::

    SET(${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT ON)

  in the `<projectDir>/ProjectName.cmake`_ file.

.. _${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES:

**${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES**

  If ``${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES`` is ``ON``, then the
  data-structures needed to generate ``Makefile.export.<PackageName>`` and
  ``<PackageName>Config.cmake`` are created.  These data structures are also
  needed in order to generate export makefiles on demand using the function
  `TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES()`_.  The default in
  TriBITS is to turn this ``ON`` automatically by default if
  ``${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES`` or
  ``${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES`` are ``ON``.  Else, by
  default, TriBITS sets this to ``OFF``.  The only reason for the project to
  override the default is to set it to ``ON`` as with::

    SET(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT ON)

  is so that the necessary data-structures are generated in order to use the
  function `TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES()`_.
  
.. _${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS:

**${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS**

  If ``${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`` is set to ``ON``, then
  any defined libraries or header files that are listed in calls to
  `TRIBITS_ADD_LIBRARY()`_ or `TRIBITS_INSTALL_HEADERS()`_ will be installed
  (unless options are passed into `TRIBITS_ADD_LIBRARY()`_ that disable
  installs).  If set to ``OFF``, then headers and libraries will *not* be
  installed by default and only ``INSTALLABLE`` executables added with
  `TRIBITS_ADD_EXECUTABLE()`_ will be installed.  However, as described in
  `TribitsBuildReference`_, shared libraries will always be installed if
  enabled since they are needed by the installed executables.
  
  For a TriBITS project that is primarily delivering libraries
  (e.g. Trilinos), then it makes sense to leave the TriBITS default which is
  ``ON`` or explicitly set::
  
    SET(${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT  ON)
  
  For a TriBITS project that is primarily delivering executables (e.g. VERA),
  then it makes sense to set the default to::
  
    SET(${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT  OFF)


.. _${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS:

**${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS**

  Determines if all of the libraries listed in ``<tplName>_LIBRARY_NAMES`` for
  a given TPL must be found for each enabled TPL.  By default, this is
  ``FALSE`` which means that the determination if all of the listed libs for a
  TPL should be found is determined by the ``MUST_FIND_ALL_LIBS`` option to
  the `TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()`_ function in the TPL
  find module.  To change the default for this, set::

    SET(${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS_DEFAULT  TRUE)

  in the `<projectDir>/ProjectName.cmake`_ file.


.. _${PROJECT_NAME}_REQUIRES_PYTHON:

**${PROJECT_NAME}_REQUIRES_PYTHON**

  If the TriBITS project requires Python, set::

    SET(${PROJECT_NAME}_REQUIRES_PYTHON  TRUE)

  in the `<projectDir>/ProjectName.cmake`_ file (See `Python Support`_).  The
  default is implicitly ``FALSE``.


.. _${PROJECT_NAME}_SET_INSTALL_RPATH:

**${PROJECT_NAME}_SET_INSTALL_RPATH**

  The cache variable ``${PROJECT_NAME}_SET_INSTALL_RPATH`` is used to define
  the default RPATH mode for the TriBITS project (see `Setting install RPATH`_
  for details).  The TriBITS default is to set this to ``TRUE`` but the
  TriBITS project can be set the default to ``FALSE`` by setting::

    SET(${PROJECT_NAME}_SET_INSTALL_RPATH_DEFAULT FALSE)

  in the project's `<projectDir>/ProjectName.cmake`_ file (see `RPATH
  Handling`_).


.. _${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME:

**${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME**

  The cache variable ``${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME``
  determines if the start and end date/time for each advanced test (i.e. added
  with `TRIBITS_ADD_ADVANCED_TEST()`_) is printed or not with each test.  The
  TriBITS default is ``OFF`` but a TriBITS project can change this default by
  setting::

    SET(${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME_DEFAULT ON)

  The implementation of this feature currently uses ``EXECUTE_PROCESS(date)``
  and therefore will only work on Linux/Unix/Mac systems and not Windows
  systems.

  NOTE: In a future version of CTest, this option may turn on start and end
  date/time for regular tests added with `TRIBITS_ADD_TEST()`_ (which uses a
  raw command with ``ADD_TEST()``).

.. _${PROJECT_NAME}_SKIP_EXTRAREPOS_FILE:

**${PROJECT_NAME}_SKIP_EXTRAREPOS_FILE**

  The cache variable ``${PROJECT_NAME}_SKIP_EXTRAREPOS_FILE`` is set in the
  `<projectDir>/ProjectName.cmake`_ file as::

    SET(${PROJECT_NAME}_SKIP_EXTRAREPOS_FILE TRUE)

  for projects that don't have a
  `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file.  This variable needs
  to be set when using the CTest driver script and does not need to be set for
  the basic configure and build process.

.. _${PROJECT_NAME}_TEST_CATEGORIES:
.. _${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT:

**${PROJECT_NAME}_TEST_CATEGORIES**

  The cache variable ``${PROJECT_NAME}_TEST_CATEGORIES`` determines what tests
  defined using `TRIBITS_ADD_TEST()`_ and `TRIBITS_ADD_ADVANCED_TEST()`_ will
  be added for ``ctest`` to run (see `Test Test Category`_).  The TriBITS
  default is ``NIGHTLY`` for a standard local build.  The `checkin-test.py`_
  tool sets this to ``BASIC`` by default.  A TriBITS project can override the
  default for a basic configure using, for example::

    SET(${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT BASIC)

  The justification for having the default `Test Test Category`_ be
  ``NIGHTLY`` instead of ``BASIC`` is that when someone is enabling a package
  to develop on it or install it, we want them by default to be seeing the
  full version of the test suite (shy of the `Test Test Category HEAVY`_
  tests which can be very expensive) for the packages they are explicitly
  enabling.  Typically they will not be enabling forward/`downstream`_
  dependent packages so the cost of running the test suite should not be too
  prohibitive.  This all depends on how good of a job the development teams do
  in making their test suites run fast and keeping the cost of running the
  tests down.  See the section `TriBITS Automated Testing`_ for a more
  detailed discussion.

.. _${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS:
.. _${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS_DEFAULT:

**${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS**

  If ``${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS`` is set to ``TRUE``, then the
  ``SYSTEM`` flag will be passed into the ``INCLUDE_DIRECTORIES()`` command
  for TPL include directories for every TPL for every package, by default.  On
  some systems this will result in include directories being passed to the
  compiler with ``-isystem`` instead of ``-I``.  This helps to avoid compiler
  warning coming from TPL header files for C and C++.  However, with CMake
  version 3.2 and less, this also results in ``-isystem`` being passed to the
  Fortran compiler (e.g. gfortran) as well.  This breaks the reading of
  Fortran module files (perhaps a bug in gfortran).  Because of this issue
  with Fortran, the TriBITS default for this option is set to ``FALSE`` but a
  project can override the default using::

    SET(${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS_DEFAULT  TRUE)

  (This would be a good default if the project has not Fortran files or has
  not Fortran files that use modules provided by TPLs).

  However, if a package or subpackage sets::

    SET(${PACKAGE_NAME}_SKIP_TPL_SYSTEM_INCLUDE_DIRS  TRUE)

  in its ``CMakeLists.txt`` files before the ``TRIBITS_ADD_LIBRARY()`` or
  ``TRIBITS_ADD_EXECUTABLE()`` commands are called in that package, then
  ``SYSTEM`` will **not** be passed into ``INCLUDE_DIRECTORIES()`` for TPL
  include dirs.  This is how some TriBITS packages with Fortran files that use
  Fortran modules avoid passing in ``-isystem`` to the Fortran compiles and
  thereby avoid the defect with gfortran described above.  If CMake version
  3.3 or greater is used, this variable is not required.

  NOTE: Currently, a TriBITS SE package must have a direct dependency on a TPL
  to have ``-isystem`` added to a TPL's include directories on the compile
  lines for that package.  That is, the TPL must be listed in the
  ``LIB_REQUIRED_TPLS`` or ``LIB_OPTIONAL_TPLS`` arguments passed into the
  `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_ function in the SE package's
  `<packageDir>/cmake/Dependencies.cmake`_ file.  In addition, to have
  ``-isystem`` added to the include directories for a TPL when compiling the
  tests for an SE package, it must be listed in the ``TEST_REQUIRED_TPLS`` or
  ``TEST_OPTIONAL_TPLS`` arguments.  This is a limitation of the TriBITS
  implementation that will be removed in a future version of TriBITS.

.. _${PROJECT_NAME}_TRACE_ADD_TEST:
.. _${PROJECT_NAME}_TRACE_ADD_TEST_DEFAULT:

**${PROJECT_NAME}_TRACE_ADD_TEST**

  If ``${PROJECT_NAME}_TRACE_ADD_TEST`` is set to ``TRUE``, then a single line
  will be printed for each call to `TRIBITS_ADD_TEST()`_ and
  `TRIBITS_ADD_ADVANCED_TEST()`_ for if the test is added or not and if not
  then why.  The default is set based on the value of
  ``${PROJECT_NAME}_VERBOSE_CONFIGURE`` but a project can override the default
  by setting::

    SET(${PROJECT_NAME}_TRACE_ADD_TEST_DEFAULT  TRUE)

.. _${PROJECT_NAME}_USE_GNUINSTALLDIRS:

**${PROJECT_NAME}_USE_GNUINSTALLDIRS**

  If ``${PROJECT_NAME}_USE_GNUINSTALLDIRS`` is set to ``TRUE``, then the
  default install paths will be determined by the standard CMake module
  ``GNUInstallDirs``.  Otherwise, platform independent install paths are used
  by default.

  A project can use the paths given the cmake module ``GNUInstallDirs`` by
  default by setting::
  
    SET(${PROJECT_NAME}_USE_GNUINSTALLDIRS_DEFAULT FALSE)

  in the project's top-level `<projectDir>/CMakeLists.txt`_ file or its
  `<projectDir>/ProjectName.cmake`_ file.  The default is ``FALSE``.

.. _${PROJECT_NAME}_USES_PYTHON:

**${PROJECT_NAME}_USES_PYTHON**

  If the TriBITS project can use Python, but does not require it, set::

    SET(${PROJECT_NAME}_USES_PYTHON  TRUE)

  in the `<projectDir>/ProjectName.cmake`_ file (see `Python Support`_).  The
  default for a TriBITS project is implicitly ``TRUE``.  To explicitly state
  that Python is never needed, set::

    SET(${PROJECT_NAME}_USES_PYTHON  FALSE)

.. _DART_TESTING_TIMEOUT:

**DART_TESTING_TIMEOUT**

  The cache variable ``DART_TESTING_TIMEOUT`` is a built-in CMake variable
  that provides a default timeout for all tests (see `Setting test timeouts at
  configure time`_).  By default, TriBITS defines this to be 1500 seconds
  (which is also the raw CMake default) but the project can change this
  default, from 1500 to 300 for example, by setting the following in the
  project's `<projectDir>/ProjectName.cmake`_ or
  `<projectDir>/CMakeLists.txt`_ file::

    SET(DART_TESTING_TIMEOUT_DEFAULT 300)

.. _CMAKE_INSTALL_RPATH_USE_LINK_PATH:

**CMAKE_INSTALL_RPATH_USE_LINK_PATH**

  The cache variable ``CMAKE_INSTALL_RPATH_USE_LINK_PATH`` is a built-in CMake
  variable that determines if the paths for external libraries (i.e. from
  TPLs) is put into the installed library RPATHS (see `RPATH Handling`_).
  TriBITS sets the default for this to ``TRUE`` but a project can change the
  default back to ``FALSE`` by setting the following in the project's
  `<projectDir>/ProjectName.cmake`_ file::

    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH_DEFAULT FALSE)

.. _MPI_EXEC_MAX_NUMPROCS:

**MPI_EXEC_MAX_NUMPROCS**

  The variable ``MPI_EXEC_MAX_NUMPROCS`` gives the maximum number of processes
  for an MPI test that will be allowed as defined by `TRIBITS_ADD_TEST()`_ and
  `TRIBITS_ADD_ADVANCED_TEST()`_.  The TriBITS default is set to be ``4`` (for
  no good reason really but it needs to stay that way for backward
  compatibility).  This default can be changed by setting::

    SET(MPI_EXEC_MAX_NUMPROCS_DEFAULT <newDefaultMax>)

  While this default can be changed for the project as a whole on all
  platforms, it is likely better to change this default on a
  machine-by-machine basis to correspond to the load that can be accommodated
  by a given machine (or class of machines).  For example if a given machine
  has 64 cores, a reasonable number for ``MPI_EXEC_MAX_NUMPROCS_DEFAULT`` is
  64.

.. _PythonInterp_FIND_VERSION:

**PythonInterp_FIND_VERSION**

  Determines the version of Python that is looked for.  TriBITS requires at
  least version "2.4".  A particular TriBITS project can require a higher
  version of TriBITS and this is set using, for example:

    SET(PythonInterp_FIND_VERSION_DEFAULT "2.6.6")

  in the `<projectDir>/ProjectName.cmake`_ file (See `Python Support`_).  The
  default is version "2.4".  The user can force a more recent version of
  Python by configuring with, for example::

    -D PythonInterp_FIND_VERSION="2.7.3"


TriBITS Macros and Functions
----------------------------

The following subsections give detailed documentation for the CMake macros and
functions that make up the core TriBITS system.  These are what are used by
TriBITS project developers in their ``CMakeLists.txt`` and other files.  All
of these functions and macros should be automatically available when
processing the project's and package's variables files if used properly.
Therefore, no explicit ``INCLUDE()`` statements should be needed other than
the initial include of the ``TriBITS.cmake`` file in the top-level
`<projectDir>/CMakeLists.txt`_ file so the command `TRIBITS_PROJECT()`_ can be
executed.

.. include:: TribitsMacroFunctionDoc.rst


General Utility Macros and Functions
------------------------------------

The following subsections give detailed documentation for some CMake macros
and functions which are *not* a core part of the TriBITS system but are
included in the TriBITS source tree, are used inside of the TriBITS system,
and are provided as a convenience to TriBITS project developers.  One will see
many of these functions and macros used throughout the implementation of
TriBITS and even in the ``CMakeLists.txt`` files for different projects that
use TriBITS.

These macros and functions are *not* prefixed with ``TRIBITS_``.  However,
there is really not a large risk to defining and using these non-namespaces
utility functions and macros.  It turns out that CMake allows one to redefine
any macro or function, even built-in ones, inside of one's project.
Therefore, even if CMake did add new commands that clashed with these names,
there would be no conflict.  When overriding a built-in command,
e.g. ``some_builtin_command()``, one can always access the original built-in
command as ``_some_builtin_command()``.


.. include:: UtilsMacroFunctionDoc.rst

.. ToDo: Edited and spell-checked through here on 5/15/2014.


FAQ
===

**Q:** Why does not TriBITS just use the standard CMake ``Find<PACKAGE_NAME>.cmake``
modules and the standard ``FIND_PACKAGE()`` function to find TPLs?

**A:** The different "standard" CMake ``Find<PACKAGE_NAME>.cmake`` modules do not
have a standard set of outputs and therefore, can't be handled in a uniform
way.  For example, 


Appendix
========


History of TriBITS
------------------

TriBITS started development in November 2007 as a set of helper macros to
provide a CMake build system for a small subset of packages in Trilinos.  The
initial goal was to support a native Windows build (using Visual C++) to
compile and install these few Trilinos packages on Windows for usage by
another project (the Sandia Titan project which included VTK).  At that time,
Trilinos was using a highly customized and augmented autotools build system.
Initially, this CMake system was just a set of macros to streamline creating
executables and tests.  Some of the conventions started in that early effort
(e.g. naming conventions of variables and macros where functions use upper
case like old FORTRAN and variables are mixed case) were continued in later
efforts and are reflected in the current implementation.  Then, stating in
early 2008, a more detailed evaluation was performed to see if Trilinos should
switch over to CMake as the default (and soon only) supported build and test
system (see "Why CMake?" in `TriBITS Overview`_).  This lead to the initial
implementation of a scalable package-based architecture (PackageArch) for the
Trilinos CMake project in late 2008.  This Trilinos CMake PackageArch system
evolved over the next few years with development in the system slowing
into 2010.  This Trilinos CMake build system was then adopted as the build
infrastructure for the CASL VERA effort in 2011 where CASL VERA packages were
treated as add-on Trilinos packages (see Section `Multi-Repository Support`_).
Over the next year, there was significant development of the system to support
larger multi-repo projects in support of CASL VERA.  That lead to the decision
to formally generalize the Trilinos CMake PackageArch build system outside of
Trilinos and the name TriBITS was formally adopted in November 2011.  Work to
refactor the Trilinos CMake system into a general reusable stand-alone
CMake-based build system started in October 2011 and an initial implementation
was complete in December 2011 when it was used for the CASL VERA build system.
In early 2012, the ORNL CASL-related projects Denovo and SCALE (see [`SCALE,
2011`_]) adopted TriBITS as their native development build systems.  Shortly
after, TriBITS was adopted as the native build system for the CASL-related
University of Michigan code MPACT.  In addition to being used in CASL, all of
these codes also had a significant life outside of CASL.  Because they used
the same TriBITS build system, it proved relatively easy to keep these various
codes integrated together in the CASL VERA code meta-build.  At the same time,
TriBITS well served the independent development teams and non-CASL projects
independent from CASL VERA.  Since the initial extraction of TriBITS from
Trilinos, the TriBITS system was further extended and refined, driven by CASL
VERA development and expansion.  Independently, an early version of TriBITS
from 2012 was adopted by the LiveV
project (see [`LiveV`_]) which was forked and extended
independently.

Why a TriBITS Package is not a CMake Package
--------------------------------------------

Note that a `TriBITS Package`_ is not the same thing as a "Package" in raw
CMake terminology.  In raw CMake, a "Package" is some externally provided bit
of software or other utility for which the current CMake project has an
optional or required dependency (see `CMake: How to Find Libraries
<http://www.cmake.org/Wiki/CMake:How_To_Find_Libraries>`_).  Therefore, a raw
CMake "Package" actually maps to a `TriBITS TPL`_.  A raw CMake "Package"
(e.g. Boost, CUDA, etc.)  can be found using a standard CMake find module
``Find<rawPackageName>.cmake`` using the built-in CMake command
``FIND_PACKAGE(<rawPackageName>)``.  It is unfortunate that the TriBITS and
the raw CMake definitions of the term "Package" are not exactly the same.
However, the term "Package" was coined by the Trilinos project long ago before
CMake was adopted as the Trilinos build system and Trilinos' definition of
"Package" (going back to 1998) pre-dates the development of CMake (see
`History of CMake <http://en.wikipedia.org/wiki/CMake#History>`_) and
therefore Trilinos dictated the terminology of TriBITS and the definition of
the term "Package" in the TriBITS system.  However, note that both meanings of
the term "Package" are consistent with the more general software engineering
definition of a "Package" according to `Software Engineering Packaging
Principles`_.


Design Considerations for TriBITS
---------------------------------

Some of the basic requirements and design goals for TriBITS are outlined in the `TriBITS Overview`_ document.

.. ToDo: Discuss design requirements for TriBITS in more detail.

As stated in `TriBITS Dependency Handling Behaviors`_, `No circular
dependencies of any kind are allowed`_.  That is, no TriBITS SE package (or
its tests) can declare a dependency on a `downstream`_ SE package, period!  To
some, this might seem over constraining but adding support for circular
dependencies to the TriBITS system would add significant complexity and
space/time overhead and is a bad idea from a basic software engineering
perspective (see the *ADP (Acyclic Dependencies Principle)* in `Software
Engineering Packaging Principles`_).  From a versioning, building, and
change-prorogation perspective, any packages involved in a circular dependency
would need to be treated as a single software engineering package anyway so
TriBITS forces development teams to glob all of this stuff together into a
single `TriBITS SE Package`_ when cycles in software exist.  There are
numerous wonderful ways to break circular dependencies between packages that
are proven and well established in the SE community (for example, see [`Agile
Software Development, 2003`_]).

.. ToDo: Discuss why it is a good idea to explicitly list packages instead of
.. just searching for them.  Hint: Think performance and circular
.. dependencies!


clone_extra_repos.py --help
---------------------------

Below is a snapshot of the output from ``clone_extra_repos.py --help``.  For
more details on the usage of ``clone_extra_repos.py``, see `Multi-Repository
Support`_ and `Multi-Repository Development Workflow`_.

.. include:: clone_extra_repos-help.txt
   :literal:

gitdist documentation
---------------------

The sections below show snapshots of the output from the `gitdist`_ tool from
`gitdist --help`_ and ``gitdist --dist-help=<topic>``:

* `gitdist --help`_
* `gitdist --dist-help=overview`_
* `gitdist --dist-help=repo-selection-and-setup`_
* `gitdist --dist-help=dist-repo-status`_
* `gitdist --dist-help=repo-versions`_
* `gitdist --dist-help=aliases`_
* `gitdist --dist-help=default-branch`_
* `gitdist --dist-help=move-to-base-dir`_
* `gitdist --dist-help=usage-tips`_
* `gitdist --dist-help=script-dependencies`_
* `gitdist --dist-help=all`_

For more details on the usage of ``gitdist``, see `Multi-Repository Support`_
and `Multi-Repository Development Workflow`_.


gitdist --help
++++++++++++++

.. include:: gitdist-help.txt
   :literal:

gitdist --dist-help=overview
++++++++++++++++++++++++++++

.. include:: gitdist-dist-help-overview.txt
   :literal:

gitdist --dist-help=repo-selection-and-setup
++++++++++++++++++++++++++++++++++++++++++++

.. include:: gitdist-dist-help-repo-selection-and-setup.txt
   :literal:

gitdist --dist-help=dist-repo-status
++++++++++++++++++++++++++++++++++++

.. include:: gitdist-dist-help-dist-repo-status.txt
   :literal:

gitdist --dist-help=repo-versions
+++++++++++++++++++++++++++++++++

.. include:: gitdist-dist-help-repo-versions.txt
   :literal:

gitdist --dist-help=aliases
+++++++++++++++++++++++++++

.. include:: gitdist-dist-help-aliases.txt
   :literal:

gitdist --dist-help=default-branch
++++++++++++++++++++++++++++++++++

.. include:: gitdist-dist-help-default-branch.txt
   :literal:

gitdist --dist-help=move-to-base-dir
++++++++++++++++++++++++++++++++++++

.. include:: gitdist-dist-help-move-to-base-dir.txt
   :literal:

gitdist --dist-help=usage-tips
++++++++++++++++++++++++++++++

.. include:: gitdist-dist-help-usage-tips.txt
   :literal:

gitdist --dist-help=script-dependencies
+++++++++++++++++++++++++++++++++++++++

.. include:: gitdist-dist-help-script-dependencies.txt
   :literal:

gitdist --dist-help=all
+++++++++++++++++++++++++++++++++++++++

.. include:: gitdist-dist-help-all.txt
   :literal:


.. _snapshot-dir.py:

snapshot-dir.py --help
----------------------

Below is a snapshot of the output from ``snapshot-dir.py --help``.  For more
details on the usage of ``snapshot-dir.py``, specifically for snapshotting the
`<projectDir>/cmake/tribits/`_ directory, see `TriBITS directory
snapshotting`_.

.. include:: snapshot-dir-help.txt
   :literal:


checkin-test.py --help
----------------------

Below is a snapshot of the output from ``checkin-test.py --help``.  This
``--help`` output contains a lot of information about the recommended
development workflow (mostly related to pushing commits) and outlines a number
of different use cases for using the tool.

.. include:: checkin-test-help.txt
   :literal:


.. _is_checkin_tested_commit.py:

is_checkin_tested_commit.py --help
----------------------------------

Below is a snapshot of the output from ``is_checkin_tested_commit.py --help``.
For more details see `Using Git Bisect with checkin-test.py workflows`_.

.. include:: is_checkin_tested_commit.txt
   :literal:


.. _get-tribits-packages-from-files-list.py:

get-tribits-packages-from-files-list.py --help
----------------------------------------------

Below is a snapshot of the output from
``get-tribits-packages-from-files-list.py --help``.  For more details see
`TriBITS Project Dependencies XML file and tools`_.

.. include:: get-tribits-packages-from-files-list.txt
   :literal:


.. _get-tribits-packages-from-last-tests-failed.py:

get-tribits-packages-from-last-tests-failed.py --help
-----------------------------------------------------

Below is a snapshot of the output from
``get-tribits-packages-from-last-tests-failed.py --help``.  For more details
see `TriBITS Project Dependencies XML file and tools`_.

.. include:: get-tribits-packages-from-last-tests-failed.txt
   :literal:


.. _filter-packages-list.py:

filter-packages-list.py --help
------------------------------

Below is a snapshot of the output from ``filter-packages-list.py --help``.
For more details see `TriBITS Project Dependencies XML file and tools`_.

.. include:: filter-packages-list.txt
   :literal:


install_devtools.py --help
--------------------------

Below is a snapshot of the output from ``install_devtools.py --help``.

.. include:: install_devtools-help.txt
   :literal:



.. ***
.. *** Common references
.. **

.. Common references to TribitsBuildReference document

.. NOTE: These references
.. are for when published in the same directory using public_docs.sh

.. _<Project>BuildReference: TribitsBuildReference.html

.. _TribitsBuildReference: `<Project>BuildReference`_

.. _Selecting the list of packages to enable: TribitsBuildReference.html#selecting-the-list-of-packages-to-enable

.. _Enabling extra repositories with add-on packages: TribitsBuildReference.html#enabling-extra-repositories-with-add-on-packages

.. _Getting set up to use CMake: TribitsBuildReference.html#getting-set-up-to-use-cmake

.. _Setting install RPATH: TribitsBuildReference.html#setting-install-rpath

.. _Dashboard Submissions: TribitsBuildReference.html#dashboard-submissions

.. _<Project>_VERBOSE_CONFIGURE: TribitsBuildReference.html#project-verbose-configure

.. _<Project>_ENABLE_DEVELOPMENT_MODE: TribitsBuildReference.html#project-enable-development-mode

.. _<Project>_CONFIGURE_OPTIONS_FILE: TribitsBuildReference.html#project-configure-options-file

.. _<Project>_EXTRAREPOS_FILE: TribitsBuildReference.html#project-extrarepos-file

.. _${PROJECT_NAME}_EXTRAREPOS_FILE: `<Project>_EXTRAREPOS_FILE`_

.. _<Project>_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE: TribitsBuildReference.html#project-enable-known-external-repos-type

.. _<Project>_PRE_REPOSITORIES: TribitsBuildReference.html#project-pre-repositories

.. _<Project>_EXTRA_REPOSITORIES: TribitsBuildReference.html#project-extra-repositories

.. _<Project>_GENERATE_REPO_VERSION_FILE: TribitsBuildReference.html#generating-a-project-repo-version-file

.. _Creating a tarball of the source tree: TribitsBuildReference.html#creating-a-tarball-of-the-source-tree

.. _Enabling support for an optional Third-Party Library (TPL): TribitsBuildReference.html#enabling-support-for-an-optional-third-party-library-tpl

.. _${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE: TribitsBuildReference.html#project-configure-options-file

.. _Outputting package dependency information: TribitsBuildReference.html#outputting-package-dependency-information

.. _${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE: TribitsBuildReference.html#outputting-package-dependency-information

.. _${PROJECT_NAME}_TRACE_FILE_PROCESSING: TribitsBuildReference.html#project-trace-file-processing

.. _Setting test timeouts at configure time: TribitsBuildReference.html#dart-testing-timeout

.. _${PROJECT_NAME}_SCALE_TEST_TIMEOUT: TribitsBuildReference.html#project-scale-test-timeout-testing-timeout

.. _Overriding test timeouts: TribitsBuildReference.html#overriding-test-timeouts

.. _make dashboard: TribitsBuildReference.html#dashboard-submissions

.. _Setting the install prefix at configure time: TribitsBuildReference.html#setting-the-install-prefix-at-configure-time

.. _TRIBITS_2ND_CTEST_DROP_SITE: TribitsBuildReference.html#tribits-2nd-ctest-drop-site

.. _TRIBITS_2ND_CTEST_DROP_LOCATION: TribitsBuildReference.html#tribits-2nd-ctest-drop-location

.. _Enabling extra repositories through a file: TribitsBuildReference.html#enabling-extra-repositories-through-a-file

.. Common references to the TribitsOverview document

.. _TriBITS Overview: ../overview/TribitsOverview.pdf


.. Common references to the TribitsLifecycleModel document

.. _TriBITS Lifecycle Model: ../lifecycle_model/TribitsLifecycleModel.pdf

..  Other common references

.. _CTest documentation: http://www.cmake.org/Wiki/CMake/Testing_With_CTest

.. Common references to raw CMake commands:

.. _CONFIGURE_FILE(): https://cmake.org/cmake/help/v3.10/command/configure_file.html

.. Other references

.. _Software Framework: https://en.wikipedia.org/wiki/Software_framework

..  LocalWords:  TribitsOverview TribitsLifecycleModel Lifecycle CMakeLists
..  LocalWords:  TriBITS Subpackage subpackage Subpackages subpackages TPL TPLs Kitware
..  LocalWords:  CMake cmake CTest ctest CDash CPack WithSubpackages WithSubpackagesA
..  LocalWords:  WithSubpackagesB WithSubpackagesC executables FOREACH ENDFOREACH

