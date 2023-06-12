CMake Language Overview and Gotchas
-----------------------------------

TriBITS removes a lot of the boiler plate code needed to write a CMake
project.  As a result, many people can come into a project that uses TriBITS
and quickly start to contribute by adding new source files, adding new
libraries, adding new tests, and even adding new TriBITS packages and external
packages/TPLs; all without really having learned anything about CMake.  Often
one can use existing example CMake code as a guide and be successful using
basic functionality. As long as nothing out of the ordinary happens, many
people can get along just fine in this mode for a time.

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

  some_func(A B C)
  some_func("A" "B" "C")
  some_func("A;B;C")

However, the above is *not* the same as::

  some_func("A B C")

which just passes in a single element with value [``A B C``].  Raw quotes in
CMake basically escape the interpretation of space characters as array element
boundaries.  Quotes around arguments with no spaces does nothing (as seen
above, except for the interpretation as variable names in an ``if()``
statement).  In order to get a quote char [``"``] into string, you must escape
it as::

  some_func(\"A\")

which passes an array with the single argument [``\"A\"``].

Variables are set using the built-in CMake ``set()`` command that just takes
string arguments like::

  set(SOME_VARIABLE "some_value")

In CMake, the above is identical, in every way, to::

  set(SOME_VARIABLE some_value)
  set("SOME_VARIABLE";"some_value")
  set("SOME_VARIABLE;some_value")

The function ``set()`` simply interprets the first argument to as the name of
a variable to set in the local scope.  Many other built-in and user-defined
CMake functions work the same way.  That is, some of the string arguments are
interpreted as the names of variables.  There is no special language feature
that interprets them as variables (except in an ``if()`` statement).

However, CMake appears to parse arguments differently for built-in CMake
control structure functions like ``foreach()`` and ``if()`` and does not just
interpret them as a string array.  For example::

  foreach (SOME_VAR "a;b;c")
    message("SOME_VAR='${SOME_VAR}'")
  endforeach()

prints ```SOME_VAR='a;b;c'`` instead of printing ``SOME_VAR='a'`` followed by
``SOME_VAR='b'``, etc., as you would otherwise expect.  Therefore, this simple
rule for the handling of function arguments as string arrays does not hold for
CMake logic control commands.  Just follow the CMake documentation for these
control structures (i.e. see ``cmake --help-command if`` and ``cmake
--help-command foreach``).

CMake offers a rich assortment of built-in commands for doing all sorts of
things.  Two of these are the built-in ``macro()`` and the ``function()``
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

   set(SOME_VAR a ${SOME_UNDEFINED_VAR} c)

(where ``SOME_UNDEFINED_VAR`` is an undefined variable) produces
``SOME_VAR='a;c'`` and *not* ``'a;;c'``!  The same thing occurs when an empty
variable is de-references such as with::

   set(EMPTY_VAR "")
   set(SOME_VAR a ${EMPTY_VAR} c)

which produces ``SOME_VAR='a;c'`` and *not* ``'a;;c'``.  In order to always
produce an element in the array even if the variable is empty, one must quote
the argument as with::

   set(EMPTY_VAR "")
   set(SOME_VAR a "${EMPTY_VAR}" c)

which produces ``SOME_VAR='a;;c'``, or three elements as one might assume.

This is a common error that people make when they call CMake functions
(built-in or TriBITS-defined) involving variables that might be undefined or
empty.  For example, for the macro::

   macro(some_macro  A_ARG  B_ARG  C_ARG)
      ...
   endmacro()

if someone tries to call it with (misspelled variable?)::

  some_macro(a ${SOME_OHTER_VAR} c)

and if ``SOME_OHTER_VAR=""`` or if it is undefined, then CMake will error out
with the error message saying that the macro ``some_macro()`` takes 3
arguments but only 2 were provided.  If a variable might be empty but that is
still a valid argument to the command, then it must be quoted as::

  some_macro(a "${SOME_OHTER_VAR}" c)

Related to this problem is that if you misspell the name of a variable in a
CMake ``if()`` statement like::

   if (SOME_VARBLE)
     ...
   endif()

then it will always be false and the code inside the if statement will never
be executed!  To avoid this problem, use the utility function
`assert_defined()`_ as::

   assert_defined(SOME_VARBLE)
   if (SOME_VARBLE)
     ...
   endif()

In this case, the misspelled variable would be caught.

While on the subject of ``if()`` statements, CMake has a strange convention.
When you say::

  if (SOME_VAR)
    do_something()
  endif()

then ``SOME_VAR`` is interpreted as a variable and will be considered true and
``do_something()`` will be called if ``${SOME_VAR}`` does *not* evaluate to
``0``, ``OFF``, ``NO``, ``FALSE``, ``N``, ``IGNORE``, ``""``, or ends in the
suffix ``-NOTFOUND``.  How about that for a true/false rule!  To be safe, use
``ON/OFF`` and ``TRUE/FALSE`` pairs for setting variables.  Look up native
CMake documentation on ``if()`` for all the interesting details and all the
magical things it can do.

**WARNING:** If you mistype ``"ON"`` as ``"NO"``, it evaluates to
``FALSE``/``OFF``!  (That is a fun defect to track down!)

CMake language behavior with respect to case sensitivity is also strange:

* Calls of built-in and user-defined macros and functions is *case
  insensitive*!  That is ``set(...)``, ``set(...)``, ``set()``, and all other
  combinations of upper and lower case characters for 'S', 'E', 'T' all call
  the built-in ``set()`` function.  The convention in TriBITS is to use
  ``lower_case_with_underscores()`` for functions and macros.

* However, the names of CMake (local or cache/global) variables are *case
  sensitive*!  That is, ``SOME_VAR`` and ``some_var`` are *different*
  variables.  Built-in CMake variables tend use all caps with underscores
  (e.g. ``CMAKE_CURRENT_SOURCE_DIR``) but other built-in CMake variables tend
  to use mixed case with underscores (e.g. ``CMAKE_Fortran_FLAGS``).  TriBITS
  tends to use a similar naming convention where project-level and cache
  variables have mostly upper-case letters except for parts that are proper
  nouns like the project, package or external package/TPL name
  (e.g. ``TribitsExProj_TRIBITS_DIR``, ``TriBITS_SOURCE_DIR``,
  ``Boost_INCLUDE_DIRS``).  Local variables and function/macro parameters can
  use camelCase or lower_case_with_underscores.

I don't know of any other programming language that uses different case
sensitivity rules for variables and functions.  However, because we must parse
macro and function arguments when writing user-defined macros and functions,
it is a good thing that CMake variables are case sensitive.  Case insensitivity
would make it much harder and more expensive to parse argument lists that take
keyword-based arguments.

Other mistakes that people make result from not understanding how CMake scopes
variables and other entities.  CMake defines a global scope (i.e. "cache"
variables) and several nested local scopes that are created by
``add_subdirectory()`` and entering functions.  See `dual_scope_set()`_ for a
short discussion of these scoping rules.  And it is not just variables that
can have local and global scoping rules.  Other entities, like defines set
with the built-in command ``add_definitions()`` only apply to the local scope
and child scopes.  That means that if you call ``add_definitions()`` to set a
define that affects the meaning of a header-file in C or C++, for example,
that definition will *not* carry over to a peer subdirectory and those
definitions will not be set (see warning in `Miscellaneous Notes
(tribits_add_library())`_).
