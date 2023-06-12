TriBITS System Data Structures
------------------------------

This section describes the global CMake variables that make up the
data-structures and the macros/functions that create them that define the
TriBITS package dependency system.  All of these variables all exist at the
base project-level ``CMakeLists.txt`` file and are typically not cache
variables (and therefore are recomputed on every reconfigure and can therefore
accommodate changing enables/disables without a full reconfigure from
scratch).  These variables define a graph of external packages/TPLs
(i.e. pre-built and found out on the system) and internal packages
(i.e. buildable CMake packages).  This information is meant for maintainers of
the TriBITS system itself and should **not** need to be known by `TriBITS
Project Developers`_ or even `TriBITS Project Architects`_.


TriBITS naming conventions
++++++++++++++++++++++++++

Before describing the TriBITS package architecture data-structures and the
macros/functions that create and manipulate those data-structures in detail,
first we define some naming conventions for TriBITS macros/function and
variables.

First, the term "Package" is overloaded in the context of TriBITS (and in
CMake for that matter).  In the context of TriBITS, here are the **different
types of entities called a "Package"**:

* **TriBITS External packages/TPLs** listed in the `<repoDir>/TPLsList.cmake`_
  file and require a `<tplDefsDir>/FindTPL<tplName>.cmake`_ file.

* **TriBITS Top-level internal packages** directly listed in the
  `<repoDir>/PackagesList.cmake`_ file

* **TriBITS Subpackages** for a top-level package (defined by the
  ``SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS`` argument in
  `tribits_package_define_dependencies()`_ command called in the
  `<packageDir>/cmake/Dependencies.cmake`_ file)

* **Native CMake packages** which are found using ``find_package()`` and
  require a ``<Package>Config.cmake`` package config file or a
  ``Find<Package>.cmake`` find module file

* **CMake Packaging** of a CMake project which means breaking it into
  components and provide various source and binary distributions.  (A TriBITS
  top-level package maps to a CMake packaging "component".)

Also, there are the collection of all three of the TriBITS-related "packages".

To try to avoid ambiguity, we define the following identifiers that appear in
TriBITS variable, macro, and function names:

* **TPL** or **EXTERNAL_PACKAGE**: Variable or function specifically for
  processing TriITS external packages/TPLs.

* **TOPLEVEL_PACKAGE**: Variable or a function specifically for a top-level
  package (or sometimes if there is no ambiguity just **PACKAGE**).

* **SUBPACKAGE**: Variable or a function specifically for a subpackage (and
  **SUBPACKAGE_FULLNAME** if it is the full name of a subpackage which
  includes the top-level parent package name
  ``${PARENT_PACKAGE_NAME}${SUBPACKAGE}``).

* **PACKAGE**: Variable or function that could pertain to an external
  package/TPL, or a top-level package, or a subpackage.

TriBITS uses the follow general case naming conventions for variables, macros,
functions and module files:

* ``ProperNoun_UPPER_CASE`` is generally used for global and cache variables.
  The proper nouns using CamelCase include the names of TriBITS entities like
  Projects, Repositories, and Packages.  UPPER_CASE is used for non-proper
  nouns. (E.g. ``Trilinos_SOURCE_DIR``)

* ``camelCaseVariableName`` or ``lower_case_variable_name`` is generally used
  for local variable names.

* ``tribits_lower_case_name()`` is generally used for TriBITS functions and
  macros. (E.g. ``trilinos_package_define_dependencies()``)

* ``TribitsModuleName.cmake`` is generally used for TriBITS CMake module file
  names. (E.g. ``TribitsAddTest.cmake``)


User-level TriBITS Project, Repository, Package and Subpackage Variables
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The user-level variables that define a TriBITS Project, Repository, Package
and Subpackage are listed in:

* `TriBITS Project Core Variables`_
* `TriBITS Repository Core Variables`_
* `TriBITS Package Core Variables`_
* `TriBITS Subpackage Core Variables`_
* `TriBITS External Package/TPL Core Variables`_

These are variables that can be accessed by `TriBITS Project Developers`_ but
are also used in the internal implementation of TriBITS functionality.


Lists of external and internal packages
+++++++++++++++++++++++++++++++++++++++

List of non-cache top-level project variables:

* `${PROJECT_NAME}_DEFINED_TPLS`_: List of all defined external packages/TPLs

  * Size `${PROJECT_NAME}_NUM_DEFINED_TPLS`_

* `${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES`_: List of all defined
  internal top-level packages

  * Size `${PROJECT_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES`_

* `${PROJECT_NAME}_DEFINED_TOPLEVEL_PACKAGES`_: List of all defined external
  packages/TPLs and top-level internal packages

  * Size `${PROJECT_NAME}_NUM_DEFINED_TOPLEVEL_PACKAGES`_

* `${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES`_: List of all defined internal
  top-level packages and subpackages

  * Size `${PROJECT_NAME}_NUM_DEFINED_INTERNAL_PACKAGES`_

* `${PROJECT_NAME}_DEFINED_PACKAGES`_: List of all defined external
  packages/TPLs, internal top-level packages, and subpackages (size
  `${PROJECT_NAME}_NUM_DEFINED_PACKAGES`_)

  * Size ${PROJECT_NAME}_NUM_DEFINED_PACKAGES

* `${PROJECT_NAME}_ENABLED_PACKAGES`_: Subset of all enabled packages from
  ``${PROJECT_NAME}_DEFINED_PACKAGES``

  * Size `${PROJECT_NAME}_NUM_ENABLED_PACKAGES`_

All of the above list variables are sorted in a valid dependency ordering in
that any upstream dependent packages are listed before a given package in
these lists.  After these variables have been set in the macro
`tribits_read_all_project_deps_files_create_deps_graph()`_, they should
considered to be constant and **not** modified.

These variables are described in more detail below.

.. _${PROJECT_NAME}_DEFINED_TPLS:

The original list of all defined external packages (TPLs) read from the
processed `<repoDir>/TPLsList.cmake`_ files is given in the list variable::

  ${PROJECT_NAME}_DEFINED_TPLS

.. _${PROJECT_NAME}_NUM_DEFINED_TPLS:

with size::

  ${PROJECT_NAME}_NUM_DEFINED_TPLS

.. _${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES:

The original list of all defined internal top-level packages read in from the
processed `<repoDir>/PackagesList.cmake`_ files is given in the list
variable::

  ${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES

.. _${PROJECT_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES:

with size::

  ${PROJECT_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES

In this list, a defined internal TriBITS Package (i.e. a package that can be
built from source) will have `${PACKAGE_NAME}_SOURCE_DIR`_ ``!= ""`` while a
defined external package/TPL will have a non-empty `${PACKAGE_NAME}_FINDMOD`_
``!= ""``.

.. _${PROJECT_NAME}_DEFINED_TOPLEVEL_PACKAGES:

The full list of defined external packages/TPLs and top-level internal
packages (i.e. TriBITS top-level packages) (**not** including subpackages) is
stored in the project-level non-cache list variable::

  ${PROJECT_NAME}_DEFINED_TOPLEVEL_PACKAGES

.. _${PROJECT_NAME}_NUM_DEFINED_TOPLEVEL_PACKAGES:

with size::

  ${PROJECT_NAME}_NUM_DEFINED_TOPLEVEL_PACKAGES

The first set of elements in this list are the defined external packages/TPLs
that are read in from the `<repoDir>/TPLsList.cmake`_ files from each
processed TriBITS repository, in order.  This is followed by the set of
internal packages (TriBITS packages) that are defined in the
`<repoDir>/PackagesList.cmake`_ files from each processed TriBITS repository,
read in order.  This list does **not** include any subpackages.

Note that some of the packages listed in
`${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES`_ may actually be treated
as external packages and not build from source code and instead will be found
on the system as pre-built/pre-installed packages using
``find_package(<PackageName>)``.  The final decision for if a package is
treated as an internal or external package is determined by the variable::

  ${PACKAGE_NAME}_PACKAGE_BUILD_STATUS=[INTERNAL|EXTERNAL]

which gets set using various criteria as described in section `Determining if
a package is internal or external`_.  This variable determines what
pre-built/pre-installed packages must be found out on the system if enabled
and what internal packages need to be built if enabled.

The set of external packages, internal top-level packages, and internal
sub-packages are just called the list of "Packages".  When the term "Packages"
is used without an adjective, it is usually meant in this more general
context.

.. _${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES:

The set of all of the defined internal top-level packages **and subpackages**
is given by the non-cache project-level list variable::

  ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES

.. _${PROJECT_NAME}_NUM_DEFINED_INTERNAL_PACKAGES:

with the size::

  ${PROJECT_NAME}_NUM_DEFINED_INTERNAL_PACKAGES

.. _${PROJECT_NAME}_DEFINED_PACKAGES:

The set of all of the defined external packages/TPLs, internal top-level
packages **and subpackages** is given by the non-cache project-level list
variable::

  ${PROJECT_NAME}_DEFINED_PACKAGES

.. _${PROJECT_NAME}_NUM_DEFINED_PACKAGES:

with the size::

  ${PROJECT_NAME}_NUM_DEFINED_PACKAGES

These data-structures as well as the package dependencies graph is built up in
the macro `tribits_read_all_project_deps_files_create_deps_graph()`_ with the
call graph described in the section `Function call tree for constructing
package dependency graph`_.  These data-structures don't consider what
packages are actually enabled or disabled.

The enable/disable logic (given an initial set of enables and disables) is
applied in the macro `tribits_adjust_package_enables()`_.  Once all of this
logic has been applied, several lists of enabled and non-enabled packages are
computed.

The list of enabled internal **top-level** packages is given in the non-cache
project-level list variable::

  ${PROJECT_NAME}_ENABLED_INTERNAL_TOPLEVEL_PACKAGES

with size::

  ${PROJECT_NAME}_NUM_ENABLED_INTERNAL_TOPLEVEL_PACKAGES

The list of enabled external packages/TPLs and internal **top-level** packages
is given in the non-cache project-level list variable::

  ${PROJECT_NAME}_ENABLED_TOPLEVEL_PACKAGES

with size::

  ${PROJECT_NAME}_NUM_ENABLED_TOPLEVEL_PACKAGES

.. _${PROJECT_NAME}_ENABLED_PACKAGES:

The list of enabled external packages/TPLs, internal **top-level and
subpackages** is given in the non-cache project-level list variable::

  ${PROJECT_NAME}_ENABLED_PACKAGES

.. _${PROJECT_NAME}_NUM_ENABLED_PACKAGES:

with size::

  ${PROJECT_NAME}_NUM_ENABLED_PACKAGES


Variables defining the package dependencies graph
+++++++++++++++++++++++++++++++++++++++++++++++++

TriBITS sets up the following project-level non-cache variables that define
the dependencies for each external package/TPL and internal package:

  .. _${PACKAGE_NAME}_LIB_DEFINED_DEPENDENCIES:

  ``${PACKAGE_NAME}_LIB_DEFINED_DEPENDENCIES``

    The list of all **defined direct** required and optional upstream external
    package/TPL and internal package dependencies, regardless if they are
    enabled or not.  To determine if a given direct upstream package
    ``<depPkg>`` in this list is enabled or not for this package
    ``${PACKAGE_NAME}``, check the value of
    ``${PACKAGE_NAME}_ENABLE_<depPkg>``.  NOTE: The variables
    ``${PACKAGE_NAME}_ENABLE_<depPkg>`` will be set even for required upstream
    packages to allow for uniform loops involving required and optional
    upstream dependencies.  (And for a parent package with subpackages, it is
    possible for a required subpackage to **not** be enabled and for
    ``${PACKAGE_NAME}_ENABLE_<depPkg>`` to be ``OFF`` as explained in
    `Subpackage enable does not auto-enable the parent package`_.)  This list
    will be set regardless of if the package ``${PACKAGE_NAME}`` is enabled or
    not.

  .. _${PACKAGE_NAME}_LIB_ENABLED_DEPENDENCIES:

  ``${PACKAGE_NAME}_LIB_ENABLED_DEPENDENCIES``

    List of all **enabled direct** required and optional upstream external
    package/TPL and internal package dependencies.  This is strict subset of
    `${PACKAGE_NAME}_LIB_DEFINED_DEPENDENCIES`_ (i.e. all of the ``<depPkg>``
    items in this list will have ``${PACKAGE_NAME}_ENABLE_<depPkg>`` set to
    ``ON``).

  .. _${PACKAGE_NAME}_LIB_DEP_REQUIRED_<depPkg>:

  ``${PACKAGE_NAME}_LIB_DEP_REQUIRED_<depPkg>``

    Is ``TRUE`` if the entry ``<depPkg>`` in
    `${PACKAGE_NAME}_LIB_DEFINED_DEPENDENCIES`_ or
    `${PACKAGE_NAME}_LIB_ENABLED_DEPENDENCIES`_ is a required LIB dependency
    and is ``FALSE`` if it is only an optional LIB dependency.

  .. _${PACKAGE_NAME}_TEST_DEFINED_DEPENDENCIES:

  ``${PACKAGE_NAME}_TEST_DEFINED_DEPENDENCIES``

    This list of all **define direct** extra package test required and
    optional upstream external package/TPL and internal package dependencies.
    This list is set regardless if the package ``${PACKAGE_NAME}`` is enabled
    or not.

  .. _${PACKAGE_NAME}_TEST_ENABLED_DEPENDENCIES:

  ``${PACKAGE_NAME}_TEST_ENABLED_DEPENDENCIES``

    The list of all **enabled direct** extra required and optional upstream
    external package/TPL and internal package dependencies.  This is a strict
    subset of `${PACKAGE_NAME}_TEST_DEFINED_DEPENDENCIES`_.

  .. _${PACKAGE_NAME}_TEST_DEP_REQUIRED_<depPkg>:

  ``${PACKAGE_NAME}_TEST_DEP_REQUIRED_<depPkg>``

    Is ``TRUE`` if the entry ``<depPkg>`` in
    `${PACKAGE_NAME}_TEST_DEFINED_DEPENDENCIES`_ or
    `${PACKAGE_NAME}_TEST_ENABLED_DEPENDENCIES`_ is a required TEST dependency
    and is ``FALSE`` if it is only an optional TEST dependency.  For the sake
    of simplicity and generality,
    ``${PACKAGE_NAME}_TEST_DEP_REQUIRED_<depPkg>`` will also be set to
    ``TRUE`` or ``FALSE`` for ``<depPkg>`` in the lists
    `${PACKAGE_NAME}_LIB_DEFINED_DEPENDENCIES`_ or
    `${PACKAGE_NAME}_LIB_ENABLED_DEPENDENCIES`_ because a LIB dependency is
    also implicitly a TEST dependency.

NOTE: The same upstream package ``<depPkg>`` can be included in both the lists
`${PACKAGE_NAME}_LIB_DEFINED_DEPENDENCIES`_ and
`${PACKAGE_NAME}_TEST_DEFINED_DEPENDENCIES`_ if ``<depPkg>`` is optional in
the former but required in the latter (which is a valid situation if you think
about it as a package that may be optional for the lib(s) of a package is
required by the tests for a package).  (Otherwise, duplicate entries will be
removed from the list ``${PACKAGE_NAME}_TEST_DEFINED_DEPENDENCIES``.)

NOTE: Having flat lists containing both optional and required dependencies
with the bool variables ``${PACKAGE_NAME}_[LIB|TEST]_DEP_REQUIRED_<depPkg>``
defining which entries are required or optional is modeled after the CMake
standard for handing the ``COMPONENTS`` and ``OPTIONAL_COMPONENTS`` arguments
to ``find_package()`` in that it passes that info to the
``<Package>Config.cmake`` file as the single list variable
``${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS`` and the bool vars
``${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED_<comp>``.

Given the above upstream dependency list variables, the following derived list
variables are then constructed which provide navigation from a package to its
downstream/forward dependent packages:

  .. _${PACKAGE_NAME}_FORWARD_LIB_DEFINED_DEPENDENCIES:

  ``${PACKAGE_NAME}_FORWARD_LIB_DEFINED_DEPENDENCIES``

    For a given package ``${PACKAGE_NAME}``, lists the names of all of the
    forward packages ``<fwdDepPkg>`` that list this package in their
    ``<fwdDepPkg>_LIB_DEFINED_PACKAGES`` variables.

  .. _${PACKAGE_NAME}_FORWARD_TEST_DEFINED_DEPENDENCIES:

  ``${PACKAGE_NAME}_FORWARD_TEST_DEFINED_DEPENDENCIES``

    For a given package ``${PACKAGE_NAME}``, lists the names of all of the
    forward packages ``<fwdDepPkg>`` that list this package in their
    ``<fwdDepPkg>_TEST_DEFINED_PACKAGES`` variables.


User enable/disable cache variables
+++++++++++++++++++++++++++++++++++

The following variables can be set by the user to determine what packages get
enabled or disabled::
  
  ${PROJECT_NAME}_ENABLE_ALL_PACKAGES
  
  ${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES
  
  ${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES

  ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}
  
  ${PROJECT_NAME}_ENABLE_TESTS
  
  ${PROJECT_NAME}_ENABLE_EXAMPLES
  
  ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE_NAME}
  
  ${PACKAGE_NAME}_ENABLE_TESTS
  
  ${PACKAGE_NAME}_ENABLE_EXAMPLES

according to the rules described in `Package Dependencies and Enable/Disable
Logic`_.


Determining if a package is internal or external
++++++++++++++++++++++++++++++++++++++++++++++++

As mentioned above, some subset of initially internal packages listed in
`${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES`_ (which all have
``${PACKAGE_NAME}_SOURCE_DIR != ""``) may be chosen to be external packages.
Packages that could be built internally may be chosen to be treated as
external packages (and therefore located on the system using
``find_package()``) by setting::

  -D TPL_ENABLE_<PackageTreatedAsExternal>=ON

.. _${PACKAGE_NAME}_PACKAGE_BUILD_STATUS:
.. _${TPL_NAME}_PACKAGE_BUILD_STATUS:

The final status of whether a package is treated as an internal package or an
external package is provided by the variable::

  ${PACKAGE_NAME}_PACKAGE_BUILD_STATUS=[INTERNAL|EXTERNAL]

(NOT: The value of ``${PACKAGE_NAME}_PACKAGE_BUILD_STATUS`` is only changed
after all of the enable/disable dependency logic is complete.)

As a result, every other package upstream from any of these
``<PackageTreatedAsExternal>`` packages must therefore also be treated as
external packages automatically and will have
``${PACKAGE_NAME}_PACKAGE_BUILD_STATUS=EXTERNAL`` set accordingly.  Also, if
any subpackage is determined to be EXTERNAL, then the parent package of that
subpackage and every other peer subpackage will also be set to EXTERNAL.


Processing of external packages/TPLs and TriBITS-compliant external packages
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The processing of external packages/TPLs is influenced by whether the external
package is a regular TriBITS TPL (i.e with a ``FindTPL<tplName>.cmake``
modules) or is a TriBITS-compliant external package.  Here, a
**TriBITS-Compliant External Package** has a ``<tplName>Config.cmake`` file
that satisfies the following properties:

* Has the target ``<tplName>::all_libs`` which is a fully specified modern
  CMake target.
* Calls ``find_dependency()`` or the equivalent for all upstream packages that
  it depends on.
* Every upstream dependent package ``<UpstreamTpl>`` has the target
  ``<UpstreamTpl>::all_libs``.  (But a non-fully TriBITS-compliant external
  package need not define this for all of its upstream dependencies.)

That means that when calling ``find_package()`` for a fully TriBITS-compliant
external package, there is no need to worry about finding any of its upstream
dependent external packages.  That means that any external packages/TPLs
defined a TriBITS project which is upstream from a TriBITS-compliant external
package will be uniquely defined by calling ``find_package()`` on the most
downstream TriBITS-compliant external package that depends on it.  Therefore,
defining the external packages and their targets in this set of external
packages just involves calling ``find_package()`` on the terminal
TriBITS-compliant external packages (i.e. TriBITS-compliant external packages
that don't have any downstream dependencies that are external packages).  Then
the remaining subset of external packages/TPLs that don't have a downstream
TriBITS-compliant external package dependency will be defined as usual.
(However, as mentioned above, some of these are not fully TriBITS compliant and
don't fully define the ``<UpstreamTpl>::all_libs`` for all of their upstream
dependencies (see below).)

By having all fully TriBITS-compliant external packages, an external
dependency is never found more than once.

The variables that are set internally to define these different subsets of
external packages/TPLs are:

* ``<Package>_IS_TRIBITS_COMPLIANT``: Set the ``TRUE`` if the package
  ``<Package>`` provides the ``<Package>::all_libs`` by just calling
  ``find_package(<Package> CONFIG REQUIRED)``.

* ``<Package>_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE``: Set to
  ``TRUE`` if the external package/TPL should be (or was after the fact)
  defined by a downstream TriBITS-compliant external package.

An external package with ``<Package>_IS_TRIBITS_COMPLIANT=TRUE`` **AND**
``<Package>_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE=FALSE`` is the
one for which ``find_package(<Package> CONFIG REQUIRED)`` will be called and
does not have any downstream packages that are being treated as external
packages.  Also, ``find_package(<Package> CONFIG REQUIRED)`` will be called on
TriBITS-compliant external packages if ``<Package>::all_libs`` was not defined
by a downstream non fully TriBITS-compliant external package.

The variable ``<Package>_IS_TRIBITS_COMPLIANT`` is set right when the packages
are initially defined by reading in the various input files.  That is, all
initially internal packages that are listed in a
`<repoDir>/PackagesList.cmake`_ file will have
``<Package>_IS_TRIBITS_COMPLIANT=TRUE`` set while all external packages/TPLs
listed in a `<repoDir>/TPLsList.cmake`_ file will have
``<Package>_IS_TRIBITS_COMPLIANT=FALSE`` set (except for those tagged with
``TRIBITS_PKG`` which will have ``<Package>_IS_TRIBITS_COMPLIANT=FALSE`` set).

The processing of external packages/TPLs is done in two loops:

* The first loop over external packages/TPLs will be those external
  packages/TPLs that have ``<Package>_IS_TRIBITS_COMPLIANT=TRUE`` **OR**
  ``<Package>_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE=TRUE``.  And we
  only call ``find_package()`` for those TriBITS-compliant external packages
  that have ``<Package>_IS_TRIBITS_COMPLIANT=TRUE`` **AND** don't have
  ``<Package>::all_libs`` already defined.  This is a reverse loop to give an
  opportunity for downstream TriBITS-compliant external packages to define
  their upstream external packages/TPLs.  NOTE: If
  ``<Package>_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE`` was set to
  ``TRUE`` before this loop starts, it will be set to ``FALSE`` for
  non-TriBITS-compliant external packages
  (i.e. ``<Package>_IS_TRIBITS_COMPLIANT=FALSE``).

* The second loop processes remaining external packages/TPLs that where not
  defined by a downstream TriBITS-compliant external package in the first
  loop.  These are all TriBITS TPLs for which
  ``<Package>_IS_TRIBITS_COMPLIANT=FALSE`` **AND**
  ``<Package>_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE=FALSE`` **AND**
  for which ``<Package>_FINDMOD`` is not empty and is not ``TRIBITS_PKG``.

For more details, see the implementation in `tribits_process_enabled_tpls()`_.


Other package-related variables
+++++++++++++++++++++++++++++++

The following global internal cache variables are used to provide more
information about a given internal package:

  ``${PACKAGE_NAME}_LIBRARIES``

    Defines list of *only* the libraries associated with the given
    (sub)package and does *not* list libraries in upstream packages.  Linkages
    to upstream packages is taken care of with calls to
    target_link_libraries(...) and the dependency management system in CMake
    takes care of adding these to various link lines as needed (this is what
    CMake does well).  However, when a package has no libraries of its own
    (which is often the case for packages that have subpackages, for example),
    then this list of libraries will contain the libraries to the direct
    dependent upstream packages in order to allow the chain of dependencies to
    be handled correctly in downstream packages and executables in the same
    package.  In this case, ${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES will be
    false.  The primary purpose of this variable is to passe to
    target_link_libraries(...) by downstream libraries and executables.

  ``${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES``

    Will be true if a package has native libraries.  Otherwise, it will be
    false.  This information is used to build export makefiles to avoid
    duplicate libraries on the link line.

  ``${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES``

    Lists out, in order, all of the enabled upstream packages that the given
    package depends on and support that package is enabled in the given
    package.  This is only computed if
    ${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES=ON.  NOTE: This list
    does *not* include the package itself.  This list is created after all of
    the enable/disable logic is applied.

  ``${PARENT_PACKAGE_NAME}_LIB_TARGETS``

    Lists all of the library targets for this package only that are as part of
    this package added by the `tribits_add_library()`_ function.  This is used
    to define a target called ${PACKAGE_NAME}_libs that is then used by
    `tribits_ctest_driver()`_ in the package-by-package mode.  If a package
    has no libraries, then the library targets for all of the immediate
    upstream direct dependent packages will be added.  This is needed for the
    chain of dependencies to work correctly.  Note that subpackages don't have
    this variable defined for them.

  ``${PARENT_PACKAGE_NAME}_ALL_TARGETS``

    Lists all of the targets associated with this package.  This includes all
    libraries and tests added with `tribits_add_library()`_ and
    `tribits_add_executable()`_.  If this package has no targets (no libraries
    or executables) this this will have the dependency only on
    ${PARENT_PACKAGE_NAME}_libs.  Note that subpackages don't have this
    variable defined for them.


Function call tree for constructing package dependency graph
------------------------------------------------------------

Below is the CMake macro and function call graph for constructing the packages
lists and dependency data-structures described above.

| `tribits_read_all_project_deps_files_create_deps_graph()`_
|   `tribits_read_defined_external_and_internal_toplevel_packages_lists()`_
|     Foreach ``<repoDir>`` in ``${PROJECT_NAME}_ALL_REPOSITORIES``:
|       ``include(`` `<repoDir>/TPLsList.cmake`_ ``)``
|       `tribits_process_tpls_lists()`_
|       ``include(`` `<repoDir>/PackagesList.cmake`_ ``)``
|       `tribits_process_packages_and_dirs_lists()`_
|   `tribits_read_deps_files_create_deps_graph()`_
|     `tribits_process_all_repository_deps_setup_files()`_
|       Foreach ``<repoDir>`` in ``${PROJECT_NAME}_ALL_REPOSITORIES``:
|         ``include(`` `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ ``)``
|     `tribits_process_project_dependency_setup_file()`_
|       ``include(``  `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_ ``)``
|     `tribits_read_all_package_deps_files_create_deps_graph()`_
|       Foreach  ``EXTERNAL_PACKAGE`` in ``${PROJECT_NAME}_DEFINED_TPLS``:
|         `tribits_read_external_package_deps_files_add_to_graph()`_
|       Foreach ``TOPLEVEL_PACKAGE`` in ``${PACKAGE_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES``:
|         `tribits_read_toplevel_package_deps_files_add_to_graph()`_
|           `tribits_prep_to_read_dependencies()`_
|           ``include(`` `<packageDir>/cmake/Dependencies.cmake`_ ``)``
|           `tribits_assert_read_dependency_vars()`_
|           `tribits_save_off_dependency_vars()`_
|           `tribits_parse_subpackages_append_packages_add_options()`_
|           `tribits_read_package_subpackage_deps_files_add_to_graph()`_
|             Foreach ``SUBPACKAGE``:
|               `tribits_read_subpackage_deps_file_add_to_graph()`_
|                 `tribits_prep_to_read_dependencies()`_
|                 ``include(`` `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_ ``)``
|                 `tribits_assert_read_dependency_vars()`_
|                 `tribits_process_package_dependencies_lists()`_
|                   See same call stack for this macro as shown below
|           `tribits_read_back_dependencies_vars()`_
|           `tribits_process_package_dependencies_lists()`_
|             `tribits_set_dep_packages()`_
|               `tribits_abort_on_self_dep()`_
|               `tribits_abort_on_missing_package()`_
|             `tribits_append_forward_dep_packages()`_
|               `tribits_abort_on_missing_package()`_

..  LocalWords:  acyclic TriBITS SUBPACKAGES CTEST subpackages buildable TPLs TPLS
..  LocalWords:  Subpackage CMake CMakeLists
