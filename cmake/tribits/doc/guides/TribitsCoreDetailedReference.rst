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
``tribits_define_global_options_and_define_extra_repos()`` in the TriBITS
CMake code file ``TribitsGlobalMacros.cmake`` which gets called inside of the
`tribits_project()`_ macro.  That function and that file are the definitive
source the options that a TriBITS project takes and what the default values
are but we strive to document them here as well.  Many of these global options
(i.e. cache variables) such as ``${PROJECT_NAME}_<SOME_OPTION>`` allow the
project to define a default by setting a local variable
``${PROJECT_NAME}_<SOME_OPTION>_DEFAULT`` as::

  set(${PROJECT_NAME}_<SOME_OPTION>_DEFAULT <someDefault>)

either in its top-level ``CMakeLists.txt`` file or in its
``ProjectName.cmake`` file (depends on what variable it is as to where it
should be set).  If ``${PROJECT_NAME}_<SOME_OPTION>_DEFAULT`` is not set by
the project, then TriBITS provides a reasonable default value.  The TriBITS
code that uses these defaults for this looks like::

  if ("${${PROJECT_NAME}_<SOME_OPTION>_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_<SOME_OPTION>_DEFAULT <someDefault>)
  endif()

  advanced_set( ${PROJECT_NAME}_<SOME_OPTION>
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

* `${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE`_
* `${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES`_
* `${PROJECT_NAME}_C_Standard`_
* `${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS`_
* `${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_APPEND`_
* `${PROJECT_NAME}_CPACK_SOURCE_GENERATOR`_
* `${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE`_
* `${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES`_
* `${PROJECT_NAME}_ELEVATE_ST_TO_PT`_
* `${PROJECT_NAME}_ENABLE_CPACK_PACKAGING`_
* `${PROJECT_NAME}_ENABLE_CXX`_
* `${PROJECT_NAME}_ENABLE_C`_
* `${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE`_
* `${PROJECT_NAME}_ENABLE_Fortran`_
* `${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES`_
* `${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE`_
* `${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION`_
* `${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES`_
* `${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES`_
* `${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE`_
* `${PROJECT_NAME}_IMPORTED_NO_SYSTEM`_
* `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_
* `${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE`_
* `${PROJECT_NAME}_MAKE_INSTALL_GROUP_WRITABLE`_
* `${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE`_
* `${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS`_
* `${PROJECT_NAME}_REQUIRES_PYTHON`_
* `${PROJECT_NAME}_SET_INSTALL_RPATH`_
* `${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME`_
* `${PROJECT_NAME}_SKIP_INSTALL_PROJECT_CMAKE_CONFIG_FILES`_
* `${PROJECT_NAME}_TEST_CATEGORIES`_
* `${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS`_
* `${PROJECT_NAME}_TRACE_ADD_TEST`_
* `${PROJECT_NAME}_USE_GNUINSTALLDIRS`_
* `${PROJECT_NAME}_USES_PYTHON`_
* `DART_TESTING_TIMEOUT`_
* `CMAKE_INSTALL_RPATH_USE_LINK_PATH`_
* `MPI_EXEC_MAX_NUMPROCS`_
* `PythonInterp_FIND_VERSION`_
* `TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE`_

These options are described below.

.. _${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE:

**${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE**

  The CMake cache variable ``${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE`` is
  used to define how some invalid TriBITS usage checks are handled.  The valid
  values include 'FATAL_ERROR', 'SEND_ERROR', 'WARNING', and 'IGNORE'.  The
  default value is 'FATAL_ERROR' for a project when
  ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE=ON``, which is best for
  development mode for a project that currently has no invalid usage patterns.
  The default is 'IGNORE' when
  ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE=OFF``.  But a project with some
  existing invalid usage patterns might want to set, for example, a default of
  'WARNING' in order to allow for a smooth upgrade of TriBITS.  To do so,
  set::

    set(${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE_DEFAULT WARNING)

  in the project's base `<projectDir>/ProjectName.cmake`_ file.


.. _${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES:

**${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES**

  To set ``${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES`` a different default,
  set::

    set(${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_DEFAULT  <newDefault>)

  in the project's base `<projectDir>/ProjectName.cmake`_ file, where
  ``<newDefault>`` can be ``FATAL_ERROR``, ``SEND_ERROR``, ``WARNING``,
  ``NOTICE`` or ``IGNORE``

  Otherwise, the default is ``WARNING`` when
  ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE`` is ``ON`` and if ``IGNORE`` if
  ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE`` is ``OFF``.


.. _${PROJECT_NAME}_C_Standard:

**${PROJECT_NAME}_C_Standard**

  The variable ``${PROJECT_NAME}_C_Standard`` is used define the C standard
  pass to the compiler in ``--std=<cstd>`` for GCC builds of the project.
  TriBITS sets the default as ``c99`` but the project can set a new default in
  the project's base `<projectDir>/CMakeLists.txt`_ file with, for example::

    set(${PROJECT_NAME}_C_Standard_DEFAULT c11)

.. _${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS:

**${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS**

  The variable ``${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS`` determines how
  unparsed and otherwise ignored arguments are handled in TriBITS functions
  that are called by the client TriBITS projects.  These are arguments that
  are left over from parsing input options to functions and macros that take
  both positional arguments and keyword arguments/options handled with the
  ``cmake_parse_arguments()`` function.  For example, for the a TriBITS
  function declared like::

    tribits_copy_files_to_binary_dir(
      <targetName>
      [SOURCE_FILES <file1> <file2> ...]
      [SOURCE_DIR <sourceDir>]
      ...
      )

  the arguments ``SOURCE_FILES <file1> <file2> ...`` and those that follow are
  parsed by the ``cmake_parse_arguments()`` function while the argument
  ``<targetName>`` is a positional argument.  The problem is that any
  arguments passed between the first ``<targetName>`` argument and the
  specified keyword arguments like ``SOURCE_FILES`` and ``SOURCE_DIR`` are
  returned as unparsed arguments and are basically ignored (which is what
  happened in earlier versions of TriBITS).  For example, calling the function
  as::

    tribits_copy_files_to_binary_dir( FooTestCopyFiles
      ThisArgumentIsNotParsedAndIsIgnored
      SOURCE_FILES file1.cpp file2.cpp ...
      ...
      )

  would result in the unparsed argument
  ``ThisArgumentIsNotParsedAndIsIgnored``.

  The value of ``${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS`` determines how
  that ignored argument is handled.  If the value is ``WARNING``, then it will
  just result in a ``message(WARNING ...)`` command that states the warning
  but configure is allowed to be completed. This would be the right value to
  allow an old TriBITS project to keep configuring until the warnings can be
  cleaned up.  If the value is ``SEND_ERROR``, then ``message(SEND_ERROR
  ...)`` is called.  This will result in the configure failing but will allow
  configure to continue until the end (or a ``FATAL_ERROR`` is raised).  This
  would be the right value when trying to upgrade a TriBITS project where you
  wanted to see all of the warnings when upgrading TriBITS (so you could fix
  them all in one shot).  Finally, the value of ``FATAL_ERROR`` will result in
  ``message(FATAL_ERROR ...)`` being called which will halt configure right
  away.  This is the best value when developing on a TriBITS project that is
  already clean but you want to catch new developer-inserted errors right
  away.

  The default value for ``${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS`` is
  ``WARNING``, so that it will be backward compatible for TriBITS projects
  that might have previously undetected unparased and therefore ignored
  argument .  However, a project can change the default by setting, for
  example::

    set(${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS_DEFAULT FATAL_ERROR)

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

    set(${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_APPEND_DEFAULT
      "${CMAKE_CURRENT_LIST_DIR}/cmake/StdDevEnvs.cmake")

  **before** the `tribits_project()`_ command.  By including this file(s)
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

    set(${PROJECT_NAME}_CPACK_SOURCE_GENERATOR_DEFAULT "TGZ;TBZ2")

  This variable should generally be set in the file::

     <projectDir>/cmake/CallbackDefineProjectPackaging.cmake

  instead of in the base-level ``CMakeLists.txt`` file so that it goes along
  with rest of the project-specific CPack packaging options.

.. _${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE:

**${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE**

  The variable ``${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE`` determines if the
  CTest driver scripts using `tribits_ctest_driver()`_ configure, build, test
  and submit results to CDash all-at-once for all of the packages being tested
  or if instead is done package-by-package.  Currently, the default is set to
  ``FALSE`` for the package-by-package mode (for historical reasons) but the
  default can be set to ``TRUE`` by setting:

    set(${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE_DEFAULT "TRUE")

  in the project's `<projectDir>/ProjectName.cmake`_ file.  (This default must
  be changed in the ``<projectDir>/ProjectName.cmake`` file and **NOT** the
  `<projectDir>/CMakeLists.txt`_ file because the latter is not directly
  processed in CTest -S driver scripts using ``tribits_ctest_driver()``.)

  In general, a project should change the default to ``TRUE`` when using a
  newer CDash installation with CDash versions 3.0+ that can accommodate the
  results coming from ctest -S and display them package-by-package very
  nicely.  Otherwise, most projects are better off with package-by-package
  mode since it results in nicer display on CDash for older CDash versions.

.. _${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES:

**${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES**

  If `${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON`_ (the TriBITS
  default value), then any explicitly enabled packages that have disabled
  `upstream`_ required packages or TPLs will be disabled.  If
  `${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=OFF`_, then an
  configure error will occur.  For more details also see
  `TribitsBuildReference`_ and `Disables trump enables where there is a
  conflict`_.  A project can define a different default value by setting::
  
    set(${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT FALSE)

.. _${PROJECT_NAME}_ELEVATE_ST_TO_PT:

**${PROJECT_NAME}_ELEVATE_ST_TO_PT**

  If ``${PROJECT_NAME}_ELEVATE_ST_TO_PT`` is set to ``ON``, then all ``ST``
  packages will be elevated to ``PT`` packages.  The TriBITS default is
  obviously ``OFF``.  The default can be changed by setting::

    set(${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT ON)

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

    set(${PROJECT_NAME}_ENABLE_CPACK_PACKAGING_DEFAULT ON)

.. _${PROJECT_NAME}_ENABLE_CXX:
  
**${PROJECT_NAME}_ENABLE_CXX**
  
  If ``${PROJECT_NAME}_ENABLE_CXX`` is ``ON``, then C++ language support for
  the project will be enabled and the C++ compiler must be found.  By default,
  TriBITS sets this to ``ON`` for all systems.  A project never requires C++
  can set this to off by default by setting:
  
    set(${PROJECT_NAME}_ENABLE_CXX_DEFAULT FALSE)

.. _${PROJECT_NAME}_ENABLE_C:
  
**${PROJECT_NAME}_ENABLE_C**
  
  If ``${PROJECT_NAME}_ENABLE_C`` is ``ON``, then C language support for the
  project will be enabled and the C compiler must be found.  By default,
  TriBITS sets this to ``ON`` for all systems.  A project never requires C can
  set this to off by default by setting:
  
    set(${PROJECT_NAME}_ENABLE_C_DEFAULT FALSE)
  
  If a project does not have any native C code a good default would be::
  
    set(${PROJECT_NAME}_ENABLE_C_DEFAULT FALSE)
  
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
  development mode.  For example,
  ``${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES`` is set to ``FATAL_ERROR`` by
  default in development mode but is set to ``IGNORE`` by default in release
  mode.  In addition, strong compiler warnings are enabled by default in
  development mode but are disabled by default in release mode.  This variable
  also affects the behavior of `tribits_set_st_for_dev_mode()`_.
 
.. _${PROJECT_NAME}_ENABLE_Fortran:
  
**${PROJECT_NAME}_ENABLE_Fortran**
  
  If ``${PROJECT_NAME}_ENABLE_Fortran`` is ``ON``, then Fortran support for
  the project will be enabled and the Fortran compiler(s) must be found.  By
  default, TriBITS sets this to ``ON`` .
  
  If a project does not have any native Fortran code a good default would be::
  
    set(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT OFF)

  This default can be set in `<projectDir>/ProjectName.cmake`_ or `<projectDir>/CMakeLists.txt`_.

  On WIN32 systems, the default for ``${PROJECT_NAME}_ENABLE_Fortran_DEFAULT``
  is set to ``OFF`` since it can be difficult to get a Fortran compiler for
  native Windows.

  Given that a native Fortran compiler is not supported by default on Windows
  and on most Mac OSX systems, projects that have optional Fortran code may
  decide to set the default depending on the platform by setting, for example::

    if ( (WIN32 AND NOT CYGWIN) OR (CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin") )
      message(STATUS "Warning: Setting ${PROJECT_NAME}_ENABLE_Fortran=OFF by default"
       " because this is Windows (not cygwin) and we assume to not have Fortran!")
      set(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT OFF)
    else()
      set(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT ON)
    endif()
  
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
  include directories and libraries.  The TriBITS default is ``OFF`` but a
  project can change the default by setting, for example::

    set(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES_DEFAULT ON)

  A project would want to leave off the creation and installation of
  ``<PackageName>Config.cmake`` files if it was only installing and providing
  executables (see `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_).
  However, if it is wanting to provide libraries for other projects to use,
  then it should turn on the default generation of these files.

.. _${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE:

**${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE**

  If ``${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE`` is ``ON``, then packages
  and subpackages marked as ``ST`` in the `<repoDir>/PackagesList.cmake`_ file
  will be implicitly enabled along with the ``PT`` packages.  Additional code
  and tests may also be enabled using this option.  The TriBITS default is
  ``OFF`` but this can be changed by setting::

    set(${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE_DEFAULT ON)

  in the `<projectDir>/ProjectName.cmake`_ file.

.. _${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION:

**${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION**

  If ``${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION`` is
  ``TRUE``, then the directories for subpackages that are not enabled are left
  out of the source tarball.  This reduces the size of the tarball as much as
  possible but does require that the TriBITS packages and subpackages be
  properly set up to allow disabled subpackages from being excluded.  The
  TriBITS default is ``TRUE`` but this can be changed by setting::

    set(${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION_DEFAULT FALSE)

.. _${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES:

**${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES**

  If ``${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES`` is ``ON``, then the
  data-structures needed to generate ``<PackageName>Config.cmake`` files are
  created.  These data structures are also needed in order to generate export
  makefiles on demand using the function
  `tribits_write_flexible_package_client_export_files()`_.  The default in
  TriBITS is to turn this ``ON`` automatically by default if
  ``${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES`` is ``ON``.  Else, by
  default, TriBITS sets this to ``OFF``.  The only reason for the project to
  override the default is to set it to ``ON`` as with::

    set(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT ON)

  is so that the necessary data-structures are generated in order to use the
  function `tribits_write_flexible_package_client_export_files()`_.

.. _${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES:

**${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES**

  If ``${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES`` is ``ON``, then the files
  ``VersionDate.cmake`` and ``<RepoName>_version_date.h`` will get generated
  and the generated file ``<RepoName>_version_date.h`` will get installed for
  each TriBITS version-controlled repository when the local directories are
  git repositories.  The default is ``OFF`` but the project can change that by
  setting::

    set(${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES ON)

  in the `<projectDir>/ProjectName.cmake`_ file.

.. _${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE:

**${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE**

  If ``${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE`` is ``ON``, then the file
  ``<Project>RepoVersion.txt`` will get generated as a byproduct of
  configuring with CMake.  See `Multi-Repository Support`_ and
  `<Project>_GENERATE_REPO_VERSION_FILE`_.  The default is ``OFF`` but the
  project can change that by setting::

    set(${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT ON)

  in the `<projectDir>/ProjectName.cmake`_ file.

  Note that if a ``git`` exectauble cannot be found at configure time, then
  the default ``${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT`` will be
  overridden to ``OFF``.  But if the user sets
  ``${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE=ON`` in the cache and ``git``
  can't be found, then an configure-time error will occur.


.. _${PROJECT_NAME}_IMPORTED_NO_SYSTEM:

**${PROJECT_NAME}_IMPORTED_NO_SYSTEM**

  By default, include directories from IMPORTED library targets from the
  TriBITS project's installed ``<Package>Config.cmake`` files will be
  considered ``SYSTEM`` headers and therefore be included on the compile lines
  of downstream CMake projects with ``-isystem`` with most compilers.
  However, if ``${PROJECT_NAME}_IMPORTED_NO_SYSTEM`` is set to ``ON`` (only
  supported for CMake versions 3.23 or greater), then all of the IMPORTED
  library targets exported into the set of installed ``<Package>Config.cmake``
  files will have the ``IMPORTED_NO_SYSTEM`` property set.  This will cause
  downstream customer CMake projects to apply the include directories from
  these IMPORTED library targets as non-system include directories.  On most
  compilers, that means that the include directories will be listed on the
  compile lines with ``-I`` instead of with ``-isystem``.  (See more details
  in the TriBITS Build Reference for ``<Project>_IMPORTED_NO_SYSTEM``.)

  The default value set by TriBITS itself is ``OFF`` but a TriBITS project can
  change the default value to ``ON`` by adding::

    if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.23)
      set(${PROJECT_NAME}_IMPORTED_NO_SYSTEM_DEFAULT ON)
    endif()

  in the `<projectDir>/ProjectName.cmake`_ file.  (NOTE: The above ``if()``
  statement ensures that a configure error will not occur if a version of
  CMake less than 3.23 is used.  But if the TriBITS project minimum CMake
  version is 3.23 or greater, then the above ``if()`` statement guard can be
  removed.)


.. _${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS:

**${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS**

  If ``${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`` is set to ``ON``, then
  any defined libraries or header files that are listed in calls to
  `tribits_add_library()`_ or `tribits_install_headers()`_ will be installed
  (unless options are passed into `tribits_add_library()`_ that disable
  installs).  If set to ``OFF``, then headers and libraries will *not* be
  installed by default and only ``INSTALLABLE`` executables added with
  `tribits_add_executable()`_ will be installed.  However, as described in
  `TribitsBuildReference`_, shared libraries will always be installed if
  enabled since they are needed by the installed executables.
  
  For a TriBITS project that is primarily delivering libraries
  (e.g. Trilinos), then it makes sense to leave the TriBITS default which is
  ``ON`` or explicitly set::
  
    set(${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT  ON)
  
  For a TriBITS project that is primarily delivering executables (e.g. VERA),
  then it makes sense to set the default to::
  
    set(${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT  OFF)

.. _${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE:

.. _${PROJECT_NAME}_MAKE_INSTALL_GROUP_WRITABLE:

.. _${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE:

**${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE**
**${PROJECT_NAME}_MAKE_INSTALL_GROUP_WRITABLE**
**${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE**

  Determines the permissions for directories and files created during the
  execution of the of the ``install`` and ``isntall_package_by_package``
  targets.

  To make the created directories by only group readable for the project by
  default, set::

    set(${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE_DEFAULT TRUE)

  To make the created directories by only group writable (and readable) for
  the project by default, set::

    set(${PROJECT_NAME}_MAKE_INSTALL_WORLD_WRITABLE_DEFAULT TRUE)

  To make the created directories by world readable for the project by
  default, set::

    set(${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE_DEFAULT TRUE)

  On non-Windows systems, these set permissions for all files and directories
  from the the user-set base directory
  ``${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR`` on down.
  For more details see `Installation considerations`_.

  These defaults can be set in the `<projectDir>/ProjectName.cmake`_ file.

.. _${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS:

**${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS**

  Determines if all of the libraries listed in ``<tplName>_LIBRARY_NAMES`` for
  a given TPL must be found for each enabled TPL.  By default, this is
  ``FALSE`` which means that the determination if all of the listed libs for a
  TPL should be found is determined by the ``MUST_FIND_ALL_LIBS`` option to
  the `tribits_tpl_find_include_dirs_and_libraries()`_ function in the TPL
  find module.  To change the default for this, set::

    set(${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS_DEFAULT  TRUE)

  in the `<projectDir>/ProjectName.cmake`_ file.


.. _${PROJECT_NAME}_REQUIRES_PYTHON:

**${PROJECT_NAME}_REQUIRES_PYTHON**

  If the TriBITS project requires Python, set::

    set(${PROJECT_NAME}_REQUIRES_PYTHON  TRUE)

  in the `<projectDir>/ProjectName.cmake`_ file (See `Python Support`_).  The
  default is implicitly ``FALSE``.


.. _${PROJECT_NAME}_SET_INSTALL_RPATH:

**${PROJECT_NAME}_SET_INSTALL_RPATH**

  The cache variable ``${PROJECT_NAME}_SET_INSTALL_RPATH`` is used to define
  the default RPATH mode for the TriBITS project (see `Setting install RPATH`_
  for details).  The TriBITS default is to set this to ``TRUE`` but the
  TriBITS project can be set the default to ``FALSE`` by setting::

    set(${PROJECT_NAME}_SET_INSTALL_RPATH_DEFAULT FALSE)

  in the project's `<projectDir>/ProjectName.cmake`_ file (see `RPATH
  Handling`_).


.. _${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME:

**${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME**

  The cache variable ``${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME``
  determines if the start and end date/time for each advanced test (i.e. added
  with `tribits_add_advanced_test()`_) is printed or not with each test.  If
  set to ``TRUE`` this also causes in the timing for each ``TEST_<IDX>`` block
  to be printed as well.  The TriBITS default is ``OFF`` but a TriBITS project
  can change this default by setting::

    set(${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME_DEFAULT ON)

  The implementation of this feature currently uses ``execute_process(date)``
  and therefore will work on many (but perhaps not all) Linux/Unix/Mac systems
  and not on Windows systems.

  NOTE: In a future version of CTest, this option may turn on start and end
  date/time for regular tests added with `tribits_add_test()`_ (which uses a
  raw command with ``add_test()``).

.. _${PROJECT_NAME}_SKIP_INSTALL_PROJECT_CMAKE_CONFIG_FILES:

**${PROJECT_NAME}_SKIP_INSTALL_PROJECT_CMAKE_CONFIG_FILES**

  To change the default value of the
  ``${PROJECT_NAME}_SKIP_INSTALL_PROJECT_CMAKE_CONFIG_FILES`` to ``TRUE``, for
  example, for a TriBITS project, set::

    set(${PROJECT_NAME}_SKIP_INSTALL_PROJECT_CMAKE_CONFIG_FILES_DEFAULT TRUE)

  in the project's `<projectDir>/CMakeLists.txt`_ or
  `<projectDir>/ProjectName.cmake`_ files.

.. _${PROJECT_NAME}_SKIP_EXTRAREPOS_FILE:

**${PROJECT_NAME}_SKIP_EXTRAREPOS_FILE**

  The cache variable ``${PROJECT_NAME}_SKIP_EXTRAREPOS_FILE`` is set in the
  `<projectDir>/ProjectName.cmake`_ file as::

    set(${PROJECT_NAME}_SKIP_EXTRAREPOS_FILE TRUE)

  for projects that don't have a
  `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file.  This variable needs
  to be set when using the CTest driver script and does not need to be set for
  the basic configure and build process.

.. _${PROJECT_NAME}_TEST_CATEGORIES:
.. _${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT:

**${PROJECT_NAME}_TEST_CATEGORIES**

  The cache variable ``${PROJECT_NAME}_TEST_CATEGORIES`` determines what tests
  defined using `tribits_add_test()`_ and `tribits_add_advanced_test()`_ will
  be added for ``ctest`` to run (see `Test Test Category`_).  The TriBITS
  default is ``NIGHTLY`` for a standard local build.  The `checkin-test.py`_
  tool sets this to ``BASIC`` by default.  A TriBITS project can override the
  default for a basic configure using, for example::

    set(${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT BASIC)

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
  ``SYSTEM`` flag will be passed into the ``include_directories()`` command
  for TPL include directories for every TPL for every package, by default.  On
  some systems this will result in include directories being passed to the
  compiler with ``-isystem`` instead of ``-I``.  This helps to avoid compiler
  warning coming from TPL header files for C and C++.  However, with CMake
  version 3.2 and less, this also results in ``-isystem`` being passed to the
  Fortran compiler (e.g. gfortran) as well.  This breaks the reading of
  Fortran module files (perhaps a bug in gfortran).  Because of this issue
  with Fortran, the TriBITS default for this option is set to ``FALSE`` but a
  project can override the default using::

    set(${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS_DEFAULT  TRUE)

  (This would be a good default if the project has not Fortran files or has
  not Fortran files that use modules provided by TPLs).

  However, if a package or subpackage sets::

    set(${PACKAGE_NAME}_SKIP_TPL_SYSTEM_INCLUDE_DIRS  TRUE)

  in its ``CMakeLists.txt`` files before the ``tribits_add_library()`` or
  ``tribits_add_executable()`` commands are called in that package, then
  ``SYSTEM`` will **not** be passed into ``include_directories()`` for TPL
  include dirs.  This is how some TriBITS packages with Fortran files that use
  Fortran modules avoid passing in ``-isystem`` to the Fortran compiles and
  thereby avoid the defect with gfortran described above.  If CMake version
  3.3 or greater is used, this variable is not required.

  NOTE: Currently, a TriBITS package must have a direct dependency on a TPL
  to have ``-isystem`` added to a TPL's include directories on the compile
  lines for that package.  That is, the TPL must be listed in the
  ``LIB_REQUIRED_TPLS`` or ``LIB_OPTIONAL_TPLS`` arguments passed into the
  `tribits_package_define_dependencies()`_ function in the package's
  `<packageDir>/cmake/Dependencies.cmake`_ file.  In addition, to have
  ``-isystem`` added to the include directories for a TPL when compiling the
  tests for an package, it must be listed in the ``TEST_REQUIRED_TPLS`` or
  ``TEST_OPTIONAL_TPLS`` arguments.  This is a limitation of the TriBITS
  implementation that will be removed in a future version of TriBITS.

.. _${PROJECT_NAME}_TRACE_ADD_TEST:
.. _${PROJECT_NAME}_TRACE_ADD_TEST_DEFAULT:

**${PROJECT_NAME}_TRACE_ADD_TEST**

  If ``${PROJECT_NAME}_TRACE_ADD_TEST`` is set to ``TRUE``, then a single line
  will be printed for each call to `tribits_add_test()`_ and
  `tribits_add_advanced_test()`_ for if the test is added or not and if not
  then why.  The default is set based on the value of
  ``${PROJECT_NAME}_VERBOSE_CONFIGURE`` but a project can override the default
  by setting::

    set(${PROJECT_NAME}_TRACE_ADD_TEST_DEFAULT  TRUE)

.. _${PROJECT_NAME}_USE_GNUINSTALLDIRS:

**${PROJECT_NAME}_USE_GNUINSTALLDIRS**

  If ``${PROJECT_NAME}_USE_GNUINSTALLDIRS`` is set to ``TRUE``, then the
  default install paths will be determined by the standard CMake module
  ``GNUInstallDirs``.  Otherwise, platform independent install paths are used
  by default.

  A project can use the paths given the cmake module ``GNUInstallDirs`` by
  default by setting::
  
    set(${PROJECT_NAME}_USE_GNUINSTALLDIRS_DEFAULT TRUE)

  in the project's top-level `<projectDir>/CMakeLists.txt`_ file or its
  `<projectDir>/ProjectName.cmake`_ file.  The default is ``FALSE``.

.. _${PROJECT_NAME}_USES_PYTHON:

**${PROJECT_NAME}_USES_PYTHON**

  If the TriBITS project can use Python, but does not require it, set::

    set(${PROJECT_NAME}_USES_PYTHON  TRUE)

  in the `<projectDir>/ProjectName.cmake`_ file (see `Python Support`_).  The
  default for a TriBITS project is implicitly ``TRUE``.  To explicitly state
  that Python is never needed, set::

    set(${PROJECT_NAME}_USES_PYTHON  FALSE)

.. _DART_TESTING_TIMEOUT:

**DART_TESTING_TIMEOUT**

  The cache variable ``DART_TESTING_TIMEOUT`` is a built-in CMake variable
  that provides a default timeout for all tests (see `Setting test timeouts at
  configure time`_).  By default, TriBITS defines this to be 1500 seconds
  (which is also the raw CMake default) but the project can change this
  default, from 1500 to 300 for example, by setting the following in the
  project's `<projectDir>/ProjectName.cmake`_ or
  `<projectDir>/CMakeLists.txt`_ file::

    set(DART_TESTING_TIMEOUT_DEFAULT 300)

.. _CMAKE_INSTALL_RPATH_USE_LINK_PATH:

**CMAKE_INSTALL_RPATH_USE_LINK_PATH**

  The cache variable ``CMAKE_INSTALL_RPATH_USE_LINK_PATH`` is a built-in CMake
  variable that determines if the paths for external libraries (i.e. from
  TPLs) is put into the installed library RPATHS (see `RPATH Handling`_).
  TriBITS sets the default for this to ``TRUE`` but a project can change the
  default back to ``FALSE`` by setting the following in the project's
  `<projectDir>/ProjectName.cmake`_ file::

    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH_DEFAULT FALSE)

.. _MPI_EXEC_MAX_NUMPROCS:

**MPI_EXEC_MAX_NUMPROCS**

  The variable ``MPI_EXEC_MAX_NUMPROCS`` gives the maximum number of processes
  for an MPI test that will be allowed as defined by `tribits_add_test()`_ and
  `tribits_add_advanced_test()`_.  The TriBITS default is set to be ``4`` (for
  no good reason really but it needs to stay that way for backward
  compatibility).  This default can be changed by setting::

    set(MPI_EXEC_MAX_NUMPROCS_DEFAULT <newDefaultMax>)

  While this default can be changed for the project as a whole on all
  platforms, it is likely better to change this default on a
  machine-by-machine basis to correspond to the load that can be accommodated
  by a given machine (or class of machines).  For example if a given machine
  has 64 cores, a reasonable number for ``MPI_EXEC_MAX_NUMPROCS_DEFAULT`` is
  64.

.. _PythonInterp_FIND_VERSION:

**PythonInterp_FIND_VERSION**

  Determines the version of Python that is looked for.  TriBITS requires at
  least version "2.7".  A particular TriBITS project can require a higher
  version of TriBITS and this is set using, for example:

    set(PythonInterp_FIND_VERSION_DEFAULT "3.5.2")

  in the `<projectDir>/ProjectName.cmake`_ file (See `Python Support`_).  The
  default is version "2.7".  The user can force a more recent version of
  Python by configuring with, for example::

    -D PythonInterp_FIND_VERSION="3.6.2"

.. _TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE:

**TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE**

  Determines how the function `tribits_deprecated()`_ behaves.  To change the
  default behavor, such as call ``message(FATAL_ERROR ...)``, set::

    set(TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE_DEFAULT  FATAL_ERROR)

  in the project's `<projectDir>/ProjectName.cmake`_ file, or
  `<projectDir>/CMakeLists.txt`_ file, or on the individual package basis in
  its `<packageDir>/CMakeLists.txt`_ file.


TriBITS Macros and Functions
----------------------------

The following subsections give detailed documentation for the CMake macros and
functions that make up the core TriBITS system.  These are what are used by
TriBITS project developers in their ``CMakeLists.txt`` and other files.  All
of these functions and macros should be automatically available when
processing the project's and package's variables files if used properly.
Therefore, no explicit ``include()`` statements should be needed other than
the initial include of the ``TriBITS.cmake`` file in the top-level
`<projectDir>/CMakeLists.txt`_ file so the command `tribits_project()`_ can be
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


.. *** NOTE: This file is symlinked into users_guide/ and maintainers_guide/
.. *** and the symlinked file gets included in the top-level *.rst document.
.. *** This is so that the includes for the generated files
.. *** TribitsMacroFunctionDoc.rst and UtilsMacroFunctionDoc.rst in those
.. *** directories work correctly.
