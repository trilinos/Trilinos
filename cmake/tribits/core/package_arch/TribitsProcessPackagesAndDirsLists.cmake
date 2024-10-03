# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include(SetCacheOnOffEmpty)
include(MultilineSet)
include(AdvancedOption)
include(AdvancedSet)
include(AssertDefined)
include(PrintVar)
include(MessageWrapper)
include(TribitsListHelpers)


#
# Helper functions
#


# @MACRO: tribits_repository_define_packages()
#
# Define the set of packages for a given `TriBITS Repository`_.  This macro is
# typically called from inside of a `<repoDir>/PackagesList.cmake`_ file for a
# given TriBITS repo.
#
# Usage::
#
#    tribits_repository_define_packages(
#       <pkg0>  <pkg0_dir>  <pkg0_classif>
#       <pkg1>  <pkg1_dir>  <pkg1_classif>
#       ...
#       )
#
# This macro sets up a 2D array of ``NumPackages`` by ``NumColumns`` listing
# out the packages for a TriBITS repository.  Each row (with 3 column entries)
# specifies a package which contains the columns (ordered 0-2):
#
# 0. **PACKAGE** (``<pkgi>``): The name of the TriBITS package.  This name
#    must be unique across all other TriBITS packages in this or any other
#    TriBITS repo that might be combined into a single TriBITS project
#    meta-build (see `Globally unique TriBITS package names`_).  The name
#    should be a valid identifier (e.g. matches the regex
#    ``[a-zA-Z_][a-zA-Z0-9_]*``).  The package names tend to use mixed case
#    (e.g. ```SomePackge`` not ``SOMEPACKGE``).
#
# 1. **DIR** (``<pkgi_dir>``): The relative directory for the package
#    ``<packageDir>``.  This directory is relative to the TriBITS repository
#    base directory ``<repoDir>``.  Under this directory will be a
#    package-specific ``cmake/`` directory with the file
#    `<packageDir>/cmake/Dependencies.cmake`_ and a base-level
#    `<packageDir>/CMakeLists.txt`_ file.  The entire contents of the package
#    including all of the source code and all of the tests should be contained
#    under this directory.  The TriBITS testing infrastructure relies on the
#    mapping of changed files to these base directories when deciding what
#    packages are modified and need to be retested (along with downstream
#    packages).  For details, see `checkin-test.py`_.
#
# 2. **CLASSIFICATION** (``<pkgi_classif>``): Gives the `Package Test
#    Group`_ `PT`_, `ST`_, or `EX`_ and the maturity level ``EP``, ``RS``,
#    ``PG``, ``PM``, ``GRS``, ``GPG``, ``GPM``, ``UM``.  These are separated
#    by a coma with no space in between such as ``"RS,PT"`` for a "Research
#    Stable", "Primary Tested" package.  No spaces are allowed so that CMake
#    treats this a one field in the array.  The maturity level can be left off
#    in which case it is assumed to be ``UM`` for "Unspecified Maturity".
#    This classification for individual packages can be changed to ``EX`` for
#    specific platforms by calling `tribits_disable_package_on_platforms()`_.
#
# **IMPORTANT:** The packages must be listed in increasing order of package
# dependencies.  That is `No circular dependencies of any kind are allowed`_
# (see the *ADP (Acyclic Dependencies Principle)* in `Software Engineering
# Packaging Principles`_).  Package ``i`` can only list dependencies (in
# `<packageDir>/cmake/Dependencies.cmake`_) for packages listed before this
# package in this list (or in upstream TriBITS repositories).  This avoids an
# expensive package sorting algorithm and makes it easy to flag packages with
# circular dependencies or misspelling of package names.
#
# NOTE: For some rare use cases, the package directory ``<pkgi_dir>`` is
# allowed to be specified as an absolute directory but this absolute directory
# must be a subdirectory of the project source base directory given by
# `PROJECT_SOURCE_DIR`_.  If not, ``message(FATAL_ERROR ...)`` is called and
# processing stops immediately.
#
# NOTE: This macro just sets the variable::
#
#   ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
#
# in the current scope.  The advantages of using this macro instead of
# directly setting this variable are that the macro:
#
# * Asserts that the variable ``REPOSITORY_NAME`` is defined and set
#
# * Avoids having to hard-code the assumed repository name
#   ``${REPOSITORY_NAME}``.  This provides more flexibility for how other
#   TriBITS projects choose to name a given TriBITS repo (i.e. the name of
#   repo subdirs).
#
# * Avoid misspelling the name of the variable
#   ``${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS``.  If one
#   misspells the name of the macro, it is an immediate error in CMake.
#
macro(tribits_repository_define_packages)
  assert_defined(REPOSITORY_NAME)
  set(${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS "${ARGN}")
endmacro()


# @FUNCTION: tribits_allow_missing_external_packages()
#
# Allow listed packages to be missing and automatically excluded from the
# package dependency data-structures.
#
# Usage::
#
#   tribits_allow_missing_external_packages(<pkg0> <plg1> ...)
#
# If the missing upstream package ``<pkgi>`` is optional, then the effect
# will be to simply ignore the missing package (i.e. it will never be added to
# package's list and not added to dependency data-structures) and remove it
# from the dependency lists for downstream packages that have an optional
# dependency on the missing upstream package ``<pkgi>``.  However, all
# downstream packages that have a required dependency on the missing
# upstream package ``<pkgi>`` will be hard disabled,
# i.e. ``${PROJECT_NAME}_ENABLE_{CURRENT_PACKAGE}=OFF`` and a note on the
# disable will be printed.
# 
# **WARNING**: This macro just sets the cache variable
# ``<pkgi>_ALLOW_MISSING_EXTERNAL_PACKAGE=TRUE`` for each package
# ``<pkgi>``.  Therefore, using this function effectively turns off error
# checking for misspelled package names so it is important to only use it when
# it absolutely is needed (use cases mentioned below).  Also note that missing
# packages are silently ignored by default.  Therefore, when doing development
# involving these packages, it is usually a good idea to set::
#
#   -D<pkgi>_ALLOW_MISSING_EXTERNAL_PACKAGE=FALSE
#
# so that it will catch errors in the misspelling of package names or source
# directories.  However, notes on what missing packages are being ignored can
# printed by configuring with::
#
#   -D <Project>_WARN_ABOUT_MISSING_EXTERNAL_PACKAGES=TRUE
#
# This macro is typically called in one of two different contexts:
#
# * Called from `<packageDir>/cmake/Dependencies.cmake`_ when the upstream
#   package ``<pkgi>`` is defined in an optional upstream `TriBITS
#   Repository`_.  This allows the downstream repos and packages to still be
#   enabled (assuming they don't have required dependencies on the missing
#   packages) when one or more upstream repos are missing.
#
# * Called from `<repoDir>/PackagesList.cmake`_ when the package ``<pkgi>``
#   might be defined in an optional non-TriBITS repo (see `How to insert a
#   package into an upstream repo`_).
#
# For some meta-projects that composes packages from may different TriBITS
# repositories, one might need to also call this function from the file
# `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_.
#
function(tribits_allow_missing_external_packages)
  foreach(TRIBITS_PACKAGE ${ARGN})
    advanced_set(${TRIBITS_PACKAGE}_ALLOW_MISSING_EXTERNAL_PACKAGE TRUE
      CACHE BOOL
      "Default set by tribits_allow_missing_external_packages()!"
      )
  endforeach()
endfunction()


#
# Below, we change the value of user cache values like
# ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME},
# ${PACKAGE_NAME}_ENABLE_TESTS, and ${PACKAGE_NAME}_ENABLE_EXAMPLES by
# just setting them to regular variables that live in the top scope
# instead of putting them back in the cache.  That means that they are
# used as global variables but we don't want to disturb the cache
# since that would change the behavior for future invocations of cmake
# (which is very confusing).  Because of this, these macros must all
# be called from the top-level ${PROJECT_NAME} CMakeLists.txt file and
# macros must call macros as not to change the variable scope.
#
# I had to do it this way in order to be able to get the right behavior which
# is:
#
# 1) Override the value of these variables in all CMake processing
#
# 2) Avoid changing the user cache values because that would be confusing and
# would make it hard to change package enables/disable later without blowing
# away the cache
#


# Macro that sets up standard user options a package
#
# On completion, the following variables are set:
#
# * ${PACKAGE_NAME_IN}_TESTGROUP: Set to PT, ST, or EX
#
macro(tribits_insert_standard_package_options  PACKAGE_NAME_IN  PACKAGE_TESTGROUP_IN)

  if (TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS_DEBUG)
    message("TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS: ${PACKAGE_NAME_IN} ${PACKAGE_TESTGROUP_IN}")
    print_var(${PACKAGE_NAME_IN}_TESTGROUP)
  endif()

  set(PACKAGE_TESTGROUP_LOCAL ${PACKAGE_TESTGROUP_IN})

  # ${PROJECT_NAME}_ELEVATE_ST_TO_PT is deprecated but allowed for backward compatibility
  if (${PROJECT_NAME}_ELEVATE_SS_TO_PS)
    set(${PROJECT_NAME}_ELEVATE_ST_TO_PT ON)
  endif()

  if (${PACKAGE_TESTGROUP_IN} STREQUAL PT OR ${PACKAGE_TESTGROUP_IN} STREQUAL ST)
    if (TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS_DEBUG)
      message("-- " "PT or ST")
      print_var(${PROJECT_NAME}_ELEVATE_ST_TO_PT)
    endif()
    if (${PROJECT_NAME}_ELEVATE_ST_TO_PT)
      set(PACKAGE_TESTGROUP_LOCAL PT)
    endif()
    if (TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS_DEBUG)
      print_var(PACKAGE_TESTGROUP_LOCAL)
    endif()
    set(PACKAGE_ENABLE "")
  elseif (${PACKAGE_TESTGROUP_IN} STREQUAL EX)
    if (TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS_DEBUG)
      message("-- " "EX")
    endif()
    set(PACKAGE_ENABLE OFF)
  else()
    message(FATAL_ERROR "Error the package classification '${PACKAGE_TESTGROUP_IN}'"
      " for the package ${PACKAGE_NAME_IN} is not a valid classification." )
  endif()

  if (TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS_DEBUG)
    print_var(PACKAGE_ENABLE)
    print_var(${PACKAGE_NAME_IN}_TESTGROUP)
  endif()

  if ("${${PACKAGE_NAME_IN}_TESTGROUP}" STREQUAL "") # Allow testing override
    if (TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS_DEBUG)
      message("-- " "Setting classification to ${PACKAGE_TESTGROUP_LOCAL}")
      print_var(PACKAGE_TESTGROUP_LOCAL)
    endif()
    set(${PACKAGE_NAME_IN}_TESTGROUP "${PACKAGE_TESTGROUP_LOCAL}")
  endif()

  if (TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS_DEBUG)
    print_var(${PACKAGE_NAME_IN}_TESTGROUP)
  endif()

  multiline_set(DOCSTR
    "Enable the package ${PACKAGE_NAME_IN}.  Set to 'ON', 'OFF', or leave"
    " empty to allow for other logic to decide."
    )
  set_cache_on_off_empty( ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME_IN}
    "${PACKAGE_ENABLE}" ${DOCSTR} )

endmacro()


# Function that determines if a package is a primary meta-project package
#  according to the variables
#  ${PARENT_REPO_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES[_EXCEPT].
#
function(tribits_is_primary_meta_project_package  PACKAGE_NAME_IN
  IS_PRIMARY_META_PROJECT_PACKAGE_OUT
  )

  set(IS_PRIMARY_META_PROJECT_PACKAGE TRUE)

  assert_defined(${PACKAGE_NAME_IN}_PARENT_REPOSITORY)
  set(PARENT_REPO_NAME ${${PACKAGE_NAME_IN}_PARENT_REPOSITORY})
  if (${PARENT_REPO_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES)
    find_list_element(
      ${PARENT_REPO_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES_EXCEPT
      ${PACKAGE_NAME_IN}  PACKAGE_EXCEPTED
      )
    if (PACKAGE_EXCEPTED)
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("-- "
          "NOTE: ${PACKAGE_NAME_IN} is classified as a primary meta-project packages even"
          " though ${PARENT_REPO_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES=ON "
          " because the package is included in the list ${PARENT_REPO_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES_EXCEPT!")
      endif()
    else()
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("-- "
          "NOTE: ${PACKAGE_NAME_IN} is not as a primary meta-project packages"
          " because ${PARENT_REPO_NAME}_NO_PRIMARY_META_PROJECT_PACKAGES=ON "
          )
      endif()
      set(IS_PRIMARY_META_PROJECT_PACKAGE FALSE)
    endif()
  endif()

  set(${IS_PRIMARY_META_PROJECT_PACKAGE_OUT} ${IS_PRIMARY_META_PROJECT_PACKAGE}
    PARENT_SCOPE )

endfunction()


# Function that determines if it is okay to allow an implicit enable of an
# upstream package given the disable of a downstream package that depends on
# it.
#
function(tribits_implicit_package_enable_is_allowed  upstreamPackageName
    packageName  implictPackageEnableAllowedOut
  )

  if (${packageName}_PACKAGE_BUILD_STATUS  STREQUAL  "EXTERNAL")
    set(implicitPackageEnableAllowed  FALSE)
  elseif (${packageName}_TESTGROUP  STREQUAL  "PT")
    set(implicitPackageEnableAllowed TRUE)
  elseif (${packageName}_TESTGROUP  STREQUAL  "ST"
      AND ${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE
    )
    set(implicitPackageEnableAllowed  TRUE)
  else()
    if (upstreamPackageName)
      message("-- " "NOTE: Not Setting ${PROJECT_NAME}_ENABLE_${packageName}=ON"
        " even though ${upstreamPackageName} has an optional dependence on"
        " ${packageName} because ${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE=OFF" )
    endif()
    set(implicitPackageEnableAllowed  FALSE)
  endif()

  set(${implictPackageEnableAllowedOut} "${implicitPackageEnableAllowed}"
    PARENT_SCOPE )

endfunction()


# @MACRO: tribits_process_packages_and_dirs_lists()
#
# Usage::
#
#   tribits_process_packages_and_dirs_lists()
#
# Macro that processes the list variable::
#
#    ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
#
# from a `<repoDir>/PackagesList.cmake`_ file that just got read in and
# creates/updates the top-level non-cache variables:
#
#   * `${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES`_
#   * `${PROJECT_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES`_
#   * ``${PROJECT_NAME}_LAST_DEFINED_INTERNAL_TOPLEVEL_PACKAGE_IDX``
#
# For each of the listed top-level (parent) packages ${PACKAGE_NAME}, it also
# sets up constant variables defined in `TriBITS Package Top-Level Local
# Variables`_ like:
#
# * `${PACKAGE_NAME}_SOURCE_DIR`_
# * `${PACKAGE_NAME}_REL_SOURCE_DIR`_
# * `${PACKAGE_NAME}_PARENT_PACKAGE`_ (to empty "")
# * `${PACKAGE_NAME}_PARENT_REPOSITORY`_ (to empty "")
# * `${PACKAGE_NAME}_TESTGROUP`_
# * `${PACKAGE_NAME}_PACKAGE_BUILD_STATUS`_ (to ``INTERNAL``)
#
# and sets up some standard enable/disable vars with default values as defined
# in `TriBITS Package Cache Variables`_ like:
#
# * `${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}`_
#
# NOTE: Set ``TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE=TRUE`` to see
# really verbose debug output from this macro.
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_process_packages_and_dirs_lists  REPOSITORY_NAME  REPOSITORY_DIR)

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS:  '${REPOSITORY_NAME}'  '${REPOSITORY_DIR}'")
  endif()

  #
  # Separate out separate lists of package names and directories
  #

  # Get the total number of packages defined

  assert_defined(${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS)
  if (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
    print_var(${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS)
  endif()
  list(LENGTH ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
    ${REPOSITORY_NAME}_NUM_PACKAGES_AND_FIELDS )
  math(EXPR ${REPOSITORY_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES
    "${${REPOSITORY_NAME}_NUM_PACKAGES_AND_FIELDS}/${PLH_NUM_FIELDS_PER_PACKAGE}")
  if (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
    print_var(${REPOSITORY_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES)
  endif()
  math(EXPR ${REPOSITORY_NAME}_LAST_DEFINED_INTERNAL_TOPLEVEL_PACKAGE_IDX
    "${${REPOSITORY_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES}-1")

  # Process each of the packages defined

  if (${REPOSITORY_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES GREATER 0)

    foreach(PACKAGE_IDX  RANGE
        ${${REPOSITORY_NAME}_LAST_DEFINED_INTERNAL_TOPLEVEL_PACKAGE_IDX}
      )

      if (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        message("")
        print_var(${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES)
      endif()

      math(EXPR PACKAGE_NAME_IDX "${PACKAGE_IDX}*${PLH_NUM_FIELDS_PER_PACKAGE}+0")
      list(GET ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
        ${PACKAGE_NAME_IDX} TRIBITS_PACKAGE )
      if (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        print_var(TRIBITS_PACKAGE)
      endif()

      math(EXPR PACKAGE_DIR_IDX
        "${PACKAGE_IDX}*${PLH_NUM_FIELDS_PER_PACKAGE}+${PLH_NUM_PACKAGE_DIR_OFFSET}")
      list(GET ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
        ${PACKAGE_DIR_IDX} PACKAGE_DIR )
      if (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        print_var(PACKAGE_DIR)
      endif()

      math(EXPR PACKAGE_CLASSIFICATION_IDX
        "${PACKAGE_IDX}*${PLH_NUM_FIELDS_PER_PACKAGE}+${PLH_NUM_PACKAGE_CLASSIFICATION_OFFSET}")
      list(GET ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
        ${PACKAGE_CLASSIFICATION_IDX} PACKAGE_CLASSIFICATION )
      if (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        print_var(PACKAGE_CLASSIFICATION)
      endif()

      # ToDo: Parse out TESTGROUP and MATURITYLEVEL (Trilinos #6042)
      set(PACKAGE_TESTGROUP ${PACKAGE_CLASSIFICATION})

      tribits_update_ps_pt_ss_st(Package ${TRIBITS_PACKAGE} PACKAGE_TESTGROUP)

      if (${TRIBITS_PACKAGE}_SOURCE_DIR_OVERRIDE)
        message("-- "
          "NOTE: ${TRIBITS_PACKAGE}_SOURCE_DIR_OVERRIDE='${${TRIBITS_PACKAGE}_SOURCE_DIR_OVERRIDE}'"
          " is overriding default path '${PACKAGE_DIR}'") 
        set(PACKAGE_DIR "${${TRIBITS_PACKAGE}_SOURCE_DIR_OVERRIDE}")
      endif()

      if (IS_ABSOLUTE "${PACKAGE_DIR}")

        set(PACKAGE_ABS_DIR "${PACKAGE_DIR}")

        string(LENGTH "${PROJECT_SOURCE_DIR}" PROJECT_SOURCE_DIR_LEN)
        string(LENGTH "${PACKAGE_ABS_DIR}" PACKAGE_ABS_DIR_LEN)

        # See if the package dir is under the project dir

        set(PACKAGE_ABS_DIR_UNDER_PROJECT_SOURCE_DIR TRUE)

        if (PACKAGE_ABS_DIR_UNDER_PROJECT_SOURCE_DIR)
          # Determine package abs dir is too short to be under project
          if (PACKAGE_ABS_DIR_LEN LESS PROJECT_SOURCE_DIR_LEN)
            set(PACKAGE_ABS_DIR_UNDER_PROJECT_SOURCE_DIR FALSE)
          endif()
        endif()

        if (PACKAGE_ABS_DIR_UNDER_PROJECT_SOURCE_DIR)
          # Determine if the package abs base dir base is the project dir
          string(SUBSTRING "${PACKAGE_ABS_DIR}" 0 ${PROJECT_SOURCE_DIR_LEN}
            PROJECT_SOURCE_DIR_BASE_MATCH)
          print_var(PROJECT_SOURCE_DIR_BASE_MATCH)
          if (NOT PROJECT_SOURCE_DIR_BASE_MATCH STREQUAL "${PROJECT_SOURCE_DIR}")
            set(PACKAGE_ABS_DIR_UNDER_PROJECT_SOURCE_DIR FALSE)
          endif()
        endif()

        if (PACKAGE_ABS_DIR_UNDER_PROJECT_SOURCE_DIR)
          # Get the path of the package dir under the project dir
          math(EXPR PACKAGE_REL_DIR_BEGIN "${PROJECT_SOURCE_DIR_LEN}+1")
          string(SUBSTRING "${PACKAGE_ABS_DIR}" ${PACKAGE_REL_DIR_BEGIN} -1
            REPOSITORY_AND_PACKAGE_DIR)
        else()
          message_wrapper(FATAL_ERROR
            "Error: The package '${TRIBITS_PACKAGE}' was given an absolute directory '${PACKAGE_ABS_DIR}' which is *not* under the project's source directory '${PROJECT_SOURCE_DIR}/'!")
          set(REPOSITORY_AND_PACKAGE_DIR "ERROR-BAD-PACKAGE-ABS-DIR")
          # ToDo: We could just put in a relative path but that requires
          # knowing the common path between the two directory paths but CMake
          # does not give an easy way to determine that.  I would have to
          # write that function myself.
        endif()

      else()

         # PACKAGE_DIR is a relative path

        if ("${REPOSITORY_DIR}" STREQUAL ".")
          set(REPOSITORY_AND_PACKAGE_DIR "${PACKAGE_DIR}")
        elseif("${PACKAGE_DIR}" STREQUAL ".")
          set(REPOSITORY_AND_PACKAGE_DIR "${REPOSITORY_DIR}")
        else()
          set(REPOSITORY_AND_PACKAGE_DIR "${REPOSITORY_DIR}/${PACKAGE_DIR}")
        endif()
        set(PACKAGE_ABS_DIR "${PROJECT_SOURCE_DIR}/${REPOSITORY_AND_PACKAGE_DIR}")

      endif()

      set(packageDependenciesFile "${PACKAGE_ABS_DIR}/cmake/Dependencies.cmake")
      if (EXISTS "${packageDependenciesFile}")
        set(PACKAGE_EXISTS TRUE)
      else()
        set(PACKAGE_EXISTS FALSE)
      endif()

      if (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        print_var(PROJECT_SOURCE_DIR)
        print_var(REPOSITORY_AND_PACKAGE_DIR)
        print_var(PACKAGE_ABS_DIR)
        print_var(PACKAGE_EXISTS)
        print_var(${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES)
        print_var(${TRIBITS_PACKAGE}_ALLOW_MISSING_EXTERNAL_PACKAGE)
      endif()

      if (${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES  IN_LIST
          ${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_ERROR_VALUES_LIST
        AND NOT PACKAGE_EXISTS
        AND NOT ${TRIBITS_PACKAGE}_ALLOW_MISSING_EXTERNAL_PACKAGE
        )
        message(
          "\n***"
          "\n*** Error, the package ${TRIBITS_PACKAGE} dependencies file"
	    " '${packageDependenciesFile}' does *NOT* exist!"
          "\n***\n" )
        message(FATAL_ERROR "Stopping due to above error!")
      elseif((NOT PACKAGE_EXISTS) AND (EXISTS "${PACKAGE_ABS_DIR}")
          AND (${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES STREQUAL "WARNING")
        )
        message(WARNING "${TRIBITS_PACKAGE}: Package base directory '${PACKAGE_ABS_DIR}'"
	  " exists but the dependencies file '${packageDependenciesFile}' does *NOT*"
	  " exist!  Package is being ignored anyway!")
      endif()

      if (PACKAGE_EXISTS OR ${PROJECT_NAME}_IGNORE_PACKAGE_EXISTS_CHECK)
        list(APPEND ${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES ${TRIBITS_PACKAGE})
        set(${TRIBITS_PACKAGE}_SOURCE_DIR
          "${PROJECT_SOURCE_DIR}/${REPOSITORY_AND_PACKAGE_DIR}")
        set(${TRIBITS_PACKAGE}_REL_SOURCE_DIR
          "${REPOSITORY_AND_PACKAGE_DIR}")
        set(${TRIBITS_PACKAGE}_PARENT_PACKAGE "")
        set(${TRIBITS_PACKAGE}_PARENT_REPOSITORY ${REPOSITORY_NAME})
        tribits_insert_standard_package_options(${TRIBITS_PACKAGE}  ${PACKAGE_TESTGROUP})
        set(${TRIBITS_PACKAGE}_PACKAGE_BUILD_STATUS INTERNAL)
        set(${TRIBITS_PACKAGE}_IS_TRIBITS_COMPLIANT TRUE)
      else()
        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          message(
            "\n***"
            "\n*** NOTE: Excluding package ${TRIBITS_PACKAGE} because ${PACKAGE_ABS_DIR}"
              " does not exist!"
            "\n***\n" )
        endif()
      endif()
      # NOTE: The variable ${PROJECT_NAME}_IGNORE_PACKAGE_EXISTS_CHECK only
      # gets set to TRUE for some unit tests.  Otherwise, in every legitimate
      # usage of this macro it is always FALSE.

      if (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE
          OR  ${PROJECT_NAME}_VERBOSE_CONFIGURE
        )
        print_var(${TRIBITS_PACKAGE}_SOURCE_DIR)
        print_var(${TRIBITS_PACKAGE}_REL_SOURCE_DIR)
        print_var(${TRIBITS_PACKAGE}_PARENT_PACKAGE)
        print_var(${TRIBITS_PACKAGE}_PARENT_REPOSITORY)
        print_var(${TRIBITS_PACKAGE}_PACKAGE_BUILD_STATUS)
        print_var(${TRIBITS_PACKAGE}_IS_TRIBITS_COMPLIANT)
      endif()

      if (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        print_var(${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES)
      endif()

    endforeach()

    # Get the actual number of packages that actually exist

    list(LENGTH ${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES
      ${PROJECT_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES )
    math(EXPR ${PROJECT_NAME}_LAST_DEFINED_INTERNAL_TOPLEVEL_PACKAGE_IDX
      "${${PROJECT_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES}-1")

  else()

    set(${REPOSITORY_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES 0)

  endif()

  message("-- After reading above PackagesList.cmake file: "
    "${PROJECT_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES"
    "='${${PROJECT_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES}'")

  # Print the final set of packages in debug mode

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    print_var(${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES)
  endif()

endmacro()
