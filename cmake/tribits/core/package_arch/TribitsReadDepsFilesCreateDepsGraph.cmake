# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(TribitsPackageDefineDependencies)
include(TribitsPackageDependencies)
include(SetDefault)
include(DualScopeSet)

# @MACRO: tribits_read_deps_files_create_deps_graph()
#
# Usage::
#
#   tribits_read_deps_files_create_deps_graph()
#
# This macro reads of all the package dependencies and builds the package
# dependency graph.  This first executes the logic in the files
# `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ (for each TriBITS repo)
# and `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_ and then reads in
# all of the `<packageDir>/cmake/Dependencies.cmake`_ and
# `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_ files and builds the
# package dependency graph variables.
#
# This macro reads from the variables:
#
#   * `${PROJECT_NAME}_ALL_REPOSITORIES`_
#   * `${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES`_
#
# and writes to the variables:
#
#   * `${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES`_
#
# as well creates the package dependency variables described in `Variables
# defining the package dependencies graph`_ that defines the directed acyclic
# dependency (DAG) package dependency graph (with navigation up and down the
# graph).
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_read_deps_files_create_deps_graph)

  message("")
  message("Processing Project, Repository, and Package dependency files and building internal dependencies graph ...")
  message("")

  tribits_process_all_repository_deps_setup_files()

  tribits_process_project_dependency_setup_file()

  tribits_read_all_package_deps_files_create_deps_graph()

endmacro()


# @MACRO: tribits_process_all_repository_deps_setup_files()
#
# Process any dependency logic at the repo level by loading
# `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ files.
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_process_all_repository_deps_setup_files)
  foreach(TIBITS_REPO  IN LISTS  ${PROJECT_NAME}_ALL_REPOSITORIES)
    tribits_get_repo_name_dir(${TIBITS_REPO}  REPO_NAME  REPO_DIR)
    tribits_set_base_repo_dir(${PROJECT_SOURCE_DIR}  ${REPO_DIR}  BASE_REPO_DIR)
    tribits_get_repo_name(${TIBITS_REPO} REPOSITORY_NAME)
    set(REPO_DEPENDENCIES_SETUP_FILE
      "${BASE_REPO_DIR}/cmake/RepositoryDependenciesSetup.cmake")
    if (EXISTS ${REPO_DEPENDENCIES_SETUP_FILE})
      tribits_trace_file_processing(REPOSITORY  INCLUDE
        "${REPO_DEPENDENCIES_SETUP_FILE}")
      include(${REPO_DEPENDENCIES_SETUP_FILE})
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        print_var(${REPO_NAME}_REPOSITORY_EMAIL_URL_ADDRESS_BASE)
        print_var(${REPO_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESS)
      endif()
    else()
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("-- "
          "The ${REPO_NAME} file ${REPO_DEPENDENCIES_SETUP_FILE} does not exist! ...")
      endif()
    endif()
  endforeach()
endmacro()


# @MACRO: tribits_process_project_dependency_setup_file()
#
# Process any dependency logic at the project level by loading the
# `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_ file
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_process_project_dependency_setup_file)
  set(PROJECT_DEPENDENCIES_SETUP_FILE
    "${PROJECT_SOURCE_DIR}/cmake/ProjectDependenciesSetup.cmake")
  if (EXISTS ${PROJECT_DEPENDENCIES_SETUP_FILE})
    tribits_trace_file_processing(PROJECT  INCLUDE
      "${PROJECT_DEPENDENCIES_SETUP_FILE}")
    include(${PROJECT_DEPENDENCIES_SETUP_FILE})
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      print_var(${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESS_BASE)
      print_var(${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESS)
    endif()
  else()
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- "
        "The ${PROJECT_NAME} file ${PROJECT_DEPENDENCIES_SETUP_FILE} does not exist! ...")
    endif()
  endif()
endmacro()


# @MACRO: tribits_read_all_package_deps_files_create_deps_graph()
#
# Usage::
#
#   tribits_read_all_package_deps_files_create_deps_graph()
#
# This macro reads in all of the `<packageDir>/cmake/Dependencies.cmake`_ and
# `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_ files for top-level
# packages and subpackages, respectively, and builds the package dependency
# graph variables.
#
# This macro reads from the variables:
#
#   * `${PROJECT_NAME}_ALL_REPOSITORIES`_
#   * `${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES`_
#
# And writes to the variable:
#
#   * `${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES`_
#
# as well creates the package dependency variables described in `Variables
# defining the package dependencies graph`_ that defines the directed acyclic
# dependency (DAG) package dependency graph (with navigation up and down the
# graph).
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_read_all_package_deps_files_create_deps_graph)

  foreach(tribitsExternalPkg  IN LISTS  ${PROJECT_NAME}_DEFINED_TPLS)
    tribits_read_external_package_deps_files_add_to_graph(${tribitsExternalPkg})
  endforeach()

  set(${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES "") # Packages and subpackages

  foreach(TRIBITS_PACKAGE  IN LISTS ${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES)
    tribits_read_toplevel_package_deps_files_add_to_graph(${TRIBITS_PACKAGE})
  endforeach()

  list(LENGTH ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES
    ${PROJECT_NAME}_NUM_DEFINED_INTERNAL_PACKAGES)
  print_var(${PROJECT_NAME}_NUM_DEFINED_INTERNAL_PACKAGES)

endmacro()


# @MACRO: tribits_read_external_package_deps_files_add_to_graph()
#
# Reads in dependencies for the external packages/TPL ``<tplName>`` and
# creates the package dependency graph entries for it.
#
# Usage::
#
#   tribits_read_external_package_deps_files_add_to_graph(<tplName>)
#
# This reads in the file ``${<tplName>_DEPENDENCIES_FILE}`` and sets the
# variable::
#
#   <tplName>_LIB_DEFINED_DEPENDENCIES
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_read_external_package_deps_files_add_to_graph  tplName)
  # Set up empty lists for forward dependencies
  set(${tplName}_FORWARD_LIB_DEFINED_DEPENDENCIES "")
  set(${tplName}_FORWARD_TEST_DEFINED_DEPENDENCIES "")
  # Read in and process the external package/TPL dependency file
  if (IS_ABSOLUTE "${${tplName}_DEPENDENCIES_FILE}")
    set(absTplDepsFile "${${tplName}_DEPENDENCIES_FILE}")
  else()
    set(absTplDepsFile "${${PROJECT_NAME}_SOURCE_DIR}/${${tplName}_DEPENDENCIES_FILE}")
  endif()
  if (EXISTS "${absTplDepsFile}")
    tribits_trace_file_processing(TPL  INCLUDE  "${absTplDepsFile}")
    include(${absTplDepsFile})
    foreach(depPkg  IN LISTS  ${tplName}_LIB_DEFINED_DEPENDENCIES)
      global_set(${tplName}_LIB_DEP_REQUIRED_${depPkg}  FALSE)
    endforeach()
    tribits_append_forward_dep_packages(${tplName}  LIB)
  endif()
endmacro()


# @MACRO: tribits_read_toplevel_package_deps_files_add_to_graph()
#
# Usage::
#
#   tribits_read_toplevel_package_deps_files_add_to_graph(<packageName>)
#
# Macro that reads in package dependencies for a top-level package from the
# file `<packageDir>/cmake/Dependencies.cmake`_ and appends the forward
# dependencies list vars for packages already read in for this package
# ``<packageName>`` (see `Variables defining the package dependencies
# graph`_).
#
# It also appends the list variable:
#
#   * `${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES`_
#
# Also, the subpackage dependencies under this top-level package are read in
# order and then this top-level package is appended and dependencies are
# created for them.
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_read_toplevel_package_deps_files_add_to_graph  PACKAGE_NAME)

  # A) Get ready to read in the contents of this package's Dependencies.cmake
  # file

  tribits_prep_to_read_dependencies(${PACKAGE_NAME})

  # Listing of subpackages
  set(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS) # Allow to be empty

  # B) Read in this package's Dependencies file and save off read dependency vars.

  set(PACKAGE_DEPENDENCIES_FILE
    "${PROJECT_SOURCE_DIR}/${${PACKAGE_NAME}_REL_SOURCE_DIR}/cmake/Dependencies.cmake")

  tribits_trace_file_processing(PACKAGE  INCLUDE  "${PACKAGE_DEPENDENCIES_FILE}")
  include(${PACKAGE_DEPENDENCIES_FILE})

  tribits_assert_read_dependency_vars(${PACKAGE_NAME})

  tribits_save_off_dependency_vars(PARENTPACK)

  # B.1) Set up the mail addresses (one regression email list for the package
  # and all subpackages)

  tribits_set_package_regression_email_list(${PACKAGE_NAME})

  # B.2) Process this package's subpackages first *before* finishing this packages!

  tribits_parse_subpackages_append_packages_add_options(${PACKAGE_NAME})

  tribits_read_package_subpackage_deps_files_add_to_graph(${PACKAGE_NAME})

  # C) Finish processing this package's dependencies into dependency graph vars
  #
  # NOTE: The subpackages for this package are automatically treated as
  # optional or required library dependent packages for this outer package!

  tribits_read_back_dependencies_vars(PARENTPACK)

  # Append the subpackages to the dependencies list if this top-level package
  set(SUBPACKAGE_IDX 0)
  foreach(TRIBITS_SUBPACKAGE   IN LISTS  ${PACKAGE_NAME}_SUBPACKAGES)
    set(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
    list(GET ${PACKAGE_NAME}_SUBPACKAGE_OPTREQ ${SUBPACKAGE_IDX} SUBPACKAGE_OPTREQ)
    list(APPEND LIB_${SUBPACKAGE_OPTREQ}_DEP_PACKAGES ${SUBPACKAGE_FULLNAME})
    math(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")
  endforeach()

  # Append this package to list of packages *after* subpackages are added!
  list(APPEND ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES ${PACKAGE_NAME})

  # Process this parent package's dependency lists!
  tribits_process_package_dependencies_lists(${PACKAGE_NAME})

endmacro()


################################################################################
# Helper macros for reading in and processing dependency lists from a single
# Dependencies.cmake file.
################################################################################


# @MACRO: tribits_prep_to_read_dependencies()
#
# Usage::
#
#   tribits_prep_to_read_dependencies(<packageName>)
#
# Macro that sets to undefined all of the variables that must be set by the
# `tribits_package_define_dependencies()`_ macro.
#
# It also sets to empty the forward dependency list vars:
#
#  * `${PACKAGE_NAME}_FORWARD_LIB_DEFINED_DEPENDENCIES`_
#  * `${PACKAGE_NAME}_FORWARD_TEST_DEFINED_DEPENDENCIES`_
#
# for each of the forward/downstream package/dependency in `Variables defining
# the package dependencies graph`_.
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_prep_to_read_dependencies  PACKAGE_NAME_IN)

  # Initial vars that must be set in the Dependencies.cmake file
  tribits_declare_undefined(LIB_REQUIRED_DEP_PACKAGES)
  tribits_declare_undefined(LIB_OPTIONAL_DEP_PACKAGES)
  tribits_declare_undefined(TEST_REQUIRED_DEP_PACKAGES)
  tribits_declare_undefined(TEST_OPTIONAL_DEP_PACKAGES)

  tribits_declare_undefined(LIB_REQUIRED_DEP_TPLS)
  tribits_declare_undefined(LIB_OPTIONAL_DEP_TPLS)
  tribits_declare_undefined(TEST_REQUIRED_DEP_TPLS)
  tribits_declare_undefined(TEST_OPTIONAL_DEP_TPLS)

  set(REGRESSION_EMAIL_LIST "") # Allow to be empty

  # Initialize other vars
  set(${PACKAGE_NAME_IN}_FORWARD_LIB_DEFINED_DEPENDENCIES "")
  set(${PACKAGE_NAME_IN}_FORWARD_TEST_DEFINED_DEPENDENCIES "")

endmacro()


# @MACRO: tribits_assert_read_dependency_vars()
#
# Usage::
#
#   tribits_assert_read_dependency_vars(<packageName>)
#
# Assert that all of the required variables set by the function
# `tribits_package_define_dependencies()`_ in the file
# `<packageDir>/cmake/Dependencies.cmake`_ have been set.
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_assert_read_dependency_vars  PACKAGE_NAME)

  tribits_assert_defined_package_var(LIB_REQUIRED_DEP_PACKAGES ${PACKAGE_NAME})
  tribits_assert_defined_package_var(LIB_OPTIONAL_DEP_PACKAGES ${PACKAGE_NAME})
  tribits_assert_defined_package_var(TEST_REQUIRED_DEP_PACKAGES ${PACKAGE_NAME})
  tribits_assert_defined_package_var(TEST_OPTIONAL_DEP_PACKAGES ${PACKAGE_NAME})

  tribits_assert_defined_package_var(LIB_REQUIRED_DEP_TPLS ${PACKAGE_NAME})
  tribits_assert_defined_package_var(LIB_OPTIONAL_DEP_TPLS ${PACKAGE_NAME})
  tribits_assert_defined_package_var(TEST_REQUIRED_DEP_TPLS ${PACKAGE_NAME})
  tribits_assert_defined_package_var(TEST_OPTIONAL_DEP_TPLS ${PACKAGE_NAME})

endmacro()


# @MACRO: tribits_save_off_dependency_vars()
#
# Usage::
#
#   tribits_save_off_dependency_vars(<postfix>)
#
# Saves off package dependency variables with variable suffix ``_<postfix>``.
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_save_off_dependency_vars  POSTFIX)

  set(LIB_REQUIRED_DEP_PACKAGES_${POSTFIX} ${LIB_REQUIRED_DEP_PACKAGES})
  set(LIB_OPTIONAL_DEP_PACKAGES_${POSTFIX} ${LIB_OPTIONAL_DEP_PACKAGES})
  set(TEST_REQUIRED_DEP_PACKAGES_${POSTFIX} ${TEST_REQUIRED_DEP_PACKAGES})
  set(TEST_OPTIONAL_DEP_PACKAGES_${POSTFIX} ${TEST_OPTIONAL_DEP_PACKAGES})

  set(LIB_REQUIRED_DEP_TPLS_${POSTFIX} ${LIB_REQUIRED_DEP_TPLS})
  set(LIB_OPTIONAL_DEP_TPLS_${POSTFIX} ${LIB_OPTIONAL_DEP_TPLS})
  set(TEST_REQUIRED_DEP_TPLS_${POSTFIX} ${TEST_REQUIRED_DEP_TPLS})
  set(TEST_OPTIONAL_DEP_TPLS_${POSTFIX} ${TEST_OPTIONAL_DEP_TPLS})

endmacro()


# @MACRO: tribits_read_back_dependencies_vars()
#
# Usage::
#
#   tribits_read_back_dependencies_vars(<postfix>)
#
# Read back the local package dependency vars from the saved-off vars with
# suffix ``_<postfix>``.
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_read_back_dependencies_vars  POSTFIX)

  set(LIB_REQUIRED_DEP_PACKAGES ${LIB_REQUIRED_DEP_PACKAGES_${POSTFIX}})
  set(LIB_OPTIONAL_DEP_PACKAGES ${LIB_OPTIONAL_DEP_PACKAGES_${POSTFIX}})
  set(TEST_REQUIRED_DEP_PACKAGES ${TEST_REQUIRED_DEP_PACKAGES_${POSTFIX}})
  set(TEST_OPTIONAL_DEP_PACKAGES ${TEST_OPTIONAL_DEP_PACKAGES_${POSTFIX}})

  set(LIB_REQUIRED_DEP_TPLS ${LIB_REQUIRED_DEP_TPLS_${POSTFIX}})
  set(LIB_OPTIONAL_DEP_TPLS ${LIB_OPTIONAL_DEP_TPLS_${POSTFIX}})
  set(TEST_REQUIRED_DEP_TPLS ${TEST_REQUIRED_DEP_TPLS_${POSTFIX}})
  set(TEST_OPTIONAL_DEP_TPLS ${TEST_OPTIONAL_DEP_TPLS_${POSTFIX}})

endmacro()


# @MACRO: tribits_process_package_dependencies_lists()
#
# Usage::
#
#   tribits_process_package_dependencies_lists(<packageName>)
#
# Sets up the upstream/backward and downstream/forward package dependency list
# variables for ``<packageName>`` described in `Variables defining the package
# dependencies graph`_.  Note that the downstream/forward dependencies of
# upstream packages for this package ``<packageName>`` are built up
# incrementally.  (The forward dependency list vars are initialized to empty
# in `tribits_prep_to_read_dependencies()`_.)
#
# See `Function call tree for constructing package dependency graph`_.
# 
macro(tribits_process_package_dependencies_lists  packageName)

  # Initialize backward dep vars
  set(${packageName}_LIB_DEFINED_DEPENDENCIES "")
  set(${packageName}_TEST_DEFINED_DEPENDENCIES "")

  # Append the XXX_TPLS list on the end of the XXX_PACKAGES list
  list(APPEND LIB_REQUIRED_DEP_PACKAGES ${LIB_REQUIRED_DEP_TPLS})
  list(APPEND LIB_OPTIONAL_DEP_PACKAGES ${LIB_OPTIONAL_DEP_TPLS})
  list(APPEND TEST_REQUIRED_DEP_PACKAGES ${TEST_REQUIRED_DEP_TPLS})
  list(APPEND TEST_OPTIONAL_DEP_PACKAGES ${TEST_OPTIONAL_DEP_TPLS})
  set(LIB_REQUIRED_DEP_TPLS "")
  set(LIB_OPTIONAL_DEP_TPLS "")
  set(TEST_REQUIRED_DEP_TPLS "")
  set(TEST_OPTIONAL_DEP_TPLS "")

  # Fill the backward dependency vars
  tribits_set_dep_packages(${packageName} LIB  REQUIRED  PACKAGES)
  tribits_set_dep_packages(${packageName} LIB  OPTIONAL  PACKAGES)
  tribits_set_dep_packages(${packageName} TEST  REQUIRED  PACKAGES)
  tribits_set_dep_packages(${packageName} TEST  OPTIONAL  PACKAGES)

  # Fill forward deps lists #63
  tribits_append_forward_dep_packages(${packageName}  LIB)
  tribits_append_forward_dep_packages(${packageName}  TEST)

endmacro()


# @MACRO: tribits_set_dep_packages()
#
# Macro set up backward package dependency lists for a given package given the
# vars read in from the macro `tribits_package_define_dependencies()`_.
#
# Usage::
#
#   tribits_set_dep_packages(<packageName> <testOrLib> <requiredOrOptional> <pkgsOrTpls>)
#
# where:
#
# * ``<testOrLib>``: ``LIB`` or ``TEST``
# * ``<requiredOrOptional>``: ``REQUIRED`` or ``OPTIONAL``
# * ``<pkgsOrTpls>``: ``PACKAGES`` or ``TPLS``
#
#
# Sets the upstream/backward dependency variables defined in the section
# `Variables defining the package dependencies graph`_.
#
# This also handles the several types of issues:
#
# * A package declaring a dependency on itself
#   (`tribits_abort_on_self_dep()`_).
#
# * A missing upstream dependent package (either error out with
#   `tribits_abort_on_missing_package()`_ or allow to be missing and disable
#   this package if this is a required dependency).
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_set_dep_packages  packageName  testOrLib  requiredOrOptional  pkgsOrTpls)

  set(inputListType  ${testOrLib}_${requiredOrOptional}_DEP_${pkgsOrTpls})
  set(packageEnableVar  ${PROJECT_NAME}_ENABLE_${packageName})

  foreach(depPkg  IN LISTS  ${inputListType})
    if (${depPkg} STREQUAL ${packageName})
      tribits_abort_on_self_dep("${packageName}" "${inputListType}")
    endif()
    tribits_is_pkg_defined(${depPkg} depPkgIsDefined)
    if (depPkgIsDefined)
      list(APPEND ${packageName}_${testOrLib}_DEFINED_DEPENDENCIES ${depPkg})
      if ("${requiredOrOptional}"  STREQUAL  "REQUIRED")
        global_set(${packageName}_${testOrLib}_DEP_REQUIRED_${depPkg}  TRUE)
      elseif ("${requiredOrOptional}"  STREQUAL  "OPTIONAL")
        global_set(${packageName}_${testOrLib}_DEP_REQUIRED_${depPkg}  FALSE)
      else()
        message(FATAL_ERROR
          "Invalid value for requiredOrOptional='${requiredOrOptional}'!")
      endif()
    else()
      tribits_set_dep_packages__handle_undefined_pkg(${packageName} ${depPkg}
        ${requiredOrOptional} ${pkgsOrTpls} ${packageEnableVar})
    endif()
  endforeach()

endmacro()


# Determine if a (internal or external) package is defined or not
#
function(tribits_is_pkg_defined  depPkg    depPkgIsDefinedOut)
  set(depPkgIsDefined  FALSE)
  if (${depPkg}_SOURCE_DIR)
    set(depPkgIsDefined  TRUE)
  elseif(${depPkg}_FINDMOD)
    set(depPkgIsDefined  TRUE)
  endif()
  set(${depPkgIsDefinedOut} ${depPkgIsDefined} PARENT_SCOPE)
endfunction()


# Implementation macro for tribits_set_dep_packages() to deal with a package
# that is not defined by TriBITS.
#
# ToDo #63: This may need to be modified when dealing with TriBITS-compliant
# packages already installed out on the system.  We may need a mode where we
# don't assert packages that are not defined but instead just assume they are
# TriBITS-compliant packages already installed.
#
macro(tribits_set_dep_packages__handle_undefined_pkg  packageName  depPkg
    requiredOrOptional  pkgsOrTpls  packageEnableVar
  )
  # Determine if it is allowed for this depPkg to not be defined
  set(errorOutForUndefinedDepPkg  TRUE)
  if (${depPkg}_ALLOW_MISSING_EXTERNAL_PACKAGE)
    set(errorOutForUndefinedDepPkg  FALSE)
  elseif (NOT  ${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES  IN_LIST
      ${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_ERROR_VALUES_LIST
    )
    set(errorOutForUndefinedDepPkg  FALSE)
  endif()
  # Produce error or deal with allowed undefined ${depPkg}
  if (errorOutForUndefinedDepPkg)
    tribits_abort_on_missing_package("${depPkg}" "${packageName}")
  else()
    if (${depPkg}_ALLOW_MISSING_EXTERNAL_PACKAGE)
      if (${PROJECT_NAME}_WARN_ABOUT_MISSING_EXTERNAL_PACKAGES)
        message_wrapper("NOTE: ${depPkg} is being ignored since its directory"
          " is missing and ${depPkg}_ALLOW_MISSING_EXTERNAL_PACKAGE ="
          " ${${depPkg}_ALLOW_MISSING_EXTERNAL_PACKAGE}!")
      endif()
      if ("${requiredOrOptional}" STREQUAL "REQUIRED")
        message_wrapper("NOTE: Setting ${packageEnableVar}=OFF because"
          " package ${packageName} has a required dependency on missing"
          " package ${depPkg}!")
        set(${packageEnableVar} OFF)
      endif()
    endif()
    # Must set enable vars for missing package to off so that logic in
    # existing downstream packages that key off of these vars will still
    # work.
    set(${PROJECT_NAME}_ENABLE_${depPkg} OFF)
    set(${packageName}_ENABLE_${depPkg} OFF)
  endif()
endmacro()


# @MACRO: tribits_append_forward_dep_packages()
#
# Appends forward/downstream package dependency lists for the upstream
# dependent package list provided.
#
# Usage::
#
#   tribits_append_forward_dep_packages(<packageName> <listType>)
#
# In particular, it appends the var::
#
#    <packageName>_FORWARD_<listType>
#
# for one of the vars listed in `Variables defining the package dependencies
# graph`_.
#
# This function is called multiple times to build up the forward package
# dependencies for a given ``<packageName>`` by the downstream packages that
# declare dependencies on it.
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_append_forward_dep_packages  packageName  libOrTest)

  foreach(depPkg  IN LISTS  ${packageName}_${libOrTest}_DEFINED_DEPENDENCIES)
    set(fwdDepPkgListName ${depPkg}_FORWARD_${libOrTest}_DEFINED_DEPENDENCIES)
    if (DEFINED ${fwdDepPkgListName})
      list(APPEND ${fwdDepPkgListName} ${packageName})
    else()
      if (${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES  IN_LIST
          ${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_ERROR_VALUES_LIST
        )
        tribits_abort_on_missing_package(${depPkg} ${packageName})
      else()
        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          message(
            "\n***"
            "\n*** NOTE: The package ${depPkg} has forward dependent package"
              " ${packageName}, but that dependency is being ignored because the package"
              " ${depPkg} is missing!"
            "\n***\n" )
        endif()
      endif()
    endif()
  endforeach()

endmacro()


# @MACRO: tribits_set_package_regression_email_list()
#
# Usage::
#
#  tribits_set_package_regression_email_list(<packageName>)
#
# Macro that sets a pacakge's regression email address
# ``${PACKAGE_NAME}_REGRESSION_EMAIL_LIST`` as described in ???.
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_set_package_regression_email_list  PACKAGE_NAME)

  # Lower-case package name To be used with auto email naming based on base email address
  string(TOLOWER "${PACKAGE_NAME}" LPACKAGE)
  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    print_var(REGRESSION_EMAIL_LIST)
  endif()

  tribits_get_repo_name(${${PACKAGE_NAME}_PARENT_REPOSITORY} REPOSITORY_NAME)

  if(${REPOSITORY_NAME}_REPOSITORY_OVERRIDE_PACKAGE_EMAIL_LIST)
    set(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST
      ${${REPOSITORY_NAME}_REPOSITORY_OVERRIDE_PACKAGE_EMAIL_LIST})
  elseif (REGRESSION_EMAIL_LIST)
    set(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST ${REGRESSION_EMAIL_LIST})
  elseif (${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESS_BASE)
    set(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST
      "${LPACKAGE}-regression@${${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESS_BASE}")
  elseif (${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESS)
    set(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST
      "${${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESS}")
  elseif (${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESS_BASE)
    set(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST
      "${LPACKAGE}-regression@${${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESS_BASE}")
  elseif (${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESS)
    set(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST
      "${${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESS}")
  else()
    set(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST "")
  endif()

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    print_var(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST)
  endif()

endmacro()


# @FUNCTION: tribits_abort_on_missing_package()
#
# Usage::
#
#   tribits_abort_on_missing_package(<depPkg>  <packageName>)
#
# Function that creates error message about missing/misspelled package.  This
# error message also suggests that the package might be defining an upstream
# dependency on a downstream dependency (i.e. a circular dependency).
#
# See `Function call tree for constructing package dependency graph`_.
#
function(tribits_abort_on_missing_package   depPkg  packageName)
  multiline_set(ERRMSG
    "Error, the package '${depPkg}' is listed as a dependency of the package"
    " '${packageName}' but the package '${depPkg}' is either not defined or"
    " is listed later in the package order."
    "  This may also be an attempt to create a circular dependency between"
    " the packages '${depPkg}' and '${packageName}' (which is not allowed)."
    "  Check the spelling of '${depPkg}' or see how it is listed in"
    " a call to tribits_repository_define_packages() in relation to"
    " '${packageName}'."
    "  To ignore/disable the undefined package '${depPkg}', set the cache"
    " variable ${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES=IGNORE."
    )
  message(${${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES} "${ERRMSG}")
endfunction()


# @FUNCTION: tribits_abort_on_self_dep()
#
# Usage::
#
#   tribits_abort_on_self_dep(<packageName> <depPkgListName>)
#
# Prints a fatal error message for an attempt for a self dependency
# declaration and which list it comes from.
#
# See `Function call tree for constructing package dependency graph`_.
#
function(tribits_abort_on_self_dep  PACKAGE_NAME  DEP_PKG_LIST_NAME)
  multiline_set(ERRMSG
    "Error, the package '${PACKAGE_NAME}' is listed as a dependency of itself"
    " in the list '${DEP_PKG_LIST_NAME}'!")
  message(FATAL_ERROR ${ERRMSG})
endfunction()


################################################################################
# Macros/functions for processing dependency info for subpackages
################################################################################


# @MACRO: tribits_parse_subpackages_append_packages_add_options()
#
# Usage::
#
#   tribits_parse_subpackages_append_packages_add_options(<parentPackageName>)
#
# Macro that parses the read-in variable
# ``SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS`` set by the macro
# `tribits_package_define_dependencies()`_ , adds subpackages to the list of
# defined packages, and defines user cache var options for those subpackages.
#
# This sets the list variables for the parent package
# ``<parentPackageName>``::
#
#   <parentPackageName>_SUBPACKAGES
#   <parentPackageName>_SUBPACKAGE_DIRS
#   <parentPackageName>_SUBPACKAGE_OPTREQ
#
# For each subpackage ``<subpackageFullName>``, this sets::
#
#   <subpackageFullName>_SOURCE_DIR
#   <subpackageFullName>_REL_SOURCE_DIR
#   <subpackageFullName>_PARENT_PACKAGE
#   <subpackageFullName>_PARENT_REPOSITORY
#
# And it appends each subpackage to the list variable:
#
#   * `${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES`_
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_parse_subpackages_append_packages_add_options  parentPackageName)

  # Fields in the list var SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  set(SPDC_SP_NAME_OFFSET 0)
  set(SPDC_SP_DIR_OFFSET 1)
  set(SPDC_SP_CLASSIFICATION_OFFSET 2)
  set(SPDC_SP_OPTREQ_OFFSET 3)
  set(SPDC_NUM_FIELDS 4) # Number of the above files

  set(${parentPackageName}_SUBPACKAGES "")
  set(${parentPackageName}_SUBPACKAGE_DIRS "")
  set(${parentPackageName}_SUBPACKAGE_OPTREQ "")

  if (SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS)

    list(LENGTH  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS  SPDC_TOTAL_LENGTH)
    math(EXPR  numSubpackages  "${SPDC_TOTAL_LENGTH}/${SPDC_NUM_FIELDS}")
    math(EXPR  subpackagesLastIdx  "${numSubpackages}-1")

    foreach(SUBPACKAGE_IDX  RANGE  ${subpackagesLastIdx})

      # subpkgName
      math(EXPR  subpkgNameIdx
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_NAME_OFFSET}")
      list(GET  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${subpkgNameIdx} subpkgName)

      # subpkgFullname
      set(subpkgFullname ${parentPackageName}${subpkgName})

      # subpkgDir
      math(EXPR  subpkgDirIdx
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_DIR_OFFSET}")
      list(GET  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${subpkgDirIdx} subpkgDir)

      # subpkgClassification
      math(EXPR  subpkgClassificationIdx
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_CLASSIFICATION_OFFSET}")
      list(GET  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${subpkgClassificationIdx}
        subpkgClassification )

      # ToDo: Parse out TESTGROUP and MATURITYLEVEL (Trilinos #6042)
      set(subpkgTestgroup ${subpkgClassification})

      tribits_update_ps_pt_ss_st(Subpackage ${subpkgFullname} subpkgTestgroup)

      # subpkgOptreq
      math(EXPR  subpkgOptreqIdx
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_OPTREQ_OFFSET}")
      list(GET  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${subpkgOptreqIdx}
        subpkgOptreq )

      # Determine if this subpackage exists
      set(subpkgFullSrcDir
        "${PROJECT_SOURCE_DIR}/${${parentPackageName}_REL_SOURCE_DIR}/${subpkgDir}")
      if (EXISTS "${subpkgFullSrcDir}")
         set(subpkgExists TRUE)
      else()
         set(subpkgExists FALSE)
      endif()

      if (NOT  subpkgExists  AND  ${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES  IN_LIST
          ${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_ERROR_VALUES_LIST
         )
         message(SEND_ERROR "ERROR: Subpackage dir '${subpkgFullSrcDir}' is missing!")
      endif()

      # Append to lists and global variables

      if (subpkgExists)

        list(APPEND ${parentPackageName}_SUBPACKAGES ${subpkgName})
        list(APPEND ${parentPackageName}_SUBPACKAGE_DIRS ${subpkgDir})
        list(APPEND ${parentPackageName}_SUBPACKAGE_OPTREQ ${subpkgOptreq})
        list(APPEND ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES ${subpkgFullname})
        set(${subpkgFullname}_SOURCE_DIR "${subpkgFullSrcDir}")
        set(${subpkgFullname}_REL_SOURCE_DIR
          "${${parentPackageName}_REL_SOURCE_DIR}/${subpkgDir}")
        set(${subpkgFullname}_PARENT_PACKAGE ${parentPackageName})
        set(${subpkgFullname}_PARENT_REPOSITORY
          ${${parentPackageName}_PARENT_REPOSITORY})

        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          print_var(${subpkgFullname}_PARENT_PACKAGE)
          print_var(${subpkgFullname}_PARENT_REPOSITORY)
        endif()

        set(${subpkgFullname}_PACKAGE_BUILD_STATUS INTERNAL)

        # Set up the input options for this subpackage
        tribits_insert_standard_package_options(${subpkgFullname} ${subpkgTestgroup})

      endif()

    endforeach()

  endif()

endmacro()


# @MACRO: tribits_read_package_subpackage_deps_files_add_to_graph()
#
# Usage::
#
#   tribits_read_package_subpackage_deps_files_add_to_graph(<toplevelPackageName>)
#
# Read in subpackages dependencies files and add to dependencies graph
# variables.
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_read_package_subpackage_deps_files_add_to_graph  PACKAGE_NAME)

  set(SUBPACKAGE_IDX 0)
  foreach(TRIBITS_SUBPACKAGE  IN LISTS  ${PACKAGE_NAME}_SUBPACKAGES)
    list(GET ${PACKAGE_NAME}_SUBPACKAGE_DIRS ${SUBPACKAGE_IDX} SUBPACKAGE_DIR)
    tribits_read_subpackage_deps_file_add_to_graph(${TRIBITS_PACKAGE}
      ${TRIBITS_SUBPACKAGE}  ${SUBPACKAGE_DIR})
    math(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")
  endforeach()

endmacro()


# @MACRO: tribits_read_subpackage_deps_file_add_to_graph()
#
# Usage::
#
#   tribits_read_subpackage_deps_file_add_to_graph(<toplevelPackageName>
#     <subpackageName> <subpackageDir>)
#
# Macro that reads in a single subpackage dependencies file
# `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_ and sets up the
# dependency structure for it.
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_read_subpackage_deps_file_add_to_graph  PACKAGE_NAME
  SUBPACKAGE_NAME  SUBPACKAGE_DIR
  )

  set(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${SUBPACKAGE_NAME})

  #
  # A) Get ready to read in the contents of this this subpackage's
  # Dependencies.cmake file
  #

  tribits_prep_to_read_dependencies(${SUBPACKAGE_FULLNAME})

  # NOTE: Subpackages use the regression email list from the parent package.

  # NOTE: Subpackages are not allowed to have subpackages!
  set(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS)

  #
  # B) Read in this subpackage's Dependencies file
  #

  set(SUBPACKAGE_FULL_DIR "${${PACKAGE_NAME}_REL_SOURCE_DIR}/${SUBPACKAGE_DIR}")

  set(SUBPACKAGE_ABS_DIR "${PROJECT_SOURCE_DIR}/${SUBPACKAGE_FULL_DIR}")
  set(SUBPAKCAGE_DEPENDENCIES_FILE "${SUBPACKAGE_ABS_DIR}/cmake/Dependencies.cmake")

  if (EXISTS ${SUBPAKCAGE_DEPENDENCIES_FILE})
    set(SUBPACKAGE_EXISTS TRUE)
  else()
    set(SUBPACKAGE_EXISTS FALSE)
  endif()

  if (SUBPACKAGE_EXISTS  OR  ${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES  IN_LIST
      ${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_ERROR_VALUES_LIST
    )

    tribits_trace_file_processing(PACKAGE  INCLUDE  "${SUBPAKCAGE_DEPENDENCIES_FILE}")
    include(${SUBPAKCAGE_DEPENDENCIES_FILE})

    tribits_assert_read_dependency_vars(${SUBPACKAGE_FULLNAME})

    tribits_process_package_dependencies_lists(${SUBPACKAGE_FULLNAME})

    set(${SUBPACKAGE_FULLNAME}_IS_TRIBITS_COMPLIANT TRUE)

    set(${SUBPACKAGE_FULLNAME}_REGRESSION_EMAIL_LIST
      ${${PACKAGE_NAME}_REGRESSION_EMAIL_LIST})

  endif()

endmacro()


#
# Private utility functions
#


# Function that sets a variable to DECLARED-UNDEFINED
#
function(tribits_declare_undefined  VAR_NAME)
  set(${VAR_NAME}  DECLARED-UNDEFINED  PARENT_SCOPE)
endfunction()


# Function that asserts that a package dependency variable is defined
# correctly
#
function(tribits_assert_defined_package_var  PACKAGE_VAR  PACKAGE_NAME)
  if (${PACKAGE_VAR} STREQUAL DECLARED-UNDEFINED)
    message(FATAL_ERROR
      "Error, the package variable ${PACKAGE_VAR} was not defined correctly for package ${PACKAGE_NAME}!"
      )
  endif()
endfunction()
