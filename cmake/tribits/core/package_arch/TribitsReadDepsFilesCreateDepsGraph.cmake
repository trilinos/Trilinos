# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

include(TribitsPackageDefineDependencies)
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
# This macro reads from the variables::
#
#   ${PROJECT_NAME}_ALL_REPOSITORIES (old)
#   ${PROJECT_NAME}_PACKAGES (old)
#
# and writes to the variable::
#
#   ${PROJECT_NAME}_SE_PACKAGES (old)
#
# as well creates the package dependency variables described in `List
# variables defining the package dependencies graph`_ that defines the
# directed acyclic dependency (DAG) package dependency graph (with navigation
# up and down the graph).
#
# See `Function call tree for constructing package dependency graph`_
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
# See `Function call tree for constructing package dependency graph`_
#
macro(tribits_process_all_repository_deps_setup_files)
  foreach(TIBITS_REPO ${${PROJECT_NAME}_ALL_REPOSITORIES})
    tribits_get_repo_name_dir(${TIBITS_REPO}  REPO_NAME  REPO_DIR)
    tribits_set_base_repo_dir(${PROJECT_SOURCE_DIR}  ${REPO_DIR}  BASE_REPO_DIR)
    tribits_get_repo_name(${TIBITS_REPO} REPOSITORY_NAME)
    #print_var(TIBITS_REPO)
    #print_var(REPO_NAME)
    #print_var(REPO_DIR)
    #print_var(REPOSITORY_NAME)
    set(REPO_DEPENDENCIES_SETUP_FILE
      "${BASE_REPO_DIR}/cmake/RepositoryDependenciesSetup.cmake")
    #print_var(REPO_DEPENDENCIES_SETUP_FILE)
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
# See `Function call tree for constructing package dependency graph`_
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
# This macro reads from the variables::
#
#   ${PROJECT_NAME}_ALL_REPOSITORIES
#   ${PROJECT_NAME}_PACKAGES (old)
#
# And writes to the variable::
#
#   ${PROJECT_NAME}_SE_PACKAGES (old)
#
# as well creates the package dependency variables described in `List
# variables defining the package dependencies graph`_ that defines the
# directed acyclic dependency (DAG) package dependency graph (with navigation
# up and down the graph).
#
# See `Function call tree for constructing package dependency graph`_
#
macro(tribits_read_all_package_deps_files_create_deps_graph)

  foreach(tribitsExternalPkg  IN  LISTS  ${PROJECT_NAME}_DEFINED_TPLS)
    tribits_read_external_package_deps_files_add_to_graph(${tribitsExternalPkg})
  endforeach()

  set(${PROJECT_NAME}_SE_PACKAGES "") # Packages and subpackages

  foreach(TRIBITS_PACKAGE  IN  LISTS ${PROJECT_NAME}_PACKAGES)
    tribits_read_toplevel_package_deps_files_add_to_graph(${TRIBITS_PACKAGE}
      ${${TRIBITS_PACKAGE}_REL_SOURCE_DIR})
  endforeach()

  # Create a reverse SE packages list for later use
  set(${PROJECT_NAME}_REVERSE_SE_PACKAGES  ${${PROJECT_NAME}_SE_PACKAGES})
  if (${PROJECT_NAME}_REVERSE_SE_PACKAGES)
    list(REVERSE  ${PROJECT_NAME}_REVERSE_SE_PACKAGES)
  endif()

  list(LENGTH ${PROJECT_NAME}_SE_PACKAGES ${PROJECT_NAME}_NUM_SE_PACKAGES)
  print_var(${PROJECT_NAME}_NUM_SE_PACKAGES)

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
# varaible::
#
#   <tplName>_LIB_ALL_DEPENDENCIES
#
# See `Function call tree for constructing package dependency graph`_
#
macro(tribits_read_external_package_deps_files_add_to_graph  tplName)
  if (IS_ABSOLUTE "${${tplName}_DEPENDENCIES_FILE}")
    set(absTplDepsFile "${${tplName}_DEPENDENCIES_FILE}")
  else()
    set(absTplDepsFile "${${PROJECT_NAME}_SOURCE_DIR}/${${tplName}_DEPENDENCIES_FILE}")
  endif()
  if (EXISTS "${absTplDepsFile}")
    tribits_trace_file_processing(TPL  INCLUDE  "${absTplDepsFile}")
    include(${absTplDepsFile})
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
# ``<packageName>``.
#
# Modifies the global variables::
#
#   ${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES
#   ${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES
#   ${PACKAGE_NAME}_TEST_REQUIRED_DEP_PACKAGES
#   ${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES
#   ${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES
#   ${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES
#   ${PACKAGE_NAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES
#   ${PACKAGE_NAME}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES
#
# It also appends the list variable::
#
#   ${PROJECT_NAME}_SE_PACKAGES (old)
#
# as for the subpackage dependencies under this top-level package are read in
# order and then this top-level package is appended and dependencies are
# dependencies are created for them.
#
# See `Function call tree for constructing package dependency graph`_
#
macro(tribits_read_toplevel_package_deps_files_add_to_graph  PACKAGE_NAME)

  # A) Get ready to read in the contents of this this pakages's Dependencies.cmake file

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

  tribits_parse_subpackages_append_se_packages_add_options(${PACKAGE_NAME})

  tribits_read_package_subpackage_deps_files_add_to_graph(${PACKAGE_NAME})

  # C) Finish processing this package's dependencies into dependency graph vars
  #
  # NOTE: The subpackages for this package are automatically treated as
  # optional or required library dependent packages for this outer package!

  tribits_read_back_dependencies_vars(PARENTPACK)

  # Append the subpackages to the dependencies list if this top-level package
  set(SUBPACKAGE_IDX 0)
  foreach(TRIBITS_SUBPACKAGE ${${PACKAGE_NAME}_SUBPACKAGES})
    set(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
    list(GET ${PACKAGE_NAME}_SUBPACKAGE_OPTREQ ${SUBPACKAGE_IDX} SUBPACKAGE_OPTREQ)
    list(APPEND LIB_${SUBPACKAGE_OPTREQ}_DEP_PACKAGES ${SUBPACKAGE_FULLNAME})
    math(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")
  endforeach()

  # Append this package to list of SE packages *after* subpackages are added!
  list(APPEND ${PROJECT_NAME}_SE_PACKAGES ${PACKAGE_NAME})

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
# It also sets to empty the forward dependency list vars::
#
#    <packageName>_FORWARD_<listType>
#
# for each of the forward/downstream in `List variables defining the package
# dependencies graph`_.
#
# See `Function call tree for constructing package dependency graph`_
#
macro(tribits_prep_to_read_dependencies  PACKAGE_NAME_IN)

  tribits_declare_undefined(LIB_REQUIRED_DEP_PACKAGES)
  tribits_declare_undefined(LIB_OPTIONAL_DEP_PACKAGES)
  tribits_declare_undefined(TEST_REQUIRED_DEP_PACKAGES)
  tribits_declare_undefined(TEST_OPTIONAL_DEP_PACKAGES)

  tribits_declare_undefined(LIB_REQUIRED_DEP_TPLS "")
  tribits_declare_undefined(LIB_OPTIONAL_DEP_TPLS "")
  tribits_declare_undefined(TEST_REQUIRED_DEP_TPLS "")
  tribits_declare_undefined(TEST_OPTIONAL_DEP_TPLS "")

  set(REGRESSION_EMAIL_LIST "") # Allow to be empty

  set(${PACKAGE_NAME_IN}_FORWARD_LIB_REQUIRED_DEP_PACKAGES "")
  set(${PACKAGE_NAME_IN}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES "")
  set(${PACKAGE_NAME_IN}_FORWARD_TEST_REQUIRED_DEP_PACKAGES "")
  set(${PACKAGE_NAME_IN}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES "")

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
# See `Function call tree for constructing package dependency graph`_
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
# See `Function call tree for constructing package dependency graph`_
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
# See `Function call tree for constructing package dependency graph`_
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
# Sets up the upstream and downstream/forward package dependency list
# variables for ``<packageName>`` described in `List variables defining the
# package dependencies graph`_.  Note that the downstream/forward dependencies
# of upstream packages on this package ``<packageName>`` are built up
# incrementally.
#
# See `Function call tree for constructing package dependency graph`_
# 
macro(tribits_process_package_dependencies_lists  PACKAGE_NAME)

  tribits_set_dep_packages(${PACKAGE_NAME} LIB  REQUIRED)
  tribits_set_dep_packages(${PACKAGE_NAME} LIB  OPTIONAL)
  tribits_set_dep_packages(${PACKAGE_NAME} TEST  REQUIRED)
  tribits_set_dep_packages(${PACKAGE_NAME} TEST  OPTIONAL)

  set(${PACKAGE_NAME}_LIB_REQUIRED_DEP_TPLS ${LIB_REQUIRED_DEP_TPLS})
  set(${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS ${LIB_OPTIONAL_DEP_TPLS})
  set(${PACKAGE_NAME}_TEST_REQUIRED_DEP_TPLS ${TEST_REQUIRED_DEP_TPLS})
  set(${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS ${TEST_OPTIONAL_DEP_TPLS})

  tribits_append_forward_dep_packages(${PACKAGE_NAME} LIB_REQUIRED_DEP_PACKAGES)
  tribits_append_forward_dep_packages(${PACKAGE_NAME} LIB_OPTIONAL_DEP_PACKAGES)
  tribits_append_forward_dep_packages(${PACKAGE_NAME} TEST_REQUIRED_DEP_PACKAGES)
  tribits_append_forward_dep_packages(${PACKAGE_NAME} TEST_OPTIONAL_DEP_PACKAGES)

endmacro()


# @FUNCTION: tribits_set_dep_packages()
#
# Usage::
#
#   tribits_set_dep_packages(<packageName>  LIB|TEST  REQUIRED|OPTIONAL)
#
# Function that helps to set up backward package dependency lists for a given
# package given the vars read in from the macro
# `tribits_package_define_dependencies()`_.
#
# Sets the upstream/backward dependency variables defined in the section `List
# variables defining the package dependencies graph`_.
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
# See `Function call tree for constructing package dependency graph`_
#
function(tribits_set_dep_packages  PACKAGE_NAME   LIB_OR_TEST  REQUIRED_OR_OPTIONAL)

  if (TRIBITS_SET_DEP_PACKAGES_DEBUG_DUMP)
    message("\nTRIBITS_SET_DEP_PACKAGES:  ${PACKAGE_NAME}  ${LIB_OR_TEST}  ${REQUIRED_OR_OPTIONAL})")
  endif()

  set(LIST_TYPE  ${LIB_OR_TEST}_${REQUIRED_OR_OPTIONAL}_DEP_PACKAGES)
  set(PACKAGE_DEPS_LIST "")
  set(SE_PACKAGE_ENABLE_VAR  ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

  foreach(DEP_PKG ${${LIST_TYPE}})
    if (TRIBITS_SET_DEP_PACKAGES_DEBUG_DUMP)
      print_var(DEP_PKG)
    endif()
    if (${DEP_PKG} STREQUAL ${PACKAGE_NAME})
      tribits_abort_on_self_dep("${PACKAGE_NAME}" "${LIST_TYPE}")
    endif()
    if (${DEP_PKG}_SOURCE_DIR)
      set(DEP_PKG_DEFINED_AND_EXISTS TRUE)
    else()
      set(DEP_PKG_DEFINED_AND_EXISTS FALSE)
    endif()
    if (TRIBITS_SET_DEP_PACKAGES_DEBUG_DUMP)
      print_var(DEP_PKG_DEFINED_AND_EXISTS)
    endif()
    if (DEP_PKG_DEFINED_AND_EXISTS)
      list(APPEND PACKAGE_DEPS_LIST ${DEP_PKG})
    else()
      if (${PROJECT_NAME}_ASSERT_MISSING_PACKAGES
          AND NOT ${DEP_PKG}_ALLOW_MISSING_EXTERNAL_PACKAGE
        )
        tribits_abort_on_missing_package(
          "${DEP_PKG}" "${PACKAGE_NAME}" "${PROJECT_NAME}_SE_PACKAGES")
      else()
        if (${DEP_PKG}_ALLOW_MISSING_EXTERNAL_PACKAGE)
          if (${PROJECT_NAME}_WARN_ABOUT_MISSING_EXTERNAL_PACKAGES)
            message_wrapper("NOTE: ${DEP_PKG} is being ignored since its directory"
              " is missing and ${DEP_PKG}_ALLOW_MISSING_EXTERNAL_PACKAGE ="
              " ${${DEP_PKG}_ALLOW_MISSING_EXTERNAL_PACKAGE}!")
          endif()
          if (REQUIRED_OR_OPTIONAL STREQUAL "REQUIRED")
            message_wrapper("NOTE: Setting ${SE_PACKAGE_ENABLE_VAR}=OFF because"
              " package ${PACKAGE_NAME} has a required dependency on missing"
              " package ${DEP_PKG}!")
            dual_scope_set(${SE_PACKAGE_ENABLE_VAR} OFF)
          endif()
        endif()
        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          message(
            "\n***"
            "\n*** NOTE: The package ${DEP_PKG} which is a dependent package of"
              " ${PACKAGE_NAME} being ignored because ${DEP_PKG} is missing!"
            "\n***\n" )
        endif()
        # Must set enable vars for missing package to off so that logic in
        # existing downstream packages that key off of these vars will still
        # work.
        dual_scope_set(${PROJECT_NAME}_ENABLE_${DEP_PKG} OFF)
        dual_scope_set(${PACKAGE_NAME}_ENABLE_${DEP_PKG} OFF)
      endif()
    endif()
  endforeach()

  #print_var(PACKAGE_DEPS_LIST)

  global_set(${PACKAGE_NAME}_${LIST_TYPE} ${PACKAGE_DEPS_LIST})

endfunction()


# @FUNCTION: tribits_append_forward_dep_packages()
#
# Usage: tribits_append_forward_dep_packages(<packageName> <listType>)
#
# Function that helps to set up forward package dependency lists for an
# upstream package given that a downstream package declared a dependency on
# it.  In particular, it appends the var::
#
#    <packageName>_FORWARD_<listType>
#
# for one of the vars listed in `List variables defining the package
# dependencies graph`_.
#
# This function is called multiple times to build up the forward package
# dependencies for a given ``<packageName>`` by the downstream packages that
# declare dependencies on it.
#
# See `Function call tree for constructing package dependency graph`_
#
function(tribits_append_forward_dep_packages PACKAGE_NAME LIST_TYPE)

  set(DEP_PKG_LIST_NAME "${PACKAGE_NAME}_${LIST_TYPE}")

  assert_defined(${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)
  foreach(DEP_PKG ${${DEP_PKG_LIST_NAME}})
    set(FWD_DEP_PKG_LIST_NAME "${DEP_PKG}_FORWARD_${LIST_TYPE}")
    if (NOT DEFINED ${FWD_DEP_PKG_LIST_NAME})
      if (${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)
        tribits_abort_on_missing_package(${DEP_PKG} ${PACKAGE_NAME} ${DEP_PKG_LIST_NAME})
      else()
        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          message(
            "\n***"
            "\n*** NOTE: The package ${DEP_PKG} has forward dependent package"
              " ${PACKAGE_NAME}, but that dependency is being ignored because the package"
              " ${DEP_PKG} is missing!"
            "\n***\n" )
        endif()
      endif()
    else()
      set(${FWD_DEP_PKG_LIST_NAME} ${${FWD_DEP_PKG_LIST_NAME}} ${PACKAGE_NAME}
        PARENT_SCOPE)
    endif()
  endforeach()

endfunction()


# @MACRO: tribits_set_package_regression_email_list()
#
# Usage::
#
#  tribits_set_package_regression_email_list(<packageName>)
#
# Macro that sets a pacakge's regression email address
# ``${PACKAGE_NAME}_REGRESSION_EMAIL_LIST`` as described in ???.
#
# See `Function call tree for constructing package dependency graph`_
#
macro(tribits_set_package_regression_email_list PACKAGE_NAME)

  # Lower-case package name To be used with auto email naming based on base email address
  string(TOLOWER "${PACKAGE_NAME}" LPACKAGE)
  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    print_var(REGRESSION_EMAIL_LIST)
  endif()

  tribits_get_repo_name(${${PACKAGE_NAME}_PARENT_REPOSITORY} REPOSITORY_NAME)
  #print_var(REPOSITORY_NAME)

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
#   tribits_abort_on_missing_package(<depPkg>  <packageName> <depPkgListName>)
#
# Function that creates error message about missing/misspelled package.  This
# error message also suggests that the package might be defining an upstream
# dependency on a downstream dependency (i.e. a circular dependency).
#
# See `Function call tree for constructing package dependency graph`_
#
function(tribits_abort_on_missing_package   DEP_PKG  PACKAGE_NAME  DEP_PKG_LIST_NAME)
  multiline_set(ERRMSG
    "Error, the package '${DEP_PKG}' is listed as a dependency of the package"
    " '${PACKAGE_NAME}' is in the list '${DEP_PKG_LIST_NAME}' but the package"
    " '${DEP_PKG}' is either not defined or is listed later in the package order."
    "  This may also be an attempt to create a circular dependency between"
    " the packages '${DEP_PKG}' and '${PACKAGE_NAME}' (which is not allowed)."
    "  Check the spelling of '${DEP_PKG}' or see how it is listed in"
    " a call to tribits_repository_define_packages() in relation to"
    " '${PACKAGE_NAME}'.")
  message(FATAL_ERROR ${ERRMSG})
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
# See `Function call tree for constructing package dependency graph`_
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


# @MACRO: tribits_parse_subpackages_append_se_packages_add_options()
#
# Usage::
#
#   tribits_parse_subpackages_append_se_packages_add_options(<toplevelPackageName>)
#
# Macro that parses the read-in variable
# ``SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS`` set by the macro
# `tribits_package_define_dependencies()`_ , add subpackages to the list of
# defined packages, and define user cache var options for those subpackages.
#
# This sets the list variables for the parent package ``<toplevelPackageName>``::
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
# And it appends for each subpackage to variable::
#
#   ${PROJECT_NAME}_SE_PACKAGES (old)
#
# See `Function call tree for constructing package dependency graph`_
#
macro(tribits_parse_subpackages_append_se_packages_add_options
  PACKAGE_NAME
  )

  #message("TRIBITS_PARSE_SUBPACKAGES_APPEND_SE_PACKAGES_ADD_OPTIONS: ${PACKAGE_NAME}")

  # Structure of SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  set(SPDC_SP_NAME_OFFSET 0)
  set(SPDC_SP_DIR_OFFSET 1)
  set(SPDC_SP_CLASSIFICATION_OFFSET 2)
  set(SPDC_SP_OPTREQ_OFFSET 3)
  set(SPDC_NUM_FIELDS 4)

  set(${PACKAGE_NAME}_SUBPACKAGES "")
  set(${PACKAGE_NAME}_SUBPACKAGE_DIRS "")
  set(${PACKAGE_NAME}_SUBPACKAGE_OPTREQ "")

  if (SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS)

    list(LENGTH SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS SPDC_TOTAL_LENGTH)
    math(EXPR NUM_SUBPACKAGES "${SPDC_TOTAL_LENGTH}/${SPDC_NUM_FIELDS}")
    math(EXPR SUBPACKAGES_LAST_IDX "${NUM_SUBPACKAGES}-1")

    foreach(SUBPACKAGE_IDX RANGE ${SUBPACKAGES_LAST_IDX})

      #message("")
      #print_var(SUBPACKAGE_IDX)

      # SUBPACKAGE_NAME
      math(EXPR SUBPACKAGE_NAME_IDX
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_NAME_OFFSET}")
      list(GET SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${SUBPACKAGE_NAME_IDX}
        SUBPACKAGE_NAME )

      set(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${SUBPACKAGE_NAME})

      # SUBPACKAGE_DIR
      math(EXPR SUBPACKAGE_DIR_IDX
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_DIR_OFFSET}")
      list(GET SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${SUBPACKAGE_DIR_IDX}
        SUBPACKAGE_DIR )

      # SUBPACKAGE_CLASSIFICATION
      math(EXPR SUBPACKAGE_CLASSIFICATION_IDX
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_CLASSIFICATION_OFFSET}")
      list(GET SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${SUBPACKAGE_CLASSIFICATION_IDX}
        SUBPACKAGE_CLASSIFICATION )

      # ToDo: Parse out TESTGROUP and MATURITYLEVEL (Trilinos #6042)
      set(SUBPACKAGE_TESTGROUP ${SUBPACKAGE_CLASSIFICATION})

      tribits_update_ps_pt_ss_st(Subpackage ${SUBPACKAGE_FULLNAME} SUBPACKAGE_TESTGROUP)

      # SUBPACKAGE_OPTREQ
      math(EXPR SUBPACKAGE_OPTREQ_IDX
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_OPTREQ_OFFSET}")
      list(GET SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${SUBPACKAGE_OPTREQ_IDX}
        SUBPACKAGE_OPTREQ )

      # Determine if this subpackage exists
      set(SUBPACKAGE_FULL_SOURCE_DIR
        ${PROJECT_SOURCE_DIR}/${${PACKAGE_NAME}_REL_SOURCE_DIR}/${SUBPACKAGE_DIR})
      if (EXISTS ${SUBPACKAGE_FULL_SOURCE_DIR})
         set(SUBPACKAGE_EXISTS TRUE)
      else()
         set(SUBPACKAGE_EXISTS FALSE)
      endif()

      if (NOT SUBPACKAGE_EXISTS AND ${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)
         message(SEND_ERROR "ERROR: Subpackage dir '${SUBPACKAGE_FULL_SOURCE_DIR}'"
           " is missing!")
      endif()

      # Append to lists and global variables

      if (SUBPACKAGE_EXISTS)

        list(APPEND ${PACKAGE_NAME}_SUBPACKAGES ${SUBPACKAGE_NAME})
        list(APPEND ${PACKAGE_NAME}_SUBPACKAGE_DIRS ${SUBPACKAGE_DIR})
        list(APPEND ${PACKAGE_NAME}_SUBPACKAGE_OPTREQ ${SUBPACKAGE_OPTREQ})
        list(APPEND ${PROJECT_NAME}_SE_PACKAGES ${SUBPACKAGE_FULLNAME})
        set(${SUBPACKAGE_FULLNAME}_SOURCE_DIR "${SUBPACKAGE_FULL_SOURCE_DIR}")
        set(${SUBPACKAGE_FULLNAME}_REL_SOURCE_DIR
          "${${PACKAGE_NAME}_REL_SOURCE_DIR}/${SUBPACKAGE_DIR}")
        set(${SUBPACKAGE_FULLNAME}_PARENT_PACKAGE ${PACKAGE_NAME})
        set(${SUBPACKAGE_FULLNAME}_PARENT_REPOSITORY
          ${${PACKAGE_NAME}_PARENT_REPOSITORY})

        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          print_var(${SUBPACKAGE_FULLNAME}_PARENT_PACKAGE)
          print_var(${SUBPACKAGE_FULLNAME}_PARENT_REPOSITORY)
        endif()

        set(${SUBPACKAGE_FULLNAME}_PACKAGE_BUILD_STATUS INTERNAL)

        # Set up the input options for this subpackage
        tribits_insert_standard_package_options(${SUBPACKAGE_FULLNAME}
          ${SUBPACKAGE_TESTGROUP})

      endif()

    endforeach()

  endif()

  #print_var(${PACKAGE_NAME}_SUBPACKAGES)
  #print_var(${PACKAGE_NAME}_SUBPACKAGE_OPTREQ)

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
# See `Function call tree for constructing package dependency graph`_
#
macro(tribits_read_package_subpackage_deps_files_add_to_graph  PACKAGE_NAME)

  #message("TRIBITS_READ_PACKAGE_SUBPACKAGE_DEPS_FILES_ADD_TO_GRAPH: ${PACKAGE_NAME}")

  #print_var(${PROJECT_NAME}_SE_PACKAGES)

  set(SUBPACKAGE_IDX 0)
  foreach(TRIBITS_SUBPACKAGE ${${PACKAGE_NAME}_SUBPACKAGES})
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
# See `Function call tree for constructing package dependency graph`_
#
macro(tribits_read_subpackage_deps_file_add_to_graph  PACKAGE_NAME
  SUBPACKAGE_NAME  SUBPACKAGE_DIR
  )

  #message("TRIBITS_READ_SUBPACKAGE_DEPS_FILE_ADD_TO_GRAPH: ${PACKAGE_NAME} ${SUBPACKAGE_NAME} ${SUBPACKAGE_DIR}")

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

  if (SUBPACKAGE_EXISTS OR ${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)

    tribits_trace_file_processing(PACKAGE  INCLUDE  "${SUBPAKCAGE_DEPENDENCIES_FILE}")
    include(${SUBPAKCAGE_DEPENDENCIES_FILE})

    tribits_assert_read_dependency_vars(${SUBPACKAGE_FULLNAME})

    tribits_process_package_dependencies_lists(${SUBPACKAGE_FULLNAME})

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
