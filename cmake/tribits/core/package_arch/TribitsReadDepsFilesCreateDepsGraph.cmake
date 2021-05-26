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

INCLUDE(TribitsPackageDefineDependencies)
INCLUDE(SetDefault)
INCLUDE(DualScopeSet)

# @MACRO: TRIBITS_READ_DEPS_FILES_CREATE_DEPS_GRAPH()
#
# Usage::
#
#   TRIBITS_READ_DEPS_FILES_CREATE_DEPS_GRAPH()
#
# This macro reads of all the package dependencies and builds the package
# dependency graph.  This first executes the logic in the files
# `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ (for each TriBITS repo)
# and `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_ and then reads in
# all of the `<packageDir>/cmake/Dependencies.cmake`_ and
# `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_ files and builds the
# package depenency graph varibles.
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
# as well creates the package dependency varaibles described in `List
# variables defining the package dependencies graph`_ that defines the
# directed acyclic depenency (DAG) package dependency graph (with navigation
# up and down the graph).
#
# See `Function call tree for constructing package dependency graph`_
#
MACRO(TRIBITS_READ_DEPS_FILES_CREATE_DEPS_GRAPH)

  MESSAGE("")
  MESSAGE("Processing Project, Repository, and Package dependency files and building internal dependencies graph ...")
  MESSAGE("")

  TRIBITS_PROCESS_ALL_REPOSITORY_DEPS_SETUP_FILES()

  TRIBITS_PROCESS_PROJECT_DEPENDENCY_SETUP_FILE()

  TRIBITS_READ_ALL_PACKAGE_DEPS_FILES_CREATE_DEPS_GRAPH()

ENDMACRO()


# @MACRO: TRIBITS_PROCESS_ALL_REPOSITORY_DEPS_SETUP_FILES()
#
# Process any dependency logic at the repo level by loading
# `<repoDir>/cmake/RepositoryDependenciesSetup.cmake`_ files.
#
# See `Function call tree for constructing package dependency graph`_
#
MACRO(TRIBITS_PROCESS_ALL_REPOSITORY_DEPS_SETUP_FILES)
  FOREACH(TIBITS_REPO ${${PROJECT_NAME}_ALL_REPOSITORIES})
    TRIBITS_GET_REPO_NAME_DIR(${TIBITS_REPO}  REPO_NAME  REPO_DIR)
    TRIBITS_SET_BASE_REPO_DIR(${PROJECT_SOURCE_DIR}  ${REPO_DIR}  BASE_REPO_DIR)
    TRIBITS_GET_REPO_NAME(${TIBITS_REPO} REPOSITORY_NAME)
    #PRINT_VAR(TIBITS_REPO)
    #PRINT_VAR(REPO_NAME)
    #PRINT_VAR(REPO_DIR)
    #PRINT_VAR(REPOSITORY_NAME)
    SET(REPO_DEPENDENCIES_SETUP_FILE
      "${BASE_REPO_DIR}/cmake/RepositoryDependenciesSetup.cmake")
    #PRINT_VAR(REPO_DEPENDENCIES_SETUP_FILE)
    IF (EXISTS ${REPO_DEPENDENCIES_SETUP_FILE})
      TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  INCLUDE
        "${REPO_DEPENDENCIES_SETUP_FILE}")
      INCLUDE(${REPO_DEPENDENCIES_SETUP_FILE})
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        PRINT_VAR(${REPO_NAME}_REPOSITORY_EMAIL_URL_ADDRESS_BASE)
        PRINT_VAR(${REPO_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESS)
      ENDIF()
    ELSE()
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("-- "
          "The ${REPO_NAME} file ${REPO_DEPENDENCIES_SETUP_FILE} does not exist! ...")
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDMACRO()


# @MACRO: TRIBITS_PROCESS_PROJECT_DEPENDENCY_SETUP_FILE()
#
# Process any dependency logic at the project level by loading the
# `<projectDir>/cmake/ProjectDependenciesSetup.cmake`_ file
#
# See `Function call tree for constructing package dependency graph`_
#
MACRO(TRIBITS_PROCESS_PROJECT_DEPENDENCY_SETUP_FILE)
  SET(PROJECT_DEPENDENCIES_SETUP_FILE
    "${PROJECT_SOURCE_DIR}/cmake/ProjectDependenciesSetup.cmake")
  IF (EXISTS ${PROJECT_DEPENDENCIES_SETUP_FILE})
    TRIBITS_TRACE_FILE_PROCESSING(PROJECT  INCLUDE
      "${PROJECT_DEPENDENCIES_SETUP_FILE}")
    INCLUDE(${PROJECT_DEPENDENCIES_SETUP_FILE})
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESS_BASE)
      PRINT_VAR(${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESS)
    ENDIF()
  ELSE()
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("-- "
        "The ${PROJECT_NAME} file ${PROJECT_DEPENDENCIES_SETUP_FILE} does not exist! ...")
    ENDIF()
  ENDIF()
ENDMACRO()


# @MACRO: TRIBITS_READ_ALL_PACKAGE_DEPS_FILES_CREATE_DEPS_GRAPH()
#
# Usage::
#
#   TRIBITS_READ_ALL_PACKAGE_DEPS_FILES_CREATE_DEPS_GRAPH()
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
MACRO(TRIBITS_READ_ALL_PACKAGE_DEPS_FILES_CREATE_DEPS_GRAPH)

  SET(${PROJECT_NAME}_SE_PACKAGES) # Packages and subpackages

  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES})
    TRIBITS_READ_TOPLEVEL_PACKAGE_DEPS_FILES_ADD_TO_GRAPH(${TRIBITS_PACKAGE}
      ${${TRIBITS_PACKAGE}_REL_SOURCE_DIR})
  ENDFOREACH()

  # Create a reverse SE packages list for later use
  SET(${PROJECT_NAME}_REVERSE_SE_PACKAGES ${${PROJECT_NAME}_SE_PACKAGES})
  IF (${PROJECT_NAME}_REVERSE_SE_PACKAGES)
    LIST(REVERSE ${PROJECT_NAME}_REVERSE_SE_PACKAGES)
  ENDIF()

  LIST(LENGTH ${PROJECT_NAME}_SE_PACKAGES ${PROJECT_NAME}_NUM_SE_PACKAGES)
  PRINT_VAR(${PROJECT_NAME}_NUM_SE_PACKAGES)

ENDMACRO()


# @MACRO: TRIBITS_READ_TOPLEVEL_PACKAGE_DEPS_FILES_ADD_TO_GRAPH()
#
# Usage::
#
#  TRIBITS_READ_TOPLEVEL_PACKAGE_DEPS_FILES_ADD_TO_GRAPH(<packageName>)
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
# It also appends the list varaible::
#
#   ${PROJECT_NAME}_SE_PACKAGES (old)
#
# as for the subpackage dependencies under this top-level package are read in
# order and then this top-level package is appended and dependencies are
# dependencies are created for them.
#
# See `Function call tree for constructing package dependency graph`_
#
MACRO(TRIBITS_READ_TOPLEVEL_PACKAGE_DEPS_FILES_ADD_TO_GRAPH  PACKAGE_NAME)

  # A) Get ready to read in the contents of this this pakages's Dependencies.cmake file

  TRIBITS_PREP_TO_READ_DEPENDENCIES(${PACKAGE_NAME})

  # Listing of subpakages
  SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS) # Allow to be empty

  # B) Read in this package's Dependencies file and save off read dependency vars.

  SET(PACKAGE_DEPENDENCIES_FILE
    "${PROJECT_SOURCE_DIR}/${${PACKAGE_NAME}_REL_SOURCE_DIR}/cmake/Dependencies.cmake")

  TRIBITS_TRACE_FILE_PROCESSING(PACKAGE  INCLUDE  "${PACKAGE_DEPENDENCIES_FILE}")
  INCLUDE(${PACKAGE_DEPENDENCIES_FILE})

  TRIBITS_ASSERT_READ_DEPENDENCY_VARS(${PACKAGE_NAME})

  TRIBITS_SAVE_OFF_DEPENDENCIES_VARS(PARENTPACK)

  # B.1) Set up the mail addresses (one regression email list for the package
  # and all subpackages)

  TRIBITS_SET_PACAKGE_REGRESSION_EMAIL_LIST(${PACKAGE_NAME})

  # B.2) Process this package's subpackages first *before* finishing this packages!

  TRIBITS_PARSE_SUBPACKAGES_APPEND_SE_PACKAGES_ADD_OPTIONS(${PACKAGE_NAME})

  TRIBITS_READ_PACKAGE_SUBPACKAGE_DEPS_FILES_ADD_TO_GRAPH(${PACKAGE_NAME})

  # C) Finish processing this package's dependencies into dependency graph vars
  #
  # NOTE: The subpackages for this package are automatically treated as
  # optional or required library dependent packages for this outer package!

  TRIBITS_READ_BACK_DEPENDENCIES_VARS(PARENTPACK)

  # Append the subpackages to the dependencies list if this top-level package
  SET(SUBPACKAGE_IDX 0)
  FOREACH(TRIBITS_SUBPACKAGE ${${PACKAGE_NAME}_SUBPACKAGES})
    SET(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
    LIST(GET ${PACKAGE_NAME}_SUBPACKAGE_OPTREQ ${SUBPACKAGE_IDX} SUBPACKAGE_OPTREQ)
    LIST(APPEND LIB_${SUBPACKAGE_OPTREQ}_DEP_PACKAGES ${SUBPACKAGE_FULLNAME})
    MATH(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")
  ENDFOREACH()

  # Append this package to list of SE packages *after* subpackages are added!
  LIST(APPEND ${PROJECT_NAME}_SE_PACKAGES ${PACKAGE_NAME})

  # Process this parent package's dependency lists!
  TRIBITS_PROCESS_PACKAGE_DEPENDENCIES_LISTS(${PACKAGE_NAME})

ENDMACRO()


################################################################################
# Helper macros for reading in and processing dependency lists from a single
# Dependencies.cmake file.
################################################################################


# @MACRO: TRIBITS_PREP_TO_READ_DEPENDENCIES()
#
# Usage::
#
#   TRIBITS_PREP_TO_READ_DEPENDENCIES(<packageName>)
#
# Macro that sets to undefined all of the variables that must be set by the
# `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_ macro.
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
MACRO(TRIBITS_PREP_TO_READ_DEPENDENCIES  PACKAGE_NAME_IN)

  TRIBITS_DECLARE_UNDEFINED(LIB_REQUIRED_DEP_PACKAGES)
  TRIBITS_DECLARE_UNDEFINED(LIB_OPTIONAL_DEP_PACKAGES)
  TRIBITS_DECLARE_UNDEFINED(TEST_REQUIRED_DEP_PACKAGES)
  TRIBITS_DECLARE_UNDEFINED(TEST_OPTIONAL_DEP_PACKAGES)

  TRIBITS_DECLARE_UNDEFINED(LIB_REQUIRED_DEP_TPLS "")
  TRIBITS_DECLARE_UNDEFINED(LIB_OPTIONAL_DEP_TPLS "")
  TRIBITS_DECLARE_UNDEFINED(TEST_REQUIRED_DEP_TPLS "")
  TRIBITS_DECLARE_UNDEFINED(TEST_OPTIONAL_DEP_TPLS "")

  SET(REGRESSION_EMAIL_LIST "") # Allow to be empty

  SET(${PACKAGE_NAME_IN}_FORWARD_LIB_REQUIRED_DEP_PACKAGES "")
  SET(${PACKAGE_NAME_IN}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES "")
  SET(${PACKAGE_NAME_IN}_FORWARD_TEST_REQUIRED_DEP_PACKAGES "")
  SET(${PACKAGE_NAME_IN}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES "")

ENDMACRO()


# @MACRO: TRIBITS_ASSERT_READ_DEPENDENCY_VARS()
#
# Usage::
#
#   TRIBITS_ASSERT_READ_DEPENDENCY_VARS(<packageName>)
#
# Assert that all of the required variables set by the function
# `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_ in the file
# `<packageDir>/cmake/Dependencies.cmake`_ have been set.
#
# See `Function call tree for constructing package dependency graph`_
#
MACRO(TRIBITS_ASSERT_READ_DEPENDENCY_VARS  PACKAGE_NAME)

  TRIBITS_ASSERT_DEFINED_PACKAGE_VAR(LIB_REQUIRED_DEP_PACKAGES ${PACKAGE_NAME})
  TRIBITS_ASSERT_DEFINED_PACKAGE_VAR(LIB_OPTIONAL_DEP_PACKAGES ${PACKAGE_NAME})
  TRIBITS_ASSERT_DEFINED_PACKAGE_VAR(TEST_REQUIRED_DEP_PACKAGES ${PACKAGE_NAME})
  TRIBITS_ASSERT_DEFINED_PACKAGE_VAR(TEST_OPTIONAL_DEP_PACKAGES ${PACKAGE_NAME})

  TRIBITS_ASSERT_DEFINED_PACKAGE_VAR(LIB_REQUIRED_DEP_TPLS ${PACKAGE_NAME})
  TRIBITS_ASSERT_DEFINED_PACKAGE_VAR(LIB_OPTIONAL_DEP_TPLS ${PACKAGE_NAME})
  TRIBITS_ASSERT_DEFINED_PACKAGE_VAR(TEST_REQUIRED_DEP_TPLS ${PACKAGE_NAME})
  TRIBITS_ASSERT_DEFINED_PACKAGE_VAR(TEST_OPTIONAL_DEP_TPLS ${PACKAGE_NAME})

ENDMACRO()


# @MACRO: TRIBITS_SAVE_OFF_DEPENDENCIES_VARS()
#
# Usage::
#
#   TRIBITS_SAVE_OFF_DEPENDENCIES_VARS(<postfix>)
#
# Saves off package depeneency varaibles with variable suffix ``_<postfix>``.
#
# See `Function call tree for constructing package dependency graph`_
#
MACRO(TRIBITS_SAVE_OFF_DEPENDENCIES_VARS  POSTFIX)

  SET(LIB_REQUIRED_DEP_PACKAGES_${POSTFIX} ${LIB_REQUIRED_DEP_PACKAGES})
  SET(LIB_OPTIONAL_DEP_PACKAGES_${POSTFIX} ${LIB_OPTIONAL_DEP_PACKAGES})
  SET(TEST_REQUIRED_DEP_PACKAGES_${POSTFIX} ${TEST_REQUIRED_DEP_PACKAGES})
  SET(TEST_OPTIONAL_DEP_PACKAGES_${POSTFIX} ${TEST_OPTIONAL_DEP_PACKAGES})

  SET(LIB_REQUIRED_DEP_TPLS_${POSTFIX} ${LIB_REQUIRED_DEP_TPLS})
  SET(LIB_OPTIONAL_DEP_TPLS_${POSTFIX} ${LIB_OPTIONAL_DEP_TPLS})
  SET(TEST_REQUIRED_DEP_TPLS_${POSTFIX} ${TEST_REQUIRED_DEP_TPLS})
  SET(TEST_OPTIONAL_DEP_TPLS_${POSTFIX} ${TEST_OPTIONAL_DEP_TPLS})

ENDMACRO()


# @MACRO: TRIBITS_READ_BACK_DEPENDENCIES_VARS()
#
# Usage::
#
#   TRIBITS_READ_BACK_DEPENDENCIES_VARS(<postfix>)
#
# Read back the local package dependency vars from the saved-off vars with
# suffix ``_<postfix>``.
#
# See `Function call tree for constructing package dependency graph`_
#
MACRO(TRIBITS_READ_BACK_DEPENDENCIES_VARS  POSTFIX)

  SET(LIB_REQUIRED_DEP_PACKAGES ${LIB_REQUIRED_DEP_PACKAGES_${POSTFIX}})
  SET(LIB_OPTIONAL_DEP_PACKAGES ${LIB_OPTIONAL_DEP_PACKAGES_${POSTFIX}})
  SET(TEST_REQUIRED_DEP_PACKAGES ${TEST_REQUIRED_DEP_PACKAGES_${POSTFIX}})
  SET(TEST_OPTIONAL_DEP_PACKAGES ${TEST_OPTIONAL_DEP_PACKAGES_${POSTFIX}})

  SET(LIB_REQUIRED_DEP_TPLS ${LIB_REQUIRED_DEP_TPLS_${POSTFIX}})
  SET(LIB_OPTIONAL_DEP_TPLS ${LIB_OPTIONAL_DEP_TPLS_${POSTFIX}})
  SET(TEST_REQUIRED_DEP_TPLS ${TEST_REQUIRED_DEP_TPLS_${POSTFIX}})
  SET(TEST_OPTIONAL_DEP_TPLS ${TEST_OPTIONAL_DEP_TPLS_${POSTFIX}})

ENDMACRO()


# @MACRO: TRIBITS_PROCESS_PACKAGE_DEPENDENCIES_LISTS()
#
# Usage::
#
#   TRIBITS_PROCESS_PACKAGE_DEPENDENCIES_LISTS(<packageName>)
#
# Sets up the upstsream and downstream/forward package dependency list
# varaibles for ``<packageName>`` descrdibed in `List variables defining the
# package dependencies graph`_.  Note that the downstream/forward dependencies
# of upstream packages on this package ``<packageName>`` are built up
# incrimentally.
#
# See `Function call tree for constructing package dependency graph`_
# 
MACRO(TRIBITS_PROCESS_PACKAGE_DEPENDENCIES_LISTS  PACKAGE_NAME)

  TRIBITS_SET_DEP_PACKAGES(${PACKAGE_NAME} LIB  REQUIRED)
  TRIBITS_SET_DEP_PACKAGES(${PACKAGE_NAME} LIB  OPTIONAL)
  TRIBITS_SET_DEP_PACKAGES(${PACKAGE_NAME} TEST  REQUIRED)
  TRIBITS_SET_DEP_PACKAGES(${PACKAGE_NAME} TEST  OPTIONAL)

  SET(${PACKAGE_NAME}_LIB_REQUIRED_DEP_TPLS ${LIB_REQUIRED_DEP_TPLS})
  SET(${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS ${LIB_OPTIONAL_DEP_TPLS})
  SET(${PACKAGE_NAME}_TEST_REQUIRED_DEP_TPLS ${TEST_REQUIRED_DEP_TPLS})
  SET(${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS ${TEST_OPTIONAL_DEP_TPLS})

  TRIBITS_APPEND_FORWARD_DEP_PACKAGES(${PACKAGE_NAME} LIB_REQUIRED_DEP_PACKAGES)
  TRIBITS_APPEND_FORWARD_DEP_PACKAGES(${PACKAGE_NAME} LIB_OPTIONAL_DEP_PACKAGES)
  TRIBITS_APPEND_FORWARD_DEP_PACKAGES(${PACKAGE_NAME} TEST_REQUIRED_DEP_PACKAGES)
  TRIBITS_APPEND_FORWARD_DEP_PACKAGES(${PACKAGE_NAME} TEST_OPTIONAL_DEP_PACKAGES)

ENDMACRO()


# @FUNCTION: TRIBITS_SET_DEP_PACKAGES()
#
# Usage::
#
#   TRIBITS_SET_DEP_PACKAGES(<packageName>  LIB|TEST  REQUIRED|OPTIONAL)
#
# Function that helps to set up backward package dependency lists for a given
# package given the vars read in from the macro
# `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_.
#
# Sets the upstream/backward dependency variables defined in the section `List
# variables defining the package dependencies graph`_.
#
# This also handles the several types of issues:
#
# * A package declaring a dependency on itself
#   (`TRIBITS_ABORT_ON_SELF_DEP()`_).
#
# * A missing upstream dependent package (either error out with
#   `TRIBITS_ABORT_ON_MISSING_PACKAGE()`_ or allow to be missing and disable
#   this package if this is a required dependency).
#
# See `Function call tree for constructing package dependency graph`_
#
FUNCTION(TRIBITS_SET_DEP_PACKAGES  PACKAGE_NAME   LIB_OR_TEST  REQUIRED_OR_OPTIONAL)

  IF (TRIBITS_SET_DEP_PACKAGES_DEBUG_DUMP)
    MESSAGE("\nTRIBITS_SET_DEP_PACKAGES:  ${PACKAGE_NAME}  ${LIB_OR_TEST}  ${REQUIRED_OR_OPTIONAL})")
  ENDIF()

  SET(LIST_TYPE  ${LIB_OR_TEST}_${REQUIRED_OR_OPTIONAL}_DEP_PACKAGES)
  SET(PACKAGE_DEPS_LIST)
  SET(SE_PACKAGE_ENABLE_VAR  ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

  FOREACH(DEP_PKG ${${LIST_TYPE}})
    IF (TRIBITS_SET_DEP_PACKAGES_DEBUG_DUMP)
      PRINT_VAR(DEP_PKG)
    ENDIF()
    IF (${DEP_PKG} STREQUAL ${PACKAGE_NAME})
      TRIBITS_ABORT_ON_SELF_DEP("${PACKAGE_NAME}" "${LIST_TYPE}")
    ENDIF()
    IF (${DEP_PKG}_SOURCE_DIR)
      SET(DEP_PKG_DEFINED_AND_EXISTS TRUE)
    ELSE()
      SET(DEP_PKG_DEFINED_AND_EXISTS FALSE)
    ENDIF()
    IF (TRIBITS_SET_DEP_PACKAGES_DEBUG_DUMP)
      PRINT_VAR(DEP_PKG_DEFINED_AND_EXISTS)
    ENDIF()
    IF (DEP_PKG_DEFINED_AND_EXISTS)
      LIST(APPEND PACKAGE_DEPS_LIST ${DEP_PKG})
    ELSE()
      IF (${PROJECT_NAME}_ASSERT_MISSING_PACKAGES
          AND NOT ${DEP_PKG}_ALLOW_MISSING_EXTERNAL_PACKAGE
        )
        TRIBITS_ABORT_ON_MISSING_PACKAGE(
          "${DEP_PKG}" "${PACKAGE_NAME}" "${PROJECT_NAME}_SE_PACKAGES")
      ELSE()
        IF (${DEP_PKG}_ALLOW_MISSING_EXTERNAL_PACKAGE)
          IF (${PROJECT_NAME}_WARN_ABOUT_MISSING_EXTERNAL_PACKAGES)
            MESSAGE_WRAPPER("NOTE: ${DEP_PKG} is being ignored since its directory"
              " is missing and ${DEP_PKG}_ALLOW_MISSING_EXTERNAL_PACKAGE ="
              " ${${DEP_PKG}_ALLOW_MISSING_EXTERNAL_PACKAGE}!")
          ENDIF()
          IF (REQUIRED_OR_OPTIONAL STREQUAL "REQUIRED")
            MESSAGE_WRAPPER("NOTE: Setting ${SE_PACKAGE_ENABLE_VAR}=OFF because"
              " package ${PACKAGE_NAME} has a required dependency on missing"
              " package ${DEP_PKG}!")
            DUAL_SCOPE_SET(${SE_PACKAGE_ENABLE_VAR} OFF)
          ENDIF()
        ENDIF()
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          MESSAGE(
            "\n***"
            "\n*** NOTE: The package ${DEP_PKG} which is a dependent package of"
              " ${PACKAGE_NAME} being ignored because ${DEP_PKG} is missing!"
            "\n***\n" )
        ENDIF()
        # Must set enable vars for missing package to off so that logic in
        # existing downstream packages that key off of these vars will still
        # work.
        DUAL_SCOPE_SET(${PROJECT_NAME}_ENABLE_${DEP_PKG} OFF)
        DUAL_SCOPE_SET(${PACKAGE_NAME}_ENABLE_${DEP_PKG} OFF)
      ENDIF()
    ENDIF()
  ENDFOREACH()

  #PRINT_VAR(PACKAGE_DEPS_LIST)

  GLOBAL_SET(${PACKAGE_NAME}_${LIST_TYPE} ${PACKAGE_DEPS_LIST})

ENDFUNCTION()


# @FUNCTION: TRIBITS_APPEND_FORWARD_DEP_PACKAGES()
#
# Usage: TRIBITS_APPEND_FORWARD_DEP_PACKAGES(<packageName> <listType>)
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
FUNCTION(TRIBITS_APPEND_FORWARD_DEP_PACKAGES PACKAGE_NAME LIST_TYPE)

  SET(DEP_PKG_LIST_NAME "${PACKAGE_NAME}_${LIST_TYPE}")

  ASSERT_DEFINED(${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)
  FOREACH(DEP_PKG ${${DEP_PKG_LIST_NAME}})
    SET(FWD_DEP_PKG_LIST_NAME "${DEP_PKG}_FORWARD_${LIST_TYPE}")
    IF (NOT DEFINED ${FWD_DEP_PKG_LIST_NAME})
      IF (${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)
        TRIBITS_ABORT_ON_MISSING_PACKAGE(${DEP_PKG} ${PACKAGE_NAME} ${DEP_PKG_LIST_NAME})
      ELSE()
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          MESSAGE(
            "\n***"
            "\n*** NOTE: The package ${DEP_PKG} has forward dependent package"
              " ${PACKAGE_NAME}, but that dependency is being ignored because the package"
              " ${DEP_PKG} is missing!"
            "\n***\n" )
        ENDIF()
      ENDIF()
    ELSE()
      SET(${FWD_DEP_PKG_LIST_NAME} ${${FWD_DEP_PKG_LIST_NAME}} ${PACKAGE_NAME}
        PARENT_SCOPE)
    ENDIF()
  ENDFOREACH()

ENDFUNCTION()


# @MACRO: TRIBITS_SET_PACAKGE_REGRESSION_EMAIL_LIST()
#
# Usage::
#
#  TRIBITS_SET_PACAKGE_REGRESSION_EMAIL_LIST(<packageName>)
#
# Macro that sets a pacakge's regression email address
# ``${PACKAGE_NAME}_REGRESSION_EMAIL_LIST`` as described in ???.
#
# See `Function call tree for constructing package dependency graph`_
#
MACRO(TRIBITS_SET_PACAKGE_REGRESSION_EMAIL_LIST PACKAGE_NAME)

  # Lower-case package name To be used with auto email naming based on base email address
  STRING(TOLOWER "${PACKAGE_NAME}" LPACKAGE)
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(REGRESSION_EMAIL_LIST)
  ENDIF()

  TRIBITS_GET_REPO_NAME(${${PACKAGE_NAME}_PARENT_REPOSITORY} REPOSITORY_NAME)
  #PRINT_VAR(REPOSITORY_NAME)

  IF(${REPOSITORY_NAME}_REPOSITORY_OVERRIDE_PACKAGE_EMAIL_LIST)
    SET(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST
      ${${REPOSITORY_NAME}_REPOSITORY_OVERRIDE_PACKAGE_EMAIL_LIST})
  ELSEIF (REGRESSION_EMAIL_LIST)
    SET(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST ${REGRESSION_EMAIL_LIST})
  ELSEIF (${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESS_BASE)
    SET(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST
      "${LPACKAGE}-regression@${${REPOSITORY_NAME}_REPOSITORY_EMAIL_URL_ADDRESS_BASE}")
  ELSEIF (${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESS)
    SET(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST
      "${${REPOSITORY_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESS}")
  ELSEIF (${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESS_BASE)
    SET(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST
      "${LPACKAGE}-regression@${${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESS_BASE}")
  ELSEIF (${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESS)
    SET(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST
      "${${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESS}")
  ELSE()
    SET(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST "")
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PACKAGE_NAME}_REGRESSION_EMAIL_LIST)
  ENDIF()

ENDMACRO()


# @FUNCTION: TRIBITS_ABORT_ON_MISSING_PACKAGE()
#
# Usage::
#
#   TRIBITS_ABORT_ON_MISSING_PACKAGE(<depPkg>  <packageName> <depPkgListName>)
#
# Function that creates error message about missing/misspelled package.  This
# error message also suggests that the package might be defining an upstream
# dependency on a downstream dependency (i.e. a circular dependency).
#
# See `Function call tree for constructing package dependency graph`_
#
FUNCTION(TRIBITS_ABORT_ON_MISSING_PACKAGE   DEP_PKG  PACKAGE_NAME  DEP_PKG_LIST_NAME)
  MULTILINE_SET(ERRMSG
    "Error, the package '${DEP_PKG}' is listed as a dependency of the package"
    " '${PACKAGE_NAME}' is in the list '${DEP_PKG_LIST_NAME}' but the package"
    " '${DEP_PKG}' is either not defined or is listed later in the package order."
    "  This may also be an attempt to create a circular dependency between"
    " the packages '${DEP_PKG}' and '${PACKAGE_NAME}' (which is not allowed)."
    "  Check the spelling of '${DEP_PKG}' or see how it is listed in"
    " a call to TRIBITS_REPOSITORY_DEFINE_PACKAGES() in relation to"
    " '${PACKAGE_NAME}'.")
  MESSAGE(FATAL_ERROR ${ERRMSG})
ENDFUNCTION()


# @FUNCTION: TRIBITS_ABORT_ON_SELF_DEP()
#
# Usage::
#
#   TRIBITS_ABORT_ON_SELF_DEP(<packageName> <depPkgListName>)
#
# Prints a fatal error message for an attempt for a self dependency
# declaration and which list it comes from.
#
# See `Function call tree for constructing package dependency graph`_
#
FUNCTION(TRIBITS_ABORT_ON_SELF_DEP  PACKAGE_NAME  DEP_PKG_LIST_NAME)
  MULTILINE_SET(ERRMSG
    "Error, the package '${PACKAGE_NAME}' is listed as a dependency of itself"
    " in the list '${DEP_PKG_LIST_NAME}'!")
  MESSAGE(FATAL_ERROR ${ERRMSG})
ENDFUNCTION()


################################################################################
# Macros/functions for processing dependency info for subpackages
################################################################################


# @MACRO: TRIBITS_PARSE_SUBPACKAGES_APPEND_SE_PACKAGES_ADD_OPTIONS()
#
# Usage::
#
#   TRIBITS_PARSE_SUBPACKAGES_APPEND_SE_PACKAGES_ADD_OPTIONS(<toplevelPackageName>)
#
# Macro that parses the read-in variable
# ``SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS`` set by the macro
# `TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()`_ , add subpackages to the list of
# defined packages, and define user cache var options for those subpackages.
#
# This sets the list varaibles for the parent package ``<toplevelPackageName>``::
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
# And it appends for each subpackage to varaible::
#
#   ${PROJECT_NAME}_SE_PACKAGES (old)
#
# See `Function call tree for constructing package dependency graph`_
#
MACRO(TRIBITS_PARSE_SUBPACKAGES_APPEND_SE_PACKAGES_ADD_OPTIONS
  PACKAGE_NAME
  )

  #MESSAGE("TRIBITS_PARSE_SUBPACKAGES_APPEND_SE_PACKAGES_ADD_OPTIONS: ${PACKAGE_NAME}")

  # Structure of SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  SET(SPDC_SP_NAME_OFFSET 0)
  SET(SPDC_SP_DIR_OFFSET 1)
  SET(SPDC_SP_CLASSIFICATION_OFFSET 2)
  SET(SPDC_SP_OPTREQ_OFFSET 3)
  SET(SPDC_NUM_FIELDS 4)

  SET(${PACKAGE_NAME}_SUBPACKAGES "")
  SET(${PACKAGE_NAME}_SUBPACKAGE_DIRS "")
  SET(${PACKAGE_NAME}_SUBPACKAGE_OPTREQ "")

  IF (SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS)

    LIST(LENGTH SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS SPDC_TOTAL_LENGTH)
    MATH(EXPR NUM_SUBPACKAGES "${SPDC_TOTAL_LENGTH}/${SPDC_NUM_FIELDS}")
    MATH(EXPR SUBPACKAGES_LAST_IDX "${NUM_SUBPACKAGES}-1")

    FOREACH(SUBPACKAGE_IDX RANGE ${SUBPACKAGES_LAST_IDX})

      #MESSAGE("")
      #PRINT_VAR(SUBPACKAGE_IDX)

      # SUBPACKAGE_NAME
      MATH(EXPR SUBPACKAGE_NAME_IDX
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_NAME_OFFSET}")
      LIST(GET SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${SUBPACKAGE_NAME_IDX}
        SUBPACKAGE_NAME )

      SET(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${SUBPACKAGE_NAME})

      # SUBPACKAGE_DIR
      MATH(EXPR SUBPACKAGE_DIR_IDX
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_DIR_OFFSET}")
      LIST(GET SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${SUBPACKAGE_DIR_IDX}
        SUBPACKAGE_DIR )

      # SUBPACKAGE_CLASSIFICATION
      MATH(EXPR SUBPACKAGE_CLASSIFICATION_IDX
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_CLASSIFICATION_OFFSET}")
      LIST(GET SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${SUBPACKAGE_CLASSIFICATION_IDX}
        SUBPACKAGE_CLASSIFICATION )

      # ToDo: Parse out TESTGROUP and MATURITYLEVEL (Trilinos #6042)
      SET(SUBPACKAGE_TESTGROUP ${SUBPACKAGE_CLASSIFICATION})

      TRIBITS_UPDATE_PS_PT_SS_ST(Subpackage ${SUBPACKAGE_FULLNAME} SUBPACKAGE_TESTGROUP)

      # SUBPACKAGE_OPTREQ
      MATH(EXPR SUBPACKAGE_OPTREQ_IDX
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_OPTREQ_OFFSET}")
      LIST(GET SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${SUBPACKAGE_OPTREQ_IDX}
        SUBPACKAGE_OPTREQ )

      # Determine if this subpackage exists
      SET(SUBPACKAGE_FULL_SOURCE_DIR
        ${PROJECT_SOURCE_DIR}/${${PACKAGE_NAME}_REL_SOURCE_DIR}/${SUBPACKAGE_DIR})
      IF (EXISTS ${SUBPACKAGE_FULL_SOURCE_DIR})
         SET(SUBPACKAGE_EXISTS TRUE)
      ELSE()
         SET(SUBPACKAGE_EXISTS FALSE)
      ENDIF()

      IF (NOT SUBPACKAGE_EXISTS AND ${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)
         MESSAGE(SEND_ERROR "ERROR: Subpackage dir '${SUBPACKAGE_FULL_SOURCE_DIR}'"
           " is missing!")
      ENDIF()

      # Append to lists and global variables

      IF (SUBPACKAGE_EXISTS)

        LIST(APPEND ${PACKAGE_NAME}_SUBPACKAGES ${SUBPACKAGE_NAME})
        LIST(APPEND ${PACKAGE_NAME}_SUBPACKAGE_DIRS ${SUBPACKAGE_DIR})
        LIST(APPEND ${PACKAGE_NAME}_SUBPACKAGE_OPTREQ ${SUBPACKAGE_OPTREQ})
        LIST(APPEND ${PROJECT_NAME}_SE_PACKAGES ${SUBPACKAGE_FULLNAME})
        SET(${SUBPACKAGE_FULLNAME}_SOURCE_DIR "${SUBPACKAGE_FULL_SOURCE_DIR}")
        SET(${SUBPACKAGE_FULLNAME}_REL_SOURCE_DIR
          "${${PACKAGE_NAME}_REL_SOURCE_DIR}/${SUBPACKAGE_DIR}")
        SET(${SUBPACKAGE_FULLNAME}_PARENT_PACKAGE ${PACKAGE_NAME})
        SET(${SUBPACKAGE_FULLNAME}_PARENT_REPOSITORY
          ${${PACKAGE_NAME}_PARENT_REPOSITORY})

        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          PRINT_VAR(${SUBPACKAGE_FULLNAME}_PARENT_PACKAGE)
          PRINT_VAR(${SUBPACKAGE_FULLNAME}_PARENT_REPOSITORY)
        ENDIF()

        # Set up the input options for this subpackage
        TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS(${SUBPACKAGE_FULLNAME}
          ${SUBPACKAGE_TESTGROUP})

      ENDIF()

    ENDFOREACH()

  ENDIF()

  #PRINT_VAR(${PACKAGE_NAME}_SUBPACKAGES)
  #PRINT_VAR(${PACKAGE_NAME}_SUBPACKAGE_OPTREQ)

ENDMACRO()


# @MACRO: TRIBITS_READ_PACKAGE_SUBPACKAGE_DEPS_FILES_ADD_TO_GRAPH()
#
# Usage::
#
#   TRIBITS_READ_PACKAGE_SUBPACKAGE_DEPS_FILES_ADD_TO_GRAPH(<toplevelPackageName>)
#
# Read in subpackages dependencies files and add to dependencies graph
# variables.
#
# See `Function call tree for constructing package dependency graph`_
#
MACRO(TRIBITS_READ_PACKAGE_SUBPACKAGE_DEPS_FILES_ADD_TO_GRAPH  PACKAGE_NAME)

  #MESSAGE("TRIBITS_READ_PACKAGE_SUBPACKAGE_DEPS_FILES_ADD_TO_GRAPH: ${PACKAGE_NAME}")

  #PRINT_VAR(${PROJECT_NAME}_SE_PACKAGES)

  SET(SUBPACKAGE_IDX 0)
  FOREACH(TRIBITS_SUBPACKAGE ${${PACKAGE_NAME}_SUBPACKAGES})
    LIST(GET ${PACKAGE_NAME}_SUBPACKAGE_DIRS ${SUBPACKAGE_IDX} SUBPACKAGE_DIR)
    TRIBITS_READ_SUBPACKAGE_DEPS_FILE_ADD_TO_GRAPH(${TRIBITS_PACKAGE}
      ${TRIBITS_SUBPACKAGE}  ${SUBPACKAGE_DIR})
    MATH(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")
  ENDFOREACH()

ENDMACRO()


# @MACRO: TRIBITS_READ_SUBPACKAGE_DEPS_FILE_ADD_TO_GRAPH()
#
# Usage::
#
#   TRIBITS_READ_SUBPACKAGE_DEPS_FILE_ADD_TO_GRAPH(<toplevelPackageName>
#     <subpackageName> <subpackageDir>)
#
# Macro that reads in a single subpackage dependencies file
# `<packageDir>/<spkgDir>/cmake/Dependencies.cmake`_ and sets up the
# dependency structure for it.
#
# See `Function call tree for constructing package dependency graph`_
#
MACRO(TRIBITS_READ_SUBPACKAGE_DEPS_FILE_ADD_TO_GRAPH  PACKAGE_NAME
  SUBPACKAGE_NAME  SUBPACKAGE_DIR
  )

  #MESSAGE("TRIBITS_READ_SUBPACKAGE_DEPS_FILE_ADD_TO_GRAPH: ${PACKAGE_NAME} ${SUBPACKAGE_NAME} ${SUBPACKAGE_DIR}")

  SET(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${SUBPACKAGE_NAME})

  #
  # A) Get ready to read in the contents of this this subpakages's
  # Dependencies.cmake file
  #

  TRIBITS_PREP_TO_READ_DEPENDENCIES(${SUBPACKAGE_FULLNAME})

  # NOTE: Subpackages use the regression email list from the parent package.

  # NOTE: Subpackages are not allowed to have subpackages!
  SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS)

  #
  # B) Read in this subpackage's Dependencies file
  #

  SET(SUBPACKAGE_FULL_DIR "${${PACKAGE_NAME}_REL_SOURCE_DIR}/${SUBPACKAGE_DIR}")

  SET(SUBPACKAGE_ABS_DIR "${PROJECT_SOURCE_DIR}/${SUBPACKAGE_FULL_DIR}")
  SET(SUBPAKCAGE_DEPENDENCIES_FILE "${SUBPACKAGE_ABS_DIR}/cmake/Dependencies.cmake")

  IF (EXISTS ${SUBPAKCAGE_DEPENDENCIES_FILE})
    SET(SUBPACKAGE_EXISTS TRUE)
  ELSE()
    SET(SUBPACKAGE_EXISTS FALSE)
  ENDIF()

  IF (SUBPACKAGE_EXISTS OR ${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)

    TRIBITS_TRACE_FILE_PROCESSING(PACKAGE  INCLUDE  "${SUBPAKCAGE_DEPENDENCIES_FILE}")
    INCLUDE(${SUBPAKCAGE_DEPENDENCIES_FILE})

    TRIBITS_ASSERT_READ_DEPENDENCY_VARS(${SUBPACKAGE_FULLNAME})

    TRIBITS_PROCESS_PACKAGE_DEPENDENCIES_LISTS(${SUBPACKAGE_FULLNAME})

    SET(${SUBPACKAGE_FULLNAME}_REGRESSION_EMAIL_LIST
      ${${PACKAGE_NAME}_REGRESSION_EMAIL_LIST})

  ENDIF()

ENDMACRO()


#
# Private utility functions
#


# Function that sets a varaible to DECLARED-UNDEFINED
#
FUNCTION(TRIBITS_DECLARE_UNDEFINED  VAR_NAME)
  SET(${VAR_NAME}  DECLARED-UNDEFINED  PARENT_SCOPE)
ENDFUNCTION()


# Function that asserts that a package dependency variable is defined
# correctly
#
FUNCTION(TRIBITS_ASSERT_DEFINED_PACKAGE_VAR  PACKAGE_VAR  PACKAGE_NAME)
  IF (${PACKAGE_VAR} STREQUAL DECLARED-UNDEFINED)
    MESSAGE(FATAL_ERROR
      "Error, the package variable ${PACKAGE_VAR} was not defined correctly for package ${PACKAGE_NAME}!"
      )
  ENDIF()
ENDFUNCTION()
