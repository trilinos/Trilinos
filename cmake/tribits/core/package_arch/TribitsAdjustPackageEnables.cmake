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


INCLUDE(TribitsProcessPackagesAndDirsLists)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsGeneralMacros)

INCLUDE(AdvancedOption)
INCLUDE(AdvancedSet)
INCLUDE(AppendStringVar)
INCLUDE(CMakeBuildTypesList)
INCLUDE(FindListElement)
INCLUDE(GlobalNullSet)
INCLUDE(PrintNonemptyVar)
INCLUDE(PrintNonemptyVarWithSpaces)
INCLUDE(PrintVar)
INCLUDE(RemoveGlobalDuplicates)
INCLUDE(SetDefault)
INCLUDE(MessageWrapper)
INCLUDE(DualScopeSet)
INCLUDE(CMakeParseArguments)


#
# @FUNCTION: TRIBITS_SET_ST_FOR_DEV_MODE()
#
# Function that allows packages to easily make a feature ``ST`` for
# development builds and ``PT`` for release builds by default.
#
# Usage::
#
#   TRIBITS_SET_ST_FOR_DEV_MODE(<outputVar>)
#
# This function is typically called in a package's top-level
# `<packageDir>/CMakeLists.txt`_ file before defining other options for the
# package.  The output variable ``${<outputVar>}`` is set to ``ON`` or ``OFF``
# based on the configure state.  In development mode
# (i.e. ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE==ON``), ``${<outputVar>}``
# will be set to ``ON`` only if ``ST`` code is enabled
# (i.e. ``${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE==ON``), otherwise it is
# set to ``OFF``. In release mode
# (i.e. ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE==OFF``), ``${<outputVar>}``
# is always set to ``ON``.  This allows some parts of a TriBITS package to be
# considered ``ST`` for development mode (thereby reducing testing time by not
# enabling the dependent features/tests), while still having important
# functionality available to users by default in a release of the package.
#
FUNCTION(TRIBITS_SET_ST_FOR_DEV_MODE  OUTPUT_VAR)
  IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
    SET(OUTPUT_VAL ${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE})
  ELSE()
    SET(OUTPUT_VAL ON)
  ENDIF()
  SET(${OUTPUT_VAR} ${OUTPUT_VAL} PARENT_SCOPE)
ENDFUNCTION()


# For backward compatibility
MACRO(TRIBITS_SET_SS_FOR_DEV_MODE  OUTPUT_VAR)
  MESSAGE(WARNING
    "WARNING: TRIBITS_SET_SS_FOR_DEV_MODE() is deprecated,"
    " use TRIBITS_SET_ST_FOR_DEV_MODE() instead!")
  TRIBITS_SET_ST_FOR_DEV_MODE(${OUTPUT_VAR})
ENDMACRO()


#
# Set the combined directory name taking into account '.' repos.
#
FUNCTION(TRIBITS_GET_REPO_NAME  REPO_DIR  REPO_NAME_OUT)
  IF (REPO_DIR STREQUAL ".")
    SET(REPO_NAME ${PROJECT_NAME})
  ELSE()
    SET(REPO_NAME ${REPO_DIR})
  ENDIF()
  SET(${REPO_NAME_OUT} "${REPO_NAME}" PARENT_SCOPE)
ENDFUNCTION()

#
# Set the combined directory name taking into account '.' repos.
#
FUNCTION(TRIBITS_SET_BASE_REPO_DIR  BASE_DIR  REPO_DIR  BASE_REPO_DIR_OUT)
  IF (REPO_DIR STREQUAL ".")
    SET(REPO_DIR_STR "")
  ELSE()
    SET(REPO_DIR_STR "/${REPO_DIR}")
  ENDIF()
  SET(${BASE_REPO_DIR_OUT} "${BASE_DIR}${REPO_DIR_STR}" PARENT_SCOPE)
ENDFUNCTION()


#
# Function that creates error message about missing/misspelled package.
#

FUNCTION(TRIBITS_ABORT_ON_MISSING_PACKAGE   DEP_PKG  PACKAGE_NAME  DEP_PKG_LIST_NAME)
  MULTILINE_SET(ERRMSG
    "Error, the package '${DEP_PKG}' is listed as a dependency of the package"
    " '${PACKAGE_NAME}' is in the list '${DEP_PKG_LIST_NAME}' but the package"
    " '${DEP_PKG}' is either not defined or is listed later in the package order."
    "  This may also be an attempt to create a cicular dependency between"
    " the packages '${DEP_PKG}' and '${PACKAGE_NAME}' (which is not allowed)."
    "  Check the spelling of '${DEP_PKG}' or see how it is listed in"
    " ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS in relationship to"
    " '${PACKAGE_NAME}'.")
  MESSAGE(FATAL_ERROR ${ERRMSG})
ENDFUNCTION()


FUNCTION(TRIBITS_ABORT_ON_SELF_DEP  PACKAGE_NAME  DEP_PKG_LIST_NAME)
  MULTILINE_SET(ERRMSG
    "Error, the package '${PACKAGE_NAME}' is listed as a dependency itself"
    " in the list '${DEP_PKG_LIST_NAME}'!")
  MESSAGE(FATAL_ERROR ${ERRMSG})
ENDFUNCTION()


#
# Function that helps to set up backward package dependency lists
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


#
# Macro that helps to set up forward package dependency lists
#

FUNCTION(TRIBITS_APPEND_FORWARD_DEP_PACKAGES PACKAGE_NAME LIST_TYPE)

  #MESSAGE("\nPACKAGE_ARCH_APPEND_FORWARD_DEP_PACKAGES: ${PACKAGE_NAME} ${LIST_TYPE}")

  SET(DEP_PKG_LIST_NAME "${PACKAGE_NAME}_${LIST_TYPE}")

  #MESSAGE("DEP_PKG_LIST_NAME = ${DEP_PKG_LIST_NAME}")
  #MESSAGE("${DEP_PKG_LIST_NAME} = ${${DEP_PKG_LIST_NAME}}")

  ASSERT_DEFINED(${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)
  FOREACH(DEP_PKG ${${DEP_PKG_LIST_NAME}})
    #MESSAGE("DEP_PKG = ${DEP_PKG}")
    SET(FWD_DEP_PKG_LIST_NAME "${DEP_PKG}_FORWARD_${LIST_TYPE}")
    #MESSAGE("FWD_DEP_PKG_LIST_NAME = ${FWD_DEP_PKG_LIST_NAME}")
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
      SET(${FWD_DEP_PKG_LIST_NAME} ${${FWD_DEP_PKG_LIST_NAME}} ${PACKAGE_NAME} PARENT_SCOPE)
    ENDIF()
  ENDFOREACH()

ENDFUNCTION()


################################################################################
#
# Helper macros for TRIBITS_READ_PACKAGE_DEPENDENCIES()


MACRO(TRIBITS_PREP_TO_READ_DEPENDENCIES)

  DECLARE_UNDEFINED(LIB_REQUIRED_DEP_PACKAGES)
  DECLARE_UNDEFINED(LIB_OPTIONAL_DEP_PACKAGES)
  DECLARE_UNDEFINED(TEST_REQUIRED_DEP_PACKAGES)
  DECLARE_UNDEFINED(TEST_OPTIONAL_DEP_PACKAGES)

  DECLARE_UNDEFINED(LIB_REQUIRED_DEP_TPLS "")
  DECLARE_UNDEFINED(LIB_OPTIONAL_DEP_TPLS "")
  DECLARE_UNDEFINED(TEST_REQUIRED_DEP_TPLS "")
  DECLARE_UNDEFINED(TEST_OPTIONAL_DEP_TPLS "")

ENDMACRO()


#
# @MACRO: TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()
#
# Define the dependencies for a given `TriBITS SE Package`_ (i.e. a top-level
# `TriBITS Package`_ or a `TriBITS Subpackage`_) in the package's
# `<packageDir>/cmake/Dependencies.cmake`_ file.
#
# Usage::
#
#   TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
#      [LIB_REQUIRED_PACKAGES <pkg1> <pkg2> ...]
#      [LIB_OPTIONAL_PACKAGES <pkg1> <pkg2> ...]
#      [TEST_REQUIRED_PACKAGES <pkg1> <pkg2> ...]
#      [TEST_OPTIONAL_PACKAGES <pkg1> <pkg2> ...]
#      [LIB_REQUIRED_TPLS <tpl1> <tpl2> ...]
#      [LIB_OPTIONAL_TPLS <tpl1> <tpl2> ...]
#      [TEST_REQUIRED_TPLS <tpl1> <tpl2> ...]
#      [TEST_OPTIONAL_TPLS <tpl1> <tpl2> ...]
#      [REGRESSION_EMAIL_LIST  <regression-email-address>]
#      [SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
#        <spkg1_name>  <spkg1_dir>  <spkg1_classifications>  <spkg1_optreq>
#        <spkg2_name>  <spkg2_dir>  <spkg2_classifications>  <spkg2_optreq>
#        ...
#        ]
#      )
#
# Every argument in this macro is optional (that is, an SE package can have no
# upstream dependencies).  The arguments that apply to all SE packages are:
#
#   ``LIB_REQUIRED_PACKAGES``
#
#     List of required upstream SE packages that must be enabled in order to
#     build and use the libraries (or capabilities) in this SE package.
#
#   ``LIB_OPTIONAL_PACKAGES``
#
#     List of additional optional upstream SE packages that can be used in
#     this SE package if enabled.  These upstream SE packages need not be
#     enabled in order to enable this SE package but not enabling one or more
#     of these optional upstream SE packages will result in diminished
#     capabilities of this SE package.
#
#   ``TEST_REQUIRED_PACKAGES``
#
#     List of additional upstream SE packages that must be enabled in order to
#     build and/or run the tests and/or examples in this SE package.  If any
#     of these upstream SE packages are not enabled, then there will be no
#     tests or examples defined or run for this SE package.
#
#   ``TEST_OPTIONAL_PACKAGES``
#
#     List of additional optional upstream SE packages that can be used by the
#     tests in this SE package.  These upstream SE packages need not be
#     enabled in order to run some basic tests or examples for this SE
#     package.  Typically, extra tests that depend on optional test SE
#     packages involve integration testing of some type.
#
#   ``LIB_REQUIRED_TPLS``
#
#     List of required upstream TPLs that must be enabled in order to build
#     and use the libraries (or capabilities) in this SE package.
#
#   ``LIB_OPTIONAL_TPLS``
#
#     List of additional optional upstream TPLs that can be used in this SE
#     package if enabled.  These upstream TPLs need not be enabled in order to
#     use this SE package but not enabling one or more of these optional
#     upstream TPLs will result in diminished capabilities of this SE package.
#
#   ``TEST_REQUIRED_TPLS``
#
#     List of additional upstream TPLs that must be enabled in order to build
#     and/or run the tests and/or examples in this SE package.  If any of
#     these upstream TPLs are not enabled, then there will be no tests or
#     examples defined or run for this SE package.
#
#   ``TEST_OPTIONAL_TPLS``
#
#     List of additional optional upstream TPLs that can be used by the tests
#     in this SE package.  These upstream TPLs need not be enabled in order to
#     run basic tests for this SE package.  Typically, extra tests that depend
#     on optional TPLs involve integration testing or some additional testing
#     of some type.
#
# Only upstream SE packages can be listed (as defined by the order the SE
# packages are listed in `TRIBITS_REPOSITORY_DEFINE_PACKAGES()`_ in the
# `<repoDir>/PackagesList.cmake`_ file).  Otherwise an error will occur and
# processing will stop.  Misspelled SE package names are caught as well.
#
# Only direct SE package dependencies need to be listed.  Indirect SE package
# dependencies are automatically handled.  For example, if this SE package
# directly depends on SE package ``PKG2`` which depends on SE package ``PKG1``
# (but this SE package does not directly depend on anything in ``PKG1``) then
# this SE package only needs to list a dependency on ``PKG2``, not ``PKG1``.
# The dependency on ``PKG1`` will be taken care of automatically by the
# TriBITS dependency management system.
#
# However, currently, all TPL dependencies must be listed, even the indirect
# ones.  This is a requirement that will be dropped in a future version of
# TriBITS.
#
# The SE packages listed in ``LIB_REQUIRED_PACKAGES`` are implicitly also
# dependencies in ``TEST_REQUIRED_PACKAGES``.  Likewise
# ``LIB_OPTIONAL_PACKAGES`` are implicitly also dependencies in
# ``TEST_OPTIONAL_PACKAGES``.  Same goes for TPL dependencies.
#
# The upstream dependencies within a single list do not need to be listed in
# any order.  For example if ``PKG2`` depends on ``PKG1``, and this given SE
# package depends on both, then one can list::
#
#   LIB_REQUIRED_PACKAGES PKG2 PKG1
#
# or::
#
#   "LIB_REQUIRED_PACKAGES PKG1 PKG2
#
# Likewise the order that dependent TPLs are listed is not significant.
#
# If some upstream SE packages are allowed to be missing, this can be specified
# by calling the macro `TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES()`_.
#
# A top-level `TriBITS Package`_ can also be broken down into `TriBITS
# Subpackages`_.  In this case, the following argument must be passed in:
#
#   .. _SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS:
#
#   ``SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS``
#
#     2D array with rows listing the subpackages where each row has the
#     columns:
#
#     * **SUBPACKAGE** (Column 0): The name of the subpackage ``<spkg_name>``.
#       The full SE package name is ``${PARENT_PACKAGE_NAME}<spkg_name>``.
#       The full SE package name is what is used in listing dependencies in
#       other SE packages.
#
#     * **DIRS** (Column 1): The subdirectory ``<spkg_dir>`` relative to the
#       parent package's base directory.  All of the contents of the
#       subpackage should be under this subdirectory.  This is assumed by the
#       TriBITS testing support software when mapping modified files to SE
#       packages that need to be tested (see `checkin-test.py`_).
#
#     * **CLASSIFICATIONS** (Column 2): The `Test Test Category`_ `PT`_,
#       `ST`_, `EX`_ and the maturity level ``EP``, ``RS``, ``PG``, ``PM``,
#       ``GRS``, ``GPG``, ``GPM``, and ``UM``, separated by a coma ',' with no
#       spaces in between (e.g. ``"PT,GPM"``).  These have exactly the same
#       meaning as for full packages (see
#       `TRIBITS_REPOSITORY_DEFINE_PACKAGES()`_).
#
#     * **OPTREQ** (Column 3): Determines if the outer parent package has an
#       ``OPTIONAL`` or ``REQUIRED`` dependence on this subpackage.
#
# Other variables that this macro handles:
#
#   ``REGRESSION_EMAIL_LIST``
#
#     The email list that is used to send CDash error messages.  If this
#     argument is missing, then the email list that CDash errors go to is
#     determined by other means (see `CDash regression email addresses`_).
#
# NOTE: All this macro really does is to just define the variables:
#
# * ``LIB_REQUIRED_DEP_PACKAGES``
# * ``LIB_OPTIONAL_DEP_PACKAGES``
# * ``TEST_REQUIRED_DEP_PACKAGES``
# * ``TEST_OPTIONAL_DEP_PACKAGES``
# * ``LIB_REQUIRED_DEP_TPLS``
# * ``LIB_OPTIONAL_DEP_TPLS``
# * ``TEST_REQUIRED_DEP_TPLS``
# * ``TEST_OPTIONAL_DEP_TPLS``
# * ``REGRESSION_EMAIL_LIST``
# * ``SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS``
#
# which are then read by the TriBITS cmake code to build the SE package
# dependency graph.  The advantage of using this macro instead of just
# directly setting the variables is that an SE package only needs to list
# dependencies that exist.  Otherwise, the ``Dependencies.cmake`` file will
# need to set all of the above local variables, even those that are empty.
# This is an error checking property of the TriBITS system to avoid misspelling
# the names of these variables.
#
MACRO(TRIBITS_PACKAGE_DEFINE_DEPENDENCIES)

  CMAKE_PARSE_ARGUMENTS(
     #prefix
     PARSE
     #options
     ""
     #one_value_keywords
     ""
     #multi_value_keywords
     "LIB_REQUIRED_PACKAGES;LIB_OPTIONAL_PACKAGES;TEST_REQUIRED_PACKAGES;TEST_OPTIONAL_PACKAGES;LIB_REQUIRED_TPLS;LIB_OPTIONAL_TPLS;TEST_REQUIRED_TPLS;TEST_OPTIONAL_TPLS;REGRESSION_EMAIL_LIST;SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS"
     ${ARGN}
     )

  TRIBITS_CHECK_FOR_UNPARSED_ARGUMENTS()

  SET(LIB_REQUIRED_DEP_PACKAGES ${PARSE_LIB_REQUIRED_PACKAGES})
  SET(LIB_OPTIONAL_DEP_PACKAGES ${PARSE_LIB_OPTIONAL_PACKAGES})
  SET(TEST_REQUIRED_DEP_PACKAGES ${PARSE_TEST_REQUIRED_PACKAGES})
  SET(TEST_OPTIONAL_DEP_PACKAGES ${PARSE_TEST_OPTIONAL_PACKAGES})
  SET(LIB_REQUIRED_DEP_TPLS ${PARSE_LIB_REQUIRED_TPLS})
  SET(LIB_OPTIONAL_DEP_TPLS ${PARSE_LIB_OPTIONAL_TPLS})
  SET(TEST_REQUIRED_DEP_TPLS ${PARSE_TEST_REQUIRED_TPLS})
  SET(TEST_OPTIONAL_DEP_TPLS ${PARSE_TEST_OPTIONAL_TPLS})
  SET(REGRESSION_EMAIL_LIST ${PARSE_REGRESSION_EMAIL_LIST})
  SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    ${PARSE_SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS})

  # ToDo:
  # * Assert that REGRESSION_EMAIL_LIST has only one entry
  # * Assert that SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS is divisible
  #   by the number of columns!

ENDMACRO()


MACRO(TRIBITS_SAVE_OFF_DEPENENCIES_VARS  POSTFIX)

  SET(LIB_REQUIRED_DEP_PACKAGES_${POSTFIX} ${LIB_REQUIRED_DEP_PACKAGES})
  SET(LIB_OPTIONAL_DEP_PACKAGES_${POSTFIX} ${LIB_OPTIONAL_DEP_PACKAGES})
  SET(TEST_REQUIRED_DEP_PACKAGES_${POSTFIX} ${TEST_REQUIRED_DEP_PACKAGES})
  SET(TEST_OPTIONAL_DEP_PACKAGES_${POSTFIX} ${TEST_OPTIONAL_DEP_PACKAGES})

  SET(LIB_REQUIRED_DEP_TPLS_${POSTFIX} ${LIB_REQUIRED_DEP_TPLS})
  SET(LIB_OPTIONAL_DEP_TPLS_${POSTFIX} ${LIB_OPTIONAL_DEP_TPLS})
  SET(TEST_REQUIRED_DEP_TPLS_${POSTFIX} ${TEST_REQUIRED_DEP_TPLS})
  SET(TEST_OPTIONAL_DEP_TPLS_${POSTFIX} ${TEST_OPTIONAL_DEP_TPLS})

ENDMACRO()


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


#
# Parse the read-in varaible SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS, add
# read subpackages to list of defined SE packages, and define options.
#
# NOTE: Directly reads varaibles ${PACKAGE_NAME} and
# ${SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS} defined in
# TRIBITS_READ_PACKAGE_DEPENDENCIES
#
MACRO(TRIBITS_PARSE_SUBPACKAGES_AND_APPEND_SE_PACKAGES_AND_ADD_OPTIONS  PACKAGE_NAME
    PACKAGE_DIR)

  #MESSAGE("TRIBITS_PARSE_SUBPACKAGES_AND_APPEND_SE_PACKAGES_AND_ADD_OPTIONS: ${PACKAGE_NAME}")

  # Structure of SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  SET(SPDC_SP_NAME_OFFSET 0)
  SET(SPDC_SP_DIR_OFFSET 1)
  SET(SPDC_SP_CLASSIFICATION_OFFSET 2)
  SET(SPDC_SP_OPTREQ_OFFSET 3)
  SET(SPDC_NUM_FIELDS 4)

  SET(${PACKAGE_NAME}_SUBPACKAGES)
  SET(${PACKAGE_NAME}_SUBPACKAGE_DIRS)
  SET(${PACKAGE_NAME}_SUBPACKAGE_OPTREQ)

  IF (SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS)

    LIST(LENGTH SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS SPDC_TOTAL_LENGTH)
    MATH(EXPR NUM_SUBPACKAGES "${SPDC_TOTAL_LENGTH}/${SPDC_NUM_FIELDS}")
    #PRINT_VAR(NUM_SUBPACKAGES)
    MATH(EXPR SUBPACKAGES_LAST_IDX "${NUM_SUBPACKAGES}-1")
    #PRINT_VAR(SUBPACKAGES_LAST_IDX)

    FOREACH(SUBPACKAGE_IDX RANGE ${SUBPACKAGES_LAST_IDX})

      #MESSAGE("")
      #PRINT_VAR(SUBPACKAGE_IDX)

      # SUBPACKAGE_NAME
      MATH(EXPR SUBPACKAGE_NAME_IDX "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_NAME_OFFSET}")
      #PRINT_VAR(SUBPACKAGE_NAME_IDX)
      LIST(GET SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${SUBPACKAGE_NAME_IDX} SUBPACKAGE_NAME )
      #PRINT_VAR(SUBPACKAGE_NAME)

      SET(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${SUBPACKAGE_NAME})

      # SUBPACKAGE_DIR
      MATH(EXPR SUBPACKAGE_DIR_IDX "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_DIR_OFFSET}")
      #PRINT_VAR(SUBPACKAGE_DIR_IDX)
      LIST(GET SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${SUBPACKAGE_DIR_IDX} SUBPACKAGE_DIR )
      #PRINT_VAR(SUBPACKAGE_DIR)

      # SUBPACKAGE_CLASSIFICATION
      MATH(EXPR SUBPACKAGE_CLASSIFICATION_IDX
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_CLASSIFICATION_OFFSET}")
      #PRINT_VAR(SUBPACKAGE_CLASSIFICATION_IDX)
      LIST(GET SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${SUBPACKAGE_CLASSIFICATION_IDX}
        SUBPACKAGE_CLASSIFICATION )
      #PRINT_VAR(SUBPACKAGE_CLASSIFICATION)

      # ToDo: Parse out TESTGROUP and MATURITYLEVEL (Trilinos #6042)
      SET(SUBPACKAGE_TESTGROUP ${SUBPACKAGE_CLASSIFICATION})

      TRIBITS_UPDATE_PS_PT_SS_ST(Subpackage ${SUBPACKAGE_FULLNAME} SUBPACKAGE_TESTGROUP)

      # SUBPACKAGE_OPTREQ
      MATH(EXPR SUBPACKAGE_OPTREQ_IDX
        "${SUBPACKAGE_IDX}*${SPDC_NUM_FIELDS}+${SPDC_SP_OPTREQ_OFFSET}")
      #PRINT_VAR(SUBPACKAGE_OPTREQ_IDX)
      LIST(GET SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS ${SUBPACKAGE_OPTREQ_IDX}
        SUBPACKAGE_OPTREQ )
      #PRINT_VAR(SUBPACKAGE_OPTREQ)

      # Determine if this subpackage exists
      SET(SUBPACKAGE_FULL_SOURCE_DIR ${PROJECT_SOURCE_DIR}/${PACKAGE_DIR}/${SUBPACKAGE_DIR})
      IF (EXISTS ${SUBPACKAGE_FULL_SOURCE_DIR})
         SET(SUBPACKAGE_EXISTS TRUE)
      ELSE()
         SET(SUBPACKAGE_EXISTS FALSE)
      ENDIF()
      #PRINT_VAR(SUBPACKAGE_FULL_SOURCE_DIR)
      #PRINT_VAR(SUBPACKAGE_EXISTS)

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
        SET(${SUBPACKAGE_FULLNAME}_PARENT_PACKAGE ${PACKAGE_NAME})
        SET(${SUBPACKAGE_FULLNAME}_PARENT_REPOSITORY ${${PACKAGE_NAME}_PARENT_REPOSITORY})

        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          PRINT_VAR(${SUBPACKAGE_FULLNAME}_PARENT_PACKAGE)
          PRINT_VAR(${SUBPACKAGE_FULLNAME}_PARENT_REPOSITORY)
        ENDIF()

        # Set up the input options for this subpackage
        TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS(${SUBPACKAGE_FULLNAME}
          ${SUBPACKAGE_TESTGROUP})

        #PRINT_VAR(${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})

      ENDIF()

    ENDFOREACH()

  ENDIF()

  #PRINT_VAR(${PACKAGE_NAME}_SUBPACKAGES)
  #PRINT_VAR(${PACKAGE_NAME}_SUBPACKAGE_OPTREQ)

ENDMACRO()


#
# Macro that reads in a single subpackage dependencies file and sets up
# the dependency structure for it.
#

MACRO(TRIBITS_READ_SUBPACKAGE_DEPENDENCIES  PACKAGE_NAME  PACKAGE_DIR
  SUBPACKAGE_NAME  SUBPACKAGE_DIR)

  #MESSAGE("TRIBITS_READ_SUBPACKAGE_DEPENDENCIES: ${PACKAGE_NAME} ${PACKAGE_DIR} ${SUBPACKAGE_NAME} ${SUBPACKAGE_DIR}")

  SET(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${SUBPACKAGE_NAME})

  #
  # A) Get ready to read in the contents of this this subpakages's Dependencies.cmake file
  #

  SET(${SUBPACKAGE_FULLNAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES "")
  SET(${SUBPACKAGE_FULLNAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES "")
  SET(${SUBPACKAGE_FULLNAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES "")
  SET(${SUBPACKAGE_FULLNAME}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES "")

  TRIBITS_PREP_TO_READ_DEPENDENCIES()

  # NOTE: Subpackages use the regression email list from the parent package.

  # NOTE: Subpackages are not allowed to have subpackages!
  SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS)

  #
  # B) Read in this subpackage's Dependencies file
  #

  SET(SUBPACKAGE_FULL_DIR "${PACKAGE_DIR}/${SUBPACKAGE_DIR}")
  LIST(APPEND ${PROJECT_NAME}_SE_PACKAGE_DIRS ${SUBPACKAGE_FULL_DIR})

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

    #
    # C) Finish processing this subpackage's dependencies into dependency graph vars
    #

    TRIBITS_PROCESS_PACKAGE_DEPENDENCIES_LISTS(${SUBPACKAGE_FULLNAME})

    #PRINT_VAR(${SUBPACKAGE_FULLNAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES)
    #PRINT_VAR(${SUBPACKAGE_FULLNAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES)
    #PRINT_VAR(${SUBPACKAGE_FULLNAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES)
    #PRINT_VAR(${SUBPACKAGE_FULLNAME}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES)

    #
    # D) Set the email addresses for the subpackage to the parent package's
    #

    SET(${SUBPACKAGE_FULLNAME}_REGRESSION_EMAIL_LIST ${${PACKAGE_NAME}_REGRESSION_EMAIL_LIST})

  ENDIF()

ENDMACRO()


#
# Read in subpackages dependencies files and add to dependencies graph variables
#
MACRO(TRIBITS_READ_ALL_PACKAGE_SUBPACKAGE_DEPENDENCIES  PACKAGE_NAME  PACKAGE_DIR)

  #MESSAGE("TRIBITS_READ_ALL_PACKAGE_SUBPACKAGE_DEPENDENCIES: ${PACKAGE_NAME} ${PACKAGE_DIR}")

  #PRINT_VAR(${PROJECT_NAME}_SE_PACKAGES)

  SET(SUBPACKAGE_IDX 0)
  FOREACH(TRIBITS_SUBPACKAGE ${${PACKAGE_NAME}_SUBPACKAGES})
    LIST(GET ${PACKAGE_NAME}_SUBPACKAGE_DIRS ${SUBPACKAGE_IDX} SUBPACKAGE_DIR)
    TRIBITS_READ_SUBPACKAGE_DEPENDENCIES(${TRIBITS_PACKAGE}  ${PACKAGE_DIR}
      ${TRIBITS_SUBPACKAGE}  ${SUBPACKAGE_DIR})
    MATH(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")
  ENDFOREACH()

  LIST(APPEND ${PROJECT_NAME}_SE_PACKAGE_DIRS ${PACKAGE_DIR})

ENDMACRO()


#
# Macro that reads in package dependencies for a package and sets forward
# dependencies for packages already read in.
#
# Modifies the global variables:
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
MACRO(TRIBITS_READ_PACKAGE_DEPENDENCIES  PACKAGE_NAME  PACKAGE_DIR)

  #
  # A) Get ready to read in the contents of this this pakages's Dependencies.cmake file
  #

  SET(${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES "")
  SET(${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES "")
  SET(${PACKAGE_NAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES "")
  SET(${PACKAGE_NAME}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES "")

  TRIBITS_PREP_TO_READ_DEPENDENCIES()

  # Set one regression email list for the package and all subpackages!
  SET(REGRESSION_EMAIL_LIST "") # Allow to be empty

  # Listing of subpakages
  SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS) # Allow to be empty

  #
  # B) Read in this package's Dependencies file and save off read dependency vars.
  #

  SET(PAKCAGE_DEPENDENCIES_FILE
    "${PROJECT_SOURCE_DIR}/${PACKAGE_DIR}/cmake/Dependencies.cmake")

  TRIBITS_TRACE_FILE_PROCESSING(PACKAGE  INCLUDE  "${PAKCAGE_DEPENDENCIES_FILE}")
  INCLUDE(${PAKCAGE_DEPENDENCIES_FILE})

  TRIBITS_ASSERT_READ_DEPENDENCY_VARS(${PACKAGE_NAME})

  TRIBITS_SAVE_OFF_DEPENENCIES_VARS(PARENTPACK)

  #
  # B.1) Set up the mail addresses
  #

  # ToDo: Move this above so that it will be handled as part of subpackage
  # processing?

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

  #
  # B.2) Process this package's subpackages first *before* finishing this packages!
  #

  TRIBITS_PARSE_SUBPACKAGES_AND_APPEND_SE_PACKAGES_AND_ADD_OPTIONS(${PACKAGE_NAME}
     ${PACKAGE_DIR})

  TRIBITS_READ_ALL_PACKAGE_SUBPACKAGE_DEPENDENCIES(${PACKAGE_NAME} ${PACKAGE_DIR})

  #
  # C) Finish processing this package's dependencies into dependency graph vars
  #
  # NOTE: The subpackages for this package are automatically treated as
  # optional or required library dependent packages for this outer package!
  #

  TRIBITS_READ_BACK_DEPENDENCIES_VARS(PARENTPACK)

  # Append the subpackages to the dependencies list
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


#
# Get the REPO_NAME and REPO_DIR given the REPO
#

FUNCTION(TRIBITS_GET_REPO_NAME_DIR  REPO_IN  REPO_NAME_OUT  REPO_DIR_OUT)
  #MESSAGE("TRIBITS_GET_REPO_NAME_DIR:  '${REPO_IN}'  '${REPO_NAME_OUT}'  '${REPO_DIR_OUT}'")
  # This list of repositories is the list of directories!
  SET(REPO_DIR ${REPO_IN})
  # Get the Repository name
  IF (REPO_IN STREQUAL ".")
    # The Project and the Reposiotry are one and the same
    SET(REPO_NAME ${PROJECT_NAME})
  ELSE()
    # The Repository name is the same as the repository directory
    SET(REPO_NAME ${REPO_IN})
  ENDIF()
  SET(${REPO_NAME_OUT} ${REPO_NAME} PARENT_SCOPE)
  SET(${REPO_DIR_OUT} ${REPO_DIR} PARENT_SCOPE)
ENDFUNCTION()


#
# Macro that reads all the package dependencies and builds dependency graph
#
# Reads from the variables:
#   ${PROJECT_NAME}_ALL_REPOSITORIES
#   ${PROJECT_NAME}_PACKAGES
#
# Writes to:
#   ${PROJECT_NAME}_SE_PACKAGES
#   ${PROJECT_NAME}_SE_PACKAGES_DIRS	
#
MACRO(TRIBITS_READ_ALL_PACKAGE_DEPENDENCIES)

  MESSAGE("")
  MESSAGE("Processing Project, Repository, and Package dependency files and building internal dependencies graph ...")
  MESSAGE("")

  #
  # A) First, process the Repository and Project dependency files
  #

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
        MESSAGE("-- " "The ${REPO_NAME} file ${REPO_DEPENDENCIES_SETUP_FILE} does not exist! ...")
      ENDIF()
    ENDIF()
  ENDFOREACH()

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
      MESSAGE("-- " "The ${PROJECT_NAME} file ${PROJECT_DEPENDENCIES_SETUP_FILE} does not exist! ...")
    ENDIF()
  ENDIF()

  #
  # B) Process the package dependency files, yielding the list of subpackages as well
  #

  SET(${PROJECT_NAME}_SE_PACKAGES) # Packages and subpackages
  SET(${PROJECT_NAME}_SE_PACKAGE_DIRS)

  SET(PACKAGE_IDX 0)
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES})
    LIST(GET ${PROJECT_NAME}_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
    TRIBITS_READ_PACKAGE_DEPENDENCIES(${TRIBITS_PACKAGE} ${PACKAGE_DIR})
    #TRIBITS_ADD_OPTIONAL_PACKAGE_ENABLES(${TRIBITS_PACKAGE})
    MATH(EXPR PACKAGE_IDX "${PACKAGE_IDX}+1")
  ENDFOREACH()

  # Create a reverse se packages list for later use
  SET(${PROJECT_NAME}_REVERSE_SE_PACKAGES ${${PROJECT_NAME}_SE_PACKAGES})
  IF (${PROJECT_NAME}_REVERSE_SE_PACKAGES)
    LIST(REVERSE ${PROJECT_NAME}_REVERSE_SE_PACKAGES)
  ENDIF()

  LIST(LENGTH ${PROJECT_NAME}_SE_PACKAGES ${PROJECT_NAME}_NUM_SE_PACKAGES)
  PRINT_VAR(${PROJECT_NAME}_NUM_SE_PACKAGES)
  #PRINT_VAR(${PROJECT_NAME}_SE_PACKAGES)

  FOREACH(TPL ${${PROJECT_NAME}_TPLS})
    IF (TPL_TENTATIVE_ENABLE_${TPL})
      MESSAGE("-- Tentatively enabling TPL '${TPL}'")
      #PRINT_VAR(TPL_ENABLE_${TPL})
    ENDIF()
  ENDFOREACH()

  ADVANCED_OPTION(${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES
    "Dump the package dependency information."
    "${${PROJECT_NAME}_VERBOSE_CONFIGURE}" )

  ADVANCED_OPTION(${PROJECT_NAME}_DUMP_FORWARD_PACKAGE_DEPENDENCIES
    "Dump the package forwrad dependency information."
    "${${PROJECT_NAME}_VERBOSE_CONFIGURE}" )

  IF (${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES)
    MESSAGE("")
    MESSAGE("Printing package dependencies ...")
    MESSAGE("")
    PRINT_NONEMPTY_VAR_WITH_SPACES(${PROJECT_NAME}_PACKAGES  DUMMY_OUT)
    MESSAGE("")
    PRINT_NONEMPTY_VAR_WITH_SPACES(${PROJECT_NAME}_SE_PACKAGES  DUMMY_OUT)
    MESSAGE("")
    FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_SE_PACKAGES})
      TRIBITS_PRINT_PACKAGE_DEPENDENCIES(${TRIBITS_PACKAGE})
      MESSAGE("")
    ENDFOREACH()
  ENDIF()

ENDMACRO()


#
# Function that sets a varaible to DECLARED-UNDEFINED
#

FUNCTION(DECLARE_UNDEFINED VAR)
  SET(${VAR} DECLARED-UNDEFINED PARENT_SCOPE)
ENDFUNCTION()


#
# Function that asserts that a package dependency variable is defined
# correctly
#

FUNCTION(TRIBITS_ASSERT_DEFINED_PACKAGE_VAR PACKAGE_VAR PACKAGE_NAME)
  IF (${PACKAGE_VAR} STREQUAL DECLARED-UNDEFINED)
    MESSAGE(FATAL_ERROR
      "Error, the package variable ${PACKAGE_VAR} was not defined correctly for package ${PACKAGE_NAME}!"
      )
  ENDIF()
ENDFUNCTION()


#
# Private helper macros
#


FUNCTION(TRIBITS_PRIVATE_PRINT_DISABLE
  ENABLE_BEING_DISABLED_VAR_NAME  PACKAGE_WITH_SOMETHING_BEING_DISABLED
  DEP_TYPE_STR  THING_DISALBED_TYPE  THING_DISABLED_NAME
  )
  #PRINT_VAR(${ENABLE_BEING_DISABLED_VAR_NAME})
  IF (${ENABLE_BEING_DISABLED_VAR_NAME})
    IF (${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES)
      MESSAGE(
        " ***\n"
        " *** NOTE: Setting ${ENABLE_BEING_DISABLED_VAR_NAME}=OFF"
        " which was '${${ENABLE_BEING_DISABLED_VAR_NAME}}' because"
        " ${PACKAGE_WITH_SOMETHING_BEING_DISABLED} has"
        " a required ${DEP_TYPE_STR} dependence on disabled"
        " ${THING_DISALBED_TYPE} ${THING_DISABLED_NAME}"
        " but ${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON!\n"
        " ***\n"
        )
    ELSE()
      MESSAGE(FATAL_ERROR
        " ***\n"
        " *** ERROR: Setting ${ENABLE_BEING_DISABLED_VAR_NAME}=OFF"
        " which was '${${ENABLE_BEING_DISABLED_VAR_NAME}}' because"
        " ${PACKAGE_WITH_SOMETHING_BEING_DISABLED} has"
        " a required ${DEP_TYPE_STR} dependence on disabled"
        " ${THING_DISALBED_TYPE} ${THING_DISABLED_NAME}!\n"
        " ***\n"
        )
    ENDIF()
  ELSE()
    MESSAGE("-- "
      "Setting ${ENABLE_BEING_DISABLED_VAR_NAME}=OFF"
      " because ${PACKAGE_WITH_SOMETHING_BEING_DISABLED} has a required ${DEP_TYPE_STR}"
      " dependence on disabled ${THING_DISALBED_TYPE} ${THING_DISABLED_NAME}")
  ENDIF()
ENDFUNCTION()


MACRO(TRIBITS_PRIVATE_DISABLE_TPL_REQUIRED_PACKAGE_ENABLE
  TPL_NAME  PACKAGE_NAME  LIBRARY_DEP
  )

  #MESSAGE("TRIBITS_PRIVATE_DISABLE_TPL_REQUIRED_PACKAGE_ENABLE"
  #  " ${TPL_NAME} ${PACKAGE_NAME} ${LIBRARY_DEP}")

  # Only turn off PACKAGE_NAME libraries and test/eamples if it
  # is currently enabled or could be enabled.

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
  IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}
     OR ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} STREQUAL ""
     )

    IF ("${LIBRARY_DEP}" STREQUAL "TRUE")

      TRIBITS_PRIVATE_PRINT_DISABLE(
        ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} ${PACKAGE_NAME} "library"
        "TPL" ${TPL_NAME}
        )

      SET(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} OFF)

    ELSE()

      SET(DEP_TYPE_STR "test/example")

      IF (${PACKAGE_NAME}_ENABLE_TESTS
        OR "${${PACKAGE_NAME}_ENABLE_TESTS}" STREQUAL ""
        )
        TRIBITS_PRIVATE_PRINT_DISABLE(
          ${PACKAGE_NAME}_ENABLE_TESTS ${PACKAGE_NAME} "${DEP_TYPE_STR}"
          "TPL" ${TPL_NAME}
          )
        SET(${PACKAGE_NAME}_ENABLE_TESTS OFF)
      ENDIF()

      IF (${PACKAGE_NAME}_ENABLE_EXAMPLES
        OR "${${PACKAGE_NAME}_ENABLE_EXAMPLES}" STREQUAL ""
        )
        TRIBITS_PRIVATE_PRINT_DISABLE(
          ${PACKAGE_NAME}_ENABLE_EXAMPLES ${PACKAGE_NAME} "${DEP_TYPE_STR}"
          "TPL" ${TPL_NAME}
          )
        SET(${PACKAGE_NAME}_ENABLE_EXAMPLES OFF)
      ENDIF()

      # NOTE: We can't assert that ${PACKAGE_NAME}_ENABLE_TESTS or
      # ${PACKAGE_NAME}_ENABLE_EXAMPLES exists yet because
      # TRIBITS_ADD_OPTIONAL_PACKAGE_ENABLES() which defines them is not
      # called until after the final package enables are set.

    ENDIF()

  ENDIF()

ENDMACRO()


FUNCTION(TRIBITS_PRIVATE_PRINT_DISABLE_REQUIRED_PACKAGE_ENABLE
  PACKAGE_NAME  PACKAGE_ENABLE_SOMETHING_VAR_NAME  FORWARD_DEP_PACKAGE_NAME
  DEP_TYPE_STR
  )
  TRIBITS_PRIVATE_PRINT_DISABLE(
    ${PACKAGE_ENABLE_SOMETHING_VAR_NAME} ${FORWARD_DEP_PACKAGE_NAME}
    "${DEP_TYPE_STR}" "package" ${PACKAGE_NAME} )
ENDFUNCTION()


MACRO(TRIBITS_PRIVATE_DISABLE_REQUIRED_PACKAGE_ENABLES
  FORWARD_DEP_PACKAGE_NAME PACKAGE_NAME LIBRARY_DEP
  )

  #MESSAGE("TRIBITS_PRIVATE_DISABLE_REQUIRED_PACKAGE_ENABLES"
  #  " ${FORWARD_DEP_PACKAGE_NAME} ${LIBRARY_DEP}")

  # Only turn off FORWARD_DEP_PACKAGE libraries and test/eamples if it
  # is currently enabled or could be enabled

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
  IF (${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME}
     OR ${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME} STREQUAL ""
     )

    IF ("${LIBRARY_DEP}" STREQUAL "TRUE")

      TRIBITS_PRIVATE_PRINT_DISABLE_REQUIRED_PACKAGE_ENABLE(
        ${PACKAGE_NAME} ${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME}
        ${FORWARD_DEP_PACKAGE_NAME} "library" )

      SET(${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME} OFF)

    ELSE()

      SET(DEP_TYPE_STR "test/example")

      IF (${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS
        OR "${${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS}" STREQUAL ""
        )
        TRIBITS_PRIVATE_PRINT_DISABLE_REQUIRED_PACKAGE_ENABLE(
          ${PACKAGE_NAME} ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS
          ${FORWARD_DEP_PACKAGE_NAME} "${DEP_TYPE_STR}" )
        SET(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS OFF)
      ENDIF()

      IF (${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES
        OR "${${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES}" STREQUAL ""
        )
        TRIBITS_PRIVATE_PRINT_DISABLE_REQUIRED_PACKAGE_ENABLE(
          ${PACKAGE_NAME} ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES
          ${FORWARD_DEP_PACKAGE_NAME} "${DEP_TYPE_STR}" )
        SET(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES OFF)
      ENDIF()

    ENDIF()

  ENDIF()

ENDMACRO()


MACRO(TRIBITS_PRIVATE_DISABLE_OPTIONAL_PACKAGE_ENABLES
  FORWARD_DEP_PACKAGE_NAME PACKAGE_NAME
  )

  #MESSAGE("TRIBITS_PRIVATE_DISABLE_OPTIONAL_PACKAGE_ENABLES"
  #  " ${FORWARD_DEP_PACKAGE_NAME} ${PACKAGE_NAME}")
  #MESSAGE("-- " "${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME} = ${${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME}}")

  #ASSERT_DEFINED(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME})
  IF (${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME} OR "${${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME}}" STREQUAL "")
    # Always disable the conditional enable but only print the message if the package is enabled.
    #MESSAGE("--  Disasble ${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME} ...")
    IF (${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
      IF (${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME})  # is explicity try already!
        MESSAGE("-- "
          "NOTE: Setting ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME}=OFF"
          " which was ${${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME}}"
          " because ${FORWARD_DEP_PACKAGE_NAME} has an optional library dependence"
          " on disabled package ${PACKAGE_NAME}")
      ELSE()  # Not explicitly set
        MESSAGE("-- "
          "Setting ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME}=OFF"
          " because ${FORWARD_DEP_PACKAGE_NAME} has an optional library dependence"
          " on disabled package ${PACKAGE_NAME}")
      ENDIF()
    ENDIF()
    SET(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME} OFF)
  ENDIF()

ENDMACRO()


#
# Macro that disabled a packages if its required upstream TPL is disabled..
#
MACRO(TRIBITS_DISABLE_PACKAGE_IF_TPL_DISABLED  TRIBITS_PACKAGE)

  FOREACH(TPL_NAME ${${TRIBITS_PACKAGE}_LIB_REQUIRED_DEP_TPLS})
    IF ( (NOT TPL_ENABLE_${TPL_NAME}) AND
      (NOT "${TPL_ENABLE_${TPL_NAME}}" STREQUAL "")
      )
      TRIBITS_PRIVATE_DISABLE_TPL_REQUIRED_PACKAGE_ENABLE(
        ${TPL_NAME}  ${TRIBITS_PACKAGE}  TRUE )
    ENDIF()
  ENDFOREACH()

  FOREACH(TPL_NAME ${${TRIBITS_PACKAGE}_TEST_REQUIRED_DEP_TPLS})
    IF ( (NOT TPL_ENABLE_${TPL_NAME}) AND
      (NOT "${TPL_ENABLE_${TPL_NAME}}" STREQUAL "")
      )
      TRIBITS_PRIVATE_DISABLE_TPL_REQUIRED_PACKAGE_ENABLE(
        ${TPL_NAME}  ${TRIBITS_PACKAGE}  FALSE )
    ENDIF()
  ENDFOREACH()

ENDMACRO()


#
# Macro that disables all of the subpackages of a parent package.
#
MACRO(TRIBITS_DISABLE_PARENTS_SUBPACKAGES PARENT_PACKAGE_NAME)

  #MESSAGE("TRIBITS_DISABLE_PARENTS_SUBPACKAGES: ${PARENT_PACKAGE_NAME}")

  #PRINT_VAR(${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME})

  IF(NOT ${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME}
    AND NOT ${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME} STREQUAL ""
    )

    SET(SUBPACKAGE_IDX 0)
    FOREACH(TRIBITS_SUBPACKAGE ${${PARENT_PACKAGE_NAME}_SUBPACKAGES})

      SET(SUBPACKAGE_NAME ${TRIBITS_SUBPACKAGE})
      SET(SUBPACKAGE_FULLNAME ${PARENT_PACKAGE_NAME}${TRIBITS_SUBPACKAGE})

      #PRINT_VAR(${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
      IF (NOT ${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME} STREQUAL "OFF")
        SET(ENABLE_BEING_DISABLED_VAR_NAME ${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
        MESSAGE("-- "
          "Setting subpackage enable ${ENABLE_BEING_DISABLED_VAR_NAME}=OFF"
          " because parent package ${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME}=OFF")
        SET(${ENABLE_BEING_DISABLED_VAR_NAME} OFF)
      ENDIF()

      MATH(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")

    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Macro that enables all of the subpackages of a parent package.
#
MACRO(TRIBITS_ENABLE_PARENTS_SUBPACKAGES PARENT_PACKAGE_NAME)

  #MESSAGE("TRIBITS_ENABLE_PARENTS_SUBPACKAGES: ${PARENT_PACKAGE_NAME}")

  #PRINT_VAR(${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME})

  IF(${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME})

    SET(SUBPACKAGE_IDX 0)
    FOREACH(TRIBITS_SUBPACKAGE ${${PARENT_PACKAGE_NAME}_SUBPACKAGES})

      SET(SUBPACKAGE_NAME ${TRIBITS_SUBPACKAGE})
      SET(SUBPACKAGE_FULLNAME ${PARENT_PACKAGE_NAME}${TRIBITS_SUBPACKAGE})

      #PRINT_VAR(${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})

      IF (NOT ${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME} AND
        NOT "${${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME}}" STREQUAL ""
        )
        # The subpackage is already disabled and is not just empty!
      ELSEIF (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
        # The subpackage is already enabled so there is no reason to enable it!
      ELSE()
        # The subpackage is not hard off or on so turn it on by default
        TRIBITS_IMPLICIT_PACKAGE_ENABLE_IS_ALLOWED( "" ${SUBPACKAGE_FULLNAME}
          SUBPACKAGE_ALLOW_IMPLICIT_ENABLE)
        IF (SUBPACKAGE_ALLOW_IMPLICIT_ENABLE)
          SET(ENABLE_BEING_ENABLED_VAR_NAME ${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
          MESSAGE("-- "
            "Setting subpackage enable ${ENABLE_BEING_ENABLED_VAR_NAME}=ON"
            " because parent package ${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME}=ON")
          SET(${ENABLE_BEING_ENABLED_VAR_NAME} ON)
        ENDIF()
      ENDIF()

      MATH(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")

    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Function that disables all forward packages that depend on the given packages
#

MACRO(TRIBITS_DISABLE_FORWARD_REQUIRED_DEP_PACKAGES PACKAGE_NAME)

  #MESSAGE("TRIBITS_DISABLE_FORWARD_REQUIRED_DEP_PACKAGES: ${PACKAGE_NAME}")

  IF (
     (NOT ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
     AND
     (NOT "${${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}}" STREQUAL "")
     )

    FOREACH(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES})
      TRIBITS_PRIVATE_DISABLE_REQUIRED_PACKAGE_ENABLES(${FWD_DEP_PKG} ${PACKAGE_NAME} TRUE)
    ENDFOREACH()

    FOREACH(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES})
      TRIBITS_PRIVATE_DISABLE_OPTIONAL_PACKAGE_ENABLES(${FWD_DEP_PKG} ${PACKAGE_NAME})
    ENDFOREACH()

    FOREACH(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES})
      TRIBITS_PRIVATE_DISABLE_REQUIRED_PACKAGE_ENABLES(${FWD_DEP_PKG} ${PACKAGE_NAME} FALSE)
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Macro that prints out dependencies for a package
#
# Does not modify the global state.
#

MACRO(TRIBITS_PRINT_PACKAGE_DEPENDENCIES PACKAGE_NAME)

  SET(PRINTED_VAR)

  PRINT_NONEMPTY_VAR_WITH_SPACES(${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES PRINTED_VAR)
  PRINT_NONEMPTY_VAR_WITH_SPACES(${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES PRINTED_VAR)
  PRINT_NONEMPTY_VAR_WITH_SPACES(${PACKAGE_NAME}_TEST_REQUIRED_DEP_PACKAGES PRINTED_VAR)
  PRINT_NONEMPTY_VAR_WITH_SPACES(${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES PRINTED_VAR)

  IF (${PROJECT_NAME}_DUMP_FORWARD_PACKAGE_DEPENDENCIES)
    PRINT_NONEMPTY_VAR_WITH_SPACES(${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES
      PRINTED_VAR)
    PRINT_NONEMPTY_VAR_WITH_SPACES(${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES
      PRINTED_VAR)
    PRINT_NONEMPTY_VAR_WITH_SPACES(${PACKAGE_NAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES
      PRINTED_VAR)
    PRINT_NONEMPTY_VAR_WITH_SPACES(${PACKAGE_NAME}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES
      PRINTED_VAR)
  ENDIF()

  PRINT_NONEMPTY_VAR_WITH_SPACES(${PACKAGE_NAME}_LIB_REQUIRED_DEP_TPLS PRINTED_VAR)
  PRINT_NONEMPTY_VAR_WITH_SPACES(${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS PRINTED_VAR)
  PRINT_NONEMPTY_VAR_WITH_SPACES(${PACKAGE_NAME}_TEST_REQUIRED_DEP_TPLS PRINTED_VAR)
  PRINT_NONEMPTY_VAR_WITH_SPACES(${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS PRINTED_VAR)

  IF (NOT PRINTED_VAR)
    MESSAGE("-- ${PACKAGE_NAME}: No dependencies!")
  ENDIF()


ENDMACRO()


#
# Private helper macros
#


MACRO(TRIBITS_PRIVATE_ADD_OPTIONAL_PACKAGE_ENABLE PACKAGE_NAME  OPTIONAL_DEP_PACKAGE
  TYPE  SET_AS_CACHE_IN
  )

  #MESSAGE("\nPACKAGE_ARCH_PRIVATE_ADD_OPTIONAL_PACKAGE_ENABLE: ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE}")

  IF (SET_AS_CACHE_IN)

    MULTILINE_SET(DOCSTR
      "Enable optional ${TYPE} support in the package ${PACKAGE_NAME}"
      " for the package ${OPTIONAL_DEP_PACKAGE}."
      "  Set to 'ON', 'OFF', or leave empty"
      " to allow for other logic to decide."
      )

    SET_CACHE_ON_OFF_EMPTY( ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE} ""
      ${DOCSTR} )

  ELSE()

    IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})
      SET( ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE} "" )
    ENDIF()

  ENDIF()

ENDMACRO()


MACRO(TRIBITS_PRIVATE_ADD_OPTIONAL_TPL_ENABLE PACKAGE_NAME OPTIONAL_DEP_TPL
  TYPE  SET_AS_CACHE_IN )

  IF (SET_AS_CACHE_IN)

    MULTILINE_SET(DOCSTR
      "Enable optional ${TYPE} support in the package ${PACKAGE_NAME}"
      " for the TPL ${OPTIONAL_DEP_TPL}."
      "  Set to 'ON', 'OFF', or leave empty"
      " to allow for other logic to decide."
      )

    SET_CACHE_ON_OFF_EMPTY( ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL} ""
      ${DOCSTR} )

  ELSE()

    IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL})
      SET( ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL} "" )
    ENDIF()

  ENDIF()

ENDMACRO()


#
# Macro that enables optional package interdependencies
#

MACRO(TRIBITS_ADD_OPTIONAL_PACKAGE_ENABLES PACKAGE_NAME)

  #MESSAGE("\nPACKAGE_ARCH_ADD_OPTIONAL_PACKAGE_ENABLES: ${PACKAGE_NAME}")

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
  SET(SET_AS_CACHE ${${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}})

  IF (SET_AS_CACHE)

    MULTILINE_SET(DOCSTR
      "Build tests for the package ${PACKAGE_NAME}.  Set to 'ON', 'OFF', or leave empty ''"
       " to allow for other logic to decide."
       )
    SET_CACHE_ON_OFF_EMPTY( ${PACKAGE_NAME}_ENABLE_TESTS "" ${DOCSTR} )

    MULTILINE_SET(DOCSTR
      "Build examples for the package ${PACKAGE_NAME}.  Set to 'ON', 'OFF', or leave empty ''"
       " to allow for other logic to decide."
       )
    SET_CACHE_ON_OFF_EMPTY( ${PACKAGE_NAME}_ENABLE_EXAMPLES "" ${DOCSTR} )

  ELSE()

    IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_TESTS)
      SET( ${PACKAGE_NAME}_ENABLE_TESTS "" )
    ENDIF()
    IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_EXAMPLES)
      SET( ${PACKAGE_NAME}_ENABLE_EXAMPLES "" )
    ENDIF()

  ENDIF()

  FOREACH(OPTIONAL_DEP_PACKAGE ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES})
    TRIBITS_PRIVATE_ADD_OPTIONAL_PACKAGE_ENABLE(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE} "library" "${SET_AS_CACHE}" )
  ENDFOREACH()

  FOREACH(OPTIONAL_DEP_PACKAGE ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES})
    TRIBITS_PRIVATE_ADD_OPTIONAL_PACKAGE_ENABLE(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE} "test" "${SET_AS_CACHE}" )
  ENDFOREACH()

  FOREACH(OPTIONAL_DEP_TPL ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS})
    TRIBITS_PRIVATE_ADD_OPTIONAL_TPL_ENABLE(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL} "library" "${SET_AS_CACHE}" )
  ENDFOREACH()

  FOREACH(OPTIONAL_DEP_TPL ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS})
    TRIBITS_PRIVATE_ADD_OPTIONAL_TPL_ENABLE(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL} "test" "${SET_AS_CACHE}" )
  ENDFOREACH()

ENDMACRO()


#
# Private helper macros
#


#
# Enable optional intra-package support for enabled target package
# ${PACKAGE_NAME} (i.e. ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} is assumed to
# be TRUE before calling this macro.
#
MACRO(TRIBITS_PRIVATE_POSTPROCESS_OPTIONAL_PACKAGE_ENABLE PACKAGE_NAME OPTIONAL_DEP_PACKAGE)

  #MESSAGE("TRIBITS_PRIVATE_POSTPROCESS_OPTIONAL_PACKAGE_ENABLE: ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE}")
  #PRINT_VAR(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})
  #PRINT_VAR(${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})

  IF (${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE} AND ${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})
    MESSAGE("-- " "NOTE:"
      " ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}}"
      " is already set!")
  ELSEIF ("${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}}" STREQUAL "")
    IF (${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})
      MESSAGE("-- " "Setting ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=ON"
       " since ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=ON AND"
       " ${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=ON")
      SET(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE} ON)
    ELSE()
      MESSAGE("-- " "NOT setting ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=ON"
       " since ${OPTIONAL_DEP_PACKAGE} is NOT enabled at this point!")
    ENDIF()
  ELSEIF (NOT "${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}}" STREQUAL ""
    AND NOT ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}
    AND ${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}
    )
    MESSAGE("-- " "NOTE: ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}}"
     " is already set so not enabling even though ${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=${${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}} is set!")
  ENDIF()

  STRING(TOUPPER ${PACKAGE_NAME} PACKAGE_NAME_UPPER)
  STRING(TOUPPER ${OPTIONAL_DEP_PACKAGE} OPTIONAL_DEP_PACKAGE_UPPER)
  SET(MACRO_DEFINE_NAME HAVE_${PACKAGE_NAME_UPPER}_${OPTIONAL_DEP_PACKAGE_UPPER})

  IF(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})
    SET(${MACRO_DEFINE_NAME} ON)
  ELSE()
    SET(${MACRO_DEFINE_NAME} OFF)
  ENDIF()

ENDMACRO()


#
# Enable optional intra-package support for enabled target package
# ${PACKAGE_NAME} (i.e. ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} is assumed to
# be TRUE before calling this macro.
#
MACRO(TRIBITS_PRIVATE_POSTPROCESS_OPTIONAL_TPL_ENABLE PACKAGE_NAME OPTIONAL_DEP_TPL)

  #MESSAGE("TRIBITS_PRIVATE_POSTPROCESS_OPTIONAL_TPL_ENABLE: ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL}")

  IF (${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL} AND TPL_ENABLE_${OPTIONAL_DEP_TPL})
    MESSAGE("-- " "NOTE:"
      " ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}=${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}}"
      " is already set!")
  ELSEIF (
    (NOT TPL_ENABLE_${OPTIONAL_DEP_TPL})
    AND
    (NOT "${TPL_ENABLE_${OPTIONAL_DEP_TPL}}" STREQUAL "")
    AND
    ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}
    )
    MESSAGE(
      "\n***"
      "\n*** NOTE: Setting ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}=OFF"
      " which was ON since TPL_ENABLE_${OPTIONAL_DEP_TPL}=OFF"
      "\n***\n"
      )
    SET(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL} OFF)
  ELSEIF ("${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}}" STREQUAL ""
    AND TPL_ENABLE_${OPTIONAL_DEP_TPL}
    )
    MESSAGE("-- " "Setting ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}=ON"
      " since TPL_ENABLE_${OPTIONAL_DEP_TPL}=ON")
    SET(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL} ON)
  ELSEIF (NOT "${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}}" STREQUAL ""
    AND NOT ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}
    AND TPL_ENABLE_${OPTIONAL_DEP_TPL}
    )
    MESSAGE("-- " "NOTE: ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}=${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}}"
      " is already set so not enabling even though TPL_ENABLE_${OPTIONAL_DEP_TPL}=${TPL_ENABLE_${OPTIONAL_DEP_TPL}} is set!")
  ENDIF()

  STRING(TOUPPER ${PACKAGE_NAME} PACKAGE_NAME_UPPER)
  STRING(TOUPPER ${OPTIONAL_DEP_TPL} OPTIONAL_DEP_TPL_UPPER)
  SET(MACRO_DEFINE_NAME HAVE_${PACKAGE_NAME_UPPER}_${OPTIONAL_DEP_TPL_UPPER})

  IF (${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL})
    SET(${MACRO_DEFINE_NAME} ON)
  ELSE()
    SET(${MACRO_DEFINE_NAME} OFF)
  ENDIF()

ENDMACRO()


#
# Macro that post-processes optional dependancies after all other
# dependencies have been worked out
#

MACRO(TRIBITS_POSTPROCESS_OPTIONAL_PACKAGE_ENABLES PACKAGE_NAME)

  #MESSAGE("\nPACKAGE_ARCH_POSTPROCESS_OPTIONAL_PACKAGE_ENABLES: ${PACKAGE_NAME}")

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
  IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    FOREACH(OPTIONAL_DEP_PACKAGE ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES})
      TRIBITS_PRIVATE_POSTPROCESS_OPTIONAL_PACKAGE_ENABLE(
        ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE} )
    ENDFOREACH()

    FOREACH(OPTIONAL_DEP_PACKAGE ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES})
      TRIBITS_PRIVATE_POSTPROCESS_OPTIONAL_PACKAGE_ENABLE(
        ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE} )
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Macro that post-processes final package enables for packages with subpackage
# enables.
#

MACRO(TRIBITS_POSTPROCESS_PACKAGE_WITH_SUBPACKAGES_ENABLES  PACKAGE_NAME)
  #MESSAGE("TRIBITS_POSTPROCESS_PACKAGE_WITH_SUBPACKAGES_ENABLES  '${PACKAGE_NAME}'")
  FOREACH(TRIBITS_SUBPACKAGE ${${PACKAGE_NAME}_SUBPACKAGES})
    SET(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
    #PRINT_VAR(${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
    #PRINT_VAR(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
    #PRINT_VAR(${SUBPACKAGE_FULLNAME}_ENABLE_TESTS)
    #PRINT_VAR(${PACKAGE_NAME}_ENABLE_TESTS)
    #PRINT_VAR(${SUBPACKAGE_FULLNAME}_ENABLE_EXAMPLES)
    #PRINT_VAR(${PACKAGE_NAME}_ENABLE_EXAMPLES)
    IF (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME}
        AND NOT ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}
      )
      MESSAGE("-- "
        "Setting ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=ON"
        " because ${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME}=ON")
      SET(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} ON)
      IF (${SUBPACKAGE_FULLNAME}_ENABLE_TESTS
        AND "${${PACKAGE_NAME}_ENABLE_TESTS}" STREQUAL ""
        )
        MESSAGE("-- "
          "Setting ${PACKAGE_NAME}_ENABLE_TESTS=ON"
          " because ${SUBPACKAGE_FULLNAME}_ENABLE_TESTS=ON")
        SET(${PACKAGE_NAME}_ENABLE_TESTS ON)
      ENDIF()
      IF (${SUBPACKAGE_FULLNAME}_ENABLE_EXAMPLES
        AND "${${PACKAGE_NAME}_ENABLE_EXAMPLES}" STREQUAL ""
        )
        MESSAGE("-- "
          "Setting ${PACKAGE_NAME}_ENABLE_EXAMPLES=ON"
          " because ${SUBPACKAGE_FULLNAME}_ENABLE_EXAMPLES=ON")
        SET(${PACKAGE_NAME}_ENABLE_EXAMPLES ON)
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDMACRO()


#
# Macro that post-processes optional package TPL based on if the TPL
# has been enabled or not
#

MACRO(TRIBITS_POSTPROCESS_OPTIONAL_TPL_ENABLES PACKAGE_NAME)

  #MESSAGE("\nPACKAGE_ARCH_ADD_OPTIONAL_TPL_ENABLES: ${PACKAGE_NAME}")
  
  IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    FOREACH(OPTIONAL_DEP_TPL ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS})
      TRIBITS_PRIVATE_POSTPROCESS_OPTIONAL_TPL_ENABLE(
        ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL} )
    ENDFOREACH()

    FOREACH(OPTIONAL_DEP_TPL ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS})
      TRIBITS_PRIVATE_POSTPROCESS_OPTIONAL_TPL_ENABLE(
        ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL} )
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Set an individual package variable enable based on the global value
#

MACRO(TRIBITS_SET_ALL_PACKAGES_PACKAGE_ENABLE_VARIABLE   PACKAGE_ARCH_VAR   PACKAGE_VAR)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("")
    MESSAGE("TRIBITS_SET_ALL_PACKAGES_PACKAGE_ENABLE_VARIABLE:")
    MESSAGE("-- " "${PACKAGE_ARCH_VAR} = ${${PACKAGE_ARCH_VAR}}")
    MESSAGE("-- " "${PACKAGE_VAR} = ${${PACKAGE_VAR}}")
  ENDIF()

  IF ("${${PACKAGE_VAR}}" STREQUAL "")
    IF (${PACKAGE_ARCH_VAR})
      MESSAGE("-- " "Setting ${PACKAGE_VAR}=ON")
      SET(${PACKAGE_VAR} ON)
    ELSEIF (
      (NOT ${PACKAGE_ARCH_VAR})
      AND
      (NOT "${PACKAGE_ARCH_VAR}" STREQUAL "")
      )
      MESSAGE("-- " "Setting ${PACKAGE_VAR}=OFF")
      SET(${PACKAGE_VAR} OFF)
    ELSE()
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("-- " "ELSE")
        # Otherwise, we will leave it up the the individual package
        # to decide?
      ENDIF()
    ENDIF()
  ELSE()
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("-- " "${PACKAGE_VAR} NOT DEFAULT")
    ENDIF()
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("-- " "${PACKAGE_VAR} = ${${PACKAGE_VAR}}")
  ENDIF()

ENDMACRO()


#
# Macro used to set ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} based on
# ${PROJECT_NAME}_ENABLE_ALL_PACKAGES
#

MACRO(TRIBITS_APPLY_ALL_PACKAGE_ENABLES  PACKAGE_NAME)
  TRIBITS_IS_PRIMARY_META_PROJECT_PACKAGE(${PACKAGE_NAME}  PACKAGE_IS_PMPP)
  TRIBITS_IMPLICIT_PACKAGE_ENABLE_IS_ALLOWED( "" ${PACKAGE_NAME}
    PROCESS_PACKAGE_ENABLE )
  IF (PACKAGE_IS_PMPP  AND  PROCESS_PACKAGE_ENABLE)
    TRIBITS_SET_ALL_PACKAGES_PACKAGE_ENABLE_VARIABLE(
      ${PROJECT_NAME}_ENABLE_ALL_PACKAGES ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} )
  ENDIF()
ENDMACRO()


#
# Macro used to set ${TRIBITS_PACKAGE)_ENABLE_TESTS and ${TRIBITS_PACKAGE)_ENABLE_EXAMPLES
# based on ${PROJECT_NAME}_ENABLE_ALL_PACKAGES
#

MACRO(TRIBITS_APPLY_TEST_EXAMPLE_ENABLES PACKAGE_NAME)
  IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
    TRIBITS_IS_PRIMARY_META_PROJECT_PACKAGE(${PACKAGE_NAME}  PACKAGE_IS_PMPP)
    IF (PACKAGE_IS_PMPP)
      TRIBITS_SET_ALL_PACKAGES_PACKAGE_ENABLE_VARIABLE(
        ${PROJECT_NAME}_ENABLE_TESTS ${PACKAGE_NAME}_ENABLE_TESTS )
      TRIBITS_SET_ALL_PACKAGES_PACKAGE_ENABLE_VARIABLE(
        ${PROJECT_NAME}_ENABLE_EXAMPLES ${PACKAGE_NAME}_ENABLE_EXAMPLES )
    ENDIF()
  ENDIF()
ENDMACRO()


#
# Private helper macro
#

MACRO(TRIBITS_PRIVATE_ENABLE_FORWARD_PACKAGE  FORWARD_DEP_PACKAGE_NAME  PACKAGE_NAME)
  TRIBITS_IMPLICIT_PACKAGE_ENABLE_IS_ALLOWED( "" ${FORWARD_DEP_PACKAGE_NAME}
    ALLOW_PACKAGE_ENABLE )
  #MESSAGE("TRIBITS_PRIVATE_ENABLE_FORWARD_PACKAGE: "
  #  "${FORWARD_DEP_PACKAGE_NAME} ${PACKAGE_NAME} ${ALLOW_PACKAGE_ENABLE}")
  # Enable the forward package if it is not already set to ON or OFF
  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
  IF(${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME} STREQUAL ""
    AND ALLOW_PACKAGE_ENABLE
    )
    MESSAGE("-- " "Setting ${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME}=ON"
      " because ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=ON")
    ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
    SET(${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME} ON)
  ENDIF()
ENDMACRO()


#
# Macro used to set ${PROJECT_NAME}_ENABLE_${FWD_PACKAGE_NAME)=ON for all optional
# and required forward library dependencies of the package ${PACKAGE_NAME}
#

MACRO(TRIBITS_ENABLE_FORWARD_LIB_PACKAGE_ENABLES PACKAGE_NAME)

  #MESSAGE("\nPACKAGE_ARCH_ENABLE_FORWARD_PACKAGE_ENABLES ${PACKAGE_NAME}")
  #PRINT_VAR(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

  # Enable the forward packages if this package is enabled
  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
  IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    FOREACH(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES})
      TRIBITS_PRIVATE_ENABLE_FORWARD_PACKAGE(${FWD_DEP_PKG} ${PACKAGE_NAME})
    ENDFOREACH()

    FOREACH(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES})
      TRIBITS_PRIVATE_ENABLE_FORWARD_PACKAGE(${FWD_DEP_PKG} ${PACKAGE_NAME})
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Macro used to set ${PROJECT_NAME}_ENABLE_${FWD_PACKAGE_NAME)=ON for all optional
# and required forward test/example dependencies of the package ${PACKAGE_NAME}
#

MACRO(TRIBITS_ENABLE_FORWARD_TEST_PACKAGE_ENABLES PACKAGE_NAME)

  #MESSAGE("\nPACKAGE_ARCH_ENABLE_FORWARD_PACKAGE_ENABLES ${PACKAGE_NAME}")
  #MESSAGE("-- " "${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=${${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}}")

  # Enable the forward packages if this package is enabled
  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
  IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    FOREACH(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES})
      TRIBITS_PRIVATE_ENABLE_FORWARD_PACKAGE(${FWD_DEP_PKG} ${PACKAGE_NAME})
    ENDFOREACH()

    FOREACH(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES})
      TRIBITS_PRIVATE_ENABLE_FORWARD_PACKAGE(${FWD_DEP_PKG} ${PACKAGE_NAME})
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Private helper macros
#

MACRO(TRIBITS_PRIVATE_ENABLE_DEP_PACKAGE  PACKAGE_NAME  DEP_PACKAGE_NAME
  OPTREQ_IN
  )

  #MESSAGE("TRIBITS_PRIVATE_ENABLE_DEP_PACKAGE:  '${PACKAGE_NAME}'  '${DEP_PACKAGE_NAME}'  '${OPTREQ_IN}'")

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME})
  #PRINT_VAR(${PACKAGE_NAME}_ENABLE_${DEP_PACKAGE_NAME})

  IF (${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME})

    #MESSAGE("The package is already enabled so there is nothing to enable!")

  ELSEIF (${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME} STREQUAL "")

    SET(TRIBITS_PRIVATE_ENABLE_DEP_PACKAGES_ENABLE_PACKAGE "")

    IF ("${OPTREQ_IN}" STREQUAL "REQUIRED")

      #MESSAGE("Always enable the upstream dependency if it is required")

      MESSAGE("-- " "Setting ${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME}=ON"
        " because ${PACKAGE_NAME} has a required dependence on ${DEP_PACKAGE_NAME}")

      SET(TRIBITS_PRIVATE_ENABLE_DEP_PACKAGES_ENABLE_PACKAGE ON)

    ELSEIF (${PACKAGE_NAME}_ENABLE_${DEP_PACKAGE_NAME})

      # Enable the upstream package if the user directly specified the
      # optional package enable reguardless if it is PT or ST or even EX.

      MESSAGE("-- " "Setting ${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME}=ON"
        " because ${PACKAGE_NAME}_ENABLE_${DEP_PACKAGE_NAME}=ON")

      SET(TRIBITS_PRIVATE_ENABLE_DEP_PACKAGES_ENABLE_PACKAGE ON)

    ELSEIF (${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES)

      # Enable the package if there is an optional dependence and we are asked
      # to enabled optional dependencies.

      TRIBITS_IMPLICIT_PACKAGE_ENABLE_IS_ALLOWED(${PACKAGE_NAME} ${DEP_PACKAGE_NAME}
        ALLOW_IMPLICIT_ENABLE)
      IF (ALLOW_IMPLICIT_ENABLE)
        MESSAGE("-- " "Setting ${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME}=ON"
          " because ${PACKAGE_NAME} has an optional dependence on ${DEP_PACKAGE_NAME}")
        SET(TRIBITS_PRIVATE_ENABLE_DEP_PACKAGES_ENABLE_PACKAGE ON)
      ENDIF()

    ENDIF()

    # Enable the upstream package
    IF (TRIBITS_PRIVATE_ENABLE_DEP_PACKAGES_ENABLE_PACKAGE)
      ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME})
      SET(${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME} ON)
    ENDIF()

  ENDIF()

ENDMACRO()


MACRO(TRIBITS_PRIVATE_ENABLE_DEP_TPL  PACKAGE_NAME  DEP_TPL_NAME)
  ASSERT_DEFINED(TPL_ENABLE_${DEP_TPL_NAME})
  IF(TPL_ENABLE_${DEP_TPL_NAME} STREQUAL "")
    MESSAGE("-- " "Setting TPL_ENABLE_${DEP_TPL_NAME}=ON because"
      " it is required by the enabled package ${PACKAGE_NAME}")
    ASSERT_DEFINED(TPL_ENABLE_${DEP_TPL_NAME})
    SET(TPL_ENABLE_${DEP_TPL_NAME} ON)
    SET(TPL_${DEP_TPL_NAME}_ENABLING_PKG  ${PACKAGE_NAME})
  ENDIF()
ENDMACRO()


MACRO(TRIBITS_PRIVATE_ENABLE_OPTIONAL_DEP_TPL PACKAGE_NAME DEP_TPL_NAME)
  #ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_${DEP_TPL_NAME})
  IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}
    AND ${PACKAGE_NAME}_ENABLE_${DEP_TPL_NAME}
    AND TPL_ENABLE_${DEP_TPL_NAME} STREQUAL ""
    )
    MESSAGE("-- " "Setting TPL_ENABLE_${DEP_TPL_NAME}=ON because"
      " ${PACKAGE_NAME}_ENABLE_${DEP_TPL_NAME}=ON")
    ASSERT_DEFINED(TPL_ENABLE_${DEP_TPL_NAME})
    SET(TPL_ENABLE_${DEP_TPL_NAME} ON)
  ENDIF()
ENDMACRO()


#
# Macro that enables the optional TPLs for given package
#

MACRO(TRIBITS_ENABLE_OPTIONAL_TPLS PACKAGE_NAME)

  #MESSAGE("TRIBITS_ENABLE_OPTIONAL_TPLS: ${PACKAGE_NAME}")
  #MESSAGE("-- " "${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=${${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}}")

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

  IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    FOREACH(DEP_TPL ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS})
      TRIBITS_PRIVATE_ENABLE_OPTIONAL_DEP_TPL(${PACKAGE_NAME} ${DEP_TPL})
    ENDFOREACH()

    FOREACH(DEP_TPL ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS})
      TRIBITS_PRIVATE_ENABLE_OPTIONAL_DEP_TPL(${PACKAGE_NAME} ${DEP_TPL})
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Macro that enables upstream (required and optional) SE packages given SE
# package
#
# Here I have to enable the required packages too or the logic just does no
# work as expected.
#
MACRO(TRIBITS_ENABLE_UPSTREAM_SE_PACKAGES PACKAGE_NAME)

  #MESSAGE("TRIBITS_ENABLE_UPSTREAM_SE_PACKAGES: ${PACKAGE_NAME}")
  #MESSAGE("-- " "${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=${${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}}")

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

  IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    FOREACH(DEP_PKG ${${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES})
      TRIBITS_PRIVATE_ENABLE_DEP_PACKAGE(${PACKAGE_NAME} ${DEP_PKG} REQUIRED)
    ENDFOREACH()

    FOREACH(DEP_PKG ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES})
      TRIBITS_PRIVATE_ENABLE_DEP_PACKAGE(${PACKAGE_NAME} ${DEP_PKG} OPTIONAL)
    ENDFOREACH()

    FOREACH(DEP_PKG ${${PACKAGE_NAME}_TEST_REQUIRED_DEP_PACKAGES})
      TRIBITS_PRIVATE_ENABLE_DEP_PACKAGE(${PACKAGE_NAME} ${DEP_PKG} REQUIRED)
    ENDFOREACH()

    FOREACH(DEP_PKG ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES})
      TRIBITS_PRIVATE_ENABLE_DEP_PACKAGE(${PACKAGE_NAME} ${DEP_PKG} OPTIONAL)
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Macro that sets the required TPLs for given package
#

MACRO(TRIBITS_ENABLE_REQUIRED_TPLS PACKAGE_NAME)

  #MESSAGE("PACKAGE_ARCH_ENABLE_REQUIRED_TPL_ENABLES: ${PACKAGE_NAME}")
  #MESSAGE("-- " "${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=${${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}}")

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

  IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    FOREACH(DEP_TPL ${${PACKAGE_NAME}_LIB_REQUIRED_DEP_TPLS})
      TRIBITS_PRIVATE_ENABLE_DEP_TPL(${PACKAGE_NAME} ${DEP_TPL})
    ENDFOREACH()

    FOREACH(DEP_TPL ${${PACKAGE_NAME}_TEST_REQUIRED_DEP_TPLS})
      TRIBITS_PRIVATE_ENABLE_DEP_TPL(${PACKAGE_NAME} ${DEP_TPL})
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Get the list of explicitly enabled entries
#
# These is the list of entires in ${LISTVAR} for which:
#
#   IF (${ENABLED_PREFIX}_ENABLE_{ENTRY})
#
# evaluates to true.
#
FUNCTION(TRIBITS_GET_ENABLED_LIST  LISTVAR  ENABLED_PREFIX  
  ENABLED_LIST_OUT_OUT  NUM_ENABLED_OUT_OUT
  )
  SET(ENABLED_LIST_OUT)
  FOREACH(ENTITY ${${LISTVAR}})
    SET(ENTITY_NAME ${ENABLED_PREFIX}_ENABLE_${ENTITY})
    ASSERT_DEFINED(${ENTITY_NAME})
    SET(INCLUDE_ENTITY FALSE)
    IF (${ENTITY_NAME})
      LIST(APPEND  ENABLED_LIST_OUT  ${ENTITY})
    ENDIF()
  ENDFOREACH()
  LIST(LENGTH  ENABLED_LIST_OUT  NUM_ENABLED_OUT)
  SET(${ENABLED_LIST_OUT_OUT} ${ENABLED_LIST_OUT} PARENT_SCOPE)
  IF (NUM_ENABLED_OUT_OUT)
    SET(${NUM_ENABLED_OUT_OUT} ${NUM_ENABLED_OUT} PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()


#
# Get the list non-disabled entries
#
# These is the list of entires in ${LISTVAR} for which:
#
#   IF (
#     (${ENABLED_PREFIX}_ENABLE_{ENTRY})
#     OR
#     (${ENABLED_PREFIX}_ENABLE_{ENTRY} STREQUAL "" )
#     )
#
# evaluates to true.
#
FUNCTION(TRIBITS_GET_NONDISABLED_LIST  LISTVAR  ENABLED_PREFIX  
  NONDISABLED_LIST_OUT_OUT  NUM_NONDISABLED_OUT_OUT
  )
  SET(NONDISABLED_LIST_OUT)
  FOREACH(ENTITY ${${LISTVAR}})
    SET(ENTITY_NAME ${ENABLED_PREFIX}_ENABLE_${ENTITY})
    ASSERT_DEFINED(${ENTITY_NAME})
    SET(INCLUDE_ENTITY FALSE)
    IF (${ENTITY_NAME} OR ${ENTITY_NAME} STREQUAL "")
      LIST(APPEND  NONDISABLED_LIST_OUT  ${ENTITY})
    ENDIF()
  ENDFOREACH()
  LIST(LENGTH  NONDISABLED_LIST_OUT  NUM_NONDISABLED_OUT)
  SET(${NONDISABLED_LIST_OUT_OUT} ${NONDISABLED_LIST_OUT} PARENT_SCOPE)
  IF (NUM_NONDISABLED_OUT_OUT)
    SET(${NUM_NONDISABLED_OUT_OUT} ${NUM_NONDISABLED_OUT} PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()


#
# Get the list of explicitly disabled entries
#
# These is the list of entires in ${LISTVAR} for which:
#
#   IF (
#     (NOT ${ENABLED_PREFIX}_ENABLE_{ENTRY})
#     AND
#     (NOT ${ENABLED_PREFIX}_ENABLE_{ENTRY} STREQUAL "" )
#     )
#
# evaluates to true.
#
FUNCTION(TRIBITS_GET_DISABLED_LIST  LISTVAR  ENABLED_PREFIX  
  DISABLED_LIST_OUT_OUT  NUM_DISABLED_OUT_OUT
  )
  SET(DISABLED_LIST_OUT)
  FOREACH(ENTITY ${${LISTVAR}})
    SET(ENTITY_NAME ${ENABLED_PREFIX}_ENABLE_${ENTITY})
    ASSERT_DEFINED(${ENTITY_NAME})
    SET(INCLUDE_ENTITY FALSE)
    IF ( (NOT ${ENTITY_NAME}) AND (NOT ${ENTITY_NAME} STREQUAL "") )
      LIST(APPEND  DISABLED_LIST_OUT  ${ENTITY})
    ENDIF()
  ENDFOREACH()
  LIST(LENGTH  DISABLED_LIST_OUT  NUM_DISABLED_OUT)
  SET(${DISABLED_LIST_OUT_OUT} ${DISABLED_LIST_OUT} PARENT_SCOPE)
  IF (NUM_DISABLED_OUT_OUT)
    SET(${NUM_DISABLED_OUT_OUT} ${NUM_DISABLED_OUT} PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()


#
# Get the list of non-enabled entries
#
# These is the list of entires in ${LISTVAR} for which:
#
#   IF (NOT ${ENABLED_PREFIX}_ENABLE_{ENTRY})

# evaluates to true.
#
FUNCTION(TRIBITS_GET_NONENABLED_LIST  LISTVAR  ENABLED_PREFIX  
  NONENABLED_LIST_OUT_OUT  NUM_NONENABLED_OUT_OUT
  )
  SET(NONENABLED_LIST_OUT)
  FOREACH(ENTITY ${${LISTVAR}})
    SET(ENTITY_NAME ${ENABLED_PREFIX}_ENABLE_${ENTITY})
    ASSERT_DEFINED(${ENTITY_NAME})
    SET(INCLUDE_ENTITY FALSE)
    IF (NOT ${ENTITY_NAME}) # Note that empty "" is also false!
      LIST(APPEND  NONENABLED_LIST_OUT  ${ENTITY})
    ENDIF()
  ENDFOREACH()
  LIST(LENGTH  NONENABLED_LIST_OUT  NUM_NONENABLED_OUT)
  SET(${NONENABLED_LIST_OUT_OUT} ${NONENABLED_LIST_OUT} PARENT_SCOPE)
  IF (NUM_NONENABLED_OUT_OUT)
    SET(${NUM_NONENABLED_OUT_OUT} ${NUM_NONENABLED_OUT} PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()


#
# Macro that sets up the basic lists of enabled packages and SE packages.
#
MACRO(TRIBITS_SET_UP_ENABLED_LISTS_AND_SE_PKG_IDX)

  # ${PROJECT_NAME}_ENABLED_PACKAGES
  TRIBITS_GET_ENABLED_LIST( ${PROJECT_NAME}_PACKAGES  ${PROJECT_NAME}
    ${PROJECT_NAME}_ENABLED_PACKAGES  ${PROJECT_NAME}_NUM_ENABLED_PACKAGES)

  # ${PROJECT_NAME}_ENABLED_SE_PACKAGES
  TRIBITS_GET_ENABLED_LIST( ${PROJECT_NAME}_SE_PACKAGES  ${PROJECT_NAME}
    ${PROJECT_NAME}_ENABLED_SE_PACKAGES  ${PROJECT_NAME}_NUM_ENABLED_SE_PACKAGES)

  # ${PROJECT_NAME}_REVERSE_ENABLED_SE_PACKAGES
  SET(${PROJECT_NAME}_REVERSE_ENABLED_SE_PACKAGES
    "${${PROJECT_NAME}_ENABLED_SE_PACKAGES}")
  LIST(REVERSE ${PROJECT_NAME}_REVERSE_ENABLED_SE_PACKAGES)

  # ${PACKAGE_NAME}_SE_PKG_IDX
  SET(SE_PKG_IDX 0)
  FOREACH(TRIBITS_SE_PACKAGE ${${PROJECT_NAME}_ENABLED_SE_PACKAGES})
    SET(${TRIBITS_SE_PACKAGE}_SE_PKG_IDX ${SE_PKG_IDX})
    MATH(EXPR  SE_PKG_IDX  "${SE_PKG_IDX} + 1")
  ENDFOREACH()

  # ${PROJECT_NAME}_ENABLED_TPLS
  TRIBITS_GET_ENABLED_LIST( ${PROJECT_NAME}_TPLS  TPL
    ${PROJECT_NAME}_ENABLED_TPLS  ${PROJECT_NAME}_NUM_ENABLED_TPLS)

  # ${PROJECT_NAME}_REVERSE_ENABLED_TPLS
  SET(${PROJECT_NAME}_REVERSE_ENABLED_TPLS
    "${${PROJECT_NAME}_ENABLED_TPLS}")
  LIST(REVERSE ${PROJECT_NAME}_REVERSE_ENABLED_TPLS)

ENDMACRO()


#
# Macro that adjusts all of the package enables from what the user input
# to the final set that will be used to enable packages
#
MACRO(TRIBITS_ADJUST_PACKAGE_ENABLES)

  IF (${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES)
    MESSAGE("")
    MESSAGE("Setting to empty '' all enabled packages on reqeust ...")
    MESSAGE("")
    FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_SE_PACKAGES})
      IF (${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
        SET_CACHE_ON_OFF_EMPTY(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} ""
          "Forced to empty '' by ${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES=ON" FORCE)
        SET(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} "")
      ENDIF()
      #PRINT_VAR(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
      # NOTE: Above, we don't want to set to empty those packages that have hard
      # disables because this will mess up the logic in later invocations.
    ENDFOREACH()
    ADVANCED_SET(${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES OFF CACHE BOOL
      "Forced to FALSE after use" FORCE)
  ENDIF()

  #
  # A) Sweep forward through and apply all disables first!
  #

  TRIBITS_GET_NONDISABLED_LIST( ${PROJECT_NAME}_SE_PACKAGES  ${PROJECT_NAME}
    ${PROJECT_NAME}_NOTDISABLED_SE_PACKAGES "")

  MESSAGE("")
  MESSAGE("Disabling all packages that have a required dependency"
    " on disabled TPLs and optional package TPL support based on TPL_ENABLE_<TPL>=OFF ...")
  MESSAGE("")
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_NOTDISABLED_SE_PACKAGES})
    TRIBITS_DISABLE_PACKAGE_IF_TPL_DISABLED(${TRIBITS_PACKAGE})
  ENDFOREACH()

  MESSAGE("")
  MESSAGE("Disabling subpackages for hard disables of parent packages"
    " due to ${PROJECT_NAME}_ENABLE_<PARENT_PACKAGE>=OFF ...")
  MESSAGE("")
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_SE_PACKAGES})
    TRIBITS_DISABLE_PARENTS_SUBPACKAGES(${TRIBITS_PACKAGE})
  ENDFOREACH()

  MESSAGE("")
  MESSAGE("Disabling forward required SE packages and optional intra-package"
    " support that have a dependancy on disabled SE packages"
    " ${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=OFF ...")
  MESSAGE("")
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_SE_PACKAGES})
    TRIBITS_DISABLE_FORWARD_REQUIRED_DEP_PACKAGES(${TRIBITS_PACKAGE})
  ENDFOREACH()

  TRIBITS_GET_NONDISABLED_LIST( ${PROJECT_NAME}_NOTDISABLED_SE_PACKAGES  ${PROJECT_NAME}
    ${PROJECT_NAME}_NOTDISABLED_SE_PACKAGES "")

  SET(${PROJECT_NAME}_REVERSE_NOTDISABLED_SE_PACKAGES
    "${${PROJECT_NAME}_NOTDISABLED_SE_PACKAGES}")
  LIST(REVERSE ${PROJECT_NAME}_REVERSE_NOTDISABLED_SE_PACKAGES)

  #
  # B) Apply all forward enables
  #

  MESSAGE("")
  MESSAGE("Enabling subpackages for hard enables of parent packages"
    " due to ${PROJECT_NAME}_ENABLE_<PARENT_PACKAGE>=ON ...")
  MESSAGE("")
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_NOTDISABLED_SE_PACKAGES})
    TRIBITS_ENABLE_PARENTS_SUBPACKAGES(${TRIBITS_PACKAGE})
  ENDFOREACH()

  IF (${PROJECT_NAME}_ENABLE_ALL_PACKAGES)
    MESSAGE("")
    MESSAGE("Enabling all SE packages that are not currently disabled because of"
      " ${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON"
      " (${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE=${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE})"
      " ...")
    MESSAGE("")
    FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_NOTDISABLED_SE_PACKAGES})
      TRIBITS_APPLY_ALL_PACKAGE_ENABLES(${TRIBITS_PACKAGE})
    ENDFOREACH()
  ENDIF()

  IF (${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES)
    MESSAGE("")
    MESSAGE("Sweep forward enabling all forward library dependent packages because"
      " ${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON ...")
    MESSAGE("")
    FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_NOTDISABLED_SE_PACKAGES})
      TRIBITS_ENABLE_FORWARD_LIB_PACKAGE_ENABLES(${TRIBITS_PACKAGE})
    ENDFOREACH()
    MESSAGE("")
    MESSAGE("Sweep backward enabling all forward test dependent packages because"
      " ${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON ...")
    MESSAGE("")
    FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_REVERSE_NOTDISABLED_SE_PACKAGES})
      TRIBITS_ENABLE_FORWARD_TEST_PACKAGE_ENABLES(${TRIBITS_PACKAGE})
    ENDFOREACH()
    # NOTE: Above, we want to sweep backward to enable test-dependent packages
    # because we don't want to enable package Z just because package Y was enabled
    # because it had a test-only dependency on package X.  Sweeping backwards through
    # the packages makes sure this does not happen.
    SET(${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES ON)
  ENDIF()

  TRIBITS_GET_ENABLED_LIST( ${PROJECT_NAME}_NOTDISABLED_SE_PACKAGES  ${PROJECT_NAME}
    ${PROJECT_NAME}_ENABLED_SE_PACKAGES  "")

  #
  # C) Enable tests for currently enabled SE packages
  #

  IF (${PROJECT_NAME}_ENABLE_TESTS OR ${PROJECT_NAME}_ENABLE_EXAMPLES)
    MESSAGE("")
    MESSAGE("Enabling all tests and/or examples that have not been"
      " explicitly disabled because ${PROJECT_NAME}_ENABLE_[TESTS,EXAMPLES]=ON ...")
    MESSAGE("")
    FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_SE_PACKAGES})
      TRIBITS_APPLY_TEST_EXAMPLE_ENABLES(${TRIBITS_PACKAGE})
    ENDFOREACH()
  ENDIF()
  # NOTE: Above, we enable tests and examples here, before the remaining required
  # packages so that we don't enable tests that don't need to be enabled based
  # on the use of the option ${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES.

  #
  # D) Sweep backwards and enable upstream required and optional SE packages
  #

  IF (${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES)
    SET(EXTRA_MSG_STR " (and optional since ${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES=ON)")
  ELSE()
    SET(EXTRA_MSG_STR "")
  ENDIF()

  MESSAGE("")
  MESSAGE("Enabling all required${EXTRA_MSG_STR} upstream SE packages for current set of"
    " enabled packages"
    " (${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE=${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE})"
    " ...")
  MESSAGE("")
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_REVERSE_NOTDISABLED_SE_PACKAGES})
    TRIBITS_ENABLE_UPSTREAM_SE_PACKAGES(${TRIBITS_PACKAGE})
  ENDFOREACH()
  # NOTE: Above, we have to loop through the packages backward to enable all
  # the packages that feed into these packages.  This has to include *all*
  # upstream SE package enables including required SE packages, optional SE
  # packages (when ${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES), and SE
  # packages

  TRIBITS_GET_ENABLED_LIST( ${PROJECT_NAME}_NOTDISABLED_SE_PACKAGES  ${PROJECT_NAME}
    ${PROJECT_NAME}_ENABLED_SE_PACKAGES  "")

  MESSAGE("")
  MESSAGE("Enabling all optional intra-package enables <TRIBITS_PACKAGE>_ENABLE_<DEPPACKAGE>"
    " that are not currently disabled if both sets of packages are enabled ...")
  MESSAGE("")
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_SE_PACKAGES})
    TRIBITS_POSTPROCESS_OPTIONAL_PACKAGE_ENABLES(${TRIBITS_PACKAGE})
  ENDFOREACH()

  #
  # E) Enable TPLs
  #

  MESSAGE("")
  MESSAGE("Enabling all remaining required TPLs for current set of"
    " enabled packages ...")
  MESSAGE("")
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_SE_PACKAGES})
    TRIBITS_ENABLE_REQUIRED_TPLS(${TRIBITS_PACKAGE})
  ENDFOREACH()

  MESSAGE("")
  MESSAGE("Enabling all optional package TPL support"
    " <TRIBITS_PACKAGE>_ENABLE_<DEPTPL> not currently disabled for"
    " enabled TPLs ...")
  MESSAGE("")
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_SE_PACKAGES})
    TRIBITS_POSTPROCESS_OPTIONAL_TPL_ENABLES(${TRIBITS_PACKAGE})
  ENDFOREACH()

  MESSAGE("")
  MESSAGE("Enabling TPLs based on <TRIBITS_PACKAGE>_ENABLE_<TPL>=ON if TPL is not explicitly disabled ...")
  MESSAGE("")
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_SE_PACKAGES})
    TRIBITS_ENABLE_OPTIONAL_TPLS(${TRIBITS_PACKAGE})
  ENDFOREACH()
  # NOTE: We need to do this after the above optional package TPL support
  # logic so that the TPL will be turned on for this package only as requested
  # in bug 4298.

  #
  # F) Set user cache variables for current set of enabled SE packages
  #

  MESSAGE("")
  MESSAGE("Set cache entries for optional packages/TPLs and tests/examples for packages actually enabled ...")
  MESSAGE("")
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_SE_PACKAGES})
    TRIBITS_ADD_OPTIONAL_PACKAGE_ENABLES(${TRIBITS_PACKAGE})
  ENDFOREACH()

  #
  # G) Turn on parent packages where at least one subpackage has been enabled
  #

  MESSAGE("")
  MESSAGE("Enabling the shell of non-enabled parent packages (mostly for show) that have at least one subpackage enabled ...")
  MESSAGE("")
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES})
    TRIBITS_POSTPROCESS_PACKAGE_WITH_SUBPACKAGES_ENABLES(${TRIBITS_PACKAGE})
  ENDFOREACH()
  # NOTE: The above ensures that loops involving the parent package will
  # process the parent package but doing this last ensures that no downstream
  # dependencies will be enabled.

  TRIBITS_SET_UP_ENABLED_LISTS_AND_SE_PKG_IDX()

ENDMACRO()


#
# Function that sets up the full package dependencies for each enabled
# package.
#
# This is needed in several different parts of the TriBITS implementation.
#
FUNCTION(TRIBITS_PACKAGE_SET_FULL_ENABLED_DEP_PACKAGES  PACKAGE_NAME)

  FOREACH(DEP_PKG ${${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES})
    IF (${PROJECT_NAME}_ENABLE_${DEP_PKG})
      LIST(APPEND  PACKAGE_FULL_DEPS_LIST  ${DEP_PKG})
    ENDIF()
    # NOTE: This if() should not be needed but this is a safeguard
  ENDFOREACH()

  FOREACH(DEP_PKG ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES})
    IF (${PACKAGE_NAME}_ENABLE_${DEP_PKG})
      LIST(APPEND  PACKAGE_FULL_DEPS_LIST  ${DEP_PKG})
    ENDIF()
  ENDFOREACH()

  IF(PACKAGE_FULL_DEPS_LIST)
    LIST(REMOVE_DUPLICATES  PACKAGE_FULL_DEPS_LIST)

    FOREACH(DEP_PACKAGE  ${PACKAGE_FULL_DEPS_LIST})
      LIST(APPEND PACKAGE_FULL_DEPS_LIST  ${${DEP_PACKAGE}_FULL_ENABLED_DEP_PACKAGES})
    ENDFOREACH()

    LIST(REMOVE_DUPLICATES PACKAGE_FULL_DEPS_LIST)
  ENDIF()

  SET(ORDERED_PACKAGE_FULL_DEPS_LIST)

  FOREACH(DEP_PACKAGE  ${PACKAGE_FULL_DEPS_LIST})

    #PRINT_VAR(${DEP_PACKAGE}_SE_PKG_IDX)
    SET(DEP_PACKAGE_VALUE  ${${DEP_PACKAGE}_SE_PKG_IDX})

    SET(SORTED_INDEX 0)
    SET(INSERTED_DEP_PACKAGE FALSE)

    FOREACH(SORTED_PACKAGE  ${ORDERED_PACKAGE_FULL_DEPS_LIST})

      #PRINT_VAR(${SORTED_PACKAGE}_SE_PKG_IDX)
      SET(SORTED_PACKAGE_VALUE  ${${SORTED_PACKAGE}_SE_PKG_IDX})

      IF (${DEP_PACKAGE_VALUE} GREATER ${SORTED_PACKAGE_VALUE})
        LIST(INSERT  ORDERED_PACKAGE_FULL_DEPS_LIST  ${SORTED_INDEX}  ${DEP_PACKAGE})
        SET(INSERTED_DEP_PACKAGE TRUE)
        BREAK()
      ENDIF()

      MATH(EXPR SORTED_INDEX ${SORTED_INDEX}+1)

    ENDFOREACH()

    IF(NOT INSERTED_DEP_PACKAGE)
      LIST(APPEND  ORDERED_PACKAGE_FULL_DEPS_LIST  ${DEP_PACKAGE})
    ENDIF()

  ENDFOREACH()

  GLOBAL_SET(${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES
    ${ORDERED_PACKAGE_FULL_DEPS_LIST})

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES)
  ENDIF()

ENDFUNCTION()


#
# Function that creates enable-only dependency data-structures
#
FUNCTION(TRIBITS_SET_UP_ENABLED_ONLY_DEPENDENCIES)

  SET(GENERATE_EXPORT_DEPENDENCIES ${${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES})
  SET(LAST_EXPORT_SE_PACKAGE)

  IF (GENERATE_EXPORT_DEPENDENCIES
      AND ${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES
    )
    # Find the last enabled SE package for which an export file is requested.
    SET(LAST_SE_PKG_IDX -1)
    SET(LAST_SE_PKG)
    FOREACH(SE_PKG ${${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES})
      #PRINT_VAR(SE_PKG)
      SET(SE_PKG_IDX ${${SE_PKG}_SE_PKG_IDX})
      #PRINT_VAR(SE_PKG_IDX)
      IF (SE_PKG_IDX)
        # The listed package is enabled so we will consider it
        IF (SE_PKG_IDX GREATER ${LAST_SE_PKG_IDX})
          SET(LAST_SE_PKG_IDX ${SE_PKG_IDX})
          SET(LAST_SE_PKG ${SE_PKG})
         #PRINT_VAR(LAST_SE_PKG_IDX)
         #PRINT_VAR(LAST_SE_PKG)
        ENDIF()
      ENDIF()
    ENDFOREACH()
    IF (LAST_SE_PKG)
      # At least one listed package was enabled
      SET(LAST_EXPORT_SE_PACKAGE ${LAST_SE_PKG})
    ELSE()
      # None of the listed packages were enabled so don't bother generating
      # any export dependencies
      SET(GENERATE_EXPORT_DEPENDENCIES FALSE)
    ENDIF()

  ENDIF()

  IF (GENERATE_EXPORT_DEPENDENCIES)

    IF (LAST_EXPORT_SE_PACKAGE)
      MESSAGE("\nSetting up export dependencies up through ${LAST_EXPORT_SE_PACKAGE} ...\n")
    ELSE()
      MESSAGE("\nSetting up export dependencies for all enabled SE packages ...\n")
    ENDIF()

    FOREACH(TRIBITS_SE_PACKAGE ${${PROJECT_NAME}_ENABLED_SE_PACKAGES})
      TRIBITS_PACKAGE_SET_FULL_ENABLED_DEP_PACKAGES(${TRIBITS_SE_PACKAGE})
      IF (${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES)
        SET(PRINTED_VAR FALSE)
        PRINT_NONEMPTY_VAR_WITH_SPACES(${TRIBITS_SE_PACKAGE}_FULL_ENABLED_DEP_PACKAGES
          PRINTED_VAR)
        IF (NOT PRINTED_VAR)
          MESSAGE("-- ${TRIBITS_SE_PACKAGE}: No library dependencies!")
        ENDIF()
      ENDIF()
      IF ("${LAST_EXPORT_SE_PACKAGE}" STREQUAL ${TRIBITS_SE_PACKAGE})
        BREAK()
      ENDIF()
    ENDFOREACH()

  ENDIF()

ENDFUNCTION()
