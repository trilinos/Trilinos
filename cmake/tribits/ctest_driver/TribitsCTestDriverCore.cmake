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

#
# Tribits platform-independent test driver.
#
#

MESSAGE("")
MESSAGE("*******************************")
MESSAGE("*** TribitsCTestDriverCore ***")
MESSAGE("*******************************")
MESSAGE("")


CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11 FATAL_ERROR)

#
# Get the basic variables that define the project and the build
#

#
# Set TRIBITS_PROJECT_ROOT
#
# We must locate the source code root directory before processing this
# file. Once the root directory is located, we can make good guesses
# at other properties, but this file is a CTest script that is always
# executed outside of the context of the project's CMake build.
#
# We allow the environment variable TRIBITS_PROJECT_ROOT to locate the
# root directory. If the variable doesn't exist, we fall back on the
# default convention.
#MESSAGE("TRIBITS_PROJECT_ROOT (before env) = '${TRIBITS_PROJECT_ROOT}'")
IF (NOT TRIBITS_PROJECT_ROOT)
  SET(TRIBITS_PROJECT_ROOT "$ENV{TRIBITS_PROJECT_ROOT}")
ENDIF()
#MESSAGE("TRIBITS_PROJECT_ROOT (after env) = '${TRIBITS_PROJECT_ROOT}'")
IF (NOT TRIBITS_PROJECT_ROOT)
  # Fall back on the default convention, in which this file is located at:
  #   <root>/cmake/tribits/ctest.
  GET_FILENAME_COMPONENT(CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
  SET(TRIBITS_PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../..")
ENDIF()
GET_FILENAME_COMPONENT(TRIBITS_PROJECT_ROOT "${TRIBITS_PROJECT_ROOT}" ABSOLUTE)
MESSAGE("TRIBITS_PROJECT_ROOT = '${TRIBITS_PROJECT_ROOT}'")

#
# Read in PROJECT_NAME
#
# Assert that the ProjectName.cmake file exists.
SET(TRIBITS_PROJECT_NAME_INCLUDE "${TRIBITS_PROJECT_ROOT}/ProjectName.cmake")
IF(NOT EXISTS "${TRIBITS_PROJECT_NAME_INCLUDE}")
  MESSAGE(FATAL_ERROR
    "Could not locate ProjectName.cmake.\n"
    "  TRIBITS_PROJECT_ROOT = ${TRIBITS_PROJECT_ROOT}\n"
    "  Set the TRIBITS_PROJECT_ROOT environment variable "
    "to point at the source root.")
ENDIF()
# Include the ProjectName.cmake file and get PROJECT_NAME
INCLUDE(${TRIBITS_PROJECT_NAME_INCLUDE})
IF(NOT PROJECT_NAME)
  MESSAGE(FATAL_ERROR
    "The project name has not been set!"
    "  It should be set in ${TRIBITS_PROJECT_ROOT}/ProjectName.cmake.")
ENDIF()
MESSAGE("PROJECT_NAME = ${PROJECT_NAME}")

#
# Set ${PROJECT_NAME}_TRIBITS_DIR
#
IF (NOT ${PROJECT_NAME}_TRIBITS_DIR)
  SET(${PROJECT_NAME}_TRIBITS_DIR "$ENV{${PROJECT_NAME}_TRIBITS_DIR}")
ENDIF()
IF (NOT ${PROJECT_NAME}_TRIBITS_DIR)
  SET(${PROJECT_NAME}_TRIBITS_DIR "${TRIBITS_PROJECT_ROOT}/cmake//tribits")
ENDIF()
MESSAGE("${PROJECT_NAME}_TRIBITS_DIR = ${${PROJECT_NAME}_TRIBITS_DIR}")

#
# Set default for CTEST_SOURCE_DIRECTORY
#
IF ("${CTEST_SOURCE_DIRECTORY}" STREQUAL "")
  MESSAGE("Set default for CTEST_SOURCE_DIRECTORY to TRIBITS_PROJECT_ROOT='TRIBITS_PROJECT_ROOT='${TRIBITS_PROJECT_ROOT}'")
  SET(CTEST_SOURCE_DIRECTORY ${TRIBITS_PROJECT_ROOT})
ENDIF()

#
# Set default for CTEST_BINARY_DIRECTORY
#
IF ("${CTEST_BINARY_DIRECTORY}" STREQUAL "")
  MESSAGE("Set defualt for CTEST_BINARY_DIRECTORY to $PWD/BUILD='$ENV{PWD}/BUILD'")
  SET(CTEST_BINARY_DIRECTORY $ENV{PWD}/BUILD)
ENDIF()

#
# Set CMAKE_MODULE_PATH
#
SET( CMAKE_MODULE_PATH
  "${TRIBITS_PROJECT_ROOT}"
  "${TRIBITS_PROJECT_ROOT}/cmake"
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch"
  "${${PROJECT_NAME}_TRIBITS_DIR}/ci_support"
  )

INCLUDE(TribitsConstants)
TRIBITS_ASESRT_MINIMUM_CMAKE_VERSION()
INCLUDE(TribitsCMakePolicies)

INCLUDE(PrintVar)
INCLUDE(MultilineSet)
INCLUDE(SetDefaultAndFromEnv)
INCLUDE(AssertDefined)
INCLUDE(AppendSet)
INCLUDE(AppendStringVar)
INCLUDE(TribitsGlobalMacros)
INCLUDE(TribitsStripCommentsFromCMakeCacheFile)

INCLUDE(${CMAKE_CURRENT_LIST_DIR}/TribitsUpdateExtraRepo.cmake)

# Need to include the project's version file to get some Git and CDash
# settings specific to the given version
TRIBITS_PROJECT_READ_VERSION_FILE(${TRIBITS_PROJECT_ROOT})

INCLUDE(TribitsFindPythonInterp)
TRIBITS_FIND_PYTHON()
MESSAGE("PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")

#############################
### Do some initial setup ###
#############################


# Get the host type

IF(WIN32)
  SET(HOST_TYPE $ENV{OS})
ELSE()
  FIND_PROGRAM(UNAME_EXE NAMES uname)
  EXECUTE_PROCESS(
    COMMAND ${UNAME_EXE}
    OUTPUT_VARIABLE HOST_TYPE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
ENDIF()


# Find git

FIND_PROGRAM(GIT_EXE NAMES ${GIT_NAME})
MESSAGE("GIT_EXE=${GIT_EXE}")

IF(NOT GIT_EXE)
  QUEUE_ERROR("error: could not find git: GIT_EXE='${GIT_EXE}'")
ENDIF()
IF(NOT EXISTS "${GIT_EXE}")
  QUEUE_ERROR("error: GIT_EXE='${GIT_EXE}' does not exist")
ENDIF()

# Find svn

FIND_PROGRAM(SVN_EXE NAMES svn)
MESSAGE("SVN_EXE=${SVN_EXE}")
# If we don't find svn, no big deal unless we have an SVN repo!

# Get the host name

SITE_NAME(CTEST_SITE_DEFAULT)


########################
### Functions/Macros ###
########################


#
# Wrapper used for unit testing purposes
#

MACRO(EXTRAREPO_EXECUTE_PROCESS_WRAPPER)
  IF (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    EXECUTE_PROCESS(${ARGN}
      RESULT_VARIABLE  EXTRAREPO_EXECUTE_PROCESS_WRAPPER_RTN_VAL)
    IF (NOT EXTRAREPO_EXECUTE_PROCESS_WRAPPER_RTN_VAL STREQUAL "0")
      MESSAGE(SEND_ERROR
        "Error: EXECUTE_PROCESS(${ARGN}) returned"
	" '${EXTRAREPO_EXECUTE_PROCESS_WRAPPER_RTN_VAL}'")
    ENDIF()
  ELSE()
    MESSAGE("EXECUTE_PROCESS(${ARGN})")
  ENDIF()
ENDMACRO()

#
# Wrapper for getting the tracking branch
#
FUNCTION(EXTRAREPO_GET_TRACKING_BRANCH  EXTRAREPO_SRC_DIR  TRACKING_BRANCH_OUT)
  IF (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    EXECUTE_PROCESS(
      COMMAND "${GIT_EXE}" rev-parse --abbrev-ref --symbolic-full-name @{u}
      WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE  EP_RTN
      OUTPUT_VARIABLE  TRACKING_BRANCH
      )
    IF (NOT EP_RTN STREQUAL "0")
      MESSAGE(SEND_ERROR "Error: obtaining tracking branch for repo"
        " '${EXTRAREPO_SRC_DIR}' failed!" )
    ENDIF()
  ELSE()
    SET(TRACKING_BRANCH "tracking/branch")
  ENDIF()
  SET(${TRACKING_BRANCH_OUT}  ${TRACKING_BRANCH}  PARENT_SCOPE)
ENDFUNCTION()


#
# Update or clone a single extra repo
#

FUNCTION(CLONE_OR_UPDATE_EXTRAREPO  EXTRAREPO_NAME_IN  EXTRAREPO_DIR_IN
  EXTRAREPO_REPOTYPE_IN  EXTRAREPO_REPOURL_IN  CHANGES_STR_VAR_INOUT
  )

  #MESSAGE("CLONE_OR_UPDATE_EXTRAREPO: ${EXTRAREPO_NAME_IN} ${EXTRAREPO_REPOURL_IN}")

  SET(CHANGES_STR "${${CHANGES_STR_VAR_INOUT}}")

  SET(EXTRAREPO_SRC_DIR "${${PROJECT_NAME}_SOURCE_DIRECTORY}/${EXTRAREPO_DIR_IN}")
  #PRINT_VAR(EXTRAREPO_SRC_DIR)

  SET(EXTRAREPO_CLONE_OUT_FILE "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.clone.out")
  SET(EXTRAREPO_CHECKOUT_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.checkout.out")

  SET(GIT_LOG_SHORT "${GIT_EXE}" log
    "--pretty=format:%h:  %s%nAuthor: %an <%ae>%nDate:   %ad%n"
    --name-status -C)

  IF (NOT EXISTS "${EXTRAREPO_SRC_DIR}")

    MESSAGE("\n${EXTRAREPO_NAME_IN}: Doing initial ${EXTRAREPO_REPOTYPE_IN}"
      " clone/checkout from URL '${EXTRAREPO_REPOURL_IN}' to dir '${EXTRAREPO_DIR_IN}' ...")

    # Determine the commands to clone/update and find the version
    IF (${EXTRAREPO_REPOTYPE_IN} STREQUAL GIT)
      SET(CLONE_CMND_ARGS
        COMMAND "${GIT_EXE}" clone "${EXTRAREPO_REPOURL}" ${EXTRAREPO_DIR_IN}
        WORKING_DIRECTORY "${${PROJECT_NAME}_SOURCE_DIRECTORY}"
        OUTPUT_FILE "${EXTRAREPO_CLONE_OUT_FILE}" )
      SET(CLONE_VERSION_CMND_ARGS
        COMMAND ${GIT_LOG_SHORT} -1)
    ELSEIF (${EXTRAREPO_REPOTYPE_IN} STREQUAL SVN)
      IF (NOT SVN_EXE)
        MESSAGE(SEND_ERROR "Error, could not find SVN executable!")
      ENDIF()
      SET(CLONE_CMND_ARGS
        COMMAND "${SVN_EXE}" checkout "${EXTRAREPO_REPOURL}" ${EXTRAREPO_DIR_IN}
        WORKING_DIRECTORY "${${PROJECT_NAME}_SOURCE_DIRECTORY}"
        OUTPUT_FILE "${EXTRAREPO_CHECKOUT_OUT_FILE}" )
      SET(CLONE_VERSION_CMND_ARGS "echo") # ToDo: Define this for SVN
    ELSE()
      MESSAGE(SEND_ERROR
	"Error, Invalid EXTRAREPO_REPOTYPE_IN='${EXTRAREPO_REPOTYPE_IN}'!")
    ENDIF()

    # Do the clone/update
    EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${CLONE_CMND_ARGS})

    # Determine the cloned/checkout version
    EXTRAREPO_EXECUTE_PROCESS_WRAPPER(
      ${CLONE_VERSION_CMND_ARGS}
      OUTPUT_VARIABLE  CLONE_OR_PULL_OUTPUT
      WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}")

    SET(VERSION_STATUS_TYPE "cloned version")

  ELSE()

    MESSAGE("\n${EXTRAREPO_NAME_IN}: Doing ${EXTRAREPO_REPOTYPE_IN} update"
      " from URL '${EXTRAREPO_REPOURL_IN}' to dir '${EXTRAREPO_SRC_DIR}' ...")

    # Pull/update changes
    IF (${EXTRAREPO_REPOTYPE_IN} STREQUAL GIT)
      EXTRAREPO_CLEAN_FETCH_RESET("${GIT_EXE}" "${EXTRAREPO_SRC_DIR}")
      SET(PULL_DIFF_CMND_ARGS
        COMMAND ${GIT_LOG_SHORT} HEAD ^ORIG_HEAD)
    ELSEIF (${EXTRAREPO_REPOTYPE_IN} STREQUAL SVN)
      SET(PULL_CMND_ARGS
        COMMAND "${SVN_EXE}" update
        WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
        OUTPUT_FILE "${EXTRAREPO_UPDATE_OUT_FILE}" )
      EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${PULL_CMND_ARGS})
      SET(PULL_DIFF_CMND_ARGS "echo") # ToDo: Determine this for SVN
    ELSE()
      MESSAGE(SEND_ERROR
	"Error, Invalid EXTRAREPO_REPOTYPE_IN='${EXTRAREPO_REPOTYPE_IN}'!")
    ENDIF()

    # Determine what changed
    EXTRAREPO_EXECUTE_PROCESS_WRAPPER(
      ${PULL_DIFF_CMND_ARGS}
      OUTPUT_VARIABLE  CLONE_OR_PULL_OUTPUT
      WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}")

    SET(VERSION_STATUS_TYPE "updates")

  ENDIF()

  SET(CHANGES_STR
    "${CHANGES_STR}\n***\n*** ${EXTRAREPO_NAME_IN} ${VERSION_STATUS_TYPE}:\n***\n\n${CLONE_OR_PULL_OUTPUT}\n\n")

  SET(${CHANGES_STR_VAR_INOUT} "${CHANGES_STR}" PARENT_SCOPE)

ENDFUNCTION()


#
# Clone or update the whole list of extra repositories
#

FUNCTION(CLONE_OR_UPDATE_ALL_EXTRAREPOS  CHANGES_STR_VAR_OUT)

  SET(CHANGES_STR "")

  SET(EXTRAREPO_IDX 0)
  FOREACH(EXTRAREPO_NAME ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES})
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS ${EXTRAREPO_IDX} EXTRAREPO_DIR )
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES ${EXTRAREPO_IDX} EXTRAREPO_REPOTYPE )
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS ${EXTRAREPO_IDX} EXTRAREPO_REPOURL )
    CLONE_OR_UPDATE_EXTRAREPO(${EXTRAREPO_NAME} ${EXTRAREPO_DIR}
      ${EXTRAREPO_REPOTYPE} ${EXTRAREPO_REPOURL} CHANGES_STR)
    MATH(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")
  ENDFOREACH()

  SET(${CHANGES_STR_VAR_OUT} "${CHANGES_STR}" PARENT_SCOPE)

ENDFUNCTION()


#
# Select the set of extra repositories
#

MACRO(TRIBITS_SETUP_EXTRAREPOS)

  IF (EXISTS "${${PROJECT_NAME}_EXTRAREPOS_FILE}" )
    # Repos many not already exist because we have not cloned them yet!
    SET(${PROJECT_NAME}_CHECK_EXTRAREPOS_EXIST FALSE)
    TRIBITS_GET_AND_PROCESS_EXTRA_REPOSITORIES_LISTS()
  ELSE()
    MESSAGE("${${PROJECT_NAME}_EXTRAREPOS_FILE} does not exist,"
       " skipping extra repositories.")
  ENDIF()

ENDMACRO()


#
# Select the list of packages
#
# OUTPUT: Sets ${PROJECT_NAME}_DEFAULT_PACKAGES
#
# NOTE: This macro is used to clean up the main TRIBITS_CTEST_DRIVER()
# macro.
#

MACRO(TRIBITS_SETUP_PACKAGES)

  # Here, we must point into the source tree just cloned (or updated)
  # and not the "driver" source dir tree for two reasons.  First, the
  # list of core packages may be more recent in what was checked out.
  # Second, the extra repos do not even exist in the "driver" source
  # tree.

  SET(${PROJECT_NAME}_ASSERT_MISSING_PACKAGES FALSE)
  SET(${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES FALSE)
  IF (CTEST_GENERATE_OUTER_DEPS_XML_OUTPUT_FILE)
    SET(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE
       "${PROJECT_BINARY_DIR}/${${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME}")
  ELSE()
    SET(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
  ENDIF()
  IF (CTEST_SUBMIT_CDASH_SUBPROJECTS_DEPS_FILE)
    SET(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE
      "${PROJECT_BINARY_DIR}/${${PROJECT_NAME}_CDASH_SUBPROJECT_DEPS_XML_FILE_NAME}" )
  ELSE()
    SET(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE)
  ENDIF()
  SET(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE)

  # Don't ignore missing repos by default.  This will allow processing to
  # continue but this outer CTest script will fail (thereby sending a CDash
  # email from the TDD system).  However, when we configure actual packages,
  # we do set this to TRUE so that the package configures will not fail due to
  # missing extra repositories.
  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES FALSE)
  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_PRE_REPOSITORIES "")
  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_EXTRA_REPOSITORIES "")
  SPLIT("${${PROJECT_NAME}_PRE_REPOSITORIES}"  ","  ${PROJECT_NAME}_PRE_REPOSITORIES)
  SPLIT("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  ","  ${PROJECT_NAME}_EXTRA_REPOSITORIES)

  TRIBITS_READ_IN_NATIVE_REPOSITORIES()
  TRIBITS_COMBINE_NATIVE_AND_EXTRA_REPOS()
  TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()

  # When we get here, we will have the basic dependency structure set up
  # with only defaults set

  # Set this to "" so that it can be defined in ENABLE_MODIFIED_PACKAGES_ONLY()
  SET(${PROJECT_NAME}_ENABLE_ALL_PACKAGES "")

ENDMACRO()


MACRO(ENABLE_PACKAGE_IF_NOT_EXPLICITLY_EXCLUDED  TRIBITS_PACKAGE)
  IF ("${${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}}" STREQUAL "")
    MESSAGE("Enabling explicitly set package ${TRIBITS_PACKAGE} ...")
    SET(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} ON)
  ELSEIF(NOT ${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
    IF (${TRIBITS_PACKAGE}_EXPLICITY_EXCLUDED)
      MESSAGE("NOT enabling explicitly set package ${TRIBITS_PACKAGE} since it was explicitly excluded!")
    ELSE()
       MESSAGE("Enabling explicitly set package ${TRIBITS_PACKAGE} which was default or otherwise disabed!")
      SET(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} ON)
    ENDIF()
  ELSE()
    MESSAGE("Explicitly set package ${TRIBITS_PACKAGE} is already enabled?")
  ENDIF()
ENDMACRO()



#
# Select packages set by the input
#

MACRO(ENABLE_USER_SELECTED_PACKAGES)

  # 1) Set the enables for packages already set with
  # ${PROJECT_NAME}_PACKAGES_USER_SELECTED

  IF (NOT ${PROJECT_NAME}_PACKAGES_USER_SELECTED)
    SET(${PROJECT_NAME}_ENABLE_ALL_PACKAGES ON)
  ELSE()
    FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES_USER_SELECTED})
      ENABLE_PACKAGE_IF_NOT_EXPLICITLY_EXCLUDED(${TRIBITS_PACKAGE})
    ENDFOREACH()
  ENDIF()

  # 2) Set extra package enables from ${PROJECT_NAME}_ADDITIONAL_PACKAGES

  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_ADDITIONAL_PACKAGES})
    ENABLE_PACKAGE_IF_NOT_EXPLICITLY_EXCLUDED(${TRIBITS_PACKAGE})
  ENDFOREACH()

ENDMACRO()


#
# Extract the list of changed files for the main repo on put into an
# modified files file.
#

MACRO(TRIBITS_GET_MODIFIED_FILES  WORKING_DIR_IN  MODIFIED_FILES_FILE_NAME_IN)
  SET(CMND_ARGS
    COMMAND "${GIT_EXE}" diff --name-only ORIG_HEAD..HEAD
    WORKING_DIRECTORY "${WORKING_DIR_IN}"
    OUTPUT_FILE ${MODIFIED_FILES_FILE_NAME_IN}
    #OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  IF (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    EXECUTE_PROCESS(${CMND_ARGS})
  ELSE()
    MESSAGE("EXECUTE_PROCESS(${CMND_ARGS})")
  ENDIF()
ENDMACRO()


#
# Select only packages that are modified or failed in the last CI iteration
#

MACRO(ENABLE_ONLY_MODIFIED_PACKAGES)

  #
  # A) Get the list of changed packages
  #

  SET(MODIFIED_FILES_FILE_NAME "${CTEST_BINARY_DIRECTORY}/modifiedFiles.txt")

  # A.1) Get changes from main ${PROJECT_NAME} repo

  TRIBITS_GET_MODIFIED_FILES("${CTEST_SOURCE_DIRECTORY}" "${MODIFIED_FILES_FILE_NAME}")

  # A.2) Get changes from extra repos

  SET(EXTRAREPO_IDX 0)
  FOREACH(EXTRAREPO_NAME ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES})

    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS
       ${EXTRAREPO_IDX} EXTRAREPO_DIR )
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS
      ${EXTRAREPO_IDX} EXTRAREPO_PACKSTAT )

    # For now, only look for changes if it has packages.  Later, we need to
    # generalize this for the general extra repo case with deeper directory
    # and other than GIT (e.g. SVN with Dakota).  For example, we would like
    # to pick up changes to Dakota and therefore enable TriKota to build.
    IF (EXTRAREPO_PACKSTAT STREQUAL HASPACKAGES)

      SET(EXTRAREPO_SRC_DIR "${CTEST_SOURCE_DIRECTORY}/${EXTRAREPO_DIR}")
      SET(EXTRAREPO_MODIFIED_FILES_FILE_NAME
        "${CTEST_BINARY_DIRECTORY}/modifiedFiles.${EXTRAREPO_NAME}.txt")

      TRIBITS_GET_MODIFIED_FILES("${EXTRAREPO_SRC_DIR}"
        "${EXTRAREPO_MODIFIED_FILES_FILE_NAME}")

      FILE(STRINGS ${EXTRAREPO_MODIFIED_FILES_FILE_NAME} EXTRAREPO_MODIFIED_FILES_STR)
      SET(EXTRAREPO_FILES_STR "")
      FOREACH(STR_LINE ${EXTRAREPO_MODIFIED_FILES_STR})
        APPEND_STRING_VAR(EXTRAREPO_FILES_STR "${EXTRAREPO_DIR}/${STR_LINE}\n")
      ENDFOREACH()
      FILE(APPEND "${MODIFIED_FILES_FILE_NAME}" ${EXTRAREPO_FILES_STR})

    ENDIF()

    MATH(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")

  ENDFOREACH()

  # A.3) Get the names of the modified packages

  IF (NOT PYTHON_EXECUTABLE)
    MESSAGE(FATAL_ERROR "Error, Python must be enabled to map from modified"
      " files to packages!")
  ENDIF()

  IF (EXISTS "${MODIFIED_FILES_FILE_NAME}")
    EXECUTE_PROCESS(
      COMMAND ${PYTHON_EXECUTABLE}
        ${${PROJECT_NAME}_TRIBITS_DIR}/ci_support/get-tribits-packages-from-files-list.py
        --files-list-file=${MODIFIED_FILES_FILE_NAME}
        --deps-xml-file=${CTEST_BINARY_DIRECTORY}/${${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME}
      OUTPUT_VARIABLE MODIFIED_PACKAGES_LIST
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  ELSE()
    SET(MODIFIED_PACKAGES_LIST)
  ENDIF()

  PRINT_VAR(MODIFIED_PACKAGES_LIST)

  #
  # B) Get the list of packages that failed last CI iteration
  #

  # NOTE: It is critical to enable and test packages until they pass.  If you
  # don't do this, then the package will not show as updated in the above
  # logic.  In this case only downstream packages will get enabled.  If the
  # failing packages break the downstream packages, this will be bad (for lots
  # of reasons).  Therefore, we must enable failing packages from the last CI
  # iteration and keep enabling and testing them until they do pass!

  IF (EXISTS "${FAILED_PACKAGES_FILE_NAME}")
    FILE(READ "${FAILED_PACKAGES_FILE_NAME}" FAILING_PACKAGES_LIST)
    STRING(STRIP "${FAILING_PACKAGES_LIST}" FAILING_PACKAGES_LIST)
    PRINT_VAR(FAILING_PACKAGES_LIST)
  ENDIF()

  #
  # C) Enable the changed and previously failing packages
  #

  FOREACH(TRIBITS_PACKAGE ${MODIFIED_PACKAGES_LIST})
    #PRINT_VAR(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
    ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
    IF ("${${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}}" STREQUAL "")
      IF (
        ${TRIBITS_PACKAGE}_TESTGROUP STREQUAL "PT"
        OR
        (
          ${TRIBITS_PACKAGE}_TESTGROUP STREQUAL "ST"
           AND
           ${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE
           )
        )
        MESSAGE("Enabling modified package: ${TRIBITS_PACKAGE}")
        SET(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} ON)
      ELSE()
        MESSAGE("NOT enabling modified ST package: ${TRIBITS_PACKAGE}")
      ENDIF()
    ELSE()
      MESSAGE("Not enabling explicitly disabled modified package: ${TRIBITS_PACKAGE}")
    ENDIF()
  ENDFOREACH()

  FOREACH(TRIBITS_PACKAGE ${FAILING_PACKAGES_LIST})
    IF ("${${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}}" STREQUAL "")
      IF (
        ${TRIBITS_PACKAGE}_TESTGROUP STREQUAL "PT"
        OR
        (
          ${TRIBITS_PACKAGE}_TESTGROUP STREQUAL "ST"
           AND
           ${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE
           )
        )
        MESSAGE("Enabling previously failing package: ${TRIBITS_PACKAGE}")
        SET(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} ON)
      ELSE()
        MESSAGE("NOT enabling previously failing ST package: ${TRIBITS_PACKAGE}")
      ENDIF()
    ELSE()
      MESSAGE("Not enabling explicitly disabled previously"
        " failing package: ${TRIBITS_PACKAGE}")
    ENDIF()
  ENDFOREACH()

  #
  # D) Print the final status
  #

  TRIBITS_PRINT_ENABLED_SE_PACKAGE_LIST(
    "\nDirectly modified or failing non-disabled packages that need to be tested" ON FALSE)

ENDMACRO()


#
# Exclude disabled packages from ${PROJECT_NAME}_EXCLUDE_PACKAGES
#
# NOTE: These disables need to dominate over the above enables so this code is
# after all the enable code has run
#

MACRO(DISABLE_EXCLUDED_PACKAGES)
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_EXCLUDE_PACKAGES})
    MESSAGE("Disabling excluded package ${TRIBITS_PACKAGE} ...")
    SET(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} OFF)
    SET(${TRIBITS_PACKAGE}_EXPLICITY_EXCLUDED TRUE)
  ENDFOREACH()
ENDMACRO()


#
# Remove packages that are only implicitly enabled but don't have tests
# enabled.
#
#
MACRO(SELECT_FINAL_SET_OF_PACKAGES_TO_PROCESS)

  SET(${PROJECT_NAME}_PACKAGES_TO_PROCESS)

  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES})

    SET(PROCESS_THE_PACKAGE FALSE)

    IF (${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}
      AND ${TRIBITS_PACKAGE}_ENABLE_TESTS
      )
      SET(PROCESS_THE_PACKAGE  TRUE)
    ELSEIF (${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}
      AND CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES
      )
      SET(PROCESS_THE_PACKAGE  TRUE)
    ENDIF()

    IF(PROCESS_THE_PACKAGE)
      APPEND_SET(${PROJECT_NAME}_PACKAGES_TO_PROCESS  ${TRIBITS_PACKAGE})
    ENDIF()

  ENDFOREACH()

  SET(${PROJECT_NAME}_PACKAGES ${${PROJECT_NAME}_PACKAGES_TO_PROCESS})

ENDMACRO()


#
# Select the default generator.
#

MACRO(SELECT_DEFAULT_GENERATOR)
  # When the build tree is known and exists, use
  # its generator.
  SET(DEFAULT_GENERATOR "DID NOT SET!")
  IF(EXISTS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
    FILE(STRINGS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt"
      line REGEX "^CMAKE_GENERATOR:" LIMIT_COUNT 1)
    IF("${line}" MATCHES "=(.+)$")
      SET(DEFAULT_GENERATOR "${CMAKE_MATCH_1}")
    ENDIF()
  ELSE()
    SET(DEFAULT_GENERATOR "Unix Makefiles")
  ENDIF()
ENDMACRO()


#
# Call INITIALIZE_ERROR_QUEUE once at the top of TRIBITS_CTEST_DRIVER
#

MACRO(INITIALIZE_ERROR_QUEUE)
  SET(TRIBITS_CTEST_DRIVER_ERROR_QUEUE "")
ENDMACRO()


#
# QUEUE_ERROR should be called only for errors that are not already reported to
# the dashboard in some other way. For example, if calling ctest_submit fails,
# then that failure does NOT show up on the dashboard, so it is appropriate to
# call QUEUE_ERROR for that case. For a build error or test failure, it is NOT
# appropriate to call QUEUE_ERROR because those already show up on the
# dashboard (assuming a good ctest_submit...)
#
# When adding more callers of QUEUE_ERROR, just make sure that it does not
# duplicate an existing/reported dashboard failure.
#

MACRO(QUEUE_ERROR err_msg)
  SET(TRIBITS_CTEST_DRIVER_ERROR_QUEUE
    ${TRIBITS_CTEST_DRIVER_ERROR_QUEUE} "${err_msg}")
ENDMACRO()


#
# Call REPORT_QUEUED_ERRORS once at the bottom of TRIBITS_CTEST_DRIVER
#

MACRO(REPORT_QUEUED_ERRORS)
  IF("${TRIBITS_CTEST_DRIVER_ERROR_QUEUE}" STREQUAL "")
    MESSAGE("TRIBITS_CTEST_DRIVER_ERROR_QUEUE is empty. All is well.")
  ELSE()
    MESSAGE("error: TRIBITS_CTEST_DRIVER_ERROR_QUEUE reports the following error message queue:")
    FOREACH(err_msg ${TRIBITS_CTEST_DRIVER_ERROR_QUEUE})
      MESSAGE("${err_msg}")
    ENDFOREACH()
  ENDIF()
ENDMACRO()


#
# Override CTEST_SUBMIT to drive multiple submits and to detect failed
# submissions and track them as queued errors.
#

MACRO(TRIBITS_CTEST_SUBMIT)

  # Cache the original CTEST_DROP_SITE and CTEST_DROP_LOCATION
  IF ("${TRIBITS_CTEST_DROP_SITE_ORIG}" STREQUAL "")
    SET(TRIBITS_CTEST_DROP_SITE_ORIG ${CTEST_DROP_SITE})
    IF (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
      PRINT_VAR(TRIBITS_CTEST_DROP_SITE_ORIG)
    ENDIF()
  ENDIF()
  IF ("${TRIBITS_CTEST_DROP_LOCATION_ORIG}" STREQUAL "")
    SET(TRIBITS_CTEST_DROP_LOCATION_ORIG ${CTEST_DROP_LOCATION})
    IF (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
      PRINT_VAR(TRIBITS_CTEST_DROP_LOCATION_ORIG)
    ENDIF()
  ENDIF()

  # Do the first submit
  SET(CTEST_DROP_SITE ${TRIBITS_CTEST_DROP_SITE_ORIG})
  SET(CTEST_DROP_LOCATION ${TRIBITS_CTEST_DROP_LOCATION_ORIG})
  IF (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
    PRINT_VAR(CTEST_DROP_SITE)
    PRINT_VAR(CTEST_DROP_LOCATION)
  ENDIF()

  TRIBITS_CTEST_SUBMIT_DRIVER(${ARGN})

  # Do the second submit if requested!
  IF (TRIBITS_2ND_CTEST_DROP_SITE OR TRIBITS_2ND_CTEST_DROP_LOCATION)

    MESSAGE("\nDoing submit to second CDash site ...\n")

    IF (NOT "${TRIBITS_2ND_CTEST_DROP_SITE}" STREQUAL "")
      IF (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
        PRINT_VAR(TRIBITS_2ND_CTEST_DROP_SITE)
      ENDIF()
      SET(CTEST_DROP_SITE ${TRIBITS_2ND_CTEST_DROP_SITE})
    ENDIF()

    IF (NOT "${TRIBITS_2ND_CTEST_DROP_LOCATION}" STREQUAL "")
      IF (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
        PRINT_VAR(TRIBITS_2ND_CTEST_DROP_LOCATION)
      ENDIF()
      SET(CTEST_DROP_LOCATION ${TRIBITS_2ND_CTEST_DROP_LOCATION})
    ENDIF()

    TRIBITS_CTEST_SUBMIT_DRIVER(${ARGN})

  ENDIF()

ENDMACRO()


MACRO(TRIBITS_CTEST_SUBMIT_DRIVER)

  # If using a recent enough ctest with RETRY_COUNT, use it to overcome
  # failed submits:
  SET(retry_args "")
  SET(retry_args RETRY_COUNT 25 RETRY_DELAY 120)
  MESSAGE("info: using retry_args='${retry_args}' for _ctest_submit call")

  # Call the original CTEST_SUBMIT and pay attention to its RETURN_VALUE:
  CTEST_SUBMIT(${ARGN} ${retry_args} RETURN_VALUE rv)

  IF(NOT "${rv}" STREQUAL "0")
    QUEUE_ERROR("error: ctest_submit failed: rv='${rv}' ARGN='${ARGN}' retry_args='${retry_args}'")
  ENDIF()

ENDMACRO()


#
# Wrapper for CTEST_UPDATE(...) for unit testing
#

MACRO(CTEST_UPDATE_WRAPPER)
  IF (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    CTEST_UPDATE(${ARGN})
  ELSE()
    MESSAGE("CTEST_UPDATE(${ARGN})")
    SET(UPDATE_RETURN_VAL ${CTEST_UPDATE_RETURN_VAL})
  ENDIF()
ENDMACRO()


#
# @FUNCTION: TRIBITS_CTEST_DRIVER()
#
# Platform-independent package-by-package CTest/CDash driver (run by ``ctest``
# **NOT** ``cmake``).
#
# Usage::
#
#   TRIBITS_CTEST_DRIVER()
#
# This driver code that is platform independent.  This script drives the
# testing process by doing a version control (VC) source update on all of the
# VC repos and then configuring and building the top-level TriBITS packages
# one at a time, in order.  This function gets called from inside of a
# platform and build-specific ``ctest -S`` driver script.
#
# To understand this script, one must understand that it gets run in several
# different modes:
#
# **Mode 1**: Run where there are already existing source and binary
# directories (``CTEST_DASHBOARD_ROOT`` is set empty before call).  This is
# for when the ctest driver script is run on an existing source and binary
# tree.  In this case, there is one project source tree and
# ``CTEST_SOURCE_DIRECTORY`` and ``CTEST_BINARY_DIRECTORY`` must be set by the
# user before calling this function.  This is used to test a local build and
# post to CDash.
#
# **Mode 2**: A new binary directory is created and new sources are cloned (or
# updated) in a driver directory (``CTEST_DASHBOARD_ROOT`` is set before
# call).  In this case, there are always two (partial) project source tree's,
# i) a "driver" skeleton source tree (typically embedded with TriBITS
# directory) that bootstraps the testing process, and ii) a true full "source"
# that is (optionally) cloned and/or updated.
#
# There are a few different directory locations are significant for this
# script:
#
#   ``TRIBITS_PROJECT_ROOT``
#
#     The root directory to an existing source tree where the project's
#     `<projectDir>/ProjectName.cmake`_ (defining ``PROJECT_NAME`` variable)
#     and ``Version.cmake`` file's can be found.
#
#   ``${PROJECT_NAME}_TRIBITS_DIR``
#
#     The base directory for the TriBITS system's various CMake modules,
#     python scripts, and other files.  By default this is assumed to be in
#     the source tree under ``${TRIBITS_PROJECT_ROOT}`` (see below) but it can
#     be overridden to point to any location.
#
#   ``CTEST_DASHBOARD_ROOT``
#
#     If set, this is the base directory where this script runs that clones
#     the sources for the project.  If this directory does not exist, it will
#     be created.  If empty, then has no effect on the script.
#
#   ``CTEST_SOURCE_DIRECTORY``
#
#     Determines the location of the sources that are used to define packages,
#     dependencies and configure and build the software.  This is a variable
#     that CTest directly reads and must therefore be set. This is used to set
#     `PROJECT_SOURCE_DIR`_ which is used by the TriBITS system.  If
#     ``CTEST_DASHBOARD_ROOT`` is set, then this is hard-coded internally to
#     ``${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}``.
#
#   ``CTEST_BINARY_DIRECTORY``
#
#     Determines the location of the binary tree where output from CMake/CTest
#     is put.  This is used to set to `PROJECT_BINARY_DIR`_ which is used by
#     the TriBITS system.  If ``CTEST_DASHBOARD_ROOT`` is set, then this is
#     hard-coded internally to ``${CTEST_DASHBOARD_ROOT}/BUILD``.
#
# ToDo: Document input variables that have defaults, to be set before, and can
# be overridden from the env.
#
# ToDo: Finish Documentation!
#
FUNCTION(TRIBITS_CTEST_DRIVER)

  MESSAGE("")
  MESSAGE("******************************")
  MESSAGE("*** TRIBITS_CTEST_DRIVER() ***")
  MESSAGE("******************************")
  MESSAGE("")

  INITIALIZE_ERROR_QUEUE()

  # The name of the source directory. Defaults to project name, but
  # can be overridden by the environment for cases in which the source
  # directory name does not match the project name.
  SET_DEFAULT_AND_FROM_ENV(CTEST_SOURCE_NAME ${PROJECT_NAME})

  MESSAGE(
    "\n***"
    "\n*** Setting input options to default and reading from env ..."
    "\n***\n")

  SET_DEFAULT_AND_FROM_ENV( CTEST_CONFIGURATION_UNIT_TESTING OFF )

  # The type of test (e.g. Nightly, Experimental, Continuous)
  SET_DEFAULT_AND_FROM_ENV( CTEST_TEST_TYPE Experimental )

  # The default track to send the build to. This can be changed to send
  # the data to a different nightly grouping on the dashboard.
  # If the test type is set to Experimental though the track is forced
  # to "Experimental" this is so that we can have experimental tests
  # on branches.
  IF(${PROJECT_NAME}_TESTING_TRACK)
    SET(${PROJECT_NAME}_TRACK_DEFAULT ${${PROJECT_NAME}_TESTING_TRACK})
  ELSE()
    SET(${PROJECT_NAME}_TRACK_DEFAULT "")
  ENDIF()
  print_var(${PROJECT_NAME}_TRACK_DEFAULT)
  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_TRACK "${${PROJECT_NAME}_TRACK_DEFAULT}")
  IF(CTEST_TEST_TYPE STREQUAL "Experimental" OR CTEST_TEST_TYPE STREQUAL "EXPERIMENTAL")
    SET(${PROJECT_NAME}_TRACK "Experimental")
    MESSAGE("-- Test type is Experimental. Forcing ${PROJECT_NAME}_TRACK to Experimental")
    PRINT_VAR(${PROJECT_NAME}_TRACK)
  ENDIF()

  # The name of the site in the dashboard (almost never need to override this)
  SET_DEFAULT_AND_FROM_ENV( CTEST_SITE ${CTEST_SITE_DEFAULT} )

  # The root of the dasbhoard where ${PROJECT_NAME} will be cloned and the
  # BUILD directory will be create (only override for separate testing)
  SET_DEFAULT_AND_FROM_ENV( CTEST_DASHBOARD_ROOT "" )

  # The build type (e.g. DEBUG, RELEASE, NONE)
  SET_DEFAULT_AND_FROM_ENV( BUILD_TYPE NONE )

  # Verobse configure or now
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_VERBOSE_CONFIGURE OFF )

  # Set the default compiler version
  SET_DEFAULT_AND_FROM_ENV(COMPILER_VERSION UNKNOWN)

  # The name of the build that appears in the dashbaord
  SET_DEFAULT_AND_FROM_ENV( CTEST_BUILD_NAME
    "${HOST_TYPE}-${COMPILER_VERSION}-${BUILD_DIR_NAME}" )

  # Remove the entire build directory if it exists or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE )

  # Remove an existing CMakeCache.txt file or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_WIPE_CACHE TRUE )

  # Select a default generator.
  SELECT_DEFAULT_GENERATOR()
  SET_DEFAULT_AND_FROM_ENV( CTEST_CMAKE_GENERATOR ${DEFAULT_GENERATOR})

  # Do the Git updates or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_UPDATES TRUE )

  # Generate the XML dependency output files or not in the inner CMake
  # configure.  There is really no reason to do this.  This option is
  # maintained for backward compatibility.
  SET_DEFAULT_AND_FROM_ENV( CTEST_GENERATE_DEPS_XML_OUTPUT_FILE FALSE )

  # Flags used on git when doing a Git update
  SET_DEFAULT_AND_FROM_ENV( CTEST_UPDATE_ARGS "")

  # Flags used on update when doing a Git update
  SET_DEFAULT_AND_FROM_ENV( CTEST_UPDATE_OPTIONS "")

  # Flags passed to 'make' assume gnumake with unix makefiles
  IF("${CTEST_CMAKE_GENERATOR}" MATCHES "Unix Makefiles")
    SET_DEFAULT_AND_FROM_ENV( CTEST_BUILD_FLAGS "-j2")
  ELSE()
    SET_DEFAULT_AND_FROM_ENV( CTEST_BUILD_FLAGS "")
  ENDIF()

  # Do the build or use an existing build
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_BUILD TRUE )

  # Do the tests or not (Note: must be true for coverage testing)
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_TEST TRUE )

  # Maximum number of procs an mpi test can request (if more are requested,
  # the test will be skipped).  Value of 0 means no override (determined
  # internally).
  SET_DEFAULT_AND_FROM_ENV( MPI_EXEC_MAX_NUMPROCS 0 )

  # How many tests ctest will spawn simultaneously
  SET_DEFAULT_AND_FROM_ENV( CTEST_PARALLEL_LEVEL 1 )

  # Turn off or change warnings-as-errors flag(s) (i.e. -Werror)
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS "" )

  # Do coverage testing or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_COVERAGE_TESTING FALSE )

  # Command to run to get coverage results
  SET_DEFAULT_AND_FROM_ENV( CTEST_COVERAGE_COMMAND gcov )

  # Do memory testing (i.e. valgrind) or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_MEMORY_TESTING FALSE )

  # Command used to perform the memory testing (i.e. valgrind)
  FIND_PROGRAM(VALGRIND_EXE NAMES valgrind)
  PRINT_VAR(VALGRIND_EXE)
  IF (VALGRIND_EXE)
    SET(CTEST_MEMORYCHECK_COMMAND_DEFAULT "${VALGRIND_EXE}")
  ELSE()
    SET(CTEST_MEMORYCHECK_COMMAND_DEFAULT)
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV( CTEST_MEMORYCHECK_COMMAND "${CTEST_MEMORYCHECK_COMMAND_DEFAULT}" )

  # Set the default options
  SET_DEFAULT_AND_FROM_ENV( CTEST_MEMORYCHECK_COMMAND_OPTIONS "")

  # Generate the basic package dependencies XML file in the outer CTest
  # program.  This XML file is used to match up modified files with with
  # changed TriBITS packages.  This file only needs to be generated in CI
  # iterations and is not needed in Nightly testing.  Turning off its
  # generation can also speed up local manual testing for large projects with
  # lots of TriBITS packges.
  SET_DEFAULT_AND_FROM_ENV( CTEST_GENERATE_OUTER_DEPS_XML_OUTPUT_FILE TRUE )

  # Generate and submit the CDash subprojects XML file
  SET_DEFAULT_AND_FROM_ENV( CTEST_SUBMIT_CDASH_SUBPROJECTS_DEPS_FILE TRUE )

  # Submit the results to the dashboard or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_SUBMIT TRUE )

  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE OFF )

  # List of additional packges that will be enabled over the current set
  # of all packagess (that would be set by ${PROJECT_NAME}_ENABLE_ALL_PACKAGES).
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_ADDITIONAL_PACKAGES "" )

  # List of packages to not directly process.  NOTE: Listing these packages
  # here will *not* disable the package in the CMake build system.  To do
  # that, you will have to disable them in the variable
  # EXTRA_CONFIGURE_OPTIONS (set in your driver script.
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_EXCLUDE_PACKAGES "" )

  IF(${PROJECT_NAME}_REPOSITORY_BRANCH)
    SET(${PROJECT_NAME}_BRANCH_DEFAULT ${${PROJECT_NAME}_REPOSITORY_BRANCH})
  ELSE()
    SET(${PROJECT_NAME}_BRANCH_DEFAULT "")
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_BRANCH "${${PROJECT_NAME}_BRANCH_DEFAULT}" )

  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE "${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT}" )

  IF(CTEST_TEST_TYPE STREQUAL "Nightly")
    SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_REPOSITORY_LOCATION
       "${${PROJECT_NAME}_REPOSITORY_LOCATION_NIGHTLY_DEFAULT}")
  ELSE()
    SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_REPOSITORY_LOCATION
       "${${PROJECT_NAME}_REPOSITORY_LOCATION_DEFAULT}")
  ENDIF()

  # Select the ${PROJECT_NAME} packages to enable (empty means to select all
  # available).  This will override any disabled packages but not those
  # disabled by ${PROJECT_NAME}_EXCLUDE_PACKAGES.
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_PACKAGES "" )
  SET(${PROJECT_NAME}_PACKAGES_USER_SELECTED ${${PROJECT_NAME}_PACKAGES})
  SPLIT("${${PROJECT_NAME}_PACKAGES_USER_SELECTED}" ","
    ${PROJECT_NAME}_PACKAGES_USER_SELECTED)
  SET(${PROJECT_NAME}_PACKAGES "")
  # Note: above, we have to keep the name ${PROJECT_NAME}_PACKAGES to maintain
  # backward compatibility of this CTest script but we want to let
  # ${PROJECT_NAME}_PACKAGES always be the full set of packages as defined by
  # the basic readin process.

  # Set the file that the extra repos will be read from
  #
  # NOTE: Here, we have no choice but to point into the "driver"
  # ${PROJECT_NAME} source treee because the local ${PROJECT_NAME} sources have not
  # even been checked out yet!  Unless, of course, we are unit testing
  # in which case we will use whatever has been passed in.

  IF (${PROJECT_NAME}_SKIP_EXTRAREPOS_FILE)
    SET(${PROJECT_NAME}_EXTRAREPOS_FILE_DEFAULT)
  ELSE()
    SET(${PROJECT_NAME}_EXTRAREPOS_FILE_DEFAULT
      "${TRIBITS_PROJECT_ROOT}/cmake/${${PROJECT_NAME}_EXTRA_EXTERNAL_REPOS_FILE_NAME}")
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_EXTRAREPOS_FILE
    "${${PROJECT_NAME}_EXTRAREPOS_FILE_DEFAULT}")

  # Select the set of extra external repos to add in packages.
  # These are the same types as CTEST_TEST_TYPE (e.g. 'Continuous' and
  # 'Nightly').  This is set by default to ${CTEST_TEST_TYPE} can be
  # overridden independent of ${CTEST_TEST_TYPE} also.
  #
  # If in release mode generally we do not want any external repositories
  # even though the CTEST_TEST_TYPE is set to "Nightly" for most release
  # builds.
  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
  IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
    SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_DEFAULT ${CTEST_TEST_TYPE})
  ELSE()
    SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_DEFAULT "None")
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE
     "${${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_DEFAULT}" )

  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_PRE_REPOSITORIES "")
  SPLIT("${${PROJECT_NAME}_PRE_REPOSITORIES}"  "," ${PROJECT_NAME}_PRE_REPOSITORIES)

  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_EXTRA_REPOSITORIES "")
  SPLIT("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  "," ${PROJECT_NAME}_EXTRA_REPOSITORIES)

  # Set as part of CI testing in order to only enable modified packages
  SET_DEFAULT_AND_FROM_ENV( CTEST_ENABLE_MODIFIED_PACKAGES_ONLY OFF )

  # Set if implicitly enabled packages should be explicitly processes
  IF (CTEST_ENABLE_MODIFIED_PACKAGES_ONLY AND NOT CTEST_START_WITH_EMPTY_BINARY_DIRECTORY)
    SET( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES_DEFAULT FALSE )
  ELSE()
    SET( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES_DEFAULT TRUE )
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES
    ${CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES_DEFAULT})

  # Set if we should disable enabled fwd packages based on disabled required deps.
  # To make testing robust, we need to do this.
  IF ("${${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT ON)
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES
    ${${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT})

  # Set second drop site and location
  SET_DEFAULT_AND_FROM_ENV( TRIBITS_2ND_CTEST_DROP_SITE "" )
  SET_DEFAULT_AND_FROM_ENV( TRIBITS_2ND_CTEST_DROP_LOCATION "" )

  MESSAGE(
    "\n***"
    "\n*** Setting unit testing input options to default and reading from env ..."
    "\n***\n")

  SET_DEFAULT_AND_FROM_ENV( CTEST_DEPENDENCY_HANDLING_UNIT_TESTING FALSE )

  SET_DEFAULT_AND_FROM_ENV( CTEST_UPDATE_RETURN_VAL 0 )

  IF (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    SET(GIT_EXE /somebasedir/git)
    SET(SVN_EXE /someotherbasedir/svn)
  ENDIF()


  MESSAGE(
    "\n***"
    "\n*** Misc setup ..."
    "\n***\n")

  #
  # Setup and create the base dashboard directory if it is not created yet.
  #

  # NOTE: This is only used in general testing dashbaoard mode, not in local
  # experimental testing mode.

  IF (CTEST_DASHBOARD_ROOT)
    SET( CTEST_BINARY_NAME BUILD )
    SET( CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
    SET( CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")
    IF (NOT EXISTS "${CTEST_DASHBOARD_ROOT}")
      MESSAGE("Creating the dashboard root directory \"${CTEST_DASHBOARD_ROOT}\" ...")
      FILE(MAKE_DIRECTORY "${CTEST_DASHBOARD_ROOT}")
    ENDIF()
  ENDIF()
  PRINT_VAR(CTEST_SOURCE_DIRECTORY)
  PRINT_VAR(CTEST_BINARY_DIRECTORY)

  SET(PROJECT_SOURCE_DIR "${CTEST_SOURCE_DIRECTORY}")
  SET(PROJECT_BINARY_DIR "${CTEST_BINARY_DIRECTORY}")
  PRINT_VAR(PROJECT_SOURCE_DIR)
  PRINT_VAR(PROJECT_BINARY_DIR)

  # Set override hook for unit testing
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_SOURCE_DIRECTORY ${CTEST_SOURCE_DIRECTORY} )

  # Must be set here after CTEST_BINARY_DIRECTORY is set!
  SET(FAILED_PACKAGES_FILE_NAME "${CTEST_BINARY_DIRECTORY}/failedPackages.txt")

  #
  # Some platform-independent setup
  #

  INCLUDE("${TRIBITS_PROJECT_ROOT}/CTestConfig.cmake")
  SET(CMAKE_CACHE_CLEAN_FILE "${CTEST_BINARY_DIRECTORY}/CMakeCache.clean.txt")
  SET(CTEST_USE_LAUNCHERS 1)

  # For coverage dashboards, send results to specialized dashboard if
  # requested
  IF (CTEST_DO_COVERAGE_TESTING)
    # Allow override of CDash drop site but use standard by default
    SET_DEFAULT(CTEST_DROP_SITE_COVERAGE_DEFAULT
      ${CTEST_DROP_SITE})
    SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_SITE_COVERAGE
      "${CTEST_DROP_SITE_COVERAGE_DEFAULT}")
    SET(CTEST_DROP_SITE "${CTEST_DROP_SITE_COVERAGE}" )
    # Allow override of CDash drop location but use standard by default
    SET_DEFAULT(CTEST_DROP_LOCATION_COVERAGE_DEFAULT
      ${CTEST_DROP_LOCATION})
    SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_LOCATION_COVERAGE
      "${CTEST_DROP_LOCATION_COVERAGE_DEFAULT}")
    SET(CTEST_DROP_LOCATION "${CTEST_DROP_LOCATION_COVERAGE}" )

    # NOTE: You must set these down below the include of
    # CTestConfig.cmake so that CTEST_DROP_SITE and CTEST_DROP_LOCATION read
    # from that file will set the defaults for the coverage options.

  ENDIF()

  #
  # Setup for the VC update
  #

  IF (CTEST_DO_UPDATES)

    SET(UPDATE_TYPE "git")
    MESSAGE("UPDATE_TYPE = '${UPDATE_TYPE}'")

    SET(CTEST_UPDATE_COMMAND "${GIT_EXE}")
    MESSAGE("CTEST_UPDATE_COMMAND='${CTEST_UPDATE_COMMAND}'")
    IF(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
      MESSAGE("${CTEST_SOURCE_DIRECTORY} does not exist so setting up for an initial checkout")
      SET( CTEST_CHECKOUT_COMMAND
        "\"${GIT_EXE}\" clone ${CTEST_UPDATE_ARGS} ${${PROJECT_NAME}_REPOSITORY_LOCATION}" )
      MESSAGE("CTEST_CHECKOUT_COMMAND='${CTEST_CHECKOUT_COMMAND}'")
    ELSE()
      MESSAGE("${CTEST_SOURCE_DIRECTORY} exists so skipping the initial checkout.")
    ENDIF()
  ENDIF()

  #
  # Write a few variables to the global level to make cmake happy.
  #
  # If these don't get set in the base CTest script scope, CTest returns an
  # error!
  #

  SET(CTEST_SOURCE_DIRECTORY ${CTEST_SOURCE_DIRECTORY} CACHE INTERNAL "")
  SET(CTEST_BINARY_DIRECTORY ${CTEST_BINARY_DIRECTORY} CACHE INTERNAL "")
  IF ("${CTEST_COMMAND}" STREQUAL "")
    SET(CTEST_COMMAND ctest)
  ENDIF()
  SET(CTEST_COMMAND ${CTEST_COMMAND} CACHE INTERNAL "")

  #
  # Empty out the binary directory
  #

  IF (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY)
    MESSAGE("\nCleaning out binary directory '${CTEST_BINARY_DIRECTORY}' ...")
    CTEST_EMPTY_BINARY_DIRECTORY("${CTEST_BINARY_DIRECTORY}")
  ENDIF()
  # NOTE: The above command will *not* delete the build directory unless there
  # is a CMakeLists.txt file in this directory.  I think Kitware put in this
  # check to avoid accidentally deleting the wrong directory by accident.
  # Also note that you have to delete the build directory before any commands
  # are run that would write files to them (and many of the steps in this
  # process do write files to the binary directory other than just CMake).

  MESSAGE(
    "\n***"
    "\n*** Read in the set of extra repos ..."
    "\n***\n")

  TRIBITS_SETUP_EXTRAREPOS()

  # NOTE: You have to set up the set of extra repos before you can read the
  # Dependencies.cmake files since the extra repos must be cloned first.


  MESSAGE(
    "\n***"
    "\n*** Start up a new dashboard ..."
    "\n***\n")

  PRINT_VAR(CTEST_TEST_TYPE)
  PRINT_VAR(${PROJECT_NAME}_TRACK)

  IF(${PROJECT_NAME}_TRACK)
    CTEST_START(${CTEST_TEST_TYPE} TRACK ${${PROJECT_NAME}_TRACK})
  ELSE()
    CTEST_START(${CTEST_TEST_TYPE})
  ENDIF()


  MESSAGE(
    "\n***"
    "\n*** Update the source code repositories ..."
    "\n***\n")

  SET(UPDATE_FAILED TRUE)

  IF (CTEST_DO_UPDATES)

    IF(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
      QUEUE_ERROR("error: source directory does not exist just prior to CTEST_UPDATE call -- initial checkout did not work")
      REPORT_QUEUED_ERRORS()
      RETURN()
    ENDIF()

    MESSAGE("\nDoing GIT update of '${CTEST_SOURCE_DIRECTORY}' ...")
    CTEST_UPDATE_WRAPPER( SOURCE "${CTEST_SOURCE_DIRECTORY}"
      RETURN_VALUE  UPDATE_RETURN_VAL)
    MESSAGE("CTEST_UPDATE(...) returned '${UPDATE_RETURN_VAL}'")

    CLONE_OR_UPDATE_ALL_EXTRAREPOS(EXTRA_REPO_CHANGES_STR)

    IF (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
      PRINT_VAR(EXTRA_REPO_CHANGES_STR)
    ENDIF()

    FILE(WRITE "${CTEST_BINARY_DIRECTORY}/Updates.txt"
     ${EXTRA_REPO_CHANGES_STR})
    # NOTE: Above, we are making the 'Updates.txt' file come first because
    # this will be the first Notes file shown on CDash.

    # Setting branch switch to success in case we are not doing a switch to a
    # different branch.

    SET(GIT_CHECKOUT_RETURN_VAL "0")

    IF(${PROJECT_NAME}_BRANCH AND NOT "${UPDATE_RETURN_VAL}" LESS "0" AND NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)

      MESSAGE("Doing switch to branch ${${PROJECT_NAME}_BRANCH}")

      EXECUTE_PROCESS(COMMAND ${GIT_EXE} checkout ${${PROJECT_NAME}_BRANCH}
        WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
        RESULT_VARIABLE GIT_CHECKOUT_RETURN_VAL
        OUTPUT_VARIABLE BRANCH_OUTPUT
        ERROR_VARIABLE  BRANCH_ERROR
      )

      IF(NOT "${GIT_CHECKOUT_RETURN_VAL}" EQUAL "0")
        EXECUTE_PROCESS(COMMAND ${GIT_EXE} checkout --track origin/${${PROJECT_NAME}_BRANCH}
          WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
          RESULT_VARIABLE GIT_CHECKOUT_RETURN_VAL
          OUTPUT_VARIABLE BRANCH_OUTPUT
          ERROR_VARIABLE  BRANCH_ERROR
        )
      ENDIF()

      IF(NOT "${GIT_CHECKOUT_RETURN_VAL}" EQUAL "0")
        MESSAGE("Switch to branch ${${PROJECT_NAME}_BRANCH} failed with"
          " error code ${GIT_CHECKOUT_RETURN_VAL}")
        QUEUE_ERROR("Switch to branch ${${PROJECT_NAME}_BRANCH} failed with"
          " error code ${GIT_CHECKOUT_RETURN_VAL}")
      ENDIF()
      #Apparently the successful branch switch is also written to stderr.
      MESSAGE("${BRANCH_ERROR}")

    ENDIF()

    IF ("${UPDATE_RETURN_VAL}" LESS "0" OR NOT "${GIT_CHECKOUT_RETURN_VAL}" EQUAL "0")
      SET(UPDATE_FAILED TRUE)
    ELSE()
      SET(UPDATE_FAILED FALSE)
    ENDIF()

  ELSE()

     MESSAGE("Skipping the update by request!")

     SET(UPDATE_FAILED FALSE)

  ENDIF()


  MESSAGE(
    "\n***"
    "\n*** Read in the set of packages and their dependencies ..."
    "\n***\n")

  # NOTE: You must read the Dependencies.cmake files *after* you have cloned
  # (or updated) all of the code!

  TRIBITS_SETUP_PACKAGES()

  SET(CDASH_SUBPROJECT_XML_FILE
    "${CTEST_BINARY_DIRECTORY}/${${PROJECT_NAME}_CDASH_SUBPROJECT_DEPS_XML_FILE_NAME}")
  PRINT_VAR(CDASH_SUBPROJECT_XML_FILE)

  MESSAGE(
    "\n***"
    "\n*** Disabling packages based on what was set in ${PROJECT_NAME}_EXCLUDE_PACKAGES ..."
    "\n***\n")

  DISABLE_EXCLUDED_PACKAGES()

  IF (NOT CTEST_ENABLE_MODIFIED_PACKAGES_ONLY)
    MESSAGE(
      "\n***"
      "\n*** Determining what packages to enable based what was set in ${PROJECT_NAME}_PACKAGES ..."
      "\n***\n")
    ENABLE_USER_SELECTED_PACKAGES()
  ELSE()
    MESSAGE(
      "\n***"
      "\n*** Determining what packages to enable based on what changed ..."
      "\n***\n")
    ENABLE_ONLY_MODIFIED_PACKAGES()
    SET(${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES ON)
  ENDIF()


  MESSAGE(
    "\n***"
    "\n*** Adjust the package dependencies to enable upstream and"
    " (optionally) downstream packages ..."
    "\n***"
    )

  SET(${PROJECT_NAME}_ENABLE_TESTS ON)
  SET(${PROJECT_NAME}_ENABLE_EXAMPLES ON)
  SET(${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES ON)
  SET(DO_PROCESS_MPI_ENABLES FALSE) # Should not be needed but CMake is messing up
  TRIBITS_ADJUST_AND_PRINT_PACKAGE_DEPENDENCIES() # Sets ${PROJECT_NAME}_NUM_ENABLED_PACKAGES

  SELECT_FINAL_SET_OF_PACKAGES_TO_PROCESS()

  TRIBITS_PRINT_ENABLED_PACKAGE_LIST(
    "\nFinal set of packages to be explicitly processed by CTest/CDash" ON FALSE)

  MESSAGE(
    "\n***"
    "\n*** Determine if to go ahead with configure, build, test ..."
    "\n***")

  IF (CTEST_ENABLE_MODIFIED_PACKAGES_ONLY
    AND ${PROJECT_NAME}_NUM_ENABLED_PACKAGES GREATER 0
    AND MODIFIED_PACKAGES_LIST
    )
    MESSAGE("\nMODIFIED_PACKAGES_LIST='${MODIFIED_PACKAGES_LIST}'"
      ":  Found modified packages, processing enabled packages!\n")
  ELSE()
    MESSAGE(
      "\nCTEST_ENABLE_MODIFIED_PACKAGES_ONLY=${CTEST_ENABLE_MODIFIED_PACKAGES_ONLY}"
      "  Running in regular mode, processing all enabled packages!\n")
  ENDIF()

  IF (${PROJECT_NAME}_NUM_ENABLED_PACKAGES GREATER 0)
    MESSAGE(
      "\n${PROJECT_NAME}_NUM_ENABLED_PACKAGES=${${PROJECT_NAME}_NUM_ENABLED_PACKAGES}:"
      "  Configuring packages!\n")
  ELSE()
    MESSAGE(
      "\n${PROJECT_NAME}_NUM_ENABLED_PACKAGES=${${PROJECT_NAME}_NUM_ENABLED_PACKAGES}:"
      "  Exiting the script!\n")
    REPORT_QUEUED_ERRORS()
    RETURN()
  ENDIF()

  #
  # Delete the CMakeCache.txt file and the CMakeFiles directory for a clean
  # reconfigure.
  #

  IF (CTEST_WIPE_CACHE)
    SET(CACHE_FILE_NAME "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
    IF (EXISTS "${CACHE_FILE_NAME}")
      MESSAGE("Removing existing cache file '${CACHE_FILE_NAME}' ...")
      FILE(REMOVE "${CACHE_FILE_NAME}")
    ENDIF()
    SET(CMAKE_FILES_DIR "${CTEST_BINARY_DIRECTORY}/CMakeFiles/")
    IF (EXISTS "${CMAKE_FILES_DIR}")
      MESSAGE("Removing existing '${CMAKE_FILES_DIR}' ...")
      FILE(REMOVE_RECURSE "${CMAKE_FILES_DIR}")
    ENDIF()
  ENDIF()
  # NOTE: Above, we have to delete the CMakeCache.txt file only after we are
  # sure we are going to be configuring packages.  There must be a
  # CMakeCache.txt file present in the binary directory or the
  # CTEST_EMPTY_BINARY_DIRECTORY() command will *not* actually delete the
  # build directory!  Also, with updated versions of CMake (2.8.10 and above)
  # you have to delete the CMakeFiles directory in addition to the
  # CMakeCache.txt file or it will not configure correctly (due to Fortran/C
  # linkage tests for one).

  MESSAGE(
    "\n***"
    "\n*** Uploading update, notes, and the subproject dependencies XML files ..."
    "\n***\n"
    )

  IF (EXISTS ${CTEST_BINARY_DIRECTORY}/Updates.txt)
    SET(CTEST_NOTES_FILES_WO_CACHE
      "${CTEST_BINARY_DIRECTORY}/Updates.txt;${CTEST_NOTES_FILES}")
  ELSE()
    SET(CTEST_NOTES_FILES_WO_CACHE "${CTEST_NOTES_FILES}")
  ENDIF()
  PRINT_VAR(CTEST_NOTES_FILES_WO_CACHE)

  # Note: We must only do the submit after we have decided if there are any
  # packages to enable or not and otherwise exit the script!

  IF (UPDATE_FAILED)
    MESSAGE("The VC update failed so submitting update and stopping ...")
    IF (CTEST_DO_SUBMIT)
      TRIBITS_CTEST_SUBMIT( PARTS update notes )
    ENDIF()
    REPORT_QUEUED_ERRORS()
    RETURN()
  ENDIF()

  IF (CTEST_DO_SUBMIT AND EXISTS ${CDASH_SUBPROJECT_XML_FILE})
    TRIBITS_CTEST_SUBMIT( FILES ${CDASH_SUBPROJECT_XML_FILE})
    MESSAGE("\nSubmitted subproject dependencies XML file!")
  ELSE()
    MESSAGE("\nSkipping submitted subproject dependencies XML file on request!")
  ENDIF()


  MESSAGE(
    "\n***"
    "\n*** Loop through ${PROJECT_NAME} packages to configure, build, and test ..."
    "\n***")

  SET(${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE)
  SET(${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES)
  SET(${PROJECT_NAME}_FAILED_PACKAGES)
  SET(PACKAGE_IDX 0)

  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES})

    MESSAGE("")
    MESSAGE("${PACKAGE_IDX}) Processing current package ${TRIBITS_PACKAGE}:"
      " libs='${${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}}',"
      " tests='${${TRIBITS_PACKAGE}_ENABLE_TESTS}'")
    MESSAGE("")

    SET_PROPERTY(GLOBAL PROPERTY SubProject ${TRIBITS_PACKAGE})
    SET_PROPERTY(GLOBAL PROPERTY Label ${TRIBITS_PACKAGE})

    #
    # A) Configure the package and its dependent packages
    #

    MESSAGE("Configuring TRIBITS_PACKAGE='${TRIBITS_PACKAGE}'")

    # Create CONFIGURE_OPTIONS for this TRIBITS_PACKAGE
    SET( CONFIGURE_OPTIONS
	"-D${PROJECT_NAME}_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}"
      "-DCTEST_USE_LAUNCHERS:BOOL=${CTEST_USE_LAUNCHERS}"
      "-D${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON"
      "-D${PROJECT_NAME}_ENABLE_TESTS:BOOL=${${TRIBITS_PACKAGE}_ENABLE_TESTS}"
      "-D${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS:STRING=${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS}"
      "-D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON"
      "-D${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=${${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES}"
      )
    IF (NOT CTEST_GENERATE_DEPS_XML_OUTPUT_FILE)
      LIST(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE:FILEPATH=")
    ENDIF()
    IF (NOT "${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE}" STREQUAL "")
      LIST(APPEND CONFIGURE_OPTIONS
        "-D${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE:BOOL=${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE}")
    ENDIF()
    IF (NOT MPI_EXEC_MAX_NUMPROCS STREQUAL 0)
      LIST(APPEND CONFIGURE_OPTIONS
        "-DMPI_EXEC_MAX_NUMPROCS:STRING=${MPI_EXEC_MAX_NUMPROCS}")
    ENDIF()
    IF (CTEST_DO_COVERAGE_TESTING)
      LIST(APPEND CONFIGURE_OPTIONS
        "-D${PROJECT_NAME}_ENABLE_COVERAGE_TESTING:BOOL=ON")
    ENDIF()
    LIST(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_EXTRAREPOS_FILE:STRING=${${PROJECT_NAME}_EXTRAREPOS_FILE}")
    LIST(APPEND CONFIGURE_OPTIONS # See TRIBITS_SETUP_PACKAGES
      "-D${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES:BOOL=ON")
    LIST(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=${${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE}")
    IF (DEFINED ${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE)
      LIST(APPEND CONFIGURE_OPTIONS
        "-D${PROJECT_NAME}_ENABLE_${${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE}:BOOL=")
      SET(${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE)
    ENDIF()
    FOREACH(FAILED_PACKAGE ${${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES})
      LIST(APPEND CONFIGURE_OPTIONS
        "-D${PROJECT_NAME}_ENABLE_${FAILED_PACKAGE}:BOOL=OFF")
    ENDFOREACH()
    SET(CONFIGURE_OPTIONS ${CONFIGURE_OPTIONS}
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS} ${EXTRA_CONFIGURE_OPTIONS})
    LIST(APPEND CONFIGURE_OPTIONS # Package enable must be at the very end to override other stuff!
       "-D${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}:BOOL=ON" )
    MESSAGE("\nCONFIGURE_OPTIONS = '${CONFIGURE_OPTIONS}'")

    # Remember this package so we can set its enable to "" next time
    SET(${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE "${TRIBITS_PACKAGE}")

    #
    # B) Configure the package and its dependent packages
    #

    IF (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)

      CTEST_CONFIGURE(
        BUILD "${CTEST_BINARY_DIRECTORY}"
        OPTIONS "${CONFIGURE_OPTIONS}" # New option!
        RETURN_VALUE CONFIGURE_RETURN_VAL
        )

      MESSAGE("Generating the file '${CMAKE_CACHE_CLEAN_FILE}' ...")
      TRIBITS_STRIP_COMMENTS_FROM_CMAKE_CACHE_FILE(
        "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt"
        "${CMAKE_CACHE_CLEAN_FILE}"
        )

      # If the configure failed, add the package to the list
      # of failed packages
      IF (NOT "${CONFIGURE_RETURN_VAL}" EQUAL "0")
        MESSAGE("${TRIBITS_PACKAGE} FAILED to configure")
        LIST(APPEND ${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES ${TRIBITS_PACKAGE})
        LIST(APPEND ${PROJECT_NAME}_FAILED_PACKAGES ${TRIBITS_PACKAGE})
      ELSE()
        # load target properties and test keywords
        CTEST_READ_CUSTOM_FILES(BUILD "${CTEST_BINARY_DIRECTORY}")
        # Overridde from this file!
        INCLUDE("${TRIBITS_PROJECT_ROOT}/CTestConfig.cmake")
      ENDIF()

      IF (EXISTS ${CMAKE_CACHE_CLEAN_FILE})
        SET(CTEST_NOTES_FILES "${CTEST_NOTES_FILES_WO_CACHE};${CMAKE_CACHE_CLEAN_FILE}")
      ELSE()
        SET(CTEST_NOTES_FILES "${CTEST_NOTES_FILES_WO_CACHE}")
      ENDIF()
      PRINT_VAR(CTEST_NOTES_FILES)

      # Submit configure results and the notes to the dashboard
      IF (CTEST_DO_SUBMIT)
        MESSAGE("\nSubmitting configure and notes ...")
        TRIBITS_CTEST_SUBMIT( PARTS configure notes )
      ENDIF()

    ENDIF()

    #
    # C) If configure passed then try the build.  Otherwise, move on to
    # to the next package.
    #

    IF ("${CONFIGURE_RETURN_VAL}" EQUAL "0" AND
      NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING AND
      NOT CTEST_CONFIGURATION_UNIT_TESTING
      )

      # Start by trying to build just the libraries for the current package

      SET( CTEST_BUILD_TARGET ${TRIBITS_PACKAGE}_libs )
      MESSAGE("\nBuilding target: '${CTEST_BUILD_TARGET}' ...\n")
      CTEST_BUILD(
        BUILD "${CTEST_BINARY_DIRECTORY}"
        RETURN_VALUE  BUILD_LIBS_RETURN_VAL
        NUMBER_ERRORS  BUILD_LIBS_NUM_ERRORS
        APPEND
        )
      MESSAGE("Build return: RETURN_VALUE=${BUILD_LIBS_RETURN_VAL},"
        " NUMBER_ERRORS=${BUILD_LIBS_NUM_ERRORS}")

      # Determine if the build failed or not.

      SET(BUILD_LIBS_SUCCESS FALSE)
      IF ("${BUILD_LIBS_NUM_ERRORS}" EQUAL "0")
        SET(BUILD_LIBS_SUCCESS TRUE)
      ENDIF()
      # Above: Since make -i is used BUILD_LIBS_RETURN_VAL might be 0, but
      # if there are errors the build should fail, so both
      # BUILD_LIBS_RETURN_VAL and BUILD_LIBS_NUM_ERRORS should be 0 for a
      # good build and for the all target to be built.

      # Submit the library build results to the dashboard

      IF (CTEST_DO_SUBMIT)
        TRIBITS_CTEST_SUBMIT( PARTS build )
      ENDIF()

      # If the build of the libraries passed, then go on the build
      # the tests/examples and run them.

      IF (BUILD_LIBS_SUCCESS)

        SET(BUILD_OR_TEST_FAILED FALSE)

        # Build the ALL target, but append the results to the last build.xml
        SET(CTEST_BUILD_TARGET)
        MESSAGE("\nBuild ALL target for '${TRIBITS_PACKAGE}' ...\n")
        CTEST_BUILD(
          BUILD "${CTEST_BINARY_DIRECTORY}"
          RETURN_VALUE  BUILD_ALL_RETURN_VAL
          NUMBER_ERRORS  BUILD_ALL_NUM_ERRORS
          APPEND
          )
        MESSAGE("Build all: BUILD_ALL_NUM_ERRORS='${BUILD_ALL_NUM_ERRORS}',"
          "BUILD_ALL_RETURN_VAL='${BUILD_ALL_RETURN_VAL}'" )

        IF (NOT "${BUILD_LIBS_NUM_ERRORS}" EQUAL "0")
          SET(BUILD_OR_TEST_FAILED TRUE)
        ENDIF()

        # Submit the build for all target
        IF (CTEST_DO_SUBMIT)
          TRIBITS_CTEST_SUBMIT( PARTS build )
        ENDIF()

        IF (CTEST_DO_TEST)
          # Remove the LastTestsFailed log so we can detect if there are any
          # failed tests.
          SET(TEST_TMP_DIR "${CTEST_BINARY_DIRECTORY}/Testing/Temporary")
          FILE(GLOB logfiles "${TEST_TMP_DIR}/LastTestsFailed*.log")
          FOREACH(logfile ${logfiles})
            FILE(REMOVE "${logfile}")
          ENDFOREACH()
          # Run the tests that match the ${TRIBITS_PACKAGE} name
          MESSAGE("\nRunning test for package '${TRIBITS_PACKAGE}'"
            " (parallel level ${CTEST_PARALLEL_LEVEL}) ...\n")
          CTEST_TEST(
            BUILD "${CTEST_BINARY_DIRECTORY}"
            PARALLEL_LEVEL "${CTEST_PARALLEL_LEVEL}"
            INCLUDE_LABEL "^${TRIBITS_PACKAGE}$"
            #NUMBER_FAILED  TEST_NUM_FAILED
            )
          # See if a 'LastTestsFailed*.log' file exists to determine if there
          # are failed tests
          FILE(GLOB FAILED_TEST_LOG_FILE "${TEST_TMP_DIR}/LastTestsFailed*.log")
          PRINT_VAR(FAILED_TEST_LOG_FILE)
          IF (FAILED_TEST_LOG_FILE)
            SET(BUILD_OR_TEST_FAILED TRUE)
          ENDIF()
          # 2009/12/05: ToDo: We need to add an argument to CTEST_TEST(...)
          # called something like 'NUMBER_FAILED numFailedTests' to allow us
          # to detect when the tests have filed.
          #IF (TEST_NUM_FAILED GREATER 0)
          #  SET(BUILD_OR_TEST_FAILED TRUE)
          #ENDIF()
          IF (CTEST_DO_SUBMIT)
            TRIBITS_CTEST_SUBMIT( PARTS Test )
          ENDIF()
        ENDIF()

        IF (CTEST_DO_COVERAGE_TESTING)
          MESSAGE("\nRunning coverage for package '${TRIBITS_PACKAGE}' ...\n")
          CTEST_COVERAGE(
            BUILD "${CTEST_BINARY_DIRECTORY}"
            LABELS ${TRIBITS_PACKAGE} ${TRIBITS_PACKAGE}Libs ${TRIBITS_PACKAGE}Exes
            )
          IF (CTEST_DO_SUBMIT)
            TRIBITS_CTEST_SUBMIT( PARTS Coverage )
          ENDIF()
        ENDIF()

        IF (CTEST_DO_MEMORY_TESTING)
          MESSAGE("\nRunning memory testing for package '${TRIBITS_PACKAGE}' ...\n")
          PRINT_VAR(CTEST_MEMORYCHECK_COMMAND)
	    PRINT_VAR(CTEST_MEMORYCHECK_COMMAND_OPTIONS)
          PRINT_VAR(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE)
          CTEST_MEMCHECK(
            BUILD "${CTEST_BINARY_DIRECTORY}"
            PARALLEL_LEVEL "${CTEST_PARALLEL_LEVEL}"
            INCLUDE_LABEL "^${TRIBITS_PACKAGE}$")
          IF (CTEST_DO_SUBMIT)
            TRIBITS_CTEST_SUBMIT( PARTS MemCheck )
          ENDIF()
        ENDIF()

        IF (BUILD_OR_TEST_FAILED)
          LIST(APPEND ${PROJECT_NAME}_FAILED_PACKAGES ${TRIBITS_PACKAGE})
        ENDIF()

      ELSE()

        MESSAGE("FAILED library build for package '${TRIBITS_PACKAGE}'")
        LIST(APPEND ${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES ${TRIBITS_PACKAGE})
        LIST(APPEND ${PROJECT_NAME}_FAILED_PACKAGES ${TRIBITS_PACKAGE})

      ENDIF()

    ENDIF()

    IF (CTEST_DO_SUBMIT)
      MESSAGE("\nSubmit the update file that will trigger the notification email ...\n")
      TRIBITS_CTEST_SUBMIT( PARTS update )
    ENDIF()

    MATH(EXPR PACKAGE_IDX "${PACKAGE_IDX}+1")

  ENDFOREACH(TRIBITS_PACKAGE)

  IF(${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES)
    MESSAGE(
      "\nFinal set packages that failed to configure or have the libraries build:"
      " '${${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES}'")
  ENDIF()

  IF(${PROJECT_NAME}_FAILED_PACKAGES)
    MESSAGE(
      "\nFinal set packages that had any failures: '${${PROJECT_NAME}_FAILED_PACKAGES}'")
  ENDIF()

  # Write a file listing the packages that failed.  This will be read in on the next CI
  # iteration since these packages must be enabled
  FILE(WRITE "${FAILED_PACKAGES_FILE_NAME}" "${${PROJECT_NAME}_FAILED_PACKAGES}\n")

  # This is no longer necessary with CMake 2.8.1
  #MESSAGE("\nKill all hanging Zoltan processes ...")
  #EXECUTE_PROCESS(COMMAND killall -s 9 zdrive.exe)

  MESSAGE("\nDone with the incremental building and testing of ${PROJECT_NAME} packages!\n")

  REPORT_QUEUED_ERRORS()

  STRING(REPLACE "submit.php"
       "index.php" CDASH_PROJECT_LOCATION
       "${CTEST_DROP_LOCATION}")
  MULTILINE_SET( SEE_CDASH_LINK_STR
    "\nSee results for:\n"
    "  Site: ${CTEST_SITE}\n"
    "  Build Name: ${CTEST_BUILD_NAME}\n"
    "at:\n"
    "  http://${CTEST_DROP_SITE}${CDASH_PROJECT_LOCATION}&display=project\n"
    )
  # ToDo: Above: We would love to provide the buildID to point to the exact
  # URL but CTest does not currently give that to you.

  IF ("${${PROJECT_NAME}_FAILED_PACKAGES}" STREQUAL "")
    MESSAGE(
      "${SEE_CDASH_LINK_STR}\n"
      "TRIBITS_CTEST_DRIVER: OVERALL: ALL PASSSED\n")
  ELSE()
    MESSAGE(FATAL_ERROR
      "${SEE_CDASH_LINK_STR}\n"
      "TRIBITS_CTEST_DRIVER: OVERALL: ALL FAILED\n")
    # NOTE: FATAL_ERROR is needed so that the ctest -S script returns != 0
    # Also, it is critical to dislplay the "See results" in this
    # MESSAGE(FATAL_ERROR ...) command in order for it to be printed last.
    # Otherwise, if you run with ctest -V -S, then the ouptut from
    # CTEST_TEST() will be printed last.
  ENDIF()

ENDFUNCTION()
