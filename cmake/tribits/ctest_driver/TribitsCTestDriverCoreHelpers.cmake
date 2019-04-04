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
# Update an existing git repo
#
FUNCTION(TRIBITS_UPDATE_GIT_EXTRAREPO  GIT_EXE  EXTRAREPO_SRC_DIR)

  SET(EXTRAREPO_FETCH_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.fetch.out")
  SET(EXTRAREPO_CLEAN_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.clean.out")
  SET(EXTRAREPO_RESET_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.reset.out")
  SET(EXTRAREPO_SET_BRANCH_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.set_branch.out")

  SET(FETCH_CMND_ARGS
    COMMAND "${GIT_EXE}" fetch ${${PROJECT_NAME}_GIT_REPOSITORY_REMOTE}
    TIMEOUT 600 # seconds
    WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
    OUTPUT_FILE "${EXTRAREPO_FETCH_OUT_FILE}" )
  SET(CLEAN_CMND_ARGS
    COMMAND "${GIT_EXE}" clean -fdx
    WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
    OUTPUT_FILE "${EXTRAREPO_CLEAN_OUT_FILE}" )
  SET(RESET_CMND_ARGS
    COMMAND "${GIT_EXE}" reset --hard HEAD
    WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
    OUTPUT_FILE "${EXTRAREPO_RESET_OUT_FILE}" )
  IF (${PROJECT_NAME}_EXTRAREPOS_BRANCH)
    SET(SET_BRANCH_CMND_ARGS
      COMMAND "${GIT_EXE}" checkout -B ${${PROJECT_NAME}_EXTRAREPOS_BRANCH}
        --track ${${PROJECT_NAME}_GIT_REPOSITORY_REMOTE}/${${PROJECT_NAME}_EXTRAREPOS_BRANCH}
      WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
      OUTPUT_FILE "${EXTRAREPO_SET_BRANCH_OUT_FILE}" )
  ELSE ()
    SET(SET_BRANCH_CMND_ARGS
      COMMAND "${GIT_EXE}" reset --hard "@{u}"
      WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
      OUTPUT_FILE "${EXTRAREPO_SET_BRANCH_OUT_FILE}" )
  ENDIF()

  EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${FETCH_CMND_ARGS})
  EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${CLEAN_CMND_ARGS})
  EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${RESET_CMND_ARGS})
  EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${SET_BRANCH_CMND_ARGS})

ENDFUNCTION()


#
# Update or clone a single extra repo
#
FUNCTION(TRIBITS_CLONE_OR_UPDATE_EXTRAREPO  EXTRAREPO_NAME_IN  EXTRAREPO_DIR_IN
  EXTRAREPO_REPOTYPE_IN  EXTRAREPO_REPOURL_IN
  )

  #MESSAGE("TRIBITS_CLONE_OR_UPDATE_EXTRAREPO: ${EXTRAREPO_NAME_IN} ${EXTRAREPO_REPOURL_IN}")

  SET(EXTRAREPO_SRC_DIR "${${PROJECT_NAME}_SOURCE_DIRECTORY}/${EXTRAREPO_DIR_IN}")
  #PRINT_VAR(EXTRAREPO_SRC_DIR)

  SET(EXTRAREPO_CLONE_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.clone.out")
  SET(EXTRAREPO_CHECKOUT_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.checkout.out")

  IF (NOT EXISTS "${EXTRAREPO_SRC_DIR}")

    MESSAGE("\n${EXTRAREPO_NAME_IN}: Doing initial ${EXTRAREPO_REPOTYPE_IN}"
      " clone/checkout from URL '${EXTRAREPO_REPOURL_IN}' to dir '${EXTRAREPO_DIR_IN}' ...")

    # Set the command to clone
    IF (${EXTRAREPO_REPOTYPE_IN} STREQUAL GIT)
      IF (${PROJECT_NAME}_EXTRAREPOS_BRANCH) 
        SET(CHECKOUT_BRANCH_ARG -b ${${PROJECT_NAME}_EXTRAREPOS_BRANCH})
      ELSE()
        SET(CHECKOUT_BRANCH_ARG)
      ENDIF()
      SET(CLONE_CMND_ARGS
        COMMAND "${GIT_EXECUTABLE}" clone
        ${CHECKOUT_BRANCH_ARG} -o ${${PROJECT_NAME}_GIT_REPOSITORY_REMOTE}
	"${EXTRAREPO_REPOURL}" ${EXTRAREPO_DIR_IN}
        WORKING_DIRECTORY "${${PROJECT_NAME}_SOURCE_DIRECTORY}"
        OUTPUT_FILE "${EXTRAREPO_CLONE_OUT_FILE}" )
    ELSE()
      MESSAGE(SEND_ERROR
	"Error, Invalid EXTRAREPO_REPOTYPE_IN='${EXTRAREPO_REPOTYPE_IN}'!")
    ENDIF()

    # Do the clone
    EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${CLONE_CMND_ARGS})

  ELSE()

    MESSAGE("\n${EXTRAREPO_NAME_IN}: Doing ${EXTRAREPO_REPOTYPE_IN} update"
      " from URL '${EXTRAREPO_REPOURL_IN}' to dir '${EXTRAREPO_SRC_DIR}' ...")

  ENDIF()

  IF (${EXTRAREPO_REPOTYPE_IN} STREQUAL GIT)
    # Always update the git repo, even after a clone.  See
    # TRIBITS_CTEST_DRIVER() documentation.
    TRIBITS_UPDATE_GIT_EXTRAREPO("${GIT_EXECUTABLE}" "${EXTRAREPO_SRC_DIR}")
  ELSE()
    MESSAGE(SEND_ERROR
      "Error, Invalid EXTRAREPO_REPOTYPE_IN='${EXTRAREPO_REPOTYPE_IN}'!")
  ENDIF()

ENDFUNCTION()


#
# Clone or update all of the extra repos and put them on the right branch.
#
# NOTE: The base repo is cloned by CTEST_START() and updated by CTEST_UPDATE()
# before calling this function.  This function only operates on the extra
# repos.
#
FUNCTION(TRIBITS_CLONE_OR_UPDATE_EXTRA_REPOS  CTEST_UPDATE_RETURN_VAL
  UPDATE_FAILED_VAR_OUT
  )

  SET(UPDATE_FAILED FALSE)

  IF (${PROJECT_NAME}_EXTRAREPOS_BRANCH)
    MESSAGE("For extra repos, doing switch to branch ${${PROJECT_NAME}_EXTRAREPOS_BRANCH}")
  ENDIF()

  SET(EXTRAREPO_IDX 0)
  FOREACH(EXTRAREPO_NAME ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES})
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS ${EXTRAREPO_IDX}
      EXTRAREPO_DIR )
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES ${EXTRAREPO_IDX}
      EXTRAREPO_REPOTYPE )
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS ${EXTRAREPO_IDX}
      EXTRAREPO_REPOURL )
    TRIBITS_CLONE_OR_UPDATE_EXTRAREPO( ${EXTRAREPO_NAME} ${EXTRAREPO_DIR}
      ${EXTRAREPO_REPOTYPE} ${EXTRAREPO_REPOURL} )
    # ToDo: Detect and return failure in cloning or updating extra repos!
    MATH(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")
  ENDFOREACH()

  SET(${UPDATE_FAILED_VAR_OUT} ${UPDATE_FAILED} PARENT_SCOPE)

ENDFUNCTION()


#
# Create the Updates.txt file
#
FUNCTION(TRIBITS_CREATE_REPO_UPDATES_FILE)
  EXTRAREPO_EXECUTE_PROCESS_WRAPPER(
    COMMAND ${PYTHON_EXECUTABLE}
      ${GITDIST_EXE} --dist-no-color
      log "--pretty=format:%h:  %s%nAuthor: %an <%ae>%nDate:   %ad%n"
      --name-status -C ORIG_HEAD..HEAD
    WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
    OUTPUT_FILE "${CTEST_BINARY_DIRECTORY}/Updates.txt"
    )
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
    MESSAGE("Setting ${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON since"
      " ${PROJECT_NAME}_PACKAGES_USER_SELECTED='${${PROJECT_NAME}_PACKAGES_USER_SELECTED}'")
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
    COMMAND "${GIT_EXECUTABLE}" diff --name-only ORIG_HEAD..HEAD
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

    # For now, only look for changes if it has packages.  Later, we may need
    # to generalize this for the general extra repo case with deeper directory
    # and other VC systems than GIT.
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
        --project-dir=${TRIBITS_PROJECT_ROOT}
        --deps-xml-file=${CTEST_BINARY_DIRECTORY}/${${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME}
      OUTPUT_VARIABLE MODIFIED_PACKAGES_LIST
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  ELSE()
    SET(MODIFIED_PACKAGES_LIST)
  ENDIF()

  SPLIT("${MODIFIED_PACKAGES_LIST}" "," MODIFIED_PACKAGES_LIST)
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

  IF (FAILING_PACKAGES_LIST STREQUAL "ALL_PACKAGES")
    MESSAGE("Enabling previously failing ALL_PACKAGES")
    SET(${PROJECT_NAME}_ENABLE_ALL_PACKAGES ON)
  ELSE()
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
  ENDIF()

  #
  # D) Print the final status
  #

  IF (${PROJECT_NAME}_ENABLE_ALL_PACKAGES)
    IF (NOT ${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE)
      MESSAGE(FATAL_ERROR
	"Error, failing 'ALL_PACKAGES' only allowed with all-at-once mode!")
    ENDIF()
    MESSAGE("\nDirectly modified or failing non-disabled packages that need"
      " to be tested:  ALL_PACKAGES")
  ELSE()
    TRIBITS_PRINT_ENABLED_SE_PACKAGE_LIST(
      "\nDirectly modified or failing non-disabled packages that need to be tested"
      ON FALSE )
  ENDIF()

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
MACRO(SELECT_FINAL_SET_OF_PACKAGES_TO_DIRECTLY_TEST)

  SET(${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST)

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
      APPEND_SET(${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST  ${TRIBITS_PACKAGE})
    ENDIF()

  ENDFOREACH()

ENDMACRO()


#
# Set mapping of labels to subprojects (i.e. TriBITS packages) for CDash.
#
# NOTE: Unlike for the inner CMake configure, only subprojects that are
# explicitly tested will be marked as a CDash subproject.  This limits the
# rows in CDash.  This does not seem to be a problem for when running ctest
# locally.  When run locally, ctest will just report aggregated times for
# subprojects that have 1 or more tests.  Not true for CDash.
#
MACRO(TRIBITS_CTEST_DRIVER_SET_LABELS_TO_SUBPROJECTS_MAPPING)
  SET(CTEST_LABELS_FOR_SUBPROJECTS)
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST})
    LIST(APPEND CTEST_LABELS_FOR_SUBPROJECTS ${TRIBITS_PACKAGE})
  ENDFOREACH()
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
# Call REPORT_QUEUED_ERRORS() once at the bottom of TRIBITS_CTEST_DRIVER()
#
MACRO(REPORT_QUEUED_ERRORS)
  IF ("${TRIBITS_CTEST_DRIVER_ERROR_QUEUE}" STREQUAL "")
    MESSAGE("TRIBITS_CTEST_DRIVER_ERROR_QUEUE is empty. All is well.")
  ELSE()
    MESSAGE("ERROR: TRIBITS_CTEST_DRIVER_ERROR_QUEUE reports the following error message queue:")
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
  SET(retry_args
    RETRY_COUNT ${CTEST_SUBMIT_RETRY_COUNT}
    RETRY_DELAY ${CTEST_SUBMIT_RETRY_DELAY})
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
# Helper macros to pass through common CMake configure arguments used by both
# package-by-pacakge approach and all-at-once approach
#
MACRO(TRIBITS_FWD_CMAKE_CONFIG_ARGS_0)
  SET( CONFIGURE_OPTIONS
    "-D${PROJECT_NAME}_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}"
    "-DCTEST_USE_LAUNCHERS:BOOL=${CTEST_USE_LAUNCHERS}"
    "-D${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON"
"-D${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS:STRING=${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS}"
    "-D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON"
    "-D${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=${${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES}"
    )
  IF (NOT CTEST_GENERATE_DEPS_XML_OUTPUT_FILE)
    LIST(APPEND CONFIGURE_OPTIONS
    "-D${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE:FILEPATH=")
  ENDIF()
  IF (NOT "${${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES}" STREQUAL "")
    LIST(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES:BOOL=${${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES}")
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
  IF (${PROJECT_NAME}_EXTRAREPOS_FILE STREQUAL "NONE")
    SET(EXTRAREOS_FILE_PASSED "")
  ELSE()
    SET(EXTRAREOS_FILE_PASSED "${${PROJECT_NAME}_EXTRAREPOS_FILE}")
  ENDIF()
  LIST(APPEND CONFIGURE_OPTIONS
    "-D${PROJECT_NAME}_EXTRAREPOS_FILE:STRING=${EXTRAREOS_FILE_PASSED}")
  LIST(APPEND CONFIGURE_OPTIONS # See TRIBITS_SETUP_PACKAGES
    "-D${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES:BOOL=ON")
  LIST(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=${${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE}")
ENDMACRO()


MACRO(TRIBITS_FWD_CMAKE_CONFIG_ARGS_1)
  SET(CONFIGURE_OPTIONS ${CONFIGURE_OPTIONS}
    ${EXTRA_SYSTEM_CONFIGURE_OPTIONS} ${EXTRA_CONFIGURE_OPTIONS})
ENDMACRO()


# Remove the all of the LastTestsFailed*.log files so we can determine if any
# tests have failed.
#
MACRO(TRIBITS_REMOVE_LAST_TEST_FAILED_LOG_FILE)
  # Remove the LastTestsFailed log so we can detect if there are any
  # failed tests.
  SET(TEST_TMP_DIR "${CTEST_BINARY_DIRECTORY}/Testing/Temporary")
  SET(LAST_TESTS_FILED_LOG_FILE_GLOB "${TEST_TMP_DIR}/LastTestsFailed*.log")
  FILE(GLOB logfiles "${LAST_TESTS_FILED_LOG_FILE_GLOB}")
  FOREACH(logfile ${logfiles})
    FILE(REMOVE "${logfile}")
  ENDFOREACH()
ENDMACRO()


# Sets the var FAILED_TEST_LOG_FILE if the file is found
MACRO(TRIBITS_FIND_LAST_TEST_FAILED_LOG_FILE)
  FILE(GLOB FAILED_TEST_LOG_FILE "${LAST_TESTS_FILED_LOG_FILE_GLOB}")
ENDMACRO()


# Get names of failed packages from failed tests
FUNCTION(TRIBITS_GET_FAILED_PACKAGES_FROM_FAILED_TESTS
   LAST_TESTS_FAILED_FILE  FAILED_PACKAGES_OUT
   )
  EXECUTE_PROCESS(
    COMMAND ${PYTHON_EXECUTABLE}
      "${${PROJECT_NAME}_TRIBITS_DIR}/ci_support/get-tribits-packages-from-last-tests-failed.py"
      "--deps-xml-file=${CTEST_BINARY_DIRECTORY}/${${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME}"
      "--last-tests-failed-file=${LAST_TESTS_FAILED_FILE}"
          OUTPUT_VARIABLE  FAILED_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  SPLIT("${FAILED_PACKAGES}" "," FAILED_PACKAGES)
  SET(${FAILED_PACKAGES_OUT} "${FAILED_PACKAGES}" PARENT_SCOPE)
ENDFUNCTION()


#
# Drive the configure, build, test, and submit package-by-package
#
# Sets ${PROJECT_NAME}_FAILED_PACKAGES as an indication if there are any
# failures.
#
MACRO(TRIBITS_CTEST_PACKAGE_BY_PACKAGE)

  MESSAGE(
    "\n***"
    "\n*** Loop through ${PROJECT_NAME} packages to configure, build, and test ..."
    "\n***")

  SET(${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE)
  SET(${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES)
  SET(PACKAGE_IDX 0)

  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST})

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
    TRIBITS_FWD_CMAKE_CONFIG_ARGS_0()
    LIST(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_ENABLE_TESTS:BOOL=${${TRIBITS_PACKAGE}_ENABLE_TESTS}")
    IF (DEFINED ${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE)
      LIST(APPEND CONFIGURE_OPTIONS
        "-D${PROJECT_NAME}_ENABLE_${${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE}:BOOL=")
      SET(${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE)
    ENDIF()
    FOREACH(FAILED_PACKAGE ${${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES})
      LIST(APPEND CONFIGURE_OPTIONS
        "-D${PROJECT_NAME}_ENABLE_${FAILED_PACKAGE}:BOOL=OFF")
    ENDFOREACH()
    TRIBITS_FWD_CMAKE_CONFIG_ARGS_1()
    LIST(APPEND CONFIGURE_OPTIONS # Package enable must be at the very end to override other stuff!
       "-D${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}:BOOL=ON" )
    MESSAGE("\nCONFIGURE_OPTIONS = '${CONFIGURE_OPTIONS}'")

    # Remember this package so we can set its enable to "" next time
    SET(${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE "${TRIBITS_PACKAGE}")

    #
    # B) Configure the package and its dependent packages
    #

    SET(PBP_CONFIGURE_PASSED TRUE)

    IF (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)

      MESSAGE("${TRIBITS_PACKAGE}: Skipping configure due"
	" to running in unit testing mode!")

    ELSE()

      #
      # We always have to configure if we are going to do anything for the
      # package.  We just want submit configure results to CDash if we are not
      # asked to configure!
      #

      SET(PBP_CONFIGURE_PASSED FALSE)

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
      IF ("${CONFIGURE_RETURN_VAL}" EQUAL "0")
        MESSAGE("\n${TRIBITS_PACKAGE}: Configure passed!\n")
        SET(PBP_CONFIGURE_PASSED TRUE)
        # load target properties and test keywords
        CTEST_READ_CUSTOM_FILES(BUILD "${CTEST_BINARY_DIRECTORY}")
        # Overridde from this file!
        INCLUDE("${TRIBITS_PROJECT_ROOT}/CTestConfig.cmake")
      ELSE()
        MESSAGE("\n${TRIBITS_PACKAGE} FAILED to configure!\n")
      ENDIF()

      IF (EXISTS ${CMAKE_CACHE_CLEAN_FILE})
        SET(CTEST_NOTES_FILES "${CTEST_NOTES_FILES_WO_CACHE};${CMAKE_CACHE_CLEAN_FILE}")
      ELSE()
        SET(CTEST_NOTES_FILES "${CTEST_NOTES_FILES_WO_CACHE}")
      ENDIF()

      PRINT_VAR(CTEST_NOTES_FILES)

      IF (NOT CTEST_DO_CONFIGURE AND CTEST_DO_SUBMIT)
        MESSAGE("${TRIBITS_PACKAGE}: Skipping submitting configure"
	  " and notes due to CTEST_DO_CONFIGURE='${CTEST_DO_CONFIGURE}'!")
      ELSEIF (CTEST_DO_SUBMIT)
        MESSAGE("\nSubmitting configure and notes ...")
        TRIBITS_CTEST_SUBMIT( PARTS configure notes )
      ENDIF()

    ENDIF()

    #
    # C) Build the library and then ALL
    #

    SET(PBP_BUILD_PASSED TRUE)
    SET(PBP_BUILD_LIBS_PASSED TRUE)

    PRINT_VAR(PBP_CONFIGURE_PASSED)

    IF ( NOT PBP_CONFIGURE_PASSED AND CTEST_DO_BUILD )

      MESSAGE("\n${TRIBITS_PACKAGE}: Skipping build due"
	" to configure failing!")

      SET(PBP_BUILD_PASSED FALSE)
      SET(PBP_BUILD_LIBS_PASSED FALSE)

    ELSEIF (NOT CTEST_DO_BUILD)

      MESSAGE("\n${TRIBITS_PACKAGE}: Skipping build due"
	" to CTEST_DO_BUILD='${CTEST_DO_BUILD}'!")

    ELSEIF (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING OR
      CTEST_CONFIGURATION_UNIT_TESTING
      )

      MESSAGE("\n${TRIBITS_PACKAGE}: Skipping build due"
	" to running in unit testing mode!")

    ELSE()

      # Start by trying to build just the libraries for the current package

      SET( CTEST_BUILD_TARGET ${TRIBITS_PACKAGE}_libs )
      MESSAGE("\nBuilding target: '${CTEST_BUILD_TARGET}' ...\n")
      SET(PBP_BUILD_LIBS_PASSED FALSE)
      CTEST_BUILD(
        BUILD "${CTEST_BINARY_DIRECTORY}"
        RETURN_VALUE  BUILD_LIBS_RETURN_VAL
        NUMBER_ERRORS  BUILD_LIBS_NUM_ERRORS
        APPEND
        )
      MESSAGE("Build return: RETURN_VALUE=${BUILD_LIBS_RETURN_VAL},"
        " NUMBER_ERRORS=${BUILD_LIBS_NUM_ERRORS}")

      # Determine if the build failed or not.

      IF ("${BUILD_LIBS_NUM_ERRORS}" EQUAL "0")
        MESSAGE("\n${TRIBITS_PACKAGE}: Libs build passed!")
        SET(PBP_BUILD_LIBS_PASSED TRUE)
      ELSE()
        MESSAGE("\nFAILED library build for package '${TRIBITS_PACKAGE}'!")
        SET(PBP_BUILD_PASSED FALSE)
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

      IF (PBP_BUILD_LIBS_PASSED)

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

        IF (NOT "${BUILD_ALL_NUM_ERRORS}" EQUAL "0")
          MESSAGE("${TRIBITS_PACKAGE}: All build FAILED!")
          SET(PBP_BUILD_PASSED FALSE)
        ELSE()
          MESSAGE("${TRIBITS_PACKAGE}: All build passed!")
        ENDIF()

        # Submit the build for all target
        IF (CTEST_DO_SUBMIT)
          TRIBITS_CTEST_SUBMIT( PARTS build )
        ENDIF()

      ENDIF()

    ENDIF()

    #
    # D) Run the tests
    #

    SET(PBP_TESTS_PASSED TRUE)

    IF (NOT PBP_BUILD_LIBS_PASSED AND CTEST_DO_TEST)

      MESSAGE("\n${TRIBITS_PACKAGE}: Skipping tests since libray build failed!\n")

      SET(PBP_TESTS_PASSED FALSE)

    ELSEIF (NOT CTEST_DO_TEST)

      MESSAGE("\n${TRIBITS_PACKAGE}: Skipping running tests due"
        " to CTEST_DO_TEST='${CTEST_DO_TEST}'!")

    ELSE()
      
      #
      # D.1) Run the regular tests
      #

      SET(PBP_TESTS_PASSED FALSE)

      # Run the tests that match the ${TRIBITS_PACKAGE} name
      MESSAGE("\nRunning test for package '${TRIBITS_PACKAGE}'"
        " (parallel level ${CTEST_PARALLEL_LEVEL}) ...\n")
      TRIBITS_REMOVE_LAST_TEST_FAILED_LOG_FILE()
      CTEST_TEST(
        BUILD "${CTEST_BINARY_DIRECTORY}"
        PARALLEL_LEVEL "${CTEST_PARALLEL_LEVEL}"
        INCLUDE_LABEL "^${TRIBITS_PACKAGE}$"
          )
      # See if a 'LastTestsFailed*.log' file exists to determine if there are
      # failed tests
      TRIBITS_FIND_LAST_TEST_FAILED_LOG_FILE()
      IF (FAILED_TEST_LOG_FILE)
        MESSAGE("\n${TRIBITS_PACKAGE}: File '${FAILED_TEST_LOG_FILE}'"
          " exists so there were failed tests!")
      ELSE()
        MESSAGE("\n${TRIBITS_PACKAGE}: File '${FAILED_TEST_LOG_FILE}'"
	  " does NOT exist so all tests passed!")
        SET(PBP_TESTS_PASSED TRUE)
      ENDIF()
      # 2009/12/05: ToDo: We need to add an argument to CTEST_TEST(...)
      # called something like 'NUMBER_FAILED numFailedTests' to allow us to
      # detect when the tests have filed.
      IF (CTEST_DO_SUBMIT)
        TRIBITS_CTEST_SUBMIT( PARTS Test )
      ENDIF()

      #
      # D.2) Collect coverage results
      #

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

    ENDIF()

    #
    # E) Run memory testing
    #

    IF (NOT PBP_BUILD_LIBS_PASSED AND CTEST_DO_MEMORY_TESTING)

      MESSAGE("\n${TRIBITS_PACKAGE}: Skipping running memory checking"
	 "tests since libray build failed!\n")

    ELSEIF (NOT CTEST_DO_MEMORY_TESTING)

      MESSAGE("\n${TRIBITS_PACKAGE}: Skipping running memory checking tests due"
        " to CTEST_DO_MEMORY_TESTING='${CTEST_DO_MEMORY_TESTING}'!")

    ELSE()

      MESSAGE("\nRunning memory testing for package '${TRIBITS_PACKAGE}' ...\n")

      PRINT_VAR(CTEST_MEMORYCHECK_COMMAND)
      PRINT_VAR(CTEST_MEMORYCHECK_COMMAND_OPTIONS)
      PRINT_VAR(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE)

      CTEST_MEMCHECK(
        BUILD "${CTEST_BINARY_DIRECTORY}"
        PARALLEL_LEVEL "${CTEST_PARALLEL_LEVEL}"
        INCLUDE_LABEL "^${TRIBITS_PACKAGE}$"
	)
      # ToDo: Determine if memory testing passed or not and affect overall
      # pass/fail!

      IF (CTEST_DO_SUBMIT)
        TRIBITS_CTEST_SUBMIT( PARTS MemCheck )
      ENDIF()

    ENDIF()

    #
    # F) Record if this package failed the build or any tests
    #

    IF (NOT PBP_CONFIGURE_PASSED OR NOT PBP_BUILD_LIBS_PASSED)
      LIST(APPEND ${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES ${TRIBITS_PACKAGE})
    ENDIF()

    IF (NOT PBP_BUILD_PASSED OR NOT PBP_TESTS_PASSED)
      LIST(APPEND ${PROJECT_NAME}_FAILED_PACKAGES ${TRIBITS_PACKAGE})
    ENDIF()

    #
    # G) Do submit of update
    #

    IF (CTEST_DO_SUBMIT)
      MESSAGE("\nSubmit the update file that will trigger the notification email ...\n")
      TRIBITS_CTEST_SUBMIT( PARTS update )
    ENDIF()

    MATH(EXPR PACKAGE_IDX "${PACKAGE_IDX}+1")

  ENDFOREACH(TRIBITS_PACKAGE)

  IF (${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES)
    MESSAGE(
      "\nFinal set packages that failed to configure or have the libraries build:"
      " '${${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES}'")
  ENDIF()

  MESSAGE("\nDone with the incremental building and testing of"
    " ${PROJECT_NAME} packages!\n")

ENDMACRO()


#
# Drive the configure, build, test, and submit all at once for all of the
# enabled packages.
#
# Sets ${PROJECT_NAME}_FAILED_PACKAGES as an indication if there are any
# failures.
#
MACRO(TRIBITS_CTEST_ALL_AT_ONCE)

  MESSAGE(
    "\n***"
    "\n*** Configure, build, test and submit results all-at-once for all enabled packages ..."
    "\n***")

  SET(AAO_CONFIGURE_FAILED FALSE)
  SET(AAO_BUILD_FAILED FALSE)
  SET(AAO_INSTALL_FAILED FALSE)

  #
  # A) Define mapping from labels to subprojects and gather configure arguments
  #

  TRIBITS_CTEST_DRIVER_SET_LABELS_TO_SUBPROJECTS_MAPPING()
  PRINT_VAR(CTEST_LABELS_FOR_SUBPROJECTS)

  MESSAGE("")
  MESSAGE("Configuring ...")
  MESSAGE("")

  # Create CONFIGURE_OPTIONS
  TRIBITS_FWD_CMAKE_CONFIG_ARGS_0()
  IF (${PROJECT_NAME}_ENABLE_ALL_PACKAGES)
    LIST(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON" )
    FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_EXCLUDE_PACKAGES})
      LIST(APPEND CONFIGURE_OPTIONS
        "-D${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}=OFF" )
    ENDFOREACH()
    # NOTE: Above we have to explicitly set disables for the excluded pacakges
    # since we are pssing in ${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON.  This is
    # effectively the "black-listing" approach.
  ELSE()
    FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST})
      LIST(APPEND CONFIGURE_OPTIONS
         "-D${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}=ON" )
    ENDFOREACH()
    # NOTE: Above we don't have to consider the packages excluded in
    # ${PROJECT_NAME}_EXCLUDE_PACKAGES because they are not enabled ad this
    # point and therefore no in ${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST.
    # This is effectively the "white-listing" approach.
  ENDIF()
  LIST(APPEND CONFIGURE_OPTIONS
    "-D${PROJECT_NAME}_ENABLE_TESTS:BOOL=ON")
  TRIBITS_FWD_CMAKE_CONFIG_ARGS_1()
  MESSAGE("\nCONFIGURE_OPTIONS = '${CONFIGURE_OPTIONS}'")

  #
  # B) Configure the package and its dependent packages
  #

  IF (NOT CTEST_DO_CONFIGURE)

    MESSAGE("\nSkipping configure due to CTEST_DO_CONFIGURE='${CTEST_DO_CONFIGURE}'!\n")
    SET(AAO_CONFIGURE_PASSED TRUE)
    # Just assume configure passeed for the purpose of running the build.

  ELSEIF (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)

    MESSAGE("Skipping actual ctest_configure() because"
      " CTEST_DEPENDENCY_HANDLING_UNIT_TESTING='${CTEST_DEPENDENCY_HANDLING_UNIT_TESTING}'!"
      )
    SET(AAO_CONFIGURE_PASSED TRUE)

  ELSE()

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
    
    IF (NOT "${CONFIGURE_RETURN_VAL}" EQUAL "0")
      MESSAGE("Configure FAILED!")
      SET(AAO_CONFIGURE_PASSED FALSE)
      SET(AAO_CONFIGURE_FAILED TRUE)
    ELSE()
      MESSAGE("Configure PASSED!")
      SET(AAO_CONFIGURE_PASSED TRUE)
      # Load target properties and test keywords
      CTEST_READ_CUSTOM_FILES(BUILD "${CTEST_BINARY_DIRECTORY}")
      # Overridde from this file!
      INCLUDE("${TRIBITS_PROJECT_ROOT}/CTestConfig.cmake")
    ENDIF()
  
    SET(CTEST_NOTES_FILES "${CTEST_NOTES_FILES_WO_CACHE}")
  
    IF (EXISTS ${CMAKE_CACHE_CLEAN_FILE})
      LIST(APPEND CTEST_NOTES_FILES "${CMAKE_CACHE_CLEAN_FILE}")
    ENDIF()
  
    IF (EXISTS "${REPO_VERSION_FILE}")
      SET(CTEST_NOTES_FILES "${REPO_VERSION_FILE};${CTEST_NOTES_FILES}")
    ENDIF()
  
    PRINT_VAR(CTEST_NOTES_FILES)
  
    # Submit configure results and the notes to the dashboard
    IF (CTEST_DO_SUBMIT)
      MESSAGE("\nSubmitting update, configure and notes ...")
      TRIBITS_CTEST_SUBMIT( PARTS update configure notes )
    ENDIF()

  ENDIF()

  #
  # C) Do the build
  #

  IF (NOT CTEST_DO_BUILD)

    MESSAGE("\nSkipping build due to CTEST_DO_BUILD='${CTEST_DO_BUILD}'!\n")

  ELSEIF (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING AND AAO_CONFIGURE_PASSED)

    MESSAGE("Skipping build because"
      " CTEST_DEPENDENCY_HANDLING_UNIT_TESTING='${CTEST_DEPENDENCY_HANDLING_UNIT_TESTING}'!"
      )

  ELSEIF (AAO_CONFIGURE_PASSED)
  
    MESSAGE("")
    MESSAGE("Building all targets ...")
    MESSAGE("")
  
    CTEST_BUILD(
      BUILD "${CTEST_BINARY_DIRECTORY}"
      RETURN_VALUE  BUILD_ALL_RETURN_VAL
      NUMBER_ERRORS  BUILD_ALL_NUM_ERRORS
      )
    MESSAGE("Build output: BUILD_ALL_NUM_ERRORS='${BUILD_ALL_NUM_ERRORS}',"
      "BUILD_ALL_RETURN_VAL='${BUILD_ALL_RETURN_VAL}'" )
  
    IF (NOT "${BUILD_ALL_NUM_ERRORS}" EQUAL "0")
      MESSAGE("Build FAILED!")
      SET(AAO_BUILD_FAILED TRUE)
    ELSE()
      MESSAGE("Build PASSED!")
    ENDIF()
  
    # Submit the build for all target
    IF (CTEST_DO_SUBMIT)
      TRIBITS_CTEST_SUBMIT( PARTS build )
    ENDIF()

    IF (CTEST_DO_INSTALL)

      MESSAGE("")
      MESSAGE("Installing (i.e. building target 'install') ...")
      MESSAGE("")

      CTEST_BUILD(
        BUILD "${CTEST_BINARY_DIRECTORY}"
        TARGET install
        RETURN_VALUE  BUILD_INSTALL_RETURN_VAL
        NUMBER_ERRORS  BUILD_INSTALL_NUM_ERRORS
        )
      MESSAGE("Build install output:"
        " BUILD_INSTALL_NUM_ERRORS='${BUILD_INSTALL_NUM_ERRORS}',"
        "BUILD_INSTALL_RETURN_VAL='${BUILD_INSTALL_RETURN_VAL}'" )

      IF (NOT "${BUILD_INSTALL_NUM_ERRORS}" EQUAL "0")
        MESSAGE("Install FAILED!")
        SET(AAO_INSTALL_FAILED TRUE)
      ELSE()
        MESSAGE("Install PASSED!")
      ENDIF()

      # Submit the build for all target
      IF (CTEST_DO_SUBMIT)
        TRIBITS_CTEST_SUBMIT( PARTS build )
      ENDIF()

    ENDIF()

  ELSE()
  
    MESSAGE("")
    MESSAGE("Skipping build because configure failed!")
    MESSAGE("")
  
  ENDIF()

  #
  # D) Run tests
  #

  IF (NOT CTEST_DO_TEST)
  
    MESSAGE("")
    MESSAGE("Skipping tests because CTEST_DO_TEST='${CTEST_DO_TEST}'!")
    MESSAGE("")

  ELSEIF (NOT AAO_CONFIGURE_PASSED)
  
    MESSAGE("")
    MESSAGE("Skipping tests because configure failed!")
    MESSAGE("")

  ELSEIF (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING AND AAO_CONFIGURE_PASSED)

    MESSAGE("Skipping testing because"
      " CTEST_DEPENDENCY_HANDLING_UNIT_TESTING='${CTEST_DEPENDENCY_HANDLING_UNIT_TESTING}'!"
      )

  ELSE()

    # NOTE: We always run the tests if the configure passed no matter if there
    # are build failures because the only way that we can detect what packages
    # have build failures is to see what packages have test failures.

    TRIBITS_REMOVE_LAST_TEST_FAILED_LOG_FILE()

    # Run the tests
    MESSAGE("")
    MESSAGE("\nRunning tests (parallel level ${CTEST_PARALLEL_LEVEL}) ...\n")
    MESSAGE("")

    CTEST_TEST(
      BUILD "${CTEST_BINARY_DIRECTORY}"
      PARALLEL_LEVEL "${CTEST_PARALLEL_LEVEL}"
      )

    # See if a 'LastTestsFailed*.log' file exists to determine if there are
    # failed tests.
    TRIBITS_FIND_LAST_TEST_FAILED_LOG_FILE()
    IF (FAILED_TEST_LOG_FILE)
      MESSAGE("File '${FAILED_TEST_LOG_FILE}' exists so there were non-passing tests!")
    ELSE()
      MESSAGE("File '${FAILED_TEST_LOG_FILE}' does NOT exist so all tests passed!")
    ENDIF()

    IF (CTEST_DO_SUBMIT)
      TRIBITS_CTEST_SUBMIT( PARTS Test )
    ENDIF()

  ENDIF()

  #
  # E) Gather coverage results
  #

  IF (NOT CTEST_DO_COVERAGE_TESTING)
  
    MESSAGE("")
    MESSAGE("Skipping converage tests because CTEST_DO_COVERAGE_TESTING='${CTEST_DO_COVERAGE_TESTING}'!")
    MESSAGE("")

  ELSEIF (NOT AAO_CONFIGURE_PASSED)
  
    MESSAGE("")
    MESSAGE("Skipping coverage tests because configure failed!")
    MESSAGE("")

  ELSEIF (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING AND AAO_CONFIGURE_PASSED)

    MESSAGE("Skipping coverage testing because"
      " CTEST_DEPENDENCY_HANDLING_UNIT_TESTING='${CTEST_DEPENDENCY_HANDLING_UNIT_TESTING}'!"
      )

  ELSE()
    
    # NOTE: We always gather the coverage results if the configure passed
    # independent if there was any build or test failures.  The coverage stats
    # may not be very valid if there are build or test failures but there is
    # no harm and showing the coverage based on tests that actually run (even
    # if they fail).

    MESSAGE("\nGathering coverage results ...\n")
    CTEST_COVERAGE(
      BUILD "${CTEST_BINARY_DIRECTORY}"
      )
    IF (CTEST_DO_SUBMIT)
      TRIBITS_CTEST_SUBMIT( PARTS Coverage )
    ENDIF()

  ENDIF()

  #
  # F) Do memory testing
  #

  IF (NOT CTEST_DO_MEMORY_TESTING)
  
    MESSAGE("")
    MESSAGE("Skipping memory testing because CTEST_DO_MEMORY_TESTING='${CTEST_DO_MEMORY_TESTING}'!")
    MESSAGE("")

  ELSEIF (NOT AAO_CONFIGURE_PASSED)
  
    MESSAGE("")
    MESSAGE("Skipping memory tests because configure failed!")
    MESSAGE("")

  ELSEIF (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING AND AAO_CONFIGURE_PASSED)

    MESSAGE("Skipping memory testing because"
      " CTEST_DEPENDENCY_HANDLING_UNIT_TESTING='${CTEST_DEPENDENCY_HANDLING_UNIT_TESTING}'!"
      )

  ELSE()
    
    # NOTE: We always gather the memory results if the configure passed
    # independent if there was any build or test failures.  The memory stats
    # may not be very valid if there are build or test failures but there is
    # no harm and showing the memory based on tests that actually run (even
    # if they fail).

    MESSAGE("\nRunning memory tests ...\n")
    PRINT_VAR(CTEST_MEMORYCHECK_COMMAND)
    PRINT_VAR(CTEST_MEMORYCHECK_COMMAND_OPTIONS)
    PRINT_VAR(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE)
    CTEST_MEMCHECK(
      BUILD "${CTEST_BINARY_DIRECTORY}"
      )
    IF (CTEST_DO_SUBMIT)
      TRIBITS_CTEST_SUBMIT( PARTS MemCheck )
    ENDIF()

  ENDIF()

  #
  # G) Determine final pass/fail by gathering list of failing packages
  #

  IF (AAO_CONFIGURE_FAILED OR AAO_BUILD_FAILED OR AAO_INSTALL_FAILED)
    IF (${PROJECT_NAME}_ENABLE_ALL_PACKAGES)
      # Special value "ALL_PACAKGES" so that it will trigger enabling all
      # packages on the next CI iteration!
      SET(${PROJECT_NAME}_FAILED_PACKAGES  ALL_PACKAGES)
    ELSE()
      # Specific packages were selected to be tested so fail all of them!
      SET(${PROJECT_NAME}_FAILED_PACKAGES  ${${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST})
    ENDIF()
    # NOTE: With the all-at-once appraoch, there is no way to determine which
    # packages have build or install failures given the current CTEST_BUILD()
    # command.  And since some build targets don't get used in tests, we can't
    # look at what packages have test failures in order to know that a build
    # failure will cause a test failure.  And in the case of install failures,
    # those will never cause test failures.  Therefore, if there are any build
    # or install failures, we just have to assume that any tested package
    # could have failed.  Hense, we set the above just like for a (global)
    # configure failures.  Perhaps we couild read the generated *.xml files to
    # figure that out but that is not worth the work righ now.  The only bad
    # consequence of this is that a CI build would end up building and testing
    # every package even if only one dowstream package had a build failure,
    # for example.  That is just one of the downsides of the all-at-once
    # appraoch vs the package-by-package appraoch.
  ELSEIF (FAILED_TEST_LOG_FILE)
    TRIBITS_GET_FAILED_PACKAGES_FROM_FAILED_TESTS("${FAILED_TEST_LOG_FILE}"
      ${PROJECT_NAME}_FAILED_PACKAGES )
  ELSE()
    # If no tests failed, then there are no failed packages!
    SET(${PROJECT_NAME}_FAILED_PACKAGES)
  ENDIF()
  # ToDo: Optionally determine pass/fail based 

  MESSAGE("\nDone with the all-at-once configure, build, test, ans submit of ${PROJECT_NAME} packages!\n")

ENDMACRO()
