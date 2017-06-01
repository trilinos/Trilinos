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
# Function for getting the tracking branch
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
    # For unit testing purpose, just return generic tacking branch string
    SET(TRACKING_BRANCH "tracking/branch")
  ENDIF()
  SET(${TRACKING_BRANCH_OUT}  ${TRACKING_BRANCH}  PARENT_SCOPE)
ENDFUNCTION()


#
# Update an existing git repo
#
FUNCTION(TRIBITS_UPDATE_GIT_EXTRAREPO  GIT_EXE  EXTRAREPO_SRC_DIR)

  SET(EXTRAREPO_CLEAN_OUT_FILE "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.clean.out")
  SET(EXTRAREPO_RESET_OUT_FILE "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.reset.out")
  SET(EXTRAREPO_FETCH_OUT_FILE "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.fetch.out")
  SET(EXTRAREPO_SET_BRANCH_OUT_FILE "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.set_branch.out")

  EXTRAREPO_GET_TRACKING_BRANCH("${EXTRAREPO_SRC_DIR}"
    EXTRAREPO_TRACKING_BRANCH)
  #PRINT_VAR(EXTRAREPO_TRACKING_BRANCH)
  SET(CLEAN_CMND_ARGS
    COMMAND "${GIT_EXE}" clean -fdx
    WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
    OUTPUT_FILE "${EXTRAREPO_CLEAN_OUT_FILE}" )
  SET(RESET_CMND_ARGS
    COMMAND "${GIT_EXE}" reset --hard HEAD
    WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
    OUTPUT_FILE "${EXTRAREPO_RESET_OUT_FILE}" )
  SET(FETCH_CMND_ARGS
    COMMAND "${GIT_EXE}" fetch origin
    TIMEOUT 600 # seconds
    WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
    OUTPUT_FILE "${EXTRAREPO_FETCH_OUT_FILE}" )
  IF (${PROJECT_NAME}_EXTRAREPOS_BRANCH)
    SET(SET_BRANCH_CMND_ARGS
      COMMAND "${GIT_EXE}" checkout -B ${${PROJECT_NAME}_EXTRAREPOS_BRANCH}
        --track origin/${${PROJECT_NAME}_EXTRAREPOS_BRANCH}
      WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
      OUTPUT_FILE "${EXTRAREPO_SET_BRANCH_OUT_FILE}" )
  ELSE ()
    SET(SET_BRANCH_CMND_ARGS
      COMMAND "${GIT_EXE}" reset --hard "${EXTRAREPO_TRACKING_BRANCH}"
      WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
      OUTPUT_FILE "${EXTRAREPO_SET_BRANCH_OUT_FILE}" )
  ENDIF()

  EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${CLEAN_CMND_ARGS})
  EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${RESET_CMND_ARGS})
  EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${FETCH_CMND_ARGS})
  EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${SET_BRANCH_CMND_ARGS})

ENDFUNCTION()



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
# Update or clone a single extra repo
#

FUNCTION(TRIBITS_CLONE_OR_UPDATE_EXTRAREPO  EXTRAREPO_NAME_IN  EXTRAREPO_DIR_IN
  EXTRAREPO_REPOTYPE_IN  EXTRAREPO_REPOURL_IN
  )

  #MESSAGE("TRIBITS_CLONE_OR_UPDATE_EXTRAREPO: ${EXTRAREPO_NAME_IN} ${EXTRAREPO_REPOURL_IN}")

  SET(EXTRAREPO_SRC_DIR "${${PROJECT_NAME}_SOURCE_DIRECTORY}/${EXTRAREPO_DIR_IN}")
  #PRINT_VAR(EXTRAREPO_SRC_DIR)

  SET(EXTRAREPO_CLONE_OUT_FILE "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.clone.out")
  SET(EXTRAREPO_CHECKOUT_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.checkout.out")

  IF (NOT EXISTS "${EXTRAREPO_SRC_DIR}")

    MESSAGE("\n${EXTRAREPO_NAME_IN}: Doing initial ${EXTRAREPO_REPOTYPE_IN}"
      " clone/checkout from URL '${EXTRAREPO_REPOURL_IN}' to dir '${EXTRAREPO_DIR_IN}' ...")

    # Set the command to clone
    IF (${EXTRAREPO_REPOTYPE_IN} STREQUAL GIT)
      SET(CLONE_CMND_ARGS
        COMMAND "${GIT_EXE}" clone "${EXTRAREPO_REPOURL}" ${EXTRAREPO_DIR_IN}
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
    TRIBITS_UPDATE_GIT_EXTRAREPO("${GIT_EXE}" "${EXTRAREPO_SRC_DIR}")
  ELSE()
    MESSAGE(SEND_ERROR
      "Error, Invalid EXTRAREPO_REPOTYPE_IN='${EXTRAREPO_REPOTYPE_IN}'!")
  ENDIF()

ENDFUNCTION()


#
# Update the branch of the base git repo
#
FUNCTION(TRIBITS_SET_BASE_REPO_BRANCH  CTEST_UPDATE_RETURN_VAL
  UPDATE_FAILED_VAR_OUT
  )

  SET(GIT_CHECKOUT_RETURN_VAL "0")

  IF (${PROJECT_NAME}_BRANCH AND NOT "${CTEST_UPDATE_RETURN_VAL}" LESS "0")

    MESSAGE("For base repo, doing switch to branch ${${PROJECT_NAME}_BRANCH}")

    SET(EXECUTE_PROCESS_COMMANDS_ARGS
      COMMAND ${GIT_EXE} checkout
        -B ${${PROJECT_NAME}_BRANCH} --track origin/${${PROJECT_NAME}_BRANCH}
      WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
      RESULT_VARIABLE GIT_CHECKOUT_RETURN_VAL
      OUTPUT_VARIABLE BRANCH_OUTPUT
      ERROR_VARIABLE  BRANCH_ERROR
      )
     # NOTE: Above will work smoothly even if the local branch already
     # exists and/or is already on that branch.  This command does not move
     # ORIG_HEAD so it will not mess up the pull and update that CTest did
     # for the base repo.

    IF (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
      EXECUTE_PROCESS(${EXECUTE_PROCESS_COMMANDS_ARGS})
    ELSE()
      MESSAGE("EXECUTE_PROCESS(${EXECUTE_PROCESS_COMMANDS_ARGS})")
      SET(GIT_CHECKOUT_RETURN_VAL 0)
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

  IF ("${CTEST_UPDATE_RETURN_VAL}" LESS "0" OR NOT "${GIT_CHECKOUT_RETURN_VAL}" EQUAL "0")
    SET(${UPDATE_FAILED_VAR_OUT} TRUE PARENT_SCOPE)
  ELSE()
    SET(${UPDATE_FAILED_VAR_OUT} FALSE PARENT_SCOPE)
  ENDIF()

ENDFUNCTION()


#
# Clone or update all of the repos and put them on right branch
#
# NOTE: The base repo is cloned and updated by CTEST_UPDATE() before calling
# this function.  This function only puts the base repo on the right branch.
#

FUNCTION(TRIBITS_CLONE_OR_UPDATE_ALL_REPOS  CTEST_UPDATE_RETURN_VAL
  UPDATE_FAILED_VAR_OUT
  )

  SET(UPDATE_FAILED FALSE)

  # A) Put the base repo on the right branch

  TRIBITS_SET_BASE_REPO_BRANCH(${CTEST_UPDATE_RETURN_VAL}  BASE_REPO_UPDATE_FAILED)
  IF (BASE_REPO_UPDATE_FAILED)
    SET(UPDATE_FAILED TRUE)
  ENDIF()

  # B) Clone and update the extra repos

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
    # ToDo: Detect and return failure
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
