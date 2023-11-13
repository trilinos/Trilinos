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


include(TribitsGetCDashUrlsInsideCTestS)


# Wrapper used for unit testing purposes
#
macro(extrarepo_execute_process_wrapper)
  if (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    execute_process(${ARGN}
      RESULT_VARIABLE  EXTRAREPO_EXECUTE_PROCESS_WRAPPER_RTN_VAL)
    if (NOT EXTRAREPO_EXECUTE_PROCESS_WRAPPER_RTN_VAL STREQUAL "0")
      message(SEND_ERROR
        "Error: execute_process(${ARGN}) returned"
        " '${EXTRAREPO_EXECUTE_PROCESS_WRAPPER_RTN_VAL}'")
    endif()
  else()
    message("execute_process(${ARGN})")
  endif()
endmacro()


# Update an existing git repo
#
function(tribits_update_git_extrarepo  GIT_EXE  EXTRAREPO_SRC_DIR)

  set(EXTRAREPO_FETCH_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.fetch.out")
  set(EXTRAREPO_CLEAN_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.clean.out")
  set(EXTRAREPO_RESET_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.reset.out")
  set(EXTRAREPO_SET_BRANCH_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.set_branch.out")

  set(FETCH_CMND_ARGS
    COMMAND "${GIT_EXE}" fetch ${${PROJECT_NAME}_GIT_REPOSITORY_REMOTE}
    TIMEOUT 600 # seconds
    WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
    OUTPUT_FILE "${EXTRAREPO_FETCH_OUT_FILE}" )
  set(CLEAN_CMND_ARGS
    COMMAND "${GIT_EXE}" clean -fdx
    WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
    OUTPUT_FILE "${EXTRAREPO_CLEAN_OUT_FILE}" )
  set(RESET_CMND_ARGS
    COMMAND "${GIT_EXE}" reset --hard HEAD
    WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
    OUTPUT_FILE "${EXTRAREPO_RESET_OUT_FILE}" )
  if (${PROJECT_NAME}_EXTRAREPOS_BRANCH)
    set(SET_BRANCH_CMND_ARGS
      COMMAND "${GIT_EXE}" checkout -B ${${PROJECT_NAME}_EXTRAREPOS_BRANCH}
        --track ${${PROJECT_NAME}_GIT_REPOSITORY_REMOTE}/${${PROJECT_NAME}_EXTRAREPOS_BRANCH}
      WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
      OUTPUT_FILE "${EXTRAREPO_SET_BRANCH_OUT_FILE}" )
  ELSE ()
    set(SET_BRANCH_CMND_ARGS
      COMMAND "${GIT_EXE}" reset --hard "@{u}"
      WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
      OUTPUT_FILE "${EXTRAREPO_SET_BRANCH_OUT_FILE}" )
  endif()

  extrarepo_execute_process_wrapper(${FETCH_CMND_ARGS})
  extrarepo_execute_process_wrapper(${CLEAN_CMND_ARGS})
  extrarepo_execute_process_wrapper(${RESET_CMND_ARGS})
  extrarepo_execute_process_wrapper(${SET_BRANCH_CMND_ARGS})

endfunction()


# Update or clone a single extra repo
#
function(tribits_clone_or_update_extrarepo  EXTRAREPO_NAME_IN  EXTRAREPO_DIR_IN
  EXTRAREPO_REPOTYPE_IN  EXTRAREPO_REPOURL_IN
  )

  #message("TRIBITS_CLONE_OR_UPDATE_EXTRAREPO: ${EXTRAREPO_NAME_IN} ${EXTRAREPO_REPOURL_IN}")

  set(EXTRAREPO_SRC_DIR "${${PROJECT_NAME}_SOURCE_DIRECTORY}/${EXTRAREPO_DIR_IN}")
  #print_var(EXTRAREPO_SRC_DIR)

  set(EXTRAREPO_CLONE_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.clone.out")
  set(EXTRAREPO_CHECKOUT_OUT_FILE
    "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.checkout.out")

  if (NOT EXISTS "${EXTRAREPO_SRC_DIR}")

    message("\n${EXTRAREPO_NAME_IN}: Doing initial ${EXTRAREPO_REPOTYPE_IN}"
      " clone/checkout from URL '${EXTRAREPO_REPOURL_IN}' to dir '${EXTRAREPO_DIR_IN}' ...")

    # Set the command to clone
    if (${EXTRAREPO_REPOTYPE_IN} STREQUAL GIT)
      if (${PROJECT_NAME}_EXTRAREPOS_BRANCH) 
        set(CHECKOUT_BRANCH_ARG -b ${${PROJECT_NAME}_EXTRAREPOS_BRANCH})
      else()
        set(CHECKOUT_BRANCH_ARG)
      endif()
      set(CLONE_CMND_ARGS
        COMMAND "${GIT_EXECUTABLE}" clone
        ${CHECKOUT_BRANCH_ARG} -o ${${PROJECT_NAME}_GIT_REPOSITORY_REMOTE}
        "${EXTRAREPO_REPOURL}" ${EXTRAREPO_DIR_IN}
        WORKING_DIRECTORY "${${PROJECT_NAME}_SOURCE_DIRECTORY}"
        OUTPUT_FILE "${EXTRAREPO_CLONE_OUT_FILE}" )
    else()
      message(SEND_ERROR
        "Error, Invalid EXTRAREPO_REPOTYPE_IN='${EXTRAREPO_REPOTYPE_IN}'!")
    endif()

    # Do the clone
    extrarepo_execute_process_wrapper(${CLONE_CMND_ARGS})

  else()

    message("\n${EXTRAREPO_NAME_IN}: Doing ${EXTRAREPO_REPOTYPE_IN} update"
      " from URL '${EXTRAREPO_REPOURL_IN}' to dir '${EXTRAREPO_SRC_DIR}' ...")

  endif()

  if (${EXTRAREPO_REPOTYPE_IN} STREQUAL GIT)
    # Always update the git repo, even after a clone.  See
    # tribits_ctest_driver() documentation.
    tribits_update_git_extrarepo("${GIT_EXECUTABLE}" "${EXTRAREPO_SRC_DIR}")
  else()
    message(SEND_ERROR
      "Error, Invalid EXTRAREPO_REPOTYPE_IN='${EXTRAREPO_REPOTYPE_IN}'!")
  endif()

endfunction()


# Clone or update all of the extra repos and put them on the right branch.
#
# NOTE: The base repo is cloned by ctest_start() and updated by ctest_update()
# before calling this function.  This function only operates on the extra
# repos.
#
function(tribits_clone_or_update_extra_repos  CTEST_UPDATE_RETURN_VAL
  UPDATE_FAILED_VAR_OUT
  )

  set(UPDATE_FAILED FALSE)

  if (${PROJECT_NAME}_EXTRAREPOS_BRANCH)
    message("For extra repos, doing switch to branch ${${PROJECT_NAME}_EXTRAREPOS_BRANCH}")
  endif()

  set(EXTRAREPO_IDX 0)
  foreach(EXTRAREPO_NAME ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES})
    list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS ${EXTRAREPO_IDX}
      EXTRAREPO_DIR )
    list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES ${EXTRAREPO_IDX}
      EXTRAREPO_REPOTYPE )
    list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS ${EXTRAREPO_IDX}
      EXTRAREPO_REPOURL )
    tribits_clone_or_update_extrarepo( ${EXTRAREPO_NAME} ${EXTRAREPO_DIR}
      ${EXTRAREPO_REPOTYPE} ${EXTRAREPO_REPOURL} )
    # ToDo: Detect and return failure in cloning or updating extra repos!
    math(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")
  endforeach()

  set(${UPDATE_FAILED_VAR_OUT} ${UPDATE_FAILED} PARENT_SCOPE)

endfunction()


# Create the Updates.txt file
#
function(tribits_create_repo_updates_file)
  extrarepo_execute_process_wrapper(
    COMMAND ${PYTHON_EXECUTABLE}
      ${GITDIST_EXE} --dist-no-color
      log "--pretty=format:%h:  %s%nAuthor: %an <%ae>%nDate:   %ad%n"
      --name-status -C ORIG_HEAD..HEAD
    WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
    OUTPUT_FILE "${CTEST_BINARY_DIRECTORY}/Updates.txt"
    )
endfunction()


# Select the set of extra repositories
#
macro(tribits_setup_extrarepos)

  if (EXISTS "${${PROJECT_NAME}_EXTRAREPOS_FILE}" )
    # Repos many not already exist because we have not cloned them yet!
    set(${PROJECT_NAME}_CHECK_EXTRAREPOS_EXIST FALSE)
    tribits_get_and_process_extra_repositories_lists()
  else()
    message("${${PROJECT_NAME}_EXTRAREPOS_FILE} does not exist,"
       " skipping extra repositories.")
  endif()

endmacro()


# Select the list of packages
#
# OUTPUT: Sets ${PROJECT_NAME}_DEFAULT_PACKAGES
#
# NOTE: This macro is used to clean up the main tribits_ctest_driver()
# macro.
#
macro(tribits_setup_packages)

  include(TribitsPrintDependencyInfo)
  include(TribitsWriteXmlDependenciesFiles)

  # Here, we must point into the source tree just cloned (or updated)
  # and not the "driver" source dir tree for two reasons.  First, the
  # list of core packages may be more recent in what was checked out.
  # Second, the extra repos do not even exist in the "driver" source
  # tree.

  set(${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES  OFF)
  set(${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES  FALSE)
  if (CTEST_GENERATE_OUTER_DEPS_XML_OUTPUT_FILE)
    set(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE
       "${PROJECT_BINARY_DIR}/${${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME}")
  else()
    set(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
  endif()
  if (CTEST_SUBMIT_CDASH_SUBPROJECTS_DEPS_FILE)
    set(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE
      "${PROJECT_BINARY_DIR}/${${PROJECT_NAME}_CDASH_SUBPROJECT_DEPS_XML_FILE_NAME}" )
  else()
    set(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE)
  endif()
  set(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE)

  # Don't ignore missing repos by default.  This will allow processing to
  # continue but this outer CTest script will fail (thereby sending a CDash
  # email from the TDD system).  However, when we configure actual packages,
  # we do set this to TRUE so that the package configures will not fail due to
  # missing extra repositories.
  set_default_and_from_env(${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES FALSE)
  set_default_and_from_env(${PROJECT_NAME}_PRE_REPOSITORIES "")
  set_default_and_from_env(${PROJECT_NAME}_EXTRA_REPOSITORIES "")
  split("${${PROJECT_NAME}_PRE_REPOSITORIES}"  ","  ${PROJECT_NAME}_PRE_REPOSITORIES)
  split("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  ","  ${PROJECT_NAME}_EXTRA_REPOSITORIES)

  tribits_read_in_native_repositories()
  tribits_combine_native_and_extra_repos()
  tribits_read_all_project_deps_files_create_deps_graph()
  tribits_print_initial_dependency_info()
  tribits_write_xml_dependency_files()

  # When we get here, we will have the basic dependency structure set up
  # with only defaults set

  # Set this to "" so that it can be defined in enable_modified_packages_only()
  set(${PROJECT_NAME}_ENABLE_ALL_PACKAGES "")

endmacro()


macro(enable_package_if_not_explicitly_excluded  TRIBITS_PACKAGE)
  if ("${${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}}" STREQUAL "")
    message("Enabling explicitly set package ${TRIBITS_PACKAGE} ...")
    set(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} ON)
  elseif(NOT ${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
    if (${TRIBITS_PACKAGE}_EXPLICITY_EXCLUDED)
      message("NOT enabling explicitly set package ${TRIBITS_PACKAGE} since it was explicitly excluded!")
    else()
       message("Enabling explicitly set package ${TRIBITS_PACKAGE} which was default or otherwise disabed!")
      set(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} ON)
    endif()
  else()
    message("Explicitly set package ${TRIBITS_PACKAGE} is already enabled?")
  endif()
endmacro()


# Select packages set by the input
#
macro(enable_user_selected_packages)

  # 1) Set the enables for packages

  if (${PROJECT_NAME}_PACKAGE_ENABLES_FILE)
    message("Setting package enables specified in file"
      " '${${PROJECT_NAME}_PACKAGE_ENABLES_FILE}'")
    include(${${PROJECT_NAME}_PACKAGE_ENABLES_FILE})
  elseif (NOT "${${PROJECT_NAME}_PACKAGES_USER_SELECTED}" STREQUAL "")
    foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES_USER_SELECTED})
      enable_package_if_not_explicitly_excluded(${TRIBITS_PACKAGE})
    endforeach()
  else()
    message("Setting ${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON since"
      " ${PROJECT_NAME}_PACKAGES_USER_SELECTED='${${PROJECT_NAME}_PACKAGES_USER_SELECTED}'")
    set(${PROJECT_NAME}_ENABLE_ALL_PACKAGES ON)
  endif()

  # 2) Set extra package enables from ${PROJECT_NAME}_ADDITIONAL_PACKAGES

  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_ADDITIONAL_PACKAGES})
    enable_package_if_not_explicitly_excluded(${TRIBITS_PACKAGE})
  endforeach()

endmacro()


# Extract the list of changed files for the main repo on put into an
# modified files file.
#
macro(tribits_get_modified_files  WORKING_DIR_IN  MODIFIED_FILES_FILE_NAME_IN)
  set(CMND_ARGS
    COMMAND "${GIT_EXECUTABLE}" diff --name-only ORIG_HEAD..HEAD
    WORKING_DIRECTORY "${WORKING_DIR_IN}"
    OUTPUT_FILE ${MODIFIED_FILES_FILE_NAME_IN}
    #OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    execute_process(${CMND_ARGS})
  else()
    message("execute_process(${CMND_ARGS})")
  endif()
endmacro()


# Select only packages that are modified or failed in the last CI iteration
#
macro(enable_only_modified_packages)

  #
  # A) Get the list of changed packages
  #

  set(MODIFIED_FILES_FILE_NAME "${CTEST_BINARY_DIRECTORY}/modifiedFiles.txt")

  # A.1) Get changes from main ${PROJECT_NAME} repo

  tribits_get_modified_files("${CTEST_SOURCE_DIRECTORY}" "${MODIFIED_FILES_FILE_NAME}")

  # A.2) Get changes from extra repos

  set(EXTRAREPO_IDX 0)
  foreach(EXTRAREPO_NAME ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES})

    list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS
       ${EXTRAREPO_IDX} EXTRAREPO_DIR )
    list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS
      ${EXTRAREPO_IDX} EXTRAREPO_PACKSTAT )

    # For now, only look for changes if it has packages.  Later, we may need
    # to generalize this for the general extra repo case with deeper directory
    # and other VC systems than GIT.
    if (EXTRAREPO_PACKSTAT STREQUAL HASPACKAGES)

      set(EXTRAREPO_SRC_DIR "${CTEST_SOURCE_DIRECTORY}/${EXTRAREPO_DIR}")
      set(EXTRAREPO_MODIFIED_FILES_FILE_NAME
        "${CTEST_BINARY_DIRECTORY}/modifiedFiles.${EXTRAREPO_NAME}.txt")

      tribits_get_modified_files("${EXTRAREPO_SRC_DIR}"
        "${EXTRAREPO_MODIFIED_FILES_FILE_NAME}")

      file(STRINGS ${EXTRAREPO_MODIFIED_FILES_FILE_NAME} EXTRAREPO_MODIFIED_FILES_STR)
      set(EXTRAREPO_FILES_STR "")
      foreach(STR_LINE ${EXTRAREPO_MODIFIED_FILES_STR})
        string(APPEND EXTRAREPO_FILES_STR "${EXTRAREPO_DIR}/${STR_LINE}\n")
      endforeach()
      file(APPEND "${MODIFIED_FILES_FILE_NAME}" ${EXTRAREPO_FILES_STR})

    endif()

    math(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")

  endforeach()

  # A.3) Get the names of the modified packages

  if (NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR "Error, Python must be enabled to map from modified"
      " files to packages!")
  endif()

  if (EXISTS "${MODIFIED_FILES_FILE_NAME}")
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE}
        ${${PROJECT_NAME}_TRIBITS_DIR}/ci_support/get-tribits-packages-from-files-list.py
        --files-list-file=${MODIFIED_FILES_FILE_NAME}
        --project-dir=${TRIBITS_PROJECT_ROOT}
        --deps-xml-file=${CTEST_BINARY_DIRECTORY}/${${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME}
      OUTPUT_VARIABLE MODIFIED_PACKAGES_LIST
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  else()
    set(MODIFIED_PACKAGES_LIST)
  endif()

  split("${MODIFIED_PACKAGES_LIST}" "," MODIFIED_PACKAGES_LIST)
  print_var(MODIFIED_PACKAGES_LIST)

  #
  # B) Get the list of packages that failed last CI iteration
  #

  # NOTE: It is critical to enable and test packages until they pass.  If you
  # don't do this, then the package will not show as updated in the above
  # logic.  In this case only downstream packages will get enabled.  If the
  # failing packages break the downstream packages, this will be bad (for lots
  # of reasons).  Therefore, we must enable failing packages from the last CI
  # iteration and keep enabling and testing them until they do pass!

  if (EXISTS "${FAILED_PACKAGES_FILE_NAME}")
    file(READ "${FAILED_PACKAGES_FILE_NAME}" FAILING_PACKAGES_LIST)
    string(STRIP "${FAILING_PACKAGES_LIST}" FAILING_PACKAGES_LIST)
    print_var(FAILING_PACKAGES_LIST)
  endif()

  #
  # C) Enable the changed and previously failing packages
  #

  foreach(TRIBITS_PACKAGE ${MODIFIED_PACKAGES_LIST})
    #print_var(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
    assert_defined(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
    if ("${${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}}" STREQUAL "")
      if (
        ${TRIBITS_PACKAGE} STREQUAL "ALL_PACKAGES"
        OR
        ${TRIBITS_PACKAGE}_TESTGROUP STREQUAL "PT"
        OR
        (
          ${TRIBITS_PACKAGE}_TESTGROUP STREQUAL "ST"
           AND
           ${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE
           )
        )
        message("Enabling modified package: ${TRIBITS_PACKAGE}")
        set(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} ON)
      else()
        message("NOT enabling modified ST package: ${TRIBITS_PACKAGE}")
      endif()
    else()
      message("Not enabling explicitly disabled modified package: ${TRIBITS_PACKAGE}")
    endif()
  endforeach()

  if (FAILING_PACKAGES_LIST STREQUAL "ALL_PACKAGES")
    message("Enabling previously failing ALL_PACKAGES")
    set(${PROJECT_NAME}_ENABLE_ALL_PACKAGES ON)
  else()
    foreach(TRIBITS_PACKAGE ${FAILING_PACKAGES_LIST})
      if ("${${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}}" STREQUAL "")
        if (
          ${TRIBITS_PACKAGE}_TESTGROUP STREQUAL "PT"
          OR
          (
            ${TRIBITS_PACKAGE}_TESTGROUP STREQUAL "ST"
             AND
             ${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE
             )
          )
          message("Enabling previously failing package: ${TRIBITS_PACKAGE}")
          set(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} ON)
        else()
          message("NOT enabling previously failing ST package: ${TRIBITS_PACKAGE}")
        endif()
      else()
        message("Not enabling explicitly disabled previously"
          " failing package: ${TRIBITS_PACKAGE}")
      endif()
    endforeach()
  endif()

  #
  # D) Print the final status
  #

  if (${PROJECT_NAME}_ENABLE_ALL_PACKAGES)
    if (NOT ${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE)
      message(FATAL_ERROR
        "Error, failing 'ALL_PACKAGES' only allowed with all-at-once mode!")
    endif()
    message("\nDirectly modified or failing non-disabled packages that need"
      " to be tested:  ALL_PACKAGES")
  else()
    tribits_print_package_list_enable_status(
      "\nDirectly modified or failing non-disabled packages that need to be tested"
      INTERNAL ON NONEMPTY )
  endif()

endmacro()


# Exclude disabled packages from ${PROJECT_NAME}_EXCLUDE_PACKAGES
#
# NOTE: These disables need to dominate over the above enables so this code is
# after all the enable code has run
#
macro(disable_excluded_packages)
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_EXCLUDE_PACKAGES})
    message("Disabling excluded package ${TRIBITS_PACKAGE} ...")
    set(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} OFF)
    set(${TRIBITS_PACKAGE}_EXPLICITY_EXCLUDED TRUE)
  endforeach()
endmacro()


# Remove packages that are only implicitly enabled but don't have tests
# enabled.
#
macro(select_final_set_of_packages_to_directly_test)

  set(${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST)

  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES})

    set(PROCESS_THE_PACKAGE FALSE)

    if (${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}
      AND ${TRIBITS_PACKAGE}_ENABLE_TESTS
      )
      set(PROCESS_THE_PACKAGE  TRUE)
    elseif (${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}
      AND CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES
      )
      set(PROCESS_THE_PACKAGE  TRUE)
    endif()

    if(PROCESS_THE_PACKAGE)
      append_set(${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST  ${TRIBITS_PACKAGE})
    endif()

  endforeach()

endmacro()


# Set mapping of labels to subprojects (i.e. TriBITS packages) for CDash.
#
# NOTE: Unlike for the inner CMake configure, only subprojects that are
# explicitly tested will be marked as a CDash subproject.  This limits the
# rows in CDash.  This does not seem to be a problem for when running ctest
# locally.  When run locally, ctest will just report aggregated times for
# subprojects that have 1 or more tests.  Not true for CDash.
#
macro(tribits_ctest_driver_set_labels_to_subprojects_mapping)
  set(CTEST_LABELS_FOR_SUBPROJECTS)
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST})
    list(APPEND CTEST_LABELS_FOR_SUBPROJECTS ${TRIBITS_PACKAGE})
  endforeach()
endmacro()


# Select the default generator.
#
macro(select_default_generator)
  # When the build tree is known and exists, use
  # its generator.
  set(DEFAULT_GENERATOR "DID NOT SET!")
  if(EXISTS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
    file(STRINGS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt"
      line REGEX "^CMAKE_GENERATOR:" LIMIT_COUNT 1)
    if("${line}" MATCHES "=(.+)$")
      set(DEFAULT_GENERATOR "${CMAKE_MATCH_1}")
    endif()
  else()
    set(DEFAULT_GENERATOR "Unix Makefiles")
  endif()
endmacro()


# Call INITIALIZE_ERROR_QUEUE once at the top of TRIBITS_CTEST_DRIVER
#
macro(initialize_error_queue)
  set(TRIBITS_CTEST_DRIVER_ERROR_QUEUE "")
endmacro()


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
macro(queue_error err_msg)
  set(TRIBITS_CTEST_DRIVER_ERROR_QUEUE
    ${TRIBITS_CTEST_DRIVER_ERROR_QUEUE} "${err_msg}")
endmacro()


# Call report_queued_errors() once at the bottom of tribits_ctest_driver()
#
macro(report_queued_errors)
  if ("${TRIBITS_CTEST_DRIVER_ERROR_QUEUE}" STREQUAL "")
    message("TRIBITS_CTEST_DRIVER_ERROR_QUEUE is empty. All is well.")
  else()
    message("ERROR: TRIBITS_CTEST_DRIVER_ERROR_QUEUE reports the following error message queue:")
    foreach(err_msg ${TRIBITS_CTEST_DRIVER_ERROR_QUEUE})
      message("${err_msg}")
    endforeach()
  endif()
endmacro()


# Setup for tracking if a configure is being attempted to keep memory if it
# will pass or not.
#
# This will wrilte files in the directory ${CTEST_BINARY_DIRECTORY} to keep
# track of this across multiple ctest -S script invocations.
#
macro(tribits_remember_if_configure_attempted)

  # Must always define these files names as they they are used in functions
  # called later in the same ctest -S invocation!
  set(CONFIGURE_ATTEMPTED_FILE
    "${CTEST_BINARY_DIRECTORY}/ConfigureAttempted.txt")
  set(CONFIGURE_PASSED_FILE
    "${CTEST_BINARY_DIRECTORY}/ConfigurePasssed.txt")

  if (CTEST_DO_CONFIGURE)
    file(WRITE "${CONFIGURE_ATTEMPTED_FILE}" "Attempting configure")
    if (EXISTS "${CONFIGURE_PASSED_FILE}")
      file(REMOVE "${CONFIGURE_PASSED_FILE}")
    endif()
  elseif(CTEST_DO_NEW_START)
    if (EXISTS "${CONFIGURE_ATTEMPTED_FILE}")
      file(REMOVE "${CONFIGURE_ATTEMPTED_FILE}")
    endif()
    if (EXISTS "${CONFIGURE_PASSED_FILE}")
      file(REMOVE "${CONFIGURE_PASSED_FILE}")
    endif()
  endif()

endmacro()
# NOTE: Above, this is made a macro because it defines the vars
# CONFIGURE_ATTEMPTED_FILE and CONFIGURE_PASSED_FILE at the top function
# scope.  This is needed so the below functions will see them set.


# Determine if a past configure was attempted but did not pass
#
function(tribits_previous_configure_attempted_but_not_passsed
  PREVIOUS_CONFIGURE_ATTEMPTED_BUT_NOT_PASSSED_VAR_OUT
  )

  #print_var(CONFIGURE_ATTEMPTED_FILE)
  #print_var(CONFIGURE_PASSED_FILE)

  if(
    (EXISTS "${CONFIGURE_ATTEMPTED_FILE}")
    AND
    (NOT EXISTS "${CONFIGURE_PASSED_FILE}")
    )
    set(PREVIOUS_CONFIGURE_ATTEMPTED_BUT_NOT_PASSSED TRUE)
  else()
    set(PREVIOUS_CONFIGURE_ATTEMPTED_BUT_NOT_PASSSED FALSE)
  endif()

  set(${PREVIOUS_CONFIGURE_ATTEMPTED_BUT_NOT_PASSSED_VAR_OUT}
    ${PREVIOUS_CONFIGURE_ATTEMPTED_BUT_NOT_PASSSED} PARENT_SCOPE)

endfunction()


# Remember that the configure passed for later ctest -S invocations
#
function(tribits_remember_configure_passed)
  file(WRITE "${CONFIGURE_PASSED_FILE}" "Configure Passed!")
endfunction()


# Override CTEST_SUBMIT to drive multiple submits and to detect failed
# submissions and track them as queued errors.
#
macro(tribits_ctest_submit)

  # Cache the original CTEST_DROP_SITE and CTEST_DROP_LOCATION
  if ("${TRIBITS_CTEST_DROP_SITE_ORIG}" STREQUAL "")
    set(TRIBITS_CTEST_DROP_SITE_ORIG ${CTEST_DROP_SITE})
    if (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
      print_var(TRIBITS_CTEST_DROP_SITE_ORIG)
    endif()
  endif()
  if ("${TRIBITS_CTEST_DROP_LOCATION_ORIG}" STREQUAL "")
    set(TRIBITS_CTEST_DROP_LOCATION_ORIG ${CTEST_DROP_LOCATION})
    if (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
      print_var(TRIBITS_CTEST_DROP_LOCATION_ORIG)
    endif()
  endif()

  # Do the first submit
  set(CTEST_DROP_SITE ${TRIBITS_CTEST_DROP_SITE_ORIG})
  set(CTEST_DROP_LOCATION ${TRIBITS_CTEST_DROP_LOCATION_ORIG})
  if (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
    print_var(CTEST_DROP_SITE)
    print_var(CTEST_DROP_LOCATION)
  endif()

  tribits_ctest_submit_driver(${ARGN})

  # Do the second submit if requested!
  if (TRIBITS_2ND_CTEST_DROP_SITE OR TRIBITS_2ND_CTEST_DROP_LOCATION)

    message("\nDoing submit to second CDash site ...\n")

    if (NOT "${TRIBITS_2ND_CTEST_DROP_SITE}" STREQUAL "")
      if (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
        print_var(TRIBITS_2ND_CTEST_DROP_SITE)
      endif()
      set(CTEST_DROP_SITE ${TRIBITS_2ND_CTEST_DROP_SITE})
    endif()

    if (NOT "${TRIBITS_2ND_CTEST_DROP_LOCATION}" STREQUAL "")
      if (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
        print_var(TRIBITS_2ND_CTEST_DROP_LOCATION)
      endif()
      set(CTEST_DROP_LOCATION ${TRIBITS_2ND_CTEST_DROP_LOCATION})
    endif()

    tribits_ctest_submit_driver(${ARGN})

  endif()

endmacro()


macro(tribits_ctest_submit_driver)

  # If using a recent enough ctest with RETRY_COUNT, use it to overcome
  # failed submits:
  set(retry_args "")
  set(retry_args
    RETRY_COUNT ${CTEST_SUBMIT_RETRY_COUNT}
    RETRY_DELAY ${CTEST_SUBMIT_RETRY_DELAY})
  message("info: using retry_args='${retry_args}' for _ctest_submit call")

  # Call the original CTEST_SUBMIT and pay attention to its RETURN_VALUE:
  ctest_submit(${ARGN} ${retry_args} RETURN_VALUE rv)

  if(NOT "${rv}" STREQUAL "0")
    queue_error("error: ctest_submit failed: rv='${rv}' ARGN='${ARGN}' retry_args='${retry_args}'")
  endif()

endmacro()


# Wrapper for ctest_update(...) for unit testing
#
macro(ctest_update_wrapper)
  if (NOT CTEST_UPDATE_UNIT_TESTING_MODE)
    ctest_update(${ARGN})
  else()
    message("ctest_update(${ARGN})")
    set(UPDATE_RETURN_VAL ${CTEST_UPDATE_RETURN_VAL})
  endif()
endmacro()


# Helper macros to pass through common CMake configure arguments used by both
# package-by-package approach and all-at-once approach
#
macro(tribits_fwd_cmake_config_args_0)
  set( CONFIGURE_OPTIONS
    "-D${PROJECT_NAME}_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}"
    "-DCTEST_USE_LAUNCHERS:BOOL=${CTEST_USE_LAUNCHERS}"
    "-D${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON"
"-D${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS:STRING=${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS}"
    "-D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON"
    "-D${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=${${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES}"
    )
  if (NOT CTEST_GENERATE_DEPS_XML_OUTPUT_FILE)
    list(APPEND CONFIGURE_OPTIONS
    "-D${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE:FILEPATH=")
  endif()
  if (NOT "${${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES}" STREQUAL "")
    list(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES:BOOL=${${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES}")
  endif()
  if (NOT "${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE}" STREQUAL "")
    list(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE:BOOL=${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE}")
  endif()
  if (NOT MPI_EXEC_MAX_NUMPROCS STREQUAL 0)
    list(APPEND CONFIGURE_OPTIONS
      "-DMPI_EXEC_MAX_NUMPROCS:STRING=${MPI_EXEC_MAX_NUMPROCS}")
  endif()
  if (${PROJECT_NAME}_SKIP_CTEST_ADD_TEST)
    list(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_SKIP_CTEST_ADD_TEST:BOOL=${${PROJECT_NAME}_SKIP_CTEST_ADD_TEST}")
  endif()
  if (CTEST_DO_COVERAGE_TESTING)
    list(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_ENABLE_COVERAGE_TESTING:BOOL=ON")
  endif()
  if (${PROJECT_NAME}_EXTRAREPOS_FILE STREQUAL "NONE")
    set(EXTRAREOS_FILE_PASSED "")
  else()
    set(EXTRAREOS_FILE_PASSED "${${PROJECT_NAME}_EXTRAREPOS_FILE}")
  endif()
  list(APPEND CONFIGURE_OPTIONS
    "-D${PROJECT_NAME}_EXTRAREPOS_FILE:STRING=${EXTRAREOS_FILE_PASSED}")
  list(APPEND CONFIGURE_OPTIONS # See TRIBITS_SETUP_PACKAGES
    "-D${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES:BOOL=ON")
  list(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=${${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE}")
  if (CTEST_DO_INSTALL)
    list(APPEND CONFIGURE_OPTIONS
      "-DCMAKE_SKIP_INSTALL_ALL_DEPENDENCY=ON")
  endif()
endmacro()


macro(tribits_fwd_cmake_config_args_1)
  set(CONFIGURE_OPTIONS ${CONFIGURE_OPTIONS}
    ${EXTRA_SYSTEM_CONFIGURE_OPTIONS} ${EXTRA_CONFIGURE_OPTIONS}
    ${${PROJECT_NAME}_EXTRA_CONFIGURE_OPTIONS} )
endmacro()


# Remove the all of the LastTestsFailed*.log files so we can determine if any
# tests have failed.
#
macro(tribits_remove_last_test_failed_log_file)
  # Remove the LastTestsFailed log so we can detect if there are any
  # failed tests.
  set(TEST_TMP_DIR "${CTEST_BINARY_DIRECTORY}/Testing/Temporary")
  set(LAST_TESTS_FILED_LOG_FILE_GLOB "${TEST_TMP_DIR}/LastTestsFailed*.log")
  file(GLOB logfiles "${LAST_TESTS_FILED_LOG_FILE_GLOB}")
  foreach(logfile ${logfiles})
    file(REMOVE "${logfile}")
  endforeach()
endmacro()


# Sets the var FAILED_TEST_LOG_FILE if the file is found
macro(tribits_find_last_test_failed_log_file)
  file(GLOB FAILED_TEST_LOG_FILE "${LAST_TESTS_FILED_LOG_FILE_GLOB}")
endmacro()


# Get names of failed packages from failed tests
function(tribits_get_failed_packages_from_failed_tests
   LAST_TESTS_FAILED_FILE  FAILED_PACKAGES_OUT
   )
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE}
      "${${PROJECT_NAME}_TRIBITS_DIR}/ci_support/get-tribits-packages-from-last-tests-failed.py"
      "--deps-xml-file=${CTEST_BINARY_DIRECTORY}/${${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME}"
      "--last-tests-failed-file=${LAST_TESTS_FAILED_FILE}"
          OUTPUT_VARIABLE  FAILED_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  split("${FAILED_PACKAGES}" "," FAILED_PACKAGES)
  set(${FAILED_PACKAGES_OUT} "${FAILED_PACKAGES}" PARENT_SCOPE)
endfunction()


# Drive the configure, build, test, and submit package-by-package
#
# Sets ${PROJECT_NAME}_FAILED_PACKAGES as an indication if there are any
# failures.
#
macro(tribits_ctest_package_by_package)

  message(
    "\n***"
    "\n*** Loop through ${PROJECT_NAME} packages to configure, build, and test ..."
    "\n***")

  set(${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE)
  set(${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES)
  set(PACKAGE_IDX 0)

  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST})

    message("")
    message("${PACKAGE_IDX}) Processing current package ${TRIBITS_PACKAGE}:"
      " libs='${${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}}',"
      " tests='${${TRIBITS_PACKAGE}_ENABLE_TESTS}'")
    message("")

    set_property(GLOBAL PROPERTY SubProject ${TRIBITS_PACKAGE})
    set_property(GLOBAL PROPERTY Label ${TRIBITS_PACKAGE})

    #
    # A) Configure the package and its dependent packages
    #

    message("Configuring TRIBITS_PACKAGE='${TRIBITS_PACKAGE}'")

    # Create CONFIGURE_OPTIONS for this TRIBITS_PACKAGE
    tribits_fwd_cmake_config_args_0()
    list(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_ENABLE_TESTS:BOOL=${${TRIBITS_PACKAGE}_ENABLE_TESTS}")
    if (DEFINED ${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE)
      list(APPEND CONFIGURE_OPTIONS
        "-D${PROJECT_NAME}_ENABLE_${${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE}:BOOL=")
      set(${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE)
    endif()
    list(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_DEFINE_MISSING_PACKAGE_LIBS_TARGETS=ON")
    foreach(FAILED_PACKAGE ${${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES})
      list(APPEND CONFIGURE_OPTIONS
        "-D${PROJECT_NAME}_ENABLE_${FAILED_PACKAGE}:BOOL=OFF")
    endforeach()
    tribits_fwd_cmake_config_args_1()
    list(APPEND CONFIGURE_OPTIONS # Package enable must be at the very end to override other stuff!
       "-D${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}:BOOL=ON" )
    message("\nCONFIGURE_OPTIONS = '${CONFIGURE_OPTIONS}'")

    # Remember this package so we can set its enable to "" next time
    set(${PROJECT_NAME}_LAST_CONFIGURED_PACKAGE "${TRIBITS_PACKAGE}")

    #
    # B) Configure the package and its dependent packages
    #

    set(PBP_CONFIGURE_PASSED TRUE)

    if (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)

      message("${TRIBITS_PACKAGE}: Skipping configure due"
        " to running in unit testing mode!")

    else()

      #
      # We always have to configure if we are going to do anything for the
      # package.  We just want submit configure results to CDash if we are not
      # asked to configure!
      #

      set(PBP_CONFIGURE_PASSED FALSE)

      ctest_configure(
        BUILD "${CTEST_BINARY_DIRECTORY}"
        OPTIONS "${CONFIGURE_OPTIONS}" # New option!
        RETURN_VALUE CONFIGURE_RETURN_VAL
        )

      message("Generating the file '${CMAKE_CACHE_CLEAN_FILE}' ...")
      tribits_strip_comments_from_cmake_cache_file(
        "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt"
        "${CMAKE_CACHE_CLEAN_FILE}"
        )

      # If the configure failed, add the package to the list
      # of failed packages
      if ("${CONFIGURE_RETURN_VAL}" EQUAL "0")
        message("\n${TRIBITS_PACKAGE}: Configure passed!\n")
        set(PBP_CONFIGURE_PASSED TRUE)
        # load target properties and test keywords
        ctest_read_custom_files(BUILD "${CTEST_BINARY_DIRECTORY}")
        # Overridde from this file!
        include("${TRIBITS_PROJECT_ROOT}/CTestConfig.cmake")
      else()
        message("\n${TRIBITS_PACKAGE} FAILED to configure!\n")
      endif()

      if (EXISTS ${CMAKE_CACHE_CLEAN_FILE})
        set(CTEST_NOTES_FILES "${CTEST_NOTES_FILES_WO_CACHE};${CMAKE_CACHE_CLEAN_FILE}")
      else()
        set(CTEST_NOTES_FILES "${CTEST_NOTES_FILES_WO_CACHE}")
      endif()

      print_var(CTEST_NOTES_FILES)

      if (NOT CTEST_DO_CONFIGURE AND CTEST_DO_SUBMIT)
        message("${TRIBITS_PACKAGE}: Skipping submitting configure"
          " and notes due to CTEST_DO_CONFIGURE='${CTEST_DO_CONFIGURE}'!")
      elseif (CTEST_DO_SUBMIT)
        message("\nSubmitting configure and notes ...")
        tribits_ctest_submit( PARTS configure notes )
      endif()

    endif()

    # Print out values read from project CTestCustom.cmake file!
    print_var(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE)
    print_var(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE)

    #
    # C) Build the library and then ALL
    #

    set(PBP_BUILD_PASSED TRUE)
    set(PBP_BUILD_LIBS_PASSED TRUE)

    print_var(PBP_CONFIGURE_PASSED)

    if ( NOT PBP_CONFIGURE_PASSED AND CTEST_DO_BUILD )

      message("\n${TRIBITS_PACKAGE}: Skipping build due"
        " to configure failing!")

      set(PBP_BUILD_PASSED FALSE)
      set(PBP_BUILD_LIBS_PASSED FALSE)

    elseif (NOT CTEST_DO_BUILD)

      message("\n${TRIBITS_PACKAGE}: Skipping build due"
        " to CTEST_DO_BUILD='${CTEST_DO_BUILD}'!")

    elseif (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING OR
      CTEST_CONFIGURATION_UNIT_TESTING
      )

      message("\n${TRIBITS_PACKAGE}: Skipping build due"
        " to running in unit testing mode!")

    else()

      # Start by trying to build just the libraries for the current package

      set( CTEST_BUILD_TARGET ${TRIBITS_PACKAGE}_libs )
      message("\nBuilding target: '${CTEST_BUILD_TARGET}' ...\n")
      set(PBP_BUILD_LIBS_PASSED FALSE)
      ctest_build(
        BUILD "${CTEST_BINARY_DIRECTORY}"
        RETURN_VALUE  BUILD_LIBS_RETURN_VAL
        NUMBER_ERRORS  BUILD_LIBS_NUM_ERRORS
        APPEND
        )
      message("Build return: RETURN_VALUE=${BUILD_LIBS_RETURN_VAL},"
        " NUMBER_ERRORS=${BUILD_LIBS_NUM_ERRORS}")

      # Determine if the build failed or not.

      if ("${BUILD_LIBS_NUM_ERRORS}" EQUAL "0")
        message("\n${TRIBITS_PACKAGE}: Libs build passed!")
        set(PBP_BUILD_LIBS_PASSED TRUE)
      else()
        message("\nFAILED library build for package '${TRIBITS_PACKAGE}'!")
        set(PBP_BUILD_PASSED FALSE)
      endif()
      # Above: Since make -i is used BUILD_LIBS_RETURN_VAL might be 0, but
      # if there are errors the build should fail, so both
      # BUILD_LIBS_RETURN_VAL and BUILD_LIBS_NUM_ERRORS should be 0 for a
      # good build and for the all target to be built.

      # Submit the library build results to the dashboard
      if (CTEST_DO_SUBMIT)
        tribits_ctest_submit( PARTS build )
      endif()

      # If the build of the libraries passed, then go on the build
      # the tests/examples and run them.

      if (PBP_BUILD_LIBS_PASSED)

        # Build the ALL target, but append the results to the last build.xml
        set(CTEST_BUILD_TARGET)
        message("\nBuild ALL target for '${TRIBITS_PACKAGE}' ...\n")
        ctest_build(
          BUILD "${CTEST_BINARY_DIRECTORY}"
          RETURN_VALUE  BUILD_ALL_RETURN_VAL
          NUMBER_ERRORS  BUILD_ALL_NUM_ERRORS
          APPEND
          )
        message("Build all: BUILD_ALL_NUM_ERRORS='${BUILD_ALL_NUM_ERRORS}',"
          "BUILD_ALL_RETURN_VAL='${BUILD_ALL_RETURN_VAL}'" )

        if (NOT "${BUILD_ALL_NUM_ERRORS}" EQUAL "0")
          message("${TRIBITS_PACKAGE}: All build FAILED!")
          set(PBP_BUILD_PASSED FALSE)
        else()
          message("${TRIBITS_PACKAGE}: All build passed!")
        endif()

        # Submit the build for all target
        if (CTEST_DO_SUBMIT)
          tribits_ctest_submit( PARTS build )
        endif()

      endif()

    endif()

    #
    # D) Run the tests
    #

    set(PBP_TESTS_PASSED TRUE)

    if (NOT PBP_BUILD_LIBS_PASSED AND CTEST_DO_TEST)

      message("\n${TRIBITS_PACKAGE}: Skipping tests since library build failed!\n")

      set(PBP_TESTS_PASSED FALSE)

    elseif (NOT CTEST_DO_TEST)

      message("\n${TRIBITS_PACKAGE}: Skipping running tests due"
        " to CTEST_DO_TEST='${CTEST_DO_TEST}'!")

    else()
      
      #
      # D.1) Run the regular tests
      #

      set(PBP_TESTS_PASSED FALSE)

      # Run the tests that match the ${TRIBITS_PACKAGE} name
      message("\nRunning test for package '${TRIBITS_PACKAGE}'"
        " (parallel level ${CTEST_PARALLEL_LEVEL}) ...\n")
      tribits_remove_last_test_failed_log_file()
      ctest_test(
        BUILD "${CTEST_BINARY_DIRECTORY}"
        PARALLEL_LEVEL "${CTEST_PARALLEL_LEVEL}"
        INCLUDE_LABEL "^${TRIBITS_PACKAGE}$"
          )
      # See if a 'LastTestsFailed*.log' file exists to determine if there are
      # failed tests
      tribits_find_last_test_failed_log_file()
      if (FAILED_TEST_LOG_FILE)
        message("\n${TRIBITS_PACKAGE}: File '${FAILED_TEST_LOG_FILE}'"
          " exists so there were failed tests!")
      else()
        message("\n${TRIBITS_PACKAGE}: File '${FAILED_TEST_LOG_FILE}'"
          " does NOT exist so all tests passed!")
        set(PBP_TESTS_PASSED TRUE)
      endif()
      # 2009/12/05: ToDo: We need to add an argument to ctest_test(...)
      # called something like 'NUMBER_FAILED numFailedTests' to allow us to
      # detect when the tests have filed.
      if (CTEST_DO_SUBMIT)
        tribits_ctest_submit( PARTS Test )
      endif()

      #
      # D.2) Collect coverage results
      #

      if (CTEST_DO_COVERAGE_TESTING)

        message("\nRunning coverage for package '${TRIBITS_PACKAGE}' ...\n")

        ctest_coverage(
          BUILD "${CTEST_BINARY_DIRECTORY}"
          LABELS ${TRIBITS_PACKAGE} ${TRIBITS_PACKAGE}Libs ${TRIBITS_PACKAGE}Exes
          )

        if (CTEST_DO_SUBMIT)
          tribits_ctest_submit( PARTS Coverage )
        endif()

      endif()

    endif()

    #
    # E) Run memory testing
    #

    if (NOT PBP_BUILD_LIBS_PASSED AND CTEST_DO_MEMORY_TESTING)

      message("\n${TRIBITS_PACKAGE}: Skipping running memory checking"
         "tests since library build failed!\n")

    elseif (NOT CTEST_DO_MEMORY_TESTING)

      message("\n${TRIBITS_PACKAGE}: Skipping running memory checking tests due"
        " to CTEST_DO_MEMORY_TESTING='${CTEST_DO_MEMORY_TESTING}'!")

    else()

      message("\nRunning memory testing for package '${TRIBITS_PACKAGE}' ...\n")

      print_var(CTEST_MEMORYCHECK_COMMAND)
      print_var(CTEST_MEMORYCHECK_COMMAND_OPTIONS)
      print_var(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE)

      ctest_memcheck(
        BUILD "${CTEST_BINARY_DIRECTORY}"
        PARALLEL_LEVEL "${CTEST_PARALLEL_LEVEL}"
        INCLUDE_LABEL "^${TRIBITS_PACKAGE}$"
        )
      # ToDo: Determine if memory testing passed or not and affect overall
      # pass/fail!

      if (CTEST_DO_SUBMIT)
        tribits_ctest_submit( PARTS MemCheck )
      endif()

    endif()

    #
    # F) Record if this package failed the build or any tests
    #

    if (NOT PBP_CONFIGURE_PASSED OR NOT PBP_BUILD_LIBS_PASSED)
      list(APPEND ${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES ${TRIBITS_PACKAGE})
    endif()

    if (NOT PBP_BUILD_PASSED OR NOT PBP_TESTS_PASSED)
      list(APPEND ${PROJECT_NAME}_FAILED_PACKAGES ${TRIBITS_PACKAGE})
    endif()

    #
    # G) Do submit of update
    #

    if (CTEST_DO_SUBMIT)
      message("\nSubmit the update file that will trigger the notification email ...\n")
      tribits_ctest_submit( PARTS update )
    endif()

    math(EXPR PACKAGE_IDX "${PACKAGE_IDX}+1")

  endforeach(TRIBITS_PACKAGE)

  if (${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES)
    message(
      "\nFinal set packages that failed to configure or have the libraries build:"
      " '${${PROJECT_NAME}_FAILED_LIB_BUILD_PACKAGES}'")
  endif()

  message("\nDone with the incremental building and testing of"
    " ${PROJECT_NAME} packages!\n")

endmacro()
# NOTE: Above, the option
# ${PROJECT_NAME}_DEFINE_MISSING_PACKAGE_LIBS_TARGETS=ON is passed down
# through to the inner CMake TriBITS configure to trigger the creation of
# dummy targets <PackageName>_libs for all the packages for the case where a
# package is disabled due to a disabled upstream package and
# ${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON but the target
# <thePackage>_libs is attempted to be built anyway and we expect it to build
# nothing and result in no error.  (The outer ctest -S driver is not smart
# enough to know all the lgoic for if a package will actually be enabled or
# not.  That is the job of the inner TriBITS dependency logic and
# ${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON.) Otherwise, with
# CMake 3.19+, cmake_build() catches errors in undefined global build targets
# like this and reports them correctly.  This workaround allows the
# package-by-package mode to gracefully disable downstream packages that can't
# be enabled due to the disable of a broken upstream packages.  See the test
# TriBITS_CTestDriver_PBP_ST_BreakConfigureRequiredPkg that exercises this use
# case.


# Drive the configure, build, test, and submit all at once for all of the
# enabled packages.
#
# Sets ${PROJECT_NAME}_FAILED_PACKAGES as an indication if there are any
# failures.
#
macro(tribits_ctest_all_at_once)

  message(
    "\n***"
    "\n*** Configure, build, test and submit results all-at-once for all enabled packages ..."
    "\n***")

  set(AAO_CONFIGURE_FAILED FALSE)
  set(AAO_BUILD_FAILED FALSE)
  set(AAO_INSTALL_FAILED FALSE)

  #
  # A) Define mapping from labels to subprojects and gather configure arguments
  #

  tribits_ctest_driver_set_labels_to_subprojects_mapping()
  print_var(CTEST_LABELS_FOR_SUBPROJECTS)

  message("")
  message("Configuring ...")
  message("")

  # Create CONFIGURE_OPTIONS
  tribits_fwd_cmake_config_args_0()
  if (NOT "${${PROJECT_NAME}_PACKAGE_ENABLES_FILE}" STREQUAL "")
    # NOTE: For now, the user is expected to pass through this file in the
    # inner CMake cache var ${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE!  We should
    # fix this in the future but that is what it is for now.
  elseif (${PROJECT_NAME}_ENABLE_ALL_PACKAGES)
    list(APPEND CONFIGURE_OPTIONS
      "-D${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON" )
    foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_EXCLUDE_PACKAGES})
      list(APPEND CONFIGURE_OPTIONS
        "-D${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}=OFF" )
    endforeach()
    # NOTE: Above we have to explicitly set disables for the excluded packages
    # since we are pssing in ${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON.  This is
    # effectively the "black-listing" approach.
  else()
    foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST})
      list(APPEND CONFIGURE_OPTIONS
         "-D${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}=ON" )
    endforeach()
    # NOTE: Above we don't have to consider the packages excluded in
    # ${PROJECT_NAME}_EXCLUDE_PACKAGES because they are not enabled ad this
    # point and therefore no in ${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST.
    # This is effectively the "white-listing" approach.
  endif()
  list(APPEND CONFIGURE_OPTIONS
    "-D${PROJECT_NAME}_ENABLE_TESTS:BOOL=${${PROJECT_NAME}_INNER_ENABLE_TESTS}")
  tribits_fwd_cmake_config_args_1()
  message("\nCONFIGURE_OPTIONS = '${CONFIGURE_OPTIONS}'")

  #
  # B) Configure the package and its dependent packages
  #

  tribits_previous_configure_attempted_but_not_passsed(
    PREVIOUS_CONFIGURE_ATTEMPTED_BUT_NOT_PASSSED)
  #print_var(PREVIOUS_CONFIGURE_ATTEMPTED_BUT_NOT_PASSSED)

  if ((NOT CTEST_DO_CONFIGURE) AND PREVIOUS_CONFIGURE_ATTEMPTED_BUT_NOT_PASSSED)

    message(
      "\nSkipping configure due to CTEST_DO_CONFIGURE='${CTEST_DO_CONFIGURE}'!\n"
      "\nHOWEVER: A configure was previously attempted but did not pass so consider configure FAILED!")
    set(AAO_CONFIGURE_PASSED FALSE)
    set(AAO_CONFIGURE_FAILED TRUE)

  elseif (NOT CTEST_DO_CONFIGURE)

    message("\nSkipping configure due to CTEST_DO_CONFIGURE='${CTEST_DO_CONFIGURE}'!\n")
    set(AAO_CONFIGURE_PASSED TRUE)
    # Just assume configure passeed for the purpose of running the build.

  elseif (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)

    message("Skipping actual ctest_configure() because"
      " CTEST_DEPENDENCY_HANDLING_UNIT_TESTING='${CTEST_DEPENDENCY_HANDLING_UNIT_TESTING}'!"
      )
    set(AAO_CONFIGURE_PASSED TRUE)

  else()

    ctest_configure(
      BUILD "${CTEST_BINARY_DIRECTORY}"
      OPTIONS "${CONFIGURE_OPTIONS}" # New option!
      RETURN_VALUE CONFIGURE_RETURN_VAL
      )
  
    message("Generating the file '${CMAKE_CACHE_CLEAN_FILE}' ...")
    tribits_strip_comments_from_cmake_cache_file(
      "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt"
      "${CMAKE_CACHE_CLEAN_FILE}"
      )
    
    if (NOT "${CONFIGURE_RETURN_VAL}" EQUAL "0")
      message("Configure FAILED!")
      set(AAO_CONFIGURE_PASSED FALSE)
      set(AAO_CONFIGURE_FAILED TRUE)
    else()
      message("Configure PASSED!")
      set(AAO_CONFIGURE_PASSED TRUE)
    endif()

    if (AAO_CONFIGURE_PASSED)
      tribits_remember_configure_passed()
    endif()
  
    set(CTEST_NOTES_FILES "${CTEST_NOTES_FILES_WO_CACHE}")
  
    if (EXISTS ${CMAKE_CACHE_CLEAN_FILE})
      list(APPEND CTEST_NOTES_FILES "${CMAKE_CACHE_CLEAN_FILE}")
    endif()
  
    if (EXISTS "${REPO_VERSION_FILE}")
      set(CTEST_NOTES_FILES "${REPO_VERSION_FILE};${CTEST_NOTES_FILES}")
    endif()
  
    print_var(CTEST_NOTES_FILES)
  
    # Submit configure results and the notes to the dashboard
    if (CTEST_DO_SUBMIT)
      message("\nSubmitting update, configure and notes ...")
      tribits_ctest_submit( PARTS update configure notes )
    endif()

  endif()

  # Read in configured CTestCustom.cmake
  ctest_read_custom_files(BUILD "${CTEST_BINARY_DIRECTORY}")
  # NOTE: Above, it is safe to call ctest_read_custom_files() even if the
  # configure failed and the file CTestCustom.cmake does exist.  In this case,
  # CTest will just do nothing.

  # Overridde any values by loading <projectDir>/CTestConfig.cmake
  include("${TRIBITS_PROJECT_ROOT}/CTestConfig.cmake")

  # Print out values read from project CTestCustom.cmake file
  print_var(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE)
  print_var(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE)

  #
  # C) Do the build
  #

  if (NOT CTEST_DO_BUILD)

    message("\nSkipping build due to CTEST_DO_BUILD='${CTEST_DO_BUILD}'!\n")

  elseif (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING AND AAO_CONFIGURE_PASSED)

    message("Skipping build because"
      " CTEST_DEPENDENCY_HANDLING_UNIT_TESTING='${CTEST_DEPENDENCY_HANDLING_UNIT_TESTING}'!"
      )

  elseif (AAO_CONFIGURE_PASSED)
  
    message("")
    message("Building all targets ...")
    message("")
  
    ctest_build(
      BUILD "${CTEST_BINARY_DIRECTORY}"
      RETURN_VALUE  BUILD_ALL_RETURN_VAL
      NUMBER_ERRORS  BUILD_ALL_NUM_ERRORS
      )
    message("Build output: BUILD_ALL_NUM_ERRORS='${BUILD_ALL_NUM_ERRORS}',"
      "BUILD_ALL_RETURN_VAL='${BUILD_ALL_RETURN_VAL}'" )
  
    if (NOT "${BUILD_ALL_NUM_ERRORS}" EQUAL "0")
      message("Build FAILED!")
      set(AAO_BUILD_FAILED TRUE)
    else()
      message("Build PASSED!")
    endif()
  
    # Submit the build for all target
    if (CTEST_DO_SUBMIT)
      tribits_ctest_submit( PARTS build )
    endif()

    if (CTEST_DO_INSTALL)

      message("")
      message("Installing (i.e. building target 'install_package_by_package') ...")
      message("")

      ctest_build(
        BUILD "${CTEST_BINARY_DIRECTORY}"
        TARGET install_package_by_package
        RETURN_VALUE  BUILD_INSTALL_RETURN_VAL
        NUMBER_ERRORS  BUILD_INSTALL_NUM_ERRORS
        )
      message("Build install output:"
        " BUILD_INSTALL_NUM_ERRORS='${BUILD_INSTALL_NUM_ERRORS}',"
        "BUILD_INSTALL_RETURN_VAL='${BUILD_INSTALL_RETURN_VAL}'" )

      if (NOT "${BUILD_INSTALL_NUM_ERRORS}" EQUAL "0")
        message("Install FAILED!")
        set(AAO_INSTALL_FAILED TRUE)
      else()
        message("Install PASSED!")
      endif()

      # Submit the build for all target
      if (CTEST_DO_SUBMIT)
        tribits_ctest_submit( PARTS build )
      endif()

    endif()

  else()
  
    message("")
    message("Skipping build because configure failed!")
    message("")
  
  endif()

  #
  # D) Run tests
  #

  if (NOT CTEST_DO_TEST)
  
    message("")
    message("Skipping tests because CTEST_DO_TEST='${CTEST_DO_TEST}'!")
    message("")

  elseif (NOT AAO_CONFIGURE_PASSED)
  
    message("")
    message("Skipping tests because configure failed!")
    message("")

  elseif (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING AND AAO_CONFIGURE_PASSED)

    message("Skipping testing because"
      " CTEST_DEPENDENCY_HANDLING_UNIT_TESTING='${CTEST_DEPENDENCY_HANDLING_UNIT_TESTING}'!"
      )

  else()

    # NOTE: We always run the tests if the configure passed no matter if there
    # are build failures because the only way that we can detect what packages
    # have build failures is to see what packages have test failures.

    tribits_remove_last_test_failed_log_file()

    # Run the tests
    message("")
    message("\nRunning tests (parallel level ${CTEST_PARALLEL_LEVEL}) ...\n")
    message("")

    ctest_test(
      BUILD "${CTEST_BINARY_DIRECTORY}"
      PARALLEL_LEVEL "${CTEST_PARALLEL_LEVEL}"
      )

    # See if a 'LastTestsFailed*.log' file exists to determine if there are
    # failed tests.
    tribits_find_last_test_failed_log_file()
    if (FAILED_TEST_LOG_FILE)
      message("File '${FAILED_TEST_LOG_FILE}' exists so there were non-passing tests!")
    else()
      message("File '${FAILED_TEST_LOG_FILE}' does NOT exist so all tests passed!")
    endif()

    if (CTEST_DO_SUBMIT)
      tribits_ctest_submit( PARTS Test )
    endif()

  endif()

  #
  # E) Gather coverage results
  #

  if (NOT CTEST_DO_COVERAGE_TESTING)
  
    message("")
    message("Skipping converage tests because CTEST_DO_COVERAGE_TESTING='${CTEST_DO_COVERAGE_TESTING}'!")
    message("")

  elseif (NOT AAO_CONFIGURE_PASSED)
  
    message("")
    message("Skipping coverage tests because configure failed!")
    message("")

  elseif (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING AND AAO_CONFIGURE_PASSED)

    message("Skipping coverage testing because"
      " CTEST_DEPENDENCY_HANDLING_UNIT_TESTING='${CTEST_DEPENDENCY_HANDLING_UNIT_TESTING}'!"
      )

  else()
    
    # NOTE: We always gather the coverage results if the configure passed
    # independent if there was any build or test failures.  The coverage stats
    # may not be very valid if there are build or test failures but there is
    # no harm and showing the coverage based on tests that actually run (even
    # if they fail).

    message("\nGathering coverage results ...\n")
    ctest_coverage(
      BUILD "${CTEST_BINARY_DIRECTORY}"
      )
    if (CTEST_DO_SUBMIT)
      tribits_ctest_submit( PARTS Coverage )
    endif()

  endif()

  #
  # F) Do memory testing
  #

  if (NOT CTEST_DO_MEMORY_TESTING)
  
    message("")
    message("Skipping memory testing because CTEST_DO_MEMORY_TESTING='${CTEST_DO_MEMORY_TESTING}'!")
    message("")

  elseif (NOT AAO_CONFIGURE_PASSED)
  
    message("")
    message("Skipping memory tests because configure failed!")
    message("")

  elseif (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING AND AAO_CONFIGURE_PASSED)

    message("Skipping memory testing because"
      " CTEST_DEPENDENCY_HANDLING_UNIT_TESTING='${CTEST_DEPENDENCY_HANDLING_UNIT_TESTING}'!"
      )

  else()
    
    # NOTE: We always gather the memory results if the configure passed
    # independent if there was any build or test failures.  The memory stats
    # may not be very valid if there are build or test failures but there is
    # no harm and showing the memory based on tests that actually run (even
    # if they fail).

    message("\nRunning memory tests ...\n")
    print_var(CTEST_MEMORYCHECK_COMMAND)
    print_var(CTEST_MEMORYCHECK_COMMAND_OPTIONS)
    print_var(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE)
    ctest_memcheck(
      BUILD "${CTEST_BINARY_DIRECTORY}"
      )
    if (CTEST_DO_SUBMIT)
      tribits_ctest_submit( PARTS MemCheck )
    endif()

  endif()

  #
  # G) Determine final pass/fail by gathering list of failing packages
  #

  if (AAO_CONFIGURE_FAILED OR AAO_BUILD_FAILED OR AAO_INSTALL_FAILED)
    if (${PROJECT_NAME}_ENABLE_ALL_PACKAGES)
      # Special value "ALL_PACKAGES" so that it will trigger enabling all
      # packages on the next CI iteration!
      set(${PROJECT_NAME}_FAILED_PACKAGES  ALL_PACKAGES)
    else()
      # Specific packages were selected to be tested so fail all of them!
      set(${PROJECT_NAME}_FAILED_PACKAGES  ${${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST})
    endif()
    # NOTE: With the all-at-once approach, there is no way to determine which
    # packages have build or install failures given the current ctest_build()
    # command.  And since some build targets don't get used in tests, we can't
    # look at what packages have test failures in order to know that a build
    # failure will cause a test failure.  And in the case of install failures,
    # those will never cause test failures.  Therefore, if there are any build
    # or install failures, we just have to assume that any tested package
    # could have failed.  Hence, we set the above just like for a (global)
    # configure failures.  Perhaps we could read the generated *.xml files to
    # figure that out but that is not worth the work right now.  The only bad
    # consequence of this is that a CI build would end up building and testing
    # every package even if only one downstream package had a build failure,
    # for example.  That is just one of the downsides of the all-at-once
    # approach vs. the package-by-package approach.
  elseif (FAILED_TEST_LOG_FILE)
    tribits_get_failed_packages_from_failed_tests("${FAILED_TEST_LOG_FILE}"
      ${PROJECT_NAME}_FAILED_PACKAGES )
  else()
    # If no tests failed, then there are no failed packages!
    set(${PROJECT_NAME}_FAILED_PACKAGES)
  endif()
  # ToDo: Optionally determine pass/fail based 

  message("\nDone with the all-at-once configure, build, test, and submit of ${PROJECT_NAME} packages!\n")

endmacro()
