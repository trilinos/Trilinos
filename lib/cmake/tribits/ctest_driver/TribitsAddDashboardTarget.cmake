# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


################################################################################
#
# This file gets included in the main TriBITS framework.  It is put here to
# reduce the size of the tribits/core/ directory.
#
################################################################################

include("${CMAKE_CURRENT_LIST_DIR}/../core/utils/TribitsGitRepoVersionInfo.cmake")

#
# Macro that drives a experimental 'dashboard' target
#

macro(tribits_add_dashboard_target)

  if (NOT (WIN32 AND NOT CYGWIN))

    if ("${${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE_DEFAULT}" STREQUAL "")
      set(${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE_DEFAULT FALSE)
    endif()
    advanced_set(${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE
      ${${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE_DEFAULT}
      CACHE BOOL
      "If set to TRUE, use all-at-once mode for configure, build, test, and submit.  Otherwise, use package-by-package mode.")

    advanced_set(${PROJECT_NAME}_DASHBOARD_CTEST_ARGS "-V" CACHE STRING
      "Extra arguments to pass to CTest when calling 'ctest -S' to run the 'dashboard' make target." )

    advanced_set(CTEST_BUILD_FLAGS "" CACHE STRING
      "Sets CTEST_BUILD_FLAGS on the env before invoking 'ctest -S'." )

    advanced_set(CTEST_PARALLEL_LEVEL "" CACHE STRING
      "Sets CTEST_PARALLEL_LEVEL on the env before invoking 'ctest -S'." )

    # H.1) Enable all packages that are enabled and have tests enabled

    set(${PROJECT_NAME}_ENABLED_PACKAGES_LIST)
    set(${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST)
    foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES})
      if (${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} AND ${TRIBITS_PACKAGE}_ENABLE_TESTS)
        if (${PROJECT_NAME}_ENABLED_PACKAGES_LIST)
          set(${PROJECT_NAME}_ENABLED_PACKAGES_LIST
            "${${PROJECT_NAME}_ENABLED_PACKAGES_LIST}\;${TRIBITS_PACKAGE}")
        else()
          set(${PROJECT_NAME}_ENABLED_PACKAGES_LIST "${TRIBITS_PACKAGE}")
        endif()
        set(${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST
          ${${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST} -D${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}=ON)
      endif()
    endforeach()
    #print_var(${PROJECT_NAME}_ENABLED_PACKAGES_LIST)

    set(EXPR_CMND_ARGS)

    # Hard override options used by basic build and tests
    append_set(EXPR_CMND_ARGS "TRIBITS_PROJECT_ROOT=${${PROJECT_NAME}_SOURCE_DIR}")
    append_set(EXPR_CMND_ARGS "${PROJECT_NAME}_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}")
    append_set(EXPR_CMND_ARGS "${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS='${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS}'")
    append_set(EXPR_CMND_ARGS "${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE=${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE}")

    # Determine if base repo is a git repo (by seeing if SHA1 can be extracted)
    tribits_git_repo_sha1("${PROJECT_SOURCE_DIR}" projectGitRepoSha1
      FAILURE_MESSAGE_OUT  projectGitRepoSha1FailureMsg)
    if (projectGitRepoSha1 STREQUAL "")
      append_set(EXPR_CMND_ARGS "CTEST_DO_UPDATES=OFF")
    endif()

    # Conditionally override options used only for the 'dashboard' target.
    # These options have no use in a a basic build/test so we don't want to
    # interfere with options users might set on the env.
    if (NOT "${${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE=${${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE}")
    endif()
    if (NOT "${CTEST_BUILD_NAME}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "CTEST_BUILD_NAME=${CTEST_BUILD_NAME}")
    endif()
    if (NOT "${${PROJECT_NAME}_ENABLE_COVERAGE_TESTING}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "CTEST_DO_COVERAGE_TESTING=${${PROJECT_NAME}_ENABLE_COVERAGE_TESTING}")
    endif()
    if (NOT "${CTEST_BUILD_FLAGS}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "CTEST_BUILD_FLAGS=${CTEST_BUILD_FLAGS}")
    endif()
    if (NOT "${CTEST_PARALLEL_LEVEL}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "CTEST_PARALLEL_LEVEL=${CTEST_PARALLEL_LEVEL}")
    endif()
    if (NOT "${CTEST_DO_SUBMIT}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "CTEST_DO_SUBMIT=${CTEST_DO_SUBMIT}")
    endif()
    if (NOT "${CTEST_DROP_METHOD}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "CTEST_DROP_METHOD=${CTEST_DROP_METHOD}")
    endif()
    if (NOT $"{CTEST_DROP_SITE}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "CTEST_DROP_SITE=${CTEST_DROP_SITE}")
    endif()
    if (NOT "${CTEST_DROP_LOCATION}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "CTEST_DROP_LOCATION=${CTEST_DROP_LOCATION}")
    endif()
    if (NOT "${CTEST_DROP_SITE_COVERAGE}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "CTEST_DROP_SITE_COVERAGE=${CTEST_DROP_SITE_COVERAGE}")
    endif()
    if (NOT "${CTEST_DROP_LOCATION_COVERAGE}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "CTEST_DROP_LOCATION_COVERAGE=${CTEST_DROP_LOCATION_COVERAGE}")
    endif()
    if (NOT "${TRIBITS_2ND_CTEST_DROP_LOCATION}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "TRIBITS_2ND_CTEST_DROP_LOCATION=${TRIBITS_2ND_CTEST_DROP_LOCATION}")
    endif()
    if (NOT "${TRIBITS_2ND_CTEST_DROP_SITE}" STREQUAL "")
      append_set(EXPR_CMND_ARGS "TRIBITS_2ND_CTEST_DROP_SITE=${TRIBITS_2ND_CTEST_DROP_SITE}")
    endif()

    #print_var(${PROJECT_NAME}_EXTRA_REPOSITORIES)
    append_set(EXPR_CMND_ARGS
      ${PROJECT_NAME}_EXTRAREPOS_FILE=${${PROJECT_NAME}_EXTRAREPOS_FILE})
    append_set(EXPR_CMND_ARGS
      ${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=${${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE})
    append_set(EXPR_CMND_ARGS
      ${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES=${${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES})
    join(${PROJECT_NAME}_EXTRA_REPOSITORIES_JOINED "," FALSE
      ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
    append_set(EXPR_CMND_ARGS
      ${PROJECT_NAME}_EXTRA_REPOSITORIES=${${PROJECT_NAME}_EXTRA_REPOSITORIES_JOINED})

    #print_var(EXPR_CMND_ARGS)

    # H.2) Add the custom target to enable all the packages with tests enabled

    if (${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE)

      set(RUNNING_EXP_DASHBOARD_MSG_HEADER
        "Running all-at-once experimental dashboard"
        )

      set(DASHBOARD_TARGET_PRE_CTEST_DRIVER_CMNDS)

      set(DASHBOARD_TARGET_CTEST_DRIVER_CMND_NUM)

      set(DASHBOARD_TARGET_POST_CTEST_DRIVER_CMNDS)

    else()

      set(RUNNING_EXP_DASHBOARD_MSG_HEADER
        "Running package-by-package experimental dashboard"
        )

      set(DASHBOARD_TARGET_PRE_CTEST_DRIVER_CMNDS
        COMMAND echo
        COMMAND echo "***"
        COMMAND echo "*** A) Clean out the list of packages"
        COMMAND echo "***"
        COMMAND echo
        COMMAND echo Running: ${CMAKE_COMMAND} -D${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES:BOOL=TRUE
          -D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=OFF ${PROJECT_SOURCE_DIR}
        COMMAND echo
        COMMAND ${CMAKE_COMMAND} -D${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES:BOOL=TRUE
          -D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=OFF ${PROJECT_SOURCE_DIR}
        # NOTE: Above, if ${PROJECT_NAME}_ENABLE_ALL_PACKAGES was set in CMakeCache.txt, then setting
        # -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=OFF will turn it off in the cache.  Note that it will
        # never be turned on again which means that the list of packages will be set explicitly below.
        )

      set(DASHBOARD_TARGET_CTEST_DRIVER_CMND_NUM "B) ")

      set(DASHBOARD_TARGET_POST_CTEST_DRIVER_CMNDS
        COMMAND echo
        COMMAND echo "***"
        COMMAND echo "*** C) Clean out the list of packages again to clean the cache file"
        COMMAND echo "***"
        COMMAND echo
        COMMAND echo Running: ${CMAKE_COMMAND} -D${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES:BOOL=TRUE
          -D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=OFF ${PROJECT_SOURCE_DIR}
        COMMAND echo
        COMMAND ${CMAKE_COMMAND} -D${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES:BOOL=TRUE
          -D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=OFF ${PROJECT_SOURCE_DIR}
  
        COMMAND echo
        COMMAND echo "***"
        COMMAND echo "*** D) Reconfigure with the original package list"
        COMMAND echo "***"
        COMMAND echo
        COMMAND echo Running: ${CMAKE_COMMAND} ${${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST}
          -D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON ${PROJECT_SOURCE_DIR}
        COMMAND echo
        COMMAND ${CMAKE_COMMAND} ${${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST}
          -D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON ${PROJECT_SOURCE_DIR}
  
        COMMAND echo
        COMMAND echo "See the results at http://${CTEST_DROP_SITE}${CTEST_DROP_LOCATION}&display=project\#Experimental"
        COMMAND echo
        )

    endif()

    add_custom_target( dashboard

      USES_TERMINAL  # Allow real-time STDOUT with ninja target

      VERBATIM  # Recommended

      # WARNING: The echoed command and the actual commands are duplicated!  You have to reproduce them!

      COMMAND echo
      COMMAND echo "**************************************************"
      COMMAND echo "*** ${RUNNING_EXP_DASHBOARD_MSG_HEADER} ***"
      COMMAND echo "**************************************************"
      COMMAND echo
      COMMAND echo ${PROJECT_NAME}_ENABLED_PACKAGES_LIST=${${PROJECT_NAME}_ENABLED_PACKAGES_LIST}
      COMMAND echo

      ${DASHBOARD_TARGET_PRE_CTEST_DRIVER_CMNDS}

      COMMAND echo
      COMMAND echo "***"
      COMMAND echo "*** ${DASHBOARD_TARGET_CTEST_DRIVER_CMND_NUM}Run the dashboard command setting the list of packages"
      COMMAND echo "***"
      COMMAND echo
      COMMAND echo Running: env ${EXPR_CMND_ARGS}
        ${PROJECT_NAME}_PACKAGES=${${PROJECT_NAME}_ENABLED_PACKAGES_LIST}
        PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
        ${CMAKE_CTEST_COMMAND} ${${PROJECT_NAME}_DASHBOARD_CTEST_ARGS} -S
          ${${PROJECT_NAME}_TRIBITS_DIR}/ctest_driver/experimental_build_test.cmake
      COMMAND echo
      COMMAND env ${EXPR_CMND_ARGS}
        ${PROJECT_NAME}_PACKAGES=${${PROJECT_NAME}_ENABLED_PACKAGES_LIST}
        PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
        ${CMAKE_CTEST_COMMAND} ${${PROJECT_NAME}_DASHBOARD_CTEST_ARGS} -S
          ${${PROJECT_NAME}_TRIBITS_DIR}/ctest_driver/experimental_build_test.cmake || echo
      # 2009/07/05: rabartl: Above, I added the ending '|| echo' to always make
      # the command pass so that 'make' will not stop and avoid this last command
      # to set back the enabled packages.

      ${DASHBOARD_TARGET_POST_CTEST_DRIVER_CMNDS}

      )

  endif()

endmacro()
