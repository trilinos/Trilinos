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
# This file gets included in the main TriBITS framework.  It is put here to
# reduce the size of the tribits/core/ directory.
#


#
# Macro that drives a experimental 'dashboard' target
#

MACRO(TRIBITS_ADD_DASHBOARD_TARGET)

  IF (NOT (WIN32 AND NOT CYGWIN))

    ADVANCED_SET(${PROJECT_NAME}_DASHBOARD_CTEST_ARGS "-V" CACHE STRING
      "Extra arguments to pass to CTest when calling 'ctest -S' to run the 'dashboard' make target." )

    ADVANCED_SET(CTEST_BUILD_FLAGS "" CACHE STRING
      "Sets CTEST_BUILD_FLAGS on the env before invoking 'ctest -S'." )

    ADVANCED_SET(CTEST_PARALLEL_LEVEL "" CACHE STRING
      "Sets CTEST_PARALLEL_LEVEL on the env before invoking 'ctest -S'." )

    # H.1) Enable all packages that are enabled and have tests enabled

    SET(${PROJECT_NAME}_ENABLED_PACKAGES_LIST)
    SET(${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST)
    FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_SE_PACKAGES})
      IF (${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} AND ${TRIBITS_PACKAGE}_ENABLE_TESTS)
        IF (${PROJECT_NAME}_ENABLED_PACKAGES_LIST)
          SET(${PROJECT_NAME}_ENABLED_PACKAGES_LIST
            "${${PROJECT_NAME}_ENABLED_PACKAGES_LIST}\;${TRIBITS_PACKAGE}")
        ELSE()
          SET(${PROJECT_NAME}_ENABLED_PACKAGES_LIST "${TRIBITS_PACKAGE}")
        ENDIF()
        SET(${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST
          ${${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST} -D${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE}=ON)
      ENDIF()
    ENDFOREACH()
    #PRINT_VAR(${PROJECT_NAME}_ENABLED_PACKAGES_LIST)

    SET(EXPR_CMND_ARGS)

    # Hard override options used by basic build and tests
    APPEND_SET(EXPR_CMND_ARGS "TRIBITS_PROJECT_ROOT=${${PROJECT_NAME}_SOURCE_DIR}")
    APPEND_SET(EXPR_CMND_ARGS "${PROJECT_NAME}_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}")
    APPEND_SET(EXPR_CMND_ARGS "${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS='${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS}'")
    APPEND_SET(EXPR_CMND_ARGS "${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE='${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE}'")

    # Conditionally override options used only for testing.  These options
    # have no use in a a basic build/test so we don't want to interfere with
    # options users might set on the env.
    IF (${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE)
      APPEND_SET(EXPR_CMND_ARGS "${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE=TRUE")
    ENDIF()
    IF (${PROJECT_NAME}_CTEST_USE_NEW_AAO_FEATURES)
      APPEND_SET(EXPR_CMND_ARGS "${PROJECT_NAME}_CTEST_USE_NEW_AAO_FEATURES=TRUE")
    ENDIF()
    IF (${PROJECT_NAME}_ENABLE_COVERAGE_TESTING)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_DO_COVERAGE_TESTING=TRUE")
    ENDIF()
    IF (CTEST_BUILD_FLAGS)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_BUILD_FLAGS='${CTEST_BUILD_FLAGS}'")
    ENDIF()
    IF (CTEST_PARALLEL_LEVEL)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_PARALLEL_LEVEL=${CTEST_PARALLEL_LEVEL}")
    ENDIF()
    IF (NOT "${CTEST_DO_SUBMIT}" STREQUAL "")
      APPEND_SET(EXPR_CMND_ARGS "CTEST_DO_SUBMIT=${CTEST_DO_SUBMIT}")
    ENDIF()
    IF (CTEST_DROP_SITE)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_DROP_SITE=${CTEST_DROP_SITE}")
    ENDIF()
    IF (CTEST_DROP_LOCATION)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_DROP_LOCATION=${CTEST_DROP_LOCATION}")
    ENDIF()
    IF (CTEST_DROP_SITE_COVERAGE)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_DROP_SITE_COVERAGE=${CTEST_DROP_SITE_COVERAGE}")
    ENDIF()
    IF (CTEST_DROP_LOCATION_COVERAGE)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_DROP_LOCATION_COVERAGE=${CTEST_DROP_LOCATION_COVERAGE}")
    ENDIF()

    #PRINT_VAR(${PROJECT_NAME}_EXTRA_REPOSITORIES)
    APPEND_SET(EXPR_CMND_ARGS
      ${PROJECT_NAME}_EXTRAREPOS_FILE=${${PROJECT_NAME}_EXTRAREPOS_FILE})
    APPEND_SET(EXPR_CMND_ARGS
      ${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=${${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE})
    APPEND_SET(EXPR_CMND_ARGS
      ${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES=${${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES})
    JOIN(${PROJECT_NAME}_EXTRA_REPOSITORIES_JOINED "," FALSE
      ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
    APPEND_SET(EXPR_CMND_ARGS
      ${PROJECT_NAME}_EXTRA_REPOSITORIES=${${PROJECT_NAME}_EXTRA_REPOSITORIES_JOINED})

    #PRINT_VAR(EXPR_CMND_ARGS)

    # H.2) Add the custom target to enable all the packages with tests enabled

    IF (${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE)

      SET(RUNNING_EXP_DASHBOARD_MSG_HEADER
        "Running all-at-once experimental dashboard"
        )

      SET(DASHBOARD_TARGET_PRE_CTEST_DRIVER_CMNDS)

      SET(DASHBOARD_TARGET_CTEST_DRIVER_CMND_NUM)

      SET(DASHBOARD_TARGET_POST_CTEST_DRIVER_CMNDS)

    ELSE()

      SET(RUNNING_EXP_DASHBOARD_MSG_HEADER
        "Running package-by-package experimental dashboard"
        )

      SET(DASHBOARD_TARGET_PRE_CTEST_DRIVER_CMNDS
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

      SET(DASHBOARD_TARGET_CTEST_DRIVER_CMND_NUM "B) ")

      SET(DASHBOARD_TARGET_POST_CTEST_DRIVER_CMNDS
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

    ENDIF()

    ADD_CUSTOM_TARGET( dashboard

      VERBATIM

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

      COMMAND echo
      COMMAND echo "See the results at http://${CTEST_DROP_SITE}${CTEST_DROP_LOCATION}&display=project\#Experimental"
      COMMAND echo

      )

  ENDIF()

ENDMACRO()
