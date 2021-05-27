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


# Standard TriBITS system includes
INCLUDE(TribitsConstants)
INCLUDE(TribitsProcessExtraRepositoriesList)
INCLUDE(TribitsProcessPackagesAndDirsLists)
INCLUDE(TribitsProcessTplsLists)
INCLUDE(TribitsAdjustPackageEnables)

# Standard TriBITS utilities includes
INCLUDE(TimingUtils)


# @MACRO: TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()
#
# Usage::
#
#   TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()
#
# Macro run at the top project-level scope that reads in packages and TPLs,
# process dependencies, and (optimally) writes XML files of dependency
# information.
#
# On output, this creates all of the package lists and dependency
# data-structures described in `TriBITS System Data Structures and
# Functions`_.
#
MACRO(TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML)

  #
  # A) Read in list of packages and package dependencies
  #

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(SET_UP_DEPENDENCIES_TIME_START_SECONDS)
  ENDIF()

  TRIBITS_READ_DEFINED_EXTERNAL_AND_INTENRAL_TOPLEVEL_PACKAGES_LISTS()

  TRIBITS_READ_PROJECT_AND_PACKAGE_DEPENDENCIES_CREATE_GRAPH_PRINT_DEPS()

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(SET_UP_DEPENDENCIES_TIME_STOP_SECONDS)
    TIMER_PRINT_REL_TIME(${SET_UP_DEPENDENCIES_TIME_START_SECONDS}
      ${SET_UP_DEPENDENCIES_TIME_STOP_SECONDS}
      "\nTotal time to read in and process all package dependencies")
  ENDIF()

  #
  # B) Dump dependnecy info as XML files if asked
  #

  TRIBITS_WRITE_XML_DEPENDENCY_FILES_IF_SUPPORTED()

ENDMACRO()


# @MACRO: TRIBITS_READ_DEFINED_EXTERNAL_AND_INTENRAL_TOPLEVEL_PACKAGES_LISTS()
#
# Usage::
#
#   TRIBITS_READ_DEFINED_EXTERNAL_AND_INTENRAL_TOPLEVEL_PACKAGES_LISTS()
#
# Macro run at the top project-level scope that reads in the contents of all
# of the `<repoDir>/TPLsList.cmake`_ and `<repoDir>/PackagesList.cmake`_ files
# to get the list of defined external packages (TPLs) and internal top-level
# (TriBITS) packages.
#
# On output, this produces the list varaibles::
#
#   ${PROJECT_NAME}_DEFINED_TPLS
#   ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES
#   ${PROJECT_NAME}_ALL_DEFINED_TOPLEVEL_PACKAGES
#   ${PROJECT_NAME}_NUM_ALL_DEFINED_TOPLEVEL_PACKAGES
#
#   ${PROJECT_NAME}_NUM_DEFINED_TPLS
#   ${PROJECT_NAME}_NUM_DEFINED_INTERNAL_PACKAGES
#   ${PROJECT_NAME}_NUM_ALL_DEFINED_TOPLEVEL_PACKAGES
#
#   ${PROJECT_NAME}_PACKAGES (old)
#   ${PROJECT_NAME}_TPLS (old)
#
# and related varaibles.
#
# This includes the files:
#
#  * `<repoDir>/TPLsList.cmake`_ 
#  * `<repoDir>/PackagesList.cmake`_
#
# and calls the macros:
#
#  * `TRIBITS_PROCESS_TPLS_LISTS()`_
#  * `TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS()`_
#
# which set their varaibles.
#
MACRO(TRIBITS_READ_DEFINED_EXTERNAL_AND_INTENRAL_TOPLEVEL_PACKAGES_LISTS)

  TRIBITS_SET_ALL_EXTRA_REPOSITORIES()

  # Set to empty
  SET(${PROJECT_NAME}_PACKAGES)
  SET(${PROJECT_NAME}_TPLS)

  #
  # A) Read list of packages and TPLs from 'PRE' extra repos
  #

  SET(READ_PRE_OR_POST_EXRAREPOS  PRE)
  TRIBITS_READ_EXTRA_REPOSITORIES_LISTS()

  #
  # B) Read list of packages and TPLs from native repos
  #

  FOREACH(NATIVE_REPO ${${PROJECT_NAME}_NATIVE_REPOSITORIES})

    TRIBITS_GET_REPO_NAME_DIR(${NATIVE_REPO}  NATIVE_REPO_NAME  NATIVE_REPO_DIR)
    #PRINT_VAR(NATIVE_REPO_NAME)
    #PRINT_VAR(NATIVE_REPO_DIR)

    # Need to make sure this gets set because logic in Dependencies.cmake files
    # looks for the presents of this variable.
    TRIBITS_SET_BASE_REPO_DIR(${PROJECT_SOURCE_DIR} ${NATIVE_REPO_DIR}
      ${NATIVE_REPO_NAME}_SOURCE_DIR)
    #PRINT_VAR(${NATIVE_REPO_NAME}_SOURCE_DIR)

    IF (NATIVE_REPO STREQUAL ".")
      SET(REPOSITORY_NAME ${PROJECT_NAME})
    ELSE()
      SET(REPOSITORY_NAME ${NATIVE_REPO_NAME})
    ENDIF()

    #
    # B.1) Define the lists of all ${NATIVE_REPO_NAME} native packages and TPLs
    #

    # B.1.a) Read in the list of TPLs for this repo

    IF (${NATIVE_REPO_NAME}_TPLS_FILE_OVERRIDE)
      IF (IS_ABSOLUTE "${${NATIVE_REPO_NAME}_TPLS_FILE_OVERRIDE}")
        SET(${NATIVE_REPO_NAME}_TPLS_FILE
          "${${NATIVE_REPO_NAME}_TPLS_FILE_OVERRIDE}")
      ELSE()
        SET(${NATIVE_REPO_NAME}_TPLS_FILE
          "${${NATIVE_REPO_NAME}_SOURCE_DIR}/${${NATIVE_REPO_NAME}_TPLS_FILE_OVERRIDE}")
      ENDIF()
    ELSE()
      SET(${NATIVE_REPO_NAME}_TPLS_FILE
        "${${NATIVE_REPO_NAME}_SOURCE_DIR}/${${PROJECT_NAME}_TPLS_FILE_NAME}")
    ENDIF()

    MESSAGE("")
    MESSAGE("Reading list of native TPLs from ${${NATIVE_REPO_NAME}_TPLS_FILE}")
    MESSAGE("")

    TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  INCLUDE
      "${${NATIVE_REPO_NAME}_TPLS_FILE}")
    INCLUDE(${${NATIVE_REPO_NAME}_TPLS_FILE})
    TRIBITS_PROCESS_TPLS_LISTS(${NATIVE_REPO_NAME}  ${NATIVE_REPO_DIR})

    # B.1.b) Read in list of packages for this repo

    IF (${NATIVE_REPO_NAME}_PACKAGES_FILE_OVERRIDE)
      IF (IS_ABSOLUTE "${${NATIVE_REPO_NAME}_PACKAGES_FILE_OVERRIDE}")
        SET(${NATIVE_REPO_NAME}_PACKAGES_FILE
          "${${NATIVE_REPO_NAME}_PACKAGES_FILE_OVERRIDE}")
      ELSE()
        SET(${NATIVE_REPO_NAME}_PACKAGES_FILE
          "${${NATIVE_REPO_NAME}_SOURCE_DIR}/${${NATIVE_REPO_NAME}_PACKAGES_FILE_OVERRIDE}")
      ENDIF()
    ELSE()
      SET(${NATIVE_REPO_NAME}_PACKAGES_FILE
        "${${NATIVE_REPO_NAME}_SOURCE_DIR}/${${PROJECT_NAME}_PACKAGES_FILE_NAME}")
    ENDIF()

    MESSAGE("")
    MESSAGE("Reading list of native packages from ${${NATIVE_REPO_NAME}_PACKAGES_FILE}")
    MESSAGE("")

    TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  INCLUDE
      "${${NATIVE_REPO_NAME}_PACKAGES_FILE}")
    INCLUDE(${${NATIVE_REPO_NAME}_PACKAGES_FILE})

    TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${NATIVE_REPO_NAME} ${NATIVE_REPO_DIR})

  ENDFOREACH()

  #
  # C) Read list of packages and TPLs from 'POST' extra repos
  #

  SET(READ_PRE_OR_POST_EXRAREPOS  POST)
  TRIBITS_READ_EXTRA_REPOSITORIES_LISTS()

  #
  # D) Set names of new vars (#63)
  #
  SET(${PROJECT_NAME}_DEFINED_TPLS ${${PROJECT_NAME}_TPLS})
  LIST(LENGTH ${PROJECT_NAME}_DEFINED_TPLS ${PROJECT_NAME}_NUM_DEFINED_TPLS)
  SET(${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES ${${PROJECT_NAME}_PACKAGES})
  LIST(LENGTH ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES
    ${PROJECT_NAME}_NUM_DEFINED_INTERNAL_PACKAGES)
  SET(${PROJECT_NAME}_ALL_DEFINED_TOPLEVEL_PACKAGES
    ${${PROJECT_NAME}_DEFINED_TPLS} ${${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES})
  LIST(LENGTH ${PROJECT_NAME}_ALL_DEFINED_TOPLEVEL_PACKAGES
    ${PROJECT_NAME}_NUM_ALL_DEFINED_TOPLEVEL_PACKAGES)

ENDMACRO()


# Function that writes XML dependnecy files if support for that exists in this
# installation of TriBITS.
#
FUNCTION(TRIBITS_WRITE_XML_DEPENDENCY_FILES_IF_SUPPORTED)
  SET(TRIBITS_PROJECT_CI_SUPPORT_DIR
     "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CI_SUPPORT_DIR}")
  SET(TRIBITS_DUMP_XML_DEPS_MODULE
   "${TRIBITS_PROJECT_CI_SUPPORT_DIR}/TribitsDumpXmlDependenciesFiles.cmake")
  IF (EXISTS "${TRIBITS_DUMP_XML_DEPS_MODULE}")
    INCLUDE(${TRIBITS_DUMP_XML_DEPS_MODULE})
    TRIBITS_WRITE_XML_DEPENDENCY_FILES()
  ENDIF()
ENDFUNCTION()


# Macro that sets ${PROJECT_NAME}_ALL_REPOSITORIES from
# ${PROJECT_NAME}_PRE_REPOSITORIES and ${PROJECT_NAME}_EXTRA_REPOSITORIES if
# it is not alrady set.  Also, it replaces ',' with ';' in the latter.
#
# This function is needed in use cases where extra repos are used where the
# extra repos are not read in through an ExtraRepositoriesList.cmake file and
# instead are directly passed in by the user.
#
MACRO(TRIBITS_SET_ALL_EXTRA_REPOSITORIES)
  IF ("${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES}"   STREQUAL  "")
    # Allow list to be seprated by ',' instead of just by ';'.  This is needed
    # by the unit test driver code
    SPLIT("${${PROJECT_NAME}_PRE_REPOSITORIES}"  ","
      ${PROJECT_NAME}_PRE_REPOSITORIES)
    SPLIT("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  ","
      ${PROJECT_NAME}_EXTRA_REPOSITORIES)
    SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES
      ${${PROJECT_NAME}_PRE_REPOSITORIES}  ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
  ENDIF()
ENDMACRO()


# Macro that processes the list of package and TPLs for the set of 'PRE' or
# 'POST' extra repos.
#
MACRO(TRIBITS_READ_EXTRA_REPOSITORIES_LISTS)

  LIST(LENGTH  ${PROJECT_NAME}_PRE_REPOSITORIES  PRE_EXTRAREPOS_LEN)
  LIST(LENGTH  ${PROJECT_NAME}_EXTRA_REPOSITORIES  POST_EXTRAREPOS_LEN)
  MATH(EXPR  ALL_EXTRAREPOS_LEN  "${PRE_EXTRAREPOS_LEN} + ${POST_EXTRAREPOS_LEN}")

  # See if processing 'PRE' or 'POST' extra repos
  IF (READ_PRE_OR_POST_EXRAREPOS  STREQUAL  "PRE")
    SET(EXTRAREPO_IDX_START  0)
    SET(EXTRAREPO_IDX_END  ${PRE_EXTRAREPOS_LEN})
  ELSEIF (READ_PRE_OR_POST_EXRAREPOS  STREQUAL  "POST")
    SET(EXTRAREPO_IDX_START  ${PRE_EXTRAREPOS_LEN})
    SET(EXTRAREPO_IDX_END  ${ALL_EXTRAREPOS_LEN})
  ELSE()
    MESSAGE(FATAL_ERROR "Invalid value for READ_PRE_OR_POST_EXRAREPOS='${READ_PRE_OR_POST_EXRAREPOS}' ")
  ENDIF()
  # NOTE: For some reason, we can't pass this argument to the function and
  # have it read.  Instead, we have to pass it a local variable.  I will never
  # understand CMake.

  SET(EXTRAREPO_IDX  ${EXTRAREPO_IDX_START})
  WHILE(EXTRAREPO_IDX  LESS  EXTRAREPO_IDX_END)

    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES  ${EXTRAREPO_IDX}  EXTRA_REPO )
    SET(REPOSITORY_NAME  ${EXTRA_REPO})

    #PRINT_VAR(EXTRA_REPO)
    #PRINT_VAR(EXTRAREPO_IDX)
    #PRINT_VAR(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS)

    # Need to make sure this gets set because logic in Dependencies.cmake files
    # looks for the presents of this variable.
    SET(${EXTRA_REPO}_SOURCE_DIR "${PROJECT_SOURCE_DIR}/${EXTRA_REPO}")
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(${EXTRA_REPO}_SOURCE_DIR)
    ENDIF()
    # ToDo: TriBITS:73: Get ${EXTRA_REPO}_SOURCE_DIR from
    # ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIR when it exists.

    SET(EXTRAREPO_PACKSTAT "")
    IF (${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS)
      LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS ${EXTRAREPO_IDX}
        EXTRAREPO_PACKSTAT )
    ENDIF()

    IF (EXTRAREPO_PACKSTAT STREQUAL NOPACKAGES)

      MESSAGE("")
      MESSAGE("Skipping reading packages and TPLs for ${READ_PRE_OR_POST_EXRAREPOS} extra repo ${EXTRA_REPO} because marked NOPACKAGES ... ")
      MESSAGE("")
      # ToDo: TriBITS:73: Don't print the above message by default.  It is
      # just clutter.

    ELSE()

      # Read in the add-on TPLs from the extra repo

      SET(${EXTRA_REPO}_TPLS_FILE
        "${${EXTRA_REPO}_SOURCE_DIR}/${${PROJECT_NAME}_EXTRA_TPLS_FILE_NAME}")
      PRINT_VAR(${EXTRA_REPO}_SOURCE_DIR)
      PRINT_VAR(${PROJECT_NAME}_EXTRA_TPLS_FILE_NAME)
      PRINT_VAR(${EXTRA_REPO}_TPLS_FILE)

      MESSAGE("")
      MESSAGE("Reading list of ${READ_PRE_OR_POST_EXRAREPOS} extra TPLs from ${${EXTRA_REPO}_TPLS_FILE} ... ")
      MESSAGE("")

      IF (NOT EXISTS "${${EXTRA_REPO}_TPLS_FILE}")
        IF (${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES)
          MESSAGE("-- "
            "NOTE: Ignoring missing ${READ_PRE_OR_POST_EXRAREPOS} extra repo '${EXTRA_REPO}' TPLs list file '${${EXTRA_REPO}_TPLS_FILE}' on request!" )
        ELSE()
          MESSAGE( SEND_ERROR
            "ERROR: Skipping missing ${READ_PRE_OR_POST_EXRAREPOS} extra repo '${EXTRA_REPO}' TPLs list file '${${EXTRA_REPO}_TPLS_FILE}'!")
        ENDIF()
      ELSE()
        TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  INCLUDE  "${${EXTRA_REPO}_TPLS_FILE}")
        INCLUDE("${${EXTRA_REPO}_TPLS_FILE}")
        SET(APPEND_TO_TPLS_LIST  TRUE)
        TRIBITS_PROCESS_TPLS_LISTS(${EXTRA_REPO}  ${EXTRA_REPO})
      ENDIF()

      # Read in the add-on packages from the extra repo

      #PRINT_VAR(${EXTRA_REPO}_PACKAGES_LIST_FILE)
      IF (${EXTRA_REPO}_PACKAGES_LIST_FILE)
        SET(EXTRAREPO_PACKAGES_FILE
          "${PROJECT_SOURCE_DIR}/${${EXTRA_REPO}_PACKAGES_LIST_FILE}")
      ELSE()
        SET(EXTRAREPO_PACKAGES_FILE
          "${${EXTRA_REPO}_SOURCE_DIR}/${${PROJECT_NAME}_EXTRA_PACKAGES_FILE_NAME}")
      ENDIF()

      MESSAGE("")
      MESSAGE("Reading list of ${READ_PRE_OR_POST_EXRAREPOS} extra packages from ${EXTRAREPO_PACKAGES_FILE} ... ")
      MESSAGE("")

      IF (NOT EXISTS "${EXTRAREPO_PACKAGES_FILE}")
        IF (${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES)
          MESSAGE("-- "
            "NOTE: Ignoring missing ${READ_PRE_OR_POST_EXRAREPOS} extra repo '${EXTRA_REPO}' packages list file '${EXTRAREPO_PACKAGES_FILE}' on request!")
        ELSE()
          MESSAGE( SEND_ERROR
            "ERROR: Skipping missing ${READ_PRE_OR_POST_EXRAREPOS} extra repo '${EXTRA_REPO}' packages list file '${EXTRAREPO_PACKAGES_FILE}'!")
          # ToDo: TriBITS:73: Change to FATAL_ERROR to abort early
        ENDIF()
      ELSE()
        TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  INCLUDE  "${EXTRAREPO_PACKAGES_FILE}")
        INCLUDE("${EXTRAREPO_PACKAGES_FILE}")
        SET(APPEND_TO_PACKAGES_LIST  TRUE)
        TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${EXTRA_REPO} ${EXTRA_REPO})
      ENDIF()

    ENDIF()

    MATH(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")

  ENDWHILE()

ENDMACRO()
