# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# Standard TriBITS system includes
include("${CMAKE_CURRENT_LIST_DIR}/../common/TribitsConstants.cmake")
include(TribitsProcessExtraRepositoriesList)
include(TribitsProcessPackagesAndDirsLists)
include(TribitsProcessTplsLists)
include(TribitsReadDepsFilesCreateDepsGraph)
include(TribitsConfigureTiming)

# Standard TriBITS utilities includes
include(TimingUtils)


# @MACRO: tribits_read_all_project_deps_files_create_deps_graph()
#
# Usage::
#
#   tribits_read_all_project_deps_files_create_deps_graph()
#
# Macro run at the top project-level scope that reads the lists of packages
# and TPLs and creates the packages dependency graph.
#
# On output, this creates all of the package lists and dependency
# data-structures described in the section `TriBITS System Data Structures`_
# and more specifically the sections:
#
# * `Lists of external and internal packages`_
# * `Variables defining the package dependencies graph`_
# * `TriBITS Package Top-Level Local Variables`_
# * `TriBITS Subpackage Top-Level Local Variables`_
# * `TriBITS Package Cache Variables`_
#
# See `Function call tree for constructing package dependency graph`_.
#
macro(tribits_read_all_project_deps_files_create_deps_graph)

  tribits_config_code_start_timer(SET_UP_DEPENDENCIES_TIME_START_SECONDS)

  tribits_read_defined_external_and_internal_toplevel_packages_lists()

  tribits_read_deps_files_create_deps_graph()

  # ${PROJECT_NAME}_DEFINED_PACKAGES
  set(${PROJECT_NAME}_DEFINED_PACKAGES
    ${${PROJECT_NAME}_DEFINED_TPLS}
    ${${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES})
  list(LENGTH ${PROJECT_NAME}_DEFINED_PACKAGES
    ${PROJECT_NAME}_NUM_DEFINED_PACKAGES)

  tribits_config_code_stop_timer(SET_UP_DEPENDENCIES_TIME_START_SECONDS
    "\nTotal time to read in all dependencies files and build dependencies graph")

endmacro()


# @MACRO: tribits_read_defined_external_and_internal_toplevel_packages_lists()
#
# Usage::
#
#   tribits_read_defined_external_and_internal_toplevel_packages_lists()
#
# Macro run at the top project-level scope that reads in the contents of all
# of the `<repoDir>/TPLsList.cmake`_ and `<repoDir>/PackagesList.cmake`_ files
# to get the list of defined external packages (TPLs) and internal top-level
# (TriBITS) packages.
#
# On output, this produces the local variables:
#
#   * `${PROJECT_NAME}_DEFINED_TPLS`_
#   * `${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES`_
#   * `${PROJECT_NAME}_DEFINED_TOPLEVEL_PACKAGES`_
#
# and the length vars for these:
#
#   * `${PROJECT_NAME}_NUM_DEFINED_TPLS`_
#   * `${PROJECT_NAME}_NUM_DEFINED_INTERNAL_TOPLEVEL_PACKAGES`_
#   * `${PROJECT_NAME}_NUM_DEFINED_TOPLEVEL_PACKAGES`_
#
# This includes the files:
#
#  * `<repoDir>/TPLsList.cmake`_ 
#  * `<repoDir>/PackagesList.cmake`_
#
# and calls the macros:
#
#  * `tribits_process_tpls_lists()`_
#  * `tribits_process_packages_and_dirs_lists()`_
#
# which set their variables.
#
# See `Function call tree for constructing package dependency graph`_
#
macro(tribits_read_defined_external_and_internal_toplevel_packages_lists)

  tribits_set_all_extra_repositories()

  # Set package list vars to empty
  set(${PROJECT_NAME}_DEFINED_TPLS "")
  set(${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES "")

  #
  # A) Read list of external packages/TPLs and top-level internal packages
  # from 'PRE' extra repos
  #

  set(READ_PRE_OR_POST_EXRAREPOS  PRE)
  tribits_read_extra_repositories_lists()

  #
  # B) Read list of external packages/TPLs and top-level internal packages
  # from the native repos
  #

  foreach(NATIVE_REPO ${${PROJECT_NAME}_NATIVE_REPOSITORIES})

    tribits_get_repo_name_dir(${NATIVE_REPO}  NATIVE_REPO_NAME  NATIVE_REPO_DIR)
    #print_var(NATIVE_REPO_NAME)
    #print_var(NATIVE_REPO_DIR)

    # Need to make sure this gets set because logic in Dependencies.cmake files
    # looks for the presents of this variable.
    tribits_set_base_repo_dir(${PROJECT_SOURCE_DIR} ${NATIVE_REPO_DIR}
      ${NATIVE_REPO_NAME}_SOURCE_DIR)
    #print_var(${NATIVE_REPO_NAME}_SOURCE_DIR)

    if (NATIVE_REPO STREQUAL ".")
      set(REPOSITORY_NAME ${PROJECT_NAME})
    else()
      set(REPOSITORY_NAME ${NATIVE_REPO_NAME})
    endif()

    #
    # B.1) Define the lists of all ${NATIVE_REPO_NAME} native packages and TPLs
    #

    # B.1.a) Read in the list of TPLs for this repo

    if (${NATIVE_REPO_NAME}_TPLS_FILE_OVERRIDE)
      if (IS_ABSOLUTE "${${NATIVE_REPO_NAME}_TPLS_FILE_OVERRIDE}")
        set(${NATIVE_REPO_NAME}_TPLS_FILE
          "${${NATIVE_REPO_NAME}_TPLS_FILE_OVERRIDE}")
      else()
        set(${NATIVE_REPO_NAME}_TPLS_FILE
          "${${NATIVE_REPO_NAME}_SOURCE_DIR}/${${NATIVE_REPO_NAME}_TPLS_FILE_OVERRIDE}")
      endif()
    else()
      set(${NATIVE_REPO_NAME}_TPLS_FILE
        "${${NATIVE_REPO_NAME}_SOURCE_DIR}/${${PROJECT_NAME}_TPLS_FILE_NAME}")
    endif()

    message("")
    message("Reading list of native TPLs from ${${NATIVE_REPO_NAME}_TPLS_FILE}")
    message("")

    tribits_trace_file_processing(REPOSITORY  INCLUDE
      "${${NATIVE_REPO_NAME}_TPLS_FILE}")
    include(${${NATIVE_REPO_NAME}_TPLS_FILE})
    tribits_process_tpls_lists(${NATIVE_REPO_NAME}  ${NATIVE_REPO_DIR})

    # B.1.b) Read in list of packages for this repo

    if (${NATIVE_REPO_NAME}_PACKAGES_FILE_OVERRIDE)
      if (IS_ABSOLUTE "${${NATIVE_REPO_NAME}_PACKAGES_FILE_OVERRIDE}")
        set(${NATIVE_REPO_NAME}_PACKAGES_FILE
          "${${NATIVE_REPO_NAME}_PACKAGES_FILE_OVERRIDE}")
      else()
        set(${NATIVE_REPO_NAME}_PACKAGES_FILE
          "${${NATIVE_REPO_NAME}_SOURCE_DIR}/${${NATIVE_REPO_NAME}_PACKAGES_FILE_OVERRIDE}")
      endif()
    else()
      set(${NATIVE_REPO_NAME}_PACKAGES_FILE
        "${${NATIVE_REPO_NAME}_SOURCE_DIR}/${${PROJECT_NAME}_PACKAGES_FILE_NAME}")
    endif()

    message("")
    message("Reading list of native packages from ${${NATIVE_REPO_NAME}_PACKAGES_FILE}")
    message("")

    tribits_trace_file_processing(REPOSITORY  INCLUDE
      "${${NATIVE_REPO_NAME}_PACKAGES_FILE}")
    include(${${NATIVE_REPO_NAME}_PACKAGES_FILE})

    tribits_process_packages_and_dirs_lists(${NATIVE_REPO_NAME} ${NATIVE_REPO_DIR})

  endforeach()

  #
  # C) Read list of external packages/TPLs and top-level internal packages
  # from 'POST' extra repos
  #

  set(READ_PRE_OR_POST_EXRAREPOS  POST)
  tribits_read_extra_repositories_lists()

  #
  # D) Compute lengths and other combined quantities
  #

  # ${PROJECT_NAME}_NUM_DEFINED_TPLS
  list(LENGTH ${PROJECT_NAME}_DEFINED_TPLS ${PROJECT_NAME}_NUM_DEFINED_TPLS)

  # ${PROJECT_NAME}_DEFINED_TOPLEVEL_PACKAGES
  set(${PROJECT_NAME}_DEFINED_TOPLEVEL_PACKAGES
    ${${PROJECT_NAME}_DEFINED_TPLS}
    ${${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES})
  list(LENGTH ${PROJECT_NAME}_DEFINED_TOPLEVEL_PACKAGES
    ${PROJECT_NAME}_NUM_DEFINED_TOPLEVEL_PACKAGES)

endmacro()


# @FUNCTION: tribits_write_xml_dependency_files_if_supported()
#
# Usage::
#
#   tribits_write_xml_dependency_files_if_supported()
#
# Function that writes XML dependency files if support for that exists in this
# instance of TriBITs.
#
# See `Function call tree for constructing package dependency graph`_
#
function(tribits_write_xml_dependency_files_if_supported)
  set(TRIBITS_PROJECT_CI_SUPPORT_DIR
     "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CI_SUPPORT_DIR}")
  set(TRIBITS_DUMP_XML_DEPS_MODULE
   "${TRIBITS_PROJECT_CI_SUPPORT_DIR}/TribitsWriteXmlDependenciesFiles.cmake")
  if (EXISTS "${TRIBITS_DUMP_XML_DEPS_MODULE}")
    include(${TRIBITS_DUMP_XML_DEPS_MODULE})
    tribits_write_xml_dependency_files()
  endif()
endfunction()


# Macro that sets ${PROJECT_NAME}_ALL_REPOSITORIES from
# ${PROJECT_NAME}_PRE_REPOSITORIES and ${PROJECT_NAME}_EXTRA_REPOSITORIES if
# it is not already set.  Also, it replaces ',' with ';' in the latter.
#
# This function is needed in use cases where extra repos are used where the
# extra repos are not read in through an ExtraRepositoriesList.cmake file and
# instead are directly passed in by the user.
#
macro(tribits_set_all_extra_repositories)
  if ("${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES}"   STREQUAL  "")
    # Allow list to be separated by ',' instead of just by ';'.  This is needed
    # by the unit test driver code
    split("${${PROJECT_NAME}_PRE_REPOSITORIES}"  ","
      ${PROJECT_NAME}_PRE_REPOSITORIES)
    split("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  ","
      ${PROJECT_NAME}_EXTRA_REPOSITORIES)
    set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES
      ${${PROJECT_NAME}_PRE_REPOSITORIES}  ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
  endif()
endmacro()


# Macro that processes the list of package and TPLs for the set of 'PRE' or
# 'POST' extra repos.
#
macro(tribits_read_extra_repositories_lists)

  list(LENGTH  ${PROJECT_NAME}_PRE_REPOSITORIES  PRE_EXTRAREPOS_LEN)
  list(LENGTH  ${PROJECT_NAME}_EXTRA_REPOSITORIES  POST_EXTRAREPOS_LEN)
  math(EXPR  ALL_EXTRAREPOS_LEN  "${PRE_EXTRAREPOS_LEN} + ${POST_EXTRAREPOS_LEN}")

  # See if processing 'PRE' or 'POST' extra repos
  if (READ_PRE_OR_POST_EXRAREPOS  STREQUAL  "PRE")
    set(EXTRAREPO_IDX_START  0)
    set(EXTRAREPO_IDX_END  ${PRE_EXTRAREPOS_LEN})
  elseif (READ_PRE_OR_POST_EXRAREPOS  STREQUAL  "POST")
    set(EXTRAREPO_IDX_START  ${PRE_EXTRAREPOS_LEN})
    set(EXTRAREPO_IDX_END  ${ALL_EXTRAREPOS_LEN})
  else()
    message(FATAL_ERROR "Invalid value for READ_PRE_OR_POST_EXRAREPOS='${READ_PRE_OR_POST_EXRAREPOS}' ")
  endif()
  # NOTE: For some reason, we can't pass this argument to the function and
  # have it read.  Instead, we have to pass it a local variable.  I will never
  # understand CMake.

  set(EXTRAREPO_IDX  ${EXTRAREPO_IDX_START})
  while(EXTRAREPO_IDX  LESS  EXTRAREPO_IDX_END)

    list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES  ${EXTRAREPO_IDX}  EXTRA_REPO )
    set(REPOSITORY_NAME  ${EXTRA_REPO})

    #print_var(EXTRA_REPO)
    #print_var(EXTRAREPO_IDX)
    #print_var(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS)

    # Need to make sure this gets set because logic in Dependencies.cmake files
    # looks for the presents of this variable.
    set(${EXTRA_REPO}_SOURCE_DIR "${PROJECT_SOURCE_DIR}/${EXTRA_REPO}")
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      print_var(${EXTRA_REPO}_SOURCE_DIR)
    endif()
    # ToDo: TriBITS:73: Get ${EXTRA_REPO}_SOURCE_DIR from
    # ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIR when it exists.

    set(EXTRAREPO_PACKSTAT "")
    if (${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS)
      list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS ${EXTRAREPO_IDX}
        EXTRAREPO_PACKSTAT )
    endif()

    if (EXTRAREPO_PACKSTAT STREQUAL NOPACKAGES)

      message("")
      message("Skipping reading packages and TPLs for ${READ_PRE_OR_POST_EXRAREPOS} extra repo ${EXTRA_REPO} because marked NOPACKAGES ... ")
      message("")
      # ToDo: TriBITS:73: Don't print the above message by default.  It is
      # just clutter.

    else()

      # Read in the add-on TPLs from the extra repo

      set(${EXTRA_REPO}_TPLS_FILE
        "${${EXTRA_REPO}_SOURCE_DIR}/${${PROJECT_NAME}_EXTRA_TPLS_FILE_NAME}")
      print_var(${EXTRA_REPO}_SOURCE_DIR)
      print_var(${PROJECT_NAME}_EXTRA_TPLS_FILE_NAME)
      print_var(${EXTRA_REPO}_TPLS_FILE)

      message("")
      message("Reading list of ${READ_PRE_OR_POST_EXRAREPOS} extra TPLs from ${${EXTRA_REPO}_TPLS_FILE} ... ")
      message("")

      if (NOT EXISTS "${${EXTRA_REPO}_TPLS_FILE}")
        if (${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES)
          message("-- "
            "NOTE: Ignoring missing ${READ_PRE_OR_POST_EXRAREPOS} extra repo '${EXTRA_REPO}' TPLs list file '${${EXTRA_REPO}_TPLS_FILE}' on request!" )
        else()
          message( SEND_ERROR
            "ERROR: Skipping missing ${READ_PRE_OR_POST_EXRAREPOS} extra repo '${EXTRA_REPO}' TPLs list file '${${EXTRA_REPO}_TPLS_FILE}'!")
        endif()
      else()
        tribits_trace_file_processing(REPOSITORY  INCLUDE  "${${EXTRA_REPO}_TPLS_FILE}")
        include("${${EXTRA_REPO}_TPLS_FILE}")
        set(APPEND_TO_TPLS_LIST  TRUE)
        tribits_process_tpls_lists(${EXTRA_REPO}  ${EXTRA_REPO})
      endif()

      # Read in the add-on packages from the extra repo

      #print_var(${EXTRA_REPO}_PACKAGES_LIST_FILE)
      if (${EXTRA_REPO}_PACKAGES_LIST_FILE)
        set(EXTRAREPO_PACKAGES_FILE
          "${PROJECT_SOURCE_DIR}/${${EXTRA_REPO}_PACKAGES_LIST_FILE}")
      else()
        set(EXTRAREPO_PACKAGES_FILE
          "${${EXTRA_REPO}_SOURCE_DIR}/${${PROJECT_NAME}_EXTRA_PACKAGES_FILE_NAME}")
      endif()

      message("")
      message("Reading list of ${READ_PRE_OR_POST_EXRAREPOS} extra packages from ${EXTRAREPO_PACKAGES_FILE} ... ")
      message("")

      if (NOT EXISTS "${EXTRAREPO_PACKAGES_FILE}")
        if (${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES)
          message("-- "
            "NOTE: Ignoring missing ${READ_PRE_OR_POST_EXRAREPOS} extra repo '${EXTRA_REPO}' packages list file '${EXTRAREPO_PACKAGES_FILE}' on request!")
        else()
          message( SEND_ERROR
            "ERROR: Skipping missing ${READ_PRE_OR_POST_EXRAREPOS} extra repo '${EXTRA_REPO}' packages list file '${EXTRAREPO_PACKAGES_FILE}'!")
          # ToDo: TriBITS:73: Change to FATAL_ERROR to abort early
        endif()
      else()
        tribits_trace_file_processing(REPOSITORY  INCLUDE  "${EXTRAREPO_PACKAGES_FILE}")
        include("${EXTRAREPO_PACKAGES_FILE}")
        set(APPEND_TO_PACKAGES_LIST  TRUE)
        tribits_process_packages_and_dirs_lists(${EXTRA_REPO} ${EXTRA_REPO})
      endif()

    endif()

    math(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")

  endwhile()

endmacro()
