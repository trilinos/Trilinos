# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include("${CMAKE_CURRENT_LIST_DIR}/../common/TribitsConstants.cmake")
include(TribitsListHelpers)

include(PrintVar)
include(Split)


# @MACRO: tribits_repository_define_tpls()
#
# Define the list of `TriBITS External Packages/TPLs`_ for a given `TriBITS
# Repository`_ which includes the external package/TPL name, TriBITS TPL find
# module, and classification .  This macro is typically called from inside of
# a TriBITS Repository's `<repoDir>/TPLsList.cmake`_ file.
#
# Usage::
#
#   tribits_repository_define_tpls(
#     <tpl0_name>   <tpl0_findmod>  <tpl0_classif>
#     <tpl1_name>   <tpl1_findmod>  <tpl1_classif>
#     ...
#     )
#
# This macro sets up a 2D array of ``NumTPLS`` by ``NumColumns`` listing out
# the `TriBITS TPLs`_ for a `TriBITS Repository`_.  Each row (with 3 entries)
# specifies a different TriBITS exernal package/TPL which contains the columns
# (ordered 0-2):
#
# 0. **TPL** (``<tpli_name>``): The name of the TriBITS external package/TPL
#    ``<tplName>``.  This name must be unique across all other TriBITS TPLs in
#    this or any other TriBITS repo that might be combined into a single
#    TriBITS project meta-build (see `Globally unique TriBITS TPL names`_).
#    However, a TPL can be redefined from an upstream repo (see below).  The
#    name should be a valid identifier (e.g. matches the regex
#    ``[a-zA-Z_][a-zA-Z0-9_]*``).  TPL names typically use mixed case
#    (e.g. ``SomeTpl``, not ``SOMETPL``).
#
# 1. **FINDMOD** (``<tpli_findmod>``): For a TriBITS external package/TPL that
#    **is not** a `TriBITS-compliant external package`_, this is set to the
#    relative or absolute path for the TriBITS TPL find module, usually with
#    the name `FindTPL<tplName>.cmake`_.  If it is a relative path, it is
#    considered relative to the repository base directory ``<repoDir>``.  If
#    just the base path for the find module is given, ending with ``"/"``
#    (e.g. ``"cmake/tpls/"``), then the find module will be assumed to be
#    under that this directory with the standard name
#    ``FindTPL<tplName>.cmake``.  (See `Creating the FindTPL<tplName>.cmake
#    file`_.)  However, if the external package **is** a `TriBITS-compliant
#    external package`_, provide the value ``TRIBITS_PKG`` instead and no
#    ``FindTPL<tplName>.cmake`` file is needed.  If a
#    ``FindTPL<tplName>Dependencies.cmake`` file is needed in this case, then
#    provide the path to that file (relative or absolute, directory or file
#    path) using ``TRIBITS_PKG:<path-to-dir-or-file>``.  This field is used to
#    set the variables `<tplName>_FINDMOD`_ and
#    `<tplName>_DEPENDENCIES_FILE`_.
#
# 2. **CLASSIFICATION** (``<pkgi_classif>``): Gives the `Package Test Group`_
#    `PT`_, `ST`_, or `EX`_ and the maturity level ``EP``, ``RS``, ``PG``,
#    ``PM``, ``GRS``, ``GPG``, ``GPM``, ``UM``.  These are separated by a coma
#    with no space in between such as ``"RS,PT"`` for a "Research Stable",
#    "Primary Tested" package.  No spaces are allowed so that CMake treats
#    this a one field in the array.  The maturity level can be left off in
#    which case it is assumed to be ``UM`` for "Unspecified Maturity".  This
#    field is used to set the variable `<tplName>_TESTGROUP`_ and
#    ``<tplName>_MATURITY_LEVEL``.
#
# A TPL defined in a upstream repo can be listed again in a downstream repo,
# which allows redefining the find module that is used to specify the external
# package/TPL.  This allows downstream repos to add additional requirements
# for a given TPL (i.e. add more libraries, headers, etc.).  However, the
# downstream repo's find module file must find the TPL components that are
# fully compatible with the upstream defined find module in terms of what it
# provides for packages in the upstream repos.
#
macro(tribits_repository_define_tpls)
  assert_defined(REPOSITORY_NAME)
  set(${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS "${ARGN}")
endmacro()


# @MACRO: tribits_process_tpls_lists()
#
# This macro that processes the project-level variable::
#
#   ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS
#
# This updates the project-level variables:
#
#   * `${PROJECT_NAME}_DEFINED_TPLS`_
#   * `${PROJECT_NAME}_NUM_DEFINED_TPLS`_
#
# For each TPL, it also sets the variables:
#
#   * `${TPL_NAME}_FINDMOD`_
#   * `${TPL_NAME}_TESTGROUP`_
#   * `${TPL_NAME}_DEPENDENCIES_FILE`_
#   * `${TPL_NAME}_TPLS_LIST_FILE`_
#   * `${TPL_NAME}_PACKAGE_BUILD_STATUS`_ (to ``EXTERNAL``)
#
# See `Function call tree for constructing package dependency graph`_
#
macro(tribits_process_tpls_lists  REPOSITORY_NAME  REPOSITORY_DIR)

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("TRIBITS_PROCESS_TPLS_LISTS:  '${REPOSITORY_NAME}'  '${REPOSITORY_DIR}'")
  endif()

  #set(TRIBITS_PROCESS_TPLS_LISTS_DEBUG ON)
  set(TRIBITS_PROCESS_TPLS_LISTS_DEBUG OFF)

  set(TPL_NAME_OFFSET 0)
  set(TPL_FINDMOD_OFFSET 1)
  set(TPL_CLASSIFICATION_OFFSET 2)
  set(TPL_NUM_COLUMNS 3)

  list(LENGTH ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS
    ${REPOSITORY_NAME}_CURR_NUM_TPLS_FULL)
  math(EXPR ${REPOSITORY_NAME}_CURR_NUM_TPLS
    "${${REPOSITORY_NAME}_CURR_NUM_TPLS_FULL}/${TPL_NUM_COLUMNS}")

  if (${REPOSITORY_NAME}_CURR_NUM_TPLS GREATER 0)

    math(EXPR ${REPOSITORY_NAME}_LAST_TPL_IDX
      "${${REPOSITORY_NAME}_CURR_NUM_TPLS}-1")

    foreach(TPL_IDX RANGE ${${REPOSITORY_NAME}_LAST_TPL_IDX})

      if (TRIBITS_PROCESS_TPLS_LISTS_DEBUG)
        print_var(TPL_IDX)
      endif()

      # Get fields for this TPL

      math(EXPR TPL_NAME_IDX
        "${TPL_IDX}*${TPL_NUM_COLUMNS}+${TPL_NAME_OFFSET}")
      list(GET ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS ${TPL_NAME_IDX}
        TPL_NAME)
      if (TRIBITS_PROCESS_TPLS_LISTS_DEBUG)
        print_var(TPL_NAME)
      endif()

      math(EXPR TPL_FINDMOD_IDX
        "${TPL_IDX}*${TPL_NUM_COLUMNS}+${TPL_FINDMOD_OFFSET}")
      list(GET ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS ${TPL_FINDMOD_IDX}
        TPL_FINDMOD)
      if (TRIBITS_PROCESS_TPLS_LISTS_DEBUG)
        print_var(TPL_FINDMOD)
      endif()

      math(EXPR TPL_CLASSIFICATION_IDX
        "${TPL_IDX}*${TPL_NUM_COLUMNS}+${TPL_CLASSIFICATION_OFFSET}")
      list(GET ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS ${TPL_CLASSIFICATION_IDX}
        TPL_CLASSIFICATION)
      if (TRIBITS_PROCESS_TPLS_LISTS_DEBUG)
        print_var(TPL_CLASSIFICATION)
      endif()

      # ToDo: Parse out TESTGROUP and MATURITYLEVEL (Trilinos #6042)
      set(TPL_TESTGROUP ${TPL_CLASSIFICATION})

      tribits_update_ps_pt_ss_st(TPL  ${TPL_NAME}  TPL_TESTGROUP)

      # Update TPLS list (unless the TPL already exists)

      if (${TPL_NAME}_FINDMOD)
        # If the variable ${TPL_NAME}_FINDMOD already exists, then this TPL
        # has already been defined in a previous repository.  In this case, we
        # will just leave the TPL in its current position.
        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          message("-- " "NOTE: The TPL ${TPL_NAME} has already been defined so leaving it"
            " in the same location and not adding it again!")
        endif()
      else()
        list(APPEND ${PROJECT_NAME}_DEFINED_TPLS ${TPL_NAME})
      endif()

      # Set ${TPL_NAME}_TESTGROUP

      if (TPL_TESTGROUP STREQUAL PT
        OR TPL_TESTGROUP STREQUAL ST
        OR TPL_TESTGROUP STREQUAL TT
        OR TPL_TESTGROUP STREQUAL EX
        )
      else()
        message(FATAL_ERROR "Error the TPL classification '${TPL_TESTGROUP}'"
          " for the TPL ${TPL_NAME} is not a valid classification." )
      endif()

      if (NOT ${TPL_NAME}_TESTGROUP) # Allow for testing override
        set(${TPL_NAME}_TESTGROUP ${TPL_TESTGROUP})
      endif()

      # Set ${TPL_NAME}_FINDMOD

      if ("${REPOSITORY_DIR}" STREQUAL "." OR IS_ABSOLUTE ${TPL_FINDMOD})
        set(REPOSITORY_DIR_AND_SEP "")
      else()
        set(REPOSITORY_DIR_AND_SEP "${REPOSITORY_DIR}/")
      endif()

      set(TPL_FINDMOD "${REPOSITORY_DIR_AND_SEP}${TPL_FINDMOD}")

      set(TPL_FINDMOD_STD_NAME "FindTPL${TPL_NAME}.cmake")

      if (TPL_FINDMOD)
        string(REGEX MATCH ".+/$" FINDMOD_IS_DIR "${TPL_FINDMOD}")
        if (FINDMOD_IS_DIR)
          set(${TPL_NAME}_FINDMOD "${TPL_FINDMOD}${TPL_FINDMOD_STD_NAME}")
        else()
          set(${TPL_NAME}_FINDMOD ${TPL_FINDMOD})
        endif()
      else()
        set(${TPL_NAME}_FINDMOD ${TPL_FINDMOD_STD_NAME})
      endif()

      # Set ${TPL_NAME}_DEPENDENCIES_FILE

      # Breakdown ${${TPL_NAME}_FINDMOD} into parts <dir>/<base>.<ext>
      get_filename_component(tpl_findmod_dir ${${TPL_NAME}_FINDMOD} DIRECTORY)
      get_filename_component(tpl_findmod_base ${${TPL_NAME}_FINDMOD} NAME_WLE)
      get_filename_component(tpl_findmod_ext ${${TPL_NAME}_FINDMOD} LAST_EXT)

      set(${TPL_NAME}_DEPENDENCIES_FILE
        "${tpl_findmod_dir}/${tpl_findmod_base}Dependencies${tpl_findmod_ext}")

      # Set ${TPL_NAME}_TPLS_LIST_FILE

      assert_defined(${REPOSITORY_NAME}_TPLS_FILE)
      set(${TPL_NAME}_TPLS_LIST_FILE ${${REPOSITORY_NAME}_TPLS_FILE})

      # Set ${TPL_NAME}_PACKAGE_BUILD_STATUS

      set(${TPL_NAME}_PACKAGE_BUILD_STATUS  EXTERNAL)

      # Set ${TPL_NAME}_IS_TRIBITS_COMPLIANT

      if (${TPL_NAME}_FINDMOD STREQUAL "TRIBITS_PKG")
        set(${TPL_NAME}_IS_TRIBITS_COMPLIANT TRUE)
      else()
        set(${TPL_NAME}_IS_TRIBITS_COMPLIANT FALSE)
      endif()

      # Print variables/properties for the TPL

      if (${PROJECT_NAME}_VERBOSE_CONFIGURE  OR  TRIBITS_PROCESS_TPLS_LISTS_VERBOSE)
        print_var(${TPL_NAME}_TESTGROUP)
        print_var(${TPL_NAME}_FINDMOD)
        print_var(${TPL_NAME}_DEPENDENCIES_FILE)
        print_var(${TPL_NAME}_TPLS_LIST_FILE)
        print_var(${TPL_NAME}_PACKAGE_BUILD_STATUS)
        print_var(${TPL_NAME}_IS_TRIBITS_COMPLIANT)
      endif()

      # Set cache var TPL_ENABLE_${TPL_NAME} with default ""

      multiline_set(DOCSTR
        "Enable support for the TPL ${TPL_NAME} in all supported ${PROJECT_NAME} packages."
        "  This can be set to 'ON', 'OFF', or left empty ''."
        )
      set_cache_on_off_empty( TPL_ENABLE_${TPL_NAME} "" ${DOCSTR} )

      # 2008/11/25: rabartl: Above, we use the prefix TPL_ instead of
      # ${PROJECT_NAME}_ in order to make it clear that external TPLs are
      # different from packages so users don't get confused and think that the
      # project actually includes the source code and internal build for some
      # TPL when it does not!

    endforeach()

  endif()

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    print_var(${PROJECT_NAME}_DEFINED_TPLS)
  endif()

  # Get and print length
  list(LENGTH ${PROJECT_NAME}_DEFINED_TPLS ${PROJECT_NAME}_NUM_DEFINED_TPLS)
  message("-- After reading above TPLsList.cmake file: "
    "${PROJECT_NAME}_NUM_DEFINED_TPLS='${${PROJECT_NAME}_NUM_DEFINED_TPLS}'")

endmacro()
