# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include(SetCacheOnOffEmpty)
include(MultilineSet)
include(AdvancedOption)
include(AssertDefined)
include(Split)
include(MessageWrapper)
include(TribitsSortListAccordingToMasterList)


# @MACRO: tribits_project_define_extra_repositories()
#
# Declare a set of extra repositories for the `TriBITS Project`_ (i.e. in the
# project's `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file).
#
# Usage::
#
#   tribits_project_define_extra_repositories(
#     <repo0_name> <repo0_dir> <repo0_vctype> <repo0_url> <repo0_packstat> <repo0_classif>
#     <repo1_name> <repo1_dir> <repo1_vctype> <repo1_url> <rep10_packstat> <repo1_classif>
#     ...
#    )
#
# This macro takes in a 2D array with 6 columns, where each row defines an
# extra repository.  The 6 columns (ordered 0-5) are:
#
# 0. **REPO_NAME** (``<repoi_name>``): The name given to the repository
#    ``REPOSITORY_NAME``.
#
# 1. **REPO_DIR** (``<repoi_dir>``): The relative directory for the repository
#    under the project directory ``${PROJECT_SOURCE_DIR}`` (or
#    ``<projectDir>``).  If this is set to empty quoted string ``""``, then
#    the relative directory name is assumed to be same as the repository name
#    ``<repoi_name>``.  NOTE: If the repo is a `TriBITS Repository`_ (see
#    ``REPO_PACKSTAT`` below) then, currently, the repo dir must be the same
#    as the package name.
#
# 2. **REPO_VCTYPE** (``<repoi_vctype>``): The version control (VC) type of
#    the repo.  Value choices include ``GIT``, ``SVN`` (i.e. Subversion) or
#    empty ``""`` (in which case ``<repoi_url>`` must be empty as well).
#    WARNING: Only VC repos of type ``GIT`` can fully participate in the
#    TriBITS development tool workflows.  The other VC types are only
#    supported for basic cloning and updating using `tribits_ctest_driver()`_
#    scripts.
#
# 3. **REPO_URL** (``<repoi_url>``): The URL of the VC repo.  This info is
#    used to initially obtain the repo source code using the VC tool listed in
#    ``<repoi_vctype>``.  If the repo does not need to be cloned for the
#    needed use cases, then this can be the empty quoted string ``""``.  Also,
#    this field must be the empty string ``""`` if ``<repoi_vctype>`` is empty
#    ``""``.
#
# 4. **REPO_PACKSTAT** (``<repoi_packstat>``): Determines if the repository is
#    a `TriBITS Repository`_ and contains any TriBITS packages (or if it just
#    provides directories and files) and if its packages are listed before or
#    after the project's own native packages.  If this is a TriBITS Repository
#    (and therefore contains `TriBITS Repository Core Files`_) then this field
#    must contain the keyword ``HASPACKAGES`` or left empty.  If the listed
#    repository is **not** a TriBITS repository, and just provides directories
#    and files, then this field must contain the keyword ``NOPACKAGES``.  The
#    default is assumed to be ``HASPACKAGES`` if neither of these keywords are
#    provided.  If the keyword ``PRE`` is provided, then the TriBITS packages
#    in this repo come before the project's native packages.  If the keyword
#    ``POST`` is provided then the packages are listed after the project's
#    native packages. The default is assumed to be ``POST`` if neither of
#    these keywords are provided.  The keywords must be separated by a comma
#    with no spaces such as with "``PRE,HASPACKAGES``",
#    "``HASPACKAGES,POST``", "``POST,NOPACKAGES``", etc.  If no keywords are
#    provided, then the empty string "" must be used (which defaults to
#    ``"HASPACKAGES,POST"``).
#
# 5. **REPO_CLASSIFICATION** (``<repoi_classif>``): Gives the `Repository Test
#    Classification`_ which also happens to be the CTest/CDash testing mode
#    and the default dashboard track.  Valid values are ``Continuous``,
#    ``Nightly``, and ``Experimental``.  See `Repository Test Classification`_
#    for a detailed description.
#
# This command is used to put together one or more VC and/or TriBITS
# repositories to construct a composite `TriBITS Project`_.  The option
# `<Project>_EXTRAREPOS_FILE`_ is used to point to files that call this macro.
#
# Repositories with ``<repoi_packstat>=NOPACKAGES`` are **not** TriBITS
# Repositories and are technically not considered at all during the basic
# configuration of the a TriBITS project.  They are only listed in this file
# so that they can be used in the version control logic for tools that perform
# version control with the repositories (such as getting git versions,
# cloning, updating, looking for changed files, etc.).  For example, a
# non-TriBITS repo can be used to grab a set of directories and files that
# fill in the definition of a package in an upstream repository (see `How to
# insert a package into an upstream repo`_).  Also, non-TriBITS repos can be
# used to provide extra test data for a given package or a set of packages so
# that extra tests can be run.
#
# Repositories with ``<repoi_repotype>=''`` are not VC repos.  This can be
# used, for example, to represent the project's native repos or it can be used
# to point to a TriBITS repository that was cloned in an early listed VC repo.
#
# NOTE: These repositories must be listed in the order of package
# dependencies.  That is, all of the packages listed in repository ``i`` must
# have upstream TPL and package dependencies listed before this package in
# this repository or in upstream repositories ``i-1``, ``i-2``, etc.
#
# NOTE: This module just sets the local variable::
#
#  ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
#
# in the current scope.  The advantages of using this macro instead of
# directly setting this variable are that the macro:
#
# * Asserts that the variable ``PROJECT_NAME`` is defined and set.
#
# * Avoids misspelling the name of the variable
#   ``${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY``.  If
#   one misspells the name of a macro, it is an immediate error in CMake.  A
#   misspelled set variable is just ignored.
#
# * The variable name can change in the future as an implementation detail.
#
macro(tribits_project_define_extra_repositories)
  assert_defined(PROJECT_NAME)
  if ("${ARGN}" STREQUAL "")
    set(${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY)
  else()
    set(${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY  "${ARGN}")
  endif()
endmacro()


#
# Field offsets
#

set(ERP_REPO_NAME_OFFSET 0)
set(ERP_REPO_DIR_OFFSET 1)
set(ERP_REPO_VCTYPE_OFFSET 2)
set(ERP_REPO_REPOURL_OFFSET 3)
set(ERP_REPO_PACKSTAT_OFFSET 4)
set(ERP_REPO_CLASSIFICATION_OFFSET 5)

set(ERP_NUM_FIELDS_PER_REPO 6)

#
# Dump the list of extra repos in verbose mode
#
function(tribits_dump_extra_repositories_list)
  if (${PROJECT_NAME}_VERBOSE_CONFIGURE OR TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
    print_var(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT)
    print_var(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS)
    print_var(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES)
    print_var(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS)
    print_var(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS)
    print_var(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS)
    print_var(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES)
    print_var(${PROJECT_NAME}_PRE_REPOSITORIES_DEFAULT)
    print_var(${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT)
  endif()
endfunction()


#
# Function that parses the PACKSTAT field in the array
# ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY and returns
# HASPKGS and PREPOST.
#
function(tribits_parse_extrarepo_packstat  PACKSTAT_IN
  HASPKGS_OUT  PREPOST_OUT
  )

  #print_var(PACKSTAT_IN)
  split("${PACKSTAT_IN}"  ","  PACKSTAT_IN_ARRAY)
 
  # Set the defaults
  set(HASPKGS  "HASPACKAGES")
  set(PREPOST  "POST")

  foreach(PACKSTAT_ELE  ${PACKSTAT_IN_ARRAY})
    #print_var(PACKSTAT_ELE)
    string(STRIP "${PACKSTAT_ELE}" PACKSTAT_ELE)
    #print_var(PACKSTAT_ELE)
    if (PACKSTAT_ELE  STREQUAL  "HASPACKAGES")
      set(HASPKGS  "HASPACKAGES")
    elseif (PACKSTAT_ELE  STREQUAL  "NOPACKAGES")
      set(HASPKGS  "NOPACKAGES")
    elseif (PACKSTAT_ELE  STREQUAL  "PRE")
      set(PREPOST  "PRE")
    elseif (PACKSTAT_ELE  STREQUAL  "POST")
      set(PREPOST  "POST")
    else()
      message_wrapper(FATAL_ERROR  "Error, the value of 'PACKSTAT' element"
        " '${PACKSTAT_ELE}' is not valid!  Valid choices are '' (empty),"
        " 'HASPACKAGES', 'NOPACKAGES', 'PRE', and 'POST'.  The defaults if all"
        " fields are empty are 'HASPACKAGES' and 'POST'") 
    endif()
  endforeach()
  # NOTE: In the above foreach(PACKSTAT_ELE ${PACKSTAT_IN_ARRAY}) loop, empty
  # elements are skipped!

  # Set the output arguments
  set(${HASPKGS_OUT}  "${HASPKGS}"  PARENT_SCOPE)
  set(${PREPOST_OUT}  "${PREPOST}"  PARENT_SCOPE)

endfunction()

#
# Macro that processes the list variable contents in
# ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY into separate arrays:
#
#   ${PROJECT_NAME}_PRE_REPOSITORIES_DEFAULT
#   ${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES
#
# The macro responds to ${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE
# to match the categories.
#
macro(tribits_process_extrarepos_lists)

  if (
      ("${${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY}" STREQUAL "")
      AND
      (NOT "${${PROJECT_NAME}_EXTRAREPOS_DIR_REPOTYPE_REPOURL_PACKSTAT_CATEGORY}" STREQUAL "")
    )
    message(WARNING "Warning! Usage of the variable"
     " '${PROJECT_NAME}_EXTRAREPOS_DIR_REPOTYPE_REPOURL_PACKSTAT_CATEGORY'"
     " is deprecated.  Please use the macro tribits_project_define_extra_repositories()"
     " instead!")
    set( ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      "${${PROJECT_NAME}_EXTRAREPOS_DIR_REPOTYPE_REPOURL_PACKSTAT_CATEGORY}" )
  endif()

  # A) Get the total number of extrarepos defined

  if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
    print_var(${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY)
  endif()
  assert_defined(${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY)
  list(LENGTH ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
    ${PROJECT_NAME}_NUM_EXTRAREPOS_AND_FIELDS )
  math(EXPR ${PROJECT_NAME}_NUM_EXTRAREPOS
    "${${PROJECT_NAME}_NUM_EXTRAREPOS_AND_FIELDS}/${ERP_NUM_FIELDS_PER_REPO}")
  if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
    print_var(${PROJECT_NAME}_NUM_EXTRAREPOS)
  endif()
  math(EXPR ${PROJECT_NAME}_LAST_EXTRAREPO_IDX "${${PROJECT_NAME}_NUM_EXTRAREPOS}-1")
  if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
    print_var(${PROJECT_NAME}_LAST_EXTRAREPO_IDX)
  endif()

  # B) Process the list of extra repos

  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES)
  set(${PROJECT_NAME}_PRE_REPOSITORIES_DEFAULT)
  set(${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT)

  set(PROCESSED_POST_EXTRAREPO  FALSE)

  foreach(EXTRAREPO_IDX RANGE ${${PROJECT_NAME}_LAST_EXTRAREPO_IDX})

    # B.1) Extract the fields for the current extrarepo row

    # NAME
    math(EXPR EXTRAREPO_NAME_IDX
      "${EXTRAREPO_IDX}*${ERP_NUM_FIELDS_PER_REPO}+${ERP_REPO_NAME_OFFSET}")
    if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      print_var(EXTRAREPO_NAME_IDX)
    endif()
    list(GET ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      ${EXTRAREPO_NAME_IDX} EXTRAREPO_NAME )
    if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      print_var(EXTRAREPO_NAME)
    endif()

    # DIR
    math(EXPR EXTRAREPO_DIR_IDX
      "${EXTRAREPO_IDX}*${ERP_NUM_FIELDS_PER_REPO}+${ERP_REPO_DIR_OFFSET}")
    if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      print_var(EXTRAREPO_DIR_IDX)
    endif()
    list(GET ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      ${EXTRAREPO_DIR_IDX} EXTRAREPO_DIR )
    if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      print_var(EXTRAREPO_DIR)
    endif()
    if (EXTRAREPO_DIR STREQUAL "")
      set(EXTRAREPO_DIR ${EXTRAREPO_NAME})
    endif()
    if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      print_var(EXTRAREPO_DIR)
    endif()

    # VCTYPE
    math(EXPR EXTRAREPO_VCTYPE_IDX
      "${EXTRAREPO_IDX}*${ERP_NUM_FIELDS_PER_REPO}+${ERP_REPO_VCTYPE_OFFSET}")
    if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      print_var(EXTRAREPO_VCTYPE_IDX)
    endif()
    list(GET ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      ${EXTRAREPO_VCTYPE_IDX} EXTRAREPO_VCTYPE )
    if (EXTRAREPO_VCTYPE STREQUAL GIT
      OR EXTRAREPO_VCTYPE STREQUAL SVN
      )
      # Okay
    elseif (EXTRAREPO_VCTYPE  STREQUAL  HG)
      # not quite okay
      message(WARNING "Warning: the repo ${EXTRAREPO_NAME} is a Mercurial repo: these are tolerated, but not fully supported.")
    elseif (EXTRAREPO_VCTYPE  STREQUAL  "")
      # We are okay with no VC type
    else()
      message(FATAL_ERROR "Error, the repo type of '${EXTRAREPO_VCTYPE}' for"
        " extra repo ${EXTRAREPO_NAME} is *not* valid.  Valid choices are 'GIT', 'HG' and 'SVN'!")
    endif()
    if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      print_var(EXTRAREPO_VCTYPE)
    endif()

    # REPOURL
    math(EXPR EXTRAREPO_REPOURL_IDX
      "${EXTRAREPO_IDX}*${ERP_NUM_FIELDS_PER_REPO}+${ERP_REPO_REPOURL_OFFSET}")
    if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      print_var(EXTRAREPO_REPOURL_IDX)
    endif()
    list(GET ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      ${EXTRAREPO_REPOURL_IDX} EXTRAREPO_REPOURL )
    if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      print_var(EXTRAREPO_REPOURL)
    endif()

    # PACKSTAT (PACKSTAT and PREPOST)
    math(EXPR EXTRAREPO_PACKSTAT_IDX
      "${EXTRAREPO_IDX}*${ERP_NUM_FIELDS_PER_REPO}+${ERP_REPO_PACKSTAT_OFFSET}")
    list(GET ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      ${EXTRAREPO_PACKSTAT_IDX} EXTRAREPO_PACKSTAT )
    if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      print_var(EXTRAREPO_PACKSTAT)
    endif()
    tribits_parse_extrarepo_packstat("${EXTRAREPO_PACKSTAT}"
      EXTRAREPO_HASPKGS  EXTRAREPO_PREPOST )

    # CLASSIFICATION
    math(EXPR EXTRAREPO_CLASSIFICATION_IDX
      "${EXTRAREPO_IDX}*${ERP_NUM_FIELDS_PER_REPO}+${ERP_REPO_CLASSIFICATION_OFFSET}")
    list(GET ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      ${EXTRAREPO_CLASSIFICATION_IDX} EXTRAREPO_CLASSIFICATION )
    if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      print_var(EXTRAREPO_CLASSIFICATION)
    endif()

    # Assert that PRE repos never come after a POST repo
    if (NOT  PROCESSED_POST_EXTRAREPO  AND  EXTRAREPO_PREPOST  STREQUAL  "POST")
      set(PROCESSED_POST_EXTRAREPO  TRUE)
    elseif (PROCESSED_POST_EXTRAREPO  AND  EXTRAREPO_PREPOST  STREQUAL  "PRE")
      message_wrapper(FATAL_ERROR  "Error, the 'PRE' extra repo '${EXTRAREPO_NAME}'"
        " specified in the PACKSTAT field '${EXTRAREPO_PACKSTAT}' came directly after"
        " a 'POST' extra repo!  All 'PRE' extra repos must be listed before all"
        " 'POST' extra repos!"
        ) 
    endif()

    # B.2) Unconditionally add the extrarepo to the list

    list(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT ${EXTRAREPO_NAME})
    list(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS ${EXTRAREPO_DIR})
    list(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES ${EXTRAREPO_VCTYPE})
    list(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS ${EXTRAREPO_REPOURL})
    list(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS ${EXTRAREPO_HASPKGS})
    list(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS ${EXTRAREPO_PREPOST})
    list(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES
      ${EXTRAREPO_CLASSIFICATION})
    if (EXTRAREPO_PREPOST  STREQUAL  "PRE")
      list(APPEND ${PROJECT_NAME}_PRE_REPOSITORIES_DEFAULT ${EXTRAREPO_NAME})
    elseif (EXTRAREPO_PREPOST  STREQUAL  "POST")
      list(APPEND ${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT ${EXTRAREPO_NAME})
    else()
      message(FATAL_ERROR "Error, bad value for EXTRAREPO_PREPOST!")
    endif()

  endforeach()

  # C) Get the actual number of active extra repos

  list(LENGTH ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT ${PROJECT_NAME}_NUM_EXTRAREPOS )
  if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
    print_var(${PROJECT_NAME}_NUM_EXTRAREPOS)
  endif()
  math(EXPR ${PROJECT_NAME}_LAST_EXTRAREPO_IDX "${${PROJECT_NAME}_NUM_EXTRAREPOS}-1")

  # D) Print the final set of extrarepos in verbose mode

  tribits_dump_extra_repositories_list()

endmacro()


#
# Assert the existence and the order of the list of extra repositories in
# ${PROJECT_NAME}_PRE_REPOSITORIES listed in
# ${PROJECT_NAME}_EXTRA_REPOSITORIES according to the list read in from the
# extra repos file as determined by the variable
# ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT.
#
function(tribits_extra_repositories_assert_subset_and_order_wrt_file)
  set(ALL_EXTRA_REPOSITORIES_IN
    ${${PROJECT_NAME}_PRE_REPOSITORIES}  ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
  set(ALL_EXTRA_REPOSITORIES_SORTED ${ALL_EXTRA_REPOSITORIES_IN})
  tribits_sort_list_according_to_master_list(
    "${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT}" ALL_EXTRA_REPOSITORIES_SORTED)
  if (NOT "${ALL_EXTRA_REPOSITORIES_IN}" STREQUAL "${ALL_EXTRA_REPOSITORIES_SORTED}")
    message(FATAL_ERROR
      "ERROR!  The list of extra repos passed in '${ALL_EXTRA_REPOSITORIES_IN}'"
      " is not a subset and in the same order as read in from extra repos file"
      " '${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT}'"
      )
  endif()
endfunction()


#
# Filter out or assert msising repos read from an extra repos.
#
# This function keys off of the variables:
#
#  ${PROJECT_NAME}_PRE_REPOSITORIES
#  ${PROJECT_NAME}_EXTRA_REPOSITORIES
#  ${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES
#
# and the variables:
#
#   ${PROJECT_NAME}_PRE_REPOSITORIES_DEFAULT
#   ${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES
#
# which contain the extra repos read from the extra repos file.
#
# If ${PROJECT_NAME}_PRE_REPOSITORIES or ${PROJECT_NAME}_EXTRA_REPOSITORIES
# are non-empty (it is assumed that the extra repos listed are a subset of
# ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT, which can be asserted so by
# calling tribits_extra_repositories_assert_subset_and_order_wrt_file()), then
# the set of repos and the associated data will be filtered based on
# ${PROJECT_NAME}_PRE_REPOSITORIES and ${PROJECT_NAME}_EXTRA_REPOSITORIES.
#
# If ${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES==TRUE, then the set of
# repos will be filtered based on what repos are present.  If
# ${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES==FALSE, then all of the
# repos must exist or message(FATAL_ERROR ...) is called and will fail the
# configure.
#
# On output the variables:
#
#   ${PROJECT_NAME}_PRE_REPOSITORIES
#   ${PROJECT_NAME}_EXTRA_REPOSITORIES
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES
#
# will be set according to the logic described above and the other extra repo
# variables will be filtered in a consistent way.
#
function(tribits_filter_or_assert_extra_repos)

  # Get the list of repos to filter

  if ("${${PROJECT_NAME}_PRE_REPOSITORIES}"  STREQUAL "")
    set(${PROJECT_NAME}_PRE_REPOSITORIES_IN  ${${PROJECT_NAME}_PRE_REPOSITORIES_DEFAULT})
  else()
    set(${PROJECT_NAME}_PRE_REPOSITORIES_IN  ${${PROJECT_NAME}_PRE_REPOSITORIES})
  endif()

  if ("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  STREQUAL "")
    set(${PROJECT_NAME}_EXTRA_REPOSITORIES_IN  ${${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT})
  else()
    set(${PROJECT_NAME}_EXTRA_REPOSITORIES_IN  ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
  endif()

  set(ALL_EXTRA_REPOSITORIES_IN
    ${${PROJECT_NAME}_PRE_REPOSITORIES_IN} ${${PROJECT_NAME}_EXTRA_REPOSITORIES_IN})

  # Get out of function if there are no pre-extra or post-extra repos
  if ("${ALL_EXTRA_REPOSITORIES_IN}"  STREQUAL  "")
    return()
  endif()

  # A) Loop through and copy info for existing repos to temp arrays

  set(PRE_REPOSITORIES_TMP)
  set(EXTRA_REPOSITORIES_TMP)
  set(ALL_EXTRA_REPOSITORIES_TMP)
  set(ALL_EXTRA_REPOSITORIES_DIRS_TMP)
  set(ALL_EXTRA_REPOSITORIES_VCTYPES_TMP)
  set(ALL_EXTRA_REPOSITORIES_REPOURLS_TMP)
  set(ALL_EXTRA_REPOSITORIES_HASPKGS_TMP)
  set(ALL_EXTRA_REPOSITORIES_PREPOSTS_TMP)
  set(ALL_EXTRA_REPOSITORIES_CATEGORIES_TMP)

  # Set-up for filtering based on ALL_EXTRA_REPOSITORIES_IN != ""
  list(LENGTH  ALL_EXTRA_REPOSITORIES_IN  ALL_EXTRA_REPOSITORIES_IN_LEN)
  #print_var(ALL_EXTRA_REPOSITORIES_IN_LEN)
  set(EXTRAREPO_IN_IDX  0)
  list(GET  ALL_EXTRA_REPOSITORIES_IN  ${EXTRAREPO_IN_IDX}  EXTRAREPO_IN)

  # Loop over full list of extra repos from extra repos file
  set(EXTRAREPO_IDX 0)
  foreach(EXTRAREPO_NAME ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT})

    #print_var(EXTRAREPO_NAME)
    #print_var(EXTRAREPO_IN)

    # A.1) Extract the data for current extra repo from file
    list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS ${EXTRAREPO_IDX}
      EXTRAREPO_DIR )
    list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES ${EXTRAREPO_IDX}
      EXTRAREPO_VCTYPE )
    list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS ${EXTRAREPO_IDX}
      EXTRAREPO_REPOURL )
    list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS ${EXTRAREPO_IDX}
      EXTRAREPO_HASPKGS )
    list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS ${EXTRAREPO_IDX}
      EXTRAREPO_PREPOST )
    list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES ${EXTRAREPO_IDX}
      EXTRAREPO_CATEGORY )

    # A.2) Determine if to add the extra repo EXTRAREPO_NAME
    if (EXTRAREPO_IN_IDX  EQUAL  ALL_EXTRA_REPOSITORIES_IN_LEN)
      # All of the extra repos in ALL_EXTRA_REPOSITORIES_IN have already
      # been processed.
      set(ADD_EXTRAREPO  FALSE)
    elseif (EXTRAREPO_IN  STREQUAL  EXTRAREPO_NAME)
      # We have a match, add the extra repo!
      set(ADD_EXTRAREPO  TRUE)
      # Update EXTRAREPO_IN to look for next!
      math(EXPR  EXTRAREPO_IN_IDX  "${EXTRAREPO_IN_IDX}+1")
      if (EXTRAREPO_IN_IDX  LESS  ALL_EXTRA_REPOSITORIES_IN_LEN)
        list(GET  ALL_EXTRA_REPOSITORIES_IN  ${EXTRAREPO_IN_IDX}  EXTRAREPO_IN)
      else()
        # We have found the last repo already so move on
        set(EXTRAREPO_IN  "")
      endif()
    else()
      # We are not at the end of the list in ALL_EXTRA_REPOSITORIES_IN yet
      # and have not reached the next entry in the list so don't add.
      set(ADD_EXTRAREPO  FALSE)
    endif()

    #print_var(ADD_EXTRAREPO)

    # A.3) Determine the match of the category

    if (ADD_EXTRAREPO)

      set(ADD_EXTRAREPO  FALSE)

      #assert_defined(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE)
      if (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
        print_var(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE)
      endif()
  
      if (${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE STREQUAL "Continuous" AND
          EXTRAREPO_CATEGORY STREQUAL "Continuous"
        )
        set(ADD_EXTRAREPO TRUE)
      elseif (${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE STREQUAL "Nightly" AND
          (EXTRAREPO_CATEGORY STREQUAL "Continuous"
             OR EXTRAREPO_CATEGORY STREQUAL "Nightly")
        )
        set(ADD_EXTRAREPO TRUE)
      endif()

    endif()

    # A.3) Determine if the repo exists
    set(EXTRAREPO_EXISTS  TRUE)
    if (ADD_EXTRAREPO  AND  ${PROJECT_NAME}_CHECK_EXTRAREPOS_EXIST
      AND  NOT  UNITTEST_SKIP_FILTER_OR_ASSERT_EXTRA_REPOS
      )

      assert_defined(PROJECT_SOURCE_DIR)
      set(EXTRAREPO_SOURCE_DIR  "${PROJECT_SOURCE_DIR}/${EXTRAREPO_DIR}")
      if (EXISTS  "${EXTRAREPO_SOURCE_DIR}")
        set(EXTRAREPO_EXISTS  TRUE)
      else()
        set(EXTRAREPO_EXISTS  FALSE)
      endif()
      #print_var(EXTRAREPO_EXISTS)

      if (NOT  EXTRAREPO_EXISTS)
        set(ADD_EXTRAREPO  FALSE)
        if (${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES)
          message("-- "
            "NOTE: Ignoring missing extra repo '${EXTRAREPO_NAME}'"
            " as requested since ${EXTRAREPO_SOURCE_DIR} does not exist" )
        else()
          message( FATAL_ERROR
            "ERROR!  Skipping missing extra repo '${EXTRAREPO_NAME}'"
            " since ${EXTRAREPO_SOURCE_DIR} does not exist!\n")
        endif()
      endif()

    endif()

    #print_var(ADD_EXTRAREPO)

    # A.4) Conditionally copy the info for the extra repo
    if (ADD_EXTRAREPO)
      message("-- " "Adding ${EXTRAREPO_PREPOST} extra ${EXTRAREPO_CATEGORY} repository ${EXTRAREPO_NAME} ...")
      if (EXTRAREPO_PREPOST  STREQUAL  "PRE")
        list(APPEND  PRE_REPOSITORIES_TMP ${EXTRAREPO_NAME})
      elseif (EXTRAREPO_PREPOST  STREQUAL  "POST")
        list(APPEND  EXTRA_REPOSITORIES_TMP ${EXTRAREPO_NAME})
      endif()
      list(APPEND ALL_EXTRA_REPOSITORIES_TMP ${EXTRAREPO_NAME})
      list(APPEND ALL_EXTRA_REPOSITORIES_DIRS_TMP ${EXTRAREPO_DIR})
      list(APPEND ALL_EXTRA_REPOSITORIES_VCTYPES_TMP ${EXTRAREPO_VCTYPE})
      list(APPEND ALL_EXTRA_REPOSITORIES_REPOURLS_TMP ${EXTRAREPO_REPOURL})
      list(APPEND ALL_EXTRA_REPOSITORIES_HASPKGS_TMP ${EXTRAREPO_HASPKGS})
      list(APPEND ALL_EXTRA_REPOSITORIES_PREPOSTS_TMP ${EXTRAREPO_PREPOST})
      list(APPEND ALL_EXTRA_REPOSITORIES_CATEGORIES_TMP ${EXTRAREPO_CATEGORY})
    else()
      list(APPEND CPACK_SOURCE_IGNORE_FILES
        "${${PROJECT_NAME}_SOURCE_DIR}/${EXTRAREPO_DIR}/")
      if (EXTRAREPO_EXISTS)
        message("-- " "*NOT* adding ${EXTRAREPO_PREPOST} extra ${EXTRAREPO_CATEGORY} repository ${EXTRAREPO_NAME}!")
      endif()
   endif()

    math(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")

  endforeach()

  # B) Copy over extra repos arrays with filtered arrays
  set(${PROJECT_NAME}_PRE_REPOSITORIES
    ${PRE_REPOSITORIES_TMP}  PARENT_SCOPE)
  set(${PROJECT_NAME}_EXTRA_REPOSITORIES
    ${EXTRA_REPOSITORIES_TMP}  PARENT_SCOPE)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES
    ${ALL_EXTRA_REPOSITORIES_TMP}  PARENT_SCOPE)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS
    ${ALL_EXTRA_REPOSITORIES_DIRS_TMP}  PARENT_SCOPE)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES
    ${ALL_EXTRA_REPOSITORIES_VCTYPES_TMP}  PARENT_SCOPE)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS
    ${ALL_EXTRA_REPOSITORIES_REPOURLS_TMP}  PARENT_SCOPE)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS
    ${ALL_EXTRA_REPOSITORIES_HASPKGS_TMP}  PARENT_SCOPE)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS
    ${ALL_EXTRA_REPOSITORIES_PREPOSTS_TMP}  PARENT_SCOPE)
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES
    ${ALL_EXTRA_REPOSITORIES_CATEGORIES_TMP}  PARENT_SCOPE)

  tribits_dump_extra_repositories_list()

endfunction()


#
# Macro that reads extra repos file, processes the list of extra repos, etc.
#
# On input, the following variables are read:
#
#   ${PROJECT_NAME}_EXTRAREPOS_FILE
#   ${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES
#   ${PROJECT_NAME}_PRE_REPOSITORIES
#   ${PROJECT_NAME}_EXTRA_REPOSITORIES
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES
#
# On output, the following variables are set:
#
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES
#
macro(tribits_get_and_process_extra_repositories_lists)

  if (${PROJECT_NAME}_EXTRAREPOS_FILE AND ${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE)

    #
    # A) Read in the extra repos list variable and process the list
    #

    message("")
    message("Reading the list of extra repositories from ${${PROJECT_NAME}_EXTRAREPOS_FILE}")
    message("")

    include(${${PROJECT_NAME}_EXTRAREPOS_FILE})

    tribits_process_extrarepos_lists()
    # Above sets ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT (and DIRS
    # ... CATEGORIES)

    #
    # B) Sort and assert the list of extra repos according to the list read into the file
    #
    if (
      (${PROJECT_NAME}_PRE_REPOSITORIES  OR  ${PROJECT_NAME}_EXTRA_REPOSITORIES)
      AND
      ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT
      )
      tribits_extra_repositories_assert_subset_and_order_wrt_file()
    endif()

    #
    # C) Filter out the missing extra repos or assert errors
    #
    if (${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES)
      set(SELECT_REPOS_MSG_SUFFIX "(ignoring missing repos)")
    else()
      set(SELECT_REPOS_MSG_SUFFIX "(asserting all selected repos exist)")
    endif()
    message("")
    message("Selecting the set of '${${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE}' extra repos ${SELECT_REPOS_MSG_SUFFIX} ...")
    message("")
    tribits_filter_or_assert_extra_repos()

  else()

    if (${PROJECT_NAME}_PRE_REPOSITORIES OR ${PROJECT_NAME}_EXTRA_REPOSITORIES)
      message("")
      if (${PROJECT_NAME}_PRE_REPOSITORIES)
        message("Processing list of PRE extra repos from ${PROJECT_NAME}_PRE_REPOSITORIES"
          "='${${PROJECT_NAME}_PRE_REPOSITORIES}' ...")
      endif()
      if (${PROJECT_NAME}_EXTRA_REPOSITORIES)
        message("Processing list of POST extra repos from ${PROJECT_NAME}_EXTRA_REPOSITORIES"
          "='${${PROJECT_NAME}_EXTRA_REPOSITORIES}' ...")
      endif()
      message("")
    endif()

  endif()

endmacro()


#
# Extract the final name of the extra repo
#
function(tribits_get_extrarepo_base_name  EXTRAREPO_NAME EXTRAREPO_NAME_OUT)
  get_filename_component(EXTRAREPO_NAME "${EXTRAREPO_NAME}" NAME)
  set(${EXTRAREPO_NAME_OUT} "${EXTRAREPO_NAME}" PARENT_SCOPE)
endfunction()
