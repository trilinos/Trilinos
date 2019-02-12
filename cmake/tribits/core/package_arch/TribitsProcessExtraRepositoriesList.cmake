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


INCLUDE(SetCacheOnOffEmpty)
INCLUDE(MultilineSet)
INCLUDE(AdvancedOption)
INCLUDE(AssertDefined)
INCLUDE(Split)
INCLUDE(MessageWrapper)
INCLUDE(TribitsSortListAccordingToMasterList)


#
# @MACRO: TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES()
#
# Declare a set of extra repositories for the `TriBITS Project`_ (i.e. in the
# project's `<projectDir>/cmake/ExtraRepositoriesList.cmake`_ file).
#
# Usage::
#
#   TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES(
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
#    supported for basic cloning and updating using `TRIBITS_CTEST_DRIVER()`_
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
# have upstream TPL and SE package dependencies listed before this package in
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
MACRO(TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES)
  ASSERT_DEFINED(PROJECT_NAME)
  IF ("${ARGN}" STREQUAL "")
    SET(${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY)
  ELSE()
    SET(${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY  "${ARGN}")
  ENDIF()
ENDMACRO()


#
# Field offsets
#

SET(ERP_REPO_NAME_OFFSET 0)
SET(ERP_REPO_DIR_OFFSET 1)
SET(ERP_REPO_VCTYPE_OFFSET 2)
SET(ERP_REPO_REPOURL_OFFSET 3)
SET(ERP_REPO_PACKSTAT_OFFSET 4)
SET(ERP_REPO_CLASSIFICATION_OFFSET 5)

SET(ERP_NUM_FIELDS_PER_REPO 6)

#
# Dump the list of extra repos in verbose mode
#
FUNCTION(TRIBITS_DUMP_EXTRA_REPOSITORIES_LIST)
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE OR TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
    PRINT_VAR(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT)
    PRINT_VAR(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS)
    PRINT_VAR(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES)
    PRINT_VAR(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS)
    PRINT_VAR(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS)
    PRINT_VAR(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS)
    PRINT_VAR(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES)
    PRINT_VAR(${PROJECT_NAME}_PRE_REPOSITORIES_DEFAULT)
    PRINT_VAR(${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT)
  ENDIF()
ENDFUNCTION()


#
# Function that parses the PACKSTAT field in the array
# ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY and returns
# HASPKGS and PREPOST.
#
FUNCTION(TRIBITS_PARSE_EXTRAREPO_PACKSTAT  PACKSTAT_IN
  HASPKGS_OUT  PREPOST_OUT
  )

  #PRINT_VAR(PACKSTAT_IN)
  SPLIT("${PACKSTAT_IN}"  ","  PACKSTAT_IN_ARRAY)
 
  # Set the defaults
  SET(HASPKGS  "HASPACKAGES")
  SET(PREPOST  "POST")

  FOREACH(PACKSTAT_ELE  ${PACKSTAT_IN_ARRAY})
    #PRINT_VAR(PACKSTAT_ELE)
    STRING(STRIP "${PACKSTAT_ELE}" PACKSTAT_ELE)
    #PRINT_VAR(PACKSTAT_ELE)
    IF (PACKSTAT_ELE  STREQUAL  "HASPACKAGES")
      SET(HASPKGS  "HASPACKAGES")
    ELSEIF (PACKSTAT_ELE  STREQUAL  "NOPACKAGES")
      SET(HASPKGS  "NOPACKAGES")
    ELSEIF (PACKSTAT_ELE  STREQUAL  "PRE")
      SET(PREPOST  "PRE")
    ELSEIF (PACKSTAT_ELE  STREQUAL  "POST")
      SET(PREPOST  "POST")
    ELSE()
      MESSAGE_WRAPPER(FATAL_ERROR  "Error, the value of 'PACKSTAT' element"
        " '${PACKSTAT_ELE}' is not valid!  Valid choices are '' (empty),"
        " 'HASPACKAGES', 'NOPACKAGES', 'PRE', and 'POST'.  The defaults if all"
        " fields are empty are 'HASPACKAGES' and 'POST'") 
    ENDIF()
  ENDFOREACH()
  # NOTE: In the above FOREACH(PACKSTAT_ELE ${PACKSTAT_IN_ARRAY}) loop, empty
  # elements are skipped!

  # Set the output arguments
  SET(${HASPKGS_OUT}  "${HASPKGS}"  PARENT_SCOPE)
  SET(${PREPOST_OUT}  "${PREPOST}"  PARENT_SCOPE)

ENDFUNCTION()

#
# Macro that processes the list varaible contents in
# ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY into sperate arrays:
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
MACRO(TRIBITS_PROCESS_EXTRAREPOS_LISTS)

  IF (
      ("${${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY}" STREQUAL "")
      AND
      (NOT "${${PROJECT_NAME}_EXTRAREPOS_DIR_REPOTYPE_REPOURL_PACKSTAT_CATEGORY}" STREQUAL "")
    )
    MESSAGE(WARNING "Warning! Usage of the variable"
     " '${PROJECT_NAME}_EXTRAREPOS_DIR_REPOTYPE_REPOURL_PACKSTAT_CATEGORY'"
     " is deprecated.  Please use the macro TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES()"
     " instead!")
    SET( ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      "${${PROJECT_NAME}_EXTRAREPOS_DIR_REPOTYPE_REPOURL_PACKSTAT_CATEGORY}" )
  ENDIF()

  # A) Get the total number of extrarepos defined

  IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
    PRINT_VAR(${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY)
  ENDIF()
  ASSERT_DEFINED(${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY)
  LIST(LENGTH ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
    ${PROJECT_NAME}_NUM_EXTRAREPOS_AND_FIELDS )
  MATH(EXPR ${PROJECT_NAME}_NUM_EXTRAREPOS
    "${${PROJECT_NAME}_NUM_EXTRAREPOS_AND_FIELDS}/${ERP_NUM_FIELDS_PER_REPO}")
  IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
    PRINT_VAR(${PROJECT_NAME}_NUM_EXTRAREPOS)
  ENDIF()
  MATH(EXPR ${PROJECT_NAME}_LAST_EXTRAREPO_IDX "${${PROJECT_NAME}_NUM_EXTRAREPOS}-1")
  IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
    PRINT_VAR(${PROJECT_NAME}_LAST_EXTRAREPO_IDX)
  ENDIF()

  # B) Process the list of extra repos

  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES)
  SET(${PROJECT_NAME}_PRE_REPOSITORIES_DEFAULT)
  SET(${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT)

  SET(PROCESSED_POST_EXTRAREPO  FALSE)

  FOREACH(EXTRAREPO_IDX RANGE ${${PROJECT_NAME}_LAST_EXTRAREPO_IDX})

    # B.1) Extract the fields for the current extrarepo row

    # NAME
    MATH(EXPR EXTRAREPO_NAME_IDX
      "${EXTRAREPO_IDX}*${ERP_NUM_FIELDS_PER_REPO}+${ERP_REPO_NAME_OFFSET}")
    IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      PRINT_VAR(EXTRAREPO_NAME_IDX)
    ENDIF()
    LIST(GET ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      ${EXTRAREPO_NAME_IDX} EXTRAREPO_NAME )
    IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      PRINT_VAR(EXTRAREPO_NAME)
    ENDIF()

    # DIR
    MATH(EXPR EXTRAREPO_DIR_IDX
      "${EXTRAREPO_IDX}*${ERP_NUM_FIELDS_PER_REPO}+${ERP_REPO_DIR_OFFSET}")
    IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      PRINT_VAR(EXTRAREPO_DIR_IDX)
    ENDIF()
    LIST(GET ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      ${EXTRAREPO_DIR_IDX} EXTRAREPO_DIR )
    IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      PRINT_VAR(EXTRAREPO_DIR)
    ENDIF()
    IF (EXTRAREPO_DIR STREQUAL "")
      SET(EXTRAREPO_DIR ${EXTRAREPO_NAME})
    ENDIF()
    IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      PRINT_VAR(EXTRAREPO_DIR)
    ENDIF()

    # VCTYPE
    MATH(EXPR EXTRAREPO_VCTYPE_IDX
      "${EXTRAREPO_IDX}*${ERP_NUM_FIELDS_PER_REPO}+${ERP_REPO_VCTYPE_OFFSET}")
    IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      PRINT_VAR(EXTRAREPO_VCTYPE_IDX)
    ENDIF()
    LIST(GET ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      ${EXTRAREPO_VCTYPE_IDX} EXTRAREPO_VCTYPE )
    IF (EXTRAREPO_VCTYPE STREQUAL GIT
      OR EXTRAREPO_VCTYPE STREQUAL SVN
      )
      # Okay
    ELSEIF (EXTRAREPO_VCTYPE  STREQUAL  HG)
      # not quite okay
      MESSAGE(WARNING "Warning: the repo ${EXTRAREPO_NAME} is a Mercurial repo: these are tolerated, but not fully supported.")
    ELSEIF (EXTRAREPO_VCTYPE  STREQUAL  "")
      # We are okay with no VC type
    ELSE()
      MESSAGE(FATAL_ERROR "Error, the repo type of '${EXTRAREPO_VCTYPE}' for"
        " extra repo ${EXTRAREPO_NAME} is *not* valid.  Valid choices are 'GIT', 'HG' and 'SVN'!")
    ENDIF()
    IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      PRINT_VAR(EXTRAREPO_VCTYPE)
    ENDIF()

    # REPOURL
    MATH(EXPR EXTRAREPO_REPOURL_IDX
      "${EXTRAREPO_IDX}*${ERP_NUM_FIELDS_PER_REPO}+${ERP_REPO_REPOURL_OFFSET}")
    IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      PRINT_VAR(EXTRAREPO_REPOURL_IDX)
    ENDIF()
    LIST(GET ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      ${EXTRAREPO_REPOURL_IDX} EXTRAREPO_REPOURL )
    IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      PRINT_VAR(EXTRAREPO_REPOURL)
    ENDIF()

    # PACKSTAT (PACKSTAT and PREPOST)
    MATH(EXPR EXTRAREPO_PACKSTAT_IDX
      "${EXTRAREPO_IDX}*${ERP_NUM_FIELDS_PER_REPO}+${ERP_REPO_PACKSTAT_OFFSET}")
    LIST(GET ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      ${EXTRAREPO_PACKSTAT_IDX} EXTRAREPO_PACKSTAT )
    IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      PRINT_VAR(EXTRAREPO_PACKSTAT)
    ENDIF()
    TRIBITS_PARSE_EXTRAREPO_PACKSTAT("${EXTRAREPO_PACKSTAT}"
      EXTRAREPO_HASPKGS  EXTRAREPO_PREPOST )

    # CLASSIFICATION
    MATH(EXPR EXTRAREPO_CLASSIFICATION_IDX
      "${EXTRAREPO_IDX}*${ERP_NUM_FIELDS_PER_REPO}+${ERP_REPO_CLASSIFICATION_OFFSET}")
    LIST(GET ${PROJECT_NAME}_EXTRAREPOS_DIR_VCTYPE_REPOURL_PACKSTAT_CATEGORY
      ${EXTRAREPO_CLASSIFICATION_IDX} EXTRAREPO_CLASSIFICATION )
    IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
      PRINT_VAR(EXTRAREPO_CLASSIFICATION)
    ENDIF()

    # Assert that PRE repos never come after a POST repo
    IF (NOT  PROCESSED_POST_EXTRAREPO  AND  EXTRAREPO_PREPOST  STREQUAL  "POST")
      SET(PROCESSED_POST_EXTRAREPO  TRUE)
    ELSEIF (PROCESSED_POST_EXTRAREPO  AND  EXTRAREPO_PREPOST  STREQUAL  "PRE")
      MESSAGE_WRAPPER(FATAL_ERROR  "Error, the 'PRE' extra repo '${EXTRAREPO_NAME}'"
        " specified in the PACKSTAT field '${EXTRAREPO_PACKSTAT}' came directly after"
        " a 'POST' extra repo!  All 'PRE' extra repos must be listed before all"
        " 'POST' extra repos!"
        ) 
    ENDIF()

    # B.2) Unconditionally add the extrarepo to the list

    LIST(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT ${EXTRAREPO_NAME})
    LIST(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS ${EXTRAREPO_DIR})
    LIST(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES ${EXTRAREPO_VCTYPE})
    LIST(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS ${EXTRAREPO_REPOURL})
    LIST(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS ${EXTRAREPO_HASPKGS})
    LIST(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS ${EXTRAREPO_PREPOST})
    LIST(APPEND ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES
      ${EXTRAREPO_CLASSIFICATION})
    IF (EXTRAREPO_PREPOST  STREQUAL  "PRE")
      LIST(APPEND ${PROJECT_NAME}_PRE_REPOSITORIES_DEFAULT ${EXTRAREPO_NAME})
    ELSEIF (EXTRAREPO_PREPOST  STREQUAL  "POST")
      LIST(APPEND ${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT ${EXTRAREPO_NAME})
    ELSE()
      MESSAGE(FATAL_ERROR "Error, bad value for EXTRAREPO_PREPOST!")
    ENDIF()

  ENDFOREACH()

  # C) Get the actual number of active extra repos

  LIST(LENGTH ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT ${PROJECT_NAME}_NUM_EXTRAREPOS )
  IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
    PRINT_VAR(${PROJECT_NAME}_NUM_EXTRAREPOS)
  ENDIF()
  MATH(EXPR ${PROJECT_NAME}_LAST_EXTRAREPO_IDX "${${PROJECT_NAME}_NUM_EXTRAREPOS}-1")

  # D) Print the final set of extrarepos in verbose mode

  TRIBITS_DUMP_EXTRA_REPOSITORIES_LIST()

ENDMACRO()


#
# Assert the existance and the order of the list of extra repositories in
# ${PROJECT_NAME}_PRE_REPOSITORIES listed in
# ${PROJECT_NAME}_EXTRA_REPOSITORIES according to the list read in from the
# extra repos file as determined by the varaible
# ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT.
#
FUNCTION(TRIBITS_EXTRA_REPOSITORIES_ASSERT_SUBSET_AND_ORDER_WRT_FILE)
  SET(ALL_EXTRA_REPOSITORIES_IN
    ${${PROJECT_NAME}_PRE_REPOSITORIES}  ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
  SET(ALL_EXTRA_REPOSITORIES_SORTED ${ALL_EXTRA_REPOSITORIES_IN})
  TRIBITS_SORT_LIST_ACCORDING_TO_MASTER_LIST(
    "${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT}" ALL_EXTRA_REPOSITORIES_SORTED)
  IF (NOT "${ALL_EXTRA_REPOSITORIES_IN}" STREQUAL "${ALL_EXTRA_REPOSITORIES_SORTED}")
    MESSAGE(FATAL_ERROR
      "ERROR!  The list of extra repos passed in '${ALL_EXTRA_REPOSITORIES_IN}'"
      " is not a subset and in the same order as read in from extra repos file"
      " '${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT}'"
      )
  ENDIF()
ENDFUNCTION()


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
# calling TRIBITS_EXTRA_REPOSITORIES_ASSERT_SUBSET_AND_ORDER_WRT_FILE()), then
# the set of repos and the associated data will be filtered based on
# ${PROJECT_NAME}_PRE_REPOSITORIES and ${PROJECT_NAME}_EXTRA_REPOSITORIES.
#
# If ${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES==TRUE, then the set of
# repos will be filtered based on what repos are present.  If
# ${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES==FALSE, then all of the
# repos must exist or MESSSAGE(FATAL_ERROR ...) is called and will fail the
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
FUNCTION(TRIBITS_FILTER_OR_ASSERT_EXTRA_REPOS)

  # Get the list of repos to filter

  IF ("${${PROJECT_NAME}_PRE_REPOSITORIES}"  STREQUAL "")
    SET(${PROJECT_NAME}_PRE_REPOSITORIES_IN  ${${PROJECT_NAME}_PRE_REPOSITORIES_DEFAULT})
  ELSE()
    SET(${PROJECT_NAME}_PRE_REPOSITORIES_IN  ${${PROJECT_NAME}_PRE_REPOSITORIES})
  ENDIF()

  IF ("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  STREQUAL "")
    SET(${PROJECT_NAME}_EXTRA_REPOSITORIES_IN  ${${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT})
  ELSE()
    SET(${PROJECT_NAME}_EXTRA_REPOSITORIES_IN  ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
  ENDIF()

  SET(ALL_EXTRA_REPOSITORIES_IN
    ${${PROJECT_NAME}_PRE_REPOSITORIES_IN} ${${PROJECT_NAME}_EXTRA_REPOSITORIES_IN})

  # Get out of function if there are no pre-extra or post-extra repos
  IF ("${ALL_EXTRA_REPOSITORIES_IN}"  STREQUAL  "")
    RETURN()
  ENDIF()

  # A) Loop through and copy info for existing repos to temp arrays

  SET(PRE_REPOSITORIES_TMP)
  SET(EXTRA_REPOSITORIES_TMP)
  SET(ALL_EXTRA_REPOSITORIES_TMP)
  SET(ALL_EXTRA_REPOSITORIES_DIRS_TMP)
  SET(ALL_EXTRA_REPOSITORIES_VCTYPES_TMP)
  SET(ALL_EXTRA_REPOSITORIES_REPOURLS_TMP)
  SET(ALL_EXTRA_REPOSITORIES_HASPKGS_TMP)
  SET(ALL_EXTRA_REPOSITORIES_PREPOSTS_TMP)
  SET(ALL_EXTRA_REPOSITORIES_CATEGORIES_TMP)

  # Set-up for filtering based on ALL_EXTRA_REPOSITORIES_IN != ""
  LIST(LENGTH  ALL_EXTRA_REPOSITORIES_IN  ALL_EXTRA_REPOSITORIES_IN_LEN)
  #PRINT_VAR(ALL_EXTRA_REPOSITORIES_IN_LEN)
  SET(EXTRAREPO_IN_IDX  0)
  LIST(GET  ALL_EXTRA_REPOSITORIES_IN  ${EXTRAREPO_IN_IDX}  EXTRAREPO_IN)

  # Loop over full list of extra repos from extra repos file
  SET(EXTRAREPO_IDX 0)
  FOREACH(EXTRAREPO_NAME ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT})

    #PRINT_VAR(EXTRAREPO_NAME)
    #PRINT_VAR(EXTRAREPO_IN)

    # A.1) Extract the data for current extra repo from file
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS ${EXTRAREPO_IDX}
      EXTRAREPO_DIR )
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES ${EXTRAREPO_IDX}
      EXTRAREPO_VCTYPE )
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS ${EXTRAREPO_IDX}
      EXTRAREPO_REPOURL )
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS ${EXTRAREPO_IDX}
      EXTRAREPO_HASPKGS )
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS ${EXTRAREPO_IDX}
      EXTRAREPO_PREPOST )
    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES ${EXTRAREPO_IDX}
      EXTRAREPO_CATEGORY )

    # A.2) Determine if to add the extra repo EXTRAREPO_NAME
    IF (EXTRAREPO_IN_IDX  EQUAL  ALL_EXTRA_REPOSITORIES_IN_LEN)
      # All of the extra repos in ALL_EXTRA_REPOSITORIES_IN have already
      # been processed.
      SET(ADD_EXTRAREPO  FALSE)
    ELSEIF (EXTRAREPO_IN  STREQUAL  EXTRAREPO_NAME)
      # We have a match, add the extra repo!
      SET(ADD_EXTRAREPO  TRUE)
      # Update EXTRAREPO_IN to look for next!
      MATH(EXPR  EXTRAREPO_IN_IDX  "${EXTRAREPO_IN_IDX}+1")
      IF (EXTRAREPO_IN_IDX  LESS  ALL_EXTRA_REPOSITORIES_IN_LEN)
        LIST(GET  ALL_EXTRA_REPOSITORIES_IN  ${EXTRAREPO_IN_IDX}  EXTRAREPO_IN)
      ELSE()
        # We have found the last repo already so move on
        SET(EXTRAREPO_IN  "")
      ENDIF()
    ELSE()
      # We are not at the end of the list in ALL_EXTRA_REPOSITORIES_IN yet
      # and have not reached the next entry in the list so don't add.
      SET(ADD_EXTRAREPO  FALSE)
    ENDIF()

    #PRINT_VAR(ADD_EXTRAREPO)

    # A.3) Determine the match of the category

    IF (ADD_EXTRAREPO)

      SET(ADD_EXTRAREPO  FALSE)

      #ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE)
      IF (TRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG)
        PRINT_VAR(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE)
      ENDIF()
  
      IF (${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE STREQUAL "Continuous" AND
          EXTRAREPO_CATEGORY STREQUAL "Continuous"
        )
        SET(ADD_EXTRAREPO TRUE)
      ELSEIF (${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE STREQUAL "Nightly" AND
          (EXTRAREPO_CATEGORY STREQUAL "Continuous"
             OR EXTRAREPO_CATEGORY STREQUAL "Nightly")
        )
        SET(ADD_EXTRAREPO TRUE)
      ENDIF()

    ENDIF()

    # A.3) Determine if the repo exists
    SET(EXTRAREPO_EXISTS  TRUE)
    IF (ADD_EXTRAREPO  AND  ${PROJECT_NAME}_CHECK_EXTRAREPOS_EXIST
      AND  NOT  UNITTEST_SKIP_FILTER_OR_ASSERT_EXTRA_REPOS
      )

      ASSERT_DEFINED(PROJECT_SOURCE_DIR)
      SET(EXTRAREPO_SOURCE_DIR  "${PROJECT_SOURCE_DIR}/${EXTRAREPO_DIR}")
      IF (EXISTS  "${EXTRAREPO_SOURCE_DIR}")
        SET(EXTRAREPO_EXISTS  TRUE)
      ELSE()
        SET(EXTRAREPO_EXISTS  FALSE)
      ENDIF()
      #PRINT_VAR(EXTRAREPO_EXISTS)

      IF (NOT  EXTRAREPO_EXISTS)
        SET(ADD_EXTRAREPO  FALSE)
        IF (${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES)
          MESSAGE("-- "
            "NOTE: Ignoring missing extra repo '${EXTRAREPO_NAME}'"
            " as requested since ${EXTRAREPO_SOURCE_DIR} does not exist" )
        ELSE()
          MESSAGE( FATAL_ERROR
            "ERROR!  Skipping missing extra repo '${EXTRAREPO_NAME}'"
            " since ${EXTRAREPO_SOURCE_DIR} does not exist!\n")
        ENDIF()
      ENDIF()

    ENDIF()

    #PRINT_VAR(ADD_EXTRAREPO)

    # A.4) Conditionally copy the info for the extra repo
    IF (ADD_EXTRAREPO)
      MESSAGE("-- " "Adding ${EXTRAREPO_PREPOST} extra ${EXTRAREPO_CATEGORY} repository ${EXTRAREPO_NAME} ...")
      IF (EXTRAREPO_PREPOST  STREQUAL  "PRE")
        LIST(APPEND  PRE_REPOSITORIES_TMP ${EXTRAREPO_NAME})
      ELSEIF (EXTRAREPO_PREPOST  STREQUAL  "POST")
        LIST(APPEND  EXTRA_REPOSITORIES_TMP ${EXTRAREPO_NAME})
      ENDIF()
      LIST(APPEND ALL_EXTRA_REPOSITORIES_TMP ${EXTRAREPO_NAME})
      LIST(APPEND ALL_EXTRA_REPOSITORIES_DIRS_TMP ${EXTRAREPO_DIR})
      LIST(APPEND ALL_EXTRA_REPOSITORIES_VCTYPES_TMP ${EXTRAREPO_VCTYPE})
      LIST(APPEND ALL_EXTRA_REPOSITORIES_REPOURLS_TMP ${EXTRAREPO_REPOURL})
      LIST(APPEND ALL_EXTRA_REPOSITORIES_HASPKGS_TMP ${EXTRAREPO_HASPKGS})
      LIST(APPEND ALL_EXTRA_REPOSITORIES_PREPOSTS_TMP ${EXTRAREPO_PREPOST})
      LIST(APPEND ALL_EXTRA_REPOSITORIES_CATEGORIES_TMP ${EXTRAREPO_CATEGORY})
    ELSE()
      LIST(APPEND CPACK_SOURCE_IGNORE_FILES
        "${${PROJECT_NAME}_SOURCE_DIR}/${EXTRAREPO_DIR}/")
      IF (EXTRAREPO_EXISTS)
        MESSAGE("-- " "*NOT* adding ${EXTRAREPO_PREPOST} extra ${EXTRAREPO_CATEGORY} repository ${EXTRAREPO_NAME}!")
      ENDIF()
   ENDIF()

    MATH(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")

  ENDFOREACH()

  # B) Copy over extra repos arrays with filtered arrays
  SET(${PROJECT_NAME}_PRE_REPOSITORIES
    ${PRE_REPOSITORIES_TMP}  PARENT_SCOPE)
  SET(${PROJECT_NAME}_EXTRA_REPOSITORIES
    ${EXTRA_REPOSITORIES_TMP}  PARENT_SCOPE)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES
    ${ALL_EXTRA_REPOSITORIES_TMP}  PARENT_SCOPE)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS
    ${ALL_EXTRA_REPOSITORIES_DIRS_TMP}  PARENT_SCOPE)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES
    ${ALL_EXTRA_REPOSITORIES_VCTYPES_TMP}  PARENT_SCOPE)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS
    ${ALL_EXTRA_REPOSITORIES_REPOURLS_TMP}  PARENT_SCOPE)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS
    ${ALL_EXTRA_REPOSITORIES_HASPKGS_TMP}  PARENT_SCOPE)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS
    ${ALL_EXTRA_REPOSITORIES_PREPOSTS_TMP}  PARENT_SCOPE)
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES
    ${ALL_EXTRA_REPOSITORIES_CATEGORIES_TMP}  PARENT_SCOPE)

  TRIBITS_DUMP_EXTRA_REPOSITORIES_LIST()

ENDFUNCTION()


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
# On output, the following varaibles are set:
#
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS
#   ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES
#
MACRO(TRIBITS_GET_AND_PROCESS_EXTRA_REPOSITORIES_LISTS)

  IF (${PROJECT_NAME}_EXTRAREPOS_FILE AND ${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE)

    #
    # A) Read in the extra repos list variable and process the list
    #

    MESSAGE("")
    MESSAGE("Reading the list of extra repositories from ${${PROJECT_NAME}_EXTRAREPOS_FILE}")
    MESSAGE("")

    INCLUDE(${${PROJECT_NAME}_EXTRAREPOS_FILE})

    TRIBITS_PROCESS_EXTRAREPOS_LISTS()
    # Above sets ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT (and DIRS
    # ... CATEGORIES)

    #
    # B) Sort and assert the list of extra repos according to the list read into the file
    #
    IF (
      (${PROJECT_NAME}_PRE_REPOSITORIES  OR  ${PROJECT_NAME}_EXTRA_REPOSITORIES)
      AND
      ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT
      )
      TRIBITS_EXTRA_REPOSITORIES_ASSERT_SUBSET_AND_ORDER_WRT_FILE()
    ENDIF()

    #
    # C) Filter out the missing extra repos or assert errors
    #
    IF (${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES)
      SET(SELECT_REPOS_MSG_SUFFIX "(ignoring missing repos)")
    ELSE()
      SET(SELECT_REPOS_MSG_SUFFIX "(asserting all selected repos exist)")
    ENDIF()
    MESSAGE("")
    MESSAGE("Selecting the set of '${${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE}' extra repos ${SELECT_REPOS_MSG_SUFFIX} ...")
    MESSAGE("")
    TRIBITS_FILTER_OR_ASSERT_EXTRA_REPOS()

  ELSE()

    IF (${PROJECT_NAME}_PRE_REPOSITORIES OR ${PROJECT_NAME}_EXTRA_REPOSITORIES)
      MESSAGE("")
      IF (${PROJECT_NAME}_PRE_REPOSITORIES)
        MESSAGE("Processing list of PRE extra repos from ${PROJECT_NAME}_PRE_REPOSITORIES"
          "='${${PROJECT_NAME}_PRE_REPOSITORIES}' ...")
      ENDIF()
      IF (${PROJECT_NAME}_EXTRA_REPOSITORIES)
        MESSAGE("Processing list of POST extra repos from ${PROJECT_NAME}_EXTRA_REPOSITORIES"
          "='${${PROJECT_NAME}_EXTRA_REPOSITORIES}' ...")
      ENDIF()
      MESSAGE("")
    ENDIF()

  ENDIF()

ENDMACRO()


#
# Extract the final name of the extra repo
#
FUNCTION(TRIBITS_GET_EXTRAREPO_BASE_NAME  EXTRAREPO_NAME EXTRAREPO_NAME_OUT)
  GET_FILENAME_COMPONENT(EXTRAREPO_NAME "${EXTRAREPO_NAME}" NAME)
  SET(${EXTRAREPO_NAME_OUT} "${EXTRAREPO_NAME}" PARENT_SCOPE)
ENDFUNCTION()
