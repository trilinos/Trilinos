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
# Tribits platform-independent test driver.
#
#

MESSAGE("")
MESSAGE("*******************************")
MESSAGE("*** TribitsCTestDriverCore ***")
MESSAGE("*******************************")
MESSAGE("")


CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11 FATAL_ERROR)

#
# Get the basic variables that define the project and the build
#

#
# Set TRIBITS_PROJECT_ROOT
#
# We must locate the source code root directory before processing this
# file. Once the root directory is located, we can make good guesses
# at other properties, but this file is a CTest script that is always
# executed outside of the context of the project's CMake build.
#
# We allow the environment variable TRIBITS_PROJECT_ROOT to locate the
# root directory. If the variable doesn't exist, we fall back on the
# default convention.
#MESSAGE("TRIBITS_PROJECT_ROOT (before env) = '${TRIBITS_PROJECT_ROOT}'")
IF (NOT TRIBITS_PROJECT_ROOT)
  SET(TRIBITS_PROJECT_ROOT "$ENV{TRIBITS_PROJECT_ROOT}")
ENDIF()
#MESSAGE("TRIBITS_PROJECT_ROOT (after env) = '${TRIBITS_PROJECT_ROOT}'")
IF (NOT TRIBITS_PROJECT_ROOT)
  # Fall back on the default convention, in which this file is located at:
  #   <root>/cmake/tribits/ctest.
  GET_FILENAME_COMPONENT(CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
  SET(TRIBITS_PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../..")
ENDIF()
GET_FILENAME_COMPONENT(TRIBITS_PROJECT_ROOT "${TRIBITS_PROJECT_ROOT}" ABSOLUTE)
MESSAGE("TRIBITS_PROJECT_ROOT = '${TRIBITS_PROJECT_ROOT}'")

#
# Read in PROJECT_NAME
#
# Assert that the ProjectName.cmake file exists.
SET(TRIBITS_PROJECT_NAME_INCLUDE "${TRIBITS_PROJECT_ROOT}/ProjectName.cmake")
IF(NOT EXISTS "${TRIBITS_PROJECT_NAME_INCLUDE}")
  MESSAGE(FATAL_ERROR
    "Could not locate ProjectName.cmake.\n"
    "  TRIBITS_PROJECT_ROOT = ${TRIBITS_PROJECT_ROOT}\n"
    "  Set the TRIBITS_PROJECT_ROOT environment variable "
    "to point at the source root.")
ENDIF()
# Include the ProjectName.cmake file and get PROJECT_NAME
INCLUDE(${TRIBITS_PROJECT_NAME_INCLUDE})
IF(NOT PROJECT_NAME)
  MESSAGE(FATAL_ERROR
    "The project name has not been set!"
    "  It should be set in ${TRIBITS_PROJECT_ROOT}/ProjectName.cmake.")
ENDIF()
MESSAGE("PROJECT_NAME = ${PROJECT_NAME}")

#
# Set ${PROJECT_NAME}_TRIBITS_DIR
#
IF (NOT "$ENV{${PROJECT_NAME}_TRIBITS_DIR}" STREQUAL "")
  SET(${PROJECT_NAME}_TRIBITS_DIR "$ENV{${PROJECT_NAME}_TRIBITS_DIR}")
ENDIF()
IF ("${${PROJECT_NAME}_TRIBITS_DIR}" STREQUAL "")
  SET(${PROJECT_NAME}_TRIBITS_DIR "${TRIBITS_PROJECT_ROOT}/cmake/tribits")
ENDIF()
MESSAGE("${PROJECT_NAME}_TRIBITS_DIR = ${${PROJECT_NAME}_TRIBITS_DIR}")

#
# Set default for CTEST_SOURCE_DIRECTORY
#
IF ("${CTEST_SOURCE_DIRECTORY}" STREQUAL "")
  MESSAGE("Set default for CTEST_SOURCE_DIRECTORY to TRIBITS_PROJECT_ROOT='TRIBITS_PROJECT_ROOT='${TRIBITS_PROJECT_ROOT}'")
  SET(CTEST_SOURCE_DIRECTORY ${TRIBITS_PROJECT_ROOT})
ENDIF()

#
# Set default for CTEST_BINARY_DIRECTORY
#
IF ("${CTEST_BINARY_DIRECTORY}" STREQUAL "")
  MESSAGE("Set defualt for CTEST_BINARY_DIRECTORY to $PWD/BUILD='$ENV{PWD}/BUILD'")
  SET(CTEST_BINARY_DIRECTORY $ENV{PWD}/BUILD)
ENDIF()

#
# Set CMAKE_MODULE_PATH
#
SET( CMAKE_MODULE_PATH
  "${TRIBITS_PROJECT_ROOT}"
  "${TRIBITS_PROJECT_ROOT}/cmake"
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch"
  "${${PROJECT_NAME}_TRIBITS_DIR}/ci_support"
  )

INCLUDE(TribitsConstants)
TRIBITS_ASESRT_MINIMUM_CMAKE_VERSION()
INCLUDE(TribitsCMakePolicies)

INCLUDE(PrintVar)
INCLUDE(MultilineSet)
INCLUDE(SetDefaultAndFromEnv)
INCLUDE(AssertDefined)
INCLUDE(AppendSet)
INCLUDE(AppendStringVar)
INCLUDE(TribitsGlobalMacros)
INCLUDE(TribitsStripCommentsFromCMakeCacheFile)

# Need to include the project's version file to get some Git and CDash
# settings specific to the given version
TRIBITS_PROJECT_READ_VERSION_FILE(${TRIBITS_PROJECT_ROOT})

INCLUDE(TribitsFindPythonInterp)
TRIBITS_FIND_PYTHON()
MESSAGE("PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")

#############################
### Do some initial setup ###
#############################


# Get the host type

IF(WIN32)
  SET(HOST_TYPE $ENV{OS})
ELSE()
  FIND_PROGRAM(UNAME_EXE NAMES uname)
  EXECUTE_PROCESS(
    COMMAND ${UNAME_EXE}
    OUTPUT_VARIABLE HOST_TYPE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
ENDIF()


# Find git

FIND_PROGRAM(GIT_EXE NAMES ${GIT_NAME})
MESSAGE("GIT_EXE=${GIT_EXE}")

IF(NOT GIT_EXE)
  QUEUE_ERROR("error: could not find git: GIT_EXE='${GIT_EXE}'")
ENDIF()
IF(NOT EXISTS "${GIT_EXE}")
  QUEUE_ERROR("error: GIT_EXE='${GIT_EXE}' does not exist")
ENDIF()


# Find gitdist

SET(GITDIST_EXE "${${PROJECT_NAME}_TRIBITS_DIR}/python_utils/gitdist")


# Get the host name

SITE_NAME(CTEST_SITE_DEFAULT)


# Get helper functions

INCLUDE(${CMAKE_CURRENT_LIST_DIR}/TribitsCTestDriverCoreHelpers.cmake)


#
# @FUNCTION: TRIBITS_CTEST_DRIVER()
#
# Platform-independent CTest/CDash driver
#
# Usage::
#
#   TRIBITS_CTEST_DRIVER()
#
# This driver code is platform independent.  CTest -S scripts use this
# function to drive the testing process for a TriBITS projects by doing a
# version control (VC) source update on all of the VC repos and then
# configuring, building, testing, and submitting results to CDash for the
# TriBITS packages in a project either in package-by-package mode (the
# default), or in an all-at-once mode (see
# ``${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE``).
#
# **Source and Binary Directory Locations:**
#
# To understand how to set the source and binary directories, one must
# understand that it gets run in one of two different modes:
#
# **Mode 1**: Run where there are already existing source and binary
# directories (i.e. ``CTEST_DASHBOARD_ROOT`` is set empty before call).  In
# this case, ``CTEST_SOURCE_DIRECTORY`` and ``CTEST_BINARY_DIRECTORY`` must be
# set by the user before calling this function.  This is used to test a local
# build and post to CDash (see the custom ``dashboard`` target).
#
# **Mode 2**: A new binary directory is created and new sources are cloned (or
# updated) in a driver directory (i.e. ``CTEST_DASHBOARD_ROOT`` is set before
# call).  In this case, there are always two (partial) project source tree's,
# a) a "driver" skeleton source tree (typically embedded with TriBITS
# directory) that bootstraps the testing process, and b) a true full "source"
# that is (optionally) cloned and/or updated and tested..
#
# There are a few different directory locations are significant for this
# script:
#
#   ``TRIBITS_PROJECT_ROOT``
#
#     The root directory to an existing source tree where the project's
#     `<projectDir>/ProjectName.cmake`_ (defining the ``PROJECT_NAME``
#     variable) and ``Version.cmake`` files can be found.
#
#   ``${PROJECT_NAME}_TRIBITS_DIR``
#
#     The base directory for the TriBITS system's various CMake modules,
#     python scripts, and other files.  By default this is assumed to be in
#     the source tree under ``${TRIBITS_PROJECT_ROOT}`` (see below) but it can
#     be overridden to point to any location.
#
#   ``CTEST_DASHBOARD_ROOT``
#
#     If set, this is the base directory where this script runs that clones
#     the sources for the project.  If this directory does not exist, it will
#     be created.  If provided as the special value "PWD", then the present
#     working directory is used.  If empty, then this var has no effect.
#
#   ``CTEST_SOURCE_DIRECTORY``
#
#     Determines the location of the sources that are used to define packages,
#     dependencies and configure, build, and test the software.  This is a
#     variable that CTest directly reads and must therefore be set. This is
#     used to set `PROJECT_SOURCE_DIR`_ which is used by the TriBITS system.
#     If ``CTEST_DASHBOARD_ROOT`` is set, then this is hard-coded internally
#     to ``${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}``.
#
#   ``CTEST_BINARY_DIRECTORY``
#
#     Determines the location of the binary tree where output from CMake/CTest
#     is put.  This is used to set to `PROJECT_BINARY_DIR`_ which is used by
#     the TriBITS system.  If ``CTEST_DASHBOARD_ROOT`` is set, then this is
#     hard-coded internally to ``${CTEST_DASHBOARD_ROOT}/BUILD``.
#
# **Determining What Packages Get Tested:**
#
# Before any testing is done, the set of packages to be tested is determined
# are determined.  By default, the set of packages to be tested is determined
# by the var:
#
#   ``${PROJECT_NAME}_PACKAGES``
#
#     Determines the specific set of packages to test.  If left at the default
#     value of empty "", then `${PROJECT_NAME}_ENABLE_ALL_PACKAGES`_ is set to
#     ``ON`` and that enables packages as described in
#     `<Project>_ENABLE_ALL_PACKAGES enables all PT (cond. ST) SE packages`_.
#     This variable can also be specified and read from from the env and can
#     use `,` to separate package names instead of ';'.
#
# ToDo: Document other variables that determine what packages get tested.
#
# The other mode is to only test the packages that have changes since the last
# time this build was run and testing packages that previously failed.  That
# mode is turned on by the var:
#
#   ``CTEST_ENABLE_MODIFIED_PACKAGES_ONLY``
#
#     If ``TRUE``, then only packages that have changes pulled from the git
#     repos since the last time the build ran will be tested (in addition to
#     packages that failed in the last build).  If ``FALSE``, the set of
#     packages to be tested is determined by ``${PROJECT_NAME}_PACKAGES`` and
#     other variables as described above.
#
# ToDo: Document other input variables that have defaults, to be set before,
# and can be overridden from the env.
#
# **All-at-once vs. package-by-package mode:**
#
# This function supports driving the configure, build, testing, and submitting
# to CDash of the packages in the TriBITS project either all-at-once or
# package-by-package, based on the value:
#
#   ``${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE``
#
#     If ``TRUE``, then single calls to ``CTEST_CONFIGURE()``,
#     ``CTEST_BUILD()`` and ``CTEST_TEST()` are made for all of the packages
#     to be tested all at once with ``CTEST_SUBMIT()`` called after each of
#     these.  If ``FALSE`` then ``CTEST_CONFIGURE()``, ``CTEST_BUILD()`` and
#     ``CTEST_TEST()` and ``CTEST_SUBMIT()`` are called in a loop, once for
#     each package to be explicitly tested.
#
# Both the all-at-once mode and the package-by-package mode should produce
# equivalent builds of the project and submits to CDash (for correctly
# constructed TriBITS projects and package).  But the package-by-package mode
# will disable failing packages when processing downstream packages, and
# therefore reduce the propagation of failures to downstream packages.
#
# For versions of CMake 3.10.0 and above and newer versions of CDash, the
# all-at-once mode will break down the build and test results on a
# package-by-package basis on CDash.  For older versions of CMake or CDash, it
# will not break down results on a package-by-package basis on CDash and all
# of the build warnings and errors and test will be all globed together.
#
# **Repository Updates:**
#
# Like the rest of TriBITS, ``ctest -S`` scripts written using this function
# supports a collection of extra repositories in addition to the base git
# repository.  The basic clone of the extra repositories requires all repos to
# use the git version control system.
#
# Whether the repos are updated (or left as is) is determined by the var:
#
#  ``CTEST_DO_UPDATES``
#
#    If set to ``TRUE``, then each of the git repos will be cloned if they are
#    missing and if already present will be updated as described below (and
#    will wipe out any local changes).  If ``FALSE``, then the git repos will
#    be left alone (and must therefore already be cloned and updated at the
#    desired state).
#
# **WARNING:** If you don't want local changes in your git repos to get blown
# away, then set ``CTEST_DO_UPDATES`` to ``FALSE``!
#
# CTest itself is used for handling the cloning and the pull of the base
# repository by calling ``CTEST_UPDATE()``.  The other extra git repositories
# are completely handled by the code in this function as described below.
#
# After the base repository is cloned for the first time (by calling
# ``CTEST_UPDATE()``), the extra repositories are cloned using the following
# command::
#
#   git clone <extrarepo_url>
# 
# Therefore, by default, whatever the default branch is set to on clone in the
# base repos, that is the branch that will be used.  Also, future repository
# updates will be done on those branches.
#
# However, if ``${PROJECT_NAME}_BRANCH`` is set to non-empty, then that branch
# will be checked out in all of the repositories.  For the base repository,
# after the clone or update is performed in the call to ``CTEST_UPDATE()``,
# then that branch is checked out using the command::
#
#   $ git checkout -B ${${PROJECT_NAME}_BRANCH} \
#       --track origin/${${PROJECT_NAME}_BRANCH}`
#
# That command is robust and will pass even if the current branch is already a
# tracking branch for ``origin/${${PROJECT_NAME}_BRANCH}``.
#
# If ``${PROJECT_NAME}_BRANCH`` is empty, then each extra repository is
# updated (even after the initial clone) using the commands::
#
#   $ git clean -fdx         # Remove untracked ignored files
#   $ git reset --hard HEAD  # Removed untracked and modified tracked files
#   $ git fetch origin       # Get updated commits
#   $ git reset --hard @{u}  # Deal with forced push
#
# The command ``git clone -fdx`` removes any untracked ignored files that may
# have been created since the last update (either by the build process or by
# someone messing around in that local git repository).  The command ``git
# reset --hard HEAD`` removes any untracked non-ignored files, any modified
# tracked files, and sets ``ORIG_HEAD`` to the current ``HEAD``.  This sets
# ``ORIG_HEAD`` after the initial clone (as ``ORIG_HEAD`` is not set after a
# ``git clone``).  This allows using the range ``ORIG_HEAD..HEAD`` with git
# diff and log commands even after the initial clone.  The ``git fetch``
# command followed by the ``git reset --hard @{u}`` command is used to update
# the local repo to match the remote tracking branch instead of ``git pull` or
# ``git fetch ; git merge @{u}``.  This is done to deal with a possible forced
# push of the remote tracking branch.  Using ``git fetch ; git reset --hard
# @{u}`` ensures that the local branch is exactly the same the remote tracking
# branch no matter what.
#
# If ``${PROJECT_NAME}_BRANCH`` is non-empty, then each extra repository is
# updated (even after the initial clone) using the commands::
#
#   $ git clean -fdx         # Remove untracked ignored files
#   $ git reset --hard HEAD  # Clean files and set ORIG_HEAD to HEAD
#   $ git fetch origin       # Get updated commits
#   $ git checkout -B ${${PROJECT_NAME}_BRANCH} \
#       --track origin/${${PROJECT_NAME}_BRANCH}  # Put on tracking branch
#
# These are the same commands as for the case where ``${PROJECT_NAME}_BRANCH``
# is empty except for the last command which does a checkout of the tracking
# branch ``${PROJECT_NAME}_BRANCH``.  In this case, the ``git reset --hard
# HEAD`` serves an additional purpose.  It sets ``ORIG_HEAD`` to the current
# ``HEAD`` before the update of the branch.  This is important because the
# command ``git checkout -B ...`` does not move ``ORIG_HEAD`` so this reset is
# needed so that the git range ``ORIG_HEAD..HEAD`` gives the changes since the
# last update.  So for an update where no new commits are pulled,
# ``ORIG_HEAD..HEAD`` will return no commits.
#
# Note that the repository updating approach using non-empty
# ``${PROJECT_NAME}_BRANCH`` is more robust, because it can recover from a
# state where someone may have put A repo on a detached head or checked out a
# different branch.  One of these repos might get into this state when a
# person is messing around in the Nightly build and source directories to try
# to figure out what happened and forgot to put the repos back on the correct
# tracking branch.  Therefore, it is recommended to always set
# ``${PROJECT_NAME}_BRANCH`` to a non-null value like ``master`` for git
# repos.
#
# ToDo: Finish Documentation!
#
FUNCTION(TRIBITS_CTEST_DRIVER)

  MESSAGE("")
  MESSAGE("******************************")
  MESSAGE("*** TRIBITS_CTEST_DRIVER() ***")
  MESSAGE("******************************")
  MESSAGE("")

  INITIALIZE_ERROR_QUEUE()

  # The name of the source directory. Defaults to project name, but
  # can be overridden by the environment for cases in which the source
  # directory name does not match the project name.
  SET_DEFAULT_AND_FROM_ENV(CTEST_SOURCE_NAME ${PROJECT_NAME})

  MESSAGE(
    "\n***"
    "\n*** Setting input options to default and reading from env ..."
    "\n***\n")

  SET_DEFAULT_AND_FROM_ENV( CTEST_CONFIGURATION_UNIT_TESTING OFF )

  # The type of test (e.g. Nightly, Experimental, Continuous)
  SET_DEFAULT_AND_FROM_ENV( CTEST_TEST_TYPE Experimental )

  # The default track to send the build to. This can be changed to send
  # the data to a different nightly grouping on the dashboard.
  # If the test type is set to Experimental though the track is forced
  # to "Experimental" this is so that we can have experimental tests
  # on branches.
  IF(${PROJECT_NAME}_TESTING_TRACK)
    SET(${PROJECT_NAME}_TRACK_DEFAULT ${${PROJECT_NAME}_TESTING_TRACK})
  ELSE()
    SET(${PROJECT_NAME}_TRACK_DEFAULT "")
  ENDIF()
  print_var(${PROJECT_NAME}_TRACK_DEFAULT)
  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_TRACK "${${PROJECT_NAME}_TRACK_DEFAULT}")
  IF(CTEST_TEST_TYPE STREQUAL "Experimental" OR CTEST_TEST_TYPE STREQUAL "EXPERIMENTAL")
    SET(${PROJECT_NAME}_TRACK "Experimental")
    MESSAGE("-- Test type is Experimental. Forcing ${PROJECT_NAME}_TRACK to Experimental")
    PRINT_VAR(${PROJECT_NAME}_TRACK)
  ENDIF()

  # The name of the site in the dashboard (almost never need to override this)
  SET_DEFAULT_AND_FROM_ENV( CTEST_SITE ${CTEST_SITE_DEFAULT} )

  # The root of the dasbhoard where ${PROJECT_NAME} will be cloned and the
  # BUILD directory will be create (only override for separate testing)
  SET_DEFAULT_AND_FROM_ENV( CTEST_DASHBOARD_ROOT "" )
  IF (CTEST_DASHBOARD_ROOT STREQUAL "PWD")
    SET(CTEST_DASHBOARD_ROOT ${CMAKE_CURRENT_BINARY_DIR})
    PRINT_VAR(CTEST_DASHBOARD_ROOT)
  ENDIF()

  # The build type (e.g. DEBUG, RELEASE, NONE)
  SET_DEFAULT_AND_FROM_ENV( BUILD_TYPE NONE )

  # Verobse configure or now
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_VERBOSE_CONFIGURE OFF )

  # Set the default compiler version
  SET_DEFAULT_AND_FROM_ENV(COMPILER_VERSION UNKNOWN)

  # The name of the build that appears in the dashboard
  SET_DEFAULT_AND_FROM_ENV( CTEST_BUILD_NAME
    "${HOST_TYPE}-${COMPILER_VERSION}-${BUILD_DIR_NAME}" )

  # Remove the entire build directory if it exists or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE )

  # Remove an existing CMakeCache.txt file or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_WIPE_CACHE TRUE )

  # Select a default generator.
  SELECT_DEFAULT_GENERATOR()
  SET_DEFAULT_AND_FROM_ENV( CTEST_CMAKE_GENERATOR ${DEFAULT_GENERATOR})

  # Do the Git updates or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_UPDATES TRUE )

  # Generate the XML dependency output files or not in the inner CMake
  # configure.  There is really no reason to do this.  This option is
  # maintained for backward compatibility.
  SET_DEFAULT_AND_FROM_ENV( CTEST_GENERATE_DEPS_XML_OUTPUT_FILE FALSE )

  # Flags used on git when doing a Git update
  SET_DEFAULT_AND_FROM_ENV( CTEST_UPDATE_ARGS "")

  # Flags used on update when doing a Git update
  SET_DEFAULT_AND_FROM_ENV( CTEST_UPDATE_OPTIONS "")

  # If doing all-at-one approach, use new CMkae/CTest/CDash features to allow
  # it to split out results into different rows on CDash like the
  # package-by-packages appraoch.
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_CTEST_USE_NEW_AAO_FEATURES  FALSE )
 
  # Do all-at-once configure, build, test and submit (or package-by-package)
  IF (${PROJECT_NAME}_CTEST_USE_NEW_AAO_FEATURES)
    SET(${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE_DEFAULT TRUE)
  ELSE()
    SET(${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE_DEFAULT FALSE)
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE
    ${${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE_DEFAULT} )

  # Flags passed to 'make' assume gnumake with unix makefiles
  IF("${CTEST_CMAKE_GENERATOR}" MATCHES "Unix Makefiles")
    SET_DEFAULT_AND_FROM_ENV( CTEST_BUILD_FLAGS "-j2")
  ELSE()
    SET_DEFAULT_AND_FROM_ENV( CTEST_BUILD_FLAGS "")
  ENDIF()

  # Do the build or use an existing build
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_BUILD TRUE )

  # Do the tests or not (Note: must be true for coverage testing)
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_TEST TRUE )

  # Maximum number of procs an mpi test can request (if more are requested,
  # the test will be skipped).  Value of 0 means no override (determined
  # internally).
  SET_DEFAULT_AND_FROM_ENV( MPI_EXEC_MAX_NUMPROCS 0 )

  # How many tests ctest will spawn simultaneously
  SET_DEFAULT_AND_FROM_ENV( CTEST_PARALLEL_LEVEL 1 )

  # Turn off or change warnings-as-errors flag(s) (i.e. -Werror)
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS "" )

  # Do coverage testing or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_COVERAGE_TESTING FALSE )

  # Command to run to get coverage results
  SET_DEFAULT_AND_FROM_ENV( CTEST_COVERAGE_COMMAND gcov )

  # Do memory testing (i.e. valgrind) or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_MEMORY_TESTING FALSE )

  # Command used to perform the memory testing (i.e. valgrind)
  FIND_PROGRAM(VALGRIND_EXE NAMES valgrind)
  PRINT_VAR(VALGRIND_EXE)
  IF (VALGRIND_EXE)
    SET(CTEST_MEMORYCHECK_COMMAND_DEFAULT "${VALGRIND_EXE}")
  ELSE()
    SET(CTEST_MEMORYCHECK_COMMAND_DEFAULT)
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV( CTEST_MEMORYCHECK_COMMAND "${CTEST_MEMORYCHECK_COMMAND_DEFAULT}" )

  # Set the default options
  SET_DEFAULT_AND_FROM_ENV( CTEST_MEMORYCHECK_COMMAND_OPTIONS "")

  # Generate the basic package dependencies XML file in the outer CTest
  # program.  This XML file is used to match up modified files with with
  # changed TriBITS packages.  This file only needs to be generated in CI
  # iterations and is not needed in Nightly testing.  Turning off its
  # generation can also speed up local manual testing for large projects with
  # lots of TriBITS packges.
  SET_DEFAULT_AND_FROM_ENV( CTEST_GENERATE_OUTER_DEPS_XML_OUTPUT_FILE TRUE )

  # Generate and submit the CDash subprojects XML file
  SET_DEFAULT_AND_FROM_ENV( CTEST_SUBMIT_CDASH_SUBPROJECTS_DEPS_FILE TRUE )

  # Submit the results to the dashboard or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_SUBMIT TRUE )

  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE OFF )

  # List of additional packges that will be enabled over the current set
  # of all packagess (that would be set by ${PROJECT_NAME}_ENABLE_ALL_PACKAGES).
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_ADDITIONAL_PACKAGES "" )

  # List of packages to not directly process.  NOTE: Listing these packages
  # here will *not* disable the package in the CMake build system.  To do
  # that, you will have to disable them in the variable
  # EXTRA_CONFIGURE_OPTIONS (set in your driver script.
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_EXCLUDE_PACKAGES "" )

  IF(${PROJECT_NAME}_REPOSITORY_BRANCH)
    SET(${PROJECT_NAME}_BRANCH_DEFAULT ${${PROJECT_NAME}_REPOSITORY_BRANCH})
  ELSE()
    SET(${PROJECT_NAME}_BRANCH_DEFAULT "")
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_BRANCH "${${PROJECT_NAME}_BRANCH_DEFAULT}" )

  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_EXTRAREPOS_BRANCH "${${PROJECT_NAME}_BRANCH}" )

  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE "${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT}" )

  IF(CTEST_TEST_TYPE STREQUAL "Nightly")
    SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_REPOSITORY_LOCATION
       "${${PROJECT_NAME}_REPOSITORY_LOCATION_NIGHTLY_DEFAULT}")
  ELSE()
    SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_REPOSITORY_LOCATION
       "${${PROJECT_NAME}_REPOSITORY_LOCATION_DEFAULT}")
  ENDIF()

  # Select the ${PROJECT_NAME} packages to enable (empty means to select all
  # available).  This will override any disabled packages but not those
  # disabled by ${PROJECT_NAME}_EXCLUDE_PACKAGES.
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_PACKAGES "" )
  SET(${PROJECT_NAME}_PACKAGES_USER_SELECTED ${${PROJECT_NAME}_PACKAGES})
  SPLIT("${${PROJECT_NAME}_PACKAGES_USER_SELECTED}" ","
    ${PROJECT_NAME}_PACKAGES_USER_SELECTED)
  SET(${PROJECT_NAME}_PACKAGES "")
  # Note: above, we have to keep the name ${PROJECT_NAME}_PACKAGES to maintain
  # backward compatibility of this CTest script but we want to let
  # ${PROJECT_NAME}_PACKAGES always be the full set of packages as defined by
  # the basic readin process.

  # Set the file that the extra repos will be read from
  #
  # NOTE: Here, we have no choice but to point into the "driver"
  # ${PROJECT_NAME} source treee because the local ${PROJECT_NAME} sources have not
  # even been checked out yet!  Unless, of course, we are unit testing
  # in which case we will use whatever has been passed in.

  IF (${PROJECT_NAME}_SKIP_EXTRAREPOS_FILE)
    SET(${PROJECT_NAME}_EXTRAREPOS_FILE_DEFAULT)
  ELSE()
    SET(${PROJECT_NAME}_EXTRAREPOS_FILE_DEFAULT
      "${TRIBITS_PROJECT_ROOT}/cmake/${${PROJECT_NAME}_EXTRA_EXTERNAL_REPOS_FILE_NAME}")
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_EXTRAREPOS_FILE
    "${${PROJECT_NAME}_EXTRAREPOS_FILE_DEFAULT}")

  # Select the set of extra external repos to add in packages.
  # These are the same types as CTEST_TEST_TYPE (e.g. 'Continuous' and
  # 'Nightly').  This is set by default to ${CTEST_TEST_TYPE} can be
  # overridden independent of ${CTEST_TEST_TYPE} also.
  #
  # If in release mode generally we do not want any external repositories
  # even though the CTEST_TEST_TYPE is set to "Nightly" for most release
  # builds.
  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
  IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
    SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_DEFAULT ${CTEST_TEST_TYPE})
  ELSE()
    SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_DEFAULT "None")
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE
     "${${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_DEFAULT}" )

  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_PRE_REPOSITORIES "")
  SPLIT("${${PROJECT_NAME}_PRE_REPOSITORIES}"  "," ${PROJECT_NAME}_PRE_REPOSITORIES)

  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_EXTRA_REPOSITORIES "")
  SPLIT("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  "," ${PROJECT_NAME}_EXTRA_REPOSITORIES)

  # Set as part of CI testing in order to only enable modified packages
  SET_DEFAULT_AND_FROM_ENV( CTEST_ENABLE_MODIFIED_PACKAGES_ONLY OFF )

  # Set if implicitly enabled packages should be explicitly processed in
  # package-by-package mode.
  IF (CTEST_ENABLE_MODIFIED_PACKAGES_ONLY AND NOT CTEST_START_WITH_EMPTY_BINARY_DIRECTORY)
    SET( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES_DEFAULT FALSE )
  ELSE()
    SET( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES_DEFAULT TRUE )
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES
    ${CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES_DEFAULT})

  # Set if we should disable enabled fwd packages based on disabled required deps.
  # To make testing robust, we need to do this.
  IF ("${${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT ON)
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES
    ${${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT})

  # Set second drop site and location
  SET_DEFAULT_AND_FROM_ENV( TRIBITS_2ND_CTEST_DROP_SITE "" )
  SET_DEFAULT_AND_FROM_ENV( TRIBITS_2ND_CTEST_DROP_LOCATION "" )

  MESSAGE(
    "\n***"
    "\n*** Setting unit testing input options to default and reading from env ..."
    "\n***\n")

  SET_DEFAULT_AND_FROM_ENV( CTEST_DEPENDENCY_HANDLING_UNIT_TESTING FALSE )

  SET_DEFAULT_AND_FROM_ENV( CTEST_UPDATE_RETURN_VAL 0 )

  IF (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    SET(GIT_EXE /somebasedir/git)
  ENDIF()


  MESSAGE(
    "\n***"
    "\n*** Misc setup ..."
    "\n***\n")

  #
  # Setup and create the base dashboard directory if it is not created yet.
  #

  # NOTE: This is only used in general testing dashbaoard mode, not in local
  # experimental testing mode.

  IF (CTEST_DASHBOARD_ROOT)
    SET( CTEST_BINARY_NAME BUILD )
    SET( CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
    SET( CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")
    IF (NOT EXISTS "${CTEST_DASHBOARD_ROOT}")
      MESSAGE("Creating the dashboard root directory \"${CTEST_DASHBOARD_ROOT}\" ...")
      FILE(MAKE_DIRECTORY "${CTEST_DASHBOARD_ROOT}")
    ENDIF()
  ENDIF()
  PRINT_VAR(CTEST_SOURCE_DIRECTORY)
  PRINT_VAR(CTEST_BINARY_DIRECTORY)

  SET(PROJECT_SOURCE_DIR "${CTEST_SOURCE_DIRECTORY}")
  SET(PROJECT_BINARY_DIR "${CTEST_BINARY_DIRECTORY}")
  PRINT_VAR(PROJECT_SOURCE_DIR)
  PRINT_VAR(PROJECT_BINARY_DIR)

  # Set override hook for unit testing
  SET_DEFAULT_AND_FROM_ENV( ${PROJECT_NAME}_SOURCE_DIRECTORY ${CTEST_SOURCE_DIRECTORY} )

  # Must be set here after CTEST_BINARY_DIRECTORY is set!
  SET(FAILED_PACKAGES_FILE_NAME "${CTEST_BINARY_DIRECTORY}/failedPackages.txt")

  #
  # Some platform-independent setup
  #

  INCLUDE("${TRIBITS_PROJECT_ROOT}/CTestConfig.cmake")
  SET(CTEST_USE_LAUNCHERS 1)

  # For coverage dashboards, send results to specialized dashboard if
  # requested
  IF (CTEST_DO_COVERAGE_TESTING)
    # Allow override of CDash drop site but use standard by default
    SET_DEFAULT(CTEST_DROP_SITE_COVERAGE_DEFAULT
      ${CTEST_DROP_SITE})
    SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_SITE_COVERAGE
      "${CTEST_DROP_SITE_COVERAGE_DEFAULT}")
    SET(CTEST_DROP_SITE "${CTEST_DROP_SITE_COVERAGE}" )
    # Allow override of CDash drop location but use standard by default
    SET_DEFAULT(CTEST_DROP_LOCATION_COVERAGE_DEFAULT
      ${CTEST_DROP_LOCATION})
    SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_LOCATION_COVERAGE
      "${CTEST_DROP_LOCATION_COVERAGE_DEFAULT}")
    SET(CTEST_DROP_LOCATION "${CTEST_DROP_LOCATION_COVERAGE}" )

    # NOTE: You must set these down below the include of
    # CTestConfig.cmake so that CTEST_DROP_SITE and CTEST_DROP_LOCATION read
    # from that file will set the defaults for the coverage options.

  ENDIF()

  #
  # Setup for the VC update
  #

  SET(CREATE_VC_UPDATE_FILE FALSE)

  IF (CTEST_DO_UPDATES)

    SET(UPDATE_TYPE "git")
    MESSAGE("UPDATE_TYPE = '${UPDATE_TYPE}'")

    SET(CTEST_UPDATE_COMMAND "${GIT_EXE}")
    MESSAGE("CTEST_UPDATE_COMMAND='${CTEST_UPDATE_COMMAND}'")

    IF(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
      MESSAGE("${CTEST_SOURCE_DIRECTORY} does not exist so setting up for an initial checkout")
      SET( CTEST_CHECKOUT_COMMAND
        "\"${GIT_EXE}\" clone ${CTEST_UPDATE_ARGS} ${${PROJECT_NAME}_REPOSITORY_LOCATION}" )
      MESSAGE("CTEST_CHECKOUT_COMMAND='${CTEST_CHECKOUT_COMMAND}'")
    ELSE()
      MESSAGE("${CTEST_SOURCE_DIRECTORY} exists so skipping the initial checkout.")
      SET(CREATE_VC_UPDATE_FILE TRUE)
    ENDIF()

  ENDIF()

  #
  # Write a few variables to the global level to make cmake happy.
  #
  # If these don't get set in the base CTest script scope, CTest returns an
  # error!
  #

  SET(CTEST_SOURCE_DIRECTORY ${CTEST_SOURCE_DIRECTORY} CACHE INTERNAL "")
  SET(CTEST_BINARY_DIRECTORY ${CTEST_BINARY_DIRECTORY} CACHE INTERNAL "")
  IF ("${CTEST_COMMAND}" STREQUAL "")
    SET(CTEST_COMMAND ctest)
  ENDIF()
  SET(CTEST_COMMAND ${CTEST_COMMAND} CACHE INTERNAL "")

  #
  # Empty out the binary directory
  #

  IF (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY)
    MESSAGE("\nCleaning out binary directory '${CTEST_BINARY_DIRECTORY}' ...")
    CTEST_EMPTY_BINARY_DIRECTORY("${CTEST_BINARY_DIRECTORY}")
  ENDIF()
  # NOTE: The above command will *not* delete the build directory unless there
  # is a CMakeLists.txt file in this directory.  I think Kitware put in this
  # check to avoid accidentally deleting the wrong directory by accident.
  # Also note that you have to delete the build directory before any commands
  # are run that would write files to them (and many of the steps in this
  # process do write files to the binary directory other than just CMake).

  MESSAGE(
    "\n***"
    "\n*** Read in the set of extra repos ..."
    "\n***\n")

  TRIBITS_SETUP_EXTRAREPOS()

  # NOTE: You have to set up the set of extra repos before you can read the
  # Dependencies.cmake files since the extra repos must be cloned first.


  MESSAGE(
    "\n***"
    "\n*** Start up a new dashboard ..."
    "\n***\n")

  PRINT_VAR(CTEST_TEST_TYPE)
  PRINT_VAR(${PROJECT_NAME}_TRACK)

  IF(${PROJECT_NAME}_TRACK)
    CTEST_START(${CTEST_TEST_TYPE} TRACK ${${PROJECT_NAME}_TRACK})
  ELSE()
    CTEST_START(${CTEST_TEST_TYPE})
  ENDIF()
  # NOTE: If the soruce directory does not yet exist, then CTEST_START() will
  # clone it!


  MESSAGE(
    "\n***"
    "\n*** Update the source code repositories ..."
    "\n***\n")

  SET(UPDATE_FAILED FALSE)

  IF (CTEST_DO_UPDATES)

    IF(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
      QUEUE_ERROR("error: source directory does not exist just prior to CTEST_UPDATE call -- initial checkout did not work")
      REPORT_QUEUED_ERRORS()
      RETURN()
    ENDIF()

    MESSAGE("\nCalling CTEST_UPDATE() to update base source repo '${CTEST_SOURCE_DIRECTORY}' ...")
    CTEST_UPDATE_WRAPPER( SOURCE "${CTEST_SOURCE_DIRECTORY}"
      RETURN_VALUE  CTEST_UPDATE_RETURN_VAL)
    MESSAGE("CTEST_UPDATE(...) returned '${CTEST_UPDATE_RETURN_VAL}'")

    TRIBITS_CLONE_OR_UPDATE_ALL_REPOS(${CTEST_UPDATE_RETURN_VAL}  LOC_UPDATE_FAILED)
    IF (LOC_UPDATE_FAILED)
      SET(UPDATE_FAILED TRUE)
    ENDIF()

    IF (CREATE_VC_UPDATE_FILE)
      TRIBITS_CREATE_REPO_UPDATES_FILE()
      # NOTE: We can only create the Updates.txt file using `gitdist
      # ... ORIG_HEAD..HEAD` when doing an update and not after the initial
      # clone.  That is because ORIG_HEAD will not exist for the base git repo
      # after the initial clone.
    ENDIF()

  ELSE()

     MESSAGE("Skipping the update by request!")

  ENDIF()


  MESSAGE(
    "\n***"
    "\n*** Read in the set of packages and their dependencies ..."
    "\n***\n")

  # NOTE: You must read the Dependencies.cmake files *after* you have cloned
  # (or updated) all of the code!

  TRIBITS_SETUP_PACKAGES()

  SET(CDASH_SUBPROJECT_XML_FILE
    "${CTEST_BINARY_DIRECTORY}/${${PROJECT_NAME}_CDASH_SUBPROJECT_DEPS_XML_FILE_NAME}")
  PRINT_VAR(CDASH_SUBPROJECT_XML_FILE)

  MESSAGE(
    "\n***"
    "\n*** Disabling packages based on what was set in ${PROJECT_NAME}_EXCLUDE_PACKAGES ..."
    "\n***\n")

  DISABLE_EXCLUDED_PACKAGES()

  IF (NOT CTEST_ENABLE_MODIFIED_PACKAGES_ONLY)
    MESSAGE(
      "\n***"
      "\n*** Determining what packages to enable based what was set in ${PROJECT_NAME}_PACKAGES by the user ..."
      "\n***\n")
    ENABLE_USER_SELECTED_PACKAGES()
  ELSE()
    MESSAGE(
      "\n***"
      "\n*** Determining what packages to enable based on what changed (and failed last CI iteration) ..."
      "\n***\n")
    ENABLE_ONLY_MODIFIED_PACKAGES()
    SET(${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES ON)
  ENDIF()


  MESSAGE(
    "\n***"
    "\n*** Adjust the package dependencies to enable upstream and"
    " (optionally) downstream packages ..."
    "\n***"
    )

  SET(${PROJECT_NAME}_ENABLE_TESTS ON)
  SET(${PROJECT_NAME}_ENABLE_EXAMPLES ON)
  SET(${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES ON)
  SET(DO_PROCESS_MPI_ENABLES FALSE) # Should not be needed but CMake is messing up
  TRIBITS_ADJUST_AND_PRINT_PACKAGE_DEPENDENCIES()
  # Above sets ${PROJECT_NAME}_NUM_ENABLED_PACKAGES

  SELECT_FINAL_SET_OF_PACKAGES_TO_DIRECTLY_TEST()
  # Above sets ${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST

  TRIBITS_PRINT_ENABLED_PACKAGES_LIST_FROM_VAR( ${PROJECT_NAME}_PACKAGES_TO_DIRECTLY_TEST
    "\nFinal set of packages to be explicitly processed by CTest/CDash" ON FALSE)

  MESSAGE(
    "\n***"
    "\n*** Determine if to go ahead with configure, build, test ..."
    "\n***")

  IF (CTEST_ENABLE_MODIFIED_PACKAGES_ONLY
    AND ${PROJECT_NAME}_NUM_ENABLED_PACKAGES GREATER 0
    AND MODIFIED_PACKAGES_LIST
    )
    MESSAGE("\nMODIFIED_PACKAGES_LIST='${MODIFIED_PACKAGES_LIST}'"
      ":  Found modified packages, processing enabled packages!\n")
  ELSE()
    MESSAGE(
      "\nCTEST_ENABLE_MODIFIED_PACKAGES_ONLY=${CTEST_ENABLE_MODIFIED_PACKAGES_ONLY}"
      "  Running in regular mode, processing all enabled packages!\n")
  ENDIF()

  IF (${PROJECT_NAME}_NUM_ENABLED_PACKAGES GREATER 0)
    MESSAGE(
      "\n${PROJECT_NAME}_NUM_ENABLED_PACKAGES=${${PROJECT_NAME}_NUM_ENABLED_PACKAGES}:"
      "  Configuring packages!\n")
  ELSE()
    MESSAGE(
      "\n${PROJECT_NAME}_NUM_ENABLED_PACKAGES=${${PROJECT_NAME}_NUM_ENABLED_PACKAGES}:"
      "  Exiting the script!\n")
    REPORT_QUEUED_ERRORS()
    RETURN()
  ENDIF()

  #
  # Delete the CMakeCache.txt file and the CMakeFiles directory for a clean
  # reconfigure.
  #

  IF (CTEST_WIPE_CACHE)
    SET(CACHE_FILE_NAME "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
    IF (EXISTS "${CACHE_FILE_NAME}")
      MESSAGE("Removing existing cache file '${CACHE_FILE_NAME}' ...")
      FILE(REMOVE "${CACHE_FILE_NAME}")
    ENDIF()
    SET(CMAKE_FILES_DIR "${CTEST_BINARY_DIRECTORY}/CMakeFiles/")
    IF (EXISTS "${CMAKE_FILES_DIR}")
      MESSAGE("Removing existing '${CMAKE_FILES_DIR}' ...")
      FILE(REMOVE_RECURSE "${CMAKE_FILES_DIR}")
    ENDIF()
  ENDIF()
  # NOTE: Above, we have to delete the CMakeCache.txt file only after we are
  # sure we are going to be configuring packages.  There must be a
  # CMakeCache.txt file present in the binary directory or the
  # CTEST_EMPTY_BINARY_DIRECTORY() command will *not* actually delete the
  # build directory!  Also, with updated versions of CMake (2.8.10 and above)
  # you have to delete the CMakeFiles directory in addition to the
  # CMakeCache.txt file or it will not configure correctly (due to Fortran/C
  # linkage tests for one).

  MESSAGE(
    "\n***"
    "\n*** Uploading update, notes, and the subproject dependencies XML files ..."
    "\n***\n"
    )

  IF (EXISTS ${CTEST_BINARY_DIRECTORY}/Updates.txt)
    SET(CTEST_NOTES_FILES_WO_CACHE
      "${CTEST_BINARY_DIRECTORY}/Updates.txt;${CTEST_NOTES_FILES}")
  ELSE()
    SET(CTEST_NOTES_FILES_WO_CACHE "${CTEST_NOTES_FILES}")
  ENDIF()
  PRINT_VAR(CTEST_NOTES_FILES_WO_CACHE)

  # Note: We must only do the submit after we have decided if there are any
  # packages to enable or not and otherwise exit the script!

  IF (UPDATE_FAILED)
    MESSAGE("The VC update failed so submitting update and stopping ...")
    IF (CTEST_DO_SUBMIT)
      TRIBITS_CTEST_SUBMIT( PARTS update notes )
    ENDIF()
    REPORT_QUEUED_ERRORS()
    RETURN()
  ENDIF()

  IF (CTEST_DO_SUBMIT AND EXISTS ${CDASH_SUBPROJECT_XML_FILE})
    TRIBITS_CTEST_SUBMIT( FILES ${CDASH_SUBPROJECT_XML_FILE})
    MESSAGE("\nSubmitted subproject dependencies XML file!")
  ELSE()
    MESSAGE("\nSkipping submitted subproject dependencies XML file on request!")
  ENDIF()

  MESSAGE(
    "\n***"
    "\n*** Configure, build, test, and submit results for ${PROJECT_NAME} packages:"
    "\n***")

  SET(CMAKE_CACHE_CLEAN_FILE "${CTEST_BINARY_DIRECTORY}/CMakeCache.clean.txt")
  SET(${PROJECT_NAME}_FAILED_PACKAGES)

  IF (${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE)

    TRIBITS_CTEST_ALL_AT_ONCE()

  ELSE()

    TRIBITS_CTEST_PACKAGE_BY_PACKAGE()

  ENDIF()

  IF(${PROJECT_NAME}_FAILED_PACKAGES)
    MESSAGE(
      "\nFinal set packages that had any failures: '${${PROJECT_NAME}_FAILED_PACKAGES}'")
  ENDIF()

  # Write a file listing the packages that failed.  This will be read in on the next CI
  # iteration since these packages must be enabled
  FILE(WRITE "${FAILED_PACKAGES_FILE_NAME}" "${${PROJECT_NAME}_FAILED_PACKAGES}\n")

  REPORT_QUEUED_ERRORS()

  STRING(REPLACE "submit.php"
       "index.php" CDASH_PROJECT_LOCATION
       "${CTEST_DROP_LOCATION}")
  MULTILINE_SET( SEE_CDASH_LINK_STR
    "\nSee results for:\n"
    "  Site: ${CTEST_SITE}\n"
    "  Build Name: ${CTEST_BUILD_NAME}\n"
    "at:\n"
    "  http://${CTEST_DROP_SITE}${CDASH_PROJECT_LOCATION}&display=project\n"
    )
  # ToDo: Above: We would love to provide the buildID to point to the exact
  # URL but CTest does not currently give that to you.

  IF ("${${PROJECT_NAME}_FAILED_PACKAGES}" STREQUAL "")
    MESSAGE(
      "${SEE_CDASH_LINK_STR}\n"
      "TRIBITS_CTEST_DRIVER: OVERALL: ALL PASSSED\n")
  ELSE()
    # ToDo: Find out why other breaking tests don't fail when FATAL_ERROR is
    # removed!
    MESSAGE(FATAL_ERROR
      "${SEE_CDASH_LINK_STR}\n"
      "TRIBITS_CTEST_DRIVER: OVERALL: ALL FAILED\n")
    # NOTE: FATAL_ERROR is needed so that the ctest -S script returns != 0
    # Also, it is critical to dislplay the "See results" in this
    # MESSAGE(FATAL_ERROR ...) command in order for it to be printed last.
    # Otherwise, if you run with ctest -V -S, then the ouptut from
    # CTEST_TEST() will be printed last :-(
  ENDIF()

ENDFUNCTION()
